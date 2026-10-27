"""
Script Name: trends.py
Purpose:     Statistics for RSV composition trends — separating real within-lineage evolution
             from lineage replacement (turnover), estimating lineage growth rates, and
             running robustness checks.
Inputs:      A features DataFrame (one row per sequence) with columns:
             Subtype, Lineage, year, region, and composition columns (gc_content, oe_CG, ...).
Outputs:     tidy pandas DataFrames of model results; no file I/O here.
Notes:       The key idea is the within/between decomposition. A rising overall trend can come
             from (a) genomes changing inside each lineage, or (b) high-value lineages replacing
             low-value ones. Adding lineage fixed effects absorbs (b), so the remaining year
             coefficient is the within-lineage slope.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf


def decompose_trend(df, feature, group="Lineage", year="year", covariates=()):
    """
    Split a composition trend into total vs within-group slope.

    total   : slope of feature ~ year across all sequences (what we saw before)
    within  : slope of feature ~ year + C(group), i.e. holding lineage constant
    between : total - within, the part explained by lineage turnover

    `covariates` (e.g. ("region",)) are added as categorical controls to both models,
    which is how we adjust for sampling bias.
    Returns a one-row DataFrame.
    """
    d = df[[feature, year, group] + list(covariates)].dropna().copy()
    ctrl = "".join(f" + C({c})" for c in covariates)
    tot = smf.ols(f"{feature} ~ {year}{ctrl}", data=d).fit()
    win = smf.ols(f"{feature} ~ {year} + C({group}){ctrl}", data=d).fit()
    return pd.DataFrame([{
        "feature": feature, "n": len(d), "n_groups": d[group].nunique(),
        "slope_total": tot.params[year], "se_total": tot.bse[year], "p_total": tot.pvalues[year],
        "slope_within": win.params[year], "se_within": win.bse[year], "p_within": win.pvalues[year],
        "slope_between": tot.params[year] - win.params[year],
        "pct_from_turnover": 100 * (1 - win.params[year] / tot.params[year]) if tot.params[year] != 0 else np.nan,
    }])


def per_lineage_slopes(df, feature, year="year", min_n=100, min_years=4):
    """
    Fit feature ~ year separately inside each lineage. Shows whether individual lineages
    are evolving, and in which direction. Lineages with too few sequences/years are skipped.
    """
    out = []
    for (sub, lin), g in df.groupby(["Subtype", "Lineage"]):
        g = g[[feature, year]].dropna()
        if len(g) < min_n or g[year].nunique() < min_years:
            continue
        m = smf.ols(f"{feature} ~ {year}", data=g).fit()
        out.append({"Subtype": sub, "Lineage": lin, "n": len(g),
                    "years": g[year].nunique(),
                    "slope": m.params[year], "se": m.bse[year],
                    "p": m.pvalues[year], "mean": g[feature].mean()})
    res = pd.DataFrame(out)
    if len(res):
        res["sig"] = res["p"] < 0.05
        res = res.sort_values(["Subtype", "slope"])
    return res


def lineage_growth_rate(df, year="year", min_n=50, min_years=3):
    """
    Relative growth rate of each lineage: logistic regression of
    "is this sequence lineage L" against time, fitted within its own subtype.
    The year coefficient is the logit-per-year growth rate — positive = gaining share.
    This is the continuous fitness proxy replacing the binary expanding/fading split.
    """
    out = []
    for sub, s in df.groupby("Subtype"):
        s = s[[year, "Lineage"]].dropna()
        for lin, g in s.groupby("Lineage"):
            if len(g) < min_n or g[year].nunique() < min_years:
                continue
            y = (s["Lineage"] == lin).astype(int)
            X = sm.add_constant(s[[year]].astype(float))
            try:
                m = sm.Logit(y, X).fit(disp=0)
            except Exception:
                continue
            out.append({"Subtype": sub, "Lineage": lin, "n": len(g),
                        "growth_rate": m.params[year], "se": m.bse[year],
                        "p": m.pvalues[year]})
    return pd.DataFrame(out)


def growth_vs_composition(df, growth, features):
    """
    Join each lineage's growth rate to its mean composition, then correlate.
    Answers: do faster-growing lineages have systematically different composition?
    Returns (per-lineage table, correlation table).
    """
    from scipy.stats import spearmanr
    comp = df.groupby(["Subtype", "Lineage"])[list(features)].mean().reset_index()
    tab = growth.merge(comp, on=["Subtype", "Lineage"])
    rows = []
    for sub, g in tab.groupby("Subtype"):
        for f in features:
            if g[f].notna().sum() >= 5:
                r, p = spearmanr(g["growth_rate"], g[f])
                rows.append({"Subtype": sub, "feature": f, "n_lineages": len(g),
                             "rho": r, "p": p})
    return tab, pd.DataFrame(rows)


def robustness(df, feature, group="Lineage", year="year"):
    """
    Re-run the within/between decomposition under several sensitivity scenarios
    (step 5): full data, excluding the COVID bottleneck, region-balanced, and
    region-adjusted. Returns one row per scenario per subtype.
    """
    out = []
    for sub, s in df.groupby("Subtype"):
        scenarios = {
            "full": s,
            "excl_2020_2021": s[~s[year].isin([2020, 2021])],
            "region_adjusted": s,          # same data, region added as covariate below
            "region_balanced": _balance(s, "region", year),
        }
        for name, d in scenarios.items():
            if len(d) < 200:
                continue
            cov = ("region",) if name == "region_adjusted" else ()
            try:
                r = decompose_trend(d, feature, group=group, year=year, covariates=cov)
            except Exception:
                continue
            r.insert(0, "scenario", name)
            r.insert(0, "Subtype", sub)
            out.append(r)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def _balance(d, col, year, seed=0):
    """Downsample so no region dominates within a year (caps each region-year at the median)."""
    rng = np.random.default_rng(seed)
    keep = []
    for y, g in d.groupby(year):
        counts = g[col].value_counts()
        cap = int(max(1, counts.median()))
        for r, gg in g.groupby(col):
            idx = gg.index.values
            if len(idx) > cap:
                idx = rng.choice(idx, cap, replace=False)
            keep.extend(idx)
    return d.loc[sorted(keep)]


def founder_composition(df, features, n_founder=30, min_n=50, min_years=3, year="year"):
    """
    Tree-free analogue of the paper's "earliest viruses" analysis.

    For each lineage, take its earliest `n_founder` sequences and average their composition.
    This is the lineage's founding genome signature. Returned alongside how successful the
    lineage later became (peak share, span in years, total count), so we can ask whether
    founding composition predicts eventual success.
    """
    out = []
    for (sub, lin), g in df.groupby(["Subtype", "Lineage"]):
        g = g.dropna(subset=[year]).sort_values("date" if "date" in g else year)
        if len(g) < min_n or g[year].nunique() < min_years:
            continue
        head = g.head(n_founder)
        rec = {"Subtype": sub, "Lineage": lin, "n_total": len(g),
               "first_year": g[year].min(), "last_year": g[year].max(),
               "span_years": g[year].max() - g[year].min(), "n_founder": len(head)}
        for f in features:
            rec[f"founder_{f}"] = head[f].mean()
            rec[f"overall_{f}"] = g[f].mean()
        out.append(rec)
    res = pd.DataFrame(out)
    if len(res):
        # peak share: the largest fraction of its subtype-year this lineage ever reached
        tot = df.groupby(["Subtype", year]).size().rename("yr_total")
        cnt = df.groupby(["Subtype", "Lineage", year]).size().rename("n").reset_index()
        cnt = cnt.join(tot, on=["Subtype", year])
        cnt["share"] = cnt["n"] / cnt["yr_total"]
        peak = cnt.groupby(["Subtype", "Lineage"])["share"].max().rename("peak_share")
        res = res.merge(peak.reset_index(), on=["Subtype", "Lineage"], how="left")
    return res


def founder_vs_success(founder, growth, features):
    """
    Correlate each lineage's founding composition with its later success
    (growth rate, peak share, span). Mirrors the paper's claim that persistent
    lineages already differed in composition at their origin.
    """
    from scipy.stats import spearmanr
    t = founder.merge(growth[["Subtype", "Lineage", "growth_rate"]],
                      on=["Subtype", "Lineage"], how="inner")
    rows = []
    for sub, g in t.groupby("Subtype"):
        for f in features:
            for outcome in ["growth_rate", "peak_share", "span_years"]:
                x, y = g[f"founder_{f}"], g[outcome]
                ok = x.notna() & y.notna()
                if ok.sum() >= 6:
                    r, p = spearmanr(x[ok], y[ok])
                    rows.append({"Subtype": sub, "feature": f, "outcome": outcome,
                                 "n_lineages": int(ok.sum()), "rho": r, "p": p})
    res = pd.DataFrame(rows)
    if len(res):   # Benjamini-Hochberg correction, since we test many combinations
        from statsmodels.stats.multitest import multipletests
        res["p_adj"] = multipletests(res["p"], method="fdr_bh")[1]
        res = res.sort_values("p")
    return t, res
