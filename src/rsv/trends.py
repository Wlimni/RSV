"""
Script Name: trends.py
Purpose:     Tell apart two reasons a composition average can move over time: the genomes
             inside each lineage actually changing, or one lineage simply replacing another.
Inputs:      A features DataFrame (one row per sequence) with columns
             Subtype, Lineage, year, region, and composition columns (gc_content, oe_CG, ...).
Outputs:     tidy pandas DataFrames of model results; no file I/O here.
Notes:       Adding lineage as a fixed effect absorbs the replacement part, so whatever year
             trend survives is change happening inside lineages.
"""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf


def decompose_trend(df, feature, group="Lineage", year="year", covariates=()):
    """
    Fit the year trend twice and compare.

    slope_total  : feature ~ year, across all sequences (the plain trend)
    slope_within : feature ~ year + C(lineage), holding lineage constant
    Opposite signs mean the plain trend is driven by lineage replacement.

    `covariates` (e.g. ("region",)) adds categorical controls to both fits, which is how
    we check a trend is not just a shift in which regions were sequencing.
    """
    d = df[[feature, year, group] + list(covariates)].dropna().copy()
    ctrl = "".join(f" + C({c})" for c in covariates)
    tot = smf.ols(f"{feature} ~ {year}{ctrl}", data=d).fit()
    win = smf.ols(f"{feature} ~ {year} + C({group}){ctrl}", data=d).fit()
    return pd.DataFrame([{
        "feature": feature, "n": len(d), "n_groups": d[group].nunique(),
        "slope_total": tot.params[year], "se_total": tot.bse[year], "p_total": tot.pvalues[year],
        "slope_within": win.params[year], "se_within": win.bse[year], "p_within": win.pvalues[year],
        "pct_from_turnover": 100 * (1 - win.params[year] / tot.params[year])
                             if tot.params[year] != 0 else np.nan,
    }])


def robustness(df, feature, group="Lineage", year="year"):
    """
    Repeat the decomposition under sensitivity scenarios, because sequencing effort is
    uneven: Europe dominates, and 2020-21 was disrupted by COVID-19.
    Returns one row per scenario per subtype.
    """
    out = []
    for sub, s in df.groupby("Subtype"):
        scenarios = {
            "full": s,
            "excl_2020_2021": s[~s[year].isin([2020, 2021])],
            "region_adjusted": s,          # same rows; region enters as a covariate below
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
    """Downsample so no one region dominates within a year (caps each region-year at the median)."""
    rng = np.random.default_rng(seed)
    keep = []
    for _, g in d.groupby(year):
        cap = int(max(1, g[col].value_counts().median()))
        for _, gg in g.groupby(col):
            idx = gg.index.values
            if len(idx) > cap:
                idx = rng.choice(idx, cap, replace=False)
            keep.extend(idx)
    return d.loc[sorted(keep)]
