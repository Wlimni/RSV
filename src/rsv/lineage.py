"""
Script Name: lineage.py
Purpose:     Fit composition trends over time separately inside each RSV lineage, following
             the design of Ye et al. 2026 (their Figs 2 and 3).
Inputs:      A features DataFrame with columns Subtype, Lineage, date, and composition columns.
Outputs:     Tidy pandas DataFrames of regression slopes and correlations. No file I/O here.
Notes:       Fitting inside each lineage, rather than pooling everything, is what keeps a
             trend from being produced by one lineage replacing another.

             RSV lineage labels are hierarchical and dot-separated (A.D.1.5 sits inside
             A.D.1 inside A.D). Full labels are too fine for this: most span only ~6 years
             with a few hundred sequences, which makes every slope noisy. Truncating to
             three levels gives about 8 well-populated lineages spanning 2015-2026, which
             is close to the 12 lineages the paper worked with.
"""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf


def truncate_lineage(labels, level=3):
    """Cut a hierarchical lineage label down to `level` dot-separated parts."""
    return labels.str.split(".").str[:level].str.join(".")


def decimal_date(dates):
    """
    Convert dates to a decimal year (2017-07-02 -> about 2017.5).
    The paper regresses against fraction-of-year rather than whole years, which uses the
    full precision of the collection date.
    """
    d = pd.to_datetime(dates, errors="coerce")
    year_start = pd.to_datetime(d.dt.year.astype("Int64").astype(str) + "-01-01", errors="coerce")
    year_end = pd.to_datetime((d.dt.year.astype("Int64") + 1).astype(str) + "-01-01", errors="coerce")
    return d.dt.year + (d - year_start).dt.days / (year_end - year_start).dt.days


def prepare(df, level=3, min_years=5, min_n=200, date_col="date",
            start=2001.0, end=2026.5, trim=5):
    """
    Add the truncated lineage label and decimal date, restrict to the usable date window,
    trim isolated dates at each end of every lineage, then keep only lineages with enough
    span and enough sequences to support a regression.
    Returns (filtered DataFrame, summary table of the lineages kept).

    Window. The default runs from 2001 to the end of June 2026, which is as wide as the
    data usefully allows. Sequencing effort grows steeply over this period, but uneven
    sample sizes across time do not bias a regression slope; they only make the later years
    more precisely estimated. The later months of 2026 are different and are excluded:
    comparing like months against 2025 shows August 2026 holding 18 genomes against 399 a
    year earlier, which is deposition lag rather than season, so those genomes are a
    non-random subset contributed only by fast-depositing laboratories.

    Trimming. Opening the window back to 2001 lets a lineage's span be set by one or two
    very old genomes: lineage A.D.3 would otherwise reach back 6.9 extra years on the
    strength of a single 2010 sequence. A point that isolated has enormous leverage on a
    fitted line. So each lineage is cut back to the stretch of years over which it was
    actually being sampled: leading and trailing years holding fewer than `trim` genomes
    are dropped. This removes stragglers without penalising large lineages, which a
    percentile-based trim would, and keeps genuinely well-sampled early data such as
    lineage A.D's 327 genomes in its first three years.
    """
    d = df.copy()
    d["lineage_grp"] = truncate_lineage(d["Lineage"].astype(str), level)
    d["dec_date"] = decimal_date(d[date_col])
    d = d[d["dec_date"].notna()]
    d = d[(d["dec_date"] >= start) & (d["dec_date"] < end)]

    if trim:
        keep_idx = []
        for _, s in d.groupby(["Subtype", "lineage_grp"]):
            keep_idx.extend(_trim_sparse_ends(s, trim).index)
        d = d.loc[sorted(keep_idx)]

    g = d.groupby(["Subtype", "lineage_grp"]).agg(
        n=("dec_date", "size"),
        years=("dec_date", lambda s: s.astype(int).nunique()),
        first=("dec_date", "min"), last=("dec_date", "max"))
    keep = g[(g["years"] >= min_years) & (g["n"] >= min_n)].reset_index()
    d = d.merge(keep[["Subtype", "lineage_grp"]], on=["Subtype", "lineage_grp"])
    return d, keep.sort_values(["Subtype", "n"], ascending=[True, False])


def _trim_sparse_ends(s, min_per_year):
    """
    Cut a lineage back to the run of years it was actually sampled in: drop leading and
    trailing calendar years holding fewer than `min_per_year` genomes. Years inside the
    run are left alone, so a quiet season in the middle does not split the series.
    """
    yrs = s["dec_date"].astype(int)
    counts = yrs.value_counts().sort_index()
    ok = counts[counts >= min_per_year].index
    if len(ok) == 0:
        return s.iloc[0:0]
    return s[(yrs >= ok.min()) & (yrs <= ok.max())]


def lineage_slopes(df, features, by=("Subtype", "lineage_grp"), covariates=()):
    """
    Regress each feature against decimal collection date inside each lineage.

    Returns one row per lineage per feature with the slope (change per year), its 95%
    confidence interval, R squared and the sample size. This is the paper's Fig 2a and
    Fig 3a/3d in tabular form.

    `covariates` adds categorical controls, normally ("region",). This matters because
    where RSV was sequenced changed a great deal over the study period, and it changed
    even inside single lineages: A.D goes from 71% Oceania early to 36% Europe late.
    Fitting inside a lineage removes lineage replacement but not that geographic drift,
    so a slope that survives region adjustment is the more trustworthy one.
    """
    out = []
    cols = ["dec_date"] + list(covariates)
    ctrl = "".join(f" + C({c})" for c in covariates)
    for keys, g in df.groupby(list(by)):
        for f in features:
            s = g[[f] + cols].dropna()
            if len(s) < 30 or s["dec_date"].nunique() < 5:
                continue
            if any(s[c].nunique() < 2 for c in covariates):
                ctrl_use = ""            # covariate has no variation here; skip it
            else:
                ctrl_use = ctrl
            m = smf.ols(f"Q('{f}') ~ dec_date{ctrl_use}", data=s).fit()
            lo, hi = m.conf_int().loc["dec_date"]
            rec = dict(zip(by, keys if isinstance(keys, tuple) else (keys,)))
            rec.update({"feature": f, "n": len(s), "slope": m.params["dec_date"],
                        "ci_lo": lo, "ci_hi": hi, "p": m.pvalues["dec_date"],
                        "r2": m.rsquared, "mean": s[f].mean()})
            out.append(rec)
    res = pd.DataFrame(out)
    if len(res):
        res["significant"] = res["p"] < 0.05
    return res


def gc_vs_dinucleotide_correlation(df, dinuc_cols, gc_col="gc_content",
                                   by=("Subtype", "lineage_grp")):
    """
    Correlate each dinucleotide frequency with GC content inside each lineage.

    This asks whether the dinucleotide changes are simply a by-product of GC content
    moving, or something separate. The paper's Fig 3b.
    """
    from scipy.stats import pearsonr
    out = []
    for keys, g in df.groupby(list(by)):
        for c in dinuc_cols:
            s = g[[c, gc_col]].dropna()
            if len(s) < 30:
                continue
            r, p = pearsonr(s[c], s[gc_col])
            rec = dict(zip(by, keys if isinstance(keys, tuple) else (keys,)))
            rec.update({"feature": c, "n": len(s), "r": r, "p": p})
            out.append(rec)
    return pd.DataFrame(out)
