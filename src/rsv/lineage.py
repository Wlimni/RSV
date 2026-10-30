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
            start=2015.0, end=2026.0):
    """
    Add the truncated lineage label and decimal date, restrict to the usable date window,
    then keep only lineages with enough span and enough sequences to support a regression.
    Returns (filtered DataFrame, summary table of the lineages kept).

    The default window is 2015-01-01 to 2025-12-31. The lower bound is set because RSV-A
    has essentially no sequences before 2010 and everything before 2015 is about 3% of the
    data. The upper bound drops 2026, where deposition lag leaves the later months nearly
    empty (18 genomes in September 2026 against an all-year September average near 1,400),
    so that year would contribute a biased, early-clustered tail to every regression.
    """
    d = df.copy()
    d["lineage_grp"] = truncate_lineage(d["Lineage"].astype(str), level)
    d["dec_date"] = decimal_date(d[date_col])
    d = d[d["dec_date"].notna()]
    d = d[(d["dec_date"] >= start) & (d["dec_date"] < end)]

    g = d.groupby(["Subtype", "lineage_grp"]).agg(
        n=("dec_date", "size"),
        years=("dec_date", lambda s: s.astype(int).nunique()),
        first=("dec_date", "min"), last=("dec_date", "max"))
    keep = g[(g["years"] >= min_years) & (g["n"] >= min_n)].reset_index()
    d = d.merge(keep[["Subtype", "lineage_grp"]], on=["Subtype", "lineage_grp"])
    return d, keep.sort_values(["Subtype", "n"], ascending=[True, False])


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
