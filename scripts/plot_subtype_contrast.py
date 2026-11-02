"""
Script Name: plot_subtype_contrast.py
Purpose:     Compare CpG between RSV-A and RSV-B across all ten genes, on three scales:
             CpG O/E (is CpG avoided relative to the gene's own base composition?),
             CpG count, and GC count. Counts answer "how much is there", the percentage
             printed on each bar answers "how dense is it", and the fold change says how
             far apart the subtypes are.
Inputs:      data/processed/rsv_gene_features.tsv
Outputs:     results/figures/fig6_subtype_contrast.png
             results/gene_cpg_subtype_test.tsv
Notes:       CpG % uses OVERLAPPING dinucleotide positions, so the denominator is
             (length - 1), not length/2 -- a CpG can start at any base, and reading the
             gene as non-overlapping pairs would miss roughly half of them. Counts are
             recovered as frequency x positions, which matches a direct count except
             where an ambiguous base interrupts a pair.

             The two count panels use a log y axis: L is thirty times the length of SH,
             so a linear axis would flatten every short gene, including NS1 and NS2.

             With n = 15,458 vs 17,580 every gene reaches p < 1e-300, so p cannot rank
             anything. The figure reports fold change; the table additionally carries
             Cliff's delta with a DeLong 95% CI, the common-language effect size, and the
             mean difference with a Welch CI.
             Run: .venv/bin/python scripts/plot_subtype_contrast.py
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("results/figures", exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
PAL = {"A": "#2C7FB8", "B": "#C62824"}
GENES = ["NS1", "NS2", "N", "P", "M", "SH", "G", "F", "M2", "L"]   # genome order, 3' to 5'

# (column, panel label, y axis label, column holding the % printed on each bar, log axis)
METRICS = [("oe_CG",    "CpG O/E",  "CpG observed / expected",   None,      False),
           ("n_CG",     "CpG count", "CpG pairs per gene",       "pct_CG",  True),
           ("n_GC_base", "GC count", "G + C bases per gene",     "pct_GC_base", True)]


def cliffs_delta(a, b):
    """
    Cliff's delta with a DeLong 95% CI.

    delta = P(a > b) - P(a < b), from each value's "placement" -- the fraction of the
    other group it exceeds -- which is exact and cheap at n ~ 15,000, where the full
    pairwise matrix would be 270 million comparisons per gene.
    """
    a, b = np.sort(np.asarray(a)), np.sort(np.asarray(b))
    m, n = len(a), len(b)
    v10 = (np.searchsorted(b, a, "left") + np.searchsorted(b, a, "right")) / (2 * n)
    v01 = (np.searchsorted(a, b, "left") + np.searchsorted(a, b, "right")) / (2 * m)
    auc = v10.mean()                                  # P(a > b), ties split
    se = np.sqrt(v10.var(ddof=1) / m + v01.var(ddof=1) / n)
    return 2 * auc - 1, 2 * (auc - 1.96 * se) - 1, 2 * (auc + 1.96 * se) - 1, auc


g = pd.read_csv("data/processed/rsv_gene_features.tsv", sep="\t")
g = g[g.gene.isin(GENES)].copy()
# freq_CG is already count / (valid dinucleotide positions), i.e. the (length - 1) basis
g["pct_CG"] = g.freq_CG * 100
g["pct_GC_base"] = g.gc_content * 100
g["n_CG"] = g.freq_CG * (g.valid_len - 1)
g["n_GC_base"] = g.gc_content * g.valid_len

rows = []
for gn in GENES:
    s = g[g.gene == gn]
    for col, short, _, pct_col, _ in METRICS:
        a = s.loc[s.Subtype == "A", col].dropna().to_numpy()
        b = s.loc[s.Subtype == "B", col].dropna().to_numpy()
        d, lo, hi, auc = cliffs_delta(a, b)
        diff = b.mean() - a.mean()
        se = np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
        rows.append({
            "gene": gn, "metric": short,
            "A_mean": a.mean(), "B_mean": b.mean(),
            "A_pct": s.loc[s.Subtype == "A", pct_col].mean() if pct_col else np.nan,
            "B_pct": s.loc[s.Subtype == "B", pct_col].mean() if pct_col else np.nan,
            "fold_B_over_A": b.mean() / a.mean() if a.mean() else np.nan,
            "diff_B_minus_A": diff,
            "diff_lo": diff - 1.96 * se, "diff_hi": diff + 1.96 * se,
            "cliffs_delta": d, "delta_lo": lo, "delta_hi": hi,
            # common-language effect size: chance a random B genome exceeds a random A one
            "P_B_gt_A": 1 - auc,
            "p": mannwhitneyu(a, b, alternative="two-sided")[1],
            "n_A": len(a), "n_B": len(b)})
t = pd.DataFrame(rows)
t.to_csv("results/gene_cpg_subtype_test.tsv", sep="\t", index=False)

x = np.arange(len(GENES))
w = 0.38
fig, axes = plt.subplots(3, 1, figsize=(12, 10.6), sharex=True)

for ax, (col, short, ylab, pct_col, log) in zip(axes, METRICS):
    sub = t[t.metric == short].set_index("gene").loc[GENES]
    for i, s in enumerate("AB"):
        ax.bar(x + (i - .5) * w, sub[f"{s}_mean"], w, color=PAL[s],
               label=f"RSV-{s}", edgecolor="0.25", linewidth=.4)
    ax.set_ylabel(ylab)
    ax.set_title(short, fontsize=11.5, weight="bold", loc="left")
    if log:
        ax.set_yscale("log")
        vals = sub[["A_mean", "B_mean"]].to_numpy()
        # headroom only where fold-change labels sit above the bars
        ax.set_ylim(vals.min() * .45, vals.max() * (9 if short != "GC count" else 2.6))
        # the percentage each count represents, so density is readable off the same bar
        for xi, gn in zip(x, GENES):
            r = sub.loc[gn]
            for i, s in enumerate("AB"):
                ax.text(xi + (i - .5) * w, r[f"{s}_mean"] * 1.12, f"{r[f'{s}_pct']:.2f}%",
                        ha="center", fontsize=6.9, color=PAL[s], rotation=90, va="bottom")
    else:
        ax.set_ylim(0, 1.1)
    # fold change above each gene: the most direct statement of the subtype gap
    if short != "GC count":
        hi = sub[["A_mean", "B_mean"]].max(axis=1)
        for xi, gn in zip(x, GENES):
            r, y = sub.loc[gn], hi[gn]
            ax.text(xi, y * 3.4 if log else y + .045, f"{r.fold_B_over_A:.2f}x",
                    ha="center", fontsize=8.2, weight="bold",
                    color=PAL["B"] if r.fold_B_over_A > 1 else PAL["A"])

# O/E has an absolute reference: 1.0 means "as expected from this gene's own base content"
axes[0].axhline(1, ls="--", lw=.9, color="0.35")
axes[0].text(-.45, 1.015, "1.0 = CpG as expected from this gene's C and G content",
             fontsize=8.2, color="0.35", ha="left")
axes[0].legend(frameon=False, loc="upper right", ncol=2, bbox_to_anchor=(1, .93))
axes[2].set_xticks(x); axes[2].set_xticklabels(GENES)
axes[2].set_xlabel("Gene (genome order, 3' to 5')")

fig.suptitle("CpG differs between the subtypes in every gene, and NS1/NS2 are the extreme case",
             fontsize=13, weight="bold", x=.008, y=.992, ha="left")
fig.text(.008, .006,
         "Fold change = RSV-B mean / RSV-A mean; red where B is higher, blue where A is higher. "
         "Shown in place of p, which is < 1e-300 for every gene at n = 15,458 vs 17,580.\n"
         "Counts are on a log axis because L is thirty times the length of SH. The percentage on each bar is that "
         "count as a share of the gene: CpG over overlapping dinucleotide\npositions (denominator length - 1), GC over bases. "
         "Cliff's delta, 95% CIs and mean differences: results/gene_cpg_subtype_test.tsv.",
         fontsize=8.2, color="0.35", va="bottom", linespacing=1.5)
fig.tight_layout(rect=[0, .062, 1, .975])
fig.savefig("results/figures/fig6_subtype_contrast.png", dpi=200)
print("wrote results/figures/fig6_subtype_contrast.png")
for short in ["CpG O/E", "CpG count", "GC count"]:
    print(f"\n--- {short} ---")
    print(t[t.metric == short][["gene", "A_mean", "B_mean", "A_pct", "B_pct",
                                "fold_B_over_A", "cliffs_delta"]].round(3).to_string(index=False))
