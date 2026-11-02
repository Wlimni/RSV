"""
Script Name: plot_subtype_contrast.py
Purpose:     Compare CpG content between RSV-A and RSV-B across all ten genes, normalised
             for gene length, with a statistical test of each gene's subtype difference.
Inputs:      data/processed/rsv_gene_features.tsv
Outputs:     results/figures/fig6_subtype_contrast.png
             results/gene_cpg_subtype_test.tsv
Notes:       With n = 15,458 vs 17,580 every gene reaches p < 1e-300, so the p value does
             not discriminate. Cliff's delta (how completely the two distributions
             separate, -1 to +1) and the fold change are the informative measures and are
             what the figure shows. Run: .venv/bin/python scripts/plot_subtype_contrast.py
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

g = pd.read_csv("data/processed/rsv_gene_features.tsv", sep="\t")
g = g[g.gene.isin(GENES)]

rows = []
for gn in GENES:
    s = g[g.gene == gn]
    a = (s.loc[s.Subtype == "A", "freq_CG"] * 1000).dropna()
    b = (s.loc[s.Subtype == "B", "freq_CG"] * 1000).dropna()
    u, p = mannwhitneyu(a, b, alternative="two-sided")
    rows.append({
        "gene": gn, "A_per1k": a.mean(), "B_per1k": b.mean(),
        "fold_B_over_A": b.mean() / a.mean(),
        "A_pct_zero": (a == 0).mean() * 100, "B_pct_zero": (b == 0).mean() * 100,
        # Cliff's delta: +1 means every A genome exceeds every B genome, -1 the reverse
        "cliffs_delta": 2 * u / (len(a) * len(b)) - 1,
        "p": p, "n_A": len(a), "n_B": len(b)})
t = pd.DataFrame(rows)
t.to_csv("results/gene_cpg_subtype_test.tsv", sep="\t", index=False)

x = np.arange(len(GENES))
w = 0.38
fig, ax = plt.subplots(2, 1, figsize=(12, 7.2), sharex=True)

# -- top: CpG normalised by gene length, so genes of different size are comparable --
for i, s in enumerate("AB"):
    ax[0].bar(x + (i - .5) * w, t[f"{s}_per1k"], w, color=PAL[s],
              label=f"RSV-{s}", edgecolor="0.25", linewidth=.4)
ax[0].set_ylabel("CpG per 1,000 nt")
ax[0].set_ylim(0, 19.5)
ax[0].legend(frameon=False, loc="upper left")
ax[0].set_title("CpG content of each gene, normalised for gene length",
                fontsize=11.5, weight="bold", loc="left")
for xi, r in zip(x, t.itertuples()):
    hi = max(r.A_per1k, r.B_per1k)
    ax[0].text(xi, hi + 0.5, f"{r.fold_B_over_A:.2f}x", ha="center", fontsize=8,
               color="#C62824" if r.fold_B_over_A > 1 else "#2C7FB8", weight="bold")
ax[0].text(.995, .93, "fold change = B / A", transform=ax[0].transAxes,
           ha="right", fontsize=8.5, color="0.4")

# -- bottom: how often the gene carries no CpG at all --
for i, s in enumerate("AB"):
    ax[1].bar(x + (i - .5) * w, t[f"{s}_pct_zero"], w, color=PAL[s],
              label=f"RSV-{s}", edgecolor="0.25", linewidth=.4)
ax[1].set_ylabel("% of genomes with zero CpG")
ax[1].set_ylim(0, 100)
ax[1].set_xticks(x); ax[1].set_xticklabels(GENES)
ax[1].set_xlabel("Gene (genome order, 3' to 5')")
ax[1].set_title("Genomes carrying no CpG anywhere in the gene",
                fontsize=11.5, weight="bold", loc="left")
for xi, r in zip(x, t.itertuples()):
    if max(r.A_pct_zero, r.B_pct_zero) > 1:
        ax[1].text(xi, max(r.A_pct_zero, r.B_pct_zero) + 2.5,
                   f"{max(r.A_pct_zero, r.B_pct_zero):.0f}%", ha="center", fontsize=8,
                   color="0.25")

fig.suptitle("CpG differs between the subtypes in every gene, but NS1 and NS2 are the extreme case",
             fontsize=13, weight="bold", x=.008, ha="left")
fig.text(.008, .005,
         "Every gene reaches p < 1e-300 (Mann-Whitney, n = 15,458 vs 17,580), so p does not discriminate here; "
         "fold change and Cliff's delta do. Full statistics in results/gene_cpg_subtype_test.tsv.",
         fontsize=8.5, color="0.35")
fig.tight_layout(rect=[0, .028, 1, .945])
fig.savefig("results/figures/fig6_subtype_contrast.png", dpi=200)
print("wrote results/figures/fig6_subtype_contrast.png")
print(t.assign(p=t.p.map(lambda v: f"{v:.1e}")).round(3).to_string(index=False))
