"""
Script Name: plot_subtype_contrast.py
Purpose:     One summary figure of the project's clearest finding: RSV-A has almost no CpG
             in its two interferon-antagonist genes (NS1, NS2) while RSV-B does, with GC
             content matched and the other genes near-identical between subtypes.
Inputs:      data/processed/rsv_gene_features.tsv
Outputs:     results/figures/fig6_subtype_contrast.png
Notes:       Descriptive. Run: .venv/bin/python scripts/plot_subtype_contrast.py
"""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

os.makedirs("results/figures", exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
PAL = {"A": "#2C7FB8", "B": "#C62824"}
GENES = ["NS1", "NS2", "N", "L"]

g = pd.read_csv("data/processed/rsv_gene_features.tsv", sep="\t")
g = g[g.gene.isin(GENES)]

stat = (g.groupby(["gene", "Subtype"])
          .apply(lambda s: pd.Series({
              "cpg": (s.freq_CG * s.valid_len).mean(),
              "zero": (s.freq_CG == 0).mean() * 100,
              "gc": s.gc_content.mean() * 100}), include_groups=False)
          .reset_index())

fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.4))
x = np.arange(len(GENES))
w = 0.38

def grouped(a, col, ylab, title):
    for i, s in enumerate("AB"):
        v = [stat[(stat.gene == gn) & (stat.Subtype == s)][col].iloc[0] for gn in GENES]
        b = a.bar(x + (i - .5) * w, v, w, color=PAL[s], label=f"RSV-{s}",
                  edgecolor="0.25", linewidth=.4)
        a.bar_label(b, fmt="%.2f", fontsize=8.5, padding=2)
    a.set_xticks(x); a.set_xticklabels(GENES)
    a.set_xlabel("gene"); a.set_ylabel(ylab)
    a.set_title(title, fontsize=11.5, weight="bold")

# log scale, otherwise the L gene's 40 CpG flattens the NS1/NS2 bars to nothing
grouped(ax[0], "cpg", "mean CpG count per gene (log scale)", "CpG content")
ax[0].set_yscale("log")
ax[0].set_ylim(0.08, 400)
ax[0].legend(frameon=False, loc="upper left")

grouped(ax[1], "zero", "% of genomes", "Genomes with zero CpG")
ax[1].set_ylim(0, 112)

grouped(ax[2], "gc", "mean GC content (%)", "GC content, for control")
ax[2].set_ylim(0, 52)

fig.suptitle("RSV-A has lost CpG from its two interferon-antagonist genes; RSV-B has not",
             fontsize=13, weight="bold", x=.012, ha="left")
fig.text(.012, .005,
         "NS1 and NS2 are the interferon antagonists. N and L are shown as controls. "
         "GC content is matched between subtypes, so the CpG difference is not a base-composition effect.",
         fontsize=9, color="0.35")
fig.tight_layout(rect=[0, .035, 1, .93])
fig.savefig("results/figures/fig6_subtype_contrast.png", dpi=200)
print("wrote results/figures/fig6_subtype_contrast.png")
print(stat.round(3).to_string(index=False))
