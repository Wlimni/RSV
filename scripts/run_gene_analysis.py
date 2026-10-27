"""
Script Name: run_gene_analysis.py
Purpose:     Per-gene version of the within/between-lineage trend decomposition, to find
             which RSV genes drive the composition signal (expected candidates: G and F).
Inputs:      data/processed/rsv_gene_features.tsv
Outputs:     results/gene_decomposition.tsv, results/gene_summary.tsv
Notes:       Run after build_gene_features.py.
             Run: .venv/bin/python scripts/run_gene_analysis.py
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import pandas as pd
from statsmodels.stats.multitest import multipletests
from rsv import trends

FEATURES = ["gc_content", "oe_CG", "oe_GG"]
GENE_ORDER = ["NS1", "NS2", "N", "P", "M", "SH", "G", "F", "M2", "L"]

g = pd.read_csv("data/processed/rsv_gene_features.tsv", sep="\t")
g = g[(g.year >= 2015) & (g.year <= 2026)]
print(f"{len(g)} gene-rows, {g['Accession ID'].nunique()} sequences, {g.gene.nunique()} genes")

# baseline composition per gene (how different are the genes from each other?)
summ = (g.groupby(["Subtype", "gene"])
          .agg(n=("gc_content", "size"), gc=("gc_content", "mean"),
               oe_CG=("oe_CG", "mean"), oe_GG=("oe_GG", "mean"),
               length=("valid_len", "median")).reset_index())
summ.to_csv("results/gene_summary.tsv", sep="\t", index=False)
print("\n=== composition by gene ===")
print(summ.round(4).to_string(index=False))

# within/between decomposition, per gene per subtype
rows = []
for (sub, gene), s in g.groupby(["Subtype", "gene"]):
    if len(s) < 500 or s.Lineage.nunique() < 3:
        continue
    for f in FEATURES:
        try:
            r = trends.decompose_trend(s, f)
        except Exception:
            continue
        r.insert(0, "gene", gene)
        r.insert(0, "Subtype", sub)
        rows.append(r)
dec = pd.concat(rows, ignore_index=True)
dec["p_within_adj"] = multipletests(dec["p_within"], method="fdr_bh")[1]
dec["gene"] = pd.Categorical(dec["gene"], GENE_ORDER, ordered=True)
dec = dec.sort_values(["Subtype", "feature", "gene"])
dec.to_csv("results/gene_decomposition.tsv", sep="\t", index=False)

print("\n=== per-gene CpG O/E trend (within-lineage, BH-corrected) ===")
sel = dec[dec.feature == "oe_CG"][["Subtype", "gene", "n", "slope_total",
                                   "slope_within", "p_within_adj"]]
print(sel.round(6).to_string(index=False))
print("\nwrote results/gene_decomposition.tsv, results/gene_summary.tsv")
