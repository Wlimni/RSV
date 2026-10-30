"""
Script Name: run_gene_analysis.py
Purpose:     Summarise genome composition gene by gene. Descriptive only - the trend
             fitting happens in run_lineage_analysis.py.
Inputs:      data/processed/rsv_gene_features.tsv
Outputs:     results/gene_summary.tsv
Notes:       Run after build_gene_features.py.
             Run: .venv/bin/python scripts/run_gene_analysis.py
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import pandas as pd

GENE_ORDER = ["NS1", "NS2", "N", "P", "M", "SH", "G", "F", "M2", "L"]
os.makedirs("results", exist_ok=True)

g = pd.read_csv("data/processed/rsv_gene_features.tsv", sep="\t")
g = g[g.gene != "GENOME"]
print(f"{len(g)} gene-rows, {g['Accession ID'].nunique()} genomes, {g.gene.nunique()} genes")

summ = (g.groupby(["Subtype", "gene"])
          .agg(n=("gc_content", "size"), length=("valid_len", "median"),
               gc=("gc_content", "mean"), gc3=("gc3", "mean"), gc12=("gc12", "mean"),
               per1k_CG=("per1k_CG", "mean"), per1k_GG=("per1k_GG", "mean"),
               oe_CG=("oe_CG", "mean"), oe_GG=("oe_GG", "mean"),
               cpg_per_gene=("freq_CG", "mean"),
               pct_zero_cpg=("freq_CG", lambda s: (s == 0).mean() * 100)).reset_index())
summ["cpg_per_gene"] = summ["cpg_per_gene"] * summ["length"]
summ["gene"] = pd.Categorical(summ["gene"], GENE_ORDER, ordered=True)
summ = summ.sort_values(["Subtype", "gene"])
summ.to_csv("results/gene_summary.tsv", sep="\t", index=False)

print("\n=== Composition by gene ===")
print(summ.round(4).to_string(index=False))
print("\nwrote results/gene_summary.tsv")
