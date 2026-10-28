"""
Script Name: run_gene_analysis.py
Purpose:     Summarise genome composition gene by gene, and measure how depleted of CpG
             each gene is. Descriptive only - no modelling.
Inputs:      data/processed/rsv_gene_features.tsv
Outputs:     results/gene_summary.tsv, results/cpg_depletion_by_gene.tsv
Notes:       Run after build_gene_features.py.
             Run: .venv/bin/python scripts/run_gene_analysis.py
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import pandas as pd

GENE_ORDER = ["NS1", "NS2", "N", "P", "M", "SH", "G", "F", "M2", "L"]
os.makedirs("results", exist_ok=True)

g = pd.read_csv("data/processed/rsv_gene_features.tsv", sep="\t")
print(f"{len(g)} gene-rows, {g['Accession ID'].nunique()} genomes, {g.gene.nunique()} genes")

summ = (g.groupby(["Subtype", "gene"])
          .agg(n=("gc_content", "size"), length=("valid_len", "median"),
               gc=("gc_content", "mean"),
               freq_CG=("freq_CG", "mean"), freq_GC=("freq_GC", "mean"),
               freq_GG=("freq_GG", "mean"),
               oe_CG=("oe_CG", "mean"), oe_GG=("oe_GG", "mean")).reset_index())
summ["gene"] = pd.Categorical(summ["gene"], GENE_ORDER, ordered=True)
summ = summ.sort_values(["Subtype", "gene"])
summ.to_csv("results/gene_summary.tsv", sep="\t", index=False)
print("\n=== Composition by gene ===")
print(summ.round(4).to_string(index=False))

# how much CpG each gene actually carries, and how often it carries none at all
dep = (g.groupby(["Subtype", "gene"])
         .apply(lambda s: pd.Series({
             "cpg_per_gene": (s.freq_CG * s.valid_len).mean(),
             "pct_zero_cpg": (s.freq_CG == 0).mean() * 100,
             "oe_CG": s.oe_CG.mean(),
             "gc": s.gc_content.mean()}), include_groups=False)
         .reset_index())
dep["gene"] = pd.Categorical(dep["gene"], GENE_ORDER, ordered=True)
dep = dep.sort_values(["Subtype", "gene"])
dep.to_csv("results/cpg_depletion_by_gene.tsv", sep="\t", index=False)
print("\n=== CpG depletion by gene ===")
print(dep.round(3).to_string(index=False))
print("\nwrote results/gene_summary.tsv, results/cpg_depletion_by_gene.tsv")
