"""
Script Name: run_codon_analysis.py
Purpose:     Test whether CpG depletion in each RSV gene is forced by the protein or is
             selection on the RNA itself. For every gene we compare the real CpG count with
             the minimum and maximum achievable by synonymous recoding of the same protein.
Inputs:      data/full/RSV_full.fasta, data/processed/rsv_features.tsv,
             data/reference/NC_038235.1.{fasta,gb} (A), NC_001781.1.{fasta,gb} (B)
Outputs:     results/codon_cpg_analysis.tsv (per subtype per gene),
             results/codon_cpg_per_sequence.tsv (aggregated, no accessions)
Notes:       Uses a year-stratified subsample because the recoding search is per sequence.
             frac_of_range near 0 means every removable CpG has been removed.
             Run: .venv/bin/python scripts/run_codon_analysis.py [--n 600]
"""
import os, sys, argparse
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np, pandas as pd
from rsv import genes, io, codon

ap = argparse.ArgumentParser()
ap.add_argument("--n", type=int, default=600, help="sequences sampled per subtype")
args = ap.parse_args()

REF = {"A": "data/reference/NC_038235.1", "B": "data/reference/NC_001781.1"}
GENE_ORDER = ["NS1", "NS2", "N", "P", "M", "SH", "G", "F", "M2", "L"]

meta = pd.read_csv("data/processed/rsv_features.tsv", sep="\t")
meta = meta[(meta.year >= 2015) & (meta.year <= 2026)]
seqs = io.read_fasta("data/full/RSV_full.fasta")

rows = []
for sub, s in meta.groupby("Subtype"):
    # year-stratified sample so no single season dominates
    per_year = max(1, args.n // s.year.nunique())
    accs = []
    for _, g in s.groupby("year"):
        accs += g.sample(min(len(g), per_year), random_state=0)["Accession ID"].tolist()
    print(f"RSV-{sub}: {len(accs)} sequences sampled")
    ref_row, aligned = genes.align_to_reference(REF[sub] + ".fasta",
                                                {a: seqs[a] for a in accs})
    cds = genes.read_cds(REF[sub] + ".gb")
    per_gene = genes.slice_genes(ref_row, aligned, cds)
    for gene, d in per_gene.items():
        for acc, gs in d.items():
            if len(gs) < 60:
                continue
            r = codon.cpg_position_in_range(gs)
            r.update({"Subtype": sub, "gene": gene})
            rows.append(r)
        print(f"  {gene} done", flush=True)

df = pd.DataFrame(rows)
summ = (df.groupby(["Subtype", "gene"])
          .agg(n=("obs_cpg", "size"),
               obs=("obs_cpg", "mean"),
               min_syn=("min_syn_cpg", "mean"),
               max_syn=("max_syn_cpg", "mean"),
               frac_of_range=("frac_of_range", "mean"),
               pct_at_min=("obs_cpg", "size")).reset_index())
# share of sequences sitting exactly at the synonymous minimum
atmin = (df.assign(at_min=df.obs_cpg <= df.min_syn_cpg)
           .groupby(["Subtype", "gene"]).at_min.mean().mul(100).rename("pct_at_min"))
summ = summ.drop(columns=["pct_at_min"]).merge(atmin.reset_index(), on=["Subtype", "gene"])
summ["removable"] = summ["max_syn"] - summ["min_syn"]
summ["gene"] = pd.Categorical(summ["gene"], GENE_ORDER, ordered=True)
summ = summ.sort_values(["Subtype", "gene"])
summ.to_csv("results/codon_cpg_analysis.tsv", sep="\t", index=False)

print("\n=== CpG vs what the protein allows ===")
print("obs = real CpG count | min_syn/max_syn = achievable by synonymous recoding")
print("frac_of_range: 0 = all removable CpG removed, 1 = as much CpG as protein allows\n")
print(summ.round(3).to_string(index=False))
print("\nwrote results/codon_cpg_analysis.tsv")
