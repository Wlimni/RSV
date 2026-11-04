"""
Script Name: build_codon_data.py
Purpose:     Track 1 data for notebook 05. Draws a lineage-balanced sample of genomes, aligns
             them to the subtype reference, and writes each gene's in-frame coding sequence
             as FASTA so codon-level composition and HyPhy can use them.
Inputs:      data/full/RSV_full.fasta, data/processed/rsv_features.tsv,
             data/reference/NC_038235.1.{gb,fasta} (A), NC_001781.1.{gb,fasta} (B)
Outputs:     data/processed/cds/{A,B}_{gene}.fasta  (accession-named, ungapped, in frame)
             data/processed/cds/sample.tsv          (the sampled genomes and their lineage)
             data/processed/cds/frame_report.tsv    (genes kept / dropped per subtype)
Notes:       Same 9 lineages and date window as the trend analysis (lineage.prepare), at most
             --per-lineage genomes each so one big lineage cannot dominate. A gene is kept only
             if it starts with ATG, has no more than 1% N and, after cutting at its first
             in-frame stop, is at least 90% of the reference length. That drops start-truncated
             genes, which are off frame inside their own slice. G needs the stop cut because
             its stop codon moves between lineages.
             Run: .venv/bin/python scripts/build_codon_data.py [--per-lineage 250]
"""
import os, sys, argparse
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import pandas as pd
from Bio.Seq import Seq
from rsv import genes, io, lineage

ap = argparse.ArgumentParser()
ap.add_argument("--per-lineage", type=int, default=250)
ap.add_argument("--seed", type=int, default=0)
args = ap.parse_args()

REF = {"A": ("data/reference/NC_038235.1.fasta", "data/reference/NC_038235.1.gb"),
       "B": ("data/reference/NC_001781.1.fasta", "data/reference/NC_001781.1.gb")}
GENE_ORDER = ["NS1", "NS2", "N", "P", "M", "SH", "G", "F", "M2-1", "L"]
OUT = "data/processed/cds"
os.makedirs(OUT, exist_ok=True)

feat = pd.read_csv("data/processed/rsv_features.tsv", sep="\t", parse_dates=["date"])
d, _ = lineage.prepare(feat)
d = d.sample(frac=1, random_state=args.seed).groupby(["Subtype", "lineage_grp"]).head(args.per_lineage)
d[["Accession ID", "Subtype", "lineage_grp", "Lineage", "dec_date", "region"]].to_csv(
    f"{OUT}/sample.tsv", sep="\t", index=False)
print(d.groupby(["Subtype", "lineage_grp"]).size())

seqs = io.read_fasta("data/full/RSV_full.fasta")
report = []
for sub, grp in d.groupby("Subtype"):
    ref_fa, ref_gb = REF[sub]
    cds = {k: v for k, v in genes.read_cds(ref_gb).items() if k in GENE_ORDER}
    accs = grp["Accession ID"].tolist()
    ref_row, aligned = genes.align_to_reference(ref_fa, {a: seqs[a] for a in accs})
    per_gene = genes.slice_genes(ref_row, aligned, cds)
    for gene, dd in per_gene.items():
        kept = {}
        ref_len = cds[gene][1] - cds[gene][0] + 1
        for acc, s in dd.items():
            if not s.startswith("ATG") or sum(c not in "ACGT" for c in s) > 0.01 * len(s):
                continue                 # start-truncated genes are off frame inside the slice
            s = s[:len(s) // 3 * 3]
            prot = str(Seq(s.replace("N", "A")).translate())
            if "*" in prot:              # cut at the first stop: G's stop moves between lineages
                s = s[:3 * (prot.index("*") + 1)]
            if len(s) < 0.9 * ref_len or (len(s) % 3):
                continue
            kept[acc] = s
        with open(f"{OUT}/{sub}_{gene}.fasta", "w") as fh:
            for acc, s in kept.items():
                fh.write(f">{acc}\n{s}\n")
        report.append({"Subtype": sub, "gene": gene, "sampled": len(dd), "kept": len(kept)})
        print(sub, gene, len(kept), "/", len(dd), flush=True)
pd.DataFrame(report).to_csv(f"{OUT}/frame_report.tsv", sep="\t", index=False)
