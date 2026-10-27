"""
Script Name: build_gene_features.py
Purpose:     Step 2 of the plan — per-gene composition. Aligns every genome to its subtype's
             annotated reference, cuts out the 10 RSV genes, and computes GC content and
             dinucleotide O/E ratios for each gene separately.
Inputs:      data/full/RSV_full.fasta, data/processed/rsv_features.tsv (for metadata),
             data/reference/NC_038235.1.{gb,fasta} (A), data/reference/NC_001781.1.{gb,fasta} (B)
Outputs:     data/processed/rsv_gene_features.tsv  (one row per sequence per gene)
Notes:       Alignment keeps insertions, so the ON1 (A) and BA (B) G-gene duplications are
             retained. Sequences are processed in chunks to bound memory.
             Run: .venv/bin/python scripts/build_gene_features.py [--chunk 2000] [--max-per-subtype N]
"""
import os, sys, argparse, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import pandas as pd
from rsv import genes, io, composition

ap = argparse.ArgumentParser()
ap.add_argument("--chunk", type=int, default=2000)
ap.add_argument("--max-per-subtype", type=int, default=0, help="0 = all")
args = ap.parse_args()

REF = {"A": ("data/reference/NC_038235.1.fasta", "data/reference/NC_038235.1.gb"),
       "B": ("data/reference/NC_001781.1.fasta", "data/reference/NC_001781.1.gb")}

meta = pd.read_csv("data/processed/rsv_features.tsv", sep="\t",
                   usecols=["Accession ID", "Subtype", "Lineage", "year", "region"])
seqs = io.read_fasta("data/full/RSV_full.fasta")
print(f"loaded {len(seqs)} sequences")

rows = []
for sub, grp in meta.groupby("Subtype"):
    accs = grp["Accession ID"].tolist()
    if args.max_per_subtype:
        accs = accs[:args.max_per_subtype]
    ref_fa, ref_gb = REF[sub]
    cds = genes.read_cds(ref_gb)
    print(f"\nRSV-{sub}: {len(accs)} sequences, {len(cds)} genes")
    t0 = time.time()
    for i in range(0, len(accs), args.chunk):
        batch = {a: seqs[a] for a in accs[i:i + args.chunk] if a in seqs}
        ref_row, aligned = genes.align_to_reference(ref_fa, batch)
        per_gene = genes.slice_genes(ref_row, aligned, cds)
        for gene, d in per_gene.items():
            for acc, s in d.items():
                if len(s) < 50:          # skip genes missing from a partial genome
                    continue
                f = composition.all_features(s)
                f.update({"Accession ID": acc, "Subtype": sub, "gene": gene})
                rows.append(f)
        done = min(i + args.chunk, len(accs))
        el = time.time() - t0
        print(f"  {done}/{len(accs)}  ({el:.0f}s elapsed, {el/done*len(accs)-el:.0f}s left)", flush=True)

out = pd.DataFrame(rows)
out = out.merge(meta, on=["Accession ID", "Subtype"], how="left")
os.makedirs("data/processed", exist_ok=True)
out.to_csv("data/processed/rsv_gene_features.tsv", sep="\t", index=False)
print(f"\nwrote data/processed/rsv_gene_features.tsv: {len(out)} rows "
      f"({out['Accession ID'].nunique()} sequences x {out.gene.nunique()} genes)")
print(out.groupby(["Subtype", "gene"]).size().unstack(fill_value=0).to_string())
