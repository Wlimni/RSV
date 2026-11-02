"""
Script Name: build_gene_features.py
Purpose:     Step 2 of the plan — per-gene composition. Aligns every genome to its subtype's
             annotated reference, cuts out the 10 RSV genes, and computes GC content and
             dinucleotide O/E ratios for each gene separately.
Inputs:      data/full/RSV_full.fasta, data/processed/rsv_features.tsv (for metadata),
             data/reference/NC_038235.1.{gb,fasta} (A), data/reference/NC_001781.1.{gb,fasta} (B)
Outputs:     data/processed/rsv_gene_features.tsv (one row per sequence per gene, plus a
             GENOME row per sequence holding the ten coding regions concatenated)
             Genome GC3/GC12 pool each gene's own in-frame site counts rather than being
             measured on the concatenation, which would be out of frame downstream of any
             gene whose length is not a multiple of 3.
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
# M2 carries two overlapping ORFs. Only M2-1 is analysed: M2-2 overlaps it by 35 nt,
# so including both would count that stretch twice in the concatenated genome.
GENE_ORDER = ["NS1", "NS2", "N", "P", "M", "SH", "G", "F", "M2-1", "L"]

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
    cds = {k: v for k, v in genes.read_cds(ref_gb).items() if k in GENE_ORDER}
    print(f"\nRSV-{sub}: {len(accs)} sequences, {len(cds)} genes")
    t0 = time.time()
    for i in range(0, len(accs), args.chunk):
        batch = {a: seqs[a] for a in accs[i:i + args.chunk] if a in seqs}
        ref_row, aligned = genes.align_to_reference(ref_fa, batch)
        per_gene = genes.slice_genes(ref_row, aligned, cds)
        per_acc_codon = {}           # accession -> list of (gc3, n3, gc12, n12) per gene
        for gene, d in per_gene.items():
            for acc, s in d.items():
                if len(s) < 50:          # skip genes missing from a partial genome
                    continue
                f = composition.all_features(s, coding=True)
                f.update({"Accession ID": acc, "Subtype": sub, "gene": gene})
                rows.append(f)
                per_acc_codon.setdefault(acc, []).append(f)
        # Genome-level row: the ten coding regions concatenated in genome order.
        # The paper defines "genomic" GC, GC3 and GC12 over the coding regions rather
        # than the raw genome, so the frame stays well defined.
        for acc in aligned:
            cat = "".join(per_gene[g].get(acc, "") for g in GENE_ORDER if g in per_gene)
            if len(cat) < 1000:
                continue
            f = composition.all_features(cat, coding=True)
            # GC3/GC12 must NOT come from the concatenation: a gene whose length is not a
            # multiple of 3 throws every downstream gene out of frame (about 4% of
            # genomes). Pool each gene's own in-frame counts instead.
            parts = per_acc_codon.get(acc, [])
            n3 = sum(x["n3"] for x in parts if x.get("gc3") is not None)
            n12 = sum(x["n12"] for x in parts if x.get("gc12") is not None)
            f["gc3"] = (sum(x["gc3"] * x["n3"] for x in parts if x.get("gc3") is not None) / n3
                        if n3 else None)
            f["gc12"] = (sum(x["gc12"] * x["n12"] for x in parts if x.get("gc12") is not None) / n12
                         if n12 else None)
            f.update({"Accession ID": acc, "Subtype": sub, "gene": "GENOME",
                      "n3": n3, "n12": n12})
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
