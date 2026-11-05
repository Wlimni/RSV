"""
Script Name: build_niv_features.py
Purpose:     Nipah virus (NiV): align every genome and fragment to the RefSeq reference, cut out the
             five annotated genes (N, M, F, G, L) and compute GC, dinucleotide O/E and GC3/GC12.
Inputs:      data/processed/NiV/NiV.fasta + NiV_metadata.tsv (from merge_gisaid.py),
             data/reference/NC_002728.1.{fasta,gb} (Malaysia strain, public NCBI RefSeq)
Outputs:     data/processed/NiV/NiV_gene_features.tsv   one row per sequence per gene (and GENOME)
Notes:       Most bat sequences are short fragments (median length 1.2 kb), so a gene is analysed only
             when the sequence covers at least 90% of the reference CDS. GC3/GC12 additionally need
             the gene to start with ATG and be a whole number of codons (read to its first stop);
             otherwise they are left empty. The GENOME row needs a full-length sequence (>= 17 kb).
             P/V/W/C overlap one another and are not analysed.
             Run: .venv/bin/python scripts/build_niv_features.py
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import pandas as pd
from Bio.Seq import Seq
from rsv import genes, io, composition

REF_FA, REF_GB = "data/reference/NC_002728.1.fasta", "data/reference/NC_002728.1.gb"
cds = genes.read_cds(REF_GB)
seqs = io.read_fasta("data/processed/NiV/NiV.fasta")
meta = pd.read_csv("data/processed/NiV/NiV_metadata.tsv", sep="\t", dtype=str)
print(len(seqs), "sequences;", list(cds))

ref_row, aligned = genes.align_to_reference(REF_FA, seqs)
per_gene = genes.slice_genes(ref_row, aligned, cds)

rows = []
for gene, d in per_gene.items():
    ref_len = cds[gene][1] - cds[gene][0] + 1
    for acc, s in d.items():
        if len(s) < 0.9 * ref_len:
            continue
        f = composition.all_features(s)
        f.update(gc3=None, gc12=None)
        if s.startswith("ATG") and set(s) <= set("ACGT"):
            t = s[:len(s) // 3 * 3]
            prot = str(Seq(t).translate())
            if "*" in prot:
                t = t[:3 * (prot.index("*") + 1)]
            if len(t) >= 0.9 * ref_len:
                f.update({k: v for k, v in composition.gc_by_codon_position(t).items()
                          if k in ("gc3", "gc12")})
        f.update({"Accession ID": acc, "gene": gene, "ref_cov": len(s) / ref_len})
        rows.append(f)
for acc, row in aligned.items():
    s = row.replace("-", "")
    if len(s) >= 17000:
        f = composition.all_features(s)
        f.update({"Accession ID": acc, "gene": "GENOME", "ref_cov": len(s) / len(ref_row.replace('-', ''))})
        rows.append(f)

out = pd.DataFrame(rows).merge(meta[["Accession ID", "Host", "Clade", "Lineage", "Location",
                                      "Collection date", "Sequence Length"]], on="Accession ID", how="left")
out.to_csv("data/processed/NiV/NiV_gene_features.tsv", sep="\t", index=False)
print(out.groupby("gene").size())
