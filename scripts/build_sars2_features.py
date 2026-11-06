"""
Script Name: build_sars2_features.py
Purpose:     SARS-CoV-2 genome-level composition for each genome: GC, observed/expected for CpG, UpA
             (TpA), GpC and TpC, length and ambiguous fraction.
Inputs:      data/processed/Cov_19/Cov_19.fasta + Cov_19_metadata.tsv (merge_gisaid_augur.py)
Outputs:     data/processed/Cov_19/Cov_19_features.tsv  one row per genome
Notes:       Whole-genome (about 29.8 kb), written in DNA letters (T for U). Ambiguous bases are
             ignored in all counts. Run: .venv/bin/python scripts/build_sars2_features.py
"""
import numpy as np, pandas as pd

LUT = np.full(256, 4, np.uint8)
for i, b in enumerate(b"ACGT"):
    LUT[b] = i

rows, acc, buf = [], None, []


def flush():
    if acc is None:
        return
    a = LUT[np.frombuffer("".join(buf).encode(), np.uint8)]
    valid = a < 4
    nv = int(valid.sum())
    base = np.bincount(a[valid], minlength=4) / nv
    pair = valid[:-1] & valid[1:]
    dn = np.bincount(a[:-1][pair] * 4 + a[1:][pair], minlength=16) / max(pair.sum(), 1)
    oe = lambda x, y: dn[x * 4 + y] / (base[x] * base[y])
    rows.append({"gisaid_epi_isl": acc, "length": len(a), "ambig_frac": 1 - nv / len(a),
                 "gc_content": base[1] + base[2], "oe_CG": oe(1, 2), "oe_TA": oe(3, 0),
                 "oe_GC": oe(2, 1), "oe_TC": oe(3, 1)})


for line in open("data/processed/Cov_19/Cov_19.fasta"):
    line = line.strip()
    if line.startswith(">"):
        flush(); acc, buf = line[1:], []
    else:
        buf.append(line)
flush()
f = pd.DataFrame(rows)
meta = pd.read_csv("data/processed/Cov_19/Cov_19_metadata.tsv", sep="\t", dtype=str)
f = f.merge(meta[["gisaid_epi_isl", "strain", "date", "country", "division", "host", "pangolin_lineage",
                  "GISAID_clade", "search", "originating_lab", "submitting_lab"]], on="gisaid_epi_isl")
f.to_csv("data/processed/Cov_19/Cov_19_features.tsv", sep="\t", index=False)
print(len(f), "genomes")
print(f.groupby("host")[["gc_content", "oe_CG", "oe_TA", "oe_TC", "ambig_frac"]].mean().round(4))
