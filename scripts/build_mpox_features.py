"""
Script Name: build_mpox_features.py
Purpose:     Mpox (MPXV) genome-level composition for every sequence, plus a rough clade call (I or II)
             by k-mer similarity, because the GISAID export has no clade column.
Inputs:      data/processed/Mpox/Mpox.fasta + Mpox_metadata.tsv (from merge_gisaid.py),
             data/reference/NC_003310.1.fasta (clade I, Zaire-96-I-16), NC_063383.1.fasta (clade IIb)
Outputs:     data/processed/Mpox/Mpox_features.tsv  one row per sequence
Notes:       Whole-genome only. Mpox is a ~197 kb double-stranded DNA virus with about 180 genes, so
             the per-gene pipeline used for RSV and NiV is not repeated here. Fields: length,
             ambiguous fraction, GC, and observed/expected for CpG, TpA, TpC and GpC. TpC matters
             because the host APOBEC3 enzymes edit C to T mostly after a T (TC -> TT), which is the
             main mutation signature of human-adapted MPXV.
             Clade call: share of sampled 21-mers found only in the clade I reference versus only in
             the clade IIb reference. Rough: it separates I from II, it does not split IIa/IIb.
             Run: .venv/bin/python scripts/build_mpox_features.py [Mpox|MpoxI]
"""
import sys
import numpy as np, pandas as pd

VIRUS = sys.argv[1] if len(sys.argv) > 1 else "Mpox"      # data folder name: Mpox or MpoxI

K, STEP = 21, 53
LUT = np.full(256, 4, np.uint8)
for i, b in enumerate(b"ACGT"):
    LUT[b] = i


def read_fa(path):
    h, s = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if h:
                    yield h, "".join(s).upper()
                h, s = line[1:].strip(), []
            else:
                s.append(line.strip())
    if h:
        yield h, "".join(s).upper()


def kmers(seq, step=1):
    return {seq[i:i + K] for i in range(0, len(seq) - K, step)}


ref_i = next(read_fa("data/reference/NC_003310.1.fasta"))[1]
ref_ii = next(read_fa("data/reference/NC_063383.1.fasta"))[1]
ki, kii = kmers(ref_i), kmers(ref_ii)
only_i, only_ii = ki - kii, kii - ki
print("reference k-mers:", len(only_i), "only in I,", len(only_ii), "only in II")

rows = []
for n, (h, s) in enumerate(read_fa(f"data/processed/{VIRUS}/{VIRUS}.fasta")):
    acc = h.split("|")[1]
    a = LUT[np.frombuffer(s.encode(), np.uint8)]
    valid = a < 4
    nv = int(valid.sum())
    if nv == 0:
        continue
    base = np.bincount(a[valid], minlength=4) / nv
    pair = valid[:-1] & valid[1:]
    code = a[:-1][pair] * 4 + a[1:][pair]
    dn = np.bincount(code, minlength=16) / max(len(code), 1)
    def oe(x, y):                               # indices: A0 C1 G2 T3
        return dn[x * 4 + y] / (base[x] * base[y]) if base[x] * base[y] else np.nan
    sample = {s[i:i + K] for i in range(0, len(s) - K, STEP)}
    rows.append({"Accession ID": acc, "length": len(s), "ambig_frac": 1 - nv / len(s),
                 "gc_content": base[1] + base[2], "oe_CG": oe(1, 2), "oe_TA": oe(3, 0),
                 "oe_TC": oe(3, 1), "oe_GC": oe(2, 1),
                 "hit_I": len(sample & only_i), "hit_II": len(sample & only_ii)})
    if n % 1000 == 0:
        print(n, flush=True)

f = pd.DataFrame(rows)
tot = f.hit_I + f.hit_II
f["clade_call"] = np.where(tot < 200, "unassigned", np.where(f.hit_I > f.hit_II, "I", "II"))
meta = pd.read_csv(f"data/processed/{VIRUS}/{VIRUS}_metadata.tsv", sep="\t", dtype=str)
f = f.merge(meta[["Accession ID", "Host", "Collection date", "Location", "Specimen",
                  "Additional host information", "Patient status"]], on="Accession ID", how="left")
f.to_csv(f"data/processed/{VIRUS}/{VIRUS}_features.tsv", sep="\t", index=False)
print(len(f), "sequences")
print(f.clade_call.value_counts())
