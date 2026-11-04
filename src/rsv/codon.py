"""
Script Name: codon.py
Purpose:     Codon-level composition for in-frame coding sequences: GC3, effective number of
             codons (ENC), and CpG split by where it sits relative to the codon boundary.
Inputs:      Lists of in-frame coding sequences (stop codon included, ACGT only counted).
Outputs:     Plain dicts / numpy arrays. No file I/O here.
Notes:       CpG falls into three classes by codon position. Which class a CpG is in decides
             how much the protein constrains it:
               12  C at codon position 1, G at 2  (Arg CGN; the amino acid fixes it)
               23  C at 2, G at 3                 (Ser/Pro/Thr/Ala NCG; G3 is free, silent)
               31  C at 3, G at 1 of the next codon (C3 is usually free, G1 is not)
             So a CpG shortage in classes 23 and 31 but not 12 means the pressure acts on the
             choice of synonymous codon rather than on the protein.
             Counts are pooled across the genomes of a lineage before the O/E is taken, which
             keeps per-genome noise (a gene holds only a handful of CpGs) out of the ratio.
"""
import numpy as np
from Bio.Data import CodonTable

_IDX = np.full(256, -1, np.int8)
for i, b in enumerate("ACGT"):
    _IDX[ord(b)] = i
BASES = "ACGT"
CLASSES = ("12", "23", "31")

_FWD = CodonTable.unambiguous_dna_by_id[1].forward_table
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
AA_OF = {c: _FWD[c] for c in CODONS if c in _FWD}               # stops are absent
FAMILIES = {}
for c, aa in AA_OF.items():
    FAMILIES.setdefault(aa, []).append(CODONS.index(c))


def tally(seqs):
    """
    Pool codon counts, per-position base counts and CpG observations over a list of
    in-frame sequences. The final (stop) codon of each is dropped.
    Returns dict: codons[64], base[3,4], cpg_obs[3], pairs[3].
    """
    codons = np.zeros(64)
    base = np.zeros((3, 4))
    obs = np.zeros(3)
    pairs = np.zeros(3)
    for s in seqs:
        a = _IDX[np.frombuffer(s[:-3].encode(), np.uint8)]
        k = len(a) // 3
        a = a[:3 * k].reshape(k, 3)
        for j in range(3):
            v = a[:, j][a[:, j] >= 0]
            base[j] += np.bincount(v, minlength=4)
        ok = (a >= 0).all(axis=1)
        codons += np.bincount((a[ok, 0] * 16 + a[ok, 1] * 4 + a[ok, 2]).astype(int), minlength=64)
        C, G = 1, 2
        # class 12 and 23: both bases inside one codon
        for n, (x, y) in enumerate(((0, 1), (1, 2))):
            m = (a[:, x] >= 0) & (a[:, y] >= 0)
            pairs[n] += m.sum()
            obs[n] += ((a[:, x] == C) & (a[:, y] == G) & m).sum()
        # class 31: last base of a codon with the first base of the next
        m = (a[:-1, 2] >= 0) & (a[1:, 0] >= 0)
        pairs[2] += m.sum()
        obs[2] += ((a[:-1, 2] == C) & (a[1:, 0] == G) & m).sum()
    return {"codons": codons, "base": base, "cpg_obs": obs, "pairs": pairs}


def cpg_oe(t):
    """O/E for each CpG class, expected from the base frequencies at the two positions."""
    p = t["base"] / t["base"].sum(axis=1, keepdims=True)
    exp = np.array([p[0, 1] * p[1, 2], p[1, 1] * p[2, 2], p[2, 1] * p[0, 2]]) * t["pairs"]
    return dict(zip(CLASSES, t["cpg_obs"] / exp)), dict(zip(CLASSES, t["cpg_obs"]))


def gc3(t):
    b = t["base"][2]
    return (b[1] + b[2]) / b.sum()


def enc(t):
    """
    Wright's (1990) effective number of codons, 20 (extreme bias) to 61 (none).
    Homozygosity F is taken per amino acid and averaged within each degeneracy class.
    """
    c = t["codons"]
    F = {}
    for aa, idxs in FAMILIES.items():
        n = c[idxs].sum()
        if n > 1:
            p = c[idxs] / n
            F[aa] = ((n * (p ** 2).sum()) - 1) / (n - 1)
    out = 2.0       # Met and Trp, one codon each
    for size, aas in ((2, [a for a, i in FAMILIES.items() if len(i) == 2]),
                      (3, [a for a, i in FAMILIES.items() if len(i) == 3]),
                      (4, [a for a, i in FAMILIES.items() if len(i) == 4]),
                      (6, [a for a, i in FAMILIES.items() if len(i) == 6])):
        vals = [F[a] for a in aas if a in F]
        if vals and np.mean(vals) > 0:
            out += len(aas) / np.mean(vals)
    return min(out, 61.0)


def aa_conditioned_cpg_share(t):
    """
    Share of codons that carry a CpG inside the codon (class 12 or 23), given the amino
    acid. Compares codon *choice* with the protein held fixed: it is computed only over
    amino acids that have both a CpG-containing and a CpG-free codon.
    Returns (CpG codons, eligible codons).
    """
    c = t["codons"]
    has = np.array([("CG" in cod) for cod in CODONS])
    k = n = 0.0
    for aa, idxs in FAMILIES.items():
        sel = has[idxs]
        if sel.any() and not sel.all():
            k += c[idxs][sel].sum()
            n += c[idxs].sum()
    return k, n


# --- amino-acid-controlled test of the codon-junction CpG (class 31) ---------------------
# A junction CpG needs a C at codon position 3 followed by a G at position 1 of the next
# codon. G-initial codons belong only to Gly, Ala, Val, Asp and Glu, so when the next amino
# acid is one of those its first base is a G whatever codon is chosen. That lets the protein
# be held fixed: compare how often a codon ends in C before a G-initial amino acid versus
# before any other. With no CpG pressure the two rates match.
_G_INIT_AA = set("GAVDE")
_ENDS_C = np.array([c[2] == "C" for c in CODONS])
_G_INIT = np.array([AA_OF.get(c) in _G_INIT_AA for c in CODONS])
# codons whose amino acid offers both a C-ending and a non-C-ending synonym (a real choice)
_HAS_CHOICE = np.zeros(64, bool)
for _aa, _idx in FAMILIES.items():
    if _ENDS_C[_idx].any() and not _ENDS_C[_idx].all():
        _HAS_CHOICE[_idx] = True
_VALID_CODON = np.array([c in AA_OF for c in CODONS])


def junction_c3_rates(seqs):
    """
    Returns (c_before_G, n_before_G, c_before_other, n_before_other): among codons with a
    synonymous choice at position 3, how many end in C when the next codon must start with
    G, and when it need not.
    """
    cg = ng = co = no = 0
    for s in seqs:
        a = _IDX[np.frombuffer(s[:-3].encode(), np.uint8)]
        k = len(a) // 3
        a = a[:3 * k].reshape(k, 3)
        ok = (a >= 0).all(axis=1)
        idx = np.where(ok, a[:, 0] * 16 + a[:, 1] * 4 + a[:, 2], 0).astype(int)
        cur, nxt = idx[:-1], idx[1:]
        m = ok[:-1] & ok[1:] & _VALID_CODON[cur] & _VALID_CODON[nxt] & _HAS_CHOICE[cur]
        gi = _G_INIT[nxt]
        ec = _ENDS_C[cur]
        cg += (m & gi & ec).sum(); ng += (m & gi).sum()
        co += (m & ~gi & ec).sum(); no += (m & ~gi).sum()
    return cg, ng, co, no
