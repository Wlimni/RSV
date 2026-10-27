"""
Script Name: codon.py
Purpose:     Decide whether CpG depletion in a coding gene reflects selection on the RNA
             itself or just the constraints of the protein it encodes.
Inputs:      A coding sequence (in frame, start of CDS at position 0).
Outputs:     Observed CpG count, plus the minimum and maximum CpG achievable by synonymous
             recoding of the same protein, and a position of the observed value in that range.
Notes:       The min/max are found by dynamic programming over codons. CpG can form inside a
             codon or across a codon boundary, so the DP state is the last base of the codon
             chosen so far. If the observed count sits at the synonymous minimum, the virus has
             removed every CpG it could remove without changing the protein, which points to
             selection acting on the RNA (e.g. by the antiviral protein ZAP) rather than on the
             protein sequence.
"""
from collections import defaultdict

BASES = "TCAG"
CODON_TABLE = {}
AAS = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
_i = 0
for _b1 in BASES:
    for _b2 in BASES:
        for _b3 in BASES:
            CODON_TABLE[_b1 + _b2 + _b3] = AAS[_i]
            _i += 1

# amino acid -> list of codons encoding it
SYN = defaultdict(list)
for _c, _a in CODON_TABLE.items():
    SYN[_a].append(_c)


def translate(seq):
    """Translate an in-frame sequence; returns protein string ('X' for incomplete/ambiguous)."""
    out = []
    for i in range(0, len(seq) - len(seq) % 3, 3):
        out.append(CODON_TABLE.get(seq[i:i + 3], "X"))
    return "".join(out)


def count_cpg(seq):
    return sum(1 for a, b in zip(seq, seq[1:]) if a == "C" and b == "G")


def _codon_cpg(codon):
    """CpG entirely inside this codon."""
    return (codon[0] == "C" and codon[1] == "G") + (codon[1] == "C" and codon[2] == "G")


def synonymous_cpg_range(seq):
    """
    Minimum and maximum CpG count over all sequences encoding the same protein.

    DP over codons; the state is the last base of the codon just placed, because a CpG can
    straddle the boundary (previous codon ends in C, next begins with G).
    Returns (observed, min_possible, max_possible, n_codons_used).
    Codons containing ambiguous bases are kept fixed (they contribute as-is and reset state).
    """
    n = len(seq) - len(seq) % 3
    codons = [seq[i:i + 3] for i in range(0, n, 3)]
    obs = count_cpg(seq)

    # DP tables: state = last base of placed codon
    best_min = {b: 0 for b in "ACGT"}
    best_max = {b: 0 for b in "ACGT"}
    started = False
    used = 0

    for cod in codons:
        aa = CODON_TABLE.get(cod)
        options = SYN[aa] if aa and aa != "*" and cod in CODON_TABLE else [cod]
        if aa == "*":
            options = [cod]            # keep stop codon fixed
        if cod not in CODON_TABLE:
            options = [cod]            # ambiguous codon: leave as is
        else:
            used += 1
        new_min, new_max = {}, {}
        for opt in options:
            inside = _codon_cpg(opt)
            last = opt[2]
            if not started:
                cand_lo = cand_hi = inside
            else:
                cand_lo = min(best_min[p] + inside + (p == "C" and opt[0] == "G")
                              for p in best_min)
                cand_hi = max(best_max[p] + inside + (p == "C" and opt[0] == "G")
                              for p in best_max)
            if last not in new_min or cand_lo < new_min[last]:
                new_min[last] = cand_lo
            if last not in new_max or cand_hi > new_max[last]:
                new_max[last] = cand_hi
        # states not reachable this round get a large/small sentinel
        for b in "ACGT":
            if b not in new_min:
                new_min[b] = min(new_min.values()) + 10_000
            if b not in new_max:
                new_max[b] = max(new_max.values()) - 10_000
        best_min, best_max = new_min, new_max
        started = True

    if not started:
        return obs, 0, 0, 0
    return obs, min(best_min.values()), max(best_max.values()), used


def cpg_position_in_range(seq):
    """
    Where the real sequence sits between the synonymous minimum and maximum, as a fraction.
    0.0 = every removable CpG has been removed (strong RNA-level selection against CpG)
    1.0 = as many CpG as the protein allows
    Returns dict with the raw counts too. None if the range is degenerate.
    """
    obs, lo, hi, used = synonymous_cpg_range(seq)
    frac = None if hi == lo else (obs - lo) / (hi - lo)
    return {"obs_cpg": obs, "min_syn_cpg": lo, "max_syn_cpg": hi,
            "codons": used, "frac_of_range": frac}
