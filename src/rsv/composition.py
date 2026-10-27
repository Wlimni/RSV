"""
Script Name: composition.py
Purpose:     Core nucleotide-composition features from the paper's method — GC content,
             the four GC-related dinucleotide frequencies (CpG, GpC, GpG, CpC) and their
             observed/expected (O/E) ratios. Called by scripts and notebooks; no I/O here.
Inputs:      A DNA sequence string (A/C/G/T, plus ambiguity codes which are ignored).
Outputs:     Python dicts / floats of composition metrics.
Notes:       O/E ratio for dinucleotide XY = freq(XY) / (freq(X) * freq(Y)), computed over
             valid ACGT positions only so ambiguous bases (N, gaps) do not distort counts.
"""
from collections import Counter

DINUCS = ("CG", "GC", "GG", "CC")  # the four GC-related dinucleotides used in the paper


def base_counts(seq):
    """Count A/C/G/T in an upper-case sequence, ignoring N, gaps and other codes."""
    c = Counter(seq.upper())
    return {b: c.get(b, 0) for b in "ACGT"}


def gc_content(seq):
    """Fraction of valid ACGT bases that are G or C. Returns None if no valid bases."""
    n = base_counts(seq)
    total = sum(n.values())
    if total == 0:
        return None
    return (n["G"] + n["C"]) / total


def dinucleotide_freqs(seq):
    """
    Observed frequency of each of the four dinucleotides, counted only where both
    positions are valid ACGT (any pair touching an ambiguous base is skipped).
    Frequency = count(XY) / (number of valid ACGT dinucleotide positions).
    """
    s = seq.upper()
    valid = set("ACGT")
    pair_total = 0
    obs = Counter()
    for a, b in zip(s, s[1:]):
        if a in valid and b in valid:
            pair_total += 1
            obs[a + b] += 1
    if pair_total == 0:
        return {d: None for d in DINUCS}
    return {d: obs[d] / pair_total for d in DINUCS}


def dinucleotide_oe(seq):
    """
    Observed/expected ratio for each dinucleotide: freq(XY) / (freq(X) * freq(Y)),
    where single-base freqs are taken over the same valid ACGT positions.
    Value ~1 means as-expected; <1 means suppressed (e.g. CpG suppression).
    """
    n = base_counts(seq)
    total = sum(n.values())
    freqs = dinucleotide_freqs(seq)
    if total == 0:
        return {d: None for d in DINUCS}
    bf = {b: n[b] / total for b in "ACGT"}
    out = {}
    for d in DINUCS:
        exp = bf[d[0]] * bf[d[1]]
        out[d] = (freqs[d] / exp) if (exp > 0 and freqs[d] is not None) else None
    return out


def all_features(seq):
    """Full feature dict for one sequence: GC content + 4 dinuc freqs + 4 O/E ratios."""
    feats = {"gc_content": gc_content(seq), "valid_len": sum(base_counts(seq).values())}
    for d, v in dinucleotide_freqs(seq).items():
        feats[f"freq_{d}"] = v
    for d, v in dinucleotide_oe(seq).items():
        feats[f"oe_{d}"] = v
    return feats
