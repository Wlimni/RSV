"""
Script Name: composition.py
Purpose:     Core nucleotide-composition features — GC content, the four GC-related
             dinucleotide frequencies (CpG, GpC, GpG, CpC) with their observed/expected
             ratios, GC at separate codon positions (GC3 / GC12), and UpA O/E.
             Called by scripts and notebooks; no I/O here.
Inputs:      A DNA sequence string (A/C/G/T, plus ambiguity codes which are ignored).
             Codon-position features additionally require the sequence to be in frame.
Outputs:     Python dicts / floats of composition metrics.
Notes:       O/E ratio for dinucleotide XY = freq(XY) / (freq(X) * freq(Y)), computed over
             valid ACGT positions only so ambiguous bases (N, gaps) do not distort counts.

             GC3 / GC12 follow Ye et al. 2026 (Methods): GC measured at the third codon
             site, versus the first and second. Third positions are mostly synonymous, so
             a faster change there points to pressure on the RNA rather than the protein.
             That paper reports GC3 falling 5.65x faster than GC12.

             UpA (TpA in DNA letters) is not in that paper. It is included because the
             RNA-virus literature treats UpA as CpG's companion motif, both suppressed and
             both linked to innate immune recognition (Greenbaum et al. 2008, PLoS Pathog
             4:e1000079). It is less well characterised than CpG, so treat it as
             exploratory rather than as an established result.
"""
from collections import Counter

DINUCS = ("CG", "GC", "GG", "CC")  # the four GC-related dinucleotides used in the paper
EXTRA_DINUCS = ("TA",)             # UpA, written TpA in DNA letters
ALL_DINUCS = DINUCS + EXTRA_DINUCS


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
        return {d: None for d in ALL_DINUCS}
    return {d: obs[d] / pair_total for d in ALL_DINUCS}


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
        return {d: None for d in ALL_DINUCS}
    bf = {b: n[b] / total for b in "ACGT"}
    out = {}
    for d in ALL_DINUCS:
        exp = bf[d[0]] * bf[d[1]]
        out[d] = (freqs[d] / exp) if (exp > 0 and freqs[d] is not None) else None
    return out


def gc_by_codon_position(seq):
    """
    GC content split by position within the codon, for an in-frame coding sequence.

    gc3  : GC at the third codon site, where most changes are synonymous
    gc12 : GC at the first and second sites together, where most changes alter the protein

    Comparing how fast these two move separates pressure on the RNA from pressure on the
    protein, in the same way that dS and dN are compared. Ye et al. 2026 report GC3
    changing 5.65x faster than GC12 in persistent mammalian influenza lineages.
    Returns None for a site class with no valid bases.
    """
    s = seq.upper()
    n = len(s) - len(s) % 3
    valid = set("ACGT")
    gc3 = tot3 = gc12 = tot12 = 0
    for i in range(0, n, 3):
        for j, ch in enumerate(s[i:i + 3]):
            if ch not in valid:
                continue
            is_gc = ch in "GC"
            if j == 2:
                tot3 += 1
                gc3 += is_gc
            else:
                tot12 += 1
                gc12 += is_gc
    return {"gc3": gc3 / tot3 if tot3 else None,
            "gc12": gc12 / tot12 if tot12 else None}


def all_features(seq, coding=False):
    """
    Full feature dict for one sequence.

    Always: GC content, dinucleotide frequencies (as a fraction and per 1,000 positions,
    the unit the paper reports) and O/E ratios, for CpG, GpC, GpG, CpC and UpA.
    With coding=True the sequence is treated as in-frame and gc3 / gc12 are added.
    """
    feats = {"gc_content": gc_content(seq), "valid_len": sum(base_counts(seq).values())}
    for d, v in dinucleotide_freqs(seq).items():
        feats[f"freq_{d}"] = v
        feats[f"per1k_{d}"] = v * 1000 if v is not None else None
    for d, v in dinucleotide_oe(seq).items():
        feats[f"oe_{d}"] = v
    if coding:
        feats.update(gc_by_codon_position(seq))
    return feats
