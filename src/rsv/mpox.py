"""
Script Name: mpox.py
Purpose:     Mutation-spectrum helpers for mpox (MPXV): parse Nextclade substitution strings and
             count APOBEC3-type changes. Called by the mpox analysis scripts; no I/O except
             reading the reference.
Inputs:      Nextclade substitution cells ("C123T,G456A"), the Nextclade reference FASTA.
Outputs:     Counts (total, APOBEC3-type, net G+C loss).
Notes:       APOBEC3-type = C->T where the reference base before it is T (TC -> TT), or G->A where
             the reference base after it is A (GA -> AA, the same event read on the other strand).
             Net G+C loss = (C->T + G->A) - (T->C + A->G). Context is taken from the reference, so
             for "private" mutations (relative to the nearest tree node) it is an approximation.
"""
import re

_SUB = re.compile(r"^([ACGT])(\d+)([ACGT])$")


def load_reference(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()


def parse(cell):
    """'C123T,G456A' -> list of (ref, pos, alt); indels and ambiguity codes are skipped."""
    if not isinstance(cell, str):
        return []
    out = []
    for x in cell.split(","):
        m = _SUB.match(x)
        if m:
            out.append((m.group(1), int(m.group(2)), m.group(3)))
    return out


def classify(subs, ref):
    """Return (number of substitutions, APOBEC3-type count, net G+C loss)."""
    apo = loss = 0
    for r, p, a in subs:
        if r == "C" and a == "T":
            loss += 1
            apo += p >= 2 and ref[p - 2] == "T"
        elif r == "G" and a == "A":
            loss += 1
            apo += p < len(ref) and ref[p] == "A"
        elif (r == "T" and a == "C") or (r == "A" and a == "G"):
            loss -= 1
    return len(subs), apo, loss
