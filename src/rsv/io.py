"""
Script Name: io.py
Purpose:     Load merged RSV FASTA + metadata, apply QC filters, and expose tidy tables.
Inputs:      Paths to RSV_full.fasta and RSV_full_metadata.tsv (from merge_gisaid_batches.py).
Outputs:     dict of {accession: sequence}; pandas DataFrame of metadata with parsed fields.
Notes:       Filtering is kept here so scripts and notebooks share one definition of "clean".
"""
import re
import pandas as pd


def read_fasta(path):
    """Return {accession: sequence}. Accession is the 2nd '|'-field of the GISAID header."""
    seqs, acc, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if acc:
                    seqs[acc] = "".join(buf)
                acc = line[1:].split("|")[1]
                buf = []
            else:
                buf.append(line)
    if acc:
        seqs[acc] = "".join(buf)
    return seqs


def load_metadata(path):
    """Load metadata TSV and add parsed columns: date, year, region, country."""
    m = pd.read_csv(path, sep="\t", dtype=str)
    m["Sequence Length"] = pd.to_numeric(m["Sequence Length"], errors="coerce")
    # keep only full YYYY-MM-DD collection dates; partial dates become NaT
    m["date"] = pd.to_datetime(
        m["Collection date"].where(m["Collection date"].str.match(r"^\d{4}-\d{2}-\d{2}$")),
        errors="coerce",
    )
    m["year"] = m["date"].dt.year
    loc = m["Location"].fillna("").str.split(" / ")
    m["region"] = loc.str[0].str.strip()
    m["country"] = loc.str[1].str.strip()
    return m


def ambiguous_fraction(seq):
    """Fraction of positions that are not A/C/G/T (N, gaps, ambiguity codes)."""
    s = seq.upper()
    if not s:
        return 1.0
    good = sum(s.count(b) for b in "ACGT")
    return 1 - good / len(s)


def qc_filter(meta, seqs, min_len=14900, max_len=15400, max_ambig=0.01, require_date=True):
    """
    Apply QC to the metadata table against the sequence dict.
    Keeps rows whose sequence length is in [min_len, max_len], ambiguous fraction
    <= max_ambig, and (optionally) a full collection date. Returns filtered meta
    plus a per-step count report (dict) for the QC notebook.
    """
    report = {"input": len(meta)}
    m = meta[meta["Accession ID"].isin(seqs)].copy()
    report["with_sequence"] = len(m)

    m["ambig_frac"] = m["Accession ID"].map(lambda a: ambiguous_fraction(seqs[a]))
    m["seq_len_actual"] = m["Accession ID"].map(lambda a: len(seqs[a]))

    m = m[m["seq_len_actual"].between(min_len, max_len)]
    report["pass_length"] = len(m)
    m = m[m["ambig_frac"] <= max_ambig]
    report["pass_ambiguity"] = len(m)
    if require_date:
        m = m[m["date"].notna()]
        report["pass_date"] = len(m)
    return m, report
