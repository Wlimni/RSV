"""
Script Name: clusters.py
Purpose:     Cluster genomes by pairwise substitution differences, for the sporadic-versus-sustained
             analyses (mpox, SARS-CoV-2). Shared helpers; no file I/O except none.
Inputs:      Nextclade substitution strings, missing/deletion ranges, alignment start and end.
Outputs:     Integer codes, boolean masks, a distance matrix, a cluster table.
Notes:       Distance between two genomes = number of positions where their substitutions (against
             the Nextclade reference) differ, counted only at positions readable in both, so
             stretches of N never look like differences. Clusters = connected components of
             "distance <= T". Classes: sporadic = 1-2 genomes; sustained = at least 10 genomes
             spanning at least 0.5 year; the rest intermediate.
"""
import numpy as np, pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from . import mpox

ALT = {"A": 0, "C": 1, "G": 2, "T": 3}


def ranges_to_mask(cells, start, end, length):
    """Boolean array, True = unreadable (missing ranges, deletions, outside the alignment)."""
    m = np.zeros(length + 2, bool)
    m[: int(start)] = True
    m[int(end) + 1:] = True
    for cell in cells:
        if isinstance(cell, str):
            for part in cell.split(","):
                a, _, b = part.partition("-")
                try:
                    m[int(a): int(b or a) + 1] = True
                except ValueError:
                    pass
    return m


def encode(subs):
    """SNP list -> sorted int codes position*4 + alt."""
    return np.array(sorted({p * 4 + ALT[a] for r, p, a in subs}), dtype=np.int64)


def distances(codes, masks):
    n = len(codes)
    D = np.zeros((n, n), np.int32)
    for i in range(n):
        for j in range(i + 1, n):
            x = np.setxor1d(codes[i], codes[j], assume_unique=True)
            if len(x):
                pos = x >> 2
                keep = ~(masks[i][pos] | masks[j][pos])
                D[i, j] = D[j, i] = len(np.unique(pos[keep]))
    return D


def cluster_table(df, D, T, date_col="dec_date", country_col="country", min_span=0.5):
    """Add cluster id, size, span in years and class to a copy of df (row order = D order)."""
    _, lab = connected_components(csr_matrix(D <= T), directed=False)
    df = df.copy()
    df["cluster"] = lab
    g = df.groupby("cluster").agg(size=("cluster", "size"), first=(date_col, "min"), last=(date_col, "max"),
                                  countries=(country_col, "nunique"))
    g["span_years"] = (g["last"] - g["first"]).fillna(0)
    g["class"] = np.where(g["size"] <= 2, "sporadic",
                  np.where((g["size"] >= 10) & (g.span_years >= min_span), "sustained", "intermediate"))
    return df.merge(g[["size", "span_years", "countries", "class"]].rename(columns={"countries": "cluster_countries"}),
                    left_on="cluster", right_index=True)
