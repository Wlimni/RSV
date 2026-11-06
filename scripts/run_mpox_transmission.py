"""
Script Name: run_mpox_transmission.py
Purpose:     Sporadic versus sustained transmission in mpox, defined from the genomes themselves and
             not from the mutation spectrum, then compared on the APOBEC3 share of recent mutations.
             Clusters of near-identical genomes are built from pairwise substitution differences;
             a cluster's size and time span decide whether it looks sporadic or sustained.
Inputs:      data/processed/MpoxI/MpoxI_analysis.tsv    (run_mpoxi_analysis.py)
             data/processed/Mpox/Mpox_lineage_features.tsv (run_mpox_lineage_analysis.py)
             data/processed/{MpoxI,Mpox}/nextclade/nextclade.tsv, dataset/reference.fasta
Outputs:     results/mpox_transmission_clusters.tsv   one row per cluster (T = 5)
             results/mpox_transmission_share.tsv      APOBEC3 share by class, clade, threshold
             results/mpox_transmission_sensitivity.tsv
Notes:       Pairwise distance = number of substitution differences (against the Nextclade
             reference) at positions readable in both genomes, so gaps do not look like
             differences. Clusters = connected components of "distance <= T" (T = 3, 5, 10).
             Classes: sporadic = cluster of 1-2 genomes; intermediate = 3-9; sustained = 10 or
             more genomes spanning at least 0.5 year. Genomes with no usable date count toward
             size but not span.
             Clade IIb is a random subsample (all lineage-A genomes plus 1,300 others) because
             the pairwise step grows with the square of the sample.
             APOBEC3 share is the pooled share of "recent" mutations (against the nearest
             Nextclade tree node), genomes with <= 1% ambiguous bases and <= 30 recent mutations.
             Intervals resample whole clusters, not genomes, because genomes in a cluster are
             related.
             Run: .venv/bin/python scripts/run_mpox_transmission.py
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np, pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from rsv import mpox

ref = mpox.load_reference("data/processed/Mpox/nextclade/dataset/reference.fasta")
L = len(ref)
ALT = {"A": 0, "C": 1, "G": 2, "T": 3}


def ranges_to_mask(cells, start, end):
    """Boolean array (True = unreadable): missing ranges, deletions, and sites outside the alignment."""
    m = np.zeros(L + 2, bool)
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


def load(df, ncpath):
    nc = pd.read_csv(ncpath, sep="\t", low_memory=False)
    nc["Accession ID"] = nc.seqName.str.split("|").str[1]
    nc = nc.set_index("Accession ID")
    codes, masks = [], []
    for acc in df["Accession ID"]:
        r = nc.loc[acc]
        codes.append(encode(mpox.parse(r["substitutions"])))
        masks.append(ranges_to_mask([r["missing"], r["deletions"]], r["alignmentStart"], r["alignmentEnd"]))
    return codes, masks


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


def cluster_table(df, D, T):
    n_comp, lab = connected_components(csr_matrix(D <= T), directed=False)
    df = df.copy()
    df["cluster"] = lab
    g = df.groupby("cluster").agg(size=("cluster", "size"), first=("dec_date", "min"), last=("dec_date", "max"),
                                  countries=("country", "nunique"))
    g["span_years"] = (g["last"] - g["first"]).fillna(0)
    g["class"] = np.where(g["size"] <= 2, "sporadic",
                  np.where((g["size"] >= 10) & (g.span_years >= 0.5), "sustained", "intermediate"))
    return df.merge(g[["size", "span_years", "class"]], left_on="cluster", right_index=True)


def boot(x, n=2000, seed=0):
    """Pooled APOBEC3 share; resample whole clusters."""
    rng = np.random.default_rng(seed)
    cl = x.groupby("cluster").agg(a=("priv_apobec", "sum"), t=("priv", "sum"))
    a, t = cl.a.values, cl.t.values
    idx = rng.integers(0, len(cl), (n, len(cl)))
    s = a[idx].sum(1) / np.maximum(t[idx].sum(1), 1)
    return a.sum() / max(t.sum(), 1), np.percentile(s, 2.5), np.percentile(s, 97.5), len(cl)


# ---------------------------------------------------------------- datasets
c1 = pd.read_csv("data/processed/MpoxI/MpoxI_analysis.tsv", sep="\t", low_memory=False)
c1 = c1[(c1.ambig_frac <= 0.05) & c1.group.isin(["Ia no outbreak", "Ia sh2024", "Ib sh2023"])].copy()
c1["dataset"] = "clade I"
c1["priv_apobec"] = c1.priv_apobec

b = pd.read_csv("data/processed/Mpox/Mpox_lineage_features.tsv", sep="\t", low_memory=False)
b = b[(b.Host == "Human") & (b.ambig_frac <= 0.05) & b.nc_lineage.notna()].copy()
fam_a = b[b.nc_lineage.str.startswith("A")]
rest = b[~b.nc_lineage.str.startswith("A")].sample(1300, random_state=0)
b = pd.concat([fam_a, rest]).copy()
b = b.rename(columns={"priv_subs": "priv"})
b["group"] = np.where(b.nc_lineage.str.startswith("A"), "IIb lineage A", "IIb B-J")
b["dataset"] = "clade IIb"
b["dec_date"] = b.dec_date

sets = {"clade I": (c1, "data/processed/MpoxI/nextclade/nextclade.tsv"),
        "clade IIb": (b, "data/processed/Mpox/nextclade/nextclade.tsv")}
all_rows, cluster_rows, sens = [], [], []
for name, (df, path) in sets.items():
    df = df.reset_index(drop=True)
    print(name, len(df), "genomes", flush=True)
    codes, masks = load(df, path)
    D = distances(codes, masks)
    for T in (3, 5, 10):
        t = cluster_table(df, D, T)
        print(f"  T={T}: clusters={t.cluster.nunique()}", t.drop_duplicates('cluster')["class"].value_counts().to_dict(), flush=True)
        for filt_name, filt in (("<=1% N, <=30 recent", (t.ambig_frac <= 0.01) & (t.priv <= 30)),
                                ("<=5% N, <=30 recent", (t.ambig_frac <= 0.05) & (t.priv <= 30)),
                                ("<=1% N, <=15 recent", (t.ambig_frac <= 0.01) & (t.priv <= 15))):
            for cls, x in t[filt].groupby("class"):
                s, lo, hi, ncl = boot(x)
                sens.append({"dataset": name, "T": T, "filter": filt_name, "class": cls, "genomes": len(x),
                             "clusters": ncl, "recent mutations": int(x.priv.sum()), "APOBEC3 share": s, "ci_lo": lo, "ci_hi": hi})
        if T == 5:
            t["dataset"] = name
            all_rows.append(t[["Accession ID", "dataset", "group", "cluster", "size", "span_years", "class",
                               "priv", "priv_apobec", "ambig_frac", "country", "year", "date"]
                              if "year" in t else ["Accession ID", "dataset", "group", "cluster", "size", "span_years", "class",
                                                    "priv", "priv_apobec", "ambig_frac", "country"]])
            cl = (t.groupby("cluster").agg(dataset=("dataset", "first"), size=("size", "first"), span_years=("span_years", "first"),
                                           klass=("class", "first"), countries=("country", "nunique"),
                                           groups=("group", lambda s: ",".join(sorted(set(s))))))
            cluster_rows.append(cl.reset_index().rename(columns={"klass": "class"}))
sens = pd.DataFrame(sens)
sens.to_csv("results/mpox_transmission_sensitivity.tsv", sep="\t", index=False)
main = sens[(sens["T"] == 5) & (sens["filter"] == "<=1% N, <=30 recent")]
main.to_csv("results/mpox_transmission_share.tsv", sep="\t", index=False)
pd.concat(cluster_rows).to_csv("results/mpox_transmission_clusters.tsv", sep="\t", index=False)
pd.concat(all_rows).to_csv("data/processed/Mpox/Mpox_transmission_genomes.tsv", sep="\t", index=False)
print(main.round(3).to_string())
