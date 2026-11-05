"""
Script Name: run_mpoxi_analysis.py
Purpose:     Clade I mpox (Ia and Ib) pipeline: APOBEC3 share of recent mutations by group, cluster
             size as an independent sporadic-versus-sustained marker, the mutation clock inside
             the Ib outbreak, and composition by group. Clade IIb is run alongside for comparison.
Inputs:      data/processed/MpoxI/MpoxI_features.tsv   (build_mpox_features.py MpoxI)
             data/processed/MpoxI/nextclade/nextclade.tsv (nextclade run, all-clades dataset)
             data/processed/Mpox/Mpox_lineage_features.tsv (clade IIb, run_mpox_lineage_analysis.py)
Outputs:     data/processed/MpoxI/MpoxI_analysis.tsv   one row per genome (not pushed)
             results/mpoxi_group_summary.tsv, results/mpoxi_apobec_share.tsv,
             results/mpoxi_cluster_share.tsv, results/mpoxi_clock.tsv
Notes:       Groups come from Nextclade: "Ib sh2023" and "Ia sh2024" are the outbreaks it names
             (sustained human-to-human spread); "Ia no outbreak" is everything else in clade Ia.
             Those outbreak labels are assigned by the lineage nomenclature, which is about
             sustained human transmission, so they are not independent of the APOBEC3 signature.
             The cluster-size marker (how many other genomes lie within a few substitutions) does
             not use the mutation spectrum and is the independent check.
             Quality filters for the spectrum: at most 1% ambiguous bases, and at most 30 recent
             mutations (more than that points to sequencing error or an unusual genome).
             Run: .venv/bin/python scripts/run_mpoxi_analysis.py
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np, pandas as pd
from rsv import lineage, mpox

ref = mpox.load_reference("data/processed/Mpox/nextclade/dataset/reference.fasta")
nc = pd.read_csv("data/processed/MpoxI/nextclade/nextclade.tsv", sep="\t", low_memory=False)
nc["Accession ID"] = nc.seqName.str.split("|").str[1]
keep = ["Accession ID", "clade", "outbreak", "lineage", "substitutions",
        "privateNucMutations.labeledSubstitutions", "privateNucMutations.unlabeledSubstitutions",
        "founderMuts['outbreak'].substitutions"]
f = pd.read_csv("data/processed/MpoxI/MpoxI_features.tsv", sep="\t", dtype={"Collection date": str})
d = f.merge(nc[keep], on="Accession ID")
d["country"] = d.Location.str.split(" / ").str[1]
cd = d["Collection date"]
d["date"] = pd.to_datetime(cd.where(cd.str.match(r"^\d{4}-\d{2}-\d{2}$"))
                           .fillna(cd.where(cd.str.match(r"^\d{4}-\d{2}$")).map(lambda s: f"{s}-15" if isinstance(s, str) else s)),
                           errors="coerce")
d["year"] = pd.to_numeric(cd.str[:4], errors="coerce")
d["dec_date"] = lineage.decimal_date(d["date"])
d["group"] = np.where(d.outbreak == "sh2023", "Ib sh2023",
              np.where(d.outbreak == "sh2024", "Ia sh2024", np.where(d.clade == "Ia", "Ia no outbreak", d.clade)))
d = d[d.group != "IIb"]                      # one stray IIb genome in the Ib download

# ---- recent-mutation spectrum
priv = lambda r: mpox.parse(r["privateNucMutations.labeledSubstitutions"]) + mpox.parse(r["privateNucMutations.unlabeledSubstitutions"])
d[["priv", "priv_apobec", "priv_loss"]] = d.apply(lambda r: mpox.classify(priv(r), ref), axis=1, result_type="expand")
d["good"] = (d.ambig_frac <= 0.01) & (d.priv <= 30)

# ---- clade IIb for comparison, same filters
b = pd.read_csv("data/processed/Mpox/Mpox_lineage_features.tsv", sep="\t", low_memory=False)
b = b[(b.Host == "Human") & (b.ambig_frac <= 0.01) & (b.priv_subs <= 30) & b.nc_lineage.notna()].copy()
b["group"] = np.where(b.nc_lineage.str.startswith("A"), "IIb lineage A (pre-2022)", "IIb B-J (2022 outbreak)")
b = b.rename(columns={"priv_subs": "priv", "priv_apobec": "priv_apobec"})


def boot_share(x, n=2000, seed=0):
    """Pooled APOBEC3 share with a genome-level bootstrap interval (optimistic: genomes are related)."""
    rng = np.random.default_rng(seed)
    a, t = x.priv_apobec.values, x.priv.values
    idx = rng.integers(0, len(x), (n, len(x)))
    s = a[idx].sum(1) / np.maximum(t[idx].sum(1), 1)
    return a.sum() / max(t.sum(), 1), np.percentile(s, 2.5), np.percentile(s, 97.5)


rows = []
for name, x in list(d[d.good].groupby("group")) + list(b.groupby("group")):
    s, lo, hi = boot_share(x)
    per = (x.priv_apobec / x.priv)[x.priv >= 3]
    rows.append({"group": name, "genomes": len(x), "recent mutations": int(x.priv.sum()),
                 "APOBEC3 share": s, "ci_lo": lo, "ci_hi": hi, "median per-genome share (>=3 mutations)": per.median()})
share = pd.DataFrame(rows)
share.to_csv("results/mpoxi_apobec_share.tsv", sep="\t", index=False)
print(share.round(3).to_string())

# ---- era inside clade Ia no outbreak: old (zoonotic-era) versus 2022 on
ia = d[(d.group == "Ia no outbreak") & d.good]
era = ia.assign(era=np.where(ia.year <= 2018, "1970-2018", "2021-2026")).groupby("era").apply(
    lambda x: pd.Series({"genomes": len(x), "recent mutations": int(x.priv.sum()),
                         "APOBEC3 share": x.priv_apobec.sum() / max(x.priv.sum(), 1)}), include_groups=False)
print(era.round(3))

# ---- cluster size: other genomes within <= T substitutions (spectrum-independent)
cl = d[d.ambig_frac <= 0.01].copy()
sets = [frozenset((p, a) for r, p, a in mpox.parse(s)) for s in cl.substitutions]
n = len(sets)
dist = np.zeros((n, n), np.int32)
for i in range(n):
    for j in range(i + 1, n):
        dist[i, j] = dist[j, i] = len(sets[i] ^ sets[j])
for T in (3, 5):
    cl[f"nbr{T}"] = (dist <= T).sum(1) - 1
cl["cluster_class"] = pd.cut(cl.nbr5, [-1, 0, 4, 10**6], labels=["singleton (0 within 5 SNPs)", "small (1-4)", "large (5+)"])
cs = (cl[cl.priv <= 30].groupby("cluster_class", observed=True)
      .apply(lambda x: pd.Series({"genomes": len(x), "recent mutations": int(x.priv.sum()),
                                  "APOBEC3 share": x.priv_apobec.sum() / max(x.priv.sum(), 1),
                                  "share Ib sh2023": (x.group == "Ib sh2023").mean()}), include_groups=False))
cs.to_csv("results/mpoxi_cluster_share.tsv", sep="\t")
print(cs.round(3))
print(pd.crosstab(cl.group, cl.cluster_class))

# ---- clock inside the Ib outbreak: APOBEC3-type changes since the outbreak founder
d = d.merge(cl[["Accession ID", "nbr3", "nbr5", "cluster_class"]], on="Accession ID", how="left")
fo = d["founderMuts['outbreak'].substitutions"].map(mpox.parse)
d[["fm_subs", "fm_apobec", "fm_loss"]] = [mpox.classify(x, ref) for x in fo]
d["Subtype"] = "MPXV-I"
d["Lineage"] = d.lineage.fillna(d.country).fillna("NA")
out = []
for grp in ("Ib sh2023", "Ia sh2024"):
    g = d[(d.group == grp) & (d.ambig_frac <= 0.01) & d.date.notna()].copy()
    g["lineage_grp"] = grp
    g["clade_fine"] = g.Lineage
    r = lineage.lineage_slopes(g, ["fm_apobec", "fm_subs", "gc_content", "oe_TC"], by=("Subtype", "lineage_grp"),
                               cluster="clade_fine", min_clusters=3)
    out.append(r)
clock = pd.concat(out)
clock.to_csv("results/mpoxi_clock.tsv", sep="\t", index=False)
print(clock[["lineage_grp", "feature", "n", "slope", "ci_lo", "ci_hi", "p", "n_clusters", "cluster_by"]].round(6).to_string())

# ---- composition by group
comp = pd.concat([d[["group", "gc_content", "oe_CG", "oe_TA", "oe_TC", "year", "country"]],
                  b[["group", "gc_content", "oe_CG", "oe_TA", "oe_TC"]]]).groupby("group").agg(
    genomes=("gc_content", "size"), gc=("gc_content", "mean"), cpg_oe=("oe_CG", "mean"),
    upa_oe=("oe_TA", "mean"), tpc_oe=("oe_TC", "mean"))
comp.to_csv("results/mpoxi_group_summary.tsv", sep="\t")
print(comp.round(4))
d.to_csv("data/processed/MpoxI/MpoxI_analysis.tsv", sep="\t", index=False)
