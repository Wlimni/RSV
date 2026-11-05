"""
Script Name: run_mpox_lineage_analysis.py
Purpose:     Mpox analysis inside lineages. Joins Nextclade lineage calls and substitutions to the
             genome composition table, counts APOBEC3-type mutations per genome, and fits time
             trends inside each lineage with cluster-robust errors (same method as RSV).
Inputs:      data/processed/Mpox/Mpox_features.tsv (build_mpox_features.py)
             data/processed/Mpox/nextclade/nextclade.tsv and dataset/reference.fasta
             (nextclade run --input-dataset nextstrain/mpox/all-clades)
Outputs:     data/processed/Mpox/Mpox_lineage_features.tsv   one row per genome (not pushed)
             results/mpox_lineage_slopes.tsv                 slope per lineage per feature
             results/mpox_spectrum_by_lineage.tsv            APOBEC3 share of recent mutations
Notes:       APOBEC3-type mutation = C->T where the reference base before it is T (TC -> TT), or
             G->A where the reference base after it is A (GA -> AA, the same event on the other
             strand). "GC loss" counts C->T and G->A minus the reverse changes, so it is the
             change in G+C count that the mutations imply.
             Total substitutions are counted against the Nextclade reference (a 2017 clade IIb
             isolate), "private" mutations against the nearest node of the Nextclade tree, which
             makes them recent mutations and comparable across lineages.
             Dates: full dates as given; YYYY-MM becomes the 15th of that month; year-only
             dates are dropped.
             Run: .venv/bin/python scripts/run_mpox_lineage_analysis.py
"""
import os, sys, re
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np, pandas as pd
from rsv import lineage

ref = "".join(l.strip() for l in open("data/processed/Mpox/nextclade/dataset/reference.fasta") if not l.startswith(">")).upper()
SUB = re.compile(r"^([ACGT])(\d+)([ACGT])$")


def parse(cell):
    """'C123T,G456A' -> list of (ref, pos, alt) for plain A/C/G/T substitutions."""
    if not isinstance(cell, str):
        return []
    out = []
    for x in cell.split(","):
        m = SUB.match(x)
        if m:
            out.append((m.group(1), int(m.group(2)), m.group(3)))
    return out


def classify(subs):
    """Counts: total, APOBEC3-type, and net G+C loss."""
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


nc = pd.read_csv("data/processed/Mpox/nextclade/nextclade.tsv", sep="\t", low_memory=False)
nc["Accession ID"] = nc.seqName.str.split("|").str[1]
rows = []
for _, r in nc.iterrows():
    n, a, g = classify(parse(r["substitutions"]))
    priv = parse(r["privateNucMutations.labeledSubstitutions"]) + parse(r["privateNucMutations.unlabeledSubstitutions"])
    pn, pa, pg = classify(priv)
    rows.append({"Accession ID": r["Accession ID"], "nc_lineage": r["lineage"], "nc_clade": r["clade"],
                 "n_subs": n, "n_apobec": a, "gc_loss": g,
                 "priv_subs": pn, "priv_apobec": pa, "priv_gc_loss": pg})
sub = pd.DataFrame(rows)

f = pd.read_csv("data/processed/Mpox/Mpox_features.tsv", sep="\t", dtype={"Collection date": str})
d = f.merge(sub, on="Accession ID")
cd = d["Collection date"]
d["date"] = pd.to_datetime(cd.where(cd.str.match(r"^\d{4}-\d{2}-\d{2}$"))
                           .fillna(cd.where(cd.str.match(r"^\d{4}-\d{2}$")).map(lambda s: f"{s}-15" if isinstance(s, str) else s)),
                           errors="coerce")
d["country"] = d.Location.str.split(" / ").str[1]
d["Subtype"] = "MPXV"
d["Lineage"] = d.nc_lineage.fillna("NA")
d["dec_date"] = lineage.decimal_date(d["date"])
d.to_csv("data/processed/Mpox/Mpox_lineage_features.tsv", sep="\t", index=False)
print(len(d), "genomes;", d.date.notna().sum(), "with usable date")

# window: through end of May 2026 (later months are still being deposited); human genomes only
h = d[(d.Host == "Human") & d.date.notna()].copy()
FEATS = ["gc_content", "oe_TC", "oe_CG", "n_apobec", "n_subs", "gc_loss"]
res_all = []
for level, min_years in ((2, 3),):
    p, kept = lineage.prepare(h, level=level, min_years=min_years, min_n=200,
                              start=2016.0, end=2026 + 151 / 365, trim=5)
    print(f"level {level}: lineages kept\\n", kept.to_string())
    r = lineage.lineage_slopes(p, FEATS, by=("Subtype", "lineage_grp"), cluster="clade_fine", min_clusters=3)
    r["level"] = level
    res_all.append(r)
res = pd.concat(res_all)
res.to_csv("results/mpox_lineage_slopes.tsv", sep="\t", index=False)
print(res[["lineage_grp", "feature", "n", "slope", "ci_lo", "ci_hi", "p", "n_clusters", "cluster_by"]].round(6).to_string())

# recent-mutation spectrum by lineage family and year
h["family"] = h.Lineage.str.split(".").str[0]
h["year"] = h.date.dt.year
spec = (h[h.priv_subs > 0].groupby(["family", "year"])
        .agg(genomes=("priv_subs", "size"), priv_subs=("priv_subs", "sum"), priv_apobec=("priv_apobec", "sum"))
        .assign(apobec_share=lambda x: x.priv_apobec / x.priv_subs).reset_index())
spec.to_csv("results/mpox_spectrum_by_lineage.tsv", sep="\t", index=False)
print(spec.round(3).to_string())
