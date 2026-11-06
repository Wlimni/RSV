"""
Script Name: run_sars2_hosts.py
Purpose:     SARS-CoV-2 in animal hosts (mink, white-tailed deer, cats): how many separate transmission
             clusters each host has, how large and long-lived they are, and what the recent mutation
             spectrum looks like. Sporadic and sustained are defined from cluster structure.
Inputs:      data/processed/Cov_19/Cov_19_features.tsv  (build_sars2_features.py)
             data/processed/Cov_19/nextclade/nextclade.tsv, dataset/reference.fasta
Outputs:     data/processed/Cov_19/Cov_19_clusters.tsv (one row per genome, not pushed)
             results/sars2_host_clusters.tsv      one row per cluster
             results/sars2_host_summary.tsv       per host: genomes, clusters, class counts
             results/sars2_spectrum_by_host.tsv   recent-mutation spectrum by host and class
Notes:       Host groups: mink = Neovison vison + Neogale vison; deer = Odocoileus virginianus +
             White-tailed deer; cat = Felis catus. Clusters use pairwise substitution differences
             (rsv.clusters), threshold 3 (main) with 2 and 5 as sensitivity; clusters are built
             inside each host group, and across hosts as a second run to see shared lineages.
             Recent mutations = substitutions against the nearest Nextclade tree node; genomes
             need <= 1% ambiguous bases and <= 30 recent mutations for the spectrum.
             Run: .venv/bin/python scripts/run_sars2_hosts.py
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np, pandas as pd
from rsv import clusters, lineage, mpox

ref = mpox.load_reference("data/processed/Cov_19/nextclade/dataset/reference.fasta")
L = len(ref)
nc = pd.read_csv("data/processed/Cov_19/nextclade/nextclade.tsv", sep="\t", low_memory=False)
nc["gisaid_epi_isl"] = nc.seqName
f = pd.read_csv("data/processed/Cov_19/Cov_19_features.tsv", sep="\t", dtype={"date": str})
d = f.merge(nc[["gisaid_epi_isl", "clade", "substitutions", "missing", "deletions", "alignmentStart", "alignmentEnd",
                "privateNucMutations.labeledSubstitutions", "privateNucMutations.unlabeledSubstitutions"]],
            on="gisaid_epi_isl")
d["hostg"] = d.host.map(lambda h: "mink" if h in ("Neovison vison", "Neogale vison")
                        else "deer" if h in ("Odocoileus virginianus", "White-tailed deer") else "cat")
d["date"] = pd.to_datetime(d.date, errors="coerce")
d["dec_date"] = lineage.decimal_date(d["date"])
d = d.reset_index(drop=True)
print(len(d), "genomes;", d.hostg.value_counts().to_dict())

priv = lambda r: mpox.parse(r["privateNucMutations.labeledSubstitutions"]) + mpox.parse(r["privateNucMutations.unlabeledSubstitutions"])
d["priv_list"] = d.apply(priv, axis=1)
d["priv"] = d.priv_list.map(len)
d["good"] = (d.ambig_frac <= 0.01) & (d.priv <= 30)

codes = [clusters.encode(mpox.parse(s)) for s in d.substitutions]
masks = [clusters.ranges_to_mask([m, dl], s, e, L) for m, dl, s, e in
         zip(d.missing, d.deletions, d.alignmentStart, d.alignmentEnd)]

out_clusters, summary, spec_rows = [], [], []
for scope, groups in (("within host", [(h, d.index[d.hostg == h]) for h in ("mink", "deer", "cat")]),
                      ("all animals together", [("all", d.index)])):
    for name, idx in groups:
        sub = d.loc[idx].reset_index(drop=True)
        D = clusters.distances([codes[i] for i in idx], [masks[i] for i in idx])
        for T in (2, 3, 5):
            t = clusters.cluster_table(sub, D, T)
            cl = t.drop_duplicates("cluster")
            row = {"scope": scope, "host": name, "T": T, "genomes": len(t), "clusters": len(cl),
                   "singletons": int((cl["size"] == 1).sum()), "largest_cluster": int(cl["size"].max()),
                   "sporadic": int((cl["class"] == "sporadic").sum()), "intermediate": int((cl["class"] == "intermediate").sum()),
                   "sustained": int((cl["class"] == "sustained").sum()),
                   "median_span_years_clusters_2plus": cl[cl["size"] >= 2].span_years.median(),
                   "max_span_years": cl.span_years.max(), "genomes_in_sustained": int((t["class"] == "sustained").sum())}
            summary.append(row)
            print(row, flush=True)
            if T == 3 and scope == "within host":
                t["scope_host"] = name
                out_clusters.append(t)
sm = pd.DataFrame(summary)
sm.to_csv("results/sars2_host_summary.tsv", sep="\t", index=False)
oc = pd.concat(out_clusters)
cl_tab = (oc.groupby(["scope_host", "cluster"]).agg(size=("size", "first"), span_years=("span_years", "first"),
          class_=("class", "first"), countries=("country", "nunique"), lineages=("pangolin_lineage", lambda s: ",".join(sorted(set(s))[:4])),
          first=("dec_date", "min"), last=("dec_date", "max")).reset_index().rename(columns={"class_": "class"}))
cl_tab.to_csv("results/sars2_host_clusters.tsv", sep="\t", index=False)
oc.drop(columns=["priv_list"]).to_csv("data/processed/Cov_19/Cov_19_clusters.tsv", sep="\t", index=False)

# recent-mutation spectrum: share of each substitution type among recent mutations
types = [a + ">" + b for a in "ACGT" for b in "ACGT" if a != b]
rows = []
for (h, cls), x in oc[oc.good].groupby(["scope_host", "class"]):
    cnt = pd.Series(0, index=types, dtype=float)
    for lst in x.priv_list:
        for r, p, a in lst:
            cnt[r + ">" + a] += 1
    n = cnt.sum()
    row = {"host": h, "class": cls, "genomes": len(x), "clusters": x.cluster.nunique(), "recent mutations": int(n)}
    row.update({k: v / n for k, v in cnt.items()})
    rows.append(row)
pd.DataFrame(rows).to_csv("results/sars2_spectrum_by_host.tsv", sep="\t", index=False)
print(pd.DataFrame(rows)[["host", "class", "genomes", "clusters", "recent mutations", "C>T", "G>T", "A>G", "T>C"]].round(3).to_string())

# ---- extra tables: baseline spectrum from the Nextclade reference tree (mostly human), spectrum by
# submitting lab, and the spectrum of recurrent mutations (seen in 2+ genomes, so not sequencing noise)
import json, re
from collections import Counter
tree = json.load(open("data/processed/Cov_19/nextclade/dataset/tree.json"))
base = Counter()


def walk(n):
    for x in n.get("branch_attrs", {}).get("mutations", {}).get("nuc", []):
        m = re.match(r"^([ACGT])(\d+)([ACGT])$", x)
        if m:
            base[m.group(1) + ">" + m.group(3)] += 1
    for c in n.get("children", []):
        walk(c)


walk(tree["tree"])
tot = sum(base.values())
pd.DataFrame({"type": types, "share": [base[t] / tot for t in types]}).assign(mutations=tot).to_csv(
    "results/sars2_baseline_spectrum.tsv", sep="\t", index=False)

good = oc[oc.good].copy()
lab_rows = []
for (h, lab), x in good.groupby(["scope_host", "submitting_lab"]):
    c = Counter(r + ">" + a for lst in x.priv_list for r, p, a in lst)
    n = sum(c.values())
    if n >= 150:
        lab_rows.append({"host": h, "lab": lab[:60], "genomes": len(x), "recent mutations": n,
                         **{t: c[t] / n for t in ("C>T", "G>T", "A>G", "T>C", "G>A")}})
pd.DataFrame(lab_rows).to_csv("results/sars2_spectrum_by_lab.tsv", sep="\t", index=False)
rec_rows = []
for h, x in good.groupby("scope_host"):
    c = Counter((p, a, r) for lst in x.priv_list for r, p, a in lst)
    rec = Counter(r + ">" + a for (p, a, r), v in c.items() if v >= 2)
    n = sum(rec.values())
    rec_rows.append({"host": h, "recurrent mutations": n, **{t: rec[t] / n for t in ("C>T", "G>T", "A>G", "T>C", "G>A")}})
pd.DataFrame(rec_rows).to_csv("results/sars2_spectrum_recurrent.tsv", sep="\t", index=False)
print(pd.DataFrame(rec_rows).round(3).to_string())
