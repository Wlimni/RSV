"""
Script Name: run_sars2_human_compare.py
Purpose:     Within-lineage comparison, animal host versus human, for SARS-CoV-2: composition and the
             recent-mutation spectrum. This is the test notebook SARS2/01 flagged as missing: the
             deer-versus-baseline comparison used a mixed-lineage reference tree, not matched human
             genomes of the same lineage, country and period.
Inputs:      data/processed/Cov_19/Cov_19_features.tsv  (build_sars2_features.py, host column, search
             column identifies which GISAID search a genome came from)
             data/processed/Cov_19/nextclade/nextclade.tsv (same dataset as run_sars2_hosts.py)
Outputs:     results/sars2_within_lineage_composition.tsv
             results/sars2_within_lineage_spectrum.tsv
             results/sars2_within_lineage_mw.tsv        Mann-Whitney test per feature per pair
Notes:       Pairs: (deer, AY.103, USA) vs (human, AY.103, USA); (deer, B.1.1.7, USA) vs (human,
             B.1.1.7, USA); (deer, B.1.2, USA) vs (human, B.1.2, USA); (mink, B.1.1.298, Denmark) vs
             (human, B.1.1.298, Denmark); (mink, B.1.8, Netherlands) vs (human, B.1.8, Netherlands);
             (mink, B.1, USA) vs (human, B.1, USA). Spectrum uses recent (private) mutations relative
             to the nearest Nextclade tree node, genomes with <= 1% ambiguous bases and <= 30 recent
             mutations, same filters as the other mpox/SARS-CoV-2 spectrum scripts.
             Run: .venv/bin/python scripts/run_sars2_human_compare.py
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
from rsv import mpox

f = pd.read_csv("data/processed/Cov_19/Cov_19_features.tsv", sep="\t", dtype={"date": str})
f["hostg"] = f.host.map(lambda h: "mink" if h in ("Neovison vison", "Neogale vison")
                        else "deer" if h in ("Odocoileus virginianus", "White-tailed deer")
                        else "human" if h == "Human" else "cat")
nc = pd.read_csv("data/processed/Cov_19/nextclade/nextclade.tsv", sep="\t", low_memory=False)
nc["gisaid_epi_isl"] = nc.seqName
priv = lambda r: mpox.parse(r["privateNucMutations.labeledSubstitutions"]) + mpox.parse(r["privateNucMutations.unlabeledSubstitutions"])
nc["priv_list"] = nc.apply(priv, axis=1)
nc["priv"] = nc.priv_list.map(len)
d = f.merge(nc[["gisaid_epi_isl", "priv_list", "priv"]], on="gisaid_epi_isl")
d["good"] = (d.ambig_frac <= 0.01) & (d.priv <= 30)

PAIRS = [("deer", "AY.103", "USA"), ("deer", "B.1.1.7", "USA"), ("deer", "B.1.2", "USA"),
         ("mink", "B.1.1.298", "Denmark"), ("mink", "B.1.8", "Netherlands"), ("mink", "B.1", "USA")]
FEATS = ["gc_content", "oe_CG", "oe_TA", "oe_TC"]
types = [a + ">" + b for a in "ACGT" for b in "ACGT" if a != b]

comp_rows, mw_rows, spec_rows = [], [], []
for hostg, lin, country in PAIRS:
    a = d[(d.hostg == hostg) & (d.pangolin_lineage == lin) & (d.country == country)]
    h = d[(d.hostg == "human") & (d.pangolin_lineage == lin) & (d.country == country)]
    label = f"{hostg} vs human, {lin}, {country}"
    for feat in FEATS:
        comp_rows.append({"pair": label, "feature": feat, "n_animal": len(a), "n_human": len(h),
                          "mean_animal": a[feat].mean(), "mean_human": h[feat].mean()})
        if len(a) >= 5 and len(h) >= 5:
            stat, p = mannwhitneyu(a[feat].dropna(), h[feat].dropna())
            mw_rows.append({"pair": label, "feature": feat, "p": p, "n_animal": len(a), "n_human": len(h)})
    for grp_name, x in (("animal", a[a.good]), ("human", h[h.good])):
        cnt = pd.Series(0, index=types, dtype=float)
        for lst in x.priv_list:
            for r, p, al in lst:
                cnt[r + ">" + al] += 1
        n = cnt.sum()
        row = {"pair": label, "which": grp_name, "genomes": len(x), "recent_mutations": int(n)}
        row.update({k: (v / n if n else np.nan) for k, v in cnt.items()})
        spec_rows.append(row)

comp = pd.DataFrame(comp_rows)
comp.to_csv("results/sars2_within_lineage_composition.tsv", sep="\t", index=False)
mw = pd.DataFrame(mw_rows)
mw.to_csv("results/sars2_within_lineage_mw.tsv", sep="\t", index=False)
spec = pd.DataFrame(spec_rows)
spec.to_csv("results/sars2_within_lineage_spectrum.tsv", sep="\t", index=False)

print(comp.round(5).to_string())
print()
print(mw.round(4).to_string())
print()
print(spec[["pair", "which", "genomes", "recent_mutations", "C>T", "G>T", "A>G", "T>C", "G>A"]].round(3).to_string())

# ---- bootstrap CI on the pooled C>T share difference (genome-resampled; genomes within a pair
# share some ancestry so this likely understates the true interval, same caveat as the mpox analysis)
def boot_share(x, n=3000, seed=0):
    rng = np.random.default_rng(seed)
    c = np.array([sum(1 for r, p, a in lst if r == "C" and a == "T") for lst in x.priv_list])
    t = x.priv.values
    idx = rng.integers(0, len(x), (n, len(x)))
    s = c[idx].sum(1) / np.maximum(t[idx].sum(1), 1)
    return c.sum() / max(t.sum(), 1), np.percentile(s, 2.5), np.percentile(s, 97.5)


ci_rows = []
for hostg, lin, country in PAIRS:
    a = d[(d.hostg == hostg) & (d.pangolin_lineage == lin) & (d.country == country) & d.good]
    h = d[(d.hostg == "human") & (d.pangolin_lineage == lin) & (d.country == country) & d.good]
    sa, la, ha = boot_share(a)
    sh, lh, hh = boot_share(h)
    ci_rows.append({"pair": f"{hostg}, {lin}, {country}", "n_animal": len(a), "animal_C>T": sa,
                    "animal_ci_lo": la, "animal_ci_hi": ha, "n_human": len(h), "human_C>T": sh,
                    "human_ci_lo": lh, "human_ci_hi": hh})
ci = pd.DataFrame(ci_rows)
ci.to_csv("results/sars2_within_lineage_ct_ci.tsv", sep="\t", index=False)
print(ci.round(3).to_string())
