"""
Script Name: make_notebooks_nipah_mpox.py
Purpose:     Generate the first-look notebooks for Nipah virus (NiV) and mpox (MPXV), one folder each
             under notebooks/. These test whether the RSV pipeline is worth running on the other
             GISAID viruses and what the host field offers.
Inputs:      none (the notebooks read data/processed/{NiV,Mpox}/ built by merge_gisaid.py,
             build_niv_features.py and build_mpox_features.py)
Outputs:     notebooks/NiV/01_host_composition.ipynb, notebooks/Mpox/01_genome_composition.ipynb
Notes:       Then execute: jupyter nbconvert --execute --inplace notebooks/NiV/*.ipynb notebooks/Mpox/*.ipynb
"""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell as md, new_code_cell as code

SETUP = '''\
import sys, os
sys.path.insert(0, os.path.abspath("../../src"))
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

os.chdir(os.path.abspath("../.."))          # repo root, so data/ and results/ paths work
pd.set_option("display.width", 130)
plt.rcParams.update({"figure.dpi": 110, "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False})
'''


def save(path, cells):
    nb = new_notebook()
    nb.cells = cells
    nb.metadata["kernelspec"] = {"name": "rsv", "display_name": "rsv (.venv)", "language": "python"}
    nbf.write(nb, path)
    print("wrote", path)


# ------------------------------------------------------------------------------- Nipah
save("notebooks/NiV/01_host_composition.ipynb", [
md("""# Nipah virus (NiV) — does genome composition differ by host?

---

### Why this virus

RSV has no host variation (all 33,040 genomes are human). Nipah jumps from **fruit bats** (the natural host) to **pigs** and **humans**, so the GISAID host field carries real information. The question here: **does a Nipah genome look different when it comes from a bat, a pig or a human?** If the paper's idea (composition shifts with the host the virus lives in) held for Nipah, bat and human genomes would differ in GC or CpG.

### The data

GISAID EpiNiV export: 343 sequences. Merged by `scripts/merge_gisaid.py`; genes cut by `scripts/build_niv_features.py` against the RefSeq Malaysia strain NC_002728.1 (public NCBI sequence)."""),
code(SETUP),

md("""## 1. Host, clade and completeness"""),
code('''m = pd.read_csv("data/processed/NiV/NiV_metadata.tsv", sep="\\t", dtype=str)
def host_group(h):
    if h == "Human": return "human"
    if "Sus" in h or "Canis" in h: return "pig/dog"
    if "Pteropus" in h or "Chiroptera" in h: return "bat"
    return "unknown"
m["hostg"] = m.Host.map(host_group)
m["L"] = pd.to_numeric(m["Sequence Length"])
m["year"] = pd.to_numeric(m["Collection date"].str[:4], errors="coerce")
m["country"] = m.Location.str.split(" / ").str[1]
m["complete"] = np.where(m.L >= 17000, "full genome (>=17 kb)", "fragment")
print(m.Host.value_counts().to_string())
print()
print(pd.crosstab(m.hostg, m.complete, margins=True))'''),
md("""**Most bat sequences are short fragments.** Only 19 of the 216 bat genomes are full-length, against 74 of 110 from humans. That matters: a fragment covers one gene region, so it can only be compared with other sequences over the **same** region."""),
code('''print(pd.crosstab(m.country, m.hostg, margins=True))
print()
print("collection year known for", m.year.notna().sum(), "of", len(m), "(missing mostly bat sequences)")
print(m.groupby("hostg").year.agg(["min", "max", "count"]))'''),

md("""## 2. Host and timing: lineages BD1, BD2, BD3

GISAID labels the Bangladesh/India clade (NiV-B) as BD1, BD2 and BD3. This is the "host plus timing" information that could separate sporadic from sustained circulation."""),
code('''m["lineage"] = m.Lineage.fillna("none")
t = m.groupby("lineage").apply(lambda g: pd.Series({
    "genomes": len(g), "from bats": (g.hostg == "bat").sum(), "from humans": (g.hostg == "human").sum(),
    "first year": g.year.min(), "last year": g.year.max()}), include_groups=False)
t'''),
code('''fig, ax = plt.subplots(figsize=(8, 3))
d = m[m.year.notna() & m.lineage.isin(["BD1", "BD2", "BD3"]) & m.hostg.isin(["bat", "human"])]
cols = {"bat": "#6A3D9A", "human": "#E08214"}
for i, lin in enumerate(["BD1", "BD2", "BD3"]):
    for h, off in (("bat", .12), ("human", -.12)):
        y = d[(d.lineage == lin) & (d.hostg == h)].year
        ax.scatter(y + np.random.default_rng(0).uniform(-.25, .25, len(y)), np.full(len(y), i + off), s=14, color=cols[h], alpha=.7,
                   label=h if i == 0 else None)
ax.set_yticks(range(3)); ax.set_yticklabels(["BD1", "BD2", "BD3"]); ax.set_xlabel("collection year")
ax.legend(frameon=False, ncol=2, loc="upper left"); plt.tight_layout()
fig.savefig("results/figures/niv_lineage_host_timeline.png", dpi=200, bbox_inches="tight"); plt.show()'''),
md("""**Reading it.** BD2 runs from 2004 to 2025 and appears in **both bats and humans** across that whole span; BD1 and BD3 are the same pattern over shorter periods. That is what a persistent bat reservoir with repeated human spillover looks like: the lineage is sustained in bats, and each human outbreak is a short dead-end. It is **consistent with** that picture, not proof: 73 of the 216 bat sequences have no collection year, and the human genomes come from only a few outbreaks."""),

md("""## 3. Composition by host

For each gene, genomes are kept only if they cover at least 90% of the reference gene (so every sequence is compared over the same region). GC3 additionally needs a gene that starts with ATG and reads in frame."""),
code('''d = pd.read_csv("data/processed/NiV/NiV_gene_features.tsv", sep="\\t")
d["hostg"] = d.Host.map(lambda h: host_group(h) if isinstance(h, str) else "unknown")
d["clade"] = d.Clade.str.extract(r"\\((.*)\\)")
order = ["N", "M", "F", "G", "L", "GENOME"]
print(pd.crosstab(d.gene, [d.clade, d.hostg]).reindex(order))'''),
code('''# NiV-B only, where bats and humans can be compared like for like
b = d[d.clade == "NiV-B"]
tab = b.groupby(["gene", "hostg"])[["gc_content", "oe_CG", "oe_TA", "gc3"]].mean().round(4)
cmp = tab.unstack("hostg")
out = pd.DataFrame({
    "GC bat": cmp[("gc_content", "bat")], "GC human": cmp[("gc_content", "human")],
    "CpG O/E bat": cmp[("oe_CG", "bat")], "CpG O/E human": cmp[("oe_CG", "human")],
    "GC3 bat": cmp[("gc3", "bat")], "GC3 human": cmp[("gc3", "human")]}).reindex(order)
out["GC diff"] = (out["GC bat"] - out["GC human"]).round(4)
out["CpG diff"] = (out["CpG O/E bat"] - out["CpG O/E human"]).round(4)
out'''),
code('''fig, axs = plt.subplots(1, 2, figsize=(10, 3.4), sharex=True)
cols = {"bat": "#6A3D9A", "human": "#E08214", "pig/dog": "#1B9E77"}
for ax, col, lab in zip(axs, ["gc_content", "oe_CG"], ["GC content", "CpG O/E"]):
    for i, g in enumerate(order):
        for h, off in (("bat", -.2), ("human", 0), ("pig/dog", .2)):
            v = d[(d.gene == g) & (d.hostg == h)][col]
            if len(v):
                ax.scatter(np.full(len(v), i + off) + np.random.default_rng(1).uniform(-.04, .04, len(v)), v, s=10, color=cols[h], alpha=.6,
                           label=h if i == 0 else None)
    ax.set_xticks(range(len(order))); ax.set_xticklabels(order); ax.set_ylabel(lab)
axs[0].legend(frameon=False, fontsize=7)
plt.tight_layout(); fig.savefig("results/figures/niv_host_composition.png", dpi=200, bbox_inches="tight"); plt.show()'''),
md("""**Reading it.** Bat and human genomes are **indistinguishable in composition**: whole-genome GC is 0.386 in bats and 0.385 in humans, CpG O/E 0.450 and 0.441. The largest differences are 0.0005 in GC and 0.014 in CpG O/E, far smaller than the gaps between genes (e.g. N is 45% GC, L is 38%).

This is the expected result for a virus that is the **same lineage** in bats and in the humans they infect: a human outbreak starts from a bat strain, and it ends before any composition change can build up. Nipah's overall CpG O/E (about 0.44) is moderately suppressed: far less than RSV-A's NS1 (0.016), similar to most RSV genes (0.2-0.4)."""),

md("""## 4. What this means

| Question | Answer |
|---|---|
| Does composition differ by host? | **No.** Bat and human GC differ by at most 0.0005 and CpG O/E by at most 0.014, per gene. |
| Is there host and timing information? | **Yes, but thin.** BD2 spans 2004-2025 in both bats and humans, which fits a persistent bat reservoir with repeated spillover. |
| Can the RSV trend pipeline run? | **No.** There are 106 full genomes, and no lineage reaches the 200 genomes over 5 years that the trend analysis needs. |

**Limits.** 343 sequences, most of them fragments; 82 have no collection year; human genomes come from a handful of outbreaks (Bangladesh alone gives 88), so they are near-copies of each other, and counting them as independent overstates the evidence. The result above is therefore a **negative**, not a failed analysis: Nipah does not show a host-linked composition signal in this data."""),
])

# ------------------------------------------------------------------------------- Mpox
save("notebooks/Mpox/01_genome_composition.ipynb", [
md("""# Mpox (MPXV) — host, composition and time

---

### Why this virus

Mpox was suggested as another virus with an animal-to-human story. It is a **double-stranded DNA virus** with a ~197 kb genome (RSV and Nipah are RNA viruses of 15-18 kb), and human-adapted mpox mutates mainly through the host **APOBEC3** enzymes, which turn C into T (usually after a T: TC to TT). That is a different mechanism from the CpG pressure studied in RSV, so the questions here are:

1. Does the host field offer a real contrast?
2. Is mpox CpG-suppressed like the RNA viruses?
3. Does the APOBEC3 signature show up as a falling GC or TpC over time?

### The data

GISAID EpiPox export: 8,396 sequences, merged by `scripts/merge_gisaid.py`; whole-genome features by `scripts/build_mpox_features.py`. No gene-level analysis: mpox has about 180 genes."""),
code(SETUP),

md("""## 1. Host and sequence quality"""),
code('''f = pd.read_csv("data/processed/Mpox/Mpox_features.tsv", sep="\\t", dtype={"Collection date": str})
def host_group(h):
    if h == "Human": return "human"
    if "Environment" in h: return "environment (wastewater)"
    if h == "unknown": return "unknown"
    return "animal"
f["hostg"] = f.Host.map(host_group)
f["qc_pass"] = f.length.between(190000, 200000) & (f.ambig_frac <= 0.01)
print(f.Host.value_counts().to_string())
print()
print(f.groupby("hostg").agg(genomes=("qc_pass", "size"), pass_qc=("qc_pass", "sum")))'''),
md("""**This settles the host question.** The "non-human" entries are 220 **wastewater** samples (not an animal), 292 **unknown**, and only 40 real animals (29 macaques, 8 unspecified monkeys, 2 chimpanzees, 1 dog). **None of the 260 wastewater or animal sequences is a full genome** (they are short fragments), so no host comparison can be made on whole-genome data. Only human (7,527) and unknown-host (184) genomes pass quality control (190-200 kb, at most 1% ambiguous bases)."""),

md("""## 2. Is mpox CpG-suppressed?"""),
code('''q = f[f.qc_pass].copy()
rows = {"Mpox (this data)": [q.gc_content.mean(), q.oe_CG.mean(), q.oe_TA.mean(), q.oe_TC.mean()]}
pd.DataFrame(rows, index=["GC content", "CpG O/E", "UpA (TpA) O/E", "TpC O/E"]).T.round(3)'''),
md("""**No.** CpG O/E is **1.04**, i.e. CpG occurs as often as chance predicts. The RNA viruses we studied are all well below 1 (RSV genes mostly 0.2-0.4, Nipah about 0.44). CpG depletion in RNA viruses is usually read as the footprint of host pressure on viral RNA (ZAP). Why a large DNA virus like mpox lacks it is **not something this data can answer** (poxviruses replicate in the cytoplasm and make their own mRNAs, so the simple explanation does not hold). What we can say is that the RSV NS1/NS2 question does not transfer to mpox, because there is no CpG gap to explain.

TpC O/E is **1.13**: the APOBEC3 target dinucleotide is slightly over-represented, not depleted."""),

md("""## 3. APOBEC3 signature over time

APOBEC3 changes C to T (and G to A on the other strand). Each change removes one G or C, so GC content should fall, and TC should become TT, so TpC O/E should fall. The size of the effect is a few mutations per year in a 197 kb genome, i.e. a GC change of about **-0.000025 to -0.00005 per year** (one lost G/C is 1/197,000 = 5e-6). We test this on human genomes collected from 2022 onwards (the 2022 outbreak onwards), with errors clustered by country, and adjusting for genome length and ambiguous bases because the covered region changes composition slightly."""),
code('''from rsv import lineage
h = q[q.Host == "Human"].copy()
cd = h["Collection date"]
full = cd.str.match(r"^\\d{4}-\\d{2}-\\d{2}$")
h["date"] = pd.to_datetime(cd.where(full), errors="coerce")
h = h[h.date.notna()]
h["dec"] = lineage.decimal_date(h.date)
h = h[h.dec >= 2022.0]
h["country"] = h.Location.str.split(" / ").str[1]
print(len(h), "human genomes with a full collection date, 2022 onwards;", h.country.nunique(), "countries")
rows = []
for y in ["gc_content", "oe_TC", "oe_CG", "oe_TA"]:
    s = h[[y, "dec", "country", "length", "ambig_frac"]].dropna()
    r = smf.ols(f"{y} ~ dec + length + ambig_frac", data=s).fit(
        cov_type="cluster", cov_kwds={"groups": s.country}, use_t=True)
    lo, hi = r.conf_int().loc["dec"]
    rows.append({"feature": y, "mean": s[y].mean(), "slope per year": r.params["dec"], "ci_lo": lo, "ci_hi": hi, "p": r.pvalues["dec"]})
res = pd.DataFrame(rows).set_index("feature")
res.to_csv("results/mpox_time_slopes.tsv", sep="\\t")
res.round(6)'''),
code('''fig, axs = plt.subplots(1, 2, figsize=(9, 3.2))
h["year"] = h.dec.astype(int)
for ax, col, lab in zip(axs, ["gc_content", "oe_TC"], ["GC content", "TpC O/E (APOBEC3 target)"]):
    g = h.groupby("year")[col].agg(["mean", "sem", "count"])
    g = g[g["count"] >= 20]
    ax.errorbar(g.index, g["mean"], yerr=1.96 * g["sem"], marker="o", color="#2C7FB8", capsize=3)
    ax.set_xlabel("collection year"); ax.set_ylabel(lab); ax.set_xticks(g.index)
plt.tight_layout(); fig.savefig("results/figures/mpox_time_trend.png", dpi=200, bbox_inches="tight"); plt.show()'''),
md("""**Reading it.**

- **GC content:** the slope is about **-0.000025 per year**, which is **right where the APOBEC3 arithmetic predicts**. The 95% interval still includes zero (p = 0.18), because the signal is tiny next to the genome-to-genome noise (SD 0.0002).
- **TpC O/E:** falls by about 0.0002 per year (p = 0.005): the dinucleotide APOBEC3 edits is slowly being used up, as expected.
- **CpG O/E** also falls (p < 0.0001), which is not what APOBEC3 predicts as directly; an overall shift in which clades are sampled is the more likely cause.

**Limits (important).** The export has **no lineage column**, so this is the *all-lineages-pooled* analysis the RSV project showed to be unreliable: a later lineage replacing an earlier one can create a trend without any lineage changing (Simpson's paradox). The result is therefore exploratory. It needs lineage labels (Nextclade or the GISAID lineage field) and a re-run inside each lineage. With 55 countries the clustered errors are reasonable, but sequencing is dominated by the USA and Germany."""),

md("""## 4. What this means

| Question | Answer |
|---|---|
| Does the host field offer a contrast? | **No.** 94% human; the rest is wastewater, unknown host, or about 40 animal fragments with no full genome. |
| Is mpox CpG-suppressed? | **No** (CpG O/E 1.04), so the RSV CpG story does not carry over. |
| Does the APOBEC3 signature appear? | **Weakly.** GC and TpC fall at the rate predicted, but only TpC is statistically clear, and lineage replacement is not excluded. |

**Verdict.** Mpox does not offer the sporadic-versus-sustained contrast that the host field would need. It is worth a second look only for the APOBEC3 signature, and that requires lineage labels first."""),
])
