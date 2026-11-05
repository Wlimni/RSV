"""
Script Name: make_notebook_feasibility.py
Purpose:     Generate the virus-suitability notebook: applies the same tests to RSV, Nipah and mpox
             to show, with numbers, which virus can go through the RSV pipeline and why.
Inputs:      none (the notebook reads data/processed/{RSV features, NiV, Mpox})
Outputs:     notebooks/Feasibility/01_virus_suitability.ipynb
Notes:       Run: .venv/bin/python scripts/make_notebook_feasibility.py
             then: jupyter nbconvert --execute --inplace notebooks/Feasibility/01_virus_suitability.ipynb
"""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell as md, new_code_cell as code

SETUP = '''\
import sys, os
sys.path.insert(0, os.path.abspath("../../src"))
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from rsv import lineage, io

os.chdir(os.path.abspath("../.."))
pd.set_option("display.width", 130)
plt.rcParams.update({"figure.dpi": 110, "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False})
'''

cells = [
md("""# Which virus can go through the RSV pipeline?

---

### The question

After the RSV trend analysis came out negative, the suggestion was to try other viruses. Three were downloaded: **RSV** (done), **Nipah (NiV)** and **mpox (MPXV)**. This notebook applies the **same tests to all three** and shows, with numbers, which ones the pipeline can run on, and which cannot and why.

### The tests

The RSV pipeline needs the following. A virus has to pass them all.

| # | Test | Why the pipeline needs it |
|---|---|---|
| 1 | **Full-length genomes** in the thousands | Composition is measured over the whole genome or whole genes; fragments cover different regions and cannot be compared. |
| 2 | **At least two hosts**, each with at least 30 full genomes | Without it there is no host contrast (the reason to leave RSV). |
| 3 | **A lineage label** | Trends are fitted *inside* a lineage, to avoid one lineage replacing another (Simpson's paradox). |
| 4 | **Lineages with at least 200 genomes over at least 5 years** | The slope needs enough points spread over time (`lineage.prepare`). |
| 5 | **Genomes that are not near-copies** | One outbreak gives many identical sequences; counting them as independent overstates the evidence. |"""),
code(SETUP),

md("""## Test 1 and 2: full genomes and hosts"""),
code('''def host_group(h):
    h = str(h)
    if h == "Human": return "human"
    if "Environment" in h: return "environment (wastewater)"
    if h in ("unknown", "nan"): return "unknown"
    if "Sus" in h or "Canis" in h: return "pig/dog"
    if "Pteropus" in h or "Chiroptera" in h: return "bat"
    return "animal"

rsv = pd.read_csv("data/full/RSV_full_metadata.tsv", sep="\\t", dtype=str)
rsv["L"] = pd.to_numeric(rsv["Sequence Length"]); rsv["hostg"] = rsv.Host.map(host_group)
niv = pd.read_csv("data/processed/NiV/NiV_metadata.tsv", sep="\\t", dtype=str)
niv["L"] = pd.to_numeric(niv["Sequence Length"]); niv["hostg"] = niv.Host.map(host_group)
mpx = pd.read_csv("data/processed/Mpox/Mpox_features.tsv", sep="\\t", dtype={"Collection date": str})
mpx["L"] = mpx.length; mpx["hostg"] = mpx.Host.map(host_group)
mpx_ok = mpx[mpx.L.between(190000, 200000) & (mpx.ambig_frac <= 0.01)]

# "full length" = within the normal size range of that virus' genome
FULL = {"RSV": (14900, 15400), "NiV": (17000, 19000)}
rsv_full = rsv[rsv.L.between(*FULL["RSV"])]
niv_full = niv[niv.L.between(*FULL["NiV"])]
tests = []
for name, allr, full in (("RSV", rsv, rsv_full), ("Nipah", niv, niv_full), ("Mpox", mpx, mpx_ok)):
    hc = full.hostg.value_counts()
    tests.append({"virus": name, "genomes downloaded": len(allr), "full-length genomes": len(full),
                  "% full-length": round(100 * len(full) / len(allr)),
                  "hosts with >=30 full genomes": ", ".join(f"{h} ({n})" for h, n in hc.items() if n >= 30 and h != "unknown") or "none",
                  "n hosts": int(((hc >= 30) & (hc.index != "unknown")).sum())})
t = pd.DataFrame(tests).set_index("virus")
t'''),
md("""**Reading it.**

- **RSV** passes test 1 (33,040) but has **one host** (human), so there is no host contrast.
- **Nipah** has 106 full genomes in total (31% of the 343), and **only humans reach 30** (74); bats have 19 and pigs/dogs 9. Most bat sequences (197 of 216) are fragments. So there is host variation, but not enough full genomes per host.
- **Mpox** has plenty of full genomes (7,711) but they are all **human**. The wastewater (220) and animal (40) entries are fragments, so no second host reaches 30."""),

md("""## Test 3 and 4: lineage labels, and the pipeline's own filter

Here the *actual* RSV filter, `lineage.prepare` (3-level label, at least 5 years, at least 200 genomes, date window), is applied to each virus. For Nipah, collection years are turned into mid-year dates because many records give only the year."""),
code('''feat = pd.read_csv("data/processed/rsv_features.tsv", sep="\\t", parse_dates=["date"])
_, kept_rsv = lineage.prepare(feat)

n = niv.copy()
n["year"] = pd.to_numeric(n["Collection date"].str[:4], errors="coerce")
n["date"] = pd.to_datetime(n.year.dropna().astype(int).astype(str) + "-07-01")
n["Subtype"] = n.Clade.str.extract(r"\\((.*)\\)")[0].fillna("NA")
n["Lineage"] = n.Lineage.fillna("none")
n = n[n.date.notna()]
_, kept_niv = lineage.prepare(n, start=1990.0, end=2027.0, trim=0)       # same rule, wider date window
_, kept_niv_full = lineage.prepare(n[n.L.between(17000, 19000)], start=1990.0, end=2027.0, trim=0)

summary = pd.DataFrame({"lineages passing the pipeline filter": [len(kept_rsv), len(kept_niv), len(kept_niv_full), 0]},
                       index=["RSV (all QC genomes)", "Nipah (all sequences)", "Nipah (full genomes only)", "Mpox (no lineage column)"])
summary'''),
code('''# what Nipah's lineages look like against the thresholds (200 genomes, 5 years)
g = n.groupby(["Subtype", "Lineage"]).agg(genomes=("year", "size"), years_with_data=("year", "nunique"),
                                          first=("year", "min"), last=("year", "max"))
full_ct = n[n.L.between(17000, 19000)].groupby(["Subtype", "Lineage"]).size().rename("full genomes")
g = g.join(full_ct).fillna(0).astype({"full genomes": int})
g["passes 200 genomes"] = g.genomes >= 200
g'''),
md("""**Reading it.**

- **RSV: 9 lineages pass**, which is why the pipeline works there.
- **Nipah: 0 pass.** The closest is **BD2: 189 sequences in total, but only 132 have a collection year (needed for a trend), and only 39 are full genomes.** That is well short of 200, over a span where most human genomes come from a few outbreaks.
- **Mpox: not testable.** The GISAID export has no lineage or clade column. This test cannot be passed until lineage labels are added (Nextclade, or the GISAID lineage field if it exists)."""),

md("""## Test 5: are the genomes independent?

Counting related genomes as independent observations overstates the evidence (this is what inflated the RSV p-values). First check the simplest case, exact duplicates, among Nipah's full genomes."""),
code('''seqs = io.read_fasta("data/processed/NiV/NiV.fasta")
f = niv_full.copy()
f["seq"] = f["Accession ID"].map(seqs)
rows = []
for h, gg in f.groupby("hostg"):
    rows.append({"host": h, "full genomes": len(gg), "distinct sequences": gg.seq.nunique(),
                 "% duplicated": round(100 * (1 - gg.seq.nunique() / len(gg)))})
pd.DataFrame(rows).set_index("host")'''),
code('''# how many human genomes come from each outbreak country
f["country"] = f.Location.str.split(" / ").str[1]
print(pd.crosstab(f.country, f.hostg, margins=True))'''),
md("""**Reading it.** Exact duplicates are rare (1 of 74 human genomes), so duplication is **not** the problem here. The dependence is by **place and outbreak**: 57 of the 74 human genomes are from Bangladesh. Genomes from one outbreak share ancestry, so the effective number of independent observations is smaller than 74 (this was not measured here; it would need a tree)."""),

md("""## Verdict

| Virus | Can the RSV pipeline run? | Why |
|---|---|---|
| **RSV** | Yes (already run) | 9 lineages pass; but one host, so no host contrast. |
| **Nipah** | **No** | 0 lineages pass (closest: BD2, 189 sequences, 132 dated, 39 full genomes); only 106 full genomes; no host has 30 full genomes except humans, and most human genomes come from one country. |
| **Mpox** | **Not as it stands** | No lineage labels; all full genomes are human; also a DNA virus with no CpG suppression (O/E 1.04), so the CpG question has nothing to explain. |

What each *did* give: Nipah, a host and timing picture (BD2 in bats and humans 2004-2025); mpox, an exploratory APOBEC3 trend (see their own notebooks).

---
## What a suitable virus must look like

The tests above give the checklist for the next download:

1. **Tens of thousands of full genomes** (not fragments).
2. **At least two hosts, each with at least 30 full genomes**, spread over several years.
3. **A lineage or clade column** (or enough sequences to build a tree).
4. **At least one lineage with 200 or more genomes over 5 or more years**.
5. **Sequences from many independent outbreaks**, not one.

## The next step

| Candidate | Tests it is likely to pass | Cheapest way to check |
|---|---|---|
| **SARS-CoV-2** (hCoV-19) | 1, 3, 4, 5 pass trivially (millions of genomes, Pango lineages). Test 2 is the question: mink and deer give a few hundred to a few thousand genomes. | Export host, date and lineage only (no sequences), then count genomes per host. |
| **H5N1** (EpiFlu) | Test 2 should pass (birds, cattle and other mammals); lineage labels depend on a clade field. | The same metadata-only export: host, date, clade or lineage. |

**Do this next:** export the metadata only for those two (host, collection date, lineage or clade, location, accession) and run this same notebook on them. It takes minutes, and it tells us which one to download in full. Two outcomes are possible: one passes all five tests and gets the full pipeline, or neither does and the project reports a negative across viruses."""),
]

nb = new_notebook()
nb.cells = cells
nb.metadata["kernelspec"] = {"name": "rsv", "display_name": "rsv (.venv)", "language": "python"}
nbf.write(nb, "notebooks/Feasibility/01_virus_suitability.ipynb")
print("wrote notebooks/Feasibility/01_virus_suitability.ipynb")
