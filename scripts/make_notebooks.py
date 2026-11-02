"""
Script Name: make_notebooks.py
Purpose:     Generate the four analysis notebooks. The notebooks are the main pipeline and
             are written to be read: each step has a markdown cell explaining what it does,
             and every section ends with a figure.
Inputs:      none
Outputs:     notebooks/01..04 *.ipynb
Notes:       Notebook content is held here as plain Python so it stays version-controlled
             and regenerable. Run: .venv/bin/python scripts/make_notebooks.py
             Then execute: jupyter nbconvert --execute --inplace notebooks/*.ipynb
"""
import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

os.makedirs("notebooks", exist_ok=True)

SETUP = '''\
import sys, os, subprocess
sys.path.insert(0, os.path.abspath("../src"))
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from IPython.display import Image, display

ROOT = os.path.abspath("..")          # repo root; the R scripts expect to run from here
pd.set_option("display.width", 130)
plt.rcParams.update({"figure.dpi": 110, "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False})
PAL = {"A": "#2C7FB8", "B": "#C62824"}    # same colours as the R figures

def run_r(script="R/figures.R"):
    """Run the R figure script from the repo root and report what it wrote."""
    r = subprocess.run(["Rscript", script], capture_output=True, text=True, cwd=ROOT)
    print(r.stdout.strip() or r.stderr.strip()[-400:])

def show(fig):
    """Display a figure from results/figures/ inline."""
    display(Image(os.path.join(ROOT, "results", "figures", fig)))
'''


def build(path, cells):
    nb = new_notebook()
    nb.cells = cells
    # name the project kernel explicitly: a bare "python3" can resolve to the Anaconda
    # base interpreter, whose NumPy 2.x breaks extensions compiled against 1.x
    nb.metadata["kernelspec"] = {"name": "rsv", "display_name": "rsv (.venv)",
                                 "language": "python"}
    nbf.write(nb, path)
    print("wrote", path)


# ----------------------------------------------------------------- 01 the data
build("notebooks/01_data.ipynb", [
    new_markdown_cell(
        "# 01 — The data\n\n"
        "---\n\n"
        "### The question\n\n"
        "Ye et al. 2026 (*Nature Microbiology*) measured GC content and GC-containing "
        "dinucleotide frequencies across 115,520 influenza A genomes, and found both "
        "**declined over time** within lineages that had achieved sustained transmission in "
        "mammals.\n\n"
        "Applying the same measurements to human RSV:\n\n"
        "> **Within each RSV lineage, how does genome composition change over time, and does "
        "the pattern differ between RSV-A and RSV-B?**\n\n"
        "This notebook covers the dataset: size, coverage in time and place, and what QC "
        "removes.\n\n"
        "---"),
    new_code_cell(SETUP),
    new_markdown_cell(
        "## 1.1  The sequences\n\n"
        "`rsv_features.tsv` is written by `scripts/build_features.py`: one row per genome, "
        "GISAID metadata alongside composition features."),
    new_code_cell(
        'd = pd.read_csv("../data/processed/rsv_features.tsv", sep="\\t")\n'
        'print("genomes :", len(d))\n'
        'print("subtypes:", d.Subtype.value_counts().to_dict())\n'
        'print("years   :", d.year.min(), "to", d.year.max())\n'
        '# accession IDs are dropped from displays so no GISAID record is stored here\n'
        'd.drop(columns=["Accession ID"]).head(3)'),
    new_markdown_cell(
        "## 1.2  Quality control\n\n"
        "An RSV genome is about 15,200 nucleotides. We keep a sequence only if it is\n\n"
        "| Check | Threshold | Why |\n"
        "|---|---|---|\n"
        "| Length | 14,900–15,400 nt | keeps complete genomes, drops fragments |\n"
        "| Ambiguous bases | at most 1% | `N` positions would distort letter counts |\n"
        "| Collection date | complete | needed to place the genome in time |\n\n"
        "The filter is deliberately gentle, because the GISAID download was already "
        "restricted to complete, high-coverage genomes."),
    new_code_cell(
        'fig, ax = plt.subplots(figsize=(7, 3))\n'
        'for s in "AB":\n'
        '    ax.hist(d.loc[d.Subtype == s, "seq_len_actual"], bins=110,\n'
        '            alpha=.6, label=f"RSV-{s}", color=PAL[s])\n'
        'for x in (14900, 15400):\n'
        '    ax.axvline(x, ls="--", lw=1, c="k")\n'
        'ax.set_xlim(14700, 15600)\n'
        'ax.set_xlabel("genome length (nt)"); ax.set_ylabel("genomes")\n'
        'ax.set_title("Genome length, dashed lines mark the filter window")\n'
        'ax.legend(frameon=False); plt.tight_layout(); plt.show()'),
    new_markdown_cell(
        "**Two genomes out of 33,040 were removed.** One RSV-B had 2.11% ambiguous bases. "
        "The other was an RSV-B at 15,854 nt — 632 nt over the reference, carrying a "
        "583-nt non-RSV tail past the 3' end at 49% GC against the virus's 34%, i.e. "
        "contamination. That second case is why the length check has an upper bound and "
        "not only a minimum."),
    new_markdown_cell(
        "## 1.3  When and where the sequences come from\n\n"
        "Two sampling features matter for everything downstream:\n\n"
        "- **Not even over time** — rises steeply from 2018, dips in 2020–21 when COVID-19 "
        "disrupted both RSV circulation and surveillance.\n"
        "- **Not even over the world** — Europe contributes far more genomes than anywhere "
        "else.\n\n"
        "So any trend over calendar time is partly a trend in *who was sequencing*. This is "
        "the main reason notebook 03 fits inside lineages rather than across the whole "
        "dataset."),
    new_code_cell(
        'fig, ax = plt.subplots(1, 2, figsize=(11, 3.2))\n'
        'yr = d.groupby(["year", "Subtype"]).size().unstack(fill_value=0)\n'
        'yr.plot.bar(stacked=True, ax=ax[0], color=[PAL["A"], PAL["B"]], width=.85,\n'
        '            legend=False)\n'
        'ax[0].set_title("genomes per year"); ax[0].set_xlabel("")\n'
        'ax[0].axvspan(yr.index.get_loc(2020) - .5, yr.index.get_loc(2021) + .5,\n'
        '              color="grey", alpha=.18)\n'
        'reg = d.groupby(["region", "Subtype"]).size().unstack(fill_value=0).sort_values("A")\n'
        'reg.plot.barh(stacked=True, ax=ax[1], color=[PAL["A"], PAL["B"]])\n'
        'ax[1].set_title("genomes per region"); ax[1].set_ylabel("")\n'
        'ax[1].legend(title="subtype", frameon=False)\n'
        'plt.tight_layout(); plt.show()'),
    new_markdown_cell(
        "*Grey band marks the COVID-19 disruption.*"),
    new_markdown_cell(
        "## 1.4  Known problems with this dataset\n\n"
        "Four, found by auditing before analysing. Two change what we do.\n\n"
        "| | Problem | Decision |\n"
        "|---|---|---|\n"
        "| **a** | 2026 tail incomplete (deposition lag) | cut the window at end of May 2026 |\n"
        "| **b** | Sequencing geography shifted, even inside lineages | region as a covariate |\n"
        "| **c** | Strong seasonality | noted, not corrected |\n"
        "| **d** | Sparse early data gives single genomes high leverage | trim sparse lineage ends |\n\n"
        "---\n\n"
        "### (a) The 2026 tail is incomplete\n\n"
        "The download was made in October 2026, but sequences take months to be deposited. "
        "Comparing 2026 month by month against the same months of 2025 separates deposition "
        "lag from ordinary seasonality."),
    new_code_cell(
        'm = {y: d[d.date.str.startswith(str(y))].date.str.slice(5, 7).value_counts().sort_index()\n'
        '     for y in (2025, 2026)}\n'
        'cmp = pd.DataFrame(m).fillna(0).astype(int)\n'
        'cmp["2026 as % of 2025"] = (cmp[2026] / cmp[2025] * 100).round(0)\n'
        'print(cmp.to_string())\n'
        'print("\\nlast collection date:", d.date.max())'),
    new_markdown_cell(
        "- Through May, 2026 holds **78–93%** of the previous year’s count — normal, and May "
        "itself is at **92%**.\n"
        "- June falls to **53%**, July **29%**, August **4%** — eighteen genomes where 2025 "
        "had 399 in the same month.\n"
        "- That is deposition lag, not season: those months are a non-random subset "
        "contributed only by laboratories that deposit quickly.\n\n"
        "**Decision: the window runs 2001 to the end of May 2026.** May is the last month "
        "still at full strength, so it is the natural cut. Everything earlier is "
        "kept — sequencing effort grows steeply over the period, but uneven sample sizes "
        "across time do not bias a slope, they only make later years more precisely "
        "estimated.\n\n"
        "---\n\n"
        "### (b) Where sequencing happens has shifted completely\n\n"
        "The more serious problem."),
    new_code_cell(
        'reg = pd.crosstab(d.year, d.region, normalize="index").mul(100).round(1)\n'
        'reg.loc[2015:2026]'),
    new_markdown_cell(
        "- Oceania supplied **38%** of genomes in 2015 and **4%** in 2026; Europe went from "
        "**18%** to **56%**. North America contributes nothing at all in 2017.\n"
        "- So calendar time and geography are **entangled** — a composition trend over time "
        "could partly be a trend in which part of the world was sequencing.\n"
        "- The drift happens *inside* individual lineages too, so fitting within lineages "
        "does not remove it. Lineage A.D goes from 71% Oceania early to 36% Europe late.\n\n"
        "But a shifting region mix only matters if composition actually differs by region. "
        "So: the composition across all years, and the same measure split by region."),
    new_code_cell(
        'w = d[d.year.between(2015, 2026)]\n'
        'fig, ax = plt.subplots(1, 2, figsize=(11, 3.2))\n'
        '\n'
        '# left: every genome pooled by year, ignoring lineage structure entirely\n'
        'for s in "AB":\n'
        '    t = w[w.Subtype == s].groupby("year").oe_CG.mean()\n'
        '    ax[0].plot(t.index, t.values, "-o", ms=4, color=PAL[s], label=f"RSV-{s}")\n'
        'ax[0].set_xlabel("year"); ax[0].set_ylabel("mean CpG O/E")\n'
        'ax[0].set_title("Composition across all years, pooled")\n'
        'ax[0].legend(frameon=False)\n'
        '\n'
        '# right: the same measure by region, to see whether region carries any signal\n'
        'r = w.pivot_table(index="region", columns="Subtype", values="oe_CG", aggfunc="mean")\n'
        'r.plot.barh(ax=ax[1], color=[PAL["A"], PAL["B"]])\n'
        'ax[1].legend(["RSV-A", "RSV-B"], frameon=False, loc="lower right")\n'
        'ax[1].set_xlim(0.230, 0.245); ax[1].set_ylabel("")\n'
        'ax[1].set_xlabel("mean CpG O/E"); ax[1].set_title("The same measure, by region")\n'
        'plt.tight_layout(); plt.show()\n'
        '\n'
        'print("spread across regions (max - min), within subtype:")\n'
        'for v in ["gc_content", "oe_CG"]:\n'
        '    for s in "AB":\n'
        '        x = w[w.Subtype == s].groupby("region")[v].mean()\n'
        '        print(f"  {v:11s} RSV-{s}: {x.max() - x.min():.5f}")'),
    new_markdown_cell(
        "Both halves matter:\n\n"
        "- **The pooled line does move** — RSV-B's CpG O/E climbs from 0.2357 in 2015 to "
        "0.2430 in 2023. Taken at face value that is a trend, and it is exactly the kind of "
        "number that section 3.1 explains you cannot trust, because lineages are replacing "
        "one another underneath it.\n"
        "- **Region carries real signal.** Within RSV-B the regional spread in CpG O/E is "
        "**0.0070**; the largest within-lineage slope anywhere in this project is 0.0011 per "
        "year. So the gap between two regions is worth roughly six years of trend. A shift "
        "in which regions were sequencing can therefore imitate a time trend.\n\n"
        "**Decision: region enters the regressions as a covariate, and results are reported "
        "both ways** (notebook 03).\n\n"
        "---\n\n"
        "### (c) Strong seasonality — noted, not corrected\n\n"
        "RSV is seasonal, so genomes are not spread evenly through the year."),
    new_code_cell(
        'mo = d.date.str.slice(5, 7).value_counts().sort_index()\n'
        'fig, ax = plt.subplots(figsize=(7, 2.6))\n'
        'ax.bar(range(1, 13), mo.values, color="#6B8FA8", edgecolor="grey", linewidth=.3)\n'
        'ax.set_xticks(range(1, 13)); ax.set_xlabel("month"); ax.set_ylabel("genomes")\n'
        'ax.set_title("Collection month, all years pooled")\n'
        'plt.tight_layout(); plt.show()'),
    new_markdown_cell(
        "November and December carry three to four times September. The two hemispheres also "
        "peak at opposite times, so season and geography are partly entangled as well. Minor "
        "for a slope measured across a decade, so we record it rather than correct it.\n\n"
        "---\n\n"
        "### (d) Sparse early data gives single genomes high leverage\n\n"
        "- Everything before 2015 is about 3% of the dataset, and **RSV-A has essentially no "
        "sequences before 2010**.\n"
        "- Keeping it is still worth doing: for some lineages it adds years of span, and a "
        "longer span makes a slope more precise.\n"
        "- But the leverage is dangerous — lineage A.D.3 would reach back **6.9 extra years "
        "on the strength of one genome from 2010**, and a point that isolated drags the "
        "fitted line hard.\n\n"
        "**Decision: each lineage is cut back to the run of years in which it was actually "
        "being sampled**, dropping leading and trailing years holding fewer than 5 genomes. "
        "The rule, and what it does to each lineage:"),
    new_code_cell(
        'f = pd.read_csv("../data/processed/rsv_features.tsv", sep="\\t")\n'
        'f["lineage_grp"] = f.Lineage.astype(str).str.split(".").str[:3].str.join(".")\n'
        'f["yr"] = f.date.str.slice(0, 4).astype(int)\n'
        '\n'
        'rows = []\n'
        'for (sub, grp), s in f.groupby(["Subtype", "lineage_grp"]):\n'
        '    counts = s.yr.value_counts().sort_index()\n'
        '    ok = counts[counts >= 5].index                      # years holding >= 5 genomes\n'
        '    if len(ok) == 0:\n'
        '        continue\n'
        '    kept = s[(s.yr >= ok.min()) & (s.yr <= ok.max())]   # keep the run between them\n'
        '    if kept.yr.nunique() >= 5 and len(kept) >= 200:     # the inclusion rule\n'
        '        rows.append({"Subtype": sub, "lineage": grp,\n'
        '                     "first yr raw": s.yr.min(), "first yr kept": kept.yr.min(),\n'
        '                     "yrs cut": kept.yr.min() - s.yr.min(),\n'
        '                     "genomes cut": len(s) - len(kept)})\n'
        'print(pd.DataFrame(rows).sort_values(["Subtype", "lineage"]).to_string(index=False))\n'
        '\n'
        '# well-sampled early data is not touched by the rule\n'
        'ad = f[(f.Subtype == "A") & (f.lineage_grp == "A.D")]\n'
        'c = ad.yr.value_counts().sort_index(); ok = c[c >= 5].index\n'
        'kept = ad[(ad.yr >= ok.min()) & (ad.yr <= ok.max())]\n'
        'dt = pd.to_datetime(kept.date)\n'
        'print("\\nA.D genomes within 3 years of its first:",\n'
        '      int((dt < dt.min() + pd.DateOffset(years=3)).sum()))'),
    new_markdown_cell(
        "**A.D.3 loses 7 years and only 1 genome** — that is the whole point of the rule: it "
        "removes span, not data. A.D keeps **361 genomes in its first three years**, so "
        "well-sampled early data is untouched.\n\n"
        "This is the same rule `lineage.prepare()` applies in notebooks 03 and 04; the kept "
        "start years here match `results/lineages.tsv` exactly. It lives in "
        "`src/rsv/lineage.py` only because those notebooks and the analysis script all need "
        "it."),
    new_markdown_cell(
        "## 1.5  Lineages\n\n"
        "GISAID assigns each genome to a **lineage**, a named branch of the RSV family tree. "
        "Labels are hierarchical and dot-separated: `A.D.1.5` sits inside `A.D.1`, which sits "
        "inside `A.D`.\n\n"
        "Lineages rise and fall over a few seasons. Notebook 03 covers why that matters and "
        "how the level is chosen."),
    new_code_cell(
        'for s in "AB":\n'
        '    top = d[d.Subtype == s].Lineage.value_counts().head(5)\n'
        '    print(f"RSV-{s} most common lineage labels")\n'
        '    print(top.to_string(), "\\n")'),
])

# ------------------------------------------------------------- 02 what we measure
build("notebooks/02_features.ipynb", [
    new_markdown_cell(
        "# 02 — What we measure\n\n"
        "---\n\n"
        "Every feature here is a property of the **letter sequence itself**, independent of "
        "what the genome codes for.\n\n"
        "| Feature | Definition |\n"
        "|---|---|\n"
        "| **GC content** | fraction of letters that are G or C |\n"
        "| **Dinucleotide frequency** | how often a given pair of adjacent letters occurs |\n"
        "| **O/E ratio** | observed pair frequency divided by what chance predicts |\n"
        "| **GC3 / GC12** | GC content at the third codon position / at the first and second |\n\n"
        "### The pairs tracked\n\n"
        "| Pair | Written | Note |\n"
        "|---|---|---|\n"
        "| C then G | CpG | the motif of interest |\n"
        "| G then C | GpC | |\n"
        "| G then G | GpG | |\n"
        "| C then C | CpC | |\n\n"
        "### Why O/E rather than raw counts\n\n"
        "A GC-rich genome contains more CpG by chance alone, so we divide by chance "
        "expectation:\n\n"
        "$$ \\mathrm{O/E}(\\mathrm{CpG}) = \\frac{f(\\mathrm{CpG})}{f(\\mathrm{C}) \\times f(\\mathrm{G})} $$\n\n"
        "Near **1** = as common as chance predicts. **Below 1** = the virus is avoiding that "
        "pair.\n\n"
        "### Why split GC by codon position\n\n"
        "- **GC3** — third position; most changes are synonymous and leave the protein "
        "unchanged.\n"
        "- **GC12** — first and second; most changes do alter the protein.\n\n"
        "Comparing how fast the two move separates pressure on the RNA from pressure on the "
        "protein — the same logic as dS versus dN. The paper reports GC3 changing **5.65× "
        "faster** than GC12 in influenza.\n\n"
        "---"),
    new_code_cell(SETUP + "from rsv import composition"),
    new_markdown_cell(
        "## 2.1  The functions, on a short example\n\n"
        "All measurements come from `src/rsv/composition.py`. Ambiguous positions are skipped "
        "rather than guessed."),
    new_code_cell(
        'demo = "ATGGCGTTAACGCGTAAATTACGATAG"\n'
        'f = composition.all_features(demo, coding=True)\n'
        'for k in ["gc_content", "gc3", "gc12", "per1k_CG", "oe_CG"]:\n'
        '    print(f"{k:12s} {f[k]:.4f}")'),
    new_markdown_cell(
        "## 2.2  The measurements across all genomes\n\n"
        "`GENOME` means the ten coding regions concatenated, which is how the paper defines "
        "genome-level GC, GC3 and GC12 — it keeps the reading frame well defined."),
    new_code_cell(
        'g = pd.read_csv("../data/processed/rsv_gene_features.tsv", sep="\\t")\n'
        'genome = g[g.gene == "GENOME"]\n'
        'cols = ["gc_content", "gc3", "gc12", "oe_CG", "oe_GC", "oe_GG", "oe_CC"]\n'
        'genome.groupby("Subtype")[cols].mean().round(4)'),
    new_markdown_cell(
        "### Reading that table\n\n"
        "- **CpG O/E sits far below 1 in both subtypes** — RSV carries only a fraction of the "
        "CpG pairs its base composition would predict. The other GC pairs sit near or above "
        "1, so the effect is specific to CpG, not a general shortage of G and C.\n"
        "- **GC3 is much lower than GC12** — typical, because third positions are free to "
        "drift when changing them does not change the protein."),
    new_code_cell(
        'fig, ax = plt.subplots(1, 2, figsize=(11, 3.4))\n'
        'for s in "AB":\n'
        '    x = genome[genome.Subtype == s]\n'
        '    ax[0].hist(x.gc_content, bins=70, alpha=.6, color=PAL[s], label=f"RSV-{s}")\n'
        'ax[0].set_xlabel("GC content"); ax[0].set_ylabel("genomes")\n'
        'ax[0].set_title("GC content"); ax[0].legend(frameon=False)\n'
        'oe = ["oe_CG", "oe_GC", "oe_GG", "oe_CC"]\n'
        'names = ["CpG", "GpC", "GpG", "CpC"]\n'
        'pos = np.arange(len(oe))\n'
        'for i, s in enumerate("AB"):\n'
        '    x = genome[genome.Subtype == s]\n'
        '    ax[1].bar(pos + (i - .5) * .38, [x[c].mean() for c in oe], width=.38,\n'
        '              color=PAL[s], label=f"RSV-{s}", edgecolor="0.3", linewidth=.3)\n'
        'ax[1].axhline(1, ls="--", lw=1, c="k")\n'
        'ax[1].set_xticks(pos); ax[1].set_xticklabels(names)\n'
        'ax[1].set_ylabel("mean O/E"); ax[1].set_title("Dinucleotide O/E (dashed line = 1)")\n'
        'ax[1].legend(frameon=False)\n'
        'plt.tight_layout(); plt.show()'),
    new_markdown_cell(
        "Only **CpG** falls clearly below the dashed line. That is the pattern reported "
        "across RNA viruses that circulate in humans."),
])

# ------------------------------------------------- 03 the core: time within lineage
build("notebooks/03_lineage_trends.ipynb", [
    new_markdown_cell(
        "# 03 — Composition over time, within lineages\n\n"
        "---\n\n"
        "**The core analysis.** Follows Figures 2 and 3 of the paper.\n\n"
        "### Why fit inside lineages\n\n"
        "A dataset-wide average can move for two different reasons:\n\n"
        "1. the genomes inside each lineage are genuinely changing, or\n"
        "2. nothing changes, but a lineage that already had a high value grows while one with "
        "a low value shrinks.\n\n"
        "The second produces a confident-looking trend with no evolution behind it. Fitting "
        "**separately inside each lineage** removes it, which is what the paper does.\n\n"
        "### Choosing the lineage level\n\n"
        "Full GISAID labels are too fine — most span about 6 years with a few hundred "
        "genomes, and a slope's uncertainty grows as both span and sample size shrink. "
        "Truncating the label to **three levels** gives 9 well-populated lineages, close to "
        "the 12 the paper used.\n\n"
        "---"),
    new_code_cell(SETUP + "from rsv import lineage"),
    new_markdown_cell(
        "## 3.1  The lineages we analyse\n\n"
        "Inclusion requires **at least 5 distinct years and at least 200 genomes**. The date "
        "window (2001 to the end of May 2026) and the sparse-end trimming are both applied inside "
        "`lineage.prepare()`, so no extra year filter is needed here."),
    new_code_cell(
        'g = pd.read_csv("../data/processed/rsv_gene_features.tsv", sep="\\t")\n'
        'dates = pd.read_csv("../data/processed/rsv_features.tsv", sep="\\t",\n'
        '                    usecols=["Accession ID", "date"])\n'
        'genome = g[g.gene == "GENOME"].merge(dates, on="Accession ID")\n'
        'gdf, keep = lineage.prepare(genome)\n'
        'print(f"{len(keep)} lineages, {len(gdf):,} genomes\\n")\n'
        'keep.round(2)'),
    new_markdown_cell(
        "Two things to note:\n\n"
        "- **The lineages do not cover the same calendar window.** B.D runs 2006–2016 while "
        "B.D.E runs 2018–2026, so they barely overlap. Comparisons between lineages are "
        "therefore partly comparisons between different periods.\n"
        "- **RSV-B contributes three lineages, and one of them (B.D) is small** at 342 "
        "genomes. Any A-versus-B statement rests on a narrow base."),
    new_markdown_cell(
        "## 3.2  Trend of each feature within each lineage\n\n"
        "Each feature is regressed against the **decimal collection date**, inside each "
        "lineage separately. The slope is the change per year."),
    new_code_cell(
        'FEATS = ["gc_content", "gc3", "gc12",\n'
        '         "per1k_CG", "per1k_GC", "per1k_GG", "per1k_CC",\n'
        '         "oe_CG", "oe_GC", "oe_GG", "oe_CC"]\n'
        'slopes = lineage.lineage_slopes(gdf, FEATS)\n'
        'slopes[slopes.feature == "gc_content"][\n'
        '    ["Subtype", "lineage_grp", "n", "slope", "ci_lo", "ci_hi", "r2", "significant"]\n'
        '].round(6)'),
    new_markdown_cell(
        "### Figure — GC content over time, one panel per lineage\n\n"
        "The equivalent of the paper's Figure 2a. Points are individual genomes; the line is "
        "the fitted trend with its confidence band."),
    new_code_cell('run_r()'),
    new_code_cell('show("fig1_lineage_trends.png")'),
    new_markdown_cell(
        "## 3.3  Where along the codon is GC changing?\n\n"
        "The sharpest test available:\n\n"
        "- pressure acting on the **RNA** → strongest effect at third codon positions, where "
        "changes are mostly synonymous\n"
        "- pressure acting on the **protein** → first and second positions move comparably"),
    new_code_cell(
        'piv = (slopes[slopes.feature.isin(["gc_content", "gc3", "gc12"])]\n'
        '       .pivot_table(index=["Subtype", "lineage_grp"], columns="feature", values="slope"))\n'
        'piv["gc3 / gc12"] = (piv.gc3 / piv.gc12).round(2)\n'
        'piv.round(6)'),
    new_code_cell('show("fig2_gc3_vs_gc12.png")'),
    new_markdown_cell(
        "### How this differs from the paper\n\n"
        "| | influenza (paper) | RSV (here) |\n"
        "|---|---|---|\n"
        "| GC3 | falls | **rises in 6 of 9 lineages** |\n"
        "| GC12 | falls | **falls in 6 of 9** |\n"
        "| relationship | same direction, GC3 5.65× faster | **opposite directions**, ratio negative |\n\n"
        "Two things follow:\n\n"
        "- Synonymous and non-synonymous sites are clearly under different pressure, which is "
        "interesting in itself.\n"
        "- This is **not** the influenza signature and should not be described as reproducing "
        "it.\n\n"
        "Note the size: slopes are around 10⁻⁴ per year, so roughly a tenth of a percentage "
        "point of GC across the whole decade. Small."),
    new_markdown_cell(
        "## 3.4  The dinucleotides\n\n"
        "The same regression applied to each dinucleotide, reported both as frequency per "
        "1,000 nucleotides (the paper's unit) and as the O/E ratio. The O/E version matters "
        "because it already corrects for GC content, so a trend there is not simply GC "
        "content moving."),
    new_code_cell(
        'cg = slopes[slopes.feature == "oe_CG"][\n'
        '    ["Subtype", "lineage_grp", "n", "slope", "ci_lo", "ci_hi", "p", "significant"]]\n'
        'cg.round(6)'),
    new_code_cell('show("fig3_dinucleotide_slopes.png")'),
    new_markdown_cell(
        "### Does geography explain these slopes?\n\n"
        "Notebook 01 showed that where RSV was sequenced shifted over the period, including "
        "inside single lineages. Repeating the fit with **region held constant** tells us "
        "whether a slope is about time or about place."),
    new_code_cell(
        'adj = lineage.lineage_slopes(gdf, ["oe_CG"], covariates=("region",))\n'
        'cmp = (slopes[slopes.feature == "oe_CG"][["Subtype","lineage_grp","n","slope","p"]]\n'
        '       .merge(adj[["Subtype","lineage_grp","slope","p"]],\n'
        '              on=["Subtype","lineage_grp"], suffixes=("", "_adj")))\n'
        'cmp["same_direction"] = (cmp.slope * cmp.slope_adj) > 0\n'
        'cmp.round(6)'),
    new_markdown_cell(
        "Most slopes keep their direction and rough size, so geography is not generating "
        "them. The ones that move most are the lineages notebook 01 flagged as having the "
        "largest geographic drift — the adjustment bites exactly where the problem is."),
    new_markdown_cell(
        "## 3.5  Is the dinucleotide signal just GC content moving?\n\n"
        "The paper's Figure 3b. A dinucleotide that correlates almost perfectly with GC "
        "content carries no information beyond GC; correlations well below 1 mean it is "
        "moving on its own."),
    new_code_cell(
        'corr = lineage.gc_vs_dinucleotide_correlation(\n'
        '    gdf, ["per1k_CG", "per1k_GC", "per1k_GG", "per1k_CC"])\n'
        'corr.groupby("feature").r.agg(["mean", "min", "max"]).round(3)'),
    new_code_cell('show("fig4_gc_correlation.png")'),
    new_markdown_cell(
        "## 3.6  Sensitivity: does the lineage level drive the result?\n\n"
        "Truncating to three levels was a choice. Repeating at level 4 gives more, smaller "
        "lineages. If direction and rough size agree, the choice is not driving the "
        "conclusion."),
    new_code_cell(
        'gdf4, keep4 = lineage.prepare(genome, level=4)\n'
        's4 = lineage.lineage_slopes(gdf4, ["gc_content", "gc3", "oe_CG"])\n'
        'cmp = pd.concat([slopes[slopes.feature.isin(["gc_content","gc3","oe_CG"])].assign(level=3),\n'
        '                 s4.assign(level=4)])\n'
        'cmp.groupby(["level", "feature"]).slope.agg(\n'
        '    median="median", lineages="size",\n'
        '    n_negative=lambda s: int((s < 0).sum())).round(6)'),
    new_markdown_cell(
        "**The CpG result does not survive this check.**\n\n"
        "| Lineage level | median slope | lineages negative |\n"
        "|---|---|---|\n"
        "| 3 | slightly positive | 3 of 9 |\n"
        "| 4 | negative | 15 of 20 |\n\n"
        "The direction flips depending on how finely lineages are defined, which means the "
        "granularity is driving it rather than the biology.\n\n"
        "GC3 is steadier — positive at both levels, though weaker at the finer one. GC "
        "content overall sits near zero either way."),
    new_markdown_cell(
        "---\n\n"
        "## 3.7  What this notebook found\n\n"
        "Stated plainly, because the result is mostly negative.\n\n"
        "| | Finding |\n"
        "|---|---|\n"
        "| **1** | **No consistent CpG trend in RSV.** Lineages disagree: A.D.1 and B.D.E rise, A.D.3 and A.D.5 fall. In the influenza paper all twelve persistent lineages moved the same way. Here they do not, and section 3.6 shows the overall direction is an artefact of lineage definition. |\n"
        "| **2** | **GC3 and GC12 move in opposite directions**, unlike influenza where both fell. Synonymous and non-synonymous sites are under different pressure in RSV, but not in the pattern the paper describes. |\n"
        "| **3** | **Geography is not the explanation.** Slopes keep their direction after region adjustment, so the inconsistency between lineages is real, not a sampling artefact. |\n"
        "| **4** | **All effects are small** — order 10⁻⁴ per year, a fraction of a percentage point across the study period. |\n\n"
        "### How to read this\n\n"
        "RSV has circulated in humans for a very long time and is already well adapted to it. "
        "The influenza result came from lineages that had recently entered a *new* host and "
        "were still adjusting. There is no equivalent transition inside human RSV, so finding "
        "no directional trend is a reasonable outcome rather than a failed analysis.\n\n"
        "The clearest signal in this project is not temporal at all — it is the difference "
        "between genes and between the two subtypes, which notebook 04 covers.")
])

# --------------------------------------------------------------------- 04 by gene
build("notebooks/04_genes.ipynb", [
    new_markdown_cell(
        "# 04 — Composition gene by gene\n\n"
        "---\n\n"
        "Everything so far treated the coding regions as one block, but genes do different "
        "jobs under different pressure. The paper does the same split across the eight "
        "influenza segments.\n\n"
        "RSV's ten genes, in order along the genome:\n\n"
        "`NS1 · NS2 · N · P · M · SH · G · F · M2 · L`\n\n"
        "Four are worth knowing:\n\n"
        "| Gene | Role |\n"
        "|---|---|\n"
        "| **NS1, NS2** | interferon antagonists — they switch off the cell's antiviral alarm |\n"
        "| **G** | attachment protein, surface-exposed, the most variable gene |\n"
        "| **F** | fusion protein, target of nirsevimab and the recent vaccines |\n"
        "| **L** | polymerase, about 43% of the genome, so it dominates any genome-wide average |\n\n"
        "---"),
    new_code_cell(SETUP),
    new_markdown_cell(
        "## 4.1  Composition of each gene"),
    new_code_cell(
        'g = pd.read_csv("../data/processed/rsv_gene_features.tsv", sep="\\t")\n'
        'GENES = ["NS1","NS2","N","P","M","SH","G","F","M2","L"]\n'
        'per_gene = g[g.gene != "GENOME"]\n'
        'tab = (per_gene.groupby(["Subtype", "gene"])\n'
        '       .agg(n=("gc_content","size"), length=("valid_len","median"),\n'
        '            gc=("gc_content","mean"), gc3=("gc3","mean"),\n'
        '            oe_CG=("oe_CG","mean"), oe_GG=("oe_GG","mean")).reset_index())\n'
        'tab["gene"] = pd.Categorical(tab["gene"], GENES, ordered=True)\n'
        'tab.sort_values(["Subtype","gene"]).round(4)'),
    new_markdown_cell(
        "### Patterns\n\n"
        "- **The G gene stands out** — the highest GC content of any gene, well above the "
        "genome-wide figure. It is the most variable gene and the most exposed to antibodies, "
        "so it is the natural region to examine on its own.\n"
        "- **CpG suppression is uneven along the genome** — in RSV-A the two interferon "
        "antagonist genes sit far below everything else."),
    new_code_cell('run_r()'),
    new_code_cell('show("fig5_gene_composition.png")'),
    new_markdown_cell(
        "## 4.2  The NS1 and NS2 contrast\n\n"
        "CpG for every gene, normalised per 1,000 nucleotides so genes of different "
        "length are comparable, with the subtype fold change."),
    new_code_cell(
        'rows = []\n'
        'for gn in GENES:\n'
        '    s = per_gene[per_gene.gene == gn]\n'
        '    a = s.loc[s.Subtype == "A", "freq_CG"] * 1000      # CpG per 1,000 nt\n'
        '    b = s.loc[s.Subtype == "B", "freq_CG"] * 1000\n'
        '    rows.append({"gene": gn, "A per 1k": a.mean(), "B per 1k": b.mean(),\n'
        '                 "fold B/A": b.mean() / a.mean(),\n'
        '                 "A % zero": (a == 0).mean() * 100,\n'
        '                 "B % zero": (b == 0).mean() * 100})\n'
        'pd.DataFrame(rows).set_index("gene").round(2)'),
    new_markdown_cell(
        "### Reading the table honestly\n\n"
        "Looking at all ten genes rather than a chosen few changes what can be claimed.\n\n"
        "- **CpG differs between the subtypes in every gene**, and in *both* directions: "
        "RSV-B is higher in NS1, NS2, P and M; RSV-A is higher in N, SH, G, F, M2 and L. "
        "A subtype difference in CpG is the rule here, not the exception.\n"
        "- **What singles out NS1 and NS2 is the size of the gap** — 12.3x and 7.4x, against "
        "2.2x for the next largest (P). They are extreme, not unique.\n"
        "- **Only NS1 and NS2 reach an absolute floor.** 84.8% of RSV-A genomes carry no CpG "
        "anywhere in NS1 and 34.2% none in NS2. No other gene does this except SH, weakly.\n"
        "- **GC content is matched** in those genes (NS1: 30.5% in A, 30.3% in B), so the "
        "difference is specific to the CpG pair rather than to base composition.\n\n"
        "A note on statistics: with 15,458 against 17,580 genomes, *every* gene returns "
        "p < 1e-300 on a Mann-Whitney test. The p value cannot discriminate at this sample "
        "size, so the fold change and Cliff's delta in "
        "`results/gene_cpg_subtype_test.tsv` are the measures that carry information.\n\n"
        "This is an observation about sequence, not a demonstrated mechanism — testing it "
        "would require laboratory work."),
    new_markdown_cell(
        "## 4.3  Which genes carry the time trend?\n\n"
        "The within-lineage regression from notebook 03, run on each gene separately. The "
        "equivalent of the paper's Figure 2b."),
    new_code_cell(
        'bg = pd.read_csv("../results/slopes_by_gene.tsv", sep="\\t")\n'
        'bg = bg[bg.feature == "oe_CG"]\n'
        'piv = bg.pivot_table(index="gene", columns=["Subtype","lineage_grp"], values="slope")\n'
        'piv = piv.reindex(GENES)\n'
        'piv.round(5)'),
    new_code_cell(
        'fig, ax = plt.subplots(figsize=(9, 3.6))\n'
        'vals = piv.values.astype(float)\n'
        'lim = np.nanmax(np.abs(vals))\n'
        'im = ax.imshow(vals, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")\n'
        'ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index)\n'
        'ax.set_xticks(range(piv.shape[1]))\n'
        'ax.set_xticklabels([f"{a}: {b}" for a, b in piv.columns], rotation=40, ha="right")\n'
        'ax.set_title("CpG O/E change per year, by gene and lineage")\n'
        'fig.colorbar(im, ax=ax, shrink=.85, label="slope per year")\n'
        'plt.tight_layout(); plt.show()'),
    new_markdown_cell(
        "Red means CpG O/E rising in that gene and lineage, blue means falling. Reading "
        "across a row shows whether a gene behaves consistently across lineages."),
])

print("done")
