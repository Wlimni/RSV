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
    nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3",
                                 "language": "python"}
    nbf.write(nb, path)
    print("wrote", path)


# ----------------------------------------------------------------- 01 the data
build("notebooks/01_data.ipynb", [
    new_markdown_cell(
        "# 01 — The data\n\n"
        "---\n\n"
        "### What this project asks\n\n"
        "Ye et al. 2026 (*Nature Microbiology*) measured GC content and GC-containing "
        "dinucleotide frequencies across 115,520 influenza A genomes, and found that these "
        "features **declined over time** within lineages that had achieved sustained "
        "transmission in mammals.\n\n"
        "We apply the same measurements to human RSV and ask:\n\n"
        "> **Within each RSV lineage, how does genome composition change over time, and does "
        "the pattern differ between RSV-A and RSV-B?**\n\n"
        "This notebook covers the dataset itself: how much there is, when and where it came "
        "from, and what quality control removes.\n\n"
        "---"),
    new_code_cell(SETUP),
    new_markdown_cell(
        "## 1.1  The sequences\n\n"
        "`rsv_features.tsv` is written by `scripts/build_features.py`: one row per genome, "
        "carrying the GISAID metadata alongside composition features."),
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
        "## 1.3  When and where the sequences come from\n\n"
        "Two features of the sampling matter for everything that follows.\n\n"
        "**Sequencing is not evenly spread over time.** It rises steeply from 2018 and dips "
        "in 2020–21, when COVID-19 disrupted both RSV circulation and surveillance.\n\n"
        "**Sequencing is not evenly spread over the world.** Europe contributes far more "
        "genomes than anywhere else.\n\n"
        "Any trend measured over calendar time is therefore partly a trend in *who was "
        "sequencing*. This is the main reason the analysis in notebook 03 is done inside "
        "lineages rather than across the whole dataset."),
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
        "Four, found by auditing the data before analysing it. Two of them change what we "
        "do.\n\n"
        "---\n\n"
        "### (a) 2026 is incomplete, so we exclude it\n\n"
        "The download was made in October 2026, but sequences take time to be deposited. "
        "The later months of 2026 are nearly empty, far below even the seasonal low."),
    new_code_cell(
        'm26 = d[d.year == 2026].date.str.slice(5, 7).value_counts().sort_index()\n'
        'allyr = d.date.str.slice(5, 7).value_counts().sort_index()\n'
        'cmp = pd.DataFrame({"2026": m26, "all years (mean/yr)": (allyr / 12).round(0)})\n'
        'print(cmp.fillna(0).astype(int).to_string())\n'
        'print("\\nlast collection date:", d.date.max())'),
    new_markdown_cell(
        "September 2026 has a couple of dozen genomes where a normal September has over a "
        "thousand. That is deposition lag, not biology. Keeping 2026 would add a tail of "
        "points bunched into the first half of the year, contributed only by laboratories "
        "that deposit quickly.\n\n"
        "**Decision: the analysis window is 2015-01-01 to 2025-12-31.**\n\n"
        "---\n\n"
        "### (b) Where sequencing happens has shifted completely\n\n"
        "This is the more serious problem."),
    new_code_cell(
        'reg = pd.crosstab(d.year, d.region, normalize="index").mul(100).round(1)\n'
        'reg.loc[2015:2026]'),
    new_markdown_cell(
        "Oceania supplied 38% of genomes in 2015 and 4% in 2026; Europe went from 18% to "
        "56%. North America contributes nothing at all in 2017.\n\n"
        "So **calendar time and geography are entangled**. A composition trend over time "
        "could partly be a trend in which part of the world was sequencing.\n\n"
        "Worse, this drift happens *inside* individual lineages too, so analysing within "
        "lineages does not remove it. Lineage A.D goes from 71% Oceania in its early period "
        "to 36% Europe in its late period.\n\n"
        "**Decision: region is included as a covariate in the regressions, and results are "
        "reported with and without it** (notebook 03).\n\n"
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
        "December and November carry three to four times September. Northern and southern "
        "hemisphere epidemics also peak at opposite times, so season and geography are "
        "partly entangled as well. For a slope measured across a decade this is minor, so "
        "we record it rather than correct it.\n\n"
        "---\n\n"
        "### (d) Before 2015 there is too little data — already excluded\n\n"
        "Everything before 2015 is about 3% of the dataset, and **RSV-A has essentially no "
        "sequences before 2010**, so the two subtypes do not even start at the same time. "
        "The 2015 lower bound handles this."),
    new_markdown_cell(
        "## 1.5  Lineages\n\n"
        "GISAID assigns each genome to a **lineage**, a named branch of the RSV family tree. "
        "The labels are hierarchical and dot-separated, so `A.D.1.5` sits inside `A.D.1`, "
        "which sits inside `A.D`.\n\n"
        "Lineages rise and fall over a few seasons. Notebook 03 explains why that matters "
        "and how we choose the level to work at."),
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
        "### GC content\n"
        "The fraction of letters that are G or C.\n\n"
        "### Dinucleotide frequency\n"
        "How often a given pair of adjacent letters occurs. We track the four GC-containing "
        "pairs used in the paper, plus one more:\n\n"
        "| Pair | Written | Note |\n"
        "|---|---|---|\n"
        "| C then G | CpG | the motif of interest |\n"
        "| G then C | GpC | |\n"
        "| G then G | GpG | |\n"
        "| C then C | CpC | |\n"
        "| T then A | UpA | **not in the paper**; added as exploratory |\n\n"
        "### Observed / expected ratio\n"
        "Counting pairs directly is misleading, because a GC-rich genome contains more CpG "
        "by chance alone. So we divide by what chance predicts:\n\n"
        "$$ \\mathrm{O/E}(\\mathrm{CpG}) = \\frac{f(\\mathrm{CpG})}{f(\\mathrm{C}) \\times f(\\mathrm{G})} $$\n\n"
        "A value near **1** means as common as chance predicts; **below 1** means the virus "
        "is avoiding that pair.\n\n"
        "### GC3 and GC12\n"
        "GC content split by position inside the codon. **GC3** is the third position, where "
        "most changes are synonymous and do not alter the protein. **GC12** is the first and "
        "second, where most changes do alter it.\n\n"
        "Comparing how fast the two move separates pressure acting on the RNA from pressure "
        "acting on the protein — the same logic as comparing dS with dN. The paper reports "
        "GC3 changing **5.65× faster** than GC12 in influenza.\n\n"
        "---"),
    new_code_cell(SETUP + "from rsv import composition"),
    new_markdown_cell(
        "## 2.1  The functions, on a short example\n\n"
        "All measurements come from `src/rsv/composition.py`. Ambiguous positions are "
        "skipped rather than guessed."),
    new_code_cell(
        'demo = "ATGGCGTTAACGCGTAAATTACGATAG"\n'
        'f = composition.all_features(demo, coding=True)\n'
        'for k in ["gc_content", "gc3", "gc12", "per1k_CG", "oe_CG", "oe_TA"]:\n'
        '    print(f"{k:12s} {f[k]:.4f}")'),
    new_markdown_cell(
        "## 2.2  The measurements across all genomes\n\n"
        "Here `GENOME` means the ten coding regions concatenated, which is how the paper "
        "defines genome-level GC, GC3 and GC12 — it keeps the reading frame well defined."),
    new_code_cell(
        'g = pd.read_csv("../data/processed/rsv_gene_features.tsv", sep="\\t")\n'
        'genome = g[g.gene == "GENOME"]\n'
        'cols = ["gc_content", "gc3", "gc12", "oe_CG", "oe_GC", "oe_GG", "oe_CC", "oe_TA"]\n'
        'genome.groupby("Subtype")[cols].mean().round(4)'),
    new_markdown_cell(
        "### Reading that table\n\n"
        "**CpG O/E sits far below 1 in both subtypes.** RSV carries only a fraction of the "
        "CpG pairs its base composition would predict. The other GC pairs sit near or above "
        "1, so the effect is specific to CpG rather than a general shortage of G and C.\n\n"
        "**GC3 is much lower than GC12**, which is typical: third positions are free to drift "
        "because changing them usually does not change the protein."),
    new_code_cell(
        'fig, ax = plt.subplots(1, 2, figsize=(11, 3.4))\n'
        'for s in "AB":\n'
        '    x = genome[genome.Subtype == s]\n'
        '    ax[0].hist(x.gc_content, bins=70, alpha=.6, color=PAL[s], label=f"RSV-{s}")\n'
        'ax[0].set_xlabel("GC content"); ax[0].set_ylabel("genomes")\n'
        'ax[0].set_title("GC content"); ax[0].legend(frameon=False)\n'
        'oe = ["oe_CG", "oe_GC", "oe_GG", "oe_CC", "oe_TA"]\n'
        'names = ["CpG", "GpC", "GpG", "CpC", "UpA"]\n'
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
        "Only **CpG**, and to a lesser extent **UpA**, fall clearly below the dashed line. "
        "That is the pattern reported across RNA viruses that circulate in humans."),
])

# ------------------------------------------------- 03 the core: time within lineage
build("notebooks/03_lineage_trends.ipynb", [
    new_markdown_cell(
        "# 03 — Composition over time, within lineages\n\n"
        "---\n\n"
        "**This is the core analysis.** It follows Figures 2 and 3 of the paper.\n\n"
        "### Why fit inside lineages\n\n"
        "An average across the whole dataset can move for two different reasons:\n\n"
        "1. the genomes inside each lineage are genuinely changing, or\n"
        "2. nothing changes, but a lineage that already had a high value grows while one "
        "with a low value shrinks.\n\n"
        "The second produces a confident-looking trend with no evolution behind it. Fitting "
        "**separately inside each lineage** removes it, which is what the paper does.\n\n"
        "### Choosing the lineage level\n\n"
        "Full GISAID labels are too fine. Most span only about 6 years with a few hundred "
        "genomes, and a regression slope's uncertainty grows as both the time span and the "
        "sample size shrink. Truncating the label to **three levels** gives roughly 8 "
        "well-populated lineages spanning 2015–2026 — close to the 12 the paper used.\n\n"
        "---"),
    new_code_cell(SETUP + "from rsv import lineage"),
    new_markdown_cell(
        "## 3.1  The lineages we analyse\n\n"
        "Inclusion requires **at least 5 distinct years and at least 200 genomes**."),
    new_code_cell(
        'g = pd.read_csv("../data/processed/rsv_gene_features.tsv", sep="\\t")\n'
        'dates = pd.read_csv("../data/processed/rsv_features.tsv", sep="\\t",\n'
        '                    usecols=["Accession ID", "date"])\n'
        'genome = g[(g.gene == "GENOME") & g.year.between(2015, 2026)].merge(dates, on="Accession ID")\n'
        'gdf, keep = lineage.prepare(genome)\n'
        'print(f"{len(keep)} lineages, {len(gdf):,} genomes\\n")\n'
        'keep.round(2)'),
    new_markdown_cell(
        "Note that all eight lineages start between 2015 and 2017 and run to 2025–2026. "
        "Because they cover essentially the same calendar window, differences between them "
        "are not confounded with *when* they happened to circulate.\n\n"
        "RSV-B contributes only two lineages. They are large, so the estimates are sound, "
        "but any A-versus-B statement rests on two groups rather than a broad survey."),
    new_markdown_cell(
        "## 3.2  Trend of each feature within each lineage\n\n"
        "Each feature is regressed against the **decimal collection date**, inside each "
        "lineage separately. The slope is the change per year."),
    new_code_cell(
        'FEATS = ["gc_content", "gc3", "gc12",\n'
        '         "per1k_CG", "per1k_GC", "per1k_GG", "per1k_CC", "per1k_TA",\n'
        '         "oe_CG", "oe_GC", "oe_GG", "oe_CC", "oe_TA"]\n'
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
        "This is the sharpest test available. If composition is being pushed by something "
        "acting on the RNA, the effect should be strongest at **third codon positions**, "
        "where changes are mostly synonymous. If it reflects protein-level constraint, the "
        "first and second positions should move comparably."),
    new_code_cell(
        'piv = (slopes[slopes.feature.isin(["gc_content", "gc3", "gc12"])]\n'
        '       .pivot_table(index=["Subtype", "lineage_grp"], columns="feature", values="slope"))\n'
        'piv["gc3 / gc12"] = (piv.gc3 / piv.gc12).round(2)\n'
        'piv.round(6)'),
    new_code_cell('show("fig2_gc3_vs_gc12.png")'),
    new_markdown_cell(
        "### What this shows, and how it differs from the paper\n\n"
        "In influenza, GC content fell at **both** site classes, with GC3 falling about "
        "5.65 times faster. That is a coherent picture: composition drifting down, fastest "
        "where changes are invisible to the protein.\n\n"
        "RSV does not look like that. **GC3 rises in 6 of 8 lineages while GC12 falls in 6 "
        "of 8** — the two site classes move in *opposite* directions, so the ratio between "
        "them is negative rather than a clean multiple.\n\n"
        "Two things follow. The synonymous and non-synonymous sites are clearly under "
        "different pressure, which is interesting in itself. But this is not the influenza "
        "signature, and it should not be described as reproducing it.\n\n"
        "Note the size: the slopes are around 10^-4 per year, so roughly a tenth of a "
        "percentage point of GC across the whole decade. Small."),
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
        "them. The two that move most are the lineages flagged in notebook 01 as having the "
        "largest geographic drift, which is what you would hope: the adjustment bites "
        "exactly where the problem is."),
    new_markdown_cell(
        "## 3.5  Is the dinucleotide signal just GC content moving?\n\n"
        "The paper's Figure 3b. If a dinucleotide's frequency correlates almost perfectly "
        "with GC content, it carries no information beyond GC. Correlations well below 1 "
        "mean the dinucleotide is moving on its own."),
    new_code_cell(
        'corr = lineage.gc_vs_dinucleotide_correlation(\n'
        '    gdf, ["per1k_CG", "per1k_GC", "per1k_GG", "per1k_CC", "per1k_TA"])\n'
        'corr.groupby("feature").r.agg(["mean", "min", "max"]).round(3)'),
    new_code_cell('show("fig4_gc_correlation.png")'),
    new_markdown_cell(
        "## 3.6  Sensitivity: does the lineage level drive the result?\n\n"
        "Truncating to three levels was a choice. Repeating at level 4 gives more, smaller "
        "lineages. If the direction and rough size of the slopes agree, the choice is not "
        "driving the conclusion."),
    new_code_cell(
        'gdf4, keep4 = lineage.prepare(genome, level=4)\n'
        's4 = lineage.lineage_slopes(gdf4, ["gc_content", "gc3", "oe_CG"])\n'
        'cmp = pd.concat([slopes[slopes.feature.isin(["gc_content","gc3","oe_CG"])].assign(level=3),\n'
        '                 s4.assign(level=4)])\n'
        'cmp.groupby(["level", "feature"]).slope.agg(\n'
        '    median="median", lineages="size",\n'
        '    n_negative=lambda s: int((s < 0).sum())).round(6)'),
    new_markdown_cell(
        "**The CpG result does not survive this check.** At three levels the median slope is "
        "slightly positive with 3 of 8 lineages negative; at four levels it is negative with "
        "14 of 19 negative. The direction flips depending on how finely lineages are "
        "defined, which means the granularity is driving it rather than the biology.\n\n"
        "GC3 is steadier: positive at both levels, though weaker at the finer one. GC "
        "content overall sits near zero either way."),
    new_markdown_cell(
        "---\n\n"
        "## 3.7  What this notebook found\n\n"
        "Stated plainly, because the result is mostly negative.\n\n"
        "**1. There is no consistent CpG trend in RSV.** Lineages disagree with each other: "
        "A.D.1 and B.D.E rise, A.D.3 and A.D.5 fall. In the influenza paper all twelve "
        "persistent lineages moved the same way. Here they do not, and the sensitivity "
        "check above shows the overall direction is an artefact of lineage definition.\n\n"
        "**2. GC3 and GC12 move in opposite directions**, unlike influenza where both fell. "
        "Synonymous and non-synonymous sites are under different pressure in RSV, but the "
        "pattern is not the one the paper describes.\n\n"
        "**3. Geography is not the explanation.** Slopes keep their direction after region "
        "adjustment, so the inconsistency between lineages is real, not a sampling artefact.\n\n"
        "**4. All effects are small** — of the order of 10^-4 per year, a fraction of a "
        "percentage point across the study period.\n\n"
        "### How to read this\n\n"
        "RSV has circulated in humans for a very long time and is already well adapted to "
        "it. The influenza result came from lineages that had recently entered a *new* host "
        "and were still adjusting. There is no equivalent transition inside human RSV, so "
        "finding no directional trend is a reasonable outcome rather than a failed "
        "analysis.\n\n"
        "The clearest signal in this project is not temporal at all. It is the difference "
        "between genes, and between the two subtypes, which notebook 04 covers.")
])

# --------------------------------------------------------------------- 04 by gene
build("notebooks/04_genes.ipynb", [
    new_markdown_cell(
        "# 04 — Composition gene by gene\n\n"
        "---\n\n"
        "Everything so far treated the coding regions as one block. But a genome is not "
        "uniform: its genes do different jobs and are under different pressure. The paper "
        "does the same split across the eight influenza segments.\n\n"
        "RSV has ten genes, in this order along the genome:\n\n"
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
        "**The G gene stands out** — the highest GC content of any gene, well above the "
        "genome-wide figure. It is the most variable gene and the most exposed to antibodies, "
        "so it is the natural region to examine on its own.\n\n"
        "**CpG suppression is uneven along the genome**, and in RSV-A the two interferon "
        "antagonist genes sit far below everything else."),
    new_code_cell('run_r()'),
    new_code_cell('show("fig5_gene_composition.png")'),
    new_markdown_cell(
        "## 4.2  The NS1 and NS2 contrast\n\n"
        "Taking the two interferon antagonist genes directly, with two ordinary genes for "
        "comparison."),
    new_code_cell(
        'sel = per_gene[per_gene.gene.isin(["NS1","NS2","N","L"])]\n'
        'cmp = (sel.groupby(["gene","Subtype"])\n'
        '       .apply(lambda s: pd.Series({\n'
        '           "CpG per gene": (s.freq_CG * s.valid_len).mean(),\n'
        '           "% genomes with zero CpG": (s.freq_CG == 0).mean() * 100,\n'
        '           "GC content": s.gc_content.mean()}), include_groups=False)\n'
        '       .round(3))\n'
        'cmp'),
    new_markdown_cell(
        "**Read the NS1 rows.** In RSV-A the gene carries almost no CpG at all, and in most "
        "RSV-A genomes there is none anywhere in it. RSV-B's copy of the same gene carries "
        "several. GC content is nearly identical between the two, so this is not a side "
        "effect of base composition — it is specific to the CpG pair.\n\n"
        "N and L, shown for comparison, behave the same way in both subtypes, so this is not "
        "a general difference between A and B.\n\n"
        "This is an observation about sequence, not a demonstrated mechanism. Testing it "
        "would require laboratory work."),
    new_markdown_cell(
        "## 4.3  Which genes carry the time trend?\n\n"
        "The same within-lineage regression from notebook 03, run on each gene separately. "
        "The equivalent of the paper's Figure 2b."),
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
