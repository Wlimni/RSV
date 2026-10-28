"""
Script Name: make_notebooks.py
Purpose:     Generate the analysis notebooks. Notebooks are the main pipeline and are written
             to be read: each step has a markdown cell explaining what it does and why, and
             each figure is produced by the matching R script and displayed inline.
Inputs:      none
Outputs:     notebooks/01..05 *.ipynb
Notes:       Notebook content is kept here as plain Python so it stays version-controlled and
             regenerable. Run: .venv/bin/python scripts/make_notebooks.py
             Then execute with: jupyter nbconvert --execute --inplace notebooks/*.ipynb
"""
import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

os.makedirs("notebooks", exist_ok=True)

# Shared first cell: imports, paths, and the helper that runs an R script and shows its figure.
SETUP = '''\
import sys, os, subprocess
sys.path.insert(0, os.path.abspath("../src"))
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from IPython.display import Image, display

ROOT = os.path.abspath("..")          # repo root; R scripts expect to run from here
pd.set_option("display.width", 120)

def run_r(script):
    """Run one of the R figure scripts from the repo root and report what it wrote."""
    r = subprocess.run(["Rscript", script], capture_output=True, text=True, cwd=ROOT)
    print(r.stdout.strip() or r.stderr.strip()[-400:])

def show(fig):
    """Display a figure from results/figures/ inline."""
    display(Image(os.path.join(ROOT, "results", "figures", fig)))
'''

DINUC_NOTE = """\
**The four features we measure**, all from the paper:

| Feature | What it counts |
|---|---|
| **GC content** | fraction of all letters that are G or C |
| **CpG** | a C immediately followed by a G |
| **GpC** | a G followed by a C |
| **GpG** | a G followed by a G |
| **CpC** | a C followed by a C |

For each dinucleotide we also compute the **O/E ratio** (observed / expected):

> O/E = frequency(XY) / (frequency(X) x frequency(Y))

O/E near **1.0** means the pair occurs about as often as chance predicts. Below **1.0** means
the virus has fewer of them than expected - the pair is *suppressed*. This matters because the
human antiviral protein ZAP recognises viral RNA at CpG sites, so viruses that circulate in
humans long-term tend to carry fewer CpGs than chance would give them.
"""


def build(path, cells):
    nb = new_notebook()
    nb.cells = cells
    nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3",
                                 "language": "python"}
    nbf.write(nb, path)
    print("wrote", path)


# ---------------------------------------------------------------- 01 data overview
build("notebooks/01_data_overview.ipynb", [
    new_markdown_cell(
        "# 01 - What data do we have?\n\n"
        "**Goal of the project:** the paper (Ye et al. 2026, *Nature Microbiology*) showed that\n"
        "influenza A lineages which became established in mammals had lower GC content and fewer\n"
        "GC-related dinucleotides than those that did not. We are testing whether similar genome\n"
        "composition features show any pattern in RSV.\n\n"
        "This notebook just looks at the dataset: how many sequences, from when, from where, and\n"
        "what the quality filter removes."),
    new_code_cell(SETUP),
    new_markdown_cell(
        "## The sequences\n\n"
        "`rsv_features.tsv` is produced by `scripts/build_features.py`. One row per genome, with\n"
        "the GISAID metadata plus the composition features already computed."),
    new_code_cell(
        'd = pd.read_csv("../data/processed/rsv_features.tsv", sep="\\t")\n'
        'print("sequences:", len(d))\n'
        'print("subtypes:", d.Subtype.value_counts().to_dict())\n'
        'print("years:", d.year.min(), "to", d.year.max())\n'
        '# accession IDs are dropped from displays so no GISAID record is stored in the notebook\n'
        'd.drop(columns=["Accession ID"]).head(3)'),
    new_markdown_cell(
        "## How the QC filter was chosen\n\n"
        "RSV genomes are about 15,200 nt. We keep sequences between **14,900 and 15,400 nt** and\n"
        "with at most **1% ambiguous bases** (N or gaps), because ambiguous letters would distort\n"
        "the dinucleotide counts. The histogram below shows why that window is sensible: almost\n"
        "everything sits in a tight peak."),
    new_code_cell(
        'fig, ax = plt.subplots(figsize=(7, 3.5))\n'
        'for s, c in zip("AB", ["#2c7fb8", "#de2d26"]):\n'
        '    ax.hist(d.loc[d.Subtype == s, "seq_len_actual"], bins=100, alpha=.55, label=f"RSV-{s}", color=c)\n'
        'ax.axvline(14900, ls="--", c="k", lw=1); ax.axvline(15400, ls="--", c="k", lw=1)\n'
        'ax.set_xlim(14700, 15600); ax.set_xlabel("sequence length (nt)"); ax.set_ylabel("sequences")\n'
        'ax.legend(); ax.set_title("Genome length, dashed lines = filter window")\n'
        'plt.tight_layout(); plt.show()'),
    new_markdown_cell(
        "## When and where the sequences come from\n\n"
        "Two things to notice, because they affect every trend later:\n\n"
        "1. Sequencing ramps up sharply after 2018, and dips in 2020-21 (COVID-19 disrupted RSV\n"
        "   circulation and surveillance).\n"
        "2. Europe contributes far more sequences than anywhere else. Any trend over time is\n"
        "   partly a trend in *who was sequencing*, not only in the virus."),
    new_code_cell(
        'fig, ax = plt.subplots(1, 2, figsize=(12, 3.5))\n'
        'yr = d.groupby(["year", "Subtype"]).size().unstack(fill_value=0)\n'
        'yr.plot.bar(stacked=True, ax=ax[0], color=["#2c7fb8", "#de2d26"], width=.85)\n'
        'ax[0].set_title("sequences per year"); ax[0].set_xlabel("")\n'
        'reg = d.groupby(["region", "Subtype"]).size().unstack(fill_value=0).sort_values("A")\n'
        'reg.plot.barh(stacked=True, ax=ax[1], color=["#2c7fb8", "#de2d26"])\n'
        'ax[1].set_title("sequences per region"); ax[1].set_ylabel("")\n'
        'plt.tight_layout(); plt.show()'),
    new_markdown_cell(
        "## Lineages\n\n"
        "GISAID assigns each genome to a **lineage** - a named branch of the RSV family tree.\n"
        "Lineages rise and fall over time, which becomes important in notebook 04."),
    new_code_cell(
        'for s in "AB":\n'
        '    top = d[d.Subtype == s].Lineage.value_counts().head(5)\n'
        '    print(f"RSV-{s} most common lineages:")\n'
        '    print(top.to_string(), "\\n")'),
])

# ------------------------------------------------------- 02 whole-genome composition
build("notebooks/02_composition.ipynb", [
    new_markdown_cell(
        "# 02 - Genome composition: the basic measurements\n\n" + DINUC_NOTE),
    new_code_cell(SETUP),
    new_code_cell(
        'd = pd.read_csv("../data/processed/rsv_features.tsv", sep="\\t")\n'
        'FEATS = ["gc_content", "oe_CG", "oe_GC", "oe_GG", "oe_CC"]\n'
        'd.groupby("Subtype")[FEATS].mean().round(4)'),
    new_markdown_cell(
        "### Reading that table\n\n"
        "**CpG O/E is about 0.24 in both subtypes.** That is strong suppression: RSV carries only\n"
        "a quarter of the CpG dinucleotides that its own base composition would predict. The other\n"
        "three dinucleotides sit near or above 1.0, so this is specific to CpG rather than a\n"
        "general property of the genome.\n\n"
        "RSV-A and RSV-B differ, but the differences in GC content are small (0.3344 vs 0.3364).\n"
        "The dinucleotide ratios separate them more clearly."),
    new_markdown_cell(
        "## Figures\n\n"
        "These are produced by `R/figures.R`, which reads the same feature table. Running it here\n"
        "keeps the notebook and the files in `results/figures/` in step."),
    new_code_cell('run_r("R/figures.R")'),
    new_code_cell('show("fig1_gc_by_subtype.png")'),
    new_code_cell('show("fig2_oe_dinuc_by_subtype.png")'),
    new_markdown_cell(
        "The dashed line in the second figure marks O/E = 1. Only CpG sits far below it. GpC, GpG\n"
        "and CpC hover around or above 1, so whatever is removing CpG is not simply removing G and\n"
        "C letters in general."),
    new_markdown_cell(
        "## Does composition change over the years?"),
    new_code_cell('show("fig3_yearly_trends.png")'),
    new_markdown_cell(
        "There is movement, and RSV-A and RSV-B seem to go in different directions after 2021.\n"
        "**Be careful with this plot.** Notebook 04 shows that most of this apparent trend is not\n"
        "the virus changing - it is one lineage replacing another. We come back to it there."),
])

# ------------------------------------------------------------------ 03 gene regions
build("notebooks/03_gene_regions.ipynb", [
    new_markdown_cell(
        "# 03 - Looking at specific regions of the genome\n\n"
        "So far every measurement used the whole genome at once. But a genome is not uniform: it\n"
        "is made of genes that do different jobs and are under different pressure. The paper did\n"
        "the same thing, measuring each influenza segment separately.\n\n"
        "RSV has **10 genes**, in this order along the genome:\n\n"
        "`NS1 - NS2 - N - P - M - SH - G - F - M2 - L`\n\n"
        "Two are worth knowing by name:\n\n"
        "- **G** (attachment glycoprotein) - sits on the virus surface, the most variable gene,\n"
        "  under constant pressure from antibodies.\n"
        "- **F** (fusion protein) - also on the surface, and the target of the new preventives\n"
        "  (nirsevimab, the RSVpreF vaccine).\n\n"
        "`scripts/build_gene_features.py` cuts every genome into these 10 genes by aligning it to\n"
        "a reference genome with MAFFT, then measures composition in each one separately."),
    new_code_cell(SETUP),
    new_code_cell(
        'g = pd.read_csv("../data/processed/rsv_gene_features.tsv", sep="\\t")\n'
        'GENES = ["NS1","NS2","N","P","M","SH","G","F","M2","L"]\n'
        'print(len(g), "rows =", g["Accession ID"].nunique(), "genomes x", g.gene.nunique(), "genes")\n'
        'g.drop(columns=["Accession ID"]).head(3)'),
    new_markdown_cell(
        "## Composition gene by gene\n\n"
        "This is the main table of this notebook. `freq_CG`, `freq_GC` and `freq_GG` are the raw\n"
        "frequencies of those dinucleotides inside each gene."),
    new_code_cell(
        'tab = (g.groupby(["Subtype","gene"])\n'
        '        .agg(n=("gc_content","size"), length=("valid_len","median"),\n'
        '             gc=("gc_content","mean"),\n'
        '             freq_CG=("freq_CG","mean"), freq_GC=("freq_GC","mean"), freq_GG=("freq_GG","mean"),\n'
        '             oe_CG=("oe_CG","mean"), oe_GG=("oe_GG","mean"))\n'
        '        .reset_index())\n'
        'tab["gene"] = pd.Categorical(tab["gene"], GENES, ordered=True)\n'
        'tab = tab.sort_values(["Subtype","gene"])\n'
        'tab.round(4)'),
    new_markdown_cell(
        "### Patterns to notice\n\n"
        "1. **The G gene stands out.** It has the highest GC content of any gene (0.43 in RSV-A vs\n"
        "   about 0.33 genome-wide) and the highest GpG O/E. It is the most variable gene, so this\n"
        "   is a region worth focusing on.\n"
        "2. **CpG suppression is not uniform along the genome.** Some genes are far more suppressed\n"
        "   than others.\n"
        "3. **NS1 and NS2 in RSV-A have almost no CpG at all** - look at how small `freq_CG` is for\n"
        "   those two rows compared with RSV-B. We look at this directly below."),
    new_code_cell('run_r("R/figures_genes.R")'),
    new_code_cell('show("fig9_gene_composition.png")'),
    new_markdown_cell(
        "## Zooming in: the G gene\n\n"
        "Since G is the most distinctive region, here is its composition on its own, over time and\n"
        "by subtype. This is the kind of focused, single-region view that is easy to interpret."),
    new_code_cell(
        'G = g[g.gene == "G"]\n'
        'yr = (G[(G.year >= 2015) & (G.year <= 2026)]\n'
        '      .groupby(["Subtype","year"])\n'
        '      .agg(n=("gc_content","size"), gc=("gc_content","mean"),\n'
        '           freq_CG=("freq_CG","mean"), freq_GC=("freq_GC","mean"), freq_GG=("freq_GG","mean"))\n'
        '      .reset_index())\n'
        'yr = yr[yr.n >= 20]\n'
        'fig, ax = plt.subplots(1, 3, figsize=(13, 3.4), sharex=True)\n'
        'for i, col in enumerate(["freq_CG","freq_GC","freq_GG"]):\n'
        '    for s, c in zip("AB", ["#2c7fb8","#de2d26"]):\n'
        '        x = yr[yr.Subtype == s]\n'
        '        ax[i].plot(x.year, x[col], "-o", color=c, ms=4, label=f"RSV-{s}")\n'
        '    ax[i].set_title(f"G gene: {col}"); ax[i].set_xlabel("year")\n'
        'ax[0].legend(); plt.tight_layout(); plt.show()'),
    new_markdown_cell(
        "## The clearest pattern in the whole project\n\n"
        "NS1 and NS2 are RSV's **interferon antagonist** genes - their job is to switch off the\n"
        "cell's antiviral alarm. In RSV-A these two genes have had almost all their CpG removed."),
    new_code_cell(
        'cmp = (g[g.gene.isin(["NS1","NS2","N","L"])]\n'
        '       .groupby(["gene","Subtype"])\n'
        '       .apply(lambda s: pd.Series({\n'
        '           "CpG per gene": (s.freq_CG * s.valid_len).mean(),\n'
        '           "% genomes with zero CpG": (s.freq_CG == 0).mean() * 100,\n'
        '           "GC content": s.gc_content.mean()}), include_groups=False)\n'
        '       .round(3))\n'
        'cmp'),
    new_markdown_cell(
        "**Read the NS1 rows.** RSV-A averages about 0.16 CpG in the whole gene and roughly 85% of\n"
        "RSV-A genomes contain **none at all**, while RSV-B has around 1.9 and almost always has\n"
        "some. GC content is nearly identical between them (0.305 vs 0.303), so this is not a side\n"
        "effect of base composition - it is specific to CpG. N and L, shown for comparison, behave\n"
        "the same in both subtypes.\n\n"
        "This is an observation, not a proven mechanism. It is consistent with ZAP pressure, but\n"
        "showing that would need laboratory work."),
    new_code_cell('show("fig11_cpg_depletion.png")'),
])

# ----------------------------------------------------------------- 04 trend caution
build("notebooks/04_trends.ipynb", [
    new_markdown_cell(
        "# 04 - Trends over time, and why the obvious reading is wrong\n\n"
        "Notebook 02 showed RSV-B's CpG O/E apparently rising year on year. The tempting conclusion\n"
        "is that the virus is evolving towards more CpG. This notebook shows that conclusion does\n"
        "not hold, using a simple idea.\n\n"
        "## The idea\n\n"
        "An average can move for two completely different reasons:\n\n"
        "1. **Within-lineage change** - each lineage's genomes are actually changing.\n"
        "2. **Replacement** - no lineage changes at all, but a lineage that already had high CpG\n"
        "   grows while a lineage with low CpG shrinks. The average rises without any evolution.\n\n"
        "To separate them we fit the trend twice: once plainly, and once while holding lineage\n"
        "constant. The second fit cannot be affected by replacement, so whatever trend survives it\n"
        "is real change. (This is a classic statistical trap called **Simpson's paradox**.)"),
    new_code_cell(SETUP + "from rsv import trends"),
    new_code_cell(
        'd = pd.read_csv("../data/processed/rsv_features.tsv", sep="\\t")\n'
        'd = d[(d.year >= 2015) & (d.year <= 2026)]\n'
        'res = pd.concat([trends.decompose_trend(s, "oe_CG").assign(Subtype=sub)\n'
        '                 for sub, s in d.groupby("Subtype")], ignore_index=True)\n'
        'res[["Subtype","feature","slope_total","slope_within","pct_from_turnover"]].round(6)'),
    new_markdown_cell(
        "### What this says\n\n"
        "For **RSV-B** the two slopes have **opposite signs**: the plain trend is positive, but\n"
        "once lineage is held constant it is negative. So the rise in notebook 02 came from\n"
        "lineage replacement. Inside lineages, CpG O/E is going *down* - which is the same\n"
        "direction the paper found for influenza.\n\n"
        "For **RSV-A** both slopes are negative, so there is no contradiction there."),
    new_code_cell('run_r("R/figures_trends.R")'),
    new_code_cell('show("fig5_within_vs_between.png")'),
    new_markdown_cell(
        "## Checking it is not an artefact of sampling\n\n"
        "Because sequencing is Europe-heavy and 2020-21 was disrupted, we repeat the same fit\n"
        "while removing those influences."),
    new_code_cell(
        'rob = trends.robustness(d, "oe_CG")\n'
        'rob[["Subtype","scenario","n","slope_within","p_within"]].round(6)'),
    new_markdown_cell(
        "**RSV-B's within-lineage decline survives every version.** RSV-A's does not: it disappears\n"
        "when the pandemic years are dropped (p = 0.97) and when regions are balanced (p = 0.16).\n"
        "So we trust the RSV-B pattern and treat the RSV-A one as a sampling artefact."),
    new_code_cell('show("fig8_robustness.png")'),
    new_markdown_cell(
        "## Summary of this notebook\n\n"
        "| Claim | Verdict |\n"
        "|---|---|\n"
        "| RSV-B CpG O/E rises over time | True of the raw average, but caused by lineage replacement |\n"
        "| RSV-B genomes are gaining CpG | **False** - within lineages it falls |\n"
        "| RSV-A shows a within-lineage trend | Not supported once sampling is controlled |\n\n"
        "The practical lesson: with this kind of data, always check whether a trend is the thing\n"
        "changing or the mix changing."),
])

# ---------------------------------------------------------------------- 05 summary
build("notebooks/05_summary.ipynb", [
    new_markdown_cell(
        "# 05 - What we found, and the limit of this dataset\n\n"
        "A short recap of the patterns, and one honest statement about what RSV data can and\n"
        "cannot answer."),
    new_code_cell(SETUP),
    new_markdown_cell(
        "## Patterns observed\n\n"
        "1. **RSV strongly suppresses CpG** genome-wide (O/E about 0.24), while the other three\n"
        "   GC dinucleotides sit near 1. So CpG is specifically avoided.\n"
        "2. **Composition varies a lot between genes.** The G gene is the most GC-rich and most\n"
        "   distinctive region; it is the natural region to focus on.\n"
        "3. **RSV-A has removed nearly all CpG from NS1 and NS2**, its two interferon-antagonist\n"
        "   genes, at the same GC content as RSV-B. This is the clearest pattern we found.\n"
        "4. **The apparent rise in CpG over time in RSV-B is lineage replacement, not evolution.**\n"
        "   Within lineages it declines, matching the paper's direction.\n"
        "5. Comparisons at the level of individual lineages were **too few to be conclusive** and\n"
        "   did not survive multiple-testing correction."),
    new_code_cell(
        'show("fig11_cpg_depletion.png")'),
    new_markdown_cell(
        "## The limit of this dataset\n\n"
        "The paper's central comparison is between viruses that achieved **sustained** transmission\n"
        "in a host and those that only caused **sporadic** infections. Those labels come from\n"
        "knowing each virus's history, not from its sequence.\n\n"
        "**Every GISAID RSV sequence is from a human, and RSV has circulated in humans for\n"
        "decades.** So this dataset contains only one of the two classes. No model can learn to\n"
        "separate sustained from sporadic here, because there are no sporadic examples.\n\n"
        "That does not make the work useless: the pipeline applies to any virus, and finding 3\n"
        "stands on its own. But answering the original question needs a dataset with both classes,\n"
        "such as a panel of paramyxoviruses where Nipah and Hendra supply the sporadic side.\n"
        "See `docs/PROJECT_BRIEF.md`."),
])

print("done")
