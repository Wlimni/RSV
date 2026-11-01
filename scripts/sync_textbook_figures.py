"""
Script Name: sync_textbook_figures.py
Purpose:     Copy the analysis figures into docs/figures/ under the short names the
             textbook HTML references, so docs/ stays self-contained and shareable.
Inputs:      results/figures/fig1..fig6 *.png
Outputs:     docs/figures/r1..r6 *.png
Notes:       Run after regenerating figures (Rscript R/figures.R, and
             scripts/plot_subtype_contrast.py). The teaching diagrams t1..t5 are made
             separately by scripts/make_teaching_figures.py and are not touched here.
"""
import shutil
from pathlib import Path

PAIRS = {
    "fig1_lineage_trends.png":     "r1.png",
    "fig2_gc3_vs_gc12.png":        "r2.png",
    "fig3_dinucleotide_slopes.png": "r3.png",
    "fig4_gc_correlation.png":     "r4.png",
    "fig5_gene_composition.png":   "r5.png",
    "fig6_subtype_contrast.png":   "r6.png",
}

src, dst = Path("results/figures"), Path("docs/figures")
dst.mkdir(parents=True, exist_ok=True)
for a, b in PAIRS.items():
    if not (src / a).exists():
        print(f"  missing, skipped: {a}")
        continue
    shutil.copy2(src / a, dst / b)
    print(f"  {a}  ->  docs/figures/{b}")
