"""
Script Name: run_trend_analysis.py
Purpose:     Run steps 1, 3 and 5 of the analysis plan — within/between-lineage decomposition
             of composition trends, lineage growth rates vs composition, and robustness checks.
Inputs:      data/processed/rsv_features.tsv
Outputs:     results/decomposition.tsv, results/per_lineage_slopes.tsv,
             results/lineage_growth.tsv, results/growth_vs_composition.tsv,
             results/robustness.tsv
Notes:       Analyses are run per subtype (A and B separately) because their trends differ.
             Run: .venv/bin/python scripts/run_trend_analysis.py
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import pandas as pd
from rsv import trends

FEATURES = ["gc_content", "oe_CG", "oe_GC", "oe_GG", "oe_CC"]
os.makedirs("results", exist_ok=True)

d = pd.read_csv("data/processed/rsv_features.tsv", sep="\t")
d = d[(d.year >= 2015) & (d.year <= 2026)].copy()
print(f"sequences 2015-2026: {len(d)}  (A={sum(d.Subtype=='A')}, B={sum(d.Subtype=='B')})")

# --- Step 1: within vs between lineage -------------------------------------------
dec = []
for sub, s in d.groupby("Subtype"):
    for f in FEATURES:
        r = trends.decompose_trend(s, f)
        r.insert(0, "Subtype", sub)
        dec.append(r)
dec = pd.concat(dec, ignore_index=True)
dec.to_csv("results/decomposition.tsv", sep="\t", index=False)
print("\n=== Step 1: trend decomposition (slope per year) ===")
print(dec[["Subtype", "feature", "slope_total", "p_total",
           "slope_within", "p_within", "pct_from_turnover"]].round(6).to_string(index=False))

slopes = trends.per_lineage_slopes(d, "oe_CG")
slopes.to_csv("results/per_lineage_slopes.tsv", sep="\t", index=False)
print(f"\nper-lineage CpG O/E slopes: {len(slopes)} lineages, "
      f"{slopes.sig.sum() if len(slopes) else 0} significant")
if len(slopes):
    print(slopes.groupby("Subtype").agg(n_lin=("slope", "size"),
          median_slope=("slope", "median"), n_pos=("slope", lambda x: (x > 0).sum())).to_string())

# --- Step 3: lineage growth rate vs composition ----------------------------------
growth = trends.lineage_growth_rate(d)
growth.to_csv("results/lineage_growth.tsv", sep="\t", index=False)
tab, corr = trends.growth_vs_composition(d, growth, FEATURES)
tab.to_csv("results/growth_vs_composition.tsv", sep="\t", index=False)
corr.to_csv("results/growth_composition_corr.tsv", sep="\t", index=False)
print(f"\n=== Step 3: lineage growth vs composition ({len(growth)} lineages) ===")
print(corr.round(4).to_string(index=False))

# --- Step 5: robustness -----------------------------------------------------------
rob = trends.robustness(d, "oe_CG")
rob.to_csv("results/robustness.tsv", sep="\t", index=False)
print("\n=== Step 5: robustness of the CpG O/E trend ===")
print(rob[["Subtype", "scenario", "n", "slope_total", "p_total",
           "slope_within", "p_within"]].round(6).to_string(index=False))
print("\nwrote results/*.tsv")
