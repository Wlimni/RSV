"""
Script Name: run_trend_analysis.py
Purpose:     Separate real within-lineage change from lineage replacement, and check the
             result holds once uneven sequencing effort is accounted for.
Inputs:      data/processed/rsv_features.tsv
Outputs:     results/decomposition.tsv, results/robustness.tsv
Notes:       RSV-A and RSV-B are analysed separately because their trends differ.
             Run: .venv/bin/python scripts/run_trend_analysis.py
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import pandas as pd
from rsv import trends

FEATURES = ["gc_content", "oe_CG", "oe_GC", "oe_GG", "oe_CC"]
os.makedirs("results", exist_ok=True)

d = pd.read_csv("data/processed/rsv_features.tsv", sep="\t")
d = d[(d.year >= 2015) & (d.year <= 2026)]
print(f"sequences 2015-2026: {len(d)}  (A={sum(d.Subtype=='A')}, B={sum(d.Subtype=='B')})")

dec = pd.concat([trends.decompose_trend(s, f).assign(Subtype=sub)
                 for sub, s in d.groupby("Subtype") for f in FEATURES],
                ignore_index=True)
dec = dec[["Subtype"] + [c for c in dec.columns if c != "Subtype"]]
dec.to_csv("results/decomposition.tsv", sep="\t", index=False)
print("\n=== Trend split into total vs within-lineage (slope per year) ===")
print("opposite signs mean the plain trend comes from lineage replacement\n")
print(dec[["Subtype", "feature", "slope_total", "p_total",
           "slope_within", "p_within"]].round(6).to_string(index=False))

rob = trends.robustness(d, "oe_CG")
rob.to_csv("results/robustness.tsv", sep="\t", index=False)
print("\n=== Does the CpG O/E result survive uneven sampling? ===")
print(rob[["Subtype", "scenario", "n", "slope_within", "p_within"]].round(6).to_string(index=False))
print("\nwrote results/decomposition.tsv, results/robustness.tsv")
