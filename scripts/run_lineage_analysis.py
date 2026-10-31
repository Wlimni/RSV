"""
Script Name: run_lineage_analysis.py
Purpose:     Fit composition trends over time inside each RSV lineage, following the design
             of Ye et al. 2026 (their Figs 2 and 3), and write the summary tables the R
             figure script reads.
Inputs:      data/processed/rsv_gene_features.tsv, data/processed/rsv_features.tsv
Outputs:     results/lineages.tsv          the lineages analysed
             results/slopes_genome.tsv     slope per lineage per feature, genome level
             results/slopes_by_gene.tsv    slope per lineage per gene
             results/gc_correlation.tsv    dinucleotide frequency vs GC content
             results/slopes_genome_region_adjusted.tsv  the same fits holding region constant
             results/sensitivity_level4.tsv  same fits at a finer lineage level
Notes:       "Genome level" means the ten coding regions concatenated, which is how the
             paper defines genomic GC, GC3 and GC12.
             The window runs 2001 to end of June 2026. The later months of 2026 are cut
             because deposition lag makes them a non-random subset, and each lineage is
             trimmed to the years it was actually sampled in.
             Run: .venv/bin/python scripts/run_lineage_analysis.py
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import pandas as pd
from rsv import lineage

GENOME_FEATURES = ["gc_content", "gc3", "gc12",
                   "per1k_CG", "per1k_GC", "per1k_GG", "per1k_CC", "per1k_TA",
                   "oe_CG", "oe_GC", "oe_GG", "oe_CC", "oe_TA"]
DINUC_FREQ = ["per1k_CG", "per1k_GC", "per1k_GG", "per1k_CC", "per1k_TA"]
os.makedirs("results", exist_ok=True)

g = pd.read_csv("data/processed/rsv_gene_features.tsv", sep="\t")
dates = pd.read_csv("data/processed/rsv_features.tsv", sep="\t",
                    usecols=["Accession ID", "date"])
g = g.merge(dates, on="Accession ID", how="left")
# prepare() applies the date window and per-lineage trimming
print(f"{len(g)} rows covering {g['Accession ID'].nunique()} genomes")

# ---- genome level (the concatenated coding regions) --------------------------------
genome = g[g.gene == "GENOME"]
gdf, keep = lineage.prepare(genome)
keep.to_csv("results/lineages.tsv", sep="\t", index=False)
print(f"\n=== Lineages analysed (level 3, >=5 years, >=200 sequences) ===")
print(keep.round(2).to_string(index=False))

slopes = lineage.lineage_slopes(gdf, GENOME_FEATURES)
slopes.to_csv("results/slopes_genome.tsv", sep="\t", index=False)

# Same fits with region held constant. Sequencing geography shifted strongly over the
# period, including inside single lineages, so a slope that survives this is the safer one.
slopes_adj = lineage.lineage_slopes(gdf, GENOME_FEATURES, covariates=("region",))
slopes_adj.to_csv("results/slopes_genome_region_adjusted.tsv", sep="\t", index=False)

# Slim per-sequence table for the scatter plots. Kept in data/ (not committed) because it
# is per-sequence GISAID-derived; the truncation logic lives only here so R cannot drift.
gdf[["Subtype", "lineage_grp", "dec_date"] + GENOME_FEATURES].to_csv(
    "data/processed/plotdata_genome.tsv", sep="\t", index=False)

print("\n=== GC content, GC3 and GC12: change per year ===")
piv = (slopes[slopes.feature.isin(["gc_content", "gc3", "gc12"])]
       .pivot_table(index=["Subtype", "lineage_grp"], columns="feature", values="slope"))
piv["gc3_vs_gc12"] = (piv["gc3"] / piv["gc12"]).round(2)
print(piv.round(6).to_string())
ratio = piv["gc3"].sum() / piv["gc12"].sum()
print(f"\noverall GC3 changes {ratio:.2f}x as fast as GC12 "
      f"(paper reports 5.65x for influenza)")

print("\n=== CpG O/E slope per lineage, before and after region adjustment ===")
cg = slopes[slopes.feature == "oe_CG"][["Subtype", "lineage_grp", "n", "slope", "p"]]
cga = slopes_adj[slopes_adj.feature == "oe_CG"][["Subtype", "lineage_grp", "slope", "p"]]
cg = cg.merge(cga, on=["Subtype", "lineage_grp"], suffixes=("", "_adj"))
cg["same_sign"] = (cg.slope * cg.slope_adj) > 0
print(cg.round(6).to_string(index=False))
print(f"\nslopes keeping their direction after region adjustment: "
      f"{cg.same_sign.sum()}/{len(cg)}")

# ---- per gene ----------------------------------------------------------------------
per_gene = g[g.gene != "GENOME"]
pg, _ = lineage.prepare(per_gene)
by_gene = lineage.lineage_slopes(pg, ["gc_content", "oe_CG"],
                                 by=("Subtype", "lineage_grp", "gene"))
by_gene.to_csv("results/slopes_by_gene.tsv", sep="\t", index=False)
print(f"\nper-gene slopes: {len(by_gene)} lineage-gene-feature fits")

# ---- is the dinucleotide signal just GC moving? -------------------------------------
corr = lineage.gc_vs_dinucleotide_correlation(gdf, DINUC_FREQ)
corr.to_csv("results/gc_correlation.tsv", sep="\t", index=False)
print("\n=== Correlation of dinucleotide frequency with GC content (mean over lineages) ===")
print(corr.groupby("feature").r.agg(["mean", "min", "max"]).round(3).to_string())

# ---- sensitivity: does a finer lineage level change the picture? --------------------
gdf4, keep4 = lineage.prepare(genome, level=4)
s4 = lineage.lineage_slopes(gdf4, ["gc_content", "gc3", "gc12", "oe_CG"])
s4.to_csv("results/sensitivity_level4.tsv", sep="\t", index=False)
print(f"\n=== Sensitivity check at level 4 ({len(keep4)} lineages) ===")
cmp = pd.concat([
    slopes[slopes.feature.isin(["gc_content", "gc3", "oe_CG"])].assign(level=3),
    s4[s4.feature.isin(["gc_content", "gc3", "oe_CG"])].assign(level=4)])
print(cmp.groupby(["level", "feature"]).slope.agg(
    median="median", n_lineages="size",
    n_negative=lambda s: (s < 0).sum()).round(6).to_string())
print("\nwrote results/*.tsv")
