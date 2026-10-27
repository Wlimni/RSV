# Findings

Human RSV-A and RSV-B, 33,038 QC-passing GISAID genomes (2001–2026; analyses restricted to
2015–2026, n = 32,059). Method adapted from Ye et al. 2026 (*Nat Microbiol*, influenza A).

Everything below is computed **separately for RSV-A and RSV-B**, because the two subtypes
behave differently. The only combined analysis is the A-versus-B contrast itself.

---

## 1. The raw time trend is mostly lineage replacement, not evolution

The first-pass result (RSV-B CpG O/E rising over time, Spearman rho = 0.82) **does not survive
decomposition**. Splitting the trend into a within-lineage part (lineage fixed effects) and a
between-lineage part reverses the sign:

| RSV-B, CpG O/E | slope per year |
|---|---|
| Total (raw) | **+0.00065** |
| Within-lineage | **−0.00052** |
| Attributable to turnover | 180% |

This is a Simpson's paradox. **Inside** lineages CpG O/E is falling — the same direction as the
declining GC-related content Ye et al. report — while the lineages that gain share happen to
carry higher CpG O/E. Two opposing forces, not one trend.

RSV-A falls both within (−0.00015) and overall (−0.00032).

*Figure:* `fig5_within_vs_between.png`, `fig6_per_lineage_slopes.png`

## 2. RSV-B's within-lineage decline is robust; RSV-A's is not

| Scenario | RSV-A within | RSV-B within |
|---|---|---|
| Full data | −0.00015 (p<0.001) | −0.00052 (p<0.001) |
| Excluding 2020–21 | +0.000001 (p = 0.97) | −0.00054 (p<0.001) |
| Region-adjusted | −0.00013 (p<0.001) | −0.00057 (p<0.001) |
| Region-balanced | +0.00005 (p = 0.16) | −0.00032 (p<0.001) |

RSV-A's signal disappears once the COVID-19 period is removed or sampling is balanced across
regions, so it is a sampling artefact. RSV-B's holds in every scenario.

*Figure:* `fig8_robustness.png`

## 3. RSV-A has almost eliminated CpG from its interferon-antagonist genes

The clearest biological result, and it is a subtype contrast rather than a time trend.

| Gene | RSV-A CpG per gene | RSV-B CpG per gene | RSV-A % sequences with zero CpG |
|---|---|---|---|
| **NS1** | **0.16** | 1.94 | **84.8%** |
| **NS2** | **0.76** | 5.64 | **34.2%** |
| N | 6.86 | 6.09 | 0% |
| L | 41.2 | 38.2 | 0% |

NS1 and NS2 are RSV's two interferon antagonists. In RSV-A they are ~7x more CpG-depleted than
in RSV-B, and most RSV-A genomes carry **no CpG at all** in NS1. GC content is essentially
identical between the subtypes for these genes (NS1: 0.305 vs 0.303), so this is specific CpG
avoidance, not a GC-content side effect. CpG is the motif recognised by the host antiviral
protein ZAP, which makes the two genes that block interferon signalling an interesting place
to find it.

This is an observation, not a demonstrated mechanism — see the caveats.

*Figure:* `fig11_cpg_depletion.png`, `fig9_gene_composition.png`

## 4. Within RSV-B, surface and internal genes move in opposite directions

Within-lineage CpG O/E slope per year (all significant after Benjamini-Hochberg correction):

| Gaining CpG | Losing CpG |
|---|---|
| G +0.0021, F +0.0019 | P −0.0024, NS2 −0.0019, M −0.0009, L −0.0009 |

G and F are the surface glycoproteins and the targets of immunity (F is the nirsevimab and
RSVpreF target). The genome-wide decline is driven by L, which is 43% of the genome.

*Figure:* `fig10_gene_trends.png`

## 5. Lineage-level analyses are underpowered (honest negative)

Neither the lineage growth-rate correlations nor the founder-composition analysis (a tree-free
analogue of the paper's "earliest viruses" result) survives multiple-testing correction:

| Test | rho | p | p (BH) |
|---|---|---|---|
| RSV-B founder CpG O/E vs growth rate | +0.69 | 0.0065 | 0.20 |
| RSV-B GpG O/E vs growth rate | +0.58 | 0.030 | 0.30 |

Only 14 RSV-B and 32 RSV-A lineages meet the inclusion thresholds. The direction is suggestive
but the sample of lineages is too small. Reported as exploratory.

*Figure:* `fig7_growth_vs_composition.png`

---

## Caveats

- P-values with n = 32,059 are not a measure of importance; the effect sizes above are small in
  absolute terms. Judge by direction and consistency, not by p.
- "Transmission" here is a proxy (lineage turnover and growth), not measured transmission. There
  is no host contrast in RSV data, unlike the avian-to-mammal comparison in the paper.
- Sampling is Europe-heavy and the 2020–21 bottleneck distorts counts; trends are shown with and
  without region adjustment for this reason.
- The ZAP interpretation of the NS1/NS2 result is a hypothesis consistent with the data, not a
  tested mechanism. It would need expression or replication experiments.
- Lineage assignment is taken from GISAID as provided; no phylogeny was reconstructed.

## What would strengthen this next

1. Compare NS1/NS2 CpG depletion against a ZAP-binding motif model, and check whether the
   depleted positions are synonymous (selection on RNA) or non-synonymous (protein constraint).
2. Build per-subtype ML trees (IQ-TREE) and reconstruct ancestral composition, replacing the
   tree-free founder proxy used here.
3. More lineages, or a finer clade definition, to give the lineage-level tests enough power.
4. Test whether the G and F CpG increase in RSV-B coincides with the 2023+ rollout of
   nirsevimab and RSVpreF.
