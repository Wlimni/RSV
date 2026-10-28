# Genomic composition and host transmission in RSV

Project brief for discussion. William Ni Hong Hung, HKU, 28 October 2026.
Shareable version with figures: see `results/figures/` or the published page.

## The question

Ye et al. 2026 (*Nature Microbiology*, doi:10.1038/s41564-025-02257-4) measured GC content
and GC-related dinucleotide frequencies (CpG, GpC, GpG, CpC) across 115,520 influenza A
genomes. Lineages achieving **sustained** mammalian transmission showed declining GC-related
content and separated from avian-only and **sporadic** mammalian viruses. The earliest
isolates of persistent lineages already carried reduced GC, which is what makes the signal
useful for risk assessment rather than merely descriptive.

The brief was to test whether comparable features associate with host transmission in
paramyxoviruses, beginning with RSV.

## What was built

GISAID EpiRSV: 8 batches merged to 33,040 genomes (15,458 A, 17,582 B, 2001-2026); 33,038
pass QC (length 14,900-15,400 nt, <=1% ambiguous, complete date).

- Whole-genome GC content, 4 dinucleotide frequencies, 4 O/E ratios.
- Per-gene composition for all 10 genes via MAFFT alignment to annotated references
  (NC_038235.1 for A, NC_001781.1 for B); 330,380 gene records. Insertions retained, so the
  ON1/BA G-gene duplications survive (validated: G recovers at 969 nt vs the 897 nt
  reference, exactly the +72 nt ON1 duplication).
- Decomposition of temporal trends into within-lineage change versus lineage replacement,
  with robustness scenarios for uneven sequencing effort.

GISAID data is not redistributed; the repo holds code, aggregate summaries and figures only.

## Findings

### 1. The apparent temporal trend is lineage replacement, not evolution

Raw CpG O/E rises in RSV-B (Spearman rho = 0.82 on annual means). Adding lineage fixed
effects reverses the sign.

| RSV-B, CpG O/E | Slope/year | p |
|---|---|---|
| Total (raw) | +0.00065 | <0.001 |
| **Within-lineage** | **-0.00052** | <0.001 |

A Simpson's paradox: lineages with higher CpG O/E are displacing those with lower. The
within-lineage direction agrees with Ye et al.

Robustness: RSV-B's within-lineage decline holds excluding 2020-21, under region adjustment
and under region balancing. **RSV-A's does not** (p = 0.97 without pandemic years; p = 0.16
region-balanced), so RSV-A's apparent signal is a sampling artefact.

Figure: `results/figures/fig5_within_vs_between.png`, `fig8_robustness.png`

### 2. RSV-A has eliminated CpG from its interferon-antagonist genes

| Gene | CpG/gene A | CpG/gene B | % A genomes with zero CpG |
|---|---|---|---|
| **NS1** | 0.16 | 1.94 | **84.8%** |
| **NS2** | 0.76 | 5.64 | **34.2%** |
| N | 6.86 | 6.09 | 0% |
| L | 41.2 | 38.2 | 0% |

NS1 GC content is 0.305 (A) vs 0.303 (B), so this is CpG-specific avoidance rather than a GC
effect. Whether the NS1 protein sequence itself constrains CpG placement has not been tested
here; that check was dropped to keep the analysis descriptive.

Figure: `results/figures/fig11_cpg_depletion.png`, `fig9_gene_composition.png`

### 3. Lineage-level tests are underpowered

Growth-rate and founder-composition correlations do not survive Benjamini-Hochberg
correction. Strongest: RSV-B founding CpG O/E vs growth rate (rho = +0.69, p = 0.0065,
p_adj = 0.20). Only 14 B and 32 A lineages qualify. Exploratory only.

## The limitation that determines what comes next

The paper's central contrast is **sustained versus sporadic transmission**, and those labels
come from host and epidemiological history, not from sequence. Every GISAID RSV record is
human, and RSV has been sustained in humans for decades. **The dataset contains one
transmission class.**

RSV alone cannot reproduce the sustained/sporadic analysis, and no model or additional RSV
sequencing changes that. The lineage turnover measure used here is a proxy for within-host
competitive success, a different question from establishment in a new host.

The pipeline is virus-agnostic and findings 2-3 stand on their own as an immune-evasion
result. But answering the original brief requires a dataset containing both classes.

## Proposed next dataset

Complete-length genomes in NCBI Nucleotide:

| Virus | Host situation | Genomes |
|---|---|---|
| **Nipah henipavirus** | bat -> human, **sporadic** | 138 |
| **Hendra henipavirus** | bat -> horse/human, **sporadic** | 92 |
| Human RSV | sustained, human | 26,399 |
| Measles morbillivirus | sustained, human | 3,527 |
| Human metapneumovirus | sustained, human | 2,349 |
| Newcastle disease virus | sustained, avian | 1,912 |
| Human parainfluenza 3 | sustained, human | 1,547 |
| Mumps orthorubulavirus | sustained, human | 1,380 |
| Human parainfluenza 1 | sustained, human | 517 |
| Canine morbillivirus | sustained, canine | 348 |
| Avian metapneumovirus | sustained, avian | 215 |
| Bovine RSV | sustained, bovine | 71 |
| Sendai virus | sustained, murine | 36 |

### The confound to design around

Pooling these into a classifier gives near-perfect accuracy that means nothing: it learns
virus identity, not transmission class. Ye et al. avoided this because every genome compared
was influenza A - one virus, many host lineages. Across paramyxovirus genera, phylogeny and
transmission class are confounded.

The design that controls for it is **matched within-genus host contrasts**:

- *Orthopneumovirus* - human RSV vs bovine RSV
- *Metapneumovirus* - human vs avian
- *Morbillivirus* - measles (human) vs canine distemper (canine)
- *Respirovirus* - human parainfluenza 1/3 vs Sendai (murine)
- *Henipavirus* - Nipah and Hendra as the sporadic class

Sample sizes are heavily imbalanced, so compare lineage- or species-level summaries rather
than per-genome classification, and downsample the sustained groups.

## Questions for discussion

1. **Is the matched-pair paramyxovirus panel the intended direction**, or was RSV meant as
   the endpoint with the transmission framing treated more loosely?
2. **Are Nipah and Hendra an acceptable sporadic class?** They are the clearest available
   case but differ in reservoir, route and pathology, not only transmission outcome.
3. **Is the NS1/NS2 CpG result worth developing separately?** It does not address host
   transmission but is the clearest pattern we found. Developing it would mean first testing
   whether the protein sequence permits CpG at those positions, then experimental work
   (restore CpG in NS1, assay replication in ZAP-competent vs ZAP-deficient cells).
4. **Should the RSV work be written up as a methods-and-negative-result component**, given
   the headline temporal trend proved to be a turnover artefact?

## Caveats

- At n = 32,059, p-values indicate only that an effect is not noise. Composition effect sizes
  are small in absolute terms; direction and consistency across scenarios are the evidence.
- The ZAP interpretation of finding 2 is a hypothesis consistent with sequence data. No
  mechanism has been tested.
- Sampling is Europe-weighted; the 2020-21 bottleneck distorts counts. Trends are reported
  with and without region adjustment.
- Lineage assignments taken as provided by GISAID. No phylogeny reconstructed; IQ-TREE
  analysis outstanding.
