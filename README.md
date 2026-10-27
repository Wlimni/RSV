# RSV genomic composition and transmission dynamics

Analysis of nucleotide-composition features (GC content and GC-related dinucleotide
frequencies) across human RSV-A and RSV-B genomes from GISAID, and how those features
relate to lineage turnover and time.

This adapts the method of Ye et al. 2026 (*Nature Microbiology*, "Genomic features
associated with sustained mammalian transmission of avian influenza A viruses",
doi:10.1038/s41564-025-02257-4) from influenza A to RSV.

## Biological question

The influenza paper separated persistent mammalian lineages from bird-only / sporadic
lineages using declining GC content and GC-dinucleotide (CpG, GpC, GpG, CpC) frequencies.
RSV has no host contrast here — every GISAID record is human — so we reframe the question
around **lineage turnover as the transmission proxy**:

- Does genome composition drift systematically over time (2015–2026)?
- Do lineages that expand and dominate (e.g. A.D.1, B.D.E) differ in composition from
  those that fade out?
- Can composition predict a sequence's lineage and collection year?

A positive signal would make composition a candidate feature for monitoring which RSV
lineages are gaining ground, and a template to later extend to other viruses (e.g. mpox).

## Data

- **Source:** GISAID EpiRSV, downloaded in 8 date-range batches (RSV-A and RSV-B),
  complete / high-coverage only, low-coverage excluded.
- **Merged set:** 33,040 genomes (15,458 A / 17,582 B), 2001–2026, mostly 2018+.
- **QC-passing:** 33,038 (length 14,900–15,400 bp, ambiguous bases ≤1%, full collection date).

**GISAID data is not redistributed in this repository.** All of `data/` is git-ignored.
To reproduce, download the batches into `data/Raw/` and regenerate (see below). A GISAID
accession list for acknowledgement is produced by the pipeline and kept local.

## Layout

```
src/rsv/          core functions (composition metrics, I/O + QC) imported by everything
scripts/          command-line steps (merge raw batches, build feature table)
notebooks/        exploration and analysis (import from src/rsv)
R/                publication-quality figures (ggplot2), fed from the feature table
data/             GISAID data + derived tables  (git-ignored)
  Raw/            raw per-batch downloads
  full/           merged RSV_full.fasta + metadata
  processed/      rsv_features.tsv (analysis-ready)
results/figures/  generated figures (committed)
docs/             notes
PLAN.md           roadmap / working plan
```

## Reproduce

```bash
pip install -r requirements.txt

# 1. merge raw batches -> data/full/
python3 -I scripts/merge_gisaid_batches.py --raw data/Raw --out data/full

# 2. QC + composition features -> data/processed/rsv_features.tsv
python3 -I scripts/build_features.py

# 3. explore / analyse
jupyter lab notebooks/

# 4. figures
Rscript R/figures.R
```

## Method notes

- **GC content** = (G+C) / valid ACGT bases.
- **Dinucleotide frequency** = count(XY) / valid ACGT dinucleotide positions, for the four
  GC-related dinucleotides CpG, GpC, GpG, CpC.
- **O/E ratio** = freq(XY) / (freq(X)·freq(Y)); <1 indicates suppression (e.g. CpG).
- Ambiguous bases and gaps are excluded from all counts.
- Whole-genome composition first; per-gene features are a planned extension (needs
  alignment/annotation).
