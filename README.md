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

## Headline results

See [docs/FINDINGS.md](docs/FINDINGS.md) for the full write-up and caveats.

1. **The raw time trend is mostly lineage replacement.** Splitting it into within- and
   between-lineage parts reverses the sign for RSV-B (total +0.00065/yr, within-lineage
   −0.00052/yr) — a Simpson's paradox. Within lineages, CpG O/E declines, matching the
   direction Ye et al. report for influenza.
2. **RSV-B's within-lineage decline is robust** to excluding 2020–21 and to region
   adjustment; **RSV-A's is not** and appears to be a sampling artefact.
3. **RSV-A has almost eliminated CpG from NS1 and NS2**, its two interferon-antagonist
   genes: 0.16 and 0.76 CpG per gene versus 1.94 and 5.64 in RSV-B, with 84.8% of RSV-A
   genomes carrying no CpG at all in NS1 — at identical GC content.
4. **In RSV-B the surface glycoproteins G and F gain CpG** within lineages while internal
   genes lose it.
5. Lineage-level growth and founder-composition tests do **not** survive multiple-testing
   correction (only 14 B / 32 A lineages) and are reported as exploratory.

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
src/rsv/          core functions imported by everything
  composition.py  GC content, dinucleotide frequencies, O/E ratios
  io.py           FASTA/metadata loading and QC filtering
  trends.py       within/between-lineage decomposition, growth rates, robustness
  genes.py        align to reference and cut out the 10 RSV genes
scripts/          command-line steps (merge, features, per-gene, analyses)
notebooks/        exploration and analysis (import from src/rsv)
R/                publication-quality figures (ggplot2)
data/             GISAID data + derived tables  (git-ignored)
  Raw/            raw per-batch downloads
  full/           merged RSV_full.fasta + metadata
  processed/      rsv_features.tsv, rsv_gene_features.tsv
  reference/      public NCBI reference genomes (committed)
results/          aggregate summary tables + figures/ (committed)
docs/FINDINGS.md  results write-up with caveats
PLAN.md           roadmap / working plan
```

## Reproduce

```bash
pip install -r requirements.txt

# 1. merge raw batches -> data/full/
python3 -I scripts/merge_gisaid_batches.py --raw data/Raw --out data/full

# 2. QC + composition features -> data/processed/rsv_features.tsv
.venv/bin/python scripts/build_features.py

# 3. per-gene features (aligns to reference with mafft; ~15 min)
.venv/bin/python scripts/build_gene_features.py

# 4. analyses -> results/*.tsv
.venv/bin/python scripts/run_trend_analysis.py
.venv/bin/python scripts/run_gene_analysis.py

# 5. explore / analyse interactively
jupyter lab notebooks/

# 6. figures -> results/figures/
Rscript R/figures.R && Rscript R/figures_trends.R && Rscript R/figures_genes.R
```

Requires `mafft` on PATH for the per-gene step. The isolated `.venv` is used because
Anaconda's NumPy 2.x conflicts with its NumPy 1.x-compiled extensions.

## Method notes

- **GC content** = (G+C) / valid ACGT bases.
- **Dinucleotide frequency** = count(XY) / valid ACGT dinucleotide positions, for the four
  GC-related dinucleotides CpG, GpC, GpG, CpC.
- **O/E ratio** = freq(XY) / (freq(X)·freq(Y)); <1 indicates suppression (e.g. CpG).
- Ambiguous bases and gaps are excluded from all counts.
- Whole-genome composition first; per-gene features are a planned extension (needs
  alignment/annotation).
