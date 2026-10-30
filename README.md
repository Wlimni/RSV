# RSV genome composition over time

Does genome composition change over time within human RSV lineages, and does the pattern
differ between RSV-A and RSV-B?

The method is adapted from Ye et al. 2026 (*Nature Microbiology*, "Genomic features
associated with sustained mammalian transmission of avian influenza A viruses",
doi:10.1038/s41564-025-02257-4), which found GC content and GC-containing dinucleotide
frequencies declining within influenza lineages that had established sustained transmission
in mammals.

## Headline results

**The influenza pattern does not transfer to RSV.** Stated plainly, because the temporal
result is mostly negative:

1. **No consistent CpG trend.** Lineages disagree in direction — A.D.1 and B.D.E rise,
   A.D.3 and A.D.5 fall. All twelve influenza lineages in the source paper moved the same
   way; these do not.
2. **The CpG result fails its own sensitivity check.** Defining lineages one level finer
   flips the median direction, so lineage granularity drives it rather than biology.
3. **GC3 and GC12 move in opposite directions** — GC3 rises in 6 of 8 lineages, GC12 falls
   in 6 of 8. Influenza had both falling, GC3 about 5.65x faster. Synonymous and
   non-synonymous sites are under different pressure in RSV, but this is not that signature.
4. **Geography does not explain it.** Slopes keep their direction after region adjustment.
5. **All effects are small**, around 1e-4 per year.

This is a reasonable biological outcome rather than a failed analysis. RSV has circulated in
humans for a very long time and is already adapted; the influenza result came from lineages
that had recently entered a *new* host and were still adjusting. No such transition exists
inside human RSV.

**The clearest signal in the project is not temporal.** In RSV-A the two interferon
antagonist genes, NS1 and NS2, carry almost no CpG at all (0.16 and 0.76 per gene, against
1.94 and 5.64 in RSV-B), at effectively identical GC content, and about 85% of RSV-A genomes
contain no CpG anywhere in NS1. Genes N and L behave the same in both subtypes, so this is
not a general A-versus-B difference. See notebook 04. This is an observation about sequence,
not a demonstrated mechanism.

## Start here

The **notebooks are the pipeline** and are written to be read in order. Each section ends
with a figure.

| Notebook | Covers |
|---|---|
| `01_data.ipynb` | the dataset, quality control, and four known problems with it |
| `02_features.ipynb` | what GC content, dinucleotide O/E, GC3 and GC12 actually measure |
| `03_lineage_trends.ipynb` | **the core analysis** — composition over time within lineages |
| `04_genes.ipynb` | composition gene by gene, and the NS1/NS2 contrast |

## Data

GISAID EpiRSV, downloaded in eight date-range batches. 33,040 genomes merged, 33,038 pass
quality control (length 14,900-15,400 nt, at most 1% ambiguous bases, complete collection
date).

The analysis window is **2015-01-01 to 2025-12-31**, giving 29,050 genomes in 8 lineages.
Two bounds, both set after auditing the data:

- **2026 is excluded** because deposition lag leaves its later months nearly empty
  (18 genomes in September 2026, against an all-year September average near 1,400).
- **Pre-2015 is excluded** because it is about 3% of the data and RSV-A has essentially no
  sequences before 2010.

Sequencing geography shifted heavily over the period, including inside single lineages, so
**region is included as a covariate** and results are reported with and without it.

**GISAID data is not redistributed here.** All of `data/` is git-ignored except the public
NCBI reference genomes. Accession IDs are stripped from notebook outputs. To reproduce,
download the batches into `data/Raw/` and run the pipeline below.

## Layout

```
notebooks/        the documented analysis - read these first
src/rsv/          core functions the notebooks import
  composition.py  GC content, dinucleotide frequencies and O/E, GC3/GC12, UpA
  io.py           loading and quality control
  genes.py        align to reference and cut out the 10 RSV genes
  lineage.py      lineage truncation, decimal dates, within-lineage regression
scripts/          command-line stages
R/figures.R       all five figures
results/          summary tables and figures (committed; no per-sequence data)
data/             GISAID data and derived tables (git-ignored)
  reference/      public NCBI reference genomes (committed)
```

## Reproduce

```bash
pip install -r requirements.txt           # needs mafft on PATH

python scripts/merge_gisaid_batches.py --raw data/Raw --out data/full
python scripts/build_features.py          # QC + whole-genome composition
python scripts/build_gene_features.py     # per-gene via alignment, ~20 min
python scripts/run_gene_analysis.py       # per-gene summary
python scripts/run_lineage_analysis.py    # the core within-lineage trends

jupyter lab notebooks/                    # notebooks call R themselves
```

The project uses an isolated `.venv`, because the Anaconda installation on the original
machine has a NumPy version conflict that breaks some scientific packages.

## Method notes

- **GC content** = (G+C) / valid ACGT bases.
- **Dinucleotide O/E** = f(XY) / (f(X) x f(Y)). Below 1 means the pair is suppressed.
  Reported for CpG, GpC, GpG, CpC, and UpA.
- **UpA is exploratory.** It is not used in the source paper; it is included because the
  RNA-virus literature treats it as CpG's companion motif (Greenbaum et al. 2008,
  *PLoS Pathogens* 4:e1000079).
- **GC3 / GC12** = GC at the third codon position, versus the first and second. Third
  positions are mostly synonymous, so comparing the two separates pressure on the RNA from
  pressure on the protein.
- **Genome level** means the ten coding regions concatenated, which is how the paper defines
  genomic GC, GC3 and GC12 and keeps the reading frame defined.
- Gene extraction aligns each genome to its subtype's annotated reference with MAFFT,
  keeping insertions, so the ON1 and BA G-gene duplications survive. Verified: the G gene
  recovers at 969 nt against an 897 nt reference, exactly the known 72 nt duplication.
- Trends are fitted **separately inside each lineage**, against decimal collection date.
  Fitting across pooled data would let one lineage replacing another masquerade as change.
- Lineage labels are truncated to three levels (`A.D.1`, `B.D.E`), requiring at least 5
  years and 200 genomes. Full labels are too fine to support a regression; level 4 is run
  as a sensitivity check.
