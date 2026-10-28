# Working plan

Scratch roadmap. The findings write-up is [docs/PROJECT_BRIEF.md](docs/PROJECT_BRIEF.md);
the analysis itself is in `notebooks/`, which is the main pipeline.

## Goal
Adapt Ye et al. 2026 (Nat Microbiol, influenza A GC/dinucleotide method) to human RSV-A/B
from GISAID, and see whether genome composition shows any pattern worth following up.

Analysis is deliberately kept descriptive: observe patterns, do not over-model them.

## Done

- [x] Merge 8 GISAID batches (A+B) -> 33,040 genomes; QC -> 33,038
- [x] `src/rsv/composition.py` — GC content, 4 dinucleotide frequencies, 4 O/E ratios
- [x] `src/rsv/io.py` — QC filter (len 14,900–15,400, ambiguity <=1%, full date)
- [x] `src/rsv/genes.py` — MAFFT alignment to annotated reference, 10-gene extraction.
      Keeps insertions, so the ON1/BA G-gene duplications survive (validated: G recovers
      at 969 nt vs the 897 nt reference, exactly the +72 nt ON1 duplication)
- [x] `src/rsv/trends.py` — within-lineage change vs lineage replacement, robustness
- [x] Per-gene features: 33,038 genomes x 10 genes = 330,380 rows
- [x] Notebooks 01–05, documented, figures shown inline
- [x] GitHub repo, no GISAID data pushed

## Findings so far

1. RSV suppresses CpG strongly genome-wide (O/E ~0.24); the other three dinucleotides sit
   near 1, so CpG specifically is avoided.
2. Composition varies a lot between genes. **G** is the most GC-rich and most distinctive
   region — the natural one to focus on.
3. **RSV-A has removed nearly all CpG from NS1 and NS2**, its interferon-antagonist genes
   (0.16 and 0.76 CpG per gene vs 1.94 and 5.64 in RSV-B; 84.8% of RSV-A genomes have none
   at all in NS1), at the same GC content. Clearest pattern found.
4. The apparent rise in CpG O/E over time in RSV-B is **lineage replacement, not evolution**.
   Within lineages it declines, matching the paper's direction. Robust for B, not for A.

## Blocked on a decision

The paper's sustained-vs-sporadic contrast needs host and epidemiological labels. Every
GISAID RSV record is human and RSV has been sustained in humans for decades, so the dataset
holds **one transmission class** and that question cannot be answered from RSV alone.

Proposed next dataset: a matched within-genus paramyxovirus panel from NCBI, with Nipah
(138) and Hendra (92) supplying the sporadic class. See the brief for the design and the
phylogenetic confound to avoid. **Waiting on the mentor's steer before building it.**

## Possible next steps

- [ ] Confirm whether "focus on a specific region" meant a gene or a geographic region.
- [ ] If genomic: go deeper on the G gene alone (CG / GC / GG over time and by lineage).
- [ ] Per-subtype ML trees (IQ-TREE, not installed) if phylogenetic context is wanted.
- [ ] Extend the workflow to mpox once RSV is settled.

## Decisions on record
- A and B analysed separately throughout; they behave differently. Combined only for the
  A-vs-B contrast.
- Never push GISAID data. All of `data/` git-ignored except `data/reference/` (public NCBI).
  Accession IDs are also stripped from notebook outputs.
- Use `.venv`, not Anaconda base (NumPy 2.x vs 1.x-compiled extensions breaks sklearn).
