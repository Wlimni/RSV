# Working plan

Scratch roadmap. Results live in [docs/FINDINGS.md](docs/FINDINGS.md).

## Goal
Adapt Ye et al. 2026 (Nat Microbiol, influenza A GC/dinucleotide method) to human RSV-A/B
from GISAID, using lineage turnover over time as the transmission proxy (no host contrast
available — every GISAID RSV record is human).

## Done

- [x] Merge 8 raw batches (A+B) -> `data/full/` (33,040 genomes)
- [x] `src/rsv/composition.py` — GC content, 4 dinucleotide freqs, 4 O/E ratios
- [x] `src/rsv/io.py` — QC filter (len 14,900–15,400, ambiguity <=1%, full date) -> 33,038
- [x] `src/rsv/trends.py` — within/between decomposition, per-lineage slopes, growth rate,
      founder composition, robustness scenarios
- [x] `src/rsv/genes.py` — mafft alignment to annotated reference, 10-gene extraction
      (keeps insertions, so ON1/BA G-gene duplications are retained — validated: G comes
      out at 969 bp vs the 897 bp reference, exactly the +72 nt ON1 duplication)
- [x] Per-gene features for all 33,038 sequences x 10 genes (330,380 rows)
- [x] Notebooks 01–06; R figures 1–11
- [x] Step 1 — within vs between lineage decomposition **(main result: Simpson's paradox)**
- [x] Step 2 — per-gene composition **(main result: NS1/NS2 CpG depletion in RSV-A)**
- [x] Step 3 — continuous lineage growth rate vs composition (underpowered)
- [x] Step 4 — founder composition vs later success, tree-free (underpowered)
- [x] Step 5 — robustness: excl. 2020–21, region-adjusted, region-balanced
- [x] Multiple-testing correction (Benjamini-Hochberg) on all lineage-level tests
- [x] GitHub repo, no GISAID data pushed

## Next

- [ ] Are the depleted CpG positions in NS1/NS2 synonymous or non-synonymous? Synonymous
      depletion would point to selection on the RNA itself (ZAP), not the protein.
- [ ] Compare NS1/NS2 depletion against a ZAP-binding motif model (CpG spacing matters,
      not just count).
- [ ] Per-subtype ML trees (IQ-TREE, not installed) + ancestral composition reconstruction,
      to replace the tree-free founder proxy.
- [ ] More lineages / finer clades to give the lineage-level tests power.
- [ ] Does the G and F CpG increase in RSV-B line up with the 2023+ nirsevimab / RSVpreF
      rollout? Needs a date-stratified before/after test.
- [ ] Extend the workflow to mpox once RSV is settled.

## Decisions on record
- **Transmission proxy:** lineage turnover + growth rate. Confirmed with William.
- **A and B analysed separately** throughout; they behave differently (B robust, A not).
  Combined analysis only for the A-vs-B contrast.
- **Never push GISAID data.** All of `data/` is git-ignored except `data/reference/`
  (public NCBI).
- **Use `.venv`,** not Anaconda base (NumPy 2.x vs 1.x-compiled extensions breaks sklearn).
