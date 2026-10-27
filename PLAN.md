# Working plan

Delete this file once the analysis stabilises; it is a scratch roadmap.

## Goal
Adapt Ye et al. 2026 (influenza A GC/dinucleotide method) to human RSV-A/B from GISAID,
using **lineage turnover over time** as the transmission proxy (no host contrast available).

## Status

- [x] Merge 8 raw batches (A+B) -> `data/full/` (33,040 genomes, QC-clean report)
- [x] Core feature module `src/rsv/composition.py` (GC content, 4 dinuc freqs, 4 O/E ratios)
- [x] I/O + QC module `src/rsv/io.py` (length 14,900–15,400, ambiguity ≤1%, full date)
- [x] `scripts/build_features.py` -> `data/processed/rsv_features.tsv` (33,038 rows)
- [x] `.gitignore` (all `data/` + PDF excluded; GISAID not redistributed)
- [ ] Notebook 01 — QC: length/ambiguity/date/geography distributions, filter justification
- [ ] Notebook 02 — composition by subtype & lineage (A vs B, dominant vs rare lineages)
- [ ] Notebook 03 — time trends 2015–2026, incl. COVID-19 disruption (2020–2021 bottleneck)
- [ ] Notebook 04 — lineage-turnover test + predictive model (composition -> lineage/year)
- [ ] R figures (ggplot2): composition distributions, time-trend lines, per-lineage panels
- [ ] Decide what aggregate results/figures are safe to commit (no per-sequence GISAID data)
- [ ] git init, remote Wlimni/RSV, first commit (code + figures only), confirm file list

## Open questions / decisions
- **Transmission proxy:** lineage turnover (expanding vs fading lineages) + drift over time.
  Confirmed with William.
- **Per-gene features:** needs alignment to a reference (RSV-A: hRSV/A ref; RSV-B ref) and
  gene coordinates. Deferred to an extension once whole-genome signal is characterised.
- **Lineage expansion metric:** define "expanding" vs "fading" by year-on-year share change
  of each lineage; threshold TBD in notebook 04.
- **COVID-19 bottleneck:** 2020–2021 shows a dip then rebound in sequence counts; check
  whether composition shifts across the bottleneck.

## Later / stretch
- Extend to mpox once the RSV workflow is solid.
- Phylogenetic context (the influenza paper used ML trees per gene for lineage assignment);
  RSV lineages come pre-assigned in GISAID, so trees are optional for now.
