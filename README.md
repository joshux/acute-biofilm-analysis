# acute-biofilm-analysis

Testing Kolpen et al. 2022 (Thorax) — that bacteria in acute lung infection carry biofilm
matrix **and** grow fast, while chronic infection carries matrix but grows slowly — on public
metatranscriptomes.

## Result — two-dimension (quadrant) test

Kolpen's two axes are **decoupled** (their words: *"the difference lies primarily in metabolic
rates, not bacterial architecture"*), so the test is per-group quadrant placement, not a
correlation. Groups: acute BALF (PRJNA1056765, Tang 2025), chronic CF sputum (PRJEB24688,
Rossi 2018), healthy non-COPD BALF (PRJNA390194, Ren 2018).

- **Growth axis (%RP):** acute is significantly **faster** than chronic in vivo —
  **10.7 % vs 4.6 %, p = 0.0015**. The metric is validated on the dataset's own lab anchors
  (exponential 13.6 % FAST, stationary 2.1 % SLOW). Chronic in vivo sits at the stationary state;
  acute sits between stationary and exponential, closer to exponential.
- **Matrix axis:** matrix machinery is transcribed in **both** acute and chronic
  (acute 92 %, chronic 100 % present) — Kolpen's architecture claim holds and the classical
  acute-planktonic model is contradicted. Acute matrix-ON is a **lower bound** (the acute arm is
  ~300× shallower; at matched depth it reaches 100 %).
- **Healthy controls:** zero bacterial signal (kitome floor) — nothing to place.

Net: acute occupies **matrix-ON × FAST** and chronic **matrix-ON × SLOW** — Kolpen's prediction,
on both axes.

> **Correction (2026-10-10).** The first pass reported "acute just as slow as chronic." That was a
> metric artifact: on a *multi-genome* panel, conserved ribosomal-protein reads map equally to
> PAO1/PA14/putida, get MAPQ 0, and are dropped at MAPQ≥20 — deflating %RP and compressing all
> groups. Measuring growth against a single reference genome (with other genera as decoys) fixes it.
> See `RESULTS-QUADRANT-CORRECTED.md` → "What went wrong".

See **`RESULTS-QUADRANT-CORRECTED.md`** for the full account. The coupling/correlation test
(`RESULTS.md`) answered a question Kolpen never asked and is superseded.

## Documents

| file | contents |
|---|---|
| `RESULTS-QUADRANT-CORRECTED.md` | **current result** — corrected growth metric, quadrant test, three groups |
| `RESULTS-QUADRANT.md` | first pass — growth metric artifact; superseded (kept for the audit trail) |
| `RESULTS.md` | superseded coupling/correlation test (the depth-artefact caution still stands) |
| `RESULTS-firstpass-superseded.md` | earlier first-pass result; its headline ratio is superseded (see RESULTS.md "Correction") |
| `RESULTS-PILOT.md` | pilot 1: dataset/depth verified; three assumptions broken (community, panel, NC arm) |
| `RESULTS-PILOT2.md` | pilot 2: multi-strain panel rebuilt and validated (6.3x mapping gain); matrix-module gap for Acinetobacter |
| `PLAN.md` | analysis plan v2 (post-pilot) |

## Method notes worth keeping

- **Decoy competition is required.** Conserved genes (`fusA`, `rpoB`, `tuf`) map
  near-identically across genera. Without competitor genomes in the index they are
  assigned high MAPQ and counted as target-organism transcripts even when they came
  from a co-infecting species. `comp_<org>.mmi` = target panel + one reference genome
  for each of 34 co-infecting genera.
- **Measure both arms in one run.** Mapping matrix and growth in separate runs
  conflates panel composition with mapping behaviour, and produced a wrong answer
  earlier in this project.
- **Select and measure from the same source.** The deposited composition table
  disagrees with the RNA reads on the dominant genus for 18/36 samples, always in the
  same direction (table says Pseudomonas, reads say Acinetobacter). The cohort is
  therefore defined by the reads.
- **Depth is the dominant confound — and it cuts differently per axis.** A ratio (growth, %RP)
  is depth-robust; a detection call (matrix systems) is not. The acute arm's shallow depth
  (median 312 target reads vs 15,171 chronic) makes its matrix-ON rate a lower bound.
- **A public dataset's own internal anchors beat cross-dataset thresholds.** PRJEB24688 ships
  in vitro exponential and stationary cultures alongside the in vivo sputum, so the FAST/SLOW
  call is made against cultures measured in the same experiment — and those samples are labelled
  `isolation_source = laboratory culture`, which must be separated from the sputum runs.

## Data

| file | contents |
|---|---|
| `results_final.csv` | per-sample table (composition, rates, balances, CIs) |
| `analysis_final.json` | full per-sample records |
| `fig_balance.svg` | balance distribution, depth confound, partial correlation |
| `rna_genus_mapped.json` | RNA-derived genus composition per sample |
| `matrix_gene_census.json` | per-gene matrix detection counts |
| `qp_panels.json` | panel composition |
| `decoy/manifest.json` | the 34 decoy genomes |
| `refs2/manifest.json`, `ref2genus.json` | multi-strain reference panel |
| `multi_panel_mapping.json` | per-sample mapping counts by genus |

## Pipeline

`build_quant_panel.py` → `fetch_decoys.py`/`retry_decoys.py` → build `comp_<org>.mmi`
→ `rna_composition.py` → `quantify3.py` (per sample) → `analyze_final.py` →
`final_stats.py` → `make_figure.py`. Reproduce block is at the end of `RESULTS.md`.
