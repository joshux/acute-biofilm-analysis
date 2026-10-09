# acute-biofilm-analysis

Testing Kolpen et al. 2022 (Thorax) — that bacteria in acute lung infection co-express biofilm
matrix and fast growth — on public BALF metatranscriptomes (PRJNA1056765, Tang et al. 2025 Sci Data).

Bacterial acute-infection arm only; COVID contrast arm deferred.

## Result

**The measurable part of Kolpen's claim is not supported.**

- Matrix transcripts **are** present in acute bacterial BALF — the complete alginate
  operon (16/16 genes), Psl (7/11) and Pel (3/7) — the first direct transcript
  detection of this machinery in acute infection.
- But there is **no matrix excess**: `B = log(matrix) − log(growth)` gives median
  +0.05 against the full growth panel and −0.41 against ribosomal proteins alone,
  so the sign depends on how growth is defined.
- The tempting Kolpen-compatible signal — matrix and growth rates rising together
  (r = +0.50) — is a **sequencing-depth artifact**: both arms track depth (r = 0.91
  and 0.77), and partialling depth out **inverts** the association to −0.76, i.e.
  the classical matrix/growth trade-off.

See **`RESULTS.md`** for the full account and limitations.

## Documents

| file | contents |
|---|---|
| `RESULTS.md` | **current result** — corrected, with method and limitations |
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
- **Depth is the dominant confound.** At low coverage almost every mapped read is a
  high-abundance growth transcript, which manufactures a positive matrix/growth
  correlation.

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
