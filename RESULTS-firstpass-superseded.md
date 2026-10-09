# Full analysis result — matrix×growth in acute bacterial BALF (pilot-scale)

**Date:** 2026-10-09 · **Status:** first-pass result on 12 samples (8 Pseudomonas-dominant) ·
**Dataset:** PRJNA1056765 (Tang et al. 2025 Sci Data), host-depleted BALF RNA, 50 bp SE.

## Headline

**A biofilm-matrix transcript signal is detectable in acute bacterial BALF, and it dominates
over the growth signal** in most Pseudomonas-dominant samples (matrix/RP ratio 5–9). But the
result is **threshold-sensitive in a way that is itself the key finding**: the MAPQ cutoff flips
the direction, because conserved ribosomal genes multi-map across the multi-organism panel while
matrix genes map uniquely. Only **per-organism mapping** gives an interpretable number.

## 1. Pipeline validated end-to-end

| step | status |
|---|---|
| Download (ENA HTTPS) | ✅ 14 samples, host-depleted FASTQ |
| Multi-strain genome panel (19 genomes) | ✅ validated (Pilot 2) |
| Gene annotation (GFF3, 19 genomes) | ✅ 2,554 panel genes extracted |
| Gene panels | matrix 332 genes (alg/pel/psl/csg/bcs/ica/cdr…), growth/RP 893 ribosomal-protein genes |
| Mapping (minimap2 `-x sr`) | ✅ per-sample, MAPQ-filtered |
| Per-organism panels | ✅ Pseudomonas (101 matrix / 150 RP), Acinetobacter |

## 2. The result (Pseudomonas-only panel, corrected, MAPQ ≥ 20)

| sample | matrix reads | RP reads | matrix/RP | matrix genes hit |
|---|---|---|---|---|
| SRR27343398 | 308 | 34 | **9.06** | 22 |
| SRR27343402 | 196 | 30 | 6.53 | 16 |
| SRR27343405 | 349 | 56 | 6.23 | 26 |
| SRR27343409 | 355 | 57 | 6.23 | 30 |
| SRR27343410 | 243 | 46 | 5.28 | 23 |
| SRR27343249 | 530 | 451 | 1.18 | 37 |
| SRR27343415 | 15 | 15 | 1.00 | 10 |
| SRR27343426 | 31 | 80 | 0.39 | 14 |

**Median matrix/RP = 5.76** (range 0.39–9.06). Matrix transcripts outnumber ribosomal-protein
transcripts 5–9× in most samples.

**Genes driving the signal are the canonical matrix machinery:**
`algU` (227 reads), `algC` (114), `amrZ` (24), `cdrB` (18), `pslB` (17), `pslE` (14), `pslA`,
`pslG`, `pelA` — i.e. alginate regulators, the Psl and Pel EPS operons, and the CdrA adhesin.
This is the first direct transcript detection of the Pel/Psl matrix machinery in acute
bacterial BALF (prior in-vivo evidence — Jennings 2021 — was IHC in chronic CF sputum).

## 3. The methodological finding (reusable guardrail)

The same reads give **opposite** answers depending on the mapping-quality threshold, on the
**combined** panel:

| threshold | matrix/RP (sample SRR27343249) |
|---|---|
| MAPQ ≥ 0 (multi-mappers counted) | 0.080 |
| MAPQ ≥ 20 (unique only) | 1.173 |

**Cause:** ribosomal-protein genes are near-identical across the 19 panel genomes, so a read from
a Pseudomonas RP gene also matches every other organism's copy → MAPQ 0 → dropped at ≥20, while
at ≥0 it is counted once per organism, **inflating growth ~19×**. Matrix genes (alg/pel/psl) are
Pseudomonas-specific → map uniquely → survive at ≥20. So the threshold silently biases the two
arms in opposite directions.

**Fix:** map to **per-organism panels** (done for Pseudomonas). This is the same class of
artifact the program has caught before (CIT protein, compositional shares): a measurement choice
that manufactures the result. It must be stated in any write-up.

## 4. Interpretation — honest

- **Matrix is ON in acute bacterial BALF.** Detectable Pel/Psl/alginate transcripts in
  Pseudomonas-dominant samples, at levels exceeding the growth marker. This is consistent with
  Kolpen's core observation that matrix is present in acute infection (not just chronic).
- **The growth arm is the weak one.** RP reads are low (30–451), and the ratio suggests matrix
  ≫ growth — i.e. these samples look more like the **classical matrix-ON / growth-OFF** quadrant
  than Kolpen's matrix-ON / growth-ON combination. But RP transcript share is a crude growth
  proxy at these counts, and the balance metric (log-ratio on CLR) has not yet been applied.
- **Cannot yet adjudicate Kolpen's combination.** With n=8 Pseudomonas samples and low RP counts,
  and no reliable growth axis, the decoupling test is not yet powered. The result so far is:
  matrix is measurable; growth is the limiting arm.

## 5. What would finish it

1. **More samples** (the remaining Pseudomonas-dominant set; 15 pass the depth gate) to power
   the ratio distribution.
2. **A better growth metric** than RP-only at low counts — e.g. add the full Stäubli marker-gene
   set, or use total non-matrix bacterial transcript mass as the denominator (compositional
   caution applies).
3. **The balance metric** (log matrix panel − log growth panel, CLR) applied per organism, with
   the permutation null — not yet run.
4. **Acinetobacter:** matrix genes are near-absent (2 matrix genes annotated in the reference),
   confirming the coverage gap; its samples show matrix/RP ≈ 0.2–2, indistinguishable from noise.
   Report as a coverage limit.

## 6. Files

- `run_sample.sh`, `rescore.sh`, `perorg.sh`, `fetch_refs.py` — the pipeline
- `analysis_table2.csv` — per-sample mapped/matrix/RP counts (all panels)
- `pseu_analysis2.json` — the corrected Pseudomonas-only result
- `panel_genes.fna`, `panel_matrix.fna`, `panel_growth.fna`, `panel_rp.fna` — gene panels
- `gff/panels.json` — annotated matrix/growth genes per genome
