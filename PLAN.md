# Analysis plan v2 — matrix×growth co-expression in acute bacterial lung infection

**Date:** 2026-10-09 (v2, post-pilot) · **Repo:** joshux/acute-biofilm-analysis · **Scope:** acute bacterial infection arm only.

> v1 plan assumed a bacterial-vs-NC contrast and a classic-CAP pathogen panel. The pilot
> (see `RESULTS-PILOT.md`, `RESULTS-PILOT2.md`) showed both were wrong. This version reflects
> what the data actually are.

## 1. Question

Kolpen et al. 2022 (Thorax 77:1015) showed by imaging that bacteria in acute lung infection form
matrix-containing aggregates WITH high ribosome content (fast growth) — the classical paradigm
says matrix production and fast growth trade off. They never measured gene expression. **Test:
in acute bacterial BALF, do matrix-biosynthesis genes and fast-growth genes co-express (Kolpen)
or trade off (classical)?**

## 2. Dataset — PRJNA1056765 (Tang et al. 2025, Sci Data 12:1919)

| group | RNA runs | notes |
|---|---|---|
| Bacterial infection | **114** | 50 bp SE, host-depleted; median 46,202 microbial reads, 10 samples ≥1M |
| Fungal / TB | 79 / 86 | secondary comparator (have bacterial signal) |
| NC (non-infectious) | 32 | **NO bacterial signal — not a usable control (pilot)** |

All runs have paired DNA (presence/absence filter). Downloaded from ENA over HTTPS; works.
Read length 50 bp SE (verified).

**Actual community (verified, 114 bacterial samples):** dominated by non-fermenting hospital
opportunists — **Acinetobacter** (top genus 38/114, mean relab 0.151), **Pseudomonas** (32/114,
0.167), Haemophilus (10/114), plus a **reagent-contaminant signature** (Comamonas, Delftia,
Ralstonia, Sphingomonas, Rhodococcus; median 0.179 of relab, >0.30 in 31/114). NOT classic CAP
pathogens.

## 3. Reference panel (v2 — multi-strain, validated)

19 genomes / 158 contigs. Multi-strain/multi-species per genus: **Acinetobacter** (A. baumannii
17978 + AB5075, A. nosocomialis, A. pittii, A. johnsonii, A. lwoffii), **Klebsiella** (MGH78578,
NTUH-K2044, HS11286, K. variicola), **Pseudomonas** (PAO1, PA14, P. putida), Stenotrophomonas
(K279a), plus S. aureus NCTC 8325, S. pneumoniae TIGR4, E. coli MG1655, H. influenzae 86-028NP.

**Validated:** multi-strain panel raised MAPQ≥20 mapping 6.3× (0.108% → 0.681%) and made genus
assignment agree with the composition table. Single-strain-per-genus was the binding constraint.

**Still to add:** contaminant genera (Comamonas, Delftia, Ralstonia, Sphingomonas) — for
*exclusion*, so kitome reads are not attributed to infection.

## 4. Gene modules

- **GROWTH:** single-copy ribosomal protein genes (rps/rpl/rpm) → **%RP within taxon** (Gifford
  2014; invariant to library size/depletion/community composition). Secondary: Stäubli 127-MG
  classifier. **Not** rRNA:DNA (Papp 2018, Steven 2017, Blazewicz 2013) or PTR/iRep (needs ≥5×
  DNA coverage).
- **MATRIX (signed, per-organism):** *Pseudomonas* pelA-G / pslA-L / algD operon / cdrA;
  Enterobacteriaceae csg / bcs / wca / pga; *Staph* icaADBC (icaR inverted) / aap / sasG / cidA;
  *S. pneumoniae* rlrA pilus (NOT cps); *H. influenzae* — no EPS exists (flag).
- **⚠ MATRIX-MODULE GAP:** the dominant organism, **Acinetobacter, has no characterized matrix
  operon** (only bap, PNAG-like, capsule). The module cannot be scored where the biomass is
  highest. This is the central design constraint of v2.

## 5. Design (v2)

**Primary:** within-bacterial, **restricted to samples whose dominant organism has a defined
matrix operon** (Pseudomonas- and Klebsiella/Enterobacteriaceae-dominant; ~32 + ~30 samples).
Compare matrix×growth co-expression across these; report Acinetobacter-dominant samples (38) as
a **coverage limit**, not a biological null.

**Secondary:** bacterial vs fungal/TB infection (both arms carry bacterial signal).

**Contaminant handling:** flag the 31/114 high-contaminant samples; add contaminant genomes to
the index and exclude their reads; report per-sample contaminant fraction as a covariate.

## 6. Pipeline

1. Download host-depleted FASTQ (ENA HTTPS). Filter rRNA/tRNA (bbduk + SILVA / SortMeRNA) —
   bacterial rRNA was NOT depleted (human-only depletion in the source study); report per-sample
   rRNA fraction.
2. Competitive mapping (minimap2 `-x sr` / BWA-MEM) to the multi-strain panel + contaminant
   genomes; MAPQ≥20; genus-level assignment; report unassigned fraction as QC.
3. Paired-DNA presence/absence filter per taxon (gene carriage vs expression).
4. Taxon-specific normalization (Klingenberg & Meinicke 2017); %RP for growth.
5. **Primary metric: single log-ratio balance** `log(matrix panel) − log(growth panel)` on
   CLR values — NOT two ratios sharing a denominator (Filzmoser 2009; Gloor 2017 induce spurious
   negative correlation). Sensitivity: proportionality (propr), SparCC, BAnOCC.
6. Permutation null vs expression-bin-matched random control panel; adjust for dominant taxon;
   pre-specified sensitivity analyses; BH correction; effect sizes with CIs.

## 7. Guardrails

Transcript ≠ matrix (psl post-transcriptionally repressed; Gannon 2024 found no pel/psl/alg DE in
aggregates). Bulk averaging erases spatial structure (a null does not falsify Kolpen; a positive
is stronger than their imaging). Growth proxy non-linear at low µ. "Acute" is our framing of
their "Bacterial infection" label. Matrix-gene in-vivo evidence is from chronic CF, extrapolated.

## 8. Verdict logic

- Matrix+ / growth+ in matrix-measurable samples → supports Kolpen's combination.
- Trade-off (negative balance) → classical paradigm holds in vivo.
- Matrix unmeasurable (Acinetobacter-dominant, or below depth gate) → reported as assay coverage
  limit, not biology.
