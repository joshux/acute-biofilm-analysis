# Analysis plan — Kolpen's matrix+growth combination in acute bacterial infection

**Date:** 2026-10-09 · **Status:** planned (research complete, no data downloaded yet) ·
**Scope:** acute BACTERICAL infection arm only (COVID contrast arm deferred).

## 1. The question

Kolpen et al. 2022 (Thorax 77:1015, [10.1136/thoraxjnl-2021-217576](https://thorax.bmj.com/content/77/10/1015)) showed by PNA-FISH/confocal imaging that bacteria in acute lung infection form biofilm aggregates WITH high ribosome content (fast growth). The classical paradigm says matrix production and fast growth are mutually exclusive (matrix is a stationary-phase, c-di-GMP program). Kolpen never measured gene expression — the co-expression question is open and this analysis tests its cross-modality implication:

**In acute bacterial pneumonia, are matrix-biosynthesis genes and fast-growth genes co-expressed in the same samples — or does matrix expression trade off against growth (the classical view)?**

Primary comparison: bacterial-infection samples vs non-infectious controls (both BALF mNGS, same platform). This is a supporting/replication analysis of Kolpen's framework, not a test of the One-Two Punch (no COVID arm in v1).

## 2. Dataset — verified facts (subagent research, 2026-10-09)

**PRJNA1056765** — Tang et al. 2025, *Sci Data* 12:1919 ([10.1038/s41597-025-06171-6](https://www.nature.com/articles/s41597-025-06171-6)). 402-patient BALF mNGS cohort (FAHZU, Hangzhou), DNA + RNA from every patient, host-depleted reads deposited.

| group | RNA runs | paired DNA |
|---|---|---|
| **Bacterial infection** | **114** | yes, all |
| Non-infectious controls (NC) | **32** | yes (pair by Patients ID, not BioSample) |
| Fungal / TB / lung cancer (contrast, optional) | 79 / 86 / 123 | yes |

- **Group mapping:** ENA Portal API returns group in `sample_title` (868 runs; NC RNA runs are mislabeled `NCn-DNA` — filter on `library_strategy=RNA-Seq`, never on title). Cross-check: Supplementary XLSX S1/S2 at [Springer static content](https://static-content.springer.com/esm/art%3A10.1038%2Fs51597-025-06171-6/MediaObjects/41597_2025_6171_MOESM1_ESM.xlsx) (has per-run `unhost_reads`, `unhost_unrRNA_reads`, `micro_reads`).
- **⚠ The binding constraint:** only host-depleted reads were deposited. Bacterial group: median **~46k microbial reads**/sample (83% ≥10k); `unhost_unrRNA_reads` median ~443k. Depth gate required (see §6).
- **Read length discrepancy:** *Sci Data* says 50-cycle SE, companion paper says 75-cycle SE — **verify from actual FASTQ before choosing mapper settings.**
- **rRNA depletion was human-only** (Ovation Trio / AnyDeplete human): bacterial rRNA is almost certainly still in the files. Good for RP-gene mRNA (not targeted), but rRNA reads must be filtered first and the per-sample rRNA fraction reported as QC.
- Pre-computed genus relative-abundance matrices on Figshare (CC BY 4.0) let us pre-select high-biomass samples BEFORE downloading FASTQ.
- License: Figshare processed data CC BY 4.0; SRA public. Citation: Tang et al. + BioProject.

## 3. Gene panels (per-organism; full catalog in GENE-PANELS.md)

**GROWTH module** (primary axis): single-copy ribosomal protein genes (rps/rpl/rpm families) → **%RP = within-taxon fraction of RP transcripts** (Gifford 2014 — invariant to library size, depletion efficiency, and community composition). Secondary: Stäubli 127 single-copy marker-gene classifier ([bioRxiv 2025.08.26.672432](https://www.biorxiv.org/content/10.1101/2025.08.26.672432v1)). **Do NOT use rRNA:DNA ratios** (Papp 2018, Steven 2017, Blazewicz 2013 — misclassify, no universal cutoff) or PTR/iRep (need ≥5× DNA coverage per genome; 50–75bp SE won't reach it).

**MATRIX module** (species-specific, signed):
- *P. aeruginosa*: pelA–G, pslA–L, algD operon (stratify — alginate only in mucoid strains), cdrA. In-vivo anchor: Jennings 2021 (Pel 5/5, Psl 3/5 in CF sputum).
- *S. aureus/epidermidis*: icaADBC (+**icaR with inverted sign**), aap/sasG, embp, cidA (eDNA); agr sign is context-dependent — exclude or report separately.
- Enterobacteriaceae: csgBAC/csgDEFG, bcsABZC, wca (colanic), pgaABCD. **csg/bcs ABSENT from Klebsiella** — presence filter mandatory.
- *S. pneumoniae*: rlrA pilus (rrgA), lytA/cbp — **NOT cps** (capsule inversely correlated with biofilm).
- *H. influenzae*: **no EPS exists** — matrix is protein/eDNA/LOS (siaB best-validated). Flag as weakest arm.

**Presence/absence filter (the paired-DNA advantage):** score matrix genes only in samples where the taxon's panel genes are detectable in that sample's DNA — otherwise gene carriage, not expression, dominates (icaADBC carried by <50% of staph isolates; csg/bcs absent from Klebsiella).

## 4. Pipeline

1. **Pilot first** (feasibility gate, §6): 8–12 bacterial-group samples pre-selected for high biomass via the Figshare genus matrix + 4 NC samples.
2. QC; verify read length; **filter rRNA/tRNA** (bbduk + SILVA SSU/LSU NR99 + Rfam tRNA) — rRNA contamination biases TPM comparisons (Mameda & Bono 2025); report per-sample rRNA fraction and test correlation with dominant genus (depletion is taxonomically biased — SAMSA 2016).
3. **Competitive mapping** of all filtered reads to ONE concatenated index: reference genomes of panel organisms + common colonisers (Moraxella, Acinetobacter, Stenotrophomonas). BWA-MEM primary (matches source pipeline; outperforms bowtie2 per Mameda & Bono 2025), MAPQ≥20, ≥95% identity, multi-mappers via `--fraction` or msamtools. Reference strains: PAO1 (PA####), S. aureus NCTC 8325 (SAOUHSC####; note rsbU defect), S. pneumoniae TIGR4, K. pneumoniae KPPR1/MGH78578, E. coli MG1655 (b####; curli cryptic — flag), H. influenzae 86-028NP (NOT Rd KW20 — lacks hmw).
4. Per-taxon, per-gene counts (featureCounts -M --fraction); **taxon-specific normalization** (Klingenberg & Meinicke 2017 — global scaling invalid for metatranscriptomes).
5. Scores: growth = %RP within taxon; matrix = taxon-scaled panel expression, restricted to DNA-confirmed taxa, signed contributors only.
6. **Primary metric — single log-ratio balance, NOT two ratios:** `balance = log(mean matrix panel) − log(mean growth panel)` on CLR-transformed values. Two scores sharing a denominator are sum-constrained and induce a spurious NEGATIVE correlation (Filzmoser 2009; Gloor 2017) — this choice determines whether the headline is real. Sensitivity: proportionality ρp via propr (Quinn 2017), SparCC/BAnOCC.

## 5. Statistics

- Primary test: balance score vs group (bacterial vs NC), adjusted for dominant taxon; permutation null vs an expression-bin-matched random control panel (10,000×; AddModuleScore-style binning) — defends against aggregation-induced spurious correlation.
- Depth gate per Lee 2025: ~10⁴ organism-assigned reads ≈ 20% DE-gene recovery; report how many samples pass at 10⁴ / 10⁵.
- Sensitivity analyses pre-specified: ±DNA presence filter; taxon-scaled vs total-bacterial denominator; excluding Klebsiella-dominant samples (no csg/bcs); excluding agr/contested-sign genes; high-depth subsample only. BH correction; effect sizes with CIs.

## 6. Interpretation guardrails (state in any write-up)

1. **Transcript ≠ matrix.** psl is translationally repressed by RsmA without transcription change (Wei & Ma 2013); pel/psl/alg NOT differentially expressed in aggregates vs planktonic (Gannon 2024); ica transcription precedes PIA (Flückiger 2005).
2. **Bulk averages erase spatial structure** — matrix operons are surface-restricted in biofilms (Heacock-Kang 2017). A null bulk result does not falsify Kolpen; a POSITIVE bulk result is stronger evidence than Kolpen's imaging.
3. **Growth proxy is non-linear at low µ** (flat below ~0.4 h⁻¹, Matamouros 2023; single-cell ribosome content not predictive, Brettner 2024).
4. "Acute" is our framing — the dataset labels the group "Bacterial infection" (BALF within 72h of presentation; our inference, not their wording).
5. Matrix-gene in-vivo evidence (Jennings) is from CHRONIC CF sputum; applying it to acute pneumonia is extrapolation — flag as the largest evidence gap.

## 7. Verdict logic

- Matrix+ / growth+ in bacterial group (vs NC) → **supports Kolpen's combination** — first transcriptomic evidence of the novel quadrant.
- Matrix+ / growth− or trade-off (negative balance structure) → classical paradigm holds in vivo.
- Matrix undetectable at depth gate → feasibility ceiling, reported honestly (the compact-DB lesson: coverage limits are findings about the assay, not the biology).

## 8. Deliverables

This plan + GENE-PANELS.md (full catalog) in `analysis/acute-biofilm/`; scripts, per-sample scores, RESULTS.md after execution. Deterministic, seeded; all inputs from PRJNA1056765 + Figshare CC BY 4.0.
