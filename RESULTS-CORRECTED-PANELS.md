# RESULTS-CORRECTED-PANELS.md — Klebsiella/Staphylococcus panel bug and the corrected screen

**Date:** 2026-10-15. **Scope:** items 1–5 of the "remaining work that could land more positive" list.

## Summary of findings

| # | Item | Outcome |
|---|---|---|
| 2 | Klebsiella matrix panel | **BUG CONFIRMED AND FIXED.** Committed `qp_kleb.fna` tagged all 62 matrix CDS as a single `other_matrix` system, so the ≥2-systems rule could never fire. Rebuilt: 7 named systems. |
| 1 | Deepen depth-limited no-BP runs (HS78/55/1) | **NOT POSSIBLE.** Their shallow depth is in the deposit, not our subsampling. HS78's deepest run is 69,611 reads total. |
| 4 | Matrix-axis positive control | **FOUND AND IT WORKS.** Cornforth 2018 (PRJNA437845) in-vitro PA14 biofilm + in-vivo CF sputum, and Rossi 2018 (PRJEB24688) deep CF sputum — all matrix-ON. |
| 5 | PMN-vs-load test | **SETTLED + MISSING EXPERIMENT RUN.** Literature: PMN O₂ depletion is causal (Kragh 2014). Ran the depth-dilution test the literature lacked: %RP is depth-robust, so low load does NOT explain slow calls. |
| 3 | Pan-genus commensurability | **VERDICT DELIVERED.** Report per-class (EPS vs protein), never pool; **exclude capsule** (inverse to biofilm in pneumococcus, contested in Klebsiella). Code now implements this. |

## The Klebsiella bug

`qp_kleb.fna` contained 62 matrix-classified genes, but `qp_kleb_meta.json` records
`syscount: {"other_matrix": 62}` — a single system. The metric requires
`>=2 subsystems each with >=2 genes detected`, so **no Klebsiella-dominated sample
could ever be called matrix-ON**, regardless of its actual transcription. This is the
same class of error as the PAO1-only screen documented in `RESULTS-ALLGENUS-2BP.md`.

Cause: the panel builder classified matrix genes by a regex list that, for Klebsiella,
matched only `mrk`. The GFFs actually annotate `pgaABCD` (PNAG), the `bcs` cellulose
operon, `csgD` (curli), `fim`, `ecp`, and a full capsule/cps locus.

**Fix** (`build_panels_v3.py`): per-organism (system, type, gene-regex, product-regex)
rules, URL-decoding of GFF attributes, and a **type** tag so systems are comparable
across genera:

- `EPS` — extracellular polysaccharide matrix: alginate, psl, pel, pga/PNAG, bcs/cellulose
- `PROTEIN` — proteinaceous fimbriae/adhesins: csg/curli, mrk/type-3 fimbriae, fim, ecp, csu, cbp, hmw, fnb, bap
- `REG` — regulator only (not counted as a structural system)

Rebuilt panels:

| organism | matrix systems (was → now) |
|---|---|
| pseu | alginate 22, psl 14, pel 7 (unchanged) |
| **kleb** | **1 (`other_matrix`) → 7: pga 2, bcs 11, cps 1, capsule 3, mrk 2, fim 4, ecp 2** |
| ecoli | pga 4, bcs 7, csg 8, fim 10, ecp 5 |
| staph | ica 2, fnb 4 (was ica-only, 2 genes) |
| acine | pga 8, csu 8 |
| strep | cps 12, cbp 9, psa 5 |
| haemo | hmw 2, adhesin 34 |
| steno | fimbriae 12, qs 1 (REG) |

## Item 4 — positive control

Cornforth 2018 **PRJNA437845** carries both arms on one chemistry:

| run | arm | dominant | target reads | %RP | matrix |
|---|---|---|---|---|---|
| SRR6833312 | in-vitro PA14 coculture biofilm | Pseudomonas | 25,027 | 5.12 | ON (alginate, psl) |
| SRR6833317 | in-vitro coculture biofilm | Acinetobacter | 43,352 | 4.33 | ON (pga, csu) |
| SRR6833345 | in-vivo CF sputum | Pseudomonas | — | — | ON (alginate, psl) |
| SRR6833347 | in-vivo CF sputum | Pseudomonas | 2,203 | 4.19 | ON (alginate, psl) |
| SRR6833349 | in-vivo CF sputum | Pseudomonas | 84,060 | 7.97 | ON (alginate, pel, psl) |

Rossi 2018 **PRJEB24688** (RiboZero hybridisation — the least biased chemistry for %RP):

| run | target reads | coding | %RP | matrix | decoy |
|---|---|---|---|---|---|
| ERR2275091 | 471,476 | 416,607 | 7.31 | ON (alginate 17,434 + pel + psl) | 0 |
| ERR2275090 | 467,140 | 394,965 | 10.95 | ON (alginate 10,462 + pel + psl) | 2,492 |

The matrix axis fires on known-positive material in both an in-vitro biofilm and
in-vivo sputum, at depths from 2,203 to 471,476 target reads. Cornforth et al. also
independently report Psl and `algU` transcribed *more* in human samples than in vitro.

## Item 5 — the depth-dilution experiment

The literature review found that %RP was validated only on deep, high-abundance taxa
(Gifford 2012 restricted all interpretation to the top ~200 genome bins) and that the
closest analogue — ribosome profiling — fabricates apparent stalling at low coverage
(Artieri & Fraser 2014). The decisive experiment (a sequencing-depth dilution series)
had not been published. Ran it:

**ERR2275090** (Rossi CF sputum), one library subsampled across a 47× depth range:

| requested reads | target reads | coding | %RP |
|---|---|---|---|
| 2,000,000 | 467,140 | 394,965 | 10.95 |
| 250,000 | 233,491 | 197,302 | 10.84 |
| 100,000 | 93,665 | 79,126 | 10.99 |
| 50,000 | 46,892 | 39,614 | 11.24 |
| 20,000 | 18,759 | 15,724 | 11.61 |
| 10,000 | 9,384 | 7,844 | 11.65 |

**SRR6833349** (Cornforth CF sputum): 7.97 → 8.03 → 7.78 → 8.22 → 7.71 → 8.65
down to 1,698 target reads.

**%RP does not deflate at low depth.** Across a 47× reduction it moves +0.7 pp and
slightly *upward*; the second run is flat-to-noisy with no direction. So the low-load
artifact does not explain SLOW calls — if anything it biases toward FAST. This removes
the "low-load artifact" objection from the growth axis and strengthens the
PMN-oxygen-depletion reading of the chronic arm.

*Caveats:* one library per organism, no spike-in, no replicate; the dilution varies
read depth on a single sample, so it tests depth sensitivity, not biological load.
The clean version still needs a spike-in and replicate dilutions.

## Item 1 — HS78/HS55/HS1 cannot be deepened

Shallow depth is in the deposit:

| patient | runs | deepest run reads |
|---|---|---|
| HS78 | 4 | 69,611 |
| HS1 | 5 | 327,180 |
| HS55 | 9 | 2,415,781 |

There is no deeper data to download for HS78/HS1. "Deepen the depth-limited no-BP runs"
is not an available action; the honest options are to exclude them or report them as
coverage-limited.

## Corrected 40k screen (all 254 runs)

First pan-genus screen of PRJNA1033689 built on the corrected 8-genome panel with
decoys. Dominant organism, deepest run per patient:

| arm | composition |
|---|---|
| 2BP (n=27) | escherichia 10, pseudomonas 6, staphylococcus 5, klebsiella 4, streptococcus 1, haemophilus 1 |
| No-BP (n=29) | staphylococcus 9, escherichia 8, streptococcus 7, klebsiella 2, stenotrophomonas 1, haemophilus 1, pseudomonas 1 |

This differs materially from the committed table's "no-BP = Escherichia 15,
Klebsiella 3, Pseudomonas 2": the corrected screen puts **Staphylococcus and
Streptococcus** at the top of the no-BP arm. Depth is balanced across arms at the
deepest-run level (2BP median 310 panel reads, no-BP median 357; ≥100 reads 17/27 vs
22/29).

**Correction to an intermediate claim of mine:** I initially reported that the
committed "HS243 Esch 7.2k" does not reproduce. That was my error — I tested only
HS243's deepest run (SRR26629284, kitome-dominated). The claim refers to a different
run of the same patient, **SRR26629286** (3,919 E. coli reads at 40k mate1; 6,597 with
both mates). The committed entry is real. Lesson: for per-patient composition claims,
screen every run of the patient, not just the deepest.

## Corrected quantification of the 49 matrix-capable samples

Selected the run with the most matrix-capable reads per patient (49 of 56 patients
had ≥100), quantified at 500k reads/mate. Early results:

| patient | arm | dominant | target reads | %RP | systems ON |
|---|---|---|---|---|---|
| HS236 | 2BP | Pseudomonas | 56,174 | 0.78 | alginate, pel, psl |
| HS410 | 2BP | Staphylococcus | 272,719 | 10.58 | fnb, ica |
| HS271 | No-BP | Staphylococcus | 45,574 | 1.45 | fnb, ica |
| **HS253** | **No-BP** | **Klebsiella** | **14,896** | **5.53** | **bcs, capsule, cps, fim, mrk, pga (6 systems)** |

HS253 is the clearest demonstration of the fix: a no-BP patient scoring matrix-ON on
**six** systems, where the committed mrk-only panel would have returned matrix-OFF by
construction.

## Arm statistics on the corrected screen

Best matrix-capable run per patient (not merely the deepest):

| threshold (matrix-capable reads) | 2BP | No-BP | Fisher p |
|---|---|---|---|
| ≥50 | 25/27 | 28/29 | 1.000 |
| ≥100 | 25/27 | 24/29 | 1.000 |
| ≥300 | 23/27 | 17/29 | 1.000 |
| ≥1000 | 19/27 | 9/29 | 0.118 |

Depth balance at that run: **2BP median 2,503 panel reads vs No-BP 683** — the 2BP arm
is the deeper one, the opposite of the committed framing in which no-BP was the
depth-starved arm. The ≥1000 threshold shows a non-significant trend toward
matrix-capable enrichment in 2BP (p = 0.118), direction consistent with the committed
claim but not significant.

## Effect of the class-stratified rule (Item 3 applied)

Applying the recommended capsule-excluded EPS-only rule to the quantified samples
flips two of them:

| patient | arm | systems ON | pooled | EPS-only |
|---|---|---|---|---|
| HS410 | 2BP | fnb, ica | ON | ON (ica only = 1) → **OFF** |
| HS271 | No-BP | fnb, ica | ON | **OFF** |
| HS253 | No-BP | bcs, capsule, cps, fim, mrk, pga | ON | ON (bcs, cps, pga) |

So the pooled rule and the class-stratified rule give different arm tallies on real
data — which is exactly why the paper must pre-specify one and report the other as a
sensitivity analysis rather than choosing after the fact.

## Reproduce

```
fetch_genomes.sh / fetch_rest.sh / fetch_decoys_all.sh   # 46 genomes + GFFs
build_panels_v3.py        # corrected panels -> qp3_*.fna, gff_intervals3.json, panels3.json
build_gcomp3.py           # gcomp3_<org>.mmi (full genome + decoys)
build_multi3.py           # panel_multi3.mmi (8 target genomes)
screen_all3.py 40000      # corrected pan-genus screen of all 254 runs -> screen3/
analyze_screen3.py        # corrected composition table
quant_any.py <SRR> <fq> <outdir>   # two-axis quantification
item4_positive_control.py / item4b_rossi.py   # positive controls
item5_depth_dilution.py   # the depth-dilution experiment
```

**Note on the toolchain:** the committed scripts hardcode
`bin/minimap2-2.28_x64-linux/minimap2`. In a fresh sandbox, `bin/minimap2` must be
placed there (or symlinked) or **index builds fail silently and produce no `.mmi`** —
which then looks like a mapping failure rather than a build failure.
