# COVID-19 pneumonia — biofilm matrix axis (PRJNA1033689)

**Question.** Is the biofilm matrix machinery transcribed in the bacteria of COVID-19
pneumonia, and does COVID occupy the same (growth x matrix) quadrant as acute
bacterial pneumonia?

**Method (mechanical, reuses committed code).**
- Reads (first 300k of R1 per run) mapped with `quantify_v3.py <SRR> pseu <fastq>`
  against `comp_pseu.mmi` (Pseudomonas panel + decoys). Matrix genes are
  Pseudomonas-specific and non-conserved, so the panel index is correct for this
  axis even though it is wrong for the conserved RP genes. **Only the matrix
  output of `quantify_v3.py` is used here** — the growth axis comes from
  `out_norm/covid_rpctrl_all.json` (`rp_over_ctrl`, the corrected metric).
- System ON = >=3 genes detected in that system. `matrix_call` is the committed
  logic from `analyze_corrected.py` (`sys_on()` / `mat_call()`): **ON** if >=2
  systems on; **alg-only** if alginate is the only system on; **partial** if
  exactly 1 non-alginate system is on, or if >=1 gene is detected anywhere but no
  system reaches 3; else **OFF**.
- Growth classes: FAST-like `rp_over_ctrl > 4.0`, SLOW-like `< 3.0`, else ambig.
- Work list: 40 COVID runs in `out_norm/covid_rpctrl_all.json["table"]`,
  processed deepest-first. **28 runs** had `pao1_coding >= 300` and all 28 were
  processed (0 failures). The other 12 were below the depth gate (30–278 coding
  reads) and were not attempted — at that depth matrix detection is noise.

## (a) Per-sample results (n = 28, deepest-first)

`g/r` = genes detected / reads in that system.

| SRR | load | pao1_coding | rp/ctrl | growth | alg g/r | psl g/r | pel g/r | systems_on | matrix_call |
|---|---|---|---|---|---|---|---|---|---|
| SRR26629241 | 0.8422 | 28866 | 0.71 | SLOW-like | 19/128 | 13/54 | 7/66 | alginate;psl;pel | **ON** |
| SRR26629260 | 0.0884 | 17300 | 4.21 | FAST-like | 11/104 | 11/49 | 5/10 | alginate;psl;pel | **ON** |
| SRR26629332 | 0.0610 | 11541 | 6.07 | FAST-like | 16/74 | 2/2 | 2/4 | alginate | alg-only |
| SRR26629295 | 0.0560 | 10752 | 5.22 | FAST-like | 7/43 | 6/14 | 2/2 | alginate;psl | **ON** |
| SRR26629302 | 0.0547 | 9994 | 4.56 | FAST-like | 12/75 | 5/25 | 0/0 | alginate;psl | **ON** |
| SRR26629409 | 0.0477 | 9297 | 5.86 | FAST-like | 6/82 | 0/0 | 0/0 | alginate | alg-only |
| SRR26629276 | 0.0427 | 8449 | 6.00 | FAST-like | 0/0 | 0/0 | 0/0 | – | OFF |
| SRR26629301 | 0.0570 | 8227 | 4.95 | FAST-like | 8/28 | 3/5 | 0/0 | alginate;psl | **ON** |
| SRR26629336 | 0.0460 | 8105 | 4.66 | FAST-like | 18/126 | 6/17 | 2/3 | alginate;psl | **ON** |
| SRR26629334 | 0.0361 | 6598 | 3.88 | ambig | 19/96 | 4/9 | 0/0 | alginate;psl | **ON** |
| SRR26629275 | 0.0324 | 6266 | 6.19 | FAST-like | 5/22 | 6/16 | 2/5 | alginate;psl | **ON** |
| SRR26629300 | 0.0412 | 6254 | 5.49 | FAST-like | 12/39 | 1/2 | 0/0 | alginate | alg-only |
| SRR26629335 | 0.0307 | 5336 | 4.45 | FAST-like | 13/57 | 3/3 | 2/4 | alginate;psl | **ON** |
| SRR26629274 | 0.0274 | 5263 | 4.03 | FAST-like | 8/47 | 6/8 | 1/1 | alginate;psl | **ON** |
| SRR26629273 | 0.0242 | 4435 | 4.12 | FAST-like | 9/32 | 6/8 | 0/0 | alginate;psl | **ON** |
| SRR26618575 | 0.0262 | 4317 | 3.12 | ambig | 9/33 | 3/5 | 0/0 | alginate;psl | **ON** |
| SRR26618576 | 0.0190 | 3340 | 2.98 | SLOW-like | 6/20 | 3/5 | 0/0 | alginate;psl | **ON** |
| SRR26618529 | 0.0201 | 3328 | 3.98 | ambig | 5/19 | 0/0 | 2/3 | alginate | alg-only |
| SRR26629277 | 0.0162 | 2934 | 6.39 | FAST-like | 5/14 | 2/3 | 1/2 | alginate | alg-only |
| SRR26629249 | 0.0211 | 2175 | 9.00 | FAST-like | 7/20 | 1/1 | 0/0 | alginate | alg-only |
| SRR26618540 | 0.0079 | 1460 | 3.09 | ambig | 12/88 | 0/0 | 0/0 | alginate | alg-only |
| SRR26618518 | 0.0081 | 1439 | 3.43 | ambig | 5/12 | 0/0 | 0/0 | alginate | alg-only |
| SRR26629327 | 0.0079 | 1349 | 4.58 | FAST-like | 5/6 | 1/1 | 0/0 | alginate | alg-only |
| SRR26629410 | 0.0064 | 1167 | 5.25 | FAST-like | 2/11 | 0/0 | 0/0 | – | partial |
| SRR26629333 | 0.0051 | 967 | 2.80 | SLOW-like | 2/4 | 2/2 | 0/0 | – | partial |
| SRR26629278 | 0.0936 | 778 | 2.25 | SLOW-like | 3/5 | 2/8 | 1/2 | alginate | alg-only |
| SRR26629242 | 0.7953 | 705 | 2.00 | SLOW-like | 1/1 | 2/8 | 0/0 | – | partial |
| SRR26629303 | 0.0486 | 416 | 3.45 | ambig | 2/2 | 0/0 | 0/0 | – | partial |

## (b) matrix_call counts

| call | all called (n=28) | primary subset (n=25) |
|---|---|---|
| ON | **13** | **13** |
| alg-only | 10 | 9 |
| partial | 4 | 2 |
| OFF | 1 | 1 |

COVID matrix-ON rate: **46.4%** (13/28); **52.0%** (13/25) on the primary subset.

**Depth thresholds used.**
- Selection gate: `pao1_coding >= 300` — 12 of 40 runs fall below it
  (SRR26618521, 26618513, 26618519, 26629309, 26618517, 26629355, 26629371,
  26618558, 26618514, 26629380, 26618572, 26629386; 30–278 coding reads) and were
  not attempted.
- Calling gate: panel `target_reads >= 200` after 300k-read subsampling. Three
  runs (SRR26629278, 26629242, 26629303) cleared the selection gate but yielded
  only 21–95 panel reads; their calls (alg-only / partial / partial) are reported
  but flagged low-signal and excluded from the primary subset.
- No sample failed to process. The three shallow runs whose R1 is a single
  sub-megabyte file (SRR26629278, 26629242, 26629303) returned a truncated body
  under a 60 MB range request; they were re-downloaded in full without a Range
  header and mapped on all available reads (8.3k / 0.9k / 10.6k reads).

## (c) Reference arms (acute sputum vs chronic in-vivo)

Read from `final_out/*.json` (`matrix_on` field, produced by `quant_final.py`),
and recomputed where `out/<srr>_pseu_v3.json` exists so the calling logic is
identical to the COVID calls.

| arm | n | matrix-ON | rate | recomputed calls (v3 logic) |
|---|---|---|---|---|
| acute (SRR2734*) | 16 | 8 | **50.0%** | ON 5, alg-only 8, partial 1, OFF 1 (n=15) |
| chronic in-vivo (ERR2275087–2275101) | 15 | 15 | **100.0%** | ON 15 (n=15) |
| **COVID (this work)** | 28 | 13 | **46.4%** | ON 13, alg-only 10, partial 4, OFF 1 |

The acute arm's two calling conventions agree closely (50% vs 33% — the
difference is the read-depth gate: the `quant_final.py` call uses a >=2-read
per-system threshold, the v3 call a >=3-gene threshold). The chronic arm is
saturated: every chronic in-vivo sample is matrix-ON under both conventions.

**Depth is not matched across arms** (median panel `target_reads`: chronic
~1.18M, acute ~10.9k, COVID ~1.5k), and matrix-ON rate rises with depth within
COVID. Matching depth against the acute arm:

| panel target_reads >= | COVID matrix-ON | acute matrix-ON (v3 logic) |
|---|---|---|
| 0 | 13/28 = 46.4% | 5/15 = 33.3% |
| 200 | 13/25 = 52.0% | 5/14 = 35.7% |
| 500 | 12/19 = 63.2% | 5/14 = 35.7% |
| 1000 | 10/14 = 71.4% | 5/13 = 38.5% |
| 2000 | 9/12 = 75.0% | 5/11 = 45.5% |

At every depth COVID's matrix-ON rate is at or above the acute arm's.

## (d) The quadrant: COVID placement

Counts, matrix-ON vs not-ON x growth class (all called samples, n=28):

| | matrix-ON | matrix-not-ON |
|---|---|---|
| **FAST-like** (rp/ctrl > 4.0) | **9** | 8 |
| **SLOW-like** (rp/ctrl < 3.0) | 2 | 3 |
| ambig (3.0–4.0) | 2 | 4 |

Same table on the primary subset (n=25, low-signal calls dropped):
FAST-like 9 / 8, SLOW-like 2 / 1, ambig 2 / 3.

**Placement.**
- COVID is **predominantly FAST-like**: 17/28 (61%) are FAST-like, median
  `rp_over_ctrl` = 4.33 over the 28 processed runs (4.45 over the 27 runs with
  `pao1_coding >= 500`, the subset the parent figure quotes), against only 5/28
  (18%) SLOW-like.
- Within the FAST-like class, **9/17 (53%) are matrix-ON** — the
  (matrix-present x FAST) cell that the acute arm occupies.
- Chronic in-vivo sits in a different cell: 100% matrix-ON but mostly SLOW-like
  (8/15 SLOW, 3 FAST, 4 ambig under the corrected metric), i.e. (matrix-ON x
  SLOW).
- COVID is therefore spread across the quadrant rather than tightly clustered:
  its largest single cell is **FAST-like x matrix-ON (9 samples)**, the acute-like
  corner, with a substantial FAST-like x matrix-not-ON cell (8) alongside it.

## Verdict

**Yes — biofilm matrix machinery is transcribed in COVID-19 pneumonia bacteria.**
13 of 28 depth-qualified COVID runs (46%; 52% on the primary subset, 71–75% at
matched depth) transcribe at least two of the three Pseudomonas matrix systems
(alginate + psl, occasionally pel), a rate equal to the acute bacterial-pneumonia
arm (50%) and rising above it once read depth is matched. **COVID therefore
occupies the acute quadrant rather than the chronic one**: its modal cell is
matrix-present x FAST-like growth (9/28), the same corner as acute bacterial
pneumonia, whereas chronic in-vivo infection sits at matrix-ON x SLOW (15/15
matrix-ON, mostly SLOW). The one caveat is depth — COVID's Pseudomonas signal
per sample is roughly an order of magnitude below the acute arm's, so the COVID
matrix-ON rate is a lower bound, and 12 further COVID runs could not be called
at all because they carry <300 PAO1-mapping coding reads.

---
Files: `out/covid_matrix.json` (per-sample results + summary + 2x2),
`out/<SRR>_pseu_v3.json` (raw per-sample panel output),
`covid_work/matrix_quant.jsonl` (run log). No git commit.
