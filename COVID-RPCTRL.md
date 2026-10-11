# COVID cohort rp/ctrl re-quantification (PRJNA1033689)

Re-quantified **all 40 Pseudomonas-positive COVID runs** from PRJNA1033689 with the
abundance-robust control-gene metric

```
rp_over_ctrl = PAO1 ribosomal-protein reads / PAO1 reads over the 23 equally-conserved
               non-RP control genes      (MAPQ>=20, rRNA/tRNA/tmRNA excluded)
```

built into the committed `norm_map.py`. Each run was quantified at constant depth
(first 200 000 reads of the R1 fastq, or all reads if the run is shallower), exactly as
in the pilot. Driver: `covid_rpctrl_pipe.py` (resumable, deepest-first); analysis:
`covid_rpctrl_analyze.py`; per-sample JSON: `out_norm/<SRR>_covid_norm.json`;
aggregate: `out_norm/covid_rpctrl_all.json`.

**n completed: 40 / 40.** 5 runs have zero control-gene reads at 200 k depth (load
<=0.003, coding <=219 reads) so `rp_over_ctrl` is undefined; 35 have a finite value and
27 of those are high-depth (coding >=500 reads). All 6 pilot samples reproduced
**bit-identically** (same 200 k subsample window).

## (a) Per-sample table (sorted by load)

| srr | read_count | load | coding | rp | ctrl | pctRP_coding | rp/ctrl | class |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| SRR26618514 | 171,015 | 0.000801 | 137 | 6 | 0 | 4.38 | — | no_ctrl |
| SRR26618558 | 175,483 | 0.000935 | 164 | 13 | 0 | 7.93 | — | no_ctrl |
| SRR26629355 | 300,706 | 0.001015 | 201 | 14 | 5 | 6.97 | 2.80 | SLOW |
| SRR26629371 | 211,301 | 0.001060 | 185 | 3 | 9 | 1.62 | 0.33 | SLOW |
| SRR26629309 | 779,312 | 0.001095 | 219 | 0 | 0 | 0.00 | — | no_ctrl |
| SRR26618517 | 2,107,859 | 0.001115 | 217 | 18 | 3 | 8.30 | 6.00 | FAST |
| SRR26618513 | 9,647,335 | 0.001170 | 231 | 13 | 11 | 5.63 | 1.18 | SLOW |
| SRR26618519 | 1,373,814 | 0.001215 | 220 | 33 | 2 | 15.00 | 16.50 | FAST |
| SRR26618521 | 7,340,428 | 0.001475 | 278 | 19 | 14 | 6.84 | 1.36 | SLOW |
| SRR26618572 | 42,352 | 0.001629 | 69 | 0 | 0 | 0.00 | — | no_ctrl |
| SRR26629380 | 32,782 | 0.002623 | 86 | 0 | 3 | 0.00 | 0.00 | SLOW |
| SRR26629386 | 10,386 | 0.003081 | 30 | 0 | 0 | 0.00 | — | no_ctrl |
| SRR26629333 | 6,888,356 | 0.005055 | 967 | 168 | 60 | 17.37 | 2.80 | SLOW |
| SRR26629410 | 2,505,552 | 0.006370 | 1167 | 147 | 28 | 12.60 | 5.25 | FAST |
| SRR26618540 | 11,469,268 | 0.007855 | 1460 | 164 | 53 | 11.23 | 3.09 | amb |
| SRR26629327 | 2,142,173 | 0.007895 | 1349 | 325 | 71 | 24.09 | 4.58 | FAST |
| SRR26618518 | 1,017,433 | 0.008080 | 1439 | 158 | 46 | 10.98 | 3.43 | amb |
| SRR26629277 | 1,678,502 | 0.016150 | 2934 | 575 | 90 | 19.60 | 6.39 | FAST |
| SRR26618576 | 6,510,113 | 0.019015 | 3340 | 185 | 62 | 5.54 | 2.98 | SLOW |
| SRR26618529 | 3,674,033 | 0.020135 | 3328 | 215 | 54 | 6.46 | 3.98 | amb |
| SRR26629249 | 4,682,128 | 0.021145 | 2175 | 288 | 32 | 13.24 | 9.00 | FAST |
| SRR26629273 | 5,302,595 | 0.024165 | 4435 | 774 | 188 | 17.45 | 4.12 | FAST |
| SRR26618575 | 7,646,598 | 0.026235 | 4317 | 250 | 80 | 5.79 | 3.13 | amb |
| SRR26629274 | 6,566,290 | 0.027375 | 5263 | 1117 | 277 | 21.22 | 4.03 | FAST |
| SRR26629335 | 1,054,090 | 0.030670 | 5336 | 922 | 207 | 17.28 | 4.45 | FAST |
| SRR26629275 | 9,420,194 | 0.032375 | 6266 | 1126 | 182 | 17.97 | 6.19 | FAST |
| SRR26629334 | 1,096,507 | 0.036130 | 6598 | 1000 | 258 | 15.16 | 3.88 | amb |
| SRR26629300 | 2,269,455 | 0.041195 | 6254 | 1351 | 246 | 21.60 | 5.49 | FAST |
| SRR26629276 | 1,346,587 | 0.042675 | 8449 | 156 | 26 | 1.85 | 6.00 | FAST |
| SRR26629336 | 2,961,210 | 0.046045 | 8105 | 1483 | 318 | 18.30 | 4.66 | FAST |
| SRR26629409 | 557,562 | 0.047710 | 9297 | 1419 | 242 | 15.26 | 5.86 | FAST |
| SRR26629303 | 10,633 | 0.048622 | 416 | 69 | 20 | 16.59 | 3.45 | amb |
| SRR26629302 | 466,656 | 0.054695 | 9994 | 2328 | 511 | 23.29 | 4.56 | FAST |
| SRR26629295 | 6,465,279 | 0.055955 | 10752 | 2882 | 552 | 26.80 | 5.22 | FAST |
| SRR26629301 | 547,612 | 0.056975 | 8227 | 1501 | 303 | 18.25 | 4.95 | FAST |
| SRR26629332 | 3,463,396 | 0.061030 | 11541 | 3721 | 613 | 32.24 | 6.07 | FAST |
| SRR26629260 | 1,666,230 | 0.088395 | 17300 | 2976 | 707 | 17.20 | 4.21 | FAST |
| SRR26629278 | 8,344 | 0.093600 | 778 | 9 | 4 | 1.16 | 2.25 | SLOW |
| SRR26629242 | 889 | 0.795276 | 705 | 4 | 2 | 0.57 | 2.00 | SLOW |
| SRR26629241 | 34,501 | 0.842207 | 28866 | 171 | 241 | 0.59 | 0.71 | SLOW |

`load = pao1_total / n_reads_analysed`; `class` from `rp_over_ctrl` using the pilot arm
anchors (FAST > 4.0, SLOW < 3.0).

## (b) Spearman correlations

| subset | metric | n | rho | p (two-sided) |
|---|---|---:|---:|---:|
| coding >=500 (primary) | rp_over_ctrl | 27 | **-0.03** | 0.88 |
| coding >=500 (primary) | pctRP_coding | 27 | **-0.03** | 0.89 |
| all finite (coding >=30) | rp_over_ctrl | 35 | +0.15 | 0.38 |
| all finite (coding >=30) | pctRP_coding | 35 | +0.24 | 0.17 |

There is **no monotone load dependence** of either metric at full n. The pilot's
apparent gradient was a two-point contrast, not a trend.

## (c) Classification (all 40 runs)

| class | n | load range | median load | rp/ctrl range |
|---|---:|---|---:|---|
| FAST-like (>4.0) | 19 | 0.00112 – 0.08840 | 0.0324 | 4.03 – 16.50 |
| ambiguous (3.0–4.0) | 6 | 0.00786 – 0.04862 | 0.0232 | 3.09 – 3.98 |
| SLOW-like (<3.0) | 10 | 0.00102 – 0.84221 | 0.0038 | 0.00 – 2.98 |
| undefined (no ctrl reads) | 5 | 0.00080 – 0.00308 | 0.0011 | — |

On the high-depth subset (coding >=500, n=27): **FAST 17, ambiguous 5, SLOW 5**
(load ranges: FAST 0.0064–0.0884; ambiguous 0.0079–0.0361; SLOW 0.0051–0.8422).

Crucially, the SLOW class is **not** a high-load class. Only 2 of the 10 SLOW calls
(SRR26629241, SRR26629242) are the near-monoculture samples; the other 8 are ordinary
low-load samples (load 0.001–0.094) scattered right through the FAST cloud.

## (d) Gradient test (load < 0.1 vs load > 0.5)

| group | n | mean rp/ctrl | median |
|---|---:|---:|---:|
| load < 0.1 | 33 | 4.49 | 4.21 |
| load > 0.5 | 2 | 1.35 | 1.35 |

Mann-Whitney two-sided: U = 6.0, z = -1.92, **p = 0.055**. Directionally as the pilot
predicted, but at n_high = 2 this is not a usable test — and the broader high group
(load > 0.2) is empty, so there is no intermediate-load evidence for a gradient.

## (e) Pilot sanity

All 6 pilot samples reproduced their earlier `out_norm/<SRR>_norm.json` values
**exactly** (n_reads_analysed, coding, rp, ctrl, pctRP, rp_over_ctrl, load all
bit-identical), confirming the same 200 k subsample window was used and that the new
`_covid_norm.json` files are directly comparable to the pilot. No discrepancies.

## Verdict

**The gradient does not hold at full n: across all 40 COVID runs, Pseudomonas
ribosomal-protein fraction is essentially independent of bacterial load (Spearman
rho = -0.03, p = 0.88 on the 27 high-depth samples), so there is no evidence that these
populations grow faster at low burden and slower at high burden.** The one reproducible
feature is a *threshold* effect, not a gradient: the two near-monoculture samples
(load 0.80–0.84) do show a collapsed rp/ctrl (0.71 and 2.00, chronic-like), but every
sample in between — from load 0.001 to 0.094 — is a mixed population whose rp/ctrl
spans the whole acute-to-chronic range (0.0 to 16.5) with no load trend, and only 2 of
the 10 SLOW-like calls are high-load. **Honest bottom line: burden predicts growth
state only in the extreme (near-monoculture high-load tail); through the bulk of the
cohort it does not, so the pilot's "load-growth gradient" was an artefact of contrasting
only the two extreme samples, and the rp/ctrl metric is not load-confounded in the way
the pilot suggested — it is simply very noisy at low Pseudomonas depth.**

## Caveats

- 5/40 runs (all load <=0.003, coding <=219 reads) have zero control-gene reads at
  200 k depth, so rp/ctrl is undefined for them; another 8 have coding <500. All
  correlation/classification claims are therefore reported on the coding>=500 subset
  (n=27) with the full finite set (n=35) as a sensitivity check.
- rp/ctrl has a coarse denominator at low depth (single-digit ctrl counts give values
  like 6.00 = 18/3, 16.50 = 33/2), so much of the low-load scatter is counting noise.
- The COVID cohort here is dominated by low-load samples; the "high-load" anchor rests
  on 2 samples, which is why the gradient test cannot be powered even if the effect is real.
