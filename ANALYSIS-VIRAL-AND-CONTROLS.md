# Non-COVID viral pneumonia + healthy controls — the two-axis analysis

**Date:** 2026-10-12. Answers two questions: (1) can the two-axis (growth × matrix) test be run on
non-COVID viral pneumonia? (2) where do healthy controls fall on those axes?

**Headline: the matrix axis is not attainable in PRJNA699976 (2 of 78 patients), the growth axis is
attainable but only in the kitome-confounded way, and the sharpest new result is that a colonisation
control places on both axes while healthy and water blanks do not.**

---

## Part 1 — PRJNA699976 at scale (110 samples screened)

The cohort: **408 pneumonia/ARI patients + 28 healthy volunteers + 4 RNase-free-water blanks**, all
BALF, all ~20M reads, pre-COVID (Wuhan). Screened 78 patients (deepest-first) + all 28 healthy + all
4 blanks.

### 1.1 The kitome floor dominates everything

| category | n | median mapped reads | dominant genus | matrix-capable reads (per 120k) |
|---|---|---|---|---|
| RNase-free water | 4 | — | **Ralstonia pickettii, 99%+** | **17.0** (sd 4.9) |
| healthy volunteers | 28 | — | *Ralstonia* 88% | **7.5** |
| pneumonia patients | 78 | — | *Ralstonia* 99% | ~background |

**The water blanks have the highest median matrix-capable read count of the three categories.** The
healthy volunteers are *not distinguishable from water*. And 60 of 78 patients are also
kitome-dominated.

### 1.2 The apparent contradiction, resolved

`HEALTHY-CONTROLS.md` reported the healthy arm as "~80% bacterial." The composition screen reports
it as kitome-dominated. Both are correct and they are not in conflict:

- **80% of *deposited* reads are bacterial** — the deposit is host-removed, so what remains is
  nearly all bacterial.
- **But those bacterial reads are almost all *Ralstonia*** — a reagent contaminant, not lung
  biology.

So the healthy arm has plenty of bacterial signal, and almost none of it is informative. This is the
same kitome floor as PRJNA390194, reached by a different route: there, host reads swamped the
sample; here, host reads were removed but the remainder is reagent.

### 1.3 Matrix-capable organisms: 2 of 78

| threshold | samples |
|---|---|
| ≥100 matrix-capable reads | **3 / 110** (2 patients, 1 "healthy") |
| ≥1000 matrix-capable reads | **0 / 110** |
| above the water-blank maximum (22 reads) | 9 / 110, of which only **2** clear it decisively |

The two decisive hits, both validated by re-probing at 400k reads (3.33×) with linear scaling:

| sample | organism | 120k reads | 400k reads | projected full depth |
|---|---|---|---|---|
| **SRR13677669** | *Acinetobacter* | 489 | 1,693 | ~1.9×10⁵ |
| **SRR13677832** | ***Pseudomonas*** | 244 | 871 | — |

**SRR13677832 is the only patient sample in the cohort whose single dominant genus is
matrix-capable.** Ranks 3–20 are all within 1–2× of the water blanks and are reagent background.

**Verdict on the matrix axis here: n=2 gives no power.** The dataset supports case studies, not a
cohort analysis. And because there is no bacterial rRNA depletion, those counts are an *upper bound*
on coding signal — whether the ~10⁵ Pseudomonas reads in SRR13677832 survive as %RP / matrix-operon
signal is the next check before either case study is usable.

### 1.4 The search landscape, corrected

The earlier "no non-COVID viral pneumonia data" verdict was **wrong** — it applied a criterion
(documented bacterial rRNA depletion) that the matrix axis does not need. Re-search found ≥15
non-COVID respiratory metatranscriptomes with ≥50 samples, several at 300–440 samples and 17–55M
reads. Top candidates not yet probed:

| dataset | what | depth | note |
|---|---|---|---|
| **PRJNA1299073** | *P. aeruginosa* pneumonia n=7 vs 8 vs 6 control | 23.8M median | **documented bacterial depletion (Inverse rRNA)** — the cleanest available matrix **positive control** |
| HRA000143 (NGDC) | 147 non-COVID CAP BALF | deep | NGDC was the earlier blind spot |
| PRJNA1165951 | 209 ARDS/non-ARDS BALF | 0.31M median | GEO-linked host data |
| PRJNA1240018 | 291 children, **sputum** | 5.6M median | best paediatric option |

What remains true: **no published non-COVID viral pneumonia cohort *reports* matrix-gene
transcription.** That is a "not yet measured," not a "cannot be measured."

---

## Part 2 — Healthy controls on the two axes

### 2.1 The PRJNA390194 result holds and is depth-independent

9 non-COPD BALF samples: **0 Pseudomonas panel reads**, 369 decoy reads. Re-probed at 20.4M reads
(SRR5677656): still only 14 PAO1 coding reads. Not a depth artifact.

**But the generalisation was wrong.** "Healthy lung = kitome floor" does not follow — PRJNA699976's
healthy arm has abundant bacterial reads (see §1.2). The correct claim is narrow:
**PRJNA390194 is a kitome floor; healthy lung is not necessarily one.**

### 2.2 The colonisation control places on BOTH axes

The most consequential finding. Running PRJDB41581 colonisation sample **DRR1001389** through the
repo's existing *S. pneumoniae* panel (TIGR4; `cps` 12 genes, `cbp` 8 genes):

| metric | value |
|---|---|
| target reads | 32,010 |
| f_rp | **10.6% → FAST** |
| matrix_ON | **TRUE** (`cps` 8/12, `cbp` 6/8) |

A **control** sample lands in the acute quadrant.

Two honest caveats:
- A second colonisation sample (DRR1001381) does **not** place — depth-limited, same failure mode as
  the acute *Pseudomonas* arm.
- Healthy SRR13677570 has 4,142 *S. pneumoniae* target reads but **0 matrix reads** — organism
  present, capsule machinery not transcribed. A genuine biological negative.

### 2.3 Why controls look "matrix-off" — the organism explanation

The controls' dominant organisms are oral streptococci and *Prevotella*. These are **not
matrix-capable in the *Pseudomonas* sense** (alginate/psl/pel), which is precisely why healthy and
colonisation samples appear matrix-negative on a *Pseudomonas* panel. Placement is
**organism-dependent**, not infection-dependent.

⚠️ **Interpretation caveat before this is used:** the streptococcal capsule (`cps`) is *inversely*
related to biofilm in the literature — a point already flagged in `ORGANISM-VIABILITY.md`. So a
`matrix_ON` call driven by `cps`/`cbp` may mean the opposite of what a *Pseudomonas* `matrix_ON`
call means. **This must be resolved before the control-placement result is claimed.**

### 2.4 Best control type

**Colonised-but-not-infected** — clinically matched to the infected arm *and* carrying bacteria:
PRJDB41581 (pneumonia 27 vs colonisation 17), PRJNA875913 (IPC 70 vs non-infectious 49). Healthy
volunteers are the wrong control for a bacterial-signal study unless paired with a kitome
denominator — which PRJNA699976 provides (4 water blanks at 52–57M reads).

---

## Part 3 — Errors found and corrected in this pass

| # | error | status |
|---|---|---|
| 1 | **`acc2genus.json` is mislabeled for 4 of 19 genomes** — *A. baumannii* 17978 → "pseudomonas", *S. aureus* NCTC8325 → "streptococcus", *S. pneumoniae* TIGR4 → "haemophilus", *H. influenzae* → "staphylococcus" | **verified by FASTA headers.** Impact tested: see below |
| 2 | PRJNA699976 `SRR13677601` listed as a patient run; it is an **RNase-free water blank** | corrected in `PRJNA699976-PROBE.md` |
| 3 | PRJDB41581 rejected on "1.1% mapping"; the committed probe was an all-zeros failed run. Re-probe: 19.7–50.8% | corrected in `RESULTS-ACUTE-SPUTUM-SEARCH.md`; prep objection stands |
| 4 | "healthy lung = kitome floor" generalisation | corrected above |

### Impact assessment of the `acc2genus.json` error

The four mislabeled genomes are **absent from every index the growth and matrix axes use**
(`compfull_pseu.mmi`, `comp_pseu.mmi`, `gcomp_pseu.mmi`). They appear only in `panel_multi.fna`
(used solely for *genus assignment* by `quant_perorg.py` / `quant_final.py`) and `compall.fna`
(used only by `rep_screen.py`).

**Decisive test:** for every sample in the acute *Pseudomonas* arm, the panel's "pseudomonas" count
was compared against the count of reads mapping to real PAO1 (the only target in
`compfull_pseu.mmi`). A sample that is actually *A. baumannii* would show a panel count with almost
no PAO1 mapping.

| result | samples |
|---|---|
| ratio PAO1/panel ≈ 3.5–6.7 (consistent; PAO1 genome ≈ 6× the RP+matrix panel) | **14 of 15** |
| ratio 0.09 — **flagged** | **SRR27343439** |

**SRR27343439** (panel 513 "pseudomonas" vs 46 real PAO1 reads) is a candidate mislabeled
*A. baumannii*. It **is** in the headline acute arm, but it is the shallowest sample there
(46 PAO1 coding reads).

**Effect on the result: negligible.** Recomputing without it:

| analysis | median acute %RP | acute vs chronic p |
|---|---|---|
| all 15 | 10.68% | 0.0015 |
| excluding SRR27343439 | 10.80% | 0.0019 |

**Verdict: the headline survives, but `acc2genus.json` must be corrected and the mislabeling
disclosed as a methods note** — a wrong genus map that silently swaps *A. baumannii* and
*Pseudomonas* is exactly the class of error this paper is about.

---

## Part 4 — What to do next, ranked

1. **Fix `acc2genus.json`** from FASTA headers and re-run any genus-dependent script. Cheap, and it
   removes a latent error.
2. **Probe PRJNA1299073** — documented bacterial depletion plus a *P. aeruginosa* pneumonia arm. This
   is the one dataset that can serve as a **matrix-axis positive control** and calibrate the metric.
3. **Resolve the `cps`-inversion caveat** before the control-placement finding is used at all.
4. **Probe PRJNA699976's two case-study samples** (SRR13677832 *Pseudomonas*, SRR13677669
   *Acinetobacter*) for whether the signal survives as %RP / matrix-operon reads.
5. **Probe HRA000143 (147 CAP BALF, NGDC)** — the largest unexamined non-COVID CAP resource, and
   NGDC was the earlier blind spot.

---

*Scripts: `prj_genus.py`, `probe_composition2.py`, `work699976/` drivers.
Data: `out/prjna699976_composition.json`, `out/DRR1001389_strep_matrixprobe.json`,
`out/SRR13677570_strep_matrixprobe.json`. Companions: `HEALTHY-CONTROLS.md`, `PRJNA699976-PROBE.md`,
`PRJNA699976-COMPOSITION.md`.*
