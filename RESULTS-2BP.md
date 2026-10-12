# 2BP enrichment test — PRJNA1033689 / GSE246795

**Date:** 2026-10-12. Test: does the Pseudomonas burden and the two-axis (growth x matrix) state
distinguish adjudicated secondary bacterial pneumonia (2BP) from no-BP in the Langelier cohort?

## The join is exact and complete

GEO `GSE246795` carries a per-sample `2bp_status: 2BP / No-BP` characteristic, and the labels are
**patient-level**: 27 2BP / 29 no-BP. GEO sample names (`MVIR1-HS161-D10ETA1-MNGS1`) are embedded
verbatim in the ENA BioSample aliases of the 254 SRA runs (PRJNA1033689), so the join is exact:
all 40 quantified runs map to labelled patients (30 runs / 8 patients 2BP; 10 runs / 9 patients
no-BP). All 177 ETA runs: 98 2BP / 79 no-BP runs, 27 / 29 patients. 76 NAS runs: 28 / 48.

## Claim A — Pseudomonas burden enriches in 2BP (SURVIVES)

At matched screening depth (full 40k reads/run, all 177 tracheal aspirates):

| threshold | 2BP (n=89 runs) | no-BP (n=61 runs) | Fisher p |
|---|---|---|---|
| >=5 reads | 30% | 16% | 0.038 |
| >=30 | 22% | 3% | 6.1e-4 |
| >=100 | 21% | 0% | 2.0e-5 |

Patient-level (deepest run each): 6/27 2BP vs 0/29 no-BP carry high-burden Pseudomonas (p=0.009).

Note the shape: at >=1 read the arms are identical (56% vs 54%) — every patient's airway has
*some* Pseudomonas-class signal. **Abundance** distinguishes the arms, not presence.

## Claim B — the matrix contrast is confounded

The 2BP matrix numbers stand (13/28 = 46% matrix-ON, 17/28 = 61% FAST). But **no no-BP run reaches
the matrix depth gate at all.** Attempting to break it by mapping the deepest no-BP runs at full
depth:

| run | patient | reads | Pseudomonas target | matrix reads |
|---|---|---|---|---|
| SRR26629355 | HS450 | 303k | 37 | 0 |
| SRR26629337 | HS158 | 434k | 0 | 0 |
| SRR26629340 | HS158 | 496k | 0 | 0 |

So "2BP 46% matrix-ON vs no-BP 0%" is near-tautological — you cannot detect a Pseudomonas operon
where there is no Pseudomonas. Report as **descriptive**, not as specificity.

## One-Two Punch impact

- The matrix-active, fast-growing populations are **present in the 2BP group**, and Pseudomonas
  burden does enrich there at matched depth.
- The two-axis **state** cannot be shown to distinguish 2BP from no-BP, because the comparator arm
  lacks the organism. The theory's *discriminative* claim is still untested by this cohort.

Consistent with — not a challenge to — Langelier et al.: their "increased bacterial RNA mass,
dominance of the culture-confirmed pathogen" is exactly what a high-burden Pseudomonas population
in the 2BP arm looks like through a species-specific lens.

## Caveats

1. **Pseudoreplication** — 5 patients contribute 4-8 runs each; run-level p-values overstate
   independence, patient-level numbers (6/27 vs 0/29) are the honest ones.
2. **Depth imbalance** — 2BP median run depth 1.03M vs no-BP 265k; matched screening controls
   screening depth but not global depth.
3. **TRA + Pseudomonas only.**
4. Growth calls depend on the artifact-corrected rp_over_ctrl.
5. **Nasal swabs null** (5/27 vs 7/45 at >=30, p=0.49) — enrichment is lower-airway, high-burden.

## Superseded in part

The "no comparator arm" conclusion in Claim B was an artifact of screening with a **PAO1-only**
index. See `RESULTS-ALLGENUS-2BP.md` for the pan-genus re-test, which finds the no-BP arm full of
matrix-capable organisms.
