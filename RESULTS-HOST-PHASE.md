# RESULTS — Host-side viral→bacterial phase test (CORRECTED)

**Date:** 2026-10-12 (corrected; supersedes the uncommitted `RESULTS-HOST-PHASE.md` of the same date)
**Repo:** joshux/acute-biofilm-analysis
**Status:** the previously reported **null (p = 0.87) was a bug**. Corrected result: the signature **does** separate 2BP from no-BP in GSE246795 (p = 0.0014). But the phase claim Josh predicted does **not** appear in unselected adult COVID.

---

## 1. What the earlier run reported, and why it was wrong

The prior run reported: 2BP +0.243 vs No-BP +0.192, MWU **p = 0.87** → "the Mick classifier does not show the COVID viral→bacterial progression."

**Root cause: a sample-label misalignment.** `GSE246795_COVID_2BP_host_gene_counts.csv.gz` orders its 56 sample columns **differently** from the series-matrix sample order — only **34 of 56** labels agree between the two orderings. Joining the counts to the `2bp_status` labels **positionally** (rather than by sample name) scrambles the labels.

| join | 2BP | No-BP | p |
|---|---|---|---|
| **by sample name (correct)** | **+0.396** | **−0.369** | **0.0014** |
| positionally (the bug) | +0.010 | −0.009 | **0.8504** |
| *committed (uncommitted) run* | *+0.243* | *+0.192* | *0.87* |

The positional-join reproduction (p = 0.8504) is a near-exact match to the committed p = 0.87 — that is the evidence the earlier number came from the misalignment. (The original script was never committed, so it cannot be inspected directly.)

**Lesson for the project: always join GEO supplementary matrices to labels by sample name, never by column position.** This class of bug produces a plausible-looking null and would have propagated into the manuscript.

## 2. Signature recovered exactly

Source: JCI 10.1172/JCI165904, Supplemental Data File 2 (`sdd4.xlsx`, sheet `AnyBacterial_vs_PurelyViral_DE`) — the Mick/Tsitsiklis **108 genes** at BH-adjusted p < 0.05:

- **83 bacterial-high** (logFC < 0): S100A8, S100A9, S100A12, CD177, MMP9, PGLYRP1, TNF, CCL4, ORM1, …
- **25 viral-high** (logFC > 0): ALOX15, ERRFI1, MSLN, TM4SF1, SCNN1A, …

Score = mean z(log₂ CPM) over the 83 minus mean z(log₂ CPM) over the 25.

## 3. Positive control passes (pipeline is sound)

In the **pediatric derivation cohort** (eranmick/pediatric-mNGS-LRTI-classifier):

| contrast | n | diff | p |
|---|---|---|---|
| Definite vs No Evidence | 117 vs 50 | +0.803 | **3.0 × 10⁻⁶** |
| bacterial vs purely-viral | 7 vs 36 | +0.960 | **0.035** |
| coinfection vs purely-viral | 71 vs 36 | +1.278 | **9.8 × 10⁻¹⁰** |

The signature behaves exactly as designed where it was derived.

## 4. It **does** transfer to adult blood

**GSE161731** (McClain/Duke, adult whole blood; 198 samples joined):

| cohort | n | median score |
|---|---|---|
| **Bacterial pneumonia** | 24 | **+1.522** |
| Influenza | 17 | −0.090 |
| Other coronavirus | 61 | −0.143 |
| COVID-19 | 77 | −0.274 |
| Healthy | 19 | −0.410 |

Bacterial vs COVID-19 **p = 1.4 × 10⁻¹²**; bacterial vs healthy **p = 2.5 × 10⁻⁸**. The pediatric signature **is** a real adult bacterial-vs-viral discriminator — that is the transferability result the project was missing.

## 5. But Josh's prediction does **not** hold in unselected COVID

> *"find a data set that didn't differentiate 2bp vs no-bp… I'm betting many would show actually bacterial instead of viral."*

**Not confirmed.** In GSE161731 — a cohort that never applied a 2BP label — COVID-19 does **not** score bacterial. It scores viral/neutral (median −0.274), statistically indistinguishable from influenza and other coronaviruses, and **COVID-19 vs healthy p = 0.43** (this cohort is largely non-hospitalized). The temporal trend runs the **wrong way**:

| COVID by time since onset | early | middle | late |
|---|---|---|---|
| overall score | −0.063 (n19) | −0.280 (n36) | −0.306 (n22) |
| bacterial program | −0.031 | −0.130 | −0.195 |
| viral program | +0.032 | +0.150 | +0.112 |

The bacterial program **falls** over disease course and the viral program **rises** — the opposite of a viral→bacterial crossover. (Hospitalized-only: early +0.259 → middle −0.111, n = 6 each.)

**Caveats on this negative:** "late" here is >21 days (convalescent), not the day 5–10 transition window; the cohort is mild, so the relevant severe window is under-sampled; and only 6 hospitalized patients have serial timepoints.

## 6. The severity confound is not resolved

The 83 "bacterial" genes are **neutrophil-degranulation** genes. The corrected 2BP > no-BP signal is driven by exactly those genes (S100A12 +1.51, CST7 +1.34, S100A8 +1.00, CD177 +0.95) with the viral program down (ALOX15 −2.06). Because 2BP patients are sicker, this signal is **not separable from severity** in GSE246795 — whose GEO metadata carries **no severity variable** at all (only age, sex, race, 2bp_status), so adjustment is impossible in it.

In GSE161731 severity contributes but cannot explain the effect: hospitalized COVID scores only +0.015 vs −0.285 for non-hospitalized (p = 0.092), against a bacterial-pneumonia arm at +1.52.

## 7. Corrected position for the manuscript

The old rule ("do **not** cite Mick as support; the test is null, p = 0.87") is now wrong and must be replaced with:

1. **Validated, not null.** Correctly applied, the 108-gene signature separates adjudicated-2BP from no-BP COVID lower-airway samples (p = 0.0014) and transfers to adult blood (bacterial pneumonia vs COVID p = 1.4 × 10⁻¹²).
2. **But it is a discriminator for a *labelled* bacterial phenotype, not evidence of progression.** In 2BP-labelled data the separation is near-definitional and severity-confounded; in unselected adult COVID the same signature reads viral/neutral, with the bacterial program *falling* over time.
3. **Net:** the host-side viral→bacterial crossover remains **undemonstrated** — but now for a precise reason (no severe, serially-sampled, steroid-annotated adult cohort in the public domain), not because the instrument failed.

## 8. Reproduction

Scripts (this analysis): `mick108_raw.json`, `host_phase.py`, `score_variants.py`, `positive_control.py`, `positive_exact.py`, `misalignment_test.py`, `final_summary.py`.

Data: GSE246795 counts + series matrix; GSE161731 counts + key; eranmick/pediatric-mNGS-LRTI-classifier host counts + metadata; JCI supplement sdd4.xlsx.

**Method rule added:** join GEO supplementary matrices to phenotype labels **by sample name**, never by column position.
