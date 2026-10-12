# RESULTS — Host-phase test on UNSELECTED COVID cohorts (no 2BP split, no timing)

**Date:** 2026-10-12
**Repo:** joshux/acute-biofilm-analysis
**Question (Josh):** take COVID data that was *never* split into 2bp / no-bp, score it on the classifier **as a whole** (ignoring timing) — does it read bacterial?

**Answer: no. And the reason is now measured, not inferred.**

---

## What "unselected" means

"Unselected COVID" = a cohort that was never stratified by secondary-bacterial-pneumonia status. GSE246795 is the opposite — physicians pre-labelled every patient 2BP or no-BP, so that cohort is *defined* by the distinction. The three cohorts below have no such label: they are just "COVID patients," scored as one group.

## Three independent cohorts, COVID arm as a whole

| cohort | tissue | COVID n | median score | vs its own comparator |
|---|---|---|---|---|
| **GSE161731** (McClain/Duke) | whole blood | 77 | **−0.274** | bacterial pneumonia +1.522, influenza −0.090, healthy −0.410 |
| **GSE157103** (Overmyer/Albany) | whole blood | 97 | **+0.001** | non-COVID ICU +0.563 (COVID *lower*) |
| **GSE243217** (Muratsu) | PBMC | 35 | **+0.175** | bacterial sepsis +0.389, healthy −0.970 |

**COVID does not read bacterial in any of them.** It sits with influenza, other coronaviruses and healthy controls — far below genuine bacterial infection. Using the bacterial-pneumonia arm as the reference for what "bacterial" actually looks like:

- **Only 3% (2/77) of COVID patients in GSE161731** score at or above the bacterial arm's 25th percentile.
- **Only 2% (2/97) of COVID patients in GSE157103** do the same.
- In GSE243217, bacterial sepsis scores *higher* than COVID (+0.389 vs +0.175), and COVID is nowhere near healthy's separation (COVID vs healthy p < 10⁻⁴).

## The decisive result: the score tracks SEVERITY, not bacteria

GSE157103 is the only one of the three with clinical severity data. Within the COVID arm, the score rises with **how sick the patient is**:

| clinical variable | association with "bacterial" score | n |
|---|---|---|
| **Mechanical ventilation** | **+0.302 (ventilated) vs −0.378 (not), p = 0.0001** | 41 vs 56 |
| **ICU admission** | **+0.162 (ICU) vs −0.335 (ward), p = 0.013** | 49 vs 48 |
| D-dimer (coagulopathy severity) | Spearman ρ = **+0.32** | 82 |
| APACHE II (physiologic severity) | ρ = +0.20 | 56 |
| CRP | ρ = +0.15 | 89 |
| **Procalcitonin** (bacterial marker) | ρ = **+0.12** | 84 |
| SOFA | ρ = +0.09 | 55 |
| **Ferritin** | ρ = **+0.01** | 92 |

This is the confound made explicit: the "bacterial program" score is **most strongly predicted by whether the patient needed a ventilator, and least predicted by the actual markers of bacterial infection** (procalcitonin, CRP). In unselected COVID it is reading inflammation and illness severity, not bacteria.

That also explains the one counter-intuitive line above — non-COVID ICU patients scoring higher than COVID. Those are non-COVID ARDS/sepsis patients, i.e. the sickest group in the cohort.

## Compartment caveat

GSE243217 is **PBMC**, which is neutrophil-depleted; the 83 "bacterial" genes are neutrophil-degranulation genes, so this cohort is a *weak* test of the bacterial arm. It is included for completeness; the two whole-blood cohorts carry the result.

## Conclusion

Combined with the corrected GSE246795 result (`RESULTS-HOST-PHASE.md`), the host-side position is now fully characterised:

1. The signature is a **valid** bacterial-vs-viral discriminator (verified in the derivation cohort, p = 3×10⁻⁶; transfers to adult blood, p = 10⁻¹²).
2. Correctly applied to **2BP-labelled** COVID airway samples it separates the groups (p = 0.0014) — but that is near-definitional and severity-confounded.
3. In **unselected** COVID — blood or airway, whole or as a group — it reads **viral/neutral**, with the bacterial program *falling* over disease course, and the residual signal tracking ventilation and ICU rather than bacterial markers.

**The viral→bacterial host-program crossover is not present in any public unselected COVID cohort tested.** Three cohorts, three compartments/designs, same answer.

## Reproduction

`host-phase/unselected_test.py` (GSE161731 anchors + GSE157103), `host-phase/g243217_score.py` (sepsis arm).
Data: GSE161731 counts+key; GSE157103 genes.ec.tsv + series matrix; GSE243217 per-sample raw counts + series matrix.
