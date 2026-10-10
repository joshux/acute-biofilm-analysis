# Writing plan — the scoped two-axis paper

**Date:** 2026-10-10
**Decision:** publishable as a **methods-forward short report**. Target *npj Biofilms and
Microbiomes* (brief communication format); fallbacks *mSystems*, *Microbiome*, *BMC Microbiology*.
Preprint to bioRxiv first.

**What this paper is:** a Pseudomonas-centred, two-axis (growth × matrix) metatranscriptomic test
of Kolpen et al. 2022's acute-biofilm model, with a validated growth metric, an independent-proxy
robustness check, and two documented artifact classes as methods contributions.

**What it is not:** a definitive "acute grows faster than chronic" comparative claim, and not a
COVID paper.

---

## Evidence stack (what the manuscript can stand on)

| # | unit | status |
|---|---|---|
| 1 | %RP growth metric validated on internal anchors (lab exp 13.6% vs stat 2.1%, p<1e-4) | committed (`RESULTS-QUADRANT-CORRECTED.md`) |
| 2 | Two-axis result: acute matrix-ON×FAST/ambig vs chronic matrix-ON×SLOW (10.7% vs 4.6%, p=0.0015) | committed (`RESULTS-GENERALISED.md`) |
| 3 | MAPQ multi-genome artifact (multi-genome panels drop ~94% of conserved RP reads) | committed (`RESULTS-QUADRANT-CORRECTED.md`) |
| 4 | Composition-table disagreement (deposited tables vs RNA reads) | committed (earlier) |
| 5 | **%TF independent-proxy concordance** (29 non-RP genes; ρ=0.91; acute vs chronic p=0.0008; depth-restricted p=0.0118) | this commit (`TF-ROBUSTNESS.md`) |
| 6 | Depth asymmetry disclosed + depth-restricted contrast survives (p=0.0142 on %RP) | committed (`RESULTS-ACUTE-SPUTUM-SEARCH.md`) |
| 7 | Negative results: no matched acute sputum dataset exists (3 cohorts screened); Acinetobacter arm is non-baumannii (correct 0-matrix); pediatric replication dead on organism (80/80 runs non-Pseudomonas); healthy controls = kitome floor | committed (`RESULTS-ACUTE-SPUTUM-SEARCH.md`, `ORGANISM-VIABILITY.md`) |

Statistics to report (fixed now, before drafting): Mann-Whitney two-sided, acute vs chronic both
unrestricted and at the ≥10k coding-read floor; chronic arm additionally at patient level
(one sample per patient, p≈0.046–0.0495 — the arm is 5 patients); %TF as pre-specified
independent proxy; Spearman depth-correlation per arm. The acute arm has no repeated patients
(16 runs → 16 BioSamples, verified in `TF-ROBUSTNESS.md` Task 0).

## Constraints carried from earlier threads

- The chronic arm is 5 patients (P30M0 contributes 6) — pseudoreplication-corrected p≈0.05. State
  it in Results, not just Limitations.
- Cross-specimen contrast (acute BALF vs chronic sputum), bounded by data availability
  (documented in `RESULTS-ACUTE-SPUTUM-SEARCH.md`).
- Kolpen 2022 was sputum+FISH only; we are the first transcript-level test — say so in Intro,
  cite Thorax 77(10):1015, DOI 10.1136/thoraxjnl-2021-217576.
- The Oct 9 `PAPER-DRAFT.md` argues the SUPERSEDED conclusion ("matrix machinery is transcribed
  but not co-expressed with growth"; "the acute-biofilm prediction is not supported") — written
  before the MAPQ correction reversed the result. **Do not edit it; supersede it.**

## Stages

**A. Foundations (done or <half day)**
- [x] A1 %TF independent-proxy robustness (`TF-ROBUSTNESS.md`) — concordant
- [ ] A2 Archive old draft: `git mv PAPER-DRAFT.md PAPER-DRAFT-SUPERSEDED-20261009.md`
- [ ] A3 bioRxiv preprint decision + AI-use disclosure statement (journal-template wording;
      verify current npj/Springer Nature policy at submission) + ORCID on record
- [ ] A4 Citation verification pass — every DOI/PMID checked against the source; no LLM-supplied
      citation enters the reference list unverified

**B. Skeleton (half day)**
- [ ] B1 Title candidates (3): lead with the validated metric + the two-axis confirmation
- [ ] B2 Claims map: 5 defensible claims → evidence units 1–7 (one-to-one; anything unmapped
      gets cut)
- [ ] B3 Objection log: 10 hostile-reviewer objections with responses (depth asymmetry, 5-patient
      chronic arm, cross-specimen, Ovation vs RiboZero chemistry mismatch, MAPQ artifact
      generalizability, decoy-index sensitivity, Pseudomonas-only generalization, p-value at
      patient level, "why not sputum", novelty vs Kolpen)
- [ ] B4 Figure plan: Fig 1 = two-axis quadrant (exists, may need %TF inset); Fig 2 = metric
      validation + artifact (wrong vs corrected side by side); Table 1 = cohort composition;
      Table 2 = %TF robustness. Keep it to 2 figures + 2 tables (brief-report format).

**C. Drafting (2–3 days, Methods first)**
- [ ] C1 Methods (full reproducibility from repo: panels, decoy index, mapping recipe, both
      quantifiers, statistics)
- [ ] C2 Results (2.1 cohort + organism viability → 2.2 matrix axis → 2.3 growth axis + anchors
      → 2.4 %TF robustness → 2.5 negative/limitation results)
- [ ] C3 Introduction (Kolpen's model; imaging-only gap; two-axis decoupling prediction)
- [ ] C4 Discussion (what confirms, what is bounded by data; artifact as methods contribution;
      specimen-availability landscape as a finding)
- [ ] C5 Abstract + cover letter

**D. Hardening (1 day)**
- [ ] D1 Hostile pre-read against the objection log (B3); fix or pre-empt each
- [ ] D2 Reproducibility spot-check: fresh clone, re-run one sample end-to-end from committed
      scripts
- [ ] D3 bioRxiv post (dated preprint = priority on the artifact finding)

**E. Submission**
- [ ] E1 *npj Biofilms and Microbiomes* brief communication; desk-reject fallback = *mSystems*
      short report, then *BMC Microbiology*

## Do not do

- No COVID extension inside this paper (that is the second paper, on PRJNA1033689, gated on this
  one being submitted).
- No "definitive biology" language anywhere — every biological sentence carries the scope clause.
- No new cohorts before drafting (the sputum search is closed; the answer was negative).
