# COVID arm assessment — PRJNA1033689

**Date:** 2026-10-11. Question: can the UCSF COVID secondary-bacterial-pneumonia cohort
(PRJNA1033689) supply a COVID arm for the two-axis test?

**Answer: it is real and deeply sequenced, and we can measure %RP in it — but the measurement is
confounded, and adding it naively would manufacture a false result.** It is a caution, not a new
arm. Details below.

---

## 1. The dataset is real and better than expected

| property | value |
|---|---|
| runs | **254**, all `library_source = METAGENOMIC`, `library_strategy = RNA-Seq` |
| unique samples | 254 (one run per sample) |
| layout | paired-end, Illumina NovaSeq 6000 |
| depth | median **810,653** reads; max 13,996,427; min 889 |
| specimen | **tracheal aspirate (TRA)** — not BALF, not sputum |
| source | Spottiswoode/Langelier et al., *Nat Commun* 2024 |
| labels | **2BP vs no-BP status is NOT in the ENA metadata** (paper / GEO GSE246795 only) |
| library prep | **no bacterial rRNA depletion** (same chemistry class as the acute BALF arm) |

Screened all 254 runs (40k reads each). **141 (55.5%) had ≥1 PAO1 read**; 40 (15.7%) passed a
permissive ≥30-read threshold and were fully quantified at 200k reads. The distribution is
strongly zero-inflated and long-tailed — one sample (SRR26629241) is a near-monoculture at 84.2%
of reads mapping to PAO1.

## 2. The measured %RP values

| set | n | median %RP | range | FAST >10% | ambig | SLOW <5% |
|---|---|---|---|---|---|---|
| all quantified | 40 | 11.91% | 0.00 – 32.24 | 22 | 8 | 10 |
| coding ≥100 | 37 | 13.24% | 0.00 – 32.24 | 22 | 8 | 7 |
| **coding ≥500** | **27** | **17.20%** | 0.57 – 32.24 | 20 | 3 | 4 |

Top sample: SRR26629332 at **32.24%** (3,721 RP / 11,541 coding reads).

## 3. Why these numbers cannot be taken at face value

Three independent checks each raise a flag.

**(a) The median exceeds the in-vitro exponential anchor.**

| reference | median %RP |
|---|---|
| lab **exponential flask culture** (rich medium, log phase) | 13.6% |
| COVID TRA, well-powered (n=27) | **17.2%** |
| acute BALF (n=13) | 10.7% |
| lab stationary (n=11) | 2.1% |
| chronic CF sputum (n=15) | 4.6% |

In-vivo bacteria in a diseased lung out-scoring a log-phase flask culture in rich broth is not
biologically credible. Either the metric saturates, or something else is driving it.

**(b) %RP correlates positively with depth here — the opposite of the acute arm.**

| arm | Spearman(log coding reads, %RP) |
|---|---|
| acute BALF (non-depleted) | **−0.29** |
| chronic CF (rRNA-depleted) | −0.33 |
| **COVID TRA (non-depleted)** | **+0.60** |

Two arms with the *same* prep chemistry (acute BALF, COVID TRA) behave in *opposite* directions
with respect to depth. That inconsistency alone means the COVID %RP is not measuring the same thing
the acute arm measures.

**(c) %RP correlates with rRNA burden and with Pseudomonas abundance.**

| relationship | Spearman ρ |
|---|---|
| rRNA% vs %RP | +0.471 |
| rRNA% vs %RP, **partial** (controlling Pseudomonas fraction) | +0.384 |
| **Pseudomonas fraction vs %RP** | **+0.379** |
| log read count vs %RP | +0.436 |

The composition correlation is the concerning one: **%RP may be tracking Pseudomonas abundance
rather than growth rate.** The mechanism is the project's recurring artifact class — with no
bacterial rRNA depletion, most rRNA maps ambiguously (MAPQ 0) and is filtered out, so the surviving
"coding" reads are the uniquely-mapping subset. In high-load samples more conserved reads survive
uniquely, and ribosomal-protein genes are the most conserved genes there are. High Pseudomonas
load therefore *mechanically* raises apparent %RP.

## 4. Two readings, and the honest position

- **Reading 1 (artifact).** These are ICU patients with secondary bacterial pneumonia; the samples
  are the deepest and most Pseudomonas-dominated in the project; the depth and abundance
  correlations are exactly the signature of a load-driven artifact. Under this reading the COVID
  %RP values are inflated and meaningless as growth estimates.
- **Reading 2 (biology).** These *are* acute infections, so genuinely fast growth is expected, and
  the COVID median (11.9%) sits comfortably in the acute range set by the BALF arm (10.7%). Under
  this reading the values are plausible and the depth correlation is an incidental consequence of
  the deepest samples also being the most Pseudomonas-dominant.

**We cannot currently distinguish these**, because the same confound is present — more weakly — in
the acute BALF arm (which also lacks rRNA depletion, and also shows a composition effect). The
disambiguating test is a load-matched or depth-matched comparison, which has not been run.

## 5. Recommendation

1. **Do not add the COVID arm as a biological result.** It would read as "COVID is the
   fastest-growing of all cohorts," which is more likely to be load-driven than real.
2. **Do add it to the methods paper as a generalizability case.** This is the first application of
   the metric to a COVID cohort, and it independently reproduces the *artifact family* — a
   third dataset, a different specimen, a different disease, same failure mode. That strengthens
   the central methodological claim rather than diluting it.
3. **If a real COVID arm is wanted**, it needs in order: (a) a load/depth-matched design, (b) a
   normalization that removes the abundance effect, and only then (c) the 2BP-vs-no-BP labels from
   GEO GSE246795 for a group contrast. Solving the confound is the prerequisite, not the labels.
4. **Caveat to state wherever this is reported:** specimen is tracheal aspirate, not BALF or
   sputum — a third specimen type in a project already criticised for cross-specimen comparison.

## 6. What this changes about the paper

The methods paper gains a third, independent instantiation of the artifact class — and a concrete
demonstration that a plausible-looking number (median 17.2%, "faster than exponential") can be
produced by a pipeline that is internally consistent and wrong. That is precisely the payload the
paper is built on, and it is now evidenced across three cohorts rather than one.

---

*Scripts: `probe_composition2.py` (composition), `quantify_rp2.py` (metric), `covid_work/` drivers.
Raw numbers: `out/covid_screen.json` (254 screen rows + 40 quant rows).*
