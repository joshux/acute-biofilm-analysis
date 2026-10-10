# Search for an acute **sputum** arm — negative result, and a depth finding that matters

**Date:** 2026-10-10
**Objective:** remove the cross-specimen weakness in the two-axis test by finding an acute
**sputum** metatranscriptome — matching Kolpen et al. 2022 (*Thorax* 77:1015), who used
expectorated sputum, not BALF. The current acute arm is BALF (PRJNA1056765); the chronic arm is
CF sputum (PRJEB24688).

**Outcome: no such dataset exists publicly.** Two candidates were screened and both fail, for
independent reasons. The search is documented here so it is not repeated.

---

## The criterion that decides everything: bacterial rRNA depletion

%RP = (ribosomal-protein transcripts) / (bacterial coding reads). This is only interpretable if
the library was built to **deplete bacterial rRNA**.

**Poly-A selection is the wrong criterion.** Bacteria have no poly-A tail, so poly-A selection
discards exactly the transcripts being measured. What is required is bacterial 16S/23S depletion
(RiboZero, QIAseq FastSelect 5S/16S/23S, or equivalent).

The two existing arms were built with **different chemistries**:

| arm | dataset | prep | bacterial rRNA depleted? |
|---|---|---|---|
| acute BALF | PRJNA1056765 | Ovation Trio RNA-Seq (AnyDeplete) | **No** — AnyDeplete targets *human* rRNA |
| chronic CF sputum | PRJEB24688 | RiboZero Epidemiology Gold + ScriptSeq v2 | **Yes** (`library_selection = Inverse rRNA`) |

This asymmetry is the reason %RP is defined with rRNA/tRNA excluded from the denominator
(`quantify_rp2.py`) — that makes the metric independent of the depletion difference, and it is why
the lab anchors validate despite the mixed prep. But it also means any new acute arm should
ideally be **RiboZero-family**, and it raises the bar on what counts as a match.

---

## Candidate 1 — PRJNA901252 (adult CAP + SCAP, sputum, China)

**Rejected on two independent grounds.**

**(a) Organism.** The Pseudomonas-centred matrix axis does not transfer.

- The paper (Han et al., Research Square preprint 2022, `10.21203/rs.3.rs-2182064/v1`; **no PMID,
  no peer-reviewed version found**) reports ***Acinetobacter baumannii*** as the dominant
  SCAP/non-survivor pathogen. *P. aeruginosa* is **not** a significant finding; the *Pseudomonas*
  **genus** is actually enriched in the *milder* CAP group, and the only species named is
  *P. fluorescens*.
- This is consistent with the organism-viability audit: the *Acinetobacter* in these cohorts are
  non-baumannii species that genuinely lack `pga`/`csu`, so a matrix axis cannot be built for them
  without a new panel.

**Empirical probe** (3 sputum metatranscriptomes: SRR22371550 9.7M reads, SRR22371655 3.3M,
SRR22371494 0.5M; 100k-read subsamples, `probe_composition2.py`):

| run | reads | PAO1 reads (any MAPQ) | PAO1 **coding** reads | decoy (non-Pseudomonas) | % reads non-Pseudomonas |
|---|---|---|---|---|---|
| SRR22371550 | 712,077 | 1,333 | **10** | 119,856 | 16.8% |
| SRR22371655 | 918,539 | 511 | **496** | 35,822 | 3.9% |
| SRR22371494 | 546,653 | 142 | **63** | 45,526 | 8.3% |

Pseudomonas is present at **0.001–0.054%** of reads while other bacteria sit at **4–17%**. There is
no Pseudomonas signal to measure. (These runs are already host-removed — the submitted filenames
end `.unHost.clean.unPlasmid.fq.gz`.)

**(b) Library prep.** ENA metadata gives `library_selection = RANDOM PCR` for the RNA runs. There
is **no documented bacterial rRNA depletion, no poly-A selection, and no wet-lab host depletion** —
host reads were removed *computationally*, after sequencing. The preprint contains **no
library-preparation methods section at all**. Under (a) this is moot, but it independently
disqualifies the dataset for %RP.

*Note on run counts:* an early ad-hoc `awk` split reported "CAP 77 / SCAP 0". That was wrong — it
matched the substring `CAP` inside `SCAP`. The correct ENA counts are **CAP sputum 30 / SCAP
sputum 51** metatranscriptomic runs.

---

## Candidate 2 — PRJEB60515 (acute pneumonia, sputum + bronchial wash, Sweden)

**Rejected on organism; retained as a methods control.**

This is the **only dataset found anywhere with documented bacterial-specific rRNA depletion**:
*"Bacterial ribosomal RNA (rRNA) was depleted from 100 ng of total RNA with QIAseq FastSelect
5S/16S/23S (Qiagen)"* — Polland et al., *Microbiol Spectr* 2023;11(5):e0163923,
PMID 37707456. ENA confirms `library_selection = Inverse rRNA` for all 15 metatranscriptomic runs,
at excellent depth (median ~49M reads, all >5M).

The problem is biological: the cohort was **selected for *H. influenzae***, so there is essentially
no Pseudomonas, *S. aureus*, *K. pneumoniae* or *A. baumannii*.

Empirical probe (ERR11532895, 100k reads): 55% of reads map to the panel, but the hits are almost
entirely decoy genomes (top: GCF_001457635.1 20,572 reads; GCF_000696675.2 9,627; GCF_019048645.1
5,481) and **PAO1 receives 0 reads**.

**Use it as a positive control for the prep chemistry** — a documented-depleted, deeply-sequenced
sputum metatranscriptome in which the %RP machinery can be shown to behave sensibly — not as a
study arm.

---

## Candidate 3 — PRJDB41581 (deep sputum, ICU pneumonia vs colonisation)

**Rejected on prep (pending nothing — it is the wrong chemistry).**

Attractive on design: 48 metatranscriptomic runs, median 18.6M reads, **sputum**, acute pneumonia
(n=27) vs **respiratory colonisation controls (n=17)** — the closest thing to a non-infected
control in an ICU.

But ENA's own `library_construction_protocol` states the RNA libraries used the **Ovation Trio
RNA-Seq Library Preparation Kit** — the same *host*-rRNA-depleting chemistry as the existing acute
BALF arm, not bacterial depletion. Probe (DRR1001378, 100k reads): only 1.1% of reads map at all
(1,116), of which 48 hit PAO1 and the rest scatter across decoys. Bacterial signal is too thin.

---

## The screening conclusion

> **There is exactly one public dataset with documented bacterial-rRNA depletion (PRJEB60515), and
> its pathogen composition is wrong for this analysis. There is exactly one dataset with the right
> specimen, disease, controls and depth (PRJDB41581), and its prep is host-rRNA depletion only.
> They are not the same dataset.**

A matched acute **sputum** arm — classic CAP pathogens, documented bacterial rRNA depletion,
adequate bacterial coding depth, raw FASTQ public — **does not currently exist in the public
archives.** Reaching a same-specimen acute-vs-chronic contrast would require generating the data.

---

## What the probes turned up that matters more than the search

Probing for a new arm exposed a real limitation in the **existing** result. It does not reverse the
finding, but it must be reported.

**The acute arm is ~92× shallower than the chronic arm at the level that the metric actually uses.**

| arm | n | median **coding** reads (the %RP denominator) | min | samples <10k coding reads |
|---|---|---|---|---|
| acute BALF | 15 | **10,295** | 38 | **7 / 15** |
| chronic in vivo | 15 | **951,671** | 216,438 | 0 / 15 |
| lab exponential | 12 | 1,119,028 | 531,838 | 0 / 12 |
| lab stationary | 11 | 1,068,273 | 613,910 | 0 / 11 |

Two reasons this is not fatal:

1. **Depth does not inflate the metric.** Spearman(log10 coding reads, %RP) is **negative** in both
   arms (acute −0.29, chronic −0.33). If shallow samples were artefactually high, the correlation
   would be positive.
2. **The contrast survives a depth floor.** Restricting the acute arm to samples with ≥10k coding
   reads (n=8) gives median **9.7%** vs chronic **4.6%**, **p=0.0142** — versus p=0.0015
   unrestricted. The result weakens but holds.

The honest statement is therefore: **acute BALF 10.7% vs chronic CF sputum 4.6% (p=0.0015; p=0.0142
restricted to acute samples with ≥10k coding reads; patient-level chronic p≈0.05).** The acute
median rests on a few thousand reads per sample, and that belongs in the limitations, not in a
footnote.

---

## Replication cohort: dead on organism, not on method

The pediatric tracheal-aspirate cohort (PRJNA875913) is **not** a usable replication arm. Of 80
scored runs, **80/80** have more decoy reads than Pseudomonas reads; median Pseudomonas coding
reads = **2**; **zero** runs reach 1,000. No label subgroup (Definite / Suspected / No Evidence)
has Pseudomonas signal.

This is consistent with everything above: the failure is **organism availability in public
respiratory metatranscriptomes**, not the method.

---

## Recommendation

1. **Do not adopt an acute sputum arm.** Neither candidate supports it.
2. **Keep the BALF acute arm**, and add the depth asymmetry above to the limitations — it is a
   genuine constraint on the headline, and reporting it pre-empts the reviewer who finds it.
3. **Use PRJEB60515 as a methods positive control** if a documented bacterial-depletion prep is
   wanted in the paper.
4. **Frame the specimen issue honestly:** the comparison is cross-specimen and cross-cohort by
   necessity. The internal lab anchors are what carry the metric's validity, and they are
   within-dataset.
5. **The two-axis result stands** — it is now bounded by *data availability*, which is a defensible
   and checkable statement, rather than by an unexamined assumption.

---

## Scripts

- `probe_composition.py` — read-composition probe with MAPQ≥20 filter.
- `probe_composition2.py` — full contingency table by MAPQ band (0 / 1–19 / 20+) × target class
  (PAO1 rRNA-tRNA / ribosomal protein / other CDS) vs decoy. Use this one: the MAPQ≥20-only view
  hides the composition, because most short-read bacterial alignments land at MAPQ 0 against a
  multi-genome decoy index.
