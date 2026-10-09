# Pilot result — feasibility of the Kolpen matrix+growth test on PRJNA1056765

**Date:** 2026-10-09 · **Status:** pilot complete, **GO with major design corrections** ·
**Samples tested:** SRR27343983 (bacterial, top-biomass), SRR34061072 (NC), plus composition
profiling of all 114 bacterial + 32 NC RNA runs.

## Headline

Access, metadata, and sequencing depth all work. But **two findings force a redesign before
the full run**: (1) the bacterial community in these samples is dominated by
*hospital/environmental opportunists* — chiefly **Acinetobacter** (top genus in 38/114 samples)
plus a strong reagent-contaminant signature — not the classic CAP pathogens the panel was built
around; and (2) mapping yield to a 7-organism reference panel is **low (~2.6% of reads)**, and
the dominant genus does not map well to its own type strain. The pilot also kills the NC
comparison arm as a bacterial-signal contrast.

---

## 1. What was verified (all green)

| check | result |
|---|---|
| ENA HTTPS download of host-depleted FASTQ | ✅ works (85 MB and 4.6 MB files pulled) |
| Read length | **50 bp single-end** — confirms *Sci Data*; the companion paper's "75-cycle" is wrong |
| Metadata completeness | ✅ 434 RNA rows, 5 groups, per-run `unhost_reads` / `unhost_unrRNA_reads` / `micro_reads` |
| Group counts | 114 bacterial, 32 NC (+79 fungal, 86 TB, 123 cancer) |
| Paired DNA | ✅ every RNA run has a DNA run (pair NC by Patients ID) |
| Depth (bacterial group, RNA) | median **46,202** micro_reads; max 5.65M; **95/114 ≥10k**, 33 ≥100k, 10 ≥1M |
| Depth (top-biomass samples) | 0.9–11.2M microbial reads — ample for a metatranscriptomic panel |

Deposited reads are the **host-depleted fraction** (`unhost_reads`), as expected. For the test
sample: 20.5M total → 3.36M unhost → 1.48M Kraken-classified microbial.

## 2. Finding 1 — the community is not what the panel assumed

Composition of the 114 bacterial-infection RNA samples (genus relative abundance, Figshare matrix):

| genus | mean relab | median relab |
|---|---|---|
| Pseudomonas | 0.167 | 0.073 |
| **Acinetobacter** | **0.151** | 0.076 |
| Delftia | 0.061 | 0.046 |
| Haemophilus | 0.056 | 0.0004 |
| Comamonas | 0.038 | 0.007 |
| Sphingomonas | 0.021 | 0.011 |
| Ralstonia | 0.021 | 0.003 |
| Stenotrophomonas | 0.019 | 0.011 |
| Klebsiella | 0.020 | 0.005 |
| Staphylococcus | 0.016 | 0.004 |
| Streptococcus | 0.005 | 0.000 |

**Top genus per sample:** Acinetobacter 38, Pseudomonas 32, Haemophilus 10, Porphyromonas 5,
Mycoplasma 4, Treponema 3, … — i.e. the cohort is dominated by **non-fermenting Gram-negative
hospital opportunists**, consistent with ventilated/ICU patients, not community-acquired CAP.

**Reagent/environmental contaminant signature:** Comamonas + Delftia + Ralstonia + Sphingomonas
+ Rhodococcus + Paraburkholderia. Median contaminant-genera fraction **0.179**; **31/114 samples
>0.30**. This is the well-known BALF mNGS "kitome" pattern and must be handled explicitly
(contaminant-aware analysis, or restriction to samples where it is low).

**Consequence:** the panel's original organisms (S. pneumoniae, S. aureus, H. influenzae) are
minor here. A panel must be rebuilt around **Acinetobacter, Pseudomonas, Stenotrophomonas,
Klebsiella** — which is a real problem, because **Acinetobacter has almost no characterized
biofilm-matrix operon** (no pel/psl/ica/csg equivalent; only bap, PNAG-like, capsule). The
matrix module is thinnest exactly where the biomass is highest.

## 3. Finding 2 — mapping yield is low

Mapped the top-biomass bacterial sample (3.36M reads) against a 7-genome panel
(PAO1, NCTC 8325, TIGR4, MG1655, MGH78578, 86-028NP, A. baumannii ATCC 17978), minimap2 `-x sr`,
MAPQ ≥ 20:

| panel | reads mapped (all) | MAPQ ≥ 20 |
|---|---|---|
| 6 genomes (no Acinetobacter) | 81,715 (2.43%) | 2,834 (0.084%) |
| **7 genomes (+ A. baumannii)** | 88,070 (2.62%) | **3,614 (0.108%)** |

Adding the dominant genus lifted MAPQ≥20 mapping by 28% — so the panel *is* the binding
constraint. But even so, only ~2.6% of reads map to any panel genome, and only ~6% of
Kraken-classified microbial reads (88k / 1.48M).

**The deeper problem:** Acinetobacter is 29% of relative abundance in this sample, yet only
**766 reads** map to *A. baumannii* ATCC 17978. The sample's Acinetobacter is evidently not that
strain/species (A. nosocomialis / A. pittii / A. johnsonii are common in this niche). **A
single reference strain per genus is insufficient** — multi-strain, multi-species references are
required, or reads must be assigned at genus level with relaxed identity thresholds.

## 4. Finding 3 — the NC arm has no bacterial signal

The NC sample (SRR34061072) yielded **2 reads** mapping to the panel. NC samples are defined by
negative BALF microbiology and negative mNGS, so they carry near-zero bacterial RNA. This means
**the bacterial-vs-NC contrast is not the right design** — there is nothing to compare against.
A viable control arm would need a different comparator (e.g. bacterial vs fungal/TB infection
samples, which have comparable bacterial loads but different clinical context), or a within-
bacterial analysis using severity/depth as the axis.

## 5. Verdict — GO, with corrections

The analysis is feasible on the data (depth is there), but the **design must change**:

1. **Rebuild the reference panel** around the actual dominant organisms: Acinetobacter (multiple
   species), Pseudomonas, Stenotrophomonas, Klebsiella, plus retain the classic pathogens.
   Use **multiple reference strains per genus**.
2. **Accept that the matrix module is thin for Acinetobacter** — its matrix biology is
   poorly characterized. Either (a) restrict the matrix analysis to samples dominated by
   organisms with a defined matrix operon (Pseudomonas, Klebsiella/Enterobacteriaceae), or
   (b) reframe the matrix module as a general EPS/capsule/adhesin signature rather than
   named operons.
3. **Drop the bacterial-vs-NC contrast.** Use a within-bacterial design (e.g. Pseudomonas-
   dominant vs Acinetobacter-dominant samples; or severity) or bacterial vs fungal/TB.
4. **Handle contamination explicitly** — flag the 31/114 high-contaminant samples, and
   consider a contaminant-aware sensitivity analysis.
5. **Mapping strategy:** genus-level assignment with multi-strain references; report the
   unassigned fraction as a QC metric (if it is ~94%, the panel is under-specified — as here).

## 6. What this pilot does NOT say

It does not falsify Kolpen. It shows the *instruments* need recalibration to the actual
community before the biology question can be asked. The depth is sufficient; the panel and the
comparison design were wrong.

## Files

- `pilot_samples.csv` — the 14 pre-selected samples (10 bacterial + 4 NC)
- composition profiling across all 114 bacterial + 32 NC RNA runs (script + output to be committed)
- panel genomes + mapping counts for SRR27343983 (to be committed)
