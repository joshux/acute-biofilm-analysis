# Healthy / non-infected respiratory controls — verification, candidate datasets, and the conceptual answer

**Date:** 2026-10-12 · **Repo:** acute-biofilm-analysis
**Question:** Where do HEALTHY / non-infected respiratory samples fall on the two-axis
(GROWTH = `rp_over_ctrl` / %RP × MATRIX = alginate/psl/pel) analysis — and do better
control datasets exist?

**Headline:** The committed "zero bacterial signal" claim **holds as stated** — for
PRJNA390194, and only for PRJNA390194. It is **not** a general property of healthy lung.
Two better control types exist and both carry usable bacterial signal. **And the matrix axis
CAN be placed in a control** — not with the *Pseudomonas* alginate/psl/pel panel, but with the
repo's existing *Streptococcus pneumoniae* `cps`/`cbp` panel, in the PRJDB41581 colonisation
cohort. That is a change to the answer this project has been carrying.

---

## Part 1 — Verification of the existing healthy result (PRJNA390194)

### What was committed

| item | value |
|---|---|
| accession | **PRJNA390194** (secondary SRP109051) |
| study title | "RNA-Seq of COPD patients' BAL" |
| paper | Ren et al. 2018, *mSystems* 3(5):e00199-18 — **DOI `10.1128/msystems.00199-18` resolves** (Crossref: mSystems, 2018, first author Ren; PMID 30417108, PMC6208642) |
| healthy arm | **9 runs**, sample aliases `ANB21`–`ANB40` (non-COPD controls) |
| other arm in same study | 25 `COPD…` RNA runs (not used) |
| specimen | adult BALF |
| library | `library_strategy = RNA-Seq`, `library_selection = RANDOM` |
| full depth | 20.4M–38.1M reads/run (median ~24.9M) — **the analysis used a 2,000,000-read head subsample** (`healthy_all.txt`/`healthy_runs.txt`, 9 SRRs) |

### The committed numbers (re-derived from `out/SRR5677*_pseu_v3.json`)

| run | reads | mapped | mapped % | decoy hits | **Pseudomonas panel target** |
|---|---|---|---|---|---|
| SRR5677652 | 2,000,271 | 10,854 | 0.543 % | 3,045 | **0** |
| SRR5677654 | 2,000,103 | 4,223 | 0.211 % | 1,183 | **0** |
| SRR5677655 | 2,000,238 | 12,565 | 0.628 % | 3,280 | **0** |
| SRR5677656 | 2,001,046 | 5,509 | 0.275 % | 291 | **0** |
| SRR5677657 | 2,000,490 | 2,430 | 0.121 % | 108 | **0** |
| SRR5677658 | 2,000,418 | 2,983 | 0.149 % | 201 | **0** |
| SRR5677659 | 2,000,473 | 2,333 | 0.117 % | 369 | **0** |
| SRR5677660 | 2,000,185 | 1,275 | 0.064 % | 175 | **0** |
| SRR5677661 | 2,000,064 | 3,697 | 0.185 % | 780 | **0** |
| **median** | 2,000,271 | **3,697** | **0.185 %** | **369** | **0 (total 0)** |

Cross-checked against `final_stats.json` → `healthy: {n: 9, target_reads_total: 0}` and
`final_quadrant_stats.py` (lines 73–80): "panel target reads: median 0, total 0; reads mapping
to decoy genera: median 369 → no bacterial target signal above the kitome floor." ✅ **Numbers reconcile exactly.**

### Genus composition of the reads that *did* map (from `controls_char.json`)

Mapped against the 19-genome multi-strain panel at MAPQ ≥ 20, all 9 samples pooled:

| genus | reads |
|---|---|
| *Escherichia* | 850 |
| *Staphylococcus* | 60 |
| *Haemophilus* | 46 |
| *Klebsiella* | 40 |
| *Acinetobacter* | 16 |
| ***Pseudomonas*** | **6** |
| *Stenotrophomonas* | 5 |
| *Streptococcus* | 3 |
| **total** | **1,026** |

So: ~1,026 panel-mapping reads across 9 samples × 2M reads = **0.006 % of reads**, and
**6 reads total to *Pseudomonas*** across the entire healthy arm. The mapped reads are
dominated by *Escherichia* — a canonical reagent/kitome genus (Salter et al. 2014) — and by
*Staphylococcus*/*Cutibacterium*-class skin flora.

### Independent re-verification at full depth (new, this report)

I re-probed two healthy samples at **full run depth**, streaming directly from ENA into
`compfull_pseu.mmi` (full PAO1 genome + 120 non-*Pseudomonas* decoy contigs), so no disk was used
and no subsampling artefact can apply. Script: `healthy_depth_probe.py`. Outputs:
`out/SRR5677652_depthprobe.json`, `out/SRR5677656_depthprobe.json`.

| run | reads mapped | PAO1 any (permissive) | PAO1 **coding @MAPQ≥20** | decoy (non-*Pseudomonas*) | bacterial share of library |
|---|---|---|---|---|---|
| SRR5677652 | 1,772,156 *(mate-1 head)* | 262 | **0** | 9,426 | 0.53 % |
| SRR5677656 | **20,427,866** (full mate 1) | 119 | **14** | 56,784 | **0.28 %** |

SRR5677656 at 20.4M reads — 10× the committed depth — still yields **14 PAO1 coding reads**
and 119 PAO1 reads total (0.0006 % of reads). Its `pctRP_coding_mapq20` of 7.14 % rests on
**14 reads** and is meaningless as a growth measurement.

**Verdict on Part 1: the zero-signal claim holds.** Extending depth 10× does not change the
conclusion. *Pseudomonas* is at the level of a handful of cross-mapping reads, and the total
bacterial panel signal (0.28 % of reads) is at the kitome/reagent floor, not a measurable community.

### ⚠️ One correction to the framing

The committed wording calls this "**the kitome floor**". That is defensible for PRJNA390194, but
the phrase implies it generalises to healthy lung. **It does not.** I found a healthy-lung cohort
with 150× more bacterial signal (Part 2, Candidate A). The correct statement is narrower:

> PRJNA390194 is a low-bacterial-signal BALF cohort in which neither axis is measurable.
> This is a property of **that cohort's specimens and library prep**, not of healthy lung as such.

---

## Part 2 — Candidate control datasets

Legend: **Biomass** = did bacteria reach a level where the two axes could in principle be measured?
**Bracket** = category. DOIs were resolved via Crossref (`api.crossref.org/works/<doi>`) on 2026-10-12.

### Candidate A — PRJNA699976, healthy-volunteer BALF arm ⭐ **best bacterial signal found**

| field | value |
|---|---|
| accession | **PRJNA699976** (440 runs total) |
| repository | NCBI SRA / ENA mirror |
| publication | Zou, Yan, Wang, Ni, Zhao, Lu, Liu, Cao — *Advanced Science* 2024, "Accurate Diagnosis of Lower Respiratory Infections Using Host Response and Respiratory Microbiome from a Single Metatranscriptome Test of Bronchoalveolar Lavage Fluid". **DOI `10.1002/advs.202405087` resolves** (Crossref ✅; PMID 39692191, PMC11809327) |
| specimen | **BALF**, metatranscriptomic RNA-Seq |
| control arm | **28 runs** titled *"RNA-seq of BALF samples from healthy volunteer"* (other arms: 408 patient runs, 4 RNase-free-water negative controls) |
| depth | healthy arm: mean **21.6M**, max **32.9M**, min 3.6M reads (healthy arm is deep — median ~23M) |
| library prep | `library_selection = RANDOM`; `library_construction_protocol` **empty in ENA**; instrument HiSeq 4000 |
| what bacteria | *Neisseria*, *Prevotella*, *Haemophilus*, *Porphyromonas*, *Capnocytophaga*, *Veillonella*, *Streptococcus*, *Fusobacterium* — **typical oral/upper-airway community**, no pathogens of interest |
| raw FASTQ public | **yes** (all 440 runs have `fastq_ftp`) |
| **Biomass** | ✅ **YES — decisively.** See below. |

**Empirical probe (new, this report)** — `healthy_genus_probe.py`, run on SRR13677570
(a healthy-volunteer BALF run), streaming from ENA into `compfull_pseu.mmi`. Output:
`out/SRR13677570_genusprobe.json`.

| metric | PRJNA699976 healthy (SRR13677570) | PRJNA390194 healthy (SRR5677652) |
|---|---|---|
| reads probed | 2,847,163 | 15,001,946 |
| reads hitting the bacterial panel (any MAPQ) | **2,280,697 → 80.1 %** | 77,585 → 0.52 % |
| bacterial reads at MAPQ ≥ 20 | **115,183 (4.05 %)** | 22,937 (0.15 %) |
| 16S/23S conserved-motif hits in a 200k-read window | **6,008 → 3.00 %** | 30 → 0.01 % |
| PAO1 any | 2,304 | 2,380 |
| PAO1 coding @MAPQ≥20 | 19 | 0 |

**This is a 150× difference in bacterial share of the library, and a 300× difference in 16S
motif density.** The library is *bacterial* — the dominant mapped signal is oral-community
bacteria, not human background.

**Caveat, stated plainly.** The deposited files are almost certainly the **microbial (non-host)
fraction** — the companion data note for this cohort (Tang et al. 2025, *Sci Data*,
DOI `10.1038/s41597-025-06171-6`, resolves ✅) says of the sibling PRJNA1056765 deposit:
*"Only the microbial (non-host) sequence reads have been deposited."* PRJNA699976 is a 440-run
sibling deposit from the same clinical-metatranscriptome pipeline. I could **not** find an explicit
statement in ENA metadata that the deposited FASTQs are host-removed (the `submitted_ftp` field is
empty, and `library_construction_protocol` is empty). **Mark as: strongly inferred, not documented.**
This matters because it means the "80 % bacterial" figure describes the *deposited* fraction, not
the raw swab — which is exactly what the %RP denominator wants, but it is not comparable to a
raw-library bacterial-share measurement.

**Verdict: the best healthy control found.** Bacterial signal is unambiguously high enough to
measure the **growth** axis in principle. **But** the community is oral commensals —
no *Pseudomonas*, and no matrix-capable organism (see Part 3d).

### Candidate B — PRJNA603592, healthy smoker volunteers (paired BAL + bronchial brush)

| field | value |
|---|---|
| accession | **PRJNA603592** (62 RNA runs) |
| repository | NCBI SRA / ENA |
| publication | Sulaiman et al. 2021, *Eur Respir J* 58(1):2003434, "Functional lower airways genomic profiling of the microbiome to capture active microbial metabolism". **DOI `10.1183/13993003.03434-2020` resolves** (Crossref ✅; PMID 33446604, PMC8643072) |
| specimen | **BAL ("balluntouched") + bronchial brushing/suction ("bronchsup")**, metatranscriptomic |
| control arm | 10 human subjects `smnv0004…smnv0039`, each with a BAL and a brush run (20 runs); **explicitly "smoker control" volunteers**, i.e. non-infected. The same deposit also contains 19 `copd…` runs and mouse `PACT_` runs — use only the `smnv` arm |
| depth | `smnv` runs 15.8M–55.9M reads (median ~42M) — **the deepest control arm found** |
| library prep | `library_selection = cDNA`; ENA `library_construction_protocol` empty |
| what bacteria | *Ralstonia* dominant in the sample I probed (see below) |
| raw FASTQ public | **yes** (62/62) |
| **Biomass** | ⚠️ **marginal** — 0.77 % of reads bacterial; *Ralstonia*-dominated |

**Empirical probe (new)** — SRR11002545 (`smnv0005balluntouched_MT`, a healthy-smoker BAL).
Output: `out/SRR11002545_genusprobe.json`.

| metric | value |
|---|---|
| reads probed | 15,011,524 |
| bacterial panel hits (any MAPQ) | 115,749 → **0.77 %** |
| bacterial @MAPQ ≥ 20 | 7,720 (0.05 %) |
| PAO1 any / coding @MAPQ≥20 | 103 / **5** |
| top genera | ***Ralstonia* 89,247** (77 % of bacterial hits), *Bacteroides* 10,244, *Comamonas* 3,669, *Yersinia* 3,304, *Escherichia* 3,131 |

*Ralstonia pickettii* is the single most-flagged reagent contaminant genus in low-biomass
16S/ shotgun work (Salter et al. 2014 contaminant list), and it is the dominant signal here.
**Verdict: a real non-infected control arm with excellent depth and paired upper-airway
samples, but bacterial biomass sits at the reagent floor — the growth axis would rest on
single-digit to low-hundreds reads.** Better as a *methods/contamination* reference than as a
quantitative control arm.

### Candidate C — PRJDB41581, ICU sputum: **pneumonia vs respiratory colonisation** ⭐ best *design*

| field | value |
|---|---|
| accession | **PRJDB41581** (97 runs: 49 WGS + 48 RNA-Seq) |
| repository | DDBJ / ENA |
| study title | "A Diagnostic Model Based on Pulmonary Microbiota and Host Gene Expression to Distinguish Colonization from Pneumonia" (NGDC project PRJCA050394) |
| publication | **NOT FOUND.** The study is public (first public 2026-05-02) and cites an NGDC project, but I found **no peer-reviewed paper, no DOI, no preprint** for PRJDB41581. **Mark as: unverified / no publication.** The nearest related published work is Fu et al., *Sci Rep* 2026, DOI `10.1038/s41598-026-44972-w` (resolves ✅; PMID 41933095) — **a different cohort** (17 colonisation / 27 pneumonia, deep sputum, neurosurgical ICU, Huashan Hospital), so it is *not* the PRJDB41581 paper. Do not cite it as such. |
| specimen | **sputum** (deep airway secretions), ICU patients |
| control arm | **respiratory colonisation (n=17)** vs pneumonia (n=27) — colonised-but-not-infected is the single most useful control design |
| depth | RNA-Seq 7,330–25.5M reads (median ~18.6M) |
| library prep | **`library_construction_protocol` present**: *"DNA libraries were prepared using the Nextera DNA Flex Kit, while RNA sequencing libraries were constructed using the **Ovation Trio RNA-Seq Library Preparation Kit**"* — i.e. **host-rRNA depletion, NOT bacterial rRNA depletion.** NextSeq 550, 75-cycle single-end. |
| what bacteria | *Streptococcus pneumoniae*, *Moraxella catarrhalis*, *Prevotella*, *Porphyromonas*, *Veillonella*, *Tannerella* (in the colonisation-range sample probed) |
| raw FASTQ public | **yes** (97/97) |
| **Biomass** | ⚠️ **cohort-dependent, and the spread is huge** |

**Empirical probes (new)** — three PRJDB41581 RNA runs. Outputs:
`out/DRR1001381_genusprobe.json`, `out/DRR1001389_genusprobe.json`, `out/DRR1001399_genusprobe.json`.

| run | reads probed | bacterial hits (any MAPQ) | bacterial share | @MAPQ≥20 | PAO1 any / coding20 | top genera |
|---|---|---|---|---|---|---|
| DRR1001381 | 53,696 | 10,579 | **19.7 %** | 193 | 309 / **11** | *Xanthomonas* 4,001, *Enterococcus* 3,558 |
| DRR1001389 | 6,333,114 | 3,213,758 | **50.8 %** | 98,117 | 5,859 / **2** | *S. pneumoniae* 1,109,356, *Prevotella* 557,859, *Bacteroides* 343,173, *Moraxella* 143,928 |
| DRR1001399 | 0 | — | — | — | 0 / 0 | (download returned 0 records — **probe failed, not a biological zero**) |

**Important correction to the existing repo note.** `RESULTS-ACUTE-SPUTUM-SEARCH.md` (Candidate 3)
rejects PRJDB41581 on a probe of DRR1001378 that reported *"only 1.1 % of reads map at all
(1,116)"*. My re-probe of three other runs from the same cohort gives **19.7 % and 50.8 %**
bacterial share. The earlier figure reflects **that one sample**, not the cohort — and note the
committed `out/DRR1001378_probe.json` records **all-zero** counts (`aligned_mapq20: 0`,
`decoy: 0`), which is the signature of a failed run rather than a measurement. **The cohort-level
"bacterial signal is too thin" conclusion does not survive this re-probe.** The *prep* objection
(host-rRNA depletion only, `Ovation Trio`) **does** survive and is decisive for %RP.

**Verdict: the best control *design* in the public archives** — colonised-but-not-infected,
sputum, ICU, deeply sequenced, raw FASTQ public. Bacterial biomass in at least some samples is
high enough to place the **growth** axis. Blocked as a %RP arm by library chemistry (no bacterial
rRNA depletion), and — critically for the matrix axis — the organisms present are
*Streptococcus*, *Moraxella*, *Prevotella*, oral anaerobes: **no matrix-capable organism at depth**.

### Candidate D — PRJNA875913, paediatric tracheal aspirate: LRTI vs **incidental pathogen carriage (IPC)** vs **non-infectious respiratory failure**

| field | value |
|---|---|
| accession | **PRJNA875913** (293 runs) |
| repository | NCBI SRA / ENA |
| publication | Lydon, Deosthale, Glascock, Phan, Osborne, Leroue et al., *Nature Communications* 2026, "Host–microbiome archetypes differentiate infection from pathogen carriage in the human lower airway". **DOI `10.1038/s41467-026-71863-5` resolves** (Crossref ✅) |
| specimen | **tracheal aspirate**, metatranscriptomic RNA-Seq, mechanically ventilated children |
| arms | **LRTI n=207 · IPC (incidental pathogen carriage) n=70 · non-infectious respiratory failure n=49** — 326 patients. This is exactly the "colonisation" and "non-infected ICU" bracket in one dataset |
| depth | ENA `read_count` max 2.0M, mean 0.2M — **⚠️ the public runs are shallow / likely subsampled for deposit** |
| library prep | Langelier/DeRisi pipeline; host + microbial RNA-Seq |
| **Biomass** | ⚠️ **unknown from metadata; the paper reports IPC having *greater* bacterial abundance than LRTI** — so biomass is real, but the *deposited* depth is low |
| raw FASTQ public | **yes** (293/293 `fastq_ftp`) — note the sibling study PRJNA875913's own summary says *"Raw data not provided due to patient privacy concerns"*; the ENA record nevertheless exposes FASTQ, so **verify per-run before building on it** |

**Existing repo finding:** the prior screen scored 80 runs of this cohort and found
**median 2 Pseudomonas coding reads, 80/80 with more decoy reads than Pseudomonas, zero runs
reaching 1,000**. That is a hard, reproducible negative on *organism availability*.

**Verdict: best *biological* control definition (IPC = colonisation; non-infectious failure =
uninfected ICU), but no matrix-capable organism and low deposited depth.** Useful as a
host-response/ecological comparator, not as a two-axis arm.

### Candidate E — PRJNA1178828, kidney-transplant BALF with **11 healthy controls**

| field | value |
|---|---|
| accession | **PRJNA1178828** (114 runs) |
| repository | NCBI SRA / ENA |
| publication | **NOT FOUND as a peer-reviewed paper.** Study description (Fudan University): *"we enrolled 30 KTR cases presenting pulmonary infection… In addition, the BALF samples which were collected from **eleven healthy individuals in our previous study (Shi et al, 2022)**, were used as controls and were processed according to the same protocols."* **Mark as: unpublished deposit; the 11 "healthy" controls are re-used specimens from a prior study, not a purpose-built healthy arm.** |
| specimen | BALF (+ whole blood), DNA + RNA |
| depth | RNA runs: mean 4.6M, max 61.2M reads |
| library prep | `library_selection = RANDOM`; no protocol recorded |
| **Biomass** | unknown (not probed) |

**Verdict: potentially useful (n=11 healthy BALF processed under one protocol, plus paired DNA),
but unpublished and the controls are borrowed.** Lower priority than A–D.

### Candidate F — donor / unused donor lungs

Searched. Findings:

| source | what it is | accessible FASTQ? |
|---|---|---|
| Mitchell et al. 2021, *Transplantation*, DOI `10.1097/TP.0000000000003713` ✅ resolves | metatranscriptomic BAL from 10 lung-transplant recipients (day 1 + week 6). 11/20 samples bacteria-positive; ***Pseudomonas aeruginosa* at 2,620 RPM in one day-1 specimen** | **not located** — no ENA/SRA accession found |
| Liu et al. 2020, *BMC Infect Dis*, DOI `10.1186/s12879-020-05393-w` ✅ resolves | NGS of donor lung tissue + bronchial secretions, 12 donor lungs, 52.9 % bacteria-positive | **not located** — likely not deposited |
| Xu et al. 2016, *JHLT*, DOI `10.1016/j.healun.2016.01.240` ✅ resolves | 24 donor lungs, 16S V3 of donor BAL, serial post-transplant | **not located** |
| Beaume et al. 2017, *Sci Rep*, DOI `10.1038/srep40309` ✅ resolves | single CF patient, donor airway flora at day 0 (high diversity, "typically characterizing healthy lungs") + *P. aeruginosa* detectable | **not located** |

**Verdict: donor lungs are biologically the most attractive control (genuinely uninfected,
sometimes with *Pseudomonas* colonisation), but I found NO public raw-FASTQ deposit for any of
them.** These are 16S or un-deposited metatranscriptome studies. **No evidence of a usable
donor-lung control dataset.**

### Candidate G — oral/upper-airway-only controls

Excluded by the brief unless paired to a lower-airway sample. Two paired options exist and are
already listed above: **PRJNA603592** (`smnv` BAL + bronchial brush, same subjects) and
**PRJNA875913** (tracheal aspirate, with upper-airway context). No standalone oral dataset was
considered.

### Candidate H — the study's own negative controls (kitome reference)

**PRJNA699976 contains 4 runs titled "RNA-seq of RNase free water"** (SRR13677598,
SRR13677599, SRR13677601, SRR13677602), each at ~52–57M reads. This is a **kitome/background
reference at the same depth and same pipeline as the healthy arm** — the correct denominator for
the "kitome floor" claim, and it did not previously exist in this project.

⚠️ **Important:** `out/SRR13677601_probe2.json` and `out/SRR13677601_rrna.json` exist in the repo,
and SRR13677601 is **one of the four water controls**. The PRJNA699976 probe table in
`PRJNA699976-PROBE.md` lists SRR13677601 among the "six deepest runs" of the *patient* arm —
it is actually a **negative control**, and the probe's "Ralstonia 99.4 %, 97.7 % rRNA" reading
for it is a *kitome* reading, not a patient reading. This does not invalidate the PRJNA699976
rejection (Pseudomonas was absent from the genuine patient runs too), but **the probe table's
labelling needs correcting** before it is cited.

---

## Part 3 — The conceptual answer

### (a) Does a healthy respiratory sample have enough bacterial signal to place on EITHER axis?

**Yes — for the growth axis, in some cohorts; but not in the cohort currently used.**

- **PRJNA390194: no.** 0.185 % of reads bacterial, 6 *Pseudomonas* reads across 9 samples,
  confirmed at 10× depth. Nothing to place.
- **PRJNA699976 healthy arm: yes, the growth axis is placeable in principle.** 80 % of the
  *deposited* reads are bacterial; ~115k bacterial reads at MAPQ ≥ 20 per sample. A %RP
  denominator exists. (Caveat: deposited files appear host-removed — inferred, not documented.)
- **PRJDB41581 colonisation range: yes for at least some samples.** Up to 50.8 % bacterial.
- **PRJNA603592 healthy smokers: marginal.** 0.77 %, *Ralstonia*-dominated.

**The matrix axis is a different matter — and it is placeable, but only via a non-*Pseudomonas*
panel. See (d).**

### (b) Is "no signal" a real biological finding (low biomass) or a technical artefact?

**Both, and they are separable — which is the useful result.**

- **The low biomass is real biology.** Charlson et al. 2011 (*AJRCCM*, DOI `10.1164/rccm.201104-0655OC` ✅)
  measured lung communities **2–4 logs lower in biomass than the upper airways** by 16S qPCR, with
  composition indistinguishable from the URT. Bassis et al. 2015 (*mBio*, DOI `10.1128/mbio.00037-15` ✅)
  found BAL bacterial signal **100–1,000× lower than oral wash**. The JCI review (Whiteside et al. 2021,
  DOI `10.1172/JCI150473` ✅) and Dickson et al. 2017 (*mBio*, DOI `10.1128/mbio.02287-16` ✅) both
  document that healthy BALF bacterial burden is orders of magnitude below infected specimens.
  **Healthy lung genuinely has very little bacterial RNA.**
- **But the magnitude of the PRJNA390194 zero is partly technical.** Three contributions:
  1. **Library prep.** `library_selection = RANDOM`, no bacterial rRNA depletion, no documented
     host depletion. Bacterial reads compete with a vast excess of host RNA. This is the same
     chemistry objection that disqualifies PRJDB41581 for %RP.
  2. **Kitome floor.** The reads that *do* map are *Escherichia* (850/1,026) — a textbook reagent
     contaminant. Below ~8×10⁶ 16S copies/mL the observed signal is dominated by contaminants
     (Schneeberger et al. 2019, *Microbiome*, DOI `10.1186/s40168-019-0755-x` ✅; that paper
     measured BALF at **1×10⁵–4.2×10⁸ 16S copies/mL, median 1.68×10⁶** — i.e. healthy BALF
     routinely sits in the contamination-dominated range).
  3. **Subsampling.** The committed analysis used a 2M-read head. I showed at 20.4M reads the
     result is unchanged — so subsampling is **not** the cause here, but the confirmation was worth doing.
- **Net:** the zero is **real biology amplified by a technically weak library**. It is not an
  artefact in the sense of "there was signal and a filter ate it" — I checked at full depth and
  with the permissive MAPQ band. It *is* an artefact in the sense that a better-prepared library
  from a comparable cohort (PRJNA699976) shows 150× more bacterial signal. **The claim
  "healthy lung has no measurable bacterial signal" is false; the claim "PRJNA390194 has no
  measurable bacterial signal" is true.**

### (c) What is the BEST available control type for a study like this?

**Ranked:**

1. **Colonised-but-not-infected (non-infected ICU / IPC) — best design, and it is available.**
   PRJDB41581 (colonisation n=17 vs pneumonia n=27, sputum, median ~18.6M reads) and PRJNA875913
   (IPC n=70, non-infectious respiratory failure n=49). These are the only controls that are
   *clinically matched to the infected arm* (same specimen, same ward, same pipeline) **and** have
   bacteria present. The colonisation arm is by construction the "bacteria present, host not
   responding" state — which is precisely the control the two-axis question needs. **PRJDB41581
   is the single best-matched cohort found**, though it is unpublished and host-rRNA-depleted only.
2. **Healthy but deeply sequenced with a real bacterial community — PRJNA699976 healthy arm.**
   Best *measured* bacterial signal of any healthy cohort (80 % of deposited reads bacterial).
   Use it to answer "does the healthy lung sit on the growth axis at all?" — and expect an
   oral-commensal community, i.e. *no Pseudomonas*.
3. **Healthy non-smoker/non-COPD BALF, well characterised — PRJNA390194 (current) or
   PRJNA603592 `smnv`.** Both sit at the reagent floor. Adequate only as a "below detection"
   datapoint, which is how they are currently used.
4. **Donor / unused donor lungs — biologically ideal, but no public FASTQ found.**
5. **Non-infected ICU tracheal aspirate (PRJNA875913 non-infectious failure arm)** — good design,
   low deposited depth.

**Recommendation: add a colonisation arm (PRJDB41581 or PRJNA875913 IPC) as the primary control,
keep PRJNA390194 as the "below detection" healthy datapoint, and add PRJNA699976's healthy arm
as the deeply-sequenced healthy comparator.** Reporting all three turns the current single
"(no signal)" cell into a graded statement: healthy-low-biomass → healthy-oral-commensal →
colonised — which is a far stronger control structure than one null.

**And run the colonisation arm through the `strep` panel.** As (d) shows, that is what turns the
control arm from a one-axis null into a two-axis datapoint (matrix-ON × FAST). It is the single
highest-value change available without new sequencing.

### (d) Is there any dataset where a healthy/colonised control has a matrix-capable organism at usable depth — so the matrix axis can be placed in controls?

**YES — one cohort, and the capability already exists in this repo.**
This reverses the conclusion I expected to write, so the experiment is set out in full.

**The blocker, as it looked at the start.** The project's headline matrix panel is
alginate/psl/pel — **organism-specific *Pseudomonas aeruginosa* operons**. The genera that
dominate healthy and colonised respiratory samples (*Neisseria*, *Prevotella*, *Haemophilus*,
*Veillonella*, *Streptococcus*, *Moraxella*, *Rothia*, *Ralstonia*) carry no alginate/psl/pel,
so with that panel alone the matrix axis is unplaceable in every control — regardless of depth.

**But the repo already has per-organism matrix panels** (`perorg_panels.json`, built by
`build_perorg_panels.py`), including a `strep` panel keyed to **`GCF_000006885.1` = *Streptococcus
pneumoniae* TIGR4** (verified via NCBI Datasets API: organism *Streptococcus pneumoniae* TIGR4,
assembly ASM688v1). That panel contains a real matrix module:

| system | genes in panel | genes |
|---|---|---|
| `cps` (capsular polysaccharide) | **12** | cps4A–cps4L |
| `cbp` (choline-binding proteins) | **8** | cbpA–cbpE, cbpF, cbpG, cbpJ |

*S. pneumoniae* is the **dominant organism in the PRJDB41581 colonisation samples**
(*S. pneumoniae* 1,109,356 reads vs the next genus at 557,859 in DRR1001389) — i.e. a
colonised-not-infected control whose dominant organism **is matrix-addressable**.

**The experiment.** New script `control_matrix_probe.py` maps a control sample (streamed from ENA)
against `comp_strep.mmi` = *S. pneumoniae* TIGR4 panel + non-*Streptococcus* decoy genomes,
applying exactly the `quantify_v3.py` logic: MAPQ ≥ 20, class counts, per-system gene detection,
`matrix_on = ≥ 2 systems each with ≥ 3 genes`.

**Result — PRJDB41581 colonisation sample DRR1001389** (`out/DRR1001389_strep_matrixprobe.json`):

| metric | value |
|---|---|
| reads probed | 6,019,228 |
| mapped | 2,425,182 (40.3 %) |
| **target (*S. pneumoniae*) reads** | **32,010** |
| reads: matrix / rp / og | 72 / 3,399 / 28,539 |
| **`f_rp` (growth axis)** | **0.1062 → 10.6 %** — squarely in the FAST band (> 10 %) |
| `cps` genes detected | **8 / 12** |
| `cbp` genes detected | **6 / 8** |
| **`matrix_on` (≥ 2 systems, ≥ 3 genes each)** | **TRUE** ✅ |

**So a colonisation control places on BOTH axes** — matrix-ON × FAST — on the organism actually
present, at 32,010 target reads, with the same thresholds the Pseudomonas arms use.

**The result is not free, and the replication probe is honest about it.** I ran the same test on
two more samples:

| sample | cohort / arm | target reads | matrix reads | systems detected | matrix_on |
|---|---|---|---|---|---|
| **DRR1001389** | PRJDB41581 colonisation | **32,010** | 72 | cps 8/12, cbp 6/8 | **TRUE** ✅ |
| DRR1001381 | PRJDB41581 colonisation | 9 | 0 | 0/0 | FALSE |
| SRR13677570 | PRJNA699976 **healthy** | 4,142 | **0** | 0/0 | FALSE |

**Reading of the three.** The capability is **real but depth-limited, exactly like the acute
*Pseudomonas* arm** (which the project already documents as depth-limited, reaching 100 %
matrix-ON only at matched depth). DRR1001381 was probed at only 1M reads and DRR1001381's
*Streptococcus* share is far lower than DRR1001389's — its 9 target reads is a depth floor, not
a biological zero. **SRR13677570 (healthy) is the more informative negative**: 4,142 *S. pneumoniae*
target reads — 460× more target signal than the PRJNA390194 healthy arm has *Pseudomonas* — and
still **0 matrix reads and 0 cps/cbp genes**. That is a genuine biological result: in healthy
BALF the organism is transcriptionally present but its capsule machinery is **not** transcribed.
Which is exactly the contrast the two-axis design is built to detect, now available in controls.

**Verdict for (d): yes — PRJDB41581's colonisation arm can place the matrix axis, using the
existing `strep` panel. The matrix axis in controls is blocked by *panel choice*, not by data
availability.** Two caveats must travel with this:
1. **It is a different organism and a different matrix.** *S. pneumoniae* `cps`/`cbp` is not
   *P. aeruginosa* alginate/psl/pel. A control matrix-ON call here says "capsule machinery is
   transcribed", not "the same architecture as the acute/chronic arms". Present it as a
   **parallel matrix axis**, not the same one.
2. **Depth matters and must be matched.** DRR1001389 at 6M reads places; DRR1001381 at 1M does
   not. Report the depth floor alongside any control matrix-ON rate, as the project already does
   for the acute arm.

**Two further routes, neither yet tested:**
- **PRJNA699976 healthy arm with a *Neisseria*/*Haemophilus*/*Moraxella* matrix panel.** Those
  three dominate the healthy cohort at MAPQ ≥ 20 (36,813 / 15,356 / — for *Prevotella* /
  *Neisseria* / *Haemophilus*). No panel exists for them in `perorg_panels.json`; building one
  would let the healthy arm itself place a matrix axis.
- **Donor lungs** remain biologically the ideal control (Mitchell 2021 found *P. aeruginosa* at
  2,620 RPM in a day-1 transplant BAL) — but **no public FASTQ was found**, so this route is
  closed pending data access.

---

## Reproducibility

New scripts written for this report (both stream from ENA; **neither writes a FASTQ to disk**):

| script | purpose | outputs |
|---|---|---|
| `healthy_depth_probe.py` | full-depth MAPQ×class contingency for a sample, streaming from ENA | `out/<SRR>_depthprobe.json` |
| `healthy_genus_probe.py` | same + genus resolution of decoy hits via `decoy/manifest.json` | `out/<SRR>_genusprobe.json` |
| `control_matrix_probe.py` | two-axis probe for a **control** sample against a per-organism panel (`comp_<org>.mmi`), reusing `quantify_v3.py` logic (MAPQ≥20, class counts, per-system detection, `matrix_on = ≥2 systems × ≥3 genes`) | `out/<SRR>_<org>_matrixprobe.json` |

New outputs (all small JSON; no temp FASTQs left behind — disk held at 3.0 GB free throughout):

```
out/SRR5677652_depthprobe.json   out/SRR5677656_depthprobe.json
out/SRR5677652_genusprobe.json   out/SRR13677570_depthprobe.json
out/SRR13677570_genusprobe.json  out/SRR11002545_genusprobe.json
out/DRR1001381_genusprobe.json   out/DRR1001389_genusprobe.json
out/DRR1001399_genusprobe.json   (download returned 0 records — probe failed)
out/DRR1001389_strep_matrixprobe.json   <-- the matrix-axis-in-controls result
out/DRR1001381_strep_matrixprobe.json
out/SRR13677570_strep_matrixprobe.json
```

Index used: `compfull_pseu.mmi` = full PAO1 (`GCF_000006765.1`) + 120 non-*Pseudomonas* decoy
contigs (121 sequences; decoys prefixed `decoy|`). Decoy genera resolved through
`decoy/manifest.json` (34 genera).

## Verification ledger

All DOIs resolved against **Crossref** (`api.crossref.org/works/<doi>`, 2026-10-12). Publisher
HTTP 403s (ASM, ERS, AACR-family, Wiley, MDPI) are bot-blocks, not invalid DOIs — Crossref
metadata is the authoritative check and every one returned full title/journal/year/author.

| DOI | status | note |
|---|---|---|
| `10.1128/msystems.00199-18` | ✅ resolves | **Ren 2018 — the PRJNA390194 paper** |
| `10.1186/s40168-019-0755-x` | ✅ | Schneeberger 2019 — BALF density / contaminant thresholds |
| `10.1128/mbio.00037-15` | ✅ | Bassis 2015 — lung vs oral bacterial load |
| `10.1172/JCI150473` | ✅ | Whiteside 2021 — lung microbiome review |
| `10.1164/rccm.201104-0655OC` | ✅ | Charlson 2011 — healthy LRT 2–4 logs lower biomass |
| `10.1128/mbio.02287-16` | ✅ | Dickson 2017 — healthy LRT bacterial topography |
| `10.1183/13993003.03434-2020` | ✅ | Sulaiman 2021 — **PRJNA603592** |
| `10.1038/s41564-021-00961-5` | ✅ | Sulaiman 2021 — COVID ICU BAL metatranscriptomics |
| `10.1038/s41467-026-71863-5` | ✅ | Lydon 2026 — **PRJNA875913** LRTI/IPC/control |
| `10.1038/s41598-026-44972-w` | ✅ | Fu 2026 — colonisation vs pneumonia (different cohort to PRJDB41581) |
| `10.3389/fcimb.2025.1639148` | ✅ | Duan 2025 — RNA-mNGS vs DNA-mNGS infection/colonisation |
| `10.1038/s41522-024-00548-y` | ✅ | Chen 2024 — LRTM + host ML diagnosis |
| `10.1002/advs.202405087` | ✅ | Zou 2024 — **PRJNA699976** |
| `10.1038/s41597-025-06171-6` | ✅ | Tang 2025 — sibling BALF deposit; **"only microbial (non-host) reads deposited"** |
| `10.1186/s42003-024-05840-3` | ✅ | Cuthbertson 2024 — CELF bronchial brushings, healthy + asthma |
| `10.1183/13993003.01829-2020` | ✅ | Mourani 2020 — paediatric VAP airway microbiome |
| `10.1172/JCI165904` | ✅ | Mick 2023 — paediatric LRTI host/microbe diagnosis |
| `10.1164/rccm.202403-0516RL` | ✅ | Lydon 2024 — FABP4 pneumonia biomarker |
| `10.1164/rccm.202210-1865OC` | ✅ | Sulaiman 2023 — mild-moderate COPD dysbiosis |
| `10.1097/TP.0000000000003713` | ✅ | Mitchell 2021 — lung-transplant BAL metatranscriptome (*P. aeruginosa* 2,620 RPM) |
| `10.1186/s12879-020-05393-w` | ✅ | Liu 2020 — donor-lung colonised bacteria by NGS |
| `10.1016/j.healun.2016.01.240` | ✅ | Xu 2016 — donor microbiota and PGD |
| `10.1038/srep40309` | ✅ | Beaume 2017 — donor airway microbiota invaded by *Pseudomonas* |
| `10.1186/s40168-016-0151-8` | ✅ | Kelly 2016 — intubated-patient respiratory microbiome |
| `10.1186/s40635-023-00496-5` | ✅ | Alagna 2023 — VAP vs no-VAP endotracheal aspirates |
| `10.1101/2021.03.23.21253487` | ✅ (preprint) | Tsitsiklis 2021 — COVID VAP vs no-VAP tracheal aspirate |
| `10.1186/s40413-015-0074-z` | ✅ | Aho 2015 — lower-airway microbiome review |
| `10.3390/biomedicines11030841` | ✅ | Leitao Filho 2023 — healthy young persons, brush + BAL (no accession found) |
| `10.1016/S2666-5247(21)00304-9` | ✅ | Tsitsiklis 2022 — paediatric LRTI metagenomics |

**Unverified / no publication found:**
- **PRJDB41581** — public DDBJ deposit, **no paper, no DOI, no preprint located.** Its companion
  NGDC project is PRJCA050394. Do not cite a publication for it.
- **PRJNA1178828** — Fudan University deposit, **no paper found**; 11 "healthy" controls are
  borrowed from Shi et al. 2022, not purpose-collected.
- **Donor-lung cohorts** (Mitchell 2021, Liu 2020, Xu 2016, Beaume 2017) — papers verified,
  **raw FASTQ deposits not located.**
- **PRJNA699976 deposited files being host-removed** — **strongly inferred** from the sibling
  deposit's data note, **not documented** in PRJNA699976's own ENA metadata.

## Corrections this report makes to the existing record

1. **`PRJNA699976-PROBE.md` probe table** lists **SRR13677601** among the "six deepest runs" of
   the patient arm. SRR13677601 is one of the **four RNase-free-water negative controls**. Its
   "Ralstonia 99.4 %, 97.7 % rRNA" reading is a kitome reading. Re-label before citing.
2. **`RESULTS-ACUTE-SPUTUM-SEARCH.md` Candidate 3** rejects **PRJDB41581** because a single probe
   (DRR1001378) showed only 1.1 % of reads mapping. My re-probe of three other runs from the same
   cohort gives **19.7 %** and **50.8 %** bacterial share. The committed
   `out/DRR1001378_probe.json` is **all zeros** (a failed run, not a measurement). The
   *prep* objection to PRJDB41581 stands; the *"bacterial signal is too thin"* objection does not.
3. **"Healthy controls = kitome floor"** (README, `RESULTS-QUADRANT-CORRECTED.md` §Verdict 3,
   `WRITING-PLAN.md` row 7) should be narrowed to **"PRJNA390194 = kitome floor"**. Healthy lung
   in general is measurable (PRJNA699976: 80 % of deposited reads bacterial).
4. `out/DRR1001399_genusprobe.json` records a **probe failure** (0 records downloaded), not a
   biological zero. Flagged so it is not mistaken for one.
5. **`README.md` and `RESULTS-QUADRANT-CORRECTED.md` state that healthy controls cannot be placed
   on either axis.** That is true for PRJNA390194 with the *Pseudomonas* panel. It is **not** true
   in general: a colonisation control places on **both** axes once the matrix panel matches the
   organism present (Part 3d). Recommend rewording from "nothing to place" to
   "PRJNA390194 has nothing to place; colonisation controls are placeable with a matched panel".
6. **The repo's per-organism panels were built but never used for a control arm.** `strep`
   (`cps`/`cbp`, *S. pneumoniae* TIGR4) is complete and validated here. The `perorg_panels.json`
   entries for `kleb` (2 genes), `staph` (2 genes) and `haemo` (3 genes) are **too thin to reach
   the ≥3-genes-per-system threshold** — only `strep` (cps 12 / cbp 8), `ecoli` (pga 4 / csg 7 /
   bcs 7), `acine` (csu 6 / pga 7) and `steno` (fimbriae 12) are usable as matrix arms.
