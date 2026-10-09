# Matrix machinery is transcribed but not co-expressed with growth in acute bacterial pneumonia

### A metatranscriptomic test of the acute-biofilm model, and a caution about depth

**Draft — 2026-10-09**

---

## Abstract

Biofilm matrix production is conventionally associated with chronic infection and
slow or arrested growth. Kolpen et al. (*Thorax* 2022) challenged this using
PNA-FISH and confocal imaging of sputum, reporting that biofilms dominate acute as
well as chronic lung infection and that the two differ in metabolic rate rather
than architecture — implying matrix production coincides with *fast* growth in
acute infection. That claim rests on imaging and has never been tested
transcriptionally. We tested it in public BALF metatranscriptomes from 114
patients with bacterial pneumonia (PRJNA1056765), quantifying matrix-machinery and
growth-machinery transcripts from a single mapping run against a per-organism
panel augmented with 34 co-infecting genera as mapping decoys.

Matrix transcripts were readily detectable: the complete alginate operon (16/16
genes), Psl (7/11) and Pel (3/7) were recovered, which to our knowledge is the
first direct transcript detection of this machinery in acute bacterial BALF.
However, matrix was not over-represented: the log-ratio balance
`B = log(matrix) − log(growth)` had median +0.05 against the full growth panel and
−0.41 against ribosomal proteins alone, so its sign depends on how growth is
defined. More importantly, the apparent co-expression of matrix and growth
(r = +0.50) proved to be a sequencing-depth artifact: both arms track depth
(r = 0.91 growth, 0.77 matrix), and partialling depth out inverted the association
to −0.76 — the classical trade-off. The transcriptional signature predicted by
the acute-biofilm model is therefore not visible in these samples. We also show
that two design choices commonly made in this kind of analysis — selecting
samples on a deposited composition table and mapping without decoy competition —
each independently produce a wrong answer, and we document both.

---

## 1. Introduction

The canonical model of biofilm infection holds that matrix-encased communities
represent a slow-growing, chronic state, and that acute infection is
characterised instead by planktonic fast growth. Kolpen et al. (2022) tested this
in sputum from patients with acute and chronic lung infection using PNA-FISH
targeting rRNA and confocal imaging, and reported that biofilm aggregates were
present in *both*, with the difference between the states lying in metabolic
activity rather than architecture. Their model therefore predicts that in acute
bacterial infection, matrix production and fast growth occur **together** — the
opposite of the classical trade-off.

This prediction is transcriptionally testable and has not been tested. The
obstacle has been data: measuring matrix-gene expression in vivo requires
respiratory specimens with both adequate microbial RNA and paired metagenomic
information to distinguish gene-absent from gene-silent. The recently published
PRJNA1056765 resource (Tang et al. 2025, *Sci Data*; 402-patient BALF
metagenomics + metatranscriptomics, host-depleted) makes this possible for the
first time at scale.

We set out to test the prediction directly, and report a negative result together
with two methodological failures that had to be fixed before the result was
trustworthy — and which, we suspect, are not unique to this study.

---

## 2. Results

### 2.1 The deposited composition table and the RNA reads disagree

The dataset ships a per-sample genus composition table. Because our measurement
is read-level, we derived genus composition independently from the RNA reads
themselves, mapping each sample to a 19-genome multi-strain panel. The two sources
agree on the dominant genus for only **18 of 36** samples examined, and every
disagreement has the same shape: the table assigns *Pseudomonas* as dominant, the
reads assign *Acinetobacter*. For `SRR27343990`, the table gives 25.3%
*Pseudomonas* / 7.6% *Acinetobacter* while the reads give 11% / 52%.

This matters because the natural workflow — pick matrix-measurable samples from
the table, then measure from the reads — selects and measures on different
objects. We therefore defined the cohort by the reads, so that selection and
measurement share one source.

### 2.2 Matrix transcripts are present in acute bacterial BALF

Matrix machinery is transcribed in these samples. Aggregating over the
10-sample cohort defined in §2.3, 40 of 101 panel matrix genes were detected
(1,723 reads), including:

- **alginate operon: 16/16** — `algA`, `algB`, `algC`, `algD`, `algE`, `algF`,
  `algG`, `algI`, `algJ`, `algK`, `algL`, `alg8`, `alg44`, `algP`, `algR`, `algU`
- **psl: 7/11** — `pslA`, `pslB`, `pslC`, `pslD`, `pslE`, `pslG`, `pslH`
- **pel: 3/7** — `pelA`, `pelB`, `pelD`
- regulators `amrZ`, `algU`, `algR`, `algB`, `cdrB`, `siaA`

The most consistently detected genes were `algC` (10/10 samples, 294 reads),
`amrZ` (10/10, 97), `algA` (9/10, 705) and `algU` (8/10, 158). Existing in-vivo
evidence for matrix-gene expression is immunohistochemical and confined to chronic
CF sputum (Jennings et al. 2021); these data show the machinery is transcribed in
acute infection, and that acute matrix production need not be inferred by analogy
from chronic disease.

### 2.3 No matrix excess, and the sign depends on how growth is defined

We computed the log-ratio balance

    B = log(R_matrix) − log(R_growth)

where `R` is length-normalised transcript abundance (reads per kb of class panel),
against two definitions of the growth arm: `B_rp` (ribosomal-protein genes only —
the canonical growth marker) and `B_full` (ribosomal proteins plus replication,
division, transcription and translation factors).

| metric | median | mean | range | n>0 | sign test | Wilcoxon z |
|---|---|---|---|---|---|---|
| `B_full` | **+0.05** | −0.14 | [−1.14, +0.73] | 5/10 | 1.00 | −0.66 |
| `B_rp` | **−0.41** | −0.53 | [−1.52, +0.27] | 3/10 | 0.34 | −1.58 |

Matrix share of target reads was median 0.416 (range 0.068–0.606). Neither
definition yields a significant departure from zero, and the two disagree in sign.
An RP-only growth arm — a natural choice, and the one we initially made — reports
a matrix *excess*; the full growth panel does not. Since the two differ only by
the inclusion of non-ribosomal growth genes, and those turn out to dominate the
growth arm (see §2.4), the discrepancy is not incidental.

### 2.4 The apparent matrix–growth co-expression is a depth artifact

Taken at face value, the data look Kolpen-compatible: across samples, matrix and
growth abundance rise together (Pearson r = **+0.50**). But both arms track
sequencing depth almost perfectly:

| relationship | Pearson r |
|---|---|
| log R_growth vs log depth | **+0.91** |
| log R_matrix vs log depth | **+0.77** |
| log R_matrix vs log R_growth (raw) | +0.50 |
| **log R_matrix vs log R_growth, depth partialled out** | **−0.76** |

Partialling depth out inverts the association. The mechanism is a ceiling: at low
coverage, nearly every mapped read belongs to a high-abundance growth transcript
(dominated by `fusA`, `nusA`, `rpoB/C`, `tuf`), so the growth rate climbs steeply
with depth while the matrix rate remains floored near the detection limit. The
raw positive correlation is therefore a property of coverage, not of regulation.
Once depth is controlled, the association is strongly **negative** — the
classical matrix/growth trade-off, not the co-expression the acute-biofilm model
predicts.

![Matrix vs growth balance](fig_balance.svg)

*Left: per-sample balance under both growth definitions. Centre: both arms scale
with sequencing depth. Right: the matrix–growth correlation inverts when depth is
removed.*

### 2.5 Decoy competition changes the growth arm

Conserved genes map near-identically across bacterial genera. With only the target
organism in the reference index, such reads have no competitor, receive high
mapping quality, and are counted as target-organism transcripts even when they
originated from a co-infecting species. In an early version of this analysis this
inflated a single conserved gene (`fusA1`) to 31% of all mapped reads in one
sample — a signature that read-level inspection showed was a genuine expressed
transcript (386 distinct sequences spanning the full gene, both strands), but
which could not be attributed to the target organism without competitors in the
index. Adding one reference genome per co-infecting genus (34 genomes) absorbed
cross-genus reads and made the growth arm interpretable.

---

## 3. Discussion

**The acute-biofilm prediction is not supported transcriptionally.** Kolpen et al.
reported that biofilm aggregates are present in acute infection and that the
acute/chronic distinction lies in metabolic rate. Our data are consistent with
the first part — matrix transcripts, including the full alginate operon, are
present — but not with the prediction that follows from the second: we find no
positive coupling between matrix and growth investment, and a strong negative
coupling once depth is controlled. In transcriptomic terms these samples sit in
the classical matrix-ON/growth-OFF quadrant.

Two caveats constrain how strongly this can be stated. First, a bulk
metatranscriptomic ratio cannot adjudicate a claim about spatial architecture and
the metabolic state of cell aggregates; a negative result means the predicted
*transcriptional signature* is absent, not that the imaging is wrong. Had the
result been positive it would have been stronger evidence than imaging alone; a
negative result is correspondingly weaker than a refutation. Second, the cohort is
small (n = 10), and it is small because depth is limiting, which is the same
variable that confounds the analysis. A designed experiment with matched depth,
or spike-in normalisation, is required for a clean test.

**Two design choices each produce a wrong answer.** We report both because we
made both.

1. *Selecting on a composition table while measuring on reads.* The deposited
   table and the RNA reads disagreed on the dominant organism in half the samples,
   always in the same direction. Selecting matrix-measurable samples from the
   table and then quantifying from the reads is a mismatch between the object
   selected and the object measured. In a dataset where the two dominant genera
   differ sharply in matrix-gene coverage, this silently changes the cohort.

2. *Mapping without decoy competition.* Conserved genes without competitors in
   the index become sinks for reads from every co-infecting species. Because such
   genes are heavily enriched for growth functions, the artifact inflates the
   growth arm specifically — and it does so in a way that is invisible without
   inspecting individual gene counts.

**A general caution.** Both failures, and the depth confound, share a structure:
a measurement choice that is defensible in isolation produces a directional
artifact in the quantity of interest. The depth confound is the most dangerous of
the three because it points in the direction of the hypothesis under test. We
recommend that analyses of this kind report the matrix/growth balance under more
than one growth definition, inspect per-gene counts rather than class totals, and
test explicitly whether the balance is a function of depth.

**Acinetobacter remains unmeasurable.** It is the RNA-dominant organism in
roughly half the bacterial samples in this dataset, and it has essentially no
characterised matrix operon — only two matrix genes could be annotated from its
genomes. The organism with the most biomass is the one that cannot be scored.
This is a coverage limit, not a biological null, and it means our result is a
statement about *Pseudomonas*, not about acute bacterial pneumonia in general.

---

## 4. Methods

**Data.** PRJNA1056765 (Tang et al. 2025, *Sci Data*). 114 bacterial-infection
BALF RNA runs with paired DNA; 50 bp single-end; host-depleted. Non-infectious
control samples contained essentially no bacterial signal (2 mapped reads) and
were not usable as a contrast, so the design is within-bacterial.

**Reads.** ENA HTTPS where available; SRA ODP with `fasterq-dump` for runs absent
from the ENA mirror.

**Panels.** Per organism (Pseudomonas, Klebsiella), gene sequences extracted from
NCBI GFF3 annotations of three strains each and classified as `matrix`
(alg/pel/psl/csg/bcs/ica/cdr/EPS-capsule; 101 and 57 genes), `rp` (ribosomal
protein; 150 and 102) or `og` (other growth functions; 241 and 282).

**Decoy competition index.** `comp_<org>.mmi` = target panel + one reference
genome for each of 34 co-infecting genera (Salmonella, Delftia, Sphingomonas,
Rhodococcus, Xanthomonas, Burkholderia, Campylobacter, Comamonas, Ralstonia,
Porphyromonas, Fusobacterium, Treponema, Neisseria, Bacteroides, Tannerella,
Prevotella, Mycoplasma, Moraxella, Haemophilus, Streptococcus, E. coli,
Acinetobacter and others), covering the taxa observed in the cohort. Only
target-tagged references are counted.

**Mapping and quantification.** minimap2 2.28, `-x sr`, MAPQ ≥ 20, `--secondary=no`,
in a **single run per sample** so that both arms share identical mapping
conditions. Abundance `R_class` = class reads / (sum of class gene lengths / 1000),
i.e. reads per kb, which removes the gene-length bias that would otherwise
dominate a panel comparison. Balance `B = log(R_matrix) − log(R_growth)`.

**Cohort.** Samples with *Pseudomonas* fraction of mapped RNA reads ≥ 0.5 and
≥ 100 target reads in the quantification panel (n = 10).

**Statistics.** Two-sided sign test and Wilcoxon signed-rank test against zero for
`B`; Pearson and Spearman correlation, and first-order partial correlation
controlling for log sequencing depth, for the coupling analysis. Read-level
permutation intervals (resampling observed reads across panel genes, 3,000
replicates) are reported per sample.

**Code and data.** All panels, scripts, per-sample records and the figure are
available at `github.com/joshux/acute-biofilm-analysis` (commit `045723a` and
subsequent). The analysis is reproducible from the block at the end of
`RESULTS.md`.

---

## Limitations

- n = 10, limited by sequencing depth — the same variable that confounds the
  balance.
- Depth is not randomly distributed with respect to matrix share, so the
  depth-controlled estimate rests on extrapolation across a narrow range.
- rRNA depletion was human-only; bacterial rRNA is likely still present, diluting
  mRNA and lowering power.
- "Acute" reflects the dataset's "Bacterial infection" label, not an independent
  clinical adjudication.
- Bulk averaging erases spatial structure, so the result cannot speak to
  architecture.
- The result concerns *Pseudomonas*; *Acinetobacter*, the other dominant organism,
  cannot be scored with current annotations.
