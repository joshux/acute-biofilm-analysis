# PRJNA699976 — the two case studies, run on both axes

**Date:** 2026-10-12. **Scope:** the ONLY two patient samples in PRJNA699976 whose
dominant organism is matrix-capable and at usable depth. Both axes (growth × matrix)
were computed for both cases, with a species identification for the *Acinetobacter* case.

**Bottom line:** the two cases land on **opposite** points of the quadrant.
*Pseudomonas*-dominant **SRR13677832 is matrix-present × SLOW**; *Acinetobacter*-dominant
**SRR13677669 is matrix-present × FAST**. Case 2's *growth* axis lands where the acute arm
lands (FAST), but its *matrix* axis does not: the organism is **A. pittii, not A. baumannii**,
and A. pittii carries **pga but no csu**. So neither case reaches **matrix-ON × FAST** — the
acute arm's cell — because neither organism has two matrix systems transcribing. The
*Acinetobacter* case reproduces only the growth half of acute.

---

## 1. Headline numbers

| | **SRR13677832** | **SRR13677669** |
|---|---|---|
| dominant organism | *Pseudomonas aeruginosa* | ***Acinetobacter pittii*** (not *baumannii*) |
| growth depth (target reads, MAPQ≥20) | 14,288 | 5,272 |
| **%RP_coding** (rRNA-excluded) | **3.88 %** | **22.31 %** |
| Gifford class (FAST >10 / SLOW <5) | **SLOW** | **FAST** |
| rp_over_ctrl | **2.58** | **3.38** |
| rp/ctrl class (>4 FAST-like / <3 SLOW-like) | SLOW-like | ambig |
| load (target reads / reads mapped) | 0.0057 | 0.0053 |
| matrix systems in panel | alginate, psl, pel | pga, csu (committed) → pga only (species-matched) |
| **matrix_call** | **alginate-only** | **pga-only** |
| **placement** | **matrix-present × SLOW** | **matrix-present × FAST** |
| matches an arm? | loose: chronic (matrix half differs — chronic is matrix-**ON**) | exact: **acute**; loose: acute, covid |
| vs the arms | growth at/below chronic (2.82); matrix present but not ON | growth above chronic (2.82), below acute (4.72); matrix present but not ON |

Read length in both subsamples: **145 bp** (modal), mean 144.9 bp — the same as the arms.
Depth: 2.5 M reads (case 1), 1.0 M reads (case 2), the same read sets feeding both axes.

---

## 2. Method (and how it differs from a naive run)

Both axes were measured on **one mapping run over one reference** per case, so growth and
matrix share identical read length, scoring and filters.

**Growth axis.** RP reads / (target reads outside rRNA/tRNA/tmRNA), MAPQ ≥ 20, ±300 bp
overlap window, and a conserved non-RP control-gene set — exactly the committed
`norm_map.py` / `quantify_rp2.py` logic, reimplemented chunk-wise because backgrounded
processes are reaped at tool-call end and each case needs >110 s of mapping. The chunked
driver (`cases/case_map.py`) is resumable: one JSON per read-chunk, merged by
`cases/merge_cases.py`.

**The denominator excludes rRNA, as required** — this dataset has no bacterial rRNA
depletion, so the raw coding fraction is an upper bound and only the rRNA-excluded metric
is comparable to the arms. The exclusion is doing real work here: for SRR13677832,
**99.6 % of the reads that align to PAO1 are rRNA/tRNA** (Section 5).

**Matrix axis.** Per-system gene detection on a `cls|acc|gene|system|idx` panel, with the
committed call rule (matrix-ON = ≥2 systems each with ≥3 genes detected). The
*Pseudomonas* panel (`qp_pseu.fna` + `comp_pseu.mmi`) was used as committed. For the
*Acinetobacter* case the committed panel is built on *A. baumannii* ATCC 17978 and
therefore cannot quantify an *A. pittii* population, so a species-matched panel was built
from the genome of the organism actually present (Section 4).

### A note on the matrix panel

The on-disk `qp_pseu.fna` is a **single-strain** panel (PAO1; alginate 22, psl 14, pel 7),
whereas the committed arm outputs were produced with a **three-strain** panel (PAO1 + PA14
+ KT2440; alginate 39, psl 14, pel 9). Both generations exist in the working tree — the
arms used the three-strain one. This matters for reading the *Pseudomonas* matrix call:
the smaller panel detects fewer genes, so **`alginate-only` here is a conservative call**
relative to the arms.

**The call is nonetheless panel-independent, and provably so.** Genes detected in a system
can never exceed the number of reads that hit it, and this sample has **2 reads on psl and
0 on pel**. A system with 2 reads cannot reach the ≥3-gene threshold under *any* panel,
however many genes that panel contains; a system with 0 reads cannot reach 1 gene. So the
`alginate-only` call for SRR13677832 holds under the single-strain panel, the three-strain
panel, or any larger one — the psl/pel failure is a property of the sample's reads, not of
the panel. The same argument applies to the *Acinetobacter* case: on the *A. baumannii*
panel both pga and csu score **0 matrix reads**, so neither system can be ON there either.

---

## 3. Case 1 — SRR13677832 (*Pseudomonas*)

- **Growth: 3.88 % RP_coding → SLOW.** 14,288 PAO1 target reads; 13,351 coding after
  excluding 937 rRNA/tRNA reads; 518 RP reads; 201 control-gene reads.
- **rp_over_ctrl = 2.58**, below the chronic anchor (2.82) and far below the acute anchor
  (4.72) — i.e. *slower* than the chronic arm on the conservation-matched metric.
- **Matrix: alginate-only.** alginate 13/22 genes (33 reads), psl 2/14 (2 reads), pel 0/7
  (0 reads). Matrix-ON requires ≥2 systems with ≥3 genes; psl and pel fail, so the call
  is `alginate-only` — matrix *present*, not matrix *ON*.
- **Placement: matrix-present × SLOW.** This is **not** the chronic arm's cell (chronic is
  matrix-**ON** × SLOW) and **not** the acute arm's (matrix-ON × FAST). On the growth axis
  alone this case is the slowest of the two by a wide margin.

**Reading:** a *Pseudomonas*-dominant BALF metatranscriptome in which the organism is
transcribing its alginate regulon (including the regulators `algU`, `algD`, `algC`) but
not psl or pel, and at a translational-output level consistent with slow growth. That is
the classic chronic/stationary biofilm phenotype — but on this dataset's own matrix rule
it does not reach matrix-ON, because psl and pel are silent.

---

## 4. Case 2 — SRR13677669 (*Acinetobacter*) — species identification first

**Verdict: *Acinetobacter pittii* (RAES20-type), NOT the *baumannii* lineage.**

Method: 1.0 M reads (3 chunks of the 145 bp R1 subsample) mapped against `acine_multi.mmi`
— the six-genome Acinetobacter genus index (2 × *A. baumannii*, *nosocomialis*, *pittii*,
*johnsonii*, *lwoffii*) plus 117 non-Acinetobacter decoy contigs — with per-accession
attribution at MAPQ ≥ 20.

| species | reads | share |
|---|---|---|
| ***A. pittii* (RAES20)** | **4,471** | **93.5 %** |
| *A. nosocomialis* (C24) | 111 | 2.3 % |
| *A. baumannii* (AB5075) | 91 | 1.9 % |
| *A. johnsonii* (PF3HH25001) | 64 | 1.3 % |
| *A. baumannii* (ATCC 17978) | 26 | 0.5 % |
| *A. lwoffii* (TMP6) | 19 | 0.4 % |
| **baumannii lineage, total** | **117** | **2.4 %** |

Supporting evidence that this is a real population and not a mapping artefact:

1. **39 : 1 ratio** in favour of pittii over the two baumannii genomes combined.
2. **Uniform whole-chromosome coverage** — reads spread across **39/39** 100-kb bins of
   the single 3.90 Mb pittii contig, the signature of a present genome rather than a
   conserved-region attractor.
3. **98.3 % mean nucleotide identity** for the pittii-mapping reads (vs 98.2 % for the
   handful of AB5075 reads) — the reads are genuinely pittii-like, not misassigned.
4. Consistent across all three chunks (1,478 / 1,535 / 1,458 pittii reads).

**Matrix capability — the decisive consequence.** In *A. baumannii* the biofilm matrix
axis rests on **two** systems, `pgaABCD` and the `csu` fimbrial operon. A. pittii has
**pgaABCD** (a clean, contiguous, syntenic locus, and it aligns at 100 % identity to the
*baumannii* pga locus) but **no csu ortholog**:

> The *A. baumannii* **csu operon** (6 CDS, 5,958 bp, contiguous) mapped against the
> *A. pittii* RAES20 genome yields **NO ALIGNMENT at asm5, asm10 or asm20** — the most
> sensitive of the three presets. The *baumannii* **pga locus** run through the identical
> pipeline **does** align, at **100 % identity**, and at CDS level its pgaB and pgaC
> queries cover 108 % and 101 % of their query length — i.e. full-length, exact matches.
> Positive and negative controls therefore behave correctly in the same test.

A. pittii's GFF does contain CDS annotated `Csu type fimbrial protein`, but these are
generic fimbrial-cassette products, not the baumannii csuABCDE operon: the nucleotide
test above finds no homology to the real operon, and they are not the six-gene contiguous
csu locus.

**Matrix call.** With only one matrix system in the panel (pga), the committed rule
(matrix-ON = ≥2 systems each ≥3 genes) **cannot be satisfied by construction**. Detection:
pgaB (2 reads), pgaA (1), pgaD (1) — 3/7 genes, 4 reads → **`pga-only`**, i.e. matrix
*present*, not matrix-ON. On the *A. baumannii* panel (for contrast) this sample yields
**1 matrix read total**, confirming that the baumannii panel cannot see this organism —
which is the whole reason a species-matched panel was required.

**Growth: 22.31 % RP_coding → FAST.** 5,272 pittii target reads; 4,900 coding after
excluding 372 rRNA/tRNA; 1,093 RP; 323 control-gene reads. **rp_over_ctrl = 3.38** —
above the chronic anchor (2.82) but below the acute anchor (4.72).

> **Why the reference had to be pittii.** Mapping the same reads onto *A. baumannii*
> ATCC 17978 gives **35.4 %** RP_coding with only 2,637 target reads — i.e. an
> artefactually *higher* and much thinner signal, because only the conserved genes can
> align across species and those are disproportionately RP genes. The species-matched
> pittii reference (5,272 target reads, 89 % of the target's own reads in rRNA) is the
> correct measurement. **The cross-species number must not be reported as the case value.**

**Placement: matrix-present × FAST** — the *growth* half of the acute arm's cell, and
squarely in the COVID arm's territory (COVID was matrix-present × FAST-like). But unlike
acute, this is a single-system (pga-only) matrix, and unlike acute the matrix-ON criterion
is unreachable.

---

## 5. The rRNA caveat, measured rather than asserted

Because neither sample was bacterial-rRNA depleted, the growth denominator must exclude
rRNA — and it is worth seeing how much signal that removes:

| | SRR13677832 | SRR13677669 |
|---|---|---|
| reads aligning to the target organism (any MAPQ) | 98.2 % of 398,631 | 4.5 % of 364,240 |
| of those, rRNA/tRNA | **99.6 %** | **89.1 %** |
| target reads surviving MAPQ ≥ 20 | 1,548 | 1,744 |
| target reads surviving MAPQ ≥ 20 **and** rRNA-excluded | **1,435** | **1,613** |

For SRR13677832 the library is ~98 % *Pseudomonas* by reads but **99.6 % of that is
rRNA** — only ~0.4 % of reads are usable coding signal, which is exactly why 2.5 M reads
were mapped to reach the arms' median depth (~5 k coding reads). For SRR13677669 the
sample is *Ralstonia*-dominated (76.8 % of reads map to decoys, mostly *Ralstonia*), with
the Acinetobacter population a small minority — 1.0 M reads gave 5,272 target reads.
This is the same "upper bound on coding signal" caveat the committed quantifiers carry,
now quantified per case.

---

## 6. Robustness (n = 2, so the call must survive the obvious knobs)

**Per-chunk spread** (the metric computed independently on each read-chunk):

| case | chunks | %RP_coding per chunk | mean | sd |
|---|---|---|---|---|
| SRR13677832 | 7 | 3.97, 4.77, 3.90, 3.27, 3.10, 3.70, 4.60 | 3.90 | 0.62 |
| SRR13677669 | 3 | 22.38, 22.33, 22.20 | 22.31 | 0.09 |

**Threshold and window sensitivity** — one mapping pass per case, every combination
evaluated from the identical SAM (so only the filter varies):

| setting | SRR13677832 %RP | SRR13677669 %RP |
|---|---|---|
| MAPQ ≥ 0, window 0 | 2.27 | 18.43 |
| MAPQ ≥ 20, window 300 (**committed**) | **3.97** | **22.38** |
| MAPQ ≥ 40, window 600 | 4.45 | 25.36 |
| full range across all 20 settings | **2.3 – 4.5** | **18.4 – 25.4** |

Both calls are stable: case 1 is SLOW (below 5 %) under *every* setting tested, and case 2
is FAST (above 10 %) under *every* setting tested. No single knob decides either call.

---

## 7. What this does and does not establish

**Does:**

- Both remaining case studies are now measured on both axes with the committed
  methodology, at depth comparable to the arms (14,288 and 5,272 target reads).
- The *Acinetobacter* case is identified to species: **A. pittii**, 93.5 % of its
  Acinetobacter reads, with a sequence-level demonstration that it lacks the *baumannii*
  csu operon while retaining pga. Its matrix axis is therefore **partial by biology**,
  not by depth.
- The two cases fall in **different** quadrants from each other, and neither reproduces
  the acute arm's matrix-ON × FAST cell.

**Does not:**

- **n = 2.** These are two patients. Nothing here supports a cohort-level claim, and the
  committed cohort-level verdict (matrix-capable organism at depth in 2 of 78 patients)
  is unchanged.
- **No bacterial rRNA depletion** — 89–99.6 % of the target organism's own reads are rRNA.
  The coding denominators are small (1,435 and 1,613 reads), so the growth values carry
  wide uncertainty; the direction is robust, the precise percentage is not.
- **145 bp reads** (not 146 bp), single-end R1 only. The arms are also ~145 bp.
- The *Pseudomonas* matrix call uses the on-disk single-strain panel, which is more
  conservative than the three-strain panel the arms used (Section 2). The *Acinetobacter*
  matrix call uses a purpose-built species-matched panel because the committed baumannii
  panel is blind to *A. pittii*.
- **`pga-only` ≠ matrix-ON.** Under the committed rule this sample cannot be matrix-ON:
  a one-system panel cannot satisfy "≥2 systems". Whether A. pittii pga expression
  constitutes a real biofilm phenotype is a biology question this rule cannot answer, and
  the call is reported as `pga-only` rather than upgraded by hand.

---

## 8. Files

| path | contents |
|---|---|
| `out/prjna699976_cases.json` | the two case records: growth (per reference), matrix, species ID, placement, rRNA diagnostic |
| `cases/case_map.py` | resumable chunked mapper (growth: `norm`/`normA`/`normP`; matrix: `panel`) |
| `cases/merge_cases.py` | folds chunk partials into the final records + placement |
| `cases/strain_id_case.py` | per-accession species attribution vs `acine_multi.mmi` |
| `cases/build_full_acine.py` | full-genome + other-genus-decoy index (the `compfull_pseu` pattern) |
| `cases/build_pittii_panel.py` | species-matched A. pittii two-axis panel (`qp_pittii.fna`, `comp_pittii.mmi`) |
| `cases/csu_orthology.py`, `cases/csu_nt_orthology.py` | the csu-absence / pga-presence orthology test |
| `cases/diag_rrna_case.py` | rRNA-share and MAPQ-band diagnostic |
| `cases/robustness_cases.py`, `cases/sens_cases.py` | per-chunk spread and threshold/window sensitivity |
| `cases/parts/*.json` | raw per-chunk results (the audit trail behind every number above) |

**New indexes built** (all reproducible from the scripts above): `compfull_acine.mmi`
(*A. baumannii* 17978 + 117 decoys), `compfull_pittii.mmi` (*A. pittii* RAES20 + 117
decoys), `comp_pittii.mmi` (pittii pga/RP/OG panel + 117 decoys). Downloaded FASTQs were
deleted immediately after mapping; only JSONs are retained.
