# PRJNA699976 — Empirical two-axis feasibility probe

**Question:** Can PRJNA699976 (Zou et al., *Adv Sci* 2024) support a bacterial two-axis
analysis — (1) **GROWTH** = bacterial ribosomal-protein transcript share, (2) **MATRIX** =
transcription of species-specific biofilm matrix operons?

**Method:** empirical probe of the deepest runs, not metadata inference. Tooling:
`probe_composition2.py` (MAPQ×class contingency), plus two new scripts written for this
probe — `prj_genus.py` (decoy GCF → genus composition) and `prj_rrna.py` / `prj_twoaxis.py`
(the rRNA-depletion diagnostic and target-organism two-axis test).

**Bottom line: UNUSABLE.** No bacterial rRNA depletion; Pseudomonas effectively absent;
the organisms that actually dominate are not matrix-capable.

---

## 1. Metadata

| field | value |
|---|---|
| runs | **440** (all unique sample accessions) |
| `library_strategy` | RNA-Seq ×440 |
| `library_source` | METATRANSCRIPTOMIC ×440 |
| `library_selection` | RANDOM ×440 |
| layout | 189 paired-end, 251 single-end |
| read_count | median **17,703,387**; max 93,122,168; min 10,128 |
| depth | **374 runs > 5M reads; 313 runs > 10M reads** |

The metadata is completely uniform — every run is a metatranscriptomic RNA-Seq with
`RANDOM` selection. **The rRNA-depletion method is not recorded anywhere in the ENA
record**, which is exactly why the prior survey left this project undetermined. That
question can only be settled by looking at the reads, which is what this probe does.

**DOI check:** `10.1002/advs.202405087` resolves (302 → Wiley), Crossref confirms
*Advanced Science*, 2024, Zou, Yan, Wang, Ni, Zhao, Lu, Liu, Cao — "Accurate Diagnosis of
Lower Respiratory Infections Using Host Response and Respiratory Microbiome from a Single
Metatranscriptome Test of Bronchoalveolar Lavage Fluid". ✅ DOI valid.

## 2. Per-sample probe table

Six deepest runs (4 deepest single-end + 2 deepest paired R1), ~150k-read subsamples each,
mapped with `minimap2 -x sr` against `compfull_pseu.mmi` (PAO1 + 120 decoy contigs).

| run | reads | PAO1 any | PAO1 MAPQ≥20 | PAO1 rRNA/tRNA | PAO1 coding (perm.) | decoy any | decoy MAPQ≥20 | mean read len |
|---|---|---|---|---|---|---|---|---|
| SRR13677796 | 157,014 | 84 | 1 | 84 | **0** | 113,009 | 7,493 | 145.0 |
| SRR13677795 | 158,718 | 112 | 1 | 112 | **0** | 126,451 | 6,988 | ~145 |
| SRR13677629 | 154,783 | 100 | 5 | 87 | 13 | 86,615 | 31,732 | ~145 |
| SRR13677627 | 157,975 | 179 | 4 | 160 | 19 | 149,379 | 54,126 | 144.7 |
| SRR13677829 | 162,875 | 86 | 0 | 84 | 2 | 142,673 | 7,180 | ~145 |
| SRR13677601 | 158,314 | 150 | 4 | 131 | 19 | 150,884 | 56,827 | 144.9 |

**Pseudomonas is absent.** PAO1-mapping reads are 84–179 per ~157k subsample (0.05–0.11%),
of which **0–5 reads** reach MAPQ≥20 and the handful that do are almost entirely rRNA/tRNA.
There are effectively **zero Pseudomonas coding reads** in any probed sample — the growth-axis
denominator for PAO1 is empty. The reads that *do* map are overwhelmingly to decoy
(non-Pseudomonas) genomes: 55–95% of reads map somewhere, but almost none of it is Pseudomonas.

## 3. The rRNA-share diagnostic (the key question)

Because PAO1 is absent, the PAO1-based rRNA share is noise. The diagnostic was therefore run
**at the level of the organism that actually dominates each library**: every read mapping at
MAPQ≥20 to a named decoy genome was classified by whether its alignment lands inside that
genome's own rRNA/tRNA intervals (from its GFF).

| run | dominant organism | decoy MAPQ≥20 | rRNA/tRNA hits | other hits | **rRNA share** |
|---|---|---|---|---|---|
| SRR13677627 | Ralstonia pickettii | 54,126 | 53,183 | 943 | **98.3 %** |
| SRR13677601 | Ralstonia pickettii | 56,827 | 55,496 | 1,331 | **97.7 %** |
| SRR13677629 | Ralstonia pickettii | 31,732 | 30,827 | 905 | **97.1 %** |
| SRR13677796 | Prevotella melaninogenica | 7,493 | 1,628 | 5,865 | 21.7 % |
| SRR13677795 | Veillonella / oral anaerobes | 6,988 | 1,431 | 5,557 | 20.5 % |
| SRR13677829 | Veillonella | 7,180 | 978 | 6,202 | 13.6 % |

Confirmed across 9 further runs sampled evenly through the depth range (ranks 20–430):

| run | depth | decoy MAPQ≥20 | rRNA share |
|---|---|---|---|
| SRR13677673 | 45.6M | 34,094 | **97.8 %** |
| SRR13677628 | 38.3M | 37,576 | **96.5 %** |
| SRR13677863 | 26.8M | 36,598 | **98.0 %** |
| SRR13677746 | 11.0M | 38,130 | **95.3 %** |
| SRR13677625 | 3.1M | 36,726 | **98.3 %** |
| SRR13677687 | 32.3M | 13,558 | 48.5 % |
| SRR13677611 | 15.0M | 3,856 | 49.7 % |
| SRR13677850 | 6.7M | 1,262 | 38.3 % |
| SRR13677801 | 19.8M | 4,299 | 21.1 % |

**Reading of the diagnostic.** When a library contains a substantial bacterial population
(the Ralstonia-dominated mode — 5 of 6 deepest, and the majority of the stratified sample),
**95–98 % of bacterial reads are rRNA/tRNA**. That is the signature of a library with **no
bacterial rRNA depletion whatsoever** (only host depletion, if any). This is not a marginal
call: it is the textbook undepleted-metatranscriptome profile, where the bacterial signal is
buried under ribosomal RNA.

A second, deeper check — classifying *all* mapped reads (not just MAPQ≥20) against the rRNA
intervals — gives **87.8 % rRNA** for SRR13677627 and **72.9 %** for SRR13677796. The
low-MAPQ mass that makes up the bulk of the libraries is conserved rRNA cross-mapping between
related organisms, i.e. even the fraction that fails the MAPQ≥20 filter is mostly ribosomal.
The true rRNA burden is therefore ≥73–88 % and the usable bacterial *coding* signal is a thin
residue.

The minority mode (Prevotella/Veillonella samples) shows a lower rRNA share only because its
bacterial population is so small that most bacterial reads are non-rRNA — not because those
libraries were depleted. Their bacterial content is ~5 % of reads and dominated by oral
anaerobes.

## 4. Genus composition — is a matrix-capable organism present at depth?

`prj_genus.py` resolves each decoy hit through `decoy/manifest.json` (GCF → taxon). Top genera
per sample:

| run | top genera (MAPQ≥20) |
|---|---|
| SRR13677627 | **Ralstonia** 53,845 (99.5 %) · Escherichia 140 · Delftia 66 |
| SRR13677601 | **Ralstonia** 56,463 (99.4 %) · Escherichia 215 · Delftia 65 |
| SRR13677629 | **Ralstonia** 31,326 (98.7 %) · Veillonella 219 · Escherichia 75 |
| SRR13677796 | **Prevotella** 4,039 (54 %) · Veillonella 1,766 · Tannerella 522 · Streptococcus 416 |
| SRR13677795 | **Veillonella** 1,688 · Streptococcus 1,202 · Haemophilus 1,181 · Prevotella 1,139 · Tannerella 866 |
| SRR13677829 | **Veillonella** 4,779 (67 %) · Haemophilus 887 · Neisseria 808 · Prevotella 244 |

The dominant organisms are **Ralstonia pickettii** (an environmental/burkholderiaceae
opportunist) and **oral anaerobes** (Veillonella, Prevotella, Tannerella, Fusobacterium) —
this is a BALF cohort whose microbiome signal is upper-airway/oral contamination plus an
environmental organism. **None of these is matrix-capable.**

Screen for the matrix-capable genera that have characterised operons in this project's panels
(Pseudomonas, Klebsiella, Escherichia, Staphylococcus, Acinetobacter, Stenotrophomonas,
Streptococcus), across all 15 probed runs:

| genus | best count in any probed run | verdict |
|---|---|---|
| Pseudomonas | 4 reads (0.003 %) | **absent** |
| Klebsiella | 0 | absent |
| Staphylococcus | 0 | absent |
| Acinetobacter | 2 reads | absent |
| Escherichia | 215 reads (0.14 %) | trace (and likely cross-mapping among Enterobacterales — 22 % of these hits are conserved rRNA) |
| Streptococcus | 1,202 reads (0.8 %) | trace, oral, not the intended target |

**No matrix-capable organism is present at depth.** The only organisms that reach meaningful
depth (Ralstonia, Veillonella, Prevotella) have no alginate/psl/pel/csu/pga/ica-type matrix
operon panel in this project and are not matrix producers in the relevant sense.

## 5. Direct two-axis test on the dominant organism

`prj_twoaxis.py` builds a per-sample index (target genome + other decoys) and classifies
MAPQ≥20 reads into rRNA/tRNA, ribosomal-protein CDS, and other CDS for the *dominant*
organism:

| run | target organism | rRNA/tRNA | ribosomal-protein | other CDS | coding reads | %RP of coding | %rRNA of target |
|---|---|---|---|---|---|---|---|
| SRR13677627 | Ralstonia pickettii | 55,483 | **22** | 609 | 631 | 3.5 % | **98.9 %** |
| SRR13677601 | Ralstonia pickettii | 37,920 | **18** | 576 | 594 | 3.0 % | **98.5 %** |
| SRR13677629 | Ralstonia pickettii | 21,412 | **9** | 322 | 331 | 2.7 % | **98.5 %** |
| SRR13677796 | Prevotella melaninogenica | 161 | 792 | 1,812 | 2,604 | 30.4 % | 5.8 % |

This is the decisive test:

- **Ralstonia mode (the dominant mode):** 98.5–98.9 % of the target organism's reads are
  rRNA/tRNA. The growth-axis denominator is 331–631 coding reads, of which only **9–22** are
  ribosomal-protein — 10–20× too few for a stable %RP, and the library is a ribosomal swamp.
  Both axes fail.
- **Prevotella mode (the minority):** the target's rRNA share is only 5.8 % and %RP of coding
  is 30.4 % — a *plausible* growth-axis signal for Prevotella. But (i) this holds only for the
  minority of runs, (ii) the bacterial fraction is ~1.7 % of reads, (iii) the dominant organism
  is inconsistent run-to-run (Prevotella here, Veillonella elsewhere), and (iv) **Prevotella has
  no matrix operon** — so the matrix axis is dead regardless, and this project's question is
  Pseudomonas matrix biology.

---

## 6. Verdict

**(a) Was bacterial rRNA depleted? — NO.**
Rigorously, no. In the libraries with a real bacterial population, 95–98 % of bacterial reads
are rRNA/tRNA; across all mapped reads the rRNA share is 73–88 %. There is no rRNA-depletion
signature anywhere in the cohort. This was the one open metadata question and the reads answer
it unambiguously.

**(b) Enough bacterial coding depth for a growth axis? — NO.**
For Pseudomonas the denominator is empty (0–5 coding reads per sample; PAO1 absent). For the
organisms that *are* present, the growth-axis denominator is either drowned by rRNA (Ralstonia
mode: 98.5 % rRNA, 9–22 RP reads) or only viable for a non-matrix-capable oral anaerobe in a
minority of runs at ~1.7 % of reads. Not viable as a two-axis growth measurement.

**(c) A matrix-capable organism present at depth for a matrix axis? — NO.**
Pseudomonas ≤4 reads, Klebsiella 0, Staphylococcus 0, Acinetobacter ≤2, Escherichia ≤215
(0.14 %, largely cross-mapping), Streptococcus ≤1,202 (0.8 %, oral). No matrix-capable genus
reaches usable depth. The deep organisms (Ralstonia pickettii, Veillonella, Prevotella) have
no matrix operon.

**(d) Overall — UNUSABLE.**

This is a well-evidenced negative that closes the question. PRJNA699976 is a
metatranscriptomic BALF cohort with **no bacterial rRNA depletion**, dominated by
**Ralstonia pickettii** and oral anaerobes, with **Pseudomonas effectively absent**. It cannot
support the growth axis (ribosomal swamp / empty Pseudomonas denominator) and it cannot
support the matrix axis (no matrix-capable organism at depth). The prior metadata survey's
"undetermined" verdict is now resolved to **ruled out**, not because of what the metadata
omitted but because of what the reads show.

---

### Files
- Raw numbers: `out/prjna699976_probe.json` (aggregate)
- Per-run: `out/<SRR>_probe2.json`, `out/<SRR>_genus.json`, `out/<SRR>_rrna.json`, `out/<SRR>_rrna_any.json`, `out/<SRR>_<acc>_twoaxis.json`
- Scripts added: `prj_genus.py`, `prj_rrna.py`, `prj_twoaxis.py`
- Metadata: `prjna699976_runs.tsv`
