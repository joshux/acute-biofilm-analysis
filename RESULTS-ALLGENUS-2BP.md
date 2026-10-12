# All-genus two-axis test against adjudicated 2BP status — PRJNA1033689

**Date:** 2026-10-12. **Question (Josh's):** the adjudicated no-BP arm was screened for
*Pseudomonas only*. If secondary bacterial infection is being *missed* clinically in some
"no-BP" patients, those patients should still show Kolpen's fast × matrix state — just carried
by a different organism. So: do no-BP patients carry matrix-capable organisms that are
**actively growing with matrix machinery on**?

**Answer: yes, they carry the organisms — and one no-BP patient lands squarely in the
fast × matrix-ON quadrant. But the arm-level contrast does not reach significance, and the
test is still confounded by depth.** Detail below; this is a real but partial result, not a
vindication of the missed-infection hypothesis at cohort scale.

---

## Why the prior 2BP test could not answer this

The earlier 2BP analysis (RESULTS-2BP.md) screened every run against a **PAO1-only** index.
Its "no-BP matrix contrast is near-tautological — no no-BP run reaches the matrix depth gate"
conclusion was therefore an artifact of the panel, not of the biology: a no-BP patient whose
airway is dominated by *E. coli*, *Klebsiella*, *S. aureus* or *Acinetobacter* would have
scored **zero** on that screen by construction.

This run re-screens **all 254 runs** against the full 8-genus panel (Pseudomonas, Klebsiella,
Escherichia, Staphylococcus, Acinetobacter, Haemophilus, Streptococcus, Stenotrophomonas),
plus 34 decoy genera for MAPQ competition.

## Methods

- Toolchain rebuilt from committed coordinates in a fresh sandbox: minimap2 2.31, 8
  `gcomp_<org>.mmi` indexes (full genome + other-genus decoys), `panel_multi.mmi`
  (46 genomes), `gff_intervals.json` intervals, `acc2genus.json` (corrected labels).
- Labels: GEO `GSE246795` series matrix carries per-sample `2bp_status: 2BP/No-BP` for all
  56 GEO samples = **27 2BP / 29 no-BP patients**. GEO sample titles
  (`MVIR1-HS161-D10ETA1-MNGS1`) are embedded verbatim in the ENA BioSample aliases of the
  254 SRA runs (PRJNA1033689) → join is **exact and complete** (178 ETA + 76 NAS runs;
  99 ETA/2BP, 79 ETA/no-BP).
- Screening: every run subsampled to 40k reads, mapped to `panel_multi.mmi` at MAPQ ≥ 20.
- Quantification: qualifying runs subsampled to a **uniform 500k reads/mate** and scored with
  the Gifford-faithful per-organism two-axis metric —
  `%RP = RP reads / coding reads` (rRNA/tRNA excluded; FAST > 10, SLOW < 5) and
  `matrix-ON = ≥2 matrix subsystems each with ≥2 genes detected`, against the sample's own
  dominant organism's full-genome index.

## Result 1 — the no-BP arm is full of matrix-capable organisms

Deepest run per patient, best axis organism at 40k screening depth:

| arm | patients with a matrix-capable dominant at depth | examples |
|---|---|---|
| 2BP | Pseudomonas 6, Escherichia 8, Klebsiella 3 | HS236 Pseudo 57k, HS231 Esch 8.8k |
| **no-BP** | **Escherichia 15, Klebsiella 3, Pseudomonas 2** | HS243 Esch 7.2k, HS253 Kleb 5.9k, HS94 Kleb 1.7k |

The no-BP arm is **not** organism-poor. It is dominated by *Escherichia* and *Klebsiella*,
neither of which the PAO1-only screen could see.

## Result 2 — two-axis quantification (30 decisive runs, uniform 500k depth)

Usable runs = matrix-capable dominant organism with ≥300 target reads (n=20); depth-matched
subset at ≥5000 target reads (n=16):

| arm | FAST × matrix-ON | FAST × off | SLOW/ambig × ON | SLOW/ambig × off |
|---|---|---|---|---|
| **2BP** (n=12) | **6** | 3 | 1 | 2 |
| **No-BP** (n=8) | **1** | 3 | 0 | 4 |

Depth-matched (≥5000 reads): 2BP 6/11 FAST×ON vs no-BP 1/5.

Fisher exact: FAST×matrix-ON **p = 0.16**; matrix-ON alone **p = 0.070**; FAST alone p = 0.36.

### The individual no-BP placements matter more than the arm tally

| patient | organism | %RP | state | matrix | note |
|---|---|---|---|---|---|
| **HS243** | *E. coli* | **16.9** | **FAST** | **ON (bcs, csg)** | lands in the acute quadrant |
| HS94 | *Klebsiella* | 23.1 | FAST | off (mrk only) | fast growth, matrix panel too thin |
| HS523 | *E. coli* | 19.4 | FAST | off | fast growth, no matrix detected |
| HS185 | *E. coli* | 10.9 | FAST | off | borderline |
| HS253 | *Klebsiella* | 5.5 | ambig | mrk only | |
| HS78/55/1 | *E. coli* | ~0 | SLOW | off | target reads < 2000 — depth-limited |

**HS243 is the finding.** A patient three physicians classified as *no* secondary bacterial
pneumonia carries a *fast-growing E. coli* population transcribing two matrix systems. That is
exactly the state the One-Two Punch predicts for a landed second punch.

## Result 3 — the honest limits

1. **The arm-level contrast is not significant** (p = 0.16 for the quadrant; 0.070 for
   matrix-ON). The direction is right — 50% of 2BP vs 12% of no-BP runs are FAST×ON — but with
   12 vs 8 usable runs this cohort cannot separate the arms at the quadrant level.
2. **Depth still confounds the no-BP arm.** Several no-BP placements (HS78, HS55, HS1) sit at
   <2000 target reads where %RP is unstable. The Klebsiella matrix panel is also thin
   (mrk only, 2 genes) so Klebsiella-dominated samples (HS94, HS253) cannot reach matrix-ON by
   construction — the same tautology the PAO1-only screen suffered, now for a different genus.
3. **E. coli's bcs/csg systems are not equivalent to Pseudomonas alginate/psl/pel.** bcs and
   csg are cellulose and curli fimbriae; a "matrix-ON" call in *E. coli* does not mean the same
   extracellular-matrix programme as a Pseudomonas alginate/psl/pel call. Cross-genus matrix
   calls are not commensurate, and the 2BP arm's matrix-ON calls are mostly Pseudomonas while
   the no-BP arm's are E. coli.
4. **Colonisation is not excluded.** *E. coli* / *Klebsiella* in tracheal aspirate at 10–20%
   of the community can be colonisation, and this assay cannot tell colonisation from infection.
   The prior thread's finding stands: microbial burden measures are flat/non-specific.

## Verdict for the theory

- **The missed-infection hypothesis is not refuted, and it gains one concrete case (HS243).**
  The no-BP arm was previously declared un-scoreable; it is scoreable, and it is not uniformly
  matrix-off.
- **It is not confirmed at cohort scale either.** The two-axis *state* remains a weak
  discriminator between adjudicated 2BP and no-BP (p = 0.16), consistent with — not a challenge
  to — the earlier conclusion that the discriminative claim is untested by this cohort.
- **The methodological point is the durable one:** an organism-specific screen (PAO1-only)
  produced a *false* "no comparator arm" conclusion. Pan-genus screening overturns it. This is a
  second instance of the paper's central artifact class — a reference-panel choice silently
  determining the result.

## Reproduce

```
screen_nobp.py   # 40k-read pan-genus screen of all 254 runs -> screen_out/
sub_worker.sh    # dl_full.sh + sub500k.sh + quant_sub.py over /tmp/decide_srrs.txt
quant_sub.py     # per-organism two-axis quantification at 500k depth -> final_sub/
```
