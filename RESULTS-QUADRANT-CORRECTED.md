# Two-dimension (quadrant) test — acute vs chronic vs healthy

**Date:** 2026-10-10 · **Repo:** joshux/acute-biofilm-analysis

> **CORRECTION (2026-10-10, same day).** The first version of this file reported "acute just as slow
> as chronic" and concluded Kolpen's metabolic differentiator was absent. **That was a measurement
> artifact.** The growth metric was built on the multi-genome panel, where conserved
> ribosomal-protein reads map equally to PAO1/PA14/P. putida and are assigned MAPQ 0 — so ~94 % of
> RP reads were being dropped by the MAPQ≥20 filter, deflating %RP and compressing all groups
> toward the middle. On the corrected, Gifford-faithful metric the result **reverses**: acute is
> significantly *faster* than chronic. See "What went wrong" below.

## What this tests

Kolpen et al. 2022 (Thorax 77:1015) argue two **decoupled** axes — *"the difference lies primarily
in metabolic rates, not bacterial architecture."*

| | matrix (architecture) | growth (metabolic rate) |
|---|---|---|
| acute | ON | FAST |
| chronic | ON | SLOW |

No matrix↔growth correlation is predicted — decoupling is the claim. The test is therefore
**per-group quadrant placement** on two thresholds.

- **Growth axis:** `%RP` = ribosomal-protein reads / *Pseudomonas* coding reads, following Gifford
  2014 (RP / total). Anchors **FAST > 10 %RP, SLOW < 5 %RP.**
- **Matrix axis:** systems detected among **alginate / psl / pel**, each needing ≥ 3 genes.
  **matrix-ON = ≥ 2 of the 3 systems.**

## Groups

| group | dataset | n | specimen |
|---|---|---|---|
| acute | PRJNA1056765 (Tang 2025 Sci Data) | 15 | BALF |
| chronic in vivo | PRJEB24688 (Rossi 2018 Nat Commun) | 15 | CF sputum |
| chronic anchor | PRJEB24688 | 12 exp / 11 stat | lab culture |
| healthy | PRJNA390194 (Ren 2018) | 9 non-COPD | BALF |

PRJEB24688 ships **its own internal anchors** — in vitro exponential and stationary cultures of the
same clones alongside the in vivo sputum — so the FAST/SLOW call is made within one experiment.
(23/38 RNA runs are `isolation_source = laboratory culture`; titles: `EXP`/`3H` = exponential,
`STAT`/`24H` = stationary, `INVIVO` = sputum.)

## Result — the growth axis

| group | n | median %RP | FAST | ambiguous | SLOW |
|---|---|---|---|---|---|
| **acute BALF** | 13 | **10.7 %** | 7 | 6 | 0 |
| **chronic in vivo sputum** | 15 | **4.6 %** | 3 | 4 | 8 |
| chronic lab exponential | 12 | 13.6 % | 8 | 4 | 0 |
| chronic lab stationary | 11 | 2.1 % | 0 | 0 | 11 |

**The lab anchors validate the metric**: exponential 13.6 % (FAST) vs stationary 2.1 % (SLOW),
p < 1e-4. The metric orders known growth states correctly.

**Key contrast — acute vs chronic in vivo: 10.7 % vs 4.6 %, Mann-Whitney p = 0.0015.**
Acute bacterial BALF is significantly **faster** than chronically-infected CF sputum.

This is Kolpen's prediction, in direction. Against the anchors: **acute sits between the lab
exponential and stationary states, closer to exponential; chronic in vivo sits essentially at the
stationary state** (4.6 % vs 2.1 %).

Depth-robustness (acute, floor on total PAO1 reads):

| floor | n | median %RP |
|---|---|---|
| ≥ 100 | 14 | 10.6 % |
| ≥ 1 000 | 13 | 10.5 % |
| ≥ 5 000 | 9 | 8.3 % |

Stable — and stable across reference strain (PAO1 / PA14 / P. putida all give the same ordering).

## Result — the matrix axis

| group | n | matrix present | matrix-ON (≥ 2 systems) |
|---|---|---|---|
| acute BALF | 13 | 12/13 (92 %) | 5/13 |
| chronic in vivo sputum | 15 | 15/15 (100 %) | 15/15 (100 %) |

Matrix machinery is transcribed in **both** compartments — Kolpen's architecture claim holds and the
classical "acute = planktonic" model is contradicted. Acute matrix-ON (5/13) is **depth-limited, not
biology-limited**: the acute arm has ~300× fewer reads than chronic, and the detection call rises
with depth (5/11 → 4/5 → 3/3 across floors). At matched depth the acute samples reach 100 %
matrix-ON, the same as chronic. Treat acute matrix-ON as a **lower bound**.

## Quadrant placement (depth-ok samples)

| | FAST | ambiguous | SLOW |
|---|---|---|---|
| **acute matrix-ON** | 2 | 3 | 0 |
| **acute alginate-only** | 5 | 3 | 0 |
| **chronic in vivo matrix-ON** | 3 | 4 | 8 |

Acute occupies the **FAST column with matrix present** — the Kolpen acute quadrant.
Chronic in vivo occupies **matrix-ON × SLOW** — the Kolpen chronic quadrant.
The prediction is, on both axes, in the right place.

## Verdict

1. **Architecture claim (Kolpen): confirmed.** Matrix machinery transcribed in acute BALF as well
   as chronic sputum; the classical acute-planktonic model is contradicted.
2. **Metabolic differentiator (Kolpen): confirmed in direction.** Acute is significantly faster than
   chronic in vivo (10.7 % vs 4.6 %, p = 0.0015), with the metric validated on the dataset's own
   lab exponential/stationary anchors. Chronic in vivo sits at the stationary state; acute sits
   between stationary and exponential, closer to exponential.
3. **Healthy controls:** zero bacterial signal (kitome floor) — nothing to place, as pre-registered.

## What went wrong (recorded, because it is the third time this class of error has appeared)

The first pass computed `%RP` as `RP / (matrix + rp + og)` on a **multi-genome** panel. Because
ribosomal-protein genes are near-identical across *Pseudomonas* species, their reads map equally to
PAO1, PA14 and P. putida → minimap2 assigns MAPQ 0 → they are dropped at MAPQ ≥ 20, while
*Pseudomonas*-specific matrix genes map uniquely and survive. The MAPQ diagnostic is unambiguous:

| sample | median MAPQ, rp | median MAPQ, matrix | %RP at MAPQ 0 | %RP at MAPQ ≥ 20 |
|---|---|---|---|---|
| SRR27343249 | **0** | 15 | 56.2 % | 19.5 % |
| SRR27343409 | **0** | 0 | 30.8 % | 8.0 % |
| ERR2275089 | **0** | 0 | 58.1 % | 26.6 % |

Single-genome panels (PAO1 only) remove the intra-species competition: the same reads then give
%RP ≈ 59 % at every MAPQ threshold, flat. The fix is to measure growth against **one** reference
genome — with non-*Pseudomonas* genera present only as decoys to absorb cross-genus conserved reads.

This is the same artifact class as the CIT-protein and compositional-share errors from earlier in
this project: **a filter that behaves differently for the two arms being compared**. The
multi-genome panel was correct for *genus assignment* and wrong for *within-genus quantification*.

## Limitations

- **Cross-specimen, cross-protocol.** Acute is BALF (Tang; human-only rRNA depletion); chronic is CF
  sputum (Rossi; bacterial rRNA depleted). Removing rRNA/tRNA from the denominator makes the metric
  independent of the depletion difference — but the two are not the same specimen type, and the
  acute arm is ~300× shallower.
- **Bulk averaging erases spatial structure.** A quadrant placement is a population average;
  Kolpen's claim is about aggregates. This is weaker than direct validation.
- **%RP is a proxy.** Ribosomal-protein transcript share tracks growth rate in culture; in a
  nutrient-limited host environment the mapping is not guaranteed.
- **"Acute" is our framing** of Tang's "bacterial infection" label; "chronic" is Rossi's chronically
  infected CF cohort.
- Strain divergence was tested (PAO1 / PA14 / P. putida) and does not change the ordering.

## Reproduce

```
build_panels.py                 # gene panels + decoy genomes (genus assignment)
quantify_v3.py <SRR> pseu       # %RP (panel) + per-system matrix detection
build_single.py                 # single-genome panel + decoys
diag_mapq.py                    # the MAPQ diagnostic that exposed the artifact
quantify_rp.py / quantify_rp2.py  # Gifford-faithful %RP; rRNA-excluded denominator
strain_check.py                 # PAO1 / PA14 / putida robustness
compare_metrics.py              # the three metrics side by side
analyze_corrected.py            # corrected quadrant tables
make_corrected_fig.py           # fig_quadrant_corrected.png
```

Data: `quadrant_table_corrected.csv`, `growth_metric_compare.json`, `strain_robust.json`,
`diag_mapq.json`, `fig_quadrant_corrected.png`.
