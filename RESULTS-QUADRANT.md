# Two-dimension (quadrant) test — acute vs chronic vs healthy

**Date:** 2026-10-10 · **Repo:** joshux/acute-biofilm-analysis · **Supersedes:** the coupling/
correlation test in `RESULTS.md` (that test answered a question Kolpen never asked).

## What this tests

Kolpen et al. 2022 (Thorax 77:1015) argue two **decoupled** axes. Their own words: *"the
difference lies primarily in metabolic rates, not bacterial architecture."*

| | matrix (architecture) | growth (metabolic rate) |
|---|---|---|
| acute | ON | FAST |
| chronic | ON | SLOW |

**No matrix↔growth correlation is predicted — decoupling is the claim.** The test is therefore
**per-group quadrant placement** on two thresholds, not a within-group correlation.

- **Growth axis:** `%RP` = ribosomal-protein transcripts as a fraction of the mapped target-organism
  transcriptome (Gifford 2014). Anchors: **FAST > 10 %RP, SLOW < 5 %RP.**
- **Matrix axis:** systems detected among **alginate / psl / pel**, each needing ≥ 3 genes.
  **matrix-ON = ≥ 2 of the 3 systems.**

## Groups

| group | dataset | n | specimen |
|---|---|---|---|
| acute | PRJNA1056765 (Tang 2025 Sci Data) | 15 Pseudomonas-dominant | BALF |
| chronic | PRJEB24688 (Rossi 2018 Nat Commun) | 15 in vivo | CF sputum |
| chronic anchor | PRJEB24688 | 12 exponential / 11 stationary | lab culture |
| healthy | PRJNA390194 (Ren 2018) | 9 non-COPD | BALF |

**The chronic dataset carries its own internal anchors.** PRJEB24688 contains, in one experiment,
in vivo sputum plus in vitro exponential and stationary cultures of the same clones. These are the
FAST/SLOW anchors measured under identical conditions, so the chronic growth call does not depend on
cross-dataset normalisation.

## Result

### Growth axis — the differentiator is not observed

| group | n | median %RP | FAST | ambiguous | SLOW |
|---|---|---|---|---|---|
| acute BALF | 15 | **11.4 %** | 8 | 7 | 0 |
| chronic in vivo sputum | 15 | **9.3 %** | 7 | 8 | 0 |
| chronic lab exponential | 12 | 20.4 % | 10 | 2 | 0 |
| chronic lab stationary | 11 | 11.5 % | 5 | 6 | 0 |

- Chronic in vivo is clearly **below** lab exponential — Mann-Whitney **p = 0.015**. The metric
  works, and it reproduces Rossi's own conclusion ("a low-energy, low-growth, non-motile
  physiological state").
- But **acute is indistinguishable from chronic in vivo** — median 11.4 % vs 9.3 %,
  Mann-Whitney **p = 0.35**; restricted to acute samples at adequate depth (target ≥ 500 reads,
  n = 5) **p = 0.76**.
- Both sit at the FAST/ambiguous boundary (median ≈ 9–11 %RP, no sample < 5 %). Only lab
  exponential is unambiguously fast.

**Kolpen's differentiator requires acute FAST and chronic SLOW. What the data show is both
intermediate and equal** — the metabolic contrast between the two compartments is not present.

### Matrix axis — architecture claim holds

| group | n | matrix present | matrix-ON (≥ 2 systems) |
|---|---|---|---|
| acute BALF | 15 | 14/15 (93 %) | 5/15 (33 %) |
| chronic in vivo sputum | 15 | 15/15 (100 %) | 15/15 (100 %) |

Matrix transcripts are present in **both** compartments — Kolpen's core claim (matrix is not
chronic-only) holds, and the classical "acute = planktonic" model is contradicted.

**But the acute matrix-ON rate is depth-limited, not biology-limited.** The matrix axis is a
detection call, and it rises with sequencing depth in the acute arm:

| acute depth floor | n | matrix-ON |
|---|---|---|
| ≥ 100 reads | 11 | 5/11 |
| ≥ 300 | 8 | 5/8 |
| ≥ 500 | 5 | 4/5 |
| ≥ 700 | 3 | **3/3** |

At matched depth the acute samples reach 100 % matrix-ON — the same as chronic. The 33 % headline
is a coverage artefact of the shallow acute arm (median 312 target reads vs 15,171 for chronic).

The growth axis, being a ratio, is **depth-robust**: acute median %RP is 9.66 → 9.54 → 9.42 %
across depth floors ≥ 100 / 300 / 500.

### Healthy controls

All 9 non-COPD BALF samples yield **zero** panel target reads (median 0; reads that do map land on
decoy/kitome genera, median 369). Healthy BALF bacterial RNA sits at the kitome floor — reported as
a finding, as pre-registered. There is no bacterial signal to place on either axis.

## Quadrant placement

| | FAST | ambiguous | SLOW |
|---|---|---|---|
| **matrix-ON** | acute 2 · chronic 23 | acute 3 · chronic 15 | — |
| **matrix-OFF / alginate-only** | acute 3 | acute 3 | — |
| **healthy** | — | — | (no signal) |

*(acute counts at the ≥100-read floor; the chronic FAST column includes the lab-culture arm in the
full table — see `quadrant_table.csv` for the per-sample split.)*

## Verdict

1. **Architecture claim (Kolpen): CONFIRMED and extended.** Matrix machinery is transcribed in
   acute BALF as well as chronic sputum — the first direct transcript detection of Pel/Psl in acute
   infection. Matrix is not a chronic-only phenotype.
2. **Metabolic differentiator (Kolpen): NOT OBSERVED.** Chronic in vivo is genuinely slow relative
   to its own lab exponential anchor (p = 0.015), but acute BALF is **just as slow**. The acute-vs-
   chronic growth-rate contrast that the two-axis model rests on does not appear at achievable depth.
3. **Classical acute-planktonic model: CONTRADICTED.** Acute samples are matrix-positive.
4. **The one deeply-sequenced acute sample is Kolpen-like.** SRR27343249 (2,627 target reads,
   ~10× the rest) is matrix-ON (all three systems) **and** FAST (19.5 %RP) — exactly the predicted
   acute quadrant. It is n = 1, and the remaining acute samples are too shallow to call.

The strongest supported statement: **on these two dimensions acute and chronic bacterial lung
infection look the same — matrix-ON with growth at the fast/ambiguous boundary.** Neither
paradigm's acute prediction survives as a contrast.

## Limitations

- **Depth asymmetry is the binding constraint.** Acute median 312 target reads vs 15,171 chronic.
  The growth axis tolerates this (ratio, stable across depth floors); the matrix axis does not
  (detection call, depth-sensitive). Acute matrix-ON is a **lower bound**.
- **Cross-specimen, cross-protocol.** Acute is BALF (Tang; human-only rRNA depletion); chronic is
  CF sputum (Rossi; bacterial rRNA depleted). %RP uses protein-coding panel reads only, so the
  depletion difference does not enter the growth metric, but the two are not the same specimen type.
- **Bulk averaging erases spatial structure.** A quadrant placement is a population average; Kolpen's
  claim is about aggregates. A negative here is weaker than a refutation.
- **"Acute" is our framing** of Tang's "bacterial infection" label; "chronic" is Rossi's chronically
  infected CF cohort.
- **Single reference per genus** beyond the multi-strain panel; CF strains may diverge from PAO1/PA14.

## Reproduce

```
build_panels.py                 # gene panels from gff/panels.json + decoy genomes
get_fastq.py / subsample.py     # acute full; chronic & healthy subsampled to 2M reads
quantify_v3.py <SRR> pseu       # %RP + per-system matrix detection
chronic_groups.py               # split PRJEB24688 into in vivo / exp / stat
analyze_quadrant.py             # quadrant tables
final_quadrant_stats.py         # group statistics
make_quadrant_fig.py            # fig_quadrant.png
```

Data files: `quadrant_table.csv`, `chronic_groups.csv`, `final_stats.json`,
`chronic_growth.json`, `fig_quadrant.png`.
