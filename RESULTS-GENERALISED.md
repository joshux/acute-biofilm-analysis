# Generalised two-axis test — all viable organisms (2026-10-10)

Supersedes the Pseudomonas-only framing. Every sample is scored against **the
organism it actually contains**, using a per-organism matrix panel.

## Method

For each sample:
1. **Dominant organism** — multi-genome panel (19 genomes), MAPQ≥20.
2. **Growth axis** `%RP` = RP reads / coding reads (rRNA/tRNA excluded), mapped to
   **that organism's own full genome** with other genera present only as decoys
   (so the target genome keeps MAPQ — the artifact that broke the first pass).
3. **Matrix axis** — per-organism matrix subsystems, ON = ≥2 subsystems each with
   ≥2 genes detected.

Matrix panels are genus-specific, from the organism inventory (`ORGANISM-VIABILITY.md`):
Pseudomonas (alginate/psl/pel/cdr), E. coli (pga/csg/bcs), Klebsiella (mrk/fim/cps),
S. aureus (ica), Acinetobacter (pga/csu), H. influenzae (hmw/adhesin/lic),
S. pneumoniae (cps/cbp), S. maltophilia (fimbriae/qs).

## Cohort coverage

84 samples scored (16 acute + 38 chronic CF + 9 healthy + others). Dominant organism:

| organism | n | usable (coding ≥ 1000) |
|---|---|---|
| Pseudomonas | 53 | **51** |
| Acinetobacter | 20 | 1 |
| Escherichia | 7 | 0 |
| Staphylococcus | 2 | 0 |
| Klebsiella | 1 | 1 |
| Haemophilus | 1 | 0 |

**Only Pseudomonas reaches usable depth.** The others map at trace levels (the
cohort's non-Pseudomonas signal is mostly kitome), so the matrix axis is in
practice Pseudomonas-only — confirmed, not assumed.

## Result — Pseudomonas, all subgroups, one metric

| group | n | median %RP | range | FAST | ambig | SLOW | matrix-ON |
|---|---|---|---|---|---|---|---|
| **acute BALF** | 13 | **10.7 %** | 5.5–21.4 | 7 | 6 | 0 | 8/13 |
| **chronic in vivo** | 15 | **4.6 %** | 2.3–19.8 | 3 | 4 | 8 | 15/15 |
| lab exponential | 12 | 13.6 % | 6.7–15.8 | 8 | 4 | 0 | 12/12 |
| lab stationary | 11 | 2.1 % | 1.0–3.8 | 0 | 0 | 11 | 11/11 |

| contrast | medians | p |
|---|---|---|
| acute vs chronic in vivo | 10.7 % vs 4.6 % | **0.0034** |
| chronic in vivo vs lab exponential | 4.6 % vs 13.6 % | 0.0025 |
| lab exponential vs lab stationary | 13.6 % vs 2.1 % | <1e-4 |

**The lab anchors validate the metric** (exponential ≫ stationary). Acute sits
between stationary and exponential, closer to exponential; chronic in vivo sits
at the stationary state.

## Quadrant placement

| acute (n=13) | FAST | ambig | SLOW |
|---|---|---|---|
| matrix-ON | 3 | 5 | 0 |
| matrix-off | 4 | 1 | 0 |

| chronic in vivo (n=15) | FAST | ambig | SLOW |
|---|---|---|---|
| matrix-ON | 3 | 4 | 8 |
| matrix-off | 0 | 0 | 0 |

Acute occupies **matrix-present × FAST/ambiguous**; chronic occupies
**matrix-ON × SLOW**. Kolpen's two predictions land in the right quadrants.

## Verdict

1. **Architecture (Kolpen): confirmed.** Matrix machinery is transcribed in acute
   as well as chronic — matrix is not a chronic-only phenotype.
2. **Metabolic differentiator (Kolpen): confirmed in direction.** Acute is
   significantly faster than chronic in vivo (p=0.0034), with the metric validated
   on the dataset's own lab exponential/stationary anchors.
3. **The classical acute-planktonic model is contradicted** — acute samples are
   matrix-positive.
4. **Generalisation failed for a mundane reason:** the cohort contains no other
   organism at usable depth. The test is Pseudomonas-specific because the *data*
   are, not because the method is.

## Limitations

- **Depth asymmetry.** Acute median ~11k coding reads vs chronic ~1M. The growth
  axis tolerates this (ratio); the matrix axis does not (detection call).
- **Cross-specimen.** Acute = BALF (human-only rRNA depletion); chronic = CF
  sputum (bacterial rRNA depleted). rRNA excluded from the denominator to
  neutralise the protocol difference, but the specimens differ.
- **Pseudoreplication in the chronic arm.** 15 samples from 5 patients; the
  sample-level p-value overstates confidence (patient-level p ≈ 0.05).
- **Bulk averaging** erases the spatial structure Kolpen's claim is about.
- **Acinetobacter is present but non-baumannii** (A. johnsonii / A. lwoffii), which
  genuinely lack pga/csu — so their matrix-off call is real but uninformative.

## Reproduce

```
build_gcomp.py       # per-organism full-genome indexes + intervals
quant_final.py       # per-sample: dominant organism -> %RP + matrix systems
aggregate_perorg.py  # cohort tables
```
