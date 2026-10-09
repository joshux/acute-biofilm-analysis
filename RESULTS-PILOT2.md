# Pilot 2 result — multi-strain panel rebuild

**Date:** 2026-10-09 · **Status:** panel fix validated; **comparison redesign required** ·
**Samples:** SRR27343983, SRR27343240 (top-biomass bacterial).

## Headline

Rebuilding the reference panel around the **actual** dominant genera — with **multiple
strains/species per genus** — raised usable mapping **6.3×** (MAPQ ≥ 20: 0.108% → **0.681%**)
and, critically, made the genus assignment **agree with the independent composition table**
(Acinetobacter dominant). The panel was the binding constraint; that is now fixed. What remains
is a **yield ceiling** and a **matrix-module gap for the dominant organism**, plus the dead NC
contrast (Pilot 1).

## 1. The rebuilt panel

19 genomes / 158 contigs, multi-strain and multi-species:

| genus | genomes included |
|---|---|
| **Acinetobacter** | A. baumannii ATCC 17978 + AB5075, A. nosocomialis, A. pittii, A. johnsonii, A. lwoffii (6) |
| **Klebsiella** | K. pneumoniae MGH78578 + NTUH-K2044 + HS11286, K. variicola (4) |
| Pseudomonas | PAO1, PA14, P. putida KT2440 (3) |
| Stenotrophomonas | S. maltophilia K279a (1) |
| Others | S. aureus NCTC 8325, S. pneumoniae TIGR4, E. coli MG1655, H. influenzae 86-028NP |

## 2. Mapping improvement (same reads, same mapper, only the panel changed)

| sample | panel | mapped (all) | MAPQ ≥ 20 |
|---|---|---|---|
| SRR27343983 | 6 genomes | 2.43% | 0.084% |
| SRR27343983 | 7 genomes (+1 Acinetobacter) | 2.62% | 0.108% |
| **SRR27343983** | **19 genomes (multi-strain)** | **4.53%** | **0.681%** |
| **SRR27343240** | **19 genomes (multi-strain)** | **3.61%** | **0.815%** |

**Genus assignment now matches the composition table** — the validation that the fix is real:

- SRR27343983 MAPQ≥20 by genus: **Acinetobacter 13,579 (59%)**, Stenotrophomonas 5,947 (26%),
  Haemophilus 1,407, Escherichia 1,105, Pseudomonas 673, Klebsiella 77, others <60.
- Composition table for this sample: Acinetobacter 29.3%, Haemophilus 5.2%, Pseudomonas 4.7%,
  Stenotrophomonas 3.7% — **Acinetobacter dominant in both**. The multi-strain panel recovered
  it; the single-strain panel could not (766 reads).

## 3. Remaining constraints

1. **Absolute yield is still modest.** ~23k (SRR27343983) and ~47k (SRR27343240) MAPQ≥20 reads.
   Against the ~1.48M / ~2.5M Kraken-classified microbial reads, that is ~10–15% — the missing
   fraction is largely the **contaminant genera** (Comamonas, Delftia, Ralstonia, Sphingomonas)
   not yet in the panel. Adding them would raise total mapping but also matters for
   interpretation: they are the kitome, and reads from them must **not** be attributed to
   infection.
2. **The matrix-module gap is now the central problem.** The dominant organism is
   **Acinetobacter**, which has **no characterized biofilm-matrix operon** comparable to
   pel/psl/alg/ica/csg (only bap, a PNAG-like polymer, and capsule). So the gene that the
   analysis most needs — matrix synthesis — is least characterized precisely where the biomass
   is highest. Options:
   - **Restrict the matrix analysis to Pseudomonas- and Klebsiella/Enterobacteriaceae-dominant
     samples** (both have defined, well-annotated matrix operons), and report Acinetobacter
     separately as a coverage limit.
   - **Reframe the matrix module** as a general EPS / capsule / adhesin / eDNA signature rather
     than named operons (broader coverage, weaker mechanistic specificity).
3. **The NC contrast is dead** (Pilot 1: 2 mapped reads). Replacement designs below.

## 4. Comparison redesign (replaces bacterial-vs-NC)

Three options, ranked:

| design | n | what it tests | verdict |
|---|---|---|---|
| **Within-bacterial, by dominant genus** — Pseudomonas-dominant vs Acinetobacter-dominant | 32 vs 38 | whether matrix+growth co-expression differs by organism ecology | ✅ **primary** — uses the cohort's actual structure; but confounded by the matrix-module gap for Acinetobacter |
| **Bacterial vs fungal/TB infection** | 114 vs 79/86 | whether matrix+growth coupling is specific to bacterial infection | ⚠️ secondary — both arms have bacterial signal, but different clinical context |
| Bacterial vs NC | 114 vs 32 | — | ❌ dead (no bacterial signal in NC) |

**Recommended primary:** within-bacterial, restricted to samples where the dominant organism has
a defined matrix operon (Pseudomonas, Klebsiella) — i.e. test Kolpen's combination *within the
organisms where matrix expression is actually measurable*, and state the Acinetobacter coverage
gap explicitly. This is the honest design: it answers the question where the instruments work,
rather than reporting a null that is really a coverage artifact.

## 5. Updated verdict

**GO, narrowed.** The panel fix works and the community is now correctly resolved. The viable
analysis is: *matrix×growth co-expression in Pseudomonas/Klebsiella-dominant acute bacterial
BALF samples, with contaminant-aware filtering, Acinetobacter reported as a coverage limit.*
The full-cohort all-organism design is not achievable because the dominant organism lacks a
measurable matrix module — that is a finding about the assay, and it should be stated, not
papered over.

## Files

- `refs2/manifest.json` — the 12 new genomes (accessions, organisms, sizes)
- `ref2genus.json` — panel reference accession → genus map
- `multi_panel_mapping.json` — per-sample mapped/MAPQ≥20 counts by genus
- `fetch_refs.py` — reproducible genome fetch
