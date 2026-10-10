# Which organisms are viable for the two-axis test? (2026-10-10)

Answers: *"what about other organisms — are any viable for this analysis?"*

Two independent lines of evidence: (1) a gene-level inventory of biofilm matrix
operons per genus; (2) the actual per-sample dominant genus in the acute cohort,
plus a direct test of the Acinetobacter arm.

## 1. Acute cohort composition (PRJNA1056765, 114 bacterial-infection RNA runs)

Source: figshare deposit `10.6084/m9.figshare.29388539` (RNA genus abundances,
relative abundance per column) + Supplementary Table S2 group labels.

| dominant genus | n |
|---|---|
| Acinetobacter | **38** |
| Pseudomonas | **32** |
| Haemophilus | 10 |
| Mycoplasma | 4 |
| Klebsiella | 2 |
| Streptococcus | 1 |
| Staphylococcus | 1 |
| other (Porphyromonas, Treponema, Neisseria, Fusobacterium, Burkholderia, Prevotella, singletons) | 26 |

Acinetobacter + Pseudomonas = 70/114 (61%). Groups: 114 bacterial / 123 lung cancer /
86 TB / 79 fungal / 32 NC. **There is no COVID group in this cohort.**

## 2. Matrix-operon inventory per genus

| genus | matrix system | genes | viable as a matrix axis? |
|---|---|---|---|
| **Pseudomonas** | alginate / Psl / Pel / CdrA | `algD-…-algA`, `pslA-L`, `pelA-G`, `cdrAB` | **Yes — gold standard** (strain-dependent: alg only in mucoid) |
| **E. coli** | PNAG / curli / cellulose / colanic acid | `pgaABCD`, `csgBAC+csgD`, `bcsRQABZC+EFG`, `wca*` | Yes (curli/cellulose cryptic at 37 °C) |
| **S. aureus** | PIA/PNAG | `icaADBC` (+`icaR`) | Yes (phase-variable; ica-independent strains exist) |
| **Klebsiella** | capsule / type-3 fimbriae | `cps` (conserved `wzi/wza/wzb/wzc/galF`), **`mrkABCDF`** | Yes — use `mrkA`, not the hypervariable K-locus middle |
| **Acinetobacter** | PNAG / Csu pili | `pgaABCD`, `csuA/B-csuE` | **Only for A. baumannii** — see §3 |
| Stenotrophomonas | none defined | `smf-1` (fimbriae), `rpfF` (QS) | Weak — no EPS operon |
| S. pneumoniae | capsule (inverse to biofilm) | `cps*` | **Backwards** — capsule down in biofilm |
| H. influenzae | none | `hmw*`, `hia`, `lic*` | Weak — matrix is host eDNA; genes phase-variable |
| M. catarrhalis | none | `uspA1`, `hag` | Weak — no EPS; poly(G) phase variation |
| M. pneumoniae | none identified | — | **Not possible** — GlcNAc polymer exists but no biosynthetic operon assigned |

## 3. The Acinetobacter arm is not viable — verified directly

I downloaded and mapped all 21 Acinetobacter-dominant acute samples (multi-strain
index, 6 Acinetobacter genomes + 34 decoys, read-level dedup).

**Result: matrix detection is zero in 21/21 samples — and that is correct, not a bug.**

Strain-level assignment of the mapped reads:

| sample | reads assigned | dominant strain |
|---|---|---|
| SRR27343978 | 5,619 | **A. johnsonii** 4,460 · A. lwoffii 793 · A. baumannii 177 |
| SRR27343990 | 3,073 | **A. lwoffii** 1,787 · A. johnsonii 1,071 · A. baumannii 100 |

These are **non-baumannii Acinetobacter** — low-virulence environmental/commensal
species. Checking their annotations: `A. johnsonii` and `A. lwoffii` carry **0**
pga/csu CDS (A. baumannii 17978 carries 13, AB5075 12, A. nosocomialis 11,
A. pittii 7). The pga/csu biofilm operons are essentially an *A. baumannii* trait.

So the "matrix-OFF" call is **real but not about biofilm regulation** — these
organisms simply do not have the machinery. Consistent with the earlier
reagent-contaminant finding: the Acinetobacter signal in this cohort looks like
environmental/kitome background rather than A. baumannii pneumonia.

## 4. Conclusion — the viable scope

**The matrix axis is measurable in ~34 of 114 bacterial samples: Pseudomonas (32)
and Klebsiella (2).** Everything else is either absent from the cohort in
matrix-measurable form or has no clean matrix operon.

- The **growth axis (%RP)** is organism-agnostic and works wherever depth allows.
- The **matrix axis** restricts the design to Pseudomonas/Klebsiella.
- Expanding to Acinetobacter does **not** extend the matrix axis — those samples
  are non-baumannii species without the operons.

Practical consequence for the paper: the two-axis test is a **Pseudomonas-centred**
result, with Klebsiella as a small secondary. That is a narrower but honest scope,
and it should be stated as such rather than presented as a general biofilm finding.

## Gene-name traps to encode in the pipeline

1. `csuA/B` is **one ORF**, not two.
2. `pgaABCD` exists in both E. coli and Acinetobacter — map per-genus, never merged.
3. `wza/wzb/wzc/wzi` are shared across Klebsiella / Acinetobacter / E. coli capsule — not genus-discriminative.
4. Klebsiella K-locus does **not** use `wzm`/`wzt`.
5. Pneumococcal `wzg/wzh/wzd/wze` = `cpsA/cpsB/cpsC/cpsD`.
6. Mycoplasma `cpsG` is a paired-plate protein, **not** a capsule gene.
7. Only 11 of 15 `psl` genes are required (`pslB/M/N/O` dispensable).
8. `bap` is largely absent from human S. aureus and often truncated in A. baumannii.
