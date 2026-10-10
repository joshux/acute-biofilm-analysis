# %TF robustness check — is the acute-vs-chronic growth contrast an RP-family artifact?

**Question.** The committed growth metric `%RP_coding` = (reads over PAO1 ribosomal-protein
genes) / (PAO1 coding reads, rRNA/tRNA excluded) is built entirely from one gene family, and the
acute arm is shallow (median ~10 k coding reads vs ~950 k in chronic). A reviewer can fairly ask
whether the contrast is a quirk of the RP family or of the shallow acute depth. This analysis
answers that with a **second, gene-family-independent proxy computed on the same alignments**.

## The independent proxy: %TF

`%TF` = reads overlapping PAO1 **cat=="growth"** intervals **minus every RP-family gene**
(`^rp[sl][A-Z0-9]+$|^rpm[A-Z0-9]+$`), divided by the **identical** denominator used by
`%RP_coding` (PAO1-mapping reads outside rRNA/tRNA). Mapper and filters are byte-for-byte those of
`quantify_rp2.py`:

```
minimap2 -x sr -t 2 -a --secondary=no compfull_pseu.mmi <fastq>
skip @ lines · skip flag&4 · skip MAPQ<20 · skip refs starting with "decoy|"
```

The resulting TF set is **29 non-RP genes** (translation/transcription + replication/division
machinery): `tufA tufB fusA-like…` — specifically `dnaA dnaB dnaE dnaG dnaK dnaN ftsA ftsI ftsK
ftsL ftsQ ftsW ftsZ gyrA gyrB infA infB infC parC parE relA rpoA rpoB rpoC rpoS spoT tsf tufA
tufB`. The 52 RP-family genes inside the panel's `growth` set were **excluded**, so `%TF` shares no
gene with `%RP`. (29 ≫ the 10-gene floor that would have aborted this task.) Quantifier:
`quantify_tf.py`; per-sample output `out/{SRR}_tf.json`; driver `run_tf.sh`.

**Denominator integrity check.** For every sample, `coding_reads` in `out/{SRR}_tf.json` is
identical to `out/{SRR}_rp2.json` (same mapper, same reads, same rRNA/tRNA exclusion). The two
metrics therefore differ *only* in the numerator gene set — exactly the comparison we want.

## Task 0 — acute-arm patient structure

ENA filereport for **PRJNA1056765** (`ena_prjna1056765.tsv`):

| | |
|---|---|
| runs in `acute_runs.txt` | 16 |
| runs resolved in ENA | 16 |
| **unique BioSamples** | **16** |
| repeated sample_accessions | none |
| runs used by the committed analysis | 15 (`SRR27343234` has 208 reads and no `%RP` value) |

**Conclusion: no repeated patients.** Each acute run is a distinct BioSample, so run-level n =
patient-level n and no clustered/patient-level statistics are required on the acute arm.

## Per-group results (medians)

| group | n | median %TF | median %RP_coding | %TF range |
|---|---|---|---|---|
| acute BALF | 15 | **5.26 %** | 10.68 % | 3.04 – 7.91 |
| chronic in vivo sputum | 15 | **2.64 %** | 4.60 % | 1.58 – 7.48 |
| chronic lab exponential | 12 | **6.42 %** | 13.64 % | 4.32 – 7.12 |
| chronic lab stationary | 11 | **2.42 %** | 2.10 % | 1.26 – 4.14 |

`%TF` reproduces the same group ordering as `%RP_coding`: lab exponential (6.42 %) ≫ lab stationary
(2.42 %), and acute (5.26 %) > chronic in vivo (2.64 %). The internal anchor holds independently:
**lab exp vs lab stat, %TF, Mann-Whitney p < 1e-4** (U=132, z=4.06).

## Correlation: %TF vs %RP_coding (Spearman)

| scope | n | ρ |
|---|---|---|
| pooled, all four groups | 53 | **0.91** |
| acute only | 15 | **0.929** |
| chronic in vivo only | 15 | 0.925 |
| chronic lab exponential | 12 | 0.741 |
| chronic lab stationary | 11 | 0.155 |

The two metrics rank samples almost identically everywhere except the lab-stationary group, whose
%RP values are squeezed into a 0.97–3.83 % band (ρ is uninformative over such a narrow range, not
evidence of disagreement).

## The decisive contrast: acute vs chronic in vivo

| metric | median acute | median chronic | Mann-Whitney p |
|---|---|---|---|
| **%TF** (independent) | 5.26 % | 2.64 % | **p = 0.0008** |
| %RP_coding (committed) | 10.68 % | 4.60 % | p = 0.0015 |

**Depth-restricted** (acute samples with ≥ 10 000 coding reads; n = 8, the shallowest excluded):

| metric | median acute (deep) | median chronic | Mann-Whitney p |
|---|---|---|---|
| **%TF** | 4.70 % | 2.64 % | **p = 0.0118** |
| %RP_coding | 9.71 % | 4.60 % | p = 0.0142 |

The contrast survives both the change of gene family **and** the depth floor.

## Verdict

**The independent proxy CONCORDS with %RP.** `%TF` — 29 non-ribosomal translation/transcription and
replication genes, on the same reads and the same rRNA-excluded denominator — recovers the same
group structure (lab exponential ≫ lab stationary, acute > chronic), ranks samples in near-perfect
agreement with `%RP_coding` (pooled ρ = 0.91; ρ = 0.93 within each arm), and gives a
**statistically significant acute-vs-chronic separation in the same direction** (p = 0.0008, versus
p = 0.0015 for `%RP_coding`). Restricting the acute arm to its eight deepest samples — removing the
shallow-depth confound a reviewer would raise — the contrast is still significant (p = 0.0118,
versus p = 0.0142 for the committed metric). The growth signal in the acute arm is therefore **not**
an artifact of the ribosomal-protein gene family or of acute-arm sequencing depth.

## Caveats (honest)

- The acute arm remains shallow in absolute terms (median ~10 k coding reads); the depth-restricted
  test uses only n = 8 acute samples, so its power is limited. The direction and significance are
  stable, but the effect size is estimated from few reads.
- The lab-stationary group's %TF and %RP_coding medians are close (2.42 % vs 2.10 %), i.e. `%TF` is
  a slightly "noisier" proxy at the low-expression end; the per-group correlation there is weak
  (ρ = 0.155) for that reason.
- Both proxies are computed against the same reference index (`compfull_pseu.mmi`) and the same
  decoy set, so they are independent in *gene family* but not in *reference/mapper*; a mapping bias
  shared by all PAO1 genes would affect both equally.

## Files

- `quantify_tf.py` — new quantifier · `run_tf.sh` — sequential driver
- `out/{SRR}_tf.json` — per-sample counts (×54)
- `tf_robustness.py`, `tf_task0_patch.py` — Task-2 analysis
- `out/tf_robustness_summary.json` — all numbers
- `ena_prjna1056765.tsv` — Task-0 ENA metadata

## Samples that could not be processed

**None.** All 54 target fastq files were present in `fastq/` and every mapping succeeded
(15 acute + 15 chronic in vivo + 12 lab exponential + 11 lab stationary = 53 in the primary panel;
`SRR27343234`, the 16th line of `acute_runs.txt`, was additionally processed — 208 reads, no
committed `%RP` value — and is excluded from the 15-sample acute arm, matching the committed
analysis).
