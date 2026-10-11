# Normalization results — does the growth contrast survive the load confound?

**Date:** 2026-10-11. Follows `COVID-ARM-ASSESSMENT.md`. Question: is %RP biased by Pseudomonas
load and/or read length, and does the acute-vs-chronic contrast survive once those are controlled?

**Answer: the contrast survives every control applied. But the load confound is real, its direction
is the OPPOSITE of what I first reported, and load cannot be fully ruled out because the two arms
occupy disjoint load ranges.**

---

## 1. Design

Every sample was **re-mapped at a constant 200,000 reads**, so sequencing depth is held fixed and
load varies only through composition. Panels: 8 acute + 8 chronic balanced (largest `pao1_coding`),
extended to all 18 acute + 17 chronic available at constant depth, plus 8 COVID. Per-sample outputs
in `out_norm/`; script `norm_map.py`.

## 2. The load effect is real — but it is not what I said

Within-arm Spearman(load, %RP) at constant depth:

| arm | n | load range | Spearman(load, %RP) |
|---|---|---|---|
| acute | 18 | 0.0001 – 0.033 | **+0.035** |
| chronic | 17 | 0.509 – 0.958 | +0.514 |
| COVID (all quantified) | 40 | 0.0008 – 0.842 | +0.377 |
| COVID, thin-signal samples dropped | 27 | — | **−0.029** |

**The acute arm shows no load dependence at all across a 250-fold load range (ρ = +0.035).** That is
the single most reassuring number here: if %RP were load-driven, the acute arm's internal variation
would track its own load. It does not.

**The COVID +0.60 I reported was a low-load-tail artifact.** It collapses to ~0 as soon as
thin-signal samples are excluded (+0.377 for all 40 → −0.029 at coding ≥500). It was never a
high-load effect.

## 3. CORRECTION to COVID-ARM-ASSESSMENT.md

That file states the mechanism as *"high Pseudomonas load mechanically raises apparent %RP"* via the
MAPQ≥20 uniqueness filter. **The direction is wrong.** Measured MAPQ-survival of RP reads relative to
other coding reads:

| arm | unique fraction RP | unique fraction other coding | RP advantage | inflation factor |
|---|---|---|---|---|
| acute (low load) | 0.669 | 0.604 | +0.065 | **1.108×** |
| chronic (high load) | 0.958 | 0.931 | +0.027 | **1.029×** |

The filter does preferentially retain RP reads — but the **RP-specific advantage is larger at LOW
load** (1.108×) than at high load (1.029×). A conservation-biased uniqueness filter therefore
inflates %RP **at low load**, the opposite of what I wrote. Consequences:

- COVID's elevation is not a high-load artifact. Its highest-%RP samples sit at load 0.008–0.06,
  and its **highest-load sample has the lowest %RP (0.58%)**.
- The COVID elevation is explained by the low-load conservation effect plus its 146 bp read length
  — not by load in the direction I claimed.

## 4. Depth and read length are not drivers

**Within-sample depth titration** (%RP at 50k / 100k / 200k reads, same sample):

| sample | arm | 50k | 100k | 200k | spread |
|---|---|---|---|---|---|
| SRR27343249 | acute | 19.28 | 20.02 | 19.50 | 0.74 pp |
| ERR2275091 | chronic | 6.90 | 6.84 | 6.85 | 0.06 pp |
| SRR26629332 | COVID | 32.46 | 32.94 | 32.24 | 0.70 pp |

Flat. Depth is not the driver — which means the cross-sample correlation had to be composition, and
it was (the low-load tail).

**Read length** is constant within arm (acute 50 bp, chronic ~74 bp, COVID ~146 bp) and differs only
between arms — a dataset property. Base-weighting is a no-op (ρ = 1.00 vs %RP; max difference
0.107 pp). **Truncating chronic and COVID reads to 50 bp to match acute leaves the contrast
intact**: acute 9.43% vs chronic 3.46%, p = 0.0117.

## 5. The acute-vs-chronic contrast under every control

Constant 200k reads, native reads, my own recomputation from `out_norm/*_norm.json`:

| metric | acute median | chronic median | ratio | Mann-Whitney p |
|---|---|---|---|---|
| **%RP (coding denominator)** | 10.94% (n=18) | 5.18% (n=17) | 2.11× | **p = 0.0006** |
| **RP / conserved-control-gene ratio** | 4.72 | 2.82 | 1.67× | **p = 0.0146** |
| %RP, base-weighted | 10.94% | 5.31% | 2.06× | p = 0.0008 |

Balanced 8 vs 8 subset (largest `pao1_coding`):

| metric | acute | chronic | ratio | p |
|---|---|---|---|---|
| %RP | 9.50% | 5.31% | 1.79× | p = 0.046 |
| RP / control ratio | 4.63 | 2.68 | 1.73× | p = 0.059 |

Full control battery:

| control | result |
|---|---|
| constant depth (200k) | **survives** (p = 0.0006) |
| base-weighted instead of read-count | **survives** (p = 0.0008) |
| MAPQ filter removed entirely | **survives**, stronger (9.73 vs 4.51, p = 0.0030) |
| read length truncated to match | **survives** (p = 0.0117) |
| conserved-control-gene normalization | **survives** on full panel (p = 0.0146); marginal on 8v8 (p = 0.059) |
| abundance matching | **not testable** |

## 6. The residual caveat, stated plainly

**Load and cohort are near-perfectly confounded.** Acute spans load 0.0001–0.033; chronic spans
0.509–0.958 — a **55× gap with zero overlap**, so no abundance-matched comparison exists at any
subset size. This is not a fixable analysis problem; it is a property of the datasets (CF sputum is
a near-monoculture of Pseudomonas; acute BALF is not).

What can be said: the acute arm shows **no internal load dependence** (ρ = +0.035) across a 250-fold
range, and the contrast survives depth, length, MAPQ-filter, and control-gene controls. That is
consistent with a genuine growth difference. It does **not** prove load is irrelevant — only that
the load effect we can measure, in the arm where we can measure it, is absent.

## 7. What this changes

1. **The headline survives, and is now defended against the obvious reviewer question.** The paper
   can state: the contrast holds at constant depth, under base-weighting, with the MAPQ filter
   removed, at matched read length, and under conserved-control-gene normalization.
2. **A third control joins the paper:** RP/control-gene ratio as an abundance- and
   conservation-robust variant. This is a new, cheap robustness check with a real precedent
   (abundance-normalized metrics).
3. **COVID stays out of the biology** — but now for a precise reason (low-load conservation artifact
   + 146 bp reads), not a hand-wave. It remains useful as a third instantiation of the artifact
   class, with the corrected mechanism.
4. **The limitation to declare:** load and cohort are confounded; no abundance-matched contrast is
   possible.

---

*Reproduce: `python3 run_all.py main|titration|trim` → `python3 final_analysis.py`.
Per-sample numbers: `out_norm/norm_results.json`. Raw per-sample: `out_norm/*_norm.json`.*
