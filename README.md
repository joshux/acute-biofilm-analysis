# acute-biofilm-analysis

Testing Kolpen et al. 2022 (Thorax) — that bacteria in acute lung infection co-express biofilm
matrix and fast growth — on public BALF metatranscriptomes (PRJNA1056765, Tang et al. 2025 Sci Data).

Bacterial acute-infection arm only; COVID contrast arm deferred.

- `PLAN.md` — analysis plan **v2** (post-pilot; multi-strain panel, within-bacterial design, contaminant handling)
- `RESULTS-PILOT.md` — pilot 1: dataset/depth verified; three assumptions broken (community, panel, NC arm)
- `RESULTS-PILOT2.md` — pilot 2: multi-strain panel rebuilt and validated (6.3x mapping gain); matrix-module gap for Acinetobacter
- `refs2/manifest.json`, `ref2genus.json` — the multi-strain reference panel
- `multi_panel_mapping.json` — per-sample mapping counts by genus
