# Existing methods-forward short reports in biomedicine — verified precedents

**Date:** 2026-10-10. Purpose: precedents + venue formats for our methods-forward short report
(two-axis metatranscriptomic test + two artifact classes). Every DOI below was resolved against
Crossref / doi.org / Europe PMC by a verification pass; items that could not be verified are
marked. This file doubles as the citation pool for the manuscript's framing paragraphs.

---

## 1. The artifact / cautionary-note genre (closest precedents)

| # | Citation | Format / length | Artifact documented |
|---|---|---|---|
| 1 | **Chalkley 2013**, *When Target–Decoy FDR Estimations Are Inaccurate and How to Spot Instances*, J Proteome Res 12:1062–1064. [10.1021/pr301063v](https://doi.org/10.1021/pr301063v) | 3-pp single-author note, no new data | Global target-decoy FDR is subset-invalid (~90% of incorrect hits track search space, not sample). **Closest formal template for our paper.** |
| 2 | **Halpin, Jangi & Street 2019**, *Multimapping confounds ribosome profiling analysis: a case-study of the Hsp90 molecular chaperone*, Proteins 87:1212–1221. [10.1002/prot.25766](https://doi.org/10.1002/prot.25766) | Case-study | Multimappers are 25–80% of Ribo-seq reads; Hsp90 expression changes **10-fold** with mapper handling; multimapping manufactures reproducible-looking false peaks. **Our MAPQ payload, almost exactly.** |
| 3 | **Paz, Warger & Taher 2024**, *Disregarding multimappers leads to biases in the functional assessment of NGS data*, BMC Genomics 25:455. [10.1186/s12864-024-10344-9](https://doi.org/10.1186/s12864-024-10344-9) | Full research article | Discarding multimappers under-quantifies 4–6% of genes non-randomly (MHC/antigen-presentation) → GSEA returns a coherent but wrong biological story. |
| 4 | **Salter et al. 2014**, *Reagent and laboratory contamination can critically impact sequence-based microbiome analyses*, BMC Biology 12:87. [10.1186/s12915-014-0087-z](https://doi.org/10.1186/s12915-014-0087-z) | ~13-pp research article | Kit DNA swamps low-biomass libraries; the canonical "protocol choice silently reverses the biology" paper. |
| 5 | **Tsou, Olesen, Alm & Snapper 2020**, *16S rRNA sequencing analysis: the devil is in the details*, Gut Microbes 11:1139–1142. [10.1080/19490976.2020.1747336](https://doi.org/10.1080/19490976.2020.1747336) | **4 pp** | Forward/reverse/merged reads + QIIME1-vs-2 give three taxonomies for one sequence; self-described "cautionary tale". |
| 6 | **Verwilt et al. 2020**, *When DNA gets in the way: a cautionary note for DNA contamination in extracellular RNA-seq studies*, PNAS 117:18934–18936. [10.1073/pnas.2001675117](https://doi.org/10.1073/pnas.2001675117) | **3-pp Letter** | >95% of SILVER-seq reads are off-exon → the "RNA" signal is cell-free DNA; collapses 44-fold after junction filtering. Title literally announces the genre. |
| 7 | **Bogdanow, Zauber & Selbach 2016**, *Systematic Errors in Peptide and Protein Identification… by Modified Peptides*, MCP 15:2791–2801. [10.1074/mcp.M115.055103](https://doi.org/10.1074/mcp.M115.055103) | Research article | Modified peptides cause 20–50% of false IDs and are the *highest-scoring* false positives — an artifact that preferentially corrupts the confident hits. |
| 8 | **Knudsen & Chalkley 2011**, *The Effect of Using an Inappropriate Protein Database…*, PLoS ONE 6:e20873. [10.1371/journal.pone.0020873](https://doi.org/10.1371/journal.pone.0020873) | ~7 pp | Wrong-scope database fabricated confident Iridovirus/Nosema IDs from bee spectra; re-analysis corrected them. |
| 9 | **Marcelino, Holmes & Sorrell 2020**, *The use of taxon-specific reference databases compromises metagenomic classification*, BMC Genomics 21:184. [10.1186/s12864-020-6592-2](https://doi.org/10.1186/s12864-020-6592-2) | Research article | A fungal DB + conserved rRNA recovered turtles and frogs from human gut; 85% of "fungal" hits were rRNA. Reference-panel design artifact — **same class as our multi-genome trap**. |
| 10 | **Walker et al. 2020**, *Non-specific amplification of human DNA is a major challenge for 16S…*, Sci Rep 10:16455. [10.1038/s41598-020-73403-7](https://doi.org/10.1038/s41598-020-73403-7) | Research article | V3–V4 primers amplify human DNA in high-host-biomass samples; self-described "methodological warning and remedy". |
| 11 | **Gomez-Alvarez, Teal & Schmidt 2009**, *Systematic artifacts in metagenomes from complex microbial communities*, ISME J 3:1314–1317. [10.1038/ismej.2009.72](https://doi.org/10.1038/ismej.2009.72) | **4 pp** | 11–35% of 454 reads are artificial replicates inflating abundance; early short-format artifact precedent. |
| 12 | **Hicks, Townes, Teng & Irizarry 2018**, *Missing data and technical variability in single-cell RNA-sequencing experiments*, Biostatistics 19:562–578. [10.1093/biostatistics/kxx053](https://doi.org/10.1093/biostatistics/kxx053) | 17 pp, statistics journal | Detection-rate differences fabricate "novel" cell groups; a stats journal publishing "this result is an artifact" as its whole payload. |

**Compositional/absolute-abundance anchors** (payload 3 context): **Vandeputte et al. 2017** Nature
551:507–513, [10.1038/nature24460](https://doi.org/10.1038/nature24460) (the Bacteroides–Prevotella
trade-off is "an artefact of relative microbiome analyses"); **Gloor et al. 2017** Front Microbiol
8:2224, [10.3389/fmicb.2017.02224](https://doi.org/10.3389/fmicb.2017.02224).

**Public-data misannotation precedents** (composition-table disagreement context): **Javed et al.
2020** Nat Commun 11:3742, [10.1038/s41467-020-17453-5](https://doi.org/10.1038/s41467-020-17453-5)
(~1% sample swaps across 8,851 ENCODE datasets); **Toker, Feng & Pavlidis 2016** F1000Research
5:2775, [10.12688/f1000research.9471.2](https://doi.org/10.12688/f1000research.9471.2).

⚠️ **CORRECTION vs our earlier notes:** the COVID-proteomics reading plan
(joshux/covid-proteomic-analysis, docs/READING-PLAN.md Block 6) cites "Knudsen & Chalkley
(over-calling mirror image)" in JASMS. **No such paper was found** — the mirror-image/palindromic
decoy-peptide problem is documented in Moosa et al. 2020, J Proteome Res 19:1029
[10.1021/acs.jproteome.9b00538](https://doi.org/10.1021/acs.jproteome.9b00538) and Lee et al. 2021,
Proteome Science 19:13 [10.1186/s12953-021-00179-7](https://doi.org/10.1186/s12953-021-00179-7).
The citation in that reading plan is likely conflated and should be corrected there. The verified
Knudsen & Chalkley paper is the 2011 PLoS ONE wrong-database study (entry 8).

---

## 2. Growth-metric lineage (must-cite for the %RP/%TF axis)

- **Gifford, Sharma, Booth & Moran 2013**, ISME J 7:281–298,
  [10.1038/ismej.2012.96](https://doi.org/10.1038/ismej.2012.96) — origin of %RP as an in-situ
  growth proxy; **Gifford, Sharma & Moran 2014** Front Microbiol 5:185,
  [10.3389/fmicb.2014.00185](https://doi.org/10.3389/fmicb.2014.00185) — operationalized across
  200 taxa. Our metric descends from these; benchmark explicitly.
- **Long, Hou, Ignacio-Espinoza & Fuhrman 2020**, ISME J 14:2680–2691,
  [10.1038/s41396-020-00773-1](https://doi.org/10.1038/s41396-020-00773-1) — "validated metric A,
  invalidated metric B" template; PTR indices fail their benchmark.

## 3. Venue formats that accept this genre (with verified examples)

| venue / format | constraints | methods-forward published examples (verified) |
|---|---|---|
| **mSphere "Observations"** | 1,200 words; ≤2 figs; ≤25 refs; abstract ≤250 | Schloss 2021 (ASVs split genomes) [10.1128/msphere.00191-21](https://doi.org/10.1128/msphere.00191-21); Williams 2018 (virome method comparison) [10.1128/msphere.00311-18](https://doi.org/10.1128/msphere.00311-18); Gao 2024 (PacBio 16S) [10.1128/msphere.00770-24](https://doi.org/10.1128/msphere.00770-24) |
| **mSystems "Observations"** | 1,200 words; ≤2 figs; ≤25 refs | Brennan 2024 (well-to-well contamination mitigation) [10.1128/msystems.00985-24](https://doi.org/10.1128/msystems.00985-24); Armstrong 2021 (UMAP artifacts) [10.1128/msystems.00691-21](https://doi.org/10.1128/msystems.00691-21) |
| **mSystems "Methods and Protocols"** | 5,000 words; must include validation vs state-of-the-art | BiomeHorizon 2022 [10.1128/msystems.01380-21](https://doi.org/10.1128/msystems.01380-21); Metapresence 2024 [10.1128/msystems.00213-24](https://doi.org/10.1128/msystems.00213-24) |
| **GigaScience "Technical Note"** | no word limit; structured abstract ≤250 | ASaiM 2018 [10.1093/gigascience/giy057](https://doi.org/10.1093/gigascience/giy057); Meta-NanoSim 2023 [10.1093/gigascience/giad013](https://doi.org/10.1093/gigascience/giad013); 16S benchmarking 2018 [10.1093/gigascience/giy054](https://doi.org/10.1093/gigascience/giy054) |
| **GigaScience "Brief Communication"** | 1,500 words; ≤10 refs; ≤2 figs | xgt 2026 [10.1093/gigascience/giag086](https://doi.org/10.1093/gigascience/giag086) |
| **Bioinformatics "Application Note"** | ~2,600 words / 4 pages | benchdamic (benchmarking) 2023 [10.1093/bioinformatics/btac778](https://doi.org/10.1093/bioinformatics/btac778); Pavian 2020 [10.1093/bioinformatics/btz715](https://doi.org/10.1093/bioinformatics/btz715) |
| **eLife short reports** | ≤1,500 words main text; 3–4 display items | RIM-Deep method 2025 [10.7554/eLife.101143](https://doi.org/10.7554/eLife.101143); SPICE pipeline 2025 [10.7554/eLife.88833](https://doi.org/10.7554/eLife.88833) |
| **npj Biofilms Microbiomes "Brief Communication"** | abstract ≤70 words; 1,000–1,500 words; ~20 refs | Rascovan 2016 [10.1038/s41522-016-0008-8](https://doi.org/10.1038/s41522-016-0008-8); Heidrich 2025 [10.1038/s41522-025-00665-2](https://doi.org/10.1038/s41522-025-00665-2) |
| **Frontiers "Brief Research Report"** | ≤4,000 words; ≤4 figs/tables | Rozas 2022 (MinION 16S vs mock community) [10.3389/fcimb.2021.806476](https://doi.org/10.3389/fcimb.2021.806476) |
| **ISME Communications "Brief Communication"** | (check current page) | Hrovat 2024 (16S region resolution) [10.1093/ismeco/ycae034](https://doi.org/10.1093/ismeco/ycae034) |

**No short format** (do not target for this paper): Microbiome (BMC), BMC Microbiology, PLOS ONE,
PeerJ (has a "Method Paper" type but full-length), JCM, J Bacteriology.

## 4. Genre conventions (from the precedents)

1. The artifact is **named in the title** ("a cautionary note", "the devil is in the details",
   "multimapping confounds…").
2. A re-analysis or controlled comparison **demonstrates** it on real public data.
3. The consequence is stated as a **changed biological conclusion**, not a technical bias.
4. A **mitigation** (filter rule, tool, checklist) is supplied.
5. Length: 3–4 printed pages is an established, respected form (Chalkley 3pp, Verwilt 3pp,
   Tsou 4pp, Gomez-Alvarez 4pp).

Our paper satisfies all five: named artifacts (MAPQ multi-genome drop; composition-table
disagreement), demonstrated on public data with a corrected result, consequences stated as a
reversed biological conclusion, mitigation (single-reference + decoy design, rRNA-excluded
denominator, pre-specified depth floors), and a validated metric as the constructive payload.
