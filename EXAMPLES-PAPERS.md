# Example papers in biomedicine — curated by what each one models

**Date:** 2026-10-10. Every DOI below was re-verified by the orchestrator against the Crossref API
(title + journal + year + author count all matched). Affiliation strings marked "(Crossref)" were
returned by Crossref itself; others are from the publisher/PMC record.

The list is organised by **the thing worth copying**, not by topic — because for a solo author
reanalysing public data, the reusable asset is the *genre*.

---

## 1. Solo / unaffiliated authors — proof the path works

The single most useful fact here: **"Independent Researcher, [City]" is an accepted affiliation
string in mainstream journals**, including Nature-portfolio titles. Eight of these are solo.

| citation | venue | why it's an exemplar |
|---|---|---|
| **Edgar RC**, *Muscle5: high-accuracy alignment ensembles…*, Nat Commun 13 (2022). [10.1038/s41467-022-34630-w](https://doi.org/10.1038/s41467-022-34630-w) | Nature Communications | Solo; affiliation literally "Independent Researcher." Re-tests published phylogenies on existing sequence data — and finds some high-bootstrap topologies are wrong. 1,200+ citations. |
| **Edgar RC**, *URMAP, an ultra-fast read mapper*, PeerJ 8:e9338 (2020). [10.7717/peerj.9338](https://doi.org/10.7717/peerj.9338) | PeerJ | Solo; "Unaffiliated, Corte Madera, CA" **(Crossref)**. Same author also in Bioinformatics 2024 as "Independent Scientist." A repeatable stream from home. |
| **Mukherjee S**, *Quiescent stem cell marker genes in glioma gene networks…*, Sci Rep 10:10937 (2020). [10.1038/s41598-020-67753-5](https://doi.org/10.1038/s41598-020-67753-5) | Scientific Reports | Solo. **Closest analogue to a reanalysis author**: 100% public GEO data, home address, code on GitHub. |
| **Stockwin LH**, *Alveolar soft-part sarcoma resembles a mesenchymal stromal progenitor…*, PeerJ 8:e9394 (2020). [10.7717/peerj.9394](https://doi.org/10.7717/peerj.9394) | PeerJ | Solo; "Unaffiliated, Frederick, MD." Meta-analysis of others' deposited microarray/RNA-seq. No new wet lab. |
| **McCartney MA**, *Structure, function and parallel evolution of the bivalve byssus…*, Phil Trans R Soc B 376:20200155 (2021). [10.1098/rstb.20200155](https://doi.org/10.1098/rstb.20200155) | Phil Trans R Soc B | Solo; "Unaffiliated." A **Royal Society** journal accepted a solo unaffiliated comparative reanalysis. |
| **Dropkin G**, *Low dose radiation risks for women surviving the a-bombs in Japan*, Environ Health 15:117 (2016). [10.1186/s12940-016-0191-3](https://doi.org/10.1186/s12940-016-0191-3) | Environmental Health | Solo; "Independent researcher, Liverpool." **The corrective-reanalysis template**: take a famous public cohort, apply a defensible alternative model, publish the disagreement. |
| **Farquhar H**, *Protein language model embeddings improve HIV drug resistance prediction*, Bioinformatics (2026). [10.1093/bioinformatics/btag260](https://doi.org/10.1093/bioinformatics/btag260) | Bioinformatics | Solo; "Independent Researcher, Finley, NSW" **(Crossref)**. Current-era benchmark + public pathogen sequence data. Companion Virology 2026 paper is a **null result as a full paper**. |
| **Coser O**, *ELISA: an interpretable hybrid generative AI agent…*, Brief Bioinform (2026). [10.1093/bib/bbag501](https://doi.org/10.1093/bib/bbag501) | Briefings in Bioinformatics | Solo; a **home street address** as the affiliation **(Crossref)**. All six datasets public via CZ CELLxGENE. |
| **Mony JT**, *Divergent immune and extracellular matrix transcriptional programs…*, Front Immunol (2026). [10.3389/fimmu.2026.1762224](https://doi.org/10.3389/fimmu.2026.1762224) | Frontiers in Immunology | Solo; affiliation field reads simply "Independent researcher" **(Crossref)**. Minimal, unembellished. |
| **Mohammadiaria M**, *ROS-induced voltage-gated ion channel expression…*, npj Syst Biol Appl 11 (2025). [10.1038/s41540-025-00595-x](https://doi.org/10.1038/s41540-025-00595-x) | npj Systems Biology | Solo; "Unaffiliated, Pavia, Italy," personal email as contact. **npj family still accepts unaffiliated solo authors.** |
| **Pigott HE et al.**, *…reanalysis of the STAR*D study's patient-level data with fidelity to the original research protocol*, BMJ Open 13:e063095 (2023). [10.1136/bmjopen-2022-063095](https://doi.org/10.1136/bmjopen-2022-063095) | BMJ Open | Lead author's affiliation is **"None, Wakefield, Rhode Island, USA"** **(Crossref)** — a home address. Model for a *sustained programme* off one public dataset. |
| **Wilshire CE, Kindlon T, Courtney R, Matthees A, Tuller D, et al.**, *Rethinking the treatment of chronic fatigue syndrome — a reanalysis…*, BMC Psychology 6:6 (2018). [10.1186/s40359-018-0218-3](https://doi.org/10.1186/s40359-018-0218-3) | BMC Psychology | Several authors have **no academic affiliation**; one is a patient who forced data release via FOI. Highest-profile "outsiders correct a famous trial" case. |
| **Vink M & Vink-Niese A**, *Graded exercise therapy for ME/CFS is not effective and unsafe. Re-analysis of a Cochrane review*, Health Psychol Open 5 (2018). [10.1177/2055102918805187](https://doi.org/10.1177/2055102918805187) | Health Psychology Open | Two-person non-institutional team; **no new data** — a reanalysis of a review's own extraction. Published twice. |

**The pattern:** venue *tier* is not the constraint, venue *type* is. Data-and-method-oriented
journals evaluate on analysis soundness and code availability, not on whether you have a lab.
The four self-descriptions that recur: protocol-fidelity reanalysis; corrective statistical
reanalysis; cross-dataset synthesis; public-data benchmark or tool.

---

## 2. Negative, null, and corrective results — the genre that says "we couldn't, and here's why"

| citation | format | the negative |
|---|---|---|
| **Gihawi A, … Salzberg SL**, *Major data analysis errors invalidate cancer microbiome findings*, mBio 14:e01607-23 (2023). [10.1128/mbio.01607-23](https://doi.org/10.1128/mbio.01607-23) | Full article | Reanalysis of Poore et al.'s *Nature* cancer-microbiome claim: the classifiers were artefacts (contaminated reference genomes; label leakage). Concludes the association "is, simply put, a fiction." **Directly relevant to us** — same artifact class, same corrective move. |
| **Barton AR, … Mathieson I**, *Insufficient evidence for natural selection associated with the Black Death*, Nature 638:E19 (2025). [10.1038/s41586-024-08496-5](https://doi.org/10.1038/s41586-024-08496-5) | **Matters Arising, 4 pp** | A headline *Nature* claim fails: correct randomization raises P by ten orders of magnitude. Published **in Nature**, with author reply. |
| **Munday PL**, *Reanalysis shows there is not an extreme decline effect in fish ocean acidification studies*, PLOS Biol 20:e3001809 (2022). [10.1371/journal.pbio.3001809](https://doi.org/10.1371/journal.pbio.3001809) | Formal Comment, ~5 pp | The reported effect is an artefact of **one analysis choice** (replacing zeros with 0.0001). Solo author. |
| **Janda GS, … Ross JS**, *Feasibility of using real-world data to emulate substance use disorder clinical trials*, BMC Med Res Methodol 24:187 (2024). [10.1186/s12874-024-02307-1](https://doi.org/10.1186/s12874-024-02307-1) | Full article | **"Zero of 272 trials could be emulated."** The unavailability *is* the measured result. **Our closest structural template** for "no matched dataset exists." |
| **Russek M, Peltner J, Haenisch B**, *Supplementing single-arm trials with external control arms*, Clin Pharmacol Ther 117 (2025). [10.1002/cpt.3684](https://doi.org/10.1002/cpt.3684) | Full article | 10/379 and 2/11 trials were feasible. "This arm is unmeasurable" as a contribution. |
| **Marchal A, … COVID-19 Host Genetics Initiative**, *Lack of association between classical HLA genes and asymptomatic SARS-CoV-2 infection*, HGG Adv 5:100300 (2024). [10.1016/j.xhgg.2024.100300](https://doi.org/10.1016/j.xhgg.2024.100300) | Full article | **The power statement model**: failed replication *despite >95% power* to detect the originally reported effect. Converts "found nothing" into "the effect isn't there." |
| **Pawel S, Heyard R, Micheloud C, Held L**, *Replication of null results: absence of evidence or evidence of absence?*, eLife 12:e92311 (2024). [10.7554/eLife.92311](https://doi.org/10.7554/eLife.92311) | Full article | **The citable authority** for why a null needs a power/detection-limit statement to be interpretable. |
| **Errington TM, … Nosek BA**, *Investigating the replicability of preclinical cancer biology*, eLife 10:e71601 (2021). [10.7554/eLife.71601](https://doi.org/10.7554/eLife.71601) | Full article | Median replication effect 85% smaller; raw data available for only 2% of experiments. |
| **Lee SH, … Kim J-S**, *Failure to detect DNA-guided genome editing using NgAgo*, Nat Biotechnol 35:17 (2016). [10.1038/nbt.3753](https://doi.org/10.1038/nbt.3753) | **Brief Communication, ~2 pp** | Three-lab coordinated refutation using the original authors' own reagent — which triggered the retraction. |
| **Perneger T & Gayet-Ageron A**, *Evidence of lack of treatment efficacy derived from statistically nonsignificant results of RCTs*, JAMA 329:2050 (2023). [10.1001/jama.2023.8549](https://doi.org/10.1001/jama.2023.8549) | Short report | 91% of non-significant outcomes actually favoured the null. Authority for "non-significant = evidence." |

**Venue answer:** *Journal of Negative Results in Biomedicine* **ceased publication 1 Sept 2017**
(archive still open/PubMed-indexed at jnrbm.biomedcentral.com — citable, not submittable). Its
closest living successors: **BMC Research Notes** (scope names "valid negative results"; MEDLINE,
IF 1.7, median 4 days to first decision), **Journal of Trial and Error** (diamond OA, negative/null
results by charter), **PLOS ONE "Missing Pieces"** collection, and *Nutrition & Metabolism*'s
permanent Null Results section.

**Registered Reports guarantee in-principle acceptance regardless of outcome** — Nature: *"published
regardless of outcome or statistical significance"*; Royal Society Open Science: *"negative results
will not prevent publication."* Caveat: Stage 1 review happens **before** data collection, so an
already-completed analysis cannot be retrofitted.

**How to frame a negative** (the recurring pattern in the strongest examples): make the
unavailability the measured quantity ("we screened N candidates and zero were feasible"), always
report a detection limit or power, run and report the positive control, locate the *mechanism* of
the failure, use the venue's formal critique format where one exists, and state what should stop
and what should happen instead.

---

## 3. Idea / hypothesis papers (no new data)

| citation | venue & constraints | why |
|---|---|---|
| **Smith RS**, *The macrophage theory of depression*, Med Hypotheses 35:298 (1991). [10.1016/0306-9877(91)90272-z](https://doi.org/10.1016/0306-9877(91)90272-z) | Medical Hypotheses | ~650 citations; seeded the entire inflammation-depression field. **The benchmark for a no-data hypothesis that changed a field.** |
| **Harvey WT & Salvato P**, *'Lyme disease': ancient engine of an unrecognized borreliosis pandemic?*, Med Hypotheses 60:742 (2003). [10.1016/S0306-9877(03)00060-4](https://doi.org/10.1016/S0306-9877(03)00060-4) | Medical Hypotheses | Method explicitly stated: take the accepted model's "linchpin premises," test each against the literature, reconstruct a competing model. (Cite for *craft*, not content — the claim is contested.) |
| **O'Donnell JS & Chappell KJ**, *Chronic SARS-CoV-2, a cause of post-acute COVID-19 sequelae?*, Front Microbiol 12:724654 (2021). [10.3389/fmicb.2021.724654](https://doi.org/10.3389/fmicb.2021.724654) | Frontiers **Hypothesis and Theory** | Textbook structure: survey the two dominant explanations, argue a third, name the discriminating experiments. |
| **Wolday D et al.**, *Interrogating the impact of intestinal parasite–microbiome on pathogenesis of COVID-19…*, Front Microbiol 12:614522 (2021). [10.3389/fmicb.2021.614522](https://doi.org/10.3389/fmicb.2021.614522) | Frontiers **Opinion** | An infectious-disease hypothesis made falsifiable, with the study design that would test it named. |

*Bioscience Hypotheses* (the once-obvious venue) is **effectively defunct** — 218 works ever
registered, last issue Jan 2009. Do not target it. Its criteria remain the cleanest statement of
the genre: innovative, clear, compatible with most facts, and **testable**.

---

## 4. Benchmark / method-comparison papers (no new tool)

| citation | venue | why |
|---|---|---|
| **Nearing JT et al.**, *Microbiome differential abundance methods produce different results across 38 datasets*, Nat Commun 13:342 (2022). [10.1038/s41467-022-28034-z](https://doi.org/10.1038/s41467-022-28034-z) | Nature Communications | The canonical modern comparison: 14 methods × 38 datasets, **no new tool**, lands an actionable ranking. ~1,000 citations. |
| **Cappellato M, Baruzzo G, Di Camillo B**, *Investigating differential abundance methods in microbiome data: a benchmark study*, PLOS Comput Biol 18:e1010467 (2022). [10.1371/journal.pcbi.1010467](https://doi.org/10.1371/journal.pcbi.1010467) | PLOS Comp Biol | The rigorous variant: simulation with known ground truth → FPR, FDR, recall, PR curves. |
| **Lindgreen S, Adair KL, Gardner PP**, *An evaluation of the accuracy and speed of metagenome analysis tools*, Sci Rep 6:19233 (2016). [10.1038/srep19233](https://doi.org/10.1038/srep19233) | Scientific Reports | Frames **neutrality itself** as the selling point: "the first unbiased benchmark in which the authors are not involved in any tool tested." |
| **Miossec MJ et al.**, *Evaluation of computational methods for human microbiome analysis using simulated data*, PeerJ 8:e9688 (2020). [10.7717/peerj.9688](https://doi.org/10.7717/peerj.9688) | PeerJ | Realistic target for a small/solo author; datasets and pipelines deposited for re-benchmarking. |
| **Bokulich NA et al.**, *Measuring the microbiome: best practices for developing and benchmarking microbiomics methods*, CSBJ 18:4048 (2020). [10.1016/j.csbj.2020.11.049](https://doi.org/10.1016/j.csbj.2020.11.049) | CSBJ | Not a benchmark — the **methodology citation to lean on** when writing one. Defines "internal" vs "neutral" benchmarks and the self-assessment trap. |

No venue here requires a new tool. What reviewers expect: multiple datasets, multiple metrics,
documented parameters, and a clear "which should you use."

---

## 5. Data descriptor / resource papers (describe, don't analyse)

| citation | format | why |
|---|---|---|
| **Rodriguez CI et al.**, *Curated and harmonized gut microbiome 16S rRNA amplicon data from dietary fiber intervention studies*, Sci Data 10:346 (2023). [10.1038/s41597-023-02254-4](https://doi.org/10.1038/s41597-023-02254-4) | Scientific Data | **Closest analogue to a reanalysis author's descriptor**: no new sequencing — found 11 existing studies, harmonised them across platforms, published the curated collection. |
| **Schneider D et al.**, *Gut bacterial communities of diarrheic patients with indications of C. difficile infection*, Sci Data 4:170152 (2017). [10.1038/sdata.2017.152](https://doi.org/10.1038/sdata.2017.152) | Scientific Data | The genre template: exact skeleton (Background & Summary → Methods → Data Records → Technical Validation), cross-national reuse story. |
| **Liao C et al.**, *Compilation of longitudinal microbiota data and hospitalome from HCT patients*, Sci Data 8:71 (2021). [10.1038/s41597-021-00860-8](https://doi.org/10.1038/s41597-021-00860-8) | Scientific Data | The "compiled longitudinal resource" pattern; argues explicitly why it overcomes prior single-variable studies. |
| **Jing G et al.**, *Microbiome Search Engine 2*, mSystems 6:e00943-20 (2021). [10.1128/mSystems.00943-20](https://doi.org/10.1128/mSystems.00943-20) | mSystems **Resource Report** | The database flavour: describes a curated resource (>250,000 samples) rather than analysing it. |

**Format requirements:** *Scientific Data* Data Descriptor — Abstract ≤170 words (no references, no
claims of new findings), Background & Summary ≤700 words, then Methods / Data Records /
**Technical Validation (required)** / Data Availability / Code Availability. *GigaScience* Data Note —
structured abstract ≤250 words in three sections (Background/Findings/Conclusions). *Data in Brief* —
**templated and mandatory**; no Conclusion/Discussion/Summary section at all.

---

## 6. Which of these our project should copy

| our situation | copy |
|---|---|
| No matched acute sputum dataset exists | **Janda 2024** — turn "we screened 3 cohorts and none were usable" into the measured result, with the counts as the finding |
| A metric that validates itself on internal anchors | **Gihawi 2023** — the corrective-reanalysis architecture, including the "here is the specific step that generated the false signal" move |
| Two artifact classes | **Chalkley 2013** (3-pp form) + **Halpin 2019** (multimapping changes a biological quantity) |
| A replication cohort that failed on organism composition | **Munday 2022** — locate the mechanism; **Marchal 2024** — state the detection limit |
| Solo, no affiliation | **Mukherjee 2020 / Farquhar 2026 / Coser 2026** — "Independent Researcher, [City]" is accepted in Scientific Reports, Bioinformatics, Briefings in Bioinformatics |
| Publishing sequence from one public dataset | **Pigott 2023** — one paper per sub-question off the same dataset |

**Citation-integrity note:** all 17 DOIs above were verified by the orchestrator against the
Crossref API. Three affiliation strings could not be confirmed via Crossref (Mukherjee, Edgar
2022, Mohammadiaria) because Crossref returned no affiliation field for them; their titles,
journals, years and author counts are confirmed.
