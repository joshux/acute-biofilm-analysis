# Paper-type decision — independent evaluation

**Date:** 2026-10-10. Question: is "methods-forward short report" the appropriate paper type?
Answer: **the methods-forward frame is right, but "short report" is wrong, and the venue changes.**

---

## 1. The finding that forces this re-evaluation

An independent prior-art audit (three citation indexes + full-text read of Rossi 2018) returned:

| claim | verdict |
|---|---|
| "First transcript-level test of Kolpen" | **supported by absence** — no citing paper and no independent study applies metatranscriptomics/RNA-seq/gene expression to Kolpen's model; Kolpen's own 2022 *Biofilm* review names the EPS-expression question as explicitly unresolved |
| The **biological** two-axis result | **confirmatory** — a re-test of a published model using borrowed proxies (Gifford %RP; alginate/psl/pel panels) |
| The **MAPQ multi-genome artifact** | **genuinely unclaimed** — abundant general multi-mapping literature exists (Halpin 2019, Deschamps-Francoeur 2020, Zhao 2023) but **nothing identifies MAPQ/multi-mapping as a growth-rate-*inverting* bias in multi-genome metatranscriptomics** |
| No prior study measures growth × matrix from the same respiratory metatranscriptome | confirmed |
| %RP has never been applied to a human clinical respiratory sample | no evidence found (all %RP literature is marine/soil/groundwater) |

**Consequence:** the paper's defensible novelty is **methodological**. The biology is the
*demonstration*, not the contribution.

## 2. Why that makes the current plan's paper type wrong

A brief communication / Observation is **defined by carrying a concise biological result**. Its
editorial criteria say so: npj — "a concise study of high quality and broad interest"; ASM
Observations — "short descriptions of research results of exceptional importance."

If we submit in that type, the required framing puts our **weakest** leg first. A reviewer reads a
1,200-word biological brief and asks "what does this add over Kolpen plus Gifford?" — and the
honest answer ("a confirmation") is not a brief-communication-grade claim.

The paper type must follow the novelty. Ours is methodological → the type must be one whose
**defining criteria are methodological**.

## 3. Verified format constraints (from the publishers' own pages)

| type | hard constraints | framing it imposes |
|---|---|---|
| **mSystems "Methods and Protocols"** | 5,000 words (excl. refs/tables/legends); *"must include validation of, or application to, a relevant and important question in microbial cell biology or ecology and provide results demonstrating its performance in comparison to existing state-of-the-art techniques"* | **methodological by definition** — this sentence describes our paper's structure |
| mSphere / mSystems "Observations" | 1,200 words; ≤2 figs; ≤25 refs; abstract ≤250 + Importance ≤150 | biological result, "exceptional importance" |
| npj Biofilms "Brief Communication" | abstract **≤70 words**; main text 1,000–1,500; ~20 refs | "concise study… of broad interest" — biology-first |
| GigaScience "Technical Note" | no word limit; **must be an open-source software tool or method**, OSI-approved licence, test data + expected outputs, reproducible examples, "rejected without review" if not met | software/tool identity |
| GigaScience "Brief Communication" | 1,500 words; **≤10 references**; ≤2 figs | methods-forward, but 10 refs is unworkable (see §5) |
| Frontiers "Brief Research Report" | ≤4,000 words; ≤4 figs/tables | neutral; viable fallback |

## 4. Ruling

**Primary target: mSystems "Methods and Protocols" (5,000 words).**

- Its *required* elements map one-to-one onto our evidence: validation (lab anchors: exp 13.6% vs
  stat 2.1%, p<1e-4; %TF independently 6.42% vs 2.42%, p=4.9e-5), application to a relevant
  question (Kolpen's model), and comparison to existing state-of-the-art (the naive multi-genome
  panel, which is precisely the artifact).
- 5,000 words carries **both** artifact classes, the validated metric, the two-axis result, and the
  five negative results — which a 1,200-word Observation cannot.
- Removes the framing mismatch: the method is the contribution by editorial definition.

**Runner-up, if a short high-visibility piece is preferred: mSphere / mSystems "Observations"
(1,200 words).** It has a clean verified track record of method-only Observations (Schloss 2021
ASV artifact; Brennan 2024 contamination mitigation) and the tightest true "short report" format.
Cost: it forces the artifact-only frame and puts all negatives in the supplement.

**Retract the earlier npj Biofilms primary recommendation.** Reason: the Brief Communication
abstract cap (70 words) plus 1,000–1,500 words cannot carry two artifacts *and* five negative
results, and its "broad interest" framing is biology-first — exactly the mismatch in §2. npj
remains a fallback only if we later decide to lead with the biology and drop the artifact material.

**Reject GigaScience Technical Note for this paper.** Verified: it requires an open-source tool with
an OSI-approved licence, test data, expected outputs, and reproducible examples — "submissions that
do not meet these requirements will be rejected without review." The repo has **no LICENSE file, no
packaging, no tests, and no tool**; our contribution is an analysis *design*, not software. Right
genre, wrong artifact. (Would become viable if we shipped the quantifiers as a packaged tool.)

**Fallback ladder:** mSystems M&P → mSphere Observations (if shortening) → Frontiers Brief Research
Report (≤4,000 words, ≤4 figs) → npj Biofilms Brief Communication (biology-led rewrite).

## 5. What would change this ruling

- **If the artifact is independently reproduced on a second dataset**, it becomes a standalone
  methods finding and GigaScience/mSphere-Observations rises — a reproduced artifact is a
  shorter, sharper paper than artifact-plus-biology.
- **If the two-axis biology is replicated in a second cohort**, the biological claim strengthens and
  a full Research Article becomes appropriate (currently it cannot: 5 CF patients, cross-specimen,
  one organism).
- **If the 10-reference cap were workable**, GigaScience Brief Communication would be the best
  short methods-forward home. It is not: we need ≥13 (Gifford ×2, Kolpen ×2, Rossi, Zhang, Halpin,
  Chalkley, Salter, Vandeputte, Long, Walter, Kopf).

## 6. Caveats

- Format constraints above are quoted from the publishers' current pages (verified 2026-10-10);
  re-verify at submission, these change.
- Editorial *bar* (acceptance likelihood) cannot be read off guidelines. The fit judgements here are
  reasoned from the stated criteria, not from insider knowledge.
- The novelty claim rests on absence across three citation indexes; it is well-supported but is not
  proof of a negative for unindexed or very recent work.
