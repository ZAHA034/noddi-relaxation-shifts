# Decision on whether to run E1 (age × scanner TE) — Section 6a of the plan · 2026-11-18 (English translation)

> English translation of `decomp/age_te/E1_GO_DECISION.md` (SHA-256 b6f49109eba8b5a4b00b0508d7e92a6f100367936a7c84127240abe646fde662), checked by script for numbers and identifiers. The Korean file is the record; where the two differ, the Korean file governs.

**Decision: do not run.** The participant information file was not opened. E1 is written as «future work» in the manuscript.

The rule is in Section 6a of the plan `ANALYSIS_PLAN.md`, and was fixed before the literature search (record outside the repository: 2026-11-18 18:22:45 +0900, plan sha256 8b19e79e…). The gist of the rule: the largest value (L_max, ms/10 years) among Q1 literature values that report as numbers the age-related change of white-matter T2 over the age range of ON-Harmony (roughly 20–60 years) is converted generously. If that optimistic power is 0.5 or more (L_max about 0.64 ms/10 years or more), run. If there is no usable literature value, do not run.

## Basis

This decision came not from «the literature value is below the threshold» but from «not a single literature value in numbers could be read». The last line of the rule (if there is no usable value, do not run) applied. Among what was read, those that stated a direction without numbers are the following.

| Study | Extent read | Population · age | What it says |
|---|---|---|---|
| Hagiwara et al., 2021, Invest Radiol (doi:10.1097/RLI.0000000000000720) | full text (PMC7864648) · supplementary material (the formula of the T2 curve, «Supplemental Digital Content 2») could not be opened | 114 healthy adults, 21–86 years, 3 T synthetic MRI | «Stable T1 and T2 values from the 20s until 60s in both GM and WM and increases from 60s» |
| Papadaki et al., 2019, J Magn Reson Imaging (doi:10.1002/jmri.26707) | abstract only | 90 healthy adults, 22–81 years, 1.5 T multi-echo spin echo | a significant increase of long T2 «first noted at 60–69 years» |
| Canales-Rodríguez et al., 2021, Neurobiol Aging (doi:10.1016/j.neurobiolaging.2021.06.002) | abstract only | 145 healthy adults, 18–60 years, multi-component T2 | the T2 of intra- and extracellular water decreases with age, more steeply in men — the size is not in the abstract |
| Knight et al., 2016, Phys Med Biol (doi:10.1088/0031-9155/61/15/5587) | full text (PMC5390949) | 37 controls, 49–87 years | white-matter T2 increases with age. The age range is outside the rule, and the slope values are only in a figure |

Two studies (Hagiwara, Papadaki) say that white-matter T2 changes little at 20–60 years. This is on the same side as «do not run». Canales-Rodríguez says that in this range the T2 of intra- and extracellular water decreases, but its size could not be read. This is the uncertainty that remains in this decision. If those numbers (or the supplementary formula of Hagiwara, or the curve coefficients of Hasan et al., 2010) are obtained, they can be entered in the same rule and the decision made again. The rule is not changed.

## Search record

PubMed (MCP) searches, 5 times, 2026-11-18. All ran to the end (Stopped: no).

1. `(T2 relaxation OR T2 relaxometry OR T2 mapping) AND white matter AND (age OR aging OR lifespan) AND (healthy adults OR adult lifespan)` — of 158 records, a list of 25 was viewed and the bibliographic records · abstracts of 15 were opened.
2. `(white matter[Title/Abstract]) AND (T2 relaxation time[Title/Abstract] OR T2 relaxometry[Title/Abstract] OR quantitative T2[Title/Abstract] OR T2 values[Title/Abstract]) AND (age-related[Title/Abstract] OR aging[Title/Abstract] OR lifespan[Title/Abstract]) AND healthy[Title/Abstract]` — 13 records, the bibliographic records · abstracts of all opened.
3. `Bouhrara[Author] AND (relaxometry OR T2) AND (age OR aging OR lifespan)` — 39 records. Only the list was obtained; none opened.
4. `(intra- and extracellular water OR long T2 OR "T2 relaxation") AND white matter AND (age OR aging) AND (myelin water fraction) AND cross-sectional AND adults` — 3 records (already opened).
5. `(synthetic MRI OR QRAPMASTER OR MR fingerprinting OR quantitative MRI) AND (T2) AND white matter AND (age-related OR aging) AND (healthy adults OR normal) AND (relaxation times)` — 18 records, abstracts of 11 opened.

Full texts were opened for two: PMC5390949 (Knight 2016) and PMC7864648 (Hagiwara 2021). PMC8333472 (Bartzokis 1994) had an empty full text.

**Blocked (web fetch, same day).**
- Hagiwara 2021, publisher side: HTTP 402. The supplementary material could not be opened.
- Canales-Rodríguez 2021: the publisher link returned only «Redirecting».
- Hasan et al., 2010, Magn Reson Med (doi:10.1002/mrm.22515; 130 people aged 15–59 years, a U-shaped curve of white-matter T2): HTTP 403.
- Papadaki 2019 · Canales-Rodríguez 2021 · Hasan 2010 have no PMC copy.

**Not searched.** Scopus, Web of Science, Google Scholar, conference abstracts, archives (not used by the rule).

## What goes in the manuscript (supplement, one sentence in the main text)

Measuring in real data the term d(R2,en − R2,in)/dz of the continuous-covariate corollary is written as future work. In a travelling-heads design like ON-Harmony, even with an age SD of 10 years and 110 observations, the minimum detectable effect at 80 % is about 2–4 ms per 10 years in extra-neurite T2 (estimate of Section 7 of the plan). Two studies say that white-matter T2 changes little at 20–60 years. One study says that in this range the T2 of intra- and extracellular water decreases, but its size could not be read. So a wider age range and a design that also measures compartmental T2 are needed. Literature citations follow the reference rules of the manuscript (Q1 · no archives). The journal grades (public SJR export copy, 2025 values): Invest Radiol · J Magn Reson Imaging · Neurobiol Aging are Q1, Phys Med Biol is Q2 (Q1 in 2024). The manuscript cites only Hagiwara et al., 2021, whose full text was read directly.

## Addendum (after the decision, same day)

The age summary of Section 7 of the plan (31.5 ± 10.5 years, 7 female) is a value carried over from an earlier work record. It was not checked again this time against the original text of the data paper (the review session gave sex as «10 versus 10»). So the manuscript does not use this summary and writes only «assuming an age SD of about 10 years». The decision does not depend on this value (not a single literature value of the rule could be read).
