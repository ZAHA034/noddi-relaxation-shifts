# Identifiers and the names used in the article

The article and its supplementary material use descriptive names; the code and the records keep the identifiers below, which this file lists beside those names (third column). Commit hashes refer to the private development repository, whose history is available to editors and reviewers on request. Section, table and figure numbers in the last column are those of the article and its supplementary material.

| Identifier | What it identifies | Name in the article | Where |
|-----|----------------------------------|------------|--------|
| V4 | The real-data test of the frozen judge on ON-Harmony traveling heads; its cells are the acquisition pair × component cells classified as identified at the reference condition | In the main text, "the calibration test of a simulation-trained posterior estimator"; in the supplementary material, "the calibration test" of the simulation-trained estimator | Section 3.4 (without the identifier); S-M7.2 and S-M9; Tables S3, S8 (a) and S12; Fig. S2 |
| Tong T-A | Direction of the center-pair effects in the (x, y) plane, by rotation test | The negative-control direction test | Section 3.4; S-M7.3 |
| Tong T-B | Share of units at retest level for x, y and ODI | The negative-control magnitude test | Section 3.4; S-M7.3 |
| Tong T-C | T-A and T-B combined on x and y (four outcome combinations) | The combined reading | Section 3.4; Table S4 |
| ODI (a) | \|ΔODI\| against the predicted relaxation shift (Kendall τ_b; descriptive) | ODI and the predicted relaxation shift | Section 3.4; S-M7.5 |
| ODI (b) | ODI pair differences against retest variation (descriptive) | ODI against single-scan retest | Section 3.4; S-M7.5; Table S7 |
| Wongkornchaovalit linearity | Share of units linear in TE in the multi-TE series | The TE-linearity check | Section 3.4; S-M7.4; Table S8 |
| C3 arm | The judge given the ROI means of NDI, FISO and ODI at both acquisitions; it gave the verdict | The main input | S-M7.2 and S-M5; `docs/supplement_details.md`, section "From S-M5 (judge study)" |
| C7 arm | The judge given, in addition, mean diffusivity, fractional anisotropy and per-shell mean signal; not used for the verdict | The extended input | S-M7.2 and S-M5; Table S1 (a) |
| L0 (levels L0–L3) | Relaxation levels of the judge's simulator: L0, no between-person variation; L1–L3, healthy variation, diffuse pathology and focal lesion | The relaxation levels no variation (L0, the shared-relaxation scenario), healthy variation (L1), diffuse pathology (L2) and focal lesion (L3) | S-M2 and S-M5 |
| F2 | The AMICO 2.1.1 NODDI fit (intrinsic parallel diffusivity 1.7 μm²/ms, isotropic diffusivity 3.0 μm²/ms), used for all judge computations and real-data analyses | AMICO, the NODDI fit (Section 3.2) | Section 3.2; S-M2 (Fitters) and S-M9.1 |
| F1 | The cuDIMOT NODDI-Watson fit, reported for stages 1–3 only | cuDIMOT | S-M2 (Fitters), S-M3 and S-M4 |
| Stage-2 grid; grid L4 | The 54 relaxation settings of stage 2 (T2,in 70, 82 or 95 ms × r 1.3, 1.6 or 1.8 × six T1 settings; free-water T2 1000 ms); called grid L4 in the FISO checks | The relaxation grid of Section 2.2 (the stage-2 grid); in the FISO checks, the sensitivity grid | Sections 2.2 and 3.4; S-M3 (The stage-2 grid) and S-M7.5 |
| Grid L1 | The stage-3 healthy between-person setting used for the FISO predictions (16 points in each of the three T1 orderings, 48 settings in all) | The primary grid | Section 3.4; S-M2 (L1) and S-M7.5 |
| R1–R26 | The rules of the real-data tests in the fixed document | Rules 1–26 of the fixed analysis document | Sections 3.3 and 3.4; S-M1; S-M7; Table S1 |
| 5f95b18 | Commit of the document fixed before computing (2026-10-28) | The fixed analysis document (2026-10-28) | Section 3.3; S-M7.1; Terms table |
| 08a189c | Commit of the processing and analysis code that the fixed document names, with the SHA-256 of 32 guarded files | The processing and analysis code bound by the fixed analysis document | Section 3.3; S-M1; S-M7.1 |
| 5632626 | Commit recording the verdicts together (2026-11-01); "after the verdict" means after it | The verdict record (2026-11-01) | Section 3.3; S-M1; Terms table |
| 85eb8b0 | Commit freezing the judge's weights (2026-10-25) | "Frozen by hash before validation"; the commit and its date are not printed | S-M5; `docs/fixed_products_and_commits.md` |
| 16793_OXF4GEP001 | The ON-Harmony session that failed twice at the same processing step and was excluded without repair | The first OXF4GEP session of one participant | Section 3.1; Table S1 (c); S-M6 |
| 16794_NOT4GEP001 | The ON-Harmony session whose reversed-phase b0 the documentation and the sidecar describe differently; excluded from the main analyses and included in a sensitivity analysis | One NOT4GEP session | Section 3.1; Tables S1 (c) and S8 (a), (c); S-M7.2 and S-M7.5 |
| osf.io/5juwv | Registration, on the Open Science Framework, of the earlier study of the same ON-Harmony data | The earlier study; the registration is given in Section 3.5 (Data and code availability) | Sections 3.3 and 3.5; S-M1; Table S2 |
| T-in, T-app, T-none | Three conventions for how relaxation enters the tortuosity relation of the signal generator: the extra-neurite perpendicular diffusivity is d∥(1 − ν) with ν the relaxation-free neurite fraction (T-in); d∥(1 − ν_app) with ν_app the apparent fraction for the acquisition and the person's relaxation (T-app); or d∥(1 − 0.60) (T-none) | True-fraction, apparent-fraction and fixed tortuosity (supplementary material) | S-M2 (generator assumptions); S-M9.4 |
| Offset-subtracted retest share | For a comparison of two acquisitions or centers, the share of units whose difference minus a common shift (the median of the differences of the other units of the same participant and pair) is within 1.96·√2·s_u (the c10 retest SD and t₀.₉₇₅(df) for the ten centers); read as at retest level when at least 0.90 of the units are | The retest rule of Section 3.4 applied after one common shift is subtracted | Section 3.3 (last paragraph); S-M9.2; Table S12 (a)–(d) |
| ΔTE* | The TE difference at which one common shift stops sufficing: the TE difference at which the expected offset-subtracted share equals 0.90, given the between-unit SD τ of the TE slopes and the retest SD | Not named in the main text | S-M9.2; S-M9.4 (E1); Table S12 (b); Fig. S3 (a) |
| S1a | Environment analysis: the stored stage-3 fits of 210 reference tissues (central relaxation and 16 healthy between-person settings, three generators, four pair classes: TE only, TR only, direction scheme only, all three), with s the median OXF1PRI retest SD of x | The reference-tissue analysis; heading "Fitter output with reference tissues" in S-M9.3 | S-M9.3 and S-M9.4; Table S12 (a); Fig. S3 (c) |
| S1b | Environment analysis: observed ROI values (multi-TE series at TE 62 ms; ON-Harmony OXF1PRI medians) mapped to true fractions at the central relaxation and fitted at every acquisition of the dataset, giving the fitter-only TE-slope SD and the offset-subtracted share of the 28 ON-Harmony acquisition pairs | The ROI-matched analysis; heading "Fitter output with ROI-matched tissue states" in S-M9.3 | S-M9.3 and S-M9.4; Fig. S3 (a) |
| S2 | Environment analysis: Proposition 4 at fitter output, with 27 ROI-matched tissue states, simulated controls and patients, sites at TE 70 and 100 ms or at the eight ON-Harmony acquisitions, and several site-term models | The Proposition 4 analysis; heading "Proposition 4 at fitter output" in S-M9.3 and S-M9.4 | S-M9.3 and S-M9.4; Fig. S3 (d) |
| K8 | Environment analysis: the fixed FISO sign rule (S-M7.5) with the algebraic Δ logit FISO replaced by noise-free fitted values | Heading "FISO signs with the fitter" in S-M9.3 and S-M9.4 | S-M9.3 and S-M9.4 |
| E1 (environment) | TE slopes, τ, β and ΔTE* of the multi-TE series | Heading "Multi-TE series" in S-M9.4 | S-M9.3 and S-M9.4; Table S12 (b); Fig. S3 (a) |
| E2 (environment) | Offset-subtracted shares of the 19 ON-Harmony scanner pairs | Heading "Scales" in S-M9.4 | S-M9.3 and S-M9.4; Table S12 (c) |
| E3 (environment) | The offset-subtracted shares across the ten centers | As for E2 | S-M9.3 and S-M9.4; Table S12 (c) |
| E4 (environment) | Scale comparison in the multi-TE series (62 against 132 ms as primary, the other 27 TE pairs), ON-Harmony and the ten centers | As for E2 | S-M9.3 and S-M9.4; Table S12 (c); Fig. S3 (b) |
| E5 (environment) | Slope structure: how the TE slopes of x depend on the level of x, and the shares expected from it at 62 against 132 ms | The slope-structure analysis; heading "Slope structure" in S-M9.5 | S-M9.5; Table S12 (d) |
| E6 (environment) | The scale comparison with shared relaxation, from the stored noise-free fits of S1b | The shared-relaxation analysis; heading "Shared relaxation" in S-M9.5 | S-M9.5; Table S12 (d) |
| S3 | Environment analysis: the raw-scale site × group term, by algebra only, for two groups that share relaxation and differ only in true NDI | The raw-scale analysis; heading "Raw-scale site × group term" in S-M9.5 | S-M9.5; Table S12 (d) |
| judge | The simulation-trained neural posterior estimator of this study (code, records and figure scripts) | The simulation-trained estimator | S-M5; Fig. S2 |
| cards | Per-acquisition-pair summaries of what the algebra predicts from the acquisition records | Algebraic pair summaries | S-M3 |
| X9; X12 | Numbers of two analyses in the stage specifications | Fitter fidelity; the ComBat example | S-M3; S-M4 |
| SG1–SG7; H1–H4; HPI | Criteria of stages 1, 2 and 3 in their specifications | Described by content in S-M4 | S-M4 |
| Q1–Q4 | The four acquisition pairs of the stage-2 grid (A–Q1 TE only, A–Q4 TR only, A–Q3 direction scheme only, A–Q2 all three) | Described by content | S-M3; S-M9.3 |
| B-iii | The retest SD of the negative control from the three c10 sessions | The c10 retest SD | S-M7.3 |
| D1–D5 (plan of 2026-11-01) | Analyses of the dated plan after the verdict (D1, bias and overconfidence; D3, joint multi-TE fit; D5, unit-level outliers) | Described by content | Tables S1 (b) and S5 (a) |
| D1–D3 (displays) | Displays after the verdict | Displays 1–3 | Table S5 (b) |
| S73 plan | The plan of the post hoc analyses committed on 2026-11-16 | The plan of 2026-11-16 | S-M9.6; Table S5 |
| D0, D1, D2 (design simulation) | Designs of the multi-site simulation: one central TE; the eight ON-Harmony acquisitions; the same TEs with short and long TEs within each vendor | One TE; ON-Harmony; Spread | S-M11; Table S13; Fig. 4 |
| S0, S1b, S1e, S2b, S2e, S2i, S3b, S3e | Scenarios of the multi-site simulation | None; ν, both; ν, extra; T2, both; T2, extra; T2, intra; ν + T2, both; ν + T2, extra | S-M11; Table S13 |
| M1, M1c, M1r, M2, M3, M4, M5a, M5b, M6a, M6b | Analyses of the multi-site simulation | ComBat; ComBat on controls; ComBat on raw NDI; site term; group × TE; meta-regression; measured relaxation without error; measured relaxation (10 % error); second TE in every participant; second TE in half | S-M11; Table S13; Fig. 4 |
| R-a, R-b, R-c, R-d | Readings of the multi-site simulation fixed before computing | The coverage, extrapolation, spread and vendor readings | S-M11; Table S5 |

In the code of the separation attempts, C3, C5 and C7 also name the coordinate sets (the main set, a set that adds mean diffusivity and fractional anisotropy, and the extended set; S-M4). The relaxation levels L0–L3 of the estimator's simulator are not the grids L1 and L4 of the FISO checks (S-M7.5).

The label E is used in five unrelated series:

- E1 and E2 of Supplementary Methods S-M8 (definitions file of 2026-10-31);
- E4 in `decomp/manuscript/figures/src/c1_split_edge.py` (rules record of 2026-11-07);
- the eddy features E1–E6 of the exploratory analysis D5 (`decomp/explore/d5/d5_features.py`);
- the simulation runs and analyses of stages 1–3 of the signal generator: runs E1, E1-TR, E2a, E2b, E2c, E3, E5, E6 and E6n (`decomp/stage1/config.py`) and analysis E4 (`decomp/stage1/analyze.py`, output `decomp/results/stage1/e4.json`). The Zenodo file `stage3_fits_E6.zip` holds the stored fits of run E6, which environment analysis S1a reads; it is not analysis E6 of the environment;
- E1–E6 of the measurement-function environment (S-M9), the only series listed in the table above.

## Labels in the translated records

The dated records in `docs/` (the fixed analysis document, the judge specification, the verdict record and the plan of the analyses after the verdict) are translations of working documents and keep their working labels. The labels mean:

| Label | Meaning |
|---|---|
| A | The analysis session, which wrote the code and the records and ran the analyses |
| Review session (reviewer) | The separate session that reviewed the plans, code and drafts before each step; "relays" are its messages |
| Author, author decision | The first author, who made the decisions recorded as such |
| Relay ZZZ-8, ZZZi, ZZZj-3, ZZZZa … ZZZZm, … | Names of the relays (review messages), in the order they were sent; a number after a hyphen numbers an item of that relay (for example ZZZZi-4 is item 4 of relay ZZZZi). "Receipt" is A's confirmation that a relay was received |
| A's proposal a1–a16 | Numbered proposals of A in the fixed analysis document, each accepted, changed or rejected as the record states |
| (ga), (na), (da), (ra) | The options (a), (b), (c) and (d) of a decision, named with the first letters of the Korean alphabet (가, 나, 다, 라); some translations write them as (a), (b) or as (A), (B) |
| ①, ②, ③ … | Item numbers inside a sentence, a list or a table cell |
| ⑨, ⑩, ⑪ as names of steps | Steps of the work plan of the measurement: ⑨ preparing the new measurement (code and loaders); ⑩ fixing the analysis document, that is, committing it with the author's permission ("⑩ permission"), after which the loader guard lets the code read images; ⑪ the measurement itself, in which every read of image data is recorded. "⑪ pre-check" is A's full re-review made before ⑩ |
| ⟦…⟧ | Working marks of the original documents, kept unchanged |
