# Multi-site design simulation — analysis plan, fixed before computing (English translation)

> English translation of `decomp/design_sim/ANALYSIS_PLAN.md` (SHA-256 9bc6c2cbcd7557650365e169994b7d22d170ad0c4da24437ae20c901674e9171), checked by script for numbers and identifiers. The Korean file is the record; where the two differ, the Korean file governs.

Written 2026-11-17. In the grades of the manuscript this is **post hoc** (defined after the verdict and fixed before computing). All inputs are assumed values, so the results are used in the manuscript together with «illustration». No verdict is reopened. No real data are read (a simulation at the level of the signal model). The main run is made only after this document is committed to the repository. If anything must change before the run, the place and the reason are written in «Change record» at the end of this document and it is committed again. Nothing is changed after the run. A pilot run to measure time is made, but its metric values are not opened.

Basis: Section 6 of review direction S81 and Section 2 of S82 (proposed by the review session, approved by the author), and four amendments by the analysis session (S1 effect size · Monte Carlo SE per replicate · record of the K-d setting · ComBat implementation and M4).

## 1. Purpose

The applied reader's question: when the compartmental relaxation of the groups differs and the TE of the sites differs, how wrong is the estimate of a fraction-only group difference by a standard ComBat analysis, are its confidence intervals right, and what does a group × TE term or a measurement of relaxation change? It adds what the simulation of Proposition 4 in the manuscript (Supplement S-M9) lacked: variation between people · empirical Bayes (EB) ComBat · unbalanced sites. Most of the results are fixed by algebra (Proposition 4, Eq. (7), Section 2.4); this simulation translates their size into the language of confidence intervals and tests.

## 2. Data-generating model (logit scale · person i, site s, ROI unit r)

- True fraction: logit ν_isr = logit ν_r + θ·g_i + b_ir. g_i = 1 patient, 0 control. Each person is scanned at one site only.
- ROI: 27 units. ν_r equally spaced from 0.45 to 0.80 (assumed; no effect on analyses on the logit scale, used only for the raw-scale analysis M1r).
- Variation of the true fraction between people: b_ir = 0.08·(√ρ·z_i + √(1 − ρ)·z_ir), z ~ N(0, 1), ρ = 0.5 (correlation between ROIs, assumed). The between-person SD of 0.08 logit is the same value as in Supplement S-M9.7 (assumed).
- Relaxation (controls): T2,in 82.5 ms, T2,en 52 ms (Section 2.4). Patients: T2 of the assigned compartment × 1.0984 (the diffuse lengthening of Section 2.4, Brandhofe et al., 2021).
- Variation of relaxation between people (factor V): for each person a lognormal multiplier on T2,in and T2,en, coefficient of variation v (log SD √ln(1 + v²)), in–en correlation 0.8 (assumed). One per person (common to all ROIs).
- T1 · proton density: common to the groups (the condition of Proposition 4). The two tissue compartments share T1 (the central value of S-M2), so the TR term of logit NDI is 0, and TR is only recorded.
- Apparent: logit ν_app,isr = logit ν_isr + TE_s·(R2,en,i − R2,in,i), R2 = 1/T2 (ms⁻¹).
- Observed: y_isr = logit ν_app,isr + γ_sr + δ_sr·ε_isr + η_s·g_i.
  - ε ~ N(0, 0.031²): the retest SD of x for OXF1PRI (manuscript 4.1).
  - γ_sr = vendor mean (Siemens 0, Philips +0.10, GE −0.10 logit; assumed) + site · ROI deviation N(0, 0.05²) (assumed; drawn anew in each replicate).
  - δ_sr = exp(N(0, 0.1²)) (assumed; in each replicate).
  - η_s: a site × group term unrelated to relaxation (factor H).
- No fitter is used (signal-model level). Fitter effects (1.00–1.30 times in Supplement S-M9) and the noise floor are not included, and are written as limits in the results paragraph.

## 3. Design factors (cells = D × S × V × H × N = 3 × 8 × 2 × 3 × 2 = 288)

Acquisitions (the 8 ON-Harmony acquisitions, card table of the fixed analysis document; TE is the middle value of the range, as in Supplement S-M9.7):

| Acquisition | Vendor | TE (ms) | TR (s) |
|---|---|---|---|
| OXF1PRI | Siemens | 92.0 | 3.60 |
| OXF2PRI | Siemens | 92.0 | 3.60 |
| OXF3TRI | Siemens | 96.4 | 3.60 |
| NOT2ING-89 | Philips | 89.0 | 4.00–4.38 |
| NOT2ING-98 | Philips | 98.0 | 4.00–4.38 |
| NOT4GEP | GE | 68.0 | 4.90 |
| OXF4GEP | GE | 69.65 | 8.00 |
| NOT1ACH | Philips | 70.6 | 9.81–10.62 |

- D site design: D1 the table above as it is (long TEs at Siemens and NOT2ING). D2 the same eight TEs rearranged so that each vendor has both short and long TEs — Siemens sites (OXF1PRI, OXF2PRI, OXF3TRI) 68.0 · 92.0 · 98.0, Philips sites (NOT2ING-89, NOT2ING-98, NOT1ACH) 69.65 · 89.0 · 96.4, GE sites (NOT4GEP, OXF4GEP) 70.6 · 92.0. D0 all eight sites at TE_c = the mean of the eight TEs, 84.45625 ms (control).
- S scenarios (ΔR2 = patients − controls, ms⁻¹; Δ = ΔR2,en − ΔR2,in):
  - S0 no difference (θ = 0, same relaxation).
  - S2b relaxation only, both tissue compartments lengthened (Δ_b = −0.000636909 ms⁻¹, the slope of Supplement S-M9.7). S2e relaxation only, extra-neurite only (Δ_e = −0.001722786). S2i relaxation only, intra-neurite only (Δ_i = +0.001085877). θ = 0.
  - S1b fraction only, θ = TE_c·Δ_b (the same apparent difference as S2b at the central TE, about −0.0538 logit). S1e fraction only, θ = TE_c·Δ_e (about −0.1455 logit).
  - S3b = the fraction of S1b + the relaxation of S2b. S3e = the fraction of S1e + the relaxation of S2e.
  - The size of S1 is matched to the apparent difference rather than to values reported for some disease in order to show side by side two causes that cannot be told apart at one site.
- V: v = 0 or 0.10.
- H: H0 none. H1 η_s ~ N(0, 0.02²) (in each replicate). H2 +0.02 logit for patients at Siemens sites (a term tied to the vendor).
- N: balanced — 20 per group at each site (40 per site, 320 in all). Unbalanced — 40 per site, with patients OXF1PRI 12 · OXF2PRI 28 · OXF3TRI 16 · NOT2ING-89 24 · NOT2ING-98 20 · NOT4GEP 20 · OXF4GEP 28 · NOT1ACH 12 (proportion 0.3–0.7). In D2 · D0 the numbers attached to the site names are used as well.
- Replicates: 1,000 per cell. Seeds: numpy SeedSequence([20261117, cell number]) — the cell number is the lexicographic order of the order above (D, S, V, H, N).

## 4. Analyses (separately for each ROI unit; 95 % intervals of the group coefficient from the t distribution)

| ID | Analysis | Reported |
|---|---|---|
| M0 | Site ignored: y ~ group (OLS) | θ̂ |
| M1 | Standard ComBat: neuroCombat 0.2.12 (EB, parametric, group kept as a categorical covariate, the 27 ROIs as features), logit scale → OLS ~ group | θ̂ |
| M1c | ComBat trained on controls only (the recommendation of Jodoin et al., 2025): neuroCombat (EB, no covariate) trained on controls only, and its estimates applied in the same way to all participants (`combat_apply` in `decomp/rules_env/analyses.py` of the repository, the same formula as neuroCombatFromTraining 0.2.12) → OLS ~ group | θ̂ |
| M1r | The same ComBat as M1 on raw-scale NDI → OLS ~ group | raw-scale θ̂ |
| M2 | y ~ group + site (fixed) | θ̂ |
| M3 | y ~ group + site + group × TE_s (TE not centred) | θ̂(TE→0) = group coefficient · θ̂(TE_c) · test of group × TE |
| M4 | Random-effects meta-regression in 2 stages: site-wise group differences d_s (OLS within site) regressed on TE_s, between-site variance τ² by the method of moments (the regression version of DerSimonian–Laird), weights 1/(v_s + τ̂²), intercept = θ̂(TE→0), interval from t(number of sites − 2) | θ̂(TE→0) · test of the slope · share with τ̂² = 0 |
| M5a · M5b | Measured relaxation: R2,in · R2,en measured for each person, the relaxation term of the equation TE_s·(R̂2,en − R̂2,in) subtracted, then M2. M5a no measurement error (reference), M5b independent lognormal error with CV 10 % on each R2 (assumed) | θ̂ |

- M4 replaces the mixed model (REML) of S81: statsmodels is not in the locked environment, and adding a package would change the deposited lock file. The aim of treating site × group variation as a random effect is the same. It is light to compute and is run in every cell.
- M5 is a reference analysis that shows the upper limit and a realistic value of «if relaxation is measured». The manuscript text does not call this analysis restore · remove · correct.

## 5. Estimand and metrics

- Estimand θ: the true fraction difference between the groups (logit, no relaxation). In S0 · S2, θ = 0. The target of M1r is, for each ROI, the true difference on the raw scale E_b[expit(logit ν_r + θ + b) − expit(logit ν_r + b)], b ~ N(0, 0.08²) (Gauss–Hermite, 40 points).
- Metrics (for each analysis · cell): bias, empirical SD, mean of the model SE, RMSE, coverage of 95 % intervals, rejection rate of H0: θ = 0 (two-sided 0.05; false positives in S0 · S2, power in S1 · S3), rejection rate of group × TE (slope) for M3 · M4, share of τ̂² = 0 for M4.
- Summary over ROIs: in each replicate the values over the 27 ROIs (error · covered or not · rejected or not) are averaged, and then averaged over replicates. The Monte Carlo SE is the SD of the per-replicate values / √1000 (reflecting the correlation between ROIs). The empirical SD is computed for each ROI over replicates and averaged. Values by ROI are kept in the deposited result files.

## 6. Readings fixed before the run (written in this frame whichever way the results go)

- R-a (D1 · H0 · V0 · balanced): in each of S2e and S2b, if the coverage of M1 is below 0.90, «a standard ComBat analysis in this design does not give valid intervals for the fraction-only difference»; if 0.90 or above, «the within-site bias was small against sampling error at this size».
- R-b (same cell): for each of b and e, if the coverage of M3 θ̂(TE→0) in S2x is 0.93 or above and its power in S1x is below 0.50, «extrapolation to TE → 0 is valid but imprecise in this design». Otherwise the values are reported as they are.
- R-c: the difference in SE and power of M3 θ̂(TE→0) between D2 and D1 is reported as design advice («TE spread within vendor»). Descriptive only.
- R-d: in H2 · S0 (D1 · V0 · balanced) the rejection rate of group × TE of M3 (false positives arising from mixed vendors) is reported. Descriptive only.
- Printing: the supplementary table carries the bias · coverage · rejection rate of M1 · M1c · M2 · M3 · M4 · M5b in the cells balanced · H0 · H2 · all D · V · S (96 cells), and all 288 cells (all analyses · metrics · values by ROI) are kept as deposited result files and tables (S81 proposed all of them as a supplementary table — divided because of the length of the supplement).

## 7. Checks that must be able to fail (run before the main run, recording the output and the failure of the planted faults)

- K-a algebra without noise: with ε, b, v, η, γ set to 0 and δ = 1, in each cell the site-wise difference of group means must equal θ + TE_s·Δ within 1e−10 (all of S1 · S2 · S3, D0 · D1 · D2). Planted fault: drop the TE term from the generator → must fail.
- K-b null calibration: in D0 · S0 · H0 · V0 · balanced, the coverage of M2 within 0.929–0.971 (0.95 ± 3 MC SE). Planted fault: halve the SE → must fail. The coverage of M1 in the same cell is reported as a result, not as a check.
- K-c analytic bias: in D0 · S2e · H0 · V0 · balanced, the bias of M2 within 3 MC SE of TE_c·Δ_e. Control: shortening T2 in the same proportion instead of lengthening it must change the sign of the bias (the same check against the original expectation must fail).
- K-d consistency with existing values: in the setting of Supplement S-M9.7 — D1 · H0 · V0 · balanced, ε = 0 · δ = 1 · site deviation of γ 0 (between-person SD 0.08 only), the middle values of the eight TEs — the rejection rate of group × TE of M3 within 3 MC SE of 0.136 in S2b and 0.634 in S2e (the stored z-test values, `ledger/derived/s75_values.json` power.by_te.mid). If not reproduced, stop and report. The code is not changed to make them agree (if the cause is on the side of the existing calculation, that fact is written down). Planted fault: effect set to 0 → the rejection rate falls to near 0.05 and must fail.
- K-e ComBat implementation: training on all participants rather than controls only and applying with `combat_apply` to the same data must equal the output of neuroCombat 0.2.12 within 1e−10 (training without covariates). Planted fault: EB shrinkage switched off (γ̂, δ̂ instead of γ*, δ*) → must fail.
- The output of the five checks (pass · fail lines, failure under the planted faults) is kept in `checks.json` in the result folder.

## 8. Computation

- CPU only (no GPU). At most 8 concurrent processes. Other users' processes and folders are not touched. PYTHONDONTWRITEBYTECODE=1.
- The locked environment of the repository (numpy · scipy · pandas · neuroCombat 0.2.12). Code: `decomp/design_sim/`. Results: `decomp/results/design_sim/` (summary by cell · values by ROI · checks · run time).
- Time is measured with a pilot run (metric values not opened).

## 9. Place in the manuscript

A new supplementary section (design · checks · summary of results) and a new supplementary table. In the main text, the power paragraph of 5.2 is replaced by a paragraph of about 120 words (the frame of R-a – R-d, grade «post hoc simulation; assumed values»), one sentence is added at the end of the model paragraph of 5.4, and the basis column of the last row of Table 2 and Box 2 (design advice) reflect it. No figure is added.

## Change record

- 2026-11-17, before the main run (after the first commit of this plan, before the metrics of any cell were computed): «Printing» of Section 6 is changed. The supplementary table carries only the cells used by R-a – R-d and the 24 cells D0 · D1 · D2 × eight scenarios × H0 · V0 · balanced, with the bias · coverage · rejection rate (12 columns) of the analyses M1 · M2 · M3 (TE → 0) · M5b. All 288 cells (all analyses · metrics · values by ROI) are kept as deposited result files and as tables in a deposited document. Reason: a supplementary table of 96 cells × 18 columns would run to several pages, lengthening again the supplement that was shortened the same day. The analyses · cells · metrics · readings · checks are not changed.
- 2026-11-17, after the first main run stopped (before the metrics of any cell were opened): the first main run (code 9b2945c) stopped at the D0 cells. In D0 the TE of the eight sites is the same, so the group × TE of M3 and the TE regression of M4 are not defined (singular design matrix). The code is changed so that M3 · M4 are not computed in D0 and are recorded as «not defined», and all 288 cells are run again with the same seeds. The checks K-0 · K-a – K-e and the readings R-a – R-d do not use M3 · M4 of D0. The 192 cells (D1 · D2) left by the stopped run were kept unopened, and are only compared by machine for equal metrics with the results of the run made again.
