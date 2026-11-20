# Code guide

One line per Python file; the tests are listed last. Files marked «locked» are byte-identical to the versions recorded in the article, so their comments stay in Korean; the English for their messages and keys is in `docs/glossary.md`. The run record in `records/` is kept as it was run, with Korean comments. In the other files, Korean comments and docstrings were replaced by English.

| File | What it does |
|---|---|
| `decomp/age_te/__init__.py` | age_te — covariate (age · sex) × scanner TE interaction (plan ANALYSIS_PLAN.md · post hoc analysis). |
| `decomp/age_te/model.py` | age_te — within-person comparison of the covariate (age · sex) × scanner TE interaction (plan decomp/age_te/ANALYSIS_PLAN.md): model v_ia = α_i + μ_a + γ · (TE_a − TE_C) · (z_i − z̄) + ε_ia and its inference, the falsification controls, the E3 tools and the rule on whether to run; reads no data. |
| `decomp/age_te/participants.py` | age_te — rules for reading the participant-information tsv and the personal-information guard of the output (plan Sections 3 · 6); reads no data, only parses text it is given. |
| `decomp/age_te/run.py` | age_te run — runs only after the plan decomp/age_te/ANALYSIS_PLAN.md is fixed by a commit (stops if the working copy of the plan differs from the version committed at HEAD); E3 without --participants, E1 only with it. |
| `decomp/design_sim/__init__.py` | design_sim — multi-site design simulation (plan decomp/design_sim/ANALYSIS_PLAN.md). |
| `decomp/design_sim/analyses.py` | design_sim.analyses — analyses M0 – M5 (ANALYSIS_PLAN.md Section 4), separately for each ROI unit, with 95 % intervals of the group coefficient from the t distribution; ComBat uses rules_env of the repository (the neuroCombat 0.2.12 wrapper `_combat` and `combat_apply`, the same formula as neuroCombatFromTraining). |
| `decomp/design_sim/checks.py` | design_sim.checks — the checks that must be able to fail, K-0, K-a – K-e (ANALYSIS_PLAN.md Section 7), run before the main run; each must pass on the clean code and fail on its planted fault, and both must hold for exit status 0. |
| `decomp/design_sim/m6.py` | M6 — a second TE within each site (plan PLAN_M6.md); the first observation is made identical to model.generate by drawing the same random numbers in the same order, and the second observation (TE_s + ΔTE) shares the person · site · ROI terms and differs only in the relaxation term and the measurement noise. |
| `decomp/design_sim/model.py` | design_sim.model — data generation of the multi-site design simulation (ANALYSIS_PLAN.md Sections 2 · 3), at the level of the signal model without a fitter: y_isr = logit ν_r + θ·g_i + b_ir + TE_s·(R2,en,i − R2,in,i) + γ_sr + δ_sr·ε_isr + η_s·g_i. |
| `decomp/design_sim/run.py` | design_sim.run — 1,000 replicates per cell and the metrics (ANALYSIS_PLAN.md Section 5); CPU only, at most 8 concurrent processes; --pilot only measures time on a few cells and replicates (no metric values written). |
| `decomp/design_sim/run_m6.py` | M6 run (plan PLAN_M6.md) — only when the plan is fixed at HEAD; 16 cells × 1,000 replicates, at most 8 concurrent processes, CPU. |
| `decomp/design_sim/summarize.py` | design_sim.summarize — collects the per-cell results of the main run (cells/cell_###.json) into the values of readings R-a – R-d (plan Section 6), the 24 cells of the supplementary table (plan «Change record») and the table of all 288 cells (deposit); no new computation. |
| `decomp/explore/b1/b1_wong_direction.py` | Computes the b1 Wong TE trajectory direction per ROI unit (post hoc, descriptive) from roi_TE*.json, stopping if it disagrees with E1 ψ. |
| `decomp/explore/b1s/b1s_sim.py` | Step 2 of b1s: simulates the Wong processing with AMICO, using the noiseless signal formula, Wong sub-01 TE92 directions and noncentral chi noise. |
| `decomp/explore/b1s/b1s_snr.py` | Step 1 of b1s: computes SNR per TE from the b0 repeats of the Wong sub-01 8 TE series, checking each file's sha256. |
| `decomp/explore/b1s/b1s_summary.py` | Summarizes b1s results (descriptive only, not a verdict): angles and ODI TE slopes, after validation gates on synthetic lines and noise-free angles. |
| `decomp/explore/b2000/b2k_core.py` | Shared mask and metric code for the b2000 noise floor, used by the measurement (b2k_read) and the known-answer test (b2k_sim). |
| `decomp/explore/b2000/b2k_headers.py` | b2k step 1 for the b2000 noise floor: reads headers, bval and sidecars of 109 V4 first sessions into a session table, without voxel values. |
| `decomp/explore/b2000/b2k_read.py` | b2k step 3 for the b2000 noise floor: reads voxels of first-session b0 and b2000 volumes and computes the ventricle mask and within-direction mean/SD (R). |
| `decomp/explore/b2000/b2k_shape.py` | Post hoc descriptive b2000 re-access (b2k): examines the cause of 0 values and the distribution shape, keeping aggregates only and leaving the M′ verdict unchanged. |
| `decomp/explore/b2000/b2k_sim.py` | b2k step 2 for the b2000 noise floor: known-answer test on an analytic head model with eddy motion and noncentral chi noise at each L. |
| `decomp/explore/b2000/b2k_sim_rigid6.py` | Variant of b2k_sim with rigid6 motion only (eddy-current term 0) for the b2000 noise-floor known-answer test, writing to plan_rigid6/. |
| `decomp/explore/b2000/b2k_summary.py` | Summarizes the b2000 noise-floor known-answer test (calibration table, power) and produces the measurement verdict from b2k results. |
| `decomp/explore/b2000/b2k_summary_rigid6.py` | Variant of b2k_summary for the rigid6 plan: summarizes the b2000 known-answer test and measurement verdict, with output folder plan_rigid6/ and its own fixed hashes. |
| `decomp/explore/b2000/cmp_plan_rigid6.py` | Compares the original b2000 result (with eddy-current term) and the plan_rigid6 result (rigid 6 only), descriptive only, without changing the verdict. |
| `decomp/explore/b2000/versions/b2k_shape_4945a692841b.py` | Earlier version of `decomp/explore/b2000/b2k_shape.py`: the version that wrote the recorded output. Post hoc descriptive b2000 re-access (b2k): examines the cause of 0 values and the distribution shape, keeping aggregates only and leaving the M′ verdict unchanged. |
| `decomp/explore/b2000/versions/b2k_summary_432783344be3.py` | Earlier version of `decomp/explore/b2000/b2k_summary.py`: the version that wrote the recorded output. Summarizes the b2000 noise-floor known-answer test (calibration table, power) and produces the measurement verdict from b2k results. |
| `decomp/explore/b2000/versions/b2k_summary_a9e2220a0f1b.py` | Earlier version of `decomp/explore/b2000/b2k_summary.py`: the version that wrote the recorded output. Summarizes the b2000 noise-floor known-answer test (calibration table, power) and produces the measurement verdict from b2k results. |
| `decomp/explore/b2000/zero_check_inline.py` | Verbatim record of the shell code that checked b2000 0 values, σ(DN) and S2 (post hoc, descriptive); this file itself was not run. |
| `decomp/explore/cause/cause1_raw.py` | cause1: exploratory raw-value decomposition of OXF1PRI~OXF2PRI coordinate differences against simulator expectations from frozen test data; not a verdict. |
| `decomp/explore/cause/cause1b_allunits.py` | cause1b: exploratory supplement repeating the raw-value comparison on all OXF1PRI~OXF2PRI ROI units regardless of the judge, plus check 2(c); not a verdict. |
| `decomp/explore/cause/cause2_posthoc_D.py` | Post hoc descriptive cause2 step: Spearman correlation of the measurement D (OXF2PRI − OXF1PRI) with signal-strength and position predictors, within and across subjects. |
| `decomp/explore/cause/cause2_trend.py` | cause2: exploratory descriptive trends of OXF1PRI~OXF2PRI ROI units against log b2000 SNR and position; not a verdict. |
| `decomp/explore/d1/d1_c3.py` | Exploratory d1 on c3: redraws main-cell posterior samples along the verdict path and computes bias and overconfidence per ROI unit, coil combination and k. |
| `decomp/explore/d1/d1_post_sym.py` | Reads stored d1 outputs to describe the four symmetric cells and the H7 selection within ROI units; exploratory, descriptive, not a verdict. |
| `decomp/explore/d2/d2_observed.py` | Computes d2 observed six-session SDs, session-pair difference SDs, centers and subject SNR from the 12 OXF1PRI repeat sessions, per subject, ROI and coordinate. |
| `decomp/explore/d2/d2_sim.py` | Simulates d2: the scatter the judge assumes for the same scanner twice (Δ = 0, L0, OXF1PRI acquisition), calling the frozen simulator unchanged. |
| `decomp/explore/d2/d2_summary.py` | Summarizes d2 descriptively: ratio of observed to simulated difference SDs per ROI unit and coordinate, with a contour over SNR and N_eff; not a verdict. |
| `decomp/explore/d3/d3_kat.py` | Runs the d3 known-answer test without measurement fits: spherical-mean formula check, implementation, failure controls, model adequacy and identification, in envs/reference. |
| `decomp/explore/d3/d3_kat2.py` | Extends the d3 known-answer test (kat2): true r up to 3.0, fitter bounds r 0.4 – 8, profile 0.5 – 6.0, plus a group test. |
| `decomp/explore/d3/d3_real.py` | Fits d3 to the Wong measurement only if the d3_kat2 gate is true: participant and whole-group fits sharing T2_in and r, with sensitivity variants. |
| `decomp/explore/d3/d3_v3.py` | Runs d3 v3: continuous r profiles, interval union over six d∥ × T2_iso variants, and a real-noise known-answer test with coverage and reading gates. |
| `decomp/explore/d3/d3_v3_check.py` | Implementation checks of d3 v3: noise-free group estimates within 0.05 of true r, and bit-identical results across two separate processes. |
| `decomp/explore/d3/d3_v3_check2.py` | Post hoc check2 of d3 v3: asks whether the true value falling far outside the interval at × 2 comes from the inner optimization. |
| `decomp/explore/d3/d3_v3_describe.py` | Counts boundary solutions and fit mismatch sizes from stored d3 v3 outputs, descriptive and defined after seeing the results, writing JSON and Markdown. |
| `decomp/explore/d3/d3_v3_run1_failed.py` | Copy of d3 v3 from run1, which failed: continuous r profiles, interval union over six variants and a real-noise known-answer test. |
| `decomp/explore/d3/d3_v3_summary.py` | Builds d3 v3 tables from stored JSON results into explore/d3/d3_v3.md, without hand-copied numbers. |
| `decomp/explore/d4/d4_2_forward.py` | Builds the d4_2 forward sensitivity table: how apparent values move with r, T2_in, T1 ratio, TE and TR, using repository formulas and no data. |
| `decomp/explore/d4/d4_2_summary.py` | Summarizes d4_2_forward.json as Markdown: tables, analytic properties and the ΔTE direction angle, without data. |
| `decomp/explore/d5/d5_analysis.py` | Exploratory descriptive d5 analysis with modes for scanner discrimination, between-scanner prediction, within-scanner outliers, linkage to D1 and a Markdown summary. |
| `decomp/explore/d5/d5_features.py` | Extracts exploratory d5 features: E1 – E6 from each session's eddy text logs, and 27 ROI units × (x, y, ODI) from roi.json. |
| `decomp/explore/d5/d5_summary.py` | Builds exploratory, descriptive d5 tables from d5_*.json and vr_gate.json into explore/d5/d5.md, without hand-copied numbers; not a verdict. |
| `decomp/explore/d5/selftest/d1_selftest_helpers.py` | Self-test of the d1 helper functions Rec, pick_rule, decomp, kstats and compare on fake draws, kept under d5/selftest. |
| `decomp/explore/d5/selftest/d5_selftest_outliers.py` | Self-test of the d5 outlier detection on simulated data: null cutoffs, flags on null data and detection of three planted outliers. |
| `decomp/explore/d5/versions/d5_analysis_d3c15e5d1281.py` | Earlier version of `decomp/explore/d5/d5_analysis.py`: the version that wrote the recorded output. Exploratory descriptive d5 analysis with modes for scanner discrimination, between-scanner prediction, within-scanner outliers, linkage to D1 and a Markdown summary. |
| `decomp/explore/ledger_derived/derive.py` | Derives manuscript values that no single source-file field gives directly, reading only committed and fixed-copy sources, and writes derived.json. |
| `decomp/explore/rv/rv_4a_coil.py` | Step 4a: extracts only ReceiveCoilName, CoilCombinationMethod and ParallelReductionFactorInPlane from dwi sidecars of the 109 V4 first sessions; exploratory, descriptive. |
| `decomp/explore/rv/rv_ad.py` | Exploratory analyses a and d: location-matched expected misspecification exceedance rate, and a repeat of the D5 between-scanner prediction using SNR-scaled retest; not a verdict. |
| `decomp/explore/rv/rv_b_family.py` | Exploratory analysis b: finds the 25° solution family on a noise-free T2_in × r × T2_iso × true d∥ grid fitted by AMICO F2. |
| `decomp/explore/rv/rv_c_level.py` | Exploratory analysis c: repeats the D2 comparison with simulation centers matched to the observed apparent NDI and FWF levels by inverting config.apparent. |
| `decomp/explore/rv/rv_e2.py` | Exploratory e2 tables: per-cell measurement coverage and thresholds, symmetric-cell thresholds from true-combination test cases, and coil-assumption sensitivity; not a verdict. |
| `decomp/explore/s73/s73_align.py` | S73 §2 (post hoc): rigid-registration check of the multi-TE series with failure controls, the fixed alignment metric recomputed, and the analyses with one registration for all series (paths B and C); it stopped at its first failure control. |
| `decomp/explore/s73/s73_fitter.py` | S73 §3 (post hoc): regresses the fitter-output shift of x, y and raw NDI on the true tissue state in each cell of S1a and compares its tissue-dependent range with the retest width; stored stage-3 fits only. |
| `decomp/explore/s73/s73_noise.py` | S73 §5 (post hoc): session SNR of the multi-TE series by participant and TE, noise-floor simulation of the ROI-matched states (noncentral χ, L = 1, 4, 16) with failure controls, and the summaries per replicate. |
| `decomp/explore/s73/s73_noise_sim.py` | S73 §5 synthesis, noise and AMICO fit in envs/reference: noiseless signals × 1000, noncentral χ noise, F2 fit, voxel means per unit and replicate. |
| `decomp/explore/s73/s73_reading.py` | S73 §5.3 and §5.4 readings from the stored results (no new computation). |
| `decomp/explore/s73/s73_tong.py` | S73 §4 (post hoc): shares at retest level of the Tong negative control by participant and center, the equivalence bound of a common center effect (noncentral t) with its simulated failure control, and the FISO diagnosis of y. |
| `decomp/explore/stored_outputs/coil_contrast.py` | Describes the coil-only contrast between OXF1PRI~X and OXF2PRI~X for six scanners X, from stored outputs only, with no tests. |
| `decomp/explore/stored_outputs/sens_table_from_stored.py` | Builds the per-cell sensitivity table from the stored V4 output v4_c3.json, with no new computation or image reading; post hoc exploratory, not a verdict. |
| `decomp/explore/stored_outputs/snr_pair_prediction.py` | Computes per-cell ranges of the noise hypothesis predicted ratio 30 ÷ SNR_pair, plus ΔSNR and SNR ratio, from stored attribute and SNR tables only. |
| `decomp/explore/stored_outputs/snr_range_counts.py` | Counts whether the 73 of 109 first sessions outside the training range (15 – 60) lie below 15 or above 60, from stored snr_table.json. |
| `decomp/explore/stored_outputs/wong_tr_header.py` | Compares Wong TR with the published Methods by reading only pixdim[4] and time units from the NIfTI headers of the sub-01 8 TE series. |
| `decomp/explore/suff/suff_c3.py` | Exploratory sufficiency test on c3: redraws OXF1PRI~OXF2PRI along the verdict path under three coil conditions (32\|64, 32\|32, 64\|64); CPU only, not a verdict. |
| `decomp/fixed/delta_tau.py` (locked) | Sets the ODI (a) equivalence margin δ_τ from the design alone (measurement fixed analysis document §3), using no ON-Harmony values. |
| `decomp/fixed/fiso_sign.py` (locked) | Computes the predicted FISO sign between scanners by algebra from acquisition records and prior ranges (measurement fixed analysis document §2-1), using no ON-Harmony values. |
| `decomp/manuscript/cards_summary/cards_summary.py` | Summarizes the S4 identifiability cards by counting and copying fields of the two committed card files, with no new inference and standard library only. |
| `decomp/manuscript/checks/check_env_numbers.py` | Recomputes every result number printed in Supplementary Methods S-M9 from the stored results of the measurement-function environment and compares it, in context, with the text, and regenerates Supplementary Table S12 and compares it line by line; it needs the text of the supplementary material, which is not part of this deposit. |
| `decomp/manuscript/checks/make_S17.py` | Prints Supplementary Table S12 (S17 in an earlier version, hence the file name) from the stored results of the measurement-function environment, computing nothing beyond minima, maxima and formatting. |
| `decomp/manuscript/figures/src/c1_mask_location.py` | Computes c1 post hoc descriptive counts of where the stored noise-floor mask lay in template space, and the projections in the left panel of a supplementary figure of earlier versions (panel 7c). |
| `decomp/manuscript/figures/src/c1_split_edge.py` | Exploratory E4 of the rules record of 2026-11-07 (not E4 of the measurement-function environment): splits categories 2 and 3 of the c1 count into the outer-edge part and the part in or next to the ventricular gaps. |
| `decomp/manuscript/figures/src/fig_data_templates_v3.py` | Provides v3 draw functions for the data figures of earlier manuscript versions and one remaining template. |
| `decomp/manuscript/figures/src/fig_real_data.py` | Draws Supplementary Fig. S1 (Wongkornchaovalit) and a supplementary figure of earlier versions (S4 cards, C3) from stored values, computing no new statistic; its Tong figure (Fig5_Tong_negative_control) belongs to earlier manuscript versions. |
| `decomp/manuscript/figures/src/fig_supp_S1.py` | Draws the order of analysis steps and fixing commits (Supplementary Fig. S1 in earlier versions; the list is in docs/fixed_products_and_commits.md), reading each value from its record and stopping on disagreement with the manuscript. |
| `decomp/manuscript/figures/src/fig_supp_S3.py` | Draws the S4 calibration audit (Supplementary Fig. S2 of earlier versions) from stored values only, including stratum coverage of central 90 % intervals against nominal 0.90. |
| `decomp/manuscript/figures/src/fig_theory_v3.py` | Draws the theory figures of earlier manuscript versions (the log-ratio plane, then called the card plane, and the site TE x group T2 interaction) and a graphical abstract, on the shared figure rules; Figures 2 and 3 of the article are drawn by make_F2.py and make_F3.py. |
| `decomp/manuscript/figures/src/figcolors.py` | Defines shared colours and colour-vision checks for the manuscript figures. |
| `decomp/manuscript/figures/src/figstyle.py` | Defines the shared figure style enforcing Elsevier artwork rules, such as widths, with an audit() check, so all manuscript figures use one rule set. |
| `decomp/manuscript/figures/src/judge_matrix_v12.py` | Draws Figure 6 of manuscript v12: the frozen judge on traveling heads as one acquisition matrix, from stored values only. |
| `decomp/manuscript/figures/src/make_F1.py` | make_F1 draws Figure 1: one brain at three echo times, displayed after the verdict. |
| `decomp/manuscript/figures/src/make_F2.py` | make_F2 draws Figure 3 (Section 2.2): compartment weighting in log-ratio coordinates, from the signal model and ledger values only, computing no statistic. |
| `decomp/manuscript/figures/src/make_F3.py` | make_F3 draws Figure 2 (Section 2.4): a site TE x group T2 interaction, from the signal model and ledger values only, computing no statistic. |
| `decomp/manuscript/figures/src/make_F4.py` | make_F4 draws Supplementary Fig. S4 (Figure 4 of earlier versions): the share of units at retest level by ROI unit in the negative control, using guards and ledger checks. |
| `decomp/manuscript/figures/src/make_F5.py` | make_F5 draws Supplementary Fig. S5 (Figure 5 of earlier versions) from stored values only: comparison across scanner models and protocols, fitting, estimating and testing nothing. |
| `decomp/manuscript/figures/src/make_F6.py` | make_F6 draws Supplementary Fig. S6 (Figure 6 of earlier versions) from stored values only: comparison across echo times in the multi-TE series of three participants. |
| `decomp/manuscript/figures/src/make_F7.py` | make_F7 draws Supplementary Fig. S2 (Figure 7 of earlier manuscript versions): the frozen judge on traveling heads, the acquisition matrix alone. |
| `decomp/manuscript/figures/src/make_F8_design_sim.py` | Figure 4 of manuscript v16z9 → v17d (Results 4.2): the multi-site design simulation (post hoc; assumed values; Supplementary Methods S-M11). Drawn from the stored cell metrics only; nothing is simulated, fitted or tested here. |
| `decomp/manuscript/figures/src/make_FS5.py` | make_FS5 draws a supplementary figure of earlier versions: where the noise-floor mask lay and the noise floor it measured. |
| `decomp/manuscript/figures/src/make_FS6.py` | make_FS6 draws Supplementary Fig. S3 (Supplementary Methods S-M9) from the stored results of the measurement-function environment; nothing is fitted or tested except the expected-share curves of panel (a), which use the analysis function itself. |
| `decomp/manuscript/recompute/a_values.py` | Recomputes the values marked ⟦A⟧ in an earlier manuscript version (its Sections 2.2, 2.4 and 2.5 and Figs. 1 and 2; numbering of that version) with repository code, on CPU, writing nothing to the repository. |
| `decomp/manuscript/recompute/tr_fiso.py` | Recomputes the two «TR on FISO» values of Section 2.4 of an earlier manuscript version (numbering of that version) with repository code, on CPU, writing nothing to the repository. |
| `decomp/remeasure/analysis.py` (locked) | Turns ROI values into scales, pair summaries, retest SDs and session SNR for the measurement. |
| `decomp/remeasure/cmds.py` (locked) | Assembles FSL commands as pure functions without running them, using the old registered pipeline's options except eddy_cpu and phase direction read from sidecars. |
| `decomp/remeasure/fit_f2.py` (locked) | fit_f2 fits F2 (AMICO 2.1.1 NODDI, same settings as S4) and writes extended coordinate maps (MD, FA, s_low, s_2000) in envs/reference. |
| `decomp/remeasure/guard.py` (locked) | Loader guard and access record: before images open, require() checks git state of the fixed analysis document and pipeline files, stopping with GuardError on mismatch. |
| `decomp/remeasure/judge_onh.py` (locked) | Verdict code for the ON-Harmony pair tests (FISO sign, FISO τ, ODI (a), ODI (b)) and s_ref_F2, run on roi.json files behind the loader guard. |
| `decomp/remeasure/judge_tong.py` (locked) | Verdict code for Tong 2020 (blinded, negative control): B-iii retest SD, T-B magnitude, T-A direction by rotation permutation and T-C 2 × 2. |
| `decomp/remeasure/judge_v4.py` (locked) | judge_v4 is the verdict code for V4, the traveling-heads measurement calibration test: applies the frozen S4 judge (weights 85eb8b0) to measurement observations, with no training. |
| `decomp/remeasure/judge_wong.py` (locked) | Verdict code for the Wongkornchaovalit 2025 multi-TE demonstration: linearity, TE-group intercept stability and ODI TE stability, on roi_<series>.json files behind the guard. |
| `decomp/remeasure/rconfig.py` (locked) | Holds values, paths and session lists of the processing pipeline, which runs on CPU only with FSL 6.0.7.23 and F2 = AMICO 2.1.1. |
| `decomp/remeasure/roi.py` (locked) | Turns standard-space maps into 48 skeleton ROI atlas labels and 27 ROI units, with voxel means over label ∩ skeleton ∩ valid mask including 0. |
| `decomp/remeasure/run.py` (locked) | Drives session processing for ON-Harmony, Tong and Wongkornchaovalit on CPU, passing guard.require() and appending to the access record whenever real-data paths are touched. |
| `decomp/remeasure/summarize_judgements.py` | Extracts JSON from the four verdict result files in the value area (V4 C3, ON-Harmony pair test, Tong, Wongkornchaovalit) and draws Markdown from it. |
| `decomp/rules_env/__init__.py` | Package marker of rules_env, the measurement-function environment (algebra, cards, synthesis and fitting, metrics, real data, analyses) described in decomp/rules_env/ANALYSIS_PLAN.md. |
| `decomp/rules_env/_paths.py` | Finds the repository root: two levels above this file in the repository layout (decomp/rules_env), or the environment variable RULES_ENV_REPO in the development layout. |
| `decomp/rules_env/analyses.py` | Runs the analyses S1, S2, K8, E1 – E4 of the environment's analysis plan (ANALYSIS_PLAN.md) and E5, E6, S3 of its §7 (R2) in the main environment, calling the AMICO fit (simrun in envs/reference) as a subprocess. |
| `decomp/rules_env/cards.py` | Defines the acquisition-condition cards (name to TE, TR, bvals, bvecs and scheme) for the ON-Harmony card units, the E5 acquisitions, the Wong multi-TE series and the Tong data. |
| `decomp/rules_env/check_synthfit.py` | Checks rules_env.synthfit in envs/reference (AMICO), outside pytest: signal equality with CardProto and Proto and with synth.noiseless, a refit of 2,000 stored E5 rows to within ±1e-5, and a speed measurement with the Task 7 extrapolation. |
| `decomp/rules_env/empirical.py` | Reads the already computed ROI values (roi.json) with the same rules as the existing verdict code and returns per-scale ROI unit values, reusing the loaders and judge functions of decomp/remeasure; it adds only the scale transforms and the per-scale retest SD, and opens no images. |
| `decomp/rules_env/metrics.py` | Computes ORS, the offset-subtracted application of the retest rule of manuscript 3.4 (ROI unit difference \|Δ\| ≤ c·√2·s and share ≥ 0.90), and the quantities derived from it (TE tolerance, τ̂, scale models, proportional bias, rank flips, S2 ratio), without any new threshold. |
| `decomp/rules_env/model.py` | Implements the algebra of manuscript formulas (2), (3), (6), (7) and Appendix A proposition 4 (weights, apparent fractions, coordinates, translation and inversion), calling the weights, apparent fraction and coordinates of decomp/stage1/config.py as they are. |
| `decomp/rules_env/regress.py` | Stop gate that recomputes the existing reported values with the new loader (manuscript 4.1, 4.2, 4.3): exits with code 1 if any value differs and computes no new quantity. |
| `decomp/rules_env/simrun.py` | Runs the simulation jobs in envs/reference (AMICO) only: reads the job bundle (json) written by analyses in the main environment as a subprocess call and writes an npz with ndi, fwf, odi, nu_app and f_app per job. |
| `decomp/rules_env/synthfit.py` | Synthesizes noise-free signals from explicit tissue and relaxation values and fits them with decomp/stage1/amico_fit.py as a subprocess (d∥ 1.7, d_iso 3.0 µm²/ms, b0 threshold 10), to be called only from envs/reference (AMICO 2.1.1). |
| `decomp/s4/fisher.py` | Computes s4 Appendix B Fisher sensitivity: the MDE^lin ratio before and after adding nuisance factors ψ to the linear baseline; not used for verdicts. |
| `decomp/s4/plan.py` (locked) | Builds the deterministic s4 case lists per dataset (θ*, conditions, nuisance variables) for training, test, oracle, reference and cards, using numpy only. |
| `decomp/s4/s4cfg.py` (locked) | s4cfg holds the s4 specification values from decomp/SPEC_S4.md §3 and samples priors, conditions and nuisance variables, using numpy only. |
| `decomp/s4/simulate.py` | Simulates one s4 case, one person's ROI at two sites: interpolated compartment signals, coil noise, AMICO F2 fit, extended coordinates and ROI mean. |
| `decomp/s4/tables.py` | Builds s4 per-design compartment signal tables and interpolates them linearly per voxel by ODI and d⊥/d∥, in envs/reference. |
| `decomp/s4/train.py` | Trains the s4 models on CPU: C3 and C7 ensembles, oracle NPE and observation density q(x \| a), with sbi 0.27.0 and torch. |
| `decomp/s4/validate.py` (locked) | Runs s4 validation, cards and maps on CPU, writing step results to decomp/results/s4/, using shrinkage and posterior SD rather than ML point estimates. |
| `decomp/stage1/amico_fit.py` | Fits F2 for stage1 with AMICO 2.1.1 NODDI on a chunk directory, writing NDI, FWF and ODI to f2_<dpar>.npz, in envs/reference. |
| `decomp/stage1/analyze.py` (locked) | Runs the stage1 decomp analysis: directions, SG1–SG6, Π maps, E4 and report inputs, per SPEC_v1, in the main environment. |
| `decomp/stage1/analyze2.py` (locked) | analyze2 in stage1 runs the decomp stage 2 analyses X1–X6 of SPEC_v2 on CPU with F2 only, writing to decomp/results/stage2/. |
| `decomp/stage1/analyze3.py` (locked) | analyze3 in stage1 runs the decomp stage 3 analyses X7–X13 of SPEC_v3 on CPU with F2 only, writing to decomp/results/stage3/. |
| `decomp/stage1/config.py` (locked) | Defines stage1 grids, protocols, generator axes and experiment rows from decomp/SPEC_v1.md §5–§7, using numpy only. |
| `decomp/stage1/f1_compare.py` | f1_compare sets F1 and F2 verdicts side by side for decomp stages 1, 2 and 3 with unchanged definitions, choosing neither, in stage1. |
| `decomp/stage1/f1_yield.py` | f1_yield runs F1 (cuDIMOT) on GPU 1 in stage1, yielding at once when another process appears and resuming later, logging per chunk. |
| `decomp/stage1/lib.py` (locked) | Pure numpy and scipy analysis functions of decomp stage1, per SPEC_v1 §8–§10. |
| `decomp/stage1/report.py` | Writes the decomp stage1 report (summary.json, REPORT.md) from directions.json, e4.json and the direction store, per SPEC_v1 §17. |
| `decomp/stage1/results_manifest.py` | Lists sha256 of committed results under decomp/results/** so stage1 analyses cannot overwrite them and tests can compare them before and after runs. |
| `decomp/stage1/selftest_synth.py` | Known-answer gate of the stage1 synthesizer, without the fitter: reproduces Δx and Δy within 1e-10 and checks T-app and T-in signals. |
| `decomp/stage1/summary3.py` | summary3 in stage1 computes the decomp stage 3 summary.json tables from x*.json, cards.json and run records, without hand copying. |
| `decomp/stage1/synth.py` (locked) | Synthesizes one chunk of decomp stage1 signals and extended coordinates, writing data.nii and the brain mask for the fitter, in envs/reference. |
| `decomp/v4_recompute/compare.py` | Compare a recomputed V4 output with the deposited one. |
| `decomp/v4_recompute/read_log.py` | Write-ahead read log: the line that names a file is on disk before the file is read. |
| `decomp/v4_recompute/run_subset.py` | Recompute part of the V4 calibration test with the unchanged judge code. |
| `dmri_cov/__init__.py` | Package marker of dmri_cov, a coverage benchmark for dMRI microstructure uncertainty (see dMRI_PLAN.md). |
| `dmri_cov/coverage/__init__.py` | Empty package marker of dmri_cov/coverage. |
| `dmri_cov/coverage/arms.py` | Defines the comparison arms and mechanically enforces what each arm may access from the target. |
| `dmri_cov/coverage/bootstrap.py` | Cluster bootstrap that resamples one side only, because calibration uncertainty is already in the interval through k and must not count twice. |
| `dmri_cov/coverage/metrics.py` | Compares arms by always reporting coverage, sharpness and vacuity together, because coverage alone misses intervals that are calibrated but uninformative. |
| `dmri_cov/coverage/sample_size.py` | Computes the traveling-subject calibration sample size and harm threshold: calibration shifts the center but adds variance, so small samples can do worse. |
| `examples/reading_rules_example.py` | Worked example applying the reading rules of Table 2 to a table of ROI values; it does not separate relaxation from tissue. |
| `ledger/derived/derive.py` | Derived ledger values (2026-11-04): values used in the manuscript that do not come directly from one field of one source file; reads only committed versions of the repository and the fixed copies of the ledger, no images. |
| `ledger/derived/r27_values.py` | R27 (v16z3): the two ratios of the retest-SD sentence of Section 4.1 (median retest SD of c10 ÷ median retest SD of OXF1PRI, for x and y); two divisions, no new data. |
| `ledger/derived/residual_s15.py` | Supplementary table of the size of the residual of Proposition 4 (derived, no data): the logit residual of Eq. (7) and the raw-scale change in the signal model of decomp/stage1/config.py, with the row of Fig. 2 as failure control. |
| `ledger/derived/rules_env_derived.py` | Derived ledger values from the outputs of the measurement-function environment (rules_env) that are not printed directly (no new computation: pooling, counting, one formula); logs each result json before reading it. |
| `ledger/derived/rules_env_r2_derived.py` | Derived ledger values from the second set of rules_env outputs (s3.json, e6.json): one maximum and one interval, no new computation; logs each result json before reading it. |
| `ledger/derived/s15_within_site.py` | Within-site column of the supplementary residual table: patient − control difference of apparent NDI at the same TE in five scenarios (Eqs. (3) · (4); derived, no data), by the same code path and inputs as residual_s15.py. |
| `ledger/derived/s2_combat_ls.py` | Derived ledger values: summary of the ComBat location-and-scale arm (combat_ls) of analysis S2 of the measurement-function environment (no new computation: minimum and maximum of a ratio and the number of cells outside 0.80–1.25). |
| `ledger/derived/s73_derived.py` | Derived ledger values from the S73 results and rules_env s1a that are not printed directly (no new computation: maxima, minima, values before rounding); logs each result json before reading it. |
| `ledger/derived/s74_clinical.py` | Derived ledger value: difference of apparent NDI within one site when only T2 lengthens (Eq. (3), TE weighting only, no new data), with the inputs of Section 2.4 and Brandhofe et al. 2021. |
| `ledger/derived/s75_values.py` | Five small values (exploratory; existing code, CPU, no new data): the FISO size of Proposition 2, tolerance bounds, noise-only expectations of the multi-TE rules, power of the group × TE term, and existing results moved to the main text. |
| `ledger/derived/s76_values.py` | The T2 path of the «Free water» paragraph of Section 5.2 (exploratory; existing code, CPU, no new data): how the diffuse tissue T2 lengthening of Section 2.4 moves apparent FISO. |
| `ledger/derived/s82_values.py` | Values printed in the manuscript from the results of the multi-site design simulation (decomp/design_sim); no new computation, per-cell metrics copied and rounded. |
| `ledger/derived/s83_values.py` | Further values from the design-simulation results (no new computation: per-cell metrics copied and counted). |
| `ledger/derived/s84_values.py` | Further values (no new simulation: stored cell metrics copied and counted; a bulk T2 check from the signal-model formula only), each with a check that must be able to fail. |
| `ledger/derived/s86_values.py` | Further values (no new simulation: stored cell metrics, signal-model formulas and arithmetic on ledger values), each computed a second way as an explicit check. |
| `ledger/derived/s87_values.py` | Known-answer test of the continuous-covariate corollary of Proposition 4, and the maximum over the unshown cells (all cells with a relaxation difference); --plant must fail. |
| `ledger/derived/s89_values.py` | Ledger values of the age × TE analysis: the power estimate before computing for E1 (not run; decomp/age_te/E1_GO_DECISION.md) and the between-person TE slopes of the multi-TE data (E3); --plant must fail. |
| `ledger/derived/s90_values.py` | New values of v17c: bias of standard ComBat in NDI units (A1), the isotropic-T2 values (A6) and the linearity margin (A10); --plant must fail. |
| `ledger/derived/s91_values.py` | Ledger values of the second TE within each site (M6, plan decomp/design_sim/PLAN_M6.md); known-answer test: M2 · M3 of the M6 cells equal the original cells within 1e−12, and a copy paired with shifted cells must fail. |
| `ledger/derived/s92_values.py` | Table of the relaxation sensitivity value applied to published clinical results (plan decomp/rsv_table/PLAN.md): d, d/TE and the T2 lengthenings, with the bulk T2 equivalent; known-answer tests (Section 2.4 example about 5 %; bulk values of s84 within 1e−12) and two planted faults. |
| `ledger/derived/s93_values.py` | Dependence of the relaxation sensitivity value on the assumed control T2 (prior of Supplementary Methods S-M2), the proton-density path and the isotropic-T2 path of apparent FISO (Section 5.2; exploratory; signal model only), with known-answer checks against Table 3 and two planted faults that must fail. |
| `prereg/scripts/b10_synth.py` | b10 measures how far fitted logit NDI shifts when TE changes 92 → 68 ms, running three estimator arms on the same synthetic signal. |
| `prereg/scripts/lockd/_harness.py` | LOCK_D verdict harness: runs dummy data end to end before unsealing, freezing the contracts on input shape, missing cells and related rules. |
| `prereg/scripts/lockd/compute_s_ref.py` | Computes and freezes s_ref, the reference retest SD s of the stratum 1 interval, from OXF1PRI repeats of 2 people × 6 sessions. |
| `records/x12_tong_driver.py` (run record) | Run record (x12) of the Tong judgment: calls, from outside, the computation intended by the command line of `decomp/remeasure/judge_tong.py`, unpacking the loader's three return values; recorded as a deviation (Supplementary Table S1). Kept as run, with Korean comments; `reproduce/tong_judge.py` does the same. |
| `reproduce/environment_recompute.py` | Recompute the measurement-function environment (Supplementary Methods S-M10) and compare with the deposit. |
| `reproduce/get_atlas.py` | Check a copy of the atlas that step 3 needs and put it in place. Steps 1 and 2 do not need it. |
| `reproduce/judge_v4_run.py` | Runs decomp/remeasure/judge_v4.py with the acquisition schemes in data/protocols, because the locked decomp/stage1/config.py fixes the scheme file location. |
| `reproduce/ledger_fields.py` | Provides the ledger field readers (json:, jsonl:, expr:, md:, lit:, none:) used by step1_check_numbers.py. |
| `reproduce/make_local_git.py` | Make the local Git history that the loader guard reads, for a copy unpacked from the Zenodo archive. |
| `reproduce/step1_check_numbers.py` | step1 of the reproduction: checks each main-text number in ledger/numbers.csv against its source file by recomputing the field. |
| `reproduce/step2_recompute.py` | step2 of the reproduction: recomputes the fixed tests from the deposited ROI tables and, for the calibration test, the frozen weights, reading no image. |
| `reproduce/tong_judge.py` | Runs the Tong verdict through its functions, because the judge_tong.py command-line entry stops with an error (deviation (iii) in Supplementary Table S1). |

## Tests

Run with `uv run --frozen python -m pytest tests` (see `README.md`).

| File | What it tests |
|---|---|
| `tests/conftest.py` | Test settings: no stale bytecode, and CPU only. |
| `tests/test_age_te.py` | age_te — quick checks of the covariate × scanner TE interaction (plan decomp/age_te/ANALYSIS_PLAN.md) on synthetic data with known answers, each with a failure control; reads no real data. |
| `tests/test_decomp_stage1.py` | Known-answer tests for decomp stage1 per SPEC_v1 §3, §7 and §8–§10, using no fitter and no real data. |
| `tests/test_design_sim.py` | design_sim — quick checks of the multi-site design simulation (plan decomp/design_sim/ANALYSIS_PLAN.md), each with a failure control. |
| `tests/test_design_sim_m6.py` | Quick checks of design_sim M6 (a second TE within each site · plan decomp/design_sim/PLAN_M6.md), each with a failure control. |
| `tests/test_fixed_design.py` | Known-answer tests of the design calculations in the measurement fixed analysis document, using acquisition records and prior ranges and no ON-Harmony values. |
| `tests/test_judge_onh.py` | Known-answer tests on fake data for the ON-Harmony pair test loader and the end markers of the Tong and Wongkornchaovalit loaders, with no real images. |
| `tests/test_public_repository.py` | Checks of this repository as published: the guard passes, the locked files match, and no server path is present. |
| `tests/test_read_log_first.py` | The read log is written before the file is read, and a log that cannot be written stops the read. |
| `tests/test_reading_rules_example.py` | Tests examples/reading_rules_example.py, pairing each check with a failure control so that every check is able to fail. |
| `tests/test_remeasure.py` | Known-answer tests of the measurement pipeline: loader guard, commands matching the old pipeline, 27 ROI units, session lists and access record, without images or FSL. |
| `tests/test_rules_env_analyses.py` | Tests the pure functions of rules_env.analyses: fiber-direction mean coordinates, S2 group differences and methods, TE slopes, batch offsets, ComBat reapplication and L1 relaxation draws. |
| `tests/test_rules_env_cards.py` | Tests that the acquisition-condition cards of rules_env.cards equal the values of the manuscript and the algebra cards (ON-Harmony pairs and TE range, E5, Wong TE grid, Tong card, design-letter map). |
| `tests/test_rules_env_empirical.py` | Tests, on fake data, that the scale transforms, retest SDs and scanner values of rules_env.empirical equal those of the existing verdict code (judge_tong, judge_onh, s_ref_f2). |
| `tests/test_rules_env_metrics.py` | Tests rules_env.metrics: the offset-subtracted retest share (ORS), TE tolerance, scale comparison, proportional bias, rank flips and the S2 reading. |
| `tests/test_rules_env_model.py` | Tests that the algebra of rules_env.model (formulas (3), (6), (7)) reproduces the values of the manuscript and the ledger, including Fig. 2, Proposition 4, the rows of Supplementary Table S11 and the inversion round trip. |
| `tests/test_s4.py` | Known-answer tests for s4 per SPEC_S4 §9 and §13, run in the main environment with numpy and torch only. |
| `tests/test_s4_fisher.py` | Known-answer tests of s4 Appendix B Fisher sensitivity: ratio ≥ 1, ratio → 1 as w → 0, ratio = 1 when J_ψ = 0. |
