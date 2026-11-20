# Details moved from the Supplementary material

This document holds text that was in the Supplementary material of the manuscript until 2026-11-16 and was moved here to shorten it: clock times, file names and provenance records, the full rows of some tables (with a summary kept in the supplement), and details of the judge study. The text is as it was printed, without change. Each section is named after the place it came from; where the Supplementary material points here, it uses that name. Sections named in the Supplementary material carry the names it uses; other sections keep the numbering of the Supplementary material at the time the text was moved (2026-11-16), so their table, figure and subsection numbers can differ from the current ones.

## From S-M2 Relaxation levels (L3)

(T2,en × 2.1–3.4, T2,in × 0.77–1.02 and the T1 of both tissue compartments × 1.17–1.56, each drawn uniformly; the T2 ranges are the lesion ball T2 (110–173 ms) divided by the T2,en center of 51.25 ms and the lesion stick T2 (63–84 ms) divided by the T2,in center of 82 ms (Lampinen et al., 2019), and contain the corresponding factors × 2.9 and × 0.87 for the zeppelin and stick values (149 and 71 ms) of Lampinen et al. (2020); the T1 factors were derived from lesion T1 (John et al., 2024) and lesion R1 (Hagiwara et al., 2017) relative to normal-appearing white matter)

## From S-M4 Separation attempts (opening summary of the three procedures)

(stage 1 projected each pair's whitened difference off the subspace spanned by nuisance directions — TE, TR, noise floor, design, and model misspecification combined with the protocol difference — and summarized the remainder as an index I, interpreted only where Π, the preservation of the tissue directions under the projection, was at least 0.5; stage 2 asked, by exhaustive search over a dense grid, whether one tissue state measured at one acquisition and a different tissue state with different relaxation measured at another could give features that agree within noise; stage 3 formed partial-identification intervals for the paired NDI and FISO differences from one noisy observation per site, admitting every relaxation state of a stated set, the 17 points of the stage-3 level L1 with CV 3 %, ρ 0 and T1 order iii)

## From S-M4 Separation attempts (opening summary of the ComBat example)

A ComBat example (X12) applied ComBat (Fortin et al., 2018) to simulated two-site data in which a tissue difference between patients and controls was either balanced across sites or confounded with site; patients and controls shared one distribution of compartmental relaxation. Its design and result follow below; it is not a ranking of methods.

## From S-M5 (judge study)

The judge study was described in subsections S-M5.1–S-M5.8 of the Supplementary material as it was on 2026-11-16; the sections below keep those subsection numbers in their names.

## From S-M5 (clock times)

No threshold, grid, rule or model was changed after the results were first opened (2026-10-26 09:04:55 KST); points the specification left open were fixed afterwards in the judge code, before the real-data rules were released (c8825e2, 2026-10-27 19:11:18 KST; S-M7.2).

## From S-M5.1 Estimation problem

The judge is summarized here; the specification in full, with all hyperparameters and data-set sizes, is in the deposit (`docs/judge_specification.md`; Korean original `decomp/SPEC_S4.md`).

For one region of interest (ROI) in one person measured at two sites A and B, the judge estimates θ* = (Δ, tissue_A): the difference Δ = (ΔNDI, ΔFISO, ΔODI) in true (relaxation-free) fractions between B and A, together with the site-A tissue (NDI_A, FISO_A, ODI_A), because the identifiability of Δ depends on the reference tissue state. The observation x is the ROI means of the F2 outputs (NDI, FISO, ODI) at both sites in the coordinates of Section 2.2 (C3, six dimensions) or, in the extended input, these plus mean diffusivity, fractional anisotropy and per-shell mean signal (C7, fourteen dimensions). The estimator is conditioned on the acquisition (TE and TR at each site, direction scheme, b0 settings and number, coil, SNR at each site and the effective ROI size N_eff), on the relaxation scenario (S-M2) and on the generator assumption k (three tortuosity conventions × three T1 orderings = nine combinations); carrying these as conditions, not priors, lets one trained estimator report every level and assumption without retraining. A shift common to all indicators is absorbed into a difference in latent means, which is why the judge states relaxation levels explicitly. Each person's compartmental relaxation times, the tissue intrinsic parallel diffusivity (five values from 1.7 to 3.0 μm²/ms, while the fitter uses 1.7 μm²/ms), within-ROI heterogeneity, fiber orientation and noise are marginalized in training. The domain is two-shell acquisitions (b ≈ 1000 and 2000 s/mm²) with the ON-Harmony direction schemes, TE 60–110 ms and TR 3.0–11.0 s; external datasets with other designs receive algebraic cards only.

## From S-M5.2 Simulator and ROI-level observation model

Tissue_A was drawn uniformly over the generator's tissue range (NDI 0.30–0.80, FISO 0.01–0.30, ODI 0.08–0.40), and a draw was redrawn whenever the site-B tissue fell outside it. Δ was drawn uniformly within ±0.10 (NDI, ODI) and ±0.05 (FISO) for reporting and over twice this width for training, so that narrower and wider priors can be obtained by importance reweighting. Relaxation levels were drawn as L0 25 %, L1 35 %, L2 20 % and L3 20 %; acquisition pairs shared TE and TR (30 %), differed in TE only (15 %) or in TR only (15 %), or were drawn independently at each site with the direction scheme (40 %); SNR was log-uniform in [15, 60] at each site. For each site, K = 256 voxels (the ROI center plus within-ROI heterogeneity) were synthesized, given noise and fitted with F2, and the ROI observation was

$$x = m + \varepsilon, \qquad \varepsilon \sim \mathcal{N}\!\left(0,\ \hat\Sigma\,(1/N_{\mathrm{eff}} - 1/K)\right),$$

where m and Σ̂ are the mean and covariance of the K fitted values and N_eff is log-uniform in [16, 200]. This follows the order of a real pipeline (fit each voxel, then average over the ROI); N_eff was not estimated from ON-Harmony data. The data sets (among them 100,000 training cases and 20,000 test cases in 20 strata) were generated once and frozen by hash.

## From S-M5.2 and S-M5.3 (longer summaries)

#### S-M5.2 Simulator and ROI-level observation model {.unnumbered}

Tissue_A was drawn uniformly over the generator's tissue range (NDI 0.30–0.80, FISO 0.01–0.30, ODI 0.08–0.40), and Δ was drawn uniformly within ±0.10 (NDI, ODI) and ±0.05 (FISO) for reporting and over twice this width for training. Relaxation levels were drawn as L0 25 %, L1 35 %, L2 20 % and L3 20 %; SNR was log-uniform in [15, 60] at each site; for each site, K = 256 voxels were synthesized, given noise and fitted with F2, and the ROI observation was formed from their mean with N_eff, the effective number of independent voxels per ROI mean, log-uniform in [16, 200] (not estimated from ON-Harmony data). The data sets (among them 100,000 training cases and 20,000 test cases in 20 strata) were generated once and frozen by hash.

#### S-M5.3 Estimator: architecture and training {.unnumbered}

The posterior was estimated by neural posterior estimation with sbi 0.27.0 (Boelts et al., 2025): a neural spline flow (Durkan et al., 2019) with five layers and 128 hidden units, and five independently seeded models per arm mixed with equal weights (Hermans et al., 2022). An oracle estimator served the reference comparison of criterion (ii), and a conditional density estimator q(x | a, k) the misspecification score of criterion (v). All twelve models were trained on CPU only in 1.61 h; all stopped by the plateau rule. Training losses were not inspected. Weights, normalization constants and mixing rules were frozen by hash before validation.

## From S-M5.3 Estimator: architecture and training

The posterior p(θ* | x, a) was estimated by neural posterior estimation with sbi 0.27.0 (Boelts et al., 2025): a neural spline flow (Durkan et al., 2019) with five layers and 128 hidden units, conditioned on (x, a) through an embedding network, and five independently seeded models per arm mixed with equal weights (Hermans et al., 2022). An oracle estimator served the reference comparison of criterion (ii), and a conditional density estimator q(x | a, k) the misspecification score of criterion (v). All twelve models were trained on CPU only in 1.61 h; all stopped by the plateau rule, after 61–123 epochs. Training losses were not inspected. Weights, normalization constants and mixing rules were frozen by hash before validation.

## From S-M5.4 Identifiability metrics (definitions, importance reweighting, design-level summaries, MDE, directions)

**Shrinkage of one observation** along a direction u:

$$c_u(x, a) = 1 - \frac{\operatorname{Var}_{\text{post}}(u^\top \Delta \mid x, a)}{\operatorname{Var}_{\text{prior}}(u^\top \Delta)},$$

computed under the reporting prior by importance reweighting. Because both priors are uniform, the weights are 0 or 1, so the effective sample size is the number of retained samples: samples outside the reporting prior are discarded and sampling continues until 2,000 remain for each generator assumption k (18,000 per card observation), drawing up to 20 times that number; an observation that falls short is marked not reweightable.

**Design-level summaries** were computed from M = 100 simulated observations per condition (acquisition records only, no image values):
- C_u(a): the median of c_u;
- σ_u(a): the median posterior standard deviation;
- the minimum detectable difference MDE_u = 2.80 σ_u (two-sided 5 % level, 80 % power, normal approximation).

**MDE as the reported quantity.** MDE is reported together with C. It refers to **one paired comparison of ROI means**, not to a comparison of group means. In group comparisons the noise component shrinks with sample size, whereas a systematic relaxation bias does not.

**Directions u:**
- the unit vectors of NDI, FISO and ODI;
- u_relax, the direction in (ΔNDI, ΔFISO) of the relaxation shift produced by the pair's TE and TR difference at the case's reference tissue (from the Jacobian of the log-ratio coordinates);
- its orthogonal complement.

## From S-M5.5 Validation criteria (i)

### (i) Calibration

**(i) Calibration (decisive).**
- Central credible intervals at 50, 80, 90 and 95 % were evaluated from 1,000 posterior samples per case, under the conditional posterior at the true k.
- A stratum × component × level passed if empirical coverage was at least nominal − 3√(p(1 − p)/n). Coverage above nominal + 0.05 was flagged as conservative.
- Failed strata are reported as "not decidable". **The judge failed overall if ten or more of the twenty strata failed.**
- Simulation-based calibration (SBC) rank histograms and TARP joint coverage (Lemos et al., 2023) were diagnostics, not criteria.

## From S-M5.5 Validation criteria (ii)

### (ii) Prior leakage

**(ii) Prior leakage — reference comparison.**
- The oracle estimator was compared with reference posteriors computed by Markov chain Monte Carlo on the exact likelihood.
- Sampler: slice sampling (Neal, 2003), 4 chains, 500 burn-in and 2,500 retained draws, R̂ ≤ 1.01, effective sample size ≥ 500, run in a reparameterized state of noise-free apparent values with the Jacobian included. The known-answer tests listed in commit b5ca002 were passed before use.
- Pass required all three:
 - (a) median |C_npe − C_ref| ≤ 0.10 for u_relax in the TR-only and TE+TR strata;
 - (b) C_ODI ≥ 0.9;
 - (c) median classifier two-sample test accuracy ≤ 0.65 (Lopez-Paz and Oquab, 2017).
- **The judge failed overall if (ii) failed.** Failure would mean the estimator reports information that the algebra does not contain.
- Importance sampling, considered first, was replaced before training because its effective sample size was 1.0–3.9 per 10⁶ draws.

## From S-M5.5 Validation criteria (v)

### (v) Injection tests

- **Injection tests, descriptive only:** four out-of-range conditions with 1,000 cases each — SNR 10; tissue intrinsic diffusivity 3.4 μm²/ms; b-value scaling ±5 % (a factor of 0.95 or 1.05 per case, applied in the synthesis only; fitting used the recorded b-values); N_eff 400. Detection rates and the coverage of undetected cases were reported.

## From S-M5.5 and S-M5.8 (validation text before the table)

#### S-M5.5 Validation criteria and predictions fixed before validation {.unnumbered}


Validation used twenty strata: five acquisition types (matched; TE only; TR only; TE and TR; TE, TR and direction scheme) × four relaxation levels. Two acquisitions have the same TR when the sums of their card TR ranges agree to 10⁻⁹, so NOT2ING-89 and NOT2ING-98, both recorded as 4.00–4.38 s, form a TE-only pair; the 28 card pairs fall 1, 3, 0, 9 and 15 into the five types. Each stratum had 1,000 new test cases.

**(i) Calibration (decisive).** Central credible intervals at 50, 80, 90 and 95 % were evaluated from 1,000 posterior samples per case at the true k. A stratum × component × level passed if empirical coverage was at least nominal − 3√(p(1 − p)/n); failed strata are reported as "not decidable". **The judge failed overall if ten or more of the twenty strata failed.** Calibration diagnostics (rank histograms, joint coverage) were not criteria.

**(ii) Prior leakage — reference comparison.** The oracle estimator was compared with reference posteriors computed by Markov chain Monte Carlo on the exact likelihood (sampler settings in the deposit). Pass required all three: (a) median |C_npe − C_ref| ≤ 0.10 for u_relax in the TR-only and TE+TR strata; (b) C_ODI ≥ 0.9; (c) median classifier two-sample test accuracy ≤ 0.65 (Lopez-Paz and Oquab, 2017). **The judge failed overall if (ii) failed.** Failure would mean the estimator reports information that the algebra does not contain. Importance sampling, considered first, was replaced before training because its effective sample size was 1.0–3.9 per 10⁶ draws.

**(iii) Prior sensitivity.** Reported as labels, not as a pass/fail criterion.
- Card posteriors were compared at half and twice the reporting-prior width; half the width is a reweighting of the reporting-prior posterior, marked "not reweightable" when fewer than 200 samples remain. A class that changed is labeled "prior-dependent", and a prior-dependent cell is not counted as identified; the identified counts of S-M9.1 rest on this rule. How the 200-sample threshold was applied at half the width, and what it means for the NDI count, is in S-M5.8.

**(iv) Faithfulness.** The ensemble, not individual members, had to pass (i). Member coverages are reported alongside.

**(v) Misspecification score.**
- Score: s(x, a) = −log[(1/9) Σ_k q(x | a, k)].
- Threshold: the 99th percentile in each in-distribution stratum. The per-stratum false-alarm rate had to be ≤ 2 %.
- In real data, an observation above threshold is "not decidable" (its handling in V4 is given in S-M7.2).
- **Injection tests, descriptive only:** four out-of-range conditions with 1,000 cases each (SNR 10; tissue intrinsic diffusivity 3.4 μm²/ms; b-value scaling of ±5 %; N_eff 400). Detection rates and the coverage of undetected cases were reported.

**Validation predictions** (P1 and P2 decisive, P3–P6 reported):

| ID | Prediction |
|----------------|------------------------------------------------------------------------------------|
| P1 | In matched pairs at L0, all three components are identified |
| P2 | At L0–L2, ODI is identified in every acquisition type |
| P3 | In TR-only pairs at L0, C along u_relax is at least 0.2 lower than in matched pairs, while C_ODI stays within 0.1 |
| P4 | In TE-only pairs at L0, C_NDI is lower than in matched pairs |
| P5 | In matched pairs, σ_NDI at L1 exceeds that at L0, and the added spread √(σ²_NDI(L1) − σ²_NDI(L0)) (L1 cases with ρ = 0) lies within ±50 % of the stage-3 estimate, 0.397 × the root mean square of s_T2 |
| P6 | At L3, the NDI component of u_relax reverses sign relative to L0 |

An interpretation of a possible P3 failure was recorded before validation.

#### S-M5.8 Validation and identifiability results {.unnumbered}


The judge was evaluated against the validation criteria fixed before validation (S-M5.5). Calibration (criterion i) passed: no stratum of 20 failed, in either arm; the reference comparison for prior leakage (criterion ii) passed: (a) 0.007 against ≤ 0.10, (b) 0.983 against ≥ 0.9, (c) 0.616 against ≤ 0.65; the misspecification score (criterion v) did not meet its criterion: the largest per-stratum false-alarm rate was 0.024 against 0.02 (range 0.006–0.024), and the threshold was not changed; injected misspecification was detected at rates of 0.014 (b-value scaling), 0.031 (intrinsic diffusivity 3.4 μm²/ms), 0.010 (N_eff 400) and 0.062 (SNR 10). Validation predictions P1 and P2 held in both arms; of P3–P6, P3 failed in both arms and P4–P6 held (Supplementary Fig. S2).

## From S-M5.6 Linear–Gaussian baseline and Fisher sensitivity

A linear–Gaussian measure C^lin^ (a local linear regression of the observation mean on θ*, its Fisher information, nuisance parameters marginalized by the Schur complement, Gaussian priors) was computed on the same simulations and compared only with the uniform-mixture posterior, which also marginalizes over k, as a measure of non-linearity. A Fisher sensitivity analysis added a global b-value scale, gray-matter partial volume and a compartmental proton-density ratio, each with five prior widths set by design; its ratios and absolute values are not reported, because stage 1 showed that linear measures overstate identifiability.

## From S-M5.8 Validation and identifiability results

### Prior-dependent NDI cells

In the other 27 pairs, NDI reached C ≥ 0.75 under the reporting prior but changed class at half the prior width (none at twice the width), so these cells are prior-dependent and not counted as identified. This means that the NDI posterior was not narrow relative to the narrower prior, not that the data carried no information about NDI: halving the prior lowers C unless the posterior narrows in proportion (that is, unless it is set by the prior), and C rose at twice the width (0.855 to 0.927 in the stable pair). The stable pair passed at half the width with C = 0.753 against the threshold of 0.75 in both arms. The specification does not fix how the 200-sample threshold is applied at half the width: the frozen code takes the median of C over the observations that keep at least 200 posterior samples, and in 223 of the 224 cells the median remaining count was below 200 (S-M5.5; Supplementary Table S1). Under the other handling the specification leaves open, treating such cells as not reweightable, the NDI count would be 28 pairs (FISO 18; 27 in C7). We report the count under the frozen code and rest no conclusion on either count.

## From S-M5.8 (identifiability details)

At the reference condition (L0, SNR 30, N_eff 50), the judge identified paired ODI differences in all 28 card pairs, NDI differences in one pair and FISO differences in none (Fig. S3a). The pair with an identified NDI difference, OXF1PRI~OXF2PRI, joins two Siemens Prisma scanners with 32- and 64-channel head coils that share TE, TR and direction scheme (Warrington et al., 2025). In the other 27 pairs, NDI reached C ≥ 0.75 under the reporting prior but changed class at half the prior width (none at twice the width), so these cells are prior-dependent and not counted as identified; this means that the NDI posterior was not narrow relative to the narrower prior, not that the data carried no information about NDI. The stable pair passed at half the width with C = 0.753 against the threshold of 0.75 in both arms. The specification does not fix how the 200-sample threshold is applied at half the width: the frozen code takes the median of C over the observations that keep at least 200 posterior samples, and in 223 of the 224 cells the median remaining count was below 200 (S-M5.5; Supplementary Table S1). Under the other handling the specification leaves open, treating such cells as not reweightable, the NDI count would be 28 pairs (FISO 18; 27 in C7). We report the count under the frozen code and rest no conclusion on either count. The reasoning on halving the prior is in the deposit (`docs/supplement_details.md`, section "From S-M5.8 Validation and identifiability results"). The C7 arm identified the same 29 cells. Two results do not depend on the half-width step. First, at the reference condition the NDI MDE was 0.056–0.068 across the 27 pairs with acquisition differences and 0.060 in the matched pair (C7: 0.054–0.070 and 0.056), so, at the level of the design and at the reference condition, an acquisition difference changed the precision of a paired NDI comparison little; for FISO and ODI, MDE was 0.034–0.042 and 0.020–0.022 across the 28 pairs and 0.040 and 0.021 in the matched pair (C3, worst of the nine k; Fig. S3b). Second, under either handling, no pair had an identified NDI or FISO difference at L3, and ODI was identified in all 28 pairs at every relaxation level. Under the frozen code, the cells identified and not prior-dependent at L1, L2 and L3 were, for NDI, 0, 2 and 0 of 28 (C3; 0, 4 and 0 in C7), every NDI cell being partial at L3, and for FISO none at any level; we do not interpret differences between these small counts. At the reference condition, 27 NDI, 28 FISO and no ODI cells were prior-dependent, and only FISO cells, 10, were assumption-dependent (C3; C7: 27, 28 and 0 prior-dependent, 1 FISO assumption-dependent). The linear-Gaussian baseline differed from the judge by a median of 0.295 in C over all card cells (C_mix − C^lin^, C3; 0.178 for C7), reported as a measure of non-linearity only (S-M5.6).

## From S-M6 ROI values (atlas releases)

The label image is the 48-label release of the JHU ICBM-DTI-81 atlas (SHA-256 fac584ec…), which is not redistributed; it is the label image of FSL `data_atlases` tag FSL6-0-1, and `reproduce/get_atlas.py --from <file>` checks its SHA-256. The copy of this atlas in the FSL installation used for processing (SHA-256 3c3f5d2f…) is a different release: it differs in 5,850 voxels, carries other label numbers for part of the external capsule and for the uncinate fasciculus and the tapetum, and has label numbers up to 50; the repository's provenance record notes that FSL releases from 6.0.5 on carry 50 labels. Only the repository copy was used.

## From S-M6 Processing (failures)

One session, 16793_OXF4GEP001, failed both attempts at the same step: merging the b0 volumes before topup stopped because the reversed-phase b0 did not match the other b0 volumes in size (fslmerge: size mismatch along a non-concatenated dimension). The provider documents this session's reversed b0 as acquired along left–right (Table S1 (c)). As the rule requires, it was excluded and not repaired. In all, 119 of the 120 ON-Harmony sessions were processed, together with all 36 Tong sessions and the data of all three Wongkornchaovalit participants. Processing was interrupted once when the shared disk filled (2026-10-28, about 80 min); the interrupted sessions were restarted and, as fixed, not counted as failures.

## From S-M7.1 Fixing and release (order of the computations)

The release commit also fixed the order of the judgment computations: V4 first, then the ON-Harmony pair tests, with the Tong and Wongkornchaovalit judgments to run as their processing finished. Within V4 the program enforces the order (session SNR table → attribute table → test-data coverage → C3 and C7 arms) and stops a step if an earlier hash does not match; across programs the order was operational, not enforced by code, and is shown by the times in the access record.

## From S-M7.1 Fixing and release (access-record file names)

(`prereg/artifacts/decomp/onharmony_image_access_log.json`; the directory name `prereg/` is kept from the earlier study's repository and holds documents fixed before computing, not a registration with a registry; external data in `external_data_access_log.json`; released in the deposit as `records/access_records.csv`, S-M1).

## From S-M7.2 V4: rules, commits and times

### Clock times

(b625516, 2026-10-26 00:17:46 KST, and b886d76, 01:00:53 KST; validation results first opened at 09:04:55 KST), and points it left open were fixed in the judge code after those results were opened and before the real-data rules were released (c8825e2, 2026-10-27 19:11:18 KST).

## From S-M7.2 (time of commit c8825e2)

The rules were applied unchanged; among the points left open that the judge code fixed (above) was that a unit of a device without a recorded coil is excluded as misspecified only when every coil combination exceeds the threshold (commit c8825e2, 2026-10-27 19:11:18 KST)

## From S-M7.5 ON-Harmony pair tests (sensitivity analyses)

### Provenance of the four ROI units

(LOCK_B §177-6, commit 8ef9044, 2026-09-22 17:02:06 KST; that study's access log for its sealed ON-Harmony image folder, written after only the files' sizes and SHA-256 had been computed, is stamped 2026-09-27 01:17:49 KST, i.e. 2026-09-26 16:17:49 UTC)

## From S-M8.1 (details)

Computed after these results had been seen, the mean x difference was 0.30 times the retest SD of a paired difference on OXF1PRI alone (σ_d = 0.043, at a higher session SNR than OXF2PRI and without a degrees-of-freedom adjustment), and the SD of the unit differences, 0.061, was 1.40 times σ_d. For NDI, over the units used for coverage, the standardized error of the posterior mean had mean −4.43 and SD 1.49 across units, against an SD of 0.94 in test data with the same coil combination. The simulator's expectation was computed at the measured session SNR of the two scanners, 18.6 and 15.6; with equal coils it depended on the regression model (+0.71, +0.45 and −0.12) and is not robust. Under equal noise-floor settings, mean standardized errors were 0.13–0.60 and their SD 0.84–0.98; the SDs were close to the test-data reference and the ODI means less so, and because the measured difference is near zero, agreement in location under equal settings follows almost by construction, but the spread also fell to the reference level.


In a normal approximation chosen after these values had been seen (first reported 2026-11-03 01:00:56 KST), which splits the coverage deficit into shift and excess spread, the shift alone would bring coverage to 0.00 and the excess spread alone to 0.73, also below the cell's threshold, so, in this approximation, the shift contributed more than the spread.



For ODI the two contributions were comparable (mean −1.18 and SD 1.44 against 1.02; 0.67 and 0.75). Both mean errors were negative (after the verdict; Table S6, No. 19): the sign rule was defined at 2026-11-02 18:02:26 KST, and the expected direction, that of frozen test data with the second scanner's coil assumed to be 64 instead of 32 channels (−3.62 and −2.65), was recorded at 18:08:15 KST, the second the resampling started (Table S6, No. 19). Two signs agree by chance with a probability near one half.

Measured difference across ROIs (exploratory):


Over the 191 units used by the judge, the measured logit-NDI difference rose slightly toward superior ROIs (Spearman ρ with the superior–inferior coordinate +0.18, subject-bootstrap 95 % interval 0.03 to 0.34; within subjects +0.22) and toward weaker b = 2000 signal (−0.19, −0.35 to −0.00); the ODI difference showed neither. Both trends lie within the scatter of the difference (SD 0.059).

## From S-M8.2 (details)

Exploratory (Table S6, No. 24): in repeated OXF1PRI sessions, the scan–rescan dispersion of x and ODI was about half of that simulated at SNR 30 (median ratios over ROI units 0.55 and 0.58), whereas the dispersion of y, which depends on the free-water level, was larger (2.20).


In the frozen test data, assuming the wrong coil at one scanner moved the mean standardized error by about ±2.6 to ±3.7 and lowered coverage from about 0.9 to about 0.3. In the real units of OXF1PRI~OXF2PRI, the response to the corresponding change differed: moving from the recorded to equal settings shifted the mean standardized error by +4.6 to +4.7 for NDI and +1.8 for ODI, against +3.0 to +3.7 and +2.7 to +2.8 expected from the test data (post hoc; Table S6, No. 25). Across the 23 ODI cells, the mean standardized errors were described by an intercept and one term per scanner (R² 0.97; 0.93 when each cell was left out in turn), although cells share scanners and are not independent. A coil-dependent noise floor and other scanner properties predict this additive structure equally, so it does not separate them.

The post hoc prediction that coil combinations including the 64-channel setting would cover zero less often than the 32-channel-only combination did not hold (Table S6, No. 17).

Apparent FISO against the judge's test data and the simulations at the OXF1PRI acquisition: S-M8.5.

## From S-M8.3 (hash and time of the definitions file)

E1 and E2 are exploratory: they were defined, with their readings, in a definitions file outside the fixed document (sha256 239f684a…, 2026-10-31 23:30:15 KST), after the main multi-TE results were known.

## From S-M8.4 Noise floor at b = 2000 s/mm² (times of the readings)

Siemens images were read on 2026-11-04 00:01:33 – 00:01:50 KST, and Philips and GE images after their known-answer tables were fixed (2026-11-04 00:21:41 KST), on 2026-11-04 00:21:48 – 00:22:03 KST.

## From S-M8.4 Noise floor at b = 2000 s/mm² (times of the re-reading rules)

with agreement tolerances defined before the re-reading (2026-11-04 00:45:19 KST) and with the truncation point, the same-voxel baseline, the rule for dropping the caveat on c₁ and the rule for the note on zeros amended at 2026-11-04 01:36:42 KST, before the re-reading began.

## From S-M9.1 and S-M9.5 (clock times and records)

#### S-M10.1 Environment and records {.unnumbered}

The analyses, rules, readings and code were committed to the development repository before any of their new quantities was computed (2026-11-14 18:22:45 KST; `decomp/rules_env/ANALYSIS_PLAN.md`), with one change made before computing (18:27:55; S-M10.3 and Table S1 (b), No. 11). The results were computed from 18:28 to 18:33 and recorded in the working folder's version history before they were read (18:33:10; committed to the development repository with this version); until then only the start and end times of each step and the sizes of the result files had been displayed, so the outcome-specific wording, written and recorded at 18:30:03 after S1a and S1b had been computed, preceded any reading of their values; every ROI file was logged in the re-read record before it was read (185 files). Main-text numbers from the environment are tied to ledger rows; the result numbers in the text of S-M10.4 and S-M10.5 were transcribed from the result files and are recomputed in context, and Table S17 regenerated line by line, by `decomp/manuscript/checks/check_env_numbers.py` (with `make_S17.py`), reading the result files listed with their SHA-256 in `decomp/results/committed_sha256.json`.

#### S-M10.5 Analyses after the scale comparison {.unnumbered}

After the scale comparison of S-M10.4 had been seen, three further analyses were defined; their inputs, readings and the wording for each outcome were committed to the development repository before any of their quantities was computed (2026-11-15 03:52:36 KST; `decomp/rules_env/ANALYSIS_PLAN.md`, section 7), computed from 03:52:56 to 03:52:59 and recorded in the working folder's version history before they were read (03:52:59).

## From Table S1 (a) deviation (iii)

### Driver hashes and times

For (iii): the command-line run stopped on 2026-10-30 at 18:09:36 KST without writing a result; the driver, whose SHA-256 is 929c03d10668d6e0788bfaf7e416f30d33f79cfdb6e6bf5b8315c08567201365 (computed on 2026-11-04 from the file `records/x12_tong_driver.py`, last modified 2026-10-31 12:56:26 KST; in the deposited copy its one server path and the author's name are replaced, and `MANIFEST.csv` lists both hashes; its first 16 hexadecimal digits match those recorded in the work log at the run), wrote the result file once, on 2026-10-31 at 12:56:44 KST (03:56:44 UTC), and the access record holds that time and the result file's SHA-256 (c0e473996ddc5a53…).

## From Table S1 (a) C7 rerun

### Clock times

(four elements clamped, none beyond float32 rounding level; reading fixed 2026-11-01 14:53:03 KST, run 2026-11-01 15:21:42 KST – 2026-11-03 00:12:50 KST, committed 2026-11-03 01:00:48 KST)

## From Table S1 (b) (deviation records in full)

**(b) Deviation records of the post hoc and exploratory analyses.**

| No. | Recorded (KST) | Change or event | Relative to the reads and computations | Verdict; where reported |
|---|-------------|------------------------------------------------------|----------------|----------|
| 1 | 2026-11-02 10:49:42 | Order. Because the thresholds matched to the applied rule had priority, the test-data side of D1 (no real data) was started before the C7 results were committed, whereas the plan had placed D1 after that commit; the real-data side of D1 stayed after it. | After the verdict; no image read | Verdict unchanged; S-M9.2 |
| 2 | 2026-11-02 18:08:15 (author decision) | Order. D1 was started while the C7 rerun was running instead of after its commit; D1 does not use C7 results (it uses the C3 arm, its stored output and the test data). | After the verdict; no image read | Verdict unchanged; S-M9.2; S-M8.1 |
| 3 | 2026-11-03 21:40:16 | Head mask. In a timing run of the known-answer simulation (one simulated OXF1PRI session, SNR 16, L 6 and 8), the head-mask threshold, 0.2 × the 99th percentile of the whole image, took the background floor as head, so that the central box grew and 819 and 1,129 mask voxels lay outside the brain. The threshold was changed to Otsu's method on the mean b0, with the rest of the mask definition unchanged; after the change, no mask voxel lay outside the region in three simulated sessions of the same scanner at L 0.5–12. | Before any raw pixel value was read (headers only) | S-M9.1; S-M8.4 |
| 4 | 2026-11-03 23:49:22; author's approval 23:59:08 | Rules before the b = 2000 s/mm² images were read, made after the known-answer results: (i) the rule defined first, on differences in median R against a minimum detectable difference, was moved to the supplement, and the reading was set by the rule on log₂ effective L; (ii) sessions with σ_q < 3 in stored units were excluded, where the definition had only flagged stored σ below 10; (iii) the calibration range was widened from the grid L 0.5–12 to a pure-noise curve over L 0.5–128, with values outside it treated as censored intervals; (iv) descriptive diagnostics were added (R with zeros included, inner against rim voxels, outlying volumes, floor level), not used for the reading. | Before any pixel was read | S-M9.1; S-M8.4; Table S7 (a) |
| 5 | 2026-11-04 00:21:41 | Order. The Philips and GE known-answer tables were fixed after the Siemens results had been seen; these scanners had been set as descriptive units in the definition of 2026-11-03 21:23:18, and no other choice was added when their tables were fixed. | Before the Philips and GE pixels were read (00:21:48) | S-M9.1; deposit (S-M1) |
| 6 | 2026-11-04 01:36:42 | Rules set before the re-reading (01:36:51): the caveat on c₁ is dropped only when both Prisma scanners agree and the paired differences satisfy \|D_skew\| ≤ ⅙ and \|D_ML\| ≤ 0.05; agreement is judged against the mean-to-SD estimate on the same voxels and values; the maximum-likelihood arm uses a likelihood left-truncated at the fifth-smallest positive value; the labels read "agrees within this sensitivity" and "differs from the central-chi shape (cause not separated by this check)"; Philips and GE receive no agreement label; the rule for the note on zeros is set; and the gap index is undefined when the expected counts of the two lowest bins are below 20. A rule that gives no value where fewer than 50 voxels are available, absent from the definition of 00:45:19, is recorded as a deviation. | Before the re-reading (01:36:51) | S-M9.1; S-M8.5 |
| 7 | 2026-11-04 10:18:45 | Known-answer simulation as defined. The definition named rigid motion from the eddy parameters of a recorded session (three translations, three rotations); the code that produced the results also applied ten eddy-current terms. The reading was made with that code and stands; the simulation was rerun as defined, with the eddy-current terms at zero, for the three Siemens scanners with the same seeds and numbers of simulations. | After the reading; no image read | Deposit (S-M1) |
| 8 | 2026-11-06, about 20:25 (disclosed in the work record of main-text Fig. 4) | Read without a prior record. While looking for stored values by ROI unit for the display of main-text Fig. 4, the contents of 54 JSON files of the analyses after the verdict (`data/remeasure/explore`) were searched for the string «tong», with file names as the only output; no file matched, no content was displayed and no value was used. No record was added afterwards, because adding one would have meant reading the 54 files again. | After the verdict; before the shares of main-text Fig. 4 were first computed (20:38:47); no image read | Verdict unchanged; this table |
| 9 | 2026-11-06, about 20:51 (disclosed in the work record of the judge-matrix figure, Fig. S5) | Read without its own record. To check a rounding tie (cell V4-19, coverage 0.3125; Table S3), the ledger copy of the verdict file was read once more; the read relies on the record of the same file with the same SHA-256 about one minute earlier and has no line of its own. | After the verdict; ledger copy only; no image read | Verdict unchanged; this table |
| 10 | 2026-11-07, with the counts of Table S7 (c) | Wording. The sentence recorded before counting for the case «mostly in the lateral ventricles» of the location rule (Table S6 (b), D3) ended «although on the slices displayed it lay at the brain surface». The revised Fig. S4 shows neither a slice nor a single participant, so the clause reads «although not in every session» (S-M8.5). The participant chosen for the earlier single-slice display had the lowest lateral-ventricle shares of its scanners' sessions (0.000, 0.012 and 0.013), which is why the slices had shown the mask at the brain surface. | After the counts; no image read | Verdict unchanged; S-M8.4, S-M8.5, Table S7 (c) |
| 11 | 2026-11-14 18:27:55 (working-folder version record 2a6c95a; committed to the development repository with this version) | Change before computing to the environment's S2 (S-M10.3): in the design committed first, the drawn relaxation of 20 simulated participants per group and site entered the change in the group difference with an error of about 0.026 on the logit scale, the size of Eq. (7) itself; the primary analysis was set to the population mode of Proposition 4 and the drawn design kept as a secondary analysis compared with the algebra of the same participants. | Before any S2 quantity was computed (S2 computed 18:29:42–18:31:37); before 18:27:55 the S2 code had run only in unit tests on algebraic and synthetic inputs | Verdict unchanged; S-M10 |
| 12 | 2026-11-14, before 14:43:18 (recorded afterwards at 14:43:18) | Read without a prior record. While planning the environment, the keys and the first 300 characters of each value of the frozen retest-SD file of the earlier study (`prereg/artifacts/unseal/s_ref.frozen.json`, including the NDI retest SD of five ROI units) were printed before the re-read record was written; recorded afterwards; sha256 unchanged (39e987ab…, the hash fixed in the analysis document). | After the verdict; used only for an order-of-magnitude remark in a planning note; no number in the paper comes from it | Verdict unchanged |
| 13 | 2026-11-14, before computing | Implementation. `neuroCombatFromTraining` (neuroCombat 0.2.12) fails with NumPy 2.5.3 (conversion of a one-element array to an integer); the same formula was implemented (`decomp/rules_env/analyses.py`, `combat_apply`) and reproduced neuroCombat's own output for the training data within 10⁻¹⁰, with location terms only and with location and scale terms. | Before computing | S-M10 |
| 14 | 2026-11-15, about 03:27 (before the analysis plan of S-M10.5 was committed) | Value seen before the plan. While the unit tests of the analyses after the scale comparison were written, a failing assertion printed one value of S3 (the raw-scale change, 0.0025, for control NDI 0.46, a group difference of −0.10 and TE 100 against 70 ms); the tests were then changed to print no values, and the event is recorded in the plan. | S3 is algebra with no data and no reading | No reading changed; S-M10.5 |
| 15 | 2026-11-15 (code review after the results; the fix committed before recomputing, 18:02:54) | Implementation error in the pooled ΔTE* of E1 (S-M10). The finite-n factor of the rule fixed before computing uses n, the units of the group from which the offset is taken (27 per participant); the pooled computation used the 81 pooled units. Recomputed with the rule: x 18.0 → 17.4 ms (90 % interval 15.6–21.4 → 15.1–20.7), y 16.8 → 16.3 ms; the per-participant values and the ON-Harmony value were not affected; the first-run values are kept in the result file `e1.json`. | After the results; reading unchanged | S-M10.4; Section 4.3 |
| 16 | 2026-11-15 (code review after the results) | Checks of the text. The repository check of S-M10 did not include the numbers of S-M10.5, and its comparison of Table S17 did not detect added rows; both checks were revised and the numbers, which had been checked in the working folder, were unchanged. The checks added to the analysis code (input hashes and all fields of the E6 job list; a stop for covariate-trained ComBat estimates) did not change any result; the E6 inputs matched the hashes fixed in the plan. | After the results | S-M10.1 |
| 17 | 2026-11-16 10:46:22 (amendment committed; author decision) | Rerun after a failed control. The registration check of the multi-TE series defined on 2026-11-16 (S-M10.6) failed its first failure control under the setting of its plan, which fixed that the analysis then stop; instead, one setting was amended (`-nosearch`), committed before a second run, and the check was run once more. The control failed again, and the analysis ended as the amendment fixed. | After the failure of the first run had been seen; first read of the second run 10:46:47 | No result added or changed; S-M10.6; Section 4.3 |

## From Table S3

### Search 3

| Date | Question | Searched | Could not be opened, or not searched | Conclusion within what was opened |
|---|-------------|-----------|--------------|---------------------|
| Search 3, 2026-10-25 | Bibliographic verification of 28 references | OpenAlex, Semantic Scholar, Europe PMC, PubMed and PMC where they loaded; publisher and repository pages | ScienceDirect, Wiley and some PubMed and PMC pages blocked | All 28 found |

### Search 7

| Date | Question | Searched | Could not be opened, or not searched | Conclusion within what was opened |
|---|-------------|-----------|--------------|---------------------|
| Search 7, 2026-10-26 | Verification of 67 references | PubMed (55 items); a preprint connector (1); publisher, repository and bibliographic pages (11) | The reference-lookup tool failed (certificate error) and judged no item | 64 verified, 3 partial, 0 mismatch or not found |

### Check of eight citations, 2026-11-03

| Date | Question | Searched | Could not be opened, or not searched | Conclusion within what was opened |
|---|-------------|-----------|--------------|---------------------|
| 2026-11-03 | Bibliographic data and content of eight citations (a check of chosen works, not a search) | PubMed; OpenAlex; seven full texts and one abstract | Supplementary files and the publisher page of Jodoin et al. (2025) | — |

### Search 10

| Date | Question | Searched | Could not be opened, or not searched | Conclusion within what was opened |
|---|-------------|-----------|--------------|---------------------|
| Search 10, 2026-11-04 | Check of the reference list (95 entries) | PubMed; page fetch; web search; a preprint connector | IEEE Xplore, DBLP, HAL, Project Euclid, OpenReview and some publisher pages blocked; 10 entries with DOIs not checked | PubMed part: 59 of 59 verified; other part: 23 verified, 0 partial, mismatch or not found, 3 unchecked |

## From Table S3 (searches of 2026-10-02 and searches 1, 2 and 4–6)

| Date | Question | Searched | Could not be opened, or not searched | Conclusion within what was opened |
|---|-------------|-----------|--------------|---------------------|
| 2026-10-02 | Whether harmonization methods had been evaluated by coverage or calibration (a search of the earlier registered study on a related question) | Web search (34 queries); Europe PMC (18 queries; keyword checks of 19 open-access full texts); PubMed (4); a preprint-server abstract search (10); OpenReview (4) | Semantic Scholar not run (HTTP 429); not searched: full ISMRM and OHBM abstract books, IEEE Xplore, ACM DL, literature not in English | — |
| Search 1, 2026-10-25 (file date) | Measurement invariance (MI) applied to brain-imaging indices; TE/TR as grouping axes | Web search; direct opening of pages | PMC, PubMed, scinapse.io and literature APIs blocked, so no systematic database search; stopped early in part | No study treating NDI or FISO as acquisition-indexed quantities, with TE/TR as the grouping variable, predicting metric invariance with scalar non-invariance on the logit scale |
| Search 2, 2026-10-25 (file date) | Relaxation and TE/TR dependence of NODDI indices; how pathology papers interpret NDI decreases and FISO increases | Web search; direct opening of pages | PMC, PubMed article pages, ScienceDirect, AHA Journals, OUP PDFs and the Europe PMC API blocked | No study in human patients running standard and multi-TE NODDI side by side; no primary study varying TR alone and measuring FISO or the free-water fraction |
| Search 4, 2026-10-26 | MI or differential item functioning (DIF) applied to imaging indices; assumptions of harmonization | Web search; page fetch; Europe PMC (three queries); Crossref | PMC, PubMed and ScienceDirect blocked; Google Scholar, PsycINFO, Scopus and Web of Science not used | No study using the configural/metric/scalar hierarchy or DIF to judge which cross-acquisition comparisons of diffusion microstructure indices are interpretable; no neuroimaging paper reading ComBat or diffusion harmonization as MI/DIF |
| Search 5, 2026-10-26 | Multi-site NODDI differences attributed to TE/TR; which comparisons remain interpretable as tissue differences | Web search; ISMRM abstract archives; publisher and repository pages | PubMed, PMC and Europe PMC blocked; Google Scholar not used; OHBM and MICCAI abstracts not searched | No paper working out which comparisons remain interpretable, formalizing that TE shifts the logit intercept but preserves the slope, or linking site TE × group T2 to ComBat's common-slope assumption (Gong et al. (2020, Eqs. 8–9) give both TE log-ratio relations; this search did not identify them as such) |
| Search 6, 2026-10-26 | Simulation-based inference conditioned on acquisition; identifiability of differences; posterior calibration | Web search; page fetch; Crossref metadata; ISMRM 2026 abstract pages | Europe PMC, IEEE Xplore and MICCAI/IPMI not queried directly; no systematic date-filtered PubMed or Europe PMC search | No test of the posterior coverage of a frozen simulation-trained estimator on traveling subjects with a true difference of zero |

## From Table S5 (a)

### Definition and first-computation times

| No. | Analysis | Grade | Defined | First computed or reported | Reported in or status |
|---|---------------------|-----------|----------|------|----------|
| 1 | Bias against spread, Prisma pair | post hoc | 2026-11-01 22:33:57 KST | 2026-11-03 01:00:56 | S-M9.2; S-M8.1 |
| 2 | Thresholds matched to the applied rule | post hoc (order with the start of the test-data draws not determinable) | 2026-11-02 10:49:42 KST | 2026-11-02 17:37:46 | S-M9.2; Fig. S5; Table S3 (b) |
| 3 | Simulator expectation of the measured difference; judge rerun with equal noise-floor settings | post hoc | 2026-11-03 02:19:46 KST | 2026-11-03 02:20:20 (expectation); 02:19:55–02:43:33 (judge) | S-M9.2; Table S16; Table S11 |
| 4 | Measured difference over all units, with a subject bootstrap | exploratory (after No. 23 had been seen) | 2026-11-03 17:26:59 KST | 2026-11-03 17:27:38 | S-M9.2; Table S16 |
| 5 | Mean x difference and unit SD relative to σ_d | exploratory | after the results had been seen | — | S-M8.1 |
| 6 | Direction of the multi-TE unit trajectories (E1) | exploratory (after the main multi-TE results) | 2026-10-31 23:30:15 KST | 2026-10-31 23:33:10 | Section 4.3; S-M8.3 |
| 7 | TE slopes: ratio to the literature prediction | exploratory | not defined in advance | 2026-11-01 | Section 4.3 |
| 8 | TE slopes of ODI (E2) | exploratory (as No. 6) | 2026-10-31 23:30:15 | 2026-10-31 23:33:10 | S-M8.3 |
| 9 | Noise floor in the data: mask, statistic R, known-answer simulation, head-mask repair | post hoc | 2026-11-03 21:23:18 KST; repair 2026-11-03 21:40:16 KST | Siemens 2026-11-04 00:01:33 – 00:01:50 KST; Philips and GE 2026-11-04 00:21:48 – 00:22:03 KST | S-M9.2; Fig. S4; S-M8.4; S-M8.5 |
| 10 | Decision rule on log₂ effective L | post hoc | 2026-11-03 23:49:22 KST | as No. 9 | S-M9.2; S-M8.5 |
| 11 | First rule, on differences in median R | post hoc | 2026-11-03 21:23:18 | as No. 9 | Replaced before any image was read; S-M8.5; Table S7 (a) |
| 12 | Re-reading: exact zeros, distribution shape, two further estimators | exploratory (after the first readings) | tolerances 2026-11-04 00:45:19 KST; amended 2026-11-04 01:36:42 KST | 2026-11-04 01:36:51 | S-M8.5; deposit |
| 13 | Known-answer test as defined (rigid motion only) | post hoc | 2026-11-04 10:18:45 | 2026-11-04 10:55 | Deposit (same reading) |
| 14 | H1′: the noise implied by the b0-based SNR exceeds the thermal noise | exploratory | 2026-11-02 11:39:21 KST | with No. 9 | Could not be tested; S-M8.4 |
| 15 | Simulations at the OXF1PRI acquisition centered on tissue reproducing the observed NDI and FISO | exploratory | 2026-11-02 10:55:57 KST | 2026-11-02 10:59:51 | S-M8.5 |
| 16 | Apparent FISO at seven scanners against the test data | exploratory | not defined in advance | 2026-11-02 11:39:21 | S-M8.5 |
| 17 | Prediction on coil combinations including the 64-channel setting | post hoc | 2026-11-02 12:27:49 | 2026-11-03 01:00:56 | Did not hold; S-M8.2 |
| 18 | Normal approximation splitting the coverage deficit | exploratory | not defined in advance | 2026-11-03 01:00:56 KST | S-M8.1 |
| 19 | Sign of the mean standardized errors | post hoc (expected direction recorded the second the resampling started) | 2026-11-02 18:02:26; direction 18:08:15 | 2026-11-03 01:00:56 | S-M8.1 |
| 20 | Within-session spread of the per-voxel ratio | exploratory | after the results had been seen | 2026-11-04 02:16:39 | S-M8.5; deposit |
| 21 | Between-session spread of R against the known-answer test | post hoc | 2026-11-03 23:49:22 | with No. 9 | Deposit |
| 22 | Synthetic test of the shape check | exploratory | not defined in advance | 2026-11-04 01:36:42 | Deposit |
| 23 | Measured logit-NDI difference by ROI position, Prisma pair | exploratory | not defined in advance | 2026-11-03 02:23:15 | S-M8.1 |
| 24 | Scan–rescan dispersion in repeated OXF1PRI sessions against simulation | exploratory (the OXF1PRI retest SD of ODI had been read, 2026-10-25) | before computing (plan of 2026-11-01 22:33:57) | 2026-11-01 23:16 | S-M8.2 |
| 25 | Response of the judge to its coil term; additive scanner model of the ODI mean errors | post hoc | 2026-11-02 18:02:26; 2026-11-03 02:19:46; 2026-11-02 10:49:42 | 2026-11-02 18:03; 2026-11-03 02:19:55; 2026-11-03 01:05 | S-M8.2 |
| 26 | Noise-free grid simulation of the multi-TE slopes | exploratory | 2026-11-02 10:55:57 | 2026-11-02 10:57:21 | S-M8.3 |
| 27 | OXF3TRI ODI cells with their recorded coil combinations | exploratory (cells grouped after the V4 results) | 2026-11-02 18:02:26 | 2026-11-02 18:03:45 | Table S11 |
| 28 | B1-S: simulation of the multi-TE processing (flattening of the TE trajectories) | exploratory | 2026-11-01 23:53:38 | 2026-11-02 00:54 | Conclusion invalidated (2026-11-02 10:51:38) |
| 29 | D3: joint multi-TE fit of the compartmental T2 ratio | exploratory (order not determinable) | plan of 2026-11-01 22:33:57 | — | Conditions not met |
| 30 | D3 extension: known-answer test with true ratios up to 3.0 | exploratory | 2026-11-01 23:53:38 | 2026-11-02 00:12 | Conditions not met; no real-data fit |
| 31 | D3-v3: continuous profile with systematic-error variants | exploratory | 2026-11-02 01:15:58 | 2026-11-02 02:47:33 | Conditions not met; conclusion invalidated (2026-11-02 10:51:38) |
| 32 | D5 (3): unit-level outliers within each scanner | exploratory (rule R21) | 2026-11-02 09:40 | 2026-11-02 10:02 | Conclusion invalidated (2026-11-02 10:33) |
| 33 | D5 (4): concentration of the judge's failures in the units of No. 32 | exploratory | 2026-11-02 09:40 | not computed | Invalidated with No. 32 |
| 34 | Environment S1a: offset-subtracted retest share at fitter output, reference tissues (stored stage-3 fits) | exploratory (closely related results had been seen) | 2026-11-14 18:22:45 (committed, 03631a0) | 2026-11-14 18:28:33–18:28:49 | S-M10; Table S17 |
| 35 | Environment S1b: fitter output with ROI-matched tissue states | exploratory (closely related results had been seen) | 2026-11-14 18:22:45 (committed, 03631a0) | 2026-11-14 18:28:49–18:29:42 | S-M10; Table S17 |
| 36 | Environment S2: Proposition 4 at fitter output (population mode after Table S1 (b), No. 11) | exploratory (closely related results had been seen) | 2026-11-14 18:22:45 (committed, 03631a0); changed 18:27:55 | 2026-11-14 18:29:42–18:31:37 | S-M10; Table S17 |
| 37 | Environment K8: FISO signs predicted with the fitter | exploratory (closely related results had been seen) | 2026-11-14 18:22:45 (committed, 03631a0) | 2026-11-14 18:31:37–18:32:34 | S-M10; Table S17 |
| 38 | Environment E1: TE-slope heterogeneity and ΔTE*, multi-TE series (pooled ΔTE* recomputed, Table S1 (b), No. 15) | exploratory (closely related results had been seen) | 2026-11-14 18:22:45 (committed, 03631a0) | 2026-11-14 18:32:34–18:32:41 | S-M10; Table S17 |
| 39 | Environment E2: offset-subtracted shares, ON-Harmony | exploratory (closely related results had been seen) | 2026-11-14 18:22:45 (committed, 03631a0) | 2026-11-14 18:32:41–18:32:42 | S-M10; Table S17 |
| 40 | Environment E3: offset-subtracted shares, ten centers | exploratory (closely related results had been seen) | 2026-11-14 18:22:45 (committed, 03631a0) | 2026-11-14 18:32:42–18:32:44 | S-M10; Table S17 |
| 41 | Environment E4: scale comparison | exploratory (closely related results had been seen) | 2026-11-14 18:22:45 (committed, 03631a0) | 2026-11-14 18:32:44–18:32:54 | S-M10; Table S17 |
| 42 | Environment E5: slope structure and the scale comparison | exploratory (defined after No. 41 had been seen) | 2026-11-15 03:52:36 (committed, 2c7d2b6) | 2026-11-15 03:52:56–03:52:57 | S-M10.5; Table S17 (d) |
| 43 | Environment E6: scale comparison with shared relaxation (stored fits) | exploratory (defined after No. 41 had been seen) | 2026-11-15 03:52:36 (committed, 2c7d2b6) | 2026-11-15 03:52:57–03:52:59 | S-M10.5; Table S17 (d) |
| 44 | Environment S3: raw-scale site × group term (algebra) | exploratory (defined after No. 41 had been seen; one value seen before, Table S1 (b), No. 14) | 2026-11-15 03:52:36 (committed, 2c7d2b6) | 2026-11-15 03:52:59 | S-M10.5; Table S17 (d) |
| 45 | Values from the signal model and the fixed rules: size of Proposition 2, tolerance of the 0.90 criterion in NDI, expected multi-TE results under repeat noise only, power of the group × TE term | exploratory (closely related values had been seen) | 2026-11-16, in a plan outside the repository, before computing | 2026-11-16 20:08:48 (first run; later runs, the last on 2026-11-17 00:56:10, changed only the recorded description of the grade) | Sections 2.2, 4.1, 4.3 and 5.2; S-M2; S-M10.7 |
| 46 | Tissue-T2 path on apparent FISO | exploratory (closely related values had been seen) | 2026-11-17, in a plan approved before computing | 2026-11-17 00:54:27 (recomputed at 01:06:59 to record decreases rather than levels; same values) | Section 5.2; S-M10.7 |

## From Table S5 (b)

### Definition and first-read times of the displays

| No. | Display | Grade | Defined | First read | Status; where reported |
|---|--------------------------------|-------|-------|----------|---|
| D1 | Fitted maps (Fig. 1): multi-TE participant sub-01, NDI, FISO and ODI of the series TE62, TE92 and TE132; one axial template slice; maps from the NODDI fit (Section 3.2) resampled with the warps and interpolation of the ROI extraction; linear gray scale per index from the 1st to the 99th percentile of the displayed maps of that index within the template brain mask on the slice. Under the same rules, ROI-unit trajectories of this participant and maps of one Tong and one ON-Harmony participant were also read and rendered; they are not displayed | post hoc display: defined after the verdict and before any map was read for it; no test | 2026-11-06 20:17:09 KST (rules); choices written 20:31 KST, before the first read | First map read 2026-11-06 20:40:25 KST; ROI tables re-read 20:41:07 KST | Displayed; no verdict changed; Fig. 1 |
| D2 | Share of units at retest level by ROI unit in the negative control: the units of the magnitude test (participant × ROI unit × center pair; same rule, retest SD and degrees of freedom) grouped by ROI unit, for x, y and ODI, drawn on the full atlas labels of one axial and one sagittal slice (chosen from the atlas and the skeleton alone; values from skeleton voxels) and, for all 27 ROI units, in a strip beside them (values in Table S14) | post hoc: defined after the verdict and before its own computation (a display; descriptive; no test and no verdict by ROI unit) | 2026-11-06 20:17:09 KST (rules); slices chosen from the atlas and the skeleton alone 20:34:29 KST, before any share was computed; display on full labels with the strip defined 2026-11-07 03:14:31 KST, before the shares were recomputed | Shares first computed 2026-11-06 20:38:47 KST; pooled, they equal the shares of the magnitude test; recomputed for the revised display 2026-11-07 03:31 KST, equal to the first | Displayed; no verdict changed; main-text Fig. 4 |
| D3 | The b = 2000 s/mm² noise floor and where the noise-floor mask of No. 9 lay (Fig. S4), for every first session on OXF1PRI, OXF2PRI and OXF3TRI that entered the scanner summary of No. 9 (nothing chosen). Left: each session's stored mask, moved to template space with the session's existing FA-to-FMRIB58 transform (nearest neighbor; no new registration) and projected along z, y and x (axial and coronal neurological, sagittal anterior on the left; one scale); color, the share of the scanner's sessions with a mask voxel on the line (zero not painted); gray, the projection of the FSL MNI152 brain mask; dashed, the projection of the Harvard-Oxford lateral-ventricle label (`HarvardOxford-sub-maxprob-thr25-1mm`, both hemispheres). Right: one dot per session at its session statistic R of No. 9 (the median of the per-voxel R in its mask), read from the stored session records; the large dot, the scanner's median of these (its R_med in No. 9, checked) with the effective L of No. 9; dashed lines, the Rician floor and the assumed L | post hoc display: defined after the verdict and before any mask, transform or template was read for it; no test; no image read for the right half | Left: 2026-11-07 03:14:31 KST (rules, before any mask, transform or template was read); right: 2026-11-07 17:13:29 KST (rules) | Left: computed 2026-11-07 03:25:33 – 03:27:09 KST; the grid checks (processing grid = raw grid; three template images on one grid) passed before any mask or transform was read; right: drawn 2026-11-07 17:18 KST from the stored session records of No. 9; the median of the dots equalled the stored R_med on each scanner | Displayed; no verdict changed; where the mask lay: Supplementary Table S7 (c); Fig. S4 |

## From Table S5 (b) (display D3 in full)

| No. | Display | Grade | Status; where reported |
|---|--------------------------------|-------|---|
| D3 | The b = 2000 s/mm² noise floor and where the noise-floor mask of No. 9 lay (Fig. S4), for every first session on OXF1PRI, OXF2PRI and OXF3TRI that entered the scanner summary of No. 9 (nothing chosen). Left: each session's stored mask, moved to template space with the session's existing FA-to-FMRIB58 transform (nearest neighbor; no new registration) and projected along z, y and x (axial and coronal neurological, sagittal anterior on the left; one scale); color, the share of the scanner's sessions with a mask voxel on the line (zero not painted); gray, the projection of the FSL MNI152 brain mask; dashed, the projection of the Harvard-Oxford lateral-ventricle label (`HarvardOxford-sub-maxprob-thr25-1mm`, both hemispheres). Right: one dot per session at its session statistic R of No. 9 (the median of the per-voxel R in its mask), read from the stored session records; the large dot, the scanner's median of these (its R_med in No. 9, checked) with the effective L of No. 9; dashed lines, the Rician floor and the assumed L | post hoc display: defined after the verdict and before any mask, transform or template was read for it; no test; no image read for the right half | Displayed; no verdict changed; where the mask lay: Supplementary Table S7 (c); Fig. S4 |

## From Table S7 (c)

### Per-scanner rows, first table

| Scanner | Sessions | (1) Lateral ventricles | (2) Outside the brain mask (mostly in its unfilled ventricular gaps) | (3) Within 5 mm of its edge (mostly next to these gaps) | (4) Elsewhere | (2) + (3) |
|---|---|-----------|-----------|-----------|-----------|-----------|
| OXF1PRI | 20 | 0.75 [0.00–0.89] | 0.07 [0.00–0.24] | 0.06 [0.03–0.21] | 0.07 [0.01–0.90] | 0.16 [0.03–0.40] |
| OXF2PRI | 20 | 0.61 [0.01–0.89] | 0.06 [0.00–0.32] | 0.06 [0.01–0.19] | 0.16 [0.02–0.93] | 0.15 [0.02–0.50] |
| OXF3TRI | 10 | 0.77 [0.01–0.89] | 0.05 [0.00–0.18] | 0.10 [0.06–0.70] | 0.05 [0.00–0.29] | 0.20 [0.07–0.70] |

### Per-scanner rows, split table

| Split of (2) + (3) | Outer edge of the brain | In or next to the ventricular gaps |
|------|---------------------------|---------------------------|
| OXF1PRI | 0.02 [0.00–0.21] | 0.14 [0.00–0.34] |
| OXF2PRI | 0.03 [0.00–0.19] | 0.12 [0.00–0.33] |
| OXF3TRI | 0.00 [0.00–0.70] | 0.15 [0.00–0.31] |

## From Table S11 (Prisma pair with equal noise-floor settings)

**Prisma pair with equal noise-floor settings.** OXF1PRI~OXF2PRI with the recorded coil combination and with equal noise-floor settings at both scanners, C3 arm, reference condition (post hoc: chosen after the verdict and defined before computing, exploratory plan §24-1, 2026-11-03 02:19:46 KST; a comparison made after the verdict, not a pass of the fixed test). Units are the units used for coverage in V4 that had posterior samples for all nine k under that setting; coverage as in S-M7.2; the mean and SD are those of the standardized error of the posterior mean, (posterior mean − true difference)/posterior SD, across units. Thresholds: the fixed threshold, and the threshold recalculated from test cases whose true coil combination is the recorded one (32 | 64); no threshold was recalculated for equal settings. The last two rows are the test-data reference for the recorded combination.

| Coil setting (OXF1PRI \| OXF2PRI) | Component | Units | Coverage | Threshold: fixed · recorded combination | Mean standardized error | SD |
|--------------|--------|-----|--------|----------|-----------|-----|
| 32 \| 64 (recorded) | NDI | 191 | 0.010 | 0.80 · 0.80 | −4.43 | 1.49 |
| 32 \| 64 (recorded) | ODI | 191 | 0.649 | 0.80 · 0.75 | −1.18 | 1.44 |
| 32 \| 32 | NDI | 191 | 0.921 | 0.80 · 0.80 | +0.25 | 0.91 |
| 32 \| 32 | ODI | 191 | 0.901 | 0.80 · 0.75 | +0.60 | 0.98 |
| 64 \| 64 | NDI | 189 | 0.931 | 0.80 · 0.80 | +0.13 | 0.97 |
| 64 \| 64 | ODI | 189 | 0.921 | 0.80 · 0.75 | +0.53 | 0.84 |
| Test data, true combination 32 \| 64 | NDI | 228 cases | 0.912 | — | +0.05 | 0.94 |
| Test data, true combination 32 \| 64 | ODI | 228 cases | 0.895 | — | +0.06 | 1.02 |

## From Table S12

### Rows of 22 ODI cells

| Cell | C3 coverage | C3 units used | C7 coverage | C7 units used | C7 units without filled samples |
|---------------------|---------|------|---------|------|--------|
| OXF1PRI~OXF3TRI · ODI | 0.568 | 162 | 0.491 | 108 | 118 |
| OXF2PRI~OXF3TRI · ODI | 0.352 | 142 | 0.469 | 98 | 93 |
| OXF1PRI~NOT1ACH · ODI | 0.528 | 161 | 0.602 | 108 | 99 |
| OXF1PRI~NOT2ING-89 · ODI | 0.085 | 47 | 0.147 | 34 | 49 |
| OXF1PRI~NOT2ING-98 · ODI | 0.264 | 197 | 0.352 | 159 | 87 |
| OXF1PRI~NOT4GEP · ODI | 0.346 | 104 | 0.293 | 75 | 58 |
| OXF1PRI~OXF4GEP · ODI | 0.483 | 87 | 0.471 | 51 | 60 |
| OXF2PRI~NOT1ACH · ODI | 0.312 | 64 | 0.360 | 25 | 97 |
| OXF2PRI~NOT2ING-89 · ODI | 0.000 | 10 | 0.000 | 3 | 41 |
| OXF2PRI~NOT2ING-98 · ODI | 0.102 | 49 | 0.000 | 19 | 101 |
| OXF2PRI~NOT4GEP · ODI | 0.447 | 47 | 0.321 | 28 | 39 |
| OXF2PRI~OXF4GEP · ODI | 0.275 | 40 | 0.389 | 18 | 40 |
| OXF3TRI~NOT1ACH · ODI | 0.229 | 70 | 0.270 | 37 | 66 |
| OXF3TRI~NOT2ING-89 · ODI | 0.000 | 13 | 0.000 | 9 | 67 |
| OXF3TRI~NOT2ING-98 · ODI | 0.038 | 52 | 0.027 | 37 | 62 |
| OXF4GEP~NOT1ACH · ODI | 0.182 | 77 | 0.060 | 50 | 96 |
| NOT2ING-89~NOT1ACH · ODI | 0.083 | 12 | 0.000 | 2 | 57 |
| NOT2ING-98~NOT1ACH · ODI | 0.189 | 164 | 0.022 | 93 | 189 |
| NOT2ING-98~NOT4GEP · ODI | 0.145 | 110 | 0.016 | 63 | 123 |
| NOT2ING-98~OXF4GEP · ODI | 0.205 | 83 | 0.054 | 56 | 115 |
| NOT4GEP~NOT1ACH · ODI | 0.188 | 96 | 0.055 | 55 | 110 |
| NOT4GEP~OXF4GEP · ODI | 0.163 | 104 | 0.026 | 78 | 84 |

## From Fig. S3 caption (NOT2ING classification, code references)

For the judge, the two NOT2ING protocols carry the same TR range (4.00–4.38 s), so NOT2ING-89~NOT2ING-98 is a TE-only pair (`decomp/s4/plan.py:19–20`; `decomp/remeasure/judge_v4.py:94–101`).
