# Implementation note for the analyses of 2026-11-16 (English translation)

> English translation of `decomp/plans/S73_IMPLEMENTATION_2026-11-16.md`, committed before any of these quantities was computed (development commit 59e5de1, 2026-11-16 05:19:34 KST). The Korean file is the record, and where the two differ the Korean file governs. The definitions it implements are in `docs/s73_post_hoc_plan.md`.

# S73 implementation note — the choices the definition left open (fixed before computing)

Definition: `decomp/plans/S73_POST_HOC_PLAN_2026-11-16.md` (previous commit). This note records only the implementation choices that document did not fix; where they disagree, the definition prevails. No quantity of §2–§5 was computed before the commit of this note. The code is in `decomp/explore/s73/` and is committed with the results. The result files are in `data/remeasure/explore/s73/` (working folder `<WORK>/remeasure/explore/s73/`).

## Common

- CPU only. Fixed seeds for each analysis (below). The failure controls run first in each analysis; if they do not behave as expected, the analysis produces no values and stops (the result file holds only the controls and «stopped»).
- Read records: for derived tables (roi.json, stored fits), the re-read records written by the existing loaders are used as they are. §2 and §5, which read images and maps, masks and registrations derived from images, read only **after the author's explicit approval**, and each file is entered in a record before it is read (if the record cannot be written, nothing is read and the run stops). Files outside the approval are not read.
- Numbers enter the manuscript through ledger rows (new identifiers from N-904, numbered by A).

## §3 Tissue dependence of the shift at fitter output

- Inputs: the stored stage-3 E6 fit table (`analyze3.load_table("E6")`, the loader of S1a) · reference tissues = `config.refs_e5()` (210) · relaxation settings = the 17 of S1a (`analyze3.lsets()["L1|cv0.03|rho0.0|iii"][0]`) · three generators · pair classes Q1–Q4 (reference acquisition A against Q). 204 cells.
- Coordinates: x = logit ν, y = logit f_iso − ln ν (the definition of `lib.coords`). Fitted values are the fiber-direction mean coordinates of the stored table (as in S1a). Raw NDI (ν_fit) is the fiber-direction mean of the fitted NDI.
- True values: the (ν, f_iso, ODI) of the reference tissues and their coordinates (x, y). Explanatory variables of (a) = (x_true, y_true, ODI_true) plus an intercept; of (b) = (ν, f_iso, ODI) plus an intercept. Least squares (`numpy.linalg.lstsq`).
- Δ_max = largest minus smallest fitted value of the regression (intercept included) over the 210 tissues. Width = 1.96·√2·s.
- s: for x, s_x_median of S1a (0.031, the same computation as the stored s1a.json); for y, the same computation (median over units of the ON-Harmony OXF1PRI retest SD) on scale y; for raw NDI, on scale raw_ndi (`empirical.s_onh`).
- Signal-model part of (b): the apparent ν′ of Eq. (3) from the relaxation setting of the cell and the TE and TR of the two acquisitions; d_raw,model = ν′_Q − ν′_A regressed on the same explanatory variables.
- (c): at each acquisition, x_fit,a regressed on x_true (α_a, β_a, largest |residual|); β_Q/β_A for each pair. Descriptive.
- The reading (definition §3.4) is applied to x (the fixed sentence concerns logit NDI). The ratios of y go in the same table; if any pair class of y exceeds 1, the supplement says so. R′ for raw NDI = the largest ratio over all cells.
- Failure controls: (i) with d = the algebraic shift of Proposition 1 (the constant ln(w_in/w_en) of the cell), all three coefficients other than the intercept |·| < 10⁻¹⁰. (ii) With 0.10·(x_true − mean) added to that d, the x_true coefficient 0.100 ± 10⁻⁶ and Δ_max ÷ width > 1. Both through the same code path as x.

## §4 Negative control

- Input: the Tong roi.json («주» (main) table) — the unit differences and thresholds of `_tong_scales` in rules_env (t₀.₉₇₅(ν)·√2·s, s = the unit retest SD from c10r1–r3) · scales x, y, odi, raw_ndi, logit_fiso. The «share» is the value without a common shift subtracted (the same as the fixed T-B: x 0.931, ODI 0.934 and y 0.844 are reproduced first, and the run stops if they differ).
- By participant: units = (center pair, ROI) of that participant. By center: the nine pairs that include the center × three participants × ROI.
- Equivalence bound: T₆(λ) is a noncentral t (scipy.stats.nct, 6 degrees of freedom). Expected share = P(|T| ≤ t₀.₉₇₅(6)). λ is solved for the share with brentq (0 ≤ λ ≤ 10). δ/s = √2·λ·(σ/s) with σ = s (the δ/s = √2·λ of the definition's table). For the three observed shares (x, ODI, y) and for 0.90.
- Conversion: δ (logit) = (δ/s) × the median over units of the c10 retest SD of the scale. For x, in NDI = δ·ν(1−ν) at ν = 0.5 and at the median Tong ROI NDI (median of the raw NDI over all sessions and units).
- Simulated failure control: 3 participants × 45 pairs × 27 ROIs = 3,645 units. For each ROI, s² ~ σ²·χ²₆/6 drawn once and shared by its 135 units. Differences ~ N(δ, 2σ²), σ = 1. 1,000 replicates, seed 20261116. Criteria: δ = 0 → mean 0.95 ± 0.005 · δ = 1.086 → 0.90 ± 0.005 · δ = 2 → < 0.80.
- Diagnosis of y: mean FISO of a unit (ROI) = mean of the raw FISO over all Tong sessions (the r1 center columns and the three c10 repeats). Share at retest level of a unit = share over all (participant, pair) units of that ROI. Kendall τ_b (scipy.stats.kendalltau, variant b) · shares of the groups of units below and above the median (all unit differences of a group pooled). Units equal to the median go «above». The same for x and logit_fiso.

## §2 Registration of the multi-TE series (after the approval of the reads)

- Series: reference = TE62 (the first 62 ms); the other nine = TE72 … TE132, TE62R2, TE132R2. Participants sub-01, 02 and 03.
- 2.1: r recomputed with the functions of the fixed check (`mean_dw_ratio` and `align_r` of `run.py`, the brain mask of the reference series). Low-b0 voxels = the lowest 1 % within the mask of the reference mean b0 (np.percentile, set from the reference only). r also without those voxels in the same mask. The source data (provider dwi) must be read.
- 2.2: inputs = the mean b0 made during processing (reference `wong/<sub>/b0_mean.nii.gz`, series `own_<series>/b0_mean.nii.gz`) and the reference brain mask. `flirt -in series -ref reference -dof 6 -cost normmi -searchcost normmi -omat` (FSL defaults otherwise). Displacement = largest |M·p − p| over the mask voxel centers p (FSL scaled mm coordinates). Voxel size = the smallest of the three pixdim of the series (stricter for "aligned"). The 2-voxel shift of control (ii) = an integer voxel shift (zero at the edge), x axis = the first image axis, the displacement divided by the pixdim of that axis.
- Maps of paths B and C: the seven native-space F2 maps and the two fit-validity maps made for each series during processing (`own_<series>/f2_<series>_*.nii.gz`). Registration to standard space = `fa2std_warp` of the reference series, brain mask = the reference mask (nn), otherwise the commands and ROI function of `to_std_and_roi`. Path C adds `--premat=<the series-to-reference matrix of 2.2>` to the same applywarp. The reference series itself is also redone with the new code and checked against the stored `roi_TE62.json` first (stop if different).
- Recomputation: for each path, a separate working folder shaped like `wong/<sub>/` (DONE.json with «따로_정합» (separately registered) = empty list; path B with the ROI tables of passing series only), and the existing code (judge_wong.judge · the ODI range ratio of S10 · the slopes of S-M8.3 · rules_env e1, e1_fix, e4, e5, e6) run with REMEASURE_WORK set to that folder. Only the registration step differs (reported).
- β of §5.4: the β by participant of E1 and its 95 % bootstrap interval (the function of e1). Condition of the first two rows = β > 0 with an interval excluding 0 in all three participants. Third row = the interval includes 0 in at least one participant and no participant has β < 0 with an interval excluding 0. Any other pattern: no sentence is chosen, and it is reported.

## §5 Noise floor of the multi-TE data (SNR after the approval of the reads)

- Tissue states: the ROI-matched units of S1b (participant × ROI, the observed values at TE62 inverted with Eq. (3) at the central relaxation) · central relaxation · generator T-in · TR 5.8 s · direction scheme = the same cards as S1b (WONG:TE<n>).
- 200 voxels: 50 for each of the four fiber directions (FDIR_E5). Signal = the noiseless formula of synthfit × 1000 (b0 = 1000).
- Noise: noncentral χ, L ∈ {1, 4, 16}. Signal S/√L in each channel, real and imaginary noise N(0, σ²), magnitude = √(Σ over channels |·|²). σ = 1000 / SNR(participant, TE). Seed = 20261116 + participant number × 1000 + L.
- SNR: `analysis.snr` (S-M6) by participant and TE. sub-01 from the stored `b1s_snr.json`. sub-02 and 03, if not recorded, need the b0 volumes of the source data (within the approval). Without the approval, the TE values of sub-01 are used for all three participants, and this is reported.
- Fit: synthfit.fit_chunk (AMICO F2 unchanged). Unit value = mean of the fitted values (NDI, FWF, ODI) of the 200 voxels → coordinates.
- Quantities (computed in each of the 50 Monte Carlo replicates; mean and range over replicates): TE slopes of x and y by unit (OLS), τ_sim = √max(0, Var over units (replicate-mean slope) − mean over units (variance between replicates)/50) and the value without that subtraction, β_sim = `beta_level` as in E1 (unit mean levels), the share by the fixed linear rule (judge_wong), the ODI range ratio (denominator = median over units of the difference of two replicates at TE62 and TE132), trajectory direction = median over units of atan2(k_y, k_x).
- Controls: (i) at SNR = ∞ (no noise, one replicate), τ_sim within ±10 % of the fitter-only value of S1b (T-in, central relaxation, x). (ii) At SNR = 5 and L = 1, some unit has a mean NDI or FISO moved from the noiseless value by more than the retest width (1.96·√2·s, s = the repeat SD of the multi-TE data).
- τ_obs = the smallest of the corrected pooled τ_x of E1 (e1_fix) of paths A, B and C.
