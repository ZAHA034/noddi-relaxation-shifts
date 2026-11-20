# Post hoc analyses defined on 2026-11-16 (English translation)

> English translation of `decomp/plans/S73_POST_HOC_PLAN_2026-11-16.md`, the definition committed before any of these quantities was computed (development commit 93c45ce, 2026-11-16 04:55:28 KST). The Korean file is the record, and where the two differ the Korean file governs. Working labels (review session, A, relay numbers, [author decision]) are explained in `docs/identifiers.md`. What was computed, and what stopped, is reported in Supplementary Methods S-M10.6. Some process notes are shortened in this translation.

# Direction S73: making the manuscript firmer — sentences to fix and four analyses fixed before computing

From: the review session (in-project reading, not external review). Passed to A after the author confirmed it.
Subject: the submission build of v16k (main text 33 pages, sha256 `1c6456bf…`) and the supplement (90 pages, `028b4cc2…`).
Marks: [author decision] / [reviewer proposal] / [author decision needed].
It does not overlap the formal items of S72; both proceed together.

## 0. Purpose and principles

- **Purpose:** [author decision] raise the rigor of the manuscript and its chance of final acceptance by carrying out §1 (sentences) and §2–§5 (four analyses) below.
- **Principles:**
  1. **Definitions first.** The definitions, readings and per-outcome sentences of §2–§5 are committed to the repository as written here before computing; the commit hash and time are recorded. All are post hoc (defined after the verdict, fixed before computing).
  2. **Results are reported as they come, even if unfavorable.** The fixed verdicts (Table 1, Supplementary Table S5) do not change; new results are placed beside them.
  3. **Failure controls.** For each check, a deliberately wrong input is entered first to show that the check catches it. If a control does not behave as expected, that analysis stops and is reported.
  4. **CPU only.** No new data are downloaded; only files already processed are used. If data would have to be downloaded or newly read, stop and ask the author.
- **Order:** §1 proceeds at once. The analyses run in the order §3 → §4 → §2 → §5; the reading of §5 uses the results of §2.

## 1. Sentences to fix (no computation) [reviewer proposal]

### 1a. An overstatement in Section 5.1

- Current sentence:
  > Between regions of one brain, the condition of shared relaxation failed (exploratory).
- Problem: Section 4.3 says the candidate causes (region-specific compartmental T2, a noise floor) were not separated. What was observed is only that one common TE shift did not hold between regions. As §2 shows, the separate registration of each series is a candidate too.
- Revised wording (provisional; the final wording follows the table of §5.4):
  > Between regions of one brain, one common TE shift did not hold (exploratory).
- In the same paragraph, the Section 4.3 sentence «Both results concern a shift shared by regions with different relaxation» becomes «…shared by regions that may differ in relaxation».
- The scale-criterion sentences of Sections 4.3 and 5.4 are conditional («regions with different relaxation») and stay.
- List every expression in the main text, the supplement, the README and the cover letter that asserts a relaxation difference between regions (for example «shared relaxation failed», «relaxation differs between regions»). Fix them only after the list has been seen.

### 1b. The registration sentence of Section 4.3

- Current sentence:
  > Every series other than the reference failed the alignment check and was registered on its own FA, so the analysis without separately registered series could not be computed
- Problem: readers take it to mean the data were misaligned. According to Supplementary Table S10 (d), however, the repeat at the same TE was below the criterion too (r 0.003–0.21): the check did not separate aligned from misaligned series.
- Revised wording:
  > The alignment check fixed before computing did not separate aligned from misaligned series: every series other than the reference fell below it, including the repeat at the reference TE, so each was registered on its own FA and the analysis without separately registered series could not be computed under that check (Supplementary Table S10 (d)).
- The result sentence of §2 follows this sentence directly.

### 1c. Start of Section 4: the tolerance of the cone measure, and direction against slope

- At the end of Section 3.4 and at the start of Section 4, «pointed within the algebraic cone» becomes «pointed within ±1° of the algebraic cone» (definition in Supplementary Methods S-M3).
- Directly after «…which it did not measure.» insert:
  > The measure counts a departure common to all tissue states, which changes intercepts only, as well as one that differs between them, which can change slopes; the measures below concern the second.
- Add the proportional bias to the parenthesis of the exploratory sentence, but **A first reports its definition**: what is regressed on what, over which set of tissues and in which units (Supplementary Table S17 (a) gives no definition).
  - If it is a quantity that measures how far the shift varies with tissue state:
    > (share 1.000; proportional bias −0.023 to +0.035; lowest rank correlation 0.997, no rank reversals)
  - Otherwise leave it out, and let the result of §3 stand in its place.
- Add a one-line definition of the proportional bias to Supplementary Table S17 (a).

### 1d. Two theory sentences

- **At the end of the «Equivalence» paragraph of Section 2.2** (after «…assessed by simulation (Sections 3.4 and 4).») insert:
  > Models without a tortuosity constraint, such as the Standard Model with free compartmental diffusivities, meet the condition of Proposition 1, so for them the translation is exact wherever the fit is unique, which degeneracies limit (Novikov et al., 2019).
  - Basis: the review session confirmed in the PMC full text that Novikov et al. use the name "Standard Model", define it as a sum of non-exchanging Gaussian compartments and discuss its degeneracies (PMC6481929).
  - New reference: Novikov, D.S., Fieremans, E., Jespersen, S.N., Kiselev, V.G., 2019. Quantifying brain microstructure with diffusion MRI: Theory and parameter estimation. NMR in Biomedicine 32, e3998. doi:10.1002/nbm.3998
  - A checks it with the bibliographic tools. PubMed gives the online date 2018-10-15 and the volume and issue 32(4).
- **In Section 6 (Conclusions), after the second sentence** insert:
  > At a single acquisition, a change that scales both tissue weights by one factor, such as a lengthening of T1 shared by both tissue compartments, is exactly equivalent to a change in free water (Proposition 2).
  - This restates the λ statement of Proposition 2: a proved result that the conclusions do not yet contain.
- [author decision needed] Whether the same sentence also goes into the abstract is the author's decision; the abstract would then have to lose as many words.

## 2. Analysis ①: reanalysis of the multi-TE data with a working registration check [reviewer proposal]

**Why:** all numbers of Section 4.3 and some of Section 5.1 and Box 1 come from the path in which each series was registered separately: the linear shares 0.815/0.901, 25.4°, ODI 6–9 times, τ, β and ΔTE* 17 ms, and E4–E6. The path originally fixed applied one registration of the reference FA to the series as aligned by the provider. That path was blocked by the check, not by the data.

### 2.1 Cause (descriptive, no reading)

- Recompute the metric of the fixed check (voxel correlation of the mean signal normalized by b0).
- Report how the correlation changes when voxels with low b0 (for example the lowest 1 % of b0 within the mask) are left out.
- This is a fact check for the change record only; it is not used to set a new criterion.

### 2.2 New registration check (fixed before computing)

- **Objects:** each series of the provider-aligned data (for each participant, the nine series other than the reference) in the native space of the reference series (the first 62 ms).
- **Metric:** register the mean b0 of each series to the mean b0 of the reference with FLIRT, 6 degrees of freedom (rigid), cost function normalized mutual information. Measure the largest displacement (mm) that the estimated transform produces within the brain mask.
- **Criterion:** "aligned" if the largest displacement is at most half the voxel size of the series.
- **Failure controls (run first):**
  - (i) the reference with itself: largest displacement at most 0.1 voxel;
  - (ii) a copy of the reference mean b0 shifted by exactly 2 voxels along x: estimated displacement 2 ± 0.25 voxels, classified "not aligned".
  - If either fails, the check is invalid: stop and report.
- **Reported only:** the result for the repeat at the same TE (62 ms). "Aligned" is expected, but it is not a condition of the reading.

### 2.3 Paths to compute (compute both and report both, whatever the results)

- **Path B (the path originally fixed):** use only the provider-aligned series that pass the new check. Apply the one FA-to-standard registration of the reference series to all series (applywarp, trilinear). The valid-mask rule (weight 0.99) and the ROI rule are those of S-M6. If too few TEs pass for the linear reading (the minimum number of TEs of the fixed rule), report "not computable".
- **Path C:** register every series to the reference with the rigid transform of 2.2, then apply the one registration of the reference series.
- **Existing path A** (each series registered separately): the values of the current manuscript, carried alongside for comparison.
- The repeat series (62 and 132 ms) go through the same paths, and the retest SD is recomputed for each path.
- Only the registration step differs between paths; report the code diff.

### 2.4 Quantities to recompute (paths A, B and C side by side)

- Linear shares of x and y: fixed rule (S-M7.4), overall and by participant.
- ODI range across TE ÷ median repeat difference: by participant (Supplementary Table S10).
- TE slopes, slope ratio against the literature prediction, median direction of the unit trajectories (S-M8.3).
- τ, β and ΔTE* with their intervals (S-M10 E1); E4, E5 and E6 (S-M10.5).
- A figure like Fig. 5 for paths B and C (supplementary figure).

### 2.5 Per-outcome sentences (fixed before computing)

- **ODI:** for each participant, compute "ODI range ratio of path B ÷ that of path A".
  - If at most 0.5 in all three participants:
    > With one registration for all series, the range of ODI across TE was X–Y times the median retest difference, against 6–9 with separately registered series (post hoc), so registration accounted for most of it.
  - If at least 0.8 in all three, the same sentence ends «…so separate registration did not account for it.»
  - Otherwise it ends «…accounted for part of it.»
  - If path B cannot be computed, apply the same rule to path C and name the path.
- **Linear shares:** the fixed results (0.815, 0.901) stay. Section 4.3 adds the values of path B (or C) as «(post hoc, one registration for all series: …)».
- **ΔTE\*:** if the value of path B (or C) lies outside the 90 % interval of path A (15.1–20.7 ms), Section 4.3 and Box 1 give both values («about 17 ms with separate registration and about Z ms with one registration»). If inside, the main text stays and the value goes in the supplement only.
- **E4–E6 and 25.4°:** side by side in the supplement; the main-text numbers do not change. If the sign of the E4 difference changes in path B (or C), Section 4.3 says so in one clause.

## 3. Analysis ②: does the shift at fitter output depend on tissue state? (analysis 1) [reviewer proposal]

**Why:** metric invariance means "the shift does not depend on tissue state". The main text now has only the direction (cone), the size after one common shift, and ranks. This answers directly the first question a referee will ask.

### 3.1 Inputs

- The stored stage-3 fits (S1a): 210 reference tissues, 17 relaxation settings, three generators, four pair classes (TE only, TR only, direction scheme only, all different).
- Record the sha256 of the input files.

### 3.2 Computation

One cell is (relaxation setting, generator, pair). The tissues in a cell share relaxation.

- **(a) Tissue dependence of the shift:**
  - Regress dᵢ = x_fit,B,i − x_fit,A,i on the true (x, y, ODI) by least squares.
  - Report the coefficients, R² and Δ_max (largest minus smallest shift predicted by the regression within the range of tissues).
  - Quantity read: Δ_max ÷ (1.96·√2·s), s = 0.031 (the median OXF1PRI retest SD, as in S1a).
  - Do the same for the y coordinate, with s the median OXF1PRI retest SD of y computed in S-M10. The claim of a parallel shift concerns both x and y.
- **(b) Contrast with raw NDI:**
  - Regress d_raw = ν_fit,B − ν_fit,A on the true (ν, FISO, ODI).
  - Quantity read: Δ_max ÷ (1.96·√2·s_raw), s_raw the raw-NDI retest SD of S-M10.
  - Report also d_raw computed in the same cell from the signal model alone (Eq. 3, without the fitter), to separate the part predicted by Proposition 3 from the part the fitter adds.
- **(c) Slopes by acquisition:**
  - At each acquisition a, regress x_fit,a on the true x to obtain α_a, β_a and the largest residual; report β_B/β_A for each pair.
  - The reference is β_B/β_A = 1, not β = 1: within one acquisition AMICO may give β ≠ 1.
  - Descriptive, no reading.

### 3.3 Failure controls (run first)

- (i) Enter the algebraic values of Proposition 1 (an exactly common shift) instead of the fitter: all coefficients |·| < 10⁻¹⁰.
- (ii) Add 0.10·(x_true − mean) to the d of (i): the x coefficient 0.100 ± 10⁻⁶, and Δ_max ÷ width > 1.

### 3.4 Reading and per-outcome sentences (fixed before computing)

- **Reading:** if Δ_max ÷ width is at most 1 in every cell of a pair class, "the tissue-dependent shift is within retest level". Report for each pair class the share of cells at most 1 and the largest ratio.
- **If every pair class is within**, insert at the start of Section 4:
  > At fitter output, the displacement of logit NDI depended on tissue state by at most R of the retest width in any pair class, and that of raw NDI by up to R′ (post hoc; Supplementary Methods S-M10), so differences between tissue states in logit NDI changed by less than retest level.
- **If any pair class exceeds**, write «…by up to R of the retest width when [pair class]…», and state the same in the evidence cell of the «Associations and ranks» row of Table 2.
- Put the logit against raw numbers in one clause in the evidence cell of the same row of Table 2.
- **Limitation sentence (fixed):** «noise-free simulation; not a test of metric invariance in people».

## 4. Analysis ③: negative control — shares by participant, equivalence bound, diagnosis of the failure in y (analysis 2) [reviewer proposal]

### 4.1 Shares by participant and by center (descriptive, no reading)

- **By participant:** with the fixed rule (t₀.₉₇₅(6)·√2·s), the share at retest level of each participant's 1,215 units (27 ROIs × 45 center pairs), for x, y and ODI and for raw NDI and logit FISO of the fixed sensitivity set.
- **By center:** the same share over the units of the nine pairs that include the center (729).
- **Main text:** add to Section 4.1 the clause «per participant, x a–b, ODI c–d, y e–f (post hoc; Supplementary Table S14)».

### 4.2 Equivalence bound (analytic, post hoc)

- **Model:** every unit carries the same center effect δ, the noise of a difference is N(δ, 2σ²), and s is estimated with 6 degrees of freedom.
  - The expected share is P(|T₆(λ)| ≤ t₀.₉₇₅(6)), with λ = δ/(√2σ).
  - Find the δ at which the expected share is 0.90, and the δ that corresponds to each observed share.
- **Values computed by the review session** (A recomputes them with the same formula and confirms agreement):

  | Share | λ | δ/s |
  |---|---|---|
  | 0.90 (threshold) | 0.768 | 1.086 |
  | 0.931 (x) | 0.476 | 0.673 |
  | 0.934 (ODI) | 0.437 | 0.618 |
  | 0.844 (y) | 1.110 | 1.570 |

  - Without an effect the expected share is 0.950.
- **Conversion:** to logit units with the median c10 retest SD of each coordinate; for x also to NDI units at NDI 0.5 and at the median ROI NDI (dν = ν(1−ν)dx).
- **Failure control (simulation, 3,645 units, 1,000 replicates):**
  - at δ = 0 the mean share is 0.95 ± 0.005;
  - at δ = 1.086 s, 0.90 ± 0.005;
  - at δ = 2 s, below 0.80 (analytic 0.777).
  - In the simulation, s is drawn once per ROI (6 degrees of freedom) and shared by that ROI's 135 units, as in the real data.
- **Main text (end of «Basis of the thresholds» in Section 3.4):**
  > Under a center effect common to all units, the expected share falls to 0.90 at about 1.1 retest SDs, and the observed shares correspond to about 0.7 (x) and 0.6 (ODI) retest SDs, or about Q in NDI (post hoc; Supplementary Methods S-M7.3); effects that differ between units are not covered by this bound.
  - The numbers are those recomputed by A.

### 4.3 Diagnosis of the failure in y (descriptive)

- For each ROI unit, pair the mean FISO of the Tong data with the share at retest level in y. Report:
  - Kendall τ_b;
  - the shares of the units with mean FISO below the median and above it;
  - the same for x and logit FISO (for comparison).
- **Per-outcome sentence (fixed):** only if the "below" share in y is lower than the "above" share by at least 0.05, insert in Section 4.1:
  > The units below retest level in y were concentrated where FISO was lowest, where the logit enlarges small absolute differences (post hoc; Supplementary Table S14).
  - Otherwise the supplement only.

## 5. Analysis ④: simulation of the noise floor of the multi-TE data [reviewer proposal]

**Why:** Section 5.3 says "SNR falls with TE; a noise floor may have contributed to the trajectories but was not separated". The τ produced by the fitter alone already exists (0.00035/ms, against 0.0029/ms observed). The same framework adds the noise floor.

### 5.1 Inputs and settings

- **Environment:** rules_env.
- **Tissue states:** the ROI-matched states of S1b (multi-TE; the observed ROI values at TE 62 ms mapped to true fractions with Eq. 3 at the central relaxation).
- **Acquisitions:** the eight TEs, TR 5.8 s, direction scheme under the same assumption as S-M10.1.
- **SNR:** the measured session SNR by participant and TE (method of S-M6).
  - If SNR by TE is not recorded, compute it with the method of S-M6 from b0 files already processed. Read no new data.
- **Noise:** noncentral χ, L ∈ {1, 4, 16}. The effective L is unknown, so all three are computed.
- **Scale:** for each unit, 200 voxels synthesized and fitted separately and averaged (as ROI means are). 50 Monte Carlo replicates.
- **Quantities:** TE slopes of x and y by unit, τ_sim (with the noise-variance correction of E1), β_sim, the share by the fixed linear rule, the ODI range ratio, the median direction of the trajectories.

### 5.2 Failure controls (run first)

- (i) SNR = ∞: τ reproduces the fitter-only value of S1b (0.00035/ms) within ±10 %.
- (ii) SNR = 5: the unit mean NDI or FISO moves by more than the retest width, which shows that this path can produce a noise-floor effect.

### 5.3 Reading (a convention fixed before computing)

- **Reference value:** τ_obs is the smallest observed τ among paths A, B and C of §2, the case most favorable to the noise floor.
- **Criterion:** use the largest τ_sim of the three L.
  - τ_sim/τ_obs ≤ 0.25: "the noise floor and fitter do not account for it";
  - 0.25–0.75: "account for part";
  - at least 0.75: "can account for most".
- **The same rule applies to:**
  - the ODI range ratio;
  - β: report its sign and size, and apply the same bounds to β_sim/β_obs.

### 5.4 Sentences for Sections 5.1 and 4.3 (fixed before computing)

| Result of §2 (β of path B or C) | Result of §5 (τ) | Sentence of Section 5.1 |
|---|---|---|
| β > 0 and the interval excludes 0 | does not account | Between regions of one brain, one common TE shift did not hold, and neither registration nor the noise floor and fitter accounted for this in reanalysis and simulation, leaving region-specific relaxation as the remaining candidate among those examined (exploratory; post hoc). |
| β > 0 and the interval excludes 0 | part or most | Between regions of one brain, one common TE shift did not hold (exploratory); in simulation the noise floor and fitter could account for [part/most] of the between-region variation (post hoc). |
| the interval includes 0 | (any) | With one registration for all series, the dependence of the TE shift on NDI was not distinguishable from zero (post hoc), so one common TE shift was not shown to fail between regions. |

- If the third row applies, the main-text sentences on E4–E6 and «about 17 ms» are rewritten beside the values of that path; the review session checks the wording after seeing the results.
- Whatever the row, «these analyses do not separate them» at the end of Section 4.3 is changed to fit the results.

## 6. Not done (outside this scope) [reviewer proposal]

- **Predicting the FISO signs again with adjusted relaxation values:** not done. There are only eight pairs and the results have been seen (unblinded), so neither outcome could carry weight.
- **To do if a referee asks:**
  - comparison with taking the logit in each voxel before averaging;
  - an analysis explaining the differences in TE slope between ROIs by the angle between fibers and B0;
  - analysis of the real data with a second fitter;
  - the precision of estimating the group difference at TE = 0 by regressing Eq. (7) on site TE.

## 7. Report requested

- The hash and time of the commit of this document as the definition (before any computation).
- §1: the sentences before and after, the list of 1a, the definition of the proportional bias of 1c, the bibliographic check of 1d.
- For each of §2–§5:
  - the failure controls (passed or stopped);
  - the sha256 of the inputs;
  - the result tables;
  - the row chosen in the per-outcome table and the sentence inserted;
  - the ledger rows (numbered by A).
- §2: the cause of 2.1, the results of the check of 2.2 (largest displacement by series), the table of paths A, B and C, the code diff.
- The change in main-text words and pages, the result of the H1 gate.
- New PDFs (reading and submission builds) with their sha256.

To be kept as before: write only "fixed before computing"; do not use restore, remove, harmonize or correct for our outputs; do not report point estimates of the judge as results; keep the marks [A proposal]/[reviewer proposal]/[author decision]. In repository text, write «author decision».
