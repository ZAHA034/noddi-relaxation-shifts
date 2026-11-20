# Dated plan of the analyses after the verdict (English translation)

> English translation of `decomp/plans/EXPLORATORY_PLAN_2026-11-01_redacted.md` (SHA-256 e83a8aadd3cf5bc5ac8e5e3a244e1b3acf59ca078db3fa238d031a037386c7c2), checked by script for numbers and identifiers. The Korean file is the record; where the two differ, the Korean file governs. Some process notes are shortened in this translation.

The deposited Korean file `decomp/plans/EXPLORATORY_PLAN_2026-11-01_redacted.md` differs from the SHA-256 given in this header only in relative path strings (status paths-neutralized in `MANIFEST.csv`); the redaction list and its checker named in the header are not part of the deposit.

<!-- Head of the redacted copy — only this block was newly written. Below this block the text is byte-identical to the original, and differs only at the places in the list below. -->
> This file is the redacted copy of the exploratory plan `EXPLORATORY_PLAN_2026-11-01.md` (working document outside the repository · sha256 83a427ef291645785f7047aa948ea48bc2cdaf31dc3798771f55ef51d14afcb0 · 248,434 bytes · read at 2026-11-04T22:48:11+09:00) (author decision 2026-11-05 #2).
> What is redacted: personal identifying expressions — the author's name, in 80 places, became «the author». Absolute paths · e-mail addresses · accounts · machine names do not occur in the original (0 places). The 13 places where the same two syllables occur inside ordinary verbs (for example the verbs meaning "to set" and "to fix") stay as they are.
> List of the redacted places: `EXPLORATORY_PLAN_2026-11-01_redactions.json` · check that rebuilds this copy from the original: `check_redaction.py`.
<!-- End of the head of the redacted copy -->
# Exploratory analysis plan — why the S4 judge failed on the measurement (ON-Harmony V4)

- **Status:** basic plan approved (conditional · relay 2026-11-01 22:05:24 KST) · fixes 2-1 – 2-7 applied (relay 2026-11-01 22:33:57 KST) · additions 1-1 – 1-3 · decisions applied (relay 2026-11-01 · this reply) · review session verdicts (11-01 j reply) 4 · 5 · 6 applied (2026-11-01 23:22:34 KST) · the D2 · D3 truth tests start alongside C7 by the author's chat approval («proceed with no. 1» · 2026-11-01, in the 23 o'clock hour) — D1 after the C7 commit · optional items on hold · s1 – s4 not approved · computation 0. At each gate the results are reported, and the next gate opens after approval (§10).
- **Basis:** relay «On the exploratory analysis plan» (2026-11-01 21:18:26 KST) · the author decision is written as PLAN §9-1 #80 (non-reading batch 68 · computation 0 before that commit).
- **Written by:** A.
- **What stays as it is:** the fixed analysis document v10, fixed before computing (`decomp/remeasure/MEASUREMENT_FIXED.md` · 5f95b18) · the verdict results (5632626 · `decomp/results/remeasure/judgements.json`) · the reading of the verdict · the frozen judge (weights 85eb8b0).

## 0. Author decision · prohibition · what is kept (relay original text as it stands)

- Author decision: «The verdict and the fixed reading stay as they are; the reasons for the failure are analysed by post hoc exploration.»
- Prohibition: «Changing the simulator or the judge, evaluating again on the same ON-Harmony data, and writing "it is calibrated". If retraining appears necessary, the plan records only "needed · not done".»
- Kept: «Training runs on CPU only. Tong's Technical Validation and verdict/ are not opened. The judge's point estimates are not used as results.»

How this plan follows these:

- The head of every output states «post hoc exploratory — not a verdict». The verdict, coverage rates and thresholds are neither recomputed nor changed.
- «Point estimates are not used in results» refers to **the judge's posterior mean** (2-5). The posterior mean is used only to compute the z of D1; the values themselves (the size of the scanner difference) are not given in tables · figures · the manuscript · talks. The per-ROI-unit arrays stay only in the value area, for reproduction.
- The r of D3 is an estimate from classical fitting, not from the judge. It is reported together with the profile interval as an exploratory result (2-5).
- Training is 0. The frozen judge and frozen simulator are called as they are. It runs with `CUDA_VISIBLE_DEVICES=`, and if CUDA is visible the code stops (the same gate as `judge_v4.real_fns`).

## 1. What is known now (new computation 0)

**Verdict results** (`onharmony/v4/v4_c3.json` · main cells 29 · cells with data 24):

- Coverage rates are 0.000 – 0.649 and the per-cell thresholds are 0.50 – 0.85.
- Gap detection occurs in 23 of 24 cells, and the misspecification excess rate is 0.16 – 0.77.
- Among the sensitivity cells, 4 cells passed at SNR 50 · N_eff 25 and 1 cell at SNR 50 · N_eff 100. In the two cells of SNR 20, passes are 0.
- The verdict is «failure» 24 · «no data» 5.

**What the simulator and the judge assumed** (from the code):

| Item | Assumption | Location |
|---|---|---|
| Noise | non-central chi · coil 32 → L = 4 · 64 → L = 8 · σ = 1/SNR on the b0-normalised signal · independent per voxel and per volume | `decomp/s4/s4cfg.py:51` · `decomp/s4/simulate.py:61` |
| ROI mean | for each site, 256 voxels drawn separately around the same centre → F2 fit → mean m · voxel covariance Σ̂ → x = m + N(0, Σ̂(1/N_eff − 1/256)) | `simulate.py:4` · `:148` |
| Training range | SNR log-U[15, 60] · N_eff log-U[16, 200] · TE 60 – 110 ms · TR 3 – 11 s | `s4cfg.py:47` · `:52–53` |
| Relaxation prior | T2_in U(70, 95) ms · T2_in/T2_en U(1.0, 1.8) · T2_iso 1000 ms · common T1 850 – 1100 ms · T1_in/T1_en ∈ {1.3, 1/1.3, 1.0} (the three scenarios of k) | `s4cfg.py:29–31` |
| b value | magnitude ±5 % (table b0.95 · b1.05) only | `simulate.py:135–137` |
| V4 conditions | main cells SNR 30 · N_eff 50 · L0 · sensitivity cells SNR {20, 50} × N_eff {25, 100} | `decomp/remeasure/judge_v4.py:8–9` |

**What already exists from the measurement:**

- Session SNR (`snr_table.json`): the definition is the same as in S4 (`decomp/remeasure/analysis.py:109–112`).
  - The 12 OXF1PRI repeat sessions (sub-13192 · sub-16974 × ses-OXF1PRI001 – 006) are at 15.4 – 19.6. All are inside the training range but lower than the 30 of the main cells.
  - Of the 109 first sessions, 73 are outside 15 – 60. **All 73 are below 15 (4.45 – 14.99) and 0 are above 60.** The 36 inside the range are at 15.02 – 20.31. So the measurement SNR of all 109 first sessions is at most 20.31 (2-7 · `x14/snr_range_counts.json` · re-read record).

## 2. D1 — bias and overconfidence (new image reads 0 · exploratory analysis)

**Checked first — it is not stored.**

- One cell of `v4_c3.json` has 28 fields (`coverage` · `coverage_by_k` · `boot95` · `by_coil_combo` · `decomposition` · `sensitivity_passed` …). There are no per-ROI-unit posterior means · SDs · samples.
- Reason: `judge_v4.cover_by_k` (`judge_v4.py:155–186`) computes only whether the 90 % interval from the posterior samples T covers 0, and discards the samples.
- So they must be drawn again.

**Start · location:** it starts after C7 has finished and been committed separately. Output goes in the new directory `data/remeasure/explore/d1/`. `v4_c3.json` is not overwritten (that path is blocked before writing).

**1. Drawing again.** With the frozen judge (`judge_v4.real_fns("c3")` · weight comparison `check_frozen`) the posterior samples of the main cells (SNR 30 · N_eff 50 · C3) are drawn again.

- The same code path `collect` → `unit_results` is used.
- The seeds are also the same. The per-pair torch seeds are `seeds_torch` of `v4_c3.json`, and the numpy seed is 20261027.
- **The settings of the original run are used as they are (2-4):** the runner `drafts/scripts/x12_judge_resume2.sh:17` gives `CUDA_VISIBLE_DEVICES=` · `S4_THREADS=16`, and `decomp/s4/validate.py:23` calls `torch.set_num_threads(16)`. Deterministic settings (`torch.use_deterministic_algorithms` etc.) are not called, because the original run did not call them. OMP · MKL variables are not set, because the original runner did not set them.
- The input is the files read at the verdict as they are: `roi.json` «주» (main) · `attr_table.json` (sha256 compared with the access record values · re-read record).

**2. Reproduction gate (no tolerance).** In all 29 cells (24 with data + 5 without data), the values drawn again must be **exactly** the same as `v4_c3.json` before going on.

- Compared: the coverage of the nine k (`coverage_by_k` · the number of covered ROI units, obtained by multiplying by the number of ROI units used) · `n_units_used` · `n_excluded_misspec` · `n_units_dropped_k` · `by_coil_combo`. For cells without data, whether «ROI units 0» is the same.
- On failure, only the number of mismatched ROI units is reported and the work stops. Because the stored values are cell summaries only, «number of mismatched ROI units» is counted as the sum of the absolute differences of covered ROI units per cell and per k · and the differences in the number of ROI units used · the misspecification count · the k-unmet count (which ROI units they are cannot be known). In that case the D1 values are not given.

**3. Stage 1 diagnosis — PIT (2-1).** For each ROI unit u · assumed k · component d, PIT = the fraction of posterior samples smaller than 0.

- For a ROI unit with several coil combinations, a PIT is computed per combination and the value farthest from 0.5 is used. This is on the same side (the least favourable side) as coverage being «covered only if all combinations cover». Under this rule the two-tail fraction is almost equal to 1 − coverage (not exactly equal, because of quantile interpolation).
- In each cell, the distribution of PIT is compared with the reference distribution (4).
- Summary: **mean of PIT** (bias · reference near 0.5) · **two-tail fraction** (PIT < 0.05 or > 0.95 · overconfidence · reference near 0.10). The median over the nine k is given, with the minimum · maximum alongside.
- ROI units · misspecification exclusions are the same as in the verdict (card pair × subject × ROI 27). As a sensitivity analysis, values that include the ROI units dropped for misspecification are given alongside.
- Intervals come from a subject-cluster bootstrap with 2,000 replicates (descriptive).

**4. Reference distribution.** On the test data L0 · of the same acquisition type, PIT and z are computed with the same procedure (nine k · same judge · relative to the true values).

- **The coil-combination rule is used in the same way (1-1).** In each cell the candidate combinations of the measurement ROI units (a scanner with a record: that coil · a scanner without: 32 · 64 — `judge_v4.coil_options`) are used as they are. The number of combinations is 1 · 2 · 4.
  - The reference uses only test cases of that cell's acquisition type whose true coil combination is among that cell's candidates. In the measurement too, the true coil is among the candidates.
  - The posterior is drawn for each candidate combination, and the PIT farthest from 0.5 is used.
  - By comparing with the PIT of the same cases computed with the true combination alone, **the share by which the two tails thicken from this rule alone** is shown within the reference (per cell).
  - The test data are drawn once each with the four combinations (32|32 · 32|64 · 64|32 · 64|64) and selected per cell — about 4 × 26 minutes.
  - **Proviso (4a · review session verdict):** the reference contains only cases where the true coil is among the candidates. Because this condition cannot be checked on real data, a PIT outside the reference cannot by itself separate overconfidence from a coil outside the candidates.
    - Possibility of checking: the judge's coil record is only OXF1PRI 32 · OXF2PRI 64 (source LOCK_B EN :2904–2905 · `decomp/s4/plan.py` ONH · the rest «unknown»). If the ON-Harmony sidecars have a receive-coil field (e.g. ReceiveCoilName · ReceiveCoilActiveElements), it can be checked. However, the repository record shows no trace of that field having been read, and whether it exists is known only by opening the sidecars. They have not been opened yet (metadata · subject to the access record · approval needed).

- The existing `test_c3.npz` holds only coverage, so the test data are also drawn again.
- If the posterior is not normal, the values of «the well-calibrated case» can differ from the theoretical values (0.5 · 0.10 · SD 1). So the comparison is with this reference in place of the theoretical values.

**5. Stage 2 — z.** z = (posterior mean − 0) / posterior SD. Summary: mean of z (bias) · SD (overconfidence). The posterior mean is used only in this computation.

**6. Variance decomposition (2-2).** In each cell, the spread of z and PIT is split in two (computed per k, then the median over the nine k).

- (a) **Between subjects within the same ROI** = the pooled SD of what is left after subtracting each ROI's mean from the subject values of that ROI → the noise side.
- (b) **Between ROI means** = the SD of the ROI means → the side of systematic bias that differs across space. As a comparison value, «the SD expected when there is no difference between ROIs» = (a) ÷ √(mean number of subjects per ROI) is written alongside.
- The reference (test data) has no ROI · subject structure, so the reference is a single overall spread. (a) is compared with it.
- **«Overconfidence ≫ reference» in the reading frame is judged only by (a).**

**7. Quantitative prediction of the noise hypothesis (2-3 · written before computing · not a verdict).**

- If an under-assumption of noise is the whole cause, the cell's (a) spread of z ÷ reference ≈ 30 ÷ SNR_pair.
- SNR_pair = √(2 ÷ (1/SNR₁² + 1/SNR₂²)) — set from the measurement SNR of the two sessions of that ROI unit (`snr_table.json` `snr`). If the two SNRs are equal it equals that value.
- For each cell, «observed ratio = (a) ÷ reference» and «predicted ratio = the median of 30 ÷ SNR_pair over that cell's ROI units» are written side by side.
- A value summarising only the ROI units whose session SNR is outside the training range (15 – 60) (at least one of the two sessions outside) is also given.
- **Written down before computing (1-3):**
  - All 109 first sessions have SNR ≤ 20.31. 73 are below 15 (4.45 – 14.99) and 0 are above 60.
  - The predicted ratio 30 ÷ SNR_pair is 1.50 – 5.26 at the subject level and 1.79 – 3.33 as cell medians (24 cells with data · main-cell excluded subjects left out). Per-cell minimum · median · maximum and the number of subjects outside the range: `x14/snr_pair_prediction.md` (json sha256 f71f46b6ec9aef6f…).
  - These values come only from the stored attribute table · SNR table (judge 0 · re-read record). The ROI units that will drop out for misspecification are not yet known, so all subjects of the cell were used.
- The prediction applies only to z. PIT is bounded to 0 – 1 and so does not scale proportionally.

**8. Descriptive account of relations (not a verdict).** Per-cell bias (PIT mean · z mean) and the (a) spread are given in tables and scatter plots against the five below. Rank correlations are written as descriptive only; p values and thresholds are not used.

- |ΔTE|
- ΔTR (sign A → B)
- Direction design same/different (`cards.json` `design_differs`)
- Coil combination (32 · 64 · unknown)
- **Contrast differing only in coil (4b · [review session proposal · added after seeing Fig. 4]):** OXF1PRI~X and OXF2PRI~X, X = NOT1ACH · NOT2ING-98 · OXF3TRI · NOT4GEP · OXF4GEP · NOT2ING-89 (six pairs).
  - Written first from stored outputs (`x14/coil_contrast.md` · before computing): subject overlap · n_units · misspecification exclusions · k exclusions · ROI units used · coverage median and the range over nine k · session SNR median on the OXF side.
    - In all six pairs, the subjects of the two cells are the same people. The numbers of ROI units used read from the figure (161/64 · 197/49 · 162/142 · 104/47 · 87/40 · 47/10) are the same as in v4_c3.json.
  - In D1, the PIT · z summaries and (a) · (b) of the two cells are described side by side for each pair. No direction is set in advance, and no test is done.
- **ΔSNR (1-2):** the difference (B − A · card direction) and the ratio (B/A) of the measurement SNR of the pair's two sessions. It is a session-level value, so the 27 ROI units of one subject have the same value.
  - For each cell, whether the PIT mean and the sign of z move with ΔSNR is described. Both between cells (cell medians) and within a cell (between subjects) are looked at.
  - The direction is not set in advance. This is not a verdict.
  - The per-cell range of ΔSNR was written before computing in `x14/snr_pair_prediction.md`.

**What is reported:** per cell, only PIT summary · z summary · (a) · (b) · number of ROI units · coverage. The judge's per-ROI-unit posterior means and intervals in physical units are not given anywhere in tables · figures · the manuscript · talks. It is marked as exploratory analysis.

**Reading frame** (written down before computing · not a verdict):

- If PIT mean ≈ reference and the two tails · spread of (a) ≫ reference → the side of having assumed too little variation · noise (§4 candidates 1 – 5 · 11 · 12). Then it is read alongside the predicted ratio of 7.
- If the PIT mean is far from the reference and (a) ≈ reference → the side of systematic bias (candidates 6 – 9). If (b) is large, that bias differs by ROI.
- If both are large it is written as mixed.

**9. Recomputing the misspecification flags per sensitivity condition (approved · in the D1 window after the C7 commit · a few minutes).**

- For each of the four sensitivity cells (SNR {20, 50} × N_eff {25, 100}), only the misspecification score is recomputed with the frozen q density (no posterior samples · same threshold · same coil-combination rule).
- For each condition, «the number of ROI units dropped for misspecification» is given exactly.
- «The number of ROI units dropped for k not met» can be known only by drawing posterior samples. Only in cells whose stored pass flag is None (ROI units used 0) is it determined as «number of ROI units − misspecification count». Other cells are written as «unknown (samples needed)».
  - **This limit is accepted (decision).** The per-condition sample computation (about 19 hours) is not done. It is decided again after seeing the misspecification counts.
- The table made from stored outputs is `x14/v4_sensitivity_table_stored.md` (receipt h §4). The «passes 0» at SNR 20 is mostly cells where no ROI units remain (16 · 19 cells).
- PIT · z of the sensitivity cells (about 19 hours) are on hold.

**The C7 arm is not included this time.**

**Entanglement:** the 24 cells are pairs of 7 scanners. The cells are not independent of one another, and ΔTE · ΔTR · design · coil coincide with the scanner as one combination per pair. So the relations of 8 cannot separate causes.

**Outputs:** `data/remeasure/explore/d1/d1_c3.json` (cell summaries · reproduction comparison) · `d1_c3_units.npz` (per-ROI-unit arrays · not used in results) · `d1_sens_flags.json` (9).

## 3. D2 — spread of same-scanner repeats vs the judge's assumption (new image reads 0)

**Observation:**

- The `roi.json` «주» of the 12 sessions (the OXF1PRI repeats above) is used. 001 is the first session, and 002 – 006 are the 10 repeats for s_ref_F2 of R17. All are in the access record.
- For each ROI unit (27) · subject, the C3 coordinates of the 6 sessions are computed. The same `lib.coords` as the judge is used.
- The SD of the 15 session-pair differences is computed per coordinate. It has the same form as the «two-session difference» the judge sees.

**How the measurement SNR 15.4 – 19.6 was obtained:**

- **Image:** the b0 volumes (b ≤ 50 · 5 for OXF1PRI) of that session's eddy output 4D (`eddy_out.nii.gz`). This 4D was deleted after processing.
- **Formula:** W = voxels in brain mask ∩ dtifit FA ≥ 0.3. For each voxel, the mean S̄ and sample variance σ² (ddof 1) of the b0 volumes are computed. SNR = mean_W(S̄) / √(mean_W(σ²)).
- **Field:** `snr` (together with `n_b0` · `n_vox` · `in_train`) of `snr.json` in the session folder. These are collected as `sessions[<sub>_<ses>].snr` of `onharmony/v4/snr_table.json`.
- **Code:** `decomp/remeasure/run.py:234–243` (`session_snr`) → `decomp/remeasure/analysis.py:109–119` (`snr`) · `rconfig.SNR_FA_MIN = 0.3`.
- **Limits:** because it is a variance across a few b0 volumes, variation beyond thermal noise (physiology · residual motion) enters. Interpolation by eddy reduces the variance. It aims at the same quantity as the SNR of S4 (σ = 1/SNR of the b0-normalised signal) but is not exactly the same.

**The spread the judge assumed:**

- With the frozen simulator (`simulate.chunk` · `envs/reference` · AMICO · CPU), «same person · same scanner twice» is drawn. Settings:
  - Both sites with OXF1PRI acquisition: sidecar TE · TR · card design · coil 32
  - Δ = 0 · L0
  - Centre = the ROI value of that subject's first session (an approximation that takes ν · f_iso · ODI as true values — noted in the limits)
  - Relaxation = prior centre (T2_in 82.5 ms · ratio 1.4 · T1 975 ms · scenario iii)
- SNR takes four values: {30 (main cells), 20 · 50 (sensitivity), the measurement SNR of that session}.
- N_eff is {16, 25, 50, 100, 200}. N_eff changes only `roi_noise`, so no refitting is done.
- For each (ROI · subject · SNR), 500 cases are drawn.

**Comparison (descriptive):**

- For each coordinate · ROI, «observed SD / simulated SD (main cells)» is computed, and the median · quartiles across the 27 ROIs are written.
- The (SNR, N_eff) contour matching the observed spread is given alongside.

**Entanglement · limits:**

- From the spread of ROI means, SNR and N_eff cannot be set separately. The variance is roughly Σ̂(SNR)/N_eff, so only one contour results.
- The observed spread includes biological variation · registration differences · scanner drift together. So it is read as an upper limit, not as the size of noise alone.
- The approximation of setting the centre to fitted values can affect Σ̂. As a sensitivity check, the centre is shifted by ±1 quartile.

**Link with D1 (descriptive):** if observed ≫ simulated, this is material indicating that the judge's posterior is narrow even when the scanner is the same. It is written side by side with the overconfidence of D1.

**Outputs:** `data/remeasure/explore/d2_repeat.json` · simulation plan · chunks (`plan_d2.npz` · `sim_d2/`).

## 4. Candidate causes — what the simulator lacked · data to test with · whether they can be separated

| # | What the simulator lacked (or had differently) | Data to test with | Can it be separated? |
|---|---|---|---|
| 1 | Biological · physiological variation between sessions of the same person (V4 sets the truth to Δ = 0) | D2 (OXF1PRI repeats) | **No** — in the ROI table it is entangled with 2 · 3 |
| 2 | Registration · ROI mask differences and partial volume (voxels differing by session · mixing of CSF · grey matter at boundaries). The simulator only draws voxels separately around the same centre and has no systematic mixing | session difference in voxel count of `roi.json` vs D1 z · D2 difference (descriptive) | **Partly** — voxel count is an indirect indicator · voxel-level checking means reading the maps (optional s1) |
| 3 | Spatially correlated noise · interpolation (eddy · registration · partial Fourier) → actual N_eff · not white noise | D2 | **No** — entangled with SNR (contour only) |
| 4 | SNR level (main cells fixed at 30 · OXF1PRI 15 – 20 · 73 of 109 first sessions outside the training range) | `snr_table` · the existing sensitivity cells · D1 by SNR group | **Partly** — SNR differs by scanner, so it coincides with the scanner |
| 5 | Shape of the noise floor (effective number of channels · reconstruction method) — the simulator has L = 4 · 8 | none (background voxels · high-b signal = image reading) | **No** — entangled with SNR · not tested this time |
| 6 | Scanner-specific systematic effects (gradient nonlinearity · b error beyond ±5 % · residual eddy currents · reconstruction filter · Gibbs · EPI distortion) | D1 per-cell z described by scanner · pair · scanner concentration (already exists) | **No** — one scanner is one combination (would need a phantom · repeats on the same scanner with changed settings) |
| 7 | Difference in diffusion times (δ · Δ) — the compartment signal tables differ only by design (b values · directions) and have no time dependence | protocol records (if the sidecars have δ · Δ — metadata reading) | **No** — entangled with design · scanner |
| 8 | Compartment relaxation outside the prior (T2 ratio outside 1.0 – 1.8 · region dependence · myelin water · orientation-dependent relaxation · T1 ratio outside the three scenarios) | D3 (Wong joint fit) · D4 (predicted shift of the cards vs D1 bias) | Within the Wong scanner **partly** (only TE changes) · in ON-Harmony **no** (ΔTE coincides with the scanner) |
| 9 | Acquisition outside the training range (TE · TR — attribute ②) | `attr_table` ② · D1 by attribute group | **Partly** — coincides with the scanner |
| 10 | Water exchange between compartments · magnetisation transfer — the simulator is static 3 compartments | form of the D3 residuals (systematic misfit across TE) | **No** — the residuals alone cannot separate mechanisms |
| 11 | Heterogeneous structure within an ROI (crossing fibres · boundaries) — in the simulator, the spread of `pixels(rng, centre)` | voxel SD within the ROI on the maps (optional s1) | **Partly** — needs read permission |
| 12 | Motion · eddy outlier replacement | eddy text logs (`eddy_movement_rms` · `eddy_outlier_report` — not images · optional s4) vs D1 z (descriptive) | **Partly** — it differs by subject even within the same scanner, so it is less entangled with the scanner |

**Exploratory record — lines added to the candidate causes of failure (review session verdict k reply 7 · not manuscript):**

- «Wong direction 25° — relaxation outside the prior (evidence from another scanner)»: B1 (= E1 recomputation) θ median 25.4°. The TE direction of the S4 generator L0 prior (T2_in 70–95 · r 1.0–1.8 · T2_iso 1000 — `decomp/s4/s4cfg.py:29`) is 48.5–90°, and even at the L1 extreme (r_B ≈ 2.13 · the truncation is z·z ≤ 4 with the two components bound together — `s4cfg._trunc_bivariate`) θ ≈ 38.3° (fixed 2026-11-02 00:41:01 KST — the 2.29 · 34.7° first written were values from wrongly taking the truncation as ±2 per component). 25° corresponds to r ≈ 2.9 (T2_in 82).
- «Fornix FWF 0.47–0.49 > simulator upper limit 0.35» (D2 centre · `s4cfg.PIX_CLIP` FISO (0, 0.35)).

> «No» in this table means «cannot be separated within these data». It does not mean that the candidate is not a cause. When the cells are filled with D1 · D2 results, «cannot be separated» also stays as it is.

## 5. D3 — Wongkornchaovalit multi-TE joint fit (compartment T2 contrast)

**Why the single-TE slope is not used:**

- The verdict (R12) used the TE slope and linearity of F2 values (logit) fitted separately per TE. This stays as it is for the verdict.
- For that slope to equal 1/T2_in − 1/T2_en, the fit per TE (d∥ · d_iso fixed) must be right. Also this approach cannot separately handle the TE weighting of the free-water compartment and the extracellular diffusion tied by tortuosity.
- So the signals of several TE are fitted together with one model.

**Model** (direction-averaged signal → independent of ODI):

E(b, TE) = Σ_c w_c(TE) A_c(b) / Σ_c w_c(TE)

- Compartment weights:
  - w_in = (1 − f) ν e^(−TE/T2_in)
  - w_en = (1 − f)(1 − ν) e^(−TE/T2_en)
  - w_iso = f e^(−TE/T2_iso)
- Compartment signals:
  - A_in = spherical mean of the stick (d∥)
  - A_en = spherical mean of the cylinder (d⊥ = d∥(1 − ν) — generators T-in · T-app · T-none are sensitivity analyses)
  - A_iso = e^(−b d_iso)
- Fixed: d∥ 1.7 µm²/ms · d_iso 3.0 µm²/ms.
- **The free-water T2 is fixed — 1000 ms.** This is the value of the S4 simulator (`s4cfg.T2ISO`) and of the manuscript («free-water T2 1000 ms»). 2000 ms is given alongside as a sensitivity analysis.
- Free: ν · f · T2_in · T2_en (or the ratio r = T2_in/T2_en).
- Sensitivity (2-6): d∥ 1.5 · 2.0 µm²/ms (main 1.7) · T2_iso 2000 ms.
- The collinearity of r and f is written alongside as a profile: how the optimum of f (and ν) moves in the profile with r fixed on a grid · the r–f correlation from the profile curvature.

**Approval conditions (relay 2026-11-01 22:05):**

- First of all, check whether TR is the same across the 8 TE series. If it differs, a T1 term enters, so stop and report.
  - The sidecar JSON is not in the received data (json 0 · fixed analysis document §2-6 also says «no JSON sidecar»).
  - **Decision (source 1 = Methods of the original paper · comparison = NIfTI header pixdim[4] only):** the Methods record (`decomp/stage1/analyze3.py:339–340` · PMC12246157) is TR 5800 ms (all 8 series) · δ/Δ 15.2/25.2 ms.
  - Header comparison (2026-11-01 · 8 series of sub-01 · first 348 B of the gzip only · voxels 0 · external data record `x14/EXTERNAL_READ_LOG.jsonl` · `x14/wong_tr_header.json`): pixdim[4] = 5.8 s is the same in all 8 series and the same as Methods.
  - So there is no T1 term between series, and D3 does not stop because of TR. δ · Δ are not in the header, so only the Methods values are used. The ScienceDB page is not read.
  - The headers of sub-02 · sub-03 are not opened separately. D3 uses only the ROI tables (`roi_TE*.json`) and does not load images, so there is no occasion to write pixdim[4] alongside (if images come to be loaded, it is written then in `EXTERNAL_READ_LOG.jsonl` alongside).
  - **Limit (5 · review session verdict · original text):** «Because TR is the same in all series, T1 is not identified from these data. The compartment fractions include the T1 weighting of TR 5800 ms, and apart from the TE/2 term of the spin-echo recovery equation this weighting is the same across TE series. D3 says nothing about T1.»
    - Size of the TE/2 term (one-line arithmetic · not a fit): the recovery factor 1 − 2e^(−(TR − TE/2)/T1) + e^(−TR/T1) changes, over TE 62 → 132 ms, 0.85233 → 0.84890 (−0.40 %) at T1_iso 3 s, 0.72061 → 0.71628 (−0.60 %) at 4.5 s, and by −0.03 % at tissue T1 1.1 s.
- The truth test is run first, and it is also shown that the test fails on deliberately wrong inputs (e.g. input with shuffled TE order · input with swapped shells).
- Results are marked as exploratory. Comparison with the 82 / 52 ms of Gong and the implied r 2.64 of E1 is descriptive only.
- It is not used on single-TE data such as ON-Harmony.

**Data:**

- `s_low` (b 700 — `rconfig.LOW_SHELL["wong"]`) and `s_2000` of `wong/sub-0N/roi_TE62.json` · `own_TE*/roi_TE*.json` «주».
- For each of subjects 3 × ROI 27, 8 TE (62 – 132 ms) × 2 shells = 16 points.
- The per-shell noise SD is computed from the repeats (TE62R2 · TE132R2). The repeats also contain the differences of separately done registration.

**Order:**

1. **Truth test.** With the generator (`decomp/stage1/synth.py` · `decomp/s4/tables.py` — including Watson dispersion), signals are made at Wong's TE · b, noise at the repeat level is added, and then r is recovered.
   - The profile width of r is given alongside.
   - Whether the direction average is independent of ODI is also checked numerically.
2. **Identification check.** If the profile of r is flat (not identified), «not identified» is reported without a measurement fit.
3. **Measurement fit.** If identified, r and the profile interval are given per ROI · subject. The relation with the S4 prior (1.0 – 1.8) is described.
   - If the residuals are systematic across TE, their form is written (material for candidate 10 · the mechanism is not separated).

**Is it possible with our code:**

- **The ROI level is possible.** The signal · relaxation weighting equations already exist (`simulate.weights` · `S.C.apparent` · `synth.ext_coords` · s4 compartment tables).
- **There is no joint fitter.** AMICO uses only a single-TE dictionary. So a new one is written.
  - `scipy.optimize.least_squares` is used (scipy 1.18.1 is in the main environment lock, so no installation is needed).
  - The truth test of 1 is included with it.
- **A voxel-level joint fit (MTE-NODDI style) is not done this time.** It would require reading the 4D images of the 30 series again and doing new registration between series.

**Limits:**

- 9/10 series were registered separately (`wong_judge.json` `own_registration` — in all three subjects, the 9 series other than TE62). Voxels differ slightly by TE.
- The provider's denoising changes the properties of the noise.
- The ROI mean of log values has a Jensen bias.
- Two things to check:
  - Whether TR is the same in each series. If it is, the T1 term is absorbed into the fractions.
  - Whether δ · Δ are the same regardless of TE.
  - Both are checked in the Methods of the original paper or in the sidecars (metadata · optional s3).

**Absolute T2:** obtaining T2_en itself needs the TE decay of b0. `b0_mean.nii.gz` (30 files) of each series would have to be newly read (optional s2), and the receive scale may differ by series. By default it is not done.

**Link with ON-Harmony:** Wong is a different scanner · different people. Even if r is outside the prior, that is only material for candidate 8, not evidence of the cause in the ON-Harmony cells.

**Outputs:** `data/remeasure/explore/d3_wong.json` · truth-test record.

## 6. D4 — TE · T1 · T2 association (a basis for the next researchers · forward model only · new data 0)

**D4.1 — predicted shift of the cards and D1 bias (done after D1).** The committed cards (`decomp/results/stage3/cards/cards.json`) have, per pair, `dTE` · `dTR_ms` · `design_differs` and the predicted shift of relaxation strata L1 · L4 (`algebra` — direction angle · cone width · dx · dy sign · maximum magnitude). These are described side by side with the sign · magnitude of the D1 per-cell bias.
   - The judge knows TE · TR and the relaxation prior as conditions. If the prior is right, the posterior should absorb this shift.
   - If the bias remains towards ΔTE it is material for candidate 8; if not, it leans to candidate 6. Both are written as descriptive only.
**D4.2 — forward sensitivity table (approved · now · data 0).** With the forward equations, a table is made of how much the apparent values move with respect to (r, T2_in, T1 ratio, TE, TR). The ranges are those of ON-Harmony and Wong.
   - In place of NDI · FWF, the apparent signal fractions ν′ · f′_iso and the card coordinates (x′ = logit ν′ · y′ = logit f′_iso − ln ν′) are used (`decomp/stage1/config.py` `weights` · `apparent` · `frac_coords` — the same equations as the oracle observation of S4). The bias of the AMICO fit does not enter.
   - In place of s_low · s_2000, the log of the direction-averaged normalised signal (spherical mean equations of stick · cylinder · free water · d⊥ from the three generators of `config.dperp`) is used.
   - «Under which acquisition differences relaxation moves the apparent values by how much» is organised.
   - Together with the identification result of D3 (in which TE range · shell r is determined), it is kept as a basis for the next researchers.

**Outputs:** `data/remeasure/explore/d4_relax.json`.

## 7. Retraining — «needed · not done»

Depending on the results, only «needed · not done» is written on the relevant line below. There is no re-evaluation on ON-Harmony.

- Retraining with a changed noise · variation model (SNR · N_eff · session variation): (written after D1 · D2)
- Retraining with a widened relaxation prior: (written after D3 · D4)
- Retraining with scanner-specific effects included: (written after D1)

## 8. Computation time (CPU only · training 0)

| Step | Estimate | Basis |
|---|---|---|
| D1 main cells drawn again (C3) | about 4.8 hours | V4 C3 run (main cell + four sensitivity cells = five cells) 10-30 19:25:56 → 10-31 19:28:09 ≈ 24.0 hours ÷ 5 |
| D1 reference distribution (test data · four coil combinations — 1-1) | about 1.7 hours | `test-cover c3` 10-30 18:18 → 18:44 (about 26 minutes) × 4 |
| D1 summaries · bootstrap | a few minutes | — |
| D2 simulation | about 1 hour | 27 ROI × 2 subjects × SNR 4 × 500 = 108,000 cases · S4 test-data simulation 1,000 cases ≈ 33 s (16 threads · `s4/logs/sim_test.jsonl`) |
| D3 truth test + fit | within 30 minutes | ROI 81 × parameters 4 |
| D4 | within 15 minutes | forward model · cards |
| D1 four sensitivity cells — **on hold** | about 19 hours | about 4.8 hours per cell × 4 |
| D1 C7 arm — **not included this time** | about 2.6 hours | expected time of the C7 rerun (about 13 hours) ÷ 5 |

- The total is about 7 hours (when run in sequence). D2 – D4 can run together with D1 (64 cores · current load about 13).
- The C7 rerun is in progress (expected to finish around 11-02 04:16). To avoid affecting C7's time, it is recommended that D1 start after C7 has finished.

## 9. Images newly read (guard · subject to the access record)

**Basic plan: 0 items.** ON-Harmony images · map NIfTI · Wong images · Tong are not read.

What is read is only derived tables already read at the verdict, with sha256 recorded, and committed records:

- `roi.json` — for D1 the pairs of the verdict, for D2 OXF1PRI 12
- `attr_table.json` · `snr_table.json`
- Wong `roi_TE*.json` 30
- S4 frozen weights · test data (value area)
- Cards

Scripts pass through the guard (`guard.require`). The sha256 of each file read is compared with the access record value.

**Re-read record (relay 2026-11-01 22:05):** every file read again, including roi.json, leaves a «re-read» line. Figure work is the same. It is first written in `drafts/x14/REREAD_LOG.jsonl` (time · path · sha256 · bytes · what · script), and moved to the canonical access record in the next approved commit.

**Optional s1 – s4 — not approved (relay 2026-11-01 22:05):**

| Symbol | What to read | Why | Record |
|---|---|---|---|
| s1 | The 7 F2 maps (NDI · FWF · ODI · MD · FA · s_low · s_2000) and mask of the 12 OXF1PRI sessions | voxel SD within the ROI (candidate 11) · spatial correlation → N_eff (candidate 3) | ON-Harmony access record |
| s2 | `b0_mean.nii.gz` of the 30 Wong series | absolute T2 (D3) | external data access record |
| s3 | Wong sidecars · TR · δ · Δ in the Methods of the original paper (not images) | checking the assumptions of D3 | external data access record |
| s4 | eddy text logs of the ON-Harmony sessions (motion · outliers — not images) | candidate 12 | value-area processing record |

**Not done:**

- ON-Harmony raw images · preprocessed DWI 4D (not in the value area — would have to be processed again)
- Wong voxel-level joint fit
- Tong's Technical Validation · verdict/

## 10. Order · stops · location of outputs

**Order (item 4 of relay 2026-11-01 22:05):**

- Now (while C7 is running · light work): ① attribution marking of batch 68 ② sensitivity table (existing outputs) ③ report on whether (F) and W1 – W3 · T1 are possible.
- After C7 has finished and been committed: ④ D1 main cells (up to the reproduction gate) ⑤ D2 ⑥ Wong joint fit ⑦ drawing W1 – W3 · T1 ⑧ V1.
- At each gate the results are reported, and the next gate opens after approval. D4 is not in this order — it is asked again at the time of approval.
- PLAN §9-1 #80 is written as non-reading batch 68 (carried out after the author has checked the original text · basis of #79).

The wording for the manuscript is set separately with the review session after the results.

**Stops:**

- If the D1 reproduction gate fails, stop.
- If CUDA is visible, stop.
- If the guard is closed, stop.
- If the input sha256 differs from the access record, stop.
- If an installation is needed, stop and ask (R8).

**Code:** kept in the repository at `decomp/explore/` (new folder). It is committed after approval, together with truth tests on synthetic data. The verdict code · fixed analysis document · the PIPELINE list of the guard are not changed.

**Outputs:** `data/remeasure/explore/` (value area). The head of every file states «post hoc exploratory — not a verdict · point estimates are not used in results».

## 11. Real-data figures for talks (item 3 of relay 2026-11-01 22:05 — drawing only · new fits · new statistics · new verdicts 0)

- **Selection rule:** written in receipt `RECEIPT_2026-11-01h.md` §2 before opening the images (22:16 KST).
  - Participant = the first participant in sorted order for each dataset (Wong sub-01 · Tong sub-1).
  - Slice = the one axial slice nearest MNI z = +20 mm · same crop.
  - Fixed colour range · viridis (differences RdBu_r) · black background.
- **Figures:**
  - (F) manuscript Fig. 3 · 4 · 5 — done (`x14/F/`).
  - (W1) one row of the Wong pipeline — **only one provider 4D series (sub-01 TE92) is newly read** (external data record · approved). The direction rule was written in receipt i §2 before opening the images: one image at b=0 · the one direction nearest the x axis at b=700 · the three directions nearest the x · y · z axes at b=2000.
    - Two greyscale sets (decision · recorded in receipt j §2 before opening): (A) all tiles at percentile 99 of all voxels of the b=0 volume (no mask) · (B) only the three b=2000 directions at percentile 99 of the voxels of those three volumes combined. The caption states the two scales.
    - Window values (6 · review session verdict): percentile 99 of b0 · percentile 99 of the b2000 group are written in the sidecar, and the caption states which window each panel uses. Where it is used is for the author to set.
  - (W2) **difference maps are left out** (9/10 series were registered separately, so registration mismatch mixes in). The maps of TE 62 · 92 · 132 are placed side by side. The full 8 TE are kept only as supplementary.
  - (W3) card-plane trajectory + wedge 48.5 – 72.1° · (T1) Tong sub-1 centres 10 + retest — all ⑦.
  - (V1) per-ROI-unit z dot plot after the D1 reproduction gate is passed (±1.96 band · posterior means and physical-unit intervals are not drawn) — ⑧. **Two cells: ODI OXF1PRI~OXF2PRI · ODI NOT2ING-98~NOT4GEP (|ΔTE| maximum 30.00).**
  - (ON) ON-Harmony images are not drawn. The License field of dataset_description.json = «CC0» (re-read record · `x14/onh_license.json`). Display of the images is for the author to set.
- **Blocked points (reported · not done):**
  - Standard-space maps are not stored (`run.py:145–161` deletes them after the ROI step). **Approved: MNI slices are made by applying the stored transforms again** (applywarp · the same command as the pipeline `decomp/remeasure/cmds.py:66–68` · maps trilinear · masks nn · reference FMRIB58_FA_1mm). It is for display and statistics are 0. The caption states «interpolated · for display».
  - The b=700 · b=2000 images of W1 are not in our value area — as approved, only one provider 4D series is read.
- **Output:**
  - For talks: PNG without text · tiles 600 × 600 px · sRGB.
  - A sidecar JSON for each file: dataset · participant · session/centre · TE · TR · index · colour range · colour map · slice coordinate · orientation · original path and sha256 · script sha256.
  - For the paper: one more set under the figstyle.py rules (Arial family · 190 mm · vector PDF · English labels).
  - Source wording: "Data: Wongkornchaovalit et al. 2025 (CC BY 4.0); processed in this study" · "Data: Tong et al. 2020 (CC BY 4.0); processed in this study".
- **Location:** scripts and outputs are kept in `drafts/x14/`. The repository commit is made after author approval.
- **Delivery:** the place from which the author downloads is `drafts/x14/download/` (PNG · SHA256SUMS.txt · zip). W1 · W2 (three TE) · W3 · T1 are placed in the same way once drawn.
- **Stop:** if a value appears inconsistent with the existing results table, report it without changing anything. If the fit map of the step to be drawn is not stored, report it without refitting.

## 12. B1 — Wong TE trajectory direction [review session proposal · post hoc analysis set after seeing the results] — **B1 is a recomputation of E1 (2026-10-31 23:33 KST) and not a new test** (review session verdict k reply 2) (definition fixed 2026-11-01 23:33:24 KST · before computing)

- **Prior check:** the same direction angle was already computed in exploratory E1.
  - 2026-10-31 23:33 KST · `drafts/manuscript_tex_2026-10-31/exploratory/RESULTS.md` · `e1_results.json` · definition `E_DEFINITIONS.md` (fixed 23:30:15).
  - ψ overall median 25.4° · inside the wedge 10/81 · implied r median 2.64.
  - So B1 is an analysis after seeing those results.
- **Data:** `judge_wong.load` (same path as manuscript Fig. S2 · E1).
  - Participants 3 × ROI 27 = 81 ROI units · eight TE (62 – 132 ms).
  - The repeats TE62R2 · TE132R2 are left out of the direction fit and used only for repeat differences.
  - roi_TE*.json has a re-read record · new images 0.
- **Coordinates:** x = logit NDI · y = logit FWF − ln NDI = ln(free water/intra) (`lib.coords`). Both axes are natural-log ratios and so have the same unit, and the angle is measured from the x axis (same as manuscript §2.2).
- **Direction:**
  - b_x · b_y = `judge_wong.fit_unit(te, v)["slope"]` (linear slope centred at 97 ms — same as E1).
  - θ_obs = atan2(b_y, b_x) (degrees · (−180°, 180°] · direction of increasing TE).
  - Comparison: for each ROI unit, if not |θ_obs − ψ_E1| ≤ 1e-9°, stop.
- **Interval:**
  - For each ROI unit, the 2×2 covariance R of the x · y straight-line fit residuals (degrees of freedom 6) → (b_x, b_y) covariance R ÷ Σ(TE − 97)².
  - 10,000 draws from the bivariate normal, θ unwrapped around θ_obs, then percentiles 2.5 · 97.5 are computed (seed 20261141).
- **Wedge:** 48.5 – 72.1°.
  - The D4.2 values as they are (`x14/D4/d4_2_forward.json` sha256 dfc2aa513b1fc7bc… · S4 prior bounds · same as the stage 2 REPORT).
  - Centre θ_c = 60.3° · half-width 11.8°.
- **Report (descriptive only · no verdict · no threshold):**
  - Per-participant θ_obs median · range.
  - The fraction of ROI units inside the wedge (point estimate). For reference, the fraction whose interval overlaps the wedge.
  - The TE slope b_⊥ = −sin θ_c · b_x + cos θ_c · b_y (per 10 ms) of the perpendicular component p = −sin θ_c · x + cos θ_c · y.
  - Comparison with same-TE repeat differences: d_p = p(R2) − p(original) (TE62 · TE132 · ROI 27) → per participant SD_rep(p) = √(mean d_p²) → ratio = |b_⊥| × 70 ms ÷ SD_rep(p). Given as per-participant median · range.
  - The repeat series were registered separately, so SD_rep contains registration differences.
- **Arithmetic check (3 · not a computation):** for any direction inside the wedge, the share of the TE shift leaking perpendicular to the reference direction is at most |sin(θ − θ_ref)|.
  - If the reference is the wedge centre (60.3°), sin 11.8° = 0.204.
  - If the reference is 57.9° (Gong 82/52 · the «58°» of manuscript §2.2), it is 14.2° (72.1 − 57.9) to the farthest end, so sin 14.2° = 0.245.
- **The disposition is not reopened (#34).** Outputs `data/remeasure/explore/b1/`.

## 13. D3 extension — true-r grid out beyond 2.64 [review session proposal · awaiting author approval] (design 2026-11-01 23:53:38 KST · before computing)

- **What to fix:** the true r of the previous truth test was only 1.0 – 1.8 (`config.R_GRID`).
  - The fitter range (r 0.5 – 5) and the profile grid (0.6 – 4.0) already extended beyond 1.8.
  - But whether r ≈ 2.6 is recovered was not tested (pointed out by the review session).
- **New grid (set by A):**
  - True r = {1.0, 1.15, 1.3, 1.45, 1.6, 1.8, 2.0, 2.3, 2.64, 3.0} × T2_in {70, 82, 95} × WM tissues 12 = 360.
    - Appended only to the `config.T2_SETS` list at run time. config.py is not changed.
    - Reason: it contains the E1 implied value 2.64 and r ≈ 2.9 corresponding to the 25° of B1, and it evenly spans the literature value 1.58 and the range between.
  - Fitter bounds r 0.4 – 8.0 · starting points {0.8, 1.2, 1.6, 2.2, 3.0, 4.0}.
  - The profile grid is log-spaced 0.5 – 6.0, 51 points. It was widened so that the interval of the true 3.0 does not touch the grid end.
- **Tests to run again:** V5 · A-1 (360 · noise 0 · |Δ| ≤ 1e-3) · A-2 (TE shuffled · shells swapped — must fail) · B (generator B · noise 0 · bias · d∥ · T1 sensitivity) · C (noise σ · 30 times per true value).
- **Group estimation (new):**
  - Participant group: the ROI 27 share (T2_in, r), and (ν, f) are separate per ROI.
  - Overall group: one r · T2_in per participant · (ν, f) per ROI unit.
  - Group truth test: generator B makes groups of ROI 27 (random tissue · noise σ). For each true r ∈ {1.3, 1.6, 2.0, 2.64, 3.0}, 30 groups are made and the coverage · width of the profile 95 % interval are computed. Then the fraction in which the reading rule below reads correctly at true 1.6 · 2.64 respectively is computed.
- **Measurement fit (after passing the truth test):** `roi_TE*.json` (re-read record) · three participant groups + one overall group · sensitivity d∥ 1.5 · 2.0 · T2_iso 1000 · 2000. The r of a single ROI unit is not interpreted.
- **Reading (before computing · review session wording as it stands · not a verdict · the disposition is not reopened #34):**
  - If the group interval contains 1.6 and does not contain 2.6: «the T2 ratio does not explain the shallow angle of B1».
  - In the opposite case (it contains 2.6 and does not contain 1.6): «the T2 ratio is larger than the literature».
  - If it contains both or neither: «verdict on hold».
- **Time:** truth test about 20 – 30 minutes (8 processes) + measurement a few minutes.
- **Approval (author · chat choice in the 00 o'clock hour of 2026-11-02 «D3 extension → measurement fit» · «B1-S» · «start D1 right after the C7 commit» — 4a not chosen).**
- **Conditions for moving on to the measurement (before computing · 2026-11-02 00:09:25 KST):** the measurement fit is run only when all three below hold. If any one fails, it is reported without the measurement.
  - A-1 passes
  - A-2 both fail
  - In the group truth test, coverage of the 95 % interval ≥ 0.85 at true 1.6 and at 2.64 respectively
- **The participant T2_in of the group truth test is 82 ms.** The overall group (three participants) uses the same framework and has no separate truth test (written as a limit).

## 14. B1-S — simulation reproducing the Wong processing as it stands [review session proposal · awaiting author approval — includes new image reads] (design 2026-11-01 23:53:38 KST · before computing)

- **Purpose:** to settle explanation (B) «noise · fit bounds flatten the trajectory».
- **Generator:** project generator `synth.noiseless`.
  - T2: 82 / 52 / 1000 (algebra 58°) · both ends of the wedge (T2_in 95 · r 1.8 → 48.5° / T2_in 70 · r 1.3 → 72.1°) · and, as a positive control for the B1 code, r 2.9 (T2_in 82 → 25.8°).
  - Tissue: WM-like ν {0.45, 0.6, 0.75} · ODI {0.08, 0.20} · FISO {0, 0.02, 0.05, 0.10} (including near 0).
- **Acquisition:** Wong as it stands — b700 · 2000 · TE 62 – 132 (8) · TR 5800 · δ/Δ 15.2/25.2 (Methods record).
  - The directions are the bvec of Wong sub-01 TE92 (metadata · newly read).
  - The protocol object is made with those directions (at run time only · config.py not changed).
- **Noise:** non-central chi. There is no record of the Wong scanner's number of channels L, so **4 · 8** are taken as sensitivity values (fixed 2026-11-02 00:16:23 KST · before computing — the 1 · 4 first written cannot be accepted by the frozen function `simulate.noisy`: the χ² degrees of freedom 2(L−1) are 0 at L = 1, so the first run stopped with «df <= 0» · results 0 · the same as the 32 · 64 channel mapping of S4 `s4cfg.COIL_L`).
  - SNR is measured per TE series from Wong's b0 repeats (4 b0 per series). The formula is the same as `analysis.snr` (brain mask ∩ FA ≥ 0.3 · using nodif_brain_mask · dti_FA of our value area).
  - **Images newly read = the 8 TE series 4D of sub-01 (only b0 volumes used · external data record).**
  - SNR × 0.5 · × 2 sensitivity.
- **Processing:** the same F2 as the Wong processing (AMICO 2.1.1 NODDI · d∥ 1.7 · d_iso 3.0 · same bounds) is used.
  - Kernels for the Wong direction scheme are made in a new study folder.
  - ROI aggregation is a voxel mean (200 voxels per ROI · voxel heterogeneity is PIX_SD of S4).
  - Coordinates `lib.coords` → B1 code (`judge_wong.fit_unit` slope → θ).
- **Validation (gate before computing):**
  - At noise 0, θ must be within ±3° of the algebraic angle (82/52/1000 → 57.9°). Otherwise stop (pipeline bug or fitter bending — report).
  - On synthetic straight lines of known direction (θ = 25° · 58° · 72°), the B1 code must hit within ±1e-9°.
- **Reading (before computing · not a verdict):** reference line 38° (the midpoint ≈ 37.7° of the observed 25 – 28° and the lower wedge end 48.5°, rounded up). At the measurement SNR:
  - If θ_sim ≤ 38°, «can be explained by (B)».
  - If θ_sim > 38°, «cannot be explained by (B)».
  - The TE slope of ODI is reported alongside.
- **Limit (original text):** Wong is data already denoised by the provider, so the original noise distribution cannot be reproduced exactly.
- **Time:** kernels a few minutes + simulation about 10 minutes + reading 8 4D files a few minutes.

## 15. D3-v3 — continuous profile · union over systematic errors · truth test with realistic noise (author approval «run it now and include it all at once» · definition 2026-11-02 01:15:58 KST · before computing)

- **Reason:** the group coverage of the D3 extension truth test (§13) fell short of the condition for moving on (1.6 0.20 · 2.64 0.50). So the measurement fit was not done.
  - Diagnosis (after seeing the results): intervals narrower than the grid spacing · about 1 % model bias · no noise floor in the test noise.
- **(i) Continuous profile:**
  - RSS(r) = min_{T2_in} min_{(ν, f)} RSS. The inner part is the same as §13.
  - r is bracketed by a coarse grid of 13 points, log-spaced over 0.5 – 6.0, and r̂ is then found by bounded Brent (xatol 1e-4).
  - The interval ends are the roots of RSS(r) − RSS_min = 3.84 (brentq). If it does not cross by the end of the range, it is written as «open».
- **(ii) Systematic error:**
  - An interval is computed for each of the six variants (d∥ {1.5, 1.7, 2.0} × T2_iso {1000, 2000}).
  - Reading uses their union (smallest lo ~ largest hi). The reading sentences are those of §13 as they stand.
- **(iii) Truth test with realistic noise:** the B1-S setup is used.
  - Wong sub-01 TE92 directions · signal equation V-S 0.0 · simulate.noisy · L4 · SNR(TE) = b1s_snr.json × {1, 2}.
  - Tissue group = WM tissues 12 (config.WM) × fibre directions 8 = 96 ROIs (voxels 200 · S4 PIX_SD · fixed seed). For each group, 27 of them are drawn (with replacement).
  - Per voxel, s_low · s_2000 are computed by synth.ext_coords (low 700) and the ROI mean is taken. The order is the same as the measurement pipeline.
  - True T2_in 82 · r ∈ {1.6, 2.64} · T2_iso 1000 · d∥ 1.7 · T-in · T1 (1100, 1.0) · TR 5800.
  - 20 groups per (true value · SNR multiple) → 80 groups.
  - Noise σ (the residual scale of the fit) is the same as §5 (Wong repeats b700 0.0128 · b2000 0.0135).
- **(iv) Conditions for moving on to the measurement:** the two below must hold at both SNR × 1 · × 2.
  - Coverage of the union interval ≥ 0.85 (true 1.6 · 2.64 respectively)
  - The fraction in which the reading rule reads correctly ≥ 0.8 (1.6 → «the T2 ratio does not explain the shallow angle of B1» · 2.64 → «the T2 ratio is larger than the literature»)
  - If any one fails, it is reported without the measurement.
- **Measurement (only when the conditions hold · same script):**
  - Wong roi_TE*.json (re-read record) · three participant groups + overall group (one r · sum of RSS).
  - Same approach as (i) · (ii). The r of a single ROI unit is not given.
- **Time:** about 30 – 60 minutes (16 processes · depending on load).

### 15-1. D3-v3 run record (addendum 2026-11-02 01:33:30 KST · the definition above is not changed)

- **Self-check (not included in the main computation · 01:19 – 01:29):** one group (true 2.64 · SNR × 2 · composition gi 0) was run with the six variants. 630 s.
  - Result: in all six variants r̂ = 6.0 (upper end of the coarse grid) · union interval [5.845, 6.0] · upper side «open» · does not cover the true value.
- **Implementation check V-N (written at the head of the script `d3_v3_check.py` at 01:33 and run):** on the same setup with only the noise left out · main variant · true 1.6 · 2.64 all |r̂ − r| / r ≤ 0.05.
  - Result: passed — r̂ 1.5898 (−0.6 %) · 2.6323 (−0.3 %). Intervals [1.544, 1.637] · [2.499, 2.778].
  - → The 6.0 of the self-check is taken to come from the noise (non-central chi), not from the fit code.
- **Main computation:** started 2026-11-02 01:32:00 (80 groups · 16 processes · `explore/d3/d3_v3.log` · `STATUS_D3V3`).

### 15-2. Failure of the first main computation and fixing the implementation (addendum 2026-11-02 02:45:12 KST · main results 0 · definitions (i) – (iv) unchanged)

- **Failure:** the 01:32 – 02:41 main computation stopped at the end. The root search for an interval end of one group (brentq) gave «f(a) and f(b) must have different signs».
  - pool.map raised the exception after finishing all groups, so no result at all could be written (no d3_v3.json · main results 0 · log `d3_v3_run1_failed.log` · the script at that time `d3_v3_run1_failed.py`).
  - The «rc=0» in the status file is a wrong value from a shell quoting mistake (in fact it ended with an exception).
- **Cause (confirmed in the code):** the inner starting point of RSS(r) was changed at each call to the solution of the immediately preceding call. So even for the same r, the value differed depending on the order of calls.
  - A point that exceeded 3.84 on the coarse grid may not exceed it when the root search recomputes it. The value at r̂ can also change on recomputation. Then the signs of the two ends become the same.
- **Fix (written before computing):**
  1. RSS(r) is made a deterministic function of log r.
     - The 13 coarse-grid points are taken in ascending order, starting from the solution at the previous point (the first point starts from (ν, f) = (0.6, 0.1) per ROI).
     - Other r start from the solution at the coarse point nearest in log distance.
     - A log r once computed reuses the remembered value. So by definition the signs of the two ends of the root search are not inconsistent.
  2. r̂ is whichever of the Brent result and the coarse-grid minimum has the smaller RSS.
  3. The outer end of the root search is the coarse point nearest from the r̂ side that exceeds 3.84 (including the points of the grid cell containing r̂).
  4. Exceptions are caught and written per group. A group with an exception is counted as «not covered · not read correctly» (conservative). Such groups are listed.
  5. Group results are written line by line to `d3_v3_groups.jsonl` as they finish. The measurement fit also catches exceptions per job. In the status file, the exit code is captured first and then written.
- **Checks before relaunching:**
  - V-N is run again with the fixed version (pass condition unchanged).
  - Determinism check V-D: with the same input of the self-check group (true 2.64 · × 2 · composition gi 0 — a group already seen), the main variant is run twice. The results must be bit-identical.
  - If both pass, relaunch.

### 15-3. D3-v3 results (addendum 2026-11-02 04:19 KST · the definition above is not changed)

- **Relaunched main computation:** started 02:47:33 · ended 03:54:02 · exit code 0 · 80 groups · exceptions 0.
  - Script sha256 6a1ff731… · outputs `explore/d3/d3_v3.json` (00c548ac…) · `d3_v3_groups.jsonl` (51057ca0…) · summary `d3_v3.md`.
- **(iv) Conditions for moving on: not met → the measurement fit was not run** (as defined, «reported without the measurement»).

| True r · SNR multiple | Union interval coverage | Fraction read correctly | r̂ median (main variant) | Union interval median | Reading |
|---|---|---|---|---|---|
| 1.6 · × 1 | 0.55 | 0.00 | 6.000 | [0.500, 6.000] | verdict on hold 20 |
| 1.6 · × 2 | 0.00 | 0.00 | 4.479 | [3.968, 5.424] | verdict on hold 20 |
| 2.64 · × 1 | 0.95 | 0.00 | 6.000 | [0.500, 6.000] | verdict on hold 20 |
| 2.64 · × 2 | 0.00 | 0.00 | 6.000 | [5.851, 6.000] | verdict on hold 20 |

- **Descriptive (after seeing the results · not set before computing · conditions and reading not changed · `d3_v3_describe.py` → `d3_v3_describe.md` · json f3cd2ce9…):**
  - × 1: almost all r̂ of the six variants reached a range end, 0.5 or 6.0 (true 1.6 119/120 · 2.64 114/120).
    - The coverages 0.55 · 0.95 come from groups in which r̂ went to different ends by variant, so that the union interval became the whole of [0.5, 6.0]. It is not coverage that singles out the true value.
  - × 2: for true 1.6, r̂ is inside the range in all 20 groups and the median is 4.479. For true 2.64, 19/20 are at the upper end 6.0.
  - Misfit (main variant RSS_min median · 432 data points): × 1 9489 · 10292 (22.0 · 23.8 per data point) · × 2 621 · 1145 (1.4 · 2.7). At noise 0 (V-N) it is 0.37 · 0.10.
- **Implementation check V-G (after seeing the results · rule written at the head of the script at 04:00 and run · `d3_v3_check2.py` sha256 4971b451… · `d3_v3_check2.json` 002ee5d1…):**
  - Question: at × 2, is the true value being pushed outside the interval due to the inner optimisation (bounded Brent of T2_in · (ν, f) starting points) staying in a local solution?
  - Method: for the 40 groups of × 2 · main variant, at r = true value and at r = r̂, a full fit over a T2_in log grid of 41 points (20 – 300 ms) × three starting points, and the neighbourhood of the best grid cell is refined by bounded Brent.
  - Rule: if ΔRSS_exh = min(true value) − min(r̂) > 3.84 holds in all groups, «true value outside the interval» is not due to the inner optimisation.
  - Result: **passed** (40/40). ΔRSS_exh minimum · median = true 1.6 585.8 · 975.9 / 2.64 273.3 · 506.7.
  - What also came to light: groups in which the grid search at r̂ is more than 3.84 lower than RSS_min of the main computation are 22/40 (true 1.6 2/20 · 2.64 20/20 · largest difference 193.6). The T2_in search of the main computation can stay in a local solution when noise is present. So the RSS values of the main computation are read as upper limits of the profile. The comparison above measured the two points with the same grid search, so the conclusion stands.
  - Optimal T2_in of the grid search (grid cell · true 82 ms): at r̂ 152 – 245 ms · at r = true value 133 – 300 ms (for true 1.6, 18/20 are at the range end 300).
    - (fixed 04:22: at first, looking only at 22 groups, it was written «at r̂ 133 – 214 · 300 in the two groups of true 1.6». It was recounted over all 40 groups.)
  - In the first run (04:01:04), the time at the head was wrongly written as «04:03». So it was run again with only the time changed (04:09:59 – 04:18:12), and the per-group results of the two runs are the same (`d3_v3_check2_run1.json`).
- **A's reading (after seeing the results · not a verdict):** under this noise model (non-central chi L4 · measurement SNR × 1 · × 2), the multi-TE group fit does not recover the T2 ratio.
  - At × 1, r̂ goes to a range end, and at × 2 it is pushed far upwards (true 1.6 → about 4.5).
  - So D3 cannot settle, on the Wong data, «whether the T2 ratio is larger than the literature».
  - The limit of B1-S (provider denoising · the original noise distribution cannot be reproduced exactly) still applies.

## 16. D1 implementation definitions — what was set when moving §2 into code (addendum 2026-11-02 04:39 KST · before computing · the definitions of §2 are not changed)

- **Start:** after the C7 results commit script (`drafts/scripts/x14_c7_commit.sh` · sha256 dbcaadd6…) has finished commit · push · guard check, it launches `x14_d1_chain.sh` (14450b4a…). If the commit stops, D1 is not launched either.
- **Script:** `drafts/x14/D1/d1_c3.py` (sha256 a0370349…). Four modes: real ‖ (ref → sens) → summary. If real or ref fails, summary is not run.
- **1 Drawing again:**
  - judge_v4 is not changed. draw_fn is wrapped so that it returns the received values as they are (computation · random number consumption change 0), and only PIT · posterior mean · SD are kept per ROI unit · coil combination · k · component (same approach as the C7 monitoring driver).
  - The pair order is card_pairs, the seeds are v4_c3.json seeds_torch (checked that they equal 6000 + pair number).
  - Functions that write V4 result files · the access record (_write_once · _write_json · _log) are blocked.
  - real has the same settings as the original runner: CUDA_VISIBLE_DEVICES= · S4_THREADS=16 · OMP · MKL · OPENBLAS variables are not set (stops if present).
- **2 Reproduction gate:**
  - With the «main» variant of the verdict (main-cell excluded subjects left out), the same functions cell_verdict · coil_counts are called again and compared exactly with the stored values.
  - Each time a pair finishes, the cells of that pair are compared. If even one mismatches, it stops on the spot.
  - Written when it stops: number of mismatched ROI units (formula of §2-2) · difference in by_coil_combo counts · list of cells not compared. D1 values are not given.
- **3 PIT:**
  - The fraction of finite samples smaller than 0 (for test data, the true value).
  - If there are several coil combinations, the PIT farthest from 0.5 is used (if equal in distance, the earlier combination).
  - **z is also taken from that combination.** §2-5 has no combination rule, so it is set here.
- **ROI units used:** the same as the verdict — leaving out main-cell excluded subjects, leaving out misspecification, and using only ROI units with samples for all nine k in all combinations. The sensitivity analysis is the value with the misspecification ROI units included.
- **Summary:** per k, PIT mean · two-tail fraction (PIT < 0.05 or > 0.95) · z mean · z SD (ddof 1) · PIT SD → median · minimum · maximum over the nine k.
- **Bootstrap:** subject clusters 2,000 replicates · numpy seed 20261102. Percentiles 2.5 · 97.5 are computed for each k median of PIT mean · two tails · z mean · z SD.
- **4 Reference distribution:**
  - 1,000 cases are used per acquisition type of the test data L0. The cells with data have four types (match · te · tetr · tetr_design — all tr cells have no data).
  - Drawn with four combinations (32|32 · 32|64 · 64|32 · 64|64) × nine k. The seed per combination is 7000 + 10 × type number + combination number.
  - Cases included per cell: the true combination is among that cell's candidates · samples for nine k are full in all candidate combinations · not misspecification (misspecification if the scores of all candidate combinations exceed the threshold).
  - On the same set of cases, the rule PIT · z and the PIT · z of the one true combination are computed. «Two-tail share added by the rule» = two tails (rule) − two tails (true combination).
- **6 Variance decomposition:**
  - (a) = pooled SD of what is left after subtracting each ROI's mean (degrees of freedom = number of ROI units − number of ROIs).
  - (b) = SD of the ROI means (ddof 1). Comparison value = (a) ÷ √(mean number of ROI units per ROI).
  - Computed per k, then the median is used. Both z · PIT are given.
- **7 Noise prediction:**
  - Observed ratio = per k, (a)_z ÷ reference (rule) z SD → median. Predicted ratio = median of 30 ÷ SNR_pair over the ROI units used.
  - ROI units outside and inside the training range (snr_table train_range) are summarised separately.
- **8 Descriptive account of relations:**
  - The card record is `decomp/results/stage3/cards/cards.json` (ON-Harmony). It was confirmed that dTE · dTR_ms · design_differs are the same as stage3_F1. dTE · dTR_ms are B − A.
  - Attributes used: |ΔTE| · ΔTR · design_differs · coil (plan.ONH · «unknown» if absent) · cell median of ΔSNR (B − A).
  - Between cells: Spearman (descriptive · no p values). For design · coil, the median is given per group.
  - The six pairs of the coil-only contrast are written side by side.
  - Within a cell: Spearman of ΔSNR and the k median of the mean z (· PIT) of the ROI units used per subject (cells with 4 or more subjects).
  - Figures (for explanation · value area): cell z mean vs |ΔTE| · ΔTR · ΔSNR · observed ratio vs predicted ratio.
- **9 Sensitivity misspecification (sens):**
  - The main-cell misspecification count is recomputed the same way and compared with the stored value (check).
  - The k-unmet count of the four sensitivity cells is given as «number of ROI units − misspecification count» only when the stored pass flag is None and there are 3 or more subjects.
  - If None and subjects < 3, «unknown (subjects < 3)»; otherwise «unknown (samples needed)».
- **Re-read record (`drafts/x14/REREAD_LOG.jsonl`):** attr_table · v4_c3 · snr_table (sha256 compared with the access record values · stops if different) · all the roi.json read by collect.
- **Outputs (value area `explore/d1/`):** d1_gate.json · d1_c3_units/ (npz per pair · not used in results — d1_c3_units.npz of §2 split per pair) · d1_ref_<type>.npz · d1_sens_flags.json · d1_c3.json · d1_c3.md · two figures · STATUS_D1.
- **Report:** per cell, only PIT · z · (a) · (b) · number of ROI units · coverage are given. Per-ROI-unit posterior means · physical-unit intervals are not given.
- **Tests (on synthetic data without the judge · before this addendum):** the wrapper's combination · k indices · PIT · combination selection (equal distance → earlier combination) · variance decomposition · reproduction comparison (same data → same · one coverage changed → mismatch 1 · one combination count changed → mismatch 1) were all right.
- **Time (rough):** real 4 – 5 hours (16 threads) ‖ ref 1.5 – 2 hours (16 threads) → sens · summary a few minutes.
- **Fix (2026-11-02 09:13 KST · before computing):** so that the C7 results commit does not stop even if the author commits non-reading batch 68 first, the commit script and the handoff script were made to accept as parent d11daaa or a commit «directly on top of it · one changed file PLAN.md · title ‘non-reading batch 68 — ’» (x14_c7_commit.sh 19143133… · x14_c7_handoff.py fa5a1f2e…). With any other HEAD it stops. The D1 definitions stay as they are.

## 17. D5 — device fingerprint · mapping between devices · within-scanner outliers [exploratory · author approval] (definition 2026-11-02 09:40 KST · before computing)

- **Approval:** author chat «let's proceed, please design it rigorously» (2026-11-02 09:35:49 KST). This is read as approving all three of A's proposals (device identification · transfer prediction · within-scanner outliers) and reading the maps of 119 sessions (widening §9 s1) · the eddy text logs (§9 s4).
- **Questions:**
  - ① Device identification — can the 7 devices be told apart from eddy records alone, or from tissue values (ROI) alone (cross-validation leaving subjects out)? → Where does the device fingerprint remain?
  - ② Prediction of the mapping between devices — if the value on device B is predicted from the same person's value on device A, does it reach the retest level? → Is the device difference a stable systematic effect?
  - ③ Within-scanner outliers · typical-class clustering — are there sessions · ROI units that are strange for their own scanner?
  - ④ Linking with D1 (after D1 has finished) — do failures (misspecification · |z| · PIT two tails) concentrate in outlier ROI units? Do pairs with a good mapping and the D1 bias move together?
- **Sessions:** rconfig.onh_sessions() 110 + sref_extra() 10 = 120. The 119 whose processing finished (excluding sub-16793_ses-OXF4GEP001 «retry also failed»).
- **Reads:**
  - Newly read — images (value area · §9 s1 extended to 119 sessions): per session f2_{NDI · FWF · ODI · MD · FA · s_low · s_2000} · nodif_brain_mask · f2_extvalid · f2_f2valid · fa2std_warp.
  - Newly read — not images (§9 s4): eddy_out.eddy_{movement_rms · outlier_map · parameters · post_eddy_shell_alignment_parameters} · bvals.
  - Re-read: roi.json «주» · attr_table.json · snr_table.json (sha256 compared with the access record values).
  - Record: each read leaves the time · path · sha256 · bytes. New reads go in `drafts/x14/ONH_NEW_READ_LOG.jsonl` (including the guard basis), re-reads in `REREAD_LOG.jsonl`. New reads are moved to the canonical access record (onharmony_image_access_log.json «열람» (reads)) in the record commit after the C7 results commit. guard.require is called before reading.
  - Nothing is written in the session folders (processing outputs). Standard-space maps are made in `explore/d5/tmp/<session>/` (<session>) and deleted after extraction.
  - Verdict outputs · verdict/ · Tong Technical Validation are not opened. GPU 0.
- **Features (fixed before computing):**
  - E — eddy · session level. Only volumes with b > 50 are «diffusion volumes».
    - E1 mean absolute RMS movement (movement_rms column 1 · mean over all volumes · mm)
    - E2 mean relative RMS (column 2 · mean from the second volume on)
    - E3 outlier slice fraction (number of entries equal to 1 in outlier_map ÷ (number of diffusion volumes × number of slices))
    - E4 maximum relative RMS
    - E5 linear eddy-current magnitude (diffusion-volume RMS of eddy_parameters columns 7 – 9)
    - E6 maximum shell-alignment shift (maximum |value| of the three post_eddy_shell_alignment shifts · mm)
  - T — tissue · session level: roi.json «주» 27 ROI units × (x · y · ODI) = 81. Coordinates computed from ROI means with lib.coords (C3) (same as the verdict). ROI units with no value in even one session are dropped from the feature set.
  - P — maps · session × ROI unit level. Voxels of the same path as the verdict (applywarp → standard space → FMRIB58 skeleton (≥ 2000) ∩ JHU labels ∩ main valid mask) are used. ROI unit = left and right label voxels combined.
    - P1 · P2 · P3 = SD (ddof 1) of x · y · ODI computed per voxel with lib.coords (C3)
    - P4 = FWF percentile 90 (proxy for partial volume · CSF mixing)
    - P5 = valid fraction (main valid voxels ÷ skeleton ∩ label voxels)
    - P6 = log number of valid voxels
  - G — registration · session level: Pearson r of skeleton voxels between standard-space FA and FMRIB58_FA.
  - Session summary P̄ = median over the 27 ROI units (P1 – P5).
- **Reproduction gate V-R (before ③):**
  - The roi.json «주» 27 ROI units × 7 values recomputed by the same path must agree with the stored values to |Δ| ≤ 1e-12, and voxel counts must be the same.
  - If even one session mismatches, ③, which uses P · G, stops and reports (① E · ① T · ② go on as they are).
- **① Device identification:**
  - Sessions: among first sessions, V4's usable (the 1 main-cell excluded session left out · included as a sensitivity analysis). Label = device 7 (the two NOT2ING protocols are one device).
  - Feature sets: E (6) · T (81) · P̄ + G (6) · E + T.
  - Model: standardise within the training fold only → shrinkage LDA (sklearn LinearDiscriminantAnalysis · solver lsqr · shrinkage auto).
  - Validation: leave one subject out (all that subject's sessions left out). Measures = balanced accuracy · confusion matrix.
  - Chance reference: percentile 95 of the distribution from running the same validation 1,000 times (seed 20261103) with device labels shuffled within subject. If the observed value exceeds it, «above chance level» (descriptive).
  - Secondary: NOT2ING-89 vs 98 (same scanner · only the protocol differs) with T. Also given with T_c, T subject-centred (the subject's mean over devices subtracted).
- **② Prediction of the mapping between devices:**
  - ROI units: V4's pairs × subjects × ROI 27 (first sessions of the attribute table · main variant · 23 pairs with data) · coordinates x · y · ODI.
  - Mappings (per pair · coordinate):
    - M0 as it is (b̂ = a)
    - **M1 shift per ROI (main):** b̂ = a + mean (b − a) of the training subjects
    - M1g one shift per pair (common across ROIs)
    - M2 linear per ROI (b̂ = α + βa · only pairs with 6 or more training subjects)
  - Validation: leave one subject out.
  - Retest scale: OXF1PRI repeats (sub-13192 · sub-16974 × ses-OXF1PRI001 – 006 · roi.json «주»). Within-subject variance pooled per ROI (degrees of freedom 10) → σ_d,r = √2 · σ_w,r (SD of the difference of two sessions).
  - Measure: RMS of the standardised error e/σ_d,r = R (per pair · coordinate · across subject × ROI). If M1 is «perfect mapping + OXF1PRI level noise», the expected R ≈ √(n/(n − 1)) (n = number of subjects of the pair).
  - Interval: within-pair subject bootstrap 2,000 times (seed 20261104) 95 %.
  - Secondary: R_between = RMS(e) ÷ between-subject SD on device B (per ROI · pooled) — relative to individual differences.
  - Reading (before computing · [A's proposal]):
    - Upper end of the 95 % interval of M1's R ≤ 1.5 × expected → «mapped at retest level»
    - Lower end > 1.5 × expected → «outside retest level»
    - Otherwise → «cannot be separated»
  - Limit: the retest scale is that of the one device OXF1PRI (SNR 15 – 20). Pairs of devices with lower SNR can have large R even with a perfect mapping. So the SNR of the two sessions is written alongside.
- **③ Within-scanner outliers · typical-class clustering:**
  - Session level (119 sessions · repeats included). Sets Q and S are computed separately.
    - Set Q = (log10 E1 · log10 E2 · logit E3′ · log10 E4), E3′ = (k + 0.5)/(N + 1).
    - Set S = (P̄1 – P̄5 · G).
  - Within-device robust standardisation: per device, per feature (v − median)/(1.4826 · MAD).
  - Pooled robust covariance: MCD over all standardised rows (sklearn MinCovDet · random_state 20261105) → robust distance d². Outlier = d² > χ²_p 0.99 quantile.
  - Stability: subject bootstrap 500 times (seed 20261106 · median · MAD · MCD recomputed). If the fraction flagged as outlier among the resamples is 0.8 or more, «stable outlier».
  - ROI units (session × ROI · about 3,200 rows): features P1 – P6. Robust standardisation within (device × ROI) → pooled MCD → χ²_6 0.99. Stability bootstrap 200 times.
  - Check (does not change results): how the main-cell excluded session (sub-16794_ses-NOT4GEP001 · synthetic b0 topup) comes out in Q · S is written.
  - Typical-class clustering (descriptive):
    - (A) between scanners — GMM (diagonal · k 1 – 8 · BIC · n_init 10 · random_state 20261107) on globally standardised E + P̄ + G. The ARI of clusters and devices is computed.
    - (B) within scanner — the same GMM (k 1 – 5) on within-device standardised Q + S. Bootstrap ARI stability (200 times). If k = 1, «one typical class · outliers only».
- **④ Linking with D1 (after D1 summary · this rule written before seeing the D1 values):**
  - ROI unit = D1 ROI unit (pair × subject × ROI · two sessions).
  - Flags:
    - S-outlier = at least one of the two sessions is a session-level «stable outlier» (Q or S)
    - U-outlier = at least one of the two sessions' rows for that ROI is a ROI-unit «stable outlier»
  - Outcomes (z · PIT only for the ROI units used in D1 · misspecification for all ROI units except main-cell excluded subjects):
    - Misspecification fraction
    - k median of |z|
    - PIT two-tail fraction (mean of the two-tail flags over the nine k)
  - Summary: per-cell (outlier − typical) difference → between-cell median · mean weighted by number of ROI units · within-cell subject bootstrap 2,000 times (seed 20261108) 95 %.
  - Pair level: Spearman of R of ② (M1 · x · y · ODI) and PIT mean · two tails · z SD of D1 (descriptive).
  - Reading (before computing · [A's proposal]):
    - If the 95 % interval of the aggregated difference clears 0 and the outlier side is worse → «failures concentrate in outlier ROI units — the motion · partial volume · registration candidates (2 · 11 · 12 of §4) rise»
    - If the interval contains 0 and failure in typical ROI units is similar to the whole cell → «failure is the same in typical ROI units too — those candidates weaken»
- **Outputs (value area `explore/d5/`):** features_eddy.json · features_roi.npz · features_map.npz · vr_gate.json · d5_1_device.json · d5_2_transfer.json · d5_3_outliers.json · d5_4_link.json · d5.md.
- **Report:** as descriptive · exploratory only. It is written in terms of mapping · prediction · difference (harmonize · correct · remove · restore are not used). The verdict · disposition stay as they are (#34).
- **Computation:** CPU only. Map extraction 8 processes, the rest 1 process (OMP 1).

### 17-1. Fix before computing — outlier threshold · stability (addendum 2026-11-02 09:47 KST · real-data computation 0 · found in a synthetic-data test)

- **What was found:** testing ③ on synthetic data (three devices · 120 rows · 4 dimensions · normal) gave 3 – 4 unplanted rows as «stable outliers».
  - Reason ①: the χ²_p 0.99 threshold is thinner than the tail of distances after small samples · per-device MAD standardisation · MCD. In a structure the same as the real one (devices 7 · sizes 30/20/20/20/10/10/9 · subjects 20 · p 4), the calibrated threshold is 29.5, more than twice χ²_4 0.99 = 13.3.
  - Reason ②: resampling subjects with replacement (bootstrap) duplicates rows and narrows the MCD covariance.
- **Fix (definition):**
  - Outlier threshold = the 0.99 quantile of the distance from running the same procedure (within-device robust standardisation → MCD) on normal null data with the same row · device · subject structure (Monte Carlo — sessions 200 times · seed 20261111 / ROI units 50 times · seed 20261113).
  - Stability = for each of 500 (ROI units 200) subject **subsamples (without replacement, 80 %)**, median · MAD · MCD are recomputed with the drawn subjects and every row is given a verdict. The threshold of that verdict is also the null 0.99 quantile of the same subsampling procedure (seeds 20261112 · 20261114). «Stable outlier» = outlier on the full data ∧ fraction flagged across subsamples ≥ 0.8.
  - The skeleton voxels of the registration indicator G are set as «skeleton ∩ individual brain mask (standard space · > 0.5)».
- **Synthetic test of the fixed version:**
  - Null 20 times: mean outliers 1.05 (expected 1.19) · mean stable outliers 0.05 · maximum 1.
  - Three planted outliers (+6 SD): all three stable outliers · others 0.
  - ② mapping functions: on synthetic data with shifts, the error of M1 · M2 is at noise level and that of M0 · M1g is the size of the shift. Inflation of the LOSO shift estimate 1.08 (expected √(n/(n−1)) = 1.04).
  - ① identification: separable synthetic data 1.00 · labels shuffled within subject 0.60.
- Scripts: d5_features.py 79f1d608… · d5_analysis.py d3c15e5d….

### 17-2. D5 ① – ③ results (addendum 2026-11-02 10:17 KST · the definition above is not changed · ④ after D1)

- **Run:** features eddy · roi 09:48 · maps 09:48:13 – 10:01:52 (8 processes) · ② 09:48 · ① 10:02 – 10:14 · ③ 10:02 – 10:05. Summary `explore/d5/d5.md` (output of d5_summary.py).
- **V-R: passed** — in all 119 sessions the recomputed roi.json «주» has maximum difference 0.0 from the stored values (bit-identical) · voxel counts the same.
- **① Device identification (balanced accuracy · chance reference 95 % 0.21 – 0.23):** T (tissue values 81) 0.841 · T_c 0.913 · E+T 0.898 · E (eddy 6) 0.585 · P̄+G 0.498.
  - T's errors concentrate in OXF1PRI ↔ OXF2PRI (same Prisma · only the coil differs · 4 · 7) and OXF3TRI (4/10).
  - Secondary: NOT2ING-89 vs 98 (same scanner · TE 89 · 98) 0.844 · chance reference 0.688 (n 20 · 89 = 4).
- **② Mapping prediction:** of pair × coordinate 69: «mapped at retest level» 2 (x · ODI of OXF1PRI~OXF2PRI) · «cannot be separated» 3 (OXF1PRI~OXF2PRI y · OXF1PRI~OXF3TRI ODI · OXF2PRI~OXF3TRI ODI) · «outside retest level» 64.
  - Median R: M0 3.84 – 4.23 → M1 2.34 – 3.28 · M2 2.33 – 2.77. M1 R_between median 0.86 – 1.15.
- **③ Outliers:**
  - Session Q (eddy quality): stable outliers 2 (sub-10975_ses-NOT2ING001 · sub-16745_ses-OXF4GEP001).
  - Session S (map summary · registration): stable outliers 0.
  - Check: the main-cell excluded session (sub-16794_ses-NOT4GEP001) is not an outlier in either Q or S (d² 3.2 · 6.7).
  - ROI units (session × ROI · 3,213 rows): outliers 1,399 · stable outliers 895 — tens of times the normal expectation (32). **Diagnosis after seeing the results (A):** rows with valid fraction exactly 1.0 are 28.6 %, so MAD 0 → the replacement scale becomes small, and there are rows with FWF percentile 90 equal to 0. So ROI-unit «stable outliers» capture distribution tails · boundary values rather than rare anomalies. It is left as the result under the definition, and this proviso is attached to the U-outlier of ④ (the definition is not changed · a fixed version only as a post hoc sensitivity analysis after an author decision).
  - Typical-class clustering: between-scanner GMM k 7 · ARI with device 0.386 / within-scanner GMM k 5 · sizes 24 · 45 · 1 · 39 · 10 · bootstrap ARI median 0.388 (unstable) · ARI with device 0.003.
- **Record fix:** the rule text of d5_3_outliers.json is out of date at «χ²_p 0.99 · bootstrap» (what the code actually did is the calibrated threshold · subsampling of §17-1). Only the text in the script was changed (d5_analysis.py 4ef68b6b… · result files unchanged · the sha256 of the version that ran is script_sha256 in the result files).
- **Self-test output kept:** `drafts/x14/D5/selftest/` (d5_selftest_outliers.py · .out — calibrated threshold 29.45 · null outliers 1.05 · stable 0.05 · the three planted all found).

### 17-3. Definition flaw found by review (axes 1 · 2 — reviews within the project · not external review) and its fix (addendum 2026-11-02 10:33 KST · the definition is not changed, only the reading is changed)

- **③ ROI-unit outliers — definition flaw (changes results):**
  - P6 (log number of valid voxels) = ln P5 (valid fraction) + ln(number of skeleton ∩ label voxels). In standard space the second term is the same across all sessions for each ROI (recounted by A: maximum spread within ROI 1.8e-15). The two features entered the same information twice.
  - Valid fraction = 1.0 in 28.6 %, and some (device × ROI) groups have MAD 0. On top of that the threshold was calibrated with a normal null.
  - → «ROI-unit stable outliers 895» is a value made by the procedure. **The U-outlier of ④ is not used (written as invalid).** ④ is given with S-outliers (session level) only.
  - Redoing it with a fixed definition (dropping P6 · handling boundary values) would come after seeing the results, so it is done only after an author decision.
- **Fix to the diagnosis of §17-2:** «there are rows with FWF percentile 90 equal to 0» is wrong. The minimum is 0.00026 (1 row below 1e-3). It was A's mistake of reading «0.000», printed to three decimals, as 0.
- **Reading proviso for ② (changes the interpretation):**
  - The retest scale is that of the one device OXF1PRI (SNR about 18). The SNR of the pairs is 8.9 – 19. And the bias from the SD estimate with degrees of freedom 10 (about × 1.1) was not included.
  - Two reviews tried an SNR adjustment after seeing the results. On axis 1 (the SNR scale of the judge's simulator), «outside retest level» fell x 22 → 4 · ODI 20 → 6, and only y stayed firm at 22 → 19. On axis 5 (SD ∝ 1/SNR · the most lenient assumption), R ÷ expected remained at x 2.0 · y 2.1 · ODI 1.6.
  - → Depending on the adjustment method, the conclusion for x · ODI shifts. **Only y is firm.** «The rest are 2 – 5 times the retest» is not written without this proviso.
- **Reading proviso for ①:** of the 17 errors in T identification, 15 came from between OXF devices with the same protocol. So «device fingerprint» is read as a «protocol + hardware» fingerprint (the two cannot be separated).

### 16-1. D1 reading proviso — noise predicted ratio (addendum 2026-11-02 10:33 KST · before the D1 computation · definition §2-7 unchanged)

- The predicted ratio 30/SNR_pair is the assumption «observed noise grows as 1/SNR, like thermal noise».
- The judge's simulator itself does not behave that way. On the D2 contour (N_eff 50), for SNR 30 → 20 the simulated SD is x × 1.79 · y × 1.01 · ODI × 1.25 (× 1.5 if thermal noise). y is dominated by voxel heterogeneity (PIX_SD).
- → In the D1 summary, 30/SNR_pair is described only as a **thermal-noise upper limit**. It is not read against the observed ratio under the name «prediction of the noise hypothesis». In particular it does not fit ODI (23 of 24 cells) and y.
- The D1 code is not changed (only the name · reading change).

### 16-2. Before the D1 summary — predictions of competing hypotheses · measuring the coil-rule asymmetry (addendum 2026-11-02 10:49:42 KST · before the D1 computation · review session reply 2026-11-02 «6-axis analysis reply» 3 · 2)

**Competing hypotheses and support · opposition in D1** ([A's proposal] · several can be true together · in D1 scanner and acquisition are entangled, so they cannot be fully separated):

| Hypothesis | Read as support | Read as opposition |
|---|---|---|
| H1 noise (observation noise · session variation assumed too small) | PIT mean ≈ reference · two tails ≫ reference · (a) z SD ≫ reference · within a cell, \|z\| larger for subjects with lower SNR_pair · ROI units outside the training range are worse | two tails unrelated to SNR · two tails ≈ reference even in low-SNR ROI units · bias (PIT mean) dominates |
| H2 location (measurement tissue in a sparse corner of the training distribution — low apparent FWF · low y) | within a cell, \|z\| · two tails larger at lower location y (rank correlation negative) · misspecification fraction substantially predicted by location (a) | location unrelated to \|z\| · two tails · test cases at the same location are well calibrated |
| H3 coil/L (32ch L4 · 64ch L8 noise model) | per-combination coverage systematically different in combinations containing 64ch · in the coil-only contrast (OXF1PRI~X vs OXF2PRI~X) the 64 side is worse even with SNR matched | per-combination coverage unrelated to combination · the differences in the contrast fully explained by SNR |
| H4 rule asymmetry (threshold with one true combination · measurement coverage with all candidates) | per-single-combination coverage ≫ coverage of the «all» rule · rule-matched threshold of e ≪ stored threshold · «two-tail share added by the rule» of the reference is large · difference between cells with one candidate (OXF1PRI~OXF2PRI) and with several | per-single-combination coverage also low · rule-matched threshold ≈ stored threshold · share added by the rule ≈ 0 |
| H5 isotropic T2 (effective T2 of the isotropic compartment shorter than the simulator's 1000 ms · more broadly, relaxation mismatch along the TE direction) | bias (z mean · PIT mean) moves with signed ΔTE (B − A) · cells with ΔTE = 0 (OXF1PRI~OXF2PRI) have no such share | bias unrelated to signed ΔTE · cells with ΔTE = 0 have the same bias |
| H6 scanner systematic (gradient nonlinearity · b error · reconstruction etc.) | bias large and (b) > expected · additive device model (cell z mean ≈ intercept + d_B − d_A) predicts cell z means well even with leave-one-cell-out · not explained by ΔTE · ΔTR · SNR | leave-one-cell-out R² of the additive device model ≈ 0 · acquisition covariates explain the bias |

- Biological variation · registration · spatially correlated noise cannot be separated from H1 in (a) of D1 (as in the §4 table).
- «Support · opposition» in the table above are descriptive directions and have no threshold. The cells are not independent of one another.

**D1 code changes (computation · random number consumption change 0 · only storage and summary added · d1_c3.py 30ab75e5…):**
- real: per pair, the per-combination coverage (covc · ROI unit × combination × k × component) and the ROI unit coordinates (X3) are also stored.
  - Reason: the per-combination coverage rate that the review session took to be «a few minutes from existing outputs» is not in the existing outputs. by_coil_combo of v4_c3.json has only per-combination misspecification counts · k-unmet counts. So it is given alongside in D1 real, which draws again by the same path as the verdict (tonight).
- Added to summary:
  - For each cell, the per-combination coverage rate (same ROI units used · one combination at a time · k median) and the coverage rate of the «all» rule are given side by side.
  - The Spearman of within-cell location (mean y of the two sessions) with \|z\| · two tails is given.
  - The between-cell Spearman of signed ΔTE (card dTE · B − A) with z mean · PIT mean is given.
  - Additive device model (ODI cells · card level 8 · with intercept · minimum-norm least squares): fit R² and leave-one-cell-out R² are given.
- ref: if another process is computing the same type (marker file), it waits until that finishes. This is because ref is launched first for e.
- New mode rulethr = **e (threshold with the rule matched)**:
  - For each cell, the measurement rule is applied to test cases of the same acquisition type with the true combination among the candidates (covered only if all candidates cover · nine k samples in all candidates · misspecification only if all candidates exceed the threshold).
  - With that rule, the threshold (percentile 5) and percentile 95 of the expected excess rate are computed. The random keys are the same as judge_v4.drawn_thresholds, and the number of clusters is the stored n_subjects.
  - It is written side by side with the value from the «one true combination» rule on the same samples · the stored verdict threshold · the measurement coverage · the measurement excess rate. The verdict is not changed.
- **Order deviation (recorded):** e has top priority, so the test-data side of D1 (ref · measurement data 0) is launched first, before the C7 commit. The measurement side (real · sens · summary) stays after the C7 commit.

## 18. Closing the exploratory list · analysis ledger (addendum 2026-11-02 10:51:38 KST · review session reply 2026-11-02 «6-axis analysis reply» 1)

- **Final list (not extended further):** D1 · D2 · a · b · c · d · per-combination coverage rate · e · (if the author approves in chat) f · 4a.
  - Planned only and not done: two-TE precision · sin²θ · g (background noise of ON-Harmony raw images · on hold) · h (D5 ③ again with the definition fixed · not done).
- **D5 new reads (maps 1,309 items · eddy text logs 595 items = 1,904 items):** the basis of approval relies on A's interpretation («is read as approving» · §17). Until the author sets directly whether it is approved, those results (D5 ① E · P̄+G · ③) are not used in the manuscript · talks. The decision is left in the deviation record. D5 ② and ① T use only re-reads of roi.json.
- **Analysis ledger:**

| Analysis | Status | What remains · reason |
|---|---|---|
| E1 / B1 (Wong direction 25.4°) | valid (observation) | the observation that the direction is outside the wedge · ψ (81 ROI units) and r 2.64 (66 ROI units) are different sets of ROI units (axis 1) |
| D2 (spread of same-scanner repeats) | valid · proviso | the y 2.2 times may be inflated by FWF level mismatch · logit amplification (axis 3) → c · ODI pooled 1.62 is because of the fornix (value without it 0.96) |
| D4 (forward table) | valid | algebra · no fit bias |
| B1-S | **conclusion invalid** | one-sided threshold on a circular angle · the basis for leaving out L = 1 was wrong · circularity of matching noise with ODI · the simulation inflates the x slope 3 – 5 times (axes 1 · 3) |
| D3 · D3 extension | conditions not met (no measurement) | remains as a record of truth-test failure |
| D3-v3 | **conclusion invalid** | truth-test noise coarser than the measurement and floor bias 5 – 30 times the fit scale, so the 3.84 interval does not hold · V-G relies on the T2_in upper bound · three definitions changed at once · T2_iso fixed (axes 1 · 3) |
| D5 ① device identification | valid · proviso · awaiting decision on new reads (E · P̄+G) | «device fingerprint» = protocol + hardware |
| D5 ② mapping prediction | valid · proviso | retest scale SNR mismatch · degrees-of-freedom bias → d · only y is firm |
| D5 ③ session outliers | valid · awaiting decision on new reads | Q stable outliers 2 · S 0 |
| D5 ③ ROI-unit outliers · ④ U-outlier | **invalid** | P6 = ln P5 + constant · point mass at valid fraction 1.0 · threshold calibrated under normality |
| D1 · per-combination coverage · e · a · b · c · d | in progress | §16 · §16-1 · §16-2 · §19 |

## 19. Review session approved analyses a · d (definition 2026-11-02 10:51:38 KST · before computing) — c · b are written separately after reading the code

- **a Location-matched expected misspecification excess rate [review session proposal · approved]:**
  - Measurement ROI units: the same collect as the verdict (roi.json «주» re-read) · 24 cells with data · main variant (main-cell excluded subjects left out). X3 = (x_A · y_A · ODI_A · x_B · y_B · ODI_B).
  - Test cases: x_c3 of the 1,000 cases of the same acquisition-type stratum (t|L0) and the misspecification flags used in computing the verdict threshold (test_cover_c3.npz «t_flag» · true condition).
  - Location matching: 6-dimensional Euclidean distance standardised by the SD of the stratum. For each measurement ROI unit, the misspecification fraction of the nearest k = 50 cases is taken as that ROI unit's location expectation (sensitivity k 25 · 100).
  - Cell expected excess rate = mean of the location expectations of the cell's ROI units. Observed = stored exceed_rate.
  - **Definition of «what % of the observed»:** per-cell ratio = expected ÷ observed → between-cell median (main). In addition, pooled ratio = Σ(cell expected × number of ROI units) ÷ Σ(observed × number of ROI units).
  - Limits: the test flags are under the true-coil condition and the measurement uses the candidate-combination minimum rule (entangled with H4). The distance does not match acquisition conditions (TE · TR · coil · SNR) (stratum = acquisition type only).
- **d D5 ② again with an SNR-scaled retest [review session proposal · approved · marked «after seeing the results»]:**
  - The ROI units · mappings (M0 · M1) · LOSO prediction errors of D5 ② are used as they are.
  - Degrees-of-freedom change: R is multiplied by √((ν − 2)/ν) (ν = pooled degrees of freedom 10 per ROI · E[1/s²] = ν/((ν − 2)σ²)).
  - SNR change: per subject · ROI σ_d = σ_w,r · √((SNR_ref/SNR_A)^(2α) + (SNR_ref/SNR_B)^(2α)). SNR_ref = median of the SNR of the 12 OXF1PRI repeat sessions (snr_table).
    - α = 1 (thermal noise · main).
    - α_c = simulator scale (from the simulated SD ratio for D2 contour N_eff 50 · SNR 30 → 20, α_c = ln(ratio)/ln 1.5 · x 1.79 · y 1.01 · ODI 1.25).
  - The reading is the same as ② (M1 R 95 % upper ≤ 1.5 × expected → mapped at retest level · lower > 1.5 × expected → outside · otherwise cannot be separated).

### 19-1. a · d results (addendum 2026-11-02 10:55:57 KST · descriptive · not a verdict · `explore/rv/`)

- **a (location-matched expected misspecification excess rate):**
  - Cell expected (k 50) 0.009 – 0.121 · observed 0.163 – 0.774.
  - Ratio median **0.17** (k 25 0.24 · k 100 0.12) · pooled ratio 0.17.
  - → Under the definition, location (6-dimensional X3 · same acquisition type) explains only a small share of the misspecification excess. The «about 60 % of the observed» of review axis 3 is not reproduced under this definition (the definition of that computation is in the axis 3 report and differs from this one).
- **d (SNR-scaled retest · after seeing the results):** SNR_ref 18.56 · α_c = x 1.43 · y 0.035 · ODI 0.54.
  - M1's R ÷ expected median (original → α 1 → α_c): x 3.07 → 1.75 → 1.42 · y 3.06 → 1.72 → 2.70 · ODI 2.27 → 1.40 → 1.60.
  - Reading counts (23 pairs · mapped / cannot be separated / outside):
    - α 1: x 1/16/6 · y 1/15/7 · ODI 11/12/0
    - α_c: x 6/13/4 · y 1/2/20 · ODI 3/15/5
  - For OXF1PRI~OXF2PRI, all three coordinates are «mapped» under both scales.
  - → «2 – 5 times the retest» mostly disappears depending on the retest scale. Under the simulator scale, only y is firmly outside.

## 19-2. Review session approved analyses b · c (definition 2026-11-02 10:55:57 KST · before computing)

- **Observed target (B1 · 81 ROI units · given before computing):** slope median b_x 7.939 · b_y 3.424 s⁻¹. The 95 % interval of a ROI-unit bootstrap 2,000 times (seed 1) is b_x [7.293, 8.802] · b_y [2.303, 4.187]. The algebraic slopes of the literature 82/52/1000 are 7.036 · 11.195.
  - Check of the review session's arithmetic (A): the (T2_in, T2_en, T2_iso) matching the two algebraic slopes form a one-parameter family. T2_in 82 → T2_iso 118.7 · T2_iso 1000 → T2_in 209.6 · T2_en 78.6 (r 2.67) · T2_iso 2000 → T2_in 234.2 · T2_en 81.8 — same as the review session's values.
- **b — solution family of 25° (noise 0 · slopes after the fit) [review session revision applied]:**
  - Grid:
    - T2_in {60, 70, 82, 100, 125, 150, 180, 210, 250} ms
    - r = T2_in/T2_en {1.0, 1.3, 1.6, 2.0, 2.4, 2.8, 3.2}
    - T2_iso {80, 100, 119, 150, 200, 300, 500, 1000, 2000} ms
    - true d∥ {1.7, 2.0, 2.26} (the fit is fixed at 1.7)
    - → 1,701 conditions.
  - Tissue: the 24 tissues of B1-S (ν {0.45, 0.6, 0.75} × ODI {0.08, 0.20} × FISO {0, 0.02, 0.05, 0.10}) · one homogeneous voxel per tissue · fibre direction = tissue number mod number of directions (same as B1-S).
  - Signal: the same equation as B1-S pix_signal, with only T2_iso and true d∥ opened as parameters (TR 5800 · T1 (1100, 1100, T1_iso) · Wong sub-01 TE92 direction scheme).
    - V-S′: it runs only if the difference from pix_signal is ≤ 1e-12 at T2_iso 1000 · d∥ 1.7.
  - Fit: simulate.amico_fit (F2 · d∥ 1.7 · d_iso 3.0) · eight TE · lib.coords (C3). The slope is the least-squares linear slope against TE (same as the slope of fit_unit).
  - For each condition, the median of the slopes of the 24 tissues is used.
  - **Match = b_x · b_y both inside the observed 95 % interval.** Outputs: list of matching conditions · a (T2_iso, T2_in) plane plot per d∥ (matching points · r marked · algebraic family curve overlaid).
  - Check (does not change results): the angle of the literature 82/52/1000 · d∥ 1.7 is compared with the 57.4° of B1-S at noise 0.
  - **Reading: only as far as «the family of solutions matching the observation».** Separating within the family is the job of f (D3 with the floor included).
- **c — D2 re-comparison with the FWF level matched [review session approval]:**
  - The same setup as D2 (plan of d2_sim · chunk · roi_noise capture · V1) is used, changing only the centre.
  - Changed centre: the true (ν, f_iso) matching the apparent (NDI, FWF) of the observed first session is solved by inverting config.apparent.
    - OXF1PRI acquisition (TE 92 · TR 3600) · relaxation at the L0 prior centre (T2 (82, 82/1.6, 1000) · T1 (1100, 1100, T1_iso)) are used.
    - ODI as it is · clipped with PIX_CLIP (same as D2).
  - The plan's random numbers · row structure · seeds are the same as D2 (common random numbers — the same cases except for the centre). For SNR scenarios, only the chunks of the two sets snr30 · measured are run (chunks 0 – 26 · 81 – 107).
  - In addition, the difference in FWF itself (absolute scale) is captured.
  - Comparison: ROI median · pooled · N_eff contour of the observed ÷ simulated SD ratio (x · y · ODI · absolute FWF) are given side by side with the original D2. Level-matching check = simulated apparent FWF median vs observed FWF.
  - Reading (before computing · [A's proposal]): if the level-matched pooled y ratio (SNR 30 · N_eff 50) is
    - ≤ 1.2 → «the y excess is explained by level mismatch · logit amplification»
    - ≥ 2.0 → «the y excess remains even with the level matched»
    - in between → «partial»

### 19-3. b results (addendum 2026-11-02 11:01:27 KST · descriptive · not a verdict · `explore/rv/rv_b_family.json` · figure rv_b_family.png)

- **Run:** 10:57:21 – 10:57:44. The first run (10:57:03) exceeded the NIfTI axis size limit (32,767) and stopped before fitting (results 0). It was run again with only the image shape split into (condition, tissue). Script sha256 8a100734….
  - The figure was drawn in the main environment from the same JSON, because the reference environment has no matplotlib.
- V-S′ 0.0. Check: the post-fit angle of the literature 82/52/1000 · d∥ 1.7 is 55.2° (B1-S noise 0 57.4° · algebra 57.9°).
- **Matching conditions 27 / 1,701** (d∥ 1.7: 11 · 2.0: 5 · 2.26: 11).
  - d∥ 1.7 (same as the fit): T2_iso 100 – 2000 · T2_in 70 – 250 · r 1.6 – 3.2. If T2_iso ≥ 1000, T2_in 180 – 250 · r 2.4 – 3.2 match. The literature level (T2_in 70 · r 1.6) matches only at T2_iso 100.
  - d∥ 2.0: T2_iso 300 – 2000 · T2_in 82 – 125 · r 2.0 – 2.4.
  - d∥ 2.26 (larger than the fit): T2_iso 150 – 2000 · T2_in 60 – 100 · r 2.0 – 3.2 (e.g. T2_in 82 · r 2.4 – 2.8 · T2_iso 1000).
  - In this grid, no condition with «r ≤ 1.6 and T2_iso ≥ 300» matches. The simulator prior (r 1.0 – 1.8 · T2_iso 1000) is outside the family.
- **Reading (as defined): it is the family of solutions matching the observation.** It is not separated within (the job of f). The two slopes alone cannot separate «isotropic T2 is short» from «T2 ratio is large». If true d∥ differs from the 1.7 of the fit, the family moves.

### 19-4. c results — the premise of level matching does not hold (addendum 2026-11-02 11:39:21 KST · descriptive · not a verdict · `explore/rv/rv_c_level.json`)

- **Run:** 10:59:51 – 11:36:32 · chunks 54 · V1 passed (wrapped result bit-identical) · V2 all true · inversion check 1.1e-16 · centre clips 0 (D2 had fornix 2). Script sha256 d7f98822….
- Observed apparent FWF median 0.071 · the true f_iso matching it median 0.034 (NDI 0.640 → true ν 0.475).
- **Level-matching check (the check line of the definition): failed.** The simulated apparent FWF median is snr30 0.230 · measured 0.257, a little over three times the observed 0.071 (median absolute difference 0.139 · 0.174).
  - → This computation did not become «a level-matched comparison». Even with the relaxation weighting inverted, the simulated apparent FWF does not come down. The root of the level difference lies not in the relaxation weighting but elsewhere (clipping of the voxel distribution · noise floor · fit bias — hypotheses · not separated).
- Numbers (as defined · read only with the proviso above) — snr30 · N_eff 50 · observed ÷ simulated · median / pooled:
  - x 0.72 / 1.00 · y 1.45 / 1.94 · ODI 0.58 / 1.52 · FWF absolute 0.69 / 1.51
  - The original D2 is x 0.55 / 0.78 · y 2.20 / 2.60 · ODI 0.58 / 1.62
  - measured: x 0.34 / 0.46 · y 1.73 / 2.24 · ODI 0.40 / 1.11 · FWF absolute 0.58 / 1.28
  - Under the reading line of the definition (y pooled 1.94) it is «partial». But since the level was not matched, it does not answer «is the y excess explained by level mismatch · logit amplification».
- **Added descriptive account (after seeing the results · new reads 0 · frozen test data):**
  - The site A apparent FWF of test data L0 5,000 cases (from x_c3, FWF = expit(y + ln NDI)) is median 0.419 · percentile 5 0.158 · percentile 1 0.075. Cases ≤ 0.10 are 2.0 %, ≤ 0.08 are 1.2 %.
  - Measurement OXF1PRI white-matter ROIs have median 0.071 and 76 % are ≤ 0.10.
  - → Real white matter lies in the 1 – 2 % tail of the simulator distribution. From c above, even changing the tissue values does not reach that tail at this SNR. It is a structural gap.

### 16-3. One more competing hypothesis — H1′ noise assumed too large (addendum 2026-11-02 11:39:21 KST · before the D1 computation · added to the §16-2 table)

- **Reason (collected observations after seeing results · not a test):**
  - In D2, the real repeat spread of x · ODI is smaller than under the SNR 30 assumption (0.55 · 0.58).
  - Wong's measured SNR is flat over TE 62 → 92. It is not a thermal-noise SNR.
  - In c, the simulator's apparent FWF is higher than observed (consistent with noise-floor bias).
  - If the measurement «SNR» (b0 temporal variance) overstates the effective noise, the simulator's SNR 30 noise may be coarser than real.
- **H1′ support:** two tails ≤ reference (less confident) or similar to the reference · PIT mean far from the reference (failure is bias) · within a cell, no relation between SNR_pair and \|z\|, or in the opposite direction.
- **H1′ opposition:** two tails ≫ reference · \|z\| larger for lower-SNR ROI units (the H1 side).
- H1 and H1′ are opposites. They are separated by the two tails · (a) z SD vs reference · SNR relation of D1.

### 16-4. One more competing hypothesis — H7 mismatch of noise distribution shape (effective L · noise floor) [review session proposal · operationalised by A] (addendum 2026-11-02 12:27:49 KST · before the D1 computation)

- **Hypothesis:**
  - The simulator uses a non-central chi with 32ch → L 4 · 64ch → L 8 (decomp/s4/s4cfg.py:51 · decomp/s4/simulate.py:61 noisy · at :126 L is chosen from the cell's coil).
  - If the real reconstruction is of the adaptive combine · SENSE family, effective L ≈ 1 (close to Rician), so the noise floor is much lower.
  - Our processing (run.py · cmds.py) has no denoising step and no Gibbs-ringing step (A searched the code · 0 items). However, the resampling (interpolation) of eddy changes the noise distribution.
- **Operationalisation (the review session's prediction changed so that it can be measured in D1):**
  - The 24 cells with data all have one or more candidate combinations containing 64ch. So «cells with a 64ch candidate vs cells without» has no contrast.
  - ROI units within a cell share the same candidate combinations, so «coil combination and \|z\| within a cell» does not vary between ROI units.
  - → **Combinations are compared with one another within the same ROI unit:** PIT bias · \|z\| · coverage · k filling (D1 stores them per combination) of combinations containing L 8 (with a 64ch side) and combinations with L 4 only.
- **Support:** within the same ROI unit, L 8 combinations have PIT means farther from the reference, larger \|z\| and lower coverage than L 4 combinations. In the per-combination coverage rates L 8 combinations are worse. That difference is larger than the between-ROI-unit difference related to SNR_pair.
- **Opposition:** no difference between combinations, or L 4 combinations are worse.
- **What cannot be separated:** «the real coil is 32, so the 64 combinations are wrong» and «effective L is smaller than both» cannot be separated by D1 alone (needs the coil name · combination method fields of 4a).
- **Direction of the already stored values (after seeing the results · descriptive · v4_c3.json by_coil_combo · re-read record):** in the 10 cells where site A has an unknown coil, cells where the misspecification count of A = 64 combinations exceeds that of A = 32 combinations 9/10 · cells with more k-unmet 7/10.

### 16-3a. Proviso on the basis of §16-3 (addendum 2026-11-02 12:27:49 KST · review session reply 3c)

- «Wong's measured SNR is flat with TE» of §16-3 may be a phenomenon arising because Wong is data denoised by the provider. Then it does not carry over to ON-Harmony. The weight of this line in the basis of H1′ is lowered.

### 18-1. Tidying the record (addendum 2026-11-02 12:27:49 KST · review session reply 1)

- **Number of D5 new reads:** 1,904 items = maps · masks · transforms 1,309 items (sessions 119 × 11: f2_NDI · FWF · ODI · MD · FA · s_low · s_2000 · nodif_brain_mask · f2_extvalid · f2_f2valid · fa2std_warp) + eddy text logs 595 items (sessions 119 × 5: movement_rms · outlier_map · parameters · post_eddy_shell_alignment_parameters · bvals).
  - The 1,309 of the last relay counted only the map side. This 1,904 adds the eddy text logs. It is a difference in scope, not a discrepancy.
  - All are our processing outputs (data/remeasure/onharmony/<session>/). Sealed raw images (source_data/onharmony_sealed) · verdict/ · blinded data are not among them (the paths of ONH_NEW_READ_LOG.jsonl compared by machine).
  - Sessions 119 = first sessions 109 + OXF1PRI repeats 10.
- **The share explained by location:** the «about 60 % of the observed» of the review axis 3 is a **withdrawn value**. The result of the defined analysis a (§19 · §19-1) is 0.17 (k 25 0.24 · k 100 0.12).

### 19-5. Scope of c · generalisation on the observed side (addendum 2026-11-02 12:27:49 KST · review session reply 2 · new reads 0)

- **Limitation:** «2 % vs 76 %» is limited to «under the noise settings of the frozen simulator». «Even changing the tissue values does not reach that tail» was confirmed only for the two noise scenarios (snr30 · measured) and the inversely solved centres (54 · true f_iso median 0.034).
- **2b — all V4 first sessions (roi.json re-read · ROI 27):**

| Scanner | Sessions | Apparent FWF median | Fraction ≤ 0.10 |
|---|---|---|---|
| NOT1ACH | 20 | 0.074 | 0.72 |
| NOT2ING | 20 | 0.106 | 0.47 |
| NOT4GEP | 10 | 0.081 | 0.79 |
| OXF1PRI | 20 | 0.074 | 0.72 |
| OXF2PRI | 20 | 0.076 | 0.70 |
| OXF3TRI | 10 | 0.075 | 0.68 |
| OXF4GEP | 9 | 0.083 | 0.71 |
| All | 109 | 0.080 | 0.67 |

  - Simulator test data L0 has median 0.419 · ≤ 0.10 is 0.020. → On all seven scanners real white matter lies in the tail of the simulator distribution.

### 18-2. Deviation record — approval of the D5 new reads (addendum 2026-11-02 13:32:44 KST · [author decision] M1 decision 1 · 2026-11-02 13:31:06 KST)

- **Maps · masks · transforms 1,309 items: retrospective approval, option (A).** The approval relied on A's interpretation, and the author confirmed that interpretation at 2026-11-02 13:31:06 KST. The author: «09:35 'let's proceed, please design it rigorously' was an approval of reading the maps of 119 sessions.»
- **eddy text logs 595 items:** the author's rule is «(A) if the 09:35 question included the eddy records, (B) if not». Facts A confirmed from the chat record:
  - At 09:27:46 (line 10519), the choices under «for you to set · scope» were «(A) eddy only / (B) eddy and maps, all 119 sessions ← recommended / (C) eddy and maps, OXF1PRI 12 sessions only».
  - At 09:30:34 (line 10526), the three analyses proposed were «① device identification: eddy features only · tissue values only» · «③ within-scanner outliers: eddy and maps».
  - At 09:34:13 (line 10541), the timetable had «① device identification — organising 119 eddy records».
  - Proviso: in «for you to set» of line 10526, the only thing named as «whether you approve» was «reading the maps (no. 3) for all 119 sessions».
  - → Because the question included the eddy records, it is recorded, following the author's rule, as **option (A), retrospective approval**. The approval relied on A's interpretation and was set by the author's rule (13:31:06) and the facts above. If the author chooses option (B) on the proviso above, this line is changed and an addition is made.
- Location of the record: this section · the new-read work record (ONH_NEW_READ_LOG.jsonl) · together with a deviation mark when moved to the canonical access record in the next record commit.

## 20. 4a — three ON-Harmony sidecar coil fields [author decision M1 decision 2 · approved 13:31:06] (definition 2026-11-02 13:33:18 KST · before reading)

- **What is read:** from the dwi sidecars `<ONH_META>/<sub>/<ses>/dwi/<sub>_<ses>_dir-AP_dwi.json` of the 109 V4 first sessions (the sessions of the attribute table) (the very files from which TE · TR were read at the verdict), only three fields are extracted: ReceiveCoilName · CoilCombinationMethod · ParallelReductionFactorInPlane. Other fields are neither extracted nor written. If a field is absent, «no field».
- **Record:** per file, time · path · sha256 · bytes in ONH_NEW_READ_LOG.jsonl (what = «4a sidecar three fields») · through the guard (guard.require).
- **Interpretation of the coil channel number (set before reading):** the first of {8, 12, 16, 20, 32, 48, 64} occurring among the integers in ReceiveCoilName is taken as the channel number. If none, «unknown (name as it is)». If not 32 · 64, «other» (outside the simulator's options).
- **Reading by result (per scanner · not a verdict · continues from plan §16-4 H7 · §16-2 H4):**
  - ① Receive coil is 32 channels → on that scanner the explanation «the 64 combinations were wrong (actually 32)» holds. That A = 64 combinations were worse in cells containing that scanner (stored value 9/10) is read as «wrong coil candidate».
  - 64 channels → the 64 combination is the right combination. The explanation «actual coil 32» does not hold on that scanner.
  - ② Combination method of the Adaptive Combine family (adaptive · ACC · «AC» combination mark in the name) → towards effective L ≈ 1 (low noise floor), and H7 is supported on that scanner.
  - Sum of Squares family → towards non-central chi (effective L ≈ effective number of channels), and H7 weakens on that scanner.
  - Both ① and ② → the two explanations remain together. Since the true combination is known, «the coverage rate of the true combination» can be given.
  - Field absent → that scanner cannot be separated (written as «no record»).
  - The parallel imaging factor is described only (spatially varying noise · g factor — not used for judgement).
- **To be given in addition (after D1):** in cells where the channel numbers of both scanners are set as 32 · 64, the coverage rate of the true combination among the D1 per-combination coverage rates is written side by side with the «all» rule coverage rate · stored threshold · rule-matched threshold of e.
- Output: `explore/rv/rv_4a_coil.json` (three fields per session · per-scanner summary).

### 20-1. 4a results (addendum 2026-11-02 13:34:30 KST · descriptive · not a verdict · `explore/rv/rv_4a_coil.json` · script sha256 96826c24… · 13:33:41 · sidecars 109 files · new-read record 109 lines)

| Scanner | Sessions | ReceiveCoilName (channels) | CoilCombinationMethod | ParallelReductionFactorInPlane |
|---|---|---|---|---|
| OXF1PRI | 20 | Head_32 (32) | Sum of Squares | no field |
| OXF2PRI | 20 | HeadNeck_64 (64) | Sum of Squares | no field |
| OXF3TRI | 10 | 32Ch_Head (32) | Sum of Squares | no field |
| NOT1ACH | 20 | no field | no field | 2 |
| NOT2ING | 20 | no field | no field | protocol 89: 2 (4) · protocol 98: 1.5 (16) |
| NOT4GEP | 10 | no field | no field | 2 |
| OXF4GEP | 9 | no field | no field | 2 |

- **Reading as defined (§20):**
  - OXF3TRI = 32 channels → ①: on OXF3TRI «the 64 combinations were wrong (actually 32)» holds. That A = 64 combinations were worse in the three cells where OXF3TRI is site A (stored values) is read as «wrong coil candidate».
  - OXF1PRI 32 · OXF2PRI 64 are the same as the judge's coil record.
  - The combination method of the three Siemens scanners is all Sum of Squares → H7 (effective L ≈ 1 through Adaptive Combine) weakens on Siemens scanners. The non-central chi form is right, and the remaining question is the size of effective L.
  - The four Philips · GE scanners (NOT1ACH · NOT2ING · NOT4GEP · OXF4GEP) lack the coil · combination method fields among the three → cannot be separated («no record»).
  - Cells whose true combination is determined (after D1, «coverage rate of the true combination»): OXF1PRI~OXF2PRI (32|64 · also one originally) · OXF1PRI~OXF3TRI (32|32) · OXF2PRI~OXF3TRI (64|32).
- **What also came to light — record change:** the two NOT2ING protocols differ not only in TE (89 · 98) but also in the in-plane parallel imaging factor (89 = 2 · 98 = 1.5). So «same scanner · only TE differs» in the D5 ① secondary analysis and §17-2 · §17-3 was wrong. It is changed to «same scanner · TE and parallel imaging factor differ». The description of the NOT2ING protocols in manuscript §3.x (mentioning only TE) is also changed when work resumes (manuscript note).

## 21. f — Wong raw data [author decision M1 decision 3 · approved 13:31:06 · starts after seeing the e · D1 results] (addendum 2026-11-02 13:34:30 KST · reads 0)

- **Step 1:** from the headers only, check (i) whether it is data before denoising (ii) whether the background (outside the brain) remains. If either is not the case, stop and report.
- **Step 2:** only the minimum files needed for estimating σ and effective L are read (which files · how many is written after step 1, before computing).
- **Step 3:** the definition of the fit with the floor included is written before computing → truth test → measurement.
- Each read is recorded in EXTERNAL_READ_LOG.jsonl. The license is checked from the license file of the dataset and reported together. GPU 0.

### 20-2. Change to the 4a reading · one line in the gap ledger · scope of f · closing branches (addendum 2026-11-02 14:59:33 KST · review session «reply on applying M1»)

- **Change to the 4a reading [review session fix]:** the 3 Siemens scanners (OXF1PRI · OXF2PRI · OXF3TRI) are Sum of Squares, so the explanation «Adaptive Combine → low noise floor» is **withdrawn** for Siemens.
  - The remaining question is the size of effective L. If there is correlation between channels, L_eff is smaller than the number of channels.
  - The four Philips · GE scanners (NOT1ACH · NOT2ING · NOT4GEP · OXF4GEP) have no coil · combination fields and are left **undecided**.
  - That it might be a SENSE-family reconstruction, judging from the parallel imaging factor 2, is written only as a **hypothesis** (inference · not confirmed).
- **One line in the gap ledger [review session proposal]:** the differences in parallel imaging factor (NOT2ING-89 2 vs 98 1.5 · Philips · GE 2 · Siemens field absent) and the resulting spatial noise amplification (g factor) are neither in the simulator nor in the attribute columns of Fig. 4.
- **Proviso on the scope of f (added to §21) [review session proposal]:** f only measures σ and L_eff of the Wong scanner. It is not carried over to the other manufacturers of ON-Harmony (Philips · GE).
- **Closing branches [review session proposal]:** after tonight's D1 overall report, no new branch other than f is opened. Whether to do the sufficiency test for reproducing the failure is set by the author and the review session after seeing that report.

### 16-5. e results — threshold with the rule matched (addendum 2026-11-02 17:49:42 KST · descriptive · not a verdict · the verdict stays as it is #34 · `explore/d1/d1_rulethr.json` sha256 252de104… · 17:37:46 – 17:37:47 · D1 ref four types 10:49 – 17:36)

- **Effect of the rule on coverage in the test data (same samples · the cell's acquisition type · true combination among the candidates):** the coverage under the one-true-combination rule (k median) is 0.88 – 0.93. Applying the measurement rule (covered only if all candidate combinations cover) it falls as follows.
  - Cells with 2 candidates: 0.28 – 0.37
  - Cells with 4 candidates: 0.075 – 0.12
  - Cells with 1 candidate (OXF1PRI~OXF2PRI): unchanged (0.90 · 0.91)
- **Threshold with the rule matched (percentile 5):** cells with 2 candidates 0.10 – 0.20 · cells with 4 candidates 0.000. The median difference from the stored threshold (0.50 – 0.85) is 0.67.
- **Compared with the measurement coverage (descriptive):**
  - Cells at or above the rule-matched threshold 21/24. Below: OXF1PRI~OXF2PRI|NDI (0.010 vs 0.80) · OXF1PRI~OXF2PRI|ODI (0.649 vs 0.75) · OXF2PRI~NOT2ING-98|ODI (0.102 vs 0.125).
  - The rule-matched threshold of the cells with 4 candidates is 0, so «exceeding» carries no information in those cells.
- **Percentile 95 of the expected excess rate (misspecification):** almost the same even with the rule matched (the misspecification rule is the candidate minimum and so lenient). Cells where the measurement excess rate exceeds it stay at 23/24.
- **A's reading (after seeing the results · not a verdict):**
  - ① Coverage test: in the 22 cells containing a scanner with no coil record, the threshold set in advance could not be exceeded even by the simulator's own data (with the rule matched, test-data coverage is 0.08 – 0.37). So «failure» in those cells carries little information about the judge. The only calibrated coverage test is OXF1PRI~OXF2PRI with one candidate, and there NDI 0.010 · ODI 0.649 fall short of the threshold (0.80 · 0.75). OXF1PRI~OXF3TRI · OXF2PRI~OXF3TRI, whose true combination was determined by 4a, are added tonight with D1's true-combination coverage rate.
  - ② Gap detection (misspecification excess): robust to the rule (23/24).
  - ③ In the test data, conditioning on a wrong coil candidate drops coverage 0.9 → 0.1 – 0.4. The judge's posterior relies heavily on the coil (→ noise floor L) assumption — the weight of H7 grows.

## 22. Handling the reply on e [review session «reply on e results»] (addendum 2026-11-02 18:02:26 KST · 2c · 4 are definitions before computing · 3 is before the D1 computation)

- **Record:** non-reading batch 68 commit = **53fba51** (2026-11-02 09:50:53 KST · author · push · one changed file PLAN.md).
- **1a analysis ledger (added to the §18 table) — a flaw in the fixed design (confirmed after seeing the results):** in the 22 cells containing a scanner with no coil record, the pass threshold set in advance could not be exceeded even by the simulator's own data. The threshold was computed with the «one true combination» rule (judge_v4.test_cover · cond = true condition), while the measurement coverage was measured with the «all candidate combinations» rule (judge_v4.cover_by_k). In e (§16-5), applying the measurement rule to the test data gave coverage 0.08 – 0.37. In the manuscript it goes at the start of the V4 results and in the abstract (after work resumes).
- **1b:** «21/24 exceeded» is not counted as a number of passes. For each cell, the measurement coverage · rule-matched threshold · original threshold are placed side by side. Cells with 4 candidates are marked «no information (threshold 0)», and cells with 2 candidates «weak (threshold 0.10 – 0.20)». The table is `explore/rv/rv_e2.json · .md`.
- **2 Symmetric cells (marked as after seeing the results):**
  - (a) By 4a, OXF3TRI = 32 channels. So OXF1PRI~OXF3TRI (32|32) · OXF2PRI~OXF3TRI (64|32) have one candidate combination. Together with OXF1PRI~OXF2PRI (32|64 · NDI · ODI) they are reported separately as «symmetric cells».
  - (b) When the original threshold was computed, the test data were an acquisition-type stratum (t|L0), not a scanner. Each case was conditioned on that case's true coil combination (four combinations mixed · about 25 % each) (judge_v4.test_cover · cond = sim["cond"]). OXF3TRI was not separately set to any coil — it is the «one true combination» coverage of test cases with four combinations mixed.
  - (c) Definition: from our ref samples (same stratum), using only cases whose true combination is the actual combination of that cell (OXF1PRI~OXF3TRI 32|32 · OXF2PRI~OXF3TRI 64|32 · OXF1PRI~OXF2PRI 32|64), the threshold (percentile 5) is computed again from the coverage · misspecification under that combination's condition. The random keys · number of clusters are the same as judge_v4.drawn_thresholds. The measurement-side «coverage rate of that combination» is compared through the per-combination coverage of D1 real (tonight).
- **4 Definition of coil assumption sensitivity (frozen test data = ref samples · new reads 0):**
  - For each stratum, each case i (true combination t) and each assumed combination a (four), using cases with samples full for all nine k, coverage (mean per k → k median) · mean of z = (posterior mean − true value)/posterior SD (k median) · PIT mean (k median) are given per NDI · ODI component.
  - Two summaries:
    - Per side (only one side wrong) — site B true → assumed (32→32 · 32→64 · 64→32 · 64→64) with A the same as true. Site A in the same way.
    - Values with the four strata combined and per-stratum values.
  - Direction asymmetry = comparing the coverage · sign and size of z mean of 32→64 and 64→32.
- **3 Prediction before D1 — OXF1PRI~OXF2PRI (symmetric cell differing only in coil · coil record right):**
  - The sign of the NDI · ODI z mean (k median) of this cell in D1 is compared with the sign of the z mean «in the match stratum · when the site B coil is wrongly assumed» from 4. The comparison is in both directions: larger than actual (32 → 64) · smaller than actual (64 → 32).
  - **H7 support:** the signs of NDI · ODI both agree with one direction (larger or smaller) → consistent with «the effective noise floor of site B differs from the assumption (in that direction)».
  - **H7 opposition:** the sign agrees with neither direction, or |z mean| < 0.1 (no bias).
  - The direction on the site A side is written for reference in the same way.
- **5 Robustness of gap detection:** the misspecification exclusion and the cell excess rate require «all candidate combinations to exceed the threshold» (minimum score · lenient). The injection test (detection 0.010 – 0.062) and the test-data expected excess rate were measured with one true coil (validate.run_inject:708 · cond = true condition). → The direction of the rule is to lower the measurement excess rate. The «outside the distribution» of the measurement 23/24 did not arise because of the rule.
- **6 C7:** at 26 hours 39 minutes · state R · CPU 1,082 % · cumulative CPU about 288 hours. There is no progress log (stdout 0 bytes · stderr only a start warning), so progress · remaining time cannot be known (no py-spy · the process is not touched). D1 does not depend on C7 results as data (C3 arm · v4_c3.json · test data only). Whether D1 starts first is for the author to set.

### 22-1. Results of the reply on e — 1b · 2c · 4 (addendum 2026-11-02 18:03:45 KST · descriptive · not a verdict · `explore/rv/rv_e2.json · .md` sha256 79e2183b… · script 755cb639… · 18:03)

- **1b:** the per-cell table is in rv_e2.md. The marks are calibrated (candidate 1) 2 cells · weak (candidate 2) 8 cells · no information (candidate 4) 14 cells. To the two candidate-2 cells containing OXF3TRI, «symmetric cell after 4a» was added.
- **2c Symmetric cells (after seeing the results):**

| Cell | True combination | Threshold of that combination only | Original threshold | Test coverage | Measurement coverage (all-candidates rule) |
|---|---|---|---|---|---|
| OXF1PRI~OXF2PRI · NDI | 32\|64 | 0.800 | 0.800 | 0.912 | 0.010 |
| OXF1PRI~OXF2PRI · ODI | 32\|64 | 0.750 | 0.800 | 0.895 | 0.649 |
| OXF1PRI~OXF3TRI · ODI | 32\|32 | 0.800 | 0.800 | 0.938 | 0.568 |
| OXF2PRI~OXF3TRI · ODI | 64\|32 | 0.700 | 0.800 | 0.900 | 0.352 |

  - The measurement «true-combination coverage rate» of the two OXF3TRI cells is added tonight from the per-combination coverage of D1 real. That value is at or above the «all candidates» coverage (0.568 · 0.352).
- **4 Coil assumption sensitivity (four strata combined · per side · only one side wrong):**
  - With the right assumption, coverage 0.90 – 0.91 · z mean ≈ 0 · PIT mean ≈ 0.49.
  - Site B assumed larger (32 → 64): NDI coverage 0.31 · z −3.20 / ODI 0.37 · z −2.56.
  - Site B assumed smaller (64 → 32): NDI 0.34 · z +3.02 / ODI 0.33 · z +2.76.
  - Site A has the opposite sign (Δ = B − A): 32 → 64 is z +3.33 · +2.82, 64 → 32 is z −3.11 · −2.78.
  - **The direction asymmetry is small** — the sizes are similar at |z| 2.6 – 3.3 and only the signs are opposite. In the match stratum, NDI's B 32 → 64 is −3.62 · 64 → 32 is +2.99, slightly asymmetric.
  - → In the frozen judge, wrongly assuming one coil drops coverage 0.9 → around 0.3 and moves the posterior by about 3 SD. The prediction of §22-3 (z sign of OXF1PRI~OXF2PRI) is read by this sign: negative → «the same form as assuming B larger (B's effective floor is lower than the L 8 assumption) or assuming A smaller» · positive → the opposite.

### 22-2. Proposal outside the list — recorded only (addendum 2026-11-02 18:03:45 KST · review session «reply on e results» 8 · not done now)

- Measuring the per-coil noise floor from the b2000 ventricular signal of the OXF1PRI · OXF2PRI processed DWI. It is opened with author approval only if D1 points to H7 (new read · the processed DWI may not be in the value area — check before opening).

### 22-3. Addendum before the D1 computation [review session «reply on handling results» 2 · [author decision] D1 starts now] (addendum 2026-11-02 18:08:15 KST · before the D1 computation)

- **2a Direction prediction (review session inference):**
  - The c results point towards the simulator's noise floor being higher than real.
  - If so, the side that over-assumed the floor more is site B = OXF2PRI with 64 channels (L 8). So the NDI · ODI z means (D1 · k median) of OXF1PRI~OXF2PRI should be **both negative** — the side of «B assumed larger (32 → 64)» in the sensitivity table, about −3.
  - If both are positive, it opposes the «floor over-assumed» version of H7. If the signs split (one negative · one positive), no judgement is made.
  - The «sign agreement» rule of §22 stays as it is. Proviso: the chance agreement probability of one sign is 50 %, and if NDI and ODI move together the probability that both agree is also close to 50 %. It is read only as weak evidence.
- **2b:** even if the true-combination coverage rates of OXF1PRI~OXF3TRI (threshold of that combination only 0.80) · OXF2PRI~OXF3TRI (0.70) exceed the threshold, that is only a symmetric comparison after the verdict and not a pass of the fixed test.
- **2c Change to the cell classification (actual values):**
  - By the number of candidates in the judge's definition: **calibrated (candidate 1) 2 · candidate 2 12 · candidate 4 10**.
  - With 4a applied: **calibrated 2 (OXF1PRI~OXF2PRI NDI · ODI) · symmetric after 4a 2 (OXF1PRI~OXF3TRI · OXF2PRI~OXF3TRI · after seeing the results) · weak 10 (candidate 2 · rule-matched threshold 0.10 – 0.20) · no information 10 (candidate 4 · threshold 0)**.
  - Of the 10 no-information cells, the OXF3TRI~X 3 cells (NOT1ACH · NOT2ING-89 · NOT2ING-98) correspond to candidate 2 after 4a. However, the rule-matched threshold was computed under the judge's definition (candidate 4).
  - **Record change:** «weak 8 · no information 14» in §22-1 and the chat report are values A miscounted (the per-cell marks in rv_e2.md were right from the start).
- **Record — change to the C7 commit script:**
  - File: drafts/scripts/x14_c7_commit.sh. Time changed: 2026-11-02 18:03:56 KST (file modification time) · sha256 after the change 46eff42d6d95cf3d….
  - What changed: the D1 launch at the end was changed to «if explore/d1/STATUS_D1 exists, do not launch (write that it has already started) · if not, launch». This is to avoid launching twice when D1 is launched first. The rest stays as it is.
- **D1 starts together with C7 [author decision]:**
  - Order deviation: instead of «after the C7 commit» of §2, it starts while C7 is running. D1 does not depend on C7 results as data (C3 arm · v4_c3.json · test data).
  - The guard (guard.check) compares only the fixed analysis document and the pipeline files. So it is not blocked by the handoff · access record changes during the C7 commit.

## 23. D1 results (addendum 2026-11-03 01:00:56 KST · descriptive · not a verdict · verdict unchanged #34 · `explore/d1/`)

- **Run:** 2026-11-02 18:08:15 – 2026-11-03 00:53:07 (real 24,277 s · ran together with C7, sharing the CPU). real · ref · sens · summary all exit 0. Script sha256 30ab75e5… · d1_c3.json sha256 a5897fcd408e03ba….
- **Reproduction gate passed:** cells 29 (24 with data + 5 without) · coverage of nine k · number of ROI units used · number of misspecification exclusions · k-unmet count · per-combination counts are exactly the same as v4_c3.json (mismatched ROI units 0 · coils 0). ROI file re-reads 109 items recorded.
- The numbers are in d1_c3.md (script output). Only what is used for reading is moved here.

### 23-1. Selection after D1 — true combination of the symmetric cells · H7 per-combination PIT · |z| (addendum 2026-11-03 01:00:56 KST · no new computation · selected from D1 output)

- **Reason:** the review session's «D1 report format» ① ② and the §16-4 H7 definition (PIT · |z| · coverage per combination within the same ROI unit) need per-combination PIT · |z|. The D1 summary gives only per-combination coverage, and gives PIT · z under the «combination farthest from the reference» rule (pick_rule). So combinations were selected from the stored ROI unit results (d1_c3_units · d1_ref_*.npz) and the same statistics were computed. The ROI units used are the same as the D1 summary.
- **Reference (test data):** the cell's acquisition type · cases whose true combination is the cell's true combination · samples for nine k full in that combination · no misspecification flag for that combination.
- **[After seeing the results] normal-approximation decomposition:** when z ~ N(μ, s²), 90 % interval coverage = Φ((1.645 − μ)/s) − Φ((−1.645 − μ)/s). «μ only» = measurement μ · reference s, «s only» = μ 0 · measurement s. Descriptive only.
- Script drafts/x14/D1/d1_post_sym.py (sha256 30b994dd17b2ab08…) · output explore/d1/d1_post_sym.json (sha256 22a46d78bf065677…).

**Symmetric cells (true combination · measurement / reference):**

| Cell | True combination | ROI units | Coverage [k range] | Symmetric threshold (true combination) · stored threshold | PIT mean | Two tails | z mean | z SD | Normal: μ only · s only |
|---|---|---|---|---|---|---|---|---|---|
| OXF1PRI~OXF2PRI · NDI | 32\|64 | 191 | 0.010 [0.010, 0.010] | 0.80 · 0.80 | 1.00 / 0.48 | 0.99 / 0.08 | −4.43 / +0.05 | 1.49 / 0.94 | 0.00 · 0.73 |
| OXF1PRI~OXF2PRI · ODI | 32\|64 | 191 | 0.649 [0.634, 0.696] | 0.75 · 0.80 | 0.75 / 0.48 | 0.35 / 0.10 | −1.18 / +0.06 | 1.44 / 1.02 | 0.67 · 0.75 |
| OXF1PRI~OXF3TRI · ODI (after 4a) | 32\|32 | 162 | 0.883 [0.833, 0.926] | 0.80 · 0.80 | 0.30 / 0.51 | 0.12 / 0.06 | +0.71 / −0.04 | 1.36 / 0.93 | 0.84 · 0.77 |
| OXF2PRI~OXF3TRI · ODI (after 4a) | 64\|32 | 142 | 0.387 [0.338, 0.486] | 0.70 · 0.80 | 0.09 / 0.47 | 0.61 / 0.10 | +1.93 / +0.10 | 1.23 / 1.02 | 0.39 · 0.82 |

- «At or above the threshold» of OXF1PRI~OXF3TRI is a symmetric comparison after the verdict and not a pass of the fixed test (§22-3 2b).
- Other combinations of the same cells: OXF1PRI~OXF3TRI 32|64 coverage 0.68 · OXF2PRI~OXF3TRI 64|64 coverage 0.93 (|z| median 0.48).

**As the H7 definition (§16-4) — combinations with L 4 only (32|32) vs combinations containing 64 (same ROI unit · coverage):**
- The 10 cells with four candidates: coverage difference of 64|64 and 32|32 (64|64 − 32|32) 0.000 – +0.115 (all 10 cells ≥ 0). Comparing the mixed combinations (32|64 · 64|32) with the mean of the two same-L combinations, in 7 cells one is lower and the other higher (the direction changes with the scanner pair), and in 3 cells (NOT2ING-98~NOT1ACH · OXF3TRI~NOT1ACH · OXF4GEP~NOT1ACH — all B = NOT1ACH) both are lower.
- OXF1PRI~X 6 cells: cells where 32|64 has higher coverage than 32|32: 4 · lower: 2 (OXF1PRI~OXF3TRI · OXF1PRI~NOT4GEP).
- → **«Combinations containing 64 are worse» is not supported.**

### 23-2. End of C7 · automatic commit (addendum 2026-11-03 01:00:56 KST)

- C7 guard rerun: 2026-11-01 15:21:44 – 2026-11-03 00:12:50 · exit 0 · 118,268 s. v4_c7.json sha256 8f5e9105… · c7_guard_record.json sha256 37ee08eb…. Guard calls 147,300 · calls with negatives 4 · clipped elements 4 · observations 4 · «unreliable» 0 · minimum discriminant −2.86e-6.
- The monitor launched the commit script at 00:13:43 (handoff 00:13:44 · quick tests 36 passed 00:13:54 · full tests in progress). D1 had already started and was not launched again.
- For the report, A read v4_c7.json · c7_guard_record.json for the first time and re-read v4_c3.json (REREAD_LOG.jsonl).
- **Per-cell descriptive account (C7_READING_FIXED · no verdict · not combined):** the misspecification excess rate is the same as C3 in all 24 cells (misspecification is computed with X3 in both arms). ROI units dropped for not filling the nine k samples: C3 1,230 · C7 2,049 (out of 7,020 ROI units). Per-cell coverage is in the receipt table.
- (Change record 2026-11-03 01:01:26 KST) The «mixed combinations» sentence of §23-1 at first said only «split in opposite directions», from rounded values. With exact values it was changed to «7 cells opposite · 3 cells both lower».

## 24. M7 decision · changes to the D1 reading [author decision 2026-11-03 · review session reply «D1 · C7 completion report»] (addendum 2026-11-03 01:52:12 KST)

- **[Author decision] (2026-11-03):**
  - M7 = option (A).
  - The sufficiency test proceeds (CPU only · new image reads 0 · only re-reads of roi.json · not a verdict).
  - b2000 ventricular floor: definition only (reads 0 · the read is set again by the author after seeing the sufficiency results).
  - f postponed.
  - Record commit now (before the sufficiency test).
- **Order:** ① record commit → ② sufficiency definition fixed (this plan · time) → computation (CPU about 30 minutes) · while it runs ③ compression applying M7 (abstract · last paragraph of the introduction) → ④ report of sufficiency results → the author sets the b2000 read. f and other new branches are not opened.
- **Accepted (review session):** D1 reproduction gate · C7 commit 06b7361 · normal decomposition of the symmetric cells (the review session recomputed and got the same values). The verdict stays as it is (#34) · 0/24 remains as the verdict.
- **Changes to the reading [review session proposal] — the reading of §23 · receipt 2026-11-03a is rewritten as follows:**
  - (a) NDI is «bias mainly, spread secondarily» (with spread alone 0.73 < threshold 0.80 · z SD 1.49 vs 0.94). For ODI the two shares are similar (0.67 · 0.75). «Failed through bias» is not written.
  - (b) OXF1PRI~OXF2PRI has ΔTE = 0. By definition, the failure of this cell is not explained by TE relaxation (H5). This logic is written as it is in the results sentence.
  - (c) The sign pattern · H3 6/6 · H6 (leave-one-cell-out R² 0.93) show «additive per-scanner bias».
    - The coil-floor hypothesis and «some other peculiarity of one 64-channel scanner» predict the same, so it is not used as evidence separating the two.
    - The 6 comparisons of H3 share the sessions of the same two scanners and are not independent, and the cells of H6 are not independent either — both are used as descriptive only.
  - (d) The post hoc observation «covers under the same-L assumption, does not cover under a different L» holds only in the two cells between Siemens scanners (OXF1PRI~OXF3TRI 32|32 0.883 > 32|64 0.68 · OXF2PRI~OXF3TRI 64|64 0.93 > 64|32 0.387).
    - In the five cells OXF1PRI~(non-Siemens scanner), only one cell (NOT4GEP) has same L better.
    - In the 10 cells with four candidates, «both mixed combinations lower» is 3 cells.
    - → It is written only as «Siemens pairs · post hoc · two cells». For non-Siemens scanners the actual coil is unknown, so this comparison cannot say what is right or wrong.
  - (e) **Post hoc fix — H7 operationalisation [review session inference]:** «combinations containing L 8 are worse» of §16-4 was, in retrospect, not the right prediction. In the difference (B − A), the floor largely cancels between same-L assumptions, so what the L-difference hypothesis predicts is «mixed combination < same-L combination». The result of the original operationalisation (not supported · §23-1) stays as it is.
  - (f) **C7:**
    - The misspecification excess rate being the same as C3 is because both arms compute it with X3 (the same by construction) → it is not used as a basis for robustness.
    - k-unmet grew 1,230 → 2,049, so the per-cell coverages of C3 · C7 come from different sets of ROI units → when compared, the number of ROI units used per cell is written alongside.
- **M7 manuscript direction [review session proposal · author decision (A)]:**
  - **Main result:** in the only comparison on which the fixed test could give a verdict (pair differing only in coil · same TE), the judge failed.
    - NDI failed mainly through bias.
    - The sign was the same as the prediction written before computing from the judge's per-coil noise floor term (chance agreement about 1/2 — stated as weak evidence).
    - The cause is not asserted («consistent with»). There is one 64-channel scanner, so it cannot be separated from other peculiarities of that scanner.
  - **Secondary results (in this order):**
    - ① Design flaw (found after the verdict) · 22/24.
    - ② Gap (23/24 · FWF 0.080 vs 0.419 · Wong 25.4° · exploratory observation E1 10-31).
    - ③ Coil sensitivity of the frozen test data (0.9 → 0.3 · about 3 SD).
    - ④ Additive per-scanner descriptive account (H6 · descriptive · not independent · the Discussion implication is marked «not tested»).
  - **Left out of the main text:** «conclusion invalid» (B1-S · D3-v3 · D5 ③ ROI-unit outliers · ④ U-outlier) and «conditions not met» (D3 · D3 extension) get one line each in the supplementary ledger. D5 goes to the supplement. D2, together with c, is one paragraph of basis for ②.
  - **Compression order:** Abstract → last paragraph of the Introduction → Results (rearranged in the above order) → Discussion. Each section is committed and its word count · page count are reported. If the guidance gives no limit, main text within 8,000 words (review session's target).
- **f:** postponed. For Wong, only the observation (25.4° · exploratory) and the solution family of b are used. The D3 family is kept in the supplementary ledger as a record of truth-test failure.

### 24-1. Definition of the sufficiency test [author decision · review session reply 4 · before computing] (addendum 2026-11-03 02:19:46 KST · after record commit 9006d63)

- **What:** the OXF1PRI~OXF2PRI main cells (SNR 30 · N_eff 50 · C3) are drawn again by the same path as the verdict (collect → cond_rows → cover_by_k), but with three coil conditions [32|64 (the recorded true combination), 32|32, 64|64].
  - The pair seed is v4_c3.json seeds_torch (6000 + pair number). CPU only · S4_THREADS=16 · no OMP variable (same as D1) · new image reads 0 · only re-reads of roi.json (the sessions of this pair only).
- **Reproduction gate:** the first combination (32|64) must be bit-identical with the D1 stored values (pit · mean · sd · okc · covc of d1_c3_units/OXF1PRI~OXF2PRI.npz). Otherwise it stops with no values.
- **ROI units used:** the same ROI units as D1 (191) ∩ samples for nine k full in that combination.
- **Threshold:** none is newly made; the symmetric threshold of that cell (NDI 0.80 · ODI 0.75) is used.
- **Reading:**
  - Support = in all four cells, two same-L combinations × NDI · ODI, coverage ≥ threshold and |z mean| < 1.
  - Weakening = for any component, in both combinations (coverage < 0.5 or |z mean| > 2).
  - Otherwise = judgement on hold.
- **Expectation (mean response in the test data · the review session's rough figures recomputed with the match stratum values of rv_e2.json · linear assumption):**
  - 32|32 (B smaller than the actual 64 = side B 64→32): NDI −4.430 + 2.988 = **−1.443** · ODI −1.179 + 2.742 = **+1.563**. If z SD stays the same, normal-approximation coverage 0.534 · 0.510.
  - 64|64 (A larger than the actual 32 = side A 32→64): NDI −4.430 + 3.686 = **−0.744** · ODI −1.179 + 2.833 = **+1.654**. Coverage 0.672 · 0.487.
  - The review session's rough figures (pooled) were +3.019 · +2.761 · +3.334 · +2.823.
  - → If the judge responds on the measurement ROI units as in the test data, the result is «judgement on hold» (ODI near the weakening boundary).
  - The coil-term shift ratio (NDI/ODI) of the test data is 1.09 · 1.30, and the measurement bias ratio is 3.76. A mismatch of one coil alone cannot explain the two components together. It is explained only if the response differs on the measurement ROI units (low FWF tail) — this test measures that.
- **Report:** coverage · z mean · z SD · PIT mean · two tails · Δz (same L − 32|64, same ROI units) of the four cells are placed next to the expected shift.
- **Limit (written in advance):** even with «support», it goes only as far as «on these ROI units the judge's coil term is large enough to cover the bias». It still cannot be separated from other peculiarities of the 64-channel scanner.
- Script drafts/x15/suff/suff_c3.py (sha256 2c31d4b82d2649fa…) · output explore/suff/.

### 24-2. Definition of ① raw-value-level decomposition [author decision · review session addition 1 · before computing · new image reads 0 · the judge's posterior not used]

- **Measurement:**
  - On the ROI units used in D1 (191 · subjects 20), the measured coordinate difference D = X_B − X_A (x = logit NDI · y · ODI) is computed. X3 stored by D1 (values computed from roi.json «주») is used, so nothing is newly read.
  - Mean · median · SD · quartiles · subject-cluster bootstrap 2,000 (seed 20261103) 95 % interval. For reference the apparent NDI difference is also given.
- **Simulator expectation (frozen test data sim_test.npz · match|L0 · x_c3 finite · forward values only):**
  - For each of the true combinations 32|64 and 32|32 (about 250 cases), per component D = a + b·θ (dNDI · dFISO · dODI) + c_A/SNR_A² + c_B/SNR_B² is fitted by least squares.
  - It is evaluated at θ = 0 · the SNR of the measurement pair (per site, the subject median of 1/SNR² of those sessions · the two sessions of the attribute table · snr_table re-read).
  - The interval is a case bootstrap 2,000.
  - Reason: the test data drew SNR 15 – 60 separately per site and the tissue difference θ is not 0 (dNDI ±0.1 · dFISO ±0.045 · dODI ±0.09 · only design variables were checked). So «tissue difference 0 · the SNR of this pair» is obtained by regression.
  - Sensitivity: (s1) with 1/SNR · (s2) θ regression intercept for cases with both SNR ≤ 30 only · (s3) evaluated with measurement SNR below 15 held at 15.
- **r = (measurement mean − expected 32|32) / (expected 32|64 − expected 32|32)** is computed for components x · ODI (y is descriptive). 0 means at the 32|32 expectation, 1 means at the 32|64 expectation.
- **Reading (before computing · boundary values are [A's proposal]) — applied in the order below:**
  1. If, for any component, the interval of Δ_L (= expected 32|64 − expected 32|32) includes 0: «cannot be separated».
  2. **(B) side:** r_x or r_ODI is outside [−0.5, 1.5], or |r_x − r_ODI| > 1.
  3. **(A) side:** both r in [−0.5, 0.5].
  4. **Agrees with the simulator assumption (L 8 vs 4):** both r in [0.5, 1.5].
  5. Otherwise «cannot be separated».
- Review session prediction: (A) the measurement is near the 32|32 expectation and far from the 32|64 expectation (both NDI · ODI) · (B) far from both expectations, with an NDI : ODI ratio not explained by L.
- **Limits (written in advance):**
  - The tissue location of the test data differs from the measurement (low FWF tail) (competing cause (D)).
  - Some measurement session SNR are outside the test range (15), so extrapolation by 1/SNR² is used (checked with s3).
- **Check of SNR direction (review session addition 0):** in the quality check of Warrington et al. 2025, the scanner with the highest diffusion SNR/CNR was the Prisma 32-channel. One line gives the session SNR median of this pair's subjects for OXF1PRI vs OXF2PRI.
- Script drafts/x15/cause/cause1_raw.py (sha256 c50c72ca5a79a19f…) · output explore/cause/.

### 24-3. Definition of ② signal strength · location trend [author decision · review session addition 2 · before computing · new image reads 0 · descriptive]

- **ROI units:** the ROI units used in D1 for OXF1PRI~OXF2PRI (191). They are clustered in 20 subjects and are not independent.
- **Predictors:**
  - **b2000 signal/noise:** the mean over the two sessions of log(SNR_session × exp(s_2000)).
    - s_2000 = log(S̄_2000/S0) of roi.json «주», and only roi.json is re-read, without re-reading the D5 maps.
    - Proviso: local b0 intensity differences due to coil sensitivity do not enter. s_2000 itself is also affected by the noise floor.
  - **Location:** the label representative voxel in the repository reference file JHU-labels.xml (not an image) → MNI mm (x = 90 − i · y = j − 126 · z = k − 72) → the label mean of the ROI unit.
    - From this, superior-inferior z (mm) and centre distance (distance from the mean coordinate of the 27 ROI units) are computed.
    - Proviso: the representative voxel is not the centroid. Location is the same per ROI and does not vary between subjects.
- **Outcome variables:** D1 ROI-unit z (32|64 · NDI · ODI · k median) · residual of ① (measurement D − 32|64 expectation · at the two session SNR of that ROI unit's subject · x · ODI).
  - Alignment to the bias side: o = outcome × sign(ROI-unit mean) — the larger o is, the larger towards the mean bias of that cell.
- **Descriptive:** per predictor, Spearman ρ(predictor, o) and the standardised coefficients · R² of least squares with all entered together (three standardised predictors + intercept) are given. The interval is a subject-cluster bootstrap 2,000 (seed 20261103).
- **Prediction before computing (review session) · reading rule:**
  - (A) The weaker the signal, the larger the bias → ρ(b2000 SNR, o) < 0.
  - (B) Nonlinearity → follows location (especially superior-inferior z) → ρ(z, o) ≠ 0.
  - (C) Follows depth (surface vs deep) → ρ(centre distance, o) ≠ 0.
  - «Agrees with the prediction» = when the interval excludes 0 and is in the predicted direction. Otherwise «cannot be separated» (or «opposite direction»).
- **Depending on the answer to question 3:** since no gradient-nonlinearity unwarping was done (§24-4), the location predictors are looked at together as main predictors.
- Script drafts/x15/cause/cause2_trend.py (sha256 e5400be7f2e44f22…) · output explore/cause/.

### 24-4. Answer to question 3 — gradient-nonlinearity unwarping (code check only · reads 0)

- **It was not done.** The processing steps are prep → topup → bet(mask) → eddy → dtifit → flirt · fnirt (→ FMRIB58) → applywarp of decomp/remeasure/run.py:216 – 225 and cmds.py.
- Gradient-nonlinearity unwarping (gradient unwarp) steps · coefficient files · related settings: 0 items found in decomp/remeasure/ (search for unwarp · gdc · coeff · nonlin).
- → As in review session addition 0, the Siemens scanners acquired with on-scanner nonlinearity unwarping turned off, so the nonlinearity share of (B) remains in the processing outputs. In ② location is looked at together as a main predictor.

### 24-5. Factual fix — OXF1PRI · OXF2PRI [review session addition 0 · Warrington et al. 2025 full text · confirmed by the review session]

- OXF1PRI and OXF2PRI are not one scanner with only the coil changed but **two different Prisma scanners** (32-channel · 64-channel).
- «Pair differing only in coil» in the plan · receipts · ledger · abstract draft is read as «same model (Prisma) · same TE · different scanner and different coil (32 · 64 channels)». It is changed in the manuscript · figures · talks.
- **Competing causes (review session):** (A) noise floor assumption (L) · (B) scanner systematic effects (gradient strength calibration · nonlinearity + head position depending on the coil) · (C) SNR size · spatial distribution · (D) the judge's response on the measurement ROI units (low FWF tail).
- **Order (review session 4):** record commit → ① · ② · sufficiency (all defined before computing) → reported together → the author sets the b2000 floor read.
  - If ① points to (A) and ② shows signal dependence, the floor measurement is of high value as a confirmation (including OXF3TRI recommended).
  - If ① points to (B), the floor measurement is not done. The manuscript writes «scanner difference · cause unknown · the noise floor term does not match in size».
- **Out of scope (only in the manuscript's limits paragraph):** fully separating (A) and (B) would need a second 64-channel scanner, or data from one scanner acquired with 32 · 64 channels swapped.
- **b2000 floor definition (author decision: definition only · reads 0):** review session reply 5 is accepted as it stands.
  - Ratio = mean/SD. Theoretical values L 1 1.91 · L 2 2.75 · L 4 3.94 · L 8 5.61 — A recomputed with the central chi formula (mean = σ√2 Γ(L + ½)/Γ(L) · variance = 2Lσ² − mean²) and confirmed they are the same.
  - **Within-direction ratio (main) [A's proposal]:** mean/SD across the b2000 directions in one voxel → median over ventricular voxels.
    - Reason: the between-voxel ratio mixes in spatial variation of σ (coil sensitivity · g factor).
    - The between-voxel ratio is given only as a sensitivity analysis.
  - **Ventricle mask [A's proposal]:** apparent FWF ≥ 0.95 ∩ brain mask, eroded by 2 voxels (reduces partial volume).
  - **Proviso:** eddy resampling · CSF pulsation · partial volume · whether the processed DWI has a scale transform (checked from the header · processing record before opening).
  - **Truth test:** in a simulation made with known L (4 · 8), the same resampling ([A's proposal] sub-voxel shift of eddy size · trilinear interpolation) is applied, and it is first checked whether the ratio separates L 4 and 8.
  - **Prediction:** under the L-difference hypothesis, the ratios of OXF2PRI and OXF1PRI are similar. If the simulator assumption (L 8 vs 4) is right, OXF2PRI is clearly higher.
  - Whether to add OXF3TRI to the targets OXF1PRI · OXF2PRI is for the author to set.

### 24-6. ① · ② results (addendum 2026-11-03 02:23:15 KST · descriptive · not a verdict · verdict unchanged #34 · explore/cause/)

- **① Raw-value-level decomposition** (02:20:20 · cause1_raw.json sha256 86b0b8d5a783abec…):
  - **Measurement difference (OXF2PRI − OXF1PRI · ROI units 191 · subjects 20):**
    - x (logit NDI) **+0.020** [+0.010, +0.032] · SD 0.059
    - ODI **+0.002** [+0.001, +0.003] · SD 0.007
    - y +0.011 [−0.019, +0.041]
    - Apparent NDI difference mean +0.005
  - **Simulator expectation (θ = 0 · evaluation SNR A 18.6 · B 15.6):**
    - 32|64 — x **+2.62** [+2.38, +2.90] · ODI **+0.137** [+0.120, +0.153]
    - 32|32 — x **+0.71** [+0.51, +0.93] · ODI **+0.034** [+0.024, +0.044] (B's SNR is lower, so a floor difference is expected even with the same coil)
  - **r:** x **−0.36** [−0.55, −0.23] · ODI **−0.30** [−0.47, −0.19] (y −0.23). In both components the Δ_L interval excludes 0.
  - **Reading (rule before computing): «(A) side — both components are near the 32|32 expectation and far from the 32|64 expectation».**
  - The three sensitivity analyses also fall in the same band:
    - s1 (1/SNR): r x −0.23 · ODI −0.21
    - s2 (SNR ≤ 30 intercept): r about +0.10 · +0.08
    - s3: the median of the evaluation SNR is above 15, so the same as the primary
  - **[A's reading · after seeing the results]**
    - The measurement difference is even closer to 0 than the 32|32 expectation (r < 0). Between these two scanners there is almost none of the difference the simulator's noise floor predicts, whether from L difference or SNR difference.
    - The r of the two components are similar (difference 0.06). The «NDI bias 3.75 times that of ODI» seen in the judge's z does not appear at the raw-value level — it reads as a difference arising in the scale of the judge's posterior.
    - Proviso: the tissue location of the test data differs from the measurement (low FWF tail), so the expected values are values at test-data tissue ((D)).
  - **SNR direction:** the session SNR median of this pair's subjects OXF1PRI 18.56 · OXF2PRI 15.59 · ΔSNR (B − A) median −2.44. The direction is the same as the Warrington et al. quality check (Prisma 32-channel has the highest diffusion SNR).
- **② Signal strength · location trend** (02:21:00 · cause2_trend.json sha256 41f2fb1bc0486627…):
  - ρ among predictors: SNR~z +0.20 · SNR~distance −0.03 · z~distance −0.52.
  - **Judge z (aligned to the bias side):**
    - NDI: ρ(b2000 SNR) −0.26 [−0.41, −0.12] · ρ(z) −0.43 [−0.55, −0.30] · ρ(distance) +0.11 [−0.04, +0.26] → by the rule (A) agrees · (B) agrees · (C) cannot be separated. Entered together β −0.17 · −0.41 · −0.15 · R² 0.17.
    - ODI: −0.35 · −0.37 · +0.35 → (A) · (B) · (C) all agree. β −0.38 · −0.17 · +0.18 · R² 0.27.
  - **① residual (measurement − 32|64 expectation):** for both x · ODI ρ(b2000 SNR) −0.40 · −0.42 → (A) agrees. Location near 0 → (B) · (C) cannot be separated.
  - **[A's check · after seeing the results] limit of the definition:** the measurement D ≈ 0, so residual ≈ −(expectation), and the expectation is a function of the subject's two session SNR.
    - So the SNR correlation of the residual is mostly mechanical (the SNR dependence of the simulator expectation itself).
    - The location correlation is masked by the between-subject variation of the expectation.
    - → The residual line is not used as independent evidence for «(A)».
  - **SNR · location dependence of the judge's z:** if the floor difference the judge assumed is larger where the signal is weak · depending on tissue, this is the same form as arises when the measurement difference is 0. It is not independent evidence.
- **Descriptive account after ② — the measurement D itself [after seeing the results · outside the definition · not a verdict]** (cause2_posthoc_D.py sha256 c0475bd1f13d38b9… · output sha256 09eb3da0d4759336…):
  - **D_x:** ρ(b2000 SNR) −0.19 [−0.35, −0.00] · ρ(z) +0.18 [+0.03, +0.34] · ρ(distance) −0.10 [−0.23, +0.05].
    - Version with the subject mean subtracted: −0.14 [−0.27, +0.01] · +0.22 [+0.05, +0.39] · −0.11.
  - **D_ODI:** all near 0 (−0.08 · −0.02 · −0.06 · within subject −0.14 · −0.05 · −0.04 · intervals all include 0).
  - → The measurement NDI difference has a weak superior-inferior trend (OXF2PRI slightly higher in superior ROIs) and a weak signal trend. Its size is variation within D_x SD 0.059, far smaller than the floor difference the simulator expected (+0.7 – +2.6).
  - It may fit a small share of (B) (nonlinearity · head position), but at this size it is not the main cause of the judge's failure.

### 24-7. Sufficiency test results (addendum 2026-11-03 02:44:28 KST · descriptive · not a verdict · verdict unchanged #34 · explore/suff/suff_c3.json sha256 2006101930e0966b…)

- **Run:** 02:19:55 – 02:43:33 (1,415 s · CPU). Reproduction gate passed (pit · mean · sd · okc · covc of 32|64 bit-identical with the D1 stored values · rows 540).
- **Results (ROI units used · coverage · z mean · z SD · PIT mean · two tails · Δz (− 32|64) · expected shift):**

| Combination · component | ROI units | Coverage / threshold | z mean | z SD | PIT | Two tails | Δz | Expected shift (test data) |
|---|---|---|---|---|---|---|---|---|
| 32\|64 NDI | 191 | 0.010 / 0.80 | −4.43 | 1.49 | 1.00 | 0.99 | — | — |
| 32\|64 ODI | 191 | 0.649 / 0.75 | −1.18 | 1.44 | 0.75 | 0.35 | — | — |
| 32\|32 NDI | 191 | **0.921** / 0.80 | +0.25 | 0.91 | 0.44 | 0.08 | +4.69 | +2.99 |
| 32\|32 ODI | 191 | **0.901** / 0.75 | +0.60 | 0.98 | 0.31 | 0.10 | +1.81 | +2.74 |
| 64\|64 NDI | 189 | **0.931** / 0.80 | +0.13 | 0.97 | 0.48 | 0.07 | +4.61 | +3.69 |
| 64\|64 ODI | 189 | **0.921** / 0.75 | +0.53 | 0.84 | 0.33 | 0.07 | +1.77 | +2.83 |

- **Reading (rule before computing): «support».** In all four cells coverage ≥ threshold and |z mean| < 1.
  - z SD 0.84 – 0.98 is at the level of the test-data reference (NDI 0.94 · ODI 1.02). Under same-L assumptions the overconfidence also disappears.
- **What differed from the expectation:**
  - The expectation written in advance (mean response in the test data · linear) was «judgement on hold» (expected coverage 0.49 – 0.67).
  - On the measurement ROI units, the judge's coil-term response is larger than expected for NDI (+4.6 – 4.7 vs +3.0 – 3.7) and smaller for ODI (+1.8 vs +2.7 – 2.8). The shift ratio NDI/ODI is 2.6 (test data 1.1 – 1.3).
  - → «It is explained only if the response differs on the measurement ROI units (low FWF tail)» of review session reply 4 did in fact happen. The one coil floor term covers the bias of both components together (ratio 3.76).
- **Added descriptive account (a report in the definition · not used for selecting ROI units):** the number of misspecification flags of that combination (among the 191 ROI units used in D1) is 32|64 0 · 32|32 1 · 64|64 64.
- **Limits (as written in advance in §24-1):**
  - «Support» goes only as far as «on these ROI units the judge's coil term is large enough to cover the bias».
  - There is one 64-channel scanner, so it cannot be separated from other peculiarities of that scanner.
- **[A's reading · after seeing the results · together with ①]**
  - In ①, the raw-value difference between the two scanners is close to 0 (x +0.020 · ODI +0.002). Giving the judge the same noise floor at both sites yields calibrated coverage.
  - The judge's failure is consistent with «the judge expects a floor difference attached to the recorded coils (32 · 64), but the measurement has no such difference» (consistent with).
  - For (B) (scanner systematic effects) to produce this form, the scanner effect would have to almost exactly cancel the real floor difference. That is possible, but there are no data to separate it.
- **Next-step rule of review session reply addition 4:** ① is on the (A) side and ② showed signal dependence (by the rule · with the proviso of a mechanical share). → The floor measurement is of high value as a confirmation (including OXF3TRI recommended). Read approval is for the author to set.

## 25. Handling the reply on narrowing the cause [author decision 2026-11-03 · review session reply «narrowing the cause · sufficiency · record · manuscript»] (addendum 2026-11-03 17:26:59 KST · 2(b) · 2(c) · b2000 are before computing · reading)

- **[Author decision] (2026-11-03):**
  - b2000 floor read approved: OXF1PRI · OXF2PRI · OXF3TRI first sessions · step by step · CPU only.
  - The sufficiency sentence goes in the abstract (wording of review session reply 3).
  - The next sections are Results → Discussion.
- **Accepted (review session):** record commit 9006d63 · definition fixed 02:19:46 · sufficiency reproduction gate · ① «(A) side» · sufficiency «support». Verdict unchanged (#34).
  - The limit of the definition of ② was a weakness in the review session's design — ② is not used as independent evidence.
  - The review session's expectation before computing («judgement on hold») was wrong — on the measurement ROI units, the coil-term response was NDI +4.6 – 4.7 · ODI +1.8, unlike the test-data mean (+3.0 – 3.7 · +2.7 – 2.8).
    - To secondary result ③ is attached «on the measurement ROI units the response differed (NDI larger · ODI smaller)».
    - One line of methods lesson in the Discussion: sensitivity measured on the test-data mean does not carry over to the measurement ROI units.
  - «22 cells … vs 0.80» in the review session's previous reply is a slip by the review session. The threshold per cell is 0.50 – 0.85.
- **Changes to the reading [review session proposal]:**
  - (a) ① and sufficiency are not two independent pieces of evidence.
    - The main cells assume SNR 30 at both sites, so with the same L, the acquisition difference the judge expects disappears in this pair with ΔTE = 0. If the measurement D ≈ 0, the location (z mean ≈ 0) almost follows.
    - The information sufficiency adds is the spread (z SD 0.84 – 0.98).
    - In the manuscript it is written as «one fact seen by two routes».
  - (b) ① used the ROI units selected by the judge (191/540). To write «a comparison that does not pass through the judge», the ROI units must also be independent of the judge → 2(b) below.
  - (c) In the test-data expectation, the number of cases dropped by the x_c3 finite filter is written per combination → 2(c) below.
  - (d) The misspecification 64/191 of 64|64 vs 1/191 of 32|32 is not interpreted (the ROI units were selected on the 32|64 basis, so it is asymmetric).
  - (e) The descriptive account after ② (the weak superior-inferior trend of D_x) is used only as one line in the supplementary material.
- **Widening the central sentence [review session proposal]:** the main result is widened from «L difference» to «the simulator's noise floor is stronger than the measurement» (consistent with · not asserted).
  - Basis: in ①, the measurement is closer to 0 even than the same-coil expectation (+0.71 · the floor arising from the measured SNR difference).
  - Explained together: the main result · secondary result ② (seven scanners FWF 0.080 vs 0.419) · D2 (repeat spread smaller than under the SNR 30 assumption). The seven-scanner evidence does not depend on the one scanner OXF2PRI.
  - Remaining limit: in this pair, «the scanner effect exactly cancels the floor difference» cannot be excluded.
  - Priority in the abstract: these three sentences > the 22-cell design flaw > the rest. In place of «same coil assumed», «equal noise-floor settings … chosen after the verdict».
- **Not to be done now [review session proposal]:** all 24 cells are not run again with the same L (that would be changing settings after the verdict to look for a pass) — only «not tested» in the Discussion.

### 25-1. Definition of 2(b) · 2(c) (before computing · new reads 0)

- **2(b) ROI units unrelated to the judge:**
  - ROI unit = every row produced by collect for OXF1PRI~OXF2PRI (ROI units with values of the seven maps in both sessions · at most 540).
  - The values are the X3 saved by D1 (all rows · values produced by collect from roi.json «주») — collect has no relation to the judge and the values are the same, so roi.json is not read again.
  - D = X_B − X_A (x · y · ODI). Mean · median · SD · subject-cluster bootstrap 2,000 (seed 20261103) 95 % interval.
  - r = (mean − expected 32|32)/(expected 32|64 − expected 32|32), computed with the expectation of ① (same 20 subjects · same evaluation points).
  - Written side by side with the values of the 191 ROI units. If there are ROI units excluded from the main cell (main_excluded), their number is written and a value with them dropped is also given as a sensitivity check.
- **2(c):** in the match|L0 stratum, for each of the true combinations 32|64 · 32|32, the number of cases · the number dropped by the x_c3 finite filter · the number flagged clipped_A · clipped_B are written (descriptive).
- Script drafts/x15/cause/cause1b_allunits.py · output explore/cause/cause1b_allunits.json.

### 25-2. b2000 floor — steps · predictions (review session reply 4 · before computing · before reading · fixed at the time of this addendum)

- **Target:** first sessions of OXF1PRI · OXF2PRI · OXF3TRI (20 + 20 + 10).
- **Step 1 (new read · small):** check from headers and processing records whether processed DWI remains · the value area · whether a scale conversion was applied.
  - If there is no processed DWI, stop here and report. Do not process again from the sealed raw images.
  - [A's check · code only] run.py:177 – 187 prune_onh keeps in the session folder only KEEP_ONH (bvals · bvecs · brain mask · dti_FA · transforms · roi.json · snr.json etc.) · f2_* · log_* · eddy text logs and deletes the rest (including eddy_out.nii.gz). → Step 1 is likely to be «does not remain». First check with a folder listing (file names · sizes · not reading image contents).
- **Step 2 (reads 0) — truth test:** put the same resampling into simulations with known L (1 · 2 · 4 · 8), and see whether the within-direction ratio separates L · how much the resampling raises the ratio.
  - If it cannot separate them, Step 3 is not done.
  - [A's proposal] Details (used only if Step 1 is passed):
    - For each voxel, noncentral chi for as many values as the b2000 directions (bvals of that session) (residual signal ν/σ 0 · 0.5) → for each volume, spline resampling of the same order as eddy (random sub-voxel shift ±0.5 voxel).
    - Within-direction mean/SD → voxel median.
    - Separates = the ratio difference between neighbouring L is greater than 4 times the Monte Carlo SD at the measurement size (number of ventricle voxels).
- **Step 3 (new read) — ventricle b2000:**
  - Main: within-direction mean/SD → effective L through the calibration table of Step 2.
  - Secondary (descriptive): ratio of σ obtained from the floor to σ implied by the session SNR (whether the measured SNR takes the noise as large · H1′).
  - Caveat: because of prescan normalise · g factor, σ varies over space.
- **Prediction before computing (review session):**
  - If the simulator assumption holds: OXF2PRI effective L ≈ 8, OXF1PRI · OXF3TRI ≈ 4.
  - If the floor is weaker than in the simulator: the effective L of the three scanners is similar and 4 or lower, or σ ratio < 1.
  - If effective L differs as the simulator assumes but measurement D ≈ 0: (B), that the scanner effect cancels the floor difference, remains.
- **Record:** each file is logged in ONH_NEW_READ_LOG and passes through the guard. It is moved to the canonical record at the next record commit.

### 25-3. Results of 2(b) · 2(c) (addendum 2026-11-03 17:28:09 KST · computed 17:27:38 · cause1b_allunits.py sha256 94ebcccfa006e78f… · output sha256 6f69f3710e885367… · descriptive · not a verdict)

- **2(b) ROI units unrelated to the judge (all rows of collect 540 · subjects 20 · main cell exclusions 0):**
  - x **+0.013** [+0.006, +0.021] · SD 0.061 · r **−0.37** [−0.55, −0.23]
  - ODI **+0.0019** [+0.0010, +0.0027] · SD 0.009 · r **−0.31** [−0.47, −0.19]
  - y +0.006 [−0.034, +0.043] · r −0.21
  - 191 ROI units: x +0.020 · ODI +0.0021 · r −0.36 · −0.30 (intervals are from this script's bootstrap).
  - → The reading band is the same («on the (A) side»). The number in the abstract is fixed at the value of the 540 ROI units (+0.01).
- **2(c):** in match|L0, for both the 32|64 241 cases · the 32|32 248 cases, cases dropped by the x_c3 finite filter 0 · clipped_A 0 · clipped_B 0. → The expectation was not underestimated. The script confirmed that the point estimate of the expectation is the same as in ①.

### 25-4. b2000 Step 1 result — no processed DWI · stopped here (addendum 2026-11-03 17:31:24 KST · new image reads 0)

- **How checked:** file listing of the 50 target session folders (OXF1PRI 20 · OXF2PRI 20 · OXF3TRI 10 · first sessions) (names · sizes only · image contents and headers not opened).
- **Result:** none of the 50 sessions has processed DWI (4D images containing b2000 volumes, such as eddy_out.nii.gz).
  - What remains is bvals · bvecs · acqparams · index · brain mask · dti_FA · fa2std transform · f2_* maps 9 · eddy text logs 13 · log_* · roi.json · snr.json · timings · DONE. As designed, prune_onh (run.py:177 – 187) deleted the rest.
- **Further search by name (under data · names only):**
  - eddy_out.nii.gz exists only at timing/sub-13192_ses-OXF1PRI001 · 002 (the output of the 2026-10-26 eddy timing measurement).
  - OXF2PRI · OXF3TRI are absent, so these cannot be used for a scanner comparison, and they were not opened.
- **→ As per review session reply 4 · author decision, stop at Step 1.** Do not process again from the sealed raw images. Step 2 (truth test) and Step 3 are not done (without Step 3 the calibration table has no use).
- **What the author · review session are to decide (written in the report):**
  - (i) Write it in the manuscript limitations without a floor measurement.
  - (ii) Read the b2000 volumes of the sealed raw images without processing and measure the within-direction ratio — new read · approval needed. With no resampling the truth test becomes simpler, but motion · mask alignment are the caveats.
  - (iii) Process again — forbidden by the author decision.

### 25-5. Schedule of the 11/8 draft · scope frozen [author decision 2026-11-03 · review session proposal] (addendum 2026-11-03 17:40:20 KST)

- **[author decision]:** 11/8 (Sunday) = complete draft for internal review. Scope frozen from today.
- **Meaning of draft (review session):**
  - All sections at the target length of the relocation map (main text about 8,400)
  - An evidence tag on every number
  - Figures are drafts with the right content
  - Supplementary material: analysis ledger · relocated method details
  - Statements · Highlights · title are drafts
  - Not done on this day: final figure format · full check of references · submission file bundle
- **Scope frozen:**
  - New analyses are only 2(b) · 2(c) · b2000.
  - If b2000 has no result by 11/5 12:00 KST, it stops and goes to a limitations paragraph.
  - f · rerun of the 24 cells · any other new branch are not opened.
- **By date (review session):**
  - 11/3: 2(b) · 2(c) → abstract revision → start of Results compression · b2000 Steps 1 – 2
  - 11/4: end of Results · draft of the main result figure (measurement D vs simulator expectation · coverage under equal floor settings) · b2000 Step 3
  - 11/5: Discussion · b2000 deadline (noon) · start of Methods compression
  - 11/6: end of Methods · rest of the Introduction · supplementary material · redrawn draft of Fig 4
  - 11/7: statements · Highlights · title draft · full cross-check of number evidence (script) → read-through by the review session
  - 11/8: read-through fixes applied · draft frozen (PDF · commit · hash)
- **Operation:**
  - Reports are batched twice a day (noon · evening KST).
  - Each report includes the main text of the sections written · word count per section (now / target) · page count.
  - Things needing an author decision are put in one batch at the top of the report.
  - The full test suite of 1,302 is run at every commit as it is, and the next section is written while the tests run.
- **Status today (11/3):**
  - 2(b) · 2(c) finished (§25-3).
  - Abstract revision S3 being committed.
  - Results compression R-a (S4) · R-b (S5) written and being committed after S3.
  - b2000 stopped at Step 1 with «no processed DWI» (§25-4). Step 2 (truth test) is a calibration table for Step 3, and without Step 3 it has no use, so it was not done.
  - The path of reading the b2000 of the sealed raw images without processing (§25-4 (ii)) is a new read, so it goes into the author decision batch. The deadline is noon on 11/5.

## 26. NeuroImage integrated format confirmed · decisions withdrawn · framework fixed until 11/8 [author decision 2026-11-03] (addendum 2026-11-03 20:54:18 KST)

- **[author decision] (2026-11-03 · chat):**
  - ① Design confirmed: the fixed FISO · ODI check in the main text · 0.20–0.77 left out of the abstract · the claim is based on coordinates.
  - ② Framework fixed until 11/8 · check the list of decisions withdrawn below.
  - ③ The lesion T2 example has no freeze exception → the existing −0.148 is used.
  - ④ b2000: the 11/8 draft uses limitations paragraph (i). Direct reading of the sealed raw images (ii) is decided again after 11/9.
  - ⑤ Access to Tong 2020 Technical Validation allowed — the access time is recorded, only comparison · citation are done, and no reanalysis.
  - ⑥ Word count: main text only (Introduction – Conclusion). Abstract · captions · tables · boxes · appendices · statements reported separately. Methods section upper limit at 11/8 2,500 (submission target 2,100).
- **Four review session errors recorded by the review session:**
  - Used 0.20–0.77 as a simulation of the relaxation mechanism (values from a design with the same relaxation · overlapping sites and populations)
  - Used 48–72° as a literature value (the simulator's prior grid · with literature T2 it is 58°)
  - Said the fixed check should be moved down to the supplement (withdrawn)
  - Widened to «the floor is stronger overall» with the same-coil expectation (+0.71) and proposed an L = 1 computation (withdrawn · what is solid is the recorded-coil expectation +1.4 – +2.6)
- **Fixes from A accepted:**
  - Abstract: +0.013 (540 ROI units) · «differed little» · 22/24 · one 64-channel scanner
  - Title (B) and Highlight 4 dropped
  - 19-pair cell in the results ladder · population bridge «derivation only»
  - Citation caveats · banned-word cleanup · disclosure of the NDI cell 0.003 difference
  - SNR 50 / 20 sensitivity used as the main basis of the noise account
    - **[A's check · just before this addendum · judgements.json V4.cells.*.sensitivity_passed]** The cells that changed are four ODI cells.
      - OXF1PRI~OXF2PRI ODI: passed at SNR 50 (N_eff 25 · 100), failed at SNR 20 · N_eff 25, no verdict possible at SNR 20 · N_eff 100
      - OXF1PRI~OXF3TRI · OXF2PRI~OXF3TRI · OXF4GEP~NOT1ACH ODI: passed only at SNR 50 · N_eff 25
    - **The NDI cell of the same pair failed in every sensitivity cell.** So «the assumed noise changes the verdict», which this fixed sensitivity supports, applies only to ODI.
- **Decisions withdrawn (framework fixed until 11/8):**
  - Main-result status of M7 (A)
  - The sufficiency sentence in the abstract (kept in Results only)
  - Fig. 3 (identifiability map → supplement)
  - «Results → Discussion» order
  - Withdrawn by the review session: the central claim «the floor is stronger overall» · L = 1 computation · «problem size» figure · 0.26 logit example · 0.20–0.77 citation · new Box 2
- **One-sentence claim (A's version adopted · final confirmation 11/7):** «Single-echo-time NODDI indices are acquisition-indexed measurements: in the signal model, relaxation weighting is a translation that sets necessary conditions for reading comparisons; in tests fixed before computing, center effects within one scanner model and protocol stayed at retest level for logit NDI and ODI but not for the free-water coordinate, physics-based predictions of cross-scanner FISO shifts were not borne out, and a simulation-trained judge failed on traveling heads—post hoc, consistent with a simulated noise floor stronger than in the data.»
- **Results ladder (evidence grade per cell):**
  - ① Same model · protocol — Tong T-B (blinded · fixed): x 0.931 · ODI 0.934 · y 0.844
  - ② Same model · protocol · different scanner — two Prismas 540 ROI units (post hoc)
  - ③ 19 pairs across protocols · scanners (fixed · unblinded · descriptive): FISO sign 2/8 · τ 0.222 · ODI (b) 2/19 · T-A no verdict possible
  - ④ Crossing TE — Wong. The fixed linearity (x 0.815 · y 0.901) is in Table 2. **[A's addition]** Since Table 2 is «all fixed tests», the main text of ④ also gives the fixed result first. Post hoc · exploratory: slope ratio x 1.13 · y 0.31 · θ 25.4° vs literature 58°
  - ⑤ Crossing populations — derivation only · no measurement (example −0.148)
  - ⑥ Simulation judge — verdict 24/24 failed · 22/24 design flaw · SNR 50/20 fixed sensitivity (caveat above) · recorded-coil expectation +1.4 – +2.6 vs measurement +0.013 · coverage under equal floor settings 0.90 – 0.93 (post hoc · in Results only)
- **Structure:**
  - Length: Introduction 900 · Theory 1,600 (+ appendix) · Methods 2,100 (11/8 upper limit 2,500) · Results 2,600 · Discussion 1,500.
  - Figures: F1 concept · F2 existing algebra figure · F3 ①② · F4 ③④ · F5 ⑥.
  - Table 1 rules · Table 2 all fixed tests.
  - Boxes: Box 1 = attached draft v2. Box 2 = reporting list + simulation bridge 3 lines (threshold as an applied rule · raw-value decomposition · same-assumption sufficiency).
- **Limitations paragraph:** one 64-channel scanner · floor not measured (SoS question up to consistent with) · L set by channel count · σ set by repeated b0 · Tong TV comparison (after ⑤).
- **Title · Highlights:** «tested on traveling heads» is not used. Two title candidates ⟦D⟧ · Highlights are based on the version whose character counts A has checked (final 11/7).
- **Citations:** as in review session reply 8. Review session addition: in Sakaie 2018 the Siemens uses adaptive combine (our Prisma uses SoS).
- **Schedule:**
  - 11/3: decision record → ladder table · abstract v0 → start of Results · 4.4.1 revision
  - 11/4: end of Results · drafts of F3 – F5
  - 11/5: Discussion · start of Methods compression
  - 11/6: end of Methods · compression of Theory · Introduction · supplement · banned words
  - 11/7: final abstract · claim · Highlights · title · statements · full number cross-check → read-through by the review session
  - 11/8: fixes applied · freeze
  - Reports at noon · evening.
- **Attached draft v2 (Box 1 · implications paragraph · three statements):** received by chat. Saved as it is in drafts/x15/review/NI_parallel_drafts_v2.md.
  - When inserted: Box 1 · implications paragraph in the Discussion (11/5), data · ethics statements in the Methods (11/5 – 11/6).
  - **[A's check]** The two commits in v2 are in the repository (5f95b18 2026-10-28 01:56:00 fixed analysis document · 5632626 2026-11-01 15:20:30 verdict result). The Tong · Wong data DOIs · licences are now on line 805 of the manuscript (figshare 10.6084/m9.figshare.8851955.v6 · Science Data Bank 10.57760/sciencedb.o00133.00038 · both CC BY 4.0).
  - **To check ⟦C⟧:** «the reverse of healthy white matter» in the implications paragraph is used only after it is confirmed against the source text of Lampinen 2019.

### 26-1. Tong 2020 Technical Validation access · comparison [author decision ⑤] (addendum 2026-11-03 21:01:10 KST · access 2026-11-03 21:00:41 KST · PMC7253426 · no reanalysis)

- **Access:** full text of the paper (PubMed MCP · PMC7253426 · doi:10.1038/s41597-020-0493-8). Numbers inside figures were not extracted and were not seen. Recorded at line 26 of x14/EXTERNAL_READ_LOG.jsonl.
- **Contents of the Technical Validation:**
  - Image quality: SNR of b0 · ghost-to-signal ratio (GSR) · background noise of all DWI (per b value) · head motion (RMS).
  - Diffusion model results: TDI maps from MRtrix3 fibre tracking (ACT · SIFT · 100 × 10⁴ (1 million) streamlines) compared visually between centres.
  - There are no inter-centre reproducibility statistics for NODDI or other indices.
- **Comparison:** our Tong T-B (x 0.931 · ODI 0.934 at retest level · y 0.844 exceeded) is neither a reanalysis of results reported by the provider nor a result in conflict with them. The provider's TV has no results at the index level.
- **Facts to use in the manuscript:**
  - All ten scanners are 3 T MAGNETOM Prisma · the same type of 64-channel head coil · the same acquisition (TR/TE 5.4 s/71 ms · GRAPPA 2 · SMS 3 · b 1000/2000/3000 with 30 directions each).
  - Software is VD13D on 9 scanners · VE11C only at centre 09 (agrees with our sensitivity excluding c09).
  - The provided data went through MRtrix3 denoising · Gibbs unringing · TOPUP · EDDY.
  - → The data of the same model · protocol cell are denoised data. Their processing differs from ON-Harmony (no denoising), so this is stated as a limitation in the noise floor discussion.
- **Citation memo:** according to the Usage Notes, 8 centres (01 – 07 · 10) were analysed in an earlier paper by the same research group (centres A – H). That paper was not opened ⟦C: check the title · whether index reproducibility is reported⟧. The relation to the Tong 2023 conference abstract in the manuscript («ODI was highly sensitive to the scanner») is also left as ⟦C⟧.

## 27. b2000 floor — direct reading of the sealed raw images [author decision 2026-11-03 · review session addendum «definition of raw-image b2000»] (addendum 2026-11-03 21:23:18 KST · before reading · before computing · fixed at the time of this addendum)

- **[author decision] (received as it is · replaces «(ii) after 11/9» of decision ④):**
  - Scope: card level 8, first sessions 109. The same as «first sessions 109» in the earlier records. It is the session list 110 minus the processing failure sub-16793 ses-OXF4GEP001, and the main-analysis exclusion sub-16794 ses-NOT4GEP001 is included. Per card level: NOT1ACH 20 · NOT2ING-89 4 · NOT2ING-98 16 · NOT4GEP 10 · OXF1PRI 20 · OXF2PRI 20 · OXF3TRI 10 · OXF4GEP 9. The code counts these numbers again before reading.
  - Files read (three per session): `<ONH_SEALED>/<sub>/<ses>/dwi/<sub>_<ses>_dir-AP_dwi.nii.gz` · `<ONH_META>/…_dir-AP_dwi.bval` · `…_dir-AP_dwi.json`. PA images · bvec · other images are not opened. The only volumes used are b ≤ 50 (b0) and |b − 2000| ≤ 100 (b2000).
  - Record: for each file, in x14/ONH_NEW_READ_LOG.jsonl, the tag «unsealed — post hoc exploratory · verdict unchanged #34» · the basis commit of guard.require() · sha256 · bytes · script. Moved to the canonical record at the next record commit.
  - CPU only (nice 19 · at most 12 concurrent processes). The time limit is 11/5 12:00 KST. If short of time, report only the three Siemens (OXF1PRI · OXF2PRI · OXF3TRI · 50 sessions).
  - Not used to re-evaluate the judge. The judge is not run again with the measured L (author decision · prohibited item).
- **Step 0 (finished · names · sizes only):** processed 4D DWI is not in the 109 session working folders (prune_onh). A name search of all of data finds only sub-13192 OXF1PRI001 · 002 (eddy_out · cudimot/data) from the 10-26 eddy timing measurement and dummy data (eddy_probe). Not used.
- **Step 1 (headers only · new read · voxel values 0):** a session table is made from the sidecars · bval · NIfTI headers.
  - Columns written: number and position of b0 · b2000 volumes · matrix · voxel size · data type · scl_slope/scl_inter · scale fields of the sidecar (Philips*Scale* etc., as they are if present).
  - Acquisition columns: PartialFourier (or a field with the same meaning) · ParallelReductionFactorInPlane · MultibandAccelerationFactor · ReceiveCoilName · CoilCombinationMethod · ImageType (Siemens DIS2D · DIS3D · ND) · on-scanner geometry correction fields (if present) · Manufacturer · model · SoftwareVersions.
  - On-scanner gradient nonlinearity correction: Warrington 2025 wrote that it is off for Siemens · on for the others. If the sidecar has the relevant field, its value is written as well.
- **Main index (per session):** for each voxel v of the ventricle mask V, r_v = mean_k M_vk / SD_k M_vk (ddof 1 · k = b2000 volumes · M = value with the header scale applied) → R = median_v r_v.
  - Voxels with a value of 0 or below in the b0 or b2000 volumes and voxels with SD = 0 are dropped and their number is written.
  - Written as well: |V| · quartiles of r_v · mean B0 of V · mean D2 · mean M².
- **Ventricle mask (raw image space · b0 and b2000 volumes only · b2000 values are not used for the mask):**
  - B0 = mean of the b0 volumes. Head H = the largest 26-connected component of {B0 > 0.2 × percentile 99 of the whole image} · holes filled.
  - Coordinates are converted to mm with the NIfTI RAS transform (sform if present). c = centroid of H, length per axis = percentile 98 − 2 of the H coordinates.
  - Central box C: |R − c_R| ≤ 0.25 L_R · |A − c_A| ≤ 0.35 L_A · −0.20 L_S ≤ S − c_S ≤ 0.40 L_S.
  - Local contrast Q = B0 / Ŝ, Ŝ = G(B0·H)/G(H) (Gaussian σ 15 mm · divides out the coil sensitivity bias). Candidate K = H ∩ C ∩ {Q ≥ 1.5}.
  - V = K eroded by 1 voxel with a 6-connected structure. If |V| < 50, «mask insufficient» is written and that session gives no value.
  - [A's fix · from the review session's definition] The review session's example «lower part of the b1000 mean» does not agree with the author decision (b0 · b2000 only) and is not used. Replacing it with the b2000 mean would create the bias of selecting and measuring voxels from the same data, so only b0 is used for the mask. Tissue mixing is reduced by erosion, and its size is measured by putting partial volume into the truth test.
- **Secondary indices (descriptive):**
  - S1 inter-voxel ratio: for each b2000 volume mean_V / SD_V → median over volumes.
  - S2 σ ratio: σ_floor / σ_SNR. σ_floor = √(mean_V M² / (2 L̂)), and L̂ is the effective L obtained through the calibration table of that card level (only for manufacturers where it can be measured).
    - σ_SNR = mean_W of raw B0 / session SNR. W = nodif_brain_mask ∩ f2_FA ≥ 0.3, and the processing grid is kept as the raw image grid as it is.
    - These three (nodif_brain_mask · f2_FA · snr.json) are re-reads and are written in REREAD_LOG.
    - Session SNR is the value obtained from dti_FA. dti_FA has never been read, so f2_FA is used in its place. As a mask made from processing results, it carries the R21 «exploration after results» tag.
- **Step 2 — truth test (voxel reads 0 · condition for Step 3 · per card level · reported grouped by manufacturer):**
  - Simulated session: an analytic head model is stamped onto the actual grid of the card level (matrix · voxel size · number of b0 · number of b2000 from Step 1) with 2 × 2 × 2 subsampling. Partial volume arises in the signal before noise.
    - Brain ellipsoid: semi-axes (68, 82, 58) mm × U[0.92, 1.08] · centre offset from the field-of-view centre by N(0, 5 mm) per axis.
    - Lateral ventricles: left and right bodies (semi-axes 6 · 28 · 9 mm · centres ±7.5 · +5 · +18 mm · diverging anteriorly by ±10°) and atria (semi-axes 7 · 12 · 8 mm · centres ±16 · −22 · +10 mm). The size factor is U[0.7, 1.5].
    - Third ventricle: semi-axes 2 · 12 · 8 mm · centre 0 · 0 · −2 mm · factor U[0.7, 1.3].
  - Signal:
    - Tissue b0 = 1 · CSF b0 = κ ~ U[1.6, 4.0].
    - Tissue b2000 = exp(−2000 · ADC_k): ADC_k = λ⊥ + (λ∥ − λ⊥)(g_k · u(x))² · λ∥ 1.7 · λ⊥ ~ U[0.3, 0.7] (× 10⁻³ mm²/s) · fibre direction u(x) rotates smoothly with location.
      - The proportion of isotropic tissue (ADC 0.8 × 10⁻³) is U[0, 0.5].
      - Directions g_k are N2 points uniform on the sphere (bvec is not read).
    - Remaining CSF signal: κ · exp(−2000 · D_c) · D_c ~ U[2.6, 3.2] × 10⁻³.
    - Pulsation: CSF signal (b0 · b2000) multiplied per volume by (1 + ε_k), ε_k ~ N(0, p²) · p ~ U[0, 0.3].
  - Motion: one actual session of the same card level is picked at random, and its eddy_parameters (translation 3 · rotation 3) are used in volume order as rigid transforms. This is a re-read of files read in D5 and is written in REREAD_LOG.
  - Noise:
    - Noncentral chi M = √((s + σz)² + σ² · Gamma((2L − 1)/2, 2)). For L = 0.5 it is |s + σz|.
    - σ = 1/SNR_t, varying smoothly by ±30 % over space (no effect on the within-voxel ratio).
    - SNR_t ~ U[minimum recorded SNR of the card level, 2 × maximum]. This is because the recorded SNR is a value measured by retest and is lower than the thermal-noise SNR (snr.json re-read).
  - L grid: 0.5 · 1 · 1.5 · 2 · 3 · 4 · 6 · 8 · 12. Noise is added separately per L to the same noiseless signal.
  - On-scanner resampling of Philips · GE: {none · linear · order-3 spline} × sub-voxel shift (per axis U[0, 0.5] voxel + a gentle location-dependent gradient) is applied to the magnitude images after noise. The three choices are calibrated separately, and the verdict uses the range combining all three. For Siemens it is none.
  - Quantization (side check): the effect at storage-unit σ of 3 · 1 is written separately. It is flagged if the measured σ_stored is below 10.
  - Size: per card level 200 simulated sessions (L 9 on the same noiseless signal · three resamplings for Philips · GE). The same mask · index code as for the measurement is used.
  - Output: per card level × L (× resampling), median of R · percentile 2.5 – 97.5 · proportion with mask insufficient · S1.
  - **Separates (measurable):** R percentile 2.5 of L 8 > R percentile 97.5 of L 4. Philips · GE are seen with the three resampling choices combined (97.5 of the highest L 4 < 2.5 of the lowest L 8).
    - If all card levels of a manufacturer fail to separate, that manufacturer is reported as «not measurable» and is not interpreted.
    - The [A's proposal] «MC SD 4 times» criterion of §25-2 is replaced by this.
  - **Minimum difference Δ_min (per scanner pair of the main comparison):**
    - Scanner value = median of the R of that scanner's sessions. A simulated scanner value is obtained by sampling with replacement from the bank as many as the actual number of sessions (sessions with a valid mask) (2,000 times · seed fixed).
    - At the same L (pooling 1 · 2 · 4 · 8), percentile 97.5 of |difference of the two scanner values| is Δ_min.
    - Power = proportion with difference ≥ Δ_min under the simulator assumption (OXF2PRI L 8 · OXF1PRI · OXF3TRI L 4). If below 0.95, the main comparison is «not measurable (the truth test cannot separate)».
  - The truth test results (calibration table · Δ_min · power) are appended to this plan and fixed before voxel reads (§27-1).
- **Step 3 — voxel reads · prediction before computing and verdict rule (review session · set in advance):**
  - Main comparison (confirmatory in character): OXF2PRI (64 channels) vs OXF1PRI · OXF3TRI (32 channels). d₁ = R(OXF2PRI) − R(OXF1PRI) · d₂ = R(OXF2PRI) − R(OXF3TRI) · d₃ = R(OXF1PRI) − R(OXF3TRI) (32 vs 32 control).
    - Simulator assumption (L 8 vs 4) side: d₁ ≥ Δ_min and d₂ ≥ Δ_min.
    - (A) side (the L difference assumed by the simulator is absent or small): |d₁| · |d₂| · |d₃| are all smaller than their respective Δ_min.
    - Anything else is «cannot separate».
  - Remaining scanners (descriptive): per-scanner R · session bootstrap 95 % interval · effective L range through the calibration table. Conclusions are not mixed between manufacturers.
  - Pulsation · motion · partial volume lower the ratio (toward small L · toward our hypothesis), and remaining signal raises the ratio. So absolute L is used only as descriptive.
  - The relation to the FWF level gap (0.080 vs 0.419) is only descriptive («consistent with» level).
- **Sensitivity (R21 «exploration after results»):** the same R is obtained with a mask of f2_FWF ≥ 0.95 on the processing grid ∩ nodif_brain_mask ∩ central box, eroded by 1 voxel.
  - The processing grid is the same as the raw image grid (no reorientation · cropping), and the reference volume of eddy is the first b0.
  - The only stored transforms are the eddy text logs. The topup field was deleted by prune_onh, so the distortion cannot be undone. So the grids are matched as they are. This point is written as well.
  - f2_FWF · nodif_brain_mask are re-reads (D5).
- **Report:** per-session R · per-scanner median and interval · per-manufacturer truth test calibration table · number of mask voxels · sensitivity · scanners not measurable and why.
- **Manuscript:** one paragraph in Results ⑥ · one supplementary table · revision of the limitations paragraph. Deep analysis across all scanners is for a follow-up paper.
- **Code · outputs:** drafts/x15/b2000/ (b2k_core.py mask · index · b2k_headers.py Step 1 · b2k_sim.py Step 2 · b2k_read.py Step 3). Outputs are in explore/b2000/. For each script, its sha256 is written in this plan before it is run.

### 27-0. Record before running Step 1 (addendum 2026-11-03 21:25:37 KST)

- b2k_core.py sha256 47191b21b572cdae… · b2k_headers.py sha256 a5398abf8adece9b… — values before running. Noise sampler check (reads 0): mean/SD at L 0.5 · 1 · 2 · 4 · 8 1.323 · 1.912 · 2.754 · 3.945 · 5.609 (theory 1.324 · 1.913 · 2.755 · 3.943 · 5.615 · samples 2 million).

### 27-1. Step 1 result · mask flaw found during the truth test — fix before reading (addendum 2026-11-03 21:40:16 KST · raw image voxel reads 0)

- **Step 1 (headers only · 21:25:37 – 21:30:18 · files 327 · record ONH_NEW_READ_LOG):** all 109 sessions have b2000 volumes 50. b0 are 5 – 7, and voxels are about 2 mm.
  - Three Siemens: the sidecar CoilCombinationMethod is «Sum of Squares» in all (Head_32 · HeadNeck_64 · 32Ch_Head). ImageType is ND (no on-scanner geometry correction), and ImageComments is «Unaliased MB3/PE3 SENSE1». Which coil combination was actually used is measured with the floor. The prediction stays as in §27.
  - Philips: SENSE acceleration 2 · 1.5. The storage scale (PhilipsRescaleSlope · ScaleSlope) differs by session, but the ratio does not change under a constant multiplier.
  - GE: ASSET 2 · Hyperband 2.
  - The on-scanner geometry correction field (NonlinearGradientCorrection) is absent from all sidecars. Philips · GE follow «on» of Warrington 2025.
- **Flaw:** in one simulated session for time measurement (OXF1PRI · SNR 16 · L 6 · 8), the head mask threshold (0.2 × percentile 99 of the whole image) took the background floor as head. So the box grew and mask voxels 819 · 1,129 went outside the brain. This is because with large L the background floor rises in the form σ√(2L).
- **Fix:** the head mask threshold is changed to Otsu(B0). The rest of the mask definition stays as it is.
  - After the fix, in three simulated sessions of the same card level, for all of L 0.5 – 12 no mask voxel went outside the region (0).
  - |V| hardly changed with L (1,517 – 1,586 · 346 – 354 · 2,541 – 2,580).
  - R was L 4 3.93 – 3.98 · L 8 5.65 – 5.72. The pure-noise ideal values at N2 50 are 3.97 · 5.65.
- **Raw image voxels have not yet been read.** This fix deals, before reading, with a flaw found by the truth test (deviation record).
- **Bank run:** b2k_core.py sha256 42a619075ca91d2d… · b2k_sim.py sha256 30047ba8a8b3fc7d… · card level order three Siemens → three Philips → two GE · 12 concurrent processes · nice 19 · seed [card level number, simulation number, 20261103].

## 28. Answer to the evening report [author decision 2026-11-03 night] · review session reply «11-03 evening report» (addendum 2026-11-03 22:53:04 KST)

- **[author decision] (received as it is):**
  - Code: GitHub is made public at submission, and the DOI of the Zenodo archived version is used.
    - Check before going public: confirm that passwords · access information · personal information · other users' names · sealed raw data · absolute paths are not in the repository history (with the 11/7 statements · planned separately).
  - Derived ROI values: made public together. The licence conditions of each dataset are followed and the source is stated.
  - Ethics sentence: ⟦the author fills in whether the university requires wording⟧.
  - [omitted]
    - [omitted]
    - [omitted]
  - Code reviewer: ⟦the author fills in⟧. If there was no human line-by-line review, it is written truthfully as «relied on truth tests · the automated test suite · in-project review».
  - The mutation test sentence is dropped (A's proposal accepted).
  - Conclusion section: about 120 words · included in the main text total 8,820 (A's proposal accepted).
- **Review session reply 1 (received · to be added):**
  - Since all ten centres are Prisma 64-channel, agreement within the same protocol is agreement between the same coil. One sentence saying that nothing is said about coil differences goes in the Discussion (11/5 Discussion).
  - It does not overlap with the result of the two Prismas (different coils) and complements it.
  - The provider's noise reduction goes into the limitations paragraph.
- **Review session reply 2 (abstract fixes ① – ④ · written as «three» but four):** put into the working copy (S8's share). The abstract is 249 words (wc_main.py · same as the review session's count). The final abstract is 11/7.
- **Review session reply 3 (main path for Methods compression · 11/5 – 11/6):**
  - Per-section budget (total about 2,350): 3.1 data · statements 350 · 3.2 processing summary 300 · 3.3 mapping from coordinates and theory → tests 150 · 3.4 tests fixed before computing 400 · 3.5 judge 450 · 3.6 post hoc analyses 350 (D1 · raw-value decomposition · equal floor settings · b2000 — with definition times) · 3.7 blinded · exposure · deviation 200 · 3.8 reproducibility 150.
  - Do not polish sentences; move them whole.
  - To move to supplementary methods: processing commands · simulator and judge details · identifiability metrics · Fisher · separation attempts 1 – 3 · fixed analysis document procedure · rule list R1 – R26.
  - Introduction 1,342 → 900 (keeping the question · contribution paragraphs) · Theory 2,351 → 1,600 (proofs to Appendix A · generality section in one or two sentences) · Discussion 2,040 → 1,500 (starting from the places for Box 1 · Box 2 · limitations about 260).
- **Review session reply 4 (to attach to the 11/4 noon report):** all of references.bib · Wong straight-line ROI unit proportions (x · y · R12) · y · ODI values and intervals of the two Prismas 540 ROI units · Table 3 body · b2000 truth test results and the §27-2 verdict rules (fixed time). Start reading with the Siemens 50 sessions, and do the rest within the approved scope and the 11/5 12:00 limit.
- **[A's fix · a fact in the report]** The evening report gave the second place where the Tong TV access time was put as «3.9». The actual place is **5.5 Limitations item 1** («Unblinded checks and disclosed exposures»). The first place is the blinded table of 3.8.

### 27-2. Siemens truth test results fixed · verdict rule addendum [A's proposal · awaiting author approval] (addendum 2026-11-03 23:49:22 KST · raw image voxel reads 0 · fixed at the time of this addendum · applies only to the three Siemens OXF1PRI · OXF2PRI · OXF3TRI)

- **Basis:** in-project review (not external peer review). It went through one review (b2k_core · b2k_sim cross-check · blocker 0) and one workflow (three reviews of code · statistics · physics → refutation checks 9 → synthesis · verdict «proceed after fixing»). This addendum takes the recommendation of that synthesis and only widens the L range to 128 (because R ≈ 16.1 of 64-channel white SoS lies on the 64 boundary).
- **(i) Fixed inputs (full sha256 · result() stops if they differ · raise, not assert):**
  - b2k_core.py 42a619075ca91d2d08f4d1ad3cefbe529defa104bf444c582c5a88a8e7deea74 (the version used by the bank as it is)
  - b2k_sim.py 30047ba8a8b3fc7d6eab301c902b972907b09acd4d7af944d6bf56baf0793d5f
  - b2k_summary.py a9e2220a0f1bb33e11ceea6337980bb448ee883e15723dd23183c5d02f9d40e1
  - b2k_read.py 947221c00ea38df143d0f5d2dce03913ceaaa01bd71dd64aefd83ffd44e21c5f
  - headers.json 806180efe1369e0b39d62067e03a388fedc32001c13ebe545cac4bc276679380
  - kat_bank_siemens.jsonl b04035ab8911fc1223a76f2e5d088a370e23b188a3297c056f3530d082d8234b (Siemens 5,400 lines = simulations 600 × L 9 · a fixed copy taken out of kat_bank_all.jsonl)
  - kat_summary_siemens.json bacdb8e59945c539a36b011cfa658c82d8caaa01e45126a643c1f93d99d3a7c7
  - pure_curve_N2_50.json 42b0027136ca2e536cd2afb23b8453f3d225de1bda931401ee107d57d948cc81
  - drafts commit: the commit containing this addendum. The hash is written in the approval request and the receipt. b2k_read opens sealed files only when x15/b2000 and this plan are the same as in that commit.
- **(ii) Truth test results (what §27 required to be fixed before reading):**
  - Calibration table (R median [percentile 2.5, 97.5] · simulations with a valid mask 185 – 196):

| L | OXF1PRI | OXF2PRI | OXF3TRI |
|---|---|---|---|
| 0.5 | 1.335 [1.323, 1.361] | 1.334 [1.319, 1.350] | 1.333 [1.315, 1.346] |
| 1 | 1.927 [1.914, 1.951] | 1.925 [1.905, 1.941] | 1.922 [1.901, 1.943] |
| 2 | 2.772 [2.752, 2.790] | 2.770 [2.744, 2.791] | 2.770 [2.737, 2.795] |
| 4 | 3.969 [3.928, 3.995] | 3.965 [3.930, 3.995] | 3.962 [3.917, 3.997] |
| 8 | 5.655 [5.598, 5.694] | 5.646 [5.593, 5.683] | 5.649 [5.568, 5.696] |
| 12 | 6.941 [6.884, 6.988] | 6.935 [6.876, 6.986] | 6.933 [6.860, 6.993] |

  - Separates: all three card levels separate (percentile 2.5 of L 8 5.568 – 5.598 > percentile 97.5 of L 4 3.995 – 3.997).
  - Δ_min (§27 definition · nominal 20/20/10 · pooling L 1 · 2 · 4 · 8):
    - OXF2PRI − OXF1PRI 0.0180 · OXF2PRI − OXF3TRI 0.0167 · OXF1PRI − OXF3TRI 0.0208
    - Power 1.00 · 1.00 (median of the L 8 vs 4 difference 1.678 · 1.684)
    - Δ_min per L (report only): L 1 0.008 – 0.012 · L 8 0.022 – 0.027
  - Quantization side check (§27 · the bank has no rounding · in-project numbers with integer rounding put into pure noise N2 50 · A recomputed L 1 · σ 3: excluded 50 % · R +0.037 · with 0 included R −0.02):
    - L 1 · σ_q 2 / 3 / 5 / 10 / 20: voxels excluded 79 / 50 / 22 / 6 / 1.6 % · R bias +0.10 / +0.043 / +0.015 / +0.004 / +0.001.
    - L 4 · σ_q 2 / 3 / 5 / 10: excluded 0 % · R −0.08 / −0.036 / −0.011 / −0.002.
    - L 0.5: at σ_q ≤ 3 almost all excluded («mask insufficient»), at σ_q 10 · 20, 87 · 64 % excluded.
- **(iii) Session values:**
  - Calibration curve g_u(L) = P(L) · b_u(L).
    - P: median voxel ratio of pure noise N2 50 (L 0.5 – 128 · log-spaced 47 points · voxels 10⁵ per point · seed 20261103).
    - b_u: bank median / P (grid L 9 points · log L interpolation · for L > 12, b_u(12)).
  - ℓ = log2 g_u⁻¹(R) (log L linear interpolation). Outside the range, the clamped end value is not used as a point and a censored interval is used.
    - R < g_u(0.5): «censored below» · ℓ ∈ (−∞, −1].
    - R > g_u(128): «censored above» · ℓ ∈ [7, ∞).
  - Quantization: σ_q = √(mean_M2 / (2 L_eff,session)). Sessions with σ_q < 3 are dropped from both ℓ and the scanner value and their number is written. σ_q < 10 is only flagged.
- **(iv) Contrasts (common subjects · first sessions):**
  - c₁ = ℓ(OXF2PRI) − ℓ(OXF1PRI) · c₂ = ℓ(OXF2PRI) − ℓ(OXF3TRI) · c₃ = ℓ(OXF1PRI) − ℓ(OXF3TRI).
  - Each c = med_s ℓ_A − med_s ℓ_B, where s are the subjects with values on both scanners. Minimum common subjects are c₁ 15 · c₂ 8 · c₃ 8.
  - Verdict interval: sampling with replacement at the subject level 2,000 times (same subject numbers on both scanners · seed [20261103, 99] · in the order c₁ → c₂ → c₃).
    - Lower end = percentile 2.5 of (med ℓ_lo(A)* − med ℓ_hi(B)*) · upper end = percentile 97.5 of (med ℓ_hi(A)* − med ℓ_lo(B)*).
    - Without censoring this is the same as the ordinary percentile interval. ±∞ is not interpolated.
  - The paired-difference sign test interval is reported only as sensitivity. If it disagrees with the verdict, that fact is written, and the verdict follows the bootstrap interval.
- **(v) Verdict (checked in order, stopping at the first that applies · ½ = log2 √2 · ⅓ is a ratio of about 1.26 times):**
  1. Too few common subjects → «cannot separate (too few pairs)».
  2. An end of the c₁ · c₂ · c₃ intervals is infinite → «cannot separate (censored outside the calibration range)». Censoring counts and directions per scanner are written.
  3. Gate G: the c₃ interval is not within the closed [−½, ½] → «32-channel control disagreement — not read as a channel-count effect». In this case, for c₁ (two Prismas · syngo E11) the label obtained with the same threshold is written only as «descriptive · one pair within Prisma».
  4. Lower ends of the c₁ · c₂ intervals are both ≥ ½ → «simulator assumption side (64/32 effective L ratio ≥ √2)» · c* = 1.
  5. c₁ · c₂ intervals are both within the open (−½, ½) → «(A) side (64/32 effective L ratio 1/√2 – √2)» · c* = 0.
  6. Otherwise → «cannot separate». If the upper end of c₁ or c₂ < −½, the flag «opposite direction (64-channel effective L is small)» is added.
  - Size (when 4 · 5): if the c₁ · c₂ intervals are both within the closed [c* − ⅓, c* + ⅓], «size agrees», otherwise «size differs».
- **(vi) Absolute level (descriptive · not a verdict · §27 «absolute L as descriptive only»):**
  - a_u = log2(L_eff,u / L_assumed,u). L_assumed is OXF1PRI · OXF3TRI 4, OXF2PRI 8.
  - The interval is the session bootstrap interval of the scanner median R (2,000 times · seed [20261103, card level order number]) carried over. Censored ends are ±∞.
  - Labels (in order): within [−½, ½] «assumed level (within 1.41 times)» · upper end < −1 «more than 2 times lower than assumed» · lower end > 1 «more than 2 times higher than assumed» · upper end < −½ «lower than assumed (1.41 – 2 times)» · lower end > ½ «higher than assumed (1.41 – 2 times)» · otherwise «on the boundary».
  - Conclusion sentence rule:
    - Even if item (v) is «simulator assumption side», if even one scanner is not at «assumed level», the same sentence states «the ratio agrees with the assumption but the effective L level differs from the assumption» and the per-scanner L_eff · intervals.
    - If item (v) is «(A) side», L_eff and a_u of the three scanners are written together.
- **(vii) §27 literal rule (supplementary · not interpreted):**
  - Δ_min and power are computed again with the actual session numbers, and the rule is computed and reported as it is.
  - «not measurable» gate: if the smaller power of the two pairs is below 0.95 or the Siemens do not separate, «not measurable (the truth test cannot separate)»; if there is a Siemens scanner with no value, «not measurable (a Siemens scanner has no value)».
  - Whether or not the spread flag is present, this verdict is not interpreted as evidence about L.
  - Reason (simulation · in-project numbers): Δ_min is 1.07 % of the expected difference. With a scanner-level R bias τ of only 0.5 %, at the same L, P((A)) 0.20 · P(simulator side) 0.13. If the measured spread is 2 times that of the simulation, P((A) | same L) 0.67.
- **(viii) Descriptive (not used for the verdict · list set in advance):**
  - Per-scanner R · session bootstrap interval · L_eff (censoring flag).
  - Spread check: the interquartile range of the scanner's session R is divided by the simulated interquartile range at the grid L nearest in log distance to that scanner's L_eff. Flagged if above 3.
  - Number of sessions flagged by σ_q and number dropped.
  - Diagnostics: R with 0 included · exclusion proportion · R of the inner part after one more erosion vs rim · R with volumes whose volume mean has |z| > 5 (MAD) dropped.
  - Floor level: mean_V M_b2000 / S_W(B0).
  - S2 · FWF sensitivity mask (R21 «exploration after results»).
  - Physical expectation («consistent with» level · synthetic numbers from the review · the actual Siemens reconstruction method was not checked):
    - SENSE1 + PF zero-filling: R ≈ 1.92 – 1.93 (L ≈ 1)
    - Per-channel POCS + SENSE1: R 1.75 – 1.92
    - Single-image POCS: R ≈ 1.59
    - homodyne (real part): R ≈ 1.36 (near the lower censoring boundary)
    - SoS zero-filling: 32-channel R ≈ 9.8 (L ≈ 24) · 64-channel R ≈ 13.4 (L ≈ 43)
    - SoS + POCS: R ≈ 7.9 · 10.9
  - Limits of interpretation (set in advance):
    - Sidecar «Sum of Squares» is the field the CMRR guide says to select in order to use SENSE1 (CMRR C2P R017pre15 p. 18 · checked in the review). From ImageComments «SENSE1» a SENSE1 combination is inferred.
    - So this test measures the reconstruction of these data, not coil physics. The manuscript wording is limited to «in these SENSE1 reconstruction data» and is not widened to 64-channel coils in general.
    - If the simulator was meant to imitate SoS data, this test cannot refute that coil dependence.
- **(ix) Manuscript:**
  - The main result is the verdict of item (v) and the conclusion sentence of item (vi). The supplementary table is item (vii).
  - This designation is set by author approval before reading. It replaces line 41 of AMEND_27-2_draft.md (sha256 548716df…) «which reading to use is set by the author».
  - The ρ/√2 size reading M of AMEND_27-2_draft.md §3 and the review session patch b2k_summary_M.patch are «superseded — not used».
- **(x) Deviation records (all four set after seeing truth test result (ii) · before reading voxels):**
  - ① The §27 verdict rule is moved down to supplementary, and the main verdict is item (v) (M′).
  - ② Exclusion of sessions with σ_q < 3 (§27 set only the flag σ_stored < 10).
  - ③ The calibration range is widened from the L 0.5 – 12 grid to the pure curve 0.5 – 128, and outside the range censored intervals are used.
  - ④ Descriptive diagnostics are added (R with 0 included · inner vs rim · volume outliers · floor level). Not used for the verdict.
- **Philips · GE:** the bank is still being appended. Before reading those voxels, the calibration table of that card level is fixed separately in the same way (fixed copy · sha · cross-check) (§27-3). Until then only the Siemens are read.
- **Tests (voxel reads 0 · fake sessions drawn from bank values · all 19 passed):**
  - H0 (L 1 · L 4) → «(A) side · size agrees»
  - H1 (8/4/4) → «simulator assumption side · size agrees»
  - Ratio 1.5 → «simulator assumption side · size differs»
  - Trio L 2 → «32-channel control disagreement»
  - Opposite direction → «cannot separate»
  - Too few pairs · σ_q 2 → «too few pairs»
  - Five censoring cases (all three below · 64 L1 and both 32 below · only 64 below · 64 above · Trio 2 sessions below) → «censored outside the calibration range»
  - 64 1 session only below → verdict unchanged
  - A truncated last line is skipped · a broken line in the middle · tampering with the calibration table · pure curve · bank copy stops

### 27-2a. §27-2 approval [author decision 2026-11-03 23:59:08 KST · answered a chat question with «approve · read now»]

- What was approved:
  - (i) §27-2 body and fixed time (23:49:22 KST · drafts commit c3cf29875cb06bf014a53110b3ede32bf005e6a8)
  - (ii) Manuscript main result = M′ verdict (item (v)) + absolute level sentence (item (vi)) · §27 literal rule = supplementary table (not interpreted)
  - (iii) Voxel reads of the Siemens 50 sessions
- Operation: to keep to at most 12 concurrent processes (§27), the Philips · GE bank (PID 107486 and workers) is paused with SIGSTOP for a while. When reading finishes it resumes with SIGCONT. With deterministic seeds the results do not change.

### 27-4. Siemens results — post hoc exploratory · verdict unchanged #34 (addendum 2026-11-04 00:04:07 KST · reads 2026-11-04 00:01:33 – 00:01:50 KST · 50/50 session values · failures 0)

- **Reads:** b2k_read.py (sha256 947221c0…) · drafts HEAD b7a72154 · concurrency 6 (bank paused with SIGSTOP · SIGCONT after reading).
  - Record: ONH_NEW_READ_LOG 150 lines (50 × NIfTI · bval · sidecar · 9 keys · «unsealed — post hoc exploratory · verdict unchanged #34»).
  - Re-reads: REREAD_LOG 200 lines (50 × f2_FWF · nodif_brain_mask · f2_FA · snr.json).
  - Outputs: read_sessions.jsonl (sha256 141f23425081a341…) · result b2k_result_siemens.json (sha256 86d848717d22d566… · produced by b2k_summary.py a9e2220a…).
- **Main result (§27-2 (v) · (vi) · rule approved by the author):** «(A) side · size agrees». Gate G was passed.
  - c₁ (Prisma 64 vs 32 · common 20 subjects) = +0.038 [−0.002, +0.076] → effective L ratio 1.03 [1.00, 1.05]
  - c₂ (Prisma 64 vs Trio 32 · common 10 subjects) = +0.209 [+0.155, +0.278] → 1.16 [1.11, 1.21]
  - c₃ (Prisma 32 vs Trio 32 · common 10 subjects) = +0.162 [+0.122, +0.233] → 1.12 [1.09, 1.17]
  - Censored 0 · sessions dropped by σ_q 0. The paired sign test interval is on the same side as well: c₁ [−0.004, +0.080].
  - Absolute level (item (vi)): effective L is OXF1PRI 1.03 [1.02, 1.05] · OXF2PRI 1.06 [1.03, 1.09] · OXF3TRI 0.93 [0.89, 0.96].
    - a_u = −1.95 · −2.92 · −2.11 (log₂), all three «more than 2 times lower than assumed» (assumed 4 · 8 · 4).
- **Scanner R (session median [bootstrap 95 %] · range):**
  - OXF1PRI 1.960 [1.949, 1.979] · 1.899 – 2.135
  - OXF2PRI 1.986 [1.957, 2.014] · 1.834 – 2.154
  - OXF3TRI 1.852 [1.814, 1.880] · 1.784 – 1.951
  - «consistent with» level with the physical expectation (SENSE1 + PF zero-filling ≈ 1.92 – 1.93).
- **§27 literal rule (supplementary · not interpreted):** «simulator assumption side (L 8 vs 4)».
  - d₁ = +0.026 ≥ Δ_min 0.018 · d₂ = +0.134 ≥ 0.017 · d₃ = +0.108. Power 1.00.
  - As pointed out before reading, the result only separates the sign (expected d for L 8 vs 4 ≈ 1.68).
- **Descriptive:**
  - Spread check: the measured interquartile range is 3.5 · 6.6 · 5.3 times the simulated, and all three were flagged. The truth test does not capture the measured spread.
  - σ_q median 136 · 415 · 13.8 (minimum 72 · 190 · 11.5). Flagged 0 · dropped 0.
  - Mask voxel median: 162 [84, 620] · 171 [96, 536] · 408 [123, 965].
  - **Proportion of voxels dropped for values 0 or below:** 10.6 % [2.5, 21.7] · 6.2 % [0.4, 18.9] · 27.4 % [23.6, 34.6].
    - There are far more 0 values than pure Rician at σ_q ≥ 11 predicts. The cause is unknown (descriptive).
    - R with 0 included is 1.945 · 1.966 · 1.822.
    - M′ computed again with R including 0 is also «(A) side · size agrees» (c₁ +0.033 [−0.003, +0.076] · effective L 1.02 · 1.04 · 0.90 · descriptive sensitivity).
  - Inner vs rim: with one more erosion most are below 50 voxels and have no value (only OXF3TRI 2 sessions 1.869 vs 1.849).
  - Volume outlier flags: OXF1PRI 1 session (1 volume) · OXF3TRI 2 sessions (3 volumes).
  - Floor level (mean_V M_b2000 / S_W): 0.069 · 0.104 · 0.076 (median).
  - FWF mask sensitivity (R21): R 1.944 · 1.930 · 1.828 · Dice median 0.46 · 0.40 · 0.47.
  - S2 σ ratio (R21): 1.22 · 1.69 · 0.76.
- **Limits of interpretation (as set in advance in §27-2 (viii)):**
  - This test measures the floor of these SENSE1 reconstruction data, not coil physics.
  - The simulator's L 4 / 8 assumption did not hold for these data (L ≈ 1 on all three scanners).
  - Not widened to 64-channel coils in general. Not used to re-evaluate the judge (prohibited item).

### 27-3. Philips · GE truth test results fixed (addendum 2026-11-04 00:21:41 KST · Philips · GE voxel reads 0 · descriptive card levels)

- **Fact about the order:** the bank started at 21:40:16 with fixed code (b2k_sim 30047ba8…) and ended at 00:20:12. Seeds are deterministic.
  - The Siemens results (§27-4) were seen before this addendum.
  - However, Philips · GE were set in advance in §27 as «descriptive card levels». The measurable rule (§27 separates · resampling choices combined) also stays as it is. This addendum adds no other choice.
- **Fixed inputs (full sha256):**
  - kat_bank_PG.jsonl 69aec62dd11d072e3529360bec550a8ac424b1e577cdb74e1df375f01af4d670 (27,000 lines = 5 card levels × simulations 200 × L 9 × resamplings 3 · taken out of kat_bank_all.jsonl · broken lines 0 · duplicates 0)
  - kat_summary_PG.json 8cdfa151a413bedb593c5db67e8babda14d8e94c4a6574b81f540cadd94f5ac1
  - b2k_summary.py 8a1e74d7b7c8fbb99f885900adaca847a8ce729ad4c328a13b77354518c6a06c
    - What changed: Philips · GE fixed copy · sha cross-check · combined L range · not-measurable flag.
    - It reproduces the Siemens results as they are (M′ · scanners · literal rule · status all the same · regression check).
  - b2k_read.py 947221c0… · b2k_core.py 42a61907… stay as they are.
- **Result — separating (percentile 2.5 of L 8 > percentile 97.5 of L 4 · three resampling choices combined):** all five card levels «do not separate».
  - Maximum percentile 97.5 of L 4: 10.12 – 10.64 (linear resampling)
  - Minimum percentile 2.5 of L 8: 5.58 – 5.61 (no resampling)
  - → **Philips «not measurable» · GE «not measurable»** (§27: not interpreted)
- **Calibration summary (R median [2.5, 97.5]):**
  - No resampling: same as Siemens (L 1 1.925 – 1.926 · L 4 3.967 – 3.970 · L 8 5.648 – 5.655).
  - Order 3: L 1 2.34 – 2.39 [2.02, 2.80] · L 4 4.72 – 4.87.
  - Linear: L 1 3.53 – 3.72 [2.53, 5.17] · L 4 7.26 – 7.64.
- **Reads (approved scope · author decision b2000 «first sessions 109»):**
  - Target 59 sessions: NOT1ACH 20 · NOT2ING-89 4 · NOT2ING-98 16 · NOT4GEP 10 · OXF4GEP 9.
  - The output is descriptive: per-scanner R · intervals · effective L per resampling choice and combined range · diagnostics. All carry the flag «not measurable — not interpreted».

### 27-5. Philips · GE results — descriptive · not measurable · not interpreted (addendum 2026-11-04 00:22:56 KST · reads 2026-11-04 00:21:48 – 00:22:03 KST · 59/59 session values · failures 0)

- **Reads · record:**
  - drafts HEAD c5ac0ccd (after §27-3 was fixed).
  - The sealed access record «b2000 Step 3» has 327 lines in all (109 sessions × 3 · 9 keys), and re-reads are 436 lines (109 × 4).
  - Outputs: read_sessions.jsonl (109 sessions · sha256 da860cd3d7a45380…) · b2k_result_all.json (sha256 a2dbfe7c8d3265f0… · b2k_summary.py 8a1e74d7…).
  - The Siemens results (§27-4) stay as they are in the overall summary as well (M′ · literal rule the same).
- **Scanner R (session median [bootstrap 95 %]) · median of voxels dropped for values 0 or below · effective L under no resampling · range combining the three choices:**
  - NOT1ACH 1.905 [1.895, 1.917] · 3.5 % · 0.98 · 0.5 – 0.99
  - NOT2ING-89 2.066 [2.016, 2.169] · 46 % · 1.14 · 0.5 – 1.25
  - NOT2ING-98 2.141 [2.109, 2.173] · 2.3 % · 1.22 · 0.5 – 1.26
  - NOT4GEP 1.618 [1.612, 1.633] · 44 % · 0.72 · 0.5 – 0.73
  - OXF4GEP 1.587 [1.575, 1.602] · 61 % · 0.69 · 0.5 – 0.71
  - (The lower end 0.5 of the combined range is lower censoring under the linear resampling choice.)
- **By the rule:** for Philips · GE the truth test cannot separate L 4 and 8, so they are «not measurable». The values above are used only as descriptive and are not interpreted. Conclusions are not mixed between manufacturers.
- **Descriptive observations (not interpretation):**
  - R of the two GE scanners (1.59 – 1.62) is lower than the Rician value (1.925). Voxels with values 0 or below are 44 – 61 %.
  - Voxels with values 0 or below in NOT2ING-89 are also 46 %.
  - R with 0 included is NOT1ACH 1.900 · NOT2ING-89 1.966 · NOT2ING-98 2.138 · NOT4GEP 1.572 · OXF4GEP 1.525.

### 27-6. 0 values · σ(DN) · S2 check [review session reply «11-04 early-morning report» 2 ① · ② · post hoc descriptive] (addendum 2026-11-04 00:40:16 KST · new image reads 0 · read_sessions.jsonl · headers.json only)

- **σ(DN)** = √(mean_M2 / (2 L_eff,session)) (L_eff from the «no resampling» calibration curve of that card level).
- **Rounding-expected 0 proportion** = the probability under Rayleigh (L 1) that at least one of the b2000 50 volumes has M < 0.5.
- Per-scanner results (session median [minimum, maximum]):

| Scanner | Data type · header scale | σ(DN) | Observed 0-value voxels | Rounding expectation | b0 white matter (DN) | S2 = σ_floor / σ_SNR |
|---|---|---|---|---|---|---|
| OXF1PRI | uint16 · none | 136 [72, 221] | 10.6 % [2.5, 21.7] | 0.0 % [0.0, 0.1] | 1,889 | 1.19 [0.68, 2.20] |
| OXF2PRI | uint16 · none | 415 [190, 1,251] | 6.2 % [0.4, 18.9] | 0.0 % | 3,740 | 1.75 [0.92, 4.30] |
| OXF3TRI | int16 · none | 13.8 [11.5, 33.0] | 27.4 % [23.6, 34.6] | 3.3 % [0.6, 4.6] | 201 | 0.78 [0.48, 1.90] |
| NOT1ACH | float32 · int16 · none | 6,330 [25, 16,702] | 3.5 % [1.8, 65.3] | 0.2 % | 86,982 | 0.94 |
| NOT2ING-89 | float32 · none | 562 | 46.1 % | 0.0 % | 6,647 | 0.96 |
| NOT2ING-98 | float32 · int16 · none | 33 [25, 934] | 2.3 % [0.0, 69.6] | 0.6 % | 413 | 0.70 |
| NOT4GEP | int16 · none | 105 | 44.0 % | 0.1 % | 1,470 | 0.97 |
| OXF4GEP | int16 · none | 158 | 60.5 % | 0.0 % | 2,007 | 1.11 |

- **Reading:**
  - σ(DN) is Trio 11.5 – 33, other scanners 25 – 16,700. It is not in the range σ 0.7 – 1.5 DN of the review session's simulation.
  - The 0 values created by integer rounding are below 1/8 of the observed (Trio 3.3 % vs 27.4 %, Prisma below 0.1 % vs 6 – 11 %).
  - So the 0 values are not the rounding of noise samples. They are likely exact 0 values arising from scanner · conversion processing (clipping · masking · slice-level or volume-level loss) (inference · not checked).
  - At σ ≥ 10 the size by which rounding moves R is ≤ 0.004 (§27-2 (ii) quantization table). A truth test with rounding at the observed σ, or integer-interval maximum likelihood, cannot change the verdict on these data.
  - The σ(DN) difference between the two Prismas (136 vs 415) is in the same direction as the storage scale difference (b0 white matter 1,889 vs 3,740). R is invariant to scale. The 0-value proportion is 10.6 % vs 6.2 %.
  - The M′ verdict computed again with R including 0 was the same (§27-4).
- **S2 (R21 «exploration after results»):**
  - If H1′ holds (the measured SNR takes the noise as large), then S2 < 1.
  - The two Prismas are 1.19 · 1.75, above 1, and Trio is 0.78. The range between sessions is wide (0.48 – 4.30).
  - Caveat: it mixes the spatial difference in σ between ventricle (centre) and white matter W (g factor · coil sensitivity). W is the processing grid f2_FA placed on the raw image grid as it is. Session SNR was obtained from data after eddy. → It is only written that H1′ is not supported (descriptive).
- **Not done (requires new image reads):** histogram of the low-value side of ventricle b2000 · slices · volumes where 0 values concentrate. This is because voxel values were not saved. The sealed raw images would have to be opened again, so an author decision is needed.
- Output: explore/b2000/zero_check.json (sha256 2e3e9d153cd35a5c…).

## 27-7. b2000 re-access — cause of 0 values · distribution shape [author decision 2026-11-04 · review session reply «b2000 confirmation answer»] (addendum 2026-11-04 00:45:19 KST · definition · before re-access · fixed at the time of this addendum)

- **[author decision] (received as it is):**
  - The scope is the same as what was already opened (109 sessions · b0 · b2000 volumes).
  - The access record carries the tag «post hoc descriptive · verdict unchanged #34».
  - CPU only · 1-hour limit.
  - Only aggregate values are saved. To avoid opening again for the next question, per-scanner low-value histogram counts and shape statistics are kept.
- **Review session reply received:**
  - The rounding hypothesis is withdrawn (review session error · the rounding truth test is not done).
  - In a manuscript note, one sentence is written: «if the 0 values are clipping of low values, R with 0 dropped is biased upward (making L large). So the nearness to Rician is not a result caused by the 0 values».
  - The §27-6 S2 wording is «this comparison cannot separate H1′» in place of «H1′ is not supported». Reason: floor σ came from raw images and session SNR from data after eddy resampling. The resampling reduces session σ and pushes the ratio above 1. The ratio of the two Prismas (1.19 vs 1.75) is kept only as descriptive.
  - Manuscript strength (before the check): c₁ is used as the post hoc main result with the caveat «if variance other than thermal noise differs per scanner, the difference may be underestimated». If the two L agree in (B) below, the caveat is dropped. The b2000 phrase of the abstract is set on 11/7.
  - **The M′ verdict rule and verdict (§27-2 · §27-4) are not changed. Everything below is post hoc descriptive.**
- **Reads:** the three files with the same sha256 as in Step 1 (dwi NIfTI · bval · sidecar).
  - The only volumes used are b0 · b2000 (gzip decompresses the whole file).
  - Passes through the guard. Record with 9 keys · the key «무엇» (what) has «sealed re-access — post hoc descriptive · verdict unchanged #34».
  - Voxel arrays are not saved.
- **Regions:**
  - V = the ventricle mask made again with the same code as in Step 3. It must be the same as V of the saved masks/<session>.npz. If it differs, «mask mismatch» is written and (B) of that session is not produced.
  - H = the head mask (Otsu) of the same code.
- **(A) Cause of 0 values (per session → per scanner):**
  1. Per volume: for each b0 · b2000 volume, the proportion of voxels with value = 0 within V · H, and the number of voxels with value < 0.
  2. Per slice: for each b2000 volume × slice (third NIfTI axis), the 0 proportion within H ∩ slice.
     - «Missing slice» = (volume, slice) with that proportion at 0.5 or above.
     - The share of 0 voxels within V lying in missing slices is given.
  3. Connectedness: for each b2000 volume, the 0 voxels within H are divided into 6-connected 3D components. The sizes of the components to which 0 voxels within V belong are counted as 1 · 2 – 9 · 10 – 99 · 100 or more. 1 – 9 is «scattered», 10 or more is «connected».
  4. Gap just above 0:
     - A histogram of the V b2000 values divided per voxel by σ̂_v is made (cell of exact 0 + [0, 2) in 40 cells of 0.05 width). σ̂_v is the value profiled at L̂ of the (B) maximum likelihood.
     - The expected counts of central chi (L̂) for the same cells are given alongside. «Gap» index = observed / expected of the lowest two cells (0 – 0.1).
     - For integer data types, counts per integer value 0 – 127 of b2000 within V are also given (checking a clipping threshold).
  - Philips · GE are described in the same way.
- **(B) Distribution shape (per session · samples = V voxels × b2000 50 volumes):**
  - Path 1 mean/SD: the same R as in Step 3 (0-exclusion rule) and R_incl0 (0 included). They must be the same as the saved values. L_ratio · L_ratio,incl0 are obtained with the card-level calibration curve.
  - Path 2 maximum likelihood:
    - The model is central chi (degrees of freedom 2L) with a common L and per-voxel σ_v. Only voxels with 40 or more values that are not 0 are used.
    - σ_v² is profiled per L (σ̂_v² = Σ m² / (2 L n_v)). log L is maximized in 1 dimension over [log 0.3, log 30] to give L_ML.
    - Interval: percentile 2.5 – 97.5 of voxel-level sampling with replacement, 200 times (seed fixed).
    - Variant (approximate): the likelihood left-truncated at τ̂, the smallest positive value within V of the session. σ̂_v uses the untruncated value plugged in. → L_ML,trunc.
  - Path 3 skewness:
    - For each voxel, the sample skewness of values that are not 0 (n ≥ 40 · with bias adjustment) is computed, and its median is used.
    - The median skewness curve of pure noise N2 50 (L 0.5 – 128 · 47 points · voxels 10⁵ · seed 20261103 · made before re-access) is inverted to give L_skew.
    - Median kurtosis is given only as descriptive.
  - **Agreement verdict (descriptive · set in advance):** decided per scanner with the median of session values.
    - If |log₂(L_ML / L_ratio)| ≤ ⅓ and |log₂(L_skew / L_ratio)| ≤ ⅓, «agree — read as pure noise».
    - Otherwise «disagree — variance other than thermal noise is mixed in». L_skew > L_ratio is read toward symmetric added variance, L_ML · L_skew < L_ratio as a different shape (descriptive).
  - **Reading for c₁ (Siemens):**
    - If both Prismas are «agree», the c₁ caveat in the manuscript is dropped.
    - If even one is «disagree», δ_u = log₂(L_ML / L_ratio) is compared between the two Prismas. c₁^ML · c₁^skew are given (same contrast as c₁ · common 20 subjects · subject sampling with replacement 2,000 · seed [20261103, 99]), and the manuscript writes c₁ together with c₁ · c₁^ML · c₁^skew, keeping the caveat.
- **Test before reading (voxel reads 0):**
  - See whether Paths 2 · 3 recover L in pure noise (L 0.5 · 1 · 2 · 4 · 8).
  - See in which direction added Gaussian variance pushes Paths 1 · 2 · 3.
  - See whether the (A) indices separate the three in synthetic data with slice loss · scattered clipping · low-value clipping.
- **Outputs:** explore/b2000/shape_sessions.jsonl (per-session aggregates) · shape_summary.json (per scanner) · skew_curve_N2_50.json.

## 29. Paper outline v2 as a reference for writing [author decision 2026-11-04 · review session reply «apply paper outline v2»] (addendum 2026-11-04 01:09:42 KST)

- **[author decision]:** «paper outline v2» is referred to for writing. It is applied together with the previous reply (writing order · devices).
  - The outline was moved to x15/review/NI_outline_v2_2026-11-04.md (sha256 926c9a7f0a06a5ee…). Only the appendix tables were moved as a summary.
- **A's check:** «the previous reply (writing order · devices)» has never come into A's conversation. It was also not found in drafts (searches for «writing order» · «term ledger» 0). → It is requested again from the author. Until it is received, A does not guess and fix a format relying on that reply (numbers · cells of the term ledger etc.).
- **Review session reply (received):**
  1. Sections · paragraphs are written as in the outline. In the working copy each sentence carries an evidence grade mark ([fixed · blinded] · [fixed · unblinded] · [fixed · descriptive] · [post hoc] · [exploratory] · [derived]). In the submitted copy these are turned into words (blinded · unblinded · descriptive · post hoc · exploratory).
  2. All numbers are taken from the number · term ledger. Numbers in the outline are checked again against the evidence field of the ledger. If they disagree, the ledger is taken as right and the review session is told.
  3. The three ⟦⟧ (b2000 shape check sentence · limitations c₁ caveat · fifth Highlight) are filled after the §27-7 shape check.
  4. If a reason arises to write differently from the outline, it is written in one line in the report before changing.
  5. Schedule: 11/4 ledger · Results · Methods · Introduction · Theory → 11/5 Discussion · Conclusion · figures · supplement → 11/6 integration · read-through by the review session → 11/7 final → 11/8 freeze.
- **What changed (where the outline replaces earlier decisions · written in the report):**
  - Length target: main text about 10,500 – 11,000 (Introduction 1,000 · Theory 1,800 · Methods 3,000 – 3,500 · Results 2,800 · Discussion 1,700 · Conclusion 120).
    - It replaces the earlier target 8,820 (§28) and the Methods upper limit 2,500 of ⑥ (author decision «refer to outline v2»).
  - Schedule: the 11/4 – 11/6 allocation of §25-5 · §26 is replaced by 5 above.
- **⟦R⟧ of the outline (for A to check with the ledger · repository):**
  - 3.2 prior grid of the wedge 48.5 – 72.1° (whether r 1.3 – 1.8 only)
  - 3.4 input values · source of −0.148
  - 4.1 why NOT3GEM was not used
  - 4.2 fitter name · ROI definition
  - 5.2 retest values on the same scale
  - F2 content
  - Numbers needing a check: 4.1 «six scanners each» · 4.4 prior · training range · 5.4 T2_iso ≈ 119 ms · 5.6 D1 (z −4.43 · 0.00 · 0.73 · 0.67 · 0.75) · equal floor |z| ≤ 0.60 · FWF ≤ 0.10 proportion 0.67 vs 0.020 · 0.23 – 0.26 at true f_iso 0.034

### 27-7a. Fix before re-access — independent review reflected (addendum 2026-11-04 01:36:42 KST · before re-access · raw image voxel reads 0)

- **Basis:** in-project review workflow (two branches of code · statistics → refutation check).
  - Four code findings were confirmed as real (lowered to major 1 · minor 2 · 1 as it is).
  - The statistics major (direction reading) was not refuted. The refutation check of the statistics blocker (sensitivity of the agreement rule) had not yet finished at the time of this addendum. The rules changed here only make it harder to drop the caveat.
- **Verdict · label fixes ([A's proposal] · confirmation requested from the review session):**
  1. **Condition for dropping the c₁ caveat:** dropped only when both Prismas are «agree», and pair difference |D_skew| ≤ ⅙, and |D_ML| ≤ 0.05.
     - D_k = c₁^k − c₁^same. c₁^same is the c₁ of mean/SD computed with the same voxels · values as Paths 2 · 3. It is computed with same-subject sampling with replacement 2,000 · seed [20261103, 99].
     - If even one exceeds, the caveat is kept, and c₁ · c₁^ML · c₁^skew are written together with no direction.
     - Reason (review synthetic test): with only scanner median ±⅓, even when variance mixed into only one Prisma lowers c₁ by −0.06 ~ −0.10 log₂, 90 – 100 % are «agree».
  2. **Baseline for agreement:** Path 1 computed with the same voxels · values, L_ratio,same. The difference from L_ratio of Step 3 (0-exclusion rule) is written as well.
  3. **Path 2:** L_ML,trunc is used for the verdict. τ̂ = the fifth smallest positive value within V. The truncated likelihood profiles σ_v again under truncation (fixed point 30 times). Untruncated L_ML is also given.
     - The voxel-resampling interval is «sampling variation only» (the bias +1.5 – 2 % from profiling σ per voxel is not captured · it applies equally to both Prismas).
  4. **Disagreement label:** «disagree — differs from central chi shape (the cause is not separated by this check)». The sign pair (δ_ML, δ_skew) is written only as candidates.
     - (−, +): symmetric added variance · partial-volume tissue signal · near-0 values
     - (≈0, +): remaining signal · direction-dependent partial volume · weak symmetric variance
     - (+, −): σ mixing · outlier volumes · volumes with tissue · variance added only in a few volumes · coherent rank-1 · low-value clipping
     - (−, −): L that differs per voxel (thermal noise)
     - The sign pair is not used to read the sign of the c₁ bias. Whether clipping occurred is written only through (A).
     - The agreement label is «agree — within this sensitivity, matches central chi shape». «read as pure noise» is not used. This is because contamination that moves L_ratio by about −7 %, like added Gaussian 0.2σ, can also be «agree».
  5. **Philips · GE:** no agreement label is attached («descriptive · resampling unknown»).
  6. **Note rule (set in advance):**
     - If 80 % or more of the 0 values within V are in components of size 10 or more, «0 values are connected regions (masking · loss) — no clipping bias in R with 0 dropped».
     - Otherwise, if gap < 0.5 and τ̂/σ̂ > 0.02, «0 values agree with low-value clipping — R with 0 dropped is biased upward (making L large)».
     - Anything else is «cause unknown — verdict with 0 included is the same».
  7. **Gap index:** if the expected counts of the lowest two cells are below 20, «not defined».
- **Code fixes:**
  - R_incl0 with 0 included is computed again, and R · R_incl0 are compared with the Step 3 saved values.
  - Failure lines carry the card level. Exact-0 cell · per-voxel 0-count histogram · the whole (volume × slice) 0-proportion matrix are kept.
  - The storage data type and nibabel scale (slope · inter) are recorded. 0 · negative · integer counts are given from stored values.
    - [A's fix · §27-6] The Philips int16 19 sessions were read by nibabel with the scale applied, so «no header scale» in §27-6 is an unconfirmed statement for those sessions. The σ(DN) of those sessions is also in scale-applied units.
  - Per-volume negatives (V · H) are given. An L_ML boundary tag is attached. Skewness · kurtosis are computed with NaN dropped.
  - Whether saved masks exist is checked before reading. The 1-hour limit is kept by the code. The skewness curve sha is cross-checked. The summary cross-checks the §27-2 fixed inputs and calibration table. Whether c1['ratio'] is the same as c₁ of M′ is written.
  - «too few voxels» below 50 voxels (C.MIN_VOX) is a rule not in §27-7. It is written as a deviation.
- **Test before reading (synthetic · voxel reads 0):**
  - In pure noise L 0.5 – 8, Paths 2 · 3 recover L.
  - The indices separate five 0 patterns (none · slice loss · low-value clipping · scattered random 0 · connected masking). Under clipping, L_ML 1.10 becomes 1.02 with the truncated likelihood.
- **sha256:** b2k_shape.py 4945a692841ba4b3fa0b767ab6929bc2239657c24b76fa7d33df8c0c64a44db2 · skew_curve_N2_50.json 2fb3f1668fe96adf9fc5986290c16a857bdaa399a4ecd4986bf3a3d97959de26

### 27-8. Re-access results — cause of 0 values · distribution shape (addendum 2026-11-04 01:38:45 KST · post hoc descriptive · verdict unchanged #34 · re-access 01:36:51 – 01:37:34 · 109/109 · failures 0)

- **Record:**
  - Sealed re-access 327 lines (109 × 3 · 9 keys · «sealed re-access — post hoc descriptive · verdict unchanged #34»)
  - drafts HEAD cdfb37e (§27-7a)
  - Outputs: shape_sessions.jsonl (sha256 abeb73916c75840a…) · shape_summary.json (3f704ec03193a25a…)
  - Voxel arrays were not saved.
- **Consistency:**
  - In all 109 sessions the mask · R · R with 0 included are the same as in Step 3.
  - c₁ computed again is the same as c₁ of M′ (+0.038 [−0.002, +0.076]).
  - Sessions where nibabel applied a scale are NOT1ACH 4 · NOT2ING-98 15. The 0 · integer counts of these sessions were computed again from stored values.
- **(A) Cause of 0 values (descriptive · cause not asserted):**
  - Three Siemens:
    - 0 values of b2000 within V are almost all isolated voxels (size-1 components: OXF1PRI 469/471 · OXF2PRI 336/345 · OXF3TRI 1,500/1,574).
    - Missing slices 0 · no 0 in b0 · negatives 0.
    - The fifth smallest positive is 2 · 4.5 · 1 DN. No low-value clipping threshold is seen.
    - Gap (observed/expected of the lowest two cells) is 1.38 · 1.37 · 0.76.
    - In integer counts, 0 is much more common than 1 (OXF1PRI 471 vs 84 · OXF2PRI 345 vs 22 · OXF3TRI 1,574 vs 3,077). Since the Rayleigh density is proportional to the value near 0, the 0 values of the two Prismas are «exact 0» far more than expected from a continuous distribution.
    - Note rule (§27-7a 6) result: all three «cause unknown — verdict with 0 included is the same».
  - Two GE: isolated 0 (sizes 1 – 9 only) · loss 0 · gap 1.30 · 1.26. 0 is 4 – 12 times more common than 1 (2,730 vs 610 · 4,526 vs 389). Low-value counts are flat toward 0 (half-normal form · agrees with R 1.59 – 1.62 · descriptive). The note is «cause unknown».
  - Philips: NOT2ING-89 has τ̂/σ̂ 0.40 · gap 0.0, so «agrees with low-value clipping». NOT1ACH · NOT2ING-98 have isolated 0 · gap 0.85 · 0.90, so «cause unknown».
- **(B) Distribution shape (Siemens · same-voxel baseline L_ratio,same):**

| Scanner | L_ratio | L_ratio,same | L_ML,trunc | L_skew | δ_ML | δ_skew | Label · sign pair |
|---|---|---|---|---|---|---|---|
| OXF1PRI | 1.033 | 1.029 | 1.042 | 1.083 | +0.018 | +0.074 | agree · (≈0, +) |
| OXF2PRI | 1.060 | 1.054 | 1.072 | 1.129 | +0.025 | +0.099 | agree · (≈0, +) |
| OXF3TRI | 0.931 | 0.930 | 0.881 | 0.810 | −0.078 | −0.199 | agree · (−, −) |

  - c₁ variants (common 20 subjects · subject sampling with replacement):
    - ratio +0.038 [−0.002, +0.076] · same +0.034 [−0.012, +0.066]
    - ML (truncated) +0.040 [−0.010, +0.134] · skew +0.058 [−0.256, +0.239]
  - Pair differences: D_ML = +0.006 [−0.033, +0.083] · D_skew = +0.024 [−0.268, +0.192]
  - **Caveat rule (§27-7a 1) result: «drop»** — both Prismas agree · |D_skew| 0.024 ≤ ⅙ · |D_ML| 0.006 ≤ 0.05.
  - Caveat: the D_skew interval is wide (a little over ±0.2). The precision of the skewness path is low. «agree» means it matches central chi shape within this sensitivity (§27-7a 4).
- **Philips · GE (descriptive · resampling unknown · no label):** L_ratio 0.98 · 1.14 · 1.22 · 0.72 · 0.69, and L_skew 0.99 · 0.76 · 0.84 · 1.02 · 0.90.

### 27-9. Review session reply «b2000 re-access result» (addendum 2026-11-04 01:45:26 KST)

- **Accepted:** three §27-7a fixes · 0-value note as a factual sentence («if clipping, R rises» does not apply to Siemens and is not used) · c₁ caveat dropped · sensitivity sentence in the limitations paragraph · §27-6 fix.
- **Descriptive to add:** that Trio's L_ML,trunc 0.88 · L_skew 0.81 are below 1 is written in the supplementary table note. «agree» stays as it is.
  - [A's fix · fact] «in a chi distribution L cannot be below 1» in the reply is not accurate. There is also chi with 1 degree of freedom (half-normal · L 0.5), and the calibration curve also starts at 0.5. 1 is a lower bound only for magnitude images of complex Gaussian noise (Rician).
  - Note sentence [A's proposal]: «On the Trio, all three estimates were below the Rician value of one (0.93, 0.88, 0.81), which magnitude images of complex noise cannot produce; this is consistent with a partly phase-constrained partial-Fourier reconstruction or with a small departure from a chi distribution, which this check cannot separate.» The Trio L of the main estimate (mean/SD) is also 0.93.
- **To write now:**
  - The b2000 paragraph of Results 5.6
  - The b2000 phrase of abstract sentence 9 (adjusted to 250 words on 11/7)
  - One more Highlight 5 candidate (character count 11/7)
- **Confirmation:** the «writing order · devices» reply (sent again) was received with the 2026-11-04 author decision (ledger format · cross-check devices · order · reporting).

### 27-10. Last refutation check (statistics blocker) · A's check (addendum 2026-11-04 02:16:39 KST · post hoc descriptive · new image reads 0 · saved values only · re-reads recorded)

- **Check verdict:** the finding was not refuted. §27-7a (01:36:42) had already taken in the fix, so the severity is blocker → major.
- **Remaining gap (synthetic test of the check):** the ⅙ rule too can miss c₁ bias of about 0.05 – 0.07 log₂ (proportion dropping the caveat 0.17 – 0.75 at added Gaussian 0.15σ · per-volume σ lognormal 0.1 · σ drift 0.3).
  - → **Manuscript wording rule (set in advance):** written only as far as «the upper end of the c₁ interval (×1.054) is smaller than the generation-difference point estimate (×1.12)». Not written as interval vs interval (c₃ lower end ×1.088).
  - The threshold ⅛ option (optional) is not made a rule since this is after reading. The observed D_skew +0.024 gives the same result under ⅛ as well (descriptive).
- **A's check (computing the Step 3 aggregate claim of the check again with the same saved values · pure-noise expectation is n_used voxels × 50 per session · that session's L · 200 times):**
  - Median observed/expected of the within-session voxel ratio interquartile range: OXF1PRI ×1.23 (above 97.5 % 11/20) · OXF2PRI ×1.62 (16/20) · OXF3TRI ×1.30 (8/10). Close to the values of the check ×1.26 · ×1.65 · ×1.23.
  - Session R outside [2.5, 97.5 %] of bank L 1: OXF1PRI above 14 · below 2 · OXF2PRI above 16 · below 3 · OXF3TRI above 1 · below 9.
  - → The floor is not pure chi noise, and the degree differs per scanner (descriptive · these values were seen before the shape check, so they are not made a new verdict criterion).
- **D descriptive (same-subject sampling with replacement · seed [20261103, 99]):**
  - OXF2PRI − OXF1PRI (20 subjects): D_ML (untruncated) +0.004 [−0.035, +0.081] · D_ML,trunc +0.006 [−0.033, +0.083] · D_skew +0.024 [−0.268, +0.192]
  - OXF1PRI − OXF3TRI (10 subjects): D_ML +0.036 [−0.028, +0.079] · D_ML,trunc +0.119 [+0.054, +0.154] · D_skew +0.236 [+0.156, +0.415]
  - → c₁ of the Prisma pair is small on every path. c₃ on the Trio side (generation difference) grows larger depending on the path. «c₁ < generation difference» holds on every path.
- **Code (small robustness):** if common subjects are below 15, the c1_caveat key is not created. It is changed to write «keep (too few pairs)» (no effect on the present result · n 20).

## 30. Records · manuscript fixes found by the number ledger cross-check (addendum 2026-11-04 02:18:36 KST · ledger workflow · basis = repository · output files)

- **Record fixes (earlier records of this plan · fixed by addendum):**
  - «The NDI cell of the same pair failed in every sensitivity cell» of §26 (20:54:18) → **failed in the three settings where a verdict is possible · SNR 20 · N_eff 100 no verdict possible (null)**. Basis: decomp/results/remeasure/judgements.json V4.cells.OXF1PRI~OXF2PRI|NDI.sensitivity_passed = {SNR20|Neff25 false · SNR20|Neff100 null · SNR50|Neff25 false · SNR50|Neff100 false}.
  - x «+0.013 [+0.006, +0.021]» of §25-3 → **[+0.005, +0.021]** (third decimal place · basis sets.all.x.boot95 = [0.00548, 0.02058]). 0.006 is the all_minus_main_excluded version (0.00580) or double rounding through 0.0055.
  - «[A's check] … T2_in 82 → T2_iso 118.7» of §19-2 (line 869) → **T2_iso ≈ 114.0 ms · T2_ex ≈ 49.7 ms** (from the observed slope medians b_y 3.424 · b_x 7.939 s⁻¹: 1/T2_iso = 1/82 ms − 3.424 s⁻¹). 118.7 is the value corresponding to b_y = b_x · tan 25.4° ≈ 3.77 and does not agree with the slope medians. It is A's check error. Basis: data/remeasure/explore/rv/rv_b_family.json obs.by.
  - «T2_iso ≈ 119 ms» of outline v2 5.4 is a value carrying over the above error. The ledger uses 114 ms (the review session was told).
- **Manuscript working copy fixes (S8's share · tag [A's fix · 2026-11-04 · ledger cross-check]):**
  - Table 3 sensitivity row: «The NDI cell failed in every setting» → «failed in the three settings that could be judged (SNR 20 with N_eff 100 could not be judged)».
  - 4.4.2: «0.006 to 0.021» → «0.005 to 0.021».
  - 3.7.2: «which kept at least nine» → «at least eight» (NOT4GEP~OXF4GEP 8 · others 9). When the Methods section is rewritten, the whole sentence is fixed.
- **Narrative · definition caveats (ledger notes · stated when used in the manuscript):**
  - The wedge 48.5 – 72.1° is a range obtained from the 2-stage grid (r 1.3 · 1.6 · 1.8 · T2_in 70 · 82 · 95 · T2_iso 1000 · T1 6 · decomp/stage1/config.py:249), not from the simulator's relaxation prior. The TE direction of the S4 prior (r 1.0 – 1.8) is 48.5 – 90°.
  - Gong T2 82 · 52: 82 is the value of participant 1, and 52 is the mean of two people (51 · 53). With the mean rule aligned, (82.5, 52) → 57.4°. With the exponential shift reported by Gong it is 46.6°. With T2_iso 2000 ms it is 59.0°.
  - 5.4: y 0.31 times is on the algebraic slope basis (fitter basis 0.34). 25.4° is the median of per-ROI-unit angles (angle of the slope medians 23.3°).
  - Rician reference: the value to compare with the measured R (50 volumes · voxel median) is N2 50 median 1.925 (not the population 1.913).
  - Equal floor settings |z|: maximum 0.6018 → written as a range like «0.13 – 0.60» in place of «≤ 0.60».
  - S2 is in the records under two definitions (session L: 1.19 · 1.75 · 0.78 / scanner L̂: 1.218 · 1.686 · 0.765). Only one is used in the manuscript.
  - The FWF gap (0.080 vs 0.419 · ≤ 0.10 proportion 0.67 vs 0.020) is [exploratory] (analysis ledger · no definition before computing). 0.23 – 0.26 at true f_iso 0.034 is [post hoc] (§19-2 analysis c).
  - 4.4.1 «the coil of five of the seven scanners was unrecorded» is on the judge input basis (the sidecar has OXF3TRI 32 channels). It is written as «five scanners for which no coil was given to the judge».

## §31 Author decision · review session reply — reports ② · ③ · progress report (received 2026-11-04 09:0x KST · recorded by A)

- [author decision] 2026-11-04:
  - Progress: run the review synthesis again → commit fixes · 46 reference items · 11 new papers · figure block as S12 → bundle into the noon report.
  - D1 – D15 as A recommends. D12 follows reply 2-⑥ (one sentence on the old 4.1 fitter output fidelity in Theory 2.5 · new ledger row · the rest in the supplement).
  - Figures: σ_d band as [post hoc], as the source says · symmetric cell in the supplement · ODI (b) sensitivity in a supplementary table · F5 two-column width (190 mm).
  - Incident 2 (email): no personal identifying information in any external request · User-Agent is a generic name · report on checking the code · settings · commit history of the repository and how to delete it · added to the pre-publication checklist.
  - Incident 1 (.pyc): PYTHONDONTWRITEBYTECODE=1 within the repository.
- Review session reply 2 (must be in S12): ① Theory 2.2 Gong Eq. 8 · 9 sentences (together with Introduction · abstract · outline) ② «design flaw found after the verdict» → the before-computing document had «conservative toward passing», and what was learned after the verdict was the size (22 cells) — 4.6 subheading · third contribution in the Introduction · abstract · Highlights ③ «the assumed noise changed the verdict» → the verdict stands · SNR-dependence flag ④ grade of the 24-cell sentence [fixed] + ODI blinded · NDI unblinded within the sentence ⑤ b2000 interpretation limit sentence · check the CMRR basis against the source text (if not possible, the factual sentence only) ⑥ one sentence on fitter output fidelity ⑦ two generalizations in the Introduction (several · within the sources we searched) ⑧ «isotropic-fraction shifts» → «FISO shifts».
- Review session reply 3: after S12, attach to the noon report the list of review synthesis fixes · table of 46 reference fixes · the changed references.bib.
- A's handling: the repository email check result and how to delete it are in drafts/x15/PRE_PUBLICATION_CHECKLIST.md (author · committer 503 commits · prereg/scripts 20 files 21 lines · rewriting history breaks the hash evidence).

### 27-11. Plan-version truth test (rigid 6 only) — deviation record and definition (definition 2026-11-04 10:18:45 KST · before computing · image voxel reads 0 · CPU only)

- **Deviation record:** §27 Step 2 (line 1602) set motion as using «the eddy_parameters (translation 3 · rotation 3) of one actual session of the same card level as rigid transforms». However, b2k_sim.py, which produced the results (header line 3 · lines 141 · 147), also included 10 eddy-current terms (eddy order-2 model · shift along the phase-encoding direction). The verdict (§27-2 M′ · §27-4 Siemens «(A) side · size agrees») was obtained with the bank · calibration table of this code and does not change. Where found: review synthesis (11-04 · decision E13) · handling: author decision 2026-11-04 «E1 – E15 as the review session judges» · review session reply E13.
- **To do (descriptive · not a verdict):** as in the plan, the bank with the eddy-current terms set to 0 is made again for the three Siemens card levels (OXF1PRI · OXF2PRI · OXF3TRI) with the same seeds · same number of simulations (200 × L 9), and with the same summary code the calibration table · effective L of measured session R · c₁ · c₂ · c₃ · M′ reading are computed again.
  - Copies: drafts/x15/b2000/b2k_sim_rigid6.py · b2k_summary_rigid6.py (the only differences from the originals are the eddy-current terms 0 · output folder · fixed hash slots).
  - Outputs: remeasure/explore/b2000/plan_rigid6/ (kat_bank_siemens.jsonl · kat_summary_siemens.json · b2k_result.json).
  - Input re-reads: headers.json · eddy_parameters (the script writes them in REREAD_LOG) · read_sessions.jsonl · pure_curve_N2_50.json (written with reread_log.py when the copies are placed).
- **Report:** how much the effective L · c₁ · c₂ · c₃ · M′ reading change relative to the original version is described. The main text is written on the side of the code that produced the results («eddy-current displacements»), and the plan-version values go in a supplementary table. Even if the reading changes, the verdict is not changed and that fact is written.

#### 27-11 result (addendum 2026-11-04 10:55 KST · descriptive · verdict not changed)

- Run: bank b2k_sim_rigid6.py bank OXF1PRI OXF2PRI OXF3TRI (8 workers · nice 19 · CPU only · 10:19 – 10:54 · 600/600) → named kat_bank_siemens.jsonl (sha256 5e01380b…) → b2k_summary_rigid6.py kat (kat_summary_siemens.json sha256 05d7a38b…) → fixed hash slots filled → result (plan_rigid6/b2k_result.json). The four input copies (read_sessions.jsonl · pure_curve_N2_50.json · kat_bank_PG.jsonl · kat_summary_PG.json) were written with reread_log.py and copied — the two for Philips · GE are because of the fixed hash check of the summary code, and their values are the same as the original version. Comparison: cmp_plan_rigid6.py → plan_rigid6/compare_orig_vs_plan.json · compare_calibration.json (original results · calibration table were read after the re-read record). Image voxel reads 0.
- Calibration table (simulated median of R per L): maximum relative difference from the original version OXF1PRI 0.016 % · OXF2PRI 0.089 % · OXF3TRI 0.23 % (L 0.5). The mask-insufficient counts are the same in all three card levels (45 · 56 · 90). L 4 and L 8 are still separated in all three card levels.
- Effective L (scanner median · 95 % interval): OXF1PRI 1.0324 [1.0218, 1.0510] (original 1.0325 [1.0219, 1.0511]) · OXF2PRI 1.0587 [1.0304, 1.0867] (original 1.0604 [1.0320, 1.0884]) · OXF3TRI 0.9269 [0.8911, 0.9535] (original 0.9309 [0.8950, 0.9575]). Size of change 0.004 or less.
- Contrasts (log₂): c₁ +0.0362 [−0.0040, +0.0739] (original +0.0384 [−0.0017, +0.0760]) · c₂ +0.2124 [+0.1585, +0.2821] (original +0.2085 [+0.1548, +0.2781]) · c₃ +0.1676 [+0.1279, +0.2387] (original +0.1616 [+0.1219, +0.2326]). As ratios c₁ 1.025 [0.997, 1.053] (original 1.027 [0.999, 1.054]).
- M′ reading: gate passed · (A) side · size agrees — same as the original version. The §27 literal rule (supplementary · not interpreted) is also the same «simulator assumption side». Spread ratio (measured IQR ÷ simulated IQR): OXF1PRI 3.47 (original 3.52) · OXF2PRI 6.16 (6.63) · OXF3TRI 6.41 (5.32) — all three scanners remain flagged (> 3).
- Manuscript: the main text stays as the code that produced the results («eddy-current displacements»), with one phrase in the deviation list of Methods 3.6 · plan-version values in the supplementary list (S14).
