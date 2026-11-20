# Amendment of the registration check of 2026-11-16 (English translation)

> English translation of `decomp/plans/S73_AMENDMENT_2026-11-16.md` (SHA-256 65aaa80837866eb9a141973a90ac346453cf6109f635990ba9e93f58a76847ad), committed before the registration was run again (development commit 14e4925, 2026-11-16 10:46:22 KST). Checked by script for numbers and identifiers; the Korean file is the record, and where the two differ the Korean file governs. It amends `docs/s73_post_hoc_plan.md` and `docs/s73_implementation.md`. What was computed, and what stopped, is reported in Supplementary Methods S-M10.6: under this setting failure control (i) failed again, so §2 ended; the figure script of item 6 was therefore not run and is not in the repository.

# S73 §2 amendment — the setting of the registration check and the remaining choices (2026-11-16 · author decision · fixed before computing)

Definition: `decomp/plans/S73_POST_HOC_PLAN_2026-11-16.md` (93c45ce) · implementation note: `decomp/plans/S73_IMPLEMENTATION_2026-11-16.md` (59e5de1). This document changes one registration command of §2, and fixes what the two earlier documents left open in §2 and in §5.3 · §5.4, which use its results. Everything else stays as in the two earlier documents. Before this commit no quantity of paths B · C was computed, and the setting below was not tried on real data (the first computation after this commit is the controls below).

## Why it is amended

- With the setting of the implementation note (`flirt -dof 6 -cost normmi -searchcost normmi`, other settings at their defaults), §2.2 failure control (i) failed: registering the reference mean b0 (2.5 mm isotropic) of each of the three participants to itself gave rotations of 0.33–0.52° · translations of 0.7–1.1 mm, and a largest displacement in the brain of 0.25–0.33 voxel (criterion ≤ 0.1). Control (ii) (a copy shifted by exactly 2 voxels along the first axis) was within its criterion at 2.11–2.21 voxels. As defined, paths B · C were not computed (`<WORK>/remeasure/explore/s73/s73_align.json` · `s73_align_control_i_diag.json`).
- The default setting starts from a wide angular search (±90°) and ended the optimization away from the identity even for identical images. The series were already brought into one space by the provider, so the registration is a small adjustment near the identity.
- Author decision (2026-11-16): change the setting and run again.

## What is changed (one registration command)

`flirt -in <mean b0 of the series> -ref <reference mean b0> -dof 6 -cost normmi -searchcost normmi -nosearch -omat <matrix>`

- `-nosearch` is added: no angular search; the optimization starts from the identity and is local only. All other FLIRT settings are at their defaults.
- The same command is used for controls (i) · (ii) and for the nine series.

## What stays the same

- Control (i) (reference against itself: largest displacement ≤ 0.1 voxel) · (ii) (copy shifted by 2 voxels: 2 ± 0.25 voxels and classified «not aligned»), the alignment criterion (largest displacement ≤ half the voxel size of the series · voxel size = the smallest of the three pixdim), the displacement computation (largest |M·p − p| over the voxel centers in the reference brain mask, in FSL scaled-mm coordinates), the cause check of 2.1, the maps · registration · ROIs of paths B · C, the quantities to compute again, the sentences for each outcome, and the order of the analyses.
- The controls are run on the reference mean b0 of each of the three participants, and all three must pass.

## Condition for ending

- **Under this setting, if either control (i) or (ii) does not behave as expected in even one of the three participants, §2 ends. The setting is not changed further.** §5.4 is then reported without choosing a row.

## Fixed at the same time

1. **Harness check (first):** run the same execution frame used for paths B · C (rules_env e1 · e1_fix · e4 · e5 · e6, changing only the working folder) on path A (the actual working folder), and check that the numbers equal the stored results (the five in `decomp/results/rules_env/`) (all common leaves equal, with only new values added · no leaf present only in the stored values). If they differ, stop. Also, the «main» table of the reference-series ROIs produced by the new code must equal the «main» table of the stored `roi_TE62.json` (as in the implementation note).
2. **Computability of path B:**
   - Linear share of the verdict code (judge_wong): computed from the series that passed. A participant has a value only if at least three of the eight TEs and at least one of the repeat pairs (TE62R2–TE62 · TE132R2–TE132) passed. «Path B linear share» is reported only if all three participants have a value.
   - ODI range ratio · ψ · E1 · E4 · E5 · E6: computed only if the eight TEs and the two repeats passed in all three participants (because ranges · slopes · τ depend on the number of TEs). Otherwise that quantity of path B is «not computable».
   - The rule of definition §2.5, «if path B cannot be computed, the same rule with path C», is followed for each quantity.
3. **τ_obs of §5.3:** the smallest of the pooled τ_x of rules_env e1 of each computed path (A and, if computed, B · C) (τ is the same as in e1_fix, which fixed the finite-n factor).
4. **β of §5.4:** if E1 of path B is computed, path B; otherwise the per-participant β of path C and its 95 % bootstrap interval (`_boot_beta` of e1). The conditions for choosing a row are as in the implementation note. Reading categories of §5: τ as in the implementation note and the definition (the largest τ_sim of the three L). The [part/most] of the second row is filled from the category of τ.
5. **Slope ratio against the literature prediction (definition §2.4):** median of the unit TE slopes (x · y) (per ms) ÷ the literature vector. Literature vector = the T2 of Gong 2020 as in Figure 5 (ledger derivation `derived.json` gong.mean_both: T2_in 82.5 · T2_en 52.0 · T2_iso 1000 ms) — b_x = 1/T2_en − 1/T2_in · b_y = 1/T2_in − 1/T2_iso (per ms).
6. **Supplementary figure (definition §2.4 «a figure like Figure 5»):** a new script `decomp/explore/s73/s73_fig_paths.py` draws it from the stored results (`s73_align.json`) only. For each of paths A · B · C: (a) the share linear in TE (overall · per participant) and 0.8, (b) the unit TE slopes (same scale) · the wedge 48.5–72.1° · the direction of the literature vector · the median direction of the units. No new computation · no new reading. If path B cannot be computed, «not computable» is written in its panel.
7. **ΔTE\* comparison:** whether the pooled ΔTE*_x of e1_fix of path B (or C) lies outside the 90 % interval of path A (15.1–20.7 ms). Sign of E4: whether the sign of the main comparison of e4 (62 against 132 ms, share of x − share of raw) differs from path A.
