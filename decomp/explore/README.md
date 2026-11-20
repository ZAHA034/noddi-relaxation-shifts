# Analyses after the verdict

These scripts are the code as it was run for the post hoc and exploratory analyses (defined in Section 3.3 of the article; Supplementary Methods S-M8 and S-M9, Supplementary Table S6). The dated plan they follow is `docs/post_verdict_plan.md`, and their stored outputs are in `data/remeasure/explore/`.

They are a record, and most of them do not run as they stand:

- working-folder and repository paths were replaced with `<WORK>` and `<REPO>` when they entered the development repository;
- many read processed images, which are not in this deposit.

Step 1 of the top `README.md` checks the main-text numbers from these analyses against the stored outputs. It does not recompute them.

| Folder | Analysis |
|---|---|
| `cause/` | Prisma pair: raw measured differences, the simulator's expectation, signal-level trends |
| `suff/` | Prisma pair with equal noise-floor settings |
| `d1/` | Bias and overconfidence of the judge; thresholds recalculated with the rule applied to real data |
| `d2/` | Repeated sessions on one scanner against the judge's simulated repeats |
| `d3/` | Known-answer tests and fits for the tissue model checks |
| `d4/` | Forward sensitivity table |
| `d5/` | Features of the judge's failures; `d5/selftest/` holds its self-tests |
| `rv/` | Receive-coil records in the sidecars; the direction of the TE effect; comparison at matched FISO level; recorded coil combinations |
| `b1/`, `b1s/` | Multi-TE series: direction of the TE effect; a simulation reproducing its processing |
| `b2000/` | Magnitude noise floor at b = 2000 s/mm²: known-answer simulation, reading, summary, shape checks |
| `stored_outputs/` | Short tables computed from stored outputs only, and the multi-TE TR check on image headers |
| `ledger_derived/` | Values derived from other recorded values for the number ledger |
| `s73/` | Analyses defined on 2026-11-16 and committed before computing (Supplementary Methods S-M10.6; `decomp/plans/S73_POST_HOC_PLAN_2026-11-16.md`): tissue dependence of the shift at fitter output, the negative control by participant and center with its equivalence bound, the noise floor in the multi-TE data, and the registration check of the multi-TE series, which stopped at its first failure control. `s73_tong.py` needs only the deposited tables and `s73_fitter.py` also the stage-3 fits from Zenodo; `s73_align.py` and the SNR step of `s73_noise.py` read processed images and source data, which are not deposited (`s73_noise_sim.py` runs in `envs/reference`) |

**Script hashes in the outputs.** Each output records the SHA-256 of the script that wrote it, as the script stood when it ran, so the recorded value does not equal the SHA-256 of the file here, whose paths were replaced and comments translated afterwards. Of the 44 recorded values in the outputs of `data/remeasure/explore/` (the outputs of `s73/` record none):

- 6 name an earlier version of a script, kept in a `versions/` folder beside the current one with that hash in its name, and 3 name the development-repository version of the current script;
- 28 name the working-folder copy of the current script before its paths were replaced; apart from string constants (paths) and docstrings, its syntax tree is identical to that of the file here;
- 7 name one of five versions that were not kept: `b2000/kat_summary_PG.json` and `b2000/plan_rigid6/kat_summary_PG.json` (930a250f7f76), `b2000/pure_curve_N2_50.json` and `b2000/plan_rigid6/pure_curve_N2_50.json` (e3387d836e36), `b2000/skew_curve_N2_50.json` (1f56402e0860), `b2000/plan_rigid6/kat_summary_siemens.json` (aed4a6c3d6b8) and `d3/d3_v3_check2_run1.json` (ade8238f8aec).
