# Compartmental relaxation can mimic neurite loss in single-echo-time NODDI: what group and multi-site comparisons can show

Code, records and derived data for the article of the same title.

- Zenodo: https://doi.org/10.5281/zenodo.24266517 (this repository together with the frozen weights, the simulated test data and the stored simulation fits)
- GitHub: https://github.com/ZAHA034/noddi-relaxation-shifts (this repository without those files)

No participant image is stored here. The source images are available from their own repositories (see [data/README.md](data/README.md)).

**Names.** The article and its supplementary material use descriptive names; the code and records of this repository keep their working identifiers. The simulation-trained neural posterior estimator of the article is "the judge" here (specification, validation and calibration test: Supplementary Methods S-M5), and the exploratory analyses of Supplementary Methods S-M9 are "the measurement-function environment" here. These and the other identifiers of the code and records (V4, T-A to T-C, C3 and C7, L0 to L3, S1a to E6, the design-simulation labels D0–D2, S0–S3e and M1–M6b, and others) are listed in [docs/identifiers.md](docs/identifiers.md) beside the names used in the article, with the working labels of the translated records (relays such as ZZZZm, steps ⑨–⑪, proposals a1–a16, options (ga) and (na)).

## What is here

**Start with `docs/`.** It holds the English versions of the documents and records:

- the fixed analysis document (`docs/fixed_analysis_document.md`), the judge specification (`docs/judge_specification.md`) and the verdict record (`docs/verdict_record.md`);
- the plan of the analyses after the verdict (`docs/post_verdict_plan.md`), the analysis plan of the measurement-function environment (`docs/environment_analysis_plan.md`), and the definition, implementation note and amendment of the analyses defined on 2026-11-16 (`docs/s73_post_hoc_plan.md`, `docs/s73_implementation.md`, `docs/s73_amendment.md`; Supplementary Methods S-M9.6);
- the plans fixed before computing of the analyses added after the verdict: the multi-site design simulation and its addendum with a second TE (`docs/design_sim_analysis_plan.md`, `docs/design_sim_m6_plan.md`; Supplementary Methods S-M11), the age × TE analysis with its decision not to run E1 (`docs/age_te_analysis_plan.md`, `docs/age_te_e1_decision.md`; S-M12), and the table that applies the relaxation sensitivity value to published clinical results (`docs/rsv_table_plan.md`; S-M13);
- the commits and SHA-256 values of the fixed products with the order of the analysis steps (`docs/fixed_products_and_commits.md`) and the computation record of the judge (`docs/judge_computation_record.md`);
- the clock times, file names, provenance records and judge-study details that the Supplementary material points to (`docs/supplement_details.md`, one section per place in the supplement);
- the identifiers (`docs/identifiers.md`), a guide to every Python file (`docs/code_guide.md`) and a glossary of the Korean words that remain in result files and code (`docs/glossary.md`);
- the table of registered predictions and recorded outcomes of the earlier study of the same data (osf.io/5juwv; `docs/earlier_study_outcomes.md`, with its Korean original beside it). This article tests none of those predictions.

| Path | Contents |
|---|---|
| `docs/` | English documents (above) |
| `reproduce/` | The reproduction scripts described below |
| `data/` | Derived ROI values and other result tables, scheme files ([data/README.md](data/README.md)) |
| `ledger/` | Every number of the main text with its source; `ledger/sources/` keeps copies of the cited sources, and working records there are cut to the cited lines or words, in their original language; `ledger/derived/` holds the scripts that computed the derived values from stored outputs, with their outputs, as they were run (comments in English; paths refer to the development layout) |
| `records/` | Literature search log, access records, the Tong run record |
| `decomp/stage1/`, `decomp/fixed/` | Signal generator, algebraic cards and predictions |
| `decomp/s4/` | The judge (the simulation-trained estimator of the article): simulation, training, validation |
| `decomp/remeasure/` | Image processing, ROI extraction, the tests fixed before computing. `MEASUREMENT_FIXED.md` is the fixed analysis document (Korean original; English in `docs/`) |
| `decomp/rules_env/` | The measurement-function environment (Supplementary Methods S-M9). `ANALYSIS_PLAN.md` is its analysis plan (Korean original; English in `docs/`) |
| `decomp/results/` | Stored results. `decomp/results/remeasure/judgements.json` is the verdict record (English in `docs/`); `decomp/results/rules_env/` holds the outputs of the environment |
| `decomp/explore/` | Analyses after the verdict, as run ([README](decomp/explore/README.md) there) |
| `decomp/design_sim/` | The multi-site design simulation (Supplementary Methods S-M11) and its addendum M6; `ANALYSIS_PLAN.md` and `PLAN_M6.md` are its plans (Korean originals; English in `docs/`). Outputs in `decomp/results/design_sim/` and `decomp/results/design_sim_m6/` |
| `decomp/age_te/` | The age × TE analysis (Supplementary Methods S-M12): plan, decision on E1 and code (Korean originals; English in `docs/`). Output in `decomp/results/age_te/`; no participant information is stored |
| `decomp/rsv_table/` | Plan of the table of published clinical results (Supplementary Methods S-M13; Korean original, English in `docs/`). The extracted values with their locations in each article, and the computed values, are in `decomp/results/rsv_table/` |
| `decomp/v4_recompute/` | Recomputes the calibration test for its main condition or chosen pairs; `CRITERIA.md` gives the comparison rules |
| `decomp/plans/` | Korean originals of the dated plan of the analyses after the verdict and of the definition and implementation note of the analyses of 2026-11-16 (English in `docs/`) |
| `decomp/manuscript/` | Figure scripts as run ([figures/README.md](decomp/manuscript/figures/README.md) maps every figure of the article to its script and file), the checks of the numbers of Supplementary Methods S-M9 and Table S12 (`checks/`), and a summary of the judge's cards |
| `dmri_cov/` | The package that `pyproject.toml` installs, from the project of the earlier study. The code of this article does not import it; its `coverage/` modules are kept because `prereg/scripts/lockd/` refers to them |
| `envs/reference/` | Locked environment of the reference fitting tools (AMICO) |
| `prereg/` | Files of the earlier registered study of the same ON-Harmony data (osf.io/5juwv) that the code reads. The folder name is fixed by the locked code. The study of this article was not itself preregistered; its rules were fixed before computing in the fixed analysis document. The folder also holds scripts of the earlier study that the code cites (`prereg/scripts/`, kept as records; those in `lockd/` import modules of that project that are not deposited, so they do not run here), the access logs of this study (`prereg/artifacts/decomp/`, empty here; see Notes) and the SHA-256 list of the Wongkornchaovalit et al. files |
| `examples/` | Short example that applies the reading rules of Table 2 to a table of ROI values |
| `tests/` | Tests of the code, the example and this repository |

### Files that are only on Zenodo

The Zenodo record has five files. Unpack the three data files in the top folder of this repository:

| File | Contents | Unpack with |
|---|---|---|
| `noddi-relaxation-shifts-v1.0.0.zip` | This repository | `unzip`, `cd noddi-relaxation-shifts-v1.0.0`, then `python3 reproduce/make_local_git.py` (see Setup) |
| `frozen_weights_and_test_data.zip` | Frozen weights and simulated test data of the judge | `unzip <download folder>/frozen_weights_and_test_data.zip -d data/s4` |
| `rules_env_work.zip` | Stored work files of the measurement-function environment: synthetic signals, job lists and fits of S1b, S2 and K8 | `unzip <download folder>/rules_env_work.zip -d data/rules_env_work` |
| `stage3_fits_E6.zip` | Stored fits of the stage-3 simulation run E6 (not analysis E6 of the environment), which analysis S1a of the environment reads | `unzip <download folder>/stage3_fits_E6.zip -d data/stage3_fits` |
| `SHA256SUMS.txt` | SHA-256 of the four archives and of every file inside the three data archives, at its unpacked path (the files of the repository archive are listed with their SHA-256 in `MANIFEST.csv`) | In the download folder: `sha256sum -c --ignore-missing SHA256SUMS.txt`. In the top folder of the repository, after unpacking: `sha256sum -c --ignore-missing <download folder>/SHA256SUMS.txt` |

The run logs of the environment and of the stage-3 fits are not included: they contain server paths and are not read by any computation.

## What the checks cover

| Check | Needs | Time | Covers |
|---|---|---|---|
| Step 1 | This repository | Seconds | The ledger: 1034 rows, 991 recomputed from their source files (677 also compared by SHA-256) and 43 without a machine check (29 values from the literature, and 14 values of the article's own text, records or plans, such as the character counts of the Highlights and the values assumed in the plans of the simulations). The ledger holds every number of the main text that carries a ledger identifier, numbers of the supplementary material, and rows of earlier drafts and of checks of the text, whose `used_in` reads «not printed in the article». The article does not print the identifiers: find a number by its printed form (column `printed_as`) and by `used_in`, which names the section, table, figure or box of the article, or the part of the supplementary material, where it is printed |
| Step 2 | This repository | Minutes | The tests of Table 1 of the article from the ROI tables: Tong et al., Wongkornchaovalit et al., the ON-Harmony scanner pairs and the session SNR table |
| Step 2 with `--judge` | Zenodo `data/s4` as well | About 24 h on CPU | The calibration test of the judge (Supplementary Methods S-M5) with its four sensitivity conditions. Not repeated in full from this repository; see below |
| Main condition of the calibration test | Zenodo `data/s4` as well | About 5 h on CPU | `decomp/v4_recompute/`. Recomputed from the deposit; see below |
| Environment | This repository; Zenodo `rules_env_work.zip` for analysis E6 and `stage3_fits_E6.zip` for analysis S1a | Seconds to minutes | The analyses of the measurement-function environment (Supplementary Methods S-M9), compared with the stored outputs; see below |
| Step 3 | Source images, FSL | Days | Image processing. Not repeated from this repository |
| Figures and analyses after the verdict | This repository | None | The code as it was run and its stored outputs, including the design simulation, the age × TE analysis and the table of published clinical results, and the scripts of `ledger/derived/`. Step 1 checks their main-text numbers against the stored outputs; no command recomputes them |
| Tests | This repository | Minutes | The code, the example and this repository |

These checks do not add up to one command that reproduces the whole article. Training the judge needs the training data, which are not deposited because of their size; the code is in `decomp/s4/`, and the deposit holds the frozen weights.

## Requirements

Tested on Linux x86-64 only. The code uses `fcntl`, so it does not run on Windows. Install [uv](https://docs.astral.sh/uv/) 0.11 (0.11.19 or later, below 0.12): both `pyproject.toml` files set `required-version = ">=0.11.19,<0.12"`, and uv stops with any other version. uv installs Python 3.12 and the locked packages, which take about 6 GB of disk because the locked PyTorch build includes CUDA libraries. No GPU is needed: the judge code stops if a GPU is visible, and the scripts below hide it.

## Setup

From GitHub, use `git clone`. From the Zenodo archive, unpack it (`unzip noddi-relaxation-shifts-v1.0.0.zip`, then `cd noddi-relaxation-shifts-v1.0.0`) and run `python3 reproduce/make_local_git.py` once: the loader guard reads the Git history (see "Notes"), and the archive has none. Then install the environment:

```
uv sync --frozen
```

Unpack the Zenodo data files you need as shown in the table above, and check them from the top folder with `sha256sum -c --ignore-missing <download folder>/SHA256SUMS.txt`.

## Reproducing the results

**Step 1. Check the numbers.**

```
uv run --frozen python reproduce/step1_check_numbers.py
```

**Step 2. Recompute the tests of Table 1** of the article from the ROI tables, on CPU, without reading any image.

```
uv run --frozen python reproduce/step2_recompute.py
```

The script runs the unchanged analysis code with its outputs in a temporary folder and compares every number of the new outputs with the deposited outputs. The locked processing code also creates an empty lock file, `dmri_work/remeasure_locks/access_log.lock`, in the parent folder of this repository, so that folder must be writable.

**The calibration test of the judge** (Supplementary Methods S-M5) runs 5 conditions on 28 pairs (140 runs) and took about 24 h in the original run. Add `--judge` to run it all. To recompute only its main condition or chosen pairs, use `decomp/v4_recompute/run_subset.py` and `compare.py`.

The main condition (SNR 30, N_eff 50) was recomputed this way from the deposit, on CPU with 16 threads on the machine of the original run, in about 5 h. All 29 cells in the three variants (87 cells) have the same JSON text as `data/remeasure/onharmony/v4/v4_c3.json` (Level A in `decomp/v4_recompute/CRITERIA.md`). `decomp/v4_recompute/results/` holds the comparison outputs, including the two failure controls. This shows that the stored output comes from this code and these inputs on that machine. It does not promise the same digits on other hardware or with another thread count; Level B states what then counts as close; Level B alone detected neither failure control (the changed seed of control A and the changed ROI value of control B both met it), and only Level A reported them. The main-condition run used `run_subset.py` as it was before the read log that is written before each input is read was added (development commit ba1fbaf added `read_log.py` and the calls to it and removed no line); the computation is the same. The four sensitivity conditions were not recomputed, so the ledger rows that rest on them are checked only against the stored output.

**The measurement-function environment** (Supplementary Methods S-M9). `decomp/rules_env/` applies the signal model of the article to stated acquisitions and relaxation, fits the synthetic signals with AMICO and reads the ROI values with the loaders of the fixed tests; its analysis plan was committed before any of its new quantities was computed (`docs/environment_analysis_plan.md`). Its outputs are in `decomp/results/rules_env/`. To recompute them and compare every number with the stored outputs:

```
uv run --frozen python reproduce/environment_recompute.py                            # E1–E6, S3 and the E1 fix; about 1 minute
uv run --frozen python reproduce/environment_recompute.py --fits                     # also refit S1b, S2 and K8 with AMICO; 5–8 minutes on 16 cores
uv run --frozen python reproduce/environment_recompute.py --stage3 data/stage3_fits  # also S1a
```

The script works in a temporary clone of this repository (run it from a clean checkout; `git` must be available), links the protocol files to the location the locked code expects and copies `data/`; it writes nothing in your working tree. Options: `--work DIR` (default `data/rules_env_work`), `--threads N` and `--keep DIR`; the exit status is 1 if any analysis fails or differs.

- E6 reads the stored S1b fits from the Zenodo file `rules_env_work.zip` (unpack it to `data/rules_env_work`); without them E6 is skipped with a message and the other analyses still run.
- `--fits` needs the locked AMICO environment: run `uv sync --frozen` in `envs/reference` once (this needs network access). New S1b fits are used for E6 only if their SHA-256 equal those of the stored fits.
- `--stage3` needs the stage-3 fits from the Zenodo file `stage3_fits_E6.zip` (unpack it to `data/stage3_fits`).

On the machine of the original run, all twelve analyses reproduced the stored outputs byte for byte, with and without `--fits`. `decomp/manuscript/checks/make_S17.py` prints Supplementary Table S12 (numbered S17 in an earlier draft, hence the file name) from these outputs. `decomp/manuscript/checks/check_env_numbers.py` compared every number printed in Supplementary Methods S-M9 with these outputs before submission; it needs the text of the supplementary material, which is not part of this deposit, and the environment variables `REMEASURE_WORK=data/remeasure` and `DECOMP_WORK=data`.

**Step 3. Process the images.** Install FSL. Put a copy of the atlas in place with `python3 reproduce/get_atlas.py --from <file>`, which checks its SHA-256 (`prereg/artifacts/atlas/README.md`). Download the source datasets and place them where `decomp/remeasure/rconfig.py` expects them, relative to the parent of this folder:

| Setting in `rconfig.py` | Contents |
|---|---|
| `ONH_SEALED` | ON-Harmony images (OpenNeuro ds004712) |
| `ONH_META` | ON-Harmony metadata files of the same dataset |
| `TONG_ZIPS` | Archives of Tong et al. (figshare) |
| `WONG` | Data of Wongkornchaovalit et al. (Science Data Bank) |

`REMEASURE_WORK` and `REMEASURE_FSLDIR` set the output and FSL folders. The NODDI fits use AMICO 2.1.1 from the separate environment `envs/reference` (`uv sync --frozen` in that folder). Then run `decomp/remeasure/run.py`. Its commands, described in Korean in that locked file, are: `run.py onharmony-batch [N]` (the 110 first sessions and 10 repeat sessions, N at a time); `run.py onharmony <sub> <ses>` (one session); `run.py tong-extract <sub>`, then `run.py tong <sub> <ses>`; and `run.py wong <sub>`. Exit codes 77 (guard stop), 75 (too little free disk) and 76 (locked: another process is working on the same folder or session) are stops, not processing failures. The processing used FSL 6.0.7.23 with `eddy_cpu`. We have not repeated this step from this repository, and values may differ at floating-point level.

**Example.** `uv run --frozen python examples/reading_rules_example.py` applies the reading rules of Table 2 to a synthetic table.

**Tests.**

```
uv run --frozen python -m pytest tests
```

The 300 tests take about 3 minutes on CPU and need no file from Zenodo. Without the atlas, one test of the label-to-unit mapping is skipped (299 pass).

## Notes

- **Locked files.** The 32 guarded files listed in `docs/fixed_products_and_commits.md`, the fixed document, the specification and the verdict record are byte-identical to the versions recorded in the article. Their SHA-256 values can be compared directly with `sha256sum`. A few of them name the folders of the development layout in their default paths; the environment variables above and `reproduce/` set the locations without changing the files. Some of them name the first author, and their comments and messages are in Korean.
- **Access logs.** The two log files under `prereg/artifacts/decomp/` start empty here, and the code appends to them when you run it, so Git shows them as changed afterwards. Supplementary Methods S-M7.1 names these two files as the access records. Their released content is `records/access_records.csv`, together with that of two logs that are not deposited: `decomp/onharmony_access_log`, a read of ON-Harmony acquisition metadata on 2026-10-20 for the simulations, and `unseal/judge_read_log` of the earlier study. It has one row per logged file: `time` (UTC; empty for the rows of logs that list files and hashes without a time of their own), `session` (the session or reader the log names), `file`, `sha256`, `record` (the log the row comes from, as a path under `prereg/artifacts/`), `entry` (the number of the log entry) and `marker` (a deviation mark, in English).
- **Other files.** `MANIFEST.csv` and `data/MANIFEST_data.csv` list every file with its SHA-256 and its relation to the development repository. `paths-neutralized`: only path strings differ (and, in two files, the author's name). `redacted for this repository`: as `paths-neutralized`, and some process notes are shortened (the copy of the dated plan). `comments-translated`: Korean comments and docstrings were replaced by English translations, and path strings were neutralized as for `paths-neutralized`; apart from docstrings and path strings, the syntax tree is identical; in `tests/`, test names were also translated and four tests that read the manuscript or the development layout were left out. Path strings in some result files (for example `decomp/results/F1_compare.json`) name locations of the development layout that are not in this deposit, and `pyproject.toml` and `dmri_cov/__init__.py` mention development documents (`dMRI_PLAN.md`, `prereg/ext_tests.json`) that are not deposited.
- **Commit order.** The loader guard (`decomp/remeasure/guard.py`) requires the fixed document to be committed after the analysis code, so this repository was committed in that order, and `reproduce/make_local_git.py` makes the same three commits for the Zenodo archive. Keep these commits when you copy or mirror the repository (for example with `git push`); a single-commit upload does not pass the guard. The dates of these commits are not evidence of when the rules were fixed. Commit identifiers in the supplementary material and in `docs/` refer to the private development repository, whose history is available to editors and reviewers on request.
- **Tong command line.** The command-line entry of `decomp/remeasure/judge_tong.py` stops with an error, a defect recorded in Supplementary Table S1. `reproduce/tong_judge.py` calls the same functions, as the original run did.
- **Language.** All documents written for readers are in English; the translations in `docs/` are English versions, and the Korean files they translate remain the record. Korean remains, unchanged, in these places: the locked files (above); the Korean originals `decomp/remeasure/MEASUREMENT_FIXED.md`, `decomp/SPEC_S4.md`, `decomp/plans/EXPLORATORY_PLAN_2026-11-01_redacted.md`, `decomp/rules_env/ANALYSIS_PLAN.md`, `decomp/plans/S73_POST_HOC_PLAN_2026-11-16.md`, `decomp/plans/S73_IMPLEMENTATION_2026-11-16.md`, `decomp/plans/S73_AMENDMENT_2026-11-16.md`, `decomp/design_sim/ANALYSIS_PLAN.md`, `decomp/design_sim/PLAN_M6.md`, `decomp/age_te/ANALYSIS_PLAN.md`, `decomp/age_te/E1_GO_DECISION.md`, `decomp/rsv_table/PLAN.md` and `docs/earlier_study_outcomes_ko.md`, each with its English version in `docs/`; keys and values of result and processing files under `decomp/results/`, `data/` and `ledger/derived/`, and run-time messages and labels written by the code, including the messages of the shell scripts, all listed with their English in `docs/glossary.md`; the column `source_field` of `ledger/selectors.csv`, which quotes the cited Korean lines and gives the reasons for rows without a machine check (the English description of each row is in `ledger/numbers.csv`); the source copies in `ledger/sources/`; the run record `records/x12_tong_driver.py`; and the test markers described in `pyproject.toml`, which is locked. Paths named in the Korean records refer to the development layout, and many of those files are not in this repository.
- **Working marks.** Marks of the form ⟦…⟧ inside dated records and code are working marks of the original documents and are kept unchanged. The other working labels of the records are explained in `docs/identifiers.md` ("Labels in the translated records").
- **Frozen weights.** The `.obj.pt` files are Python pickles, and loading a pickle can run code. Compare their SHA-256 with `docs/fixed_products_and_commits.md` before loading them.

## Licenses

Code: MIT ([LICENSE](LICENSE)). Documents, results and derived values: CC BY 4.0 ([LICENSE-docs.md](LICENSE-docs.md)), with attribution to the source datasets, which are under CC0 (ON-Harmony) and CC BY 4.0 (Tong et al.; Wongkornchaovalit et al.); see [data/README.md](data/README.md). Third-party software keeps its own licence, and the FSL atlas is not redistributed ([THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)).

## How to cite

Cite the article and this deposit: Park, J. (2026). Code, records and derived data for: Compartmental relaxation can mimic neurite loss in single-echo-time NODDI: what group and multi-site comparisons can show (Version 1.0.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.24266517. See also [CITATION.cff](CITATION.cff).
