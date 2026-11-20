# Recomputing the V4 calibration test: comparison criteria

This document was committed before any of the runs below started. The commit that adds it is the record of that order. Dates of commits in the public repository are not evidence of when rules were fixed; the commit in the development repository is.

## Purpose

To check from the public archive alone whether the V4 calibration test recomputes to the values reported in the article. The verdict and the disposition do not change, whatever this comparison shows. A disagreement is written only as a reproducibility note in the Supplementary Methods and the README.

## What is run

- On CPU only (`CUDA_VISIBLE_DEVICES` empty), with `S4_THREADS=16` as in the original run.
- From the unpacked public archive after `reproduce/make_local_git.py`, in the environment installed from `uv.lock`, with the frozen weights and test data of the deposit. No image is read.
- With the unchanged functions of `decomp/remeasure/judge_v4.py`, called by `run_subset.py`. The torch seed of a pair in the main condition is 6000 + its pair index, and the numpy seed is 20261027, as in the original run. `run_subset.py` adds only the choice of pairs and the two controls below.
- Every input file is entered in the read log before it is read.

| Order | Run | Pairs | Why |
|---:|---|---|---|
| 1 | Self-test of `compare.py` | none | The comparison must tell equal, close and different apart |
| 2 | Failure control A: torch seed + 1 | `OXF1PRI~NOT2ING-89` | A changed seed must be reported |
| 3 | Failure control B: one ROI value + 0.5 (first row of the pair, third value of the first session) | `OXF1PRI~NOT2ING-89` | A changed input must be reported |
| 4 | Two pairs | `OXF1PRI~OXF2PRI`, `OXF1PRI~NOT2ING-89` | The first holds the only NDI cell; the second is small. Together 3 cells |
| 5 | Main condition: SNR 30, N_eff 50 | All 28 pairs | 24 cells with data, 5 without |

If neither control is reported as different, the runs stop there and the comparison is not used. If control A is reported and control B is not, this limit of the comparison is stated in the report, and the runs go on.

The four sensitivity conditions are not run. The fields `sensitivity_passed` and `snr_neff_dependent` are therefore not compared, and ledger rows that rest on them are listed as not computed.

## What is compared

**Target.** `data/remeasure/onharmony/v4/v4_c3.json` of the public archive. Before any comparison, its SHA-256 must equal both the value in `data/MANIFEST_data.csv` and the value recorded in `decomp/results/remeasure/judgements.json`. The `inputs` block of a new output must equal that of the target. No file under a `verdict/` folder is opened.

**Two pairs.** Every field of the 3 cells, including coverage, coverage by mixing parameter, threshold, pass or fail, and the number of units. `boot95` is left out as well, because its random stream runs through all earlier cells.

**Main condition.** Every field of all 29 cells in all three variants (`main`, `with_16794_NOT4GEP001`, `without_snr_out_of_train`), including the NDI cell and `boot95`, and the blocks `summary`, `info_none` and `device_concentration_descriptive`. Also every ledger row whose source is a `V4.` field of `decomp/results/remeasure/judgements.json`, evaluated on the new values.

## What counts as the same

**Level A, expected.** The JSON text of every compared field is identical. Each pair has a fixed torch seed, the numpy generators are fixed, and the runs use the same machine, lock file and thread count as the original. Bit-identical results on other machines are not promised.

**Level B, applied only where Level A fails.** For each cell with data, all of these hold:
- the same `verdict` and `passed`;
- the same `threshold`, `n_units` and `n_subjects`;
- |new coverage − reported coverage| ≤ max(0.02, 1 / `n_units_used`). The second term is one unit of the cell.

**Otherwise the cell differs.** Results are reported cell by cell with the run time of each pair.
