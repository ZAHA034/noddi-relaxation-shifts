# Data

## Derived tables (`data/remeasure/`)

ROI values and other result tables computed in this study. They contain no image. Participant labels are those of the public datasets.

| Folder | Source dataset | License of the source |
|---|---|---|
| `onharmony/` | ON-Harmony, OpenNeuro ds004712, doi:10.18112/openneuro.ds004712.v2.0.1 (Warrington et al., 2025) | CC0 |
| `tong/` | Tong et al. (2020), figshare, doi:10.6084/m9.figshare.8851955.v6 | CC BY 4.0 |
| `wong/` | Wongkornchaovalit et al. (2025), Science Data Bank, doi:10.57760/sciencedb.o00133.00038 | CC BY 4.0 |

Changes made to the source data: the images were preprocessed, fitted and registered to template space, and each table holds ROI means of the fitted maps (Sections 3.1 and 3.2 of the article). Derived ROI values are shared under CC BY 4.0, with attribution to the source datasets: ON-Harmony (CC0), Tong et al. (CC BY 4.0) and Wongkornchaovalit et al. (CC BY 4.0).

`explore/` holds the results of the analyses after the verdict. `MANIFEST_data.csv` lists every file.

## Scheme files (`data/protocols/`)

b values and gradient directions of the two acquisition designs used by the simulator, from ON-Harmony.

## Frozen weights and simulated test data (`data/s4/`, Zenodo deposit)

From the top folder of this repository: `unzip <download folder>/frozen_weights_and_test_data.zip -d data/s4`. `docs/fixed_products_and_commits.md` gives the SHA-256 of the 13 frozen products; `data/MANIFEST_data.csv` and `SHA256SUMS.txt` give all 15 files of `data/s4`. Training data and fitted maps are not deposited because of their size; the code that produced them is in this repository.

## Stored work files of the measurement-function environment (`data/rules_env_work/`, Zenodo deposit)

From the top folder of this repository: `unzip <download folder>/rules_env_work.zip -d data/rules_env_work`. The 270 files are the synthetic signals, scheme files, job lists and AMICO fits of analyses S1b (98 files), S2 (122) and K8 (50) of the environment (Supplementary Methods S-M10), in the folders `s1b/`, `s2/` and `k8/`. Analysis E6 reads the stored S1b fits, and `reproduce/environment_recompute.py --fits` compares new S1b fits with them. The run logs are not included: they contain server paths and no computation reads them.

## Stored stage-3 simulation fits (`data/stage3_fits/`, Zenodo deposit)

From the top folder of this repository: `unzip <download folder>/stage3_fits_E6.zip -d data/stage3_fits`. The 186 files (about 4.1 GB unpacked) are the stored fits of the stage-3 simulation run E6 (`E6/P2a/c*/` and `E6/P2b/c*/`, each with `DONE`, `rows.npz` and `f2_1.7.npz`). Analysis S1a of the environment reads them: `reproduce/environment_recompute.py --stage3 data/stage3_fits`. Run E6 of the simulation is not analysis E6 of the environment (`docs/identifiers.md`).

## Licence of the files on Zenodo

The frozen weights, the simulated test data, the stored work files and the stored fits are results of this study and are shared under CC BY 4.0 (`LICENSE-docs.md`).
