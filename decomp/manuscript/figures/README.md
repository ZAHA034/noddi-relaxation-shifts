# Figure scripts

`src/` holds the code that drew the figures, as it was run. Most scripts read the development layout (ledger tools, read logs, working folders), and `make_F1.py` and `make_F4.py` also read processed images or FSL templates, so they do not run from this repository as they stand. They are a record of how each figure was drawn. The values they draw come from the tables in `data/`, and step 1 of the top `README.md` checks the main-text numbers against those tables. `src/fig_theory_v3.py` runs as it stands and writes its figures into this folder; it draws the first versions of the theory figures from the equations.

## Figure in the article → script → output file

Figures 1–4 are in the main text; Figures S1–S6 are in the supplementary material. The output files are the ones named in the figure links of the manuscript and the supplement; the scripts are in `src/`.

| Figure in the article | Script | Output file |
|---|---|---|
| Fig. 1 | `make_F1.py` | `F1_maps` |
| Fig. 2 | `make_F3.py` | `F3_site_TE_group_T2` |
| Fig. 3 | `make_F2.py` | `F2_logratio` |
| Fig. 4 | `make_F8_design_sim.py` | `F8_design_sim` |
| Fig. S1 | `fig_real_data.py` (Wongkornchaovalit et al. panel) | `FigS2_Wong_multiTE` |
| Fig. S2 | `make_F7.py` | `F7_judge` |
| Fig. S3 | `make_FS6.py` | `FS6_env` |
| Fig. S4 | `make_F4.py` | `F4_negctrl_roi_prisma` |
| Fig. S5 | `make_F5.py` | `F5_scanners` |
| Fig. S6 | `make_F6.py` | `F6_te` |

The other scripts in `src/` drew figures of earlier versions of the manuscript or helper panels. `make_F8_design_sim.py` reads the stored results of the design simulation in `decomp/results/design_sim/` and `decomp/results/design_sim_m6/`.

Captions and labels inside the scripts use the numbering of the version in which they were run.
