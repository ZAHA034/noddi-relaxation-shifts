# Third-party notices

The code in this repository is under the MIT License (`LICENSE`) and its documents under CC BY 4.0 (`LICENSE-docs.md`).

## JHU ICBM-DTI-81 white-matter label atlas (through FSL `data_atlases`)

The atlas is not included. Step 3 needs the version from FSL `data_atlases`, tag `FSL6-0-1`, and `reproduce/get_atlas.py --from <file>` checks a copy by its SHA-256 (`prereg/artifacts/atlas/README.md`). The atlas stays under its owners' terms; FSL is licensed for non-commercial use. Citation: Mori, S., Oishi, K., Jiang, H., Jiang, L., Li, X., Akhter, K., Hua, K., Faria, A.V., Mahmood, A., Woods, R., Toga, A.W., Pike, G.B., Neto, P.R., Evans, A., Zhang, J., Huang, H., Miller, M.I., van Zijl, P., Mazziotta, J. (2008). Stereotaxic white matter atlas based on diffusion tensor imaging in an ICBM template. NeuroImage 40(2), 570-582. doi:10.1016/j.neuroimage.2007.12.035

## Source datasets

Source datasets: see data/README.md

## Python packages

Dependencies are installed from `uv.lock` and keep their own licenses. AMICO (`dmri-amico`, in `envs/reference`) carries a non-commercial research license. neuroCombat 0.2.12 is installed from PyPI by the lock file and is not redistributed here.
