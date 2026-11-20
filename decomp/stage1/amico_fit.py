"""Fits F2 for stage1 with AMICO 2.1.1 NODDI on a chunk directory, writing NDI, FWF and ODI to f2_<dpar>.npz, in envs/reference."""
import os, shutil, sys
import numpy as np, nibabel as nib
import amico

d, dpar = sys.argv[1], float(sys.argv[2]); nth = int(sys.argv[3]) if len(sys.argv) > 3 else 16
study = os.path.join(d, f"amico_{dpar}"); subj = os.path.join(study, "subj"); os.makedirs(subj, exist_ok=True)
for f in ("DWI.scheme",):
    shutil.copy(os.path.join(d, f), subj)
for f, g in (("data.nii", "DWI.nii"), ("nodif_brain_mask.nii", "mask.nii")):
    dst = os.path.join(subj, g)
    if not os.path.exists(dst):
        os.symlink(os.path.join(d, f), dst)
amico.setup()
ae = amico.Evaluation(study, "subj")
ae.set_config("nthreads", nth); ae.set_config("BLAS_nthreads", 1)
ae.load_data(dwi_filename="DWI.nii", scheme_filename="DWI.scheme", mask_filename="mask.nii", b0_thr=10)
ae.set_model("NODDI"); ae.model.set(dPar=dpar * 1e-3, dIso=3.0e-3)
ae.generate_kernels(regenerate=True); ae.load_kernels(); ae.fit(); ae.save_results()
o = os.path.join(subj, "AMICO", "NODDI")
g = lambda nm: nib.load(os.path.join(o, f"fit_{nm}.nii.gz")).get_fdata().reshape(-1).astype(np.float32)
np.savez(os.path.join(d, f"f2_{dpar}.npz"), ndi=g("NDI"), fwf=g("FWF"), odi=g("ODI"))
shutil.rmtree(study)
print(f"F2 d∥={dpar} 끝", flush=True)
