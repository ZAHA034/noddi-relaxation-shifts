"""b10 measures how far fitted logit NDI shifts when TE changes 92 → 68 ms, running three estimator arms on the same synthetic signal."""
import os, sys, numpy as np, nibabel as nib
from amico.scheme import Scheme
from amico.synthesis import NODDIIntraCellular, NODDIExtraCellular, NODDIIsotropic

OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
# §180 control arm. If "app", ties the extracellular tortuosity to the **apparent fraction**.
# Algebraically: after b0 normalization the signal becomes **a NODDI signal with exactly (ν_app, f_iso_app, κ)**.
# → the model misspecification disappears, so the fit must recover ν_app exactly, and the shift must be an identity.
# If the deficit disappears in this arm, it is confirmed that the 9.05% in §179-4 is caused by **tortuosity coupling**.
TORT = os.environ.get("B10_TORT", "true")
assert TORT in ("true", "app"), TORT
# §181-1: separates the **tissue** diffusivity from the **fitter's fixed value**.
# The fitter always forces d∥=1.7 and d_iso=3.0 (§170-2a, cannot be changed).
# If the tissue differs from this, **one more misspecification** arises — Howard d∥ 2.0~2.5 is the literature interval (§160-2).
# B10, B16 and BT-1 were all run under TISSUE_DPAR=1.7 (matching) conditions.
TISSUE_DPAR = float(os.environ.get("B10_TISSUE_DPAR", "1.7"))   # µm²/ms
TISSUE_DISO = float(os.environ.get("B10_TISSUE_DISO", "3.0"))
# §181-3 · review no. 1: ΔR2 is run with both the **median q50** and the **truncated mean** of the recorded prior.
# 0.0086082 is a truncated mean (§61-4). q50 = 7.405e-3 (§23).
DR2_CHOICE = os.environ.get("B10_DR2", "mean")
assert DR2_CHOICE in ("mean", "q50"), DR2_CHOICE
FIT_DPAR, FIT_DISO = 0.0017, 0.003          # Fixed values of the fitter (§170-2a) — for generating the kernel
DPAR, DISO = TISSUE_DPAR*1e-3, TISSUE_DISO*1e-3   # **Tissue** values
T2IN = 82.0
DR2  = {"mean": 0.0086082, "q50": 7.405e-3}[DR2_CHOICE]   # §61-4 truncated mean / §23 q50
T2EN = 1.0/(1.0/T2IN + DR2); T2ISO = 1000.0
TE_REF, TE_TGT = 92.0, 68.0                  # OXF1PRI → NOT4GEP (§136 worst operating point, ΔTE = −24)

k_of = lambda odi: 1.0/np.tan(odi*np.pi/2.0)
lg   = lambda p: np.log(p/(1-p))

rng = np.random.default_rng(0)
def unit(n, r):
    v = r.normal(size=(n,3)); return v/np.linalg.norm(v,axis=1,keepdims=True)
rows=[np.c_[np.zeros((9,3)),np.zeros(9)]]
for b,n in ((1000.,50),(2000.,50)): rows.append(np.c_[unit(n,rng),np.full(n,b)])
table=np.vstack(rows)
sch=Scheme(table,b0_thr=10)
ic,ec,iso=NODDIIntraCellular(sch),NODDIExtraCellular(sch),NODDIIsotropic(sch)
P=ic.protocol_hr
ARGS=(P['grad_dirs'],np.squeeze(P['gradient_strength']),np.squeeze(P['delta']),np.squeeze(P['smalldel']))

# Grid: white-matter tissue inside the claimed interval of §136, [0.40, 0.75]
NDI  = np.array([0.40,0.49,0.58,0.67,0.76,0.85])   # Recorded claimed interval [0.40, 0.850] (§68-3b)
FISO = np.array([0.00,0.05,0.15,0.30])
ODI  = np.array([0.05,0.10,0.20,0.35,0.50])
truth=np.array([(v,f,o) for v in NDI for f in FISO for o in ODI],dtype=float)
drng=np.random.default_rng(11); fdirs=unit(len(truth),drng)   # Random directions (§169-3)

def _nu_app(vic, TE):
    w_in, w_en = np.exp(-TE/T2IN), np.exp(-TE/T2EN)
    return vic*w_in/(vic*w_in + (1-vic)*w_en)

def signal(vic,fiso,odi,fdir,TE):
    """T2-weighted NODDI signal, **normalized by its own b0** (the same convention as the measured DWI). ⚠ The extracellular tortuosity uses the **intrinsic ν** — the actual compartment fraction of the tissue, not the apparent fraction. This mismatch is what DI-1 measures (§148)."""
    kp=k_of(odi); dpar=DPAR*1e-6
    S_ic=ic._synth_meas_watson_SH_cyl_neuman_PGSE(np.array([dpar,0,kp]),*ARGS,fdir)
    nu_t = vic if TORT == "true" else _nu_app(vic, TE)     # ← the only difference of the control arm
    S_ec=ec._synth_meas_watson_hindered_diffusion_PGSE(np.array([dpar,dpar*(1-nu_t),kp]),*ARGS,fdir)
    w_in,w_en,w_iso=np.exp(-TE/T2IN),np.exp(-TE/T2EN),np.exp(-TE/T2ISO)
    S =(1-fiso)*(vic*w_in*S_ic+(1-vic)*w_en*S_ec)+fiso*w_iso*iso.get_signal(DISO)
    S0=(1-fiso)*(vic*w_in+(1-vic)*w_en)+fiso*w_iso
    return np.asarray(S/S0,dtype=np.float32)

def apparent(vic,fiso,TE):
    """Apparent fraction — **the target of estimation** (§142-1). It is not the true value."""
    w_in,w_en,w_iso=np.exp(-TE/T2IN),np.exp(-TE/T2EN),np.exp(-TE/T2ISO)
    A=vic*w_in+(1-vic)*w_en
    return vic*w_in/A, fiso*w_iso/(fiso*w_iso+(1-fiso)*A)

sig=[]; app=[]
for TE in (TE_REF,TE_TGT):
    for (v,f,o),fd in zip(truth,fdirs):
        sig.append(signal(v,f,o,fd,TE)); app.append(apparent(v,f,TE))
sig=np.stack(sig,0); app=np.array(app)
REP=8
vol=np.repeat(sig[:,None,None,:],REP,axis=1).repeat(REP,axis=2)

# cuDIMOT input
nib.save(nib.Nifti1Image(vol,np.eye(4)),f"{OUT}/data.nii.gz")
nib.save(nib.Nifti1Image(np.ones(vol.shape[:3],np.uint8),np.eye(4)),f"{OUT}/nodif_brain_mask.nii.gz")
np.savetxt(f"{OUT}/bvals",table[:,3][None,:],fmt="%.1f"); np.savetxt(f"{OUT}/bvecs",table[:,:3].T,fmt="%.6f")
# AMICO input (same array)
nib.save(nib.Nifti1Image(vol,np.eye(4)),f"{OUT}/DWI.nii.gz")
nib.save(nib.Nifti1Image(np.ones(vol.shape[:3],np.uint8),np.eye(4)),f"{OUT}/mask.nii.gz")
np.savetxt(f"{OUT}/DWI.scheme",table,fmt="%.6f",header="VERSION: BVECTOR",comments="")
np.savez(f"{OUT}/truth.npz",truth=truth,rep=REP,fdirs=fdirs,app=app,
         TE_REF=TE_REF,TE_TGT=TE_TGT,DR2=DR2,T2IN=T2IN,T2EN=T2EN,
         pred_shift=(TE_TGT-TE_REF)*DR2)
print(f"[비틀림={TORT} 조직 d∥={TISSUE_DPAR} d_iso={TISSUE_DISO} ΔR2={DR2_CHOICE}({DR2:.6g})]  조합 {len(truth)} × TE 2 = {len(sig)} 행, 복제 {REP}×{REP} → 복셀 {vol.size//vol.shape[-1]:,}")
print(f"TE {TE_REF} → {TE_TGT} (ΔTE {TE_TGT-TE_REF:+.1f}),  항등식 예측 이동 = {(TE_TGT-TE_REF)*DR2:+.5f} 로짓")
print(f"b0 정규화: [{vol[...,table[:,3]<10].min():.6f}, {vol[...,table[:,3]<10].max():.6f}]  NaN={np.isnan(sig).sum()}")
print(f"겉보기 NDI 범위 {app[:,0].min():.4f}~{app[:,0].max():.4f} (참 {truth[:,0].min():.2f}~{truth[:,0].max():.2f})")
