"""Known-answer gate of the stage1 synthesizer, without the fitter: reproduces Δx and Δy within 1e-10 and checks T-app and T-in signals."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
import synth as S

ok = True
for prot in C.PROTOCOLS:
    pr = S.Proto(prot)
    for ti in (C.WM[0], C.WM[-1]):
        nu, odi, fiso = C.TISSUE[ti]
        for t2i in (C.T2_SETS.index(C.T2_CENTRE), C.T2_SETS.index((70.0, 1.8, 2000.0))):
            t2 = C.t2_triplet(C.t2_of(t2i))
            a1, a2 = C.weights(68.0, t2), C.weights(98.0, t2)
            x1, y1 = C.frac_coords(*C.apparent(nu, fiso, a1)); x2, y2 = C.frac_coords(*C.apparent(nu, fiso, a2))
            e = max(abs((x2 - x1) - np.log((a2[0] / a2[1]) / (a1[0] / a1[1]))), abs((y2 - y1) - np.log((a2[2] / a2[0]) / (a1[2] / a1[0]))))
            ok &= e < 1e-10
            cell = dict(tissue=ti, t2set=t2i, TE=98.0, TR=np.nan, t1set=-1, dpar=1.7, fdir=3, gen="T-app")
            s_app, _, nu_a, f_a = S.noiseless(pr, cell)
            kp = C.odi_to_kappa(odi)
            ref = (1 - f_a) * (nu_a * pr.E_ic(1.7, kp, 3) + (1 - nu_a) * pr.E_ec(1.7, 1.7 * (1 - nu_a), kp, 3)) + f_a * pr.E_iso()
            e2 = np.abs(s_app - ref).max(); ok &= e2 < 1e-10
            s_in = S.noiseless(pr, dict(cell, gen="T-in"))[0]
            e3 = np.abs(s_in - ref).max()
            if C.t2_of(t2i)[1] != 1.0:
                ok &= e3 > 1e-6
            print(f"{prot} 조직 {ti} T2 {C.t2_of(t2i)}: (i) {e:.1e}  (ii) {e2:.1e}  (iii) T-in 차 {e3:.2e}")
# (iv)
pr = S.Proto("P2b")
Dt = np.diag([1.7e-3, 0.3e-3, 0.3e-3])
sig = np.exp(-pr.bvals * np.einsum("vi,ij,vj->v", pr.bvecs, Dt, pr.bvecs))[None, :]
md, fa = S.ext_coords(sig, pr.bvals, pr.bvecs, 1000.0)[0, :2]
ev = np.array([1.7, 0.3, 0.3]); fa_t = np.sqrt(1.5 * ((ev - ev.mean()) ** 2).sum() / (ev ** 2).sum())
e4 = max(abs(md - ev.mean()), abs(fa - fa_t)); ok &= e4 < 1e-5
print(f"(iv) MD {md:.6f} (참 {ev.mean():.6f}) · FA {fa:.6f} (참 {fa_t:.6f})")
print("관문", "통과" if ok else "⛔ 실패"); sys.exit(0 if ok else 1)
