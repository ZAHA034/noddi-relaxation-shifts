"""Recomputes the values marked ⟦A⟧ in an earlier manuscript version (its Sections 2.2, 2.4 and 2.5 and Figs. 1 and 2; numbering of that version) with repository code, on CPU, writing nothing to the repository."""
import json
import math
import sys

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, f"{sys.argv[1]}/decomp/stage1")
import config as C  # noqa: E402

T2 = (82.0, 52.0, 1000.0)


def xy(nu, f, te, t2=T2):
    a = C.weights(te, np.asarray(t2, float))                 # TE term only (example from the manuscript — TR ignored)
    n2, f2 = C.apparent(nu, f, a)
    return np.asarray(C.frac_coords(n2, f2), float).ravel()[:2]


out = {}
# 2.2 — TE direction angle · 30 ms shift
d = xy(0.6, 0.05, 98.0) - xy(0.6, 0.05, 68.0)
out["dir_deg"] = math.degrees(math.atan2(d[1], d[0]))
out["shift_30ms"] = d.tolist()
# Fig. 1 — do all nine states move with the same vector? (translation · 0.21, 0.34)
moves = [(xy(n, f, 98.0) - xy(n, f, 68.0)).tolist() for n in (0.45, 0.60, 0.75) for f in (0.02, 0.05, 0.12)]
out["fig1_moves_max_abs_dev"] = float(np.max(np.abs(np.asarray(moves) - np.asarray(moves[0]))))
# 2.2 — direction of the neurite decrease (ν = 0.6) at fixed FISO
e = 1e-6
v = np.asarray(C.frac_coords(0.6 - e, 0.05), float).ravel()[:2] - np.asarray(C.frac_coords(0.6, 0.05), float).ravel()[:2]
out["neurite_loss_deg"] = math.degrees(math.atan2(v[1], v[0]))


def nu_app(nu, te, t2in, t2en):
    return float(C.apparent(nu, 0.0, C.weights(te, np.array([t2in, t2en, 1000.0])))[0])


# 2.4 — NDI shift (TE 68 → 92 · T2_in 82 · r 1.0 – 1.8 · ν 0.6)
sh = [nu_app(0.6, 92.0, 82.0, 82.0 / r) - nu_app(0.6, 68.0, 82.0, 82.0 / r) for r in np.linspace(1.0, 1.8, 81)]
out["ndi_shift_te68_92"] = [min(sh), max(sh)]
# 2.5 — extracellular T2 at 2–3 times (ν 0.46 · TE 80 – 110)
dd = [nu_app(0.46, te, 82.0, 52.0 * k) - nu_app(0.46, te, 82.0, 52.0) for te in np.linspace(80, 110, 31) for k in np.linspace(2, 3, 21)]
out["ndi_drop_en_x2_3"] = [min(dd), max(dd)]
# 2.5 · Fig. 2 — difference between patient A (T2_en 52 → 70) and the control · change on the logit scale
dA = {te: nu_app(0.46, te, 82.0, 70.0) - nu_app(0.46, te, 82.0, 52.0) for te in (70.0, 100.0)}
lg = lambda p: math.log(p / (1 - p))
dl = {te: lg(nu_app(0.46, te, 82.0, 70.0)) - lg(nu_app(0.46, te, 82.0, 52.0)) for te in (70.0, 100.0)}
out["patientA_diff_te70_te100"] = [dA[70.0], dA[100.0]]
out["patientA_logit_change_100_vs_70"] = dl[100.0] - dl[70.0]
# Fig. 1 — mean T2 of two participants (E12 · same input as 57°)
T2_MEAN = ((82.0 + 83.0) / 2, (51.0 + 53.0) / 2, 1000.0)
dm = xy(0.6, 0.05, 98.0, T2_MEAN) - xy(0.6, 0.05, 68.0, T2_MEAN)
out["T2_mean"] = list(T2_MEAN)
out["dir_deg_mean"] = math.degrees(math.atan2(dm[1], dm[0]))
out["shift_30ms_mean"] = dm.tolist()
moves_m = [(xy(n, f, 98.0, T2_MEAN) - xy(n, f, 68.0, T2_MEAN)).tolist() for n in (0.45, 0.60, 0.75) for f in (0.02, 0.05, 0.12)]
out["fig1_moves_max_abs_dev_mean"] = float(np.max(np.abs(np.asarray(moves_m) - np.asarray(moves_m[0]))))
# 2.4 · Fig. 2 — same mean T2 (H1 · 2026-11-04 · printed digits as they are — the preceding key is referenced by a ledger row)
t2m = T2_MEAN[0]
ddm = [nu_app(0.46, te, t2m, 52.0 * k) - nu_app(0.46, te, t2m, 52.0) for te in np.linspace(80, 110, 31) for k in np.linspace(2, 3, 21)]
out["ndi_drop_en_x2_3_mean"] = [min(ddm), max(ddm)]
dAm = {te: nu_app(0.46, te, t2m, 70.0) - nu_app(0.46, te, t2m, 52.0) for te in (70.0, 100.0)}
dlm = {te: lg(nu_app(0.46, te, t2m, 70.0)) - lg(nu_app(0.46, te, t2m, 52.0)) for te in (70.0, 100.0)}
out["patientA_diff_te70_te100_mean"] = [dAm[70.0], dAm[100.0]]
out["patientA_logit_change_100_vs_70_mean"] = dlm[100.0] - dlm[70.0]
print(json.dumps(out, indent=1))
