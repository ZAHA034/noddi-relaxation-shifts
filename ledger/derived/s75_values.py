"""Five small values (exploratory; existing code, CPU, no new data): the FISO size of Proposition 2, tolerance bounds, noise-only expectations of the multi-TE rules, power of the group × TE term, and existing results moved to the main text."""
import hashlib
import json
import math
import pathlib
import platform
import re
import subprocess
import sys

import numpy as np
from scipy import stats

sys.dont_write_bytecode = True
REPO = pathlib.Path(sys.argv[1]).resolve()
HERE = pathlib.Path(__file__).resolve().parent
LEDGER = HERE.parent
WORK = LEDGER.parents[3]                                            # <data>
sha = lambda p: hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()  # noqa: E731

# ── input: repository code (read only) ──────────────────────────────────────
sys.path.insert(0, str(REPO / "decomp" / "remeasure"))
import judge_wong as JW  # noqa: E402
S4CFG = (REPO / "decomp" / "s4" / "s4cfg.py").read_text(encoding="utf-8")
CFG = (REPO / "decomp" / "stage1" / "config.py").read_text(encoding="utf-8")
AN2 = (REPO / "decomp" / "stage1" / "analyze2.py").read_text(encoding="utf-8")
assert "L2_T2, L2_T1 = 0.0984, 0.0441" in S4CFG and "L3_EN, L3_IN, L3_T1 = (2.1, 3.4), (0.77, 1.02), (1.17, 1.56)" in S4CFG
assert "T1C_RANGE, T1ISO = (850.0, 1100.0), C.T1_ISO" in S4CFG and "for tt in (850.0, 1100.0)" in CFG
for k, v in {"OXF1PRI": "(92.0, 92.0)", "OXF2PRI": "(92.0, 92.0)", "OXF3TRI": "(96.4, 96.4)", "NOT2ING-89": "(89.0, 89.0)",
             "NOT2ING-98": "(98.0, 98.0)", "NOT4GEP": "(68.0, 68.0)", "OXF4GEP": "(68.3, 71.0)", "NOT1ACH": "(70.0, 71.2)"}.items():
    assert re.search(rf'"{k}": \("[SV]", {re.escape(v)}', AN2), k
TE_RANGE = {"OXF1PRI": (92.0, 92.0), "OXF2PRI": (92.0, 92.0), "OXF3TRI": (96.4, 96.4), "NOT2ING-89": (89.0, 89.0),
            "NOT2ING-98": (98.0, 98.0), "NOT4GEP": (68.0, 68.0), "OXF4GEP": (68.3, 71.0), "NOT1ACH": (70.0, 71.2)}

# ── 1. size of Proposition 2 ─────────────────────────────────────────────
F0 = 0.080                                                          # N-509
lam = lambda TR, T1, k: (1 - math.exp(-TR / (k * T1))) / (1 - math.exp(-TR / T1))  # noqa: E731
fnew = lambda f, l: f / (f + (1 - f) * l)                           # noqa: E731
prop2 = {}
for TR in (3.0, 3.6, 5.4):
    row = {}
    for name, ks in (("diffuse", (1.0441,)), ("lesion", (1.17, 1.56))):
        vals = [fnew(F0, lam(TR * 1000, T1, k)) for T1 in (850.0, 1100.0) for k in ks]
        row[name] = dict(fiso_min=min(vals), fiso_max=max(vals), change_min=min(vals) - F0, change_max=max(vals) - F0)
    prop2[f"TR{TR}"] = row
# check: directly with Eq. (3) (config.apparent) — same tissue T1 · TR 3.6 s · T1 1100 · diffuse
sys.path.insert(0, str(REPO / "decomp" / "stage1"))
import config as C  # noqa: E402
a0 = C.weights(0.0, (82.5, 52.0, 1000.0), TR=3600.0, t1=(1100.0, 1100.0, C.T1_ISO))
a1 = C.weights(0.0, (82.5, 52.0, 1000.0), TR=3600.0, t1=(1100.0 * 1.0441, 1100.0 * 1.0441, C.T1_ISO))
_, f_a0 = C.apparent(0.46, 0.05, a0)
_, f_a1 = C.apparent(0.46, 0.05, a1)
chk = fnew(f_a0, lam(3600.0, 1100.0, 1.0441))
assert abs(chk - f_a1) < 1e-12, (chk, f_a1)

# ── 2. tolerance bounds (S73 results as they are · reread log) ────────────────────────
TONG = WORK / "remeasure" / "explore" / "s73" / "s73_tong.json"
r = subprocess.run([sys.executable, "-I", "-B", str(LEDGER / "reread_log.py"), str(TONG), "S75a §6-2 — 0.90 기준의 허용 경계를 장부로(새 계산 없음) · 2026-11-17 다시 실행(등급 문구만 바꿈)"],
                   capture_output=True, text=True)
h = r.stdout.strip().splitlines()[-1] if r.returncode == 0 and r.stdout.strip() else None
if h != sha(TONG):
    raise SystemExit("STOP: s73_tong.json 기록 실패 또는 바뀜")
eq = json.loads(TONG.read_bytes())["equivalence"]
thr = eq["rows"]["0.90 (threshold)"]
tol = dict(delta_over_s=thr["delta_over_s"], delta_logit=thr["delta_logit_or_odi"], ndi_at_half=thr["delta_ndi_at_0.5"],
           ndi_at_roi_median=thr["delta_ndi_at_roi_median"], roi_median_ndi=eq["roi_median_ndi"], retest_sd_median_x=eq["retest_sd_median"]["x"])
assert abs(tol["delta_logit"] - tol["delta_over_s"] * tol["retest_sd_median_x"]) < 1e-12

# ── 3 · 4. noise-only expectations of the multi-TE rules ───────────────────────────────
rng = np.random.default_rng(20261116)
B, NU = 4000, 27
TES = JW.TES
lin_shares, odi_ratio = [], []
for _ in range(B):
    data = {}
    for sb in ("s1", "s2", "s3"):
        ser = {}
        slope, lev = rng.normal(0, 0.003, NU), rng.normal(0.8, 0.3, NU)
        for name, te in [(f"TE{t}", t) for t in TES] + [("TE62R2", 62), ("TE132R2", 132)]:
            v = lev[:, None] + slope[:, None] * (te - 97.0) + rng.normal(0, 1.0, (NU, 3))   # noise SD 1 (the rules do not depend on scale)
            v[:, 2] = 0.3 + rng.normal(0, 1.0, NU)                                           # ODI: does not change with TE + noise
            ser[name] = {f"u{i:02d}": v[i] for i in range(NU)}
        data[sb] = ser
    res = JW.judge(data)
    lin_shares.append(res["overall"]["x"]["linear_ratio"])
    for sb, o in res["by_subject"].items():
        odi_ratio.append(o["ODI"]["te_range_median"] / o["ODI"]["repeat_abs_diff_median"])
t = np.array(TES, float) - 97.0
q = t ** 2 - np.mean(t ** 2)
sd_c_scaled = 35.0 ** 2 / math.sqrt(np.sum(q ** 2))                 # SD of |ĉ|·35² ÷ σ
analytic_known_sigma = 2 * stats.norm.cdf(2 * math.sqrt(2) / sd_c_scaled) - 1
linearity = dict(expected_share_sim=float(np.mean(lin_shares)), sim_q025=float(np.quantile(lin_shares, 0.025)),
                 sim_q975=float(np.quantile(lin_shares, 0.975)), analytic_known_sigma=analytic_known_sigma, replicates=B,
                 observed_x=0.815, observed_y=0.901)
odi = dict(expected_ratio_median=float(np.median(odi_ratio)), q025=float(np.quantile(odi_ratio, 0.025)),
           q975=float(np.quantile(odi_ratio, 0.975)), n=len(odi_ratio))

# ── 5. power of the group × TE term ────────────────────────────────────────
S15 = json.loads((HERE / "residual_s15.json").read_text(encoding="utf-8"))
beta = {k: S15["rows"][k]["dte10"]["logit"] / 10.0 for k in ("diffuse_both", "diffuse_en")}
def power(tes, b, sd, n=20):
    tes = np.asarray(tes, float)
    se = math.sqrt((2 * sd ** 2 / n) / np.sum((tes - tes.mean()) ** 2))
    z = abs(b) / se
    return float(stats.norm.cdf(z - 1.959963984540054) + stats.norm.cdf(-z - 1.959963984540054))
SDS = (0.031, 0.08, 0.125, 0.20)
tes_mid = [sum(v) / 2 for v in TE_RANGE.values()]
pw = {}
for label, tes in (("mid", tes_mid), ("low", [v[0] for v in TE_RANGE.values()]), ("high", [v[1] for v in TE_RANGE.values()])):
    pw[label] = {k: {str(sd): power(tes, b, sd) for sd in SDS} for k, b in beta.items()}
long_te = sorted(k for k, v in TE_RANGE.items() if sum(v) / 2 > 85)

# ── 6. existing results moved to the main text (no new computation · reread log): E2 · E3 · table S14 (c) ─────────────────
def logged(p, why):
    r_ = subprocess.run([sys.executable, "-I", "-B", str(LEDGER / "reread_log.py"), str(p), why], capture_output=True, text=True)
    h_ = r_.stdout.strip().splitlines()[-1] if r_.returncode == 0 and r_.stdout.strip() else None
    if h_ != sha(p):
        raise SystemExit(f"STOP: {p} 기록 실패 또는 바뀜")
    return json.loads(pathlib.Path(p).read_bytes())
ENV = LEDGER.parent / "env" / "out"
E2 = logged(ENV / "e2.json", "S75a §4c — E2 쌍별 몫의 중앙값을 장부로(새 계산 없음) · 2026-11-17 다시 실행(등급 문구만 바꿈)")
E3 = logged(ENV / "e3.json", "S75a §4b — E3 척도별 몫을 장부로(새 계산 없음) · 2026-11-17 다시 실행(등급 문구만 바꿈)")
pairs = E2["pairs"]
e2med = {sc: float(np.median([E2["by"][f"{sc}|{p}"]["ors"] for p in pairs])) for sc in ("raw_ndi", "x", "odi", "y", "logit_fiso", "raw_fiso")}
e2 = dict(n_pairs=len(pairs), median_ors=e2med, tau_b_dTE_orsx=E2["tau_b_dTE_orsx"]["tau"], p_dTE_orsx=E2["tau_b_dTE_orsx"]["p"])
e3 = {sc: dict(raw=v["raw"], ors=v["ors"]) for sc, v in E3["by"].items()}
yd = json.loads(TONG.read_bytes())["y_diagnosis"]["by_scale"]       # the same file logged and read above (hash checked)
s14c = {sc: dict(below=yd[sc]["below"]["share"], above=yd[sc]["above"]["share"], tau=yd[sc]["kendall_tau_b"], p=yd[sc]["p"]) for sc in ("y", "x")}

out = dict(
    what="Values for the v16x revision (2026-11-16; new values exploratory: computed after closely related values had been seen; existing code; no new data): size of Proposition 2, tolerance of the 0.90 share criterion, expected results of the multi-TE rules under noise only, power of the group × TE term, and existing results moved to the main text (E2 and E3, exploratory; Supplementary Table S14 (c), post hoc)",   # 2026-11-16: so that the public copy is in English (values unchanged)
    e2=e2, e3=e3, s14c=s14c,
    sources={"s73_tong.json": sha(TONG), "residual_s15.json": sha(HERE / "residual_s15.json"),
             "e2.json": sha(ENV / "e2.json"), "e3.json": sha(ENV / "e3.json"),
             "judge_wong.py": sha(REPO / "decomp/remeasure/judge_wong.py"), "s4cfg.py": sha(REPO / "decomp/s4/s4cfg.py"),
             "config.py": sha(REPO / "decomp/stage1/config.py"), "analyze2.py": sha(REPO / "decomp/stage1/analyze2.py")},
    python=platform.python_version(), numpy=np.__version__,
    prop2=prop2, tolerance=tol, linearity=linearity, odi_noise_only=odi,
    power=dict(beta_per_ms=beta, sds=list(SDS), n_per_group_site=20, tes_mid=tes_mid, by_te=pw, long_te_acquisitions=long_te,
               assumptions="z test; known variance; no site × group variation beyond what TE explains; groups share T1 and PD; no confounding by vendor"),
)
p36 = prop2["TR3.6"]
out["printed"] = dict(     # one field per value printed in the manuscript (a ledger row points to it)
    p2_L2_T1_pct=round(0.0441 * 100, 10), p2_L3_T1_min=1.17, p2_L3_T1_max=1.56, p2_TR_s=3.6, p2_f0=F0,
    p2_diffuse_min=p36["diffuse"]["fiso_min"], p2_diffuse_max=p36["diffuse"]["fiso_max"],
    p2_lesion_min=p36["lesion"]["fiso_min"], p2_lesion_max=p36["lesion"]["fiso_max"],
    tol_ndi_at_roi_median=tol["ndi_at_roi_median"], roi_median_ndi=tol["roi_median_ndi"], tol_over_s=tol["delta_over_s"],
    e3_raw_fiso_raw=e3["raw_fiso"]["raw"], e3_raw_fiso_ors=e3["raw_fiso"]["ors"],
    e3_logit_fiso_raw=e3["logit_fiso"]["raw"], e3_logit_fiso_ors=e3["logit_fiso"]["ors"],
    s14c_y_below=s14c["y"]["below"], s14c_y_above=s14c["y"]["above"], s14c_x_below=s14c["x"]["below"], s14c_x_above=s14c["x"]["above"],
    s14c_y_tau=s14c["y"]["tau"], s14c_y_p=s14c["y"]["p"],
    e2_x=e2med["x"], e2_raw_ndi=e2med["raw_ndi"], e2_odi=e2med["odi"], e2_y=e2med["y"], e2_logit_fiso=e2med["logit_fiso"],
    e2_raw_fiso=e2med["raw_fiso"], e2_tau=e2["tau_b_dTE_orsx"], e2_p=e2["p_dTE_orsx"], e2_n_pairs=e2["n_pairs"],
    lin_expected=linearity["expected_share_sim"], odi_noise_ratio=odi["expected_ratio_median"],
    pow_both_sd08=pw["mid"]["diffuse_both"]["0.08"], pow_en_sd08=pw["mid"]["diffuse_en"]["0.08"], pow_sd=0.08, pow_n=20,
)
(HERE / "s75_values.json").write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
print(json.dumps({k: out[k] for k in ("prop2", "tolerance", "linearity", "odi_noise_only")}, ensure_ascii=False, indent=1))
print(json.dumps(out["power"]["by_te"]["mid"], indent=1), out["power"]["beta_per_ms"], long_te)
