"""Ledger values of the age × TE analysis: the power estimate before computing for E1 (not run; decomp/age_te/E1_GO_DECISION.md) and the between-person TE slopes of the multi-TE data (E3); --plant must fail."""
import json, math, pathlib, sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
REPO = pathlib.Path(sys.argv[1])
E3 = pathlib.Path(sys.argv[2])
PLANT = "--plant" in sys.argv
sys.path.insert(0, str(REPO / "decomp"))
from age_te import model as M  # noqa: E402

fails = []
te = np.array(list(M.TE_ACQ.values()), float)
var_te = float(np.var(te))

# ── known-answer test: closed form of the balanced design ─────────────────────────────────
rng = np.random.default_rng(20261118)
z = {f"p{i:02d}": float(v) for i, v in enumerate(rng.normal(3.15, 1.05, 20))}
person = [p for p in z for _ in M.TE_ACQ]
acq = [a for _ in z for a in M.TE_ACQ]
if PLANT:
    person, acq = person[1:], acq[1:]
X = M.design(person, acq, z)
inv11 = float(np.linalg.pinv(X.T @ X)[0, 0])
var_z = float(np.var(list(z.values())))
closed = 1.0 / (len(person) * var_te * var_z) if not PLANT else 1.0 / (160 * var_te * var_z)
rel = abs(inv11 - closed) / closed
if rel > 1e-12:
    fails.append(f"균형 설계 닫힌 꼴 대조 실패: 상대 차 {rel:.2e}")

# ── E1 estimate before computing (plan Section 7) ───────────────────────────────────────
raw = dict(balanced_closed_form_rel_error=rel, te_sd=math.sqrt(var_te))
for tag, sig in (("s031", 0.031), ("s060", 0.06)):
    se = sig / math.sqrt(110 * var_te * 1.05 ** 2)
    mde = M.mde(se, 81)
    raw[f"se_{tag}"] = se
    raw[f"mde80_{tag}"] = mde
    raw[f"t2en_len_{tag}"] = M.t2en_change(-mde)["ms"]
    raw[f"t2en_short_{tag}"] = -M.t2en_change(mde)["ms"]
mags = [raw[k] for k in raw if k.startswith("t2en_")]
raw["t2en_mag_min"], raw["t2en_mag_max"] = min(mags), max(mags)
raw["power"] = 0.80
raw["go_threshold_ms_per_decade"] = None
lo, hi = 0.0, 5.0
for _ in range(80):                                   # threshold of the rule on whether to run (power 0.5) — bisection
    mid = (lo + hi) / 2
    lo, hi = (mid, hi) if not M.go_rule(mid)["run"] else (lo, mid)
raw["go_threshold_ms_per_decade"] = hi

# ── E3 (the committed summary as it is) ───────────────────────────────────────────
e3 = json.loads(E3.read_text(encoding="utf-8"))["E3_multite"]
for k, d in e3["by_participant"].items():
    raw[f"e3_{k}_slope_x"] = d["x"]["slope"]
    raw[f"e3_{k}_se_x"] = d["x"]["se"]
    raw[f"e3_{k}_se_retest_x"] = d["x"]["se_retest"]
    raw[f"e3_{k}_slope_odi"] = d["odi"]["slope"]
het = e3["heterogeneity"]
for co in ("x", "y", "odi"):
    raw[f"e3_het_{co}_chi2"] = het[co]["chi2"]
    raw[f"e3_het_{co}_p"] = het[co]["p"]
    raw[f"e3_het_{co}_between_sd"] = het[co]["between_sd"]
    if "retest_se" in het[co]:
        raw[f"e3_het_{co}_p_retest"] = het[co]["retest_se"]["p"]
        raw[f"e3_het_{co}_between_sd_retest"] = het[co]["retest_se"]["between_sd"]
raw["e3_assumed_sd"] = e3["assumed_between_sd_x"]
raw["e3_assumed_shift_30"] = e3["assumed_shift_sd_at_30ms"]
raw["e3_bound_indep_30"] = e3["bounding_shift_sd_at_30ms"]["independent"]
raw["e3_bound_common_30"] = e3["bounding_shift_sd_at_30ms"]["common"]
raw["e3_tau_shift_30"] = 30.0 * het["x"]["between_sd"]
raw["e3_tau_shift_30_retest"] = 30.0 * het["x"]["retest_se"]["between_sd"]
raw["e3_ratio_tau_assumed_retest"] = het["x"]["retest_se"]["between_sd"] / e3["assumed_between_sd_x"]
raw["e3_p_bound"] = 0.001                              # upper bound for «p < 0.001» in the manuscript (actual values are e3_het_x_p_retest · e3_het_odi_p)
for kk in ("e3_het_x_p_retest", "e3_het_odi_p"):
    if not raw[kk] < raw["e3_p_bound"]:
        fails.append(f"{kk} = {raw[kk]} 이 상한 0.001 보다 작지 않다 — 원고 «p < 0.001» 이 틀림")
raw["e3_ratio_tau_assumed"] = het["x"]["between_sd"] / e3["assumed_between_sd_x"]
raw["e3_weighted_mean_slope_x"] = het["x"]["weighted_mean"]
print(json.dumps(raw, indent=0))
print("실패", fails)
if not PLANT:
    (HERE / "s89_values.json").write_text(json.dumps(dict(source=E3.name, raw=raw, fails=fails), ensure_ascii=False, indent=1), encoding="utf-8")
sys.exit(1 if fails else 0)
