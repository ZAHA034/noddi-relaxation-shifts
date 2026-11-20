"""S73 §4 (post hoc): shares at retest level of the Tong negative control by participant and center, the equivalence bound of a common center effect (noncentral t) with its simulated failure control, and the FISO diagnosis of y."""
from __future__ import annotations

import json
import math
import os
import pathlib
import sys

import numpy as np
from scipy import optimize, stats

REPO = pathlib.Path(os.environ.get("S73_REPO") or pathlib.Path(__file__).resolve().parents[3])
sys.path[:0] = [str(REPO / "decomp"), str(REPO / "decomp" / "stage1")]
from rules_env import analyses as RA  # noqa: E402
from rules_env import empirical as EM  # noqa: E402

SCALES = ("x", "y", "odi", "raw_ndi", "logit_fiso")
DF = 6
C975 = float(stats.t.ppf(0.975, DF))
SEED = 20261116


def expected_share(lam: float, df: int = DF, c: float = C975) -> float:
    """Expected share under a common effect δ, P(|T_df(λ)| ≤ c), λ = δ/(√2σ)."""
    return float(stats.nct.cdf(c, df, lam) - stats.nct.cdf(-c, df, lam)) if lam > 0 else float(stats.t.cdf(c, df) - stats.t.cdf(-c, df))


def lam_for(share: float) -> float:
    return float(optimize.brentq(lambda l: expected_share(l) - share, 1e-9, 10.0, xtol=1e-12))


def simulate(delta: float, n_rep: int = 1000, n_roi: int = 27, per_roi: int = 135, seed: int = SEED) -> dict:
    """Definition §4.2 failure control — one s² ~ χ²₆/6 (σ = 1) drawn per ROI unit and shared by its units · differences ~ N(δ, 2)."""
    rng = np.random.default_rng(seed)
    sh = np.empty(n_rep)
    for i in range(n_rep):
        s = np.sqrt(rng.chisquare(DF, n_roi) / DF)
        d = rng.normal(delta, math.sqrt(2.0), (n_roi, per_roi))
        sh[i] = np.mean(np.abs(d) <= C975 * math.sqrt(2.0) * s[:, None])
    return dict(delta=delta, mean=float(sh.mean()), sd=float(sh.std(ddof=1)), analytic=expected_share(delta / math.sqrt(2.0)))


def share(keys, d, thr) -> dict:
    v = [abs(d[k]) <= thr[k] for k in keys]
    return dict(n=len(v), share=(float(np.mean(v)) if v else None))


def main(out_dir: str) -> dict:
    out = pathlib.Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    res = dict(definition="decomp/plans/S73_POST_HOC_PLAN_2026-11-16.md §4", c=C975, df=DF)
    # failure controls (simulation) first
    ctl = {k: simulate(dl) for k, dl in (("delta_0", 0.0), ("delta_1.086", 1.086), ("delta_2", 2.0))}
    ok = (abs(ctl["delta_0"]["mean"] - 0.95) <= 0.005 and abs(ctl["delta_1.086"]["mean"] - 0.90) <= 0.005 and ctl["delta_2"]["mean"] < 0.80)
    res["controls"] = dict(simulation=ctl, ok=bool(ok))
    if not ok:
        res["stopped"] = "모의 실패 대조가 기대대로 움직이지 않음"
        (out / "s73_tong.json").write_text(json.dumps(res, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
        raise SystemExit("⛔ " + res["stopped"])
    tb = RA._tong_scales(EM)
    fixed = (round(tb["x"]["raw"], 3), round(tb["odi"]["raw"], 3), round(tb["y"]["raw"], 3))
    res["fixed_tb_reproduced"] = fixed == (0.931, 0.934, 0.844)
    if not res["fixed_tb_reproduced"]:
        res["stopped"] = f"고정 T-B 비율을 재현하지 못함 {fixed}"
        (out / "s73_tong.json").write_text(json.dumps(res, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
        raise SystemExit("⛔ " + res["stopped"])
    t = EM.read_tong()
    # 4.1 by participant · by center
    res["by_participant"], res["by_center"], res["overall"] = {}, {}, {}
    for sc in SCALES:
        d, thr = tb[sc]["d"], tb[sc]["thr"]
        keys = list(d)
        res["overall"][sc] = share(keys, d, thr)
        for sb in sorted({k[0][0] for k in keys}):
            res["by_participant"].setdefault(sc, {})[sb] = share([k for k in keys if k[0][0] == sb], d, thr)
        for c in sorted({c for k in keys for c in k[0][1:]}):
            res["by_center"].setdefault(sc, {})[c] = share([k for k in keys if c in k[0][1:]], d, thr)
    # 4.2 equivalence bound
    sd_med = {}
    for sc in SCALES:
        values, rep = EM.tong_structs(t, sc)
        units = sorted({u for c in values for sb in values[c] for u in values[c][sb]})
        sd = EM.retest_sd_scalar(rep, units)
        sd_med[sc] = float(np.median([v[0] for v in sd.values() if v[0] is not None]))
    ndi_values, _ = EM.tong_structs(t, "raw_ndi")
    raw = {}
    for (sb, ses), r in t["roi"].items():
        for u, v in EM.scale_values(r, "raw_ndi").items():
            raw.setdefault(u, []).append(v)
    ndi_med = float(np.median([np.mean(v) for v in raw.values()]))
    eq = {"0.90 (threshold)": 0.90, "x": res["overall"]["x"]["share"], "odi": res["overall"]["odi"]["share"], "y": res["overall"]["y"]["share"]}
    res["equivalence"] = dict(no_effect_share=expected_share(0.0), retest_sd_median=sd_med, roi_median_ndi=ndi_med, rows={})
    for k, p in eq.items():
        lam = lam_for(p)
        row = dict(share=p, lam=lam, delta_over_s=math.sqrt(2.0) * lam)
        sc = {"0.90 (threshold)": "x", "x": "x", "odi": "odi", "y": "y"}[k]
        row["delta_logit_or_odi"] = row["delta_over_s"] * sd_med[sc]
        if sc == "x":
            row["delta_ndi_at_0.5"] = row["delta_logit_or_odi"] * 0.25
            row["delta_ndi_at_roi_median"] = row["delta_logit_or_odi"] * ndi_med * (1 - ndi_med)
        res["equivalence"]["rows"][k] = row
    # 4.3 diagnosis of y
    fiso = {}
    for (sb, ses), r in t["roi"].items():
        for u, v in EM.scale_values(r, "raw_fiso").items():
            fiso.setdefault(u, []).append(v)
    fmean = {u: float(np.mean(v)) for u, v in fiso.items()}
    med = float(np.median(list(fmean.values())))
    res["y_diagnosis"] = dict(fiso_median=med, units_below=sorted(u for u, v in fmean.items() if v < med), by_scale={})
    for sc in ("y", "x", "logit_fiso"):
        d, thr = tb[sc]["d"], tb[sc]["thr"]
        us = sorted({k[1] for k in d})
        su = {u: np.mean([abs(d[k]) <= thr[k] for k in d if k[1] == u]) for u in us}
        tau, p = stats.kendalltau([fmean[u] for u in us], [su[u] for u in us], variant="b")
        below = share([k for k in d if fmean[k[1]] < med], d, thr)
        above = share([k for k in d if fmean[k[1]] >= med], d, thr)
        res["y_diagnosis"]["by_scale"][sc] = dict(kendall_tau_b=float(tau), p=float(p), below=below, above=above,
                                                  unit_share={u: float(v) for u, v in su.items()})
    yd = res["y_diagnosis"]["by_scale"]["y"]
    res["y_diagnosis"]["sentence_in_main_text"] = bool(yd["below"]["share"] <= yd["above"]["share"] - 0.05)
    (out / "s73_tong.json").write_text(json.dumps(res, ensure_ascii=False, indent=1, default=float) + "\n", encoding="utf-8")
    return res


if __name__ == "__main__":
    r = main(sys.argv[1])
    print(json.dumps(dict(controls=r["controls"]["ok"], fixed=r["fixed_tb_reproduced"], overall=r["overall"], by_participant=r["by_participant"],
                          equivalence=r["equivalence"], y=dict(tau=r["y_diagnosis"]["by_scale"]["y"]["kendall_tau_b"],
                                                               below=r["y_diagnosis"]["by_scale"]["y"]["below"],
                                                               above=r["y_diagnosis"]["by_scale"]["y"]["above"],
                                                               sentence=r["y_diagnosis"]["sentence_in_main_text"])),
                     ensure_ascii=False, indent=1))
