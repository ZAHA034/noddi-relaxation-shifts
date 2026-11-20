"""Derived ledger values from the S73 results and rules_env s1a that are not printed directly (no new computation: maxima, minima, values before rounding); logs each result json before reading it."""
import hashlib
import json
import pathlib
import subprocess
import sys

sys.dont_write_bytecode = True
HERE = pathlib.Path(__file__).resolve().parent
LEDGER = HERE.parent
S73 = LEDGER.parents[3] / "remeasure" / "explore" / "s73"            # <WORK>/remeasure/explore/s73
ENV = LEDGER.parent / "env" / "out"


def logged(p: pathlib.Path):
    r = subprocess.run([sys.executable, "-I", "-B", str(LEDGER / "reread_log.py"), str(p), "장부 도출 값(s73_derived.py) — S73 결과 합치기 · 새 계산 없음"],
                       capture_output=True, text=True)
    h = r.stdout.strip().splitlines()[-1] if r.returncode == 0 and r.stdout.strip() else None
    raw = p.read_bytes()
    if h != hashlib.sha256(raw).hexdigest():
        raise SystemExit(f"STOP: {p} 기록 실패 또는 바뀜")
    return json.loads(raw), h


fit, h1 = logged(S73 / "s73_fitter.json")
tong, h2 = logged(S73 / "s73_tong.json")
s1a, h3 = logged(ENV / "s1a.json")
noise, h4 = logged(S73 / "s73_noise.json")
pb = [v["prop_bias"] for v in s1a["settings"].values()]
bp = tong["by_participant"]
rng = lambda sc: [min(v["share"] for v in bp[sc].values()), max(v["share"] for v in bp[sc].values())]  # noqa: E731
eq = tong["equivalence"]["rows"]
res = dict(
    sources={"s73_fitter.json": h1, "s73_tong.json": h2, "s1a.json": h3, "s73_noise.json": h4},
    s1a=dict(prop_bias_min=min(pb), prop_bias_max=max(pb)),
    fitter=dict(R_x_max=fit["reading"]["R_x_max"], R_raw_max=fit["reading"]["R_raw_max"],
                x_within_every_pair_class=fit["reading"]["x_within_retest_in_every_pair_class"],
                y_max={k: v["y"]["max_ratio"] for k, v in fit["summary"].items()},
                raw_signal_model_max=max(v["raw_ndi_signal_model"]["max_ratio"] for v in fit["summary"].values())),
    tong=dict(x_min=rng("x")[0], x_max=rng("x")[1], odi_min=rng("odi")[0], odi_max=rng("odi")[1], y_min=rng("y")[0], y_max=rng("y")[1],
              delta_over_s_threshold=eq["0.90 (threshold)"]["delta_over_s"], delta_over_s_x=eq["x"]["delta_over_s"],
              delta_over_s_odi=eq["odi"]["delta_over_s"], delta_ndi_x_at_half=eq["x"]["delta_ndi_at_0.5"],
              delta_ndi_x_at_roi_median=eq["x"]["delta_ndi_at_roi_median"],
              y_below=tong["y_diagnosis"]["by_scale"]["y"]["below"]["share"], y_above=tong["y_diagnosis"]["by_scale"]["y"]["above"]["share"]),
    noise=dict(tau_sim_max={L: v["pooled"]["tau"] for L, v in noise["by_L"].items()}),
)
(HERE / "s73_derived.json").write_text(json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")
print(json.dumps(res, ensure_ascii=False))
