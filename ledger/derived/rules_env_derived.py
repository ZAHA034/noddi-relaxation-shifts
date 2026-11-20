"""Derived ledger values from the outputs of the measurement-function environment (rules_env) that are not printed directly (no new computation: pooling, counting, one formula); logs each result json before reading it."""
import hashlib
import json
import pathlib
import subprocess
import sys

sys.dont_write_bytecode = True
HERE = pathlib.Path(__file__).resolve().parent
LEDGER = HERE.parent
OUT = LEDGER.parent / "env" / "out"
SRC = LEDGER.parent / "env" / "src"
sys.path.insert(0, str(SRC))


def logged(name):
    p = OUT / f"{name}.json"
    r = subprocess.run([sys.executable, "-I", "-B", str(LEDGER / "reread_log.py"), str(p), "장부 도출 값(rules_env_derived.py) — rules_env 결과 합치기 · 새 계산 없음"],
                       capture_output=True, text=True)
    h = r.stdout.strip().splitlines()[-1] if r.returncode == 0 and r.stdout.strip() else None
    raw = p.read_bytes()
    if h != hashlib.sha256(raw).hexdigest():
        raise SystemExit(f"STOP: {p} 기록 실패 또는 바뀜")
    return json.loads(raw), h


s1a, h1 = logged("s1a")
s2, h2 = logged("s2")
e1, h3 = logged("e1")
k8, h4 = logged("k8")
e4, h5 = logged("e4")
from rules_env import metrics as MT  # noqa: E402
import numpy as np  # noqa: E402

st = s1a["settings"]
gens = sorted({k.split("|")[0] for k in st})
loc = ("none", "ctrl_offset", "combat_loc", "combat_train_ctrl_loc")
cells = {k: v for k, v in s2["by"].items() if k.startswith("pop|x|two|") and k.split("|")[3] in ("fig3", "l2en", "l2both")
         and k.split("|")[5] in loc}
outside = {k: v for k, v in cells.items() if not (0.80 <= v["ratio"] <= 1.25)}
l2b = [v for k, v in cells.items() if k.split("|")[3] == "l2both"]
s_on = float(s1a["s_x_median"])
tau = float(e1["pooled"]["x"]["tau"])
res = dict(
    sources={"s1a.json": h1, "s2.json": h2, "e1.json": h3, "k8.json": h4, "e4.json": h5},
    s1a=dict(n_generators=len(gens), n_settings_per_pair=len(s1a["eta"]), min_ors=min(v["ors"] for v in st.values()),
             min_spearman=min(v["spearman"] for v in st.values()), max_rank_flip=max(v["rank_flip"] for v in st.values()),
             all_retest_level=all(v["retest_level"] is True for v in st.values())),
    s2=dict(two_site_dTE=100.0 - 70.0, n_cells=len(cells), n_outside=len(outside), outside=sorted(outside),
            all_ratio_at_least_one=all(v["ratio"] >= 1.0 for v in cells.values()),
            ratio_min=min(v["ratio"] for v in cells.values()), ratio_max=max(v["ratio"] for v in cells.values()),
            ratio_range=[min(v["ratio"] for v in cells.values()), max(v["ratio"] for v in cells.values())],
            outside_fitted_range=[min(v["interaction"] for v in outside.values()), max(v["interaction"] for v in outside.values())],
            l2both_eq7=l2b[0]["eq7"], outside_fitted_min=min(v["interaction"] for v in outside.values()),
            outside_fitted_max=max(v["interaction"] for v in outside.values())),
    e1=dict(dte_star_onh=MT.te_tolerance(np.full(27, s_on), tau), s_onh_x_median=s_on, tau_pooled_x=tau),
    e4=dict(wong_62_132_x_minus_raw_ci=[e4["primary"]["wong_62_132_x_minus_raw"]["lo"], e4["primary"]["wong_62_132_x_minus_raw"]["hi"]]),
    k8={g: dict(predicted=k8[f"{g}|n_fitted_predicted"], agree=k8[f"{g}|fitted_agree_observed"]) for g in gens},
)
(HERE / "rules_env_derived.json").write_text(json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")
print(json.dumps(res, ensure_ascii=False))
