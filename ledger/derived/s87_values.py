"""Known-answer test of the continuous-covariate corollary of Proposition 4, and the maximum over the unshown cells (all cells with a relaxation difference); --plant must fail."""
import glob, json, pathlib, sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
OUT = pathlib.Path(sys.argv[1])
PLANT = "--plant" in sys.argv
fails = []

# ── known-answer test of the corollary ─────────────────────────────────────────
TES = np.array([92.0, 92.0, 96.4, 89.0, 98.0, 68.0, 69.65, 70.6])   # the eight ON-Harmony acquisitions (same values as power.tes_mid in s75_values.json)
beta, c, b = 0.30, 0.012, 0.0004          # assumed values (logit / z unit, 1/ms, 1/(ms · z unit))
z = np.linspace(0.0, 1.0, 41)             # all sites have the same range of z
rng = np.random.default_rng(20261118)
tau = rng.normal(0, 0.2, len(TES))
delta = np.exp(rng.normal(0, 0.1, len(TES)))
slopes, slopes_mult, adj_slopes = [], [], []
for s, te in enumerate(TES):
    y = 1.0 + beta * z + (0.0 if PLANT else te) * (c + b * z) + tau[s]
    slopes.append(np.polyfit(z, y, 1)[0])
    slopes_mult.append(np.polyfit(z, delta[s] * (y - y.mean()) + y.mean(), 1)[0])
    adj_slopes.append(np.polyfit(z, y - tau[s], 1)[0])          # dropping the additive site term leaves the slope unchanged
slopes, slopes_mult, adj_slopes = map(np.array, (slopes, slopes_mult, adj_slopes))
want = beta + TES * b
err = float(np.max(np.abs(slopes - want)))
err_adj = float(np.max(np.abs(adj_slopes - want)))
err_mult = float(np.max(np.abs(slopes_mult - delta * want)))
pair = float((slopes[4] - slopes[5]) - (TES[4] - TES[5]) * b)       # site at 98 ms versus 68 ms
inv = float((slopes[5] - beta) / TES[5])                            # inversion: relaxation share of the observed slope / TE
if not (err < 1e-10 and err_adj < 1e-10 and err_mult < 1e-10 and abs(pair) < 1e-10 and abs(inv - b) < 1e-10):
    fails.append(f"따름정리 대조 실패: 기울기 {err:.2e} · 가산 {err_adj:.2e} · 승법 {err_mult:.2e} · 사이트 짝 {pair:.2e} · 역산 {inv - b:.2e}")
raw = dict(corollary_max_abs_error=max(err, err_adj, err_mult, abs(pair), abs(inv - b)), corollary_tes=TES.tolist(),
           corollary_inputs=dict(beta=beta, c=c, b=b, z_points=len(z)))

# ── unshown cells: all cells with a relaxation difference (S2b · S2e · S2i · S3b · S3e) ─────────
cells = {}
for p in sorted(glob.glob(str(OUT / "cells" / "cell_*.json"))):
    d = json.loads(pathlib.Path(p).read_text(encoding="utf-8"))
    k = d["cell"]
    cells[(k["D"], k["S"], float(k["V"]), k["H"], k["N"])] = d["metrics"]
if len(cells) != 288:
    raise SystemExit(f"칸 {len(cells)}")
REL = ("S2b", "S2e", "S2i", "S3b", "S3e")
shown = lambda k: k[2] == 0.0 and k[3] == "H0" and k[4] == "bal"
raw["m1_cov_rel_unseen_noH2_max"] = max(cells[k]["M1"]["coverage"] for k in cells if k[1] in REL and not shown(k) and k[3] != "H2")
raw["m1_cov_rel_H2_max"] = max(cells[k]["M1"]["coverage"] for k in cells if k[1] in REL and k[3] == "H2")
raw["m1_cov_rel_shown_max"] = max(cells[k]["M1"]["coverage"] for k in cells if k[1] in REL and shown(k))
print(json.dumps({k: v for k, v in raw.items() if k != "corollary_tes"}, indent=0))
print("실패", fails)
if not PLANT:
    (HERE / "s87_values.json").write_text(json.dumps(dict(source=OUT.name, raw=raw, fails=fails), ensure_ascii=False, indent=1), encoding="utf-8")
sys.exit(1 if fails else 0)
