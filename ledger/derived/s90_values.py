"""New values of v17c: bias of standard ComBat in NDI units (A1), the isotropic-T2 values (A6) and the linearity margin (A10); --plant must fail."""
import glob, json, math, pathlib, sys

HERE = pathlib.Path(__file__).resolve().parent
OUT, S86, RV = (pathlib.Path(a) for a in sys.argv[1:4])
PLANT = "--plant" in sys.argv
fails = []

TES = [92.0, 92.0, 96.4, 89.0, 98.0, 68.0, 69.65, 70.6]          # the eight acquisitions (same values as s87_values.py)
TE_MEAN = sum(TES) / len(TES)
NU, T2IN, T2EN = 0.46, 82.5, 52.0
expit = lambda v: 1.0 / (1.0 + math.exp(-v))
logit = lambda p: math.log(p / (1.0 - p))
l0 = logit(NU) + TE_MEAN * (1.0 / T2EN - 1.0 / T2IN)
w_in, w_en = math.exp(-TE_MEAN / T2IN), (1.0 if PLANT else math.exp(-TE_MEAN / T2EN))
nu_app_w = NU * w_in / (NU * w_in + (1 - NU) * w_en)
if abs(expit(l0) - nu_app_w) > 1e-12:
    fails.append(f"대조군 겉보기 NDI 두 길 불일치: {expit(l0):.6f} 대 {nu_app_w:.6f}")

cells = {}
for p in [OUT / "cells" / f"cell_{i:03d}.json" for i in (36, 48, 132, 144)]:          # only the cells logged for reading (cell number = lexicographic order of D · S · V · H · N)
    d = json.loads(pathlib.Path(p).read_text(encoding="utf-8"))
    k = d["cell"]
    if float(k["V"]) == 0.0 and k["H"] == "H0" and k["N"] == "bal" and k["S"] in ("S2b", "S2e") and k["D"] in ("D0", "D1"):
        cells[(k["D"], k["S"])] = d["metrics"]["M1"]["bias"]
if len(cells) != 4:
    raise SystemExit(f"칸 {len(cells)}")
raw = dict(te_mean=TE_MEAN, control_nu_app=expit(l0))
for (D, S), b in cells.items():
    raw[f"m1_bias_logit_{D}_{S}"] = b
    raw[f"m1_bias_ndi_{D}_{S}"] = expit(l0 + b) - expit(l0)
    raw[f"m1_bias_ndi_abs_{D}_{S}"] = abs(expit(l0 + b) - expit(l0))          # the manuscript gives only the size, «biased by about …»
s86 = json.loads(S86.read_text(encoding="utf-8"))["raw"]
raw["tol_one_se_per_s"] = s86["m2_se_b"] / TE_MEAN * 1000.0
raw["e_example"] = s86["e_example"]
# A6 linearity margin: from the exact value of ledger N-247 (66/81)
share_x = 0.8148148148148148
n_units = 81
raw["lin_units_x"] = round(share_x * n_units)
raw["lin_need"] = math.ceil(0.8 * n_units)
raw["lin_margin_x"] = raw["lin_units_x"] - raw["lin_need"]
rv = json.loads(RV.read_text(encoding="utf-8"))
m = rv["matches"]
iso_low_r = sorted({x["t2iso"] for x in m if x["r"] < 2.0})
r_other = min(x["r"] for x in m if x["t2iso"] not in iso_low_r) if any(x["t2iso"] not in iso_low_r for x in m) else None
if len(iso_low_r) != 1 or r_other is None or r_other < 2.0 or rv["n_match"] != len(m):
    fails.append(f"격자 조건 꼴이 보충 S-M8 문장과 다름: 비 < 2 의 등방 T2 {iso_low_r} · 나머지 비 최솟값 {r_other}")
raw["grid_iso_t2_low_ratio"] = iso_low_r[0] if iso_low_r else None
raw["grid_min_ratio_other"] = r_other
raw["grid_n_match"], raw["grid_n_conditions"] = rv["n_match"], rv["n_conditions"]
fam = json.loads((HERE / "derived.json").read_text(encoding="utf-8"))["gong"]["mean_both"]["family_T2_iso"]   # same field as ledger N-646
raw["family_t2iso"] = fam
print(json.dumps(raw, indent=0))
print("실패", fails)
if not PLANT:
    (HERE / "s90_values.json").write_text(json.dumps(dict(source=[OUT.name, S86.name, RV.name], raw=raw, fails=fails), ensure_ascii=False, indent=1), encoding="utf-8")
sys.exit(1 if fails else 0)
