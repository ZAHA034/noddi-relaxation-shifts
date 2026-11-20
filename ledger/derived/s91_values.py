"""Ledger values of the second TE within each site (M6, plan decomp/design_sim/PLAN_M6.md); known-answer test: M2 · M3 of the M6 cells equal the original cells within 1e−12, and a copy paired with shifted cells must fail."""
import glob, json, math, pathlib, sys

HERE = pathlib.Path(__file__).resolve().parent
OUT, OUT6 = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
PLANT = "--plant" in sys.argv
fails = []

m6 = {}
for p in sorted(glob.glob(str(OUT6 / "m6_cells" / "cell_*.json"))):
    d = json.loads(pathlib.Path(p).read_text(encoding="utf-8"))
    m6[d["cell"]["idx"]] = d
if len(m6) != 16:
    raise SystemExit(f"M6 칸 {len(m6)}")
idxs = sorted(m6)
for k, i in enumerate(idxs):
    j = idxs[(k + 1) % len(idxs)] if PLANT else i
    orig = json.loads((OUT / "cells" / f"cell_{j:03d}.json").read_text(encoding="utf-8"))["metrics"]
    for key in ("M2", "M3_te0"):
        for f in ("mean_se", "coverage", "reject"):
            a, b = m6[i]["metrics"][key][f], orig[key][f]
            if abs(a - b) > 1e-12:
                fails.append(f"칸 {i} {key}.{f}: M6 {a} 대 원래 {b}")
                break

by = {(d["cell"]["S"], float(d["cell"]["V"])): d["metrics"] for d in m6.values()}
REL = [(S, V) for S in ("S2b", "S2e", "S2i") for V in (0.0, 0.10)]
raw = dict(dte=m6[idxs[0]]["dte"], noise2=m6[idxs[0]]["noise2"])
for meth in ("M6a", "M6b"):
    covs = [by[c][meth]["coverage"] for c in REL]
    raw[f"{meth}_cov_min"], raw[f"{meth}_cov_max"] = min(covs), max(covs)
    raw[f"{meth}_se_ratio_V0"] = by[("S1b", 0.0)][meth]["se_ratio_to_M2"]
    raw[f"{meth}_se_ratio_V10"] = by[("S1b", 0.10)][meth]["se_ratio_to_M2"]
    raw[f"{meth}_power_S1b_V0"] = by[("S1b", 0.0)][meth]["reject"]
    raw[f"{meth}_power_S1b_V10"] = by[("S1b", 0.10)][meth]["reject"]
raw["M3_se_ratio_V0"] = by[("S1b", 0.0)]["M3_te0"]["se_ratio_to_M2"]
raw["M3_power_S1b_V0"] = by[("S1b", 0.0)]["M3_te0"]["reject"]
raw["M6a_cov_mcse_max"] = max(by[c]["M6a"]["coverage_mcse"] for c in REL)
raw["M6b_n_cells_below_093"] = sum(by[c]["M6b"]["coverage"] < 0.93 for c in REL)
raw["M6a_in_093_097"] = all(0.93 <= by[c]["M6a"]["coverage"] <= 0.97 for c in REL)
print(json.dumps(raw, indent=0))
print("실패", fails[:4], len(fails))
if not PLANT:
    (HERE / "s91_values.json").write_text(json.dumps(dict(source=[OUT.name, OUT6.name], raw=raw, fails=fails), ensure_ascii=False, indent=1), encoding="utf-8")
sys.exit(1 if fails else 0)
