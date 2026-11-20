"""Further values (no new simulation: stored cell metrics, signal-model formulas and arithmetic on ledger values), each computed a second way as an explicit check."""
import csv, glob, hashlib, json, math, pathlib, sys

from scipy import optimize, stats

HERE = pathlib.Path(__file__).resolve().parent
OUT = pathlib.Path(sys.argv[1])
cells = {}
for p in sorted(glob.glob(str(OUT / "cells" / "cell_*.json"))):
    d = json.loads(pathlib.Path(p).read_text(encoding="utf-8"))
    c = d["cell"]
    cells[(c["D"], c["S"], float(c["V"]), c["H"], c["N"])] = d["metrics"]
fails = []
if len(cells) != 288:
    raise SystemExit(f"칸 {len(cells)}")


def m(D, S, k, V=0.0, H="H0", N="bal"):
    return cells[(D, S, V, H, N)][k]


raw = {}
# (6) SE multiple of the within-site bias — the stored mean_se (checked against the formula √((0.08²+0.031²)·2/20)/√8 = 0.00959)
for S, tag in (("S2b", "b"), ("S2e", "e")):
    raw[f"m1_bias_{tag}"] = m("D1", S, "M1")["bias"]
    raw[f"m2_se_{tag}"] = m("D1", S, "M2")["mean_se"]
    raw[f"bias_over_se_{tag}"] = abs(m("D1", S, "M2")["bias"]) / m("D1", S, "M2")["mean_se"]
formula = math.sqrt((0.08 ** 2 + 0.031 ** 2) * 2 / 20) / math.sqrt(8)
raw["se_formula"] = formula
if abs(raw["m2_se_b"] - formula) > 3e-4:
    fails.append("SE 공식과 저장값이 다름")
# (8) SE of the extrapolation
raw["m3_se_b"] = m("D1", "S1b", "M3_te0")["mean_se"]          # the same cell as the SE printed by S-M11 (R-c · D1 · S1b 0.068)
raw["m3_se_e"] = m("D1", "S2e", "M3_te0")["mean_se"]
raw["m3_over_pooled"] = raw["m3_se_b"] / raw["m2_se_b"]
raw["m3_ci_ndi"] = 1.959963984540054 * raw["m3_se_b"] * 0.46 * 0.54
# (9) unshown cells
keyset = [k for k in cells if k[1] in ("S2b", "S2e") and not (k[2] == 0.0 and k[3] == "H0" and k[4] == "bal")]
raw["m1_cov_unseen_max_noH2"] = max(cells[k]["M1"]["coverage"] for k in keyset if k[3] != "H2")
raw["m1_cov_unseen_max_H2"] = max(cells[k]["M1"]["coverage"] for k in keyset if k[3] == "H2")
raw["m3_pow_V01_S1e"] = m("D1", "S1e", "M3_te0", V=0.1)["reject"]
raw["m3_pow_V01_S1b"] = m("D1", "S1b", "M3_te0", V=0.1)["reject"]
# signal-model arithmetic (the controls of 2.4)
NU, TIN, TEN, P = 0.46, 82.5, 52.0, 0.0984


def app(tin, ten, te):
    a, b = NU * math.exp(-te / tin), (1 - NU) * math.exp(-te / ten)
    return a / (a + b)


def bulk(tin, ten, tes):
    n = len(tes)
    ys = [math.log(NU * math.exp(-t / tin) + (1 - NU) * math.exp(-t / ten)) for t in tes]
    mx, my = sum(tes) / n, sum(ys) / n
    return -1 / (sum((x - mx) * (y - my) for x, y in zip(tes, ys)) / sum((x - mx) ** 2 for x in tes))


tes = list(range(10, 101, 10))
b0 = bulk(TIN, TEN, tes)
mm = optimize.brentq(lambda k: bulk(TIN, TEN * k, tes) / b0 - (1 + P), 1.0, 3.0)
raw["bulk_match_m"] = mm
raw["bulk_match_drop_70"] = app(TIN, TEN, 70) - app(TIN, TEN * mm, 70)
raw["bulk_match_drop_100"] = app(TIN, TEN, 100) - app(TIN, TEN * mm, 100)
for lo, hi, key in ((10, 120, "10_120"), (70, 100, "70_100")):
    t2 = list(range(lo, hi + 1, 10 if lo == 10 else 5))
    bb = bulk(TIN, TEN, t2)
    k2 = optimize.brentq(lambda k: bulk(TIN, TEN * k, t2) / bb - (1 + P), 1.0, 3.0)
    raw[f"bulk_match_range_{key}"] = [app(TIN, TEN, 70) - app(TIN, TEN * k2, 70), app(TIN, TEN, 100) - app(TIN, TEN * k2, 100)]
raw["intra_raise_70"] = app(TIN * (1 + P), TEN, 70) - app(TIN, TEN, 70)
raw["intra_raise_100"] = app(TIN * (1 + P), TEN, 100) - app(TIN, TEN, 100)
# inversion (TE 90 ms · deficit 0.02)
p90 = app(TIN, TEN, 90)
lg = lambda x: math.log(x / (1 - x))
raw["inv_te"], raw["inv_deficit"] = 90.0, 0.02
raw["inv_d_logit"] = lg(p90) - lg(p90 - 0.02)
t2en = optimize.brentq(lambda t: app(TIN, t, 90) - (p90 - 0.02), 52.0, 80.0)
raw["inv_t2en"] = t2en
raw["inv_pct"] = (t2en / TEN - 1) * 100
d_r2 = raw["inv_d_logit"] / 90.0                     # Eq. (4): logit difference = TE · Δ(R2,en − R2,in)
if abs((1 / (1 / TEN - d_r2)) - t2en) > 0.01:
    fails.append("역산: 식 (4) 길과 정확 해가 다름")
# (5) strict invariance · (7) TE·e
raw["strict_ratio"] = (98.0 / 68.0) ** 2
s75 = json.loads((HERE / "s75_values.json").read_text(encoding="utf-8"))
raw["te_c_ms"] = sum(s75["power"]["tes_mid"]) / 8
raw["e_example"] = 0.5
raw["te_e_logit"] = raw["te_c_ms"] / 1000 * 0.5
raw["fiso_lesion_min"], raw["fiso_lesion_max"] = s75["printed"]["p2_lesion_min"], s75["printed"]["p2_lesion_max"]
# retest SD changed (common-effect model · ledger value) — only in S-M9
led = {r["ID"]: r for r in csv.DictReader(open(HERE.parent / "numbers.csv", encoding="utf-8"))}
g = lambda i: float(led[i]["원 값"])
share, sc, so = {"x": g("N-324"), "y": g("N-326")}, {"x": g("N-991"), "y": g("N-992")}, {"x": g("N-993"), "y": g("N-994")}
resc = {}
for name, tq in (("normal", 1.959963984540054), ("t6", stats.t.ppf(0.975, 6)), ("t8", stats.t.ppf(0.975, 8))):
    for c in ("x", "y"):
        k = optimize.brentq(lambda k: stats.norm.cdf(tq - k) - stats.norm.cdf(-tq - k) - share[c], 0, 6)
        k2 = k * sc[c] / so[c]
        resc[f"{name}_{c}"] = stats.norm.cdf(tq - k2) - stats.norm.cdf(-tq - k2)
raw["rescaled_shares"] = resc
raw["rescaled_x_range"] = [min(v for k, v in resc.items() if k.endswith("_x")), max(v for k, v in resc.items() if k.endswith("_x"))]
raw["rescaled_y_range"] = [min(v for k, v in resc.items() if k.endswith("_y")), max(v for k, v in resc.items() if k.endswith("_y"))]
if not raw["rescaled_x_range"][1] < 0.90:
    fails.append("x 가 어느 문턱에서 기준을 넘음 — «would not hold» 문장이 틀림")
# table S13 (e): unshown cells (D1) — between-person relaxation variation · unbalanced · random site × group term
rows_e = []
for lab, V, H, N in (("relaxation variation (CV 0.10)", 0.1, "H0", "bal"), ("unbalanced sites", 0.0, "H0", "unbal"), ("random site × group term", 0.0, "H1", "bal")):
    for S in ("S1b", "S1e", "S2b", "S2e"):
        mm_ = cells[("D1", S, V, H, N)]
        rows_e.append(f"| {lab} | {S} | {mm_['M1']['coverage']:.2f} | {mm_['M1']['bias']:.3f} | {mm_['M3_te0']['coverage']:.2f} | {mm_['M3_te0']['reject']:.2f} | {mm_['M5b']['coverage']:.2f} |".replace("-0.000", "0.000"))
raw["table_s13_e_rows"] = rows_e
raw["m1_bias_mcse_e"] = m("D1", "S2e", "M1")["bias_mcse"]
rec = dict(source=OUT.name, raw=raw, fails=fails,
           cells_sha256=hashlib.sha256("".join(sorted(hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
                                                      for p in glob.glob(str(OUT / "cells" / "cell_*.json")))).encode()).hexdigest())
(HERE / "s86_values.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")
for k, v in raw.items():
    if not k.startswith("table"):
        print(k, v if not isinstance(v, float) else round(v, 5))
print("실패", fails)
sys.exit(1 if fails else 0)
