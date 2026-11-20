"""Further values (no new simulation: stored cell metrics copied and counted; a bulk T2 check from the signal-model formula only), each with a check that must be able to fail."""
import csv, glob, hashlib, json, math, pathlib, sys

HERE = pathlib.Path(__file__).resolve().parent
OUT = pathlib.Path(sys.argv[1])
cells = {}
for p in sorted(glob.glob(str(OUT / "cells" / "cell_*.json"))):
    d = json.loads(pathlib.Path(p).read_text(encoding="utf-8"))
    c = d["cell"]
    cells[(c["D"], c["S"], float(c["V"]), c["H"], c["N"])] = d
if len(cells) != 288:
    raise SystemExit(f"칸 수 {len(cells)} ≠ 288")


def met(D, S, k, H="H0"):
    v = cells[(D, S, 0.0, H, "bal")]["metrics"][k]
    return None if v.get("defined") is False else v


fails = []
led = {r["ID"]: r for r in csv.DictReader(open(HERE.parent / "numbers.csv", encoding="utf-8"))}
n18, n19 = int(float(led["N-018"]["원 값"])), int(float(led["N-019"]["원 값"]))
raw = dict(odi_pairs_beyond=n19 - n18)

# standard ComBat · fixed site term: maximum coverage in S2b · S2e of all designs (D0 · D1 · D2)
raw["m1_cov_s2_max"] = max(met(D, S, "M1")["coverage"] for D in ("D0", "D1", "D2") for S in ("S2b", "S2e"))
raw["m2_cov_s2_max"] = max(met(D, S, "M2")["coverage"] for D in ("D0", "D1", "D2") for S in ("S2b", "S2e"))
raw["m1_bias_d0_minus_d1_maxabs"] = max(abs(met("D0", S, k)["bias"] - met("D1", S, k)["bias"]) for S in ("S2b", "S2e") for k in ("M1", "M2"))
if raw["m1_cov_s2_max"] >= 0.005 or raw["m2_cov_s2_max"] >= 0.005:
    fails.append("S2 칸 포함률이 0.00 으로 반올림되지 않음")

# within-site bias against the 30 ms change
s75 = json.loads((HERE / "s75_values.json").read_text(encoding="utf-8"))
te_c = sum(s75["power"]["tes_mid"]) / len(s75["power"]["tes_mid"])
raw["te_c_ms"] = te_c
raw["bias_over_change30"] = te_c / 30.0
alt = abs(met("D0", "S2e", "M2")["bias"]) / abs(float(led["N-840"]["원 값"]))
raw["bias_over_change30_from_sim"] = alt
if abs(alt - raw["bias_over_change30"]) > 0.02:
    fails.append(f"비 {raw['bias_over_change30']:.3f} 와 시뮬레이션 비 {alt:.3f} 가 다름")
raw["tes_mid"] = s75["power"]["tes_mid"]

# group effect tied to the vendor (H2) · S0
for D in ("D1", "D2"):
    m3 = met(D, "S0", "M3_te0", "H2")
    raw[f"{D.lower()}_h2_s0_m3_bias"] = m3["bias"]
    raw[f"{D.lower()}_h2_s0_m3_cov"] = m3["coverage"]
    raw[f"{D.lower()}_h2_s0_slope_reject"] = met(D, "S0", "M3_slope", "H2")["reject"]

# bulk T2 check (signal model · the controls of 2.4: ν 0.46 · T2,in 82.5 · T2,en 52 ms)
NU, TIN, TEN, P = 0.46, 82.5, 52.0, 0.0984


def bulk_t2(tin, ten, tes):
    xs = [t for t in tes]
    ys = [math.log(NU * math.exp(-t / tin) + (1 - NU) * math.exp(-t / ten)) for t in xs]
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    return -1.0 / slope


def ndi_app(tin, ten, te):
    a, b = NU * math.exp(-te / tin), (1 - NU) * math.exp(-te / ten)
    return a / (a + b)


bulk = {}
for name, tes in (("te70_100", [70 + 5 * i for i in range(7)]), ("te10_120", [10 + 10 * i for i in range(12)])):
    b0 = bulk_t2(TIN, TEN, tes)
    en_only = bulk_t2(TIN, TEN * (1 + P), tes) / b0 - 1
    both = bulk_t2(TIN * (1 + P), TEN * (1 + P), tes) / b0 - 1          # recorded only (fixed range)
    ident = bulk_t2(TIN * (1 + P), TEN * (1 + P), [t * (1 + P) for t in tes]) / b0 - 1
    lo, hi = P, 5.0                                 # extra-neurite lengthening L: the value at which the bulk increase = 9.84 % (bisection)
    for _ in range(200):
        mid = (lo + hi) / 2
        if bulk_t2(TIN, TEN * (1 + mid), tes) / b0 - 1 < P:
            lo = mid
        else:
            hi = mid
    L = (lo + hi) / 2
    drops = {te: ndi_app(TIN, TEN, te) - ndi_app(TIN, TEN * (1 + L), te) for te in (70, 100)}
    drops_p = {te: ndi_app(TIN, TEN, te) - ndi_app(TIN, TEN * (1 + P), te) for te in (70, 100)}
    bulk[name] = dict(bulk_change_en_only=en_only, bulk_change_both=both, identity_scaled_range=ident, en_lengthening_matching_bulk=L,
                      ndi_drop_matched=drops, ndi_drop_en_9_84=drops_p)
    if not en_only < P:
        fails.append(f"{name}: 신경돌기 밖만 늘린 bulk 증가 {en_only:.4f} 가 9.84 % 보다 작지 않음")
    if abs(ident - P) > 1e-9:
        fails.append(f"{name}: 대조(척도 항등식) bulk 증가 {ident:.9f} ≠ 0.0984")
    if not all(drops[t] > drops_p[t] for t in drops):
        fails.append(f"{name}: 맞춘 늘림의 NDI 하강이 더 크지 않음")
for te, want in ((70, float(led["N-923"]["원 값"])), (100, float(led["N-924"]["원 값"]))):
    got = bulk["te70_100"]["ndi_drop_en_9_84"][te]
    if abs(got - want) > 5e-4:
        fails.append(f"TE {te}: 9.84 % 하강 {got:.4f} 이 장부 {want:.4f} 와 다름")
raw["bulk_t2"] = bulk

# rows of table S13 (c) · (d) (strings · put into the supplement as they are)
rows_c = []
for D in ("D0", "D1", "D2"):
    for S in ("S0", "S1b", "S1e", "S2b", "S2e", "S2i", "S3b", "S3e"):
        cells_ = []
        for k in ("M1c", "M1r", "M4_te0", "M5a"):
            v = met(D, S, k)
            cells_ += ["—", "—"] if v is None else [f"{v['coverage']:.2f}", f"{v['bias']:.3f}".replace("-0.000", "0.000")]
        rows_c.append(f"| {D} | {S} | " + " | ".join(cells_) + " |")
rows_d = []
for D in ("D1", "D2"):
    m3, sl, m1, m5 = met(D, "S0", "M3_te0", "H2"), met(D, "S0", "M3_slope", "H2"), met(D, "S0", "M1", "H2"), met(D, "S0", "M5b", "H2")
    rows_d.append(f"| {D} | S0 | {sl['reject']:.2f} | {m3['coverage']:.2f} | {m3['bias']:.3f} | {m1['coverage']:.2f} | {m1['bias']:.3f} | "
                  f"{m5['coverage']:.2f} | {m5['bias']:.3f} |")
raw["table_s13_c_rows"], raw["table_s13_d_rows"] = rows_c, rows_d

rec = dict(source=OUT.name, raw=raw, fails=fails,
           cells_sha256=hashlib.sha256("".join(sorted(hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
                                                      for p in glob.glob(str(OUT / "cells" / "cell_*.json")))).encode()).hexdigest())
(HERE / "s84_values.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")
print(json.dumps({k: v for k, v in raw.items() if not k.startswith("table")}, ensure_ascii=False, indent=0)[:2500])
print("실패", fails)
sys.exit(1 if fails else 0)
