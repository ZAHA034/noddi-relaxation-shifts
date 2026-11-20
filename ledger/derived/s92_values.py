"""Table of the relaxation sensitivity value applied to published clinical results (plan decomp/rsv_table/PLAN.md): d, d/TE and the T2 lengthenings, with the bulk T2 equivalent; known-answer tests (Section 2.4 example about 5 %; bulk values of s84 within 1e−12) and two planted faults."""
import json, math, pathlib, sys

HERE = pathlib.Path(__file__).resolve().parent
PLANT = "--plant" in sys.argv
PLANT2 = "--plant2" in sys.argv
T2IN, T2EN = 82.5, 52.0
DR2C = 1.0 / T2EN - 1.0 / T2IN
logit = lambda p: math.log(p / (1.0 - p))
fails = []


def rsv(ndi_con, ndi_pat, te):
    d = logit(ndi_pat) - logit(ndi_con)
    s = d / te                                                     # ms⁻¹ — Δ(R2,en − R2,in)
    r2en = 1.0 / T2EN + (d * te if PLANT else s)
    ext = (1.0 / r2en / T2EN - 1.0) * 100.0
    k = 1.0 / (1.0 + s / DR2C)
    both = (k - 1.0) * 100.0
    return dict(d=d, rsv_per_s=s * 1000.0, ext_pct=ext, both_pct=both)


NU = 0.46
GRIDS = {"te10_100": list(range(10, 101, 10)), "te10_120": list(range(10, 121, 10)), "te70_100": list(range(70, 101, 5))}


def bulk_t2(tin, ten, tes):
    ys = [math.log(NU * math.exp(-t / tin) + (1 - NU) * math.exp(-t / ten)) for t in tes]
    mx, my = sum(tes) / len(tes), sum(ys) / len(ys)
    return -1.0 / (sum((x - mx) * (y - my) for x, y in zip(tes, ys)) / sum((x - mx) ** 2 for x in tes))


def bulk_change_ext(ext_pct):
    """% change of bulk T2 when extra-neurite T2 is lengthened by ext_pct % (minimum · maximum over the three grids)."""
    f = 1.0 + ext_pct / 100.0
    if f <= 0:
        return float("nan"), float("nan")
    out = []
    for tes in GRIDS.values():
        b0 = bulk_t2(T2IN, T2EN, tes)
        b1 = bulk_t2(T2IN * f, T2EN, tes) if PLANT2 else bulk_t2(T2IN, T2EN * f, tes)
        out.append((b1 / b0 - 1.0) * 100.0)
    return min(out), max(out)


kb = bulk_change_ext(9.84)
S84 = json.loads((HERE / "s84_values.json").read_text(encoding="utf-8"))["raw"]["bulk_t2"]     # the values behind «3–4 %» in the note of Table S11 (same signal model · the same two grids)
for nm in ("te70_100", "te10_120"):
    f_ = 1.0 + 9.84 / 100.0
    got = (bulk_t2(T2IN * f_, T2EN, GRIDS[nm]) if PLANT2 else bulk_t2(T2IN, T2EN * f_, GRIDS[nm])) / bulk_t2(T2IN, T2EN, GRIDS[nm]) - 1.0
    if abs(got - S84[nm]["bulk_change_en_only"]) > 1e-12:
        fails.append(f"정답 시험 실패({nm}): 밖 구획 9.84 % → bulk {got:.6f} 대 s84 {S84[nm]['bulk_change_en_only']:.6f}")
if not (round(kb[0]) == 3 and round(kb[1]) == 4):
    fails.append(f"정답 시험 실패: 밖 구획 9.84 % → bulk {kb[0]:.2f}–{kb[1]:.2f} %(표 S11 주 «3–4 %» 와 정수 반올림이 달라짐)")
for nm, tes in GRIDS.items():
    p_ = 0.0984
    idt = bulk_t2(T2IN * (1 + p_), T2EN * (1 + p_), [t * (1 + p_) for t in tes]) / bulk_t2(T2IN, T2EN, tes) - 1.0
    if abs(idt - p_) > 1e-9:
        fails.append(f"척도 대조 실패({nm}): {idt}")
ka = rsv(0.618, 0.598, 90.0)
if not (4.5 <= ka["ext_pct"] <= 5.5):
    fails.append(f"정답 시험 실패: 2.4 예의 밖 구획 % = {ka['ext_pct']:.2f}(4.5–5.5 이어야 함)")
raw = dict(known_answer=ka, known_answer_bulk=list(kb))
src = pathlib.Path(sys.argv[1])
rows = json.loads(src.read_text(encoding="utf-8"))["rows"]
out = []
for r in rows:
    vals = []
    for reg in r["regions"]:
        ndi_pat = reg["ndi_pat"] if reg.get("ndi_pat") is not None else reg["ndi_con"] + reg["diff"]
        v = rsv(reg["ndi_con"], ndi_pat, r["te_ms"])
        v.update(region=reg["name"], ndi_con=reg["ndi_con"], ndi_pat=ndi_pat)
        v["bulk_ext_min"], v["bulk_ext_max"] = bulk_change_ext(v["ext_pct"])
        vals.append(v)
    rs = [v["rsv_per_s"] for v in vals]
    ex = [v["ext_pct"] for v in vals]
    bo = [v["both_pct"] for v in vals]
    out.append(dict(study=r["study"], disease=r["disease"], te_ms=r["te_ms"], n_regions=len(vals), rsv_min=min(rs), rsv_max=max(rs),
                    ext_min=min(ex), ext_max=max(ex), both_min=min(bo), both_max=max(bo), regions=vals, t2_reported=r.get("t2_reported")))
raw["rows"] = out
# rows of Table 3 (several patient groups of the same study · same disease merged into one row — plan change record 2) · flat keys read by the ledger
groups = []
for o in out:
    key = (o["study"], o["disease"])
    g = next((x for x in groups if x["key"] == key), None)
    if g is None:
        g = dict(key=key, te=o["te_ms"], rsv=[], ext=[], both=[], bulk=[], ncon=[], npat=[], ndi_con=[], ndi_pat=[], n_groups=0)
        groups.append(g)
    g["n_groups"] += 1
    srow = next(r for r in rows if r["study"] == o["study"] and r["disease"] == o["disease"] and r["te_ms"] == o["te_ms"] and len(r["regions"]) == len(o["regions"]) and r is not None and id(r) not in g.setdefault("_seen", set()))
    g["_seen"].add(id(srow)); g["ncon"].append(srow["n_con"]); g["npat"].append(srow["n_pat"])
    for v in o["regions"]:
        g["rsv"].append(v["rsv_per_s"]); g["ext"].append(v["ext_pct"]); g["both"].append(v["both_pct"])
        g["bulk"] += [v["bulk_ext_min"], v["bulk_ext_max"]]; g["ndi_con"].append(v["ndi_con"]); g["ndi_pat"].append(v["ndi_pat"])
for i, g in enumerate(groups):
    raw[f"t3_{i}_te"] = float(g["te"])
    for k in ("rsv", "ext", "both", "bulk", "ndi_con", "ndi_pat"):
        raw[f"t3_{i}_{k}_min"], raw[f"t3_{i}_{k}_max"] = min(g[k]), max(g[k])
    raw[f"t3_{i}_n"] = len(g["rsv"])
    raw[f"t3_{i}_ncon_min"], raw[f"t3_{i}_ncon_max"] = min(g["ncon"]), max(g["ncon"])
    raw[f"t3_{i}_npat_min"], raw[f"t3_{i}_npat_max"] = min(g["npat"]), max(g["npat"])
raw["t3_rows"] = [dict(i=i, study=g["key"][0], disease=g["key"][1], n_groups=g["n_groups"], n_values=len(g["rsv"])) for i, g in enumerate(groups)]
bon = json.loads(src.read_text(encoding="utf-8")).get("bonnier2014", {}).get("regions", [])
if bon:
    pct = [(b["t2_pat_ms"] / b["t2_con_ms"] - 1.0) * 100.0 for b in bon]
    raw["bonnier_pct_min"], raw["bonnier_pct_max"] = min(pct), max(pct)
print(json.dumps(dict(known_answer=ka, rows=[{k: v for k, v in o.items() if k != "regions"} for o in out]), indent=0, ensure_ascii=False)[:3000])
if PLANT or PLANT2:
    print("심은 결함 — 실패해야 함:", fails)
print("실패", fails)
if not (PLANT or PLANT2):
    (HERE / "s92_values.json").write_text(json.dumps(dict(source=src.name, raw=raw, fails=fails), ensure_ascii=False, indent=1), encoding="utf-8")
sys.exit(1 if fails else 0)
