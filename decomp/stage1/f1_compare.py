"""f1_compare sets F1 and F2 verdicts side by side for decomp stages 1, 2 and 3 with unchanged definitions, choosing neither, in stage1."""
from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C  # noqa: E402

R = C.RESULTS
GENS = C.GENERATORS


def load(stage, tag, name="summary.json"):
    return json.load(open(R / f"stage{stage}{tag}" / name, encoding="utf-8"))


def stage1(tag, F):
    s = load(1, tag)
    out = {}
    for g in GENS:
        for st in ("C3", "C5", "C7"):
            r = s["sg"][f"{F}|{g}"][st]
            out[f"{g}|{st}"] = {k: bool(r[k]["pass_"]) for k in ("SG1", "SG2", "SG3", "SG4", "SG6")} | {"SG7": bool(s["sg7"][f"{F}|{g}|{st}"]["pass_"])}
    return out


def stage2(tag):
    s = load(2, tag)
    X2 = s["X2"]
    v = lambda g, f, q, ax: X2[f"{g}|{f}|{q}"][f"{ax}_unident"]
    h1 = {}
    for g in GENS:
        for f in ("frac", "C3", "C7"):
            ndi_ok = v(g, f, "Q2", "NDI") >= 0.9
            h1[f"{g}|{f}"] = dict(Q2_FISO=v(g, f, "Q2", "FISO") >= 0.9, Q2_NDI=ndi_ok if not (g == "T-in" and f == "C7") else "따로(예외 — 수치만)",
                                  Q4_FISO=v(g, f, "Q4", "FISO") >= 0.9)
    H1 = all(all(x is True for x in d.values() if isinstance(x, bool)) for d in h1.values())
    H3 = all(v(g, f, q, ax) <= 0.10 for g in GENS for f in ("frac", "C3", "C7") for q in ("Q3", "Q0") for ax in ("NDI", "FISO"))
    H4 = all(v(g, "S", "Q2", "NDI") >= 0.9 for g in ("T-app", "T-none")) and all(v(g, "S", "Q2", "FISO") >= 0.9 for g in GENS)
    H2 = {k: bool(d["H2_pass"]) for k, d in s["X5_H2"].items()}
    X6 = {r["pair"]: (r["class_C3"], r["class_C7"]) for r in s["X6"]}
    return dict(H1=H1, H1_cells=h1, H2=all(H2.values()), H2_cells=H2, H3=H3, H4=H4, H_psi=bool(s["X3"]["H_psi_pass"]),
                H_psi_frac_ok=s["X3"]["frac_cells_ok"], X6=X6)


def stage3(tag):
    s = load(3, tag)
    return dict(HPI=bool(s["X8"]["HPI_pass"]), HPI_frac_cells_ok=s["X8"]["frac_cells_ok"], X8=s["X8"]["by_feature_same_generator"],
                X9={q: s["X9"][q]["in_fraction_cone"] for q in ("Q1", "Q2", "Q3", "Q4")}, X9_ODI={q: s["X9"][q]["ODI_within_2sd"] for q in ("Q1", "Q2", "Q3", "Q4")},
                X7_L0={f: (s["X7"]["L0"][f]["NDI"]["saturated"], s["X7"]["L0"][f]["FISO"]["saturated"], s["X7"]["L0"][f]["NDI"]["n"]) for f in ("frac", "S", "C3", "C7")},
                X7_L0_C7_NDI_unident=s["X7"]["L0"]["C7"]["NDI"]["unident_by_Q"], X10=s["X10"],
                X12={k: (v["labels"], v["ratio_median"]) for k, v in s["X12"].items() if isinstance(v, dict)}, cards=s["X11_cards"])


def same(a, b):
    return json.dumps(a, sort_keys=True, default=str) == json.dumps(b, sort_keys=True, default=str)


def main():
    two = {"F2": "", "F1": "_F1"}
    res = {F: dict(stage1=stage1(t, F), stage2=stage2(t), stage3=stage3(t)) for F, t in two.items()}
    cmp = {}
    for k in res["F2"]["stage1"]:
        for sg in res["F2"]["stage1"][k]:
            cmp[f"1|{k}|{sg}"] = dict(F2=res["F2"]["stage1"][k][sg], F1=res["F1"]["stage1"][k][sg])
    for h in ("H1", "H2", "H3", "H4", "H_psi"):
        cmp[f"2|{h}"] = dict(F2=res["F2"]["stage2"][h], F1=res["F1"]["stage2"][h])
    for p in res["F2"]["stage2"]["X6"]:
        cmp[f"2|X6|{p}"] = dict(F2=res["F2"]["stage2"]["X6"][p], F1=res["F1"]["stage2"]["X6"].get(p))
    cmp["3|HPI"] = dict(F2=res["F2"]["stage3"]["HPI"], F1=res["F1"]["stage3"]["HPI"])
    for ver in ("확정", "v2_범위"):
        cmp[f"3|X10|{ver}|all_cells_match"] = dict(F2=res["F2"]["stage3"]["X10"][ver]["all_cells_match"], F1=res["F1"]["stage3"]["X10"][ver]["all_cells_match"])
    for k, v in cmp.items():
        v["same"] = same(v["F2"], v["F1"])
    gpu = dict(GPU1_chunks=48, GPU0_chunks=110, note="덩어리별 GPU(장부 data/f1/ledger.jsonl) — 차이를 해석하지 않는다(중계 ZZZa-2)")
    doc = dict(rule="decomp/F1_RULE.md — 두 판을 나란히 · 고르지 않는다 · 같은 판정인 곳만 «적합기 수준에서»", gpu=gpu, verdicts=cmp,
               n_items=len(cmp), n_same=sum(v["same"] for v in cmp.values()), detail=res)
    path = C.guard_write(R / "F1_compare.json")
    json.dump(doc, open(path, "w"), ensure_ascii=False, indent=1, default=float)
    print("씀", path, "· 항목", doc["n_items"], "· 같음", doc["n_same"])


if __name__ == "__main__":
    main()
