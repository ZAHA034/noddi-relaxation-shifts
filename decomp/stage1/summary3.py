"""summary3 in stage1 computes the decomp stage 3 summary.json tables from x*.json, cards.json and run records, without hand copying."""
from __future__ import annotations

import collections
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C  # noqa: E402

OUT3 = C.out_dir(3)
Q = {"Q0": "같은 장비", "Q1": "TE 만", "Q2": "TE + TR + 설계", "Q3": "설계만", "Q4": "TR 만"}


def rng(xs):
    xs = [x for x in xs if x is not None]
    return [min(xs), max(xs)] if xs else None


def load(name, key):
    return json.load(open(OUT3 / name, encoding="utf-8"))[key]


def run_record():
    log = C.WORK / "logs"
    lines = (log / f"regen3{C.TAG}.log").read_text(encoding="utf-8").strip().splitlines() if (log / f"regen3{C.TAG}.log").exists() else []
    la = [float(ln.split()[1]) for ln in (log / f"regen3{C.TAG}_loadavg.log").read_text().splitlines() if ln.strip()] if (log / f"regen3{C.TAG}_loadavg.log").exists() else []
    head = subprocess.run(["git", "rev-parse", "--short=7", "HEAD"], cwd=C.RESULTS.parents[1], capture_output=True, text=True).stdout.strip()
    return dict(log=lines, loadavg_1min_range=rng(la), head_at_summary=head, gpu="none (CPU 만)", fitter=C.TAG_FITTER[C.TAG])


def main():
    x13 = load("x13.json", "X13")
    x10 = load("x10.json", "X10")
    x7 = load("x7.json", "X7")
    x8 = load("x8.json", "X8")
    x9 = load("x9.json", "X9")
    x12 = load("x12.json", "X12")
    cards = json.load(open(OUT3 / "cards" / "cards.json", encoding="utf-8"))

    cells24 = [k for k in x10["확정"][0] if k.endswith("σ")]
    s10 = {}
    for ver, rows in x10.items():
        allm = [r["pair"] for r in rows if all(r[k] == "일치" for k in cells24)]
        others = sorted(sum(r[k] == "일치" for k in cells24) for r in rows if r["pair"] not in allm)
        s10[ver] = dict(n_pairs=len(rows), n_cells=len(cells24), all_cells_match=allm, others_match_cells_range=rng(others))

    s7 = {}
    for key, v in x7.items():
        parts = key.split("|")
        lv, feat, q = parts[0], parts[-2], parts[-1]
        d = s7.setdefault(lv, {}).setdefault(feat, {})
        for ax in ("NDI", "FISO"):
            e = d.setdefault(ax, dict(saturated=0, n=0, unident_by_Q={}))
            e["n"] += 1
            e["saturated"] += bool(v[f"{ax}_saturated"])
            e["unident_by_Q"].setdefault(q, []).append(v[f"{ax}_unident"])
    for lv in s7.values():
        for feat in lv.values():
            for e in feat.values():
                e["unident_by_Q"] = {q: rng(xs) for q, xs in sorted(e["unident_by_Q"].items())}

    s8 = dict(HPI_pass=x8["HPI_pass"], frac_cells_ok=x8["frac_cells_ok"], width_ok=x8["width_ok"], by_feature_same_generator={})
    for feat in ("frac", "S", "C3", "C7"):
        sel = [v for k, v in x8["cells"].items() if k.split("|")[2] == feat and k.split("|")[0] == k.split("|")[1]]
        if sel:
            s8["by_feature_same_generator"][feat] = dict(
                n_cells=len(sel), cov_NDI=rng([v["cov_NDI"] for v in sel]), cov_FISO=rng([v["cov_FISO"] for v in sel]),
                width_NDI_median=rng([v["width_NDI_median"] for v in sel]), width_FISO_median=rng([v["width_FISO_median"] for v in sel]),
                empty_frac=rng([v["empty_frac"] for v in sel]))

    s9 = {q: dict(name=Q[q], in_fraction_cone=rng([v["xy_dir_in_fraction_cone"] for k, v in x9.items() if k.endswith("|" + q)]),
                  ODI_within_2sd=rng([v["ODI_within_2sd"] for k, v in x9.items() if k.endswith("|" + q)])) for q in ("Q1", "Q2", "Q3", "Q4")}

    s12 = {}
    for des in ("balanced", "confounded"):
        for ax in ("F", "N"):
            sel = [v for k, v in x12.items() if k.split("|")[2] == ax and k.split("|")[3] == des]
            s12[f"{ax}|{des}"] = dict(ratio_median=rng([v["ratio_median"] for v in sel]), ratio_raw_median=rng([v["ratio_raw_median"] for v in sel]),
                                     labels=dict(collections.Counter(v["label"] for v in sel)), labels_raw=dict(collections.Counter(v["label_raw"] for v in sel)),
                                     errors=sum(v["errors"] for v in sel))
    s12["Q2|confounded|raw"] = rng([v["ratio_raw_median"] for k, v in x12.items() if k.split("|")[1] == "Q2" and k.split("|")[3] == "confounded"])

    on = [c for c in cards if c["set"] == "ON-Harmony"]
    s11 = dict(n_cards=len(cards), by_set=dict(collections.Counter(c["set"] for c in cards)),
               onharmony_match_cells=sorted(collections.Counter(c["n_match_cells"] for c in on).items()) if on and "n_match_cells" in on[0] else None)

    doc = dict(spec="decomp/SPEC_v3.md(1af2a53 · §13 D-2)", code="82be38f + 출력 경로(ab73080)", tag=C.TAG, run=run_record(),
               X13=x13, X10=s10, X7=s7, X8=s8, X9=s9, X12=s12, X11_cards=s11,
               deviations=["D-1: F1 미실행 → F1 판정은 갈래 _F1 로 따로(decomp/F1_RULE.md)",
                           "D-2: 3 단계 결과 재생성 — 원 출력 미커밋 · 2026-10-23 재부팅으로 소실(SPEC_v3 §13)"] if C.TAG == "" else ["F1 판(SPEC_v3 §0 예외)"])
    path = C.guard_write(OUT3 / "summary.json")
    json.dump(doc, open(path, "w"), ensure_ascii=False, indent=1, default=float)
    print("씀", path)


if __name__ == "__main__":
    main()
