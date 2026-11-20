#!/usr/bin/env python3
"""Summarizes the S4 identifiability cards by counting and copying fields of the two committed card files, with no new inference and standard library only."""
from __future__ import annotations

import ast
import copy
import datetime
import hashlib
import itertools
import json
import math
import os
import statistics
import sys
import time
from collections import Counter

T0 = time.perf_counter()
STARTED = datetime.datetime.now().astimezone().isoformat(timespec="seconds")

REPO = "<repo>"
CARDS = {arm: f"{REPO}/decomp/results/s4/cards_{arm}.json" for arm in ("c3", "c7")}
SRC = {"plan.py": f"{REPO}/decomp/s4/plan.py", "validate.py": f"{REPO}/decomp/s4/validate.py",
       "s4cfg.py": f"{REPO}/decomp/s4/s4cfg.py", "judge_v4.py": f"{REPO}/decomp/remeasure/judge_v4.py",
       "SPEC_S4.md": f"{REPO}/decomp/SPEC_S4.md", "REPORT.md": f"{REPO}/decomp/results/s4/REPORT.md"}
HERE = os.path.dirname(os.path.abspath(__file__))
OUT_JSON = os.path.join(HERE, "cards_summary.json")
OUT_MD = os.path.join(HERE, "cards_summary.md")
COMMAND = ("cd <project-notes>/manuscript_tex_2026-10-31/cards_summary && "
           "CUDA_VISIBLE_DEVICES= PYTHONDONTWRITEBYTECODE=1 <repo>/.venv/bin/python -B cards_summary.py")

ARMS = ("c3", "c7")
COMPS = ("NDI", "FISO", "ODI")
TYPES = ("match", "te", "tr", "tetr", "tetr_design")
TYPE_KO = {"match": "일치", "te": "TE 만", "tr": "TR 만", "tetr": "TE+TR", "tetr_design": "TE+TR+설계"}
# plan.py:86 — main cells L0–L3 × (SNR 30 · N_eff 50) + sensitivity cells L0 × SNR {20, 50} × N_eff {25, 100} · label form plan.py:94
CONDS = [(lv, "SNR30", "Neff50") for lv in ("L0", "L1", "L2", "L3")] + [("L0", f"SNR{s}", f"Neff{n}") for s in (20, 50) for n in (25, 100)]
REF, SENS4 = CONDS[0], CONDS[4:]
EXPECTED_TYPE_COUNTS = {"match": 1, "te": 3, "tr": 0, "tetr": 9, "tetr_design": 15}
REPORT_TOTALS = {   # REPORT.md:118 (C3) · :119 (C7) — classification (worst k) · marker · among «식별» (identified), «사전 의존» (prior-dependent) · «식별» without a marker
    "c3": {"식별": 740, "부분": 364, "비식별": 0, "재가중 불가": 0, "사전 의존": 705, "가정 의존": 198, "식별 중 사전 의존": 408, "표지 없는 식별": 332},
    "c7": {"식별": 792, "부분": 312, "비식별": 0, "재가중 불가": 0, "사전 의존": 719, "가정 의존": 197, "식별 중 사전 의존": 444, "표지 없는 식별": 348},
}
HALF_MIN = 200                                      # SPEC_S4.md:131 · validate.py:167
PD, AD, NR = "사전 의존", "가정 의존", "재가중 불가"
FLAG_NR = {"half": "사전 민감도 ×0.5 재가중 불가", "two": "사전 민감도 ×2 재가중 불가"}   # validate.py:840
CLASSES = ("식별", "부분", "비식별", NR)
# ⑤ Text of the quoted line — stops if the line number points to different text
LINE_BIND = [
    ("SPEC_S4.md", 10, "E3 카드의 대표 수치 «MDE 대표 · C 함께»"),
    ("SPEC_S4.md", 55, "중요도 가중은 «주 폭 안이면 상수 · 밖이면 0» 이다"),
    ("SPEC_S4.md", 106, "분류는 아홉 중 가장 나쁜 것**(가장 낮은 C · 가장 큰 MDE)으로 하며 아홉의 중앙값을 함께 적는다"),
    ("SPEC_S4.md", 107, "MDE_u = 2.80 × σ_u"),
    ("SPEC_S4.md", 107, "카드의 대표 수치는 MDE · C 를 함께 적는다(저자 결정 E3)."),
    ("SPEC_S4.md", 109, "식별 C ≥ 0.75"),
    ("SPEC_S4.md", 131, "남은 표본 200 미만이면 «재가중 불가»"),
    ("SPEC_S4.md", 153, "분류(k 아홉 중 가장 나쁜 것 · 중앙값 함께)"),
    ("SPEC_S4.md", 190, "한 k 라도 «재가중 불가» 면 그 칸은 «재가중 불가»"),
    ("SPEC_S4.md", 202, "×0.5 재가중에 남은 표본 200 이상이 되도록"),
    ("plan.py", 86, "cells = [(l, 30.0, 50.0) for l in S.L_LEVELS] + [(\"L0\", s, ne) for s in (20.0, 50.0) for ne in (25.0, 100.0)]"),
    ("plan.py", 94, "tags += [f\"{a}~{b}|{l}|SNR{int(snr)}|Neff{int(ne)}\"] * M"),
    ("judge_v4.py", 94, "def acq_type(pair: str) -> str:"),
    ("judge_v4.py", 101, "return \"te\" if dte else \"tr\" if dtr else (\"tetr_design\" if a[0] != b[0] else \"match\")"),
    ("s4cfg.py", 86, "중요도 가중은 이 지시함수다"),
    ("validate.py", 34, "cards=2000"),
    ("validate.py", 160, "dict(c_half=nan, n_half=0, c_two=nan)"),
    ("validate.py", 162, "return r"),
    ("validate.py", 165, "h = S.in_main(T, 0.5)"),
    ("validate.py", 166, "r[\"n_half\"] = int(h.sum())"),
    ("validate.py", 167, "if h.sum() >= 200 else nan"),
    ("validate.py", 169, "if f2.sum() >= 200 else nan"),
    ("validate.py", 189, "np.where(np.isnan(c), \"재가중 불가\""),
    ("validate.py", 212, "np.nanmedian(v)"),
    ("validate.py", 761, "n_half = np.zeros((m_, nk), int)"),
    ("validate.py", 773, "n_half[i, kk] = r[\"n_half\"]"),
    ("validate.py", 792, "C_worst=float(C.min()), MDE_worst=float(np.nanmax(E))"),
    ("validate.py", 818, "2.80 * s"),
    ("validate.py", 822, "n=int(np.isfinite(M[\"k_c\"][m, 0, a]).sum())"),
    ("validate.py", 836, "Cp = [_nanmed(M[f\"k_c_{pr}\"][m, kk, a])"),
    ("validate.py", 838, "if any(v is None for v in Cp):"),
    ("validate.py", 840, "재가중 불가"),
    ("validate.py", 843, "if e[f\"cls_{pr}\"] != wk[\"cls\"]:"),
    ("validate.py", 844, "flags.append(\"사전 의존\")"),
    ("validate.py", 845, "e[\"n_half_median\"] = float(np.median(M[\"n_half\"][m]))"),
    ("REPORT.md", 118, "{'식별': 740, '부분': 364} · 표지 {'사전 의존': 705, '가정 의존': 198}"),
    ("REPORT.md", 119, "{'식별': 792, '부분': 312} · 표지 {'사전 의존': 719, '가정 의존': 197}"),
]


def sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def read_line(path, n):
    with open(path, encoding="utf-8") as f:
        lines = f.read().split("\n")
    return lines[n - 1]


def onh_from_plan(path):
    """The ONH of plan.py (scanner → design · TE range · TR range · coil) taken literally — ast.literal_eval (0 execution)."""
    with open(path, encoding="utf-8") as f:
        tree = ast.parse(f.read(), filename=path)
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "ONH" for t in node.targets):
            return ast.literal_eval(node.value), [node.lineno, node.end_lineno]
    raise SystemExit("⛔ plan.py 에 ONH 가 없다")


def acq_type(onh, pair):
    """Exactly as judge_v4.py:94–101 — TE · TR count as "different" only when the sums of the card ranges (= the midpoints) differ."""
    a, b = (onh[u] for u in pair.split("~"))
    dte = abs(sum(a[1]) - sum(b[1])) > 1e-9
    dtr = abs(sum(a[2]) - sum(b[2])) > 1e-9
    if dte and dtr:
        return "tetr_design" if a[0] != b[0] else "tetr"
    return "te" if dte else "tr" if dtr else ("tetr_design" if a[0] != b[0] else "match")


def type_counts(onh, pairs):
    out = {t: 0 for t in TYPES}
    for p in pairs:
        out[acq_type(onh, p)] += 1
    return out


def classify(c, hi=0.75, lo=0.25):
    """Exactly as validate.py:188–189 (a single value)."""
    if c is None or math.isnan(c):
        return NR
    return "식별" if c >= hi else ("비식별" if c < lo else "부분")


def same(a, b):
    return a is not None and b is not None and math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-15)


def argmin(v):
    return min(range(len(v)), key=lambda i: v[i])


def argmax(v):
    return max(range(len(v)), key=lambda i: v[i])


def entry_problems(e):
    """Recounts the relations among the fields of one cell × direction entry with the rules of validate.py — list of mismatches (passes if empty)."""
    C, sd, E = e["C_by_k"], e["sd_by_k"], e["MDE_by_k"]
    if len(C) != 9 or len(sd) != 9 or len(E) != 9 or any(v is None for v in C + sd + E):
        return ["k 아홉 값이 다 있지 않다"]
    p = []
    if any(not same(m, 2.80 * s) for m, s in zip(E, sd)):
        p.append("MDE_by_k ≠ 2.80 × sd_by_k (validate.py:818)")
    for fld, want in (("C_worst", min(C)), ("MDE_worst", max(E)), ("C_median", statistics.median(C)), ("MDE_median", statistics.median(E))):
        if not same(e[fld], want):
            p.append(f"{fld} (validate.py:792)")
    if e["cls_by_k"] != [classify(c) for c in C]:
        p.append("cls_by_k (validate.py:791)")
    if e["cls"] != classify(min(C)):
        p.append("cls (validate.py:793)")
    if e["cls_sens"] != classify(min(C), 0.5, 0.1):
        p.append("cls_sens (validate.py:793)")
    if e["assumption_dependent"] != (len(set(e["cls_by_k"])) > 1) or (AD in e["flags"]) != e["assumption_dependent"]:
        p.append("가정 의존 (validate.py:794 · :832–833)")
    pd = False
    for pr in ("half", "two"):
        Cp = e[f"C_{pr}_by_k"]
        nr = len(Cp) != 9 or any(v is None for v in Cp)
        want = NR if nr else classify(min(Cp))
        if e[f"cls_{pr}"] != want:
            p.append(f"cls_{pr} (validate.py:838–842)")
        if (FLAG_NR[pr] in e["flags"]) != nr:
            p.append(f"{FLAG_NR[pr]} (validate.py:840)")
        pd = pd or (not nr and want != e["cls"])
    if (PD in e["flags"]) != pd:
        p.append("사전 의존 (validate.py:843–844)")
    if not set(e["flags"]) <= {PD, AD, *FLAG_NR.values()}:
        p.append("모르는 표지")
    if not isinstance(e.get("n_half_median"), (int, float)):
        p.append("n_half_median 없음 (validate.py:845)")
    return p


def totals(d):
    t = {k: 0 for k in REPORT_TOTALS["c3"]}
    for ent in d["cells"].values():
        for e in ent.values():
            t[e["cls"]] += 1
            t[PD] += PD in e["flags"]
            t[AD] += AD in e["flags"]
            t["식별 중 사전 의존"] += e["cls"] == "식별" and PD in e["flags"]
            t["표지 없는 식별"] += e["cls"] == "식별" and PD not in e["flags"]
    return t


def failure_controls(d, onh, pairs, key):
    """⑥ Planted defects — clean items pass, and for each defect the check must catch it."""
    base = d["cells"][key]["NDI"]
    other = {c: [x for x in CLASSES[:3] if x != c][0] for c in CLASSES}
    muts = {
        "cls 를 다른 분류로": lambda e: e.__setitem__("cls", other[e["cls"]]),
        "«사전 의존» 표지 뒤집기": lambda e: e.__setitem__("flags", [f for f in e["flags"] if f != PD] if PD in e["flags"] else e["flags"] + [PD]),
        "MDE_worst × 1.01": lambda e: e.__setitem__("MDE_worst", e["MDE_worst"] * 1.01),
        "C_half_by_k 를 다른 분류 쪽 값으로": lambda e: e.__setitem__("C_half_by_k", [0.5 if e["cls_half"] == "식별" else 0.9] * 9),
    }
    res = {"깨끗한 항목 통과": not entry_problems(base)}
    for name, f in muts.items():
        e = copy.deepcopy(base)
        f(e)
        res[name + " → 잡힘"] = bool(entry_problems(e))
    onh2 = copy.deepcopy(onh)
    onh2["OXF4GEP"] = (onh2["OXF4GEP"][0], onh2["OXF4GEP"][1], onh2["NOT4GEP"][2], onh2["OXF4GEP"][3])
    res["유형 규칙: OXF4GEP 의 TR 범위 ← NOT4GEP 의 것 → 유형 수가 달라져 잡힘"] = type_counts(onh2, pairs) != EXPECTED_TYPE_COUNTS
    return res


def ck(cond):
    return "|".join(cond)


def ref_detail(e, K):
    """One entry of the reference cell — the worst k (minimum C) of the main prior · ×2 · ×0.5, with that k · classification · marker · the changed comparison."""
    C, C2, Ch, E = e["C_by_k"], e["C_two_by_k"], e["C_half_by_k"], e["MDE_by_k"]
    half_ok = len(Ch) == 9 and all(v is not None for v in Ch)
    i, i2, ie = argmin(C), argmin(C2), argmax(E)
    ih = argmin(Ch) if half_ok else None
    changed = [lab for pr, lab in (("two", "×2"), ("half", "×0.5")) if e[f"cls_{pr}"] not in (NR, e["cls"])]
    return {
        "C_worst": e["C_worst"], "k_C_worst": K[i], "C_median": e["C_median"],
        "C_two_worst": C2[i2], "k_C_two_worst": K[i2],
        "C_half_worst": Ch[ih] if half_ok else None, "k_C_half_worst": K[ih] if half_ok else None,
        "n_half_median": e["n_half_median"],
        "cls": e["cls"], "cls_two": e["cls_two"], "cls_half": e["cls_half"],
        "prior_dependent": PD in e["flags"], "assumption_dependent": AD in e["flags"],
        "class_changed_by": "both" if len(changed) == 2 else (changed[0] if changed else "none"),
        "MDE_worst": e["MDE_worst"], "k_MDE_worst": K[ie], "MDE_median": e["MDE_median"],
        "n_finite_obs_k0": e["n"],
    }


def mm(v):
    v = [x for x in v if x is not None]
    return [min(v), max(v)] if v else None


def type_summary(d, comp, ps):
    if not ps:
        return {"n_pairs": 0, "note": "쌍 0 — 값 없음"}
    ref = {p: d["cells"][f"{p}|{ck(REF)}"][comp] for p in ps}
    sens = {p: {ck(c): d["cells"][f"{p}|{ck(c)}"][comp] for c in SENS4} for p in ps}
    sv = [e for s in sens.values() for e in s.values()]
    return {
        "n_pairs": len(ps), "pairs": ps,
        "ref_MDE_worst_k_range_over_pairs": mm([e["MDE_worst"] for e in ref.values()]),
        "ref_MDE_median_k_range_over_pairs": mm([e["MDE_median"] for e in ref.values()]),
        "ref_C_worst_k_range_over_pairs": mm([e["C_worst"] for e in ref.values()]),
        "ref_C_median_k_range_over_pairs": mm([e["C_median"] for e in ref.values()]),
        "ref_cls_counts": {c: sum(e["cls"] == c for e in ref.values()) for c in CLASSES},
        "ref_identified_without_prior_dependent": sum(e["cls"] == "식별" and PD not in e["flags"] for e in ref.values()),
        "ref_assumption_dependent": sum(AD in e["flags"] for e in ref.values()),
        "ref_prior_dependent": sum(PD in e["flags"] for e in ref.values()),
        "sens4_MDE_worst_k_range_over_pairs_and_conditions": mm([e["MDE_worst"] for e in sv]),
        "sens4_MDE_median_k_range_over_pairs_and_conditions": mm([e["MDE_median"] for e in sv]),
        "per_pair_ref": {p: {k: e[k] for k in ("MDE_worst", "MDE_median", "C_worst", "C_median", "cls", "flags")} for p, e in ref.items()},
        "per_pair_sens4_MDE_worst": {p: {c: e["MDE_worst"] for c, e in s.items()} for p, s in sens.items()},
        "per_pair_sens4_MDE_median": {p: {c: e["MDE_median"] for c, e in s.items()} for p, s in sens.items()},
    }


def level_summary(d, comp, pairs, types):
    out = {}
    for c in CONDS[:4]:
        es = {p: d["cells"][f"{p}|{ck(c)}"][comp] for p in pairs}
        idf = {p: e["cls"] == "식별" for p, e in es.items()}
        idf_np = {p: e["cls"] == "식별" and PD not in e["flags"] for p, e in es.items()}
        out[c[0]] = {
            "cell": ck(c), "n_pairs": len(es),
            "cls_counts": {k: sum(e["cls"] == k for e in es.values()) for k in CLASSES},
            "identified_without_prior_dependent": sum(idf_np.values()),
            "prior_dependent": sum(PD in e["flags"] for e in es.values()),
            "assumption_dependent": sum(AD in e["flags"] for e in es.values()),
            "identified_by_type": {t: sum(idf[p] for p in pairs if types[p] == t) for t in TYPES},
            "identified_without_prior_dependent_by_type": {t: sum(idf_np[p] for p in pairs if types[p] == t) for t in TYPES},
        }
    return out


def half_nr(d):
    """Review item 1b — ×0.5 «reweighting not possible» (재가중 불가): the card marker and n_half_median < 200 are counted separately."""
    cells = d["cells"]
    nh, uniform = {}, 0
    for cl, ent in cells.items():
        vals = {e["n_half_median"] for e in ent.values()}
        uniform += len(vals) == 1
        nh[cl] = max(vals)
    lt = sorted(cl for cl, v in nh.items() if v < HALF_MIN)
    lab = lambda e, pr: e[f"cls_{pr}"] == NR or FLAG_NR[pr] in e["flags"]
    by_cond = {}
    for c in CONDS:
        v = sorted(nh[cl] for cl in cells if cl.endswith("|" + ck(c)))
        by_cond[ck(c)] = {"n_cells": len(v), "lt200": sum(x < HALF_MIN for x in v), "min": v[0], "median": statistics.median(v), "max": v[-1]}
    vals = sorted(nh.values())
    return {
        "unit": "칸 = 카드 쌍 × 조건(28 × 8 = 224) — n_half_median 은 칸마다 값 하나(관측 100 × k 9 = 900 값의 중앙값 · validate.py:845 · 그 칸의 모든 방향 항목에 같은 값)",
        "n_cells": len(cells), "n_cell_direction_entries": sum(len(ent) for ent in cells.values()),
        "cells_with_identical_n_half_median_across_directions": uniform,
        "card_label_x05_not_reweightable_cells": sum(any(lab(e, "half") for e in ent.values()) for ent in cells.values()),
        "card_label_x05_not_reweightable_entries": sum(lab(e, "half") for ent in cells.values() for e in ent.values()),
        "entries_with_None_in_C_half_by_k": sum(any(v is None for v in e["C_half_by_k"]) for ent in cells.values() for e in ent.values()),
        "card_label_x2_not_reweightable_entries": sum(lab(e, "two") for ent in cells.values() for e in ent.values()),
        "cells_n_half_median_lt200": len(lt),
        "cells_n_half_median_ge200": {cl: v for cl, v in nh.items() if v >= HALF_MIN},
        "same_count_as_cell_x_component_NDI_FISO_ODI": [sum(c in cells[cl] for cl in lt for c in COMPS), sum(c in ent for ent in cells.values() for c in COMPS)],
        "same_count_as_cell_x_direction": [sum(len(cells[cl]) for cl in lt), sum(len(ent) for ent in cells.values())],
        "n_half_median_over_cells": {"min": vals[0], "median": statistics.median(vals), "max": vals[-1]},
        "n_half_median_by_condition": by_cond,
    }


NOT_COMPUTED = [
    "×0.5 재가중 뒤 남은 표본 수의 k 별 · 관측별 값 — 카드에는 칸마다 중앙값 하나(n_half_median · validate.py:845)만 있다 → not computed",
    "C_half_by_k · C_two_by_k 의 각 값에 들어간 관측 수(200 이상 남은 관측 수) — 카드에 없다(항목의 n 은 주 사전 · 첫째 k «T-in|i» 의 유한 관측 수 · validate.py:822) → not computed",
    "×2 · ×0.5 사전 아래의 σ · MDE — 카드에는 C 만 있다(validate.py:836–837) → not computed",
    "×2 · ×0.5 의 민감도 문턱(0.5 · 0.1) 분류 — cls_sens 는 주 사전에만 있다(validate.py:793) → not computed",
    "L1 – L3 의 민감도 조건(SNR 20 · 50 × N_eff 25 · 100) — 그런 칸이 없다(plan.py:86 · 민감도 칸은 L0 만) → not computed",
    "TR 만 다른 쌍 — 카드 28 쌍 가운데 0 쌍(judge_v4.py:94–101 규칙) → 값 없음",
    "관측 단위 사후 표본 · 구간 · 포함률 — 카드에 없다(카드 밖 중간 파일은 열지 않았다) → not computed",
]


def fmt(x, nd=3):
    return "—" if x is None else f"{x:.{nd}f}"


def frange(r, nd=3):
    if r is None:
        return "—"
    a, b = fmt(r[0], nd), fmt(r[1], nd)
    return a if a == b else f"{a}–{b}"


def kd(k):
    return "—" if k is None else k.replace("|", "·")


def main():
    fails = []
    hashes_in = {"cards_c3.json": sha256(CARDS["c3"]), "cards_c7.json": sha256(CARDS["c7"])}
    hashes_in.update({n: sha256(p) for n, p in SRC.items()})
    for f, n, s in LINE_BIND:                                                           # ⑤
        if s not in read_line(SRC[f], n):
            fails.append(f"줄 대조 {f}:{n} — «{s}» 없음")
    onh, onh_lines = onh_from_plan(SRC["plan.py"])
    pairs = [f"{a}~{b}" for a, b in itertools.combinations(list(onh), 2)]               # plan.py:84–85 · judge_v4.py:90–91
    types = {p: acq_type(onh, p) for p in pairs}
    tc = type_counts(onh, pairs)
    if len(pairs) != 28:
        fails.append(f"쌍 {len(pairs)} ≠ 28")
    if tc != EXPECTED_TYPE_COUNTS:                                                      # ②
        fails.append(f"획득 유형 수 {tc} ≠ {EXPECTED_TYPE_COUNTS}")
    cards = {}
    for arm in ARMS:
        with open(CARDS[arm], encoding="utf-8") as f:
            cards[arm] = json.load(f)
    want_cells = {f"{p}|{ck(c)}" for p in pairs for c in CONDS}
    n_checked = {}
    for arm in ARMS:
        d = cards[arm]
        if set(d["cells"]) != want_cells:                                               # ①
            fails.append(f"{arm}: 칸 집합이 28 × 8 과 다르다")
            continue
        if d["n"] != 100 * len(want_cells):
            fails.append(f"{arm}: n {d['n']} ≠ 22400")
        n_checked[arm] = 0
        for cl, ent in d["cells"].items():
            if len({e["n_half_median"] for e in ent.values()}) != 1:                  # validate.py:845 — must be one value per cell
                fails.append(f"{arm} {cl}: n_half_median 이 방향마다 다르다")
            for u, e in ent.items():
                fails += [f"{arm} {cl} {u}: {m}" for m in entry_problems(e)]            # ③
                n_checked[arm] += 1
        if totals(d) != REPORT_TOTALS[arm]:                                             # ④
            fails.append(f"{arm}: 합계 {totals(d)} ≠ REPORT {REPORT_TOTALS[arm]}")
    ctrl = failure_controls(cards["c3"], onh, pairs, f"{pairs[1]}|{ck(REF)}") if not fails else {}
    if ctrl and not all(ctrl.values()):                                                 # ⑥
        fails.append(f"실패 대조가 잡지 못함: {ctrl}")
    if fails:
        raise SystemExit("⛔ 대조 실패 — 아무것도 쓰지 않았다:\n" + "\n".join(fails[:60]))

    K = cards["c3"]["k_names"]
    if cards["c7"]["k_names"] != K:
        raise SystemExit("⛔ 두 팔의 k_names 가 다르다")
    res = {
        "status": "S4 카드 두 파일의 필드를 세고 옮긴 요약 · 새 추론 0 · 프로젝트 안 작업(외부 검토 아님)",
        "sources_sha256": hashes_in, "onh_from_plan_py_lines": onh_lines, "k_names": K,
        "checks": {"cell_direction_entries_checked": n_checked, "problems": 0, "type_counts": tc,
                   "type_counts_expected_from_task": EXPECTED_TYPE_COUNTS, "report_totals_match": True,
                   "line_bindings_checked": len(LINE_BIND), "failure_controls": ctrl},
        "acquisition_types": {t: [p for p in pairs if types[p] == t] for t in TYPES},
        "A_reference_cell": {}, "B_half_width_not_reweightable": {}, "C_by_arm_component_type": {},
        "identified_counts_by_level": {}, "not_computed": NOT_COMPUTED,
    }
    for arm in ARMS:
        d = cards[arm]
        A = {"cell": ck(REF), "NDI_per_pair": {}, "counts": {}}
        for comp in COMPS:
            rows = {p: ref_detail(d["cells"][f"{p}|{ck(REF)}"][comp], K) for p in pairs}
            if comp == "NDI":
                A["NDI_per_pair"] = {p: {"type": types[p], **r} for p, r in rows.items()}
            by = Counter(r["class_changed_by"] for r in rows.values())
            A["counts"][comp] = {
                "n_pairs": len(rows), "prior_dependent": sum(r["prior_dependent"] for r in rows.values()),
                "changed_by_x2_only": by["×2"], "changed_by_x05_only": by["×0.5"], "changed_by_both": by["both"],
                "not_prior_dependent_pairs": [p for p, r in rows.items() if not r["prior_dependent"]],
                "class_patterns_main_x2_x05": dict(Counter(f"{r['cls']} / {r['cls_two']} / {r['cls_half']}" for r in rows.values())),
                "assumption_dependent": sum(r["assumption_dependent"] for r in rows.values()),
            }
        res["A_reference_cell"][arm] = A
        res["B_half_width_not_reweightable"][arm] = half_nr(d)
        res["C_by_arm_component_type"][arm] = {comp: {t: type_summary(d, comp, [p for p in pairs if types[p] == t]) for t in TYPES} for comp in COMPS}
        res["identified_counts_by_level"][arm] = {comp: level_summary(d, comp, pairs, types) for comp in COMPS}

    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=1)
        f.write("\n")
    write_md(res, pairs, types)


def write_md(res, pairs, types):
    L = []
    w = L.append
    spec131 = read_line(SRC["SPEC_S4.md"], 131)
    w("# S4 식별성 카드 요약 (cards_summary)")
    w("")
    w("커밋된 카드 두 파일의 필드를 세고 옮긴 프로젝트 안 요약이다 — 외부 검토가 아니다. 새 추론 · 새 계산 0(최소 · 최대 · 개수 · 필드 사이 대조만).")
    w("표의 수치는 소수 셋째 자리 반올림이다(전체 자릿수는 `cards_summary.json`). k 이름의 «|» 는 표에서 «·» 로 적었다(예: T-in·iii = 카드의 «T-in|iii»).")
    w("")
    w("## 0. 출처 · 해시 · 실행 · 대조")
    w("")
    w("| 파일 | sha256 |")
    w("|---|---|")
    for n, h in res["sources_sha256"].items():
        w(f"| {n} | `{h}` |")
    w(f"| cards_summary.py (이 스크립트) | `{sha256(os.path.abspath(__file__))}` |")
    w(f"| cards_summary.json (출력) | `{sha256(OUT_JSON)}` |")
    w("")
    w("- 이 파일(cards_summary.md) 자신의 sha256 은 파일 안에 적을 수 없다 — 영수증에 적는다.")
    w(f"- 카드 경로: `{CARDS['c3']}` · `{CARDS['c7']}` (읽기만)")
    w(f"- 명령: `{COMMAND}`")
    w(f"- 실행 시작 {STARTED} · 벽시계 {time.perf_counter() - T0:.3f} 초 · 인터프리터 `{sys.executable}` (Python {sys.version.split()[0]})"
      f" · 작업 폴더 `{os.getcwd()}` · CUDA_VISIBLE_DEVICES={os.environ.get('CUDA_VISIBLE_DEVICES')!r}"
      f" · PYTHONDONTWRITEBYTECODE={os.environ.get('PYTHONDONTWRITEBYTECODE')!r} · sys.flags.dont_write_bytecode={sys.flags.dont_write_bytecode}")
    c = res["checks"]
    w(f"- 대조(어긋나면 출력 없이 멈춤): 칸 × 방향 항목 C3 {c['cell_direction_entries_checked']['c3']} · C7 {c['cell_direction_entries_checked']['c7']} 개의 필드 관계가 validate.py 규칙과 같음(어긋남 0)"
      f" · 획득 유형 수 {c['type_counts']} = 과제 문안 · 팔별 합계 = REPORT.md:118 · :119 · 인용 줄 {c['line_bindings_checked']} 개 글자 대조"
      f" · 실패 대조 {sum(c['failure_controls'].values())}/{len(c['failure_controls'])} 통과(깨끗한 항목 통과 1 + 심은 결함 {len(c['failure_controls']) - 1} 개 모두 잡힘).")
    w(f"- 저장소 모듈 import 0 — 카드 장비 기록 ONH 는 plan.py:{res['onh_from_plan_py_lines'][0]}–{res['onh_from_plan_py_lines'][1]} 을 ast 로 읽었다.")
    w("")
    w("## 1. 획득 유형 (규칙 judge_v4.py:94–101 · 카드 범위 plan.py:17–22)")
    w("")
    for t in TYPES:
        ps = res["acquisition_types"][t]
        w(f"- {TYPE_KO[t]} {len(ps)} 쌍: " + (" · ".join(ps) if ps else "없음"))
    w("")
    w("## 2. 기준 칸(L0 · SNR 30 · N_eff 50) NDI — 쌍별 (검토 항목 1a)")
    w("")
    w("열: C = 그 사전 아래 k 아홉 가운데 가장 낮은 C(괄호 안 = 그 k) — 주 사전 `C_worst`(= min `C_by_k`) · ×2 min `C_two_by_k` · ×0.5 min `C_half_by_k`."
      " n½ = `n_half_median`(관측 100 × k 9 의 ×0.5 재가중 뒤 남은 표본 수의 중앙값 · k 마다 주 지지 표본 2,000 개 가운데 · 가중이 지지 지시함수라 유효 표본 수 = 남은 수 · SPEC_S4.md:55 · s4cfg.py:86 · validate.py:165)."
      " 분류 = `cls` / `cls_two` / `cls_half`. 바꾼 비교 = 주 사전의 분류와 다른 쪽.")
    for arm in ARMS:
        A = res["A_reference_cell"][arm]
        w("")
        w(f"### {arm.upper()}")
        w("")
        w("| 쌍 | 유형 | C 주 (k) | C ×2 (k) | C ×0.5 (k) | n½ | 분류 주 / ×2 / ×0.5 | 사전 의존 | 바꾼 비교 |")
        w("|---|---|---|---|---|---|---|---|---|")
        for p in pairs:
            r = A["NDI_per_pair"][p]
            w(f"| {p} | {TYPE_KO[r['type']]} | {fmt(r['C_worst'])} ({kd(r['k_C_worst'])}) | {fmt(r['C_two_worst'])} ({kd(r['k_C_two_worst'])}) | "
              f"{fmt(r['C_half_worst'])} ({kd(r['k_C_half_worst'])}) | {r['n_half_median']:g} | {r['cls']} / {r['cls_two']} / {r['cls_half']} | "
              f"{'예' if r['prior_dependent'] else '아니오'} | {r['class_changed_by'].replace('none', '—').replace('both', '둘 다')} |")
        w("")
        for comp in COMPS:
            q = A["counts"][comp]
            npd = q["not_prior_dependent_pairs"]
            npd_s = "없음" if not npd else (f"{len(npd)} 쌍 모두" if len(npd) == q["n_pairs"] else " · ".join(npd))
            w(f"- {arm.upper()} {comp}: «사전 의존» {q['prior_dependent']}/{q['n_pairs']} · 바꾼 비교 ×2 만 {q['changed_by_x2_only']} · ×0.5 만 {q['changed_by_x05_only']} · 둘 다 {q['changed_by_both']}"
              f" · «사전 의존» 아닌 쌍: {npd_s} · 분류 꼴(주 / ×2 / ×0.5): "
              + " · ".join(f"{k} {v}" for k, v in q["class_patterns_main_x2_x05"].items()))
    w("")
    w("## 3. ×0.5 사전의 «재가중 불가» (검토 항목 1b)")
    w("")
    w("사양 원문(SPEC_S4.md:131 · 파일에서 그대로 읽음):")
    w("")
    w("> " + spec131)
    w("")
    for arm in ARMS:
        b = res["B_half_width_not_reweightable"][arm]
        ge = "; ".join(f"`{cl}` {v:g}" for cl, v in b["cells_n_half_median_ge200"].items()) or "없음"
        nm = b["n_half_median_over_cells"]
        w(f"- {arm.upper()}: 카드 표지 «재가중 불가»(×0.5 · `cls_half` 또는 표지 «{FLAG_NR['half']}») {b['card_label_x05_not_reweightable_cells']} 칸"
          f" (칸 × 방향 {b['card_label_x05_not_reweightable_entries']}/{b['n_cell_direction_entries']} · `C_half_by_k` 에 None 인 항목 {b['entries_with_None_in_C_half_by_k']} · ×2 표지 {b['card_label_x2_not_reweightable_entries']})"
          f" · `n_half_median` < 200 인 칸 {b['cells_n_half_median_lt200']}/{b['n_cells']} (200 이상: {ge})"
          f" · 같은 수를 칸 × 성분(NDI · FISO · ODI)으로 {b['same_count_as_cell_x_component_NDI_FISO_ODI'][0]}/{b['same_count_as_cell_x_component_NDI_FISO_ODI'][1]}"
          f" · 칸 × 방향으로 {b['same_count_as_cell_x_direction'][0]}/{b['same_count_as_cell_x_direction'][1]}"
          f" · 칸 사이 `n_half_median` 최소 {nm['min']:g} · 중앙값 {nm['median']:g} · 최대 {nm['max']:g}"
          f" · 방향 사이 값이 같은 칸 {b['cells_with_identical_n_half_median_across_directions']}/{b['n_cells']}.")
    rc = {arm: res["B_half_width_not_reweightable"][arm]["n_half_median_by_condition"][ck(REF)] for arm in ARMS}
    x05 = {arm: {comp: res["A_reference_cell"][arm]["counts"][comp]["changed_by_x05_only"] + res["A_reference_cell"][arm]["counts"][comp]["changed_by_both"]
                 for comp in COMPS} for arm in ARMS}
    w("- 기준 칸(L0 · SNR 30 · N_eff 50): `n_half_median` < 200 인 쌍 " + " · ".join(
        f"{arm.upper()} {rc[arm]['lt200']}/{rc[arm]['n_cells']} (범위 {rc[arm]['min']:g}–{rc[arm]['max']:g})" for arm in ARMS)
      + " — ×0.5 비교가 분류를 바꾼 기준 칸 항목(" + " / ".join(f"{arm.upper()} " + " · ".join(f"{c} {x05[arm][c]}" for c in COMPS) for arm in ARMS)
      + ")에서 `C_half_by_k` 의 각 값은 그 k 의 관측 100 개 가운데 200 이상 남은 관측만의 중앙값이고, 아홉 k 를 합친 관측 × k 900 값 가운데 200 이상인 것은 절반 이하다(중앙값 < 200 에서 바로 나오는 추론 · k 별 개수는 카드에 없다).")
    w("- 세는 단위: 칸 = 카드 쌍 × 조건(28 × 8 = 224 · plan.py:84–86). `n_half_median` 은 칸마다 하나다(validate.py:761 · :773 의 관측 × k 배열 900 값의 중앙값을 그 칸의 모든 방향 항목에 같은 값으로 적음 · validate.py:845) — 성분마다 따로 세는 수가 아니다.")
    w("- 구현 줄: validate.py:165–167 (관측 × k 마다 주 지지 표본 2,000 개 — validate.py:34 — 가운데 ×0.5 지지 안의 수 `n_half` · 200 이상일 때만 c_half, 아니면 nan)"
      " · :160–162 (주 지지 표본이 목표에 못 이른 관측은 n_half = 0 · c_half = nan) · :836 + :210–212 (칸의 k 별 C_half = 관측들의 nanmedian — nan 관측은 표지 없이 빠지고, 모든 관측이 nan 일 때만 None)"
      " · :838–840 (어느 k 라도 None 이면 `cls_half` «재가중 불가» + 표지) · :841–844 (아니면 min C_half 의 분류가 주 사전과 다르면 «사전 의존») · :845 (`n_half_median`).")
    w("- ×0.5 가 재가중 불가일 때 «사전 의존»: 칸 수준(어느 k 에서 모든 관측이 200 미만)이면 ×0.5 비교를 건너뛰고 ×2 비교만으로 정하며 표지 «사전 민감도 ×0.5 재가중 불가» 를 단다(validate.py:835–844 — 비교는 else 가지에만 있다)."
      " 관측 수준(그 관측 × k 만 200 미만)이면 표지 없이 그 값을 빼고 200 이상 남은 관측의 중앙값으로 비교한다. 이 카드들에서 칸 수준 가지는 위 수대로 0 번 걸렸다.")
    w("- 사양인가 구현인가: SPEC_S4.md:131 은 문턱(200)과 이름만 정한다 — (가) 문턱을 관측 × k 에 거는지 칸에 거는지 (나) 일부 관측만 재가중 불가일 때 칸의 C_half 를 어떻게 내는지 (다) ×0.5 가 재가중 불가일 때 «사전 의존» 을 어떻게 정하는지는 적혀 있지 않다."
      " §13 4(:190)의 «한 k 라도 «재가중 불가» 면 그 칸은 «재가중 불가»» 도 칸 · k 가 언제 재가중 불가인지는 정하지 않는다. 위 세 동작은 validate.py 의 구현 선택이다."
      " §13 14(:202)는 카드 k 마다 2,000 개를 «×0.5 재가중에 남은 표본 200 이상이 되도록 — 주 지지 표본의 약 1/8 – 1/4» 로 정했다(= 250 – 500 개 기대).")
    w("")
    w("## 4. 팔 × 성분 × 획득 유형 — 기준 칸 MDE · C · 표지 · 민감도 4 칸 MDE (검토 항목 3)")
    w("")
    w("- 사양이 카드의 대표 수치로 두는 것: MDE(C 함께) — 저자 결정 E3(SPEC_S4.md:10 «MDE 대표 · C 함께» · :107 «카드의 대표 수치는 MDE · C 를 함께 적는다»)."
      " k 아홉 사이에서는 «분류는 아홉 중 가장 나쁜 것(가장 낮은 C · 가장 큰 MDE)으로 하며 아홉의 중앙값을 함께 적는다»(SPEC_S4.md:106 · :153 · 문턱 :109)."
      " 그래서 «가장 나쁜 k» 열이 분류에 쓰는 값이고 «k 중앙값» 열은 함께 적는 값이다. 사양에서 «대표» 라는 말이 둘 가운데 하나에 직접 붙어 있지는 않다.")
    w("- MDE 가장 나쁜 k = `MDE_worst`(k 아홉의 최대 · validate.py:792) · k 중앙값 = `MDE_median` · C 가장 나쁜 k = `C_worst`. 범위는 그 유형의 쌍 사이 최소–최대(쌍 하나면 값 하나)."
      " 민감도 4 칸 = L0 · SNR {20, 50} × N_eff {25, 100}(plan.py:86) — 범위는 쌍 × 4 칸 사이 최소–최대. 분류 수 = 식별 / 부분 / 비식별(가장 나쁜 k).")
    w("")
    w("| 팔 | 성분 | 유형 | 쌍 | MDE 기준 · 가장 나쁜 k | MDE 기준 · k 중앙값 | C 기준 · 가장 나쁜 k | 분류 식/부/비 | 가정 의존 | 사전 의존 | MDE 민감도 4 칸 · 가장 나쁜 k | MDE 민감도 4 칸 · k 중앙값 |")
    w("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for arm in ARMS:
        for comp in COMPS:
            for t in TYPES:
                s = res["C_by_arm_component_type"][arm][comp][t]
                if not s["n_pairs"]:
                    w(f"| {arm.upper()} | {comp} | {TYPE_KO[t]} | 0 | — | — | — | — | — | — | — | — |")
                    continue
                cc = s["ref_cls_counts"]
                w(f"| {arm.upper()} | {comp} | {TYPE_KO[t]} | {s['n_pairs']} | {frange(s['ref_MDE_worst_k_range_over_pairs'])} | {frange(s['ref_MDE_median_k_range_over_pairs'])} | "
                  f"{frange(s['ref_C_worst_k_range_over_pairs'])} | {cc['식별']}/{cc['부분']}/{cc['비식별']} | {s['ref_assumption_dependent']} | {s['ref_prior_dependent']} | "
                  f"{frange(s['sens4_MDE_worst_k_range_over_pairs_and_conditions'])} | {frange(s['sens4_MDE_median_k_range_over_pairs_and_conditions'])} |")
    w("")
    w("## 5. 수준별 식별 수 (SNR 30 · N_eff 50 · 28 쌍)")
    w("")
    w("칸 값 = `cls` «식별»(가장 나쁜 k) 수 / 그 가운데 «사전 의존» 없는 수(SPEC_S4.md:131 의 뜻으로 «식별» 이라 쓸 수 있는 수) / «부분» / «비식별».")
    w("")
    w("| 팔 | 성분 | L0 (기준 칸) | L1 | L2 | L3 |")
    w("|---|---|---|---|---|---|")
    for arm in ARMS:
        for comp in COMPS:
            lv = res["identified_counts_by_level"][arm][comp]
            cells = []
            for l in ("L0", "L1", "L2", "L3"):
                q = lv[l]
                cells.append(f"{q['cls_counts']['식별']} / {q['identified_without_prior_dependent']} / {q['cls_counts']['부분']} / {q['cls_counts']['비식별']}")
            w(f"| {arm.upper()} | {comp} | " + " | ".join(cells) + " |")
    w("")
    w("## 6. 카드에 없는 것")
    w("")
    for s in res["not_computed"]:
        w(f"- {s}")
    w("")
    with open(OUT_MD, "w", encoding="utf-8") as f:
        f.write("\n".join(L))


if __name__ == "__main__":
    main()
