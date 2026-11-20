"""실측 고정 문서 §2-1 — FISO 예측 부호 (나)(저자 결정 D-w2 · 중계 ZZZw-2). 획득 기록 · 사전 범위만 쓰는 대수 — ON-Harmony 값 0.

겉보기 f_iso: logit f′_iso = logit f_iso + log(a_iso / (ν a_in + (1 − ν) a_en)) (`decomp/stage1/config.py` apparent 와 같은 식)
→ 두 장비의 차 Δ logit f′_iso = log(a_iso^B / A^B(ν)) − log(a_iso^A / A^A(ν)) 는 f_iso 와 무관하고 ν(조직 NDI)와 이완에 기댄다.
부호 = 이완 격자(카드와 같은 L1 · L4)와 ν 사전 범위(1 단계 NDI_GRID 0.30 – 0.75 · 촘촘히) 위에서 모두 + 면 «+» · 모두 − 면 «−» · 섞이면 «±».
(가) = 카드의 dy 부호(`decomp/stage1/analyze3.py` card_algebra) — 쌍마다 나란히 적는다. 부호 예측에 넣는 쌍 = ΔTR ≠ 0 이고 (나)가 L1 · L4 에서 같은 부호(«±» 아님).

    uv run --frozen python decomp/fixed/fiso_sign.py             # 카드 28 쌍 표(마크다운) · JSON 을 표준 출력으로
    uv run --frozen python decomp/fixed/fiso_sign.py --devices   # 장비 단위 21 쌍(중계 ZZZZf-2 ①)
"""
from __future__ import annotations

import itertools
import json
import math
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "stage1"))
import analyze2 as A2  # noqa: E402
import analyze3 as A3  # noqa: E402
import config as C  # noqa: E402

NU = np.round(np.arange(C.NDI_GRID[0], C.NDI_GRID[-1] + 1e-9, 0.01), 2)      # ν 사전 범위(1 단계 격자의 끝과 끝) · 0.01 간격


def grids():
    """카드와 같은 이완 격자 — L4(사전 격자 54) · L1(중심 ± 2 SD) · `analyze3.card_algebra` 와 같은 정의."""
    return {"L4": [(C.t2_triplet(t2), C.t1_triplet(t1)) for t2 in C.ETA_T2 for t1 in C.T1_SETS],
            "L1": [(t2, t1) for n, t2, t1 in C.ETA_E6 if n.startswith("L1|cv0.03|rho0.0|")]}


def dlogit_fiso(aA, aB, nu):
    """Δ logit f′_iso(B − A) — aA · aB 는 구획 가중 (in · en · iso) · nu 는 배열."""
    nu = np.asarray(nu, float)
    return np.log(aB[2] / (nu * aB[0] + (1 - nu) * aB[1])) - np.log(aA[2] / (nu * aA[0] + (1 - nu) * aA[1]))


def sign_of(v):
    v = np.asarray(v, float)
    return "+" if (v > 0).all() else "−" if (v < 0).all() else "±"


def fiso_sign(teA, trA, teB, trB, nu=NU):
    """(나) — 격자마다 부호와 |Δ logit f′_iso| 최대."""
    out = {}
    for g, grid in grids().items():
        vals = []
        for t2, t1 in grid:
            aA = C.weights(teA, t2, trA, t1) if trA else C.weights(teA, t2)
            aB = C.weights(teB, t2, trB, t1) if trB else C.weights(teB, t2)
            vals.append(dlogit_fiso(aA, aB, nu))
        v = np.concatenate(vals)
        out[g] = dict(sign=sign_of(v), max_abs=float(np.abs(v).max()))
    return out


def pairs():
    """카드와 같은 28 쌍 — `analyze3.x11` 그대로(NOT2ING 의 TR 고정 · TE · TR 은 기록 범위의 가운데)."""
    sc = dict(A2.SCANNERS)
    for k, tr in A3.NOT2ING_FIXED.items():
        sc[k] = (sc[k][0], sc[k][1], tr)
    for a, b in itertools.combinations(list(sc), 2):
        A_, B_ = sc[a], sc[b]
        yield a, b, A_[0] != B_[0], sum(A_[1]) / 2, 1000 * sum(A_[2]) / 2, sum(B_[1]) / 2, 1000 * sum(B_[2]) / 2


def table():
    rows = []
    for a, b, dsg, teA, trA, teB, trB in pairs():
        ga = A3.card_algebra(teA, trA, teB, trB)
        gn = fiso_sign(teA, trA, teB, trB)
        na, nb = gn["L1"]["sign"], gn["L4"]["sign"]
        include = abs(trB - trA) > 1e-9 and na == nb and na != "±"
        rows.append(dict(pair=f"{a}~{b}", dTE=round(teB - teA, 2), dTR_ms=round(trB - trA, 1), design_differs=dsg,
                         dy_sign_L1=ga["L1"]["dy_sign"], dy_sign_L4=ga["L4"]["dy_sign"], fiso_sign_L1=na, fiso_sign_L4=nb,
                         fiso_max_abs_dlogit_L1=gn["L1"]["max_abs"], fiso_max_abs_dlogit_L4=gn["L4"]["max_abs"],
                         agree_with_dy=(na == ga["L1"]["dy_sign"] and nb == ga["L4"]["dy_sign"]), in_sign_prediction=include))
    return rows


def device_of(unit: str) -> str:
    return "NOT2ING" if unit.startswith("NOT2ING") else unit


def device_rows(counts: dict, min_subjects: int = 3) -> list[dict]:
    """장비 단위 21 쌍(중계 ZZZZf-2 ① · [A 제안 — NOT2ING 은 피험자마다 프로토콜 하나]) — 카드 28 행을 장비 쌍으로 모은다(설계만 · 값 0).
    counts[frozenset({a, b})] = {"n": 공통 피험자 수, "protocols": {카드 단위: 수}}(세션 목록 · 메타데이터 — `decomp/remeasure/analysis.py` common_counts).
    - NOT2ING 이 없는 쌍: 카드 행 그대로.
    - NOT2ING 쌍: 예측 크기 = 공통 피험자의 프로토콜 비율로 가중한 평균 · 부호 = 공통 피험자에게 있는 프로토콜이 모두 «넣음» 이고 같은 부호일 때만
      그 부호(그래야 피험자마다 어느 프로토콜이든 예측이 같다) · NOT2ING 안 쌍(89 ~ 98)은 장비 사이가 아니라 뺀다.
    - 공통 피험자 < min_subjects 는 «판정 불가 쌍» — 부호 예측 · τ 에서 미리 뺀다(중계 ZZZZf-2 ③)."""
    groups: dict[tuple[str, str], list[dict]] = {}
    for r in table():
        a, b = r["pair"].split("~")
        da, db = device_of(a), device_of(b)
        if da == db:
            continue
        groups.setdefault((da, db), []).append(dict(r, unit_a=a, unit_b=b))
    assert len(groups) == 21 and not any((b, a) in groups for a, b in groups)
    out = []
    for (da, db), rs in groups.items():
        c = counts[frozenset((da, db))]
        n = c["n"]
        if len(rs) == 1:
            w = {rs[0]["pair"]: 1.0}
        else:
            prot = {r["pair"]: c["protocols"].get(r["unit_a"] if r["unit_a"].startswith("NOT2ING") else r["unit_b"], 0) for r in rs}
            tot = sum(prot.values())
            w = {k: (v / tot if tot else 0.0) for k, v in prot.items()}
        used = [r for r in rs if w[r["pair"]] > 0] or rs
        signs = {r["fiso_sign_L1"] for r in used}
        include = n >= min_subjects and all(r["in_sign_prediction"] for r in used) and len(signs) == 1
        why = ("판정 불가 쌍 — 공통 피험자 < %d" % min_subjects if n < min_subjects else
               "넣음" if include else "뺌 — " + " · ".join(f"{r['pair']} " + ("ΔTR 0" if abs(r["dTR_ms"]) < 1e-9 else f"{r['fiso_sign_L1']} · {r['fiso_sign_L4']}")
                                                     for r in used))
        wt = lambda k: sum(w[r["pair"]] * r[k] for r in rs) if len(rs) > 1 else rs[0][k]  # noqa: E731
        out.append(dict(pair=f"{da}~{db}", n_subjects=n, protocols=c["protocols"], weights=w, judgeable=n >= min_subjects,
                        dTE=round(wt("dTE"), 2), dTR_ms=round(wt("dTR_ms"), 1),
                        fiso_max_abs_dlogit_L1=wt("fiso_max_abs_dlogit_L1"), fiso_max_abs_dlogit_L4=wt("fiso_max_abs_dlogit_L4"),
                        rows=[r["pair"] for r in used], in_sign_prediction=include, sign=(signs.pop() if include else None), why=why))
    return out


if __name__ == "__main__" and "--devices" in sys.argv:
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "remeasure"))
    import analysis as AN  # noqa: E402
    import rconfig as RC  # noqa: E402
    rows = device_rows(AN.common_counts(RC.onh_sessions()))
    print("| 장비 쌍(21) | 공통 피험자 | NOT2ING 프로토콜 | ΔTE ms | ΔTR ms | 최대 \\|Δ logit f′_iso\\| L1 · L4 | 부호 예측 | 까닭 |")
    print("|---|---|---|---|---|---|---|---|")
    for r in rows:
        pr = " · ".join(f"{k[-2:]} {v}" for k, v in r["protocols"].items() if k != "-") or "—"
        print(f"| {r['pair']} | {r['n_subjects']} | {pr} | {r['dTE']:+.1f} | {r['dTR_ms']:+.1f} | {r['fiso_max_abs_dlogit_L1']:.4f} · "
              f"{r['fiso_max_abs_dlogit_L4']:.4f} | {r['sign'] or '—'} | {r['why']} |")
    print(json.dumps(rows, ensure_ascii=False))
    sys.exit(0)


if __name__ == "__main__":
    rows = table()
    print("| 쌍 | ΔTE ms | ΔTR ms | 설계 다름 | (가) dy 부호 L1 · L4 | (나) f_iso 부호 L1 · L4 | (나) 최대 \\|Δ logit f′_iso\\| L1 · L4 | (가) 와 같음 | 부호 예측에 넣음 |")
    print("|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        print(f"| {r['pair']} | {r['dTE']:+.1f} | {r['dTR_ms']:+.1f} | {'예' if r['design_differs'] else '아니오'} | {r['dy_sign_L1']} · {r['dy_sign_L4']} | "
              f"{r['fiso_sign_L1']} · {r['fiso_sign_L4']} | {r['fiso_max_abs_dlogit_L1']:.4f} · {r['fiso_max_abs_dlogit_L4']:.4f} | "
              f"{'예' if r['agree_with_dy'] else '아니오'} | {'예' if r['in_sign_prediction'] else '아니오'} |")
    print(json.dumps(dict(nu_range=[float(NU[0]), float(NU[-1])], n_nu=len(NU), rows=rows), ensure_ascii=False))
