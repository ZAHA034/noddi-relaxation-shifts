"""실측 고정 문서 §2-6 — Wongkornchaovalit 2025 다중 TE 시연(재현 · 기술) 판정 코드: 직선성(명제 (가)의 실증 시험) · TE 묶음 절편 안정 · ODI 의 TE 안정성.
값 계산 0 — ⑩ 뒤에 roi_<계열>.json 들에 돈다(로더 가드). 규칙은 고정 문서 §2-6 · R12 · R23 · 검토자 ZZZZd-3 · ZZZZf-1.

- 좌표 x · y(lib.coords) · ODI · TE 8 개(62 – 132) · TE 는 가운데 97 ms 로 옮겨 맞춘다 · 반복(TE62R2 · TE132R2)은 회귀에서 빼고 재검사에만.
- 직선성: 피험자 × ROI 단위 × 좌표(x · y)마다 이차 회귀 c(TE − 97)² 의 곡률 몫 |c| · 35² ≤ 2 · √2 · s_rt 면 «직선».
  s_rt = 피험자 × 좌표마다 27 단위를 모은 재검사 SD = √(Σ d² / (2 · 차의 수)) · d = 반복 − 같은 TE 원 촬영(TE62 · TE132).
  보고: 좌표별 «직선» 비율(피험자마다 · 전체) — 0.8 이상이면 «명제와 맞음».
- TE 묶음: 낮은 둘(62 · 72) · 높은 둘(122 · 132) · 전체(8) — 묶음마다 직선의 TE → 0 절편 · 묶음 간 범위 / 단일 TE 값의 TE 간 범위(기술).
- ODI: 단위마다 8 TE 범위 · 반복 차를 나란히(회귀하지 않는다).
- 민감도: 따로 정합한 계열(정렬 r < 0.90 · DONE.json «따로_정합»)을 뺀다.

    uv run --frozen python decomp/remeasure/judge_wong.py [작업 폴더]      # ⑩ 뒤 — JSON
"""
from __future__ import annotations

import json
import math
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import judge_tong as JT  # noqa: E402 — unit_coords(같은 좌표)
import rconfig as R  # noqa: E402

TES = (62, 72, 82, 92, 102, 112, 122, 132)
TE_C, SPAN = 97.0, 35.0
LIN_RATIO = 0.8
BUNDLES = {"낮은 둘": (62, 72), "높은 둘": (122, 132), "전체": TES}
REPEATS = {"TE62R2": "TE62", "TE132R2": "TE132"}


def fit_unit(te: np.ndarray, v: np.ndarray) -> dict:
    """한 단위 · 한 좌표 — 일차(기울기 · 절편 at 97) · 이차 곡률 c · 곡률 몫 |c| · 35²(점이 셋 이상일 때)."""
    t = te - TE_C
    b1 = np.polyfit(t, v, 1)
    out = dict(slope=float(b1[0]), level97=float(b1[1]), n=len(te))
    if len(te) >= 3:
        c = np.polyfit(t, v, 2)[0]
        out.update(c=float(c), curvature=float(abs(c) * SPAN ** 2))
    return out


def intercept0(te, v) -> float:
    """TE → 0 절편(직선) — 두 점이면 두 점을 잇는 직선."""
    return float(np.polyfit(np.asarray(te, float), np.asarray(v, float), 1)[1])


def s_rt(series: dict[str, dict[str, np.ndarray]], ci: int) -> tuple[float | None, int]:
    """피험자 하나 · 좌표 하나의 재검사 SD — 27 단위의 (반복 − 원) 차를 모아 √(Σ d² / (2 · 차의 수))."""
    d = [series[r][u][ci] - series[o][u][ci] for r, o in REPEATS.items() if r in series and o in series
         for u in set(series[r]) & set(series[o])]
    return (math.sqrt(sum(x * x for x in d) / (2 * len(d))), len(d)) if d else (None, 0)


def judge(data: dict[str, dict[str, dict[str, np.ndarray]]], drop: dict[str, set] | None = None) -> dict:
    """data = {피험자: {계열(TE62 … TE132R2): {단위: (x, y, ODI)}}} · drop = {피험자: 뺄 계열}(민감도)."""
    res = dict(by_subject={}, overall={})
    lin_all = {"x": [], "y": []}
    for sb, ser in sorted(data.items()):
        ser = {k: v for k, v in ser.items() if k not in (drop or {}).get(sb, set())}
        units = sorted(set.intersection(*[set(ser[f"TE{t}"]) for t in TES if f"TE{t}" in ser])) if any(f"TE{t}" in ser for t in TES) else []
        out = {}
        for ci, co in enumerate(("x", "y")):
            srt, nd = s_rt(ser, ci)
            lin, ratios = [], []
            for u in units:
                tt = np.array([t for t in TES if f"TE{t}" in ser and u in ser[f"TE{t}"]], float)
                vv = np.array([ser[f"TE{int(t)}"][u][ci] for t in tt])
                f = fit_unit(tt, vv)
                if srt is not None and "curvature" in f:
                    lin.append(f["curvature"] <= 2 * math.sqrt(2) * srt)
                ints = [intercept0([t for t in b if f"TE{t}" in ser], [ser[f"TE{t}"][u][ci] for t in b if f"TE{t}" in ser])
                        for b in BUNDLES.values() if sum(f"TE{t}" in ser for t in b) >= 2]
                span = vv.max() - vv.min()
                if len(ints) >= 2 and span > 0:
                    ratios.append((max(ints) - min(ints)) / span)
            lin_all[co] += lin
            out[co] = dict(s_rt=srt, n_diffs=nd, n_units=len(lin), linear_ratio=(float(np.mean(lin)) if lin else None),
                           bundle_range_over_te_range=dict(median=float(np.median(ratios)) if ratios else None,
                                                           q25=float(np.quantile(ratios, 0.25)) if ratios else None,
                                                           q75=float(np.quantile(ratios, 0.75)) if ratios else None, n=len(ratios)))
        odi_rng = [float(np.ptp([ser[f"TE{t}"][u][2] for t in TES if f"TE{t}" in ser])) for u in units]
        odi_rep = [float(ser[r][u][2] - ser[o][u][2]) for r, o in REPEATS.items() if r in ser and o in ser for u in set(ser[r]) & set(ser[o])]
        out["ODI"] = dict(te_range_median=float(np.median(odi_rng)) if odi_rng else None,
                          repeat_abs_diff_median=float(np.median(np.abs(odi_rep))) if odi_rep else None)
        res["by_subject"][sb] = out
    for co, v in lin_all.items():
        r = float(np.mean(v)) if v else None
        res["overall"][co] = dict(linear_ratio=r, n=len(v), verdict=(None if r is None else "명제와 맞음" if r >= LIN_RATIO else "맞지 않음"))
    return res


def load(work: pathlib.Path = R.WORK) -> tuple[dict, dict, list[str]]:
    """⑩ 뒤 — WORK/wong/<sub>/roi_<계열>.json(따로 정합한 계열은 own_<계열>/ 아래) → (data, 따로 정합한 계열, 뺀 피험자).
    끝 표지(DONE.json — «따로_정합» 이 여기 있다)가 없는 피험자는 쓰지 않는다(⑪ 사전 점검 — 끊긴 처리에서 따로 정합 집합이 빈 것으로 읽히지 않게)."""
    import guard as G  # noqa: E402
    G.require()
    data, own, missing = {}, {}, []
    for sb in R.WONG_SUBJECTS:
        wd = work / "wong" / sb
        if not (wd / "DONE.json").exists():
            missing.append(sb)
            continue
        own[sb] = set(json.loads((wd / "DONE.json").read_text(encoding="utf-8")).get("따로_정합", []))
        for s in R.WONG_SERIES:
            p = (wd / f"own_{s}" if s in own[sb] else wd) / f"roi_{s}.json"
            if p.exists():
                data.setdefault(sb, {})[s] = JT.unit_coords(json.loads(p.read_text(encoding="utf-8")))
    return data, own, missing


if __name__ == "__main__":
    d, own, miss = load(pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else R.WORK)
    print(json.dumps(dict(main=judge(d), without_own_registration=judge(d, own), own_registration=own, missing_subjects=miss),
                     ensure_ascii=False, default=list))
