"""실측 고정 문서 §2-3 — Tong 2020(맹검 · 음성 대조) 판정 코드: B-iii 재검사 SD · T-B 크기 · T-A 방향(회전 순열) · T-C 2 × 2.
값 계산 0 — ⑩ 뒤에 roi.json 들에 돈다(실자료 경로면 로더 가드). 규칙은 고정 문서 §2-3 «정의» · R10 · a6 · a10 · 중계 ZZZZh-2 · ZZZZh-4.

- 좌표: x = logit ν · y = logit f_iso − ln ν(`decomp/stage1/lib.py` coords · 카드와 같다) · ODI 원척도. 주 = roi.json «주»(b ≤ 2000).
- B-iii: 좌표 · ROI 단위마다 피험자 안 분산(c10r1 – r3 · ddof 1) → 피험자 평균 → 제곱근 · 자유도 ν = Σ(반복 − 1).
- T-B: 단위 = 센터 쌍(45) × 피험자(3) × ROI(27) · Δ = 두 센터 값의 차(c10 = c10r1) · «재검사 수준» = |Δ| ≤ c(ν) · √2 · s ·
  c(ν) = t₀.₉₇₅(ν)(새 측정 하나에 대한 95 % 예측 구간 — 검토자 제안 ZZZZh-2) · 민감도 B-i · B-ii 는 c = 1.96(s 를 알려진 값으로 본다).
  비율 = Δ 와 s 가 모두 있는 단위 중 · s 가 없는 ROI 의 단위는 뺀다(a10) · 문턱 0.90 · 판정 불가 = 분모 < 전체의 절반 또는 s 없는 ROI > 절반.
- T-A: 센터 효과 e_c = (센터 값 − 그 피험자 · ROI 의 센터 평균)의 피험자 · ROI 중앙값(x · y) · 쌍 벡터 d = e_B − e_A ·
  «방향이 정해짐» = |d| ≥ 1.96 · √2 · s(s = B-iii 의 x · y 단위 중앙값 SD 중 큰 값 — a6 · 크기 여유(margin)이지 잡음 검정이 아니다) ·
  정해지지 않은 쌍이 절반을 넘으면 판정 불가 · 통계량 = 정해진 쌍 중 D 안 비율 · 회전 순열 3,600 각(0.05° · 항등 포함 — ZZZZh-4) ·
  p = #{θ : 비율(θ) ≥ 비율(0)} / 3,600(한쪽) · p < 0.05 면 «정렬됨».

    uv run --frozen python decomp/remeasure/judge_tong.py [작업 폴더]      # ⑩ 뒤 — JSON
"""
from __future__ import annotations

import itertools
import json
import math
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "stage1"))
import lib  # noqa: E402
import rconfig as R  # noqa: E402

COORDS = ("x", "y", "ODI")
CENTERS = tuple(f"c{i:02d}" for i in range(1, 11))                         # c10 = c10r1(센터 사이 쌍)
REPEATS = ("ses-c10r1", "ses-c10r2", "ses-c10r3")
D_WEDGES = ((48.5, 72.1), (85.3, 94.3))                                     # 선 방향(도) — TE 쐐기 ∪ TR 부채꼴(2 단계 X1 · 68e36a1)
P_RANDOM = sum(b - a for a, b in D_WEDGES) / 180.0                          # 0.181
N_ANGLES = 3600                                                             # [0°, 180°) · 0.05° 간격 · 항등 포함(ZZZZh-4)
TB_RATIO = 0.90
C_KNOWN = 1.96                                                              # s 를 알려진 값으로 볼 때(B-i · B-ii · T-A 크기 여유)
ALPHA = 0.05


def session_center(ses: str) -> str:
    return ses.replace("ses-", "")[:3]


def unit_coords(roi_res: dict, table: str = "주") -> dict[str, np.ndarray]:
    """roi.json → {ROI 단위: (x, y, ODI)} — NDI · FWF · ODI 가 모두 있는 단위만."""
    out = {}
    for u, row in roi_res[table].items():
        v = [row.get(k) for k in ("NDI", "FWF", "ODI")]
        if any(a is None or not math.isfinite(a) for a in v):
            continue
        xy, _ = lib.coords(np.array([v[0]]), np.array([v[1]]), np.array([v[2]]))
        out[u] = xy[0]
    return out


def retest_sd(rep: dict[str, list[dict[str, np.ndarray]]], units) -> dict[str, dict[str, tuple[float | None, int]]]:
    """B-iii — rep = {피험자: [반복마다 {단위: (x, y, ODI)}]} → {좌표: {단위: (s, ν)}} · 반복이 둘 미만인 피험자는 그 단위에서 뺀다."""
    out = {c: {} for c in COORDS}
    for ci, c in enumerate(COORDS):
        for u in units:
            var, df = [], 0
            for sess in rep.values():
                v = [d[u][ci] for d in sess if u in d]
                if len(v) >= 2:
                    var.append(float(np.var(v, ddof=1)))
                    df += len(v) - 1
            out[c][u] = (math.sqrt(sum(var) / len(var)), df) if var else (None, 0)
    return out


def c_of(df: int) -> float:
    """T-B 의 c = t₀.₉₇₅(ν)(ZZZZh-2) — ν = 6 이면 2.4469 · ν 가 줄면 커진다 · ν = 0 이면 정의 없음(nan)."""
    if df < 1:
        return float("nan")
    from scipy.stats import t
    return float(t.ppf(0.975, df))


def tb(values: dict, s: dict, centers=CENTERS, c_fixed: float | None = None) -> dict:
    """T-B — values[센터][피험자][단위] = (x, y, ODI) · s = retest_sd 꼴 {좌표: {단위: (s, ν)}} · c_fixed 를 주면 그 c(민감도 B-i · B-ii).
    좌표마다: 전체 단위 수 · 분모(Δ 와 s 가 모두 있는 단위) · 재검사 수준 수 · 비율 · s 없는 ROI 수와 그 단위 수 · 판정."""
    pairs = list(itertools.combinations(centers, 2))
    subs = sorted({sb for c in centers for sb in values.get(c, {})})
    units = sorted({u for c in centers for sb in values.get(c, {}) for u in values[c][sb]})
    n_total = len(pairs) * len(R.TONG_SUBJECTS) * 27
    out = {}
    for ci, co in enumerate(COORDS):
        no_s = sorted(u for u in units if s[co].get(u, (None, 0))[0] is None)
        den = hit = dropped = 0
        for a, b in pairs:
            for sb in subs:
                va, vb = values.get(a, {}).get(sb, {}), values.get(b, {}).get(sb, {})
                for u in set(va) & set(vb):
                    su, nu = s[co].get(u, (None, 0))
                    if su is None:
                        dropped += 1
                        continue
                    c = c_fixed if c_fixed is not None else c_of(nu)
                    den += 1
                    hit += abs(vb[u][ci] - va[u][ci]) <= c * math.sqrt(2) * su
        judgeable = den >= n_total / 2 and len(no_s) <= 13
        ratio = hit / den if den else None
        out[co] = dict(n_total=n_total, n_den=den, n_retest=int(hit), ratio=ratio, roi_without_s=no_s, units_dropped=dropped,
                       judgeable=judgeable, verdict=("판정 불가" if not judgeable else "재검사 수준" if ratio >= TB_RATIO else "넘음"))
    return out


def center_effects(values: dict, centers=CENTERS) -> dict[str, np.ndarray]:
    """센터 효과 e_c(x · y) — 피험자 · 단위마다 센터 값에서 그 피험자 · 단위의 센터 평균을 뺀 뒤 피험자 · 단위 중앙값."""
    keys = sorted({(sb, u) for c in centers for sb in values.get(c, {}) for u in values[c][sb]})
    out = {}
    for c in centers:
        dev = []
        for sb, u in keys:
            have = [values[cc][sb][u][:2] for cc in centers if u in values.get(cc, {}).get(sb, {})]
            if u in values.get(c, {}).get(sb, {}) and len(have) >= 2:
                dev.append(values[c][sb][u][:2] - np.mean(have, 0))
        out[c] = np.median(np.array(dev), 0) if dev else np.full(2, np.nan)
    return out


def center_inputs(values: dict, centers=CENTERS) -> dict:
    """센터 효과의 재료 보고(검토자 ZZZZi-3 · a14 덧붙임) — (피험자 · 단위) 수 · 일부 센터가 빠져 있는 센터로 평균한 수 · 센터가 둘 미만이라 뺀 수."""
    keys = sorted({(sb, u) for c in centers for sb in values.get(c, {}) for u in values[c][sb]})
    n_have = [sum(u in values.get(c, {}).get(sb, {}) for c in centers) for sb, u in keys]
    return dict(n_subject_units=len(keys), n_partial=sum(2 <= n < len(centers) for n in n_have), n_dropped=sum(n < 2 for n in n_have))


def in_d(phi_deg: np.ndarray, wedges=D_WEDGES) -> np.ndarray:
    phi = np.mod(phi_deg, 180.0)
    return np.any([(phi >= a) & (phi <= b) for a, b in wedges], axis=0)


def line_angle(d: np.ndarray) -> np.ndarray:
    """벡터 (n, 2) → 선 방향(도 · [0, 180))."""
    return np.mod(np.degrees(np.arctan2(d[:, 1], d[:, 0])), 180.0)


def ta(e: dict[str, np.ndarray], s_xy: tuple[float, float], centers=CENTERS, wedges=D_WEDGES, n_angles: int = N_ANGLES) -> dict:
    """T-A — e = center_effects · s_xy = (B-iii 의 x 단위 중앙값 SD, y 단위 중앙값 SD)(a6).
    «방향이 정해짐» = |d| ≥ 1.96 · √2 · max(s_x, s_y) · 정해지지 않은 쌍이 절반을 넘으면 판정 불가 · 회전 순열 p(ZZZZh-4)."""
    pairs = list(itertools.combinations(centers, 2))
    d = np.array([e[b] - e[a] for a, b in pairs])
    margin = C_KNOWN * math.sqrt(2) * max(s_xy)
    det = np.isfinite(d).all(1) & (np.linalg.norm(np.nan_to_num(d), axis=1) >= margin)
    n_und = int(len(pairs) - det.sum())
    res = dict(n_pairs=len(pairs), n_determined=int(det.sum()), n_undetermined=n_und, margin=margin,
               undetermined_limit=math.floor(len(pairs) / 2) + 1, p_random=sum(b - a for a, b in wedges) / 180.0)
    if n_und > len(pairs) / 2 or det.sum() == 0:
        return res | dict(judgeable=False, verdict="판정 불가", ratio=None, p=None)
    n_in, p = rotation_p(line_angle(d[det]), wedges, n_angles)                    # 회전해도 길이는 그대로 — 정해진 쌍 집합이 같다
    return res | dict(judgeable=True, ratio=float(n_in / det.sum()), n_in_d=n_in, p=p, n_angles=n_angles,
                      verdict="정렬됨" if p < ALPHA else "정렬 안 됨")


def rotation_p(phi: np.ndarray, wedges=D_WEDGES, n_angles: int = N_ANGLES) -> tuple[int, float]:
    """회전 순열(ZZZZh-4) — 선 방향 phi(도)를 [0°, 180°) 의 0.05° 격자 3,600 각(항등 포함)으로 함께 돌려 D 안 개수를 센다 ·
    p = #{θ : 개수(θ) ≥ 개수(0)} / 3,600(한쪽 · 같은 값은 ≥ 로 센다) → (항등의 D 안 개수, p)."""
    theta = np.arange(n_angles) * (180.0 / n_angles)
    counts = in_d(np.asarray(phi, float)[None, :] + theta[:, None], wedges).sum(1)
    return int(counts[0]), float((counts >= counts[0]).mean())


def standardized_wedges(s_xy: tuple[float, float], wedges=D_WEDGES) -> tuple[tuple[float, float], ...]:
    """좌표를 (x / s_x, y / s_y) 로 옮길 때 선 방향 φ 의 상(권고 ZZZZh-7 3 · 기술 민감도) — 단조라 끝점만 옮긴다."""
    sx, sy = s_xy
    f = lambda a: float(np.mod(np.degrees(np.arctan2(math.sin(math.radians(a)) / sy, math.cos(math.radians(a)) / sx)), 180.0))  # noqa: E731
    return tuple((f(a), f(b)) for a, b in wedges)


def ta_standardized(e: dict[str, np.ndarray], s_xy: tuple[float, float], centers=CENTERS) -> dict:
    """T-A 기술 민감도 — 센터 효과를 B-iii 단위 중앙값 SD 로 나눈 좌표에서(크기 여유는 1.96 · √2 — 나눈 좌표의 SD 1) · 판정은 원좌표."""
    sx, sy = s_xy
    es = {c: v / np.array([sx, sy]) for c, v in e.items()}
    return ta(es, (1.0, 1.0), centers, standardized_wedges(s_xy))


def tc(ta_res: dict, tb_res: dict) -> str:
    """T-C 2 × 2 — T-A(정렬됨 · 정렬 안 됨) × T-B «재검사 수준» = x · y 모두 0.90 이상(ODI 는 따로 기술) · 하나라도 판정 불가면 판정 불가."""
    if not ta_res["judgeable"] or not (tb_res["x"]["judgeable"] and tb_res["y"]["judgeable"]):
        return "판정 불가"
    b = "재검사 수준" if tb_res["x"]["ratio"] >= TB_RATIO and tb_res["y"]["ratio"] >= TB_RATIO else "넘음"
    return f"{ta_res['verdict']} · {b}"


def judge(values: dict, rep: dict, centers=CENTERS) -> dict:
    """values = {센터: {피험자: {단위: (x, y, ODI)}}}(c10 = c10r1) · rep = {피험자: [c10r1 · r2 · r3 의 {단위: 좌표}]} → 판정 전부."""
    units = sorted({u for c in centers for sb in values.get(c, {}) for u in values[c][sb]})
    s = retest_sd(rep, units)
    med = lambda co: float(np.median([v[0] for v in s[co].values() if v[0] is not None])) if any(v[0] is not None for v in s[co].values()) else float("nan")  # noqa: E731
    s_xy = (med("x"), med("y"))
    e = center_effects(values, centers)
    ta_res, tb_res = ta(e, s_xy, centers), tb(values, s, centers)
    return dict(centers=list(centers), s_xy_unit_median=s_xy, df_by_coord={co: sorted({v[1] for v in s[co].values()}) for co in COORDS},
                center_effect_inputs=center_inputs(values, centers),
                T_A=ta_res, T_A_standardized_descriptive=ta_standardized(e, s_xy, centers), T_B=tb_res, T_C=tc(ta_res, tb_res))


def load(work: pathlib.Path = R.WORK, table: str = "주", fname: str = "roi.json") -> tuple[dict, dict, list[str]]:
    """⑩ 뒤 — WORK/tong/<sub>_<ses>/roi.json → (values, rep, 뺀 세션). 실자료 경로면 로더 가드.
    끝 표지(DONE.json)가 없는 세션은 쓰지 않고 뺀 세션으로 돌려준다(⑪ 사전 점검 — 끊긴 처리가 조용히 섞이지 않게)."""
    import guard as G  # noqa: E402
    G.require()
    values, rep, missing = {}, {}, []
    for sb in R.TONG_SUBJECTS:
        for ses in R.TONG_SESSIONS:
            wd = work / "tong" / f"{sb}_{ses}"
            p = wd / fname
            if not (wd / "DONE.json").exists() or not p.exists():
                missing.append(f"{sb}_{ses}")
                continue
            uc = unit_coords(json.loads(p.read_text(encoding="utf-8")), table)
            if ses in REPEATS:
                rep.setdefault(sb, []).append(uc)
            if ses.endswith("r1"):
                values.setdefault(session_center(ses), {})[sb] = uc
    return values, rep, missing


if __name__ == "__main__":
    work = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else R.WORK
    v, rp, miss = load(work)
    res = dict(main=judge(v, rp), without_c09=judge(v, rp, tuple(c for c in CENTERS if c != "c09")), missing_sessions=miss)
    vs, rps = load(work, fname="roi_allshell.json")
    res["all_shells_descriptive"] = judge(vs, rps)
    print(json.dumps(res, ensure_ascii=False, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))
