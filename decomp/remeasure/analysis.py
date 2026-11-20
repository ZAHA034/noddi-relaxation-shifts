"""새 실측(⑨-a) — ROI 값 → 척도 · 쌍 요약 · 재검사 SD · 세션 SNR(중계 ZZZZf-2 ② ③ · ZZZZf-3 · ZZZZf-6). 값 계산 0 — ⑩ 뒤에 roi.json 들에 돈다.

척도(ZZZZf-2 ②): FISO 부호 · τ 의 관측 = logit f_iso — ROI 평균(roi.json «주») 뒤에 logit 을 취한다. 0 · 1 은 S4 와 같은 자름
  (`decomp/stage1/lib.py:9` CLIP = 1e-4 · `:17–19` np.clip(…, CLIP, 1 − CLIP) — S4 는 `decomp/s4/simulate.py:149` 에서 lib.coords 를 부른다).
  Tong · Wongkornchaovalit 의 x · y 도 lib.coords 그대로. 원척도 ΔFISO 는 민감도. ODI 는 원척도(coords 의 셋째 좌표와 같다).
쌍 요약(ZZZZf-2 ③): 두 장비를 모두 찍은 피험자마다 ROI 단위별 B − A → 단위마다 피험자 중앙값 → 27 단위 중앙값(부호 · τ 의 관측) · 단위별 표는 기술.
  공통 피험자 < rconfig.PAIR_MIN_SUBJECTS(3)인 쌍은 «판정 불가 쌍»(세션 목록에서 미리 센다). 실패 세션은 뺀다 · 값이 없는 ROI 단위(유효 화소 0)는
  그 피험자에서 결측으로 두고 남은 것으로 낸다.
"""
from __future__ import annotations

import itertools
import math
import pathlib
import sys
from collections import defaultdict

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "stage1"))
import lib  # noqa: E402 — S4 와 같은 좌표 · 자름
import rconfig as R  # noqa: E402

DEVICES = ("NOT1ACH", "NOT2ING", "NOT4GEP", "OXF1PRI", "OXF2PRI", "OXF3TRI", "OXF4GEP")       # decomp/fixed/delta_tau.py 와 같은 순서
CARD_UNITS = ("NOT1ACH", "NOT2ING-89", "NOT2ING-98", "NOT4GEP", "OXF1PRI", "OXF2PRI", "OXF3TRI", "OXF4GEP")   # 카드 단위(analyze2.SCANNERS · NOT2ING 두 프로토콜 따로)


def logit_fiso(fwf) -> np.ndarray:
    """ROI 평균 FISO(AMICO FWF) → logit — S4 자름(lib.CLIP). lib.coords 의 y + ln(NDI) 와 같다(정답 시험)."""
    f = np.clip(np.asarray(fwf, float), lib.CLIP, 1 - lib.CLIP)
    return np.log(f / (1 - f))


def unit_values(roi_res: dict, feature: str, scale: str = "raw", table: str = "주") -> dict[str, float]:
    """roi.json 의 27 단위 표 → {단위: 값} · scale «logit» 은 FWF(FISO)에만 · 값 없는 단위(None · 유효 화소 0)는 뺀다."""
    if scale == "logit" and feature != "FWF":
        raise ValueError("logit 척도는 FISO(FWF)에만")
    out = {}
    for u, row in roi_res[table].items():
        v = row.get(feature)
        if v is None or not math.isfinite(v):
            continue
        out[u] = float(logit_fiso(v)) if scale == "logit" else float(v)
    return out


def pair_summary(A: dict[str, dict[str, float]], B: dict[str, dict[str, float]], units, min_subjects: int = R.PAIR_MIN_SUBJECTS) -> dict:
    """A · B = {피험자: {단위: 값}}(그 장비의 첫 세션 · 실패 세션은 넣지 않는다) → 쌍 요약(B − A).
    단위마다 공통 피험자의 B − A 중앙값 → 값이 있는 단위들의 중앙값."""
    subs = sorted(set(A) & set(B))
    by_unit, n_by_unit = {}, {}
    for u in units:
        d = [B[s][u] - A[s][u] for s in subs if u in A[s] and u in B[s]]
        by_unit[u] = float(np.median(d)) if d else None
        n_by_unit[u] = len(d)
    vals = [v for v in by_unit.values() if v is not None]
    return dict(n_subjects=len(subs), judgeable=len(subs) >= min_subjects, value=float(np.median(vals)) if vals else None,
                n_units=len(vals), by_unit=by_unit, n_by_unit=n_by_unit)


def common_counts(sessions: list[tuple[str, str]]) -> dict[frozenset, dict]:
    """세션 목록(메타데이터) → 장비 쌍마다 공통 피험자 수와 NOT2ING 프로토콜별 수 — 값 0."""
    by = defaultdict(set)
    for sub, ses in sessions:
        by[R.device_of(ses)].add(sub)
    out = {}
    for a, b in itertools.combinations(DEVICES, 2):
        c = by[a] & by[b]
        prot = defaultdict(int)
        for s in c:
            prot[R.not2ing_protocol(s) if "NOT2ING" in (a, b) else "-"] += 1
        out[frozenset((a, b))] = dict(n=len(c), subjects=sorted(c), protocols=dict(sorted(prot.items())), judgeable=len(c) >= R.PAIR_MIN_SUBJECTS)
    return out


def card_common_counts(sessions: list[tuple[str, str]]) -> dict[frozenset, int]:
    """카드 28 쌍마다 공통 피험자 수(V4 칸은 카드 쌍 × 성분 — `decomp/SPEC_S4.md:228` A-1 6) — 0 이면 그 쌍의 칸은 «자료 없음»(a11). 세션 목록 · NOT2ING 프로토콜(메타데이터)만 · 값 0."""
    by = defaultdict(set)
    for sub, ses in sessions:
        d = R.device_of(ses)
        by[R.not2ing_protocol(sub) if d == "NOT2ING" else d].add(sub)
    return {frozenset((a, b)): len(by[a] & by[b]) for a, b in itertools.combinations(CARD_UNITS, 2)}


def s_ref(values: dict[str, list[dict[str, float]]], units) -> tuple[dict[str, float | None], int]:
    """재검사 SD — `prereg/scripts/lockd/compute_s_ref.py` 와 같은 식: 피험자마다 세션 사이 분산(ddof 1) → 피험자 평균 → 제곱근.
    values = {피험자: [세션마다 {단위: 값}]} · 자유도 = Σ(세션 − 1)(옛 판 2 × 5 = 10). 값이 없는 세션은 그 단위에서 뺀다."""
    out = {}
    for u in units:
        var = []
        for sess in values.values():
            v = [d[u] for d in sess if u in d]
            if len(v) >= 2:
                var.append(float(np.var(v, ddof=1)))
        out[u] = math.sqrt(sum(var) / len(var)) if var else None
    df = sum(len(s) - 1 for s in values.values())
    return out, df


def sd_multipliers(df: int) -> tuple[float, float]:
    """SD 의 95 % 구간 배수(자유도만) — √(df/χ²_{0.975, df}) · √(df/χ²_{0.025, df})(`compute_s_ref.py:11` 과 같은 식 · df 10 → 0.6987 · 1.7549).
    s_ref_F2 의 민감도 배수 — 실패 세션이 있으면 남은 세션의 자유도로 다시 낸다."""
    from scipy.stats import chi2
    return math.sqrt(df / chi2.ppf(0.975, df)), math.sqrt(df / chi2.ppf(0.025, df))


def snr(b0: np.ndarray, mask: np.ndarray, fa: np.ndarray, fa_min: float = R.SNR_FA_MIN, train=R.SNR_TRAIN) -> dict:
    """세션 SNR [A 제안 · 검토자 확인 — 중계 ZZZZf-6 1] — 백질 화소 W = 뇌 마스크 ∩ dtifit FA ≥ fa_min 에서
    S̄ = 화소마다 b0 부피들의 평균 · σ² = 화소마다 b0 부피들의 표본 분산(ddof 1) → SNR = 평균_W(S̄) / √(평균_W(σ²)).
    S4 의 정의(각 TE 에서 자기 b0 로 정규화한 신호에 σ = 1/SNR · `decomp/SPEC_S4.md:76`)와 같은 양 — 학습 범위(U_log[15, 60]) 밖이면 표시."""
    W = (mask > 0) & (fa >= fa_min)
    x = np.asarray(b0, float)[W]
    k = x.shape[-1] if x.ndim == 2 else 0
    if k < 2 or len(x) == 0:
        return dict(snr=None, n_b0=k, n_vox=int(W.sum()), in_train=None)
    val = float(x.mean(1).mean() / math.sqrt(x.var(1, ddof=1).mean()))
    return dict(snr=val, n_b0=k, n_vox=int(W.sum()), in_train=bool(train[0] <= val <= train[1]))


def retest_level(by_unit: dict[str, float | None], s: dict[str, float | None], c: float = 1.96, ratio: float = 0.90) -> dict:
    """ODI (b) 의 쌍 규칙(고정 문서 §2-2 · §3 규칙 3) — 단위마다 |Δ| ≤ c · √2 · s 면 «재검사 수준» · 값과 s 가 모두 있는 단위 중 비율 ≥ 0.90 이면
    그 쌍은 «재검사 수준»(27 단위가 모두 있으면 25 단위 이상 — 권고 ZZZZh-7 7). Δ = 쌍 요약의 단위 값(피험자 중앙값 — 계통 차의 여유이지 잡음 검정이 아니다)."""
    use = [u for u, v in by_unit.items() if v is not None and s.get(u) is not None]
    hit = [u for u in use if abs(by_unit[u]) <= c * math.sqrt(2) * s[u]]
    need = math.ceil(ratio * len(use) - 1e-9)
    return dict(n_units=len(use), n_retest=len(hit), need=need, retest_level=bool(use) and len(hit) >= need)


def clipped_units(roi_res: dict, table: str = "주") -> list[str]:
    """logit 자름 보고(권고 ZZZZh-7 5) — ROI 평균 FISO(FWF)가 S4 자름(lib.CLIP = 1e-4) 이하 · 1 − 1e-4 이상이라 logit ±9.21 로 들어간 단위."""
    return sorted(u for u, row in roi_res[table].items()
                  if row.get("FWF") is not None and math.isfinite(row["FWF"]) and not (lib.CLIP < row["FWF"] < 1 - lib.CLIP))
