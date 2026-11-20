"""실측 고정 문서 §2-5 — V4(travelling heads 실측 보정 시험 · 확인용) 판정 코드: 동결 S4 판단기(가중치 85eb8b0)를 실측 관측에 그대로 적용한다.
규칙 = `decomp/SPEC_S4.md` 부록 A(:211–219) · A-1(:220–229)을 옮긴 것이다 — 다르면 SPEC 이 정본. 학습 0 · 값 계산 0(⑩ 뒤에 돈다).

- 관측(A-1 1): 세션마다 roi.json «주» 의 F2 일곱 지도 ROI 평균 → lib.coords(C7) → 두 세션을 이어 붙인 C7(14) · C3(6)(simulate.py:149–153 과 같은 식).
  쌍의 방향은 카드 이름 «A~B» 그대로(사이트 A = 앞 장비) · 단위 = 카드 쌍 × 두 장비를 모두 찍은 피험자(첫 세션) × ROI 27 × 성분.
- 조건 a(A-1 2): TE · TR = 그 세션의 사이드카 · 설계 = 카드 기록(plan.ONH) · b0 의 b 값 · 수 = 설계 파일(s4cfg.design_b0 — 카드와 같은 기록) ·
  코일 = 기록(OXF1PRI 32 · OXF2PRI 64) · 기록 없는 장비는 32 · 64 로 각각 내고 모든 조합이 덮을 때만 «덮음»(두 장비 모두 기록이 없으면 조합 넷 ·
  아홉 k 표본도 모든 조합에서 차야 쓴다 — a15) · L0 · 주 칸 SNR 30 · N_eff 50 ·
  민감도 칸 SNR {20, 50} × N_eff {25, 100}(«식별됨» 칸 · 문턱은 주 칸 것 그대로 — A-1 2 A 읽기).
- k 아홉(A-1 3): k 마다의 조건부 사후(앙상블 다섯 · 주 사전 재가중 · 표본 1,000 — validate.draw)의 90 % 중앙 구간이 0 을 덮는가 →
  칸 덮음률 = k 마다 덮음률의 중앙값(최소 · 최대 함께 · 문턱을 사이에 두면 «가정 의존») · 아홉의 균등 혼합 덮음률은 참고.
- 오지정(A-1 4): 점수 −log[(1/9) Σ_k q(x_c3 | a, k)](validate.q_scores 와 같은 식) > 그 층의 문턱(miss.json · 층 = 쌍의 획득 유형 | L0)이면
  그 단위를 덮음률에서 뺀다 · 칸마다 초과율 · 기대 초과율 분포(시험 자료 같은 층 · 독립 묶음 수 · 2,000 번)의 95 백분위를 넘으면
  «분포 밖 — 간극 검출» · 민감도 = 뺀 단위를 넣은 덮음률 · 코일 조합이 여럿이면 조합마다 점수를 내고 모든 조합이 넘을 때만 뺀다(최솟값 —
  검토자 판정 ZZZZj-3 · 참 코일은 후보 안에 있으므로 뺀 단위 ⊆ 참 조합으로 뺄 단위 · 칸 초과율도 같은 점수라 미상 코일 칸의 간극 검출은 하한) ·
  칸 · 코일 조합마다 문턱을 넘은 단위 수와 표본이 차지 않은 단위 수를 적는다(R26 — 수만 · 검토자 제안 ZZZZi-3).
- «식별됨» 칸(A-1 5 · SPEC_S4.md:131): 카드 분류(L0 · 주 칸 · 가장 나쁜 k)가 «식별» 이고 «사전 의존» 이 아닌 (쌍 · 성분) — C3 29 칸.
  그 밖의 칸은 «덮음 — 정보 없음»(통과로 세지 않음 · (다)).
- 문턱(A-1 7): 시험 자료 L0 · 같은 획득 유형의 경우마다 같은 절차(k 아홉 · 90 % 구간 · 오지정 문 · 동결 판단기로 다시 뽑음)로 «참 Δ 덮음» →
  독립 묶음 수(그 쌍의 두 장비를 모두 찍은 피험자 수)만큼 복원 추출해 덮음률(k 아홉 중앙값)을 2,000 번 → 5 백분위 · 못 내면 0.80.
- 판정은 칸(카드 쌍 × 성분)마다 따로(A-1 7) — 칸을 모아 «유형» 비율로 내지 않는다. 자료 없는 칸(공통 피험자 0 — a11)은 «자료 없음» 으로
  따로 적고 어떤 분모에도 넣지 않는다(검토자 제안 ZZZZh-5) · 민감도 칸에서 판정이 갈리면 «SNR/N_eff 의존».
- 불확실성((라) · A-1 6): 피험자 묶음 부트스트랩 2,000 회 · 장비 집중(실패 칸이 한 장비에 몰리는가 — 장비 순열 5,040 · [A 읽기 · 기술]).
- 기술 분해(§2-5-1 · a8): 칸마다 속성 묶음(① 잡음 · ② 획득 범위 · ③ 장비 같음/다름 · ④ 설명 안 됨)별 세션 쌍 수 · 단위 수 · 덮은 수 · 비율.

- 순서(R18 · R20 · §2-5-1 — 코드가 강제): ① snr-table(세션별 SNR 표 → sha256 을 열람 기록에) ② attr-table <SNR 표 sha256>(속성 표 → sha256 을 열람 기록에)
  ③ test-cover <팔> ④ run <팔> <속성 표 sha256> — 앞 단계의 sha256 이 맞지 않으면 멈춘다 · 결과는 파일(WORK/onharmony/v4/v4_<팔>.json)과 sha256.
  SNR 표 · 속성 표 · 결과는 한 번만 쓰고(있으면 멈춤) 쓸 때 열람 기록(rconfig.ACCESS_LOG_ONH)에 sha256 · 시각을 적는다 · V4 결과 파일이 있으면
  SNR 표 · 속성 표를 만들지 않는다(반영 확인 읽기 F3).
- 재현: 사후 표본의 torch 시드를 단계 · 칸 · 쌍마다 고정하고 출력에 적는다 · 단위 결과는 한 번만 뽑고 민감도 둘(16794 넣음 · SNR 범위 밖 뺌)은 행만 고른다 ·
  문턱 · 기대 초과율은 (획득 유형 · 성분 · 묶음 수)만의 함수라 변형이 달라도 같은 입력이면 같은 값(반영 확인 읽기 F2).
- 아홉 k 모두 표본이 찬 단위만 쓴다(시험 자료 문턱과 같은 규칙) · 묶음(피험자) 3 명 미만인 칸은 «판정 불가» · 단위 수로 뽑은 문턱은 참고로 함께(A-1 7).

    CUDA_VISIBLE_DEVICES= uv run --frozen python decomp/remeasure/judge_v4.py snr-table
    CUDA_VISIBLE_DEVICES= uv run --frozen python decomp/remeasure/judge_v4.py attr-table <SNR 표 sha256>
    CUDA_VISIBLE_DEVICES= uv run --frozen python decomp/remeasure/judge_v4.py test-cover <c3|c7>
    CUDA_VISIBLE_DEVICES= uv run --frozen python decomp/remeasure/judge_v4.py run <c3|c7> <속성 표 sha256>
"""
from __future__ import annotations

import hashlib
import itertools
import json
import math
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "stage1"))
sys.path.insert(0, str(HERE.parents[0] / "s4"))
import lib  # noqa: E402
import plan as P  # noqa: E402 — 카드 기록(ONH)
import rconfig as R  # noqa: E402
import s4cfg as S  # noqa: E402

RESULTS = HERE.parents[0] / "results" / "s4"
COMPS = ("NDI", "FISO", "ODI")
MAIN_CELL = ("L0", 30.0, 50.0)
SENS_CELLS = tuple(("L0", s, n) for s in (20.0, 50.0) for n in (25.0, 100.0))
N_KEEP = 1000
N_BOOT = 2000
FALLBACK_THRESHOLD = 0.80
CARD_ORDER = tuple(P.ONH)                                                   # 카드 쌍의 방향(앞 장비 = 사이트 A)
NO_DATA = ("NOT2ING-89~NOT2ING-98", "NOT2ING-89~NOT4GEP", "NOT2ING-89~OXF4GEP", "OXF3TRI~NOT4GEP", "OXF3TRI~OXF4GEP")   # a11


# ── 동결 산출물(ZZZZh-3) ────────────────────────────────────────────────────────
def frozen_expected() -> dict[str, str]:
    """동결 산출물의 기대 sha256 — 저장소 기록(train_freeze.json · sim_hashes.json — 가드 목록 안)에서 읽는다(손으로 적지 않음)."""
    tf = json.loads((RESULTS / "train_freeze.json").read_text(encoding="utf-8"))["models"]
    sh = json.loads((RESULTS / "sim_hashes.json").read_text(encoding="utf-8"))["files"]
    out = {f"models/{m}.obj.pt": tf[m]["sha256_obj"] for m in [f"{a}_{i}" for a in ("c3", "c7") for i in range(5)] + ["q_0"]}
    out.update({f: sh[f]["sha256"] for f in ("sim_test.npz", "plan_test.npz")})
    return out


def check_frozen(work: pathlib.Path = S.WORK) -> dict[str, str]:
    """쓰기 전에 동결 산출물의 sha256 을 대조한다 — 하나라도 다르면 멈춘다(⑩ 뒤 바뀐 산출물로 판정하지 않게)."""
    exp = frozen_expected()
    bad = {f: h for f, h in exp.items() if hashlib.sha256((work / f).read_bytes()).hexdigest() != h}
    if bad:
        raise RuntimeError(f"⛔ 동결 산출물 sha256 이 다르다 — {sorted(bad)}")
    return exp


# ── 칸 · 쌍 ─────────────────────────────────────────────────────────────────────
def card_pairs() -> list[str]:
    return [f"{a}~{b}" for a, b in itertools.combinations(CARD_ORDER, 2)]


def acq_type(pair: str) -> str:
    """쌍의 획득 유형(A-1 5 · §7 층) — 카드 기록의 ΔTE · ΔTR(범위 가운데) · 설계 다름 → match · te · tr · tetr · tetr_design."""
    a, b = (P.ONH[u] for u in pair.split("~"))
    dte = abs(sum(a[1]) - sum(b[1])) > 1e-9
    dtr = abs(sum(a[2]) - sum(b[2])) > 1e-9
    if dte and dtr:
        return "tetr_design" if a[0] != b[0] else "tetr"
    return "te" if dte else "tr" if dtr else ("tetr_design" if a[0] != b[0] else "match")


def identified_cells(cards: dict) -> dict[str, list[str]]:
    """«식별됨» 칸 — 카드 분류 L0 · 주 칸 · 가장 나쁜 k 가 «식별» 이고 «사전 의존» 이 아님(SPEC_S4.md:227 · :131) → {성분: [쌍]}."""
    out = {c: [] for c in COMPS}
    for key, v in cards["cells"].items():
        pair, lv, snr, ne = key.split("|")
        if (lv, snr, ne) != ("L0", "SNR30", "Neff50"):
            continue
        for c in COMPS:
            if v[c]["cls"] == "식별" and "사전 의존" not in v[c].get("flags", []):
                out[c].append(pair)
    return {c: sorted(v) for c, v in out.items()}


# ── 조건 · 관측 ────────────────────────────────────────────────────────────────
def coil_options(unit: str) -> tuple[int, ...]:
    c = P.ONH[unit][3]
    return (c,) if c else S.COILS


def cond_rows(pair: str, te: np.ndarray, tr_ms: np.ndarray, coils: tuple[int, int], snr: float, neff: float) -> np.ndarray:
    """조건 a — te (n, 2) ms · tr_ms (n, 2) ms(사이드카) · 설계 · b0 = 카드 기록 · L0 · 생성기 칸은 0(k 는 validate.with_gs 로 채운다)."""
    a, b = pair.split("~")
    n = len(te)
    acq = dict(TE_A=te[:, 0], TE_B=te[:, 1], TR_A=tr_ms[:, 0], TR_B=tr_ms[:, 1],
               design_A=np.full(n, P.ONH[a][0]), design_B=np.full(n, P.ONH[b][0]),
               coil_A=np.full(n, coils[0]), coil_B=np.full(n, coils[1]), snr_A=np.full(n, snr), snr_B=np.full(n, snr), neff=np.full(n, neff))
    lv, p = S.sample_level(np.random.default_rng(0), n, np.full(n, "L0"))
    return S.cond_vector(acq, lv, p, np.full(n, S.GENERATORS[0]), np.full(n, S.SCENARIOS[0]))


def session_x(roi_res: dict, units) -> dict[str, np.ndarray]:
    """한 세션 roi.json «주» → {ROI 단위: C7 좌표 7}(일곱 지도가 모두 있는 단위만 · lib.coords)."""
    out = {}
    for u in units:
        row = roi_res["주"].get(u, {})
        v = [row.get(f) for f in R.FEATURES]
        if any(a is None or not math.isfinite(a) for a in v):
            continue
        y, _ = lib.coords(np.array([v[0]]), np.array([v[1]]), np.array([v[2]]), np.array([v[3:7]]), "C7")
        out[u] = y[0]
    return out


# ── 덮음 · 오지정 ───────────────────────────────────────────────────────────────
def wquantile(x: np.ndarray, w: np.ndarray, q: float) -> float:
    """가중 분위(validate.wquantile 과 같은 식)."""
    o = np.argsort(x)
    cw = np.cumsum(w[o]) - 0.5 * w[o]
    return float(np.interp(q, cw, x[o]))


def cover_by_k(draw_fn, X: np.ndarray, cond_list: list[np.ndarray], k_list, n_keep: int = N_KEEP, truth: np.ndarray | None = None):
    """단위마다 · k 마다 · 성분마다 90 % 중앙 구간이 참값(실측은 0)을 덮는가 — 코일 조합(cond_list)이 여럿이면 모두 덮을 때만 덮음.
    draw_fn(X, cond, n_keep, gs) → (T (b, n_keep, ≥3), ok (b,), acc (b,) 주 지지 비율) ·
    반환 cov (n, 9, 3) · ok (n, 9)(모든 조합에서 표본이 찼을 때) · cov_mix (n, 3) 아홉의 균등 혼합(참고 — 가중 ∝ 주 지지 비율 · validate.mix_pool) ·
    ok_by_combo (n, 조합 수)(그 조합에서 아홉 k 모두 표본이 찼는가 — R26 보고)."""
    n = len(X)
    truth = np.zeros((n, 3)) if truth is None else truth
    cov = np.ones((n, len(k_list), 3), bool)
    ok = np.ones((n, len(k_list)), bool)
    ok_by_combo = np.ones((n, len(cond_list)), bool)
    cov_mix = np.ones((n, 3), bool)
    for j, cond in enumerate(cond_list):
        Ts, accs = [], []
        for kk, gs in enumerate(k_list):
            T, okk, acc = draw_fn(X, cond, n_keep, gs)
            lo, hi = np.nanquantile(T[:, :, 0:3], 0.05, axis=1), np.nanquantile(T[:, :, 0:3], 0.95, axis=1)
            cov[:, kk] &= (lo <= truth) & (truth <= hi)
            ok[:, kk] &= np.asarray(okk, bool)
            ok_by_combo[:, j] &= np.asarray(okk, bool)
            Ts.append(T[:, :, 0:3])
            accs.append(np.asarray(acc, float))
        for i in range(n):
            if not ok[i].all():
                cov_mix[i] = False
                continue
            Pm = np.concatenate([T[i] for T in Ts])
            w = np.concatenate([np.full(len(T[i]), a[i] / len(T[i])) for T, a in zip(Ts, accs)])
            w = w / w.sum()
            for d in range(3):
                lo, hi = wquantile(Pm[:, d], w, 0.05), wquantile(Pm[:, d], w, 0.95)
                cov_mix[i, d] &= bool(lo <= truth[i, d] <= hi)
    return cov, ok, cov_mix, ok_by_combo


def miss_flags(logq_fn, X3: np.ndarray, cond_list: list[np.ndarray], k_list, threshold: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """오지정 점수(−log 아홉 k 균등 혼합 밀도) — logq_fn(X3, cond, gs) → (n,) · 코일 조합마다 내고 가장 작은 점수로 문턱과 비교(모든 조합이 넘을 때만 —
    검토자 판정 ZZZZj-3) · 둘째 = 그 점수(칸 초과율도 이 표시로 낸다) ·
    셋째 = 조합마다 문턱을 넘었는가 (n, 조합 수)(R26 보고 — 수만 · 판정에 쓰지 않음)."""
    by = np.empty((len(X3), len(cond_list)))
    for j, cond in enumerate(cond_list):
        L = np.stack([logq_fn(X3, cond, gs) for gs in k_list], 1)
        m = L.max(1, keepdims=True)
        by[:, j] = -(m[:, 0] + np.log(np.exp(L - m).mean(1)))
    sc = by.min(1)
    return sc > threshold, sc, by > threshold


def rates_by_k(cov: np.ndarray, ok: np.ndarray, use: np.ndarray) -> tuple[np.ndarray, int]:
    """한 성분 cov (n, 9) · ok (n, 9) · use (n,) → k 마다 덮음률 (9,) · 쓴 단위 수.
    아홉 k 모두 표본이 찬 단위만 쓴다(시험 자료 문턱과 같은 규칙 — 검토 F8 · 한 k 라도 «재가중 불가» 인 단위는 뺀다)."""
    m = use & ok.all(1)
    n = int(m.sum())
    return ((cov[m].sum(0) / n) if n else np.full(cov.shape[1], np.nan)), n


def threshold_from_test(cov_t: np.ndarray, ok_t: np.ndarray, flag_t: np.ndarray, n_clusters: int, rng, n_boot: int = N_BOOT) -> float | None:
    """A-1 7 — 시험 자료(한 층 · 한 성분) cov_t (m, 9) · 오지정 경우는 뺌 · n_clusters 개를 복원 추출해 k 마다 덮음률 → 중앙값 → 2,000 번 → 5 백분위."""
    use = np.flatnonzero(~flag_t & ok_t.all(1))
    if n_clusters < 1 or len(use) == 0:
        return None
    draws = rng.choice(use, (n_boot, n_clusters), replace=True)
    med = np.median(cov_t[draws].mean(1), axis=1)                                # (n_boot, 9) → k 중앙값
    return float(np.quantile(med, 0.05))


def exceed_p95(flag_t: np.ndarray, n_clusters: int, rng, n_boot: int = N_BOOT) -> float | None:
    """A-1 4 — 기대 초과율 분포(시험 자료 같은 층 · 같은 독립 묶음 수 · 2,000 번)의 95 백분위."""
    if n_clusters < 1 or len(flag_t) == 0:
        return None
    return float(np.quantile(rng.choice(flag_t.astype(float), (n_boot, n_clusters), replace=True).mean(1), 0.95))


# ── 칸 판정 ─────────────────────────────────────────────────────────────────────
BLIND = {"NDI": "비맹검 점검", "FISO": "비맹검 점검", "ODI": "맹검 확인"}                 # SPEC_S4.md 부록 A (마)


def cell_verdict(cov: np.ndarray, ok: np.ndarray, flag: np.ndarray, subjects: np.ndarray, thr: float | None, rng,
                 p95: float | None = None, n_boot: int = N_BOOT, cov_mix: np.ndarray | None = None, thr_units: float | None = None) -> dict:
    """한 칸(쌍 × 성분) — cov · ok (n, 9) · flag (n,) 오지정 · subjects (n,) 피험자(묶음) · thr 문턱(None 이면 0.80) · cov_mix (n,) 참고 ·
    thr_units = 단위 수로 뽑은 문턱(A-1 7 «참고로만»). 묶음(피험자)이 3 명 미만이면 «판정 불가»(문턱 · p95 가 퇴화 — 검토 F7 · R.PAIR_MIN_SUBJECTS)."""
    t = FALLBACK_THRESHOLD if thr is None else thr
    subs = np.unique(subjects)
    rk, n_used = rates_by_k(cov, ok, ~flag)
    med = float(np.median(rk)) if n_used else None
    rk_in, _ = rates_by_k(cov, ok, np.ones(len(cov), bool))
    few = len(subs) < R.PAIR_MIN_SUBJECTS
    boot = []
    for _ in range(n_boot if not few else 0):
        pick = rng.choice(subs, len(subs), replace=True)
        idx = np.concatenate([np.flatnonzero(subjects == s) for s in pick])
        r, n = rates_by_k(cov[idx], ok[idx], ~flag[idx])
        if n:
            boot.append(float(np.median(r)))
    exceed = float(flag.mean()) if len(flag) else None
    passed = None if (med is None or few) else bool(med >= t)
    return dict(n_units=int(len(cov)), n_units_used=n_used, n_units_dropped_k=int((~flag & ~ok.all(1)).sum()), n_subjects=int(len(subs)),
                too_few_subjects=few, n_excluded_misspec=int(flag.sum()), exceed_rate=exceed, exceed_p95=p95,
                gap_detected=(None if p95 is None or exceed is None or few else bool(exceed > p95)),
                coverage_by_k=[None if not np.isfinite(v) else float(v) for v in rk], coverage=med,
                coverage_min=float(rk.min()) if n_used else None, coverage_max=float(rk.max()) if n_used else None,
                threshold=t, threshold_fallback=thr is None, threshold_units_reference=thr_units, passed=passed,
                assumption_dependent=(None if passed is None else bool(rk.min() < t <= rk.max())),
                coverage_incl_excluded=float(np.median(rk_in)) if np.isfinite(rk_in).all() else None,
                coverage_uniform_mixture_reference=(float(cov_mix[~flag & ok.all(1)].mean()) if cov_mix is not None and n_used else None),
                boot95=[float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))] if boot else None)


def device_concentration(fail_cells: list[str], tested_cells: list[str]) -> dict:
    """[A 읽기 · 기술 · 판정에 쓰지 않음] (라) «장비 공유는 장비 단위 순열로 함께 본다» — 실패 칸이 한 장비에 몰리는가.
    칸 결과를 장비 쌍(NOT2ING 두 프로토콜은 한 장비 · 하나라도 실패면 실패)으로 모은다 · 통계량 = 장비마다 그 장비가 든 실패 쌍 수의 최댓값 ·
    영 = 장비 7 대의 순열 π 아래 쌍 (a, b) 의 결과를 (π(a), π(b)) 로 옮긴다 — 시험한 쌍 집합은 그대로 두고 옮겨 간 자리가 시험한 쌍이 아니면 버린다 ·
    p = #{π : 통계량 ≥ 항등} / 5,040(관측 쌍이 순열에 닫혀 있지 않아 근사)."""
    dev = lambda u: "NOT2ING" if u.startswith("NOT2ING") else u  # noqa: E731
    devs = sorted({dev(u) for u in CARD_ORDER})
    key = lambda a, b: frozenset((a, b))  # noqa: E731
    fail = {}
    for c in tested_cells:
        a, b = (dev(u) for u in c.split("~"))
        if a != b:
            fail[key(a, b)] = fail.get(key(a, b), False) or c in fail_cells
    tested = set(fail)

    def stat(f):
        return max([sum(v for k, v in f.items() if d in k) for d in devs] or [0])
    s0 = stat(fail)
    cnt = 0
    for perm in itertools.permutations(devs):
        m = dict(zip(devs, perm))
        f2 = {}
        for k, v in fail.items():
            k2 = key(*(m[x] for x in k))
            if k2 in tested:
                f2[k2] = v
        cnt += stat(f2) >= s0
    per_dev = {d: dict(failed=sum(v for k, v in fail.items() if d in k), tested=sum(1 for k in fail if d in k)) for d in devs}
    return dict(stat=s0, p=cnt / math.factorial(len(devs)), n_perm=math.factorial(len(devs)), per_device=per_dev)


def decomposition(cov: np.ndarray, ok: np.ndarray, flag: np.ndarray, subjects: np.ndarray, attrs: dict[str, dict]) -> dict:
    """§2-5-1 · a8 — 속성 묶음마다 세션 쌍 수 · 단위 수 · 덮은 단위 수(k 중앙값) · 비율(k 마다 낸 뒤 중앙값) · 검정 없음.
    attrs = {피험자: {"noise": bool, "range": bool, "same_maker": bool}}(속성 표 — V4 값보다 먼저 고정한 파일)."""
    groups = {"① 잡음": lambda a: a["noise"], "② 획득 범위": lambda a: a["range"],
              "③ 장비 같음": lambda a: a["same_maker"], "③ 장비 다름": lambda a: not a["same_maker"],
              "④ 설명 안 됨": lambda a: not (a["noise"] or a["range"] or not a["same_maker"])}
    out = {}
    for g, f in groups.items():
        subs = [s for s in np.unique(subjects) if f(attrs[s])]
        if not subs:
            out[g] = dict(n_session_pairs=0, verdict="해당 없음")
            continue
        m = np.isin(subjects, subs) & ~flag & ok.all(1)
        n = int(m.sum())
        cov_n = cov[m].sum(0)
        out[g] = dict(n_session_pairs=len(subs), n_units=n, n_covered=float(np.median(cov_n)) if n else 0.0,
                      ratio=float(np.median(cov_n / n)) if n else None)
    return out


def judge_cells(identified: dict[str, list[str]], unit_cov: dict, thresholds: dict, p95s: dict, rng, sens: dict | None = None,
                attrs: dict | None = None, thr_units: dict | None = None) -> dict:
    """칸마다 판정 — unit_cov[쌍] = dict(cov (n, 9, 3) · ok (n, 9) · flag (n,) · subjects (n,)) · thresholds[(쌍, 성분)] · p95s[쌍] ·
    sens[(snr, neff)][쌍] = 같은 꼴(민감도 칸) · 자료 없는 쌍(단위 0)은 «자료 없음»(분모에 넣지 않음 — ZZZZh-5) ·
    «식별됨» 이 아닌 칸은 «덮음 — 정보 없음» 으로 따로 적는다((다)) · 칸마다 획득 유형((바) «이 유형») · 맹검 표지((마))."""
    out = dict(cells={}, summary={}, info_none={})
    for c in COMPS:
        out["info_none"][c] = [p for p in card_pairs() if p not in identified.get(c, [])]
    for c, pairs in identified.items():
        ci = COMPS.index(c)
        for pair in pairs:
            u = unit_cov.get(pair)
            base = dict(acq_type=acq_type(pair), blind=BLIND[c])
            if u is None or len(u["cov"]) == 0:
                out["cells"][f"{pair}|{c}"] = dict(base, verdict="자료 없음")
                continue
            r = base | cell_verdict(u["cov"][:, :, ci], u["ok"], u["flag"], u["subjects"], thresholds.get((pair, c)), rng, p95s.get(pair),
                                    cov_mix=u["cov_mix"][:, ci] if "cov_mix" in u else None,
                                    thr_units=(thr_units or {}).get((pair, c)))
            r["verdict"] = "통과" if r["passed"] else "실패" if r["passed"] is not None else "판정 불가"
            if "combos" in u:
                r["by_coil_combo"] = coil_counts(u)
            if sens:
                sv = {}
                for key, d in sens.items():
                    su = d.get(pair)
                    if su is not None and len(su["cov"]):
                        rr = cell_verdict(su["cov"][:, :, ci], su["ok"], su["flag"], su["subjects"], thresholds.get((pair, c)), rng, None, n_boot=0)
                        sv[f"SNR{int(key[0])}|Neff{int(key[1])}"] = rr["passed"]
                r["sensitivity_passed"] = sv
                r["snr_neff_dependent"] = len({v for v in sv.values() if v is not None} | {r["passed"]}) > 1
            if attrs is not None and pair in attrs:
                r["decomposition"] = decomposition(u["cov"][:, :, ci], u["ok"], u["flag"], u["subjects"], attrs[pair])
            out["cells"][f"{pair}|{c}"] = r
    for c in COMPS:
        cells = {k: v for k, v in out["cells"].items() if k.endswith(f"|{c}")}
        out["summary"][c] = dict(tested=len(cells), no_data=sum(v["verdict"] == "자료 없음" for v in cells.values()),
                                 with_data=sum(v["verdict"] != "자료 없음" for v in cells.values()),
                                 passed=sum(v["verdict"] == "통과" for v in cells.values()),
                                 failed=sum(v["verdict"] == "실패" for v in cells.values()),
                                 not_judgeable=sum(v["verdict"] == "판정 불가" for v in cells.values()),
                                 info_none=len(out["info_none"][c]), blind=BLIND[c])
    out["device_concentration_descriptive"] = {}
    for c in COMPS:
        tested = [k.split("|")[0] for k, v in out["cells"].items() if k.endswith(f"|{c}") and v["verdict"] in ("통과", "실패")]
        failed = [k.split("|")[0] for k, v in out["cells"].items() if k.endswith(f"|{c}") and v["verdict"] == "실패"]
        if tested:
            out["device_concentration_descriptive"][c] = device_concentration(failed, tested)
    return out


# ── ⑩ 뒤 실행 — 순서: snr-table → attr-table → test-cover → run(R18 · R20 · §2-5-1 · 검토 F3) ─────────────────────
ARMS = ("c3", "c7")
V4_DIR = "v4"                                                               # 결과 · 표는 WORK/onharmony/v4/


def sha256_file(path: pathlib.Path) -> str:
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def _write_json(path: pathlib.Path, obj) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(obj, ensure_ascii=False, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)) + "\n",
                   encoding="utf-8")
    tmp.replace(path)
    return sha256_file(path)


def _write_once(path: pathlib.Path, obj) -> str:
    """한 번만 쓴다 — 이미 있으면 멈춘다(V4 값을 본 뒤 표를 다시 만들어 돌리지 않게 · 반영 확인 읽기 F3) · 다시 만들어야 하면 이탈로 기록하고 손으로 지운 뒤."""
    if path.exists():
        raise RuntimeError(f"⛔ {path} 가 이미 있다 — 다시 쓰지 않는다(다시 만들어야 하면 이탈로 기록하고 손으로 지운 뒤)")
    return _write_json(path, obj)


def _no_v4_values(root: pathlib.Path) -> None:
    """SNR 표 · 속성 표는 V4 값보다 먼저만 만든다(R18 · R20) — V4 결과 파일이 하나라도 있으면 멈춘다."""
    done = sorted(p.name for p in (root / V4_DIR).glob("v4_*.json"))
    if done:
        raise RuntimeError(f"⛔ V4 값이 이미 있다({done}) — SNR 표 · 속성 표는 V4 값보다 먼저만 만든다(R18 · R20)")


def _log(what: str, path: pathlib.Path, basis: dict | None) -> None:
    """열람 기록(rconfig.ACCESS_LOG_ONH «열람»)에 시각(UTC) · 근거 커밋 · 파일 · sha256 — run.log 와 같은 기록 · 잠금(R18 · 반영 확인 읽기 F3)."""
    import run as RUN  # noqa: E402
    RUN.log("onharmony", what, [path], basis or None, "V4")


def _require_sha(path: pathlib.Path, expect: str, what: str) -> None:
    if not path.exists() or sha256_file(path) != expect:
        raise RuntimeError(f"⛔ {what} {path} 의 sha256 이 주어진 값과 다르거나 파일이 없다 — 순서(SNR 표 → 속성 표 → V4)를 지킨다")


def _validate():
    import validate as V  # noqa: E402 — torch · sbi(메인 환경)
    return V


def real_fns(arm: str):
    """동결 판단기(앙상블 다섯 · q) — CUDA 가 보이면 멈춘다(CPU 만 · 검토 F11) · sha256 대조 뒤에 연다 · S4 실패 층이 있으면 멈춘다(F12)."""
    import torch  # noqa: E402
    if torch.cuda.is_available():
        raise SystemExit("⛔ CUDA 가 보인다 — CPU 만(PLAN §9-1 #49). CUDA_VISIBLE_DEVICES= 로 다시 돌려라")
    check_frozen()
    fails = json.loads((RESULTS / f"test_{arm}.json").read_text(encoding="utf-8")).get("n_layers_fail")
    if fails:
        raise RuntimeError(f"⛔ S4 시험 자료에 실패 층 {fails} — 그 층 유형의 칸은 «판정 불가»(§7 i) · 이 코드는 실패 층 0 을 전제한다")
    V = _validate()
    posts = V.load_posts(arm, 5)
    est = torch.load(S.WORK / "models" / "q_0.obj.pt", weights_only=False)

    def draw_fn(X, cond, n_keep, gs):
        T, _, _, ok, acc = V.draw(posts, X, cond, n_keep, gs)
        return T, ok, acc

    def logq_fn(X3, cond, gs):
        a = torch.as_tensor(V.with_gs(cond, *gs), dtype=torch.float32)
        x = torch.as_tensor(X3, dtype=torch.float32)
        with torch.no_grad():
            return est.log_prob(x[None], condition=a).reshape(-1).numpy()
    return draw_fn, logq_fn, V.GEN_SCEN


def seed_torch(seed: int) -> int:
    """사후 표본의 torch 난수를 고정한다(검토 F1 — 다시 돌려도 같은 덮음) · 쓴 시드를 돌려준다(출력에 적음)."""
    import torch  # noqa: E402
    torch.manual_seed(seed)
    return seed


def test_cover(arm: str, out: pathlib.Path) -> str:
    """시험 자료 L0 다섯 층(획득 유형)의 경우마다 k 별 «참 Δ 덮음» · 오지정 표시 → npz(문턱 · 기대 초과율의 재료 · ⑩ 뒤) · sha256 을 돌려준다."""
    draw_fn, logq_fn, K = real_fns(arm)
    V = _validate()
    sim, _ = V.load_set("test")
    miss = json.loads((RESULTS / "miss.json").read_text(encoding="utf-8"))
    res = {}
    for ti, t in enumerate(S.ACQ_TYPES):
        seed_torch(5000 + 100 * ARMS.index(arm) + ti)
        idx = np.flatnonzero(sim["layer"] == f"{t}|L0")
        X, cond, th = V.xobs(sim, arm)[idx], sim["cond"][idx], sim["theta"][idx, 0:3]
        cov, ok, _, _ = cover_by_k(draw_fn, X, [cond], K, N_KEEP, th)
        flag, _, _ = miss_flags(logq_fn, sim["x_c3"][idx], [cond], K, miss["layers"][f"{t}|L0"]["threshold"])
        res.update({f"{t}_cov": cov, f"{t}_ok": ok, f"{t}_flag": flag})
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(out, **res)
    return sha256_file(out)


USABLE = "쓸 수 있음"


def session_states(work: pathlib.Path, sessions) -> dict[str, str]:
    """세션마다 상태 — 다시도 실패(FAILED.json · batch_summary «실패_다시도_실패») · 충돌(끝 표지 + 실패 둘) · 끝 표지 없음 · 쓸 수 있음
    (검토 F10 · 요약이 없으면 멈춤 · 쌍 검정 로더 judge_onh 도 같은 규칙 — 반복 세션 포함)."""
    root = work / "onharmony"
    sp = root / "batch_summary.json"
    if not sp.exists():
        raise RuntimeError(f"⛔ {sp} 이 없다 — 일괄 처리의 실패 · 충돌 기록 없이 판정하지 않는다")
    summ = json.loads(sp.read_text(encoding="utf-8"))
    import run as RUN  # noqa: E402 — prior_failures(읽기만)
    fails = RUN.prior_failures(work)
    out = {}
    for sub, ses in sessions:
        key, wd = f"{sub}_{ses}", root / f"{sub}_{ses}"
        if key in summ.get("실패_다시도_실패", []) or (wd / "FAILED.json").exists():
            out[key] = "다시도 실패"
        elif key in summ.get("충돌", []) or ((wd / "DONE.json").exists() and len(fails.get((sub, ses), [])) >= 2):
            out[key] = "충돌"
        elif not (wd / "DONE.json").exists():
            out[key] = "끝 표지 없음"
        else:
            out[key] = USABLE
    return out


def usable_sessions(work: pathlib.Path) -> tuple[dict, dict]:
    """처리가 끝난 첫 세션(DONE.json) — 다시도 실패와 충돌은 뺀다(session_states) → ({카드 단위: {피험자: 세션}}, 뺀 세션과 까닭)."""
    st = session_states(work, R.onh_sessions())
    dropped, by_unit = {}, {}
    for sub, ses in R.onh_sessions():
        key = f"{sub}_{ses}"
        if st[key] != USABLE:
            dropped[key] = st[key]
            continue
        d = R.device_of(ses)
        by_unit.setdefault(R.not2ing_protocol(sub) if d == "NOT2ING" else d, {})[sub] = ses
    return by_unit, dropped


def snr_table(work: pathlib.Path) -> tuple[pathlib.Path, str]:
    """1 단계(⑩ 뒤 · V4 값보다 먼저) — 세션마다 snr.json 을 모은 한 파일 → 경로 · sha256(열람 기록에 sha256 · 시각을 적는다 · R18)."""
    import guard as G  # noqa: E402
    basis = G.require()
    root = work / "onharmony"
    _no_v4_values(root)
    tab = {f"{sub}_{ses}": json.loads((root / f"{sub}_{ses}" / "snr.json").read_text(encoding="utf-8"))
           for sub, ses in R.onh_sessions() + R.sref_extra() if (root / f"{sub}_{ses}" / "snr.json").exists()}
    path = root / V4_DIR / "snr_table.json"
    sha = _write_once(path, dict(what="세션별 SNR 표(R18 · analysis.snr) — V4 값보다 먼저", train_range=list(R.SNR_TRAIN), sessions=tab))
    _log("V4 1 단계 — 세션별 SNR 표(R18)", path, basis)
    return path, sha


def attr_table(work: pathlib.Path, snr_sha: str) -> tuple[pathlib.Path, str]:
    """2 단계 — 속성 표(§2-5-1 · a8 — 카드 쌍 × 피험자마다 ① 잡음 · ② 획득 범위 · ③ 제조사 같음) · SNR 표의 sha256 이 맞아야 만든다."""
    import guard as G  # noqa: E402
    basis = G.require()
    root = work / "onharmony"
    _no_v4_values(root)
    sp = root / V4_DIR / "snr_table.json"
    _require_sha(sp, snr_sha, "SNR 표")
    snr = json.loads(sp.read_text(encoding="utf-8"))["sessions"]
    by_unit, dropped = usable_sessions(work)
    tab = {}
    for pair in card_pairs():
        a, b = pair.split("~")
        for sub in sorted(set(by_unit.get(a, {})) & set(by_unit.get(b, {}))):
            ss = (by_unit[a][sub], by_unit[b][sub])
            meta = [json.loads((R.ONH_META / sub / s / "dwi" / f"{sub}_{s}_dir-AP_dwi.json").read_text(encoding="utf-8")) for s in ss]
            shells = [np.loadtxt(R.ONH_META / sub / s / "dwi" / f"{sub}_{s}_dir-AP_dwi.bval").ravel() for s in ss]
            tab.setdefault(pair, {})[sub] = dict(
                sessions=list(ss), te_ms=[m["EchoTime"] * 1000 for m in meta], tr_ms=[m["RepetitionTime"] * 1000 for m in meta],
                noise=any(not snr.get(f"{sub}_{s}", {}).get("in_train", False) for s in ss),
                range=any(not (60 <= m["EchoTime"] * 1000 <= 110 and 3.0 <= m["RepetitionTime"] <= 11.0
                               and set(np.round(bv[bv > 50] / 500) * 500) <= {1000.0, 2000.0}) for m, bv in zip(meta, shells)),
                same_maker=meta[0].get("Manufacturer") == meta[1].get("Manufacturer"),
                main_excluded=any((sub, s) in R.MAIN_EXCLUDE_ONH for s in ss))
    path = root / V4_DIR / "attr_table.json"
    sha = _write_once(path, dict(what="V4 속성 표(§2-5-1) — V4 값보다 먼저", snr_table_sha256=snr_sha, dropped_sessions=dropped, pairs=tab))
    _log("V4 2 단계 — 속성 표(§2-5-1 · R20)", path, basis)
    return path, sha


def collect(work: pathlib.Path, attrs: dict) -> dict:
    """관측 행 — 속성 표의 (쌍, 피험자)마다 두 세션 roi.json «주» → ROI 단위별 C7 14 · TE · TR(속성 표의 값) · 16794 표시."""
    import guard as G  # noqa: E402
    G.require()
    root = work / "onharmony"
    rows = {}
    for pair in card_pairs():
        out = []
        for sub, at in sorted(attrs.get(pair, {}).items()):
            ra, rb = (json.loads((root / f"{sub}_{s}" / "roi.json").read_text(encoding="utf-8")) for s in at["sessions"])
            units = list(ra["주"])
            xa, xb = session_x(ra, units), session_x(rb, units)
            for u in sorted(set(xa) & set(xb)):
                out.append(dict(sub=sub, unit=u, x=np.r_[xa[u], xb[u]], te=at["te_ms"], tr=at["tr_ms"], main_excluded=at["main_excluded"],
                                noise=at["noise"]))
        rows[pair] = out
    return rows


def unit_results(pair: str, rows: list[dict], arm: str, draw_fn, logq_fn, K, miss_thr: float, snr: float, neff: float) -> dict:
    """한 쌍의 단위 결과 — X(C3 = 사이트 A 의 앞 셋 + 사이트 B 의 앞 셋 · C7 = 일곱 + 일곱) · 조건(코일 조합마다) · k 아홉 덮음 · 오지정 ·
    조합마다의 오지정 표시 · 표본이 찼는가(R26)."""
    a, b = pair.split("~")
    combos = [(ca, cb) for ca in coil_options(a) for cb in coil_options(b)]
    names = [f"{ca}|{cb}" for ca, cb in combos]
    if not rows:
        m = len(combos)
        return dict(cov=np.zeros((0, 9, 3), bool), ok=np.zeros((0, 9), bool), cov_mix=np.zeros((0, 3), bool), flag=np.zeros(0, bool),
                    subjects=np.array([]), main_excluded=np.zeros(0, bool), noise=np.zeros(0, bool),
                    combos=names, flag_by_combo=np.zeros((0, m), bool), ok_by_combo=np.zeros((0, m), bool))
    X7 = np.array([r["x"] for r in rows])
    X3 = np.c_[X7[:, 0:3], X7[:, 7:10]]
    X = X7 if arm == "c7" else X3
    te, tr = np.array([r["te"] for r in rows], float), np.array([r["tr"] for r in rows], float)
    conds = [cond_rows(pair, te, tr, cc, snr, neff) for cc in combos]
    cov, ok, cov_mix, ok_by_combo = cover_by_k(draw_fn, X, conds, K)
    flag, _, flag_by_combo = miss_flags(logq_fn, X3, conds, K, miss_thr)
    return dict(cov=cov, ok=ok, cov_mix=cov_mix, flag=flag, subjects=np.array([r["sub"] for r in rows]),
                main_excluded=np.array([r["main_excluded"] for r in rows]), noise=np.array([r["noise"] for r in rows]),
                combos=names, flag_by_combo=flag_by_combo, ok_by_combo=ok_by_combo)


def coil_counts(u: dict) -> dict[str, dict[str, int]]:
    """R26(검토자 제안 ZZZZi-3 · a15) — 칸 · 코일 조합(«코일_A|코일_B»)마다 그 조합의 오지정 점수가 문턱을 넘은 단위 수 ·
    그 조합에서 아홉 k 중 표본이 차지 않은 k 가 있는 단위 수 · 수만 · 판정에 쓰지 않음(판정은 위 규칙 — 오지정은 모든 조합이 넘을 때만 · 덮음과 표본은 모든 조합 ·
    조합마다의 수는 최댓값으로 읽었으면 더 뺐을 단위의 재료)."""
    return {c: dict(n_misspec=int(u["flag_by_combo"][:, j].sum()), n_k_unfilled=int((~u["ok_by_combo"][:, j]).sum()))
            for j, c in enumerate(u["combos"])}


def select(u: dict, keep: np.ndarray) -> dict:
    return {k: (v[keep] if isinstance(v, np.ndarray) and len(v) == len(keep) else v) for k, v in u.items()}


def drawn_thresholds(tc: dict, arm: str, seed: int):
    """문턱(A-1 7 · 5 백분위) · 기대 초과율(A-1 4 · 95 백분위)을 (획득 유형 · 성분 · 묶음 수)만의 함수로 — 같은 입력이면 변형 · 쌍이 달라도
    같은 난수로 같은 값(주 분석과 민감도의 차이에 몬테카를로 잡음이 섞이지 않게 · 반영 확인 읽기 F2) · get(«thr» | «p95», 유형, 묶음 수, 성분 번호)."""
    cache = {}

    def get(kind: str, t: str, n: int, ci: int = -1):
        key = (kind, t, n, ci)
        if key not in cache:
            g = np.random.default_rng([seed, ARMS.index(arm), S.ACQ_TYPES.index(t), ci + 1, n, int(kind == "p95")])
            cache[key] = (exceed_p95(tc[f"{t}_flag"], n, g) if kind == "p95"
                          else threshold_from_test(tc[f"{t}_cov"][:, :, ci], tc[f"{t}_ok"], tc[f"{t}_flag"], n, g))
        return cache[key]
    return get


def run(arm: str, work: pathlib.Path, attr_sha: str, test_npz: pathlib.Path, seed: int = 20261027) -> tuple[dict, str]:
    """4 단계 — 실측 V4(주 칸 · 민감도 칸 넷 · 기술 분해 · 16794_NOT4GEP001 넣은 민감도 · SNR 범위 밖 뺀 민감도) → 결과 파일 · sha256.
    속성 표의 sha256 이 맞아야 돈다 · 단위 결과는 한 번만 뽑고 민감도 둘은 행만 고른다(검토 F2) · 칸마다 torch 시드 고정(F1) ·
    문턱 · 기대 초과율은 입력만의 함수(drawn_thresholds) · 결과 파일은 한 번만 쓰고 열람 기록에 sha256(반영 확인 읽기 F2 · F3)."""
    import guard as G  # noqa: E402
    basis = G.require()
    root = work / "onharmony"
    ap = root / V4_DIR / "attr_table.json"
    _require_sha(ap, attr_sha, "속성 표")
    at = json.loads(ap.read_text(encoding="utf-8"))
    draw_fn, logq_fn, K = real_fns(arm)
    miss = json.loads((RESULTS / "miss.json").read_text(encoding="utf-8"))
    cards = json.loads((RESULTS / f"cards_{arm}.json").read_text(encoding="utf-8"))
    ident = identified_cells(cards)
    tc = dict(np.load(test_npz))
    drawn = drawn_thresholds(tc, arm, seed)
    rows = collect(work, at["pairs"])
    seeds = {}
    cells = [MAIN_CELL] + list(SENS_CELLS)
    units = {}
    for ci_, (_, snr, neff) in enumerate(cells):
        units[(snr, neff)] = {}
        for pi, pair in enumerate(card_pairs()):
            seeds[f"SNR{int(snr)}|Neff{int(neff)}|{pair}"] = seed_torch(6000 + 1000 * ARMS.index(arm) + 100 * ci_ + pi)
            units[(snr, neff)][pair] = unit_results(pair, rows[pair], arm, draw_fn, logq_fn, K,
                                                    miss["layers"][f"{acq_type(pair)}|L0"]["threshold"], snr, neff)
    res = dict(arm=arm, seeds_torch=seeds, seed_numpy=seed, inputs=dict(
        attr_table_sha256=attr_sha, test_cover_sha256=sha256_file(test_npz),
        design_files_sha256={f"{S.DESIGNS[d]}_{e}": sha256_file(S.C.PROTO_DIR / f"{S.DESIGNS[d]}_{e}") for d in S.DESIGNS for e in ("bval", "bvec")},   # 조건의 b0 를 정하는 설계 파일(F13)
        frozen=frozen_expected()))
    variants = (("main", lambda u: ~u["main_excluded"]), ("with_16794_NOT4GEP001", lambda u: np.ones(len(u["subjects"]), bool)),
                ("without_snr_out_of_train", lambda u: ~u["main_excluded"] & ~u["noise"]))           # R11 · R18 — SNR 범위 밖 세션 쌍 뺌
    main_key = (MAIN_CELL[1], MAIN_CELL[2])
    for label, keep_fn in variants:
        rng = np.random.default_rng(seed)                                    # 칸 부트스트랩(구간 · 기술) — 변형마다 같은 출발
        main = {p: select(u, keep_fn(u)) for p, u in units[main_key].items()}
        thr, thr_u, p95 = {}, {}, {}
        for pair, u in main.items():
            n_cl, n_un = len(np.unique(u["subjects"])), len(u["subjects"])
            t = acq_type(pair)
            p95[pair] = drawn("p95", t, n_cl)
            for ci, c in enumerate(COMPS):
                thr[(pair, c)] = drawn("thr", t, n_cl, ci)
                thr_u[(pair, c)] = drawn("thr", t, n_un, ci)                 # 참고(A-1 7)
        sens = ({(s, n): {p: select(u, keep_fn(u)) for p, u in units[(s, n)].items()} for _, s, n in SENS_CELLS} if label == "main" else None)
        attrs = ({p: {s: v for s, v in d.items()} for p, d in at["pairs"].items()} if label == "main" else None)
        res[label] = judge_cells(ident, main, thr, p95, rng, sens, attrs, thr_u)
    path = root / V4_DIR / f"v4_{arm}.json"
    sha = _write_once(path, res)
    _log(f"V4 4 단계 — 결과 v4_{arm}.json", path, basis)
    return res, sha


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "snr-table":
        print(*snr_table(R.WORK))
    elif cmd == "attr-table":
        print(*attr_table(R.WORK, sys.argv[2]))
    elif cmd == "test-cover":
        arm = sys.argv[2]
        print(S.WORK / "v4" / f"test_cover_{arm}.npz", test_cover(arm, S.WORK / "v4" / f"test_cover_{arm}.npz"))
    elif cmd == "run":
        arm, attr_sha = sys.argv[2], sys.argv[3]
        _, sha = run(arm, R.WORK, attr_sha, S.WORK / "v4" / f"test_cover_{arm}.npz")
        print(R.WORK / "onharmony" / V4_DIR / f"v4_{arm}.json", sha)
    else:
        raise SystemExit(cmd)
