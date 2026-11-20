"""실측 고정 문서 §3 — ODI (a) 동등성 여유 δ_τ 를 설계만으로 정한다(중계 ZZZw-3 2 · [검토자 제안] · 계산 방법은 A 제안 — 검토자가 본다). ON-Harmony 값 0.

ODI (a): 21 장비 쌍의 Kendall τ_b(예측 이완 이동 크기 x · |관측 ODI 차| y) 의 90 % 구간이 (−δ_τ, δ_τ) 안이면 «동등».
- x: 카드의 L1 최대 이동(`analyze3.card_algebra` · 카드 28 장과 같은 쌍 정의) — NOT2ING 의 두 프로토콜은 평균 · NOT2ING 안의 쌍은 뺀다(장비 사이가 아니다).
- 영(참 τ = 0): y_ij = |u_i − u_j + e_ij| · u(장비 효과) ~ N(0, r²) · e(쌍 잡음) ~ N(0, 1) — 장비가 교환 가능한 잡음. r 은 격자로 보고 가장 큰 δ_τ 를 쓴다(보수).
- 90 % 구간(A 제안 — 순열 역산): 장비 단위 정확 순열(7! = 5,040)로 그 y 의 영분포 T 를 내고, 위치 이동을 가정해 [τ̂ − q95(T), τ̂ − q05(T)] — 참 τ = 0 이면 0 을 0.90 으로 덮는다(구성상).
- δ_τ = 참 τ = 0 에서 P(구간 ⊂ (−δ, δ)) ≥ 0.80 인 가장 작은 δ(0.01 격자). 미리 정함: δ_τ > 0.5 면 (a) 는 «기술만 — 검정력 부족».
FISO 부호 S(a9) · FISO τ(a2)도 예측 쪽 장비 순열 — 예측 부호는 설계 부호(공통 피험자 조건 없이 · 9 쌍 · 재점검 C-2) · fiso_sizes 는 영에서의 크기.

    uv run --frozen python decomp/fixed/delta_tau.py [모의 수 · 기본 1000]     # JSON
    uv run --frozen python decomp/fixed/delta_tau.py --devices                 # 판정 불가 쌍을 뺀 19 쌍 · NOT2ING 피험자 가중 x(중계 ZZZZf-2 ① ③)
    uv run --frozen python decomp/fixed/delta_tau.py --fiso-size [모의 수 · 시드 수]          # FISO 부호 S · FISO τ 의 영 크기(설계만 · 기본 4000 · 시드 1)
    uv run --frozen python decomp/fixed/delta_tau.py --fiso-size-before [모의 수 · 시드 수]   # 고치기 전 방식(지운 부호 · 항등 8 쌍 고정)의 크기 — 비교용
"""
from __future__ import annotations

import itertools
import json
import pathlib
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "stage1"))

DEVICES = ("NOT1ACH", "NOT2ING", "NOT4GEP", "OXF1PRI", "OXF2PRI", "OXF3TRI", "OXF4GEP")
PAIRS = list(itertools.combinations(range(len(DEVICES)), 2))                   # 21
R_GRID = (0.0, 0.5, 1.0, 2.0, np.inf)                                           # 장비 효과 / 쌍 잡음(∞ = 장비 효과만)
SEED = 20261026


def device_of(unit: str) -> str:
    return "NOT2ING" if unit.startswith("NOT2ING") else unit


def pair_predictions(level: str = "L1") -> np.ndarray:
    """21 장비 쌍(PAIRS 순)의 예측 이동 크기 — 카드 L1 최대 이동 · NOT2ING 두 프로토콜 평균 · NOT2ING 안 쌍 뺌."""
    import analyze3 as A3  # noqa: E402
    import fiso_sign as FS  # noqa: E402
    acc = defaultdict(list)
    for a, b, _, teA, trA, teB, trB in FS.pairs():
        da, db = device_of(a), device_of(b)
        if da == db:
            continue
        i, j = sorted((DEVICES.index(da), DEVICES.index(db)))
        acc[(i, j)].append(A3.card_algebra(teA, trA, teB, trB)[level]["max_shift"])
    assert len(acc) == len(PAIRS), len(acc)
    return np.array([np.mean(acc[p]) for p in PAIRS])


def pair_predictions_weighted(counts: dict, level: str = "L1") -> np.ndarray:
    """21 장비 쌍의 예측 이동 크기 — NOT2ING 쌍은 공통 피험자의 프로토콜 비율로 가중(중계 ZZZZf-2 ① · [A 제안] · `fiso_sign.device_rows` 와 같은 규칙).
    counts = `decomp/remeasure/analysis.py` common_counts(세션 목록 · 메타데이터)."""
    import analyze3 as A3  # noqa: E402
    import fiso_sign as FS  # noqa: E402
    acc = defaultdict(list)
    for a, b, _, teA, trA, teB, trB in FS.pairs():
        da, db = device_of(a), device_of(b)
        if da == db:
            continue
        unit = a if a.startswith("NOT2ING") else b if b.startswith("NOT2ING") else "-"
        c = counts[frozenset((da, db))]
        n = c["protocols"].get(unit, 0) if unit != "-" else 1
        i, j = sorted((DEVICES.index(da), DEVICES.index(db)))
        acc[(i, j)].append((n, A3.card_algebra(teA, trA, teB, trB)[level]["max_shift"]))
    out = []
    for p in PAIRS:
        w = [n for n, _ in acc[p]]
        out.append(sum(n * v for n, v in acc[p]) / sum(w) if sum(w) else float(np.mean([v for _, v in acc[p]])))
    return np.array(out)


def observed_finite(v, keep, what: str) -> np.ndarray:
    """관측 쌍(keep)의 값이 모두 유한해야 한다 — 판정 가능 쌍의 값이 없으면(처리 실패 등) 조용히 p = 0 이 나오지 않게 멈춘다(재점검 FISO-3).
    그런 쌍은 v5 §4 규칙대로 판정 불가로 옮겨 keep 에서 뺀 뒤 다시 부른다. 관측 없는 쌍의 값(NaN 포함)은 쓰지 않는다."""
    v = np.asarray(v, float)
    bad = [f"{DEVICES[PAIRS[k][0]]}~{DEVICES[PAIRS[k][1]]}" for k in np.asarray(keep) if not np.isfinite(v[k])]
    if bad:
        raise ValueError(f"⛔ {what}: 관측 쌍 값 없음(NaN · None) — {' · '.join(bad)} · 판정 불가로 옮겨 keep 에서 뺀 뒤 다시(v5 §4)")
    return v


def perm_taus_x(x, y, keep, idx) -> np.ndarray:
    """관측이 있는 쌍(keep)만 있을 때 — 장비 순열을 예측 x 쪽에 건다: x_π[k] = x[π(i), π(j)] · 순열마다 τ_b(x_π[keep], y[keep]).
    x 는 21 쌍 모두 정의돼 있다(설계). 쌍이 모두 있으면 perm_taus 와 같은 분포다(정답 시험)."""
    x, y = np.asarray(x, float), observed_finite(y, keep, "ODI (a) τ")
    keep = np.asarray(keep)
    i, j = np.triu_indices(len(keep), 1)
    yk = y[keep]
    sy = np.sign(yk[i] - yk[j])
    X = x[idx][:, keep]
    SX = np.sign(X[:, i] - X[:, j])
    den = np.sqrt((SX != 0).sum(1) * (sy != 0).sum())
    return np.where(den > 0, (SX @ sy) / np.where(den > 0, den, 1.0), 0.0)


def interval90_subset(x, y, keep, idx):
    """판정 불가 쌍을 뺀 관측 쌍(keep)의 순열 역산 90 % 구간 — interval90 과 같은 규칙."""
    T = perm_taus_x(x, y, keep, idx)
    t = tau_b(np.asarray(x)[keep], np.asarray(y)[keep])
    q05, q95 = np.quantile(T, [0.05, 0.95])
    return t, max(-1.0, t - q95), min(1.0, t - q05)


def design_subset(keep, x, n_sim=1000, seed=SEED, deltas=np.round(np.arange(0.0, 1.01, 0.01), 2)):
    """design 과 같은 규칙을 관측 쌍(keep)에 — 판정 불가 쌍(공통 피험자 < 3)을 미리 뺀 쌍 집합(중계 ZZZZf-2 ③)."""
    x = np.asarray(x, float)
    idx = device_perm_index()
    rng = np.random.default_rng(seed)
    res = {}
    for r in R_GRID:
        lo, hi, cov = [], [], 0
        for _ in range(n_sim):
            t, a, b = interval90_subset(x, null_y(rng, r), keep, idx)
            lo.append(a), hi.append(b)
            cov += a <= 0 <= b
        lo, hi = np.array(lo), np.array(hi)
        power = [float(((lo > -d) & (hi < d)).mean()) for d in deltas]
        ok = [d for d, p in zip(deltas, power) if p >= 0.80]
        res["inf" if np.isinf(r) else str(r)] = dict(coverage0=cov / n_sim, delta_tau=float(ok[0]) if ok else None)
    dts = [v["delta_tau"] for v in res.values()]
    dt = None if any(v is None for v in dts) else max(dts)
    return dict(x=x.tolist(), keep=[f"{DEVICES[PAIRS[k][0]]}~{DEVICES[PAIRS[k][1]]}" for k in keep], n_sim=n_sim, seed=seed, by_r=res,
                delta_tau=dt, descriptive_only=(dt is None or dt > 0.5))


def sign_vector(rows) -> np.ndarray:
    """장비 단위 부호 예측(`fiso_sign.device_rows` — 쌍 이름 «a~b» 의 부호 = b − a 의 부호) → PAIRS(장비 이름 순 i < j) 의 s_ij = (D_j − D_i) 예측 부호.
    이름 순서가 PAIRS 와 반대면 부호를 바꾼다 · 부호 예측이 없는 쌍(«—» · 판정 불가)은 0."""
    s = np.zeros(len(PAIRS))
    for r in rows:
        if not r["in_sign_prediction"]:
            continue
        a, b = r["pair"].split("~")
        v = 1.0 if r["sign"] == "+" else -1.0
        i, j = DEVICES.index(a), DEVICES.index(b)
        s[PAIRS.index((min(i, j), max(i, j)))] = v if i < j else -v
    return s


def design_signs(counts: dict) -> np.ndarray:
    """설계 부호 [A 제안 a9 — 재점검 C-2] — 공통 피험자 조건 없이 낸 장비 단위 부호 예측(`fiso_sign.device_rows(min_subjects=0)`) → PAIRS 순.
    판정 불가 쌍도 설계로는 부호가 있다(OXF3TRI~NOT4GEP − → 9 쌍) — 순열이 그 예측을 관측 쌍으로 옮기면 그 항을 더한다.
    관측 없는 쌍은 항등에서 더해지지 않으므로 항등의 통계량 · 일치 비율은 판정 가능한 8 쌍 그대로다."""
    import fiso_sign as FS  # noqa: E402
    return sign_vector(FS.device_rows(counts, min_subjects=0))


def magnitude_vector(rows, key: str = "fiso_max_abs_dlogit_L1") -> np.ndarray:
    """장비 단위 예측 크기(방향 없음) → PAIRS 순 — FISO τ 의 x(21 쌍 모두 정의)."""
    out = np.full(len(PAIRS), np.nan)
    for r in rows:
        a, b = r["pair"].split("~")
        i, j = DEVICES.index(a), DEVICES.index(b)
        out[PAIRS.index((min(i, j), max(i, j)))] = r[key]
    assert np.isfinite(out).all()
    return out


def sign_perm_stats(s, o, keep) -> np.ndarray:
    """S_π(순열 5,040 · 첫째 = 항등) — S_π = Σ_{k ∈ keep} s(π(i_k), π(j_k)) · o_k. s 는 반대칭으로 읽는다(s(b, a) = −s(a, b))."""
    s, o = np.asarray(s, float), observed_finite(o, keep, "FISO 부호")
    keep = np.asarray(keep)
    return (s[device_perm_index()] * device_perm_orient())[:, keep] @ o[keep]


def sign_perm(s, o, keep) -> dict:
    """FISO 부호 검정 [A 제안 a9 — 고정 문서 v5 §3 · 재점검 C-2] — 장비 단위 순열 5,040 을 모두 예측 쪽에 건다
    (관측 쌍이 순열에 닫혀 있지 않아 정확 검정은 아니다 — 영 크기는 fiso_sizes).
    S_π = Σ_{k ∈ keep} s(π(i_k), π(j_k)) · o_k — s · o 는 PAIRS 순의 예측 부호 · 관측 부호(D_j − D_i 의 부호 · 0 은 0).
    s = 설계 부호(design_signs · 9 쌍) · 반대칭(순열이 방향을 뒤집으면 부호가 바뀐다) · 예측이 없는 쌍(0)은 그 순열에서 더하지 않는다.
    keep = 판정 가능한 관측 쌍(19) — 순열마다 같은 관측 집합. 한쪽 p = #{π : S_π ≥ S_id} / 5,040 · 효과 크기 = 항등 배정의 일치 비율."""
    s, o = np.asarray(s, float), observed_finite(o, keep, "FISO 부호")
    S = sign_perm_stats(s, o, keep)
    both = [k for k in keep if s[k] != 0 and o[k] != 0]
    agree = sum(1 for k in both if s[k] == o[k])
    return dict(S=float(S[0]), p=float((S >= S[0]).mean()), n_perm=len(S), n_pairs=len(both), agree=agree,
                agree_rate=(agree / len(both) if both else None))


def perm_taus_sel(x, y, s, keep, idx) -> np.ndarray:
    """FISO τ [A 제안 a2 — 재점검 C-2] — 순열마다 쌍을 다시 고른다: keep_π = {k ∈ keep(관측 쌍) : s(π(i_k), π(j_k)) ≠ 0}
    (s = 설계 부호 · 크기만 보므로 방향은 쓰지 않는다) · τ_π = τ_b(x_π[keep_π], y[keep_π]) · x_π[k] = x(π(i_k), π(j_k)).
    항등의 keep_π = 관측 쌍 중 부호 예측이 있는 쌍(8). 고른 쌍이 둘 미만이면 0."""
    x, s = np.asarray(x, float), np.asarray(s, float)
    obs = np.zeros(len(PAIRS), bool)
    obs[np.asarray(keep)] = True
    y = np.where(obs, observed_finite(y, keep, "FISO τ"), 0.0)                   # 관측 없는 쌍의 y 는 쓰이지 않는다(NaN 이어도 됨)
    M = obs[None, :] & (s[idx] != 0)                                             # (순열, 21) — 순열마다 고른 쌍
    a, b = np.triu_indices(len(PAIRS), 1)
    MM = M[:, a] & M[:, b]
    X = x[idx]
    SX = np.sign(X[:, a] - X[:, b]) * MM
    SY = np.sign(y[a] - y[b])[None, :] * MM
    den = np.sqrt((SX != 0).sum(1) * (SY != 0).sum(1))
    return np.where(den > 0, (SX * SY).sum(1) / np.where(den > 0, den, 1.0), 0.0)


def fiso_tau(x, y, s, keep) -> dict:
    """FISO 크기 순서 — τ_b(x · y)(항등의 keep_π = 8 쌍) · 한쪽 p = #{π : τ_π ≥ τ_id} / 5,040(perm_taus_sel)."""
    T = perm_taus_sel(x, y, s, keep, device_perm_index())
    return dict(tau=float(T[0]), p=float((T >= T[0] - 1e-12).mean()), n_perm=len(T),
                n_pairs=int(sum(1 for k in keep if np.asarray(s, float)[k] != 0)))


def null_d(rng, r: float, spike: float | None = None) -> np.ndarray:
    """영의 쌍 차(PAIRS 순 · D_j − D_i) — u ~ N(0, r²)(장비 효과 · ∞ = 장비 효과만) + e ~ N(0, 1)(쌍 잡음).
    spike 를 주면 «장비 하나만 튐»: 무작위 장비 하나의 u = spike · 나머지 0 · 쌍 잡음 N(0, 1)."""
    if spike is not None:
        u = np.zeros(len(DEVICES))
        u[rng.integers(len(DEVICES))] = spike
    else:
        u = rng.standard_normal(len(DEVICES)) * (1.0 if np.isinf(r) else r)
    e = np.zeros(len(PAIRS)) if (spike is None and np.isinf(r)) else rng.standard_normal(len(PAIRS))
    return np.array([u[j] - u[i] + e[k] for k, (i, j) in enumerate(PAIRS)])


def fiso_sizes(s, x, keep, n_sim: int = 4000, seed: int = SEED, alpha: float = 0.05, spikes=(1.0, 3.0, 10.0),
               n_seeds: int = 1, fixed_keep: bool = False) -> dict:
    """FISO 부호 S · FISO τ 의 영 크기(설계만 · ON-Harmony 값 0 · 재점검 C-2 · FISO-1) — 영마다 한쪽 p ≤ α 인 비율.
    영 = null_d 의 r 격자(R_GRID) + «장비 하나만 튐»(spikes) · 관측 쌍 = keep · o = 쌍 차의 부호 · y = 쌍 차의 크기.
    시드 seed … seed + n_seeds − 1 에서 n_sim 씩 낸 기각률의 평균(n_seeds = 1 이면 시드 하나). fixed_keep = 고치기 전 비교용 —
    τ 의 쌍을 항등에서 고른 집합으로 모든 순열에 그대로 둔다(perm_taus_x 와 같다)."""
    s, x = np.asarray(s, float), np.asarray(x, float)
    keep = np.asarray(keep)
    idx = device_perm_index()
    SP = (s[idx] * device_perm_orient())[:, keep]                                # (순열, 관측 쌍) — S_π = SP @ o
    obs = np.zeros(len(PAIRS), bool)
    obs[keep] = True
    M = np.broadcast_to(obs & (s != 0), idx.shape) if fixed_keep else obs[None, :] & (s[idx] != 0)
    a, b = np.triu_indices(len(PAIRS), 1)
    MM = (M[:, a] & M[:, b]).astype(float)
    X = x[idx]
    SX = np.sign(X[:, a] - X[:, b]) * MM
    nx = (SX != 0).sum(1)
    nulls = [("r=" + ("inf" if np.isinf(r) else str(r)), r, None) for r in R_GRID] + [(f"spike={c}", 0.0, c) for c in spikes]
    acc = {name: [0.0, 0.0] for name, _, _ in nulls}
    for k in range(n_seeds):
        rng = np.random.default_rng(seed + k)
        for name, r, c in nulls:
            D = np.array([null_d(rng, r, c) for _ in range(n_sim)])
            S = np.sign(D[:, keep]) @ SP.T                                       # (모의, 순열)
            p_s = (S >= S[:, :1]).mean(1)
            Y = np.abs(D)
            SY = np.sign(Y[:, a] - Y[:, b])
            num = SY @ SX.T
            ny = (SY != 0).astype(float) @ MM.T
            den = np.sqrt(nx[None, :] * ny)
            T = np.where(den > 0, num / np.where(den > 0, den, 1.0), 0.0)
            p_t = (T >= T[:, :1] - 1e-12).mean(1)
            acc[name][0] += (p_s <= alpha).mean()
            acc[name][1] += (p_t <= alpha).mean()
    out = {name: dict(S=float(v[0] / n_seeds), tau=float(v[1] / n_seeds)) for name, v in acc.items()}
    return dict(n_sim=n_sim, n_seeds=n_seeds, seed=seed, alpha=alpha, se=float(np.sqrt(alpha * (1 - alpha) / (n_sim * n_seeds))),
                fixed_keep=fixed_keep, by_null=out, max_S=max(v["S"] for v in out.values()), max_tau=max(v["tau"] for v in out.values()))


def tau_b(x, y) -> float:
    """Kendall τ_b(동점 보정) — scipy.stats.kendalltau 의 기본과 같다."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    i, j = np.triu_indices(len(x), 1)
    sx, sy = np.sign(x[i] - x[j]), np.sign(y[i] - y[j])
    den = np.sqrt((sx != 0).sum() * (sy != 0).sum())
    return float((sx * sy).sum() / den) if den > 0 else 0.0


def device_perm_index() -> np.ndarray:
    """장비 순열 5,040 개 × 쌍 21 — 순열 π 아래 쌍 k = (i, j) 가 가져올 y 의 쌍 번호 (π(i), π(j))."""
    pos = {p: k for k, p in enumerate(PAIRS)}
    out = []
    for perm in itertools.permutations(range(len(DEVICES))):
        out.append([pos[tuple(sorted((perm[i], perm[j])))] for i, j in PAIRS])
    return np.array(out)


def device_perm_orient() -> np.ndarray:
    """device_perm_index 와 같은 순서 · 모양 — 순열 π 아래 쌍 k = (i, j) 의 방향: π(i) < π(j) 면 +1 · 뒤집히면 −1."""
    return np.array([[1.0 if perm[i] < perm[j] else -1.0 for i, j in PAIRS] for perm in itertools.permutations(range(len(DEVICES)))])


def perm_taus(x, y, idx) -> np.ndarray:
    """순열마다 τ_b(x, y∘π) — y 의 동점 수는 순열에 불변이라 분모는 하나."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    i, j = np.triu_indices(len(x), 1)
    sx = np.sign(x[i] - x[j])
    Y = y[idx]                                                                   # (순열, 21)
    sy = np.sign(Y[:, i] - Y[:, j])
    den = np.sqrt((sx != 0).sum() * (np.sign(y[i] - y[j]) != 0).sum())
    return (sy @ sx) / den


def interval90(x, y, idx):
    """순열 역산(위치 이동) 90 % 구간 — [τ̂ − q95(T), τ̂ − q05(T)] · [−1, 1] 로 자름."""
    T = perm_taus(x, y, idx)
    t = tau_b(x, y)
    q05, q95 = np.quantile(T, [0.05, 0.95])
    return t, max(-1.0, t - q95), min(1.0, t - q05)


def null_y(rng, r):
    u = rng.standard_normal(len(DEVICES)) * (1.0 if np.isinf(r) else r)
    e = np.zeros(len(PAIRS)) if np.isinf(r) else rng.standard_normal(len(PAIRS))
    return np.array([abs(u[i] - u[j] + e[k]) for k, (i, j) in enumerate(PAIRS)])


def design(n_sim=1000, x=None, seed=SEED, deltas=np.round(np.arange(0.0, 1.01, 0.01), 2)):
    """r 마다: 0 의 덮음률 · δ 마다 검정력 · δ_τ(검정력 ≥ 0.80 인 가장 작은 δ)."""
    x = pair_predictions() if x is None else np.asarray(x, float)
    idx = device_perm_index()
    rng = np.random.default_rng(seed)
    res = {}
    for r in R_GRID:
        lo, hi, cov = [], [], 0
        for _ in range(n_sim):
            t, a, b = interval90(x, null_y(rng, r), idx)
            lo.append(a), hi.append(b)
            cov += a <= 0 <= b
        lo, hi = np.array(lo), np.array(hi)
        power = [float(((lo > -d) & (hi < d)).mean()) for d in deltas]
        ok = [d for d, p in zip(deltas, power) if p >= 0.80]
        res["inf" if np.isinf(r) else str(r)] = dict(coverage0=cov / n_sim, delta_tau=float(ok[0]) if ok else None,
                                                     power_at=dict(zip([f"{d:.2f}" for d in deltas[::10]], power[::10])))
    dts = [v["delta_tau"] for v in res.values()]
    dt = None if any(v is None for v in dts) else max(dts)
    return dict(x=x.tolist(), pairs=[f"{DEVICES[i]}~{DEVICES[j]}" for i, j in PAIRS], n_sim=n_sim, seed=seed, by_r=res,
                delta_tau=dt, descriptive_only=(dt is None or dt > 0.5),
                rule="δ_τ = r 격자에서 가장 큰 값(보수) · δ_τ > 0.5 면 (a) 는 «기술만 — 검정력 부족»")


if __name__ == "__main__" and ("--fiso-size" in sys.argv or "--fiso-size-before" in sys.argv):
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "remeasure"))
    import analysis as AN  # noqa: E402
    import fiso_sign as FS  # noqa: E402
    import rconfig as RC  # noqa: E402
    cnt = AN.common_counts(RC.onh_sessions())
    keep = [k for k, (i, j) in enumerate(PAIRS) if cnt[frozenset((DEVICES[i], DEVICES[j]))]["judgeable"]]
    num = [int(a) for a in sys.argv[1:] if a.isdigit()]
    kw = dict(n_sim=num[0] if num else 4000, n_seeds=num[1] if len(num) > 1 else 1)
    x = magnitude_vector(FS.device_rows(cnt))
    res = (fiso_sizes(sign_vector(FS.device_rows(cnt)), x, keep, fixed_keep=True, **kw) if "--fiso-size-before" in sys.argv
           else fiso_sizes(design_signs(cnt), x, keep, **kw))
    print(json.dumps(res, ensure_ascii=False))
    sys.exit(0)

if __name__ == "__main__" and "--devices" in sys.argv:
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "remeasure"))
    import analysis as AN  # noqa: E402
    import rconfig as RC  # noqa: E402
    cnt = AN.common_counts(RC.onh_sessions())
    keep = [k for k, (i, j) in enumerate(PAIRS) if cnt[frozenset((DEVICES[i], DEVICES[j]))]["judgeable"]]
    print(json.dumps(design_subset(keep, pair_predictions_weighted(cnt), int(sys.argv[1]) if sys.argv[1:2] != ["--devices"] else 1000),
                     ensure_ascii=False))
    sys.exit(0)

if __name__ == "__main__":
    print(json.dumps(design(int(sys.argv[1]) if len(sys.argv) > 1 else 1000), ensure_ascii=False))
