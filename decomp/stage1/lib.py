"""decomp 1 단계 — 분석의 순수 함수 (SPEC_v1 §8–§10). numpy · scipy 만. 시험: tests/test_decomp_stage1.py."""
from __future__ import annotations

import math

import numpy as np
from scipy import optimize, stats

CLIP = 1e-4
TAU = math.tan(math.radians(15.0))          # §9 τ
PI_STAR = 0.5                               # §11 Π*


# ── 좌표 (§1) ─────────────────────────────────────────────────────────────
def coords(ndi, fwf, odi, ext=None, stage="C3"):
    """적합값 → (좌표 배열 (N, p), 자른 비율 (N,) bool). stage ∈ C3 · C5 · C7."""
    n = np.clip(ndi, CLIP, 1 - CLIP)
    f = np.clip(fwf, CLIP, 1 - CLIP)
    clipped = (ndi != n) | (fwf != f)
    x = np.log(n / (1 - n))
    y = np.log(f / (1 - f)) - np.log(n)
    cols = [x, y, odi]
    if stage in ("C5", "C7"):
        cols += [ext[:, 0], ext[:, 1]]
    if stage == "C7":
        cols += [ext[:, 2], ext[:, 3]]
    return np.stack(cols, -1).astype(float), clipped


STAGE_P = {"C3": 3, "C5": 5, "C7": 7}


# ── 기하 (§8) ─────────────────────────────────────────────────────────────
def inv_sqrt(S):
    """대칭 양정치 (…, p, p) → S^{-1/2}."""
    w, V = np.linalg.eigh(S)
    w = np.maximum(w, 1e-300)
    return (V * (1 / np.sqrt(w))[..., None, :]) @ np.swapaxes(V, -1, -2)


def wangle(u, v, S):
    """백색화 각(도): cos = uᵀS⁻¹v/√(uᵀS⁻¹u · vᵀS⁻¹v) — 방향 부호 무시(0–90°)."""
    Si = np.linalg.inv(S)
    c = abs(u @ Si @ v) / math.sqrt((u @ Si @ u) * (v @ Si @ v))
    return math.degrees(math.acos(min(1.0, c)))


def ls_slope(t, Y):
    """Y (T, p) 를 t (T,) 에 최소제곱 회귀한 기울기 (p,)."""
    tc = t - t.mean()
    return (tc[:, None] * (Y - Y.mean(0))).sum(0) / (tc ** 2).sum()


def representative(vecs, W12):
    """범주 대표(§9-1): 백색화(W12 = Σ^{-1/2}) · 단위화한 열의 SVD 첫 왼쪽 특이벡터와 σ2/σ1."""
    V = np.stack([W12 @ v for v in vecs], 1)
    V = V / np.linalg.norm(V, axis=0, keepdims=True)
    if V.shape[1] == 1:
        return V[:, 0], 0.0
    U, s, _ = np.linalg.svd(V, full_matrices=False)
    u = U[:, 0] * np.sign(U[:, 0] @ V.sum(1) or 1.0)
    return u, float(s[1] / s[0])


def rank_projector(Ut, tau=TAU):
    """§9-2·3: Ut (p, m) 단위 열 → (k, Q_k, P = I − Q_kQ_kᵀ, 특이값)."""
    p = Ut.shape[0]
    if Ut.shape[1] == 0:
        return 0, np.zeros((p, 0)), np.eye(p), np.zeros(0)
    U, s, _ = np.linalg.svd(Ut, full_matrices=False)
    k = int((s / s[0] >= tau).sum())
    Q = U[:, :k]
    return k, Q, np.eye(p) - Q @ Q.T, s


def span_angle(u, Others):
    """§9-5: 단위 u 와 Others (p, m) 의 열공간 사이 각(도)."""
    if Others.shape[1] == 0:
        return 90.0
    Q, _ = np.linalg.qr(Others)
    c = np.linalg.norm(Q.T @ u) / np.linalg.norm(u)
    return math.degrees(math.acos(min(1.0, c)))


def preservation(P, T):
    """§10-4: Π = σ_min(P Q_T) · Π_j · 주각(도). T (p, q) 백색화 조직 방향."""
    QT, _ = np.linalg.qr(T)
    s = np.linalg.svd(P @ QT, compute_uv=False)
    Pi = float(s.min()) if len(s) == T.shape[1] else 0.0
    Pij = np.linalg.norm(P @ T, axis=0) / np.linalg.norm(T, axis=0)
    return Pi, Pij, s


# ── 축소 추정 (§8-1) ──────────────────────────────────────────────────────
def shrink_cov(X):
    """X (…, n, p) → Schäfer–Strimmer 표적 D 꼴 축소 공분산 (…, p, p). 분산은 줄이지 않는다."""
    n = X.shape[-2]
    Xc = X - X.mean(-2, keepdims=True)
    S = np.swapaxes(Xc, -1, -2) @ Xc / (n - 1)
    sd = np.sqrt(np.maximum(np.diagonal(S, axis1=-2, axis2=-1), 1e-300))
    Z = Xc / sd[..., None, :]
    W = Z[..., :, :, None] * Z[..., :, None, :]                   # (…, n, p, p)
    Wm = W.mean(-3)
    var_r = n / (n - 1) ** 3 * ((W - Wm[..., None, :, :]) ** 2).sum(-3)
    R = Wm * n / (n - 1)
    off = ~np.eye(X.shape[-1], dtype=bool)
    num = (var_r * off).sum((-1, -2))
    den = (R ** 2 * off).sum((-1, -2))
    lam = np.clip(np.where(den > 0, num / np.maximum(den, 1e-300), 1.0), 0.0, 1.0)
    Rs = R * (1 - lam[..., None, None]) * off + np.eye(X.shape[-1])
    return Rs * sd[..., :, None] * sd[..., None, :], lam


# ── 비중심 χ² 역전 (§10-2) ────────────────────────────────────────────────
def ncx2_ci(T, df, level=0.95):
    """λ 의 등꼬리 구간. 하한 = max(0, …) · T 가 하위 (1−level)/2 미만이면 [0, 0] (공집합 없음)."""
    a = (1 - level) / 2
    f = lambda lam, q: stats.ncx2.cdf(T, df, lam) - q if lam > 0 else stats.chi2.cdf(T, df) - q
    hi_bound = max(10.0, 4 * T + 100)
    lo = 0.0 if stats.chi2.cdf(T, df) <= 1 - a else optimize.brentq(lambda l: f(l, 1 - a), 0.0, hi_bound)
    # T 가 하위 a 미만이면 상한을 정의하는 λ ≥ 0 이 없다 → [0, 0] (SPEC_v1 §10-2)
    hi = 0.0 if stats.chi2.cdf(T, df) <= a else optimize.brentq(lambda l: f(l, a), 0.0, hi_bound)
    return lo, hi


# ── 지수 (§10) ────────────────────────────────────────────────────────────
def group_basis(groups, W12, tau=TAU):
    """§9-1·2: 범주(열 묶음, 원 좌표) → 백색화 대표 벡터 열 (p, m). 범주 안 σ2/σ1 ≥ τ 이면 둘째 벡터도."""
    cols = []
    for G in groups:
        V = W12 @ G
        V = V / np.linalg.norm(V, axis=-2, keepdims=True)
        if V.shape[-1] == 1:
            cols.append(V)
            continue
        Uu, s, _ = np.linalg.svd(V, full_matrices=False)
        cols.append(Uu[..., :, :1])
        sec = (s[..., 1] / s[..., 0] >= tau)
        cols.append(Uu[..., :, 1:2] * np.asarray(sec)[..., None, None])
    if not cols:
        return np.zeros(W12.shape[:-1] + (0,))
    return np.concatenate(cols, -1)


def proj_from_basis(Ub, tau=TAU):
    """(…, p, m) 단위/0 열 → 사영 행렬 P (…, p, p) 와 계수 k (…)."""
    p = Ub.shape[-2]
    if Ub.shape[-1] == 0:
        return np.broadcast_to(np.eye(p), Ub.shape[:-2] + (p, p)).copy(), np.zeros(Ub.shape[:-2], int)
    Uu, s, _ = np.linalg.svd(Ub, full_matrices=False)
    keep = s / np.maximum(s[..., :1], 1e-300) >= tau
    Q = Uu * keep[..., None, :]
    return np.eye(p) - Q @ np.swapaxes(Q, -1, -2), keep.sum(-1)


def index_stats(D, U, rng, B=499, level=0.95, tau=TAU, groups=None, basic_ci=False):
    """한 쌍 · 한 ROI: D (n, p) 피험자 차이, U (p, m) 교란 방향(원 좌표 단위 · 열).
    반환 dict: k, Psi, Den, I(또는 None), p_value, Psi_ci, den_lower, reported, lam."""
    n, p = D.shape
    W, lam = shrink_cov(D)
    W12 = inv_sqrt(W)
    if groups is None:
        groups = [U[:, i:i + 1] for i in range(U.shape[1])]
    P, k = proj_from_basis(group_basis(groups, W12, tau), tau)
    k = int(k)
    Q = np.linalg.svd(np.eye(p) - P)[0][:, :k]
    if k >= p:
        return dict(k=k, exhausted=True, reported=False)
    Delta = D.mean(0)
    dt = W12 @ Delta
    r = P @ dt
    T = n * float(r @ r)
    Ttot = n * float(dt @ dt)
    # 매개변수 부트스트랩(§10-2): 귀무 평균 μ̂₀ = W^{1/2} Q Qᵀ W^{-1/2} Δ, 잡음 e* ~ N(0, W) · 같은 축소 추정
    W_half = np.linalg.cholesky(W)
    mu0 = np.linalg.inv(W12) @ (Q @ (Q.T @ dt))
    E = rng.standard_normal((B, n, p)) @ W_half.T
    Wb, _ = shrink_cov(E + mu0)               # 공분산은 평균 이동에 불변 — 분자 · 분모 자료가 같은 W*
    Wb12 = inv_sqrt(Wb)
    Eb = E.mean(1)                                                     # (B, p)
    mnum = np.einsum("bij,bj->bi", Wb12, mu0 + Eb)
    mden = np.einsum("bij,bj->bi", Wb12, Eb)
    Pb, _ = proj_from_basis(group_basis(groups, Wb12, tau), tau)     # §10-2 3: P*_b 를 매 b 마다 다시
    rnum = np.einsum("bij,bj->bi", Pb, mnum)
    num_b = (rnum ** 2).sum(1)
    den_b = (mden ** 2).sum(1)
    E0, D0 = float(num_b.mean()), float(den_b.mean())
    Psi = float(r @ r) - E0
    Den = float(dt @ dt) - D0
    T_b = n * (num_b)
    pval = (1 + (T_b >= T).sum()) / (B + 1)
    Tadj = T * (p - k) / max(n * E0, 1e-300)
    lo, hi = ncx2_ci(Tadj, p - k, level)
    Ttot_adj = Ttot * p / max(n * D0, 1e-300)
    den_ok = Ttot_adj > stats.chi2.ppf(1 - (1 - level) / 2, p)
    out_basic = None
    if basic_ci:                                       # SPEC_v2 §5 — 매개변수 부트스트랩 기본 구간(대립 적합: 평균 Δ̂)
        E2 = rng.standard_normal((B, n, p)) @ W_half.T + Delta
        Wa, _ = shrink_cov(E2)
        Wa12 = inv_sqrt(Wa)
        Pa, _ = proj_from_basis(group_basis(groups, Wa12, tau), tau)
        ra = np.einsum("bij,bj->bi", Pa, np.einsum("bij,bj->bi", Wa12, E2.mean(1)))
        psi_b = (ra ** 2).sum(1) - E0
        lo_b, hi_b = np.quantile(psi_b, [(1 - level) / 2, 1 - (1 - level) / 2])
        out_basic = (max(0.0, 2 * Psi - hi_b), max(0.0, 2 * Psi - lo_b))
    return dict(k=k, exhausted=False, Psi=Psi, Den=Den, I=(Psi / Den) if den_ok else None, p_value=float(pval), Psi_ci_basic=out_basic,
                Psi_ci=(lo / n, hi / n), reported=bool(den_ok), lam=float(lam), E0=E0, D0=D0, T=T, Ttot=Ttot)
