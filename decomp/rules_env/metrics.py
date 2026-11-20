"""Computes ORS, the offset-subtracted application of the retest rule of manuscript 3.4 (ROI unit difference |Δ| ≤ c·√2·s and share ≥ 0.90), and the quantities derived from it (TE tolerance, τ̂, scale models, proportional bias, rank flips, S2 ratio), without any new threshold."""
from __future__ import annotations

import math

import numpy as np
from scipy import stats
from scipy.optimize import brentq, minimize_scalar


def _ok(v) -> bool:
    return v is not None and math.isfinite(v)


def _need(n: int, ratio: float) -> int:
    return math.ceil(ratio * n - 1e-9)                    # Same integer form as analysis.retest_level


def ors(d: dict, s: dict, c: float, ratio: float = 0.90, loo: bool = True, min_units: int = 10, max_excluded: float = 0.10) -> dict:
    """Offset-subtracted retest ratio. loo=True: for each unit, subtract the median of d of the other units and judge · loo=False: no offset (the original rule)."""
    units = [u for u in d if _ok(d[u]) and _ok(s.get(u)) and s[u] > 0]
    n_exc = len(d) - len(units)
    v = np.array([d[u] for u in units], float)
    sv = np.array([s[u] for u in units], float)
    n = len(units)
    if n == 0:
        return dict(share=float("nan"), n=0, n_excluded=n_exc, offset=float("nan"), retest_level=None, pass_units={})
    if loo:
        off = np.array([np.median(np.delete(v, i)) for i in range(n)]) if n > 1 else np.zeros(1)
    else:
        off = np.zeros(n)
    hit = np.abs(v - off) <= c * math.sqrt(2.0) * sv
    level = None if (n < min_units or n_exc / max(len(d), 1) > max_excluded) else bool(hit.sum() >= _need(n, ratio))
    return dict(share=float(hit.mean()), n=n, n_excluded=n_exc, offset=float(np.median(v)) if loo else 0.0,
                retest_level=level, pass_units={u: bool(h) for u, h in zip(units, hit)})


def ors_diff(d1, s1, d2, s2, c, clusters: dict, B: int = 2000, seed: int = 0, loo: bool = True) -> dict:
    """Difference of ORS between two scales on the same units (1 − 2) and a paired bootstrap 95 % interval — if there are two or more bundles (participants), the bundles are redrawn, otherwise the units."""
    p1 = ors(d1, s1, c, loo=loo)["pass_units"]
    p2 = ors(d2, s2, c, loo=loo)["pass_units"]
    common = sorted(set(p1) & set(p2))
    diff_u = np.array([float(p1[u]) - float(p2[u]) for u in common])
    point = float(diff_u.mean()) if len(common) else float("nan")
    rng = np.random.default_rng(seed)
    cl = np.array([clusters.get(u, u) for u in common], dtype=object)
    groups = sorted(set(cl.tolist()), key=str)
    reps = np.empty(B)
    for b in range(B):
        if len(groups) >= 2:
            pick = rng.choice(len(groups), len(groups), replace=True)
            idx = np.concatenate([np.flatnonzero(cl == groups[g]) for g in pick])
        else:
            idx = rng.integers(0, len(common), len(common))
        reps[b] = diff_u[idx].mean()
    lo, hi = np.percentile(reps, [2.5, 97.5])
    return dict(diff=point, lo=float(lo), hi=float(hi), n=len(common))


def _infl(n, kind: str = "mean"):
    """Variance factor of the offset (the other n − 1 units) — mean: 1 + 1/(n − 1) (the formula of plan §1) · median: 1 + π/(2(n − 1)) (median approximation · sensitivity)."""
    n = np.asarray(n, float)
    k = 1.0 if kind == "mean" else math.pi / 2.0
    return np.where(n > 1, 1.0 + k / np.maximum(n - 1.0, 1.0), 1.0)


def _te_share(dte: float, s_u: np.ndarray, tau: float, c: float, n_group=None, infl_kind: str = "mean") -> float:
    """n_group: the number of units of the bundle for which an offset is computed per unit (a scalar or the same length as s_u) — if absent, len(s_u) (a single bundle). Code review F2 of 2026-11-15: for an array that pools several bundles, n_group must be given (ANALYSIS_PLAN §8)."""
    s_u = np.asarray(s_u, float)
    infl = _infl(len(s_u) if n_group is None else n_group, infl_kind)
    return float(np.mean(2.0 * stats.norm.cdf(c * math.sqrt(2.0) * s_u / np.sqrt(2.0 * s_u ** 2 * infl + (dte * tau) ** 2)) - 1.0))


def te_tolerance(s_u, tau: float, c: float = 1.96, ratio: float = 0.90, span: float = 70.0, n_group=None, infl_kind: str = "mean") -> float:
    """ΔTE (ms) at which the expected ORS equals ratio for offset-subtracted differences ~ N(0, 2s²(1 + 1/(n−1)) + ΔTE²τ²). inf if τ ≤ 0 or the root lies outside the span."""
    s_u = np.asarray(s_u, float)
    if not tau > 0:
        return math.inf
    f = lambda x: _te_share(x, s_u, tau, c, n_group, infl_kind) - ratio
    if f(0.0) < 0:
        return 0.0
    if f(span) > 0:
        return math.inf
    return float(brentq(f, 0.0, span, xtol=1e-6))


def tau_hat(k_hat, var_noise: float, ddof: int = 1) -> float:
    """SD of the TE slope between units after subtracting the noise variance. ddof: if the values with the bundle mean subtracted are pooled, the number of bundles (pooled variance · sensitivity F3)."""
    return math.sqrt(max(float(np.var(np.asarray(k_hat, float), ddof=ddof)) - var_noise, 0.0))


def beta_level(k_hat, xbar, var_xbar: float) -> float:
    """The share by which the TE slope varies with the level, cov(k̂, x̄)/(var(x̄) − noise variance)."""
    k, x = np.asarray(k_hat, float), np.asarray(xbar, float)
    return float(np.cov(k, x, ddof=1)[0, 1] / (np.var(x, ddof=1) - var_xbar))


def te_tolerance_interval(s_u, k_hat, var_k: float, df: int, clusters, B: int = 2000, seed: int = 0,
                          c: float = 1.96, ratio: float = 0.90, span: float = 70.0, n_group=None, infl_kind: str = "mean",
                          tau_ddof: int = 1) -> tuple[float, float, float]:
    """Point estimate and 90 % interval of ΔTE* — τ̂ by bundle bootstrap, and s is drawn together with the distribution of the sample SD (√(χ²_df/df) times) (var_k is also proportional to s²)."""
    s_u, k = np.asarray(s_u, float), np.asarray(k_hat, float)
    cl = np.asarray(clusters)
    groups = np.unique(cl)
    ng = None if n_group is None else np.broadcast_to(np.asarray(n_group, float), s_u.shape)
    point = te_tolerance(s_u, tau_hat(k, var_k, tau_ddof), c, ratio, span, ng, infl_kind)
    rng = np.random.default_rng(seed)
    reps = np.empty(B)
    for b in range(B):
        m = math.sqrt(rng.chisquare(df) / df)
        if len(groups) >= 2:
            idx = np.concatenate([np.flatnonzero(cl == g) for g in rng.choice(groups, len(groups), replace=True)])
        else:
            idx = rng.integers(0, len(k), len(k))
        reps[b] = te_tolerance(s_u * m, tau_hat(k[idx], var_k * m * m, tau_ddof), c, ratio, span, ng, infl_kind)
    reps.sort()
    lo, hi = reps[int(math.floor(0.05 * (B - 1)))], reps[int(math.ceil(0.95 * (B - 1)))]
    return point, float(lo), float(hi)


def raw_shift_pred(nu_a, c: float):
    """The shift expit(logit ν + c) − ν that the logit shift c produces on the raw scale."""
    nu = np.asarray(nu_a, float)
    return 1.0 / (1.0 + np.exp(-(np.log(nu / (1 - nu)) + c))) - nu


def scale_models(nu_a, d_raw) -> dict:
    """Two models with one parameter each: «constant raw shift» (d = a) versus «constant logit shift» (d = expit(logit ν_A + c) − ν_A)."""
    nu, d = np.asarray(nu_a, float), np.asarray(d_raw, float)
    a = float(d.mean())
    res = minimize_scalar(lambda c: float(np.sum((d - raw_shift_pred(nu, c)) ** 2)), bounds=(-5.0, 5.0), method="bounded",
                          options=dict(xatol=1e-10))
    r_raw, r_log = d - a, d - raw_shift_pred(nu, res.x)
    return dict(a=a, c=float(res.x), ss_raw=float(np.sum(r_raw ** 2)), ss_logit=float(np.sum(r_log ** 2)),
                resid_raw=r_raw, resid_logit=r_log)


def prop_bias(a, b) -> float:
    """Proportional bias (Bland & Altman 1999): the slope of an OLS of (b − a) on (a + b)/2."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.polyfit((a + b) / 2.0, b - a, 1)[0])


def spearman(a, b) -> float:
    return float(stats.spearmanr(a, b)[0])


def rank_flip_share(a, b, rc: float) -> float:
    """Among pairs of units that are farther apart than rc in A, the share whose order is reversed in B."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    i, j = np.triu_indices(len(a), 1)
    far = np.abs(a[i] - a[j]) > rc
    if not far.any():
        return float("nan")
    return float(np.mean(np.sign(a[i][far] - a[j][far]) != np.sign(b[i][far] - b[j][far])))


def s2_ratio(fit: float, alg: float) -> float:
    return float(fit / alg)


def s2_reading(ratios, lo: float = 0.80, hi: float = 1.25) -> bool:
    return bool(all(lo <= r <= hi for r in ratios))
