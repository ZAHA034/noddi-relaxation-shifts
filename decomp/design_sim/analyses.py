"""design_sim.analyses — analyses M0 – M5 (ANALYSIS_PLAN.md Section 4), separately for each ROI unit, with 95 % intervals of the group coefficient from the t distribution; ComBat uses rules_env of the repository (the neuroCombat 0.2.12 wrapper `_combat` and `combat_apply`, the same formula as neuroCombatFromTraining)."""
from __future__ import annotations

import os
import pathlib
import sys

import numpy as np
from scipy import stats

_HERE = pathlib.Path(__file__).resolve().parent
for _cand in (_HERE.parent, pathlib.Path(os.environ.get("DESIGN_SIM_REPO", "")) / "decomp"):
    if (_cand / "rules_env").is_dir():
        if str(_cand) not in sys.path:
            sys.path.insert(0, str(_cand))
        break
from rules_env import analyses as RA  # noqa: E402

ALPHA = 0.05


def ols(X: np.ndarray, Y: np.ndarray):
    """Several responses (columns) at once: coefficients (p × R), residual variance (R,), (X'X)⁻¹, degrees of freedom."""
    XtX_inv = np.linalg.inv(X.T @ X)
    B = XtX_inv @ (X.T @ Y)
    res = Y - X @ B
    df = X.shape[0] - X.shape[1]
    s2 = (res * res).sum(0) / df
    return B, s2, XtX_inv, df


def coef(B, s2, inv, j):
    return B[j], np.sqrt(s2 * inv[j, j])


def site_dummies(site: np.ndarray, k: int) -> np.ndarray:
    return (site[:, None] == np.arange(k)[None, :]).astype(float)


def _est(est, se, df):
    return dict(est=np.asarray(est, float), se=np.asarray(se, float), df=float(df))


def m0(y, g):
    X = np.column_stack([np.ones_like(g), g])
    B, s2, inv, df = ols(X, y)
    return _est(*coef(B, s2, inv, 1), df)


def m2(y, g, site, k):
    X = np.column_stack([site_dummies(site, k), g])
    B, s2, inv, df = ols(X, y)
    return _est(*coef(B, s2, inv, k), df)


def m3(y, g, site, k, te, te_c):
    """θ̂(TE→0) = group coefficient, θ̂(TE_c), test of group × TE (slope)."""
    X = np.column_stack([site_dummies(site, k), g, g * te])
    B, s2, inv, df = ols(X, y)
    jg, jt = k, k + 1
    e0, se0 = coef(B, s2, inv, jg)
    sl, sesl = coef(B, s2, inv, jt)
    ec = B[jg] + te_c * B[jt]
    sec = np.sqrt(s2 * (inv[jg, jg] + 2 * te_c * inv[jg, jt] + te_c ** 2 * inv[jt, jt]))
    return dict(te0=_est(e0, se0, df), tec=_est(ec, sec, df), slope=_est(sl, sesl, df))


def harmonized_ols(H, g):
    return m0(H, g)


def m1(y, g, site, raw=False):
    """Standard ComBat (EB, group kept) → OLS ~ group. With raw=True, on raw-scale NDI (M1r)."""
    dat = (1.0 / (1.0 + np.exp(-y))) if raw else y
    H = RA._combat(dat.T, site, [int(x) for x in g], mean_only=False)["data"].T
    return harmonized_ols(H, g)


def m1c(y, g, site):
    """Trained on controls only (EB, no covariate) → applied in the same way to all participants → OLS ~ group."""
    c = g == 0
    tr = RA._combat(y[c].T, site[c], None, mean_only=False)
    H = RA.combat_apply(y.T, site, tr["estimates"]).T
    return harmonized_ols(H, g)


def m4(y, g, site, k, te_s):
    """Random-effects meta-regression in 2 stages: site-wise differences d_s on TE_s — τ² by the method of moments (regression version of DerSimonian–Laird), interval from t(k − 2)."""
    R = y.shape[1]
    d = np.empty((k, R))
    v = np.empty((k, R))
    for s in range(k):
        m = site == s
        yp, yc = y[m & (g == 1)], y[m & (g == 0)]
        npt, nc = len(yp), len(yc)
        sp2 = (((yp - yp.mean(0)) ** 2).sum(0) + ((yc - yc.mean(0)) ** 2).sum(0)) / (npt + nc - 2)
        d[s] = yp.mean(0) - yc.mean(0)
        v[s] = sp2 * (1.0 / npt + 1.0 / nc)
    X = np.column_stack([np.ones(k), te_s])
    p = X.shape[1]
    w = 1.0 / v                                                    # (k, R)
    XtWX = np.einsum("ki,kr,kj->rij", X, w, X)
    XtWd = np.einsum("ki,kr,kr->ri", X, w, d)
    bfe = np.linalg.solve(XtWX, XtWd[..., None])[..., 0]           # (R, p)
    fit = (X @ bfe.T)                                               # (k, R)
    Q = (w * (d - fit) ** 2).sum(0)
    inv = np.linalg.inv(XtWX)
    XtW2X = np.einsum("ki,kr,kj->rij", X, w * w, X)
    c = w.sum(0) - np.einsum("rij,rji->r", inv, XtW2X)
    tau2 = np.maximum(0.0, (Q - (k - p)) / c)
    ws = 1.0 / (v + tau2[None, :])
    XtWX2 = np.einsum("ki,kr,kj->rij", X, ws, X)
    XtWd2 = np.einsum("ki,kr,kr->ri", X, ws, d)
    cov = np.linalg.inv(XtWX2)
    b = np.einsum("rij,rj->ri", cov, XtWd2)
    df = k - p
    return dict(te0=_est(b[:, 0], np.sqrt(cov[:, 0, 0]), df), slope=_est(b[:, 1], np.sqrt(cov[:, 1, 1]), df), tau0=(tau2 == 0.0))


def m5(y, g, site, k, te, r2_meas):
    """Measured relaxation: TE_s·(R̂2,en − R̂2,in) subtracted for each person, then M2."""
    m_in, m_en = r2_meas
    return m2(y - (te * (m_en - m_in))[:, None], g, site, k)


def covered(e: dict, target) -> np.ndarray:
    t = stats.t.ppf(1 - ALPHA / 2, e["df"])
    return np.abs(e["est"] - target) <= t * e["se"]


def rejected(e: dict) -> np.ndarray:
    t = stats.t.ppf(1 - ALPHA / 2, e["df"])
    return np.abs(e["est"]) > t * e["se"]
