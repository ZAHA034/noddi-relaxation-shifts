"""M6 — a second TE within each site (plan PLAN_M6.md); the first observation is made identical to model.generate by drawing the same random numbers in the same order, and the second observation (TE_s + ΔTE) shares the person · site · ROI terms and differs only in the relaxation term and the measurement noise."""
import math

import numpy as np
from scipy import stats

from design_sim import analyses as A
from design_sim import model as M

DTE = 30.0                        # ms — difference between the second TE and the site TE (plan Section 1)
T2_REF = 70.0                     # ms — approximate tissue T2 of the multiplier exp(ΔTE / T2_REF) by which the noise grows with the signal loss at the longer TE (assumed)
NOISE2 = math.exp(DTE / T2_REF)


def generate2(cell: dict, rng: np.random.Generator, rng2: np.random.Generator, dte: float = DTE, noise2: float = NOISE2) -> dict:
    """First observation y1 (= the y of model.generate) · second observation y2 · site · g · te · te_s · θ."""
    site, g = M.layout(cell["N"])
    n = len(site)
    te_s = M.site_te(cell["D"])
    te = te_s[site]
    th = M.theta_of(cell["S"])
    asg = M.SCEN[cell["S"]][1]
    fin = M.LENGTHEN if asg in ("b", "i") else 1.0
    fen = M.LENGTHEN if asg in ("b", "e") else 1.0
    t2in = np.where(g == 1, M.T2_IN_C * fin, M.T2_IN_C)
    t2en = np.where(g == 1, M.T2_EN_C * fen, M.T2_EN_C)
    v = cell["V"]
    if v > 0:
        s = math.sqrt(math.log(1.0 + v * v))
        z1 = rng.standard_normal(n)
        z2 = M.CORR_V * z1 + math.sqrt(1 - M.CORR_V ** 2) * rng.standard_normal(n)
        t2in = t2in * np.exp(s * z1 - s * s / 2)
        t2en = t2en * np.exp(s * z2 - s * s / 2)
    dr2 = 1.0 / t2en - 1.0 / t2in
    base = np.tile(np.log(M.NU_R / (1 - M.NU_R)), (n, 1)) + th * g[:, None]
    b = M.SD_B * (math.sqrt(M.RHO_ROI) * rng.standard_normal((n, 1)) + math.sqrt(1 - M.RHO_ROI) * rng.standard_normal((n, M.N_ROI)))
    gam = np.array([M.VENDOR_MEAN[M.VENDOR[a]] for a in M.ACQ])[:, None] + M.SD_GAMMA * rng.standard_normal((len(M.ACQ), M.N_ROI))
    dl = np.exp(M.SD_LOGDELTA * rng.standard_normal((len(M.ACQ), M.N_ROI)))
    e1 = rng.standard_normal((n, M.N_ROI))
    h = cell["H"]
    if h == "H1":
        eta = M.SD_ETA * rng.standard_normal(len(M.ACQ))
    elif h == "H2":
        eta = np.array([M.ETA_SIEMENS if M.VENDOR[a] == "Siemens" else 0.0 for a in M.ACQ])
    else:
        eta = np.zeros(len(M.ACQ))
    # same order of additions as model.generate (so that the floating-point results are equal)
    y1 = base + (te * dr2)[:, None]
    y1 = y1 + b
    y1 = y1 + gam[site]
    y1 = y1 + dl[site] * M.SD_EPS * e1
    y1 = y1 + (eta[site] * g)[:, None]
    e2 = rng2.standard_normal((n, M.N_ROI))
    y2 = base + ((te + dte) * dr2)[:, None]
    y2 = y2 + b
    y2 = y2 + gam[site]
    y2 = y2 + dl[site] * M.SD_EPS * noise2 * e2
    y2 = y2 + (eta[site] * g)[:, None]
    rng.standard_normal(n)                # the two sets of measured-relaxation errors (M5b) that model.generate draws at the end of a replicate — the values are not used, only the random-number stream is matched
    rng.standard_normal(n)                # (fixed 2026-11-18 · record x15/r42/RECORD_m6_rng_2026-11-18.md, a development record not in this repository)
    return dict(y1=y1, y2=y2, site=site, g=g, te=te, te_s=te_s, theta=th)


def m6a(y1, y2, g, site, k, te, dte=DTE):
    """All participants: for each person the extrapolation to TE → 0, y0 = y1 − TE·(y2 − y1)/ΔTE, with a fixed site term (M2)."""
    if dte == 0:
        raise ValueError("둘째 TE 가 사이트 TE 와 같다 — 기울기를 낼 수 없다")
    y0 = y1 - te[:, None] * (y2 - y1) / dte
    return A.m2(y0, g, site, k)


def subset_mask(site, g):
    """First half in numbering order within each site · group (plan Section 2)."""
    m = np.zeros(len(site), bool)
    for s in np.unique(site):
        for gg in (0.0, 1.0):
            idx = np.flatnonzero((site == s) & (g == gg))
            m[idx[: len(idx) // 2]] = True
    return m


def m6b(y1, y2, g, site, k, te_s, dte=DTE):
    """Half only: for each site θ̂_s = [group difference of everyone's y1] − TE_s·[group difference of the TE slopes of the half], with variance for each group s²(y1)/n + TE²·s²(slope)/m − 2·TE·cov(y1, slope)/n (the half is a subset of everyone); sites combined with inverse-variance weights · normal approximation (df = number of people − 2 × sites)."""
    if dte == 0:
        raise ValueError("둘째 TE 가 사이트 TE 와 같다 — 기울기를 낼 수 없다")
    sub = subset_mask(site, g)
    sl = (y2 - y1) / dte
    num, den = 0.0, 0.0
    for s in range(k):
        T = te_s[s]
        U, V = [], []
        for gg in (1.0, 0.0):
            a = (site == s) & (g == gg)
            bsub = a & sub
            n_, m_ = a.sum(), bsub.sum()
            u = y1[a].mean(0) - T * sl[bsub].mean(0)
            yb, sb = y1[bsub], sl[bsub]
            c = ((yb - yb.mean(0)) * (sb - sb.mean(0))).sum(0) / (m_ - 1)
            var = y1[a].var(0, ddof=1) / n_ + T ** 2 * sl[bsub].var(0, ddof=1) / m_ - 2 * T * c / n_
            U.append(u)
            V.append(var)
        th_s, v_s = U[0] - U[1], V[0] + V[1]
        w = 1.0 / v_s
        num = num + w * th_s
        den = den + w
    est, se = num / den, np.sqrt(1.0 / den)
    return dict(est=np.asarray(est, float), se=np.asarray(se, float), df=float(len(site) - 2 * k))


def one(cell: dict, rng: np.random.Generator, rng2: np.random.Generator) -> dict:
    d = generate2(cell, rng, rng2)
    y1, y2, g, site, te, te_s = d["y1"], d["y2"], d["g"], d["site"], d["te"], d["te_s"]
    k = len(M.ACQ)
    out = {"M6a": m6a(y1, y2, g, site, k, te), "M6b": m6b(y1, y2, g, site, k, te_s), "M2": A.m2(y1, g, site, k)}
    if float(np.ptp(te_s)) > 0:
        out["M3_te0"] = A.m3(y1, g, site, k, te, M.TE_C)["te0"]
    return {key: (e["est"], e["se"], A.covered(e, d["theta"]), A.rejected(e), d["theta"]) for key, e in out.items()}


def summarize(recs: list[dict]) -> dict:
    R = len(recs)
    out = {}
    for key in recs[0]:
        est = np.array([r[key][0] for r in recs])
        se = np.array([r[key][1] for r in recs])
        cov = np.array([r[key][2] for r in recs], float).mean(1)
        rej = np.array([r[key][3] for r in recs], float).mean(1)
        err = (est - np.array([r[key][4] for r in recs])[:, None]).mean(1)
        out[key] = dict(bias=float(err.mean()), bias_mcse=float(err.std(ddof=1) / math.sqrt(R)), mean_se=float(se.mean()),
                        coverage=float(cov.mean()), coverage_mcse=float(cov.std(ddof=1) / math.sqrt(R)),
                        reject=float(rej.mean()), reject_mcse=float(rej.std(ddof=1) / math.sqrt(R)))
    for key in ("M6a", "M6b", "M3_te0"):
        if key in out:
            out[key]["se_ratio_to_M2"] = out[key]["mean_se"] / out["M2"]["mean_se"]
    return out


def cells_m6() -> list[dict]:
    """The 16 target cells: D1 × eight scenarios × V (0, 0.10) × H0 × balanced (original cell numbers kept)."""
    return [c for c in M.cells() if c["D"] == "D1" and c["H"] == "H0" and c["N"] == "bal"]


def t_crit(df):
    return stats.t.ppf(1 - A.ALPHA / 2, df)
