"""age_te — within-person comparison of the covariate (age · sex) × scanner TE interaction (plan decomp/age_te/ANALYSIS_PLAN.md): model v_ia = α_i + μ_a + γ · (TE_a − TE_C) · (z_i − z̄) + ε_ia and its inference, the falsification controls, the E3 tools and the rule on whether to run; reads no data."""
import itertools
import math

import numpy as np
from scipy import stats

TE_ACQ = {"NOT1ACH": 70.6, "NOT2ING-89": 89.0, "NOT2ING-98": 98.0, "NOT4GEP": 68.0, "OXF4GEP": 69.65,
          "OXF1PRI": 92.0, "OXF2PRI": 92.0, "OXF3TRI": 96.4}          # recorded TE (middle value if a range) — same as Supplement S-M9.7
TE_C = sum(TE_ACQ.values()) / len(TE_ACQ)                              # 84.45625 ms
VENDOR = {"NOT1ACH": "Philips", "NOT2ING-89": "Philips", "NOT2ING-98": "Philips", "NOT4GEP": "GE", "OXF4GEP": "GE",
          "OXF1PRI": "Siemens", "OXF2PRI": "Siemens", "OXF3TRI": "Siemens"}
T2EN = 52.0                                                             # control extra-neurite T2 of 2.4 (ms) — only for the conversion used in interpretation
T2IN = 82.5                                                             # control intra-neurite T2 of 2.4 (ms) — only for the comparison with the assumed value (E3)


def design(person, acq, z, extra=None, te_map=None):
    """Design matrix: interaction column (first column) · person dummies · acquisition dummies (first acquisition left out) · extra columns. te_map is changed only in the falsification control (TE permutation)."""
    te_map = TE_ACQ if te_map is None else te_map
    person, acq = np.asarray(person), np.asarray(acq)
    zs = np.asarray([z[p] for p in person], float)
    zbar = float(np.mean(list(z.values())))
    te = np.asarray([te_map[a] for a in acq], float)
    cols = [(te - TE_C) * (zs - zbar)]
    ps, as_ = sorted(set(person)), sorted(set(acq))
    cols += [(person == p).astype(float) for p in ps]
    cols += [(acq == a).astype(float) for a in as_[1:]]
    if extra is not None:
        cols += list(extra)
    return np.column_stack(cols)


def fit(v, person, acq, z, extra=None, te_map=None):
    """OLS — γ (coefficient of the first column) · SE · degrees of freedom · 95 % interval · residual SD · residual sum of squares."""
    X = design(person, acq, z, extra, te_map)
    v = np.asarray(v, float)
    beta, *_ = np.linalg.lstsq(X, v, rcond=None)
    resid = v - X @ beta
    rank = np.linalg.matrix_rank(X)
    df = len(v) - rank
    s2 = float(resid @ resid) / df
    cov = s2 * np.linalg.pinv(X.T @ X)
    se = math.sqrt(cov[0, 0])
    t = stats.t.ppf(0.975, df)
    return dict(gamma=float(beta[0]), se=se, df=int(df), ci95=[float(beta[0] - t * se), float(beta[0] + t * se)],
                resid_sd=math.sqrt(s2), rss=float(resid @ resid), n_obs=int(len(v)), n_person=len(set(person)), n_acq=len(set(acq)))


def permutation_p(v, person, acq, z, n_perm=10000, seed=20261118, extra=None):
    """Two-sided p from the null distribution with the covariate shuffled between people — (1 + #|γ*| ≥ |γ|)/(n + 1)."""
    g0 = abs(fit(v, person, acq, z, extra)["gamma"])
    rng = np.random.default_rng(seed)
    keys = sorted(z)
    vals = np.array([z[k] for k in keys], float)
    hit = 0
    for _ in range(n_perm):
        zp = dict(zip(keys, rng.permutation(vals)))
        hit += abs(fit(v, person, acq, zp, extra)["gamma"]) >= g0 - 1e-15
    return (1 + hit) / (n_perm + 1)


def mde(se, df, power=0.80, alpha=0.05):
    """Minimum detectable effect (two-sided) — (t_{1−α/2} + t_{power}) · SE."""
    return float((stats.t.ppf(1 - alpha / 2, df) + stats.t.ppf(power, df)) * se)


def t2en_change(gamma_per_unit, t2en=T2EN):
    """Change of extra-neurite T2 (ms · %) when γ (logit / ms / covariate unit) is read as «intra-neurite T2 unchanged, only extra-neurite T2 changes» — for interpretation. If R2,en changes by γ, T2,en′ = 1/(1/T2,en + γ)."""
    t2p = 1.0 / (1.0 / t2en + gamma_per_unit)
    return dict(ms=float(t2p - t2en), pct=float((t2p / t2en - 1) * 100))


def slope_heterogeneity(slopes, ses):
    """TE slopes of each person and their SE → homogeneity χ² (around the weighted mean) · degrees of freedom · p."""
    s, e = np.asarray(slopes, float), np.asarray(ses, float)
    w = 1 / e ** 2
    m = float(np.sum(w * s) / np.sum(w))
    chi2 = float(np.sum(w * (s - m) ** 2))
    df = len(s) - 1
    return dict(weighted_mean=m, chi2=chi2, df=df, p=float(stats.chi2.sf(chi2, df)))


def fit_units(units, person, acq, z, coord=0):
    """Sensitivity (plan Section 5) — unit-level OLS: observation = session × unit, fixed effects = person × unit · acquisition × unit. units[k] = {unit: (x, y, ODI)}. The SE ignores the correlation between units of the same session, so it is used descriptively only."""
    pu, au, v, zs, te = [], [], [], [], []
    for us, p, a in zip(units, person, acq):
        for u, c in us.items():
            pu.append(f"{p}|{u}")
            au.append(f"{a}|{u}")
            v.append(float(c[coord]))
            zs.append(z[p])
            te.append(TE_ACQ[a])
    pu, au = np.asarray(pu), np.asarray(au)
    zbar = float(np.mean(list(z.values())))
    cols = [(np.asarray(te) - TE_C) * (np.asarray(zs) - zbar)]
    cols += [(pu == k).astype(float) for k in sorted(set(pu))]
    cols += [(au == k).astype(float) for k in sorted(set(au))]
    X = np.column_stack(cols)
    v = np.asarray(v, float)
    beta, *_ = np.linalg.lstsq(X, v, rcond=None)
    resid = v - X @ beta
    df = len(v) - np.linalg.matrix_rank(X)
    s2 = float(resid @ resid) / df
    # variance of γ: inverse of the residual sum of squares of the first column regressed on the other columns × s² (cheaper than a full pinv with many fixed effects, and equal)
    Z = X[:, 1:]
    r0 = X[:, 0] - Z @ np.linalg.lstsq(Z, X[:, 0], rcond=None)[0]
    se = math.sqrt(s2 / float(r0 @ r0))
    t = stats.t.ppf(0.975, df)
    return dict(gamma=float(beta[0]), se=se, df=int(df), ci95=[float(beta[0] - t * se), float(beta[0] + t * se)], n_obs=int(len(v)),
                note="unit-level OLS; SE ignores correlation between units of one session (descriptive)")


def te_slope(tes, v):
    """E3 — line of the value v on TE (ms) in one person's TE series → slope · SE (residual) · degrees of freedom. Slope = R2,en − R2,in (signal model, for x)."""
    t = np.asarray(tes, float)
    v = np.asarray(v, float)
    A = np.column_stack([np.ones_like(t), t - t.mean()])
    beta, *_ = np.linalg.lstsq(A, v, rcond=None)
    resid = v - A @ beta
    df = len(v) - 2
    se = math.sqrt(float(resid @ resid) / df / float(np.sum((t - t.mean()) ** 2)))
    return dict(slope=float(beta[1]), se=se, df=int(df))


def assumed_slope_sd(cv=0.10, rho=0.8, n=200000, seed=20261118):
    """Reference value of E3 — the between-person SD of the TE slope R2,en − R2,in (ms⁻¹) implied by the assumption of between-person relaxation variation of the design simulation (decomp/design_sim/ANALYSIS_PLAN.md Section 2: lognormal multipliers on T2,in · T2,en, coefficient of variation cv, in–en correlation rho). Monte Carlo (fixed seed)."""
    rng = np.random.default_rng(seed)
    s = math.sqrt(math.log(1 + cv ** 2))
    a = rng.standard_normal(n)
    b = rho * a + math.sqrt(1 - rho ** 2) * rng.standard_normal(n)
    t2in = T2IN * np.exp(s * a - s * s / 2)
    t2en = T2EN * np.exp(s * b - s * s / 2)
    return float(np.std(1 / t2en - 1 / t2in, ddof=1))


def between_sd(slopes, ses):
    """Moment estimate of the SD of the true slopes between people, τ = √max(0, sample variance − mean SE²) — descriptive (very rough with few people)."""
    s, e = np.asarray(slopes, float), np.asarray(ses, float)
    return float(math.sqrt(max(0.0, float(np.var(s, ddof=1)) - float(np.mean(e ** 2)))))


def te_permutation_p(v, person, acq, z):
    """Falsification control (plan Section 5 · change record 1) — empirical null: γ is computed again in every distinct arrangement of the TE values attached to the acquisitions, shuffled between acquisitions (only acquisitions present in the data). If covariate × acquisition deviations unrelated to TE are large, any arrangement gives a large γ, so whether the |γ| of the actual arrangement stands out in that distribution is examined. p = #{|γ*| ≥ |γ̂|} / number of arrangements (including the actual arrangement). Result: p · number of arrangements · 95th percentile of |γ*|."""
    acqs = sorted(set(acq))
    vals = [TE_ACQ[a] for a in acqs]
    g0 = abs(fit(v, person, acq, z)["gamma"])
    perms = sorted(set(itertools.permutations(vals)))
    gs = np.array([abs(fit(v, person, acq, z, te_map=dict(zip(acqs, pv)))["gamma"]) for pv in perms])
    return dict(p=float(np.mean(gs >= g0 - 1e-15)), n_assignments=len(perms), abs_gamma_q95=float(np.quantile(gs, 0.95)))


def acq_slopes_lack_of_fit(v, person, acq, z):
    """Second falsification control — model with a separate covariate slope for each acquisition (saturated) versus a model linear in TE. The difference F (degrees of freedom k − 2, residual) and p describe the size of covariate × acquisition deviations not explained by TE (descriptive)."""
    acqs = sorted(set(acq))
    zbar = float(np.mean(list(z.values())))
    zc = np.asarray([z[p] - zbar for p in person], float)
    acq = np.asarray(acq)
    extra = [(acq == a).astype(float) * zc for a in acqs]
    X_full = design(person, acq, z, extra)                           # the first column (line in TE) is contained in the acquisition-wise slopes and is redundant — the coefficients are not used, only the residual
    v = np.asarray(v, float)
    r_full = v - X_full @ np.linalg.lstsq(X_full, v, rcond=None)[0]
    df_full = len(v) - np.linalg.matrix_rank(X_full)
    lin = fit(v, person, acq, z)
    k_extra = np.linalg.matrix_rank(X_full) - np.linalg.matrix_rank(design(person, acq, z))
    rss_full = float(r_full @ r_full)
    F = ((lin["rss"] - rss_full) / k_extra) / (rss_full / df_full)
    return dict(F=float(F), df1=int(k_extra), df2=int(df_full), p=float(stats.f.sf(F, k_extra, df_full)))


def power_two_sided(effect, se, alpha=0.05):
    """Power by normal approximation Φ(|effect|/SE − z_{1−α/2}) + Φ(−|effect|/SE − z_{1−α/2}) — rule on whether to run (plan Section 6a)."""
    zc = stats.norm.ppf(1 - alpha / 2)
    r = abs(effect) / se
    return float(stats.norm.cdf(r - zc) + stats.norm.cdf(-r - zc))


def go_rule(lit_t2_change_ms_per_decade, se0=None, t2en=T2EN, w_en=0.5):
    """Rule on whether to run (plan Section 6a, fixed before computing) — reads the literature change of white-matter T2 (ms/10 years, magnitude regardless of sign) generously: single-exponential T2 ≈ signal-weighted mean of compartmental T2, with all of its change taken to be in the extra-neurite compartment (weight w_en = 0.5), so the extra-neurite T2 change = L / w_en. The power for that γ is computed with an optimistic SE (se0: σ 0.031 · balanced · plan Section 7), and if it is 0.5 or more, «run»."""
    if se0 is None:
        se0 = 0.031 / math.sqrt(110 * float(np.var(list(TE_ACQ.values()))) * 1.05 ** 2)
    e = abs(float(lit_t2_change_ms_per_decade)) / w_en
    g_len = abs(1.0 / (t2en + e) - 1.0 / t2en)
    g_sh = abs(1.0 / (t2en - e) - 1.0 / t2en)
    g = min(g_len, g_sh)                                                 # the smaller of the two directions (conservative) — the difference is within 5 %
    pw = power_two_sided(g, se0)
    return dict(t2en_change_ms=e, gamma=g, se0=se0, power=pw, run=pw >= 0.5)
