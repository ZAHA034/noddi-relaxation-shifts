"""Tests the pure functions of rules_env.analyses: fiber-direction mean coordinates, S2 group differences and methods, TE slopes, batch offsets, ComBat reapplication and L1 relaxation draws."""
import math
import pathlib
import sys

import numpy as np
import pytest

HERE = pathlib.Path(__file__).resolve().parent
for _cand in (HERE.parent, HERE.parent / "decomp"):
    if (_cand / "rules_env").is_dir():
        sys.path.insert(0, str(_cand))
        break
from rules_env import analyses as AS  # noqa: E402
from rules_env import model as M  # noqa: E402


def test_coords_mean_averages_fdir_fastest():
    ndi = np.array([0.5, 0.5, 0.6, 0.6, 0.7, 0.7, 0.7, 0.7])
    fwf = np.full(8, 0.1)
    out = AS.coords_mean(dict(ndi=ndi, fwf=fwf, odi=np.full(8, 0.2)), n_fdir=4)
    assert out["x"].shape == (2,)
    assert out["x"][0] == pytest.approx((2 * math.log(1) + 2 * math.log(0.6 / 0.4)) / 4)
    assert out["raw_ndi"][1] == pytest.approx(0.7)


def _algebraic_s2(te_sites, ctrl, pat, nu, d_nu=-0.03, n=20, tr=3600.0):
    """Algebraic apparent values instead of a fit — the position method must give equation (7) as it is."""
    X = {}
    for s, te in enumerate(te_sites):
        site = {"TE": te, "TR": tr}
        for g, grp, dn in (("c", ctrl, 0.0), ("p", pat, d_nu)):
            w = M._group_weights(site, grp)
            v = np.array([math.log(M.apparent(u + dn, 0.05, w)[0] / (1 - M.apparent(u + dn, 0.05, w)[0])) for u in nu])
            X[(s, g)] = np.repeat(v[None, :], n, 0)
    return X


def test_s2_location_methods_give_eq7_on_algebra():
    ctrl = {"t2": (82.5, 52.0, 1000.0), "t1": None}
    pat = {"t2": (82.5, 70.0, 1000.0), "t1": None}
    nu = np.linspace(0.4, 0.7, 27)
    X = _algebraic_s2((70.0, 100.0), ctrl, pat, nu)
    rng = np.random.default_rng(0)
    eq7 = M.prop4_change({"TE": 100.0}, {"TE": 70.0}, ctrl, pat)
    for method in ("none", "ctrl_offset", "combat_loc", "combat_train_ctrl_loc"):
        D = AS.s2_site_diffs(X, n_sites=2, method=method, noise_sd=np.full(27, 1e-4), rng=rng)
        assert D[1] - D[0] == pytest.approx(eq7, abs=2e-3), method


def test_s2_fraction_only_logit_has_no_interaction():
    c = {"t2": (82.0, 51.25, 1000.0), "t1": None}
    X = _algebraic_s2((70.0, 100.0), c, c, np.linspace(0.4, 0.7, 27))
    D = AS.s2_site_diffs(X, n_sites=2, method="none", noise_sd=np.zeros(27), rng=np.random.default_rng(1))
    assert D[1] - D[0] == pytest.approx(0.0, abs=1e-12)


def test_te_slopes_linear():
    te = np.array([62, 72, 82, 92, 102, 112, 122, 132], float)
    ser = {f"TE{int(t)}": {"a": 0.1 + 0.004 * t, "b": -0.2 + 0.001 * t} for t in te}
    k, xbar = AS.te_slopes(ser)
    assert k == pytest.approx({"a": 0.004, "b": 0.001}) and xbar["a"] == pytest.approx(0.1 + 0.004 * 97)


def test_grouped_loo_offsets():
    d = {("g1", "u1"): 1.0, ("g1", "u2"): 2.0, ("g1", "u3"): 3.0, ("g2", "u1"): 10.0, ("g2", "u2"): 10.0, ("g2", "u3"): 13.0}
    off = AS.grouped_loo_offsets(d)
    assert off[("g1", "u1")] == pytest.approx(2.5) and off[("g2", "u3")] == pytest.approx(10.0)


def test_l1_draws_truncated_and_seeded():
    a = AS.l1_draws((82.0, 51.25, 1000.0), (1100.0, 1100.0, 4300.0), 500, seed=3)
    b = AS.l1_draws((82.0, 51.25, 1000.0), (1100.0, 1100.0, 4300.0), 500, seed=3)
    assert a == b
    z_in = np.array([t2[0] / 82.0 - 1 for t2, _ in a]) / 0.03
    z_en = np.array([t2[1] / 51.25 - 1 for t2, _ in a]) / 0.03
    assert np.all(np.hypot(z_in, z_en) <= 2 + 1e-9) and all(t2[2] == 1000.0 and t1[2] == 4300.0 for t2, t1 in a)
    assert all(abs(t1[0] / 1100.0 - 1) <= 0.08 + 1e-9 and t1[0] == t1[1] for _, t1 in a)
    assert all(t1 is None for _, t1 in AS.l1_draws((82.5, 52.0, 1000.0), None, 5, seed=1))


def test_combat_apply_reproduces_training_output():
    # neuroCombatFromTraining (0.2.12) breaks under numpy 2 (int() of a 1-element array) — combat_apply, which uses the same formula,
    # If the trained estimate is applied again to the training data, it must equal neuroCombat's own output (no covariates · location only · location + scale)
    rng = np.random.default_rng(5)
    dat = rng.normal(0, 1, (27, 60)) + np.repeat(rng.normal(0, 1, (27, 3)), 20, 1)
    batch = np.repeat([0, 1, 2], 20)
    for mean_only in (True, False):
        tr = AS._combat(dat, batch, None, mean_only=mean_only)
        assert AS.combat_apply(dat, batch, tr["estimates"]) == pytest.approx(tr["data"], abs=1e-10)


def test_s2_subjects_population_mode_is_centre():
    gd = dict(t2=(82.0, 51.25, 1000.0), t1=(1100.0, 1100.0, 4300.0), dnu=0.0)
    assert AS.s2_subjects("pop", gd, 20, seed=1) == [((82.0, 51.25, 1000.0), (1100.0, 1100.0, 4300.0))]
    l1 = AS.s2_subjects("l1", gd, 20, seed=1)
    assert len(l1) == 20 and l1 == AS.l1_draws(gd["t2"], gd["t1"], 20, seed=1)


def test_s2_sample_algebra_cancels_draws():
    # Algebraic interaction from the relaxations of the same people — with all people's relaxations centered, it equals formula (7)
    c = dict(t2=(82.5, 52.0, 1000.0), t1=None, dnu=0.0)
    p = dict(t2=(82.5, 70.0, 1000.0), t1=None, dnu=-0.03)
    subj = {("A", "c"): [(c["t2"], None)] * 5, ("A", "p"): [(p["t2"], None)] * 5,
            ("B", "c"): [(c["t2"], None)] * 5, ("B", "p"): [(p["t2"], None)] * 5}
    sites = {"A": {"TE": 100.0, "TR": 3600.0}, "B": {"TE": 70.0, "TR": 3600.0}}
    d = AS.s2_sample_algebra(subj, sites)
    assert d["A"] - d["B"] == pytest.approx(M.prop4_change({"TE": 100.0}, {"TE": 70.0}, c, p), abs=1e-12)


# ── R2 (ANALYSIS_PLAN §7): synthetic input only ────────────────────────────────────────
@pytest.mark.parametrize("n", [2, 3, 10, 27, 28])
def test_loo_median_matches_grouped_offsets(n):
    rng = np.random.default_rng(n)
    V = rng.normal(size=(5, n))
    for r in range(5):
        d = {("g", i): float(V[r, i]) for i in range(n)}
        ref = AS.grouped_loo_offsets(d)
        assert np.allclose(AS.loo_median(V[r]), [ref[("g", i)] for i in range(n)])
    assert np.allclose(AS.loo_median(V), np.stack([AS.loo_median(V[r]) for r in range(5)]))


def test_group_share_equals_ors_rule():
    rng = np.random.default_rng(3)
    d = {"a": rng.normal(0.2, 0.05, 27), "b": rng.normal(-0.1, 0.05, 27)}
    s = {"a": 0.03, "b": 0.03}
    dd = {(g, i): float(v) for g, arr in d.items() for i, v in enumerate(arr)}
    off = AS.grouped_loo_offsets(dd)
    want = AS.MT.ors({k: v - off[k] for k, v in dd.items()}, {k: s[k[0]] for k in dd}, 1.96, loo=False)["share"]
    assert abs(float(AS.group_share(d, s)) - want) < 1e-12


def _par(beta, kbar=0.009, tau=0.003, spread=1.2, n=27, sigma=0.055):
    xb = np.linspace(-spread / 2, spread / 2, n) + 0.6
    s_raw = float(np.mean(AS._expit(xb) * (1 - AS._expit(xb)))) * sigma
    return {"p": dict(xbar=xb, kbar=kbar, beta=beta, tau=tau, var_x=sigma ** 2 / 8, sigma=sigma, s_x=sigma, s_raw=s_raw)}


def test_e5_predict_no_shift_scales_tie():
    r = AS.e5_predict(_par(0.0, kbar=0.0, tau=0.0), M_=400, seed=1)
    assert abs(r["diff"]) < 0.03


def test_e5_predict_shared_shift_favours_logit():
    r = AS.e5_predict(_par(0.0, tau=0.0), M_=400, seed=1)
    assert r["ors_x"] > 0.9 > r["ors_raw"] and r["diff"] > 0


def test_e5_predict_level_dependent_shift_favours_raw():
    xb = _par(0.0)["p"]["xbar"]
    nu_mid = float(np.mean(AS._expit(xb)))
    b = (2 * nu_mid - 1) * 0.009 * 1.6                    # The slope grows greatly with level → logit is not even, raw is more even
    r = AS.e5_predict(_par(b, tau=abs(b) * 0.35), M_=400, seed=1)
    assert r["diff"] < 0
    assert AS.e5_predict(_par(b, tau=abs(b) * 0.35), M_=400, seed=1, use_beta=False)["diff"] > r["diff"]


def test_e5_first_order_and_readings():
    f = AS.e5_first_order(0.009, 0.006, 0.65)
    assert math.isclose(f["R"], 0.006 / 0.009) and math.isclose(f["T"], 0.3) and f["raw_more_uniform"] is True
    assert AS.e5_first_order(0.009, 0.0, 0.65)["raw_more_uniform"] is False
    assert AS.e5_reading(-0.08, -0.148, -0.037) == "가" and AS.e5_reading(0.01, -0.148, -0.037) == "나"
    assert AS.e6_reading(0.95, 0.5) == "가" and AS.e6_reading(0.95, 0.92) == "나" and AS.e6_reading(0.5, 0.6) == "나" and AS.e6_reading(0.8, 0.95) == "다"


def test_e6_predict_constant_logit_shift():
    xb = np.linspace(-0.2, 1.6, 27)
    r = AS.e6_predict({"p": xb}, {"p": xb + 0.6}, {"p": 0.055}, {"p": 0.055}, {"p": 0.011}, M_=300, seed=0)
    assert r["noise_free_shift_sd"]["p"]["sd_x"] < 1e-12 and r["noise_free_shift_sd"]["p"]["sd_raw"] > 0.01
    assert r["reading"] == "가"


def test_s3_algebra_logit_zero_and_integral_identity():
    """raw change = ∫_{c_B}^{c_A} [g(x_p + c) − g(x_c + c)] dc (g = ν(1 − ν)) — the value is not printed, only true or false."""
    dr = 1 / 52.0 - 1 / 82.5
    ok_logit, ok_int = True, True
    for nc in (0.46, 0.60):
        for dn in (-0.03, -0.10, -0.20):
            v = AS.s3_algebra(nc, dn)
            ok_logit &= abs(v["logit"]) < 1e-12
            xc, xp = math.log(nc / (1 - nc)), math.log((nc + dn) / (1 - nc - dn))
            c = np.linspace(70 * dr, 100 * dr, 2001)
            g = lambda z: AS._expit(z) * (1 - AS._expit(z))
            f = g(xp + c) - g(xc + c)
            h = c[1] - c[0]
            simpson = h / 3 * (f[0] + f[-1] + 4 * f[1:-1:2].sum() + 2 * f[2:-1:2].sum())
            ok_int &= abs(v["raw"] - simpson) < 1e-10
    assert ok_logit
    assert ok_int
    assert abs(AS.s3_algebra(0.5, 0.0)["raw"]) < 1e-15


# ── Tests added after the 2026-11-15 code review (mutation tests) — synthetic input only · assertions that print no values ────────────────
MT = AS.MT


def test_ors_integer_rule_at_boundary():
    n = 20
    need = MT._need(n, 0.90)
    s = {u: 1.0 for u in range(n)}
    thr = 1.96 * math.sqrt(2.0)
    at = {u: (0.0 if u < need else 10 * thr) for u in range(n)}          # Exactly `need` pass
    below = {u: (0.0 if u < need - 1 else 10 * thr) for u in range(n)}   # need − 1 items
    ok = MT.ors(at, s, 1.96, loo=False)["retest_level"] is True and MT.ors(below, s, 1.96, loo=False)["retest_level"] is False
    assert ok


def test_ors_offset_is_leave_one_out_median():
    d = {u: 0.0 for u in range(10)}
    d[10] = 100.0                                                        # One outlier: with a median offset the other 10 pass, with a mean it fails
    r = MT.ors(d, {u: 0.1 for u in d}, 1.96, loo=True)
    assert abs(r["share"] - 10 / 11) < 1e-12


def test_te_share_finite_n_factor_and_group_size():
    s = np.full(81, 0.05)
    ok = True
    for n in (27, 81):
        want = 2 * AS.MT.stats.norm.cdf(1.96 * math.sqrt(2) * 0.05 / math.sqrt(2 * 0.05 ** 2 * (1 + 1 / (n - 1)))) - 1
        got = MT._te_share(0.0, s, 0.003, 1.96, n_group=n)
        ok &= abs(got - want) < 1e-12
    ok &= abs(MT._te_share(0.0, s, 0.003, 1.96) - MT._te_share(0.0, s, 0.003, 1.96, n_group=81)) < 1e-15
    t27, t81 = MT.te_tolerance(s, 0.003, n_group=27), MT.te_tolerance(s, 0.003)
    ok &= t27 < t81                                                       # The smaller the batch, the larger the offset noise, so the tolerance narrows
    ok &= MT.te_tolerance(s, 0.003, n_group=27, infl_kind="median") < t27  # The median factor is larger
    assert ok


def test_tau_hat_pooled_ddof():
    rng = np.random.default_rng(2)
    g = [rng.normal(m, 0.003, 27) for m in (0.008, 0.009, 0.010)]
    kc = np.concatenate([x - x.mean() for x in g])
    want = math.sqrt(sum(((x - x.mean()) ** 2).sum() for x in g) / (81 - 3))
    assert abs(MT.tau_hat(kc, 0.0, ddof=3) - want) < 1e-15 and MT.tau_hat(kc, 0.0) < MT.tau_hat(kc, 0.0, ddof=3)


def test_beta_level_noise_term():
    rng = np.random.default_rng(4)
    xt = rng.normal(0, 0.3, 20000)
    k = 0.01 + 0.006 * xt
    var_x = 0.2 ** 2
    xbar = xt + rng.normal(0, 0.2, xt.size)
    b_corr, b_naive = MT.beta_level(k, xbar, var_x), MT.beta_level(k, xbar, 0.0)
    assert abs(b_corr - 0.006) < 0.0003 and b_naive < 0.0050


def test_te_tolerance_interval_includes_sd_sampling():
    rng = np.random.default_rng(6)
    k = rng.normal(0, 0.003, 81)
    s = np.full(81, 0.05)
    cl = np.repeat([0, 1, 2], 27)
    _, lo_s, hi_s = MT.te_tolerance_interval(s, k, 0.0, 5, cl, B=400, seed=0)
    _, lo_l, hi_l = MT.te_tolerance_interval(s, k, 0.0, 10 ** 7, cl, B=400, seed=0)
    assert (hi_s - lo_s) > (hi_l - lo_l)                                   # With small degrees of freedom it is as wide as the sampling variation of s


def test_ors_diff_resamples_clusters():
    d1 = {f"{g}_{i}": (0.0 if g == 0 else 100.0) for g in range(3) for i in range(30)}
    d2 = {k: 0.0 for k in d1}
    s = {k: 1.0 for k in d1}
    cl = {k: k.split("_")[0] for k in d1}
    w_cl = MT.ors_diff(d1, s, d2, s, 1.96, cl, B=400, seed=0, loo=False)
    w_un = MT.ors_diff(d1, s, d2, s, 1.96, {}, B=400, seed=0, loo=False)
    assert (w_cl["hi"] - w_cl["lo"]) > 2 * (w_un["hi"] - w_un["lo"])


def test_combat_apply_rejects_covariate_trained_estimates():
    rng = np.random.default_rng(7)
    dat = rng.normal(0, 1, (27, 60)) + np.repeat(rng.normal(0, 1, (27, 3)), 20, 1)
    batch = np.repeat([0, 1, 2], 20)
    tr = AS._combat(dat, batch, np.tile([0, 1], 30), mean_only=True)
    with pytest.raises(ValueError):
        AS.combat_apply(dat, batch, tr["estimates"])


def test_e5_tau_e_uses_beta_squared_and_true_level_variance():
    xb = np.linspace(-1, 1, 27)
    var_x = 0.01
    vtrue = float(np.var(xb, ddof=1)) - var_x
    tau = 0.004
    ok = abs(AS.e5_tau_e(tau, 0.002, xb, var_x) - math.sqrt(tau ** 2 - 0.002 ** 2 * vtrue)) < 1e-15
    ok &= AS.e5_tau_e(tau, 0.0, xb, var_x) == tau
    ok &= AS.e5_tau_e(tau, -0.002, xb, var_x) == AS.e5_tau_e(tau, 0.002, xb, var_x)
    ok &= AS.e5_tau_e(0.001, 0.01, xb, var_x) == 0.0
    assert ok


def test_e5_predict_uses_half_te_span():
    xb = np.linspace(-0.2, 1.4, 27)
    par = lambda kbar: {"p": dict(xbar=xb, kbar=kbar, beta=0.0, tau=0.0, var_x=1e-6, sigma=1e-6, s_x=0.01, s_raw=0.002)}
    # No noise · even slope: the x difference is a single 70·k̄ → only the spread of the raw difference remains; if the half-width is not 35, the raw share changes
    r35 = AS.e5_predict(par(0.009), M_=5, seed=0)
    r70 = AS.e5_predict(par(0.018), M_=5, seed=0, half=17.5)
    assert abs(r35["ors_raw"] - r70["ors_raw"]) < 1e-12 and r35["ors_x"] == 1.0


def test_e6_noise_goes_through_inverse_logit():
    xb = np.linspace(0.2, 1.4, 27)
    nu = 1 / (1 + np.exp(-xb))
    sigma = 0.05
    s_raw = float(np.mean(nu * (1 - nu))) * sigma
    r = AS.e6_predict({"p": xb}, {"p": xb.copy()}, {"p": sigma}, {"p": sigma}, {"p": s_raw}, M_=400, seed=0)
    assert abs(r["ors_x"] - r["ors_raw"]) < 0.03 and r["ors_x"] < 0.99


def test_e5_first_order_raw_rate_sign():
    f = AS.e5_first_order(0.009, 0.006, 0.65)
    assert abs(f["raw_rate"] - (f["R"] - f["T"])) < 1e-15


def test_stop_gates_and_e6_input_checks(tmp_path):
    ok = AS.gate_close(-0.0864, -0.0864) and not AS.gate_close(-0.0864, -0.0865) and not AS.gate_close(None, 0.0)
    (tmp_path / "fits.npz").write_bytes(b"a"); (tmp_path / "jobs.json").write_bytes(b"b")
    import hashlib
    good = {"fits.npz": hashlib.sha256(b"a").hexdigest(), "jobs.json": hashlib.sha256(b"b").hexdigest()}
    ok &= AS.e6_inputs_ok(tmp_path, good) and not AS.e6_inputs_ok(tmp_path, {**good, "jobs.json": "0" * 64})
    c = [dict(nu=0.5, fiso=0.1, odi=0.2, t2=[82.0, 51.25, 1000.0], t1=None, gen="T-in", fdir=0, dpar=1.7)]
    c2 = [dict(c[0], t2=[90.0, 51.25, 1000.0])]
    ok &= AS.e6_cells_match(c, [dict(c[0])]) and not AS.e6_cells_match(c, c2) and not AS.e6_cells_match(c, c + c)
    assert ok


def test_need_integer_rule_when_ratio_times_n_is_not_whole():
    ok = MT._need(27, 0.90) == 25 and MT._need(20, 0.90) == 18 and MT._need(81, 0.90) == 73
    s = {u: 1.0 for u in range(27)}
    thr = 1.96 * math.sqrt(2.0)
    d24 = {u: (0.0 if u < 24 else 10 * thr) for u in range(27)}            # 24 / 27 = 0.889 < 0.90
    ok &= MT.ors(d24, s, 1.96, loo=False)["retest_level"] is False
    assert ok


def test_slope_level_noise_terms():
    vk, vx = AS.slope_level_noise(2.0)
    assert abs(vk - 4.0 / 4200.0) < 1e-18 and abs(vx - 4.0 / 8.0) < 1e-15


def test_e1_fix_variants_follow_section_8():
    ng = np.full(81, 27.0)
    kw = {k: AS.e1_fix_kwargs(v, ng, 3) for k, v in AS.E1_FIX_VARIANTS.items()}
    ok = kw["as_run"] == {} and set(kw["primary"]) == {"n_group"} and kw["sens_df"]["tau_ddof"] == 3
    ok &= kw["sens_median"]["infl_kind"] == "median" and set(kw["sens_both"]) == {"n_group", "infl_kind", "tau_ddof"}
    ok &= all(kw[k]["n_group"] is ng for k in ("primary", "sens_df", "sens_median", "sens_both"))
    assert ok


def test_e6_raw_noise_scales_with_nu_through_inverse_logit():
    xb = np.r_[np.full(14, -3.0), np.full(13, 3.0)]                        # ν(1 − ν) ≈ 0.045: raw noise far below a ν = 0.5 threshold
    sigma = 0.05
    r = AS.e6_predict({"p": xb}, {"p": xb.copy()}, {"p": sigma}, {"p": sigma}, {"p": 0.25 * sigma}, M_=400, seed=0)
    assert r["ors_raw"] > 0.99 and r["ors_x"] < 0.99


def test_te_tolerance_interval_scales_s_by_chi_draw_not_inverse():
    rng = np.random.default_rng(8)
    k = rng.normal(0, 0.003, 2000)
    s = np.full(2000, 0.05)
    pt, lo, hi = MT.te_tolerance_interval(s, k, 0.0, 5, np.zeros(2000), B=600, seed=0)
    assert 1.2 < hi / pt < 1.75                                             # √(χ²₀.₉₅(5)/5) ≈ 1.49; with 1/m it would be ≈ 2.1
