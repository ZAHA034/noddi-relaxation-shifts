"""Known-answer tests for s4 per SPEC_S4 §9 and §13, run in the main environment with numpy and torch only."""
from __future__ import annotations

import math
import pathlib
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "decomp" / "s4"))
import s4cfg as S  # noqa: E402


def test_prior_is_truncated_uniform_main_prior_is_part_of_training_prior():
    rng = np.random.default_rng(0)
    t = S.sample_theta(rng, 200_000)
    assert S.feasible(t).all()
    assert np.all(np.abs(t[:, 0:3]) <= S.D_TRAIN + 1e-12)
    m1, m05, m2 = S.in_main(t), S.in_main(t, 0.5), S.in_main(t, 2.0)
    assert m2.all() and (m05 <= m1).all() and 0.05 < m1.mean() < 0.5      # Half of the half-width → about 1/8 (slightly larger because of truncation)
    tm = S.sample_theta(rng, 200_000, S.D_MAIN)
    assert S.in_main(tm).all()
    # Importance weight = indicator function (§3-2): among the training-prior samples, those inside the main support have the same distribution as the main prior
    sub = t[m1]
    for j in range(6):
        a, b = np.quantile(sub[:, j], [0.25, 0.75]), np.quantile(tm[:, j], [0.25, 0.75])
        assert np.allclose(a, b, atol=0.006), j


def test_prior_variance_is_quadratic():
    rng = np.random.default_rng(1)
    t = S.sample_theta(rng, 300_000, S.D_MAIN)
    cov = np.cov(t[:, 0:3].T)
    U = rng.normal(size=(5, 3))
    U /= np.linalg.norm(U, axis=1, keepdims=True)
    v_mc = (t[:, 0:3] @ U.T).var(0)
    v_q = np.einsum("ni,ij,nj->n", U, cov, U)
    assert np.allclose(v_mc, v_q, rtol=1e-6)


def test_ROI_noise_covariance_equals_formula():
    rng = np.random.default_rng(2)
    d = 3
    A = rng.normal(size=(d, d))
    Sig = A @ A.T + 0.1 * np.eye(d)
    n = 200_000
    m = np.zeros((n, d))
    neff = np.full(n, 40.0)
    x = S.roi_noise(rng, m, np.repeat(Sig[None], n, 0), neff, K=256)
    want = Sig * (1 / 40.0 - 1 / 256)
    assert np.allclose(np.cov(x.T), want, rtol=0.03, atol=1e-6)
    assert np.allclose(S.roi_noise(rng, m[:5], np.repeat(Sig[None], 5, 0), np.full(5, 256.0), K=256), 0, atol=1e-5)  # N_eff = K → noise 0


def test_condition_vector_shape_and_one_hot_coding():
    rng = np.random.default_rng(3)
    acq = S.sample_acq(rng, 50)
    lv, p = S.sample_level(rng, 50)
    gen, scen = rng.choice(S.GENERATORS, 50), rng.choice(S.SCENARIOS, 50)
    c = S.cond_vector(acq, lv, p, gen, scen)
    assert c.shape == (50, len(S.COND_NAMES)) == (50, 34)
    for grp in ([f"L_{l}" for l in S.L_LEVELS], [f"gen_{g}" for g in S.GENERATORS], [f"scen_{s}" for s in S.SCENARIOS]):
        assert np.all(c[:, [S.COND_NAMES.index(k) for k in grp]].sum(1) == 1)
    assert np.all(acq["neff"] <= S.K_PIX)


def test_acquisition_types_respect_stratum_definition():
    rng = np.random.default_rng(4)
    kinds = np.repeat(S.ACQ_TYPES, 200)
    a = S.sample_acq(rng, len(kinds), kinds)
    k = a["kind"]
    same = lambda x, y: np.isclose(a[x], a[y])
    assert (same("TE_A", "TE_B") & same("TR_A", "TR_B") & (a["design_A"] == a["design_B"]))[k == "match"].all()
    assert (~same("TE_A", "TE_B") & same("TR_A", "TR_B"))[k == "te"].all()
    assert (same("TE_A", "TE_B") & ~same("TR_A", "TR_B"))[k == "tr"].all()
    assert (a["design_A"] != a["design_B"])[k == "tetr_design"].all()
    assert (a["design_A"] == a["design_B"])[k == "tetr"].all()


def test_relaxation_levels():
    rng = np.random.default_rng(5)
    n = 400
    lv = np.repeat(S.L_LEVELS, n // 4)
    _, p = S.sample_level(rng, n, lv)
    t2A, t1A = S.sample_relax_A(rng, n, np.full(n, "iii"))
    t2B, t1B = S.relax_B(rng, lv, p, t2A, t1A)
    L0 = lv == "L0"
    assert np.allclose(t2B[L0], t2A[L0]) and np.allclose(t1B[L0], t1A[L0])
    L3 = lv == "L3"
    assert np.allclose(t2B[L3, 1] / t2A[L3, 1], p["m_en"][L3]) and np.allclose(t2B[L3, 0] / t2A[L3, 0], p["m_in"][L3])
    L2 = lv == "L2"
    en = L2 & (p["alloc_en"] == 1)
    assert np.allclose(t2B[en, 0], t2A[en, 0]) and np.allclose(t2B[en, 1] / t2A[en, 1], 1 + S.L2_T2)
    assert np.allclose(t2B[:, 2], t2A[:, 2]) and np.allclose(t1B[:, 2], t1A[:, 2])                    # Degrees of freedom fixed
    L1 = lv == "L1"
    r = t2B[L1, 0:2] / t2A[L1, 0:2] - 1
    s = p["s_t2"][L1][:, None]
    z = np.where(s > 0, r / np.where(s > 0, s, 1), 0)
    assert (np.sqrt((z ** 2).sum(1)) <= 2.0 / math.sqrt(1 - 0.8 ** 2) + 1e-6).all()


def test_case_list_is_deterministic():
    import plan
    a, b = plan.train(30), plan.train(30)
    for k in a:
        assert np.array_equal(a[k], b[k]), k


def test_table_interpolation_exact_at_grid_points_linear_between():
    import tables as TB
    rng = np.random.default_rng(6)
    T = dict(odi=TB.ODI_T, ratio=TB.RATIO_T, ic=rng.random((2, len(TB.ODI_T), 8, 5)).astype(np.float32),
             ec=rng.random((2, len(TB.RATIO_T), len(TB.ODI_T), 8, 5)).astype(np.float32))
    di, fd = np.array([1]), np.array([3])
    j, r = 7, 11
    assert np.allclose(TB.interp_ic(T, di, np.array([TB.ODI_T[j]]), fd)[0], T["ic"][1, j, 3])
    assert np.allclose(TB.interp_ec(T, di, np.array([TB.RATIO_T[r]]), np.array([TB.ODI_T[j]]), fd)[0], T["ec"][1, r, j, 3])
    mid = (TB.ODI_T[j] + TB.ODI_T[j + 1]) / 2
    assert np.allclose(TB.interp_ic(T, di, np.array([mid]), fd)[0], (T["ic"][1, j, 3] + T["ic"][1, j + 1, 3]) / 2, atol=1e-6)


def test_noncentral_chi_noise_second_moment():
    import simulate
    rng = np.random.default_rng(7)
    n, V = 20_000, 4
    Sg = np.full((n, V), 0.6, np.float32)
    snr = np.full(n, 20.0)
    for L in (4, 8):
        x = simulate.noisy(rng, Sg, snr, np.full(n, L))
        want = 0.6 ** 2 + 2 * L * (1 / 20.0) ** 2                       # E[M²] = S² + 2Lσ²
        assert abs(np.mean(x.astype(float) ** 2) - want) / want < 0.01, L


def test_oracle_observation_noise_and_fractions():
    import simulate
    rng = np.random.default_rng(8)
    P = __import__("plan").oracle_train(2000)
    o = simulate.oracle_obs(rng, P, np.arange(2000))
    clean = simulate.oracle_obs(np.random.default_rng(9), {**P}, np.arange(2000)) - o
    assert o.shape == (2000, 6) and abs(np.std(clean) / math.sqrt(2) - simulate.ORACLE_SD) < 4e-4
    assert np.allclose(o[:, 2], P["theta"][:, 5], atol=0.03) and np.allclose(o[:, 5], P["theta"][:, 5] + P["theta"][:, 2], atol=0.03)   # ODI is independent of relaxation


def test_u_relax_absent_in_matched_pairs_TR_only_is_nearly_FISO_direction():
    pytest.importorskip("torch")
    import validate as V
    rng = np.random.default_rng(10)
    n = 300
    th = S.sample_theta(rng, n, S.D_MAIN)
    acq = S.sample_acq(rng, n, np.r_[np.repeat("match", 150), np.repeat("tr", 150)])
    P = {f"acq_{k}": v for k, v in acq.items()}
    U = V.directions(th, P, np.arange(n))
    assert np.isnan(U[:150, 3]).all()
    assert (np.abs(U[150:, 3, 1]) > 0.9).mean() > 0.9                    # Stage 2 X1: TR fan 85–94° — nearly pure FISO
    assert np.allclose(np.einsum("ni,ni->n", U[150:, 3], U[150:, 4]), 0, atol=1e-9)


def test_shrinkage_near_0_for_prior_samples_and_near_1_if_concentrated():
    pytest.importorskip("torch")
    import validate as V
    rng = np.random.default_rng(11)
    covs = dict(main=V.cov_prior(rng, n=200_000)[0:3, 0:3], half=V.cov_prior(rng, 0.5, n=200_000)[0:3, 0:3], two=V.cov_prior(rng, 2.0, n=200_000)[0:3, 0:3])
    T = S.sample_theta(rng, 10_000, S.D_MAIN)
    U = np.eye(3)
    r = V.stats_one(T, S.sample_theta(rng, 10_000, S.D_TRAIN), True, U, covs, sens=True)
    assert np.all(np.abs(r["c"]) < 0.05)
    assert np.all(np.abs(r["c_two"]) < 0.05)                                  # ×2 = posterior of the training prior without reweighting (ZZZg-3)
    assert 1000 < r["n_half"] < 2000 and np.all(np.abs(r["c_half"]) < 0.08)    # ×0.5 = additional reweighting · remaining sample ≈ 1/8
    Tc = T.copy()
    Tc[:, 0:3] = T[:, 0:3] * 0.01 + 0.0
    r = V.stats_one(Tc, Tc, True, U, covs)
    assert np.all(r["c"] > 0.99)
    assert np.isnan(V.stats_one(T, T, False, U, covs)["c"]).all()
    assert list(V.classify(np.array([0.8, 0.5, 0.1, np.nan]))) == ["식별", "부분", "비식별", "재가중 불가"]


# ── Relay ZZZg (2026-10-25 · results 0 · training 0) — conditional posterior per k · reference posterior MCMC · misspecification score ──────────────────
class _MockPost:
    """Reads the k of condition a (tortuosity × scenario one-hot) and draws samples with ΔODI = 0.01 · k index — the rest is the uniform training prior."""

    def __init__(self, seed):
        self.rng = np.random.default_rng(seed)

    def sample_batched(self, shape, x, **kw):
        import torch
        V = __import__("validate")
        c = x.numpy()[:, -len(S.COND_NAMES):]
        k = np.array([V.GEN_SCEN.index((S.GENERATORS[int(np.argmax(r[[S.COND_NAMES.index(f"gen_{g}") for g in S.GENERATORS]]))],
                                        S.SCENARIOS[int(np.argmax(r[[S.COND_NAMES.index(f"scen_{s}") for s in S.SCENARIOS]]))])) for r in c])
        t = S.sample_theta(self.rng, shape[0] * len(c), S.D_TRAIN).reshape(shape[0], len(c), 6)
        t[:, :, 2] = 0.01 * k[None, :]
        return torch.as_tensor(t, dtype=torch.float32)


def test_samples_come_from_conditional_posterior_of_true_or_fixed_k():
    pytest.importorskip("torch")
    import validate as V
    rng = np.random.default_rng(12)
    acq = S.sample_acq(rng, 6)
    lv, p = S.sample_level(rng, 6)
    gen, scen = np.array(S.GENERATORS)[[0, 1, 2, 0, 1, 2]], np.array(S.SCENARIOS)[[0, 0, 1, 1, 2, 2]]
    cond = S.cond_vector(acq, lv, p, gen, scen)
    X = np.zeros((6, 6), np.float32)
    posts = [_MockPost(1), _MockPost(2)]
    T, tag, T2, ok, acc = V.draw(posts, X, cond, 300)
    assert ok.all() and np.isfinite(T).all() and S.in_main(T.reshape(-1, 6)).all()
    want = [V.GEN_SCEN.index((g, s)) for g, s in zip(gen, scen)]
    assert np.allclose(T[:, :, 2], 0.01 * np.array(want)[:, None])                         # gs=None → the true k of the case
    assert set(np.unique(tag)) == {0, 1} and np.all((acc > 0.1) & (acc < 0.6))             # In the two models · main support ratio
    T, *_ = V.draw(posts, X, cond, 300, gs=V.GEN_SCEN[5])
    assert np.allclose(T[:, :, 2], 0.05)                                                     # gs → substituted by that k
    Pm, w = V.mix_pool([np.zeros((100, 6)), np.ones((300, 6))], [0.2, 0.6])
    assert math.isclose(w[:100].sum(), 0.25) and math.isclose(w[100:].sum(), 0.75)          # Reference mixture: weight ∝ main support ratio


def test_card_class_is_worst_k_and_split_means_assumption_dependent():
    pytest.importorskip("torch")
    import validate as V
    C = [0.9] * 8 + [0.5]
    r = V.worst_of_k(C, [0.01] * 8 + [0.03])
    assert r["cls"] == "부분" and r["C_worst"] == 0.5 and r["MDE_worst"] == 0.03 and r["C_median"] == 0.9 and r["assumption_dependent"]
    r = V.worst_of_k([0.8] * 9, [0.02] * 9)
    assert r["cls"] == "식별" and not r["assumption_dependent"]
    assert V.worst_of_k([0.8] * 8 + [None], [0.02] * 9) is None                             # If even one k is missing, «재가중 불가» («reweighting not possible»)


def test_misspecification_score_is_uniform_mixture_density_of_nine_k():
    pytest.importorskip("torch")
    import validate as V
    rng = np.random.default_rng(13)
    q = rng.uniform(0.01, 5.0, (50, 9))
    assert np.allclose(V.mix_score(np.log(q)), -np.log(q.mean(1)))
    assert np.allclose(V.mix_score(np.full((3, 9), -1000.0)), 1000.0)                       # No overflow even at a very small density


def test_split_Rhat_and_bulk_ESS_known_answers():
    pytest.importorskip("torch")
    import validate as V
    rng = np.random.default_rng(14)
    r, e = V.rhat_ess(rng.standard_normal((4, 2500)))
    assert r <= 1.01 and 0.85 * 10_000 < e < 1.15 * 10_000                                  # Independent samples: ESS ≈ total count
    rho, x = 0.5, np.zeros((4, 2500))
    eps = rng.standard_normal((4, 2500))
    for t in range(1, 2500):
        x[:, t] = rho * x[:, t - 1] + math.sqrt(1 - rho ** 2) * eps[:, t]
    r, e = V.rhat_ess(x)
    assert r <= 1.01 and abs(e / (10_000 * (1 - rho) / (1 + rho)) - 1) < 0.15                # AR(1): ESS = N(1 − ρ)/(1 + ρ)
    y = rng.standard_normal((4, 2500)) + np.array([0, 0, 0, 1.0])[:, None]
    assert V.rhat_ess(y)[0] > 1.05                                                           # If one chain is off, R̂ is large


def _ref_cases(n, levels, seed):
    import plan
    import simulate
    rng = np.random.default_rng(seed)
    th, acq, lv, p, gen, scen = plan._draw(rng, n, S.D_MAIN, np.full(n, "tetr"), np.array(levels))
    P = plan._finish(rng, th, acq, lv, p, gen, scen)
    return P, dict(oracle=simulate.oracle_obs(rng, P, np.arange(n)), cond=P["cond"])


def test_reference_posterior_Jacobian_equals_finite_differences():
    pytest.importorskip("torch")
    import validate as V
    P, sim = _ref_cases(3, ["L0", "L1", "L3"], 15)
    K = V.ref_consts(sim, P, np.arange(3))
    st = V.ref_init(np.random.default_rng(16), K)
    th, logj, ok = V.ref_theta(st, K)
    assert ok.all()
    h = 1e-6
    for i in range(3):
        J = np.zeros((6, 6))
        for d in range(6):
            sp, sm = st[i:i + 1].copy(), st[i:i + 1].copy()
            sp[0, d] += h
            sm[0, d] -= h
            Ki = {k: v[i:i + 1] for k, v in K.items()}
            J[:, d] = (V.ref_theta(sp, Ki)[0][0] - V.ref_theta(sm, Ki)[0][0]) / (2 * h)
        assert abs(np.log(abs(np.linalg.det(J))) - logj[i]) < 1e-5, i


def test_reference_MCMC_matches_fixed_relaxation_analytic_solution():
    """If the relaxation of the three compartments is the same, the apparent fraction = the true fraction → the posterior is Gaussian: A ~ N(observed_A, σ²) · Δ ~ N(observed_B − observed_A, 2σ²) · Cov(Δ, A) = −σ²."""
    pytest.importorskip("torch")
    import validate as V
    from simulate import ORACLE_SD as sd
    P, sim = _ref_cases(2, ["L0", "L0"], 17)
    th0 = np.array([[0.03, 0.01, -0.02, 0.55, 0.15, 0.24], [-0.04, -0.01, 0.03, 0.62, 0.12, 0.20]])
    obs = np.c_[th0[:, 3:6], th0[:, 3:6] + th0[:, 0:3]] + np.random.default_rng(18).normal(0, sd, (2, 6))
    sim["oracle"] = obs
    K = V.ref_consts(sim, P, np.arange(2))
    one = np.ones((2, 3))
    fixed = ((80.0 * one, 1000.0 * one), (80.0 * one, 1000.0 * one))
    TH = V.ref_mcmc(K, seed=19, fixed=fixed)
    for j in range(2):
        d = [V.rhat_ess(TH[j, :, :, k]) for k in range(6)]
        assert max(r for r, _ in d) <= 1.01 and min(e for _, e in d) >= 500
        R = TH[j].reshape(-1, 6)
        mu = np.r_[obs[j, 3:6] - obs[j, 0:3], obs[j, 0:3]]
        sdv = np.r_[np.full(3, math.sqrt(2) * sd), np.full(3, sd)]
        assert np.all(np.abs(R.mean(0) - mu) < 4 * sdv / math.sqrt(500)), j
        assert np.allclose(R.std(0), sdv, rtol=0.07), j
        assert np.allclose(np.cov(R.T)[0:3, 3:6], -sd ** 2 * np.eye(3), atol=0.15 * sd ** 2), j


def test_reference_MCMC_equals_importance_sampling_posterior_with_nuisance():
    """In the case where σ is made large (0.05) so that importance sampling is usable — the change of variables · Jacobian · relaxation (z · z′ of L0 · L1) give the same posterior as the old calculation."""
    pytest.importorskip("torch")
    import validate as V
    P, sim = _ref_cases(2, ["L0", "L1"], 20)
    K = V.ref_consts(sim, P, np.arange(2))
    TH = V.ref_mcmc(K, seed=21, sd=0.05)
    for j in range(2):
        R = TH[j].reshape(-1, 6)
        th, w, ess = V.ref_posterior_is(np.random.default_rng(22 + j), sim["cond"][j], P, j, sim["oracle"][j], n=400_000, sd=0.05)
        assert ess > 2000, ess
        m_is = w @ th
        s_is = np.sqrt(w @ (th - m_is) ** 2)
        assert np.all(np.abs(R.mean(0) - m_is) < 0.12 * s_is + 3 * s_is / math.sqrt(min(ess, 2000))), j
        assert np.allclose(R.std(0), s_is, rtol=0.12), j
