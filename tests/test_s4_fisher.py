"""Known-answer tests of s4 Appendix B Fisher sensitivity: ratio ≥ 1, ratio → 1 as w → 0, ratio = 1 when J_ψ = 0."""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "decomp" / "s4"))
import fisher as F  # noqa: E402


def _case(rng, d=6, q=6):
    Jt, Jp = rng.normal(size=(d, 6)), rng.normal(size=(d, q))
    A, B = rng.normal(size=(d, d)), rng.normal(size=(6, 6))
    R = A @ A.T / d + 0.1 * np.eye(d)
    covT = (B @ B.T / 6 + 0.1 * np.eye(6)) * 1e-2
    U = rng.normal(size=(4, 3))
    return Jt, R, covT, Jp, U / np.linalg.norm(U, axis=1, keepdims=True)


def test_ratio_is_at_least_1_Schur_complement():
    rng = np.random.default_rng(0)
    for _ in range(200):
        d = int(rng.choice([6, 14]))
        Jt, R, covT, Jp, U = _case(rng, d)
        w = np.exp(rng.uniform(np.log(1e-3), np.log(1.0), 6))
        r, sdr = F.mde_ratio(Jt, R, covT, Jp, w, U)
        assert (r >= 1 - 1e-10).all(), r
        assert (sdr <= 1 + 1e-10).all(), sdr                    # The data do not make the posterior SD of ψ larger than the prior SD


def test_ratio_is_1_as_width_goes_to_0():
    rng = np.random.default_rng(1)
    for _ in range(50):
        Jt, R, covT, Jp, U = _case(rng)
        r, _ = F.mde_ratio(Jt, R, covT, Jp, np.full(6, 1e-7), U)
        assert np.allclose(r, 1.0, atol=1e-6), r
        rs = [F.mde_ratio(Jt, R, covT, Jp, np.full(6, w), U)[0] for w in (1e-3, 1e-2, 1e-1, 1.0)]
        assert all((b >= a - 1e-10).all() for a, b in zip(rs, rs[1:]))          # The wider, the less the ratio shrinks


def test_J_psi_0_gives_ratio_1_and_psi_posterior_is_prior():
    rng = np.random.default_rng(2)
    Jt, R, covT, _, U = _case(rng)
    w = np.array([0.01, 0.02, 0.05, 0.1, 0.2, 0.3])
    r, sdr = F.mde_ratio(Jt, R, covT, np.zeros((6, 6)), w, U)
    assert np.allclose(r, 1.0, atol=1e-12) and np.allclose(sdr, 1.0, atol=1e-12)


def test_extended_posterior_equals_formula_marginalizing_nuisance_factors():
    rng = np.random.default_rng(3)
    Jt, R, covT, Jp, _ = _case(rng, 14)
    w = np.array([0.01, 0.01, 0.05, 0.05, 0.1, 0.1])
    S1, _ = F.post(Jt, R, covT, Jp, w)
    Rm = R + Jp @ np.diag(w ** 2) @ Jp.T                                          # x | θ ~ N(J_θ θ, R + J_ψ W J_ψᵀ)
    S1m = np.linalg.inv(Jt.T @ np.linalg.solve(Rm, Jt) + np.linalg.inv(covT))
    assert np.allclose(S1, S1m, rtol=1e-8, atol=1e-12)


def test_baseline_equals_validate_fisher_cell():
    pytest.importorskip("torch")
    import validate as V
    rng = np.random.default_rng(4)
    th = rng.normal(size=(100, 6)) * 0.05
    x = th @ rng.normal(size=(6, 14)) + rng.normal(size=(100, 14)) * 0.1
    covT = np.diag(np.full(6, 0.05 ** 2))
    Jt, R = F.regress(x, th)
    assert np.allclose(F.post(Jt, R, covT)[0], V.fisher_cell(x, th, covT), rtol=1e-10, atol=1e-14)


def test_common_random_differencing_recovers_linear_model_slope_exactly():
    rng = np.random.default_rng(5)
    n, d = 10, 6
    base = rng.normal(size=(n, d))
    G = rng.normal(size=(d, 6))                                                   # ∂x/∂ψ (in PSI order)
    X = []
    for _, var in F.VARIANTS:
        psi = np.zeros(6)
        for k, v in var.items():
            nm = dict(bs_A="kappa_A", bs_B="kappa_B").get(k, k)
            psi[F.PSI.index(nm)] = (v - 1.0) if k.startswith("bs_") else v          # κ = b multiplier − 1
        X.append(base + psi @ G.T)
    J = F.jpsi(np.array(X), [v for v, _ in F.VARIANTS])
    assert np.allclose(J, G, atol=1e-12)


def test_gray_matter_mix_is_unchanged_at_g_0_and_gray_matter_at_g_1():
    rng = np.random.default_rng(6)
    Sg, Sgm = rng.uniform(0.1, 1, (8, 5)).astype(np.float32), rng.uniform(0.1, 1, (8, 5)).astype(np.float32)
    s0t, s0g = rng.uniform(0.5, 1, (8, 1)), rng.uniform(0.5, 1, (8, 1))
    assert np.allclose(F.gm_mix(Sg, Sgm, s0t, s0g, 0.0), Sg)
    assert np.allclose(F.gm_mix(Sg, Sgm, s0t, s0g, 1.0), Sgm)
    g = 0.3
    want = ((1 - g) * s0t * Sg + g * s0g * Sgm) / ((1 - g) * s0t + g * s0g)          # Re-normalize the mixing of the signal before normalization
    assert np.allclose(F.gm_mix(Sg, Sgm, s0t, s0g, g), want, rtol=1e-6)


def test_b0_signal_is_fraction_sum_of_compartment_weights():
    pix = np.array([[[0.5, 0.2, 0.3], [0.6, 0.1, 0.2]]])                            # (n=1, K=2, 3)
    a = np.array([[0.4, 0.3, 0.9]])
    want = [(1 - 0.2) * (0.5 * 0.4 + 0.5 * 0.3) + 0.2 * 0.9, (1 - 0.1) * (0.6 * 0.4 + 0.4 * 0.3) + 0.1 * 0.9]
    assert np.allclose(F.s0(pix, a).ravel(), want)


def test_cells_are_map_L0_49_and_card_main_cells_28():
    maps = np.array([f"dTE{t:+.0f}|dTR{r:+.2f}|{lv}" for t in range(7) for r in range(7) for lv in ("L0", "L1", "L2", "L3") for _ in range(3)])
    cards = np.array([f"U{i}~U{j}|{lv}|SNR{s}|Neff{n}" for i in range(8) for j in range(i + 1, 8)
                      for lv, s, n in [("L0", 30, 50), ("L1", 30, 50), ("L2", 30, 50), ("L3", 30, 50), ("L0", 20, 25), ("L0", 20, 100), ("L0", 50, 25), ("L0", 50, 100)]
                      for _ in range(3)])
    cells = F.fisher_cells(maps, cards)
    assert len(cells) == 77
    assert sum(nm == "maps" for nm, _, _ in cells) == 49 and sum(nm == "cards" for nm, _, _ in cells) == 28
    assert all(cl.endswith("|L0") for nm, cl, _ in cells if nm == "maps")
    assert all(cl.endswith("|L0|SNR30|Neff50") for nm, cl, _ in cells if nm == "cards")
    assert all(len(ii) == 3 and (np.diff(ii) > 0).all() for _, _, ii in cells)


def test_width_settings_five_per_factor_and_five_for_all_three():
    cf = F.configs()
    assert len(cf) == 20
    cols, w = F.active(dict(cf)["all|4"])
    assert sorted(cols.tolist()) == list(range(6)) and np.allclose(sorted(w), sorted([0.05, 0.05, 0.20, 0.20, 0.30, 0.30]))
    cols, w = F.active(dict(cf)["G1|2"])
    assert cols.tolist() == [F.PSI.index("p_in"), F.PSI.index("p_iso")] and np.allclose(w, 0.10)
