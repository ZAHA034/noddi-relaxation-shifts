"""Quick checks of design_sim M6 (a second TE within each site · plan decomp/design_sim/PLAN_M6.md), each with a failure control."""
import math
import pathlib
import sys

import numpy as np
import pytest

HERE = pathlib.Path(__file__).resolve().parent
for _cand in (HERE.parent / "decomp", HERE.parent):
    if (_cand / "design_sim").is_dir():
        sys.path.insert(0, str(_cand))
        break
from design_sim import m6  # noqa: E402
from design_sim import model as M  # noqa: E402

K = len(M.ACQ)


def _cell(S="S2e", V=0.0, H="H0", D="D1", N="bal"):
    return next(c for c in M.cells() if (c["D"], c["S"], c["V"], c["H"], c["N"]) == (D, S, V, H, N))


@pytest.mark.parametrize("S,V,H", [("S2e", 0.0, "H0"), ("S3b", 0.10, "H1"), ("S0", 0.10, "H2")])
def test_first_observation_equals_the_original_generator(S, V, H):
    c = _cell(S, V, H)
    a = M.generate(c, np.random.default_rng(np.random.SeedSequence([M.SEED, c["idx"]])))
    b = m6.generate2(c, np.random.default_rng(np.random.SeedSequence([M.SEED, c["idx"]])), np.random.default_rng(1))
    assert np.array_equal(a["y"], b["y1"]) and np.array_equal(a["site"], b["site"]) and np.array_equal(a["g"], b["g"])


def _noise_free(S, V_person=None):
    """Makes y1, y2 directly without noise · site terms (relaxation may differ between people)."""
    c = _cell(S)
    site, g = M.layout("bal")
    te_s = M.site_te("D1")
    te = te_s[site]
    th = M.theta_of(S)
    asg = M.SCEN[S][1]
    fin = M.LENGTHEN if asg in ("b", "i") else 1.0
    fen = M.LENGTHEN if asg in ("b", "e") else 1.0
    t2in = np.where(g == 1, M.T2_IN_C * fin, M.T2_IN_C)
    t2en = np.where(g == 1, M.T2_EN_C * fen, M.T2_EN_C)
    if V_person is not None:
        t2en = t2en * V_person
    dr2 = 1 / t2en - 1 / t2in
    base = np.tile(np.log(M.NU_R / (1 - M.NU_R)), (len(site), 1)) + th * g[:, None]
    y1 = base + (te * dr2)[:, None]
    y2 = base + ((te + m6.DTE) * dr2)[:, None]
    return c, y1, y2, g, site, te, te_s, th


@pytest.mark.parametrize("S", ["S0", "S1b", "S2b", "S2e", "S2i", "S3e"])
def test_without_noise_M6a_and_M6b_recover_theta_exactly(S):
    _, y1, y2, g, site, te, te_s, th = _noise_free(S)
    assert np.allclose(m6.m6a(y1, y2, g, site, K, te)["est"], th, atol=1e-10)
    # M6b has variance 0 if the within-group variation is 0 — a very small jitter is added for the weights
    rng = np.random.default_rng(0)
    j = 1e-9 * rng.standard_normal(y1.shape)
    assert np.allclose(m6.m6b(y1 + j, y2 + j, g, site, K, te_s)["est"], th, atol=1e-6)


def test_M6a_is_exact_even_when_relaxation_differs_between_people():
    rng = np.random.default_rng(3)
    _, y1, y2, g, site, te, te_s, th = _noise_free("S2e", V_person=np.exp(0.1 * rng.standard_normal(len(M.layout("bal")[0]))))
    assert np.allclose(m6.m6a(y1, y2, g, site, K, te)["est"], th, atol=1e-10)


def test_failure_control_zero_instead_of_TE_in_the_extrapolation_is_wrong_in_relaxation_cells():
    _, y1, y2, g, site, te, te_s, th = _noise_free("S2e")
    wrong = m6.m6a(y1, y2, g, site, K, np.zeros_like(te))["est"]
    assert np.all(np.abs(wrong - th) > 0.05)


def test_failure_control_stops_when_the_second_TE_is_the_same():
    _, y1, y2, g, site, te, te_s, th = _noise_free("S2e")
    with pytest.raises(ValueError):
        m6.m6a(y1, y1, g, site, K, te, dte=0.0)
    with pytest.raises(ValueError):
        m6.m6b(y1, y1, g, site, K, te_s, dte=0.0)


def test_noise_multiplier_of_the_second_observation_and_half_mask():
    assert m6.NOISE2 == pytest.approx(math.exp(30 / 70))
    site, g = M.layout("bal")
    sub = m6.subset_mask(site, g)
    for s in range(K):
        for gg in (0.0, 1.0):
            a = (site == s) & (g == gg)
            assert sub[a].sum() == a.sum() // 2


def test_there_are_16_target_cells():
    cs = m6.cells_m6()
    assert len(cs) == 16 and {c["S"] for c in cs} == set(M.SCENARIOS) and {c["V"] for c in cs} == {0.0, 0.10}


def test_one_replicate_gives_four_methods():
    c = _cell("S1b")
    r = m6.one(c, np.random.default_rng(1), np.random.default_rng(2))
    assert set(r) == {"M6a", "M6b", "M2", "M3_te0"}


def test_SE_of_M6a_and_M6b_matches_the_spread_over_replicates():
    """Check of the variance formula — in 150 replicates of the same cell, the ratio of the standard deviation of the ROI-mean estimate to the mean SE is 0.8–1.25."""
    c = _cell("S1b", 0.10)
    rng = np.random.default_rng(np.random.SeedSequence([M.SEED, c["idx"]]))
    rng2 = np.random.default_rng(np.random.SeedSequence([m6_seed2(), c["idx"]]))
    ests, ses = {"M6a": [], "M6b": []}, {"M6a": [], "M6b": []}
    for _ in range(150):
        r = m6.one(c, rng, rng2)
        for k in ests:
            ests[k].append(r[k][0][0])
            ses[k].append(r[k][1][0])
    for k in ests:
        ratio = np.mean(ses[k]) / np.std(ests[k], ddof=1)
        assert 0.8 < ratio < 1.25, (k, ratio)


def m6_seed2():
    return 20261118


def test_second_observation_shares_the_person_and_site_common_terms(monkeypatch):
    """With measurement noise set to 0, y2 − y1 = ΔTE·(R2,en − R2,in) — for each person a constant that does not depend on the ROI (failure control: if even one common term is not shared, it differs between ROIs)."""
    monkeypatch.setattr(M, "SD_EPS", 0.0)
    c = _cell("S3e", 0.10, "H1")
    d = m6.generate2(c, np.random.default_rng(5), np.random.default_rng(6))
    diff = d["y2"] - d["y1"]
    assert np.allclose(diff, diff[:, :1], atol=1e-12)
    assert np.all(np.abs(diff[:, 0]) > 0)


def test_with_equal_relaxation_y2_minus_y1_is_the_same_for_everyone(monkeypatch):
    """S0 · V 0 · H1: all relaxation is the same, so the noise-free y2 − y1 is one constant (failure control: dropping the site × group term η from the second observation makes it differ in patients)."""
    monkeypatch.setattr(M, "SD_EPS", 0.0)
    c = _cell("S0", 0.0, "H1")
    d = m6.generate2(c, np.random.default_rng(7), np.random.default_rng(8))
    diff = d["y2"] - d["y1"]
    assert np.allclose(diff, diff[0, 0], atol=1e-12)


def test_first_observation_equals_the_original_generator_over_consecutive_replicates():
    """Random-number stream — drawing three replicates in a row with the same generator must also agree with the original generator (2026-11-18: comparing only the first replicate missed a fault that skipped the random numbers of the measurement error)."""
    c = _cell("S1b", 0.10, "H1")
    ra = np.random.default_rng(np.random.SeedSequence([M.SEED, c["idx"]]))
    rb = np.random.default_rng(np.random.SeedSequence([M.SEED, c["idx"]]))
    r2 = np.random.default_rng(9)
    for _ in range(3):
        assert np.array_equal(M.generate(c, ra)["y"], m6.generate2(c, rb, r2)["y1"])
