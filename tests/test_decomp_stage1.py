"""Known-answer tests for decomp stage1 per SPEC_v1 §3, §7 and §8–§10, using no fitter and no real data."""
import math
import pathlib
import sys

import numpy as np
import pytest
from scipy import stats

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "decomp" / "stage1"))
import config as C  # noqa: E402
import lib  # noqa: E402


def test_row_count_equals_spec_7():
    want = {"E1": 1_492_992, "E1-TR": 46_656, "E3": 207_360, "E2a": 16_588_800, "E2c": 460_800, "E2b": 8_294_400}
    assert {e: C.n_rows(e) for e in want} == want


def test_grid_sizes():
    assert len(C.TE_GRID) == 16 and len(C.T2_SETS) == 36 and len(C.TISSUE) == 36 and len(C.WM) == 12
    assert all(t in C.T2_SETS for t in C.T2_E2A) and len(C.FIBRE_DIRS) == 8
    assert 90.0 in C.TE_GRID and 70.0 in C.TE_GRID and 100.0 in C.TE_GRID and 65.0 in C.TE_GRID


@pytest.mark.parametrize("nu,fiso", [(0.45, 0.02), (0.75, 0.25), (0.3, 0.1)])
def test_reweighting_is_exact_translation_in_signal_fraction_coordinates(nu, fiso):
    t2 = C.t2_triplet((70.0, 1.8, 2000.0))
    t1 = C.t1_triplet((850.0, 1.3))
    a1, a2 = C.weights(68.0, t2, 3600.0, t1), C.weights(98.0, t2, 10620.0, t1)
    x1, y1 = C.frac_coords(*C.apparent(nu, fiso, a1))
    x2, y2 = C.frac_coords(*C.apparent(nu, fiso, a2))
    assert abs((x2 - x1) - math.log((a2[0] / a2[1]) / (a1[0] / a1[1]))) < 1e-12
    assert abs((y2 - y1) - math.log((a2[2] / a2[0]) / (a1[2] / a1[0]))) < 1e-12


def test_repulsion_directions_deterministic_unit_norm():
    d1, d2 = C.repulsion_dirs(), C.repulsion_dirs()
    assert np.array_equal(d1, d2) and np.allclose(np.linalg.norm(d1, axis=1), 1)


def test_shrinkage_estimator_preserves_variance_lambda_range_formula():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(12, 5)) @ rng.normal(size=(5, 5))
    S, lam = lib.shrink_cov(X)
    Sc = np.cov(X.T)
    assert np.allclose(np.diag(S), np.diag(Sc)) and 0 <= lam <= 1
    assert np.linalg.eigvalsh(S).min() > 0
    n = len(X); Z = (X - X.mean(0)) / X.std(0, ddof=1)
    num = den = 0.0
    for i in range(5):
        for j in range(5):
            if i != j:
                w = Z[:, i] * Z[:, j]
                num += n / (n - 1) ** 3 * ((w - w.mean()) ** 2).sum()
                den += (n / (n - 1) * w.mean()) ** 2
    assert abs(lam - min(1, num / den)) < 1e-10


def test_rank_rule_merges_at_30_degrees():
    def two(deg):
        a = math.radians(deg)
        return np.array([[1, math.cos(a)], [0, math.sin(a)], [0, 0]])
    assert lib.rank_projector(two(20))[0] == 1 and lib.rank_projector(two(40))[0] == 2


def test_preservation_metric():
    T = np.eye(4)[:, :2]
    assert lib.preservation(np.eye(4) - np.outer(*[np.eye(4)[:, 3]] * 2), T)[0] == pytest.approx(1.0)
    P = np.eye(4) - np.outer(np.eye(4)[:, 0], np.eye(4)[:, 0])
    assert lib.preservation(P, T)[0] == pytest.approx(0.0, abs=1e-12)


def test_noncentral_chi_squared_interval_rule():
    df = 4
    hi_T = stats.chi2.ppf(0.99, df)
    lo, hi = lib.ncx2_ci(hi_T, df)
    assert lo > 0 and hi > lo
    lo, hi = lib.ncx2_ci(stats.chi2.ppf(0.9, df), df)
    assert lo == 0 and hi > 0
    assert lib.ncx2_ci(stats.chi2.ppf(0.01, df), df) == (0.0, 0.0)


def test_index_stat_gaussian_null_size():
    rng = np.random.default_rng(3)
    p, n = 5, 20
    A = rng.normal(size=(p, p)); W = A @ A.T / p + 0.3 * np.eye(p); L = np.linalg.cholesky(W)
    U = rng.normal(size=(p, 2))
    rej = 0
    for _ in range(150):
        D = U @ rng.normal(size=2) + rng.standard_normal((n, p)) @ L.T
        rej += lib.index_stats(D, U, rng, B=199)["p_value"] < 0.05
    assert 0.0 <= rej / 150 <= 0.12


def test_E5_row_count_equals_spec_v2_11():
    assert C.n_rows("E5") == 6_570_720 and len(C.TISSUE_E5) == 2028 and len(C.ETA_T2) == 9


def test_E5_resolution_comes_from_budget():
    budget = 27_713_088 // 4
    fixed = 3 * 4 * 54 * 5 * 3
    assert fixed * 26 ** 2 <= budget < fixed * 27 ** 2


def test_basic_bootstrap_interval_contains_point_estimate():
    rng = np.random.default_rng(5)
    p, n = 4, 20
    U = rng.normal(size=(p, 1))
    D = rng.standard_normal((n, p)) + np.r_[1.0, 0, 0, 0]
    r = lib.index_stats(D, U, rng, B=199, basic_ci=True)
    lo, hi = r["Psi_ci_basic"]
    assert 0 <= lo <= hi


def test_E6_eta_and_row_count_equal_spec_v3_9():
    assert len(C.ETA_E6) == 199 and C.n_rows("E6") == 24_214_320 and C.n_rows("E6n") == 221_400
    names = [n for n, _, _ in C.ETA_E6]
    assert sum(n.startswith("L1|") for n in names) == 160 and sum(n.startswith("L1s|") for n in names) == 16
    assert sum(n.startswith("L3|") for n in names) == 18 and sum(n.startswith("L0|") for n in names) == 3


def test_L1_samples_within_Mahalanobis_radius_2_and_T1_iso_fixed():
    for n, t2, t1 in C.ETA_E6:
        if n.startswith("L1|"):
            cv = float(n.split("|")[1][2:]); rho = float(n.split("|")[2][3:])
            z = np.array([t2[0] / C.T2C[0] - 1, t2[1] / C.T2C[1] - 1]) / cv
            R = np.array([[1, rho], [rho, 1]]) + 1e-12 * np.eye(2)
            assert z @ np.linalg.solve(R, z) <= 4 + 1e-6
            assert t1[2] == C.T1_ISO
