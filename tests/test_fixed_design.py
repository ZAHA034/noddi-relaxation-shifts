"""Known-answer tests of the design calculations in the measurement fixed analysis document, using acquisition records and prior ranges and no ON-Harmony values."""
from __future__ import annotations

import itertools
import math
import pathlib
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "decomp" / "fixed"))
sys.path.insert(0, str(ROOT / "decomp" / "stage1"))
import config as C  # noqa: E402
import delta_tau as DT  # noqa: E402
import fiso_sign as FS  # noqa: E402


# ── D-w2 — sign of the apparent f_iso ─────────────────────────────────────────────────────
def test_apparent_f_iso_diff_equals_config_apparent_independent_of_f_iso():
    rng = np.random.default_rng(0)
    for _ in range(50):
        t2, t1 = C.t2_triplet(C.ETA_T2[rng.integers(len(C.ETA_T2))]), C.t1_triplet(C.T1_SETS[rng.integers(len(C.T1_SETS))])
        aA = C.weights(rng.uniform(60, 110), t2, rng.uniform(3000, 11000), t1)
        aB = C.weights(rng.uniform(60, 110), t2, rng.uniform(3000, 11000), t1)
        nu = rng.uniform(0.3, 0.75)
        want = []
        for f in (0.02, 0.10, 0.25):
            fa, fb = C.apparent(nu, f, aA)[1], C.apparent(nu, f, aB)[1]
            want.append(math.log(fb / (1 - fb)) - math.log(fa / (1 - fa)))
        assert np.allclose(want, FS.dlogit_fiso(aA, aB, nu), atol=1e-12)


def test_same_acquisition_is_0_TR_only_longer_looks_like_more_free_water():
    assert FS.fiso_sign(92.0, 3600.0, 92.0, 3600.0)["L1"]["max_abs"] == 0.0
    s = FS.fiso_sign(80.0, 3600.0, 80.0, 9800.0)                               # Only TR 3.6 → 9.8 s: the free water with long T1 recovers more
    assert s["L1"]["sign"] == "+" and s["L4"]["sign"] == "+"


def test_table_28_pairs_sign_prediction_only_TR_diff_same_L1_L4_sign():
    rows = FS.table()
    assert len(rows) == 28
    for r in rows:
        if r["in_sign_prediction"]:
            assert abs(r["dTR_ms"]) > 0 and r["fiso_sign_L1"] == r["fiso_sign_L4"] != "±"


# ── ZZZw-3 — δ_τ ──────────────────────────────────────────────────────────────
def test_tau_b_equals_scipy():
    stats = pytest.importorskip("scipy.stats")
    rng = np.random.default_rng(1)
    for _ in range(30):
        x, y = rng.integers(0, 5, 21).astype(float), rng.integers(0, 6, 21).astype(float)
        assert DT.tau_b(x, y) == pytest.approx(stats.kendalltau(x, y).statistic, abs=1e-12)


def test_scanner_permutations_are_5040_and_shuffle_all_pairs():
    idx = DT.device_perm_index()
    assert idx.shape == (5040, 21)
    assert (idx[0] == np.arange(21)).all()                                       # First permutation = identity
    assert all(sorted(r) == list(range(21)) for r in idx)
    assert len({tuple(r) for r in idx}) == 5040


def test_null_90_interval_covers_0_about_090_exact_match_not_equivalence():
    rng = np.random.default_rng(2)
    x = rng.permutation(np.linspace(0.01, 0.4, 21))                              # No ties
    idx = DT.device_perm_index()
    for r in (0.0, 1.0, np.inf):
        cov = np.mean([(lambda t: t[1] <= 0 <= t[2])(DT.interval90(x, DT.null_y(rng, r), idx)) for _ in range(300)])
        assert 0.85 <= cov <= 0.96, (r, cov)
    t, lo, hi = DT.interval90(x, x.copy(), idx)
    assert t == pytest.approx(1.0) and hi == pytest.approx(1.0) and not (lo > -0.5 and hi < 0.5)


def test_power_does_not_fall_with_delta_and_is_0_at_0():
    rng = np.random.default_rng(3)
    x = rng.permutation(np.linspace(0.01, 0.4, 21))
    res = DT.design(n_sim=60, x=x, seed=4)
    for v in res["by_r"].values():
        p = list(v["power_at"].values())
        assert p[0] == 0.0 and all(b >= a for a, b in zip(p, p[1:]))


def test_predicted_shift_21_pairs_same_acquisition_OXF1PRI_OXF2PRI_is_0():
    x = DT.pair_predictions()
    assert x.shape == (21,) and (x >= 0).all()
    k = DT.PAIRS.index((DT.DEVICES.index("OXF1PRI"), DT.DEVICES.index("OXF2PRI")))
    assert x[k] == 0.0


# ── Relay ZZZZf-2 ① ③ — scanner-unit pairs ─────────────────────────────────────────
def _counts():
    sys.path.insert(0, str(ROOT / "decomp" / "remeasure"))
    import analysis as AN  # noqa: E402
    import rconfig as RC  # noqa: E402
    return AN.common_counts(RC.onh_sessions())


def test_scanner_sign_prediction_8_of_21_pairs_and_two_no_verdict_pairs():
    rows = {r["pair"]: r for r in FS.device_rows(_counts())}
    assert len(rows) == 21
    inc = {p: r["sign"] for p, r in rows.items() if r["in_sign_prediction"]}
    assert inc == {"OXF1PRI~NOT4GEP": "−", "OXF1PRI~NOT1ACH": "+", "OXF2PRI~NOT4GEP": "−", "OXF2PRI~NOT1ACH": "+",
                   "NOT2ING~NOT4GEP": "−", "NOT4GEP~OXF4GEP": "+", "NOT4GEP~NOT1ACH": "+", "OXF4GEP~NOT1ACH": "+"}
    assert {p for p, r in rows.items() if not r["judgeable"]} == {"OXF3TRI~NOT4GEP", "OXF3TRI~OXF4GEP"}
    assert rows["NOT2ING~NOT4GEP"]["rows"] == ["NOT2ING-98~NOT4GEP"]                   # All 10 common subjects at TE 98
    assert not rows["OXF3TRI~NOT2ING"]["in_sign_prediction"]                            # 89 is −, 98 is + — if mixed, it is dropped
    card = {r["pair"]: r["fiso_max_abs_dlogit_L1"] for r in FS.table()}
    assert rows["OXF1PRI~NOT2ING"]["weights"] == {"OXF1PRI~NOT2ING-89": 0.2, "OXF1PRI~NOT2ING-98": 0.8}       # Common subjects 4 · 16
    assert rows["OXF1PRI~NOT2ING"]["fiso_max_abs_dlogit_L1"] == pytest.approx(0.2 * card["OXF1PRI~NOT2ING-89"] + 0.8 * card["OXF1PRI~NOT2ING-98"])


def test_subject_weighted_prediction_equals_mean_for_equal_protocol_counts():
    even = {k: dict(v, protocols={"NOT2ING-89": 1, "NOT2ING-98": 1} if "NOT2ING" in k else v["protocols"]) for k, v in _counts().items()}
    assert np.allclose(DT.pair_predictions_weighted(even), DT.pair_predictions())


def test_x_and_y_permutation_tau_distributions_match_if_all_pairs_observed():
    rng = np.random.default_rng(5)
    x, y = rng.random(21), rng.random(21)
    idx = DT.device_perm_index()
    a = np.sort(DT.perm_taus(x, y, idx))
    b = np.sort(DT.perm_taus_x(x, y, np.arange(21), idx))
    assert np.allclose(a, b)


# ── Null distribution of the FISO sign test (proposal a9 of A — confirmed as an empty definition on inspection) ─────────────────
def test_scanner_sign_uses_scanner_name_order_and_flips_when_reversed():
    s = DT.sign_vector(FS.device_rows(_counts()))
    k = DT.PAIRS.index((DT.DEVICES.index("NOT4GEP"), DT.DEVICES.index("OXF1PRI")))
    assert s[k] == 1.0                     # Card «OXF1PRI~NOT4GEP» − (= NOT4GEP − OXF1PRI < 0) → OXF1PRI − NOT4GEP > 0
    k2 = DT.PAIRS.index((DT.DEVICES.index("NOT1ACH"), DT.DEVICES.index("NOT4GEP")))
    assert s[k2] == -1.0                   # Card «NOT4GEP~NOT1ACH» + (= NOT1ACH − NOT4GEP > 0) → NOT4GEP − NOT1ACH < 0
    assert (s != 0).sum() == 8 and np.isfinite(DT.magnitude_vector(FS.device_rows(_counts()))).all()


def test_sign_permutation_identity_reproduces_observed_and_count_is_5040():
    s = np.zeros(21)
    s[0], s[5], s[10] = 1, -1, 1
    o = s.copy()                           # The observation equals the prediction
    keep = list(range(21))
    r = DT.sign_perm(s, o, keep)
    assert r["n_perm"] == 5040 and r["S"] == 3.0 and r["agree"] == 3 and r["agree_rate"] == 1.0
    assert 0 < r["p"] < 0.5
    r0 = DT.sign_perm(s, np.zeros(21), keep)
    assert r0["S"] == 0.0 and r0["p"] == 1.0 and r0["agree_rate"] is None


def test_swapping_two_scanners_in_sign_permutation_flips_that_pair():
    """Prediction · observation + only on one pair (0, 1) — +1 if π puts {0, 1} in the same order (5! = 120 ways) · −1 if reversed (120) · 0 otherwise. If the reversal is handled, p = 120/5040; if not handled, 240/5040."""
    s = np.zeros(21)
    k = DT.PAIRS.index((0, 1))
    s[k] = 1.0
    o = s.copy()
    r = DT.sign_perm(s, o, list(range(21)))
    assert r["S"] == 1.0 and r["p"] == pytest.approx(120 / 5040)


# ── Re-check C-2 · C-3 — design sign · FISO τ re-chosen at each permutation · zero size ─────────
def _keep19():
    c = _counts()
    return [k for k, (i, j) in enumerate(DT.PAIRS) if c[frozenset((DT.DEVICES[i], DT.DEVICES[j]))]["judgeable"]]


def _brute_S(s, o, keep):
    """S_π by hand for each permutation and each pair — if the direction is reversed, the sign is changed and only the pairs of keep are added."""
    pos = {pq: k for k, pq in enumerate(DT.PAIRS)}
    out = []
    for perm in itertools.permutations(range(len(DT.DEVICES))):
        tot = 0.0
        for k in keep:
            a, b = perm[DT.PAIRS[k][0]], perm[DT.PAIRS[k][1]]
            tot += (s[pos[(a, b)]] if a < b else -s[pos[(b, a)]]) * o[k]
        out.append(tot)
    return np.array(out)


def test_scanner_sign_is_kept_for_pairs_with_same_name_order():
    s = DT.sign_vector(FS.device_rows(_counts()))
    assert s[DT.PAIRS.index((DT.DEVICES.index("NOT2ING"), DT.DEVICES.index("NOT4GEP")))] == -1.0     # «NOT2ING~NOT4GEP» −
    assert s[DT.PAIRS.index((DT.DEVICES.index("NOT4GEP"), DT.DEVICES.index("OXF4GEP")))] == 1.0      # «NOT4GEP~OXF4GEP» +


def test_sign_permutation_stat_equals_brute_force_over_observed_19_pairs():
    keep = _keep19()
    assert len(keep) == 19 and 13 not in keep and 20 not in keep                 # Pairs with no verdict possible = PAIRS 13 · 20
    rng = np.random.default_rng(12)
    o = rng.choice([-1.0, 1.0], 21)                                              # Nonzero values also for the two pairs with no verdict possible — it differs if keep is ignored
    for s in (DT.design_signs(_counts()), rng.choice([-1.0, 0.0, 1.0], 21)):
        assert np.array_equal(DT.sign_perm_stats(s, o, keep), _brute_S(s, o, keep))
    assert not np.array_equal(_brute_S(s, o, keep), _brute_S(s, o, list(range(21))))


def test_design_sign_adds_no_verdict_pair_to_9_identity_stat_ratio_same():
    c = _counts()
    s8, s9 = DT.sign_vector(FS.device_rows(c)), DT.design_signs(c)
    k13 = DT.PAIRS.index((DT.DEVICES.index("NOT4GEP"), DT.DEVICES.index("OXF3TRI")))
    assert (s9 != 0).sum() == 9 and np.flatnonzero(s8 != s9).tolist() == [k13] and s9[k13] == 1.0   # «OXF3TRI~NOT4GEP» −
    keep = _keep19()
    o = np.random.default_rng(11).choice([-1.0, 1.0], 21)
    a8, a9 = DT.sign_perm_stats(s8, o, keep), DT.sign_perm_stats(s9, o, keep)
    assert a8[0] == a9[0] and not np.array_equal(a8, a9)                       # The term is added only in permutations that carry that prediction over to the observed pair
    r8, r9 = DT.sign_perm(s8, o, keep), DT.sign_perm(s9, o, keep)
    assert r8["S"] == r9["S"] and r8["agree_rate"] == r9["agree_rate"] and r9["n_pairs"] == 8


def test_FISO_tau_reselects_observed_design_sign_pairs_each_permutation():
    c = _counts()
    keep, s = _keep19(), DT.design_signs(c)
    x = DT.magnitude_vector(FS.device_rows(c))
    rng = np.random.default_rng(13)
    y = rng.random(21)
    idx = DT.device_perm_index()
    T = DT.perm_taus_sel(x, y, s, keep, idx)
    k8 = [k for k in keep if s[k] != 0]
    assert len(k8) == 8 and T[0] == pytest.approx(DT.tau_b(x[k8], y[k8]))
    for p in rng.choice(len(idx), 40, replace=False):                           # Indifferent — τ_b with the pair chosen at each permutation
        sel = [k for k in keep if s[idx[p, k]] != 0]
        assert 7 <= len(sel) <= 9 and T[p] == pytest.approx(DT.tau_b(x[idx[p, sel]], y[sel]))
    assert np.allclose(DT.perm_taus_sel(x, y, np.ones(21), np.arange(21), idx), DT.perm_taus_x(x, y, np.arange(21), idx))
    r = DT.fiso_tau(x, y, s, keep)
    assert r["n_pairs"] == 8 and r["tau"] == pytest.approx(T[0]) and r["n_perm"] == 5040


def test_FISO_sign_S_and_tau_null_size_are_near_nominal():
    """One-sided 0.05 rejection rate at the null (grid of scanner effects r · only one scanner jumps) — S ≤ 0.05 · τ allows a small anti-conservatism (≤ 0.06 · the τ of the null in which one scanner jumps in a large simulation is about 0.052 – 0.054 — because the observed pairs are not closed under permutation · recheck FISO-1). Before the fix (erased sign · fixing the 8 identity pairs), τ went up to 0.078."""
    c = _counts()
    x = DT.magnitude_vector(FS.device_rows(c))
    r = DT.fiso_sizes(DT.design_signs(c), x, _keep19(), n_sim=4000)
    assert len(r["by_null"]) == 8 and r["max_S"] <= 0.05 and r["max_tau"] <= 0.06 and r["n_seeds"] == 1, r
    old = DT.fiso_sizes(DT.sign_vector(FS.device_rows(c)), x, _keep19(), n_sim=4000, fixed_keep=True)
    assert old["max_tau"] > 0.07 and old["by_null"]["spike=3.0"]["tau"] > r["by_null"]["spike=3.0"]["tau"]


def test_stops_if_observed_pair_value_missing_ignores_NaN_of_unobserved():
    c = _counts()
    keep, s = _keep19(), DT.design_signs(c)
    x = DT.magnitude_vector(FS.device_rows(c))
    rng = np.random.default_rng(14)
    o, y = rng.choice([-1.0, 1.0], 21), rng.random(21)
    o[13] = o[20] = y[13] = y[20] = np.nan                                         # Pairs with no verdict possible (no observations) — may be used
    base = DT.sign_perm(s, o, keep)
    assert np.isfinite(base["S"]) and np.isfinite(DT.fiso_tau(x, y, s, keep)["tau"])
    o2, y2 = o.copy(), y.copy()
    o2[keep[0]] = y2[keep[0]] = np.nan                                             # There is no value for the pairs for which a verdict is possible
    for f in (lambda: DT.sign_perm(s, o2, keep), lambda: DT.sign_perm_stats(s, o2, keep), lambda: DT.fiso_tau(x, y2, s, keep),
              lambda: DT.perm_taus_x(x, y2, keep, DT.device_perm_index())):
        with pytest.raises(ValueError, match="관측 쌍 값 없음"):
            f()


def test_null_size_seeds_average_rejection_rates_SE_from_total_simulations():
    c = _counts()
    s, x, keep = DT.design_signs(c), DT.magnitude_vector(FS.device_rows(c)), _keep19()
    two = DT.fiso_sizes(s, x, keep, n_sim=200, n_seeds=2, spikes=(3.0,))
    one = [DT.fiso_sizes(s, x, keep, n_sim=200, seed=DT.SEED + k, spikes=(3.0,)) for k in (0, 1)]
    for name, v in two["by_null"].items():
        for key in ("S", "tau"):
            assert v[key] == pytest.approx((one[0]["by_null"][name][key] + one[1]["by_null"][name][key]) / 2)
    assert two["se"] == pytest.approx((0.05 * 0.95 / 400) ** 0.5) and two["n_seeds"] == 2
