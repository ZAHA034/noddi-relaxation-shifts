"""Tests rules_env.metrics: the offset-subtracted retest share (ORS), TE tolerance, scale comparison, proportional bias, rank flips and the S2 reading."""
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
from rules_env import metrics as MT  # noqa: E402
from rules_env._paths import use  # noqa: E402

use("decomp", "remeasure")
import analysis as AN  # noqa: E402


def _units(v):
    return {f"u{i:04d}": float(x) for i, x in enumerate(v)}


def test_ors_null_calibration():
    rng = np.random.default_rng(1)
    d = _units(5.0 + rng.normal(0, math.sqrt(2), 2000))
    r = MT.ors(d, {u: 1.0 for u in d}, c=1.96)
    assert 0.93 <= r["share"] <= 0.96 and r["offset"] == pytest.approx(5.0, abs=0.1) and r["retest_level"] is True


def test_ors_detects_spread():
    rng = np.random.default_rng(2)
    d = _units(rng.normal(0, math.sqrt(2), 20000) + rng.normal(0, 1.5, 20000))
    r = MT.ors(d, {u: 1.0 for u in d}, c=1.96)
    assert r["share"] == pytest.approx(0.821, abs=0.02) and r["retest_level"] is False


def test_ors_raw_matches_retest_level():
    rng = np.random.default_rng(3)
    d = _units(rng.normal(0.01, 0.02, 27))
    s = {u: float(v) for u, v in zip(d, rng.uniform(0.005, 0.02, 27))}
    r, ref = MT.ors(d, s, c=1.96, loo=False), AN.retest_level(d, s, c=1.96)
    assert (round(r["share"] * r["n"]), r["n"], r["retest_level"]) == (ref["n_retest"], ref["n_units"], ref["retest_level"])


def test_ors_undeterminable():
    assert MT.ors(_units(np.zeros(9)), {f"u{i:04d}": 1.0 for i in range(9)}, c=1.96)["retest_level"] is None
    d = _units(np.zeros(100))
    s = {u: (None if i < 11 else 1.0) for i, u in enumerate(d)}
    r = MT.ors(d, s, c=1.96)
    assert r["retest_level"] is None and r["n_excluded"] == 11


def test_ors_diff_zero_same_input():
    rng = np.random.default_rng(4)
    d = _units(rng.normal(0, 1, 60))
    s = {u: 1.0 for u in d}
    r = MT.ors_diff(d, s, d, s, c=1.96, clusters={u: i % 3 for i, u in enumerate(d)}, B=200)
    assert (r["diff"], r["lo"], r["hi"]) == (0.0, 0.0, 0.0)


def test_te_tolerance_closed_form():
    assert MT.te_tolerance(np.ones(10**6), 0.01, span=200.0) == pytest.approx(91.6, abs=0.5)


def test_te_tolerance_cap():
    assert MT.te_tolerance(np.ones(27), 0.0) == math.inf
    assert MT.te_tolerance(np.ones(27), 0.01, span=70.0) == math.inf          # Root ≈ 91 ms > 70
    assert MT.te_tolerance(np.ones(27), 0.05, span=70.0) < 70.0


def test_tau_hat_recovers():
    rng = np.random.default_rng(5)
    k = 0.02 + rng.normal(0, 0.004, 5000) + rng.normal(0, math.sqrt(1e-5), 5000)
    assert MT.tau_hat(k, 1e-5) == pytest.approx(0.004, rel=0.10)


def test_beta_level_zero_independent():
    rng = np.random.default_rng(6)
    k, x = rng.normal(0, 0.004, 5000), rng.normal(0.5, 0.3, 5000)
    assert abs(MT.beta_level(k, x, 0.0)) < 0.01 * 0.004 / 0.3 * 10      # Near 0 if independent (within sampling error)
    assert MT.beta_level(0.01 * x + k, x, 0.0) == pytest.approx(0.01, abs=0.002)


def test_te_interval_contains_point():
    rng = np.random.default_rng(7)
    k = rng.normal(0.02, 0.01, 81)
    pt, lo, hi = MT.te_tolerance_interval(np.full(81, 0.05), k, 1e-6, df=54, clusters=np.repeat([0, 1, 2], 27), B=300, span=200.0)
    assert lo <= pt <= hi


def test_scale_models_picks_logit():
    rng = np.random.default_rng(8)
    nu = rng.uniform(0.5, 0.8, 27)
    d = MT.raw_shift_pred(nu, 0.49) + rng.normal(0, 0.003, 27)
    r = MT.scale_models(nu, d)
    assert r["ss_logit"] < r["ss_raw"] and r["c"] == pytest.approx(0.49, abs=0.03)


def test_scale_models_picks_raw():
    rng = np.random.default_rng(9)
    nu = rng.uniform(0.5, 0.8, 27)
    r = MT.scale_models(nu, 0.09 + rng.normal(0, 0.003, 27))
    assert r["ss_raw"] < r["ss_logit"]


def test_raw_shift_pred_values():
    assert MT.raw_shift_pred(np.array([0.5, 0.6, 0.7, 0.8]), 0.490) == pytest.approx([0.1201, 0.1100, 0.0920, 0.0672], abs=1e-4)


def test_prop_bias_zero_for_translation():
    a = np.linspace(-0.5, 1.5, 27)
    assert MT.prop_bias(a, a + 0.3) == pytest.approx(0.0, abs=1e-12)
    assert MT.prop_bias(a, 1.2 * a) != pytest.approx(0.0, abs=1e-3)


def test_rank_flip_zero_for_translation():
    a = np.linspace(0, 1, 27)
    assert MT.rank_flip_share(a, a + 0.3, rc=0.05) == 0.0
    assert MT.rank_flip_share(a, -a, rc=0.05) == 1.0


def test_s2_reading_bounds():
    assert MT.s2_reading([0.80, 1.0, 1.25]) is True
    assert MT.s2_reading([0.799, 1.0]) is False
    assert MT.s2_ratio(-0.12, -0.148) == pytest.approx(0.12 / 0.148)
