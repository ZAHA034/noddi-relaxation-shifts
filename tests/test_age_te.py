"""age_te — quick checks of the covariate × scanner TE interaction (plan decomp/age_te/ANALYSIS_PLAN.md) on synthetic data with known answers, each with a failure control; reads no real data."""
import math
import pathlib
import sys

import numpy as np
import pytest

HERE = pathlib.Path(__file__).resolve().parent
for _cand in (HERE.parent / "decomp", HERE.parent):
    if (_cand / "age_te").is_dir():
        sys.path.insert(0, str(_cand))
        break
from age_te import model as M  # noqa: E402
from age_te.participants import assert_no_ids, parse_participants  # noqa: E402

ACQ = list(M.TE_ACQ)


def _design(seed=1, n_person=20, pool_all=None):
    """Unbalanced design similar to ON-Harmony: 3–8 acquisitions per person (only one of the two NOT2ING protocols per person). pool_all can reduce the acquisitions."""
    rng = np.random.default_rng(seed)
    acq_all = ACQ if pool_all is None else pool_all
    persons, acqs = [], []
    z = {}
    for i in range(n_person):
        p = f"p{i:02d}"
        z[p] = rng.normal(3.15, 1.05)                       # age in units of 10 years (of the form of the data-paper summary 31.5 ± 10.5 years)
        pool = [a for a in acq_all if a != ("NOT2ING-89" if i % 5 else "NOT2ING-98")]
        k = rng.integers(3, len(pool) + 1)
        for a in rng.choice(pool, size=k, replace=False):
            persons.append(p)
            acqs.append(str(a))
    return persons, acqs, z


def _gen(persons, acqs, z, gamma, noise=0.0, seed=2, beta_z=0.0, vendor_slope=None, acq_slope=None):
    rng = np.random.default_rng(seed)
    ps = sorted(set(persons))
    alpha = {p: rng.normal(0.5, 0.2) + beta_z * z[p] for p in ps}           # person level (including the covariate main effect)
    mu = {a: 0.004 * (M.TE_ACQ[a] - M.TE_C) + rng.normal(0, 0.05) for a in ACQ}   # acquisition level (including the TE main effect)
    zbar = float(np.mean(list(z.values())))
    v = []
    for p, a in zip(persons, acqs):
        val = alpha[p] + mu[a] + gamma * (M.TE_ACQ[a] - M.TE_C) * (z[p] - zbar)
        if vendor_slope:
            val += vendor_slope.get(M.VENDOR[a], 0.0) * (z[p] - zbar)
        if acq_slope:
            val += acq_slope.get(a, 0.0) * (z[p] - zbar)
        v.append(val + rng.normal(0, noise))
    return np.array(v)


def test_without_noise_the_interaction_coefficient_is_recovered_exactly():
    persons, acqs, z = _design()
    g = -0.0007
    f = M.fit(_gen(persons, acqs, z, g), persons, acqs, z)
    assert abs(f["gamma"] - g) < 1e-10


def test_person_and_acquisition_main_effects_do_not_move_the_coefficient():
    """γ stays the same even when the main effect of age (absorbed by the person fixed effects) and the main effect of TE (absorbed by the acquisition fixed effects) are large."""
    persons, acqs, z = _design()
    a = M.fit(_gen(persons, acqs, z, -0.0007), persons, acqs, z)["gamma"]
    b = M.fit(_gen(persons, acqs, z, -0.0007, beta_z=5.0), persons, acqs, z)["gamma"]
    assert abs(a - b) < 1e-10


def test_failure_control_fitting_with_shuffled_TEs_gives_a_wrong_coefficient():
    persons, acqs, z = _design()
    g = -0.0007
    v = _gen(persons, acqs, z, g)
    saved = dict(M.TE_ACQ)
    try:
        vals = list(saved.values())
        M.TE_ACQ.update(dict(zip(saved, vals[::-1])))
        wrong = M.fit(v, persons, acqs, z)["gamma"]
    finally:
        M.TE_ACQ.clear()
        M.TE_ACQ.update(saved)
    assert abs(wrong - g) > 0.2 * abs(g)


def test_with_noise_the_interval_covers_the_true_value_and_the_SE_matches_the_Monte_Carlo_spread():
    persons, acqs, z = _design()
    g, est, cover, ses = -0.0007, [], 0, []
    for r in range(300):
        f = M.fit(_gen(persons, acqs, z, g, noise=0.03, seed=100 + r), persons, acqs, z)
        est.append(f["gamma"])
        ses.append(f["se"])
        cover += f["ci95"][0] <= g <= f["ci95"][1]
    assert abs(np.mean(est) - g) < 4 * np.std(est) / math.sqrt(len(est))
    assert 0.85 < np.mean(ses) / np.std(est) < 1.15
    assert 0.91 <= cover / 300 <= 0.98


def test_under_the_null_the_t_interval_rejection_rate_is_near_5_percent():
    persons, acqs, z = _design(seed=3)
    rej = sum(not (f["ci95"][0] <= 0 <= f["ci95"][1])
              for f in (M.fit(_gen(persons, acqs, z, 0.0, noise=0.03, seed=500 + r), persons, acqs, z) for r in range(400)))
    assert 0.02 <= rej / 400 <= 0.09


def test_permutation_p_is_minimal_for_a_strong_signal_and_not_small_under_the_null():
    persons, acqs, z = _design()
    strong = _gen(persons, acqs, z, -0.003, noise=0.01)
    assert M.permutation_p(strong, persons, acqs, z, n_perm=199) == pytest.approx(1 / 200)
    null = _gen(persons, acqs, z, 0.0, noise=0.03, seed=7)
    p = M.permutation_p(null, persons, acqs, z, n_perm=199)
    assert 0.0 < p <= 1.0
    assert M.permutation_p(null, persons, acqs, z, n_perm=199) == p        # fixed seed — the same when run again


def test_permutation_p_null_rejection_rate_is_near_5_percent():
    persons, acqs, z = _design(seed=4)
    rej = sum(M.permutation_p(_gen(persons, acqs, z, 0.0, noise=0.03, seed=900 + r), persons, acqs, z, n_perm=99, seed=r) < 0.05
              for r in range(80))
    assert rej / 80 <= 0.15


def test_vendor_specific_covariate_slopes_mislead_the_main_model_but_not_the_sensitivity_model():
    """With only a vendor × age term unrelated to TE: γ of the main model comes out nonzero (failure control — why this sensitivity is needed), and the model with vendor-wise slope terms recovers γ = 0."""
    persons, acqs, z = _design()
    v = _gen(persons, acqs, z, 0.0, vendor_slope={"Siemens": 0.02, "Philips": 0.01})
    main = M.fit(v, persons, acqs, z)["gamma"]
    zbar = float(np.mean(list(z.values())))
    extra = [np.array([(M.VENDOR[a] == vv) * (z[p] - zbar) for p, a in zip(persons, acqs)], float) for vv in ("Philips", "Siemens")]
    sens = M.fit(v, persons, acqs, z, extra)["gamma"]
    assert abs(main) > 1e-4
    assert abs(sens) < 1e-10


def test_the_ROI_unit_level_model_also_recovers_exactly_without_noise():
    persons, acqs, z = _design()
    g = -0.0005
    base = _gen(persons, acqs, z, g)
    units = [{f"u{k}": (val + 0.1 * k, 0.0, 0.0) for k in range(5)} for val in base]
    f = M.fit_units(units, persons, acqs, z, coord=0)
    assert abs(f["gamma"] - g) < 1e-9


def test_extra_neurite_T2_conversion_keeps_the_reciprocal_relation():
    assert M.t2en_change(0.0) == {"ms": 0.0, "pct": 0.0}
    t2p = 55.0
    g = 1 / t2p - 1 / M.T2EN
    assert M.t2en_change(g)["ms"] == pytest.approx(t2p - M.T2EN, abs=1e-9)
    assert M.t2en_change(-0.0005)["ms"] > 0 > M.t2en_change(0.0005)["ms"]       # if R2,en falls (γ < 0), extra-neurite T2 lengthens


def test_TE_slope_and_homogeneity():
    tes = list(range(62, 133, 10))
    f = M.te_slope(tes, [0.3 + 0.002 * t for t in tes])
    assert f["slope"] == pytest.approx(0.002) and f["se"] < 1e-12
    same = M.slope_heterogeneity([0.002, 0.002, 0.002], [1e-4, 1e-4, 1e-4])
    assert same["chi2"] == pytest.approx(0.0) and same["p"] == pytest.approx(1.0)
    diff = M.slope_heterogeneity([0.001, 0.002, 0.003], [1e-4, 1e-4, 1e-4])
    assert diff["p"] < 1e-6


def test_minimum_detectable_size_is_about_2_8_times_the_SE():
    assert M.mde(1.0, 1000) == pytest.approx(1.96 + 0.8416, rel=2e-3)


def test_participant_information_parsing_rules():
    text = "participant_id\tage\tsex\textra\nsub-a\t31\tF\tq\nsub-b\t25-30\tm\tq\nsub-c\tn/a\tfemale\tq\nsub-d\t40\tn/a\tq\n"
    age, female, fmt = parse_participants(text)
    assert age == {"sub-a": 3.1, "sub-b": 2.75, "sub-d": 4.0}
    assert female == {"sub-a": 1.0, "sub-b": 0.0, "sub-c": 1.0}
    assert fmt == dict(age_column=True, sex_column=True, age_bins=1, age_missing=1, sex_missing=1)


def test_output_containing_a_participant_id_stops():
    assert_no_ids('{"n_persons": 20, "participant_1": {}}', ids=["sub-x1"])
    with pytest.raises(SystemExit):
        assert_no_ids('{"who": "sub-12345"}')
    with pytest.raises(SystemExit):
        assert_no_ids('{"who": "x1"}', ids=["x1"])


def test_assumed_SD_and_between_person_SD_estimate():
    assert M.assumed_slope_sd(cv=0.0) == pytest.approx(0.0, abs=1e-15)
    sd = M.assumed_slope_sd()
    # delta approximation: square root of 0.1² · (R2en² + R2in² − 2 · 0.8 · R2en · R2in)
    ren, rin = 1 / M.T2EN, 1 / M.T2IN
    assert sd == pytest.approx(0.1 * math.sqrt(ren ** 2 + rin ** 2 - 1.6 * ren * rin), rel=0.05)
    assert M.assumed_slope_sd(rho=0.0) > sd                              # a lower correlation widens the spread of the difference
    assert M.between_sd([0.002, 0.002, 0.002], [1e-4] * 3) == 0.0
    assert M.between_sd([0.001, 0.002, 0.003], [1e-6] * 3) == pytest.approx(0.001, rel=1e-3)


SIX = ["NOT1ACH", "NOT2ING-98", "NOT4GEP", "OXF1PRI", "OXF3TRI", "OXF4GEP"]       # six with all TEs different — 720 arrangements (test time)


def test_TE_permutation_falsification_control_flags_an_interaction_linear_in_TE():
    persons, acqs, z = _design(pool_all=SIX)
    r = M.te_permutation_p(_gen(persons, acqs, z, -0.001, noise=0.005), persons, acqs, z)
    assert r["n_assignments"] == 720
    assert r["p"] < 0.05            # arrangements that swap acquisitions with nearly the same TE (68 · 69.65 · 70.6 ms, 96.4 · 98 ms) also give similar |γ|, so the smallest p is larger than 1/720


def test_failure_control_TE_unrelated_acquisition_slopes_fool_the_person_permutation_but_not_the_TE_permutation():
    """If the age slope differs between acquisitions independently of TE (acquisition × age deviations), the person permutation p can come out small, but those deviations give a similar |γ| in any TE arrangement, so the TE permutation p is not small — why this falsification control is needed. The acquisition slopes are set out of step with the TE order."""
    persons, acqs, z = _design(pool_all=SIX)
    order = sorted(SIX, key=lambda a: M.TE_ACQ[a])
    slope = dict(zip(order, [0.03, -0.03, 0.03, -0.03, 0.03, -0.03]))         # alternating in TE order — nearly unrelated to the line in TE
    v = _gen(persons, acqs, z, 0.0, noise=0.005, acq_slope=slope)
    assert M.permutation_p(v, persons, acqs, z, n_perm=199) < 0.05
    assert M.te_permutation_p(v, persons, acqs, z)["p"] > 0.10
    lof = M.acq_slopes_lack_of_fit(v, persons, acqs, z)
    assert lof["df1"] == len(SIX) - 2 and lof["p"] < 1e-6


def test_when_acquisition_slopes_are_linear_in_TE_the_departure_test_is_not_small():
    persons, acqs, z = _design(pool_all=SIX)
    lof = M.acq_slopes_lack_of_fit(_gen(persons, acqs, z, -0.001, noise=0.03, seed=11), persons, acqs, z)
    assert lof["df1"] == len(SIX) - 2 and lof["p"] > 0.01


def test_threshold_of_the_rule_on_whether_to_run():
    assert M.go_rule(0.3)["run"] is False and M.go_rule(1.0)["run"] is True
    lo, hi = M.go_rule(0.62), M.go_rule(0.66)                              # plan Section 6a: threshold about 0.64 ms/10 years
    assert lo["run"] is False and hi["run"] is True
    assert M.go_rule(0.64)["power"] == pytest.approx(0.5, abs=0.01)
    assert M.power_two_sided(0.0, 1.0) == pytest.approx(0.05)
