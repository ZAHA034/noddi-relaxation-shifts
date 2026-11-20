"""design_sim — quick checks of the multi-site design simulation (plan decomp/design_sim/ANALYSIS_PLAN.md), each with a failure control."""
import pathlib
import sys

import numpy as np
import pytest

HERE = pathlib.Path(__file__).resolve().parent
for _cand in (HERE.parent / "decomp", HERE.parent):
    if (_cand / "design_sim").is_dir() and (_cand / "rules_env").is_dir():
        sys.path.insert(0, str(_cand))
        break
from design_sim import analyses as A  # noqa: E402
from design_sim import checks as K  # noqa: E402
from design_sim import model as M  # noqa: E402


def test_there_are_288_cells_numbered_in_order():
    cs = M.cells()
    assert len(cs) == 288 and [c["idx"] for c in cs] == list(range(288))


def test_D2_rearranges_the_D1_TEs_and_each_vendor_has_a_short_and_a_long_TE():
    assert sorted(M.TE_D2.values()) == sorted(M.TE_D1.values())
    for v in ("Siemens", "Philips", "GE"):
        tes = [M.TE_D2[a] for a in M.ACQ if M.VENDOR[a] == v]
        assert min(tes) < 75 < max(tes)


def test_failure_control_D1_does_not_give_each_vendor_a_short_and_a_long_TE():
    tes = [M.TE_D1[a] for a in M.ACQ if M.VENDOR[a] == "Siemens"]
    assert min(tes) > 75


def test_K0_delta_equals_the_stored_S_M9_7_slope_and_differs_when_the_factor_changes():
    assert K.k0()["ok"] and not K.k0(1.10)["ok"]


def test_Ka_noise_free_algebra_and_failure_of_a_generator_without_the_TE_term():
    assert K.ka()["ok"] and not K.ka(te_term=False)["ok"]


def test_Ke_ComBat_apply_equals_neuroCombat_and_differs_without_EB():
    assert K.ke()["ok"] and not K.ke(eb=False)["ok"]


def test_without_noise_M2_M3_M4_give_the_algebraic_values():
    c = dict(idx=0, D="D1", S="S3e", V=0.0, H="H0", N="unbal")
    d = M.generate(c, np.random.default_rng(0), noise=False)
    y = d["y"] + 1e-6 * np.random.default_rng(1).standard_normal(d["y"].shape)   # keeps the residual variance from being 0
    k = len(M.ACQ)
    th, dl = d["theta"], M.delta("e")
    e2 = A.m2(y, d["g"], d["site"], k)["est"]
    m3 = A.m3(y, d["g"], d["site"], k, d["te"], M.TE_C)
    m4 = A.m4(y, d["g"], d["site"], k, d["te_s"])
    assert np.allclose(m3["te0"]["est"], th, atol=1e-5) and np.allclose(m3["slope"]["est"], dl, atol=1e-7)
    assert np.allclose(m4["te0"]["est"], th, atol=1e-5) and np.allclose(m4["slope"]["est"], dl, atol=1e-7)
    assert np.all(np.abs(e2 - th) > 0.05)          # a site fixed effect alone leaves the within-site bias (the algebra of Proposition 4)


def test_failure_control_without_relaxation_assignment_M2_also_gives_theta():
    c = dict(idx=0, D="D1", S="S1e", V=0.0, H="H0", N="bal")
    d = M.generate(c, np.random.default_rng(0), noise=False)
    y = d["y"] + 1e-6 * np.random.default_rng(1).standard_normal(d["y"].shape)
    assert np.allclose(A.m2(y, d["g"], d["site"], len(M.ACQ))["est"], d["theta"], atol=1e-5)


def test_coverage_and_rejection_decisions():
    e = dict(est=np.array([0.0, 1.0]), se=np.array([1.0, 0.1]), df=1e6)
    assert A.covered(e, 0.0).tolist() == [True, False] and A.rejected(e).tolist() == [False, True]


@pytest.mark.parametrize("scen,expected", [("S1b", "b"), ("S1e", "e")])
def test_S1_size_equals_the_apparent_difference_of_the_relaxation_scenario_at_the_central_TE(scen, expected):
    assert M.theta_of(scen) == pytest.approx(M.TE_C * M.delta(expected), rel=0, abs=1e-15)


def test_in_D0_M3_and_M4_are_recorded_as_not_defined_and_do_not_stop():
    from design_sim import run as R
    c = dict(idx=0, D="D0", S="S2e", V=0.0, H="H0", N="bal")
    recs = [R.one(c, np.random.default_rng(s)) for s in (0, 1)]
    assert all(r[k] is None for r in recs for k in ("M3_te0", "M3_tec", "M4_te0", "M3_slope", "M4_slope", "M4_tau0"))
    met = R.metrics(recs)
    assert met["M3_te0"] == {"defined": False} and met["M4_tau0_share"] is None and "coverage" in met["M2"]


def test_failure_control_with_all_TEs_equal_the_M4_TE_regression_stops_on_a_singular_matrix():
    d = M.generate(dict(idx=0, D="D0", S="S0", V=0.0, H="H0", N="bal"), np.random.default_rng(0))
    with pytest.raises(np.linalg.LinAlgError):
        A.m4(d["y"], d["g"], d["site"], len(M.ACQ), d["te_s"])


def test_one_D1_replicate_gives_every_analysis():
    from design_sim import run as R
    rec = R.one(dict(idx=0, D="D1", S="S1e", V=0.1, H="H2", N="unbal"), np.random.default_rng(0))
    assert all(rec[k] is not None for k in R.KEYS + R.SLOPES) and rec["M4_tau0"] is not None
