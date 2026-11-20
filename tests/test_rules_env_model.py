"""Tests that the algebra of rules_env.model (formulas (3), (6), (7)) reproduces the values of the manuscript and the ledger, including Fig. 2, Proposition 4, the rows of Supplementary Table S11 and the inversion round trip."""
import json
import math
import pathlib
import sys

import numpy as np
import pytest

HERE = pathlib.Path(__file__).resolve().parent
for _cand in (HERE.parent, HERE.parent / "decomp"):          # drafts: x15/env/src · repository: <root>/decomp
    if (_cand / "rules_env").is_dir():
        sys.path.insert(0, str(_cand))
        break
from rules_env import model as M  # noqa: E402

CTRL = {"t2": (82.5, 52.0, 1000.0), "t1": None}               # Section 2.4 · Figure 3 control group
PAT = {"t2": (82.5, 70.0, 1000.0), "t1": None}                # Figure 3 patient (T2,en 70 ms)


def test_fig3_values():
    d70 = M.apparent(0.46, 0.0, M.weights(70.0, CTRL["t2"]))[0] - M.apparent(0.46, 0.0, M.weights(70.0, PAT["t2"]))[0]
    d100 = M.apparent(0.46, 0.0, M.weights(100.0, CTRL["t2"]))[0] - M.apparent(0.46, 0.0, M.weights(100.0, PAT["t2"]))[0]
    assert d70 == pytest.approx(0.0857477752721304, abs=1e-5)       # N-649
    assert d100 == pytest.approx(0.12025460115930686, abs=1e-5)     # N-650
    assert M.prop4_change({"TE": 100.0}, {"TE": 70.0}, CTRL, PAT) == pytest.approx(-0.14835164835164855, abs=1e-6)  # N-106


def test_prop4_t1_term_vanishes_equal_tr():
    c = {"t2": CTRL["t2"], "t1": (1100.0 * math.sqrt(1.3), 1100.0 / math.sqrt(1.3), 4300.0)}
    p = {"t2": PAT["t2"], "t1": (c["t1"][0] * 1.0441, c["t1"][1] * 1.0441 * 1.2, 4300.0)}
    te_only = M.prop4_change({"TE": 100.0}, {"TE": 70.0}, CTRL, PAT)
    same_tr = M.prop4_change({"TE": 100.0, "TR": 3600.0}, {"TE": 70.0, "TR": 3600.0}, c, p)
    diff_tr = M.prop4_change({"TE": 100.0, "TR": 3600.0}, {"TE": 70.0, "TR": 10215.0}, c, p)
    assert same_tr == pytest.approx(te_only, abs=1e-12)
    assert abs(diff_tr - te_only) > 1e-4


def test_prop4_pd_cancels():
    p = {"t2": CTRL["t2"], "t1": None, "rho": (1.2, 1.0, 1.0)}
    assert M.prop4_change({"TE": 100.0}, {"TE": 70.0}, CTRL, p) == pytest.approx(0.0, abs=1e-12)
    assert M.x_shift({"TE": 70.0}, p) - M.x_shift({"TE": 70.0}, CTRL) == pytest.approx(math.log(1.2), abs=1e-12)


def test_translation_independent_of_fractions():
    w = M.weights(92.0, (82.0, 51.25, 1000.0), 3600.0, (1254.0, 965.0, 4300.0))
    dx, dy = M.translation(w)
    for nu, f in ((0.3, 0.02), (0.6, 0.10), (0.75, 0.25)):
        x0, y0 = M.coords(nu, f)
        x1, y1 = M.coords(*M.apparent(nu, f, w))
        assert (x1 - x0, y1 - y0) == pytest.approx((dx, dy), abs=1e-12)


def test_invert_roundtrip():
    w = M.weights(68.0, (95.0, 95.0 / 1.8, 1000.0))
    for nu in np.linspace(0.1, 0.9, 9):
        for f in (0.01, 0.05, 0.2, 0.5):
            assert M.invert(*M.apparent(nu, f, w), w) == pytest.approx((nu, f), abs=1e-12)


def test_table_s15_rows():
    # Values of x15/ledger/derived/residual_s15.json (drafts): TE_B 70 · TE_A 80 / 100 · control 82.5 / 52
    en = {"t2": (82.5, 52.0 * 1.0984, 1000.0), "t1": None}
    both = {"t2": (82.5 * 1.0984, 52.0 * 1.0984, 1000.0), "t1": None}
    assert M.prop4_change({"TE": 80.0}, {"TE": 70.0}, CTRL, en) == pytest.approx(-0.017227855902291478, abs=1e-6)
    assert M.prop4_change({"TE": 100.0}, {"TE": 70.0}, CTRL, en) == pytest.approx(-0.051683567706874434, abs=1e-6)
    assert M.prop4_change({"TE": 80.0}, {"TE": 70.0}, CTRL, both) == pytest.approx(-0.0063690861214532254, abs=1e-6)


def test_load_input_requires_keys(tmp_path):
    good = {"acquisitions": [], "tissues": [], "groups": {}}
    (tmp_path / "ok.json").write_text(json.dumps(good))
    assert M.load_input(tmp_path / "ok.json")["groups"] == {}
    del good["groups"]
    (tmp_path / "bad.json").write_text(json.dumps(good))
    with pytest.raises(ValueError, match="groups"):
        M.load_input(tmp_path / "bad.json")
