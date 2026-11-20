"""Tests that the acquisition-condition cards of rules_env.cards equal the values of the manuscript and the algebra cards (ON-Harmony pairs and TE range, E5, Wong TE grid, Tong card, design-letter map)."""
import pathlib
import sys

import numpy as np
import pytest

HERE = pathlib.Path(__file__).resolve().parent
for _cand in (HERE.parent, HERE.parent / "decomp"):
    if (_cand / "rules_env").is_dir():
        sys.path.insert(0, str(_cand))
        break
from rules_env import cards as K  # noqa: E402


def test_onh_pairs_count():
    assert len(K.pairs("onh28")) == 28


def test_onh_te_range():
    tes = [K.card_meta(f"ONH:{u}")["TE"] for u in K.ONH_UNITS]
    assert (min(tes), max(tes)) == (68.0, 98.0)          # N-250


def test_scheme_differs_15():
    n = sum(K.card_meta(a)["scheme"] != K.card_meta(b)["scheme"] for a, b in K.pairs("onh28"))
    assert n == 15                                          # N-116


def test_e5_matches_config():
    q2 = K.card_meta("E5:Q2")
    assert (q2["TE"], q2["TR"], q2["scheme"]) == (70.6, 10215.0, "P2b")
    assert [p for p in K.pairs("e5")] == [("E5:A", f"E5:{q}") for q in ("Q1", "Q2", "Q3", "Q4")]


def test_design_letter_map():
    # Agreement of the arrays (bval · bvec) is checked by check_synthfit (a) in envs/reference — the P2a · P2b files are outside the repository (same place as the stage 1 code)
    a, o = K.card_meta("E5:A"), K.card_meta("ONH:OXF1PRI")
    assert (o["TE"], o["TR"], o["scheme"]) == (a["TE"], a["TR"], a["scheme"])
    assert K.card_meta("ONH:NOT1ACH")["scheme"] == "P2b" and K.card_meta("ONH:NOT1ACH")["TE"] == pytest.approx(70.6)


def test_card_meta_has_no_arrays():
    assert set(K.card_meta("ONH:NOT4GEP")) == {"name", "TE", "TR", "scheme"}


def test_wong_te_grid():
    cs = [K.card(f"WONG:TE{t}") for t in range(62, 133, 10)]
    assert [c["TE"] for c in cs] == [62.0, 72.0, 82.0, 92.0, 102.0, 112.0, 122.0, 132.0]
    assert all(c["TR"] == 5800.0 for c in cs)
    assert all(np.array_equal(c["bvals"], cs[0]["bvals"]) for c in cs)
    assert cs[0]["bvecs"].shape == (64, 3) and int((cs[0]["bvals"] < 50).sum()) == 4
    assert K.pairs("wong_end") == [("WONG:TE62", "WONG:TE132")]


def test_tong_card():
    t = K.card("TONG")
    assert (t["TE"], t["TR"]) == (71.0, 5400.0)
    assert t["bvecs"].shape == (66, 3) and int((t["bvals"] < 50).sum()) == 6


def test_unknown_card():
    with pytest.raises(KeyError):
        K.card("ONH:NOPE")


def test_onh_units_match_analysis():
    from rules_env._paths import use
    use("decomp", "remeasure")
    import analysis as AN
    assert K.ONH_UNITS == AN.CARD_UNITS


def test_onh_cards_match_fiso_sign_pairs():
    # Same TE · TR as the algebra card · FISO sign rule (fixed/fiso_sign.pairs — NOT2ING uses the confirmed TR range)
    from rules_env._paths import use
    use("decomp", "fixed")
    import fiso_sign as FS
    for a, b, dsg, teA, trA, teB, trB in FS.pairs():
        ca, cb = K.card_meta(f"ONH:{a}"), K.card_meta(f"ONH:{b}")
        assert (ca["TE"], ca["TR"], cb["TE"], cb["TR"]) == pytest.approx((teA, trA, teB, trB), abs=1e-9)
        assert (ca["scheme"] != cb["scheme"]) == dsg
