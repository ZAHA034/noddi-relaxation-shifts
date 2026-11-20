"""Tests, on fake data, that the scale transforms, retest SDs and scanner values of rules_env.empirical equal those of the existing verdict code (judge_tong, judge_onh, s_ref_f2)."""
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
from rules_env import empirical as EM  # noqa: E402

JT, JO, AN = EM.JT, EM.JO, EM.AN


def _roi(rows):
    return {"주": {u: dict(NDI=n, FWF=f, ODI=o) for u, (n, f, o) in rows.items()}}


def test_scale_values_match_unit_coords():
    roi = _roi({"a": (0.6, 0.08, 0.2), "b": (0.7, 0.0, 0.3), "c": (0.5, 0.15, 0.25)})
    uc = JT.unit_coords(roi)
    for k, ci in (("x", 0), ("y", 1), ("odi", 2)):
        got = EM.scale_values(roi, k)
        assert got == pytest.approx({u: float(v[ci]) for u, v in uc.items()}, abs=0)
    assert EM.scale_values(roi, "raw_ndi") == {"a": 0.6, "b": 0.7, "c": 0.5}
    assert EM.scale_values(roi, "logit_fiso")["a"] == pytest.approx(float(AN.logit_fiso(0.08)))


def test_scale_values_skip_incomplete():
    roi = _roi({"a": (0.6, 0.08, None), "b": (0.7, 0.1, 0.3)})
    for k in EM.SCALES:
        assert set(EM.scale_values(roi, k)) == {"b"}


def test_clipped_units():
    roi = _roi({"a": (0.6, 0.0, 0.2), "b": (1.0, 0.1, 0.3), "c": (0.6, 0.1, 0.2)})
    assert EM.clipped(roi) == {"a", "b"}


def test_retest_sd_scalar_matches_judge_tong():
    rng = np.random.default_rng(0)
    rep = {sb: [{u: np.array([rng.normal(), rng.normal(), rng.normal()]) for u in "abc"} for _ in range(3)] for sb in ("s1", "s2")}
    ref = JT.retest_sd(rep, list("abc"))["x"]
    got = EM.retest_sd_scalar({sb: [{u: float(v[u][0]) for u in v} for v in sess] for sb, sess in rep.items()}, list("abc"))
    assert got == pytest.approx(ref)


def test_s_rt_corrected_removes_mean_shift():
    rng = np.random.default_rng(1)
    units = [f"u{i}" for i in range(400)]
    orig62, orig132 = ({u: rng.normal() for u in units} for _ in range(2))
    ser = {"TE62": orig62, "TE132": orig132,
           "TE62R2": {u: orig62[u] + 0.3 + rng.normal(0, 0.05 * math.sqrt(2)) for u in units},
           "TE132R2": {u: orig132[u] - 0.2 + rng.normal(0, 0.05 * math.sqrt(2)) for u in units}}
    s, df = EM.s_rt_corrected(ser)
    assert s == pytest.approx(0.05, rel=0.1) and df == 800 - 2
    legacy, _ = JW_s_rt(ser)
    assert legacy > 3 * s


def JW_s_rt(ser):
    return EM.JW.s_rt({k: {u: np.array([v, 0.0, 0.0]) for u, v in d.items()} for k, d in ser.items()}, 0)


def test_onh_device_values_match_judge_onh():
    rois = {f"s{i}_ses-OXF1PRI001": _roi({"a": (0.6, 0.08, 0.2 + 0.01 * i), "b": (0.7, 0.1, 0.3)}) for i in range(3)}
    by_unit = {"OXF1PRI": {f"s{i}": "ses-OXF1PRI001" for i in range(3)}}
    ref = JO.device_values(by_unit, rois, "ODI", "raw", "주")
    assert EM.onh_device_values(by_unit, rois, "odi") == ref


def test_s_onh_matches_s_ref_f2_for_odi():
    R = EM.R
    rng = np.random.default_rng(2)
    rois, states = {}, {}
    for sb in R.SREF_SUBJECTS:
        for ses in R.SREF_SESSIONS:
            rois[f"{sb}_{ses}"] = _roi({u: (0.6, 0.08, 0.2 + rng.normal(0, 0.01)) for u in "abc"})
            states[f"{sb}_{ses}"] = EM.V4.USABLE
    ref_s, ref_df, _ = JO.s_ref_f2(pathlib.Path("."), states, rois, list("abc"))
    s, df = EM.s_onh(states, rois, list("abc"), "odi")
    assert (s, df) == (pytest.approx(ref_s), ref_df)
