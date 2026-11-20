"""Known-answer tests on fake data for the ON-Harmony pair test loader and the end markers of the Tong and Wongkornchaovalit loaders, with no real images."""
from __future__ import annotations

import json
import math
import pathlib
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "decomp" / "remeasure"))
import judge_onh as JO  # noqa: E402
import judge_tong as JT  # noqa: E402
import judge_wong as JW  # noqa: E402
import rconfig as RC  # noqa: E402
import roi as ROI  # noqa: E402

UNITS = ROI.label_units()
SESSIONS = RC.onh_sessions() + RC.sref_extra()


def _roi(f, o):
    row = lambda fw, od: dict(NDI=0.5, FWF=fw, ODI=od, MD=0.8, FA=0.4, s_low=0.5, s_2000=0.3)  # noqa: E731
    return {"주": {u: row(f[u], o[u]) for u in UNITS}, "민감도_옛_-M": {u: row(f[u], o[u] + 1e-3) for u in UNITS}}


def _work(tmp_path, monkeypatch, dev_f=None, dev_o=None, fail=(), clip=None, seed=0):
    """Fake working folder — 120 sessions of roi.json · DONE.json · batch summary · access record · guard basis. dev_f · dev_o = scanner effects (logit FWF · ODI)."""
    import guard as G
    import run as RUN
    rng = np.random.default_rng(seed)
    dev_f, dev_o = dev_f or {}, dev_o or {}
    subs = sorted({a for a, _ in SESSIONS})
    sf = {s: rng.normal(0, 0.3) for s in subs}
    so = {s: rng.normal(0, 0.02) for s in subs}
    bf = {u: rng.normal(-2.0, 0.5) for u in UNITS}
    bo = {u: rng.uniform(0.1, 0.4) for u in UNITS}
    root = tmp_path / "onharmony"
    for sub, ses in SESSIONS:
        d = RC.device_of(ses)
        wd = root / f"{sub}_{ses}"
        wd.mkdir(parents=True)
        f = {u: 1 / (1 + math.exp(-(bf[u] + sf[sub] + dev_f.get(d, 0.0) + rng.normal(0, 0.05)))) for u in UNITS}
        o = {u: bo[u] + so[sub] + dev_o.get(d, 0.0) + rng.normal(0, 0.002) for u in UNITS}
        if clip and (sub, ses) == clip[0]:
            f[clip[1]] = 0.0
        (wd / "roi.json").write_text(json.dumps(_roi(f, o)))
        if (sub, ses) in fail:
            (wd / "FAILED.json").write_text("{}")
        else:
            (wd / "DONE.json").write_text("{}")
    (root / "batch_summary.json").write_text(json.dumps({"실패_다시도_실패": [f"{a}_{b}" for a, b in fail], "충돌": []}))
    log = tmp_path / "access.json"
    log.write_text(json.dumps({"열람": []}))
    monkeypatch.setattr(RC, "ACCESS_LOG_ONH", log)
    monkeypatch.setattr(RUN, "LOCK", tmp_path / "access.lock")
    monkeypatch.setattr(G, "require", lambda: dict(fixed_commit="f" * 40, pipeline_commit="p" * 40))
    return log


def _pairs(res, key="main"):
    return {p["pair"]: p for p in (res[key] if key == "main" else res["sensitivity"][key])["pairs"]}


def test_pair_value_is_D_j_minus_D_i_gives_19_verdict_pairs_and_s_ref_F2(tmp_path, monkeypatch):
    log = _work(tmp_path, monkeypatch, dev_f={"NOT1ACH": 1.0})
    res, sha = JO.run(tmp_path, n_sim_delta=10, n_sim_size=10)
    P = _pairs(res)
    assert P["NOT1ACH~OXF1PRI"]["fiso"] < -0.8 and P["NOT2ING~NOT4GEP"]["fiso"] == pytest.approx(0.0, abs=0.1)   # NOT1ACH is i (first) — D_j − D_i < 0
    assert sum(p["judgeable_fiso"] for p in P.values()) == 19 and not P["NOT4GEP~OXF3TRI"]["judgeable_fiso"]
    assert res["main"]["fiso_null_size"] == JO.FISO_SIZE_FIXED and res["main"]["odi_a"]["delta_tau"] == JO.DELTA_TAU_FIXED
    assert res["main"]["odi_a"]["descriptive_only"] and res["main"]["odi_a_delta_tau_redone"] is None
    assert res["main"]["fiso_sign"]["n_perm"] == 5040 and res["main"]["fiso_sign"]["n_pairs"] == 8
    assert len(res["sessions"]["s_ref_F2_used"]) == 12 and res["s_ref_F2"]["df"] == 10 and res["s_ref_F2"]["multipliers"][1] == pytest.approx(1.7549, abs=1e-4)
    assert res["sessions"]["main_excluded"] == ["sub-16794_ses-NOT4GEP001"]
    import analysis as AN
    import delta_tau as DT
    cnt = AN.common_counts(RC.onh_sessions())
    keep = [k for k, p in enumerate(res["main"]["pairs"]) if p["judgeable_odi"]]
    y = np.array([abs(p["odi"]) if p["odi"] is not None else np.nan for p in res["main"]["pairs"]])
    T = DT.perm_taus_x(DT.pair_predictions_weighted(cnt, "L1"), y, keep, DT.device_perm_index())
    a = res["main"]["odi_a"]
    assert a["tau"] == pytest.approx(T[0]) and a["p_tau0_two_sided"] == pytest.approx(float((np.abs(T) >= abs(T[0]) - 1e-12).mean()))   # Both sides (a16)
    ent = json.loads(log.read_text())["열람"]
    assert ent[-1]["파일"][0]["sha256"] == sha
    with pytest.raises(RuntimeError, match="이미 있다"):
        JO.run(tmp_path, n_sim_delta=10, n_sim_size=10)


def test_ODI_b_breaks_only_big_scanner_effect_pairs_notes_acquisition_cell(tmp_path, monkeypatch):
    _work(tmp_path, monkeypatch, dev_o={"OXF4GEP": 0.1})
    res, _ = JO.run(tmp_path, n_sim_delta=10, n_sim_size=10)
    b = res["main"]["odi_b"]
    want = {"NOT1ACH~OXF4GEP", "NOT2ING~OXF4GEP", "NOT4GEP~OXF4GEP", "OXF1PRI~OXF4GEP", "OXF2PRI~OXF4GEP"}   # OXF3TRI~OXF4GEP are pairs for which no verdict is possible
    assert {x["pair"] for x in b["broken"]} == want and b["n_pairs"] == 19 and b["n_retest_level"] == 14
    assert all("dTE" in x and "dTR_ms" in x for x in b["broken"])
    assert set(res["main"]["odi_b_sensitivity"]) == {"옛_s_ref_5262f43", "s_ref_F2_x1.7549"}


def test_under_3_common_subjects_after_failure_no_verdict_design_redone(tmp_path, monkeypatch):
    tri = [(a, b) for a, b in RC.onh_sessions() if RC.device_of(b) == "OXF3TRI"]
    _work(tmp_path, monkeypatch, fail=tri[:8])
    res, _ = JO.run(tmp_path, n_sim_delta=10, n_sim_size=10)
    moved = res["main"]["moved_to_not_judgeable"]["fiso"]
    assert {"NOT1ACH~OXF3TRI", "OXF1PRI~OXF3TRI", "OXF2PRI~OXF3TRI"} <= set(moved)
    assert res["main"]["odi_a_delta_tau_redone"] is not None and res["main"]["fiso_null_size"] != JO.FISO_SIZE_FIXED
    assert len(res["sessions"]["dropped"]) == 8 and all(v == "다시도 실패" for v in res["sessions"]["dropped"].values())


def test_five_sens_16794_in_NOT2ING_98_only_risk_ROI_out_old_tab_raw_scale(tmp_path, monkeypatch):
    _work(tmp_path, monkeypatch, clip=(("sub-13192", "ses-NOT1ACH001"), UNITS[0]))
    res, _ = JO.run(tmp_path, n_sim_delta=10, n_sim_size=10)
    assert set(res["sensitivity"]) == {"fiso_raw_scale", "with_16794_NOT4GEP001", "not2ing_98_only", "without_risk_roi", "old_nonzero_table"}
    P, W = _pairs(res), _pairs(res, "with_16794_NOT4GEP001")
    assert W["NOT2ING~NOT4GEP"]["n_subjects_fiso"] == P["NOT2ING~NOT4GEP"]["n_subjects_fiso"] + 1          # 16794 is dropped from the main analysis
    N = _pairs(res, "not2ing_98_only")
    assert P["NOT1ACH~NOT2ING"]["n_subjects_fiso"] == 20 and N["NOT1ACH~NOT2ING"]["n_subjects_fiso"] == 16
    assert res["sensitivity"]["not2ing_98_only"]["odi_a_delta_tau_redone"] is not None                      # The prediction changes, so the design is redone
    assert all(r["n_units"] == 23 for r in res["sensitivity"]["without_risk_roi"]["odi_b"]["per_pair"].values())
    assert res["clipped"]["by_session"] == {"sub-13192_ses-NOT1ACH001": [UNITS[0]]}
    assert res["clipped"]["by_pair"]["NOT1ACH~OXF1PRI"] == 1


def test_no_batch_summary_stops_missing_repeats_lower_s_ref_F2_dof(tmp_path, monkeypatch):
    _work(tmp_path, monkeypatch, fail=[("sub-16974", "ses-OXF1PRI004")])
    res, _ = JO.run(tmp_path, n_sim_delta=10, n_sim_size=10)
    assert len(res["sessions"]["s_ref_F2_used"]) == 11 and res["s_ref_F2"]["df"] == 9
    (tmp_path / "onharmony" / "pairs" / "onh_pairs.json").unlink()
    monkeypatch.setattr(JO, "OLD_SREF_SHA", "0" * 64)
    with pytest.raises(RuntimeError, match="sha256"):                                            # Stops if the old s_ref file changes
        JO.run(tmp_path, n_sim_delta=10, n_sim_size=10)
    (tmp_path / "onharmony" / "batch_summary.json").unlink()
    with pytest.raises(RuntimeError, match="batch_summary"):
        JO.run(tmp_path, n_sim_delta=10, n_sim_size=10)


def _guard_ok(monkeypatch):
    import guard as G
    monkeypatch.setattr(G, "require", lambda: {})


def test_Tong_loader_returns_without_sessions_lacking_end_marker(tmp_path, monkeypatch):
    _guard_ok(monkeypatch)
    row = {u: dict(NDI=0.5, FWF=0.1, ODI=0.2) for u in UNITS}
    for i, ses in enumerate(("ses-c01r1", "ses-c02r1")):
        wd = tmp_path / "tong" / f"sub-1_{ses}"
        wd.mkdir(parents=True)
        (wd / "roi.json").write_text(json.dumps({"주": row}))
        if i == 0:
            (wd / "DONE.json").write_text("{}")
    values, rep, missing = JT.load(tmp_path)
    assert "sub-1" in values["c01"] and "c02" not in values and "sub-1_ses-c02r1" in missing


def test_Wong_loader_drops_unmarked_subjects_separate_reg_set_not_empty(tmp_path, monkeypatch):
    _guard_ok(monkeypatch)
    row = {"주": {u: dict(NDI=0.5, FWF=0.1, ODI=0.2) for u in UNITS}}
    for sb, done in (("sub-01", True), ("sub-02", False)):
        wd = tmp_path / "wong" / sb
        (wd / "own_TE72").mkdir(parents=True)
        (wd / "roi_TE62.json").write_text(json.dumps(row))
        (wd / "own_TE72" / "roi_TE72.json").write_text(json.dumps(row))
        if done:
            (wd / "DONE.json").write_text(json.dumps({"따로_정합": ["TE72"]}))
    data, own, missing = JW.load(tmp_path)
    assert set(data) == {"sub-01"} and own["sub-01"] == {"TE72"} and "TE72" in data["sub-01"] and missing == ["sub-02", "sub-03"]   # sub-03 does not even have a folder
