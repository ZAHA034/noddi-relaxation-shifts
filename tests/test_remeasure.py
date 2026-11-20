"""Known-answer tests of the measurement pipeline: loader guard, commands matching the old pipeline, 27 ROI units, session lists and access record, without images or FSL."""
from __future__ import annotations

import hashlib
import json
import os
import pathlib
import subprocess
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "decomp" / "remeasure"))
import cmds as K  # noqa: E402
import fit_f2 as F  # noqa: E402
import guard as G  # noqa: E402
import rconfig as R  # noqa: E402
import roi  # noqa: E402
import run  # noqa: E402

TIMING = ROOT / "prereg" / "scripts" / "lockd" / "timing_run.sh"


# ── ⑨-b loader guard ────────────────────────────────────────────────────────────
def _git(tmp, *a):
    subprocess.run(["git", "-C", str(tmp), *a], check=True, capture_output=True)


def _repo(tmp, with_fixed=False):
    _git(tmp, "init", "-q")
    _git(tmp, "config", "user.name", "t")
    _git(tmp, "config", "user.email", "t@t")
    for rel in G.PIPELINE:
        (tmp / rel).parent.mkdir(parents=True, exist_ok=True)
        (tmp / rel).write_text("x\n")
    if with_fixed:
        (tmp / G.FIXED_DOC).write_text("고정\n")
    _git(tmp, "add", "-A")
    _git(tmp, "commit", "-qm", "pipeline")


def _fix(tmp):
    (tmp / G.FIXED_DOC).write_text("고정\n")
    _git(tmp, "add", "-A")
    _git(tmp, "commit", "-qm", "fixed")


def test_guard_stops_without_fixed_analysis_document(tmp_path):
    _repo(tmp_path)
    with pytest.raises(G.GuardError, match="HEAD 에 없다"):
        G.check(tmp_path)


def test_guard_stops_if_fixed_analysis_document_working_copy_differs(tmp_path):
    _repo(tmp_path)
    _fix(tmp_path)
    (tmp_path / G.FIXED_DOC).write_text("고침\n")
    with pytest.raises(G.GuardError, match="고정 문서 .* 작업 사본"):
        G.check(tmp_path)


def test_guard_same_commit_as_pipeline_is_not_after(tmp_path):
    _repo(tmp_path, with_fixed=True)
    with pytest.raises(G.GuardError, match="뒤가 아니다"):
        G.check(tmp_path)


def test_guard_blocks_again_if_code_is_committed_after_fixing(tmp_path):
    _repo(tmp_path)
    _fix(tmp_path)
    (tmp_path / G.PIPELINE[0]).write_text("y\n")
    _git(tmp_path, "commit", "-qam", "code")
    with pytest.raises(G.GuardError, match="뒤가 아니다"):
        G.check(tmp_path)


def test_guard_does_not_run_with_uncommitted_code(tmp_path):
    _repo(tmp_path)
    _fix(tmp_path)
    (tmp_path / G.PIPELINE[-1]).write_text("y\n")
    with pytest.raises(G.GuardError, match="파이프라인 파일 .* 작업 사본"):
        G.check(tmp_path)


def test_guard_ok_returns_three_basis_commits(tmp_path):
    _repo(tmp_path)
    _fix(tmp_path)
    b = G.check(tmp_path)
    assert b["fixed_commit"] == b["head"] and b["pipeline_commit"] != b["fixed_commit"]


def test_real_repo_blocked_no_fixed_analysis_document_in_HEAD_else_basis():
    """Blocked before ⑩, open after the ⑩ commit (⑪ pre-check — so that the full test also runs after ⑩)."""
    root = pathlib.Path(__file__).resolve().parents[1]
    in_head = subprocess.run(["git", "cat-file", "-e", f"HEAD:{G.FIXED_DOC}"], cwd=root, capture_output=True).returncode == 0
    if not in_head:
        with pytest.raises(G.GuardError):
            G.require()
    else:
        b = G.require()
        assert b["fixed_commit"] and b["pipeline_commit"] and b["fixed_commit"] != b["pipeline_commit"]


# ── ⑨-a command — same as the old registered path (timing_run.sh)? ───────────────────────────
def _old(step_name):
    line = next(ln for ln in TIMING.read_text(encoding="utf-8").splitlines() if ln.startswith(f"step {step_name}"))
    return line.split(None, 2)[2]


def test_eddy_has_old_options_differing_only_in_executable_and_threads():
    old = _old("04_eddy").split()
    new = K.eddy(8)
    assert old[0] == "eddy_cuda" and pathlib.Path(new[0]).name == "eddy_cpu"
    assert new[1:-1] == old[1:] and new[-1] == "--nthr=8"


def test_topup_bet_registration_matches_old_pipeline():
    old = _old("02_topup").replace('"$FSLDIR/etc/flirtsch/b02b0.cnf"', str(R.FSLDIR / R.TOPUP_CONFIG)).split()
    new = K.topup()
    assert [pathlib.Path(new[0]).name] + new[1:] == old
    assert "bet hifi_b0_mean nodif_brain -m -f 0.25" in _old("03_mask")
    assert [pathlib.Path(a[0]).name for a in K.mask()] == ["fslmaths", "bet"] and K.mask()[1][1:] == ["hifi_b0_mean", "nodif_brain", "-m", "-f", "0.25"]
    old_reg = _old("08_reg_fa").replace("\\$FSLDIR/data/standard/FMRIB58_FA_1mm", str(R.FSLDIR / R.FA_REF))
    flirt, fnirt = K.register()
    assert " ".join([pathlib.Path(flirt[0]).name] + flirt[1:]) in old_reg
    assert " ".join([pathlib.Path(fnirt[0]).name] + fnirt[1:]) in old_reg


def test_acqparams_follows_sidecar_direction():
    assert K.acqparams("j-", "j", 0.069) == "0 -1 0 0.069\n0 1 0 0.069\n"                      # Same as the old fixed value
    assert K.acqparams("j", "j-", 0.07) == "0 1 0 0.07\n0 -1 0 0.07\n"                         # Both flipped — the same adjustment in topup and eddy
    assert K.acqparams("j-", "i", 0.04).splitlines()[1] == "1 0 0 0.04"                         # Sessions in which PA is left-right (README)
    assert K.first_b0(["5", "1000", "0"]) == 0 and K.first_b0(["1000", "50", "0"]) == 1


# ── ROI ──────────────────────────────────────────────────────────────────────


def _fake_space():
    atlas = np.zeros((4, 4, 3), int)
    atlas[0, :, :] = 7
    atlas[1, :, :] = 8
    atlas[2, :, :] = 1
    skel = np.full(atlas.shape, 3000)
    skel[:, 3, :] = 0                                                         # Outside the skeleton
    return atlas, skel


def test_label_mean_main_rule_includes_0_and_old_rule_excludes_it():
    atlas, skel = _fake_space()
    m = np.ones(atlas.shape)
    m[0, 0, :] = 0.0
    main = roi.label_means({"FWF": m}, atlas, skel)
    old = roi.label_means({"FWF": m}, atlas, skel, nonzero=True)
    assert main["n"][6] == 9 and main["FWF"][6] == pytest.approx(6 / 9)                     # Label 7: 9 voxels inside the skeleton · includes three 0 values
    assert old["FWF"][6] == pytest.approx(1.0) and old["n_FWF"][6] == 6
    brain = np.ones(atlas.shape)
    brain[0, 1, :] = 0
    assert roi.label_means({"FWF": m}, atlas, skel, roi.valid_mask(brain=brain))["n"][6] == 6


def test_valid_mask_is_brain_F2_ext_valid_intersection_keeping_inner_0():
    """Relay ZZZZe-2 (2) — outside the brain · outside the registration · fit failures are dropped by the mask, and true zeros inside the mask (FWF = 0) are included · the seven maps use the same voxels."""
    atlas, skel = _fake_space()
    fwf = np.full(atlas.shape, 0.2)
    fwf[0, 0, :] = 0.0                                                        # True 0 values inside the mask (three)
    md = np.ones(atlas.shape)
    brain, f2v, extv = (np.ones(atlas.shape) for _ in range(3))
    brain[0, 1, 0] = 0                                                        # Outside the brain 1
    f2v[0, 2, 0] = 0.5                                                        # F2 fit failures (trilinear < 0.99) mixed in: 1
    extv[0, 2, 1] = 0.0                                                       # Extended-coordinate failure 1
    v = roi.valid_mask(brain, f2v, extv)
    lm = roi.label_means({"FWF": fwf, "MD": md}, atlas, skel, v)
    assert lm["n"][6] == 6 and lm["n_FWF"][6] == lm["n_MD"][6] == 6                       # 9 − 3 · two maps, same voxels
    assert lm["FWF"][6] == pytest.approx((0.2 * 3) / 6)                                    # mean including the three 0 values (dropping the 0 pushes it up to 0.2)
    old = roi.label_means({"FWF": fwf}, atlas, skel, None, nonzero=True)
    assert old["FWF"][6] == pytest.approx(0.2) and old["n_FWF"][6] == 6                    # Sensitivity (old -M) — drops the 0
    assert roi.valid_mask() is None and roi.valid_mask(f2_valid=f2v).sum() == f2v.size - 1


def test_F2_valid_map_is_only_finite_positive_mean_b0_in_mask():
    msk = np.ones((2, 2, 1), bool)
    msk[1, 1, 0] = False
    b0 = np.array([[[100.0], [0.0]], [[50.0], [80.0]]])
    ndi = np.full((2, 2, 1), 0.5)
    ndi[1, 0, 0] = np.nan
    ok = F.f2_valid(msk, b0, [ndi, np.zeros((2, 2, 1)), np.ones((2, 2, 1))])
    assert ok.tolist() == [[[True], [False]], [[False], [False]]]                          # b0 0 · NaN · outside the mask


def test_27_ROI_units_are_simple_LR_means_and_drop_empty_labels():
    lm = {"n": np.ones(48, int), "NDI": np.full(48, np.nan)}
    lm["NDI"][6], lm["NDI"][7], lm["NDI"][0] = 0.4, 0.6, 0.5                                   # Labels 7 · 8 = corticospinal tract left and right · label 1 = middle cerebellar peduncle
    t = roi.unit_table(lm, ("NDI",))
    assert t["Corticospinal tract"]["NDI"] == pytest.approx(0.5) and t["Middle cerebellar peduncle"]["NDI"] == pytest.approx(0.5)
    lm["NDI"][7] = np.nan
    assert roi.unit_table(lm, ("NDI",))["Corticospinal tract"]["NDI"] == pytest.approx(0.4)
    assert roi.unit_table(lm, ("NDI",))["Fornix (column and body of fornix)"]["NDI"] is None


# ── Session · per dataset ─────────────────────────────────────────────────────────────
def test_ON_Harmony_first_sessions_are_110_without_NOT3GEM():
    s = R.onh_sessions()
    assert len(s) == 110 and len(set(s)) == 110 and all(b.endswith("001") and "NOT3GEM" not in b for _, b in s)
    cnt = {}
    for _, b in s:
        cnt[R.device_of(b)] = cnt.get(R.device_of(b), 0) + 1
    assert cnt == {"NOT1ACH": 20, "NOT2ING": 20, "NOT4GEP": 10, "OXF1PRI": 20, "OXF2PRI": 20, "OXF3TRI": 10, "OXF4GEP": 10}


def test_Tong_main_shells_are_b0_1000_2000():
    assert run.select_vols(np.array([0, 1000, 2000, 3000, 5, 995, 2005, 2990]), R.TONG_MAIN_BMAX).tolist() == [0, 1, 2, 4, 5, 6]


def test_Wong_has_10_series_with_repeats_and_reference_is_first_TE62():
    assert len(R.WONG_SERIES) == 10 and R.WONG_REF == "TE62" and {"TE62R2", "TE132R2"} <= set(R.WONG_SERIES)


def test_alignment_metric_is_1_for_identical_maps():
    rng = np.random.default_rng(0)
    d = rng.uniform(100, 200, (5, 5, 4, 6))
    bv = np.array([0, 0, 700, 700, 2000, 2000])
    m = run.mean_dw_ratio(d, bv)
    mask = np.ones(m.shape, bool)
    assert run.align_r(m, m, mask) == pytest.approx(1.0) and run.align_r(m, m * 2 + 1, mask) == pytest.approx(1.0)


def test_ext_coords_same_when_chunked_drop_at_or_below_0_outlier_voxels():
    rng = np.random.default_rng(1)
    sig = rng.uniform(0.1, 1, (45001, 7))
    fn = lambda s, b, g, low: np.c_[s.sum(1), s[:, 0], s[:, 1] * low, s[:, -1]]
    out, ok = F.ext_maps(sig, None, None, 2.0, fn)
    assert ok.all() and np.allclose(out, fn(sig, None, None, 2.0), atol=1e-5)
    sig[5, 3] = 0.0                                                          # voxel with 0 mixed in (mask edge)
    sig[7, 0] = 0.123456                                                     # Outlier voxel marker

    def fn2(s, b, g, low):
        if np.any(np.isclose(s[:, 0], 0.123456)):
            raise np.linalg.LinAlgError("Singular matrix")
        return fn(s, b, g, low)
    out, ok = F.ext_maps(sig, None, None, 2.0, fn2)
    assert not ok[5] and not ok[7] and ok.sum() == len(sig) - 2
    assert (out[[5, 7]] == 0).all() and np.allclose(out[ok], fn(sig[ok], None, None, 2.0), atol=1e-5)


def test_invalid_ext_coord_voxel_drops_seven_maps_main_ext_only_sens():
    atlas, skel = _fake_space()
    md = np.ones(atlas.shape)
    md[0, 0, :] = 0.0                                                         # Invalid voxel (value 0 · valid 0)
    valid = np.ones(atlas.shape)
    valid[0, 0, :] = 0.0
    lm = roi.label_means({"MD": md, "NDI": md}, atlas, skel, roi.valid_mask(ext_valid=valid))
    assert lm["MD"][6] == pytest.approx(1.0) and lm["n_MD"][6] == 6 and lm["n_NDI"][6] == 6     # Main — same voxel set
    old = roi.label_means({"MD": md, "NDI": md}, atlas, skel, None, nonzero=True, ext_valid=valid)
    assert old["n_MD"][6] == 6 and old["n_NDI"][6] == 6                                    # Sensitivity — nonzero voxels


def test_series_below_alignment_threshold_get_separate_registration():
    """Relay ZZZZe-2 — a series with r < 0.90 (including NaN) is registered separately with the FA of that series · reported · sensitivity = that series left out."""
    assert run.own_registration({"TE62": 1.0, "TE72": 0.95, "TE82": 0.89, "TE92": float("nan"), "TE102": 0.90}) == ["TE82", "TE92"]
    assert R.WONG_ALIGN_R_MIN == 0.90


def test_main_excluded_ON_Harmony_sessions_are_in_first_110_with_reasons():
    assert set(R.MAIN_EXCLUDE_ONH) == {("sub-16794", "ses-NOT4GEP001")} and set(R.MAIN_EXCLUDE_ONH) <= set(R.onh_sessions())
    assert all("README:122" in why for why in R.MAIN_EXCLUDE_ONH.values())


# ── Real data verdict · access record ───────────────────────────────────────────────────


def test_access_record_appends_in_lock_with_relative_path_and_sha256(tmp_path):
    lg = tmp_path / "log.json"
    lg.write_text(json.dumps({"문서": "x"}))
    f = tmp_path / "data" / "a.nii.gz"
    f.parent.mkdir()
    f.write_bytes(b"abc")
    basis = dict(fixed_commit="f" * 40, pipeline_commit="p" * 40)
    for _ in range(2):
        G.log_access(lg, ("Tong 2020", "열람(새 분석)"), G.access_entry("시험", [f], tmp_path, basis, "s"), tmp_path / "lk")
    d = json.loads(lg.read_text())
    e = d["Tong 2020"]["열람(새 분석)"]
    assert len(e) == 2 and e[0]["파일"][0] == {"file": "data/a.nii.gz", "sha256": hashlib.sha256(b"abc").hexdigest()}
    assert "#51" in e[0]["근거"] and not any(str(v).startswith("/") for v in e[0]["파일"][0].values())


# ── Relay ZZZZf — scale · pair summary · retest SD · SNR · session ────────────────────────────
import analysis as AN  # noqa: E402
sys.path.insert(0, str(ROOT / "decomp" / "stage1"))
import lib  # noqa: E402


def test_FISO_logit_equals_S4_clipping_and_matches_y_of_coords():
    f = np.array([0.0, 1e-6, 0.2, 0.9, 1.0])
    n = np.full(5, 0.5)
    xy, _ = lib.coords(n, f, np.zeros(5))
    assert np.allclose(AN.logit_fiso(f), xy[:, 1] + np.log(0.5)) and lib.CLIP == 1e-4
    assert AN.logit_fiso(0.0) == pytest.approx(np.log(1e-4 / (1 - 1e-4)))


def test_ROI_unit_values_drop_units_without_values_and_logit_is_FISO_only():
    res = {"주": {"A": {"FWF": 0.0, "ODI": 0.3}, "B": {"FWF": None, "ODI": 0.4}}}
    assert AN.unit_values(res, "ODI") == {"A": 0.3, "B": 0.4}
    assert list(AN.unit_values(res, "FWF", "logit")) == ["A"]
    with pytest.raises(ValueError):
        AN.unit_values(res, "ODI", "logit")


def test_pair_summary_subject_ROI_unit_medians_under_3_common_no_verdict():
    A = {"s1": {"u1": 0.0, "u2": 0.0}, "s2": {"u1": 0.0, "u2": 0.0}, "s3": {"u1": 0.0}}
    B = {"s1": {"u1": 1.0, "u2": 5.0}, "s2": {"u1": 3.0, "u2": 7.0}, "s3": {"u1": 2.0, "u2": 9.0}, "s4": {"u1": 9.0}}
    r = AN.pair_summary(A, B, ["u1", "u2", "u3"])
    assert r["n_subjects"] == 3 and r["judgeable"] and r["by_unit"] == {"u1": 2.0, "u2": 6.0, "u3": None}
    assert r["value"] == pytest.approx(4.0) and r["n_by_unit"] == {"u1": 3, "u2": 2, "u3": 0}
    assert not AN.pair_summary({"s1": {"u1": 0.0}}, B, ["u1"])["judgeable"]


def test_common_subjects_counted_from_session_list_two_pairs_have_0():
    c = AN.common_counts(R.onh_sessions())
    assert len(c) == 21
    zero = sorted("~".join(sorted(k)) for k, v in c.items() if v["n"] == 0)
    assert zero == ["NOT4GEP~OXF3TRI", "OXF3TRI~OXF4GEP"] and sum(v["judgeable"] for v in c.values()) == 19
    assert c[frozenset(("NOT2ING", "OXF1PRI"))]["protocols"] == {"NOT2ING-89": 4, "NOT2ING-98": 16}
    assert c[frozenset(("NOT2ING", "NOT4GEP"))]["protocols"] == {"NOT2ING-98": 10}


def test_retest_SD_uses_same_formula_as_compute_s_ref():
    vals = {"a": [{"u": 1.0}, {"u": 3.0}, {"u": 5.0}], "b": [{"u": 2.0}, {"u": 2.0}]}
    s, df = AN.s_ref(vals, ["u", "v"])
    assert df == 3 and s["u"] == pytest.approx(((4.0 + 0.0) / 2) ** 0.5) and s["v"] is None


def test_SNR_is_WM_voxel_b0_mean_over_between_b0_noise():
    rng = np.random.default_rng(3)
    shape = (6, 6, 4)
    b0 = 500.0 + rng.normal(0, 20.0, shape + (6,))
    mask = np.ones(shape, bool)
    fa = np.full(shape, 0.5)
    fa[0] = 0.1                                                               # Outside white matter
    r = AN.snr(b0, mask, fa)
    assert r["n_b0"] == 6 and r["n_vox"] == 5 * 6 * 4 and r["snr"] == pytest.approx(25.0, rel=0.1) and r["in_train"]
    assert AN.snr(b0[..., :1], mask, fa)["snr"] is None
    assert not AN.snr(5000.0 + rng.normal(0, 20.0, shape + (6,)), mask, fa)["in_train"]   # 250 > 60


def test_s_ref_F2_adds_10_OXF1PRI_repeat_sessions_NOT2ING_TE89_has_four():
    ex = R.sref_extra()
    assert len(ex) == 10 and all(b in {f"ses-OXF1PRI00{i}" for i in range(2, 7)} for _, b in ex)
    assert {a for a, _ in ex} == {"sub-13192", "sub-16974"} and not set(ex) & set(R.onh_sessions())
    assert R.not2ing_protocol("sub-13192") == "NOT2ING-89" and R.not2ing_protocol("sub-16745") == "NOT2ING-98" and len(R.NOT2ING_TE89) == 4
    assert R.SNR_TRAIN == (15.0, 60.0) and R.PAIR_MIN_SUBJECTS == 3


def test_guard_also_binds_analysis_and_prediction_code():
    assert {"decomp/remeasure/analysis.py", "decomp/fixed/fiso_sign.py", "decomp/fixed/delta_tau.py",
            "decomp/stage1/lib.py", "decomp/stage1/synth.py"} <= set(G.PIPELINE)


# ── A guard stop is not a processing failure (review finding) · guard list · SD multiples ─────────────
def test_after_session_action_separates_guard_stop_from_failure():
    assert run.after_exit(0, 1) == "끝" and run.after_exit(1, 1) == "다시" and run.after_exit(1, 2) == "실패"
    assert run.after_exit(G.GUARD_EXIT, 1) == "가드 멈춤" and run.after_exit(G.GUARD_EXIT, 2) == "가드 멈춤"
    assert G.GUARD_EXIT not in (0, 1, 2)


def test_step_ending_with_guard_code_reraises_GuardError(tmp_path, monkeypatch):
    monkeypatch.setattr(run, "free_gb", lambda p: 1e6)                         # Does not rely on the actual free space of the test folder
    with pytest.raises(G.GuardError):
        run.step(tmp_path, "t", ["sh", "-c", f"exit {G.GUARD_EXIT}"])
    with pytest.raises(RuntimeError) as e:
        run.step(tmp_path, "t2", ["sh", "-c", "exit 3"])
    assert not isinstance(e.value, G.GuardError)


def test_guard_also_binds_stage1_dependencies_setting_prediction_values():
    assert {"decomp/stage1/config.py", "decomp/stage1/analyze.py", "decomp/stage1/analyze2.py",
            "decomp/stage1/analyze3.py"} <= set(G.PIPELINE) and len(G.PIPELINE) == 32


def test_guard_also_binds_verdict_code_frozen_record_and_env_lock():
    """Relay ZZZZh-3 · recommendation ZZZZh-7 4 — three verdict codes · the S4 code that V4 calls · frozen record · environment lock (weights outside the repository are checked by judge_v4.check_frozen)."""
    assert {"decomp/remeasure/judge_tong.py", "decomp/remeasure/judge_wong.py", "decomp/remeasure/judge_v4.py", "decomp/remeasure/judge_onh.py",
            "decomp/s4/validate.py", "decomp/s4/s4cfg.py", "decomp/s4/plan.py", "decomp/results/s4/train_freeze.json",
            "decomp/results/s4/sim_hashes.json", "decomp/results/s4/miss.json", "decomp/results/s4/cards_c3.json",
            "decomp/results/s4/cards_c7.json", "uv.lock", "pyproject.toml", "envs/reference/uv.lock",
            "envs/reference/pyproject.toml"} <= set(G.PIPELINE)
    assert len(set(G.PIPELINE)) == len(G.PIPELINE)


def test_SD_multipliers_match_compute_s_ref_and_are_redone_by_dof():
    lo, hi = AN.sd_multipliers(10)
    assert lo == pytest.approx(0.6987, abs=1e-4) and hi == pytest.approx(1.7549, abs=1e-4)
    assert AN.sd_multipliers(9)[1] > hi


def test_common_subject_protocol_order_is_sorted():
    c = AN.common_counts(R.onh_sessions())
    assert list(c[frozenset(("NOT2ING", "OXF3TRI"))]["protocols"]) == ["NOT2ING-89", "NOT2ING-98"]


# ── Re-check C-1 · C-3 — stop path of the batch · resumed counting of a relaunched batch · card pairs ─────────
class _FakeProc:
    """Fake child session — finishes at the polls-th poll (None before that) · wait finishes immediately."""
    def __init__(self, rc, polls=1):
        self.rc, self.left = rc, polls

    def poll(self):
        self.left -= 1
        return self.rc if self.left <= 0 else None

    def wait(self):
        return self.rc


def _fake_batch(monkeypatch, work, sessions, script, guard_ok=lambda n: True):
    """Fake child session — script[(sub, ses)] = exit code per attempt or (exit code, number of polls) (if 0, DONE.json is written) · guard_ok(n) = whether the n-th guard check passes · free disk is ample · waiting uses a fake clock (30 s each)."""
    launched, calls, clock = [], {"n": 0}, {"t": 0.0}

    def require():
        calls["n"] += 1
        if not guard_ok(calls["n"]):
            raise G.GuardError("가짜 가드 닫힘")
        return {}

    def launch(sub, ses):
        launched.append((sub, ses))
        v = script[(sub, ses)].pop(0)
        rc, polls, fin = (v + (v[0] == 0,))[:3] if isinstance(v, tuple) else (v, 1, v == 0)
        if fin:                                                              # Child session that wrote the end marker (separate from the exit code — RUN5-1)
            (work / "onharmony" / f"{sub}_{ses}").mkdir(parents=True, exist_ok=True)
            (work / "onharmony" / f"{sub}_{ses}" / "DONE.json").write_text("{}")
        return _FakeProc(rc, polls)

    monkeypatch.setattr(run.R, "onh_sessions", lambda: list(sessions))
    monkeypatch.setattr(run.R, "sref_extra", lambda: [])
    monkeypatch.setattr(run.G, "require", require)
    monkeypatch.setattr(run, "_launch", launch)
    monkeypatch.setattr(run, "free_gb", lambda p: 1e6)
    monkeypatch.setattr(run, "_now", lambda: clock["t"])
    monkeypatch.setattr(run, "_sleep", lambda s: clock.update(t=clock["t"] + s))
    return launched


SA, SB, SC = ("sub-a", "ses-X001"), ("sub-b", "ses-X001"), ("sub-c", "ses-X001")


def test_batch_guard_closed_before_launch_starts_none_lists_remaining(tmp_path, monkeypatch):
    launched = _fake_batch(monkeypatch, tmp_path, [SA, SB, SC], {SA: [0], SB: [0], SC: [0]}, guard_ok=lambda n: n <= 2)
    s = run.batch_onh(1, work=tmp_path)
    assert launched == [SA] and s["실패_다시도_실패"] == [] and "가짜 가드 닫힘" in s["가드로_멈춤"]["까닭"]
    assert s["가드로_멈춤"]["세션"] == [] and s["가드로_멈춤"]["남은_세션"] == ["sub-b_ses-X001", "sub-c_ses-X001"]


def test_batch_stops_without_counting_failure_if_child_session_exits_77(tmp_path, monkeypatch):
    launched = _fake_batch(monkeypatch, tmp_path, [SA, SB], {SA: [G.GUARD_EXIT], SB: [0]})
    s = run.batch_onh(1, work=tmp_path)
    assert launched == [SA] and s["실패_다시도_실패"] == [] and s["가드로_멈춤"]["세션"] == ["sub-a_ses-X001"]
    assert s["가드로_멈춤"]["남은_세션"] == ["sub-b_ses-X001"] and run.prior_failures(tmp_path) == {}
    assert not (tmp_path / "onharmony" / "sub-a_ses-X001" / "FAILED.json").exists()


def test_relaunch_continues_fail_count_skips_retry_failed_keeps_summary(tmp_path, monkeypatch):
    launched = _fake_batch(monkeypatch, tmp_path, [SA, SB, SC], {SA: [1, 1], SB: [1, G.GUARD_EXIT], SC: [0]})
    s1 = run.batch_onh(1, work=tmp_path)
    assert launched == [SA, SB, SC, SA, SB] and s1["실패_다시도_실패"] == ["sub-a_ses-X001"]
    assert s1["가드로_멈춤"]["세션"] == ["sub-b_ses-X001"] and run.prior_failures(tmp_path) == {SA: [1, 1], SB: [1]}
    launched2 = _fake_batch(monkeypatch, tmp_path, [SA, SB, SC], {SB: [1]})
    s2 = run.batch_onh(1, work=tmp_path)
    assert launched2 == [SB]                                                    # Do not launch SA (retry failure) or SC (end), and launch SB only once
    assert s2["실패_다시도_실패"] == ["sub-a_ses-X001", "sub-b_ses-X001"] and s2["가드로_멈춤"] is None
    assert [r["이번에_다시도_실패"] for r in s2["실행"]] == [["sub-a_ses-X001"], ["sub-b_ses-X001"]]
    assert s2["실행"][0]["가드로_멈춤"]["세션"] == ["sub-b_ses-X001"]
    log = [json.loads(x) for x in (tmp_path / "onharmony" / "batch_log.jsonl").read_text().splitlines()]
    assert [r["시도"] for r in log if r["sub"] == "sub-b" and r["무엇"] == "시작"] == [1, 2, 2]     # Guard stops are not counted, and counting resumes
    with pytest.raises(RuntimeError, match="다시도 실패"):
        run.onh_session("sub-a", "ses-X001", work=tmp_path, inputs={"ap": tmp_path / "x", "pa": tmp_path / "y"})


def test_retry_failure_of_batch_cut_before_marker_gets_marker_no_relaunch(tmp_path, monkeypatch):
    (tmp_path / "onharmony").mkdir()
    (tmp_path / "onharmony" / "attempts.jsonl").write_text(
        "".join(json.dumps(dict(utc="t", sub="sub-a", ses="ses-X001", rc=rc)) + "\n" for rc in (1, G.GUARD_EXIT, -9)))
    launched = _fake_batch(monkeypatch, tmp_path, [SA], {SA: [0]})
    s = run.batch_onh(1, work=tmp_path)
    assert launched == [] and s["실패_다시도_실패"] == ["sub-a_ses-X001"]
    assert json.loads((tmp_path / "onharmony" / "sub-a_ses-X001" / "FAILED.json").read_text())["종료_코드"] == [1, -9]


def test_stopped_batch_ends_with_stop_code(monkeypatch):
    base = dict(먼저_멈춘_까닭=None, 가드로_멈춤=None, 디스크로_멈춤=None, 잠김=[])
    for extra, code in ((dict(먼저_멈춘_까닭="가드", 가드로_멈춤={"까닭": "x"}), G.GUARD_EXIT),
                        (dict(먼저_멈춘_까닭="디스크", 디스크로_멈춤={"까닭": "x"}), run.STOP_EXIT),
                        (dict(잠김=["sub-a_ses-X001"]), run.BUSY_EXIT)):
        monkeypatch.setattr(run, "batch_onh", lambda n, s=dict(base, **extra): s)
        with pytest.raises(SystemExit) as e:
            run.main(["onharmony-batch"])
        assert e.value.code == code
    monkeypatch.setattr(run, "batch_onh", lambda n: dict(base))
    run.main(["onharmony-batch"])                                               # 0 if there is no stop or lock
    assert len({G.GUARD_EXIT, run.STOP_EXIT, run.BUSY_EXIT, 0, 1}) == 5 and set(run.NOT_COUNTED) == {0, G.GUARD_EXIT, run.STOP_EXIT, run.BUSY_EXIT}


def test_card_pairs_with_0_common_subjects_are_five():
    c = AN.card_common_counts(R.onh_sessions())
    assert len(c) == 28
    zero = sorted("~".join(sorted(k)) for k, v in c.items() if v == 0)
    assert zero == ["NOT2ING-89~NOT2ING-98", "NOT2ING-89~NOT4GEP", "NOT2ING-89~OXF4GEP", "NOT4GEP~OXF3TRI", "OXF3TRI~OXF4GEP"]


# ── Reflection checks RUN-1 – RUN-4 — disk stop · cap record · one-session command · lock · collision ─────────
def test_short_disk_before_launch_stops_without_launching_or_counting(tmp_path, monkeypatch):
    launched = _fake_batch(monkeypatch, tmp_path, [SA, SB], {SA: [0], SB: [0]})
    monkeypatch.setattr(run, "free_gb", lambda p: 0.0)
    s = run.batch_onh(1, work=tmp_path)
    assert launched == [] and s["먼저_멈춘_까닭"] == "디스크" and s["디스크로_멈춤"]["남은_세션"] == ["sub-a_ses-X001", "sub-b_ses-X001"]
    assert s["실패_다시도_실패"] == [] and run.prior_failures(tmp_path) == {} and s["가드로_멈춤"] is None


def test_child_session_exit_75_is_not_counted_and_stops_batch(tmp_path, monkeypatch):
    launched = _fake_batch(monkeypatch, tmp_path, [SA, SB], {SA: [run.STOP_EXIT], SB: [0]})
    s = run.batch_onh(1, work=tmp_path)
    assert launched == [SA] and s["디스크로_멈춤"]["세션"] == ["sub-a_ses-X001"] and s["디스크로_멈춤"]["남은_세션"] == ["sub-b_ses-X001"]
    assert run.prior_failures(tmp_path) == {} and not (tmp_path / "onharmony" / "sub-a_ses-X001" / "FAILED.json").exists()


def test_step_with_short_disk_is_StopRun_not_processing_failure(tmp_path, monkeypatch):
    monkeypatch.setattr(run, "free_gb", lambda p: 0.0)
    with pytest.raises(run.StopRun):
        run.step(tmp_path, "t", ["sh", "-c", "exit 0"])
    free = iter([100.0, 0.0, 0.0])                                              # The previous check was ample · the disk filled up in the middle of the stage
    monkeypatch.setattr(run, "free_gb", lambda p: next(free))
    with pytest.raises(run.StopRun):
        run.step(tmp_path, "t2", ["sh", "-c", "exit 3"])
    monkeypatch.setattr(run, "free_gb", lambda p: 100.0)
    with pytest.raises(RuntimeError) as e:
        run.step(tmp_path, "t3", ["sh", "-c", "exit 3"])
    assert not isinstance(e.value, run.StopRun)


def test_child_session_exits_75_when_disk_fills_outside_a_step(monkeypatch):
    def boom(argv):
        raise OSError(28, "No space left on device")
    monkeypatch.setattr(run, "_main", boom)
    monkeypatch.setattr(run, "low_disk", lambda p: True)
    for worker in (True, False):                                                 # Child sessions too · other commands such as Tong · Wong too (RUN5-3)
        if worker:
            monkeypatch.setenv(run.WORKER, "1")
        else:
            monkeypatch.delenv(run.WORKER, raising=False)
        for cmd in (["onharmony", "sub-a", "ses-X001"], ["tong", "sub-1", "ses-c01"], ["wong", "sub-01"]):
            with pytest.raises(SystemExit) as e:
                run.main(cmd)
            assert e.value.code == run.STOP_EXIT
    monkeypatch.setattr(run, "low_disk", lambda p: False)
    with pytest.raises(OSError):                                                 # If the disk is ample, pass it up as is (processing failure)
        run.main(["onharmony", "sub-a", "ses-X001"])


def test_locked_session_not_counted_or_relaunched_in_this_batch(tmp_path, monkeypatch):
    launched = _fake_batch(monkeypatch, tmp_path, [SA, SB], {SA: [run.BUSY_EXIT], SB: [0]})
    s = run.batch_onh(1, work=tmp_path)
    assert launched == [SA, SB] and s["잠김"] == ["sub-a_ses-X001"] and s["먼저_멈춘_까닭"] is None
    assert run.prior_failures(tmp_path) == {} and s["실패_다시도_실패"] == []


def test_first_failure_returned_after_cap_is_added_to_left_by_cap(tmp_path, monkeypatch):
    launched = _fake_batch(monkeypatch, tmp_path, [SA, SB, SC], {SA: [(1, 3)], SB: [(0, 2)], SC: [0]})
    s = run.batch_onh(2, work=tmp_path, cap_s=45)
    assert launched == [SA, SB] and s["실행"][-1]["상한으로_남김"] == ["sub-c_ses-X001", "sub-a_ses-X001"]
    assert run.prior_failures(tmp_path) == {SA: [1]}


def test_single_session_cmd_runs_child_logs_try_two_failures_write_marker(tmp_path, monkeypatch):
    launched = _fake_batch(monkeypatch, tmp_path, [SA], {SA: [1, 1]})
    assert run.single_onh("sub-a", "ses-X001", work=tmp_path) == 1 and run.prior_failures(tmp_path) == {SA: [1]}
    assert not (tmp_path / "onharmony" / "sub-a_ses-X001" / "FAILED.json").exists()
    assert run.single_onh("sub-a", "ses-X001", work=tmp_path) == 1 and run.prior_failures(tmp_path) == {SA: [1, 1]}
    assert (tmp_path / "onharmony" / "sub-a_ses-X001" / "FAILED.json").exists() and launched == [SA, SA]
    with pytest.raises(RuntimeError, match="다시도 실패"):
        run.single_onh("sub-a", "ses-X001", work=tmp_path)
    assert run.batch_onh(1, work=tmp_path)["실패_다시도_실패"] == ["sub-a_ses-X001"] and launched == [SA, SA]
    attempts = [json.loads(x) for x in (tmp_path / "onharmony" / "attempts.jsonl").read_text().splitlines()]
    assert [a["띄운_쪽"] for a in attempts] == ["한 세션 명령", "한 세션 명령"]


def test_held_work_folder_lock_blocks_batch_and_single_session_command(tmp_path, monkeypatch):
    launched = _fake_batch(monkeypatch, tmp_path, [SA], {SA: [0, 0]})
    with run.run_lock(tmp_path):
        with pytest.raises(run.Busy):
            run.batch_onh(1, work=tmp_path)
        with pytest.raises(run.Busy):
            run.single_onh("sub-a", "ses-X001", work=tmp_path)
    assert launched == [] and run.batch_onh(1, work=tmp_path)["실패_다시도_실패"] == [] and launched == [SA]


def test_held_session_lock_makes_child_locked_passes_to_step_subprocess(tmp_path, monkeypatch):
    lock = tmp_path / "onharmony" / "locks" / "sub-a_ses-X001.lock"
    with run.flock_or_busy(lock, "다른 하위 세션"):
        with pytest.raises(run.Busy):
            run.onh_session("sub-a", "ses-X001", work=tmp_path, inputs={"ap": tmp_path / "x", "pa": tmp_path / "y"})
    monkeypatch.setattr(run, "free_gb", lambda p: 1e6)
    with run.flock_or_busy(tmp_path / "l.lock", "시험", hand_down=True):
        fd = run._HELD[-1]
        run.step(tmp_path, "fd", ["sh", "-c", "ls /proc/$$/fd"])
    assert str(fd) in (tmp_path / "log_fd.txt").read_text().split() and run._HELD == []


def test_end_marker_with_two_failures_writes_no_marker_reports_conflict(tmp_path, monkeypatch):
    (tmp_path / "onharmony" / "sub-a_ses-X001").mkdir(parents=True)
    (tmp_path / "onharmony" / "sub-a_ses-X001" / "DONE.json").write_text("{}")
    (tmp_path / "onharmony" / "attempts.jsonl").write_text(
        "".join(json.dumps(dict(utc="t", sub="sub-a", ses="ses-X001", rc=1)) + "\n" for _ in range(2)))
    launched = _fake_batch(monkeypatch, tmp_path, [SA], {SA: [0]})
    s = run.batch_onh(1, work=tmp_path)
    assert launched == [] and s["충돌"] == ["sub-a_ses-X001"] and s["실패_다시도_실패"] == []
    assert not (tmp_path / "onharmony" / "sub-a_ses-X001" / "FAILED.json").exists()


# ── Reflection check 2 — parent-side verdict · record file · Tong disk · cap after a stop · test reinforcement (RUN5-1 – RUN5-6) ────
def test_end_marker_means_done_any_exit_code_failure_plus_low_disk_is_stop(tmp_path, monkeypatch):
    launched = _fake_batch(monkeypatch, tmp_path, [SA, SB], {SA: [(120, 1, True)], SB: [1]})
    state = {"free": 1e6}
    monkeypatch.setattr(run, "free_gb", lambda p: state["free"])
    inner = run._launch

    def launch(sub, ses):
        pr = inner(sub, ses)
        if (sub, ses) == SB:
            state["free"] = 0.0                                                  # The disk fills while SB runs
        return pr
    monkeypatch.setattr(run, "_launch", launch)
    s = run.batch_onh(1, work=tmp_path)
    attempts = [json.loads(x) for x in (tmp_path / "onharmony" / "attempts.jsonl").read_text().splitlines()]
    assert [(a["sub"], a["rc"], a.get("실제_종료_코드")) for a in attempts] == [("sub-a", 0, 120), ("sub-b", run.STOP_EXIT, 1)]
    assert run.prior_failures(tmp_path) == {} and s["먼저_멈춘_까닭"] == "디스크" and launched == [SA, SB]
    assert not any((tmp_path / "onharmony" / f"{k[0]}_{k[1]}" / "FAILED.json").exists() for k in (SA, SB))


def test_half_lines_and_broken_summary_do_not_block_next_batch(tmp_path, monkeypatch):
    (tmp_path / "onharmony").mkdir()
    (tmp_path / "onharmony" / "attempts.jsonl").write_text(json.dumps(dict(utc="t", sub="sub-a", ses="ses-X001", rc=1)) + '\n{"utc": "t", "sub": "sub-a", "se')
    (tmp_path / "onharmony" / "batch_summary.json").write_text('{"실행": [')
    launched = _fake_batch(monkeypatch, tmp_path, [SA], {SA: [1]})
    s = run.batch_onh(1, work=tmp_path)
    assert launched == [SA] and s["실패_다시도_실패"] == ["sub-a_ses-X001"] and len(s["실행"]) == 1   # Previous failure 1 + this one 1 = retry failure
    lines = (tmp_path / "onharmony" / "attempts.jsonl").read_text().splitlines()
    assert json.loads(lines[-1])["rc"] == 1 and run.prior_failures(tmp_path) == {SA: [1, 1]}
    assert len(list((tmp_path / "onharmony").glob("batch_summary.broken_*.json"))) == 1 and not (tmp_path / "onharmony" / "batch_summary.tmp").exists()


def test_tong_extract_stops_without_merging_zip_when_disk_is_short(tmp_path, monkeypatch):
    monkeypatch.setattr(run, "guard_if_real", lambda paths: None)
    monkeypatch.setattr(run, "low_disk", lambda p: True)
    with pytest.raises(run.StopRun):
        run.tong_extract("sub-1", work=tmp_path)
    assert not (tmp_path / "tong" / "sub-1_merged.zip").exists()


def test_after_stop_remaining_sessions_not_moved_to_cap_even_past_cap(tmp_path, monkeypatch):
    launched = _fake_batch(monkeypatch, tmp_path, [SA, SB, SC], {SA: [(0, 5)], SB: [0], SC: [0]}, guard_ok=lambda n: n <= 2)
    s = run.batch_onh(2, work=tmp_path, cap_s=45)
    assert launched == [SA] and s["가드로_멈춤"]["남은_세션"] == ["sub-b_ses-X001", "sub-c_ses-X001"] and s["실행"][-1]["상한으로_남김"] == []


def test_after_stop_waits_for_running_sessions_and_records_their_failures(tmp_path, monkeypatch):
    launched = _fake_batch(monkeypatch, tmp_path, [SA, SB], {SA: [G.GUARD_EXIT], SB: [(1, 3)]})
    s = run.batch_onh(2, work=tmp_path)
    attempts = [json.loads(x) for x in (tmp_path / "onharmony" / "attempts.jsonl").read_text().splitlines()]
    assert launched == [SA, SB] and [(a["sub"], a["rc"]) for a in attempts] == [("sub-a", G.GUARD_EXIT), ("sub-b", 1)]
    assert run.prior_failures(tmp_path) == {SB: [1]} and s["가드로_멈춤"]["세션"] == ["sub-a_ses-X001"] and s["가드로_멈춤"]["남은_세션"] == ["sub-b_ses-X001"]


def test_first_stop_reason_is_guard_even_if_session_exits_75_after_guard(tmp_path, monkeypatch):
    _fake_batch(monkeypatch, tmp_path, [SA, SB], {SA: [run.STOP_EXIT], SB: [0]}, guard_ok=lambda n: n <= 2)
    s = run.batch_onh(2, work=tmp_path)
    assert s["먼저_멈춘_까닭"] == "가드" and s["디스크로_멈춤"]["세션"] == ["sub-a_ses-X001"] and s["가드로_멈춤"]["세션"] == []


def test_single_session_cmd_refuses_two_failures_no_marker_writes_marker(tmp_path, monkeypatch):
    (tmp_path / "onharmony").mkdir()
    (tmp_path / "onharmony" / "attempts.jsonl").write_text(
        "".join(json.dumps(dict(utc="t", sub="sub-a", ses="ses-X001", rc=1)) + "\n" for _ in range(2)))
    launched = _fake_batch(monkeypatch, tmp_path, [SA], {SA: [0]})
    with pytest.raises(RuntimeError, match="다시도 실패"):
        run.single_onh("sub-a", "ses-X001", work=tmp_path)
    assert launched == [] and (tmp_path / "onharmony" / "sub-a_ses-X001" / "FAILED.json").exists()


def test_child_session_passes_its_session_lock_to_steps(tmp_path, monkeypatch):
    seen = {}

    def body(sub, ses, work, inputs):
        seen["held"] = [os.readlink(f"/proc/self/fd/{fd}") for fd in run._HELD]
        return {}
    monkeypatch.setattr(run, "_onh_session", body)
    run.onh_session("sub-a", "ses-X001", work=tmp_path)
    assert seen["held"] == [str(tmp_path / "onharmony" / "locks" / "sub-a_ses-X001.lock")] and run._HELD == []
