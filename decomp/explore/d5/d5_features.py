"""Extracts exploratory d5 features: E1 – E6 from each session's eddy text logs, and 27 ROI units × (x, y, ODI) from roi.json."""
from __future__ import annotations

import datetime
import hashlib
import json
import math
import os
import pathlib
import shutil
import subprocess
import sys

import numpy as np

REPO = pathlib.Path("<REPO>")
X14 = pathlib.Path("<WORK>/stage1_work/drafts/x14")
for p in ("decomp/remeasure", "decomp/stage1", "decomp/s4"):
    sys.path.insert(0, str(REPO / p))
import rconfig as R  # noqa: E402
import lib  # noqa: E402

OUT = R.WORK / "explore" / "d5"
ONH = R.WORK / "onharmony"
NEWLOG, REREAD = X14 / "ONH_NEW_READ_LOG.jsonl", X14 / "REREAD_LOG.jsonl"
NOW = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
EDDY_FILES = ("eddy_out.eddy_movement_rms", "eddy_out.eddy_outlier_map", "eddy_out.eddy_parameters",
              "eddy_out.eddy_post_eddy_shell_alignment_parameters", "bvals")
MAP_IN = tuple(f"f2_{f}" for f in R.FEATURES) + ("nodif_brain_mask", "f2_extvalid", "f2_f2valid", "fa2std_warp")
VR_TOL = 1e-12


def sha(p) -> str:
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()


def rel(p) -> str:
    return os.path.relpath(p, REPO)


def entry(path, what, key, basis=None) -> dict:
    b = pathlib.Path(path).read_bytes()
    return dict(시각=NOW, utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"), 세션=key, 파일=rel(path),
                sha256=hashlib.sha256(b).hexdigest(), 바이트=len(b), 무엇=what, 근거=basis, 스크립트=__file__)


def write_log(path, rows):
    with pathlib.Path(path).open("a", encoding="utf-8") as fh:
        for r in rows:
            fh.write(json.dumps(r, ensure_ascii=False) + "\n")


def sessions() -> list[tuple[str, str]]:
    """Sessions whose processing has finished (among 110 first sessions + 10 OXF1PRI repeats, those with DONE.json and without FAILED.json)."""
    out = []
    for sub, ses in R.onh_sessions() + R.sref_extra():
        wd = ONH / f"{sub}_{ses}"
        if (wd / "DONE.json").exists() and not (wd / "FAILED.json").exists():
            out.append((sub, ses))
    return out


def unit_of(sub: str, ses: str) -> str:
    d = R.device_of(ses)
    return R.not2ing_protocol(sub) if d == "NOT2ING" else d


# ── eddy ───────────────────────────────────────────────────────────────────────
def shell_align_max(path) -> float:
    vals = []
    for ln in pathlib.Path(path).read_text(encoding="utf-8").splitlines():
        parts = ln.split()
        try:
            nums = [float(x) for x in parts]
        except ValueError:
            continue
        if len(nums) == 6:
            vals.extend(abs(v) for v in nums[:3])
    return float(max(vals)) if vals else float("nan")


def eddy_feats(wd: pathlib.Path) -> dict:
    mr = np.loadtxt(wd / "eddy_out.eddy_movement_rms", ndmin=2)
    bv = np.loadtxt(wd / "bvals").ravel()
    om = np.loadtxt(wd / "eddy_out.eddy_outlier_map", skiprows=1, ndmin=2)
    par = np.loadtxt(wd / "eddy_out.eddy_parameters", ndmin=2)
    assert len(mr) == len(bv) == len(om) == len(par), f"부피 수가 다르다 {wd.name}"
    dwi = bv > 50
    k, N = int(om[dwi].sum()), int(dwi.sum() * om.shape[1])
    return dict(E1=float(mr[:, 0].mean()), E2=float(mr[1:, 1].mean()), E3=k / N, E3_k=k, E3_N=N, E4=float(mr[1:, 1].max()),
                E5=float(np.sqrt(np.mean(par[dwi, 6:9] ** 2))), E6=shell_align_max(wd / "eddy_out.eddy_post_eddy_shell_alignment_parameters"),
                n_vol=int(len(bv)), n_dwi=int(dwi.sum()), n_slice=int(om.shape[1]))


def mode_eddy():
    OUT.mkdir(parents=True, exist_ok=True)
    import guard as G
    basis = G.require()
    rows, logs = {}, []
    for sub, ses in sessions():
        wd, key = ONH / f"{sub}_{ses}", f"{sub}_{ses}"
        logs += [entry(wd / f, "D5 eddy 글 기록(영상 아님 · §9 s4)", key, basis) for f in EDDY_FILES]
        rows[key] = dict(sub=sub, ses=ses, device=R.device_of(ses), unit=unit_of(sub, ses), **eddy_feats(wd))
    write_log(NEWLOG, logs)
    (OUT / "features_eddy.json").write_text(json.dumps(dict(made=NOW, script_sha256=sha(__file__), n=len(rows), sessions=rows), ensure_ascii=False, indent=1) + "\n",
                                            encoding="utf-8")
    print("eddy", len(rows), "세션 · 읽기 기록", len(logs))


# ── roi ────────────────────────────────────────────────────────────────────────
def roi_coords(res: dict, units: list[str]) -> np.ndarray:
    out = np.full((len(units), 3), np.nan)
    for i, u in enumerate(units):
        row = res["주"].get(u, {})
        v = [row.get(f) for f in ("NDI", "FWF", "ODI")]
        if any(a is None or not math.isfinite(a) for a in v):
            continue
        out[i] = lib.coords(np.array([v[0]]), np.array([v[1]]), np.array([v[2]]))[0][0]
    return out


def mode_roi():
    OUT.mkdir(parents=True, exist_ok=True)
    import roi as ROI
    units = ROI.unit_order()
    keys, X, logs = [], [], []
    for sub, ses in sessions():
        p = ONH / f"{sub}_{ses}" / "roi.json"
        b = p.read_bytes()
        logs.append(dict(시각=NOW, 다시읽음=str(p), sha256=hashlib.sha256(b).hexdigest(), 바이트=len(b), 무엇="D5 — roi.json «주» 좌표(① T · ②)", 스크립트=__file__))
        keys.append(f"{sub}_{ses}")
        X.append(roi_coords(json.loads(b), units))
    write_log(REREAD, logs)
    np.savez(OUT / "features_roi.npz", keys=np.array(keys), units=np.array(units), X=np.array(X), made=NOW, script_sha256=sha(__file__))
    print("roi", len(keys), "세션 ·", len(units), "단위")


# ── maps ───────────────────────────────────────────────────────────────────────
def map_worker(args):
    sub, ses = args
    import nibabel as nib
    import cmds as K
    import roi as ROI
    import run as RUN
    key, wd = f"{sub}_{ses}", ONH / f"{sub}_{ses}"
    tmp = OUT / "tmp" / key
    tmp.mkdir(parents=True, exist_ok=True)
    try:
        jobs = [(f"f2_{f}", f"std_f2_{f}", "trilinear") for f in R.FEATURES] + [("nodif_brain_mask", "std_f2_mask", "nn"),
                                                                                 ("f2_extvalid", "std_f2_extvalid", "trilinear"),
                                                                                 ("f2_f2valid", "std_f2_f2valid", "trilinear")]
        for src, out, interp in jobs:          # Same arguments as cmds.applywarp — only the input and output are absolute paths (nothing is written to the session folder)
            argv = K.applywarp(str(wd / src), str(tmp / out), str(wd / "fa2std_warp"), interp=interp)
            r = subprocess.run(["nice", "-n", str(R.NICE), *argv], env=RUN.ENV, capture_output=True, text=True)
            if r.returncode != 0:
                raise RuntimeError(f"applywarp 실패 {key} {src}: {r.stderr[-300:]}")
        std = {f: tmp / f"std_f2_{f}.nii.gz" for f in R.FEATURES}
        aux = [tmp / f"std_f2_{v}.nii.gz" for v in ("mask", "extvalid", "f2valid")]
        res = ROI.extract(std, *aux)
        stored = json.loads((wd / "roi.json").read_text(encoding="utf-8"))["주"]
        maxd, npix_same, none_same = 0.0, True, True
        for u, row in res["주"].items():
            s = stored[u]
            npix_same &= row["화소"] == s["화소"]
            for f in R.FEATURES:
                a, b = row[f], s[f]
                if (a is None) != (b is None):
                    none_same = False
                elif a is not None:
                    maxd = max(maxd, abs(a - b))
        # Voxel features
        load = lambda p: np.asarray(nib.load(str(p)).dataobj, dtype=float)  # noqa: E731
        atlas = np.asarray(nib.load(str(R.ATLAS)).dataobj).astype(int)
        skel = ROI.skeleton_mask(np.asarray(nib.load(str(R.FSLDIR / (R.SKELETON + ".nii.gz"))).dataobj))
        maps = {f: load(std[f]) for f in ("NDI", "FWF", "ODI", "FA")}
        brain, extv, f2v = (load(p) for p in aux)
        valid = ROI.valid_mask(brain, f2v, extv)
        base = skel & (atlas > 0)
        names = ROI.label_units()
        units = ROI.unit_order()
        P = np.full((len(units), 6), np.nan)
        for i, u in enumerate(units):
            labs = [j for j, nm in enumerate(names, start=1) if nm == u]
            sel_all = base & np.isin(atlas, labs)
            sel = sel_all & valid
            n_all, n_v = int(sel_all.sum()), int(sel.sum())
            if n_v >= 2:
                c, _ = lib.coords(maps["NDI"][sel], maps["FWF"][sel], maps["ODI"][sel])
                P[i, 0:3] = c.std(0, ddof=1)
                P[i, 3] = float(np.percentile(maps["FWF"][sel], 90))
            P[i, 4] = n_v / n_all if n_all else np.nan
            P[i, 5] = math.log(n_v) if n_v > 0 else np.nan
        tfa = load(R.FSLDIR / (R.FA_REF + ".nii.gz"))
        m = skel & (brain > 0.5)
        Gr = float(np.corrcoef(maps["FA"][m], tfa[m])[0, 1])
        logs = [entry(wd / f"{s}.nii.gz", "D5 지도 · 마스크 · 변환(영상 · §9 s1 을 119 세션으로)", key) for s in MAP_IN]
        logs_rr = [dict(시각=NOW, 다시읽음=str(wd / "roi.json"), sha256=sha(wd / "roi.json"), 바이트=(wd / "roi.json").stat().st_size,
                        무엇="D5 V-R — 저장된 roi.json «주» 와 대조", 스크립트=__file__)]
        return key, dict(P=P, G=Gr, vr=dict(max_abs_diff=maxd, npix_same=bool(npix_same), none_same=bool(none_same),
                                            passed=bool(npix_same and none_same and maxd <= VR_TOL))), logs, logs_rr
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def mode_maps():
    OUT.mkdir(parents=True, exist_ok=True)
    import guard as G
    import roi as ROI
    from multiprocessing import get_context
    basis = G.require()
    ss = sessions()
    with get_context("fork").Pool(8) as pool:
        out = pool.map(map_worker, ss)
    keys = [k for k, *_ in out]
    logs = [dict(e, 근거=basis) for _, _, lg, _ in out for e in lg]
    write_log(NEWLOG, logs)
    write_log(REREAD, [e for *_, rr in out for e in rr])
    vr = {k: d["vr"] for k, d, *_ in out}
    gate = dict(rule=f"다시 낸 roi.json «주» 27 × 7 이 저장값과 |Δ| ≤ {VR_TOL} · 화소 수 같음 · 값 없음 자리 같음(계획서 §17 V-R)", made=NOW,
                passed=all(v["passed"] for v in vr.values()), n=len(vr), n_failed=sum(not v["passed"] for v in vr.values()),
                max_abs_diff=max(v["max_abs_diff"] for v in vr.values()), sessions=vr)
    (OUT / "vr_gate.json").write_text(json.dumps(gate, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    np.savez(OUT / "features_map.npz", keys=np.array(keys), units=np.array(ROI.unit_order()), P=np.array([d["P"] for _, d, *_ in out]),
             G=np.array([d["G"] for _, d, *_ in out]), made=NOW, script_sha256=sha(__file__),
             P_names=np.array(["SD_x", "SD_y", "SD_ODI", "FWF_p90", "valid_frac", "log_n_valid"]))
    print("maps", len(keys), "세션 · V-R", "통과" if gate["passed"] else f"실패 {gate['n_failed']}", "· 최대 차", gate["max_abs_diff"], "· 읽기 기록", len(logs))


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    dict(eddy=mode_eddy, roi=mode_roi, maps=mode_maps).get(cmd, lambda: sys.exit("쓰임: eddy|roi|maps"))()
