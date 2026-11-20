"""S73 §2 (post hoc): rigid-registration check of the multi-TE series with failure controls, the fixed alignment metric recomputed, and the analyses with one registration for all series (paths B and C); it stopped at its first failure control."""
from __future__ import annotations

import json
import math
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

import numpy as np

REPO = pathlib.Path(os.environ.get("S73_REPO") or pathlib.Path(__file__).resolve().parents[3])
sys.path[:0] = [str(REPO / "decomp" / "remeasure"), str(REPO / "decomp" / "stage1"), str(REPO / "decomp" / "v4_recompute")]
import cmds as K  # noqa: E402
import judge_wong as JW  # noqa: E402
import rconfig as R  # noqa: E402
import roi as ROI  # noqa: E402
import run as RUN  # noqa: E402 — mean_dw_ratio · align_r (the metric of the fixed check) unchanged
from read_log import ReadLog  # noqa: E402

ENV = dict(os.environ, CUDA_VISIBLE_DEVICES="", OMP_NUM_THREADS="1", FSLDIR=str(R.FSLDIR), FSLOUTPUTTYPE="NIFTI_GZ",
           PATH=f"{R.FSLDIR / 'bin'}:{os.environ.get('PATH', '')}")
VALID = ("extvalid", "f2valid")
WEDGE = (48.5, 72.1)
LIT_T2 = (82.5, 52.0, 1000.0)          # amendment item 5 — the literature vector as in Figure 5 (Gong 2020 · ledger derivation derived.json gong.mean_both: T2_in · T2_en · T2_iso, ms)
T2IN, T2ISO = 82.0, 1000.0
LOG: ReadLog


def nib_load(p):
    import nibabel as nib
    LOG.before_read([p])
    return nib.load(str(p))


def sh(argv, cwd) -> None:
    r = subprocess.run(argv, cwd=cwd, env=ENV, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"⛔ {' '.join(map(str, argv))}\n{r.stdout[-1500:]}{r.stderr[-1500:]}")


def fsl_coords(img, ijk: np.ndarray) -> np.ndarray:
    """Voxel indices → FSL scaled mm coordinates (the first axis is flipped when the affine determinant is positive — the FLIRT convention)."""
    pix = np.asarray(img.header.get_zooms()[:3], float)
    ijk = ijk.astype(float).copy()
    if np.linalg.det(img.affine[:3, :3]) > 0:
        ijk[:, 0] = img.shape[0] - 1 - ijk[:, 0]
    return ijk * pix


def max_disp(mat: np.ndarray, pts: np.ndarray) -> tuple[float, float]:
    """Largest |M·p − p| and largest |M⁻¹·p − p| over the mask voxel centers p (mm)."""
    P = np.c_[pts, np.ones(len(pts))]
    f = np.linalg.norm((P @ mat.T)[:, :3] - pts, axis=1).max()
    b = np.linalg.norm((P @ np.linalg.inv(mat).T)[:, :3] - pts, axis=1).max()
    return float(f), float(b)


def flirt(src: pathlib.Path, ref: pathlib.Path, omat: pathlib.Path, cwd) -> np.ndarray:
    LOG.before_read([src, ref])
    sh([K.fsl("flirt"), "-in", str(src), "-ref", str(ref), "-dof", "6", "-cost", "normmi", "-searchcost", "normmi", "-nosearch", "-omat", str(omat)], cwd)   # amendment (S73_AMENDMENT): -nosearch
    return np.loadtxt(omat)


def check_series(sub: str, tmp: pathlib.Path) -> dict:
    wd = R.WORK / "wong" / sub
    ref_b0, mask_p = wd / "b0_mean.nii.gz", wd / "nodif_brain_mask.nii.gz"
    mimg = nib_load(mask_p)
    pts = fsl_coords(mimg, np.argwhere(np.asarray(mimg.dataobj) > 0))
    out = {}
    for s in R.WONG_SERIES:
        if s == R.WONG_REF:
            continue
        b0 = wd / f"own_{s}" / "b0_mean.nii.gz"
        bimg = nib_load(b0)
        vox = float(min(bimg.header.get_zooms()[:3]))
        mat = flirt(b0, ref_b0, tmp / f"{sub}_{s}.mat", tmp)
        f, b = max_disp(mat, pts)
        out[s] = dict(max_disp_mm=f, max_disp_inverse_mm=b, voxel_mm=vox, max_disp_vox=f / vox, aligned=bool(f <= 0.5 * vox),
                      mat=mat.tolist(), zooms=[float(z) for z in bimg.header.get_zooms()[:3]], shape=list(bimg.shape[:3]))
    return out


def controls(sub: str, tmp: pathlib.Path) -> dict:
    """§2.2 failure controls — (i) reference to itself ≤ 0.1 voxel · (ii) a copy shifted by exactly 2 voxels along the first axis: 2 ± 0.25 voxels and classified «not aligned»."""
    import nibabel as nib
    wd = R.WORK / "wong" / sub
    ref_b0 = wd / "b0_mean.nii.gz"
    rimg, mimg = nib_load(ref_b0), nib_load(wd / "nodif_brain_mask.nii.gz")
    pts = fsl_coords(mimg, np.argwhere(np.asarray(mimg.dataobj) > 0))
    pix = np.asarray(rimg.header.get_zooms()[:3], float)
    vmin = float(pix.min())
    m1 = flirt(ref_b0, ref_b0, tmp / f"{sub}_self.mat", tmp)
    d1 = max_disp(m1, pts)[0]
    a = np.asarray(rimg.dataobj, dtype=np.float32)
    sh2 = np.zeros_like(a)
    sh2[2:] = a[:-2]
    p2 = tmp / f"{sub}_shift2.nii.gz"
    nib.save(nib.Nifti1Image(sh2, rimg.affine, rimg.header), str(p2))
    m2 = flirt(p2, ref_b0, tmp / f"{sub}_shift2.mat", tmp)
    d2 = max_disp(m2, pts)[0]
    d2_vox = d2 / pix[0]
    ok = d1 / vmin <= 0.1 and abs(d2_vox - 2.0) <= 0.25 and not (d2 <= 0.5 * vmin)
    return dict(self_max_disp_vox=d1 / vmin, shift2_max_disp_vox=d2_vox, shift2_classified_aligned=bool(d2 <= 0.5 * vmin),
                self_mat=m1.tolist(), shift2_mat=m2.tolist(), ok=bool(ok))


def cause_check(sub: str) -> dict:
    """§2.1 — the metric of the fixed check (voxel correlation of the mean b0-normalized signal, reference mask) recomputed from the source data, and the same without the voxels in the lowest 1 % of the reference mean b0."""
    d = R.WONG / sub / "dwi"
    wd = R.WORK / "wong" / sub
    msk = np.asarray(nib_load(wd / "nodif_brain_mask.nii.gz").dataobj) > 0
    ser = {s: d / f"{sub}_acq-{s}_dwi" for s in R.WONG_SERIES}

    def ratio_and_b0(p):
        LOG.before_read([pathlib.Path(f"{p}.bval")])
        bv = np.loadtxt(f"{p}.bval").ravel()
        data = np.asarray(nib_load(pathlib.Path(f"{p}.nii.gz")).dataobj, dtype=np.float32)
        return RUN.mean_dw_ratio(data, bv), data[..., bv <= R.B0_MAX_PREP].mean(-1)

    ref_map, ref_b0 = ratio_and_b0(ser[R.WONG_REF])
    thr = float(np.percentile(ref_b0[msk], 1.0))
    keep = msk & (ref_b0 > thr)
    out = dict(b0_p1_threshold=thr, n_mask=int(msk.sum()), n_kept=int(keep.sum()), series={})
    for s, p in ser.items():
        m, _ = ratio_and_b0(p) if s != R.WONG_REF else (ref_map, None)
        out["series"][s] = dict(r_fixed=RUN.align_r(ref_map, m, msk), r_without_low_b0=RUN.align_r(ref_map, m, keep))
    return out


def roi_path(sub: str, s: str, premat: pathlib.Path | None, tmp: pathlib.Path) -> dict:
    """Native-space F2 maps of the series → (premat) → fa2std_warp of the reference series → standard space → ROI (same commands and functions as to_std_and_roi). The brain mask is that of the reference series (nn, no premat)."""
    wd = R.WORK / "wong" / sub
    sd = wd if s == R.WONG_REF else wd / f"own_{s}"
    warp = wd / "fa2std_warp.nii.gz"
    out = tmp / f"{sub}_{s}_{'C' if premat else 'B'}"
    out.mkdir(exist_ok=True)
    LOG.before_read([warp] + [sd / f"f2_{s}_{f}.nii.gz" for f in R.FEATURES + VALID] + [wd / "nodif_brain_mask.nii.gz"])
    std = {}
    for f in R.FEATURES + VALID:
        argv = K.applywarp(str(sd / f"f2_{s}_{f}"), str(out / f"std_{f}"), str(warp))
        if premat is not None:
            argv.append(f"--premat={premat}")
        sh(argv, out)
        std[f] = out / f"std_{f}.nii.gz"
    sh(K.applywarp(str(wd / "nodif_brain_mask"), str(out / "std_mask"), str(warp), interp="nn"), out)
    res = ROI.extract({f: std[f] for f in R.FEATURES}, out / "std_mask.nii.gz", std["extvalid"], std["f2valid"])
    shutil.rmtree(out)
    return res


def psi_summary(data: dict) -> dict:
    """S-M8.3 E1 — ψ = atan2(b_y, b_x) · median · 90 % interval (2,000 resamples of units within participants · seed 20261031) · inside the wedge · b_y ≤ 0 · median r."""
    te = np.array(JW.TES, float)
    per = {}
    for sb in sorted(data):
        ser = data[sb]
        if not all(f"TE{t}" in ser for t in JW.TES):
            return dict(computable=False, reason=f"{sb}: TE 여덟이 다 있지 않음")
        units = sorted(set.intersection(*[set(ser[f"TE{t}"]) for t in JW.TES]))
        rows = []
        for u in units:
            bx = JW.fit_unit(te, np.array([ser[f"TE{t}"][u][0] for t in JW.TES], float))["slope"]
            by = JW.fit_unit(te, np.array([ser[f"TE{t}"][u][1] for t in JW.TES], float))["slope"]
            psi = math.degrees(math.atan2(by, bx))
            r = None if by <= 0 else (1.0 if abs(psi - 90.0) < 1e-12 else 1 + (1 - T2IN / T2ISO) / math.tan(math.radians(psi)))
            rows.append(dict(b_x=bx, b_y=by, psi=psi, r=r))
        per[sb] = rows

    def summ(rs):
        ps = [r["psi"] for r in rs]
        rr = [r["r"] for r in rs if r["r"] is not None]
        return dict(n=len(ps), psi_median=float(np.median(ps)), in_wedge=int(sum(WEDGE[0] <= p <= WEDGE[1] for p in ps)),
                    n_by_le_0=int(sum(r["b_y"] <= 0 for r in rs)), r_median=float(np.median(rr)) if rr else None,
                    b_x_median=float(np.median([r["b_x"] for r in rs])), b_y_median=float(np.median([r["b_y"] for r in rs])))
    lit_bx, lit_by = 1.0 / LIT_T2[1] - 1.0 / LIT_T2[0], 1.0 / LIT_T2[0] - 1.0 / LIT_T2[2]       # amendment item 5: T2 of Gong 2020 · per ms
    allr = [r for v in per.values() for r in v]
    lit = dict(T2_in=LIT_T2[0], T2_en=LIT_T2[1], T2_iso=LIT_T2[2], b_x=lit_bx, b_y=lit_by,
               ratio_x=float(np.median([r["b_x"] for r in allr])) / lit_bx, ratio_y=float(np.median([r["b_y"] for r in allr])) / lit_by)
    rng = np.random.default_rng(20261031)
    arr = {sb: np.array([r["psi"] for r in per[sb]]) for sb in per}
    boot = [float(np.median(np.concatenate([a[rng.integers(0, len(a), len(a))] for a in arr.values()]))) for _ in range(2000)]
    return dict(computable=True, overall=summ([r for v in per.values() for r in v]), ci90=[float(np.percentile(boot, 5)), float(np.percentile(boot, 95))],
                by_participant={sb: summ(v) for sb, v in per.items()}, literature=lit,
                units={sb: [dict(b_x=r["b_x"], b_y=r["b_y"], psi=r["psi"]) for r in v] for sb, v in per.items()})


def env_analyses(work: pathlib.Path, out: pathlib.Path) -> dict:
    """rules_env e1 · e1_fix · e4 · e5 · e6 — the same code with only the work folder changed (same order as environment_recompute · s1b.json first)."""
    out.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(REPO / "decomp" / "results" / "rules_env" / "s1b.json", out / "s1b.json")
    env = dict(ENV, REMEASURE_WORK=str(work), PYTHONDONTWRITEBYTECODE="1", OMP_NUM_THREADS="8")
    env.pop("PYTHONPATH", None)
    rw = pathlib.Path(os.environ.get("S73_RULES_ENV_WORK", str(R.ROOT / "data" / "rules_env" / "work")))
    res = {}
    for name in ("e1", "e1_fix", "e4", "e5", "e6"):
        r = subprocess.run([sys.executable, "-B", "-m", "rules_env.analyses", name, str(out), str(rw)], cwd=REPO / "decomp", env=env,
                           capture_output=True, text=True)
        res[name] = dict(exit=r.returncode, tail=(r.stdout + r.stderr)[-600:])
        p = out / f"{name}.json"
        if r.returncode == 0 and p.exists():
            res[name]["result"] = json.loads(p.read_text(encoding="utf-8"))
    return res


def alt_tree(base: pathlib.Path, path: str, rois: dict, keep: dict) -> pathlib.Path:
    """Work folder of path B or C — wong/<sub>/{DONE.json (따로_정합 = []), roi_<series>.json} · onharmony and tong are links to the actual folders."""
    root = base / f"work_{path}"
    if root.exists():
        shutil.rmtree(root)
    (root / "wong").mkdir(parents=True)
    for d in ("onharmony", "tong"):
        (root / d).symlink_to(R.WORK / d, target_is_directory=True)
    for sub, ser in rois.items():
        wd = root / "wong" / sub
        wd.mkdir()
        done = json.loads((R.WORK / "wong" / sub / "DONE.json").read_text(encoding="utf-8"))
        done.update(따로_정합=[], S73_경로=path, S73_쓴_계열=sorted(s for s in ser if keep[sub].get(s, True)))
        (wd / "DONE.json").write_text(json.dumps(done, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
        for s, r in ser.items():
            if keep[sub].get(s, True):
                (wd / f"roi_{s}.json").write_text(json.dumps(r, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    return root


def leaves(o, pre=""):
    """Leaves of the JSON (path → value)."""
    if isinstance(o, dict):
        for k, v in o.items():
            yield from leaves(v, f"{pre}/{k}")
    elif isinstance(o, list):
        for n, v in enumerate(o):
            yield from leaves(v, f"{pre}[{n}]")
    else:
        yield pre, o


def same_numbers(new: dict, old: dict, tol: float = 1e-9) -> dict:
    """The same comparison as environment_recompute — the numbers of the common leaves agree within a relative 1e-9, and no leaf may exist on one side only."""
    a, b = dict(leaves(new)), dict(leaves(old))
    diff = [k for k in set(a) & set(b) if (isinstance(a[k], (int, float)) and isinstance(b[k], (int, float)) and not isinstance(a[k], bool)
                                          and abs(a[k] - b[k]) > tol * max(1.0, abs(b[k]))) or (not isinstance(a[k], (int, float)) and a[k] != b[k])]
    return dict(common=len(set(a) & set(b)), differ=sorted(diff)[:20], n_differ=len(diff), only_new=sorted(set(a) - set(b))[:20],
                only_stored=sorted(set(b) - set(a))[:20], ok=not diff and set(a) == set(b))


def env_summary(env: dict) -> dict:
    e1, ef, e4 = (env.get(n, {}).get("result") for n in ("e1", "e1_fix", "e4"))
    out = {}
    if e1:
        out["tau_pooled_x"] = e1["pooled"]["x"]["tau"]
        out["by_participant"] = {k.split("|")[1]: dict(tau=v["tau"], beta=v["beta"], beta_ci95=v["beta_ci95"], mean_slope=v["mean_slope"])
                                 for k, v in e1["by"].items() if k.startswith("x|")}
    if ef:
        out["dte_star_x"] = ef["pooled"]["x"]["primary"]                    # main values of the recomputation with the fix (e1_fix §8): dte_star · dte_ci90 · tau
    if e4:
        out["e4_wong_62_132_x_minus_raw"] = e4["primary"]["wong_62_132_x_minus_raw"]
    return out


def main(out_dir: str) -> dict:
    global LOG
    out = pathlib.Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    LOG = ReadLog(out / "s73_read_log.jsonl", "S73 §2 다중 TE 정합 — 정의 보충(-nosearch)에 따른 다시 실행(저자 승인 2026-11-16)")
    res = dict(definition="decomp/plans/S73_POST_HOC_PLAN_2026-11-16.md §2 · decomp/plans/S73_AMENDMENT_2026-11-16.md",
               flirt="-dof 6 -cost normmi -searchcost normmi -nosearch", subjects=list(R.WONG_SUBJECTS), reference=R.WONG_REF)
    dump = lambda: (out / "s73_align.json").write_text(json.dumps(res, ensure_ascii=False, indent=1, default=float) + "\n", encoding="utf-8")  # noqa: E731
    tmp = pathlib.Path(tempfile.mkdtemp(prefix="s73_align_"))
    try:
        # 1. Failure controls (all three participants)
        res["controls"] = {sb: controls(sb, tmp) for sb in R.WONG_SUBJECTS}
        if not all(v["ok"] for v in res["controls"].values()):
            res["stopped"] = "§2.2 실패 대조가 정의 보충의 설정에서도 기대대로 움직이지 않음 — §2 를 끝냄(설정을 더 바꾸지 않음)"
            dump()
            raise SystemExit("⛔ " + res["stopped"])
        # 2. New registration check (nine series × three participants)
        res["check"] = {sb: check_series(sb, tmp) for sb in R.WONG_SUBJECTS}
        dump()
        # 3. Cause check of 2.1 (descriptive)
        res["cause"] = {sb: cause_check(sb) for sb in R.WONG_SUBJECTS}
        dump()
        # 4. Harness check — reproduction of the reference-series ROIs and of rules_env on path A
        ref_roi = {}
        for sb in R.WONG_SUBJECTS:
            ref_roi[sb] = roi_path(sb, R.WONG_REF, None, tmp)
            stored = json.loads((R.WORK / "wong" / sb / f"roi_{R.WONG_REF}.json").read_text(encoding="utf-8"))
            res.setdefault("reference_reproduced", {})[sb] = json.dumps(ref_roi[sb]["주"], sort_keys=True) == json.dumps(stored["주"], sort_keys=True)
        if not all(res["reference_reproduced"].values()):
            res["stopped"] = "새 코드의 기준 계열 ROI 가 저장된 roi 표와 다름 — 하네스 오류"
            dump()
            raise SystemExit("⛔ " + res["stopped"])
        envA = env_analyses(R.WORK, out / "env_A")
        stored_dir = REPO / "decomp" / "results" / "rules_env"
        res["harness_A"] = {n: (same_numbers(v["result"], json.loads((stored_dir / f"{n}.json").read_text(encoding="utf-8")))
                                if "result" in v else dict(ok=False, exit=v["exit"], tail=v["tail"])) for n, v in envA.items()}
        dump()
        if not all(v["ok"] for v in res["harness_A"].values()):
            res["stopped"] = "경로 A 에서 rules_env 가 저장된 결과를 재현하지 못함 — 하네스 오류"
            dump()
            raise SystemExit("⛔ " + res["stopped"])
        # 5. ROIs of paths B · C
        rois = {"B": {sb: {R.WONG_REF: ref_roi[sb]} for sb in R.WONG_SUBJECTS}, "C": {sb: {R.WONG_REF: ref_roi[sb]} for sb in R.WONG_SUBJECTS}}
        from concurrent.futures import ThreadPoolExecutor
        jobs = [(p, sb, s) for sb in R.WONG_SUBJECTS for s in R.WONG_SERIES if s != R.WONG_REF for p in ("B", "C")]
        with ThreadPoolExecutor(max_workers=int(os.environ.get("S73_THREADS", "8"))) as ex:
            got = list(ex.map(lambda j: roi_path(j[1], j[2], tmp / f"{j[1]}_{j[2]}.mat" if j[0] == "C" else None, tmp), jobs))
        for (p, sb, s), r_ in zip(jobs, got):
            rois[p][sb][s] = r_
        keepB = {sb: {s: (s == R.WONG_REF or res["check"][sb][s]["aligned"]) for s in R.WONG_SERIES} for sb in R.WONG_SUBJECTS}
        keepC = {sb: {s: True for s in R.WONG_SERIES} for sb in R.WONG_SUBJECTS}
        # 6. Computability (amendment item 2)
        def avail(keep):
            per = {}
            for sb, k in keep.items():
                n_te = sum(k[f"TE{t}"] for t in JW.TES)
                reps = [r for r, o in JW.REPEATS.items() if k[r] and k[o]]
                per[sb] = dict(n_te=n_te, repeat_pairs=reps, linear=n_te >= 3 and len(reps) >= 1, full=all(k.values()))
            return dict(by_participant=per, linear_all=all(v["linear"] for v in per.values()), full_all=all(v["full"] for v in per.values()))
        res["computable"] = {"B": avail(keepB), "C": avail(keepC)}
        trees = {"A": R.WORK, "B": alt_tree(out, "B", rois["B"], keepB), "C": alt_tree(out, "C", rois["C"], keepC)}
        res["paths"] = {}
        for p, work in trees.items():
            data, own, missing = JW.load(work)
            j = JW.judge(data)
            comp = dict(linear_all=True, full_all=True) if p == "A" else res["computable"][p]
            pr = dict(series_used={sb: sorted(v) for sb, v in data.items()}, missing=missing)
            pr["linear"] = (dict(overall={co: j["overall"][co]["linear_ratio"] for co in ("x", "y")},
                                 by_participant={sb: {co: j["by_subject"][sb][co]["linear_ratio"] for co in ("x", "y")} for sb in j["by_subject"]})
                            if comp["linear_all"] else "not computable")
            if comp["full_all"]:
                pr["odi_range_ratio"] = {sb: j["by_subject"][sb]["ODI"]["te_range_median"] / j["by_subject"][sb]["ODI"]["repeat_abs_diff_median"]
                                         for sb in j["by_subject"]}
                pr["psi"] = psi_summary(data)
                env = envA if p == "A" else env_analyses(work, out / f"env_{p}")
                pr["env_exit"] = {n: v["exit"] for n, v in env.items()}
                pr["env"] = env_summary(env)
            else:
                pr["odi_range_ratio"] = pr["psi"] = pr["env"] = "not computable"
            pr["judge_wong"] = j
            res["paths"][p] = pr
            dump()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return res


if __name__ == "__main__":
    r = main(sys.argv[1])
    print(json.dumps(dict(controls={k: v["ok"] for k, v in r["controls"].items()}, computable=r.get("computable"),
                          aligned={sb: {s: v["aligned"] for s, v in c.items()} for sb, c in r["check"].items()}), ensure_ascii=False))
