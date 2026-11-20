"""S73 §5 (post hoc): session SNR of the multi-TE series by participant and TE, noise-floor simulation of the ROI-matched states (noncentral χ, L = 1, 4, 16) with failure controls, and the summaries per replicate."""
from __future__ import annotations

import json
import math
import os
import pathlib
import subprocess
import sys

import numpy as np

REPO = pathlib.Path(os.environ.get("S73_REPO") or pathlib.Path(__file__).resolve().parents[3])
HERE = pathlib.Path(__file__).resolve().parent
sys.path[:0] = [str(REPO / "decomp"), str(REPO / "decomp" / "remeasure"), str(REPO / "decomp" / "stage1"), str(REPO / "decomp" / "v4_recompute")]
import analysis as AN  # noqa: E402
import judge_wong as JW  # noqa: E402
import lib  # noqa: E402
import rconfig as R  # noqa: E402
from read_log import ReadLog  # noqa: E402
from rules_env import analyses as RA  # noqa: E402
from rules_env import empirical as EM  # noqa: E402
from rules_env import metrics as MT  # noqa: E402

TES = RA.WONG_TES
LS = (1, 4, 16)
N_PER_FDIR, N_REP = 50, 50
SEED0 = 20261116
B0_THR = 10.0


def snr_table(log: ReadLog) -> dict:
    import nibabel as nib
    out = {}
    for sb in R.WONG_SUBJECTS:
        wd = R.WORK / "wong" / sb
        log.before_read([wd / "nodif_brain_mask.nii.gz"])
        msk = np.asarray(nib.load(str(wd / "nodif_brain_mask.nii.gz")).dataobj) > 0
        for t in TES:
            s = f"TE{t}"
            sd = wd if s == R.WONG_REF else wd / f"own_{s}"
            p = R.WONG / sb / "dwi" / f"{sb}_acq-{s}_dwi"
            fa_p = sd / f"f2_{s}_FA.nii.gz"
            log.before_read([pathlib.Path(f"{p}.bval"), pathlib.Path(f"{p}.nii.gz"), fa_p])
            bv = np.loadtxt(f"{p}.bval").ravel()
            img = nib.load(f"{p}.nii.gz")
            b0 = np.asarray(img.dataobj, dtype=np.float32)[..., bv <= B0_THR]
            fa = np.asarray(nib.load(str(fa_p)).dataobj, dtype=float)
            out.setdefault(sb, {})[str(t)] = AN.snr(b0, msk, fa)
    return out


def run_sim(jobs: list[dict], work: pathlib.Path, tag: str, n_proc: int, nthreads: int) -> dict:
    """Splits the jobs into n_proc groups run at the same time in envs/reference and merges the result npz files."""
    work.mkdir(parents=True, exist_ok=True)
    groups = [jobs[i::n_proc] for i in range(n_proc) if jobs[i::n_proc]]
    procs = []
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", S73_REPO=str(REPO), CUDA_VISIBLE_DEVICES="")
    env.pop("VIRTUAL_ENV", None)
    for gi, g in enumerate(groups):
        jp, op = work / f"{tag}_{gi}.json", work / f"{tag}_{gi}.npz"
        jp.write_text(json.dumps({"jobs": g}), encoding="utf-8")
        log = open(work / f"{tag}_{gi}.log", "w")
        procs.append((subprocess.Popen(["uv", "run", "--frozen", "python", "-I", "-B", str(HERE / "s73_noise_sim.py"), str(jp), str(op), str(work / tag),
                                        str(nthreads)], cwd=REPO / "envs" / "reference", env=env, stdout=log, stderr=subprocess.STDOUT), op, log))
    out = {}
    for p, op, log in procs:
        rc = p.wait()
        log.close()
        if rc != 0:
            raise SystemExit(f"⛔ 시뮬레이션 실패({op.name} · 종료 {rc})")
        out.update(dict(np.load(op)))
    return out


def job(key, sb, t, units, L, sigma, n_per_fdir, n_rep, seed) -> dict:
    t2c, t1c = RA.centre_eta()
    return dict(key=key, card=f"WONG:TE{t}", units=[list(map(float, u)) for u in units], t2=list(t2c), t1=list(t1c), gen="T-in",
                fdirs=[int(f) for f in RA.FDIRS], n_per_fdir=n_per_fdir, n_rep=n_rep, L=L, sigma=sigma, seed=seed)


def coords(ndi, fwf, odi):
    xy, _ = lib.coords(np.asarray(ndi, float).ravel(), np.asarray(fwf, float).ravel(), np.asarray(odi, float).ravel())
    return xy[:, 0].reshape(np.shape(ndi)), xy[:, 1].reshape(np.shape(ndi))


def slopes(V: np.ndarray) -> np.ndarray:
    """V (…, TE) → OLS slope."""
    te = np.asarray(TES, float)
    t = te - te.mean()
    return (V * t).sum(-1) / (t ** 2).sum()


def summarize(res: dict, units_by: dict, L: int) -> dict:
    """Per Monte Carlo replicate, unit slopes and levels → τ (with and without replicate noise subtracted) · β · share by the linear rule · ODI range ratio · ψ."""
    out = dict(by_participant={})
    kc_all, var_all, xs_all = [], [], []
    for sb, keys in units_by.items():
        nd = np.stack([res[f"{sb}|TE{t}|L{L}|ndi"] for t in TES], -1)       # (replicate, unit, TE)
        fw = np.stack([res[f"{sb}|TE{t}|L{L}|fwf"] for t in TES], -1)
        od = np.stack([res[f"{sb}|TE{t}|L{L}|odi"] for t in TES], -1)
        x, y = coords(nd, fw, od)
        kx, ky = slopes(x), slopes(y)                                        # (replicate, unit)
        n_rep = kx.shape[0]
        kbar, kvar = kx.mean(0), kx.var(0, ddof=1) if n_rep > 1 else np.zeros(kx.shape[1])
        tau_raw = float(np.std(kbar, ddof=1))
        tau = float(math.sqrt(max(0.0, np.var(kbar, ddof=1) - kvar.mean() / n_rep)))
        xbar = x.mean(-1).mean(0)
        beta = MT.beta_level(kbar, xbar, float(x.mean(-1).var(0, ddof=1).mean() / n_rep) if n_rep > 1 else 0.0)
        lin, odir = {"x": [], "y": []}, []
        psi = []
        for r in range(n_rep):
            r2 = (r + 1) % n_rep
            ser = {f"TE{t}": {u: np.array([x[r, i, ti], y[r, i, ti], od[r, i, ti]]) for i, u in enumerate(keys)} for ti, t in enumerate(TES)}
            ser["TE62R2"] = {u: np.array([x[r2, i, 0], y[r2, i, 0], od[r2, i, 0]]) for i, u in enumerate(keys)}
            ser["TE132R2"] = {u: np.array([x[r2, i, -1], y[r2, i, -1], od[r2, i, -1]]) for i, u in enumerate(keys)}
            j = JW.judge({sb: ser})
            for co in ("x", "y"):
                lin[co].append(j["by_subject"][sb][co]["linear_ratio"])
            o = j["by_subject"][sb]["ODI"]
            odir.append(o["te_range_median"] / o["repeat_abs_diff_median"] if o["repeat_abs_diff_median"] else float("nan"))
            psi.append(float(np.median(np.degrees(np.arctan2(ky[r], kx[r])))))
        out["by_participant"][sb] = dict(n_units=len(keys), tau=tau, tau_uncorrected=tau_raw, beta=beta, mean_slope_x=float(kbar.mean()),
                                         linear_share_x=float(np.mean(lin["x"])) if None not in lin["x"] else None,
                                         linear_share_y=float(np.mean(lin["y"])) if None not in lin["y"] else None,
                                         odi_range_ratio=float(np.nanmean(odir)), psi_median=float(np.mean(psi)))
        kc_all += list(kbar - kbar.mean())
        out.setdefault("_kbar_all", []).extend(kbar.tolist())
        var_all += list(kvar)
        xs_all.append(n_rep)
    n_rep = xs_all[0]
    kc = np.array(kc_all)
    out["pooled"] = dict(tau=float(math.sqrt(max(0.0, np.var(kc, ddof=1) - np.mean(var_all) / n_rep))), tau_uncorrected=float(np.std(kc, ddof=1)),
                         n_units=len(kc), tau_uncentered_all=float(np.std(out.pop("_kbar_all"), ddof=1)))
    return out


def category(ratio: float) -> str:
    return "설명하지 못함" if ratio <= 0.25 else "일부 설명" if ratio < 0.75 else "대부분 설명할 수 있음"


def main(out_dir: str, n_proc: int = 6, variant: str = "memo") -> dict:
    """variant «memo» = as in the implementation note (sub-01 from the stored b1s_snr.json, sub-02 and sub-03 computed here) · «recomputed» = all three participants computed here (sensitivity)."""
    out = pathlib.Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    log = ReadLog(out / "s73_read_log.jsonl", "S73 §5 다중 TE SNR(저자 승인 2026-11-16)")
    res = dict(definition="decomp/plans/S73_POST_HOC_PLAN_2026-11-16.md §5", snr_variant=variant)
    sfx = "" if variant == "memo" else f"_snr_{variant}"
    w = EM.read_wong()
    wm = RA._wong_matched(EM, w)
    units_by = {}
    for (sb, u) in wm["units"]:
        units_by.setdefault(sb, []).append(u)
    res["units"] = {sb: len(v) for sb, v in units_by.items()}
    res["units_dropped"] = wm["dropped"]
    prev = out / "s73_noise_snr_recomputed.json"                         # if an earlier run already read these values, do not read the files again
    snr = json.loads(prev.read_text(encoding="utf-8"))["snr"] if prev.exists() else snr_table(log)
    stored = json.loads((R.WORK / "explore" / "b1s" / "b1s_snr.json").read_text(encoding="utf-8"))["te"]
    res["snr_recomputed"] = snr
    res["snr_sub01_stored_b1s"] = {t: v["snr"] for t, v in stored.items()}
    if variant == "memo":
        snr = {sb: ({t: dict(snr=stored[t]["snr"]) for t in v} if sb == "sub-01" else v) for sb, v in snr.items()}
    res["snr"] = snr
    (out / f"s73_noise{sfx}.json").write_text(json.dumps(res, ensure_ascii=False, indent=1, default=float) + "\n", encoding="utf-8")
    U = {sb: [wm["units"][(sb, u)] for u in units_by[sb]] for sb in units_by}
    U3 = {sb: [(v[0], v[2], v[1]) for v in U[sb]] for sb in U}                # (ν, ODI, f_iso)
    work = out / f"sim_work{sfx}"
    # controls (i) SNR = ∞ · (ii) SNR = 5 · L = 1
    s1b = json.loads((REPO / "decomp" / "results" / "rules_env" / "s1b.json").read_text(encoding="utf-8"))
    tau_s1b = s1b["wong"]["by"]["x|T-in|L0|iii"]["tau_fit"]
    cj = [job(f"{sb}|TE{t}|L1", sb, t, U3[sb], 1, 0.0, 1, 1, SEED0) for sb in U3 for t in TES]
    cj += [job(f"{sb}|TE{t}|SNR5", sb, t, U3[sb], 1, 1000.0 / 5.0, N_PER_FDIR, 1, SEED0 + 7) for sb in U3 for t in TES]
    cr = run_sim(cj, work, "controls", n_proc, 8)
    inf = summarize({k.replace("|L1|", "|L1|"): v for k, v in cr.items() if "SNR5" not in k}, units_by, 1)
    s_raw = {sc: {sb: EM.s_rt_corrected(EM.wong_series(w, sc)[sb])[0] for sb in units_by} for sc in ("raw_ndi", "raw_fiso")}
    moved = {}
    for sb in units_by:
        for t in TES:
            for k, sc in (("ndi", "raw_ndi"), ("fwf", "raw_fiso")):
                d = np.abs(cr[f"{sb}|TE{t}|SNR5|{k}"][0] - cr[f"{sb}|TE{t}|L1|{k}"][0])
                moved[f"{sb}|TE{t}|{k}"] = int((d > 1.96 * math.sqrt(2) * s_raw[sc][sb]).sum())
    ok1 = abs(inf["pooled"]["tau_uncentered_all"] - tau_s1b) <= 0.10 * tau_s1b          # same definition as tau_fit of S1b (81 units · not centered)
    ok2 = sum(moved.values()) > 0
    res["controls"] = dict(i=dict(tau_sim_inf=inf["pooled"]["tau_uncentered_all"], tau_sim_inf_centered=inf["pooled"]["tau"], tau_s1b=tau_s1b, ok=bool(ok1)),
                           ii=dict(units_moved_beyond_retest_width=moved, s_raw=s_raw, ok=bool(ok2)), ok=bool(ok1 and ok2))
    (out / f"s73_noise{sfx}.json").write_text(json.dumps(res, ensure_ascii=False, indent=1, default=float) + "\n", encoding="utf-8")
    if not (ok1 and ok2):
        res["stopped"] = "§5.2 실패 대조가 기대대로 움직이지 않음"
        (out / f"s73_noise{sfx}.json").write_text(json.dumps(res, ensure_ascii=False, indent=1, default=float) + "\n", encoding="utf-8")
        raise SystemExit("⛔ " + res["stopped"])
    jobs = [job(f"{sb}|TE{t}|L{L}", sb, t, U3[sb], L, 1000.0 / snr[sb][str(t)]["snr"], N_PER_FDIR, N_REP,
                SEED0 + 1000 * (i + 1) + L * 10 + ti) for i, sb in enumerate(sorted(U3)) for ti, t in enumerate(TES) for L in LS]
    sim = run_sim(jobs, work, "main", n_proc, max(4, 48 // n_proc))
    np.savez_compressed(out / f"s73_noise_sim{sfx}.npz", **sim)
    res["by_L"] = {str(L): summarize(sim, units_by, L) for L in LS}
    e1f = json.loads((REPO / "decomp" / "results" / "rules_env" / "e1_fix.json").read_text(encoding="utf-8"))
    res["reading_inputs_note"] = "τ_obs 는 정의상 §2 의 경로 A · B · C 중 가장 작은 값 — §2 가 멈춰 경로 A(저장된 e1_fix)만 있다"
    (out / f"s73_noise{sfx}.json").write_text(json.dumps(res, ensure_ascii=False, indent=1, default=float) + "\n", encoding="utf-8")
    return res


if __name__ == "__main__":
    r = main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 6, sys.argv[3] if len(sys.argv) > 3 else "memo")
    print(json.dumps(dict(controls=r["controls"], by_L={L: v["pooled"] for L, v in r["by_L"].items()}), ensure_ascii=False, indent=1, default=float))
