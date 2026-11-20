"""design_sim.run — 1,000 replicates per cell and the metrics (ANALYSIS_PLAN.md Section 5); CPU only, at most 8 concurrent processes; --pilot only measures time on a few cells and replicates (no metric values written)."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import pathlib
import sys
import time

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_k, "1")
import numpy as np  # noqa: E402

if __package__ in (None, ""):
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from design_sim import analyses as A  # noqa: E402
from design_sim import model as M  # noqa: E402

KEYS = ("M0", "M1", "M1c", "M1r", "M2", "M3_te0", "M3_tec", "M4_te0", "M5a", "M5b")
SLOPES = ("M3_slope", "M4_slope")


def one(cell: dict, rng: np.random.Generator) -> dict:
    d = M.generate(cell, rng)
    y, g, site, te, te_s = d["y"], d["g"], d["site"], d["te"], d["te_s"]
    k = len(M.ACQ)
    te_defined = float(np.ptp(te_s)) > 0.0          # D0: if all site TEs are equal, group × TE and the TE regression are not defined (plan «Change record» 2)
    out = {"M0": A.m0(y, g), "M1": A.m1(y, g, site), "M1c": A.m1c(y, g, site), "M1r": A.m1(y, g, site, raw=True),
           "M2": A.m2(y, g, site, k), "M5a": A.m5(y, g, site, k, te, d["r2"]), "M5b": A.m5(y, g, site, k, te, d["r2_meas"])}
    if te_defined:
        m3 = A.m3(y, g, site, k, te, M.TE_C)
        m4 = A.m4(y, g, site, k, te_s)
        out.update({"M3_te0": m3["te0"], "M3_tec": m3["tec"], "M3_slope": m3["slope"], "M4_te0": m4["te0"], "M4_slope": m4["slope"]})
    rec = {}
    for key in KEYS:
        if key not in out:
            rec[key] = None
            continue
        tgt = d["theta_raw"] if key == "M1r" else d["theta"]
        e = out[key]
        rec[key] = (e["est"], e["se"], A.covered(e, tgt), A.rejected(e), tgt)
    for key in SLOPES:
        rec[key] = (out[key]["est"], out[key]["se"], A.rejected(out[key])) if key in out else None
    rec["M4_tau0"] = m4["tau0"] if te_defined else None
    return rec


def metrics(recs: list[dict]) -> dict:
    R = len(recs)
    out = {}
    for key in KEYS:
        if recs[0][key] is None:
            out[key] = dict(defined=False)
            continue
        est = np.array([r[key][0] for r in recs])            # (R, ROI)
        se = np.array([r[key][1] for r in recs])
        cov = np.array([r[key][2] for r in recs], float)
        rej = np.array([r[key][3] for r in recs], float)
        tgt = np.array([np.broadcast_to(r[key][4], est.shape[1]) for r in recs])
        err = est - tgt
        er_rep, cv_rep, rj_rep = err.mean(1), cov.mean(1), rej.mean(1)
        empsd_r = est.std(0, ddof=1)
        out[key] = dict(bias=float(er_rep.mean()), bias_mcse=float(er_rep.std(ddof=1) / np.sqrt(R)),
                        empsd=float(empsd_r.mean()), empsd_mcse=float(empsd_r.mean() / np.sqrt(2 * (R - 1))),
                        mean_se=float(se.mean()), rmse=float(np.sqrt((err ** 2).mean())),
                        coverage=float(cv_rep.mean()), coverage_mcse=float(cv_rep.std(ddof=1) / np.sqrt(R)),
                        reject=float(rj_rep.mean()), reject_mcse=float(rj_rep.std(ddof=1) / np.sqrt(R)),
                        target=float(tgt.mean()),
                        by_roi=dict(bias=err.mean(0).round(6).tolist(), coverage=cov.mean(0).round(4).tolist(),
                                    reject=rej.mean(0).round(4).tolist()))
    for key in SLOPES:
        if recs[0][key] is None:
            out[key] = dict(defined=False)
            continue
        est = np.array([r[key][0] for r in recs])
        rej = np.array([r[key][2] for r in recs], float).mean(1)
        out[key] = dict(mean_est=float(est.mean()), reject=float(rej.mean()), reject_mcse=float(rej.std(ddof=1) / np.sqrt(R)))
    if recs[0]["M4_tau0"] is None:
        out["M4_tau0_share"] = None
    else:
        t0 = np.array([r["M4_tau0"] for r in recs], float).mean(1)
        out["M4_tau0_share"] = float(t0.mean())
    return out


def run_cell(args) -> dict:
    cell, reps, outdir, pilot = args
    rng = np.random.default_rng(np.random.SeedSequence([M.SEED, cell["idx"]]))
    t = time.time()
    recs = [one(cell, rng) for _ in range(reps)]
    sec = time.time() - t
    if pilot:
        return dict(idx=cell["idx"], reps=reps, seconds=sec)
    res = dict(cell=cell, reps=reps, seconds=round(sec, 2), seed=[M.SEED, cell["idx"]], theta=M.theta_of(cell["S"]),
               metrics=metrics(recs))
    p = pathlib.Path(outdir) / "cells" / f"cell_{cell['idx']:03d}.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")
    return dict(idx=cell["idx"], seconds=sec)


def code_hashes() -> dict:
    here = pathlib.Path(__file__).resolve().parent
    return {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(here.glob("*.py"))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--cells", default="0-287")
    ap.add_argument("--pilot", action="store_true")
    a = ap.parse_args()
    if a.procs > 8:
        raise SystemExit("동시 프로세스는 8 개 이하(계획 8 절)")
    lo, hi = (int(x) for x in a.cells.split("-"))
    cs = [c for c in M.cells() if lo <= c["idx"] <= hi]
    reps = 10 if a.pilot else M.REPS
    if a.pilot:
        cs = [c for c in cs if c["idx"] in (0, 143, 287)]
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    import multiprocessing as mp
    with mp.get_context("fork").Pool(a.procs) as pool:
        done = pool.map(run_cell, [(c, reps, str(out), a.pilot) for c in cs], chunksize=1)
    rec = dict(start=time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(t0)), end=time.strftime("%Y-%m-%d %H:%M:%S"),
               wall_seconds=round(time.time() - t0, 1), procs=a.procs, cells=[lo, hi], pilot=a.pilot, reps=reps,
               per_cell_seconds=[d["seconds"] for d in done], code_sha256=code_hashes(),
               numpy=np.__version__)
    name = "pilot_timing.json" if a.pilot else f"run_record_{lo:03d}_{hi:03d}.json"
    (out / name).write_text(json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")
    print(json.dumps({k: rec[k] for k in ("wall_seconds", "procs", "reps", "pilot")}))


if __name__ == "__main__":
    main()
