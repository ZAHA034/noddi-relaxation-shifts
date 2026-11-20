#!/usr/bin/env python
"""Recompute the measurement-function environment (Supplementary Methods S-M10) and compare with the deposit.

The analyses live in the package decomp/rules_env/. This script runs them in a temporary
folder with the unchanged code and compares every number in each new output with the
stored output decomp/results/rules_env/<name>.json. It needs no setup outside the
repository: it makes a temporary folder T, clones this repository into T/repo (the code
checks the Git history), links T/dmri_work/decomp_stage1/protocols to this repository's
data/protocols (the code looks for the protocol files there), copies data/ into T and
runs everything from T/repo. T is deleted at the end unless --keep is given. The files in
this working tree are never written (data/, decomp/results/ and so on). Only committed
files are cloned, so run it on a clean checkout.

Default (minutes, CPU only): e1 e2 e3 e4 e5 e6 s3 e1_fix, in this order, into one temporary
output folder. The stored s1b.json is copied there first because e1 reads fields from it,
and e1_fix reads the e1 output. Analysis e6 needs the stored synthetic-signal fits
(s1b/fits.npz and s1b/jobs.json) from the Zenodo file unpacked into data/rules_env_work/
(option --work). If that folder is missing, e6 is skipped and a message says so.

--fits   First re-run the three analyses that fit synthetic signals with AMICO
         (s1b, s2, k8). They call "uv run --frozen" in envs/reference/, so they need
         the AMICO environment of envs/reference/ (uv installed, the locked
         environment already built, or network access for uv to build it). This takes
         many processes (5–8 minutes on 16 cores). The new S1b fits are used for e6 only if their SHA-256
         equal those of the deposited fits; otherwise this is reported and e6 uses the
         deposited work folder.
--stage3 DIR
         Also run s1a. It reads the stored stage-3 simulation fits under
         DIR/E6/P2a and DIR/E6/P2b (files DONE, rows.npz and f2_1.7.npz in each c* folder).
         Unpack the Zenodo file stage3_fits_E6.zip to data/stage3_fits and pass
         --stage3 data/stage3_fits.

Report per analysis: common (numeric and other leaves present in both), differ,
only_new, only_deposited, bytes_identical. Numbers agree if they are equal or within a
relative tolerance of 1e-9. The exit status is 1 if any analysis fails or differs.

Usage: python reproduce/environment_recompute.py [--work DIR] [--fits] [--stage3 DIR]
                                                 [--threads N] [--keep DIR]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile
import time

REPO = pathlib.Path(__file__).resolve().parents[1]
DATA = REPO / "data"
DEPOSITED = REPO / "decomp" / "results" / "rules_env"
DEFAULT_ORDER = ("e1", "e2", "e3", "e4", "e5", "e6", "s3", "e1_fix")
FIT_ORDER = ("s1b", "s2", "k8")
FIT_FILES = ("fits.npz", "jobs.json")


def sha256(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def leaves(obj, prefix=""):
    """Flatten JSON into {path: value}."""
    if isinstance(obj, dict):
        out = {}
        for k, v in obj.items():
            out.update(leaves(v, f"{prefix}/{k}"))
        return out
    if isinstance(obj, list):
        out = {}
        for i, v in enumerate(obj):
            out.update(leaves(v, f"{prefix}[{i}]"))
        return out
    return {prefix: obj}


def same(a, b, tol):
    if isinstance(a, bool) or isinstance(b, bool) or not isinstance(a, (int, float)) or not isinstance(b, (int, float)):
        return a == b
    if isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b):
        return True
    return a == b or abs(a - b) <= tol * max(abs(a), abs(b), 1e-300)


def compare(new, old, tol=1e-9):
    """Compare the leaves present in both outputs. Numbers must agree within a relative tolerance."""
    a, b = leaves(new), leaves(old)
    common = sorted(set(a) & set(b))
    differ = [k for k in common if not same(a[k], b[k], tol)]
    numeric = sum(isinstance(a[k], (int, float)) and not isinstance(a[k], bool) for k in common)
    return dict(common=len(common), numeric=numeric, differ=len(differ), only_new=len(set(a) - set(b)),
                only_deposited=len(set(b) - set(a)), examples=differ[:5])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", default=str(DATA / "rules_env_work"),
                    help="folder with the stored work files (s1b/, k8/, s2/) from the Zenodo file; default data/rules_env_work")
    ap.add_argument("--fits", action="store_true",
                    help="first re-run s1b, s2 and k8 (AMICO fits through 'uv run --frozen' in envs/reference; 5–8 minutes on 16 cores)")
    ap.add_argument("--stage3", default="", metavar="DIR",
                    help="folder with the stored stage-3 simulation fits (E6/P2a, E6/P2b), i.e. the unpacked Zenodo file "
                         "stage3_fits_E6.zip (data/stage3_fits); also runs s1a")
    ap.add_argument("--threads", type=int, default=0,
                    help="limit the run to N CPU cores and N BLAS threads (default: all cores available to the process)")
    ap.add_argument("--keep", default="", help="folder to keep the recomputed outputs in (must not exist yet)")
    args = ap.parse_args()

    if args.keep:
        tmp = pathlib.Path(args.keep).resolve()
        if tmp.exists():
            raise SystemExit(f"{tmp} already exists")
        tmp.mkdir(parents=True)
    else:
        tmp = pathlib.Path(tempfile.mkdtemp(prefix="env_recompute_"))
    try:
        return body(args, tmp)
    finally:
        if not args.keep:
            shutil.rmtree(tmp, ignore_errors=True)


def body(args, tmp) -> int:
    out_dir, scratch, fresh_work = tmp / "out", tmp / "scratch", tmp / "work_new"
    out_dir.mkdir()
    scratch.mkdir()
    repo = tmp / "repo"
    r = subprocess.run(["git", "clone", "-q", str(REPO), str(repo)], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"git clone of {REPO} failed (this script needs the Git repository):\n{r.stderr[-1500:]}")
    proto = tmp / "dmri_work" / "decomp_stage1" / "protocols"
    proto.parent.mkdir(parents=True)
    proto.symlink_to(DATA / "protocols", target_is_directory=True)
    venv = REPO / "envs" / "reference" / ".venv"
    if venv.exists() and not (repo / "envs" / "reference" / ".venv").exists():
        (repo / "envs" / "reference" / ".venv").symlink_to(venv, target_is_directory=True)

    if args.threads > 0 and hasattr(os, "sched_setaffinity"):
        cpus = sorted(os.sched_getaffinity(0))[:args.threads]
        os.sched_setaffinity(0, cpus)                 # child processes inherit this
    nthr = str(args.threads) if args.threads > 0 else str(len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count())

    # Private copy of the deposited input tables. The deposited outputs of this analysis
    # are never used as inputs, except the stored s1b.json that e1 reads (seeded below).
    data_copy = tmp / "data"
    shutil.copytree(DATA, data_copy, ignore=shutil.ignore_patterns("rules_env_work", "s4", "stage3_fits"))   # stage-3 fits are read in place (--stage3)
    if args.stage3:
        s3dir = pathlib.Path(args.stage3).resolve()
        if not (s3dir / "E6").is_dir():
            raise SystemExit(f"--stage3: {s3dir}/E6 not found")
        shutil.rmtree(data_copy / "fits", ignore_errors=True)
        (data_copy / "fits").symlink_to(s3dir, target_is_directory=True)

    env = dict(os.environ, REMEASURE_WORK=str(data_copy / "remeasure"), DECOMP_WORK=str(data_copy), CUDA_VISIBLE_DEVICES="",
               PYTHONDONTWRITEBYTECODE="1", OMP_NUM_THREADS=nthr, OPENBLAS_NUM_THREADS=nthr, MKL_NUM_THREADS=nthr,
               MPLCONFIGDIR=str(tmp / "mpl"))
    env.pop("PYTHONPATH", None)

    results, notes = {}, []

    def run(name, work):
        t = time.time()
        r = subprocess.run([sys.executable, "-B", "-m", "rules_env.analyses", name, str(out_dir), str(work)],
                           cwd=repo / "decomp", env=env, capture_output=True, text=True)
        secs = time.time() - t
        if r.returncode != 0:
            sys.stderr.write(f"--- {name} failed (exit {r.returncode}); last output:\n{(r.stdout + r.stderr)[-2500:]}\n")
            log = work / name / "simrun.log"
            if log.exists():
                sys.stderr.write(f"--- tail of {log}:\n{log.read_text(errors='replace')[-1500:]}\n")
        return r.returncode, secs

    def check(name, secs, status="ok"):
        res = dict(status=status, seconds=round(secs, 1))
        new_p, old_p = out_dir / f"{name}.json", DEPOSITED / f"{name}.json"
        if status == "ok":
            if not new_p.exists():
                res["status"] = "failed"
            else:
                res.update(compare(json.loads(new_p.read_text(encoding="utf-8")), json.loads(old_p.read_text(encoding="utf-8"))))
                res["bytes_identical"] = sha256(new_p) == sha256(old_p)
        results[name] = res
        print(name, json.dumps(res, ensure_ascii=False), flush=True)

    work_dep = pathlib.Path(args.work).resolve()
    work_for_e6 = work_dep
    shutil.copyfile(DEPOSITED / "s1b.json", out_dir / "s1b.json")      # e1 adds fitter-only fields from this file

    if args.fits:
        for name in FIT_ORDER:
            rc, secs = run(name, fresh_work)
            check(name, secs, "ok" if rc == 0 else "failed")
        new = fresh_work / "s1b"
        old = work_dep / "s1b"
        if all((new / f).exists() for f in FIT_FILES) and all((old / f).exists() for f in FIT_FILES):
            if all(sha256(new / f) == sha256(old / f) for f in FIT_FILES):
                work_for_e6 = fresh_work
                notes.append("new S1b fits have the same SHA-256 as the deposited fits; e6 uses the new fits")
            else:
                notes.append("new S1b fits differ in SHA-256 from the deposited fits (numbers may still agree, see the s1b line); "
                             "e6 uses the deposited work folder")
        elif all((new / f).exists() for f in FIT_FILES):
            work_for_e6 = fresh_work
            notes.append("no deposited work folder to compare with; e6 uses the new S1b fits and checks their SHA-256 itself")
        else:
            notes.append("new S1b fits are not available; e6 uses the deposited work folder")

    for name in DEFAULT_ORDER:
        if name == "e6" and not (work_for_e6 / "s1b" / "fits.npz").exists():
            msg = (f"e6 skipped: the stored S1b fits are not found in {work_for_e6}/s1b. Unpack the Zenodo work file into "
                   f"data/rules_env_work/ or pass --work DIR.")
            print(msg, flush=True)
            notes.append(msg)
            check(name, 0.0, "skipped")
            continue
        rc, secs = run(name, work_for_e6 if name == "e6" else scratch)
        check(name, secs, "ok" if rc == 0 else "failed")

    if args.stage3:
        rc, secs = run("s1a", scratch)
        check("s1a", secs, "ok" if rc == 0 else "failed")

    for n in notes:
        print("note:", n)
    bad = {k: v for k, v in results.items() if v["status"] == "failed" or v.get("differ") or v.get("only_new") or v.get("only_deposited")}
    skipped = sorted(k for k, v in results.items() if v["status"] == "skipped")
    ident = sorted(k for k, v in results.items() if v.get("bytes_identical"))
    print("bytes identical:", ident)
    if skipped:
        print("not run:", skipped)
    print("RESULT", "all compared values agree" if not bad else f"differences or failures in: {sorted(bad)}")
    (tmp / "environment_results.json").write_text(json.dumps(dict(results=results, notes=notes), ensure_ascii=False, indent=1), encoding="utf-8")
    print("outputs kept in" if args.keep else "outputs in (deleted when this script ends; use --keep DIR to keep them)", out_dir)
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
