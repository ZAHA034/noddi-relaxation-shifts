#!/usr/bin/env python
"""Step 2 of the reproduction: recompute the fixed tests from the deposited tables.

Input: the ROI tables in data/remeasure and, with --judge, the frozen weights and test
data in data/s4 (from the Zenodo deposit). No image is read.

Without --judge this takes minutes. With --judge it also runs the full calibration test,
which took about 24 h on CPU; decomp/v4_recompute/ recomputes a part of it.

The script copies the input tables to a temporary folder without the deposited
outputs, runs the unchanged analysis code there, and compares every number in each
new output with the deposited output. It never overwrites data/.

Usage: python reproduce/step2_recompute.py [--judge] [--source-meta] [--threads N] [--keep DIR]
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
OUTPUTS = ("onharmony/pairs", "onharmony/v4", "tong/judge", "wong/judge", "explore")


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
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--judge", action="store_true", help="also run the full calibration test (about 24 h on CPU; needs data/s4)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--source-meta", action="store_true", help="with --judge: recompute the attribute table; needs the ON-Harmony sidecars")
    ap.add_argument("--keep", default="", help="folder to keep the recomputed outputs in")
    args = ap.parse_args()

    tmp = pathlib.Path(args.keep or tempfile.mkdtemp(prefix="recompute_"))
    work = tmp / "remeasure"
    if work.exists():
        raise SystemExit(f"{work} already exists")
    ignore = lambda d, names: [n for n in names if os.path.relpath(os.path.join(d, n), DATA / "remeasure") in OUTPUTS]
    shutil.copytree(DATA / "remeasure", work, ignore=ignore)

    env = dict(os.environ, REMEASURE_WORK=str(work), DECOMP_WORK=str(DATA), CUDA_VISIBLE_DEVICES="", PYTHONDONTWRITEBYTECODE="1",
               OMP_NUM_THREADS=str(args.threads), OPENBLAS_NUM_THREADS=str(args.threads), MKL_NUM_THREADS=str(args.threads),
               MPLCONFIGDIR=str(tmp / "mpl"))

    def run(*cmd):
        t = time.time()
        r = subprocess.run([sys.executable, *cmd], cwd=REPO, env=env, capture_output=True, text=True)
        if r.returncode != 0:
            sys.stderr.write(r.stderr[-3000:])
            raise SystemExit(f"failed: {' '.join(cmd)}")
        return r.stdout, time.time() - t

    results = {}

    def check(name, new_obj, deposited, secs, new_path=None):
        old_obj = json.loads(deposited.read_text(encoding="utf-8"))
        for o in (new_obj, old_obj):                  # a set in the code: its printed order varies from run to run
            if isinstance(o, dict) and "own_registration" in o:
                o["own_registration"] = {k: sorted(v) for k, v in o["own_registration"].items()}
        res = compare(new_obj, old_obj)
        res["seconds"] = round(secs, 1)
        if new_path is not None:
            res["bytes_identical"] = sha256(new_path) == sha256(deposited)
        results[name] = res
        print(name, json.dumps(res, ensure_ascii=False), flush=True)

    out, s = run("reproduce/tong_judge.py")
    check("tong", json.loads(out), DATA / "remeasure/tong/judge/tong_judge.json", s)

    out, s = run("decomp/remeasure/judge_wong.py", str(work))
    check("wong", json.loads(out), DATA / "remeasure/wong/judge/wong_judge.json", s)

    _, s = run("decomp/remeasure/judge_onh.py")
    p = work / "onharmony/pairs/onh_pairs.json"
    check("onharmony_pairs", json.loads(p.read_text(encoding="utf-8")), DATA / "remeasure/onharmony/pairs/onh_pairs.json", s, p)

    v4 = work / "onharmony/v4"
    out, s = run("decomp/remeasure/judge_v4.py", "snr-table")
    snr_sha = out.split()[-1]
    p = v4 / "snr_table.json"
    check("snr_table", json.loads(p.read_text(encoding="utf-8")), DATA / "remeasure/onharmony/v4/snr_table.json", s, p)

    if args.judge:
        missing = [f for f in ("models/c3_0.obj.pt", "sim_test.npz", "v4/test_cover_c3.npz") if not (DATA / "s4" / f).exists()]
        if missing:
            raise SystemExit(f"--judge needs data/s4 from the Zenodo deposit; missing {missing}. See README.md")
        p = v4 / "attr_table.json"
        if args.source_meta:                          # the attribute table needs the ON-Harmony sidecars (public, OpenNeuro ds004712)
            out, s = run("decomp/remeasure/judge_v4.py", "attr-table", snr_sha)
            attr_sha = out.split()[-1]
            check("attr_table", json.loads(p.read_text(encoding="utf-8")), DATA / "remeasure/onharmony/v4/attr_table.json", s, p)
        else:                                         # default: start from the deposited table
            shutil.copyfile(DATA / "remeasure/onharmony/v4/attr_table.json", p)
            attr_sha = sha256(p)
            print("attr_table: deposited table used (run with --source-meta to recompute it from the sidecars)", flush=True)
        out, s = run("reproduce/judge_v4_run.py", "run", "c3", attr_sha)
        p = v4 / "v4_c3.json"
        check("v4_c3", json.loads(p.read_text(encoding="utf-8")), DATA / "remeasure/onharmony/v4/v4_c3.json", s, p)

    bad = {k: v for k, v in results.items() if v["differ"]}
    print("RESULT", "all compared values agree" if not bad else f"values differ in: {sorted(bad)}")
    (tmp / "step2_results.json").write_text(json.dumps(results, ensure_ascii=False, indent=1), encoding="utf-8")
    print("outputs kept in", tmp)
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
