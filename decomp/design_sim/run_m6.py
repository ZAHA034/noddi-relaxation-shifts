"""M6 run (plan PLAN_M6.md) — only when the plan is fixed at HEAD; 16 cells × 1,000 replicates, at most 8 concurrent processes, CPU."""
import argparse
import hashlib
import json
import pathlib
import subprocess
import time

import numpy as np

from design_sim import m6
from design_sim import model as M

SEED2 = 20261118


def run_cell(args):
    cell, reps, outdir = args
    rng = np.random.default_rng(np.random.SeedSequence([M.SEED, cell["idx"]]))
    rng2 = np.random.default_rng(np.random.SeedSequence([SEED2, cell["idx"]]))
    t = time.time()
    recs = [m6.one(cell, rng, rng2) for _ in range(reps)]
    res = dict(cell=cell, reps=reps, seconds=round(time.time() - t, 2), seeds=[[M.SEED, cell["idx"]], [SEED2, cell["idx"]]],
               dte=m6.DTE, noise2=m6.NOISE2, theta=M.theta_of(cell["S"]), metrics=m6.summarize(recs))
    p = pathlib.Path(outdir) / "m6_cells" / f"cell_{cell['idx']:03d}.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")
    return cell["idx"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--procs", type=int, default=8)
    a = ap.parse_args()
    if a.procs > 8:
        raise SystemExit("동시 프로세스는 8 개 이하")
    here = pathlib.Path(__file__).resolve().parent
    head = subprocess.run(["git", "show", "HEAD:decomp/design_sim/PLAN_M6.md"], cwd=here.parents[1], capture_output=True)
    if head.returncode != 0 or head.stdout != (here / "PLAN_M6.md").read_bytes():
        raise SystemExit("계획이 HEAD 에 고정되어 있지 않다 — 계산하지 않는다")
    cs = m6.cells_m6()
    assert len(cs) == 16
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    import multiprocessing as mp
    with mp.get_context("fork").Pool(a.procs) as pool:
        done = pool.map(run_cell, [(c, M.REPS, str(out)) for c in cs], chunksize=1)
    rec = dict(cells=done, wall_seconds=round(time.time() - t0, 1), procs=a.procs, reps=M.REPS,
               code_sha256={f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(here.glob("*.py"))},
               plan_sha256=hashlib.sha256((here / "PLAN_M6.md").read_bytes()).hexdigest(), numpy=np.__version__)
    (out / "m6_run_record.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")
    print(json.dumps({k: rec[k] for k in ("wall_seconds", "procs", "reps")}))


if __name__ == "__main__":
    main()
