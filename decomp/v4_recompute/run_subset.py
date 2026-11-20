#!/usr/bin/env python
"""Recompute part of the V4 calibration test with the unchanged judge code. Criteria: CRITERIA.md in this folder.

    python decomp/v4_recompute/run_subset.py --pairs OXF1PRI~OXF2PRI,OXF1PRI~NOT2ING-89 --out DIR
    python decomp/v4_recompute/run_subset.py --main-condition --out DIR

Only the main condition (SNR 30, N_eff 50) is run. Results go to DIR/v4_recompute.json, in the structure of v4_c3.json without the
sensitivity fields; each finished pair is kept in DIR/units, so a stopped run continues where it stopped. Nothing under data/ is written.
The two options --seed-offset and --shift-roi exist only for the failure controls.
Every input is entered in a read log before it is read (read_log.py); if the log cannot be written, the run stops before reading anything.
"""
from __future__ import annotations

import argparse
import os
import pathlib
import pickle
import sys
import time

REPO = pathlib.Path(__file__).resolve().parents[2]
DATA = REPO / "data"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", default="", help="comma-separated pairs, for example OXF1PRI~OXF2PRI")
    ap.add_argument("--main-condition", action="store_true", help="all pairs of the main condition")
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed-offset", type=int, default=0, help="failure control: add this to each torch seed")
    ap.add_argument("--shift-roi", type=float, default=0.0, help="failure control: add this to one ROI value of the first row of each pair")
    ap.add_argument("--dry-run", action="store_true", help="load everything and stop before the first draw")
    ap.add_argument("--read-log", default="", help="read log, default DIR/read_log.jsonl of --out")
    ap.add_argument("--read-what", default="V4 recomputation", help="what the read log says this run is")
    args = ap.parse_args()
    if bool(args.pairs) == args.main_condition:
        raise SystemExit("Give either --pairs or --main-condition")
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
    from read_log import ReadLog  # noqa: E402

    log = ReadLog(args.read_log, args.read_what) if args.read_log else None
    if log is None:
        pathlib.Path(args.out).mkdir(parents=True, exist_ok=True)
        log = ReadLog(pathlib.Path(args.out) / "read_log.jsonl", args.read_what)

    os.environ.setdefault("DECOMP_WORK", str(DATA))
    os.environ.setdefault("REMEASURE_WORK", str(DATA / "remeasure"))
    os.environ["CUDA_VISIBLE_DEVICES"] = ""
    os.environ.setdefault("S4_THREADS", "16")
    sys.path[:0] = [str(REPO / "decomp" / d) for d in ("stage1", "s4", "remeasure")]
    import config as C  # noqa: E402

    C.PROTO_DIR = DATA / "protocols"                                # the same run-time location as reproduce/judge_v4_run.py
    import json  # noqa: E402

    import numpy as np  # noqa: E402
    import guard as G  # noqa: E402
    import judge_v4 as J  # noqa: E402
    import rconfig as R  # noqa: E402
    import s4cfg as S  # noqa: E402

    arm, seed = "c3", 20261027
    out = pathlib.Path(args.out)
    (out / "units").mkdir(parents=True, exist_ok=True)
    basis = G.require()
    at_path = R.WORK / "onharmony" / J.V4_DIR / "attr_table.json"
    log.before_read([at_path])
    at = json.loads(at_path.read_text(encoding="utf-8"))
    tables = sorted({R.WORK / "onharmony" / f"{sub}_{s}" / "roi.json" for d in at["pairs"].values() for sub, a in d.items() for s in a["sessions"]})
    log.before_read([*sorted((S.WORK / "models").glob("*.obj.pt")), S.WORK / "sim_test.npz", S.WORK / "plan_test.npz", S.WORK / "v4" / f"test_cover_{arm}.npz",
                     J.RESULTS / f"test_{arm}.json", J.RESULTS / "miss.json", J.RESULTS / f"cards_{arm}.json", *tables])
    draw_fn, logq_fn, K = J.real_fns(arm)
    miss = json.loads((J.RESULTS / "miss.json").read_text(encoding="utf-8"))
    cards = json.loads((J.RESULTS / f"cards_{arm}.json").read_text(encoding="utf-8"))
    test_npz = S.WORK / "v4" / f"test_cover_{arm}.npz"
    tc = dict(np.load(test_npz))
    drawn = J.drawn_thresholds(tc, arm, seed)
    rows = J.collect(R.WORK, at["pairs"])
    all_pairs = J.card_pairs()
    chosen = all_pairs if args.main_condition else args.pairs.split(",")
    unknown = sorted(set(chosen) - set(all_pairs))
    if unknown:
        raise SystemExit(f"Unknown pairs: {unknown}")
    ident = {c: [p for p in pairs if p in chosen] for c, pairs in J.identified_cells(cards).items()}
    if args.dry_run:
        print("dry run", dict(guard=basis["fixed_commit"][:12], pairs=len(chosen), cells={c: len(p) for c, p in ident.items()},
                              rows={p: len(rows[p]) for p in chosen[:3]}))
        return 0

    _, snr, neff = J.MAIN_CELL
    seeds, seconds, units = {}, {}, {}
    for pi, pair in enumerate(all_pairs):
        if pair not in chosen:
            continue
        keep = out / "units" / f"{pair}.pkl"
        if keep.exists():
            saved = pickle.loads(keep.read_bytes())
            seeds[saved["key"]], seconds[pair], units[pair] = saved["seed"], saved["seconds"], saved["units"]
            print("kept", pair, flush=True)
            continue
        row = rows[pair]
        if args.shift_roi and row:
            row[0] = dict(row[0], x=np.array(row[0]["x"], float))
            row[0]["x"][2] += args.shift_roi
        key = f"SNR{int(snr)}|Neff{int(neff)}|{pair}"
        started = time.time()
        seeds[key] = J.seed_torch(6000 + 1000 * J.ARMS.index(arm) + pi + args.seed_offset)
        units[pair] = J.unit_results(pair, row, arm, draw_fn, logq_fn, K, miss["layers"][f"{J.acq_type(pair)}|L0"]["threshold"], snr, neff)
        seconds[pair] = round(time.time() - started, 1)
        keep.write_bytes(pickle.dumps(dict(key=key, seed=seeds[key], seconds=seconds[pair], units=units[pair])))
        print(pair, "seed", seeds[key], "units", len(row), "seconds", seconds[pair], flush=True)

    res = dict(arm=arm, seeds_torch=seeds, seed_numpy=seed, inputs=dict(
        attr_table_sha256=J.sha256_file(at_path), test_cover_sha256=J.sha256_file(test_npz),
        design_files_sha256={f"{S.DESIGNS[d]}_{e}": J.sha256_file(S.C.PROTO_DIR / f"{S.DESIGNS[d]}_{e}") for d in S.DESIGNS for e in ("bval", "bvec")},
        frozen=J.frozen_expected()))
    variants = (("main", lambda u: ~u["main_excluded"]), ("with_16794_NOT4GEP001", lambda u: np.ones(len(u["subjects"]), bool)),
                ("without_snr_out_of_train", lambda u: ~u["main_excluded"] & ~u["noise"]))
    for label, keep_fn in variants:                                 # as in judge_v4.run, without the sensitivity conditions
        rng = np.random.default_rng(seed)
        main = {p: J.select(u, keep_fn(u)) for p, u in units.items()}
        thr, thr_u, p95 = {}, {}, {}
        for pair, u in main.items():
            n_cl, n_un = len(np.unique(u["subjects"])), len(u["subjects"])
            t = J.acq_type(pair)
            p95[pair] = drawn("p95", t, n_cl)
            for ci, c in enumerate(J.COMPS):
                thr[(pair, c)] = drawn("thr", t, n_cl, ci)
                thr_u[(pair, c)] = drawn("thr", t, n_un, ci)
        attrs = {p: dict(d) for p, d in at["pairs"].items()} if label == "main" else None
        res[label] = J.judge_cells(ident, main, thr, p95, rng, None, attrs, thr_u)
    res["recompute"] = dict(mode="main condition" if args.main_condition else "pairs", pairs=chosen, seconds_by_pair=seconds, total_seconds=round(sum(seconds.values()), 1),
                            threads=os.environ["S4_THREADS"], seed_offset=args.seed_offset, shift_roi=args.shift_roi, guard=basis)
    print("output sha256", J._write_json(out / "v4_recompute.json", res), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
