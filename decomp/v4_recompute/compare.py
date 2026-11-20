#!/usr/bin/env python
"""Compare a recomputed V4 output with the deposited one. Criteria: CRITERIA.md in this folder.

    python decomp/v4_recompute/compare.py --new DIR/v4_recompute.json --report DIR/compare.json
    python decomp/v4_recompute/compare.py --selftest

Every file is entered in a read log before it is read (read_log.py). If the log cannot be written, the comparison stops before reading anything.

Exit code: 0 = Level A in everything compared, 1 = Level B, 2 = different.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import pathlib
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from read_log import ReadLog  # noqa: E402

REPO = pathlib.Path(__file__).resolve().parents[2]
TARGET = REPO / "data" / "remeasure" / "onharmony" / "v4" / "v4_c3.json"
VARIANTS = ("main", "with_16794_NOT4GEP001", "without_snr_out_of_train")
BLOCKS = ("summary", "info_none", "device_concentration_descriptive")
NOT_RUN = {"sensitivity_passed", "snr_neff_dependent"}             # need the sensitivity conditions
MISSING = "<field absent>"
LOG = None


def canon(v) -> str:
    return json.dumps(v, sort_keys=True, ensure_ascii=False)


def field_diff(new: dict, old: dict, skip: set[str]) -> list[str]:
    return sorted(k for k in (set(new) | set(old)) - skip if canon(new.get(k, MISSING)) != canon(old.get(k, MISSING)))


def level_b(new: dict, old: dict) -> bool:
    """Level B of CRITERIA.md, for a cell with data."""
    if "coverage" not in old:
        return canon(new) == canon(old)
    if any(canon(new.get(k, MISSING)) != canon(old.get(k, MISSING)) for k in ("verdict", "passed", "threshold", "n_units", "n_subjects")):
        return False
    a, b = new.get("coverage"), old.get("coverage")
    if a is None or b is None:
        return a is None and b is None
    return abs(a - b) <= max(0.02, 1 / max(old.get("n_units_used") or 1, 1))


def compare_outputs(new: dict, old: dict) -> dict:
    full = new.get("recompute", {}).get("mode") == "main condition"
    skip = set(NOT_RUN) if full else NOT_RUN | {"boot95"}
    out = dict(mode=new.get("recompute", {}).get("mode"), inputs_equal=canon(new.get("inputs")) == canon(old.get("inputs")), cells=[], blocks={})
    for var in VARIANTS if full else ("main",):
        cells_new, cells_old = new[var]["cells"], old[var]["cells"]
        names = sorted(set(cells_new) | set(cells_old)) if full else sorted(cells_new)
        for name in names:
            cn, co = cells_new.get(name, {}), cells_old.get(name, {})
            diff = field_diff(cn, co, skip)
            out["cells"].append(dict(variant=var, cell=name, level="A" if not diff else "B" if cn and co and level_b(cn, co) else "differs", differing_fields=diff,
                                     verdict=[cn.get("verdict"), co.get("verdict")], coverage=[cn.get("coverage"), co.get("coverage")],
                                     threshold=[cn.get("threshold"), co.get("threshold")], n_units_used=[cn.get("n_units_used"), co.get("n_units_used")]))
        if full:
            for b in BLOCKS:
                out["blocks"][f"{var}.{b}"] = canon(new[var].get(b)) == canon(old[var].get(b))
    levels = {c["level"] for c in out["cells"]}
    exact = out["inputs_equal"] and levels <= {"A"} and all(out["blocks"].values())
    out["level"] = "A" if exact else "differs" if "differs" in levels or not out["inputs_equal"] or not all(out["blocks"].values()) else "B"   # the blocks follow from the verdicts
    out["n_cells"], out["n_cells_level_a"] = len(out["cells"]), sum(c["level"] == "A" for c in out["cells"])
    return out


def check_target() -> dict:
    """The SHA-256 of the target must equal the values recorded in the manifest and in the verdict record."""
    LOG.before_read([TARGET, REPO / "data" / "MANIFEST_data.csv", REPO / "decomp" / "results" / "remeasure" / "judgements.json"])
    sha = hashlib.sha256(TARGET.read_bytes()).hexdigest()
    rel = str(TARGET.relative_to(REPO))
    with open(REPO / "data" / "MANIFEST_data.csv", encoding="utf-8") as fh:
        rows = list(csv.reader(fh))
    in_manifest = any(sha in r and any(rel.endswith(c) or c.endswith("v4_c3.json") for c in r) for r in rows)
    verdict = json.loads((REPO / "decomp" / "results" / "remeasure" / "judgements.json").read_text(encoding="utf-8"))
    recorded = [v["sha256"] for k, v in verdict.get("입력", {}).items() if k == "v4_c3"]
    if not in_manifest or recorded != [sha]:
        raise SystemExit(f"The target does not match its recorded SHA-256 (manifest {in_manifest}, verdict record {recorded == [sha]})")
    return dict(target_sha256=sha, in_manifest=in_manifest, in_verdict_record=True)


def ledger_rows(new: dict, out_dir: pathlib.Path) -> dict:
    """Ledger rows whose source is a V4 field of the verdict record, evaluated on the new values."""
    sys.path.insert(0, str(REPO / "reproduce"))
    import ledger_fields as LF  # noqa: E402

    view = out_dir / "ledger_view.json"
    view.write_text(json.dumps({"V4": new["main"]}, ensure_ascii=False), encoding="utf-8")
    same, differ, not_computed = [], [], []
    with open(REPO / "ledger" / "selectors.csv", encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh) if r["source_file"].endswith("judgements.json") and "V4." in r["source_field"]]
    for r in rows:
        if any(k in r["source_field"] for k in NOT_RUN):
            not_computed.append(r["id"])
            continue
        try:
            value = LF.eval_field(out_dir, {"근거 파일": view.name, "근거 필드": r["source_field"], "근거 파일 sha256": "", "원 값": r["raw_value"]})
            (same if LF.same_value(value, r["raw_value"]) else differ).append(r["id"])
        except KeyError:
            not_computed.append(r["id"])
    return dict(rows=len(rows), same=same, differ=differ, not_computed=not_computed)


def selftest() -> int:
    old = json.loads(TARGET.read_text(encoding="utf-8"))
    def case(change):
        new = copy.deepcopy(old)
        new["recompute"] = dict(mode="main condition")
        change(new)
        return compare_outputs(new, old)["level"]
    name = "OXF1PRI~NOT2ING-89|ODI"
    def nudge(n):
        n["main"]["cells"][name]["coverage"] += 1e-9
    def flip(n):
        n["main"]["cells"][name]["verdict"] = "pass" if n["main"]["cells"][name]["verdict"] != "pass" else "fail"
    def far(n):
        n["main"]["cells"][name]["coverage"] += 0.5
    def drop(n):
        del n["main"]["cells"][name]
    def thr(n):
        n["main"]["cells"][name]["threshold"] += 0.01
    def inp(n):
        n["inputs"]["attr_table_sha256"] = "0" * 64
    def sens(n):
        n["main"]["cells"][name].pop("sensitivity_passed")
    def block(n):
        n["main"]["summary"]["ODI"]["failed"] += 1
    expect = [("unchanged", lambda n: None, "A"), ("sensitivity fields absent", sens, "A"), ("coverage + 1e-9", nudge, "B"), ("coverage + 0.5", far, "differs"), ("verdict changed", flip, "differs"),
              ("threshold changed", thr, "differs"), ("cell removed", drop, "differs"), ("inputs changed", inp, "differs"), ("summary changed", block, "differs")]
    bad = 0
    for label, fn, want in expect:
        got = case(fn)
        bad += got != want
        print(f"{'ok ' if got == want else 'BAD'} {label}: {got} (expected {want})")
    return 1 if bad else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--new")
    ap.add_argument("--report")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--read-log", default="", help="read log, default read_log.jsonl beside --new, or a temporary file for --selftest")
    ap.add_argument("--read-what", default="V4 comparison", help="what the read log says this run is")
    args = ap.parse_args()
    global LOG
    LOG = ReadLog(args.read_log or (pathlib.Path(args.new).resolve().parent / "read_log.jsonl" if args.new else pathlib.Path(tempfile.mkdtemp(prefix="v4_compare_")) / "read_log.jsonl"), args.read_what)
    print(check_target())
    if args.selftest:
        return selftest()
    LOG.before_read([args.new])
    new = json.loads(pathlib.Path(args.new).read_text(encoding="utf-8"))
    old = json.loads(TARGET.read_text(encoding="utf-8"))
    res = compare_outputs(new, old)
    res["target"] = check_target()
    res["recompute"] = new.get("recompute")
    if res["mode"] == "main condition":
        res["ledger"] = ledger_rows(new, pathlib.Path(args.new).resolve().parent)
    for c in res["cells"]:
        print(f"{c['level']:8} {c['variant'][:12]:12} {c['cell']:32} {str(c['verdict'][0]):8} {c['coverage'][0]} | {c['coverage'][1]}  {','.join(c['differing_fields'])}")
    print("LEVEL", res["level"], "· cells", res["n_cells"], "· Level A", res["n_cells_level_a"], "· inputs equal", res["inputs_equal"], "· blocks", res["blocks"] or "not compared")
    if "ledger" in res:
        print("ledger", {k: (v if k == "rows" else len(v)) for k, v in res["ledger"].items()}, "· differ", res["ledger"]["differ"], "· not computed", res["ledger"]["not_computed"])
    if args.report:
        pathlib.Path(args.report).write_text(json.dumps(res, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    return {"A": 0, "B": 1}.get(res["level"], 2)


if __name__ == "__main__":
    raise SystemExit(main())
