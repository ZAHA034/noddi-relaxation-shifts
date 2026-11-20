"""Run the Tong judgment through its functions.

The command-line entry of decomp/remeasure/judge_tong.py stops with an error, a defect recorded
as deviation (iii) in Supplementary Table S1: it unpacks the loader's return value into two names,
and the loader returns three. The locked file is not changed. This script does what the original
run did (records/x12_tong_driver.py): it calls the same functions and takes three return values.
"""
import json
import pathlib
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "decomp" / "remeasure"))
import judge_tong as JT  # noqa: E402
import rconfig as R  # noqa: E402

work = R.WORK
v, rp, miss = JT.load(work)
res = dict(main=JT.judge(v, rp), without_c09=JT.judge(v, rp, tuple(c for c in JT.CENTERS if c != "c09")), missing_sessions=miss)
vs, rps, _ = JT.load(work, fname="roi_allshell.json")
res["all_shells_descriptive"] = JT.judge(vs, rps)
print(json.dumps(res, ensure_ascii=False, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))
