"""Stop gate that recomputes the existing reported values with the new loader (manuscript 4.1, 4.2, 4.3): exits with code 1 if any value differs and computes no new quantity."""
from __future__ import annotations

import json
import pathlib
import sys

from . import empirical as EM

EXPECTED = dict(tong=dict(x=0.931, ODI=0.934, y=0.844),            # N-324 · N-325 · N-326
                wong=dict(x=0.815, y=0.901),                       # N-247 · N-248
                onh_odi=dict(n_retest_level=2, n_pairs=19))         # N-018 · N-019


def run(work: pathlib.Path = EM.R.WORK) -> dict:
    t = EM.read_tong(work)
    v, rep = EM.tong_structs(t)
    tb = EM.JT.judge(v, rep)["T_B"]
    w = EM.read_wong(work)
    jw = EM.JW.judge(EM.wong_series(w))["overall"]
    st = EM.onh_states(work)
    o = EM.read_onh(work, st)
    s_main, df, _ = EM.JO.s_ref_f2(work, st["states"], o["rois"], o["units"])
    os_ = EM.JO.pair_values(EM.JO.device_values(st["main_by_unit"], o["rois"], "ODI", "raw", "주"), o["units"])
    keep_o, _ = EM.JO.keep_pairs(o["cnt"], os_)
    ob = EM.JO.odi_b(os_, keep_o, s_main, o["rows_by_pair"])
    got = dict(tong={co: tb[co]["ratio"] for co in ("x", "ODI", "y")}, wong={co: jw[co]["linear_ratio"] for co in ("x", "y")},
               onh_odi=dict(n_retest_level=ob["n_retest_level"], n_pairs=ob["n_pairs"]))
    ok = (all(round(got["tong"][k], 3) == e for k, e in EXPECTED["tong"].items())
          and all(round(got["wong"][k], 3) == e for k, e in EXPECTED["wong"].items())
          and got["onh_odi"] == EXPECTED["onh_odi"])
    return dict(ok=ok, got=got, expected=EXPECTED, missing=dict(tong=t["missing"], wong=w["missing"]),
                n_sessions=dict(tong=len(t["roi"]), wong=sum(len(s) for s in w["roi"].values()), onh=len(o["rois"])), s_ref_f2_df=df)


if __name__ == "__main__":
    res = run()
    pathlib.Path(sys.argv[1]).write_text(json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")
    print(json.dumps(res, ensure_ascii=False))
    sys.exit(0 if res["ok"] else 1)
