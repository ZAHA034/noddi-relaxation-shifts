"""design_sim.summarize — collects the per-cell results of the main run (cells/cell_###.json) into the values of readings R-a – R-d (plan Section 6), the 24 cells of the supplementary table (plan «Change record») and the table of all 288 cells (deposit); no new computation."""
from __future__ import annotations

import json
import pathlib
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from design_sim import model as M  # noqa: E402

TABLE_KEYS = ("M1", "M2", "M3_te0", "M5b")


def load(out: pathlib.Path) -> dict:
    cells = {}
    for p in sorted((out / "cells").glob("cell_*.json")):
        d = json.loads(p.read_text(encoding="utf-8"))
        c = d["cell"]
        cells[(c["D"], c["S"], float(c["V"]), c["H"], c["N"])] = d
    if len(cells) != 288:
        raise SystemExit(f"칸 {len(cells)} 개 — 288 개가 아님")
    return cells


def get(cells, D, S, V=0.0, H="H0", N="bal"):
    return cells[(D, S, float(V), H, N)]["metrics"]


def readings(cells) -> dict:
    out = {}
    for x in ("e", "b"):
        m = get(cells, "D1", f"S2{x}")
        cov = m["M1"]["coverage"]
        out[f"R-a S2{x}"] = dict(M1_coverage=cov, M1_coverage_mcse=m["M1"]["coverage_mcse"], M1_bias=m["M1"]["bias"],
                                 sentence=("a standard ComBat analysis in this design does not give valid intervals for the fraction-only "
                                           "difference") if cov < 0.90 else "the within-site bias was small against sampling error at this size")
        cov3 = m["M3_te0"]["coverage"]
        pw = get(cells, "D1", f"S1{x}")["M3_te0"]["reject"]
        rule = cov3 >= 0.93 and pw < 0.50
        out[f"R-b {x}"] = dict(M3_te0_coverage_S2=cov3, M3_te0_power_S1=pw,
                               sentence="extrapolation to TE → 0 is valid but imprecise in this design" if rule else None)
    for x in ("b", "e"):
        d1, d2 = get(cells, "D1", f"S1{x}"), get(cells, "D2", f"S1{x}")
        out[f"R-c S1{x}"] = dict(D1_se=d1["M3_te0"]["mean_se"], D2_se=d2["M3_te0"]["mean_se"], D1_power=d1["M3_te0"]["reject"],
                                 D2_power=d2["M3_te0"]["reject"])
    m = get(cells, "D1", "S0", H="H2")
    out["R-d"] = dict(M3_slope_reject=m["M3_slope"]["reject"], M3_slope_reject_mcse=m["M3_slope"]["reject_mcse"])
    return out


def table24(cells) -> list:
    rows = []
    for D in M.DESIGNS:
        for S in M.SCENARIOS:
            m = get(cells, D, S)
            row = dict(D=D, S=S)
            for k in TABLE_KEYS:
                row[k] = (dict(bias=m[k]["bias"], coverage=m[k]["coverage"], reject=m[k]["reject"]) if m[k].get("defined", True)
                          else None)
            rows.append(row)
    return rows


def full_rows(cells) -> list:
    rows = []
    for key, d in sorted(cells.items(), key=lambda kv: kv[1]["cell"]["idx"]):
        m = d["metrics"]
        r = dict(idx=d["cell"]["idx"], D=key[0], S=key[1], V=key[2], H=key[3], N=key[4], theta=d["theta"])
        for k, v in m.items():
            if isinstance(v, dict) and v.get("defined", True) is False:
                r[f"{k}.defined"] = "no"
            elif isinstance(v, dict):
                for f in ("bias", "coverage", "reject", "mean_se", "empsd", "rmse", "mean_est"):
                    if f in v:
                        r[f"{k}.{f}"] = v[f]
            else:
                r[k] = v
        rows.append(r)
    return rows


def main():
    out = pathlib.Path(sys.argv[1])
    cells = load(out)
    summ = dict(readings=readings(cells), table24=table24(cells), te_c=M.TE_C)
    (out / "summary.json").write_text(json.dumps(summ, ensure_ascii=False, indent=1), encoding="utf-8")
    rows = full_rows(cells)
    cols = list(dict.fromkeys(c for r in rows for c in r))
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for r in rows:
        lines.append("| " + " | ".join(f"{r[c]:.4f}" if isinstance(r.get(c), float) else str(r.get(c, "—")) for c in cols) + " |")
    (out / "all_cells.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(summ["readings"], ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
