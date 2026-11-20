"""Supplementary Table S12 (S-M9; Table S17 and S-M10 until v16z5, 2026-11-17): printed from the stored result files of the exploratory analyses after the verdict; nothing is computed
beyond minima, maxima and formatting. Used by check_env_numbers.py, which compares every generated line with the supplement.
Run: PYTHONDONTWRITEBYTECODE=1 <repo>/.venv/bin/python -I -B make_S17.py [results folder]   (default decomp/results/rules_env)
"""
import json
import pathlib
import re
import sys

TORT = {"T-in": "true fraction", "T-app": "apparent fraction", "T-none": "fixed"}   # 2026-11-20 v17j: printed names of the stored generator keys

sys.dont_write_bytecode = True
Q = [("Q1", "TE only"), ("Q4", "TR only"), ("Q3", "Direction scheme only"), ("Q2", "TE, TR and direction scheme")]
SC = [("raw_ndi", "raw NDI"), ("x", "x"), ("raw_fiso", "raw FISO"), ("logit_fiso", "logit FISO"), ("y", "y"), ("odi", "ODI")]


def table_lines(J: dict) -> list:
    s1a, e1, e3, e4 = J["s1a"], J["e1"], J["e3"], J["e4"]
    L = ["**(a) Fitter output with reference tissues (210 tissues, 17 relaxation settings).**", "",
         "| Pair class | Tortuosity | Offset-subtracted share | Common shift (logit) | Proportional bias | Lowest rank correlation | Rank reversals |",
         "|---|---|---|---|---|---|---|"]
    for q, name in Q:
        for g in ("T-in", "T-app", "T-none"):
            r = [v for k, v in s1a["settings"].items() if k.startswith(f"{g}|{q}|")]
            o, off, pb = [x["ors"] for x in r], [x["offset"] for x in r], [x["prop_bias"] for x in r]
            sp, rf = [x["spearman"] for x in r], [x["rank_flip"] for x in r]
            L.append(f"| {name} | {TORT[g]} | {min(o):.3f} | {min(off):+.3f} to {max(off):+.3f} | {min(pb):+.4f} to {max(pb):+.4f} | {min(sp):.4f} | {max(rf):.0f} |")
    L += ["", "**(b) TE heterogeneity of the multi-TE series.**", "",
          "| Coordinate | Participant | Repeat SD | Between-unit SD of TE slope (per ms) | β (per ms) and 95 % interval | ΔTE* (ms) and 90 % interval | Share at 62 vs 132 ms, observed / predicted |",
          "|---|---|---|---|---|---|---|"]
    for co in ("x", "y"):
        for k, v in e1["by"].items():
            if not k.startswith(co + "|"):
                continue
            L.append(f"| {co} | {k.split('|')[1]} | {v['s_rt_corrected']:.4f} | {v['tau']:.5f} | {v['beta']:+.5f} ({v['beta_ci95'][0]:+.5f} to {v['beta_ci95'][1]:+.5f}) | "
                     f"{v['dte_star']:.1f} ({v['dte_ci90'][0]:.1f}–{v['dte_ci90'][1]:.1f}) | {v['ors_62_132_observed']:.3f} / {v['ors_62_132_predicted']:.3f} |")
        p = J["e1_fix"]["pooled"][co]["primary"] if "e1_fix" in J else e1["pooled"][co]   # 2026-11-15: pooled values with the error fixed (ANALYSIS_PLAN §8)
        L.append(f"| {co} | pooled | — | {p['tau']:.5f} | — | {p['dte_star']:.1f} ({p['dte_ci90'][0]:.1f}–{p['dte_ci90'][1]:.1f}) | — |")
    w = e4["wong"]["TE62-TE132"]
    L += ["", "**(c) Share of units at retest level after one common shift is subtracted, by scale.**", "",
          "| Scale | Multi-TE, 62 vs 132 ms | ON-Harmony, median of 19 pairs | Ten centers, without / with the common shift |", "|---|---|---|---|"]
    for k, nm in SC:
        L.append(f"| {nm} | {w[k]:.3f} | {e4['onh_median'][k]:.3f} | {e3['by'][k]['raw']:.3f} / {e3['by'][k]['ors']:.3f} |")
    if "e5" in J and "e6" in J and "s3" in J:
        L += table_lines_r2(J)
    L = [re.sub(r"[-−](0\.0+)(?![\d])", r"\1", ln) if ln.startswith("| ") else ln for ln in L]   # 2026-11-15 (F9): no sign on a value that rounds to zero
    return [re.sub(r"(?<![\w.])-(?=\d)", "−", ln) if ln.startswith("| ") else ln for ln in L]   # minus sign as printed


def table_lines_r2(J: dict) -> list:
    """(d) Analyses after the scale comparison (S-M9.5): E5 by participant and expected shares, E6 by generator, S3 grid."""
    e5, e6, s3 = J["e5"], J["e6"], J["s3"]
    L = ["", "**(d) Analyses after the scale comparison (S-M9.5).**", "",
         "| Slope structure, participant | Mean TE slope of x (per ms) | β (per ms) | R = β / mean slope | 2ν̄ − 1 | Raw-scale rate R − (2ν̄ − 1) |", "|---|---|---|---|---|---|"]
    for sb, v in e5["by"].items():
        f = v["first_order"]
        L.append(f"| {sb} | {v['kbar']:.5f} | {v['beta']:+.5f} | {f['R']:.2f} | {f['T']:.2f} | {f['raw_rate']:.2f} |")
    L += ["", "| Slope structure, expected at 62 vs 132 ms | Share, x | Share, raw NDI | Difference x − raw NDI | 2.5–97.5 % of replicates |", "|---|---|---|---|---|"]
    for k, nm in (("predicted", "Slopes depend on NDI (β as estimated)"), ("predicted_beta0", "Slopes independent of NDI (β = 0)")):
        v = e5[k]
        L.append(f"| {nm} | {v['ors_x']:.3f} | {v['ors_raw']:.3f} | {v['diff']:+.3f} | {v['diff_q025']:+.3f} to {v['diff_q975']:+.3f} |")
    L += ["", "| Shared relaxation, tortuosity | Share, x (central relaxation) | Share, raw NDI (central relaxation) | Share, x (17 settings) | Share, raw NDI (17 settings) | Settings with x ≥ 0.90 > raw NDI |",
          "|---|---|---|---|---|---|"]
    for g in ("T-in", "T-app", "T-none"):
        c = e6["by"][f"{g}|L0|iii"]
        r = [v for k, v in e6["by"].items() if k.startswith(g + "|")]
        xs, rs = [v["ors_x"] for v in r], [v["ors_raw"] for v in r]
        L.append(f"| {TORT[g]} | {c['ors_x']:.3f} | {c['ors_raw']:.3f} | {min(xs):.3f}–{max(xs):.3f} | {min(rs):.3f}–{max(rs):.3f} | {sum(v['reading'] == '가' for v in r)} of {len(r)} |")
    L += ["", "| Raw-scale term, control NDI | Group difference in true NDI | Raw-scale change, TE 80 vs 70 ms | TE 90 vs 70 ms | TE 100 vs 70 ms |", "|---|---|---|---|---|"]
    for nc in (0.46, 0.60):
        for dn in (-0.03, -0.10, -0.20):
            v = [s3["rows"][f"{nc}|{dn}|{ta}"]["raw"] for ta in (80, 90, 100)]
            L.append(f"| {nc:.2f} | {dn:+.2f} | " + " | ".join(f"{x:+.4f}" for x in v) + " |")
    return L


if __name__ == "__main__" and not globals().get("RES"):
    RES = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else pathlib.Path(__file__).resolve().parents[3] / "decomp" / "results" / "rules_env"
    J = {n: json.loads((RES / f"{n}.json").read_text(encoding="utf-8")) for n in ("s1a", "e1", "e1_fix", "e3", "e4", "e5", "e6", "s3") if (RES / f"{n}.json").exists()}
    print("\n".join(table_lines(J)))
