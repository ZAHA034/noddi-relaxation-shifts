"""Figure 4 of manuscript v16z9 → v17d (Results 4.2): the multi-site design simulation (post hoc; assumed values; Supplementary Methods S-M11).
Drawn from the stored cell metrics only; nothing is simulated, fitted or tested here.

(a) Coverage of the fraction-only difference by 95 % intervals where the groups differ only in compartmental T2 (S2b: both tissue
    compartments; S2e: extra-neurite only), design D1 (the eight ON-Harmony acquisitions), balanced sites, no between-person relaxation
    variation, no unrelated site x group term: M1 standard ComBat, M2 fixed site term, M3 group x TE term at zero TE, M5b measured
    relaxation with 10 % error; v17d adds M6a, a second TE 30 ms longer in every participant (decomp/design_sim/PLAN_M6.md; same cells,
    same first observation). Line at 0.95. v17j (2026-11-20): the legend gives the analysis names without the internal labels M1–M6a.
(b) Power against fraction-only differences of the same apparent sizes (S1b, S1e), M3, M5b and M6a, same cells.
(c) Share of studies in which the group x TE term was significant with a vendor-linked group effect and no relaxation difference (H2, S0):
    D1 and D2 (short and long TEs within each vendor; read after the results of (a) and (b)). Line at 0.05.
Guards (explicit failures, not asserts): 288 cells; every drawn value equals the cell metric it names; the D1 values equal the stored readings
of summary.json (R-a, R-b, R-d) to 1e-12; the M6a values equal the 16 M6 cell metrics; a planted fault (--plant) that swaps M1 and M5b must
fail the summary guard.
Reads (log each file with drafts/x15/ledger/reread_log.py before running; done 2026-11-18 03:39 KST for the 288 cells): <results>/cells/*.json,
<results>/summary.json, <m6 results>/m6_cells/*.json (v17d). Writes FIG_OUT/F8_design_sim.pdf/.png (figstyle audit must pass).
Run (CPU only):
  PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=<scratch> FIG_OUT=<out> <repo>/.venv/bin/python -I -B make_F8_design_sim.py <results folder> <m6 results folder> [--plant]
"""
import sys

sys.dont_write_bytecode = True

import glob  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
import pathlib  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import figcolors as fc  # noqa: E402
import figstyle as fs  # noqa: E402

RES = pathlib.Path(sys.argv[1])
RES6 = pathlib.Path(sys.argv[2])
PLANT = "--plant" in sys.argv
cells = {}
for p in sorted(glob.glob(str(RES / "cells" / "cell_*.json"))):
    d = json.loads(pathlib.Path(p).read_text(encoding="utf-8"))
    c = d["cell"]
    cells[(c["D"], c["S"], float(c["V"]), c["H"], c["N"])] = d["metrics"]
if len(cells) != 288:
    raise SystemExit(f"cells: {len(cells)} != 288")
summ = json.loads((RES / "summary.json").read_text(encoding="utf-8"))["readings"]
m6c = {}
for p in sorted(glob.glob(str(RES6 / "m6_cells" / "cell_*.json"))):
    d = json.loads(pathlib.Path(p).read_text(encoding="utf-8"))
    c = d["cell"]
    m6c[(c["S"], float(c["V"]))] = d["metrics"]["M6a"]
if len(m6c) != 16:
    raise SystemExit(f"M6 cells: {len(m6c)} != 16")


def m(D, S, k, f, H="H0"):
    return cells[(D, S, 0.0, H, "bal")][k][f]


AN = [("M1", "M1", "ComBat", fc.OI["verm"], "o"), ("M2", "M2", "Site term", fc.OI["orange"], "s"),
      ("M3_te0", "M3", "Group × TE", fc.OI["blue"], "^"), ("M5b", "M5b", "Measured relaxation", fc.OI["green"], "D")]
if PLANT:
    AN[0], AN[3] = (AN[3][0],) + AN[0][1:], (AN[0][0],) + AN[3][1:]
cov = {(S, lab): m("D1", S, key, "coverage") for S in ("S2b", "S2e") for key, lab, *_ in AN}
pw = {(S, lab): m("D1", S, key, "reject") for S in ("S1b", "S1e") for key, lab, *_ in AN if lab in ("M3", "M5b")}
fp = {D: m(D, "S0", "M3_slope", "reject", "H2") for D in ("D1", "D2")}
for S in ("S2b", "S2e"):
    cov[(S, "M6a")] = m6c[(S, 0.0)]["coverage"]
for S in ("S1b", "S1e"):
    pw[(S, "M6a")] = m6c[(S, 0.0)]["reject"]
bad = []
for S, rk in (("S2b", "R-a S2b"), ("S2e", "R-a S2e")):
    if abs(cov[(S, "M1")] - summ[rk]["M1_coverage"]) > 1e-12:
        bad.append(f"M1 coverage {S} {cov[(S, 'M1')]} != summary {summ[rk]['M1_coverage']}")
for S, rk in (("S2b", "R-b b"), ("S2e", "R-b e")):
    if abs(cov[(S, "M3")] - summ[rk]["M3_te0_coverage_S2"]) > 0.01:      # summary pools the S2 cells of the reading; cell value within 0.01
        bad.append(f"M3 coverage {S} {cov[(S, 'M3')]} vs summary {summ[rk]['M3_te0_coverage_S2']}")
for S, rk in (("S1b", "R-b b"), ("S1e", "R-b e")):
    if abs(pw[(S, "M3")] - summ[rk]["M3_te0_power_S1"]) > 1e-12:
        bad.append(f"M3 power {S} {pw[(S, 'M3')]} != summary {summ[rk]['M3_te0_power_S1']}")
if abs(fp["D1"] - summ["R-d"]["M3_slope_reject"]) > 1e-12:
    bad.append("R-d")
if bad:
    raise SystemExit("guard failed: " + "; ".join(bad))

fig = fs.new_fig(190.0, 52.0)
ys = {"S2b": 1.0, "S2e": 0.0, "S1b": 1.0, "S1e": 0.0}
ylab = {"S2b": "Both tissue\ncompartments", "S2e": "Extra-neurite\nonly", "S1b": "Both-compartment\nsize", "S1e": "Extra-neurite\nsize"}
axa = fs.ax_mm(fig, 27, 12, 52, 30)
axb = fs.ax_mm(fig, 103, 12, 36, 30)
axc = fs.ax_mm(fig, 156, 12, 28, 30)
off = {"M1": 0.18, "M2": 0.06, "M3": -0.06, "M5b": -0.18, "M6a": -0.29}
for key, lab, name, col, mk in AN:
    for S in ("S2b", "S2e"):
        axa.plot([cov[(S, lab)]], [ys[S] + off[lab]], ls="none", marker=mk, ms=4.2, mfc=col, mec="white", mew=0.5, color=col,
                 label=name if S == "S2b" else None, zorder=3)   # v17j: legend shows the analysis name only
    if lab in ("M3", "M5b"):
        for S in ("S1b", "S1e"):
            axb.plot([pw[(S, lab)]], [ys[S] + off[lab] * 1.5], ls="none", marker=mk, ms=4.2, mfc=col, mec="white", mew=0.5, zorder=3)
M6 = ("M6a", "Second TE", fc.OI["sky"], "v")
for S in ("S2b", "S2e"):
    axa.plot([cov[(S, "M6a")]], [ys[S] + off["M6a"]], ls="none", marker=M6[3], ms=4.2, mfc=M6[2], mec="white", mew=0.5, color=M6[2],
             label="Second TE" if S == "S2b" else None, zorder=3)
for S in ("S1b", "S1e"):
    axb.plot([pw[(S, "M6a")]], [ys[S] + off["M6a"] * 1.5], ls="none", marker=M6[3], ms=4.2, mfc=M6[2], mec="white", mew=0.5, zorder=3)
for ax, S2, xl in ((axa, ("S2b", "S2e"), "Coverage of the fraction-only difference"), (axb, ("S1b", "S1e"), "Power")):
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.5, 1.5)
    ax.set_yticks([ys[s] for s in S2])
    ax.set_yticklabels([ylab[s] for s in S2], fontsize=fs.PT_TEXT)
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_xlabel(xl, fontsize=fs.PT_TEXT)
    fs.hgrid(ax, "x")
axa.axvline(0.95, color=fc.INK2, lw=fs.LW_THIN, ls=(0, (3, 2)), zorder=1)
for i, D in enumerate(("D1", "D2")):
    axc.plot([fp[D]], [1 - i], ls="none", marker="^", ms=4.2, mfc=fc.OI["blue"], mec="white", mew=0.5, zorder=3)
axc.axvline(0.05, color=fc.INK2, lw=fs.LW_THIN, ls=(0, (3, 2)), zorder=1)
axc.set_xlim(0, 0.12)
axc.set_xticks([0, 0.05, 0.10])
axc.set_xticklabels(["0", "0.05", "0.10"])
axc.set_ylim(-0.5, 1.5)
axc.set_yticks([1, 0])
axc.set_yticklabels(["ON-Harmony TEs", "TEs within vendor"], fontsize=fs.PT_TEXT)
axc.set_xlabel("False-positive group × TE", fontsize=fs.PT_TEXT)
fs.hgrid(axc, "x")
for ax in (axa, axb, axc):
    ax.tick_params(labelsize=fs.PT_TEXT, width=fs.LW_AXIS, length=2.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_linewidth(fs.LW_AXIS)
fs.panel(axa, "a", dx_mm=-25.0)
fs.panel(axb, "b", dx_mm=-25.0)
fs.panel(axc, "c", dx_mm=-26.0)
leg = fig.legend(loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=5, frameon=False, fontsize=fs.PT_TEXT, handletextpad=0.3,
                 columnspacing=1.6)
probs = fs.audit(fig, "F8_design_sim")
if probs:
    raise SystemExit(f"audit: {probs}")
out = pathlib.Path(os.environ["FIG_OUT"])
out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / "F8_design_sim.pdf", metadata={"CreationDate": None, "ModDate": None})
fig.savefig(out / "F8_design_sim.png", dpi=300)
vals = dict(coverage=cov and {f"{k[0]}|{k[1]}": v for k, v in cov.items()}, power={f"{k[0]}|{k[1]}": v for k, v in pw.items()}, false_pos=fp)
(out / "F8_design_sim_values.json").write_text(json.dumps(vals, indent=1), encoding="utf-8")
print("wrote", out / "F8_design_sim.pdf", json.dumps(vals))
