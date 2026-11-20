"""Supplementary Fig. S6 (v16e, S-M10): exploratory analyses of the measurement-function environment after the verdict.
Drawn from the stored rules_env results only (drafts x15/env/out/*.json, recorded before they were read); nothing is fitted or tested here
except the expected-share curves of panel (a), which use the same function as the analysis (rules_env.metrics._te_share).

(a) Multi-TE series: expected share at retest level after one common shift is subtracted, against the TE difference, from each participant's
    between-unit SD of TE slopes and repeat SD (solid), and from the fitter alone (dashed, median repeat SD); dotted line 0.90; points, observed
    shares at 62 against 132 ms (E1).
(b) Multi-TE series: offset-subtracted share for x and raw NDI in each of the 28 TE pairs against the TE difference (E4).
(c) Fitter output with reference tissues: offset-subtracted share by pair class and generator (S1a).
(d) Proposition 4 at fitter output: fitted change in the group difference between TE 70 and 100 ms divided by Eq. (7), by relaxation
    scenario, four location terms (small horizontal offsets) and three generators, population mode, logit scale; band 0.80–1.25 (S2).
Run: PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=<scratch> RULES_ENV_REPO=<repo> FIG_DRAFTS=<working folder> FIG_OUT=<folder> <repo>/.venv/bin/python -I -B make_FS6.py
"""
import sys

sys.dont_write_bytecode = True

import hashlib  # noqa: E402
import importlib.util  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
import pathlib  # noqa: E402
import subprocess  # noqa: E402

os.environ["CUDA_VISIBLE_DEVICES"] = ""


def stop(msg):
    raise SystemExit(f"STOP FS6: {msg}")


HERE = pathlib.Path(__file__).resolve().parent
# Working folder (drafts): environment variable FIG_DRAFTS; if absent, the place when this file is in drafts/x15/ms/R1_files/src
DRAFTS = pathlib.Path(os.environ.get("FIG_DRAFTS") or HERE.parents[3]).resolve()
if not (DRAFTS / "x15" / "env" / "out").is_dir():
    stop("set FIG_DRAFTS to the working folder that holds x15/env/out")
RES = DRAFTS / "x15" / "env" / "out"
LEDGER = DRAFTS / "x15" / "ledger"
sys.path.insert(0, str(DRAFTS / "x15" / "env" / "src"))
if not os.environ.get("FIG_OUT"):
    stop("set FIG_OUT")
OUT = pathlib.Path(os.environ["FIG_OUT"]).resolve()

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402

_fc = HERE / "figcolors.py" if (HERE / "figcolors.py").is_file() else DRAFTS / "x15" / "ms" / "S65_files" / "src" / "figcolors.py"
_spec = importlib.util.spec_from_file_location("figcolors", _fc)
FC = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(FC)                       # palette assertions run on import
from rules_env import metrics as MT  # noqa: E402


def logged(name):
    p = RES / f"{name}.json"
    r = subprocess.run([sys.executable, "-I", "-B", str(LEDGER / "reread_log.py"), str(p), "보충 그림 S6 그리기(rules_env 결과 · 새 계산 없음)"],
                       capture_output=True, text=True)
    h = r.stdout.strip().splitlines()[-1] if r.returncode == 0 and r.stdout.strip() else None
    raw = p.read_bytes()
    if h != hashlib.sha256(raw).hexdigest():
        stop(f"{p} not logged or changed")
    return json.loads(raw)


s1a, s1b, s2, e1, e4 = (logged(n) for n in ("s1a", "s1b", "s2", "e1", "e4"))

FONT = "Liberation Sans"
if FONT not in {f.name for f in font_manager.fontManager.ttflist}:
    stop("Liberation Sans is not installed")
PT_TEXT, PT_LETTER = 7.0, 8.5
INK, INK2, MUTED, GRID, NEUTRAL = FC.INK, FC.INK2, FC.MUTED, FC.GRID, FC.NEUTRAL
AXIS = "#9a988f"
MM = 1 / 25.4
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": [FONT], "font.size": PT_TEXT, "axes.labelsize": PT_TEXT, "axes.titlesize": PT_TEXT,
    "xtick.labelsize": PT_TEXT, "ytick.labelsize": PT_TEXT, "legend.fontsize": PT_TEXT, "axes.edgecolor": AXIS, "axes.linewidth": 0.5,
    "axes.labelcolor": INK, "axes.spines.top": False, "axes.spines.right": False, "xtick.color": INK2, "ytick.color": INK2,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2.5, "ytick.major.size": 2.5, "lines.linewidth": 1.0,
    "legend.frameon": False, "mathtext.default": "regular", "pdf.fonttype": 42, "ps.fonttype": 42, "axes.unicode_minus": True,
    "figure.facecolor": "#ffffff", "axes.facecolor": "#ffffff", "savefig.facecolor": "#ffffff", "lines.solid_capstyle": "round",
})
PART = {"sub-01": FC.OI["blue"], "sub-02": FC.OI["orange"], "sub-03": FC.OI["green"]}
GEN = {"T-in": FC.OI["blue"], "T-app": FC.OI["orange"], "T-none": FC.OI["green"]}
GMARK = {"T-in": "o", "T-app": "s", "T-none": "^"}
TORT = {"T-in": "true fraction", "T-app": "apparent fraction", "T-none": "fixed"}   # v17j (2026-11-20): printed names of the stored generator keys

fig = plt.figure(figsize=(190 * MM, 128 * MM))
axs = fig.subplots(2, 2, gridspec_kw=dict(hspace=0.55, wspace=0.32, left=0.08, right=0.98, top=0.95, bottom=0.09))


def letter(ax, s):
    ax.text(-0.13, 1.06, f"({s})", transform=ax.transAxes, fontsize=PT_LETTER, fontweight="bold", color=INK, ha="left", va="bottom")


# (a) E1
ax = axs[0, 0]
d = np.linspace(0, 80, 321)
for i, (sb, col) in enumerate(PART.items()):
    v = e1["by"][f"x|{sb}"]
    s_u = np.full(v["n_units"], v["s_rt_corrected"])
    ax.plot(d, [MT._te_share(x, s_u, v["tau"], 1.96) for x in d], color=col, lw=1.0, label=sb.replace("sub-", "Participant "))
    ax.plot([70 + (i - 1) * 1.6], [v["ors_62_132_observed"]], marker="o", ms=4.5, color=col, mec="#ffffff", mew=0.6, ls="none", zorder=4)
tf = s1b["wong"]["by"]["x|T-in|L0|iii"]["tau_fit"]
s_med = float(np.median([e1["by"][f"x|{sb}"]["s_rt_corrected"] for sb in PART]))
ax.plot(d, [MT._te_share(x, np.full(27, s_med), tf, 1.96) for x in d], color=MUTED, lw=1.0, ls=(0, (3, 2)), label="Fitter alone")
ax.axhline(0.90, color=INK2, lw=0.5, ls=(0, (1, 1.5)))
ax.set_xlim(0, 80)
ax.set_ylim(0.3, 1.0)
ax.set_xlabel("TE difference (ms)")
ax.set_ylabel("Share at retest level, x")
ax.legend(loc="lower left", handlelength=1.6)
letter(ax, "a")

# (b) E4
ax = axs[0, 1]
for k, v in e4["wong"].items():
    ta, tb = (int(t.replace("TE", "")) for t in k.split("-"))
    dte = tb - ta
    ax.plot([dte - 0.8], [v["x"]], marker="o", ms=3.5, color=FC.OI["blue"], ls="none")
    ax.plot([dte + 0.8], [v["raw_ndi"]], marker="s", ms=3.5, color=FC.OI["orange"], ls="none")
ax.plot([], [], marker="o", ms=3.5, color=FC.OI["blue"], ls="none", label="x (logit NDI)")
ax.plot([], [], marker="s", ms=3.5, color=FC.OI["orange"], ls="none", label="raw NDI")
ax.axhline(0.90, color=INK2, lw=0.5, ls=(0, (1, 1.5)))
ax.set_xlim(0, 75)
ax.set_ylim(0.3, 1.0)
ax.set_xlabel("TE difference (ms)")
ax.set_ylabel("Share at retest level")
ax.legend(loc="lower left", handlelength=1.0)
letter(ax, "b")

# (c) S1a
ax = axs[1, 0]
classes = [("Q4", "TR only"), ("Q1", "TE only"), ("Q3", "Scheme only"), ("Q2", "TE, TR, scheme")]
for i, (q, nm) in enumerate(classes):
    for j, g in enumerate(GEN):
        o = [v["ors"] for k, v in s1a["settings"].items() if k.startswith(f"{g}|{q}|")]
        ax.plot([i + (j - 1) * 0.18], [min(o)], marker=GMARK[g], ms=4, color=GEN[g], ls="none")
for g in GEN:
    ax.plot([], [], marker=GMARK[g], ms=4, color=GEN[g], ls="none", label=TORT[g])
ax.axhline(0.90, color=INK2, lw=0.5, ls=(0, (1, 1.5)))
ax.set_xticks(range(len(classes)), [nm for _, nm in classes])
ax.set_xlim(-0.6, len(classes) - 0.4)
ax.set_ylim(0.80, 1.03)
ax.set_ylabel("Lowest share at retest level\n(17 relaxation settings)")
ax.legend(loc="lower right", ncol=1, handlelength=1.0, title="Tortuosity", title_fontsize=PT_TEXT)
letter(ax, "c")

# (d) S2 — ratio to Eq. (7), by scenario
ax = axs[1, 1]
scen = [("fig3", "Fig. 2 example"), ("l2en", "Diffuse,\nextra-neurite"), ("l2both", "Diffuse,\nboth compartments")]
meths = ("none", "ctrl_offset", "combat_loc", "combat_train_ctrl_loc")
ax.axhspan(0.80, 1.25, color=NEUTRAL, lw=0, zorder=0)
ax.axhline(1.0, color=INK2, lw=0.5)
for i, (sc, nm) in enumerate(scen):
    for j, g in enumerate(GEN):
        for k, m in enumerate(meths):
            v = s2["by"][f"pop|x|two|{sc}|{g}|{m}"]
            ax.plot([i + (j - 1) * 0.22 + (k - 1.5) * 0.04], [v["ratio"]], marker=GMARK[g], ms=3.5, color=GEN[g], ls="none")
for g in GEN:
    ax.plot([], [], marker=GMARK[g], ms=3.5, color=GEN[g], ls="none", label=TORT[g])
ax.set_xticks(range(len(scen)), [nm for _, nm in scen])
ax.set_xlim(-0.6, len(scen) - 0.4)
ax.set_ylim(0.75, 1.35)
ax.set_ylabel("Fitted change ÷ Eq. (7)")
ax.legend(loc="upper left", handlelength=1.0, title="Tortuosity", title_fontsize=PT_TEXT)
letter(ax, "d")

problems = [t.get_text() for t in fig.findobj(matplotlib.text.Text) if t.get_text() and t.get_fontsize() < PT_TEXT - 1e-6]
if problems:
    stop(f"text below {PT_TEXT} pt: {problems}")
OUT.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT / "FS6_env.pdf", metadata={"CreationDate": None, "ModDate": None})
fig.savefig(OUT / "FS6_env.png", dpi=300)
print("written", OUT / "FS6_env.pdf")
