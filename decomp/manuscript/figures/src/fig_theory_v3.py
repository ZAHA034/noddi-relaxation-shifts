"""Theory figures v3 — Fig. 1 (card plane), Fig. 2 (site TE x group T2), graphical abstract.
Ports drafts/figures/fig_theory.py (F2, F3) to the shared rule set (figstyle.py).

Physics and numbers are unchanged from fig_theory.py and are re-asserted against fig_theory.out:
  T2 82 / 52 / 1000 ms (intra-neurite, extra-neurite, free water; Gong et al. 2020a, 2020b)
  TE wedge 48.5-72.1 deg, TR fan 85.3-94.3 deg (stage-2 report, commit 68e36a1)
  +30 ms TE shift (0.211, 0.336), angle 57.9 deg; neurite-loss direction at nu 0.6 = 158.2 deg
  2026-11-04 (E12): Fig. 1 uses the two participants' mean T2, 82.5 / 52 ms -> shift (0.213, 0.334), angle 57.4 deg
  (recompute/a_values.json *_mean); Fig. 2 uses the same mean T2 since 2026-11-04 (H1; printed values unchanged).
  site TE 70 / 100 ms: A - control -0.086 -> -0.120; logit A - control changes by -0.148
Fig. 2 is fully analytic (no template panel).
Run: python3 -I fig_theory_v3.py
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Wedge, FancyArrowPatch, FancyBboxPatch, Rectangle  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

# ----------------------------------------------------------------------------- constants
TE_WEDGE = (48.5, 72.1)     # deg, stage-2 report (68e36a1)
TR_FAN = (85.3, 94.3)       # deg, stage-2 report (68e36a1)
T2_IN, T2_EX, T2_ISO = 82.0, 52.0, 1000.0
T2_IN_FIG1 = (82.0 + 83.0) / 2   # Figs. 1 and 2 — mean of the two participants (E12 · H1); T2_EX 52 is already the mean of 51 and 53
SITE_TE = (70.0, 100.0)

# Expected values copied from drafts/figures/fig_theory.out (old script output)
EXPECT = {
    "shift_TE30": (0.21328671328671367, 0.33363636363636373),     # Fig. 1 since 2026-11-04: recompute/a_values.json shift_30ms_mean
    "shift_angle_deg": 57.410073820719106,                         # a_values.json dir_deg_mean
    "axon_loss_angle_deg": 158.19859051364818,
    "A-control TE70": -0.0857477752721304,        # Fig. 2 since 2026-11-04 (H1): recompute/a_values.json patientA_diff_te70_te100_mean
    "A-control TE100": -0.12025460115930686,
    "B-control TE70": -0.060506259464341516,     # same model, nu 0.40 vs 0.46 at the mean T2 (decomp/stage1/config.py)
    "B-control TE100": -0.05849362556481941,
}

# ----------------------------------------------------------------------------- model
def card(nu, f):
    """x = logit(nu) = log(intra/extra); y = log(iso/intra) = logit(f) - ln(nu)."""
    return math.log(nu / (1 - nu)), math.log(f / (1 - f)) - math.log(nu)


def unit(deg):
    return np.array([math.cos(math.radians(deg)), math.sin(math.radians(deg))])


def nu_app(nu, te, t2in, t2ex):
    a = nu * np.exp(-te / t2in)
    b = (1 - nu) * np.exp(-te / t2ex)
    return a / (a + b)


def logit(p):
    return np.log(p / (1 - p))


def arrow(ax, p, q, color, lw=1.0, ms=7, z=3):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=ms, color=color,
                                 lw=lw, zorder=z, shrinkA=0, shrinkB=0, joinstyle="miter"))


def sector(ax, centre, r, a0, a1, color, alpha):
    for off in (0, 180):  # a line direction has both senses
        ax.add_patch(Wedge(centre, r, a0 + off, a1 + off, facecolor=color, alpha=alpha,
                           lw=0, zorder=1))


# ----------------------------------------------------------------------------- Fig. 1
def fig1():
    fig = fs.new_fig(190, 84)
    # (a) translation
    ax = fs.ax_mm(fig, 16, 12, 84, 64)
    nus, fis = (0.45, 0.60, 0.75), (0.02, 0.05, 0.12)
    te_ms = 30.0
    shift = te_ms * np.array([1 / T2_EX - 1 / T2_IN_FIG1, 1 / T2_IN_FIG1 - 1 / T2_ISO])
    for nu in nus:
        for f in fis:
            p = np.array(card(nu, f))
            q = p + shift
            arrow(ax, p, q, fs.TE_C, lw=0.9, ms=6, z=3)
            ax.plot(*p, "o", ms=4.0, color=fs.INK2, mec=fs.SURF, mew=0.6, zorder=4)
            ax.plot(*q, "o", ms=4.0, color=fs.TE_C, mec=fs.SURF, mew=0.6, zorder=4)
    ax.set_xlabel("x = logit NDI")
    ax.set_ylabel("y = logit FISO − ln NDI")
    ax.set_xlim(-0.35, 1.45)
    ax.set_ylim(-3.85, -0.25)
    ax.set_xticks([0.0, 0.4, 0.8, 1.2])
    fs.hgrid(ax, "both")
    ax.plot([], [], "o", ms=4.0, color=fs.INK2, mec=fs.SURF, mew=0.6, label="tissue state")
    ax.plot([], [], "o", ms=4.0, color=fs.TE_C, mec=fs.SURF, mew=0.6,
            label="same state, TE + 30 ms")
    ax.legend(loc="upper right", handletextpad=0.2, borderaxespad=0.3, labelspacing=0.35)
    fs.panel(ax, "a")
    ang = math.degrees(math.atan2(shift[1], shift[0]))

    # (b) directions
    bx = fs.ax_mm(fig, 112, 6, 74, 74)
    c = np.array([0.0, 0.0])
    sector(bx, c, 1.0, *TE_WEDGE, fs.TE_C, alpha=0.17)
    sector(bx, c, 1.0, *TR_FAN, fs.TR_C, alpha=0.32)
    bx.plot([-1.08, 1.08], [0, 0], color=fs.GRID, lw=0.5, zorder=0)
    bx.plot([0, 0], [-1.05, 1.05], color=fs.GRID, lw=0.5, zorder=0)
    for deg, lab, xy, ha, va in ((0, "R2 extra ↑", (0.60, 0.035), "left", "bottom"),
                                 (135, "R2 intra ↑", (-0.43, 0.43), "right", "bottom"),
                                 (270, "R2 free ↑", (0.04, -0.58), "left", "center")):
        arrow(bx, c, 0.56 * unit(deg), fs.MUTED, lw=0.6, ms=5, z=2)
        bx.text(*xy, lab, color=fs.INK2, ha=ha, va=va)
    nu = 0.6
    axon_deg = math.degrees(math.atan2(1 / nu, -1 / (nu * (1 - nu))))
    arrow(bx, c, 0.86 * unit(axon_deg), fs.NDI_C, lw=1.4, ms=8, z=4)
    p = 0.86 * unit(axon_deg)
    bx.text(p[0] + 0.02, p[1] - 0.09, "neurite loss", color=fs.INK, ha="center", va="top")
    arrow(bx, c, 0.86 * unit(180), fs.INK, lw=1.2, ms=7, z=4)
    bx.text(-0.86, -0.06, "extra-neurite T2 ↑", color=fs.INK, ha="left", va="top")
    arrow(bx, c, 0.86 * unit(90), fs.FISO_C, lw=1.4, ms=8, z=4)
    bx.text(-0.11, 0.90, "free water ↑ ≡ tissue T1 ↑", color=fs.INK, ha="right", va="center")
    # sector labels just outside the arcs
    te_mid = unit(sum(TE_WEDGE) / 2)
    bx.text(*(1.05 * te_mid), "TE", color=fs.TE_C, fontweight="bold", ha="left", va="bottom")
    bx.text(*(-1.05 * te_mid), "TE", color=fs.TE_C, fontweight="bold", ha="right", va="top")
    bx.text(0.07, -0.97, "TR", color=fs.INK, fontweight="bold", ha="left", va="center")
    bx.set_xlim(-1.12, 1.12)
    bx.set_ylim(-1.12, 1.12)
    bx.set_aspect("equal")
    bx.set_xticks([])
    bx.set_yticks([])
    for sp in bx.spines.values():
        sp.set_visible(False)
    fs.panel(bx, "b", dx_mm=-3.0, dy_mm=-2.0)
    fs.save(fig, "Fig1_card_plane")
    return {"shift_TE30": shift.tolist(), "shift_angle_deg": ang, "axon_loss_angle_deg": axon_deg}


# ----------------------------------------------------------------------------- Fig. 2
GROUPS = [  # label, nu, T2_in, T2_ex, line style
    ("Controls", 0.46, T2_IN_FIG1, 52.0, dict(color=fs.INK, ls="-", lw=1.2)),
    ("Relaxation-only", 0.46, T2_IN_FIG1, 70.0, dict(color=fs.INK2, ls=(0, (4.0, 2.0)), lw=1.2)),
    ("Fraction-only", 0.40, T2_IN_FIG1, 52.0, dict(color=fs.MUTED, ls="-", lw=1.2)),
]
GMARK = {"Controls": dict(marker="o", mfc=fs.INK, mec=fs.INK),
         "Relaxation-only": dict(marker="o", mfc=fs.SURF, mec=fs.INK2),
         "Fraction-only": dict(marker="D", mfc=fs.MUTED, mec=fs.MUTED)}


def bracket(ax, x, y0, y1, text, side="left", color=fs.INK2, dx=1.6):
    ax.annotate("", (x, y1), (x, y0), arrowprops=dict(arrowstyle="<->", color=color, lw=0.6,
                                                      shrinkA=0.8, shrinkB=0.8,
                                                      mutation_scale=5))
    ha = "right" if side == "left" else "left"
    xt = x - dx if side == "left" else x + dx
    return ax.text(xt, (y0 + y1) / 2, text, color=color, ha=ha, va="center")


def fig2():
    fig = fs.new_fig(140, 70)          # two panels at their former size
    te = np.linspace(0, 120, 241)
    axa = fs.ax_mm(fig, 15, 11, 52, 49)
    axb = fs.ax_mm(fig, 82, 11, 52, 49)
    out = {}
    handles = []
    for lab, nu, a, b, st in GROUPS:
        y = nu_app(nu, te, a, b)
        axa.plot(te, y, **st, zorder=3)
        axb.plot(te, logit(y), **st, zorder=3)
        out[lab] = {t: float(nu_app(nu, t, a, b)) for t in SITE_TE}
        for t in SITE_TE:
            v = nu_app(nu, t, a, b)
            axa.plot(t, v, ms=3.6, mew=0.8, zorder=4, ls="none", **GMARK[lab])
            axb.plot(t, logit(v), ms=3.6, mew=0.8, zorder=4, ls="none", **GMARK[lab])
        handles.append(Line2D([], [], ms=3.6, mew=0.8, label=lab, **st, **GMARK[lab]))
    for ax in (axa, axb):
        for i, t in enumerate(SITE_TE):
            ax.axvline(t, color=fs.AXIS, lw=0.6, zorder=1.5)
        fs.hgrid(ax, "y")
        ax.set_xlabel("TE (ms)")
        ax.set_xlim(0, 120)
        ax.set_xticks([0, 20, 40, 60, 80, 100, 120])
    axa.set_ylim(0.40, 0.68)
    axa.set_ylabel("Apparent NDI")
    axb.set_ylim(-0.45, 0.70)
    axb.set_ylabel("logit apparent NDI")
    for ax, ytop in ((axa, 0.68), (axb, 0.70)):
        for i, t in enumerate(SITE_TE):
            ax.text(t, ytop, ("site B", "site A")[i], color=fs.INK2, ha="center", va="bottom")   # as the caption and Section 2.4: TE_B 70 ms, TE_A 100 ms
    axa.legend(handles=handles, loc="upper left", handlelength=2.6, handletextpad=0.5,
               borderaxespad=0.2, labelspacing=0.35)
    # (a) A - control brackets at both sites
    c, pa, pb = out["Controls"], out["Relaxation-only"], out["Fraction-only"]
    d70 = pa[70.0] - c[70.0]
    d100 = pa[100.0] - c[100.0]
    bracket(axa, 70.0, pa[70.0], c[70.0], f"{d70:+.3f}".replace("-", "−"), side="left")
    t2 = axa.annotate("", (100.0, c[100.0]), (100.0, pa[100.0]),
                      arrowprops=dict(arrowstyle="<->", color=fs.INK2, lw=0.6, shrinkA=0.8,
                                      shrinkB=0.8, mutation_scale=5))
    axa.text(101.6, 0.553, f"{d100:+.3f}".replace("-", "−"),
             color=fs.INK2, ha="left", va="center")
    # (b) call-outs
    y0 = float(logit(0.46))
    axb.annotate("same tissue,\ndifferent T2", xy=(1.5, y0), xytext=(9, 0.45),
                 color=fs.INK2, ha="left", va="center",
                 arrowprops=dict(arrowstyle="-", color=fs.MUTED, lw=0.5, shrinkA=1, shrinkB=1))
    yb = float(logit(nu_app(0.40, 32.0, T2_IN_FIG1, 52.0)))
    axb.annotate("fraction\ndifference", xy=(32, yb), xytext=(40, -0.33),
                 color=fs.INK2, ha="left", va="center",
                 arrowprops=dict(arrowstyle="-", color=fs.MUTED, lw=0.5, shrinkA=1, shrinkB=1))
    dlog = (logit(pa[100.0]) - logit(c[100.0])) - (logit(pa[70.0]) - logit(c[70.0]))
    fs.panel(axa, "a", dy_mm=4.5)
    fs.panel(axb, "b", dy_mm=4.5)
    fs.save(fig, "Fig2_site_TE_group_T2")
    diffs = {f"A-control TE{int(t)}": pa[t] - c[t] for t in SITE_TE}
    diffs.update({f"B-control TE{int(t)}": pb[t] - c[t] for t in SITE_TE})
    diffs["logit A-control change TE70->100"] = float(dlog)
    return diffs


# ----------------------------------------------------------------------------- graphical abstract
def _box(ax, x, y, w, h, fc="#f7f6f3", ec=fs.AXIS):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.6",
                                fc=fc, ec=ec, lw=0.5, zorder=1))


def _flow(ax, x0, x1, y):
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle="-|>", mutation_scale=7,
                                 color=fs.INK2, lw=0.9, zorder=2, shrinkA=0, shrinkB=0))


def graphical_abstract():
    W, H = 130.0, 50.0
    fig = plt.figure(figsize=(W * fs.MM, H * fs.MM))
    ax = fs.canvas(fig)
    ax.set_xticks([])
    ax.set_yticks([])
    bw, bh, by, gap = 27.5, 39.0, 5.5, 5.0
    xs = [2.5 + i * (bw + gap) for i in range(4)]
    heads = ["Single-TE index", "Which comparisons?", "Identifiability map",
             "Travelling heads"]
    for x, hd in zip(xs, heads):
        _box(ax, x, by, bw, bh)
        ax.text(x + bw / 2, by + bh - 3.0, hd, ha="center", va="top", fontweight="bold",
                color=fs.INK)
    for i in range(3):
        _flow(ax, xs[i] + bw + 0.6, xs[i + 1] - 0.6, by + bh / 2)

    # 1 — index = tissue x weighting, with a mini card plane translation
    x = xs[0]
    ax.text(x + bw / 2, by + 28.5, "tissue property", ha="center", va="center", color=fs.INK)
    ax.text(x + bw / 2, by + 24.6, "×", ha="center", va="center", color=fs.INK)
    ax.text(x + bw / 2, by + 20.7, "TE, TR, T2, T1", ha="center", va="center", color=fs.TE_C)
    ax.text(x + bw / 2, by + 17.3, "weighting", ha="center", va="center", color=fs.TE_C)
    o = np.array([x + 6.0, by + 4.2])
    for k, (dx, dy) in enumerate(((0, 0), (6.5, 0.6), (1.2, 5.6), (8.0, 6.0))):
        p = o + np.array([dx, dy])
        q = p + np.array([3.4, 2.6])
        ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=4.5, color=fs.TE_C,
                                     lw=0.6, zorder=3, shrinkA=0, shrinkB=0))
        ax.plot(*p, "o", ms=2.6, color=fs.INK2, mec=fs.SURF, mew=0.4, zorder=4)
        ax.plot(*q, "o", ms=2.6, color=fs.TE_C, mec=fs.SURF, mew=0.4, zorder=4)

    # 2 — four comparison types
    x = xs[1]
    rows = ["same protocol", "different TE", "different TR", "multi-site"]
    for j, r in enumerate(rows):
        yy = by + 28.0 - j * 6.4
        ax.add_patch(Rectangle((x + 3.0, yy - 2.3), bw - 6.0, 4.6, fc=fs.SURF, ec=fs.AXIS,
                               lw=0.5, zorder=2))
        ax.text(x + bw / 2, yy, r, ha="center", va="center", color=fs.INK)
    ax.text(x + bw / 2, by + 3.4, "NDI · FISO · ODI", ha="center", va="center", color=fs.INK2)

    # 3 — identifiability map: neutral grid (no results)
    x = xs[2]
    gx0, gy0, cw, ch = x + 6.0, by + 11.5, 4.3, 4.3
    cols = [fs.NDI_C, fs.FISO_C, fs.ODI_C]
    for i in range(4):
        for j in range(3):
            ax.add_patch(Rectangle((gx0 + i * cw, gy0 + j * ch), cw, ch, fc=fs.SURF,
                                   ec=fs.AXIS, lw=0.4, zorder=2))
    for j, (cname, col) in enumerate(zip(("NDI", "FISO", "ODI"), cols)):
        ax.plot(gx0 - 1.6, gy0 + (2 - j) * ch + ch / 2, fs.CMRK[cname], ms=3.0, color=col,
                zorder=3)
    ax.text(x + bw / 2, by + 6.4, "smallest identifiable", ha="center", va="center",
            color=fs.INK2)
    ax.text(x + bw / 2, by + 3.4, "difference", ha="center", va="center", color=fs.INK2)

    # 4 — travelling heads: one head, several scanners
    x = xs[3]
    hc = np.array([x + bw / 2, by + 22.0])
    ax.add_patch(plt.Circle(hc, 3.0, fc=fs.SURF, ec=fs.INK, lw=0.7, zorder=3))
    for k, dx in enumerate((-9.0, 0.0, 9.0)):
        sc = np.array([x + bw / 2 + dx, by + 11.5])
        ax.add_patch(FancyBboxPatch((sc[0] - 3.0, sc[1] - 2.2), 6.0, 4.4,
                                    boxstyle="round,pad=0,rounding_size=1.0", fc=fs.GRID,
                                    ec=fs.INK2, lw=0.5, zorder=3))
        ax.add_patch(FancyArrowPatch(hc + np.array([dx * 0.25, -3.2]), sc + np.array([0, 2.4]),
                                     arrowstyle="-|>", mutation_scale=4.5, color=fs.INK2,
                                     lw=0.5, zorder=2, shrinkA=0, shrinkB=0))
    ax.text(x + bw / 2, by + 4.5, "true difference 0", ha="center", va="center", color=fs.INK2)

    probs = [p for p in fs.audit(fig, "GraphicalAbstract") if not p.startswith("width")]
    if probs:
        raise SystemExit(f"GraphicalAbstract: {probs}")
    base = os.path.join(fs.OUTDIR, "GraphicalAbstract")
    fig.savefig(base + ".pdf")
    fig.savefig(base + ".png", dpi=300)
    plt.close(fig)
    print("wrote", base + ".pdf/.png (130 x 50 mm; width check skipped — GA is not a column figure)")


# ----------------------------------------------------------------------------- main
def check(name, got, want):
    assert round(got, 3) == round(want, 3), (name, got, want)
    print(f"  {name:28s} {got:+.4f}  (fig_theory.out {want:+.4f})  OK")


if __name__ == "__main__":
    r1 = fig1()
    d2 = fig2()
    graphical_abstract()
    print("asserted against fig_theory.out:")
    check("shift x (TE +30 ms)", r1["shift_TE30"][0], EXPECT["shift_TE30"][0])
    check("shift y (TE +30 ms)", r1["shift_TE30"][1], EXPECT["shift_TE30"][1])
    check("shift angle (deg)", r1["shift_angle_deg"], EXPECT["shift_angle_deg"])
    check("neurite-loss angle (deg)", r1["axon_loss_angle_deg"], EXPECT["axon_loss_angle_deg"])
    for k in ("A-control TE70", "A-control TE100", "B-control TE70", "B-control TE100"):
        check(k, d2[k], EXPECT[k])
    assert round(d2["logit A-control change TE70->100"], 3) == -0.148
    print(f"  logit A-control change      {d2['logit A-control change TE70->100']:+.4f}  (caption −0.148)  OK")
