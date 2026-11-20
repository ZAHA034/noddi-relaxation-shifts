"""Shared figure style for the NODDI manuscript (NeuroImage / Elsevier artwork rules).

One rule set for every figure, so that all figures read as one system.

Elsevier artwork rules this file enforces (checked by audit()):
  * width 90 mm (single column), 140 mm (1.5 column) or 190 mm (double column); height <= 240 mm
  * text >= 7 pt at print size; no sub- or superscripts (they would fall below 6 pt):
    write T2, T1, Neff, R2 on the baseline
  * font Arial (Liberation Sans is the metric-compatible substitute when Arial is not
    installed; Elsevier accepts "fonts that look similar"); fonts embedded (TrueType, type 42)
  * line widths 0.25-1.5 pt
  * no figure title inside the figure (the title is the first sentence of the caption);
    panel letters (a), (b), ... and at most a few words per panel
  * vector PDF for submission; PNG at 300 dpi only for previews

Colour has one meaning across all figures:
  NDI  blue   #2a78d6  circle
  FISO orange #eb6834  square
  ODI  green  #1baf7a  triangle
  TE   violet #6250d6  (wedge, TE changes)
  TR   amber  #eda100  (fan, TR changes)
  subjects / devices: greys; pass vs fail is shown by FILL (filled = pass / identified,
  open = fail / not identified), never by colour.
Templates (synthetic data) carry a grey "TEMPLATE · synthetic data" mark in every data panel;
final figures drawn from repository results must not.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

# ----------------------------------------------------------------------------- palette
NDI_C, FISO_C, ODI_C = "#2a78d6", "#eb6834", "#1baf7a"
TE_C, TR_C = "#6250d6", "#eda100"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURF = "#e8e7e2", "#9a988f", "#ffffff"
ZONE_PART, ZONE_NOT = "#f5f4f1", "#e9e8e3"
GREYS = ["#0b0b0b", "#5f5d57", "#8a8880"]   # subjects: contrast on white 19.7 / 6.6 / 3.6 : 1 (2026-11-01; was #6f6d67 · #aaa8a0 at 2.4 : 1)
WM_C = "#b4b2aa"

COMPONENTS = ("NDI", "FISO", "ODI")
CCOL = {"NDI": NDI_C, "FISO": FISO_C, "ODI": ODI_C}
CMRK = {"NDI": "o", "FISO": "s", "ODI": "^"}
CMS = {"NDI": 4.0, "FISO": 3.7, "ODI": 4.3}

# ----------------------------------------------------------------------------- sizes
PT_TEXT = 7.0          # every label, tick, annotation
PT_LETTER = 8.5        # panel letters
LW_AXIS = 0.5
LW_DATA = 1.0
LW_THIN = 0.5
LW_MIN = 0.25
WIDTHS_MM = (90.0, 140.0, 190.0)
MAX_H_MM = 240.0
MM = 1 / 25.4


def _font():
    names = {f.name for f in font_manager.fontManager.ttflist}
    return "Arial" if "Arial" in names else "Liberation Sans"


FONT = _font()

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": [FONT],
    "font.size": PT_TEXT, "axes.labelsize": PT_TEXT, "axes.titlesize": PT_TEXT,
    "xtick.labelsize": PT_TEXT, "ytick.labelsize": PT_TEXT, "legend.fontsize": PT_TEXT,
    "axes.edgecolor": AXIS, "axes.linewidth": LW_AXIS, "axes.labelcolor": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.color": INK2, "ytick.color": INK2,
    "xtick.major.width": LW_AXIS, "ytick.major.width": LW_AXIS,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "xtick.minor.size": 0, "ytick.minor.size": 0,
    "xtick.major.pad": 2, "ytick.major.pad": 2, "axes.labelpad": 3,
    "lines.linewidth": LW_DATA, "patch.linewidth": LW_THIN,
    "legend.frameon": False, "mathtext.default": "regular",
    "mathtext.fontset": "custom", "mathtext.rm": FONT, "mathtext.it": FONT + ":italic",
    "mathtext.bf": FONT + ":bold", "mathtext.sf": FONT, "mathtext.cal": FONT, "mathtext.tt": FONT,
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "lines.solid_capstyle": "round", "axes.unicode_minus": True,
})

OUTDIR = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
TEMPLATE_TEXT = "TEMPLATE · synthetic data"


# ----------------------------------------------------------------------------- layout helpers
def new_fig(w_mm, h_mm):
    assert any(abs(w_mm - w) < 1e-6 for w in WIDTHS_MM), f"width {w_mm} mm not 90/140/190"
    assert h_mm <= MAX_H_MM, f"height {h_mm} mm > {MAX_H_MM}"
    return plt.figure(figsize=(w_mm * MM, h_mm * MM))


def ax_mm(fig, x, y, w, h, **kw):
    """Axes placed in millimetres from the figure's lower-left corner."""
    W, H = fig.get_size_inches() / MM
    return fig.add_axes([x / W, y / H, w / W, h / H], **kw)


def canvas(fig):
    """Full-figure axes in millimetre units (for keys, tables, connectors)."""
    W, H = fig.get_size_inches() / MM
    c = fig.add_axes([0, 0, 1, 1], zorder=-1)
    c.set_xlim(0, W)
    c.set_ylim(0, H)
    c.set_xticks([])
    c.set_yticks([])
    c.axis("off")
    return c


def panel(ax, letter, header=None, dx_mm=-7.0, dy_mm=2.0):
    """Panel letter (bold, 8.5 pt) at the axes' top-left, optional short header (<= 5 words)."""
    pts = 72.0 / 25.4
    ax.annotate(f"({letter})", (0, 1), xycoords="axes fraction",
                xytext=(dx_mm * pts, dy_mm * pts), textcoords="offset points",
                fontsize=PT_LETTER, fontweight="bold", color=INK,
                ha="left", va="bottom", annotation_clip=False)
    if header:
        assert len(header.split()) <= 6, "panel header too long — move it to the caption"
        ax.annotate(header, (0, 1), xycoords="axes fraction",
                    xytext=(dx_mm * pts + 18, dy_mm * pts + 0.6), textcoords="offset points",
                    fontsize=PT_TEXT, color=INK2, ha="left", va="bottom",
                    annotation_clip=False)


def template_mark(ax, x=0.99, y=0.01, ha="right", va="bottom"):
    t = ax.text(x, y, TEMPLATE_TEXT, transform=ax.transAxes, fontsize=PT_TEXT,
                color=WM_C, ha=ha, va=va, zorder=0.5)
    t.set_in_layout(False)
    return t


def hgrid(ax, axis="y"):
    ax.grid(True, axis=axis, color=GRID, lw=LW_THIN)
    ax.set_axisbelow(True)


def swarm(vals, dy, dx, max_off):
    """Beeswarm offsets so that no two marks closer than dy in value share a slot."""
    vals = np.asarray(vals, float)
    off = np.zeros(len(vals))
    placed = []
    cands = [0.0]
    k = 1
    while k * dx <= max_off + 1e-9:
        cands += [k * dx, -k * dx]
        k += 1
    for i in np.argsort(vals):
        for o in cands:
            if all(abs(vals[i] - v) >= dy or abs(o - q) >= dx * 0.99 for v, q in placed):
                off[i] = o
                break
        else:
            off[i] = cands[len(placed) % len(cands)]
        placed.append((vals[i], off[i]))
    return off


# ----------------------------------------------------------------------------- audit + save
def audit(fig, name):
    """Fail loudly if the figure breaks a rule. Returns a list of problems (empty = pass)."""
    problems = []
    W, H = fig.get_size_inches() / MM
    if not any(abs(W - w) < 0.05 for w in WIDTHS_MM):
        problems.append(f"width {W:.1f} mm")
    if H > MAX_H_MM + 0.05:
        problems.append(f"height {H:.1f} mm")
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    fb = fig.bbox
    all_tick_labels, drawn_tick_labels = set(), set()
    for a in fig.axes:
        if not a.axison:          # hidden axes draw no tick labels
            for axis in (a.xaxis, a.yaxis):
                for tk in axis.get_major_ticks() + axis.get_minor_ticks():
                    all_tick_labels.update({id(tk.label1), id(tk.label2)})
            continue
        for axis in (a.xaxis, a.yaxis):
            for tk in axis.get_major_ticks() + axis.get_minor_ticks():
                all_tick_labels.update({id(tk.label1), id(tk.label2)})
            for tk in axis._update_ticks():
                drawn_tick_labels.update({id(tk.label1), id(tk.label2)})
    for t in fig.findobj(matplotlib.text.Text):
        s = t.get_text()
        if not s or not t.get_visible():
            continue
        if id(t) in all_tick_labels and id(t) not in drawn_tick_labels:
            continue
        if t.get_fontsize() < PT_TEXT - 1e-6:
            problems.append(f"text {s!r} at {t.get_fontsize():.2f} pt")
        if "$" in s and ("_" in s or "^" in s):
            problems.append(f"sub/superscript in {s!r}")
        bb = t.get_window_extent(r)
        if bb.x0 < fb.x0 - 0.5 or bb.y0 < fb.y0 - 0.5 or bb.x1 > fb.x1 + 0.5 or bb.y1 > fb.y1 + 0.5:
            problems.append(f"text {s!r} outside the figure")
    for ln in fig.findobj(matplotlib.lines.Line2D):
        if ln.get_visible() and ln.get_linestyle() not in ("None", "", " ") and len(ln.get_xdata()) > 1:
            if 0 < ln.get_linewidth() < LW_MIN - 1e-6 or ln.get_linewidth() > 1.5 + 1e-6:
                problems.append(f"line width {ln.get_linewidth():.2f} pt")
    if problems:
        print(f"AUDIT FAIL {name}:")
        for p in problems:
            print("  -", p)
    else:
        print(f"AUDIT PASS {name}: {W:.0f} x {H:.0f} mm, font {FONT}")
    return problems


def save(fig, name, strict=True):
    probs = audit(fig, name)
    if strict and probs:
        raise SystemExit(f"{name}: {len(probs)} rule problem(s); not saved")
    base = os.path.join(OUTDIR, name)
    fig.savefig(base + ".pdf", metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(base + ".png", dpi=300)
    plt.close(fig)
    print("wrote", base + ".pdf/.png")
