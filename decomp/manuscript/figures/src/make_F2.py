"""make_F2 draws Figure 3 (Section 2.2): compartment weighting in log-ratio coordinates, from the signal model and ledger values only, computing no statistic."""
import sys

sys.dont_write_bytecode = True        # first statement: no __pycache__ in the repository's figures/src or the ledger folder

import ast  # noqa: E402
import csv  # noqa: E402
import hashlib  # noqa: E402
import importlib.util  # noqa: E402
import itertools  # noqa: E402
import math  # noqa: E402
import os  # noqa: E402
import pathlib  # noqa: E402
import re  # noqa: E402
import subprocess  # noqa: E402
from decimal import Decimal, ROUND_HALF_UP  # noqa: E402

os.environ["CUDA_VISIBLE_DEVICES"] = ""   # CPU only

REPO = pathlib.Path(os.environ.get("DMRI_REPO") or pathlib.Path(__file__).resolve().parents[4])
WORK = pathlib.Path(os.environ.get("DMRI_WORK") or REPO / "data")
LEDGER_DIR = WORK / "stage1_work" / "drafts" / "x15" / "ledger"
SOURCES = LEDGER_DIR / "sources"
OUT = pathlib.Path(os.environ.get("FIG_OUT") or WORK / "stage1_work" / "drafts" / "figs_2026-11-06")
OUT.mkdir(parents=True, exist_ok=True)
NAME = "F2_logratio"

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.patches import Circle, FancyArrowPatch, Wedge  # noqa: E402
from matplotlib.text import Text  # noqa: E402
from matplotlib.transforms import offset_copy  # noqa: E402

_spec = importlib.util.spec_from_file_location("ledger_check", LEDGER_DIR / "ledger_check.py")
LC = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(LC)          # read-only use: eval_field · same_value · rounding_ok · numbers_in

LEDGER = {r["ID"]: r for r in csv.DictReader((LEDGER_DIR / "numbers.csv").open(encoding="utf-8"))}
CHECKED = []                          # (ID, raw value, manuscript text, role) for F2_values.md
VISIBLE = []                          # (where, text, IDs)
REGISTERED = {}                       # id(Text) -> (where, kind) for the figure-text audit
MIN_PT = 7.0


def stop(msg):
    sys.exit(f"STOP F1: {msg}")


def fmt(v, k):
    """Round half up to k decimals, as the ledger's rounding rules do; minus sign U+2212."""
    d = Decimal(repr(float(v))).quantize(Decimal(1).scaleb(-k), rounding=ROUND_HALF_UP)
    s = f"{abs(d):.{k}f}"
    return ("−" + s) if d < 0 else s


def norm(s):
    return str(s).replace("−", "-").replace(" ", "").replace("\n", "").replace(" ", "").lower()


def row(id_):
    r = LEDGER.get(id_)
    if r is None:
        stop(f"{id_} is not in the ledger")
    return r


def reproduce(id_, role):
    """The row reproduces from its own copy (literature rows: not machine-checked) and its rounding rule holds.
    Returns the raw value(s) as floats."""
    r = row(id_)
    try:
        v = LC.eval_field(SOURCES, r)
    except Exception as ex:                                   # noqa: BLE001 — any failure to reproduce is a stop
        stop(f"{id_} does not reproduce from its copy: {type(ex).__name__}: {ex}")
    if v is not None and not LC.same_value(v, r["원 값"]):
        stop(f"{id_}: copy gives {v!r}, ledger raw value {r['원 값']!r}")
    ok, why = LC.rounding_ok(r["원 값"], r["원고 표기"], r["반올림 규칙"])
    if ok is False:
        stop(f"{id_}: manuscript text «{r['원고 표기']}» does not follow the rounding rule ({why})")
    vals = [float(p.replace(",", "").replace("−", "-")) for p in r["원 값"].split("|")]
    CHECKED.append((id_, r["원 값"], r["원고 표기"], role + ("" if v is not None else " — literature row, not recomputed")))
    return vals


def val(id_, role):
    return reproduce(id_, role)[0]


def agree(id_, computed, k=None, rel=1e-12, abs_tol=None):
    """A value recomputed here equals the ledger's raw value (to `rel`, or within `abs_tol`) and prints as its text."""
    raw = float(row(id_)["원 값"].split("|")[0].replace(",", "").replace("−", "-"))
    close = abs(computed - raw) <= abs_tol if abs_tol is not None else math.isclose(computed, raw, rel_tol=rel, abs_tol=1e-15)
    if not close:
        stop(f"{id_}: recomputed {computed!r}, ledger raw {raw!r}")
    if k is not None and fmt(computed, k) != fmt(raw, k):
        stop(f"{id_}: recomputed prints {fmt(computed, k)}, ledger {fmt(raw, k)}")


# ============================================================================ the relaxation grid of Section 2.2
def config_copy():
    r = row("N-628")
    path = SOURCES / r["근거 파일"] / r["근거 파일 sha256"][:12]
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != r["근거 파일 sha256"]:
        stop(f"{path}: sha256 differs from the ledger")
    live = REPO / r["근거 파일"]
    if live.exists() and hashlib.sha256(live.read_bytes()).hexdigest() != r["근거 파일 sha256"]:
        stop(f"repository {r['근거 파일']} differs from its ledger copy")
    return raw.decode("utf-8")


def grid_from_config():
    """Evaluate only the grid assignments of config.py (no import, no other code runs)."""
    tree = ast.parse(config_copy())
    want = ("T2IN_GRID", "T1_SETS", "T1_ISO", "ETA_T2", "ACQ_E5")
    ns = {"__builtins__": {}}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name) \
                and node.targets[0].id in want:
            ns[node.targets[0].id] = eval(compile(ast.Expression(node.value), "config.py", "eval"), ns)   # noqa: S307
    missing = [w for w in want if w not in ns]
    if missing:
        stop(f"config.py: grid assignments not found: {missing}")
    return ns


def wedge_and_fan(g):
    """As decomp/stage1/analyze2.py x1: line directions of the TE change A -> Q1 (+180 deg: increasing TE) and of the TR change
    A -> Q4 over the 9 x 6 relaxation states of the grid."""
    def weights(te, t2, tr, t1):
        return [math.exp(-te / a) * (1 - math.exp(-tr / b)) for a, b in zip(t2, t1)]
    ang = {"Q1": [], "Q4": []}
    _, teA, trA = g["ACQ_E5"]["A"]
    for (t2in, r, t2iso), (tt, ratio) in itertools.product(g["ETA_T2"], g["T1_SETS"]):
        t2 = (t2in, t2in / r, t2iso)
        t1 = (tt * math.sqrt(ratio), tt / math.sqrt(ratio), g["T1_ISO"])
        aA = weights(teA, t2, trA, t1)
        for q in ang:
            _, te, tr = g["ACQ_E5"][q]
            aB = weights(te, t2, tr, t1)
            dx = math.log((aB[0] / aB[1]) / (aA[0] / aA[1]))
            dy = math.log((aB[2] / aB[0]) / (aA[2] / aA[0]))
            ang[q].append(math.degrees(math.atan2(dy, dx)))
    return (min(ang["Q1"]) + 180, max(ang["Q1"]) + 180), (min(ang["Q4"]), max(ang["Q4"]))


def read_values():
    V = {}
    g = grid_from_config()
    # the grid equals the rows that describe it in Section 2.2
    if sorted(g["T2IN_GRID"]) != reproduce("N-628", "intra-neurite T2 of the grid"):
        stop("config.py T2IN_GRID differs from N-628")
    rs = sorted({r for _, r, _ in g["ETA_T2"]})
    if rs != [val(i, "r of the grid") for i in ("N-051", "N-052", "N-053")]:
        stop(f"config.py ETA_T2 ratios {rs} differ from N-051-N-053")
    if {iso for *_, iso in g["ETA_T2"]} != {val("N-629", "free-water T2 of the grid")}:
        stop("config.py ETA_T2 free-water T2 differs from N-629")
    if g["T1_ISO"] != val("N-634", "free-water T1 of the grid"):
        stop("config.py T1_ISO differs from N-634")
    trs = (g["ACQ_E5"]["A"][2], g["ACQ_E5"]["Q4"][2])
    if trs != (val("N-630", "short TR of the TR change"), val("N-631", "long TR of the TR change")):
        stop(f"config.py ACQ_E5 TRs {trs} differ from N-630/N-631")
    if g["ACQ_E5"]["A"][1] != g["ACQ_E5"]["Q4"][1]:
        stop("the TR change of ACQ_E5 is not at one TE")
    te_w, tr_f = wedge_and_fan(g)
    for (id_lo, id_hi), (lo, hi) in ((("N-049", "N-050"), te_w), (("N-059", "N-060"), tr_f)):
        reproduce(id_lo, "x1.json end, guard")
        reproduce(id_hi, "x1.json end, guard")
        agree(id_lo, lo, rel=1e-12)
        agree(id_hi, hi, rel=1e-12)
    for id_, (lo, hi) in (("N-206", te_w), ("N-207", tr_f)):
        lo_l, hi_l = reproduce(id_, "printed range")
        if (fmt(lo, 1), fmt(hi, 1)) != (fmt(lo_l, 1), fmt(hi_l, 1)) or row(id_)["원고 표기"] != f"{fmt(lo, 1)}–{fmt(hi, 1)}°":
            stop(f"{id_}: recomputed {fmt(lo, 1)}–{fmt(hi, 1)}°, ledger «{row(id_)['원고 표기']}»")
    V["te_wedge"], V["tr_fan"], V["grid"] = te_w, tr_f, g

    # compartmental T2 (two participants' means; free water assumed) and the 30-ms TE step
    T2in, T2en, T2iso = val("N-644", "T2 intra"), val("N-043", "T2 extra"), val("N-046", "T2 free water")
    dTE = val("N-063", "TE step")
    vec = np.array([dTE * (1 / T2en - 1 / T2in), dTE * (1 / T2in - 1 / T2iso)])
    reproduce("N-064", "x shift")
    reproduce("N-065", "y shift")
    agree("N-064", vec[0], k=2)
    agree("N-065", vec[1], k=2)
    vdeg = math.degrees(math.atan2(vec[1], vec[0]))
    reproduce("N-602", "direction of the TE step")
    agree("N-602", vdeg, k=0, rel=1e-9)
    V.update(T2=(T2in, T2en, T2iso), dTE=dTE, vec=vec, vdeg=vdeg)

    # nine tissue states
    V["nus"] = [val(i, "tissue state NDI") for i in ("N-073", "N-074", "N-075")]
    V["fisos"] = [val(i, "tissue state FISO") for i in ("N-076", "N-077", "N-078")]

    # directions of tissue and relaxation changes, from d(x, y) = (d ln w_in - d ln w_en, d ln w_iso - d ln w_in)
    def dir_of(dlnw):
        d_in, d_en, d_iso = dlnw
        return math.degrees(math.atan2(d_iso - d_in, d_in - d_en)) % 360
    dirs = {"R2en up": dir_of((0, -1, 0)), "R2in up": dir_of((-1, 0, 0)), "R2iso up": dir_of((0, 0, -1)),
            "T2en up": dir_of((0, 1, 0)), "tissue T1 up": dir_of((-1, -1, 0))}
    nu0 = val("N-067", "NDI of the neurite-loss direction")
    # neurite loss at fixed FISO: d(x, y)/d(-nu) = (-1/(nu(1-nu)), 1/nu); free water at fixed nu: (0, +)
    dirs["neurite loss"] = math.degrees(math.atan2(1 / nu0, -1 / (nu0 * (1 - nu0))))
    dirs["free water"] = math.degrees(math.atan2(1.0, 0.0))
    for k, id_ in (("R2en up", "N-068"), ("R2in up", "N-069"), ("R2iso up", "N-070"), ("T2en up", "N-071"),
                   ("tissue T1 up", "N-072"), ("free water", "N-072")):
        reproduce(id_, f"direction: {k}")
        agree(id_, dirs[k], k=0, abs_tol=1e-9)
    reproduce("N-066", "direction: neurite loss")
    agree("N-066", dirs["neurite loss"], k=0, abs_tol=1e-3)   # the ledger's value is a finite difference (a_values.py)
    V["dirs"] = dirs
    return V


# ============================================================================ drawing
FONT = "Liberation Sans"
try:
    font_manager.findfont(font_manager.FontProperties(family=FONT), fallback_to_default=False)
except ValueError:
    stop(f"font {FONT} is not installed")
OI = dict(blue="#0072B2", orange="#E69F00", orange_text="#A06E00", green="#009E73", verm="#D55E00",
          sky="#56B4E9", grey="#6E6E6E", ink="#1A1A1A")                    # Okabe–Ito, fixed to meaning (direction S39 §4)
plt.rcParams.update({
    "font.family": FONT, "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5,
    "ytick.major.size": 2.5, "pdf.fonttype": 42, "ps.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False,
    "svg.fonttype": "none",
})
MM = 1 / 25.4
W, H = 190, 84                                                      # mm, double column


def text(ax, x, y, s, ids=(), where="", **kw):
    """A figure text; numbers in it must come from the listed ledger rows."""
    t = ax.text(x, y, s, **kw)
    REGISTERED[id(t)] = (where or s, "text")
    if ids or re.search(r"\d", s):
        have = []
        for i in ids:
            if norm(row(i)["원고 표기"]) not in norm(s):
                stop(f"«{s}» does not contain the ledger text «{row(i)['원고 표기']}» of {i}")
            have += [abs(v) for v in LC.numbers_in(row(i)["원고 표기"], ms=True)]
        for v in LC.numbers_in(s):
            if not any(math.isclose(abs(v), h, rel_tol=1e-12, abs_tol=1e-12) for h in have):
                stop(f"number {v} in «{s}» has no listed ledger row")
        VISIBLE.append((where or s.replace("\n", " "), s.replace("\n", " "), tuple(ids)))
    return t


def symbol(fig, ax, x, y, base, sub, tail, ha, va, color, base_pt=8.75, sub_pt=7.0, where=""):
    """Base letter (italic), subscript (upright, 7 pt, dropped) and tail as separate texts, placed by measured widths in
    points (public API: Text, offset_copy, renderer metrics) so that the subscript stays at 7 pt."""
    rend = fig.canvas.get_renderer()
    px2pt = 72 / fig.dpi
    segs = [(base, base_pt, "italic", 0.0), (sub, sub_pt, "normal", -0.24 * base_pt), (tail, 7.0, "normal", 0.0)]
    widths = []
    for s, pt, style, _ in segs:
        w, h, d = rend.get_text_width_height_descent(s, font_manager.FontProperties(family=FONT, size=pt, style=style),
                                                     ismath=False)
        widths.append(w * px2pt)
    _, hb, db = rend.get_text_width_height_descent("Rp", font_manager.FontProperties(family=FONT, size=base_pt),
                                                   ismath=False)
    hb, db = hb * px2pt, db * px2pt
    total = sum(widths) + 0.3
    x0 = {"left": 0.0, "right": -total, "center": -total / 2}[ha]
    y0 = {"baseline": 0.0, "bottom": db, "top": -(hb - db), "center": -(hb / 2 - db)}[va]
    out, cx = [], x0
    for (s, pt, style, dy), w in zip(segs, widths):
        tr = offset_copy(ax.transData, fig=fig, x=cx, y=y0 + dy, units="points")
        t = ax.text(x, y, s, transform=tr, ha="left", va="baseline", fontsize=pt, style=style, color=color)
        REGISTERED[id(t)] = (where or base + sub + tail, "symbol")      # symbol subscripts are names, not numbers
        out.append(t)
        cx += w + (0.3 if s == base else 0.0)
    return out


def draw(V):
    fig = plt.figure(figsize=(W * MM, H * MM))
    fig.canvas.draw()                                               # a renderer for the symbol labels

    def ax_mm(l, b, w, h):
        return fig.add_axes([l / W, b / H, w / W, h / H])

    def panel_label(x_mm, y_mm, letter, title):
        for s, kw in ((letter, dict(fontsize=8, fontweight="bold")), (title, dict(fontsize=7.5))):
            t = fig.text((x_mm if s == letter else x_mm + 4.2) / W, y_mm / H, s, va="baseline", **kw)
            REGISTERED[id(t)] = (s, "title")

    vec, vdeg, dirs = V["vec"], V["vdeg"], V["dirs"]
    nus, fisos = V["nus"], V["fisos"]

    def state(nu, f):
        return np.array([math.log(nu / (1 - nu)), math.log(f / (1 - f)) - math.log(nu)])

    states = {(nu, f): state(nu, f) for nu in nus for f in fisos}

    # ------------------------------------------------------------ (a)
    XL, YL = (-0.95, 1.65), (-3.95, -0.55)
    aw = 50.0
    ah = aw * (YL[1] - YL[0]) / (XL[1] - XL[0])                     # equal aspect: angles read as in (b, c)
    axa = ax_mm(12, 9, aw, ah)
    axa.set_xlim(*XL)
    axa.set_ylim(*YL)
    axa.set_aspect("equal")
    axa.set_xticks([-0.5, 0, 0.5, 1.0, 1.5])
    axa.set_yticks([-3.5, -3, -2.5, -2, -1.5, -1])
    axa.set_xticklabels(["−0.5", "0", "0.5", "1.0", "1.5"])
    axa.set_yticklabels(["−3.5", "−3.0", "−2.5", "−2.0", "−1.5", "−1.0"])
    axa.set_xlabel("x = logit NDI", labelpad=2)
    axa.set_ylabel("y = logit FISO − ln NDI", labelpad=2)
    for f in fisos:                                                 # iso-FISO polylines: the lattice that is translated
        base = np.array([states[(nu, f)] for nu in nus])
        axa.plot(base[:, 0], base[:, 1], color="#B5B5B5", lw=0.6, zorder=1)
        axa.plot(base[:, 0] + vec[0], base[:, 1] + vec[1], color=OI["sky"], lw=0.6, zorder=1)
    for p in states.values():
        q = p + vec
        axa.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=5, lw=0.8, color=OI["blue"],
                                      shrinkA=1.6, shrinkB=1.6, zorder=3))
        axa.plot(*p, "o", ms=3.2, mfc=OI["ink"], mec="white", mew=0.4, zorder=4)
        axa.plot(*q, "o", ms=3.2, mfc=OI["blue"], mec="white", mew=0.4, zorder=4)
    for f, id_ in zip(fisos, ("N-076", "N-077", "N-078")):
        p = states[(nus[0], f)]
        text(axa, p[0] - 0.08, p[1], f"FISO {row(id_)['원고 표기']}", ids=(id_,), ha="right", va="center", color=OI["grey"])
    for nu, id_ in zip(nus, ("N-073", "N-074", "N-075")):
        p = states[(nu, fisos[0])]
        text(axa, p[0], p[1] - 0.17, f"NDI {row(id_)['원고 표기']}", ids=(id_,), ha="center", va="top", color=OI["grey"])
    p = states[(nus[0], fisos[-1])]
    q = p + vec
    axa.plot([p[0], p[0], q[0]], [p[1], q[1], q[1]], color=OI["blue"], lw=0.5, ls=(0, (2, 1.5)))
    text(axa, p[0] - 0.05, (p[1] + q[1]) / 2, f"Δy {fmt(vec[1], 2)}", ids=("N-065",), ha="right", va="center",
         color=OI["blue"])
    text(axa, (p[0] + q[0]) / 2, q[1] + 0.07, f"Δx {fmt(vec[0], 2)}", ids=("N-064",), ha="center", va="bottom",
         color=OI["blue"])
    text(axa, 1.55, -0.72, f"same states,\nTE + {row('N-063')['원고 표기']}", ids=("N-063",), ha="right", va="top",
         color=OI["blue"], linespacing=1.1)
    text(axa, 0.55, -1.70, "tissue states", ha="center", va="top", color=OI["ink"])
    panel_label(2, H - 5, "(a)", "Every state moves by the same vector")

    # ------------------------------------------------------------ compass helpers
    R, LIM = 1.0, 1.5

    def compass(ax):
        ax.set_xlim(-LIM, LIM)
        ax.set_ylim(-LIM, LIM)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.add_patch(Circle((0, 0), R, fill=False, lw=0.5, color="#C8C8C8", zorder=0))
        ax.plot([-R, R], [0, 0], color="#C8C8C8", lw=0.5, zorder=0)
        ax.plot([0, 0], [-R, R], color="#C8C8C8", lw=0.5, zorder=0)
        for a in (0, 90, 180, 270):
            t = math.radians(a)
            ax.plot([0.94 * R * math.cos(t), R * math.cos(t)], [0.94 * R * math.sin(t), R * math.sin(t)],
                    color="#9A9A9A", lw=0.6)

    def ray(ax, deg, color, lw=1.4, r=0.92, offset=0.0):
        t = math.radians(deg)
        n = np.array([-math.sin(t), math.cos(t)]) * offset
        ax.add_patch(FancyArrowPatch(tuple(n), tuple(n + r * np.array([math.cos(t), math.sin(t)])), arrowstyle="-|>",
                                     mutation_scale=7, lw=lw, color=color, shrinkA=0, shrinkB=0, zorder=5))

    def at(deg, r):
        t = math.radians(deg)
        return r * math.cos(t), r * math.sin(t)

    def sectors(ax, a0, a1, color, alpha):
        for off in (0, 180):                                        # a line direction: both senses of the change
            ax.add_patch(Wedge((0, 0), R, a0 + off, a1 + off, facecolor=color, alpha=alpha, edgecolor="none", zorder=1))

    # ------------------------------------------------------------ (b)
    cw = 54.0
    axb = ax_mm(77, 10, cw, cw)
    compass(axb)
    te_w, tr_f = V["te_wedge"], V["tr_fan"]
    sectors(axb, *te_w, OI["blue"], 0.22)
    sectors(axb, *tr_f, OI["orange"], 0.45)
    axb.plot(*zip((0, 0), at(vdeg, R)), color=OI["blue"], lw=0.6, ls=(0, (2, 1.5)), zorder=2)
    ray(axb, dirs["free water"], OI["green"])
    ray(axb, dirs["neurite loss"], OI["green"])
    text(axb, *at(vdeg, 1.06), f"{row('N-602')['원고 표기']}, as in (a)", ids=("N-602",), ha="left", va="center",
         color=OI["blue"])
    text(axb, *at(dirs["free water"], 1.07), f"free water ↑ ({row('N-072')['원고 표기']})", ids=("N-072",),
         ha="center", va="bottom", color=OI["green"])
    text(axb, *at(dirs["neurite loss"], 1.05),
         f"neurite loss\n({row('N-066')['원고 표기']} at NDI {row('N-067')['원고 표기']},\nFISO fixed)",
         ids=("N-066", "N-067"), ha="right", va="center", color=OI["green"], linespacing=1.05)
    text(axb, *at(240, 1.06), f"TE\n{row('N-206')['원고 표기']}", ids=("N-206",), ha="right", va="top", color=OI["blue"],
         linespacing=1.05)
    text(axb, *at(272, 1.07), f"TR\n{row('N-207')['원고 표기']}", ids=("N-207",), ha="left", va="top",
         color=OI["orange_text"], linespacing=1.05)
    t0 = axb.text(*at(0, 1.04), "0°", ha="left", va="center", color=OI["grey"])   # axis reference, like a tick label
    REGISTERED[id(t0)] = ("compass 0° (axis reference)", "tick")
    panel_label(72, H - 5, "(b)", "Across sites: TE and TR sectors")

    # ------------------------------------------------------------ (c)
    axc = ax_mm(135, 10, cw, cw)
    compass(axc)
    axc.add_patch(Wedge((0, 0), 0.55, dirs["R2in up"], dirs["T2en up"], facecolor=OI["verm"], alpha=0.12,
                        edgecolor="none", zorder=1))
    for k in ("R2en up", "R2in up", "R2iso up", "T2en up"):
        ray(axc, dirs[k], OI["verm"], lw=1.0)
    ray(axc, dirs["tissue T1 up"], OI["verm"], lw=1.0, offset=0.07)  # tissue T1 up, drawn beside free water
    ray(axc, dirs["free water"], OI["green"], offset=-0.07)
    ray(axc, dirs["neurite loss"], OI["green"])
    symbol(fig, axc, 0.97, 0.05, "R", "2,en", " ↑", "right", "bottom", OI["verm"])
    symbol(fig, axc, *at(dirs["R2in up"], 1.04), "R", "2,in", " ↑", "right", "bottom", OI["verm"])
    _x, _y = at(dirs["T2en up"], 1.04)
    symbol(fig, axc, _x, _y + 0.02, "T", "2,en", " ↑", "right", "bottom", OI["verm"])
    text(axc, _x, _y - 0.02, "alone", ha="right", va="top", color=OI["verm"])
    symbol(fig, axc, *at(dirs["R2iso up"], 1.05), "R", "2,iso", " ↑", "center", "top", OI["verm"])
    text(axc, -0.10, 1.07, "free water ↑", ha="right", va="bottom", color=OI["green"])
    text(axc, 0.10, 1.07, "≡ tissue T1 ↑", ha="left", va="bottom", color=OI["verm"])
    text(axc, 0.10, 1.05, "(equal T1)", ha="left", va="top", color=OI["verm"])   # label added 2026-11-16: Proposition 2 holds when the tissue compartments share T1
    text(axc, *at(dirs["neurite loss"], 1.05), "neurite\nloss", ha="right", va="center", color=OI["green"],
         linespacing=1.05)
    panel_label(133, H - 5, "(c)", "Across groups: relaxation vs tissue")
    return fig


# ============================================================================ audits
def audit(fig):
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    fb = fig.bbox
    ticks = set()
    for ax in fig.axes:
        for t in ax.get_xticklabels() + ax.get_yticklabels():
            ticks.add(id(t))
    texts, small, outside, unregistered = [], [], [], []
    for t in fig.findobj(Text):
        s = t.get_text()
        if not s.strip() or not t.get_visible():
            continue
        if id(t) in ticks:
            if t.get_fontsize() < MIN_PT:
                small.append(s)
            continue
        if t.get_fontsize() < MIN_PT - 1e-9:
            small.append(f"{s} ({t.get_fontsize()} pt)")
        bb = Text.get_window_extent(t, rend)                         # the text alone (an annotation's arrow excluded)
        if bb.x0 < fb.x0 - 0.5 or bb.x1 > fb.x1 + 0.5 or bb.y0 < fb.y0 - 0.5 or bb.y1 > fb.y1 + 0.5:
            outside.append(s)
        if id(t) not in REGISTERED and t.get_text() not in ("x = logit NDI", "y = logit FISO − ln NDI"):
            if re.search(r"\d", s):
                unregistered.append(s)
        texts.append((t, bb, REGISTERED.get(id(t), (s, "axis"))[0]))
    if small:
        stop(f"texts below {MIN_PT} pt: {small}")
    if outside:
        stop(f"texts outside the figure: {outside}")
    if unregistered:
        stop(f"figure texts with numbers that no ledger row covers: {unregistered}")
    overlaps = []
    for (t1, b1, g1), (t2, b2, g2) in itertools.combinations(texts, 2):
        if g1 == g2 and REGISTERED.get(id(t1), ("", ""))[1] == "symbol":
            continue                                                # pieces of one symbol label
        if b1.x0 < b2.x1 - 0.5 and b2.x0 < b1.x1 - 0.5 and b1.y0 < b2.y1 - 0.5 and b2.y0 < b1.y1 - 0.5:
            overlaps.append((t1.get_text(), t2.get_text()))
    if overlaps:
        stop(f"overlapping texts: {overlaps}")
    return len(texts)


# Template of the direction's caption draft (S39 §4) with ledger IDs. Two wordings differ from the draft because the caption
# must pass ledger_check --strict-all: «T2 82.5, 52 and 1000 ms» → «T2 of 82.5 ms, 52 ms and 1000 ms» (with spaces removed,
# «T2 82.5» reads as one number), and the number words without ledger rows are dropped («Nine tissue states» → «Tissue
# states»; «between two relaxation directions» → the shaded cone it means, which also gives the shading its only legend).
# v12 (review direction S41 §2–3): caption compressed; figure renumbered 1 -> 2.
CAPTION = ("![**Fig. 3. Relaxation weighting in log-ratio coordinates.** (a) Tissue states (black) and the same states at "
           "TE + {N-063} (blue) move by one vector, ({N-064}, {N-065}): a translation; axes at equal scale. (b) Directions of "
           "TE (blue) and TR (orange) differences between sites and of tissue changes (green), on the relaxation grid of "
           "Section 2.2. (c) Directions of relaxation differences between groups (vermillion)."
           "](figures/F2_logratio.pdf){#fig:fig3 .wide width=100%}")


def caption(V):
    """Fill the template; then check the caption's statements on the values and run ledger_check on it."""
    dirs, te_w, tr_f = V["dirs"], V["te_wedge"], V["tr_fan"]
    # «More free water points along a TR change and a tissue-T1 lengthening»: 90 deg inside the TR fan, equal to tissue T1
    if not (tr_f[0] <= dirs["free water"] <= tr_f[1]) or dirs["free water"] != dirs["tissue T1 up"]:
        stop("caption: free water is not along the TR fan and the tissue-T1 direction")
    # «neurite loss lies outside both sectors but between two relaxation directions»
    nl = dirs["neurite loss"] % 180
    if te_w[0] <= nl <= te_w[1] or tr_f[0] <= nl <= tr_f[1]:
        stop("caption: neurite loss lies inside a sector")
    if not (dirs["R2in up"] < dirs["neurite loss"] < dirs["T2en up"]):
        stop("caption: neurite loss is not between R2,in up and T2,en up")
    s = re.sub(r"\{(N-\d{3})\}", lambda m: f"{row(m.group(1))['원고 표기']} [{m.group(1)}]", CAPTION)
    md = OUT / "F2_caption.md"
    md.write_text(s + "\n", encoding="utf-8")
    res = subprocess.run([sys.executable, "-I", "-B", str(LEDGER_DIR / "ledger_check.py"), str(md), "--strict-all", "1"],
                         capture_output=True, text=True, env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
    first = res.stdout.splitlines()[0] if res.stdout else res.stderr[-400:]
    if not re.search(r"실패 0\b", first) or not re.search(r"쓰지 않을 말 0\b", first):
        stop(f"caption ledger_check: {first}\n{res.stdout[-1500:]}")
    words = len([w for w in re.sub(r"\[N-\d{3}\]|\]\(figures/.*$|!\[|\*\*", " ", s).split() if re.search(r"\w", w)])
    return s, first, words


def main():
    V = read_values()
    fig = draw(V)
    n_texts = audit(fig)
    cap, ledger_line, words = caption(V)
    meta = {"CreationDate": None, "Producer": None}
    fig.savefig(OUT / f"{NAME}.pdf", metadata=meta)
    fig.savefig(OUT / f"{NAME}.png", dpi=300, metadata={"Software": None})
    lines = ["# Figure 1 — values, guards and audit", "",
             f"- Recomputed on the relaxation grid (config.py copy {row('N-628')['근거 파일 sha256'][:12]}): TE wedge "
             f"{V['te_wedge'][0]!r}–{V['te_wedge'][1]!r}° (N-049/N-050), TR fan {V['tr_fan'][0]!r}–{V['tr_fan'][1]!r}° "
             f"(N-059/N-060); printed as N-206 and N-207.",
             f"- Grid: T2_in {list(V['grid']['T2IN_GRID'])} ms, r {sorted({r for _, r, _ in V['grid']['ETA_T2']})}, T2_iso "
             f"{sorted({i for *_, i in V['grid']['ETA_T2']})} ms; T1 sets {V['grid']['T1_SETS']} (geometric mean, T1_in/T1_en), "
             f"T1_iso {V['grid']['T1_ISO']} ms; acquisitions A {V['grid']['ACQ_E5']['A']}, Q1 {V['grid']['ACQ_E5']['Q1']}, "
             f"Q4 {V['grid']['ACQ_E5']['Q4']}.",
             f"- TE step vector {V['vec'].tolist()} (N-064/N-065), direction {V['vdeg']!r}° (N-602).",
             f"- Directions: " + ", ".join(f"{k} {v!r}°" for k, v in V['dirs'].items()) + ".",
             f"- Figure texts audited: {n_texts}; all ≥ {MIN_PT} pt, inside the figure, no overlap.",
             f"- Caption: {words} words; ledger_check --strict-all 1: {ledger_line}", "",
             "| ID | raw value | manuscript text | role |", "|---|---|---|---|"]
    seen = set()
    for i, raw, ms, role in CHECKED:
        if (i, role) in seen:
            continue
        seen.add((i, role))
        lines.append(f"| {i} | {raw} | {ms} | {role} |")
    lines += ["", "## Visible texts with numbers", "", "| where | text | IDs |", "|---|---|---|"]
    lines += [f"| {w} | {t} | {' · '.join(ids)} |" for w, t, ids in VISIBLE]
    (OUT / "F2_values.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"F1 written to {OUT}: {NAME}.pdf/.png, F2_caption.md ({words} words), F2_values.md; {ledger_line}")


if __name__ == "__main__":
    main()
