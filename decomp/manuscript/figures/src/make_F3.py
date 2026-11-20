"""make_F3 draws Figure 2 (Section 2.4): a site TE x group T2 interaction, from the signal model and ledger values only, computing no statistic."""
import sys

sys.dont_write_bytecode = True        # first statement: no __pycache__ in the repository's figures/src or the ledger folder

import csv  # noqa: E402
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
NAME = "F3_site_TE_group_T2"

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.patches import FancyArrowPatch  # noqa: E402
from matplotlib.text import Text  # noqa: E402

_spec = importlib.util.spec_from_file_location("ledger_check", LEDGER_DIR / "ledger_check.py")
LC = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(LC)          # read-only use: eval_field · same_value · rounding_ok · numbers_in

LEDGER = {r["ID"]: r for r in csv.DictReader((LEDGER_DIR / "numbers.csv").open(encoding="utf-8"))}
CHECKED = []                          # (ID, raw value, manuscript text, role) for F3_values.md
VISIBLE = []                          # (where, text, IDs)
REGISTERED = {}                       # id(Text) -> (where, kind) for the figure-text audit
MIN_PT = 7.0


def stop(msg):
    sys.exit(f"STOP F2: {msg}")


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


def val(id_, role):
    """The row reproduces from its own copy (literature rows: not machine-checked) and its rounding rule holds."""
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
    CHECKED.append((id_, r["원 값"], r["원고 표기"], role + ("" if v is not None else " — literature row, not recomputed")))
    return float(r["원 값"].replace("−", "-"))


def agree(id_, computed, k, rel=1e-9):
    raw = float(row(id_)["원 값"].replace("−", "-"))
    if not math.isclose(computed, raw, rel_tol=rel, abs_tol=1e-15):
        stop(f"{id_}: recomputed {computed!r}, ledger raw {raw!r}")
    if fmt(computed, k) != fmt(raw, k):
        stop(f"{id_}: recomputed prints {fmt(computed, k)}, ledger {fmt(raw, k)}")


# ============================================================================ values
def nu_app(TE, nu, T2in, T2en):
    w_in, w_en = np.exp(-TE / T2in), np.exp(-TE / T2en)
    return nu * w_in / (nu * w_in + (1 - nu) * w_en)


def logit(p):
    return np.log(p / (1 - p))


def logit_closed(TE, nu, T2in, T2en):
    return logit(nu) + TE * (1 / T2en - 1 / T2in)


def read_values():
    nu_c = val("N-093", "true NDI of controls and the relaxation-only group")
    T2in, T2en = val("N-644", "intra-neurite T2, all groups"), val("N-043", "extra-neurite T2, controls")
    T2en_r = val("N-099", "extra-neurite T2, relaxation-only group")
    nu_f = val("N-107", "true NDI of the fraction-only group")
    TE_B, TE_A = val("N-100", "TE at site B"), val("N-101", "TE at site A")
    groups = {"Controls": (nu_c, T2in, T2en), "Relaxation only": (nu_c, T2in, T2en_r), "Fraction only": (nu_f, T2in, T2en)}
    c, r, f = (groups[k] for k in ("Controls", "Relaxation only", "Fraction only"))
    dB = float(nu_app(TE_B, *r) - nu_app(TE_B, *c))
    dA = float(nu_app(TE_A, *r) - nu_app(TE_A, *c))
    val("N-649", "raw difference at site B (magnitude)")
    val("N-650", "raw difference at site A (magnitude)")
    agree("N-649", -dB, 3)
    agree("N-650", -dA, 3)
    lgB = float(logit(nu_app(TE_B, *r)) - logit(nu_app(TE_B, *c)))
    lgA = float(logit(nu_app(TE_A, *r)) - logit(nu_app(TE_A, *c)))
    inter = lgA - lgB
    closed = (TE_A - TE_B) * ((1 / r[2] - 1 / r[1]) - (1 / c[2] - 1 / c[1]))
    if abs(inter - closed) > 1e-12:
        stop(f"logit change {inter!r} differs from the closed form {closed!r}")
    val("N-106", "logit change of the difference between the sites")
    agree("N-106", inter, 3)
    gaps = [float(logit(nu_app(t, *f)) - logit(nu_app(t, *c))) for t in (0.0, TE_B, TE_A, 120.0)]
    if not np.allclose(gaps, gaps[0], rtol=0, atol=1e-12):
        stop(f"fraction-only gap is not constant in TE: {gaps}")
    te = np.linspace(0, 120, 13)
    for k, g in groups.items():
        if not np.allclose(logit(nu_app(te, *g)), logit_closed(te, *g), rtol=0, atol=1e-12):
            stop(f"{k}: logit line is not straight in TE")
    if abs(logit_closed(0.0, *c) - logit_closed(0.0, *r)) > 1e-15:
        stop("controls and relaxation-only lines do not meet at TE 0")
    return dict(groups=groups, TE_B=TE_B, TE_A=TE_A, dB=dB, dA=dA, lgB=lgB, lgA=lgA, inter=inter, gaps=gaps)


# ============================================================================ drawing
FONT = "Liberation Sans"
try:
    font_manager.findfont(font_manager.FontProperties(family=FONT), fallback_to_default=False)
except ValueError:
    stop(f"font {FONT} is not installed")
OI = dict(blue="#0072B2", green="#009E73", verm="#D55E00", grey="#6E6E6E", ink="#1A1A1A")   # Okabe–Ito, fixed to meaning
STYLE = {"Controls": dict(color=OI["ink"], ls="-"), "Relaxation only": dict(color=OI["verm"], ls=(0, (4, 1.6))),
         "Fraction only": dict(color=OI["green"], ls="-")}
plt.rcParams.update({
    "font.family": FONT, "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5,
    "ytick.major.size": 2.5, "pdf.fonttype": 42, "ps.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False,
})
MM = 1 / 25.4
W, H = 140, 64                                                     # mm, 1.5 column


def register(t, s, ids=(), where="", kind="text"):
    REGISTERED[id(t)] = (where or s, kind)
    if kind == "text" and (ids or re.search(r"\d", s)):
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


def draw(V):
    fig = plt.figure(figsize=(W * MM, H * MM))
    groups, TE_B, TE_A = V["groups"], V["TE_B"], V["TE_A"]
    c, r, f = (groups[k] for k in ("Controls", "Relaxation only", "Fraction only"))
    TE = np.linspace(0, 120, 241)
    sub = {"Controls": (f"NDI {row('N-093')['원고 표기']}", ("N-093",)),
           "Relaxation only": (f"extra-neurite\nT2 {row('N-099')['원고 표기']}", ("N-099",)),
           "Fraction only": (f"NDI {row('N-107')['원고 표기']}", ("N-107",))}

    def ax_mm(l, b, w, h):
        return fig.add_axes([l / W, b / H, w / W, h / H])

    def panel(ax, fn, ylim, yticks, ylabel, letter, title, x0):
        ax.set_xlim(0, 120)
        ax.set_ylim(*ylim)
        ax.set_xticks([0, 40, 80, 120])
        ax.set_yticks(yticks)
        ax.set_yticklabels([("−" if t < 0 else "") + f"{abs(t):g}" for t in yticks])
        ax.set_xlabel("TE (ms)", labelpad=2)
        ax.set_ylabel(ylabel, labelpad=2)
        for te, s, id_ in ((TE_B, "site B", "N-100"), (TE_A, "site A", "N-101")):
            ax.axvline(te, color="#C8C8C8", lw=0.6, zorder=0)
            lab = f"{s}\n{row(id_)['원고 표기'].replace(' ms', '')} ms"
            register(ax.text(te, ylim[1], lab, ha="center", va="bottom", color=OI["grey"], linespacing=1.05), lab, (id_,))
        for k, g in groups.items():
            ax.plot(TE, fn(TE, *g), lw=1.2, **STYLE[k], zorder=2)
            ax.plot([TE_B, TE_A], fn(np.array([TE_B, TE_A]), *g), "o", ms=3.4, mfc=STYLE[k]["color"], mec="white", mew=0.5,
                    zorder=4)
        for s, kw, dx in ((letter, dict(fontsize=8, fontweight="bold"), 0.0), (title, dict(fontsize=7.5), 4.2)):
            register(fig.text((x0 + dx) / W, (H - 4) / H, s, va="baseline", **kw), s, kind="title")

    def end_labels(ax, fn, short=()):
        for k, g in groups.items():
            s, ids = (k, ()) if k in short else (f"{k}\n{sub[k][0]}", sub[k][1])
            register(ax.annotate(s, xy=(120, fn(120.0, *g)), xytext=(3, 0), textcoords="offset points", ha="left",
                                 va="center", color=STYLE[k]["color"], linespacing=1.1, annotation_clip=False), s, ids)

    def gap_arrow(ax, te, y0, y1, color):
        ax.add_patch(FancyArrowPatch((te, y0), (te, y1), arrowstyle="-|>", mutation_scale=5, lw=0.8, color=color,
                                     shrinkA=2.2, shrinkB=2.2, zorder=3))

    # ------------------------------------------------------------ (a) raw scale
    AW, AH = 38.5, 44.0
    axa = ax_mm(12, 10, AW, AH)
    panel(axa, nu_app, (0.38, 0.70), [0.4, 0.5, 0.6, 0.7], "Apparent NDI", "(a)", "Raw scale: the gap grows with TE", 1.5)
    for te, id_ in ((TE_B, "N-649"), (TE_A, "N-650")):
        y0, y1 = nu_app(te, *c), nu_app(te, *r)
        gap_arrow(axa, te, y0, y1, OI["verm"])
        s = f"−{row(id_)['원고 표기']}"
        register(axa.text(te, y1 - 0.016, s, ha="center", va="top", color=OI["verm"]), s, (id_,))
    end_labels(axa, nu_app)

    # ------------------------------------------------------------ (b) logit scale
    axb = ax_mm(82, 10, AW, AH)
    panel(axb, logit_closed, (-0.5, 0.8), [-0.4, 0, 0.4, 0.8], "logit apparent NDI", "(b)", "Logit scale: straight lines",
          71.5)
    y0 = logit_closed(0, *c)
    axb.plot(0, y0, "o", ms=5.5, mfc="none", mec=OI["ink"], mew=0.7, zorder=5, clip_on=False)
    s = "same tissue:\nmeet at TE 0"
    register(axb.annotate(s, xy=(1.5, y0 + 0.03), xytext=(6, 0.36), ha="left", va="bottom", color=OI["ink"],
                          linespacing=1.05,
                          arrowprops=dict(arrowstyle="-", lw=0.5, color=OI["grey"], shrinkA=1, shrinkB=2)),
             s, kind="tick", where="axis origin TE 0 (exempt, like a tick label)")
    s = "parallel:\nconstant gap"
    register(axb.text(18, logit_closed(18, *f) - 0.05, s, ha="left", va="top", color=OI["green"], linespacing=1.05), s)
    for te in (TE_B, TE_A):
        gap_arrow(axb, te, logit_closed(te, *c), logit_closed(te, *r), OI["verm"])
    ghost = logit_closed(TE_A, *c) + V["lgB"]                       # site-B gap hung from the site-A control point
    xb = TE_A + 5
    axb.plot([xb, xb], [logit_closed(TE_A, *c), ghost], color=OI["verm"], lw=2.2, alpha=0.3, solid_capstyle="butt",
             zorder=3)
    axb.plot([xb, xb], [ghost, logit_closed(TE_A, *r)], color=OI["verm"], lw=2.2, solid_capstyle="butt", zorder=3)
    s = f"{row('N-106')['원고 표기']}\nvs site B"
    register(axb.annotate(s, xy=(xb + 1, (ghost + logit_closed(TE_A, *r)) / 2), xytext=(122, 0.25), ha="left", va="center",
                          color=OI["verm"], linespacing=1.05, annotation_clip=False,
                          arrowprops=dict(arrowstyle="-", lw=0.5, color=OI["verm"], shrinkA=1, shrinkB=0)), s, ("N-106",))
    # (b): the relaxation-only line is named at its end in one line (its definition stands in (a) and in the caption), so no
    # label box covers the site-B point or the arrows inside the panel
    end_labels(axb, logit_closed, short=("Relaxation only",))
    return fig


# ============================================================================ audits
def audit(fig):
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    fb = fig.bbox
    ticks = {id(t) for ax in fig.axes for t in ax.get_xticklabels() + ax.get_yticklabels()}
    axis_labels = {"TE (ms)", "Apparent NDI", "logit apparent NDI"}
    texts, small, outside, unregistered = [], [], [], []
    for t in fig.findobj(Text):
        s = t.get_text()
        if not s.strip() or not t.get_visible():
            continue
        if t.get_fontsize() < MIN_PT - 1e-9:
            small.append(f"{s} ({t.get_fontsize()} pt)")
        if id(t) in ticks:
            continue
        bb = Text.get_window_extent(t, rend)                         # the text alone (an annotation's arrow excluded)
        if bb.x0 < fb.x0 - 0.5 or bb.x1 > fb.x1 + 0.5 or bb.y0 < fb.y0 - 0.5 or bb.y1 > fb.y1 + 0.5:
            outside.append(s)
        if id(t) not in REGISTERED and s not in axis_labels and re.search(r"\d", s):
            unregistered.append(s)
        texts.append((t, bb))
    if small:
        stop(f"texts below {MIN_PT} pt: {small}")
    if outside:
        stop(f"texts outside the figure: {outside}")
    if unregistered:
        stop(f"figure texts with numbers that no ledger row covers: {unregistered}")
    overlaps = [(a.get_text(), b.get_text()) for (a, ba), (b, bb) in itertools.combinations(texts, 2)
                if ba.x0 < bb.x1 - 0.5 and bb.x0 < ba.x1 - 0.5 and ba.y0 < bb.y1 - 0.5 and bb.y0 < ba.y1 - 0.5]
    if overlaps:
        stop(f"overlapping texts: {overlaps}")
    return len(texts)


# Template of the direction's caption draft (S39 §4) with ledger IDs; «T2 82.5 and 52 ms» and «extra-neurite T2 70 ms» take
# «of» because ledger_check reads «T2 82.5» without spaces as one number. The title keeps its edit marker [.]{.aedit}.
# v12 (review direction S41 §2–3): caption compressed; figure renumbered 2 -> 3.
# v16p (S74 direction §4): renumbered 3 -> 2 (first cited in the Introduction); caption ends «(illustration; derivation in Section 2.4)».
CAPTION = ("![**Fig. 2. A site TE × group T2 interaction (illustration).** Apparent NDI of controls, a group with longer "
           "extra-neurite T2 (relaxation only) and a group with lower NDI (fraction only), at sites with TE {N-100} and "
           "{N-101}. (a) Raw scale: the relaxation-only gap widens with TE. (b) Logit scale: lines are straight; the "
           "relaxation-only gap changes by {N-106} between sites, the fraction-only gap does not (illustration; derivation in Section 2.4)."
           "](figures/F3_site_TE_group_T2.pdf){#fig:fig2 width=100%}")


def caption():
    s = re.sub(r"\{(N-\d{3})\}", lambda m: f"{row(m.group(1))['원고 표기']} [{m.group(1)}]", CAPTION)
    md = OUT / "F3_caption.md"
    md.write_text(s + "\n", encoding="utf-8")
    res = subprocess.run([sys.executable, "-I", "-B", str(LEDGER_DIR / "ledger_check.py"), str(md), "--strict-all", "1"],
                         capture_output=True, text=True, env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
    first = res.stdout.splitlines()[0] if res.stdout else res.stderr[-400:]
    # The only number allowed without a ledger ID is the axis origin in «meets the control line at TE 0» (like a tick label).
    m = re.search(r"실패 (\d+)\b", first)
    fails = [ln.strip() for ln in res.stdout.splitlines()[1:] if ln.strip().startswith("실패")]
    origin = [ln for ln in fails if re.fullmatch(r"실패 원고 \d+ 행: ID 없는 숫자 «0» \(엄격 구역\)", ln)]
    if (not m or int(m.group(1)) != len(origin) or len(origin) != s.count("at TE 0")
            or not re.search(r"쓰지 않을 말 0\b", first)):
        stop(f"caption ledger_check: {first}\n{res.stdout[-1500:]}")
    first += f" — the {len(origin)} failure(s) are the axis origin «at TE 0», allowed as a tick-like value"
    words = len([w for w in re.sub(r"\[N-\d{3}\]|\]\(figures/.*$|!\[|\*\*|\[\.\]\{\.aedit\}", " ", s).split()
                 if re.search(r"\w", w)])
    return s, first, words


def main():
    V = read_values()
    fig = draw(V)
    n_texts = audit(fig)
    cap, ledger_line, words = caption()
    fig.savefig(OUT / f"{NAME}.pdf", metadata={"CreationDate": None, "Producer": None})
    fig.savefig(OUT / f"{NAME}.png", dpi=300, metadata={"Software": None})
    lines = ["# Figure 2 — values, guards and audit", "",
             f"- Raw-scale differences (relaxation-only − control): site B {V['dB']!r}, site A {V['dA']!r} (N-649 · N-650 are the "
             f"magnitudes); logit differences B {V['lgB']!r}, A {V['lgA']!r}; change {V['inter']!r} (N-106; equals the closed "
             f"form within 1e-12).",
             f"- Fraction-only − control logit gap at TE 0, 70, 100, 120 ms: {V['gaps']!r} (constant).",
             f"- Figure texts audited: {n_texts}; all ≥ {MIN_PT} pt, inside the figure, no overlap.",
             f"- Caption: {words} words; ledger_check --strict-all 1: {ledger_line}", "",
             "| ID | raw value | manuscript text | role |", "|---|---|---|---|"]
    lines += [f"| {i} | {raw} | {ms} | {role} |" for i, raw, ms, role in CHECKED]
    lines += ["", "## Visible texts with numbers", "", "| where | text | IDs |", "|---|---|---|"]
    lines += [f"| {w} | {t} | {' · '.join(ids)} |" for w, t, ids in VISIBLE]
    (OUT / "F3_values.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"F2 written to {OUT}: {NAME}.pdf/.png, F3_caption.md ({words} words), F3_values.md; {ledger_line}")


if __name__ == "__main__":
    main()
