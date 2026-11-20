"""Data figures of the NODDI manuscript — draw functions v3 and the one template left (Fig. 4).

Status 2026-11-01: Fig. 3, Fig. 5 and Fig. S2 are drawn from stored values by fig_real_data.py with
the draw functions below (template mark off). Fig. 4 waits for the V4 judgement and is the only
figure still drawn here from synthetic data, with the template mark. Fig. 5 is drawn from the
stored Tong judge output only (tong_judge.json keeps counts and shares, not per-pair effects or
per-unit shares), so its draw function takes that structure; synth_fig5 was removed with it.
The text below describes the v3 templates as they were made.

Fig. 3  S4 identifiability map            (was F4)
Fig. 4  V4 calibration on travelling heads (was F6)
Fig. 5  Tong negative control             (was F7)
Fig. S2 Wongkornchaovalit multi-TE series  (was F8)

v3 ports fig_templates_v2.py (sha256 14a23f9a...) to the shared rule set in figstyle.py
(text >= 7 pt, no sub/superscripts, widths 90/140/190 mm, titles and explanations in the
captions, one colour meaning across figures). Definitions and thresholds are unchanged.

Structure. Each figure has one draw function that takes a data dict, and one synthetic
generator that builds that dict from a fixed seed:

    draw_fig3(d, template)   <- synth_fig3(rng)
    draw_fig4(d, template)   <- synth_fig4(rng)
    draw_fig5(d, template)   <- fig_real_data.py fig5 (stored Tong judge output; no synthetic generator)
    draw_figS1(d, template)  <- synth_figS1(rng)

The template flag (keyword, no default) adds the grey "TEMPLATE · synthetic data" mark to every
data panel; the final figures are drawn by the same draw functions from repository outputs with
the flag off.
Nothing here reads from disk or the network; no image value and no S4, V4, Tong or
Wongkornchaovalit output is used. Design-level facts are real: the 28 ON-Harmony card pairs,
their orientation, dTE, dTR and scheme difference (CARDS, copied from fixed document v10,
appendix, cards 2b6aa13 — acquisition metadata), the five acquisition types they fall into
by the S4/V4 code's same-TR rule (1 / 3 / 0 / 9 / 15, asserted; see SAME_TR) and the five V4
cells without data (NO_DATA, fixed document v10 a11, asserted). Outcomes are synthetic: C and
MDE (Fig. 3), coverage, thresholds, verdicts, flags, subject and unit counts and attribute
shares (Fig. 4), and everything in Fig. 5 and S2.

Data dicts and the repository outputs that will feed them (field names marked ⟦R⟧ are not
known to the reviewer and must be read from the files):

Fig. 3 — cards_c3.json (main) and cards_c7.json (supplement), reference condition L0 · SNR 30 · Neff 50
    pairs        list of 28 (acq_a, acq_b)                         ⟦R: card pair id field⟧
    ref_pair     ("OXF1PRI", "OXF2PRI")  same TE, TR and scheme; coil differs
    stratum      {pair: 0..4} = the five validation strata (matched · TE only · TR only ·
                 TE and TR · TE, TR and scheme)                     ⟦R: acquisition-type field⟧
    C            {comp: {pair: C}}  shrinkage, worst of nine k, C = 1 - (sigma_post / sigma_prior)^2
                                                                    ⟦R: per-card C or sigma_post field⟧
    ref_cell     (SNR 30, Neff 50); corner_cells: SNR {20, 50} x Neff {25, 100}, L0 (the only
                 cells with MDE in the cards, decomp/s4/plan.py:86)  ⟦R: cell key field⟧
    mde_ref      {comp: {pair: MDE at ref_cell}}                    MDE = 2.80 sigma_post
    mde_corner   {comp: {pair: [MDE at each corner cell]}}          ⟦R: per-card MDE field⟧
    sigma_prior  {comp: width / sqrt(12)}; reporting prior uniform +-0.10 (NDI, ODI), +-0.05 (FISO)
  Classes: identified C >= 0.75, not identified C < 0.25 (Methods 3.6.4).

Fig. 4 — V4 judge output (judge_v4.py run <arm C3> <attribute-table sha256>) ⟦R: output file⟧
    cells        24 dicts (NDI 1 + ODI 23 cells with coverage), each:
                 comp, pair, design (dTE ms, dTR s, scheme — from CARDS), units (<= 27 x subjects),
                 subjects (common subjects of the card pair)                ⟦R: field names⟧
                 k_cov  (9 coverages, one per k), coverage (median over the nine k),
                 threshold (own threshold = 5th percentile in S4 test data), fallback (True if 0.80
                 was used because the percentile cannot be computed), verdict (coverage >= threshold),
                 gap (misspecification exceedance flagged "gap detected"), c7_differs (C7 arm disagrees),
                 decomp {"noise"|"range"|"vendor"|"none": share of units not covered, or None if no such
                 session pair}  (fixed document v10 sec. 2-5-1, R20)  ⟦R: field names⟧
    nodata       5 dicts (pair, design, units 0, subjects 0): ODI cells with no subject scanned on
                 both acquisitions = NO_DATA (asserted)

Fig. 5 — Tong T-A / T-B / T-C output (fixed document v10 sec. 2-3, R10): tong_judge.json main
    n_pairs, n_within  centre pairs, and pairs whose effect length is within the margin (T_A)
    margin, limit      1.96 * sqrt(2) * max(s_x, s_y), and the count from which T-A is undeterminable
    d_share            share of axial directions in D (T_A.p_random)
    tb                 {coord: (n_retest, n_total, share)} over all units (T_B)
    tc                 (aligned, retest) cell reached, or None if undeterminable (T_C)
  Only the undeterminable branch of T-A can be drawn from stored values (asserted); per-pair effects
  and per-unit shares are not in the judge output.
  D = TE wedge 48.5-72.1 deg U TR fan 85.3-94.3 deg; expectation under random direction 0.181.
  T-B threshold 0.90 per coordinate; T-C retest cell needs x and y >= 0.90.

Fig. S2 — Wongkornchaovalit series (fixed document v10 sec. 2-6) ⟦R: output file⟧
    TE           eight TEs 62-132 ms; TE_rep repeats at 62 and 132 ms (retest only, not fitted)
    units        {"x"|"y": list of 81 dicts (3 subjects x 27 ROI units): subject, v (8 values),
                 rep (2 values), lin (slope, intercept of the linear fit centred at 97 ms),
                 ratio = |c| * 35^2 / (2 * sqrt(2) * s_rt) (curvature; linear if <= 1),
                 icpt (intercept at TE -> 0 from low TEs, high TEs, all eight), single (min, max)}}
  Share threshold 0.8.

Run: python3 -I fig_data_templates_v3.py [fig3 fig4 fig5 figS1]   (CPU only)
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs  # noqa: E402

import matplotlib  # noqa: E402
import matplotlib.ticker as mticker  # noqa: E402
from matplotlib.patches import Arc, Rectangle, Wedge  # noqa: E402

INK, INK2, MUTED, GRID, AXIS, SURF = fs.INK, fs.INK2, fs.MUTED, fs.GRID, fs.AXIS, fs.SURF
COMPONENTS, CCOL, CMRK, CMS = fs.COMPONENTS, fs.CCOL, fs.CMRK, fs.CMS
PT = fs.PT_TEXT
SEED = 20261027
ZEBRA = "#f6f6f3"
CELL_OFF, CELL_ON = "#f6f5f2", "#e6e5df"

# Design of the 28 ON-Harmony card pairs — acquisition metadata, not image values.
# Copied from fixed analysis document v10, appendix "ON-Harmony 28 장의 획득 차이" (acquisition differences of the 28 ON-Harmony cards) (cards 2b6aa13),
# FIXED_DRAFT_v10_2026-10-27.md lines 588-615: card orientation A~B, dTE (ms, 0.1), dTR (ms),
# direction scheme differs. Only coverage, verdicts, subject/unit counts and attribute shares
# below are synthetic.
CARDS = [
    ("OXF1PRI", "OXF2PRI", 0.0, 0.0, False),
    ("OXF1PRI", "OXF3TRI", 4.4, 0.0, False),
    ("OXF1PRI", "NOT2ING-89", -3.0, 417.5, True),
    ("OXF1PRI", "NOT2ING-98", 6.0, 604.5, True),
    ("OXF1PRI", "NOT4GEP", -24.0, 1300.0, True),
    ("OXF1PRI", "OXF4GEP", -22.3, 4400.0, True),
    ("OXF1PRI", "NOT1ACH", -21.4, 6615.0, True),
    ("OXF2PRI", "OXF3TRI", 4.4, 0.0, False),
    ("OXF2PRI", "NOT2ING-89", -3.0, 417.5, True),
    ("OXF2PRI", "NOT2ING-98", 6.0, 604.5, True),
    ("OXF2PRI", "NOT4GEP", -24.0, 1300.0, True),
    ("OXF2PRI", "OXF4GEP", -22.3, 4400.0, True),
    ("OXF2PRI", "NOT1ACH", -21.4, 6615.0, True),
    ("OXF3TRI", "NOT2ING-89", -7.4, 417.5, True),
    ("OXF3TRI", "NOT2ING-98", 1.6, 604.5, True),
    ("OXF3TRI", "NOT4GEP", -28.4, 1300.0, True),
    ("OXF3TRI", "OXF4GEP", -26.8, 4400.0, True),
    ("OXF3TRI", "NOT1ACH", -25.8, 6615.0, True),
    ("NOT2ING-89", "NOT2ING-98", 9.0, 187.0, False),
    ("NOT2ING-89", "NOT4GEP", -21.0, 882.5, False),
    ("NOT2ING-89", "OXF4GEP", -19.3, 3982.5, False),
    ("NOT2ING-89", "NOT1ACH", -18.4, 6197.5, False),
    ("NOT2ING-98", "NOT4GEP", -30.0, 695.5, False),
    ("NOT2ING-98", "OXF4GEP", -28.3, 3795.5, False),
    ("NOT2ING-98", "NOT1ACH", -27.4, 6010.5, False),
    ("NOT4GEP", "OXF4GEP", 1.7, 3100.0, False),
    ("NOT4GEP", "NOT1ACH", 2.6, 5315.0, False),
    ("OXF4GEP", "NOT1ACH", 0.9, 2215.0, False),
]
DESIGN = {(a, b): (dte, dtr, sch) for a, b, dte, dtr, sch in CARDS}
PAIRS = [(a, b) for a, b, *_ in CARDS]
REF_PAIR = ("OXF1PRI", "OXF2PRI")          # same TE, TR and scheme; receive coil differs
# V4 cells with no subject scanned on both acquisitions (fixed document v10 a11, sec. 2-5, R11)
NO_DATA = [("NOT2ING-89", "NOT2ING-98"), ("NOT2ING-89", "NOT4GEP"), ("NOT2ING-89", "OXF4GEP"),
           ("OXF3TRI", "NOT4GEP"), ("OXF3TRI", "OXF4GEP")]
STRATA = ["Matched", "TE only", "TR only", "TE and TR", "TE, TR and\nscheme"]
# Subjects per acquisition (Methods 3.7.2; fixed document v10 line 140): NOT2ING was run with
# TE 89 ms for 4 subjects and TE 98 ms for 16; every other acquisition has at most 20.
ACQ_SIZE = {a: 20 for p in PAIRS for a in p} | {"NOT2ING-89": 4, "NOT2ING-98": 16}


def pair_cap(p):
    return min(ACQ_SIZE[p[0]], ACQ_SIZE[p[1]])
VENDOR = {"ACH": "Philips", "ING": "Philips", "PRI": "Siemens", "TRI": "Siemens", "GEP": "GE"}

TE_WEDGE = (48.5, 72.1)
TR_FAN = (85.3, 94.3)
D_SHARE = ((TE_WEDGE[1] - TE_WEDGE[0]) + (TR_FAN[1] - TR_FAN[0])) / 180.0
assert len(PAIRS) == 28 and len(set(PAIRS)) == 28 and round(D_SHARE, 3) == 0.181
assert len({x for p in PAIRS for x in p}) == 8 and all(p in DESIGN for p in NO_DATA)

C_ID, C_NOT = 0.75, 0.25
MDE_K = 2.80
PRIOR_HALF = {"NDI": 0.10, "FISO": 0.05, "ODI": 0.10}
SIG_PRIOR = {c: 2 * h / math.sqrt(12) for c, h in PRIOR_HALF.items()}


def canvas(fig):
    """fs.canvas with its (hidden) ticks cleared, so that fs.audit does not count them."""
    c = fs.canvas(fig)
    c.set_xticks([])
    c.set_yticks([])
    return c


# Acquisition type follows the S4/V4 code (decomp/remeasure/judge_v4.py:94-101 acq_type): TR differs
# only when the card TR ranges (s, decomp/s4/plan.py:17-22 ONH) have different sums. NOT2ING-89 and
# NOT2ING-98 share (4.00, 4.38) s, so that pair is TE only, although DESIGN keeps dTR +187 ms from the
# confirmed-TR card table for display. It is the only pair where the two disagree (asserted).
CARD_TR = {"OXF1PRI": (3.60, 3.60), "OXF2PRI": (3.60, 3.60), "OXF3TRI": (3.60, 3.60),
           "NOT2ING-89": (4.00, 4.38), "NOT2ING-98": (4.00, 4.38), "NOT4GEP": (4.90, 4.90),
           "OXF4GEP": (8.00, 8.00), "NOT1ACH": (9.81, 10.62)}
SAME_TR = {p for p in PAIRS if abs(sum(CARD_TR[p[0]]) - sum(CARD_TR[p[1]])) <= 1e-9}
assert {p for p in PAIRS if (DESIGN[p][1] == 0) != (p in SAME_TR)} == {("NOT2ING-89", "NOT2ING-98")}


def acq_type(p):
    dte, _, sch = DESIGN[p]
    dtr = p not in SAME_TR                  # code rule, not the displayed dTR
    if sch:
        assert dte != 0 and dtr, p
        return 4
    return {(False, False): 0, (True, False): 1, (False, True): 2, (True, True): 3}[(dte != 0, dtr)]


N_TYPE = [sum(acq_type(p) == k for p in PAIRS) for k in range(5)]
assert N_TYPE == [1, 3, 0, 9, 15], N_TYPE
assert [p for p in PAIRS if acq_type(p) == 0] == [REF_PAIR]
assert [p for p in PAIRS if acq_type(p) == 1] == [("OXF1PRI", "OXF3TRI"), ("OXF2PRI", "OXF3TRI"),
                                                  ("NOT2ING-89", "NOT2ING-98")]


def fmt_signed(v, nd=1):
    return f"{v:+.{nd}f}".replace("-", "\u2212")


def design_cells(p):
    """(dTE ms text, dTR s text, scheme text) in the card orientation."""
    dte, dtr, sch = DESIGN[p]
    return (fmt_signed(dte) if dte else "0", fmt_signed(dtr / 1000.0) if dtr else "0",
            "yes" if sch else "no")


def pname(p, sep=" ~ "):
    return f"{p[0]}{sep}{p[1]}"


def vendor(a):
    return VENDOR[a[4:7]]


def legend_marker(cv, x, y, comp=None, filled=True, color=None, ms=None, marker=None):
    col = color or CCOL[comp]
    cv.plot([x], [y], ls="none", marker=marker or CMRK[comp], ms=ms or CMS[comp],
            mfc=col if filled else SURF, mec=col if not filled else SURF, mew=0.9 if not filled else 0.5)


# ============================================================================ Fig. 3
REF_CELL = (30.0, 50.0)                                             # SNR, Neff
CORNER_CELLS = [(20.0, 25.0), (20.0, 100.0), (50.0, 25.0), (50.0, 100.0)]
MDE_REF_RANGE = {"NDI": (0.050, 0.070), "FISO": (0.030, 0.045), "ODI": (0.020, 0.025)}   # synthetic


def synth_fig3(rng):
    """Synthetic rho = sigma / sigma_prior at the reference (C = 1 - rho^2), then synthetic MDE at the
    reference cell and at the four corner cells (corners within about +-30 % of the reference)."""
    # Pattern of the C3 cards at the reference (cards_summary.json): NDI C >= 0.75 in all 28 pairs, prior-dependent
    # except the matched pair; FISO prior-dependent in all 28, partial in 1 / 0 / 4 / 5 pairs of types
    # matched / TE / TE+TR / TE+TR+scheme; ODI identified and stable. Values of C are synthetic.
    fiso_part = {0: 1, 1: 0, 3: 4, 4: 5}
    C, prior_dep = ({c: {} for c in COMPONENTS} for _ in range(2))
    for p in PAIRS:
        t = acq_type(p)
        part = fiso_part.get(t, 0) > 0
        if part:
            fiso_part[t] -= 1
        for c in COMPONENTS:
            if c == "ODI":
                rho, pd = rng.uniform(0.12, 0.16), False
            elif c == "NDI":
                rho, pd = rng.uniform(0.35, 0.44), p != REF_PAIR
            else:
                rho, pd = (rng.uniform(0.51, 0.54) if part else rng.uniform(0.42, 0.49)), True
            C[c][p], prior_dep[c][p] = 1 - rho ** 2, pd
    mde_ref, mde_corner = ({c: {} for c in COMPONENTS} for _ in range(2))
    for p in PAIRS:
        t = acq_type(p)
        for c in COMPONENTS:
            lo, hi = MDE_REF_RANGE[c]
            v = lo + min(1.0, 0.12 * t + rng.uniform(0, 0.5)) * (hi - lo)
            f = (rng.uniform(1.22, 1.30), rng.uniform(0.96, 1.08), rng.uniform(0.90, 1.02), rng.uniform(0.70, 0.78))
            mde_ref[c][p] = v
            mde_corner[c][p] = [v * x for x in f]          # CORNER_CELLS order
    return dict(pairs=PAIRS, ref_pair=REF_PAIR, stratum={p: acq_type(p) for p in PAIRS}, C=C, prior_dep=prior_dep,
                ref_cell=REF_CELL, corner_cells=CORNER_CELLS, mde_ref=mde_ref, mde_corner=mde_corner)


def draw_fig3(d, *, template, name="Fig3_S4_identifiability"):
    """(a) shrinkage C and (b) MDE, one column per card pair shared by both panels, one row per component.
    Pairs are grouped by acquisition type and ordered within a type by |dTE|. (2026-11-01: (a) was a strip per
    acquisition type, where marks of near-equal C overlapped; it now uses the columns of (b).)"""
    pairs, ref = d["pairs"], d["ref_pair"]
    # identified = C >= 0.75 and the class did not change with the prior width (SPEC_S4.md:131 · REPORT.md:161)
    n_id = {c: sum(d["C"][c][p] >= C_ID and not d["prior_dep"][c][p] for p in pairs) for c in COMPONENTS}
    if template:
        assert n_id == {"NDI": 1, "FISO": 0, "ODI": 28}, n_id
        assert [p for p in pairs if d["C"]["NDI"][p] >= C_ID and not d["prior_dep"]["NDI"][p]] == [ref], "NDI"
    n_type = [sum(d["stratum"][p] == k for p in pairs) for k in range(5)]
    assert n_type == [1, 3, 0, 9, 15], n_type          # card design + same-TR rule of the S4/V4 code
    assert tuple(d["ref_cell"]) == REF_CELL and [tuple(x) for x in d["corner_cells"]] == CORNER_CELLS
    assert all(len(d["mde_corner"][c][p]) == 4 for c in COMPONENTS for p in pairs)

    C_AX0 = 0.5                     # (a) shows C from 0.5: no cell lies below (asserted), so the class < 0.25 is empty
    assert min(d["C"][c][p] for c in COMPONENTS for p in pairs) >= C_AX0
    # columns: one per pair, grouped by type, |dTE| within a type; an empty slot for a type without pairs
    order = sorted(pairs, key=lambda p: (d["stratum"][p], abs(DESIGN[p][0]), pname(p)))
    xpos, spans, x = {}, [], 0.0
    for k in range(5):
        ps = [p for p in order if d["stratum"][p] == k]
        x0 = x
        for p in ps:
            xpos[p] = x + 0.5
            x += 1
        x += 0 if ps else 1.0
        spans.append((x0, x, len(ps)))
        x += 0.8
    xmax = x - 0.8

    H = 169.0
    fig = fs.new_fig(190, H)
    cv = canvas(fig)
    XB0, XB1 = 27.0, 170.0          # plot columns; 170-190 mm holds the band names of (a)
    RH, RG = 15.0, 2.5              # row height and gap (mm)
    YB, YA = 36.0, 98.0             # bottom of the lowest (b) row and of the lowest (a) row
    ya = {c: YA + (2 - i) * (RH + RG) for i, c in enumerate(COMPONENTS)}
    yb = {c: YB + (2 - i) * (RH + RG) for i, c in enumerate(COMPONENTS)}
    top_a, top_b = ya["NDI"] + RH, yb["NDI"] + RH

    def frame(a):
        a.set_xlim(-0.3, xmax + 0.3)
        a.tick_params(axis="x", length=0)
        a.set_xticks([])
        for x0, x1, n in spans[1:]:
            a.axvline(x0 - 0.4, color=AXIS, lw=0.5, zorder=0.7)

    def row_label(c, y0):
        cv.text(XB0 - 13.0, y0 + RH / 2 - 1.2, c, rotation=90, ha="center", va="center", fontsize=PT, color=INK,
                fontweight="bold")
        cv.plot([XB0 - 13.0], [y0 + RH / 2 + (4.6 if c != "FISO" else 5.2)], marker=CMRK[c], color=CCOL[c],
                ms=CMS[c], mec=SURF, mew=0.5)

    # ---------------- (a) shrinkage C per pair: bands = classes at the reporting prior; open = prior-dependent
    for c in COMPONENTS:
        a = fs.ax_mm(fig, XB0, ya[c], XB1 - XB0, RH)
        frame(a)
        a.set_ylim(C_AX0, 1.0)
        a.axhspan(C_AX0, C_ID, color=fs.ZONE_PART, lw=0, zorder=0)
        a.axhline(C_ID, color=AXIS, lw=0.5, zorder=0.6)
        for p in order:
            filled = not d["prior_dep"][c][p]
            a.plot([xpos[p]], [d["C"][c][p]], ls="none", marker=CMRK[c], ms=CMS[c],
                   mfc=CCOL[c] if filled else SURF, mec=CCOL[c], mew=0.9, zorder=3, clip_on=False)
        a.set_yticks([C_AX0, C_ID, 1.0])
        a.set_yticklabels(["0.5", "0.75", "1"])
        a.spines["bottom"].set_visible(False)
        for y0, y1, t in ((C_ID, 1.0, "identified\nif filled"), (C_AX0, C_ID, "partial\nif filled")):
            a.text(1.012, (y0 + y1) / 2, t, transform=a.get_yaxis_transform(), ha="left", va="center",
                   fontsize=PT, color=INK2, linespacing=1.05, clip_on=False)
        row_label(c, ya[c])
        if c == "NDI":
            fs.panel(a, "a", "Shrinkage C per card pair", dx_mm=-25, dy_mm=14.0)
            # acquisition types over the columns (shared with (b))
            names = [s_.replace("\n", " ") for s_ in STRATA]
            for (x0, x1, n), s_ in zip(spans, names):
                a.annotate(f"{s_}\n{n} pair{'s' if n != 1 else ''}", ((x0 + x1) / 2, 1), xycoords=("data", "axes fraction"),
                           xytext=(0, 3.2), textcoords="offset points", ha="center", va="bottom", fontsize=PT,
                           color=INK2, linespacing=1.15, annotation_clip=False)
                if n:
                    a.plot([x0 + 0.1, x1 - 0.1], [1.04, 1.04], transform=a.get_xaxis_transform(), color=AXIS,
                           lw=0.5, clip_on=False)
            if template:
                fs.template_mark(a, 1.0, 1 + 14.2 / RH, ha="right", va="bottom").set_clip_on(False)
    cv.text(XB0 - 19.5, ya["FISO"] + RH / 2, "Shrinkage C (worst of nine k)", rotation=90, ha="center",
            va="center", fontsize=PT, color=INK)
    # key of (a), on its header line: fill, and the identified counts
    ka = top_a + 15.2
    kx = 92.0
    for j, c in enumerate(COMPONENTS):
        legend_marker(cv, kx - 1.0 + j * 2.6, ka, c, filled=True)
        legend_marker(cv, kx + 43.0 + j * 2.6, ka, c, filled=False)
    cv.text(kx + 6.6, ka, "class kept at all prior widths", fontsize=PT, va="center")
    cv.text(kx + 50.6, ka, "prior-dependent, not counted", fontsize=PT, va="center")
    cv.text(XB1, top_a + 10.5, f"identified (filled, C ≥ 0.75): NDI {n_id['NDI']}, FISO {n_id['FISO']}, "
            f"ODI {n_id['ODI']} of {len(pairs)} pairs", fontsize=PT, va="center", ha="right", color=INK)

    # ---------------- (b) MDE per pair: reference cell (mark) and range over the four corner cells (bar)
    hi_all = max(max(d["mde_corner"][c][p] + [d["mde_ref"][c][p]]) for c in COMPONENTS for p in pairs)
    ymax = 0.05 * math.ceil(hi_all / 0.05 + 1e-9)          # one MDE scale for the three components
    for i, c in enumerate(COMPONENTS):
        a = fs.ax_mm(fig, XB0, yb[c], XB1 - XB0, RH)
        frame(a)
        a.set_ylim(0, ymax)
        a.set_yticks(np.arange(0, ymax + 1e-9, 0.05))
        a.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))
        fs.hgrid(a)
        for p in order:
            cs = d["mde_corner"][c][p]
            a.plot([xpos[p]] * 2, [min(cs), max(cs)], color=CCOL[c], lw=1.0, zorder=3, solid_capstyle="butt")
            filled = not d["prior_dep"][c][p]                  # fill as in (a)
            a.plot([xpos[p]], [d["mde_ref"][c][p]], ls="none", marker=CMRK[c], ms=CMS[c],
                   mfc=CCOL[c] if filled else SURF, mec=CCOL[c], mew=0.9, zorder=4)
        if c == "ODI":
            a.set_xticks([xpos[p] for p in order])
            a.set_xticklabels([pname(p) for p in order], rotation=90, ha="center", va="top")
            a.tick_params(axis="x", pad=2)
        row_label(c, yb[c])
        if c == "NDI":
            fs.panel(a, "b", "Minimum detectable difference, worst k", dx_mm=-25, dy_mm=3.0)
            if template:
                fs.template_mark(a, 1.0, 1 + 3.2 / RH, ha="right", va="bottom").set_clip_on(False)
    cv.text(XB0 - 19.5, yb["FISO"] + RH / 2, "MDE (index units)", rotation=90, ha="center", va="center",
            fontsize=PT, color=INK)
    kb = top_b + 3.0 + 1.2
    cv.text(76.0, kb, "mark: SNR 30 · Neff 50 (fill as in (a))", fontsize=PT, va="center")
    cv.plot([126.0] * 2, [kb - 1.6, kb + 1.6], color=INK2, lw=1.0, solid_capstyle="butt")
    cv.text(128.2, kb, "range over SNR {20, 50} × Neff {25, 100}", fontsize=PT, va="center")
    fs.save(fig, name)


# ============================================================================ Fig. 4
def synth_fig4(rng):
    """Design columns and the five no-data cells are real (fixed document v10); coverage,
    verdicts, flags, subject and unit counts and attribute shares are synthetic."""
    nodata = sorted(NO_DATA, key=acq_type)        # grouped by type (code rule); fixed-document order within
    odi = sorted([p for p in PAIRS if p not in nodata],
                 key=lambda p: (acq_type(p), abs(DESIGN[p][0]), pname(p)))
    order = [("NDI", REF_PAIR)] + [("ODI", p) for p in odi]
    assert len(order) == 24 and sum(c == "ODI" for c, _ in order) == 23
    fallback = rng.integers(1, 24)
    gap_cells = set(rng.choice(np.arange(24), 3, replace=False).tolist())
    c7_cells = set(rng.choice(np.arange(24), 2, replace=False).tolist())
    per_pair = {}                    # subjects, units and attribute presence belong to the card pair
    for p in [REF_PAIR] + odi:
        same_vendor = vendor(p[0]) == vendor(p[1])
        subj = int(rng.integers(3, pair_cap(p) + 1))
        per_pair[p] = dict(
            subjects=subj,
            units=27 * subj - int(rng.integers(0, int(0.04 * 27 * subj) + 1)),   # a few missing ROI units
            pres={"noise": rng.random() < 0.35, "range": rng.random() < 0.25,
                  "vendor": not same_vendor,                     # attribute 3 is fixed by the devices
                  "none": same_vendor and rng.random() < 0.92})  # "none" needs the same manufacturer
    cells = []
    for i, (comp, p) in enumerate(order):
        base = float(np.clip(rng.normal(0.872, 0.042), 0.74, 0.975))
        ks = np.clip(base + rng.normal(0, 0.013, 9), 0, 1)
        thr = 0.80 if i == fallback else float(rng.uniform(0.82, 0.872))
        cov = float(np.median(ks))
        dec = {k: (float(np.clip(1 - cov + rng.normal(0.05 if k == "noise" else 0.0, 0.05), 0, 0.5))
                   if on else None) for k, on in per_pair[p]["pres"].items()}
        cells.append(dict(comp=comp, pair=p, design=design_cells(p), units=per_pair[p]["units"],
                          subjects=per_pair[p]["subjects"], k_cov=ks, coverage=cov, threshold=thr,
                          fallback=bool(i == fallback), verdict=bool(cov >= thr), gap=i in gap_cells,
                          c7_differs=i in c7_cells, decomp=dec))
    nd = [dict(pair=p, design=design_cells(p), units=0, subjects=0) for p in nodata]
    return dict(cells=cells, nodata=nd)


def draw_fig4(d, *, template, name="Fig4_V4_travelling_heads"):
    cells, nodata = d["cells"], d["nodata"]
    nd_pairs = [r["pair"] for r in nodata]
    assert sorted(nd_pairs) == sorted(NO_DATA) and nd_pairs == sorted(nd_pairs, key=acq_type), nd_pairs
    assert all(r["units"] == 0 and r["subjects"] == 0 for r in nodata)
    assert all(0 < c["units"] <= 27 * c["subjects"] for c in cells)
    assert all(c["subjects"] <= pair_cap(c["pair"]) for c in cells), "subjects above protocol size"
    assert not ({c["pair"] for c in cells} & set(NO_DATA)) and len(cells) == 24
    ndi = [c for c in cells if c["comp"] == "NDI"]
    odi = [c for c in cells if c["comp"] == "ODI"]
    pitch = 3.75
    n_lines = 1 + len(ndi) + 1 + len(odi) + 1 + len(nodata)
    H_head, H_key = 13.0, 14.0
    H = round(3 + H_head + n_lines * pitch + 2 + H_key)
    fig = fs.new_fig(190, H)
    cv = canvas(fig)

    X_PAIR, X_DTE, X_DTR, X_SCH, X_UNITS = 2.5, 45.5, 55.0, 62.5, 81.0
    PX0, PX1, CMIN, CMAX = 84.5, 122.0, 0.60, 1.00
    X_VERD, X_FLAG = 125.0, 132.5
    DEC_KEYS = ["noise", "range", "vendor", "none"]
    DEC_LAB = ["noise", "acq.\nrange", "other\nvendor", "none"]
    DX0, DW, DG = 150.5, 8.4, 1.0
    XR = DX0 + 4 * DW + 3 * DG          # right edge of the table

    def px(c):
        return PX0 + (c - CMIN) / (CMAX - CMIN) * (PX1 - PX0)

    ytop_rule = H - 3 - H_head
    yh = ytop_rule + 1.2                                  # header baseline
    for x, t, ha in ((X_PAIR, "Card pair A ~ B", "left"), (X_DTE, "ΔTE\n(ms)", "right"),
                     (X_DTR, "ΔTR\n(s)", "right"), (X_SCH, "Scheme", "center"),
                     (X_UNITS, "Units\n(subjects)", "right"), (X_VERD, "Verdict", "left")):
        cv.text(x, yh, t, fontsize=PT, color=INK2, ha=ha, va="bottom", linespacing=1.1)
    cv.text((38.0 + 67.0) / 2, yh + 7.0, "Acquisition, B − A", fontsize=PT, color=INK2,
            ha="center", va="bottom")
    cv.plot([38.0, 67.0], [yh + 6.2, yh + 6.2], color=AXIS, lw=0.5)
    cv.text((PX0 + PX1) / 2, yh + 3.4, "Coverage of zero", fontsize=PT, color=INK2, ha="center", va="bottom")
    for t in (0.6, 0.7, 0.8, 0.9, 1.0):
        cv.text(px(t), yh, f"{t:.1f}" if t < 1 else "1", fontsize=PT, color=INK2, ha="center", va="bottom")
    for k, lab in enumerate(DEC_LAB):
        cv.text(DX0 + k * (DW + DG) + DW / 2, yh, lab, fontsize=PT, color=INK2, ha="center",
                va="bottom", linespacing=1.05)
    cv.text((DX0 + XR) / 2, yh + 7.0, "Not covered, by attribute", fontsize=PT, color=INK2,
            ha="center", va="bottom")
    cv.plot([DX0, XR], [yh + 6.2, yh + 6.2], color=AXIS, lw=0.5)
    cv.plot([1.5, XR], [ytop_rule, ytop_rule], color=AXIS, lw=0.5)
    if template:
        t = cv.text(XR, H - 1.0, fs.TEMPLATE_TEXT, fontsize=PT, color=fs.WM_C, ha="right", va="top")
        t.set_in_layout(False)

    y = ytop_rule - pitch * 0.72
    row_ys = []

    def section(text, right=None):
        nonlocal y
        cv.text(X_PAIR, y, text, fontsize=PT, color=INK, fontweight="bold", va="center")
        if right:
            cv.text(X_VERD, y, right, fontsize=PT, color=INK2, va="center", ha="left")
        y -= pitch

    def design_row(r, c1, c2):
        dte, dtr, sch = r["design"]
        cv.text(X_PAIR, y, pname(r["pair"]), fontsize=PT, color=c1, va="center")
        cv.text(X_DTE, y, dte, fontsize=PT, color=c2, va="center", ha="right")
        cv.text(X_DTR, y, dtr, fontsize=PT, color=c2, va="center", ha="right")
        cv.text(X_SCH, y, sch, fontsize=PT, color=c2, va="center", ha="center")
        cv.text(X_UNITS, y, f'{r["units"]} ({r["subjects"]})', fontsize=PT, color=c2, va="center", ha="right")

    def plot_row(c, zebra):
        nonlocal y
        if zebra:
            cv.add_patch(Rectangle((1.5, y - pitch / 2), XR - 1.5, pitch, color=ZEBRA, lw=0, zorder=0))
        design_row(c, INK, INK2)
        col = CCOL[c["comp"]]
        ks = np.asarray(c["k_cov"])
        cv.plot([px(ks.min()), px(ks.max())], [y, y], color=col, lw=1.0, zorder=3, solid_capstyle="butt")
        cv.plot([px(c["threshold"])] * 2, [y - 1.25, y + 1.25], color=MUTED if c["fallback"] else INK,
                lw=1.0, zorder=4, solid_capstyle="butt")
        ok = c["verdict"]
        cv.plot([px(c["coverage"])], [y], ls="none", marker=CMRK[c["comp"]], ms=CMS[c["comp"]],
                mfc=col if ok else SURF, mec=col, mew=0.9, zorder=5)
        cv.text(X_VERD, y, "pass" if ok else "fail", fontsize=PT, va="center", color=INK,
                fontweight="normal" if ok else "bold")
        flags = [t for t, on in (("gap", c["gap"]), ("C7≠", c["c7_differs"])) if on]
        cv.text(X_FLAG, y, " ".join(flags), fontsize=PT, va="center", color=INK2)
        for k, key in enumerate(DEC_KEYS):
            xk = DX0 + k * (DW + DG)
            v = c["decomp"].get(key)
            if v is not None:
                cv.add_patch(Rectangle((xk, y - 1.3), DW, 2.6, lw=0, zorder=3,
                                       color=matplotlib.cm.Greys(0.08 + 0.82 * min(v, 0.5) / 0.5)))
            else:
                cv.text(xk + DW / 2, y, "–", fontsize=PT, color=MUTED, ha="center", va="center")
        row_ys.append(y)
        y -= pitch

    def summary(cs):
        if template:                       # no counts that read like results for real devices
            return "n pass · n fail"
        return f"{sum(c['verdict'] for c in cs)} pass · {sum(not c['verdict'] for c in cs)} fail"

    section("NDI · unblinded check", summary(ndi))
    for c in ndi:
        plot_row(c, False)
    section("ODI · blinded confirmation", summary(odi))
    for k, c in enumerate(odi):
        plot_row(c, k % 2 == 1)
    section(f"No data · no subject scanned on both acquisitions ({len(nodata)} ODI cells)")
    for r in nodata:
        design_row(r, MUTED, MUTED)
        y -= pitch
    ybot = row_ys[-1] - pitch / 2
    for t in (0.6, 0.7, 0.8, 1.0):
        cv.plot([px(t)] * 2, [ybot, ytop_rule], color=GRID, lw=0.5, zorder=1)
    cv.plot([px(0.9)] * 2, [ybot, ytop_rule], color=INK2, lw=0.6, zorder=1)
    y_rule = y + pitch / 2
    cv.plot([1.5, XR], [y_rule, y_rule], color=AXIS, lw=0.5)

    # ---------------- key
    k1, k2 = y_rule - 5.0, y_rule - 10.5
    oc = CCOL["ODI"]
    legend_marker(cv, 4, k1, "ODI", filled=True)
    cv.text(6.5, k1, "pass", fontsize=PT, va="center")
    legend_marker(cv, 17, k1, "ODI", filled=False)
    cv.text(19.5, k1, "fail", fontsize=PT, va="center")
    cv.plot([30, 34], [k1, k1], color=oc, lw=1.0, solid_capstyle="butt")
    cv.text(35.5, k1, "range over the nine k (mark = median)", fontsize=PT, va="center")
    cv.plot([4, 4], [k2 - 1.25, k2 + 1.25], color=INK, lw=1.0, solid_capstyle="butt")
    cv.text(6.5, k2, "own threshold", fontsize=PT, va="center")
    cv.plot([30, 30], [k2 - 1.25, k2 + 1.25], color=MUTED, lw=1.0, solid_capstyle="butt")
    cv.text(32.5, k2, "fallback 0.80", fontsize=PT, va="center")
    cv.plot([55, 55], [k2 - 1.25, k2 + 1.25], color=INK2, lw=0.6)
    cv.text(57.5, k2, "nominal 0.90", fontsize=PT, va="center")
    cv.text(84, k2, "gap: gap-detected flag    C7≠: C7 arm disagrees", fontsize=PT, va="center", color=INK)
    gx, gw = DX0, 4 * DW + 3 * DG
    n = 50
    for k in range(n):
        cv.add_patch(Rectangle((gx + k * gw / n, k1 - 1.2), gw / n + 0.02, 2.4, lw=0,
                               color=matplotlib.cm.Greys(0.08 + 0.82 * k / (n - 1))))
    cv.text(gx, k2, "0", fontsize=PT, color=INK2, ha="left", va="center")
    cv.text(gx + gw, k2, "≥ 0.5", fontsize=PT, color=INK2, ha="right", va="center")
    cv.text(gx + gw / 2, k2, "share", fontsize=PT, color=INK2, ha="center", va="center")
    fs.save(fig, name)


# ============================================================================ Fig. 5
def in_D(theta):
    t = np.mod(theta, 180.0)
    return ((t >= TE_WEDGE[0]) & (t <= TE_WEDGE[1])) | ((t >= TR_FAN[0]) & (t <= TR_FAN[1]))


TA_UNDET_MIN = 23   # T-A undeterminable when this many of the 45 pairs are within the margin (Methods 3.7.6)


MARGIN_C = "#ebeae5"     # inside the T-A margin: no direction is read


def draw_fig5(d, *, template, name="Fig5_Tong_negative_control"):
    """Fig. 5 from the stored Tong judge output (module docstring: d). Per-pair effects and per-unit
    shares are not stored, so (a) shows the plane, D and the margin, (b) the pairs within the margin
    against the T-A rule, (c) the T-B share over all units and (d) the T-C cell."""
    n_pairs, n_in, margin, limit = d["n_pairs"], d["n_within"], d["margin"], d["limit"]
    assert d["tc"] is None and n_in >= limit, "only the undeterminable branch of T-A is drawn from stored values"
    assert abs(d["d_share"] - D_SHARE) < 1e-9, d["d_share"]
    fig = fs.new_fig(190, 100)
    cv = canvas(fig)
    R1, H1 = 56.0, 36.0      # top row bottom, height
    R2, H2 = 10.0, 30.0      # bottom row bottom, height

    # ---------------- (a) card plane: the sectors of D and the margin
    R = 0.36
    ax = fs.ax_mm(fig, 16, R1, 2 * H1, H1)       # equal aspect: width = 2 * height
    ax.set_xlim(-R, R)
    ax.set_ylim(0, R)
    ax.set_aspect("equal")
    for (t0, t1), col, al in ((TE_WEDGE, fs.TE_C, 0.16), (TR_FAN, fs.TR_C, 0.28)):
        ax.add_patch(Wedge((0, 0), 2 * R, t0, t1, color=col, alpha=al, lw=0, zorder=0))
    ax.add_patch(Wedge((0, 0), margin, 0, 180, facecolor=MARGIN_C, edgecolor="none", alpha=0.85, zorder=1))
    ax.add_patch(Arc((0, 0), 2 * margin, 2 * margin, theta1=0, theta2=180, color=INK2, lw=0.8, zorder=2))
    tm = math.radians(np.mean(TE_WEDGE))
    ax.text(R * 0.92 / math.tan(tm) + 0.03, R * 0.92, "TE", fontsize=PT, color=fs.TE_C,
            fontweight="bold", ha="left", va="center")
    ax.text(-0.045, R * 0.92, "TR", fontsize=PT, color="#9a6a00", fontweight="bold", ha="right", va="center")
    ax.text(-0.14, 0.12, f"{n_in} of {n_pairs}\ncentre-pair effects\nlie inside", fontsize=PT,
            color=INK, ha="center", va="center", linespacing=1.15, zorder=3)
    ax.text(-0.14, 0.035, f"margin {margin:.3f}", fontsize=PT, color=INK2, ha="center",
            va="center", zorder=3)
    ax.set_xlabel("Δx = Δ logit NDI", labelpad=2)
    ax.set_ylabel("Δy = Δ log(iso/intra)", labelpad=2)
    ax.set_xticks([-0.3, -0.15, 0, 0.15, 0.3])
    ax.set_yticks([0, 0.1, 0.2, 0.3])
    fs.panel(ax, "a", "Card plane, sectors and margin", dx_mm=-12)
    if template:
        fs.template_mark(ax, 0.005, 0.995, ha="left", va="top")

    # ---------------- (b) T-A: the number of pairs within the margin, against the rule
    bb, bh = R1 + 8.0, 22.0
    bx = fs.ax_mm(fig, 112, bb, 72, bh)
    bx.axis("off")
    bx.text(0.0, 0.80, f"{n_in} of {n_pairs}", transform=bx.transAxes, fontsize=14, fontweight="bold", color=INK,
            ha="left", va="center")
    bx.text(0.0, 0.50, "centre pairs within the margin", transform=bx.transAxes, fontsize=PT, color=INK,
            ha="left", va="center")
    bx.text(0.0, 0.22, f"T-A undeterminable (rule: ≥ {limit} of {n_pairs}); no rotation test",
            transform=bx.transAxes, fontsize=PT, color=INK2, ha="left", va="center")
    fs.panel(bx, "b", "Direction (T-A)", dx_mm=-10, dy_mm=(R1 + H1) - (bb + bh) + 2.0)
    if template:
        fs.template_mark(bx, 1.0, 1.03, ha="right", va="bottom").set_clip_on(False)

    # ---------------- (c) T-B: share of units at retest level, over all units
    tx = fs.ax_mm(fig, 28, R2, 64, H2)
    pos = {"x": 0, "y": 1, "ODI": 2.45}
    for k in ("x", "y", "ODI"):
        n_rt, n_tot, sh = d["tb"][k]
        ok = sh >= 0.90                                   # filled = retest level, open = exceeds
        col, mk, ms = (CCOL["ODI"], "^", 5.0) if k == "ODI" else (INK, "D", 4.2)
        tx.plot([pos[k]], [sh], ls="none", marker=mk, ms=ms, mfc=col if ok else SURF, mec=col, mew=0.9, zorder=4)
        tx.text(pos[k] + 0.17, sh, f"{sh:.3f}\n{'retest level' if ok else 'exceeds'}", fontsize=PT, color=INK,
                ha="left", va="center", linespacing=1.1, zorder=4)
    tx.axhline(0.90, color=INK2, lw=0.6, zorder=1)
    tx.set_xlim(-0.6, 3.45)
    tx.set_ylim(0.70, 1.0)
    tx.set_xticks([pos[k] for k in ("x", "y", "ODI")])
    tx.set_xticklabels([f"{k}\n{d['tb'][k][0]} of {d['tb'][k][1]}" for k in ("x", "y", "ODI")], linespacing=1.2)
    tx.tick_params(axis="x", length=0)
    tx.axvline(1.85, color=GRID, lw=0.5)
    tx.set_yticks([0.7, 0.8, 0.9, 1.0])
    tx.set_ylabel("Share at retest level", labelpad=2)
    fs.hgrid(tx)
    fs.panel(tx, "c", "Magnitude (T-B)", dx_mm=-24)
    if template:
        fs.template_mark(tx, 0.02, 0.015, ha="left").set_zorder(5)

    # ---------------- (d) T-C cell
    gx, gy, cw, chh = 130.0, R2, 24.0, 12.5
    texts = {(False, True): "strongest\npass", (False, False): "large but\nseparable",
             (True, True): "small:\nlimitation", (True, False): "worst:\nsite term"}
    hit_cell = d["tc"]
    for (al, rt), txt in texts.items():
        cx = gx + (0 if rt else cw)
        cy = gy + (0 if al else chh)
        hit = hit_cell is not None and (al, rt) == tuple(hit_cell)
        cv.add_patch(Rectangle((cx + 0.35, cy + 0.35), cw - 0.7, chh - 0.7, lw=1.0 if hit else 0,
                               edgecolor=INK, facecolor=CELL_ON if hit else CELL_OFF, zorder=1))
        cv.text(cx + cw / 2, cy + chh / 2, txt, fontsize=PT, ha="center", va="center",
                color=INK if hit else MUTED, linespacing=1.15)
    cv.text(gx + cw / 2, gy + 2 * chh + 1.0, "retest level", fontsize=PT, ha="center", va="bottom", color=INK2)
    cv.text(gx + 1.5 * cw, gy + 2 * chh + 1.0, "exceeds", fontsize=PT, ha="center", va="bottom", color=INK2)
    cv.text(gx - 1.2, gy + chh * 1.5, "not\naligned", fontsize=PT, ha="right", va="center",
            color=INK2, linespacing=1.1)
    cv.text(gx - 1.2, gy + chh * 0.5, "aligned", fontsize=PT, ha="right", va="center", color=INK2)
    if hit_cell is None:
        cv.text(gx + cw, gy - 1.2, "No cell outlined: T-C is undeterminable\nbecause T-A is undeterminable", fontsize=PT,
                ha="center", va="top", color=INK, linespacing=1.15)
    dax = fs.ax_mm(fig, gx, gy, 2 * cw, 2 * chh)
    dax.axis("off")
    fs.panel(dax, "d", "Both (T-C)", dx_mm=-28, dy_mm=(R2 + H2) - (gy + 2 * chh) + 2.0)
    if template:
        t = cv.text(gx + cw, gy - 1.8, fs.TEMPLATE_TEXT, fontsize=PT, color=fs.WM_C, ha="center", va="top")
        t.set_in_layout(False)
    fs.save(fig, name)


# ============================================================================ Fig. S2
def synth_figS1(rng):
    TEs = np.array([62, 72, 82, 92, 102, 112, 122, 132], float)
    tc = TEs - 97.0
    rep_te = np.array([62.0, 132.0])
    out = {}
    for coord, (b_mu, b_sd, lev_sd, s_rt) in {"x": (0.0075, 0.0011, 0.35, 0.030),
                                               "y": (0.0105, 0.0018, 0.45, 0.050)}.items():
        units = []
        for si in range(3):
            for _ in range(27):
                a = rng.normal(0.2 if coord == "x" else -2.2, lev_sd)
                b = rng.normal(b_mu, b_sd)
                c = rng.normal(0, 6.5e-5 if coord == "x" else 1.0e-4)
                v = a + b * tc + c * tc ** 2 + rng.normal(0, s_rt / 1.6, len(TEs))
                rep = a + b * (rep_te - 97) + c * (rep_te - 97) ** 2 + rng.normal(0, s_rt / 1.6, 2)
                quad = np.polyfit(tc, v, 2)
                lin = np.polyfit(tc, v, 1)
                ratio = abs(quad[0]) * 35 ** 2 / (2 * math.sqrt(2) * s_rt)
                lo = np.polyfit(TEs[:2], v[:2], 1)
                hi = np.polyfit(TEs[-2:], v[-2:], 1)
                full = np.polyfit(TEs, v, 1)
                units.append(dict(subject=si, v=v, rep=rep, lin=lin, ratio=ratio,
                                  icpt=(lo[1], hi[1], full[1]), single=(v.min(), v.max())))
        out[coord] = units
    return dict(TE=TEs, TE_rep=rep_te, units=out)


def draw_figS1(d, *, template, name="FigS2_Wong_multiTE"):
    TEs, rep_te, out = np.asarray(d["TE"]), np.asarray(d["TE_rep"]), d["units"]
    tc = TEs - 97.0
    REP_X = np.array([50.0, 144.0])  # slots for the repeats at 62 and 132 ms, beside the TE axis (not TE values)
    marks = ["h", "p", "v"]          # subjects: shapes that no component uses (circle, square, triangle are NDI, FISO, ODI)
    G = fs.GREYS
    dodge = np.array([-1.6, 0.0, 1.6])
    fig = fs.new_fig(190, 116)
    cv = canvas(fig)
    R1, H1 = 69.0, 40.0
    R2, H2 = 12.0, 36.0

    def trend_panel(ax, coord, letter, header, ylabel):
        us = out[coord]
        for u in us:
            centred = u["v"] - np.polyval(u["lin"], 0.0)
            ax.plot(TEs + dodge[u["subject"]], centred, color=G[u["subject"]], lw=0.3, alpha=0.3, zorder=1)
        for si in range(3):
            su = [u for u in us if u["subject"] == si]
            cen = np.array([u["v"] - np.polyval(u["lin"], 0.0) for u in su])
            rep = np.array([u["rep"] - np.polyval(u["lin"], 0.0) for u in su])
            med = np.median(cen, axis=0)
            fit = np.polyfit(tc, med, 1)
            ax.plot(TEs + dodge[si], np.polyval(fit, tc), color=G[si], lw=1.0, zorder=3)
            ax.plot(TEs + dodge[si], med, ls="none", marker=marks[si], ms=3.6, color=G[si], mec=SURF,
                    mew=0.5, zorder=4)
            ax.plot(REP_X + dodge[si], np.median(rep, axis=0), ls="none",
                    marker="x", ms=3.6, color=G[si], mew=0.9, zorder=4)
        ax.axvline(97, color=GRID, lw=0.5, zorder=0)
        ax.set_xlim(43, 151)
        ax.set_xticks(list(TEs) + list(REP_X))
        ax.set_xticklabels([f"{t:.0f}" for t in TEs] + ["rep", "rep"])
        ax.set_xlabel("TE (ms)", labelpad=2)
        ax.set_ylabel(ylabel, labelpad=2)
        fs.hgrid(ax)
        fs.panel(ax, letter, header, dx_mm=-13)
        if template:
            fs.template_mark(ax, 0.99, 0.015).set_zorder(5)

    ax = fs.ax_mm(fig, 18, R1, 72, H1)
    trend_panel(ax, "x", "a", "x against TE", "x relative to fit at 97 ms")
    bx = fs.ax_mm(fig, 113, R1, 72, H1)
    trend_panel(bx, "y", "b", "y against TE", "y relative to fit at 97 ms")

    # ---------------- (c) intercepts at TE -> 0
    ch = (H2 - 4) / 2
    for r, coord in enumerate(("x", "y")):
        cx = fs.ax_mm(fig, 18, R2 + (1 - r) * (ch + 4), 46, ch)
        us = out[coord]
        for si in range(3):
            su = [u for u in us if u["subject"] == si]
            s_lo = np.median([u["single"][0] for u in su])
            s_hi = np.median([u["single"][1] for u in su])
            cx.add_patch(Rectangle((si - 0.13, s_lo), 0.26, s_hi - s_lo, color="#dddcd6", lw=0, zorder=1))
            ic = np.median(np.array([u["icpt"] for u in su]), axis=0)
            for k, (mk, dxo) in enumerate((("<", -0.24), (">", 0.24), ("P", 0.0))):
                cx.plot([si + dxo], [ic[k]], marker=mk, ms=3.6, color=G[si], mec=SURF, mew=0.4, zorder=3)
        cx.set_xlim(-0.6, 2.6)
        cx.set_xticks([0, 1, 2])
        cx.set_xticklabels(["S1", "S2", "S3"] if r == 1 else [])
        cx.tick_params(axis="x", length=0)
        cx.margins(y=0.16)
        cx.yaxis.set_major_locator(mticker.MaxNLocator(3))
        cx.set_ylabel(coord, labelpad=2, rotation=0, va="center", ha="right", fontweight="bold")
        fs.hgrid(cx)
        if r == 0:
            fs.panel(cx, "c", "Intercepts at TE → 0", dx_mm=-13)
        if template and r == 1:
            fs.template_mark(cx, 1.07, 0.0, ha="left", va="bottom").set_clip_on(False)
    # key for (c), to its right
    kx, ky = 70.0, R2 + H2 - 4.0
    for k, (mk, t) in enumerate((("<", "low TEs"), (">", "high TEs"), ("P", "all eight"))):
        cv.plot([kx + 1], [ky - k * 4.4], marker=mk, ms=3.6, color=INK2, mec=SURF, mew=0.4)
        cv.text(kx + 4.2, ky - k * 4.4, t, fontsize=PT, va="center")
    cv.add_patch(Rectangle((kx, ky - 3 * 4.4 - 1.0), 2.2, 2.0, color="#dddcd6", lw=0))
    cv.text(kx + 4.2, ky - 3 * 4.4, "single-TE range", fontsize=PT, va="center")

    # ---------------- (d) curvature ECDF, pooled over the 81 units; x solid, y dashed (greys are subjects)
    dx_ = fs.ax_mm(fig, 116, R2, 58, H2)
    shares, LS = {}, {"x": "-", "y": (0, (3.0, 1.6))}
    for coord in ("x", "y"):
        r = np.sort([u["ratio"] for u in out[coord]])
        F = np.arange(1, len(r) + 1) / len(r)
        dx_.step(np.concatenate([[0], r]), np.concatenate([[0], F]), where="post", color=INK, lw=1.0, ls=LS[coord])
        sh = float((r <= 1).mean())
        shares[coord] = sh
        dx_.plot([1], [sh], marker="D", ms=3.4, color=INK, mec=INK, mew=0.5, zorder=4)   # filled: share >= 0.8
    hi_c, lo_c = sorted(shares, key=shares.get, reverse=True)
    dx_.text(0.96, shares[hi_c] + 0.04, f"{hi_c}  {shares[hi_c]:.2f}", fontsize=PT, va="bottom", ha="right", color=INK)
    dx_.text(1.05, shares[lo_c] - 0.05, f"{lo_c}  {shares[lo_c]:.2f}", fontsize=PT, va="top", ha="left", color=INK)
    for k, coord in enumerate(("x", "y")):                      # line key, lower right (no data there)
        yk = 0.55 - k * 0.10
        dx_.plot([1.55, 1.75], [yk, yk], color=INK, lw=1.0, ls=LS[coord])
        dx_.text(1.80, yk, coord, fontsize=PT, va="center", ha="left", color=INK)
    dx_.axvline(1, color=INK2, lw=0.6)
    dx_.text(0.96, 0.05, "linear", fontsize=PT, color=INK2, ha="right", va="bottom")
    dx_.axhline(0.8, color=INK2, lw=0.5)
    dx_.set_xlim(0, 2.2)
    dx_.set_ylim(0, 1.02)
    dx_.set_xticks([0, 0.5, 1, 1.5, 2])
    dx_.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    dx_.set_xlabel("Curvature ratio", labelpad=2)
    dx_.set_ylabel("Share of units ≤ ratio", labelpad=2)
    for yg in (0.2, 0.4, 0.6):                                   # grid below the threshold only (labels sit above)
        dx_.axhline(yg, color=GRID, lw=fs.LW_THIN, zorder=0)
    fs.panel(dx_, "d", "Curvature", dx_mm=-16)
    if template:
        fs.template_mark(dx_, 0.99, 0.015).set_zorder(5)
    assert max(max(u["ratio"] for u in out[c]) for c in ("x", "y")) <= 2.2, "curvature ratio beyond the (d) axis"

    # ---------------- key for (a)-(c), in the empty upper left of (a); it must cover no data
    ax.set_ylim(*ax.get_ylim())                       # freeze the limits that the data set
    KEY = (0.02, 0.63, 0.69, 0.99)                    # x0, x1, y0, y1 in axes fraction
    for ln in list(ax.lines):
        if ln.get_color() == GRID:                    # the 97-ms guide stops below the key
            ln.set_ydata([0, KEY[2] - 0.02])
            continue
        xy = np.column_stack([ln.get_xdata(), ln.get_ydata()]).astype(float)
        if len(xy) > 1 and ln.get_linestyle() not in ("None", "", " ", "none"):   # lines: sample the segments
            t = np.linspace(0, 1, 25)[:, None]
            xy = np.concatenate([q0 + t * (q1 - q0) for q0, q1 in zip(xy[:-1], xy[1:])])
        pa = (ln.get_transform() - ax.transAxes).transform(xy)
        hit = (pa[:, 0] > KEY[0]) & (pa[:, 0] < KEY[1]) & (pa[:, 1] > KEY[2]) & (pa[:, 1] < KEY[3])
        assert not hit.any(), "the key of (a) would cover data"
    ax.add_patch(Rectangle((KEY[0], KEY[2]), KEY[1] - KEY[0], KEY[3] - KEY[2], transform=ax.transAxes,
                           facecolor=SURF, edgecolor="none", zorder=2.5))
    rows = (0.93, 0.83, 0.73)
    kt = dict(transform=ax.transAxes, fontsize=PT, va="center", zorder=6)
    km = dict(transform=ax.transAxes, clip_on=False, zorder=6)
    for si in range(3):
        ax.plot([0.045], [rows[si]], marker=marks[si], ms=3.8, color=G[si], mec=SURF, mew=0.5, **km)
        ax.text(0.075, rows[si], f"S{si + 1}", **kt)
    ax.plot([0.20], [rows[0]], marker="x", ms=3.6, color=INK2, mew=0.9, **km)
    ax.text(0.235, rows[0], "repeat (not fitted)", **kt)
    ax.plot([0.18, 0.22], [rows[1]] * 2, color=MUTED, lw=0.4, **km)
    ax.text(0.235, rows[1], "one ROI unit", **kt)
    ax.plot([0.18, 0.22], [rows[2]] * 2, color=INK, lw=1.0, **km)
    ax.text(0.235, rows[2], "fit to subject median", **kt)
    fs.save(fig, name)


# ============================================================================ main
# Only Fig. 4 is still drawn as a template (V4 not yet judged). Fig. 3, Fig. 5 and Fig. S2 are drawn
# from stored values by fig_real_data.py under the same file names, so their templates are not drawn
# here any more (synth_fig3 and synth_figS1 are kept for the record).
def main(which):
    for k in which:
        if k != "fig4":
            sys.exit(f"{k}: drawn from stored values by fig_real_data.py; no template is drawn")
        draw_fig4(synth_fig4(np.random.default_rng([SEED, 6])), template=True)


if __name__ == "__main__":
    main(sys.argv[1:] or ["fig4"])
