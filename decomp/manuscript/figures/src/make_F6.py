"""make_F6 draws Supplementary Fig. S6 (Figure 6 of earlier versions) from stored values only: comparison across echo times in the multi-TE series of three participants."""
import sys

sys.dont_write_bytecode = True            # no __pycache__ beside this file, in the repository or in the ledger folder

import csv  # noqa: E402
import datetime  # noqa: E402
import hashlib  # noqa: E402
import importlib.util  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import os  # noqa: E402
import pathlib  # noqa: E402
import re  # noqa: E402
import subprocess  # noqa: E402

os.environ["CUDA_VISIBLE_DEVICES"] = ""   # CPU only


def stop(msg):
    raise SystemExit(f"STOP F6: {msg}")


# Path — this file is in drafts/x15/ms/S44_files/src: the working folder data is parents[5], the repository is dmri_cov next to it (changed with DMRI_REPO).
HERE = pathlib.Path(__file__).resolve().parent
WORK = pathlib.Path(os.environ.get("DMRI_WORK") or HERE.parents[5]).resolve()
REPO = pathlib.Path(os.environ.get("DMRI_REPO") or WORK.parent / "dmri_cov").resolve()
DRAFTS = WORK / "stage1_work" / "drafts"
LEDGER = DRAFTS / "x15" / "ledger"
SOURCES = LEDGER / "sources"
X14 = DRAFTS / "x14"
REREAD_LOG = X14 / "REREAD_LOG.jsonl"
if not os.environ.get("FIG_OUT"):
    stop("set FIG_OUT (output folder)")
OUT = pathlib.Path(os.environ["FIG_OUT"]).resolve()
NOTES = pathlib.Path(os.environ.get("FIG_NOTES") or OUT).resolve()
for _p in (OUT, NOTES):
    if _p == REPO or REPO in _p.parents or _p == LEDGER or LEDGER in _p.parents:
        stop(f"{_p} is inside the repository or the ledger — outputs go elsewhere")
NAME = "F6_te"
PREVIEW = "--preview" in sys.argv
RUN_ID = "F6-" + datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ") + ("-preview" if PREVIEW else "")
WHAT = ("v13 표시 — post hoc display · displayed after the verdict · no test · 판정 불변 #34 · Fig. 6 · 그림용 기존 결과 다시 읽기"
        "(저자 결정 2026-11-08 · 방향 S45)")

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.colors import same_color  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402

_spec = importlib.util.spec_from_file_location("figcolors", HERE / "figcolors.py")
FC = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(FC)              # S42 §1.2 colours; its palette assertions run on import

_spec = importlib.util.spec_from_file_location("ledger_check", LEDGER / "ledger_check.py")
LC = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(LC)              # read-only use: same_value · eval_field · rounding_ok · numbers_in · norm_text

ROWS = {r["ID"]: r for r in csv.DictReader((LEDGER / "numbers.csv").open(encoding="utf-8"))}


# ----------------------------------------------------------------------------- logged reads (new in v13)
LOGGED = {}                               # resolved path -> dict(sha256, bytes, raw, what)
LC_LOGGED = []                            # copies logged for the ledger_check.py run on the caption
LC_REFUSED = []                           # copies that reread_log.py refuses by its path rule (not logged; see reads.md)


def sha_bytes(b):
    return hashlib.sha256(b).hexdigest()


def reread_log(path, what):
    """Append one entry to drafts/x14/REREAD_LOG.jsonl with the ledger's reread_log.py (it hashes the file and appends
    under its lock); returns the sha256 it logged."""
    r = subprocess.run([sys.executable, "-I", "-B", str(LEDGER / "reread_log.py"), str(path), f"{WHAT} · {what} · 실행 {RUN_ID}"],
                       capture_output=True, text=True)
    if r.returncode != 0 or not r.stdout.strip():
        stop(f"reread_log.py {path}: {r.stdout.strip()} {r.stderr.strip()}")
    return r.stdout.strip().splitlines()[-1].strip()


def logged_read(path, what):
    """A ledger copy or a derived file: logged first, then read once in this run; the bytes used are the bytes logged."""
    path = pathlib.Path(path).resolve()
    if "source_data" in path.parts:
        stop(f"{path}: nothing under source_data is opened")
    if path not in LOGGED:
        h = reread_log(path, what)
        raw = path.read_bytes()
        if sha_bytes(raw) != h:
            stop(f"{path} changed between the logged hash and the read")
        LOGGED[path] = dict(sha256=h, bytes=len(raw), raw=raw, what=what)
    return LOGGED[path]["raw"]


def _logged_load_src(src_dir, rel, sha=""):
    """ledger_check.load_src with the read logged (eval_field reads the ledger copies through it)."""
    base = pathlib.Path(src_dir) / rel.replace("../", "__up__/")
    p = base / sha[:12] if sha and (base / sha[:12]).is_file() else base
    return p, logged_read(p, f"장부 사본 {rel} (장부 행 대조)")


LC.load_src = _logged_load_src            # this process only; ledger_check.py itself is not changed


def log_ledger_check_copies():
    """The copies ledger_check.py --strict-all 1 compares (every ledger row's source copy), counted. Not logged: under the read
    rule of 2026-11-08 (drafts/x15/v12/READ_POLICY_2026-11-08.md; direction S46 V1) ledger_check.py comparing the ledger's copies is not a re-read of results."""
    seen = set()
    for r in ROWS.values():
        if not r["근거 파일"] or r["근거 필드"].startswith(("lit:", "none:")):
            continue
        base = SOURCES / r["근거 파일"].replace("../", "__up__/")
        sha = r["근거 파일 sha256"]
        seen.add((base / sha[:12] if sha and (base / sha[:12]).is_file() else base).resolve())
    return len(seen)


# ----------------------------------------------------------------------------- style (from the repository figstyle.py, inlined; v12 colours)
FONT = "Liberation Sans"
if FONT not in {f.name for f in font_manager.fontManager.ttflist}:
    stop("Liberation Sans is not installed")
PT_TEXT, PT_LETTER = 7.0, 8.5
LW_AXIS, LW_DATA, LW_THIN, LW_MIN = 0.5, 1.0, 0.5, 0.25
INK, INK2, MUTED, GRID, NEUTRAL = FC.INK, FC.INK2, FC.MUTED, FC.GRID, FC.NEUTRAL
AXIS, SURF = "#9a988f", "#ffffff"
MM = 1 / 25.4
WIDTHS_MM = (140.0, 190.0)
MAX_H_MM = 240.0
TAG_SOLID, TAG_DASHED = "-", (0, (2.2, 1.4))
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": [FONT],
    "font.size": PT_TEXT, "axes.labelsize": PT_TEXT, "axes.titlesize": PT_TEXT,
    "xtick.labelsize": PT_TEXT, "ytick.labelsize": PT_TEXT, "legend.fontsize": PT_TEXT,
    "axes.edgecolor": AXIS, "axes.linewidth": LW_AXIS, "axes.labelcolor": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.color": INK2, "ytick.color": INK2,
    "xtick.major.width": LW_AXIS, "ytick.major.width": LW_AXIS,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.minor.size": 0, "ytick.minor.size": 0,
    "xtick.major.pad": 2, "ytick.major.pad": 2, "axes.labelpad": 3,
    "lines.linewidth": LW_DATA, "patch.linewidth": LW_THIN, "hatch.linewidth": 0.5,
    "legend.frameon": False, "mathtext.default": "regular",
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "lines.solid_capstyle": "round", "axes.unicode_minus": True,
})


def new_fig(w_mm, h_mm):
    if not any(abs(w_mm - w) < 1e-6 for w in WIDTHS_MM) or h_mm > MAX_H_MM:
        stop(f"figure size {w_mm} x {h_mm} mm")
    return plt.figure(figsize=(w_mm * MM, h_mm * MM))


def ax_mm(fig, x, y, w, h, **kw):
    """Axes placed in millimetres from the figure's lower-left corner."""
    W, H = fig.get_size_inches() / MM
    return fig.add_axes([x / W, y / H, w / W, h / H], **kw)


def canvas(fig):
    """Full-figure axes in millimetre units (for labels, keys and frames)."""
    W, H = fig.get_size_inches() / MM
    c = fig.add_axes([0, 0, 1, 1], zorder=-1)
    c.set_xlim(0, W)
    c.set_ylim(0, H)
    c.set_xticks([])
    c.set_yticks([])
    c.axis("off")
    return c


def hgrid(ax, axis="y"):
    ax.grid(True, axis=axis, color=GRID, lw=LW_THIN)
    ax.set_axisbelow(True)


def audit(fig, name):
    """Elsevier artwork rules: width, height, text >= 7 pt inside the figure, no sub/superscripts, line widths 0.25–1.5 pt."""
    problems = []
    W, H = fig.get_size_inches() / MM
    if not any(abs(W - w) < 0.05 for w in WIDTHS_MM):
        problems.append(f"width {W:.1f} mm")
    if H > MAX_H_MM + 0.05:
        problems.append(f"height {H:.1f} mm")
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    fb = fig.bbox
    every, drawn_ = tick_label_ids(fig, drawn_only=False), tick_label_ids(fig, drawn_only=True)
    for t in fig.findobj(matplotlib.text.Text):
        s = t.get_text()
        if not s or not t.get_visible() or (id(t) in every and id(t) not in drawn_):
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
    return problems


def save(fig, name):
    probs = audit(fig, name)
    if probs:
        stop(f"{name}: {len(probs)} artwork rule problem(s); not saved")
    OUT.mkdir(parents=True, exist_ok=True)
    base = OUT / name
    fig.savefig(str(base) + ".pdf", metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(str(base) + ".png", dpi=300)
    plt.close(fig)


# ----------------------------------------------------------------------------- frozen sources + ledger guards (v11)
READ = {}                                 # (rel path, sha256) -> copy path


def frozen(rel, ledger_id):
    """Frozen copy of <rel> whose sha256 is the one recorded in ledger row <ledger_id> (read through the log)."""
    sha = ROWS[ledger_id]["근거 파일 sha256"]
    if ROWS[ledger_id]["근거 파일"] != rel:
        stop(f"{ledger_id} cites {ROWS[ledger_id]['근거 파일']}, not {rel}")
    p = SOURCES / rel.replace("../", "__up__/") / sha[:12]
    b = logged_read(p, f"장부 사본 {rel} (그림 값)")
    if sha_bytes(b) != sha:
        stop(f"sha256 of {p} differs from {ledger_id}")
    READ[(rel, sha)] = p
    return b


USED = {}                                 # ledger ID -> (value used here, what)


def led(i, value, what):
    """Stop unless <value> equals the raw value of ledger row <i> (and the row reproduces from its own source copy)."""
    r = ROWS[i]
    v = LC.eval_field(SOURCES, r)
    if v is not None and not LC.same_value(v, r["원 값"]):
        stop(f"{i}: ledger raw {r['원 값']!r} does not reproduce from its source ({v!r})")
    if not LC.same_value(value, r["원 값"]):
        stop(f"{i} ({what}): stored value {value!r} != ledger raw {r['원 값']!r}")
    ok, why = LC.rounding_ok(r["원 값"], r["원고 표기"], r["반올림 규칙"])
    if ok is False:
        stop(f"{i}: 원고 표기 {r['원고 표기']!r} is not the rounding of {r['원 값']!r} ({why})")
    USED[i] = (value, what)
    return value


def row_ok(i, what):
    """A ledger row used as a reference only: its copy has the recorded sha256, it reproduces, and it rounds as recorded."""
    r = ROWS[i]
    if not r["근거 필드"].startswith(("lit:", "none:")):
        frozen(r["근거 파일"], i)
    return led(i, r["원 값"], what)


def T(i):
    """Text of a visible number = the ledger's 원고 표기 (manuscript notation) (the row must have been checked by led())."""
    if i not in USED:
        stop(f"{i} printed before its value was checked against the ledger")
    return ROWS[i]["원고 표기"]


REG = []                                  # (Text artist, ledger IDs whose 원고 표기 (manuscript notation) it prints)


def say(ax, x, y, s, ids, **kw):
    kw.setdefault("fontsize", PT_TEXT)
    kw.setdefault("color", INK)
    t = ax.text(x, y, s, **kw)
    REG.append((t, tuple(ids)))
    return t


def same(a, b, what, tol=1e-12):
    if a is None or b is None or abs(float(a) - float(b)) > tol:
        stop(f"{what}: {a} != {b}")


DRAWN = []                                # drawn values: panel, element, file, sha12, field, value, IDs


def field_ids(rel, field):
    """Ledger rows that cite exactly this file and field."""
    return [i for i, r in ROWS.items() if r["근거 파일"] == rel and r["근거 필드"] == field]


def drawn(panel, element, rel, field, value, ids=None):
    """Log one drawn value with its source; every ledger row citing the same file and field must agree with it."""
    sha = next((s for (k, s) in READ if k == rel), None)
    if sha is None:
        stop(f"{panel} {element}: {rel} was not read through its frozen copy")
    ids = sorted(set(ids or []) | set(field_ids(rel, field)))
    for i in ids:
        if i not in USED:
            led(i, value, f"{panel}: {element}")
    DRAWN.append(dict(panel=panel, element=element, file=rel, sha=sha[:12], field=field, value=value, ids=ids))
    return value


# ----------------------------------------------------------------------------- colours of this figure (v12)
GREY = FC.OI["grey"]                      # marks without a verdict (participants)
DOT = FC.MUTED                            # the 81 unit dots of the slope panel
TE_BLUE = FC.OI["blue"]                   # TE direction, as in the theory figures (Okabe–Ito blue)
CMS = {"NDI": 4.0, "FISO": 3.7, "ODI": 4.3}

# D3 (author decision 2026-11-04 — 57): the literature direction, as in v11 and v12.
D3 = "57"
DREL = "<project-notes>/x15/ledger/derived/derived.json"
LITS = {"57": dict(dir_id="N-602", dir_rel=DREL, dir_field=("gong", "mean_both", "dir_deg"), rx="N-615", ry="N-616",
                   key="mean_both", t2_in_ids=("N-041", "N-042")),
        "58": dict(dir_id="N-048", dir_rel="decomp/manuscript/recompute/a_values.json", dir_field=("dir_deg",), rx="N-370",
                   ry="N-371", key="p1_in_mean_en", t2_in_ids=("N-041",))}
LIT = LITS[D3]


def dig(obj, path):
    for k in path:
        obj = obj[k]
    return obj


# ----------------------------------------------------------------------------- read stored values (v12 load, the multi-TE part, and load_extra's slopes)
def load():
    d = {}
    jrel = "decomp/results/remeasure/judgements.json"
    J = json.loads(frozen(jrel, "N-016"))
    if sha_bytes(logged_read(REPO / jrel, "저장소 judgements.json — 고정 사본과 sha256 대조만")) != ROWS["N-016"]["근거 파일 sha256"]:
        stop("repository judgements.json differs from the frozen copy")
    W = J["Wongkornchaovalit"]
    mrel = "decomp/remeasure/MEASUREMENT_FIXED.md"
    lines = frozen(mrel, "N-206").decode("utf-8").splitlines()
    m329 = re.search(r"— (\d+\.\d+) 이상이면 «명제와 맞음»", lines[328])          # §2-6 (line 329): reading threshold
    if not m329:
        stop("threshold line of the fixed document changed")
    d["wong_thr_read"] = float(m329.group(1))

    # TE-linearity check (fixed before computing; descriptive), overall and by participant
    wx, wy = W["overall"]["x"], W["overall"]["y"]
    for c, row_r, row_n, row_k in (("x", "N-247", "N-011", "N-012"), ("y", "N-248", "N-362", "N-013")):
        o = W["overall"][c]
        k = o["linear_ratio"] * o["n"]
        if abs(k - round(k)) > 1e-9:
            stop(f"multi-TE {c}: share x n is not a count")
        led(row_r, o["linear_ratio"], f"multi-TE {c}: share of units linear in TE")
        led(row_n, o["n"], f"multi-TE {c}: units")
        led(row_k, round(k), f"multi-TE {c}: units linear in TE")
    wthr = led("N-246", d["wong_thr_read"], "reading threshold, fixed document line 329")
    led("N-363", wthr, "reading threshold (same rule, second row)")
    BYSUB = {("sub-01", "x"): "N-365", ("sub-02", "x"): "N-366", ("sub-03", "x"): "N-364",
             ("sub-01", "y"): "N-367", ("sub-02", "y"): "N-368", ("sub-03", "y"): "N-369"}
    if sorted(W["by_subject"]) != ["sub-01", "sub-02", "sub-03"]:
        stop(f"multi-TE participants: {sorted(W['by_subject'])}")
    parts = {"x": [], "y": []}
    for (sub, c), rid in BYSUB.items():
        e = W["by_subject"][sub][c]
        k = e["linear_ratio"] * e["n_units"]
        if e["n_units"] != 27 or abs(k - round(k)) > 1e-9:
            stop(f"multi-TE {sub} {c}: {e['n_units']} units, share x n {k}")
        led(rid, e["linear_ratio"], f"multi-TE {sub} {c}: share of units linear in TE (participant mark)")
        parts[c].append(dict(sub=sub, share=e["linear_ratio"], n_units=e["n_units"], n_lin=round(k), rid=rid, at=e["linear_ratio"] >= wthr))
    for c in ("x", "y"):                                   # the participants' units add up to the overall units
        if sum(p["n_units"] for p in parts[c]) != W["overall"][c]["n"] or \
           sum(p["n_lin"] for p in parts[c]) != round(W["overall"][c]["linear_ratio"] * W["overall"][c]["n"]):
            stop(f"multi-TE {c}: participant units do not add up to the overall")
    low = min((p for c in ("x", "y") for p in parts[c]), key=lambda p: p["share"])
    below = [(c, p["sub"]) for c in ("x", "y") for p in parts[c] if not p["at"]]
    if low["rid"] != "N-364" or below != [("x", "sub-03")]:
        stop(f"multi-TE: lowest participant share is {low} and the shares below the threshold are {below}, not N-364 alone")
    d.update(wong=dict(x=wx["linear_ratio"], y=wy["linear_ratio"]), wong_thr=wthr, wong_parts=parts,
             wong_n=dict(x=round(wx["linear_ratio"] * wx["n"]), y=round(wy["linear_ratio"] * wy["n"])))

    # per-unit TE slopes and directions (exploratory)
    brel = "data/remeasure/explore/b1/b1_wong_direction.json"
    B = json.loads(frozen(brel, "N-372"))
    units = []
    for sub in sorted(B["by_participant"]):
        e = B["by_participant"][sub]
        if len(e["units"]) != e["n"]:
            stop(f"b1 {sub}: unit count")
        for u in e["units"]:
            same(math.degrees(math.atan2(u["by"], u["bx"])), u["theta"], f"b1 {sub} {u['unit']}: theta = atan2(b_y, b_x)", 1e-9)
            units.append(dict(sub=sub, unit=u["unit"], bx=u["bx"], by=u["by"], theta=u["theta"]))
    th = np.array([u["theta"] for u in units])
    if len(th) != B["overall"]["n"]:
        stop("b1: overall n")
    led("N-011", len(th), "b1: units drawn (= multi-TE units)")
    same(float(np.median(th)), B["overall"]["theta_median"], "b1: median of the drawn directions = stored median")
    led("N-372", B["overall"]["theta_median"], "b1: median direction (median of the unit angles)")
    led("N-206", B["wedge"], "b1: TE wedge used by the stored file")
    led("N-164", len(B["by_participant"]), "multi-TE participants")
    bxs, bys = np.array([u["bx"] for u in units]), np.array([u["by"] for u in units])
    d.update(units=units, theta_med=B["overall"]["theta_median"], n_by_neg=int((bys < 0).sum()), n_bx_neg=int((bxs < 0).sum()),
             in_wedge=B["overall"]["in_wedge"])

    # median TE slopes behind the slope ratios (rv_b_family obs, 1/s at four digits) against the 81 drawn units
    rrel = "data/remeasure/explore/rv/rv_b_family.json"
    obs = json.loads(frozen(rrel, "N-370"))["obs"]
    med_b1 = dict(bx=float(np.median(bxs)) * 1000.0, by=float(np.median(bys)) * 1000.0)         # 1/ms -> 1/s
    for c in ("bx", "by"):
        if abs(med_b1[c] - obs[c]) > 0.0005 + 1e-12:          # obs is the same median, stored to 3 decimals in 1/s
            stop(f"median {c} of the 81 units {med_b1[c]} is not rv obs {obs[c]} to its rounding")
    d.update(obs=obs, med_b1=med_b1, med_angle_deg=math.degrees(math.atan2(med_b1["by"], med_b1["bx"])))

    # literature slope vector (Gong 2020 T2; free-water T2 assumed) — D3 picks the T2 set, the ratios and the direction row
    G2 = json.loads(frozen(DREL, "N-615"))["gong"]          # current derived.json (with the slope ratios)
    G1 = json.loads(frozen(DREL, "N-602"))["gong"]          # earlier frozen copy (the direction row)
    for key in ("mean_both", "p1_in_mean_en"):
        for f in ("T2_in", "T2_en", "T2_iso", "dir_deg"):
            same(G1[key][f], G2[key][f], f"derived.json copies disagree on gong.{key}.{f}")
    g = G2[LIT["key"]]
    t2_in_lit = [float(ROWS[i]["원 값"]) for i in LIT["t2_in_ids"]]
    same(g["T2_in"], sum(t2_in_lit) / len(t2_in_lit), f"T2_in of gong.{LIT['key']} = mean of {LIT['t2_in_ids']}")
    for i in LIT["t2_in_ids"]:
        led(i, float(ROWS[i]["원 값"]), "Gong 2020 intra-neurite T2 (input of the literature vector)")
    for i, v in (("N-044", G2["inputs_ms"]["T2_en"][0]), ("N-045", G2["inputs_ms"]["T2_en"][1]), ("N-043", g["T2_en"])):
        led(i, v, "Gong 2020 extra-neurite T2 (input of the literature vector)")
    led("N-046", g["T2_iso"], "free-water T2 assumed in the literature vector")
    pbx = 1.0 / g["T2_en"] - 1.0 / g["T2_in"]                # per ms (derive.py formula)
    pby = 1.0 / g["T2_in"] - 1.0 / g["T2_iso"]
    dir_doc = json.loads(frozen(LIT["dir_rel"], LIT["dir_id"]))
    dir_deg = led(LIT["dir_id"], dig(dir_doc, LIT["dir_field"]), "literature direction")
    same(math.degrees(math.atan2(pby, pbx)), dir_deg, "atan2 of the literature vector = stored direction", 1e-9)
    rx = led(LIT["rx"], g["slope_ratio_x"], "median x slope / literature x slope")
    ry = led(LIT["ry"], g["slope_ratio_y"], "median y slope / literature y slope")
    same(obs["bx"] / (pbx * 1000.0), rx, "rv obs bx / literature x slope = stored ratio", 1e-12)
    same(obs["by"] / (pby * 1000.0), ry, "rv obs by / literature y slope = stored ratio", 1e-12)
    d.update(lit=dir_deg, pvec=(pbx, pby), lit_T2=(g["T2_in"], g["T2_en"], g["T2_iso"]), ratio=(rx, ry),
             ratio_b1=(med_b1["bx"] / (pbx * 1000.0), med_b1["by"] / (pby * 1000.0)))
    med = (obs["bx"] / 1000.0, obs["by"] / 1000.0)            # 1/s -> per ms: the «+» of the median slopes (v12 load_extra)
    same(med[0], rx * pbx, "median slopes x = ratio x literature x", 1e-12)
    same(med[1], ry * pby, "median slopes y = ratio x literature y", 1e-12)
    d["med_slopes"] = med

    # TE wedge of the fixed tests (relaxation grid of Section 2.2; +180 = increasing TE)
    xrel = "decomp/results/stage2/x1.json"
    q1 = json.loads(frozen(xrel, "N-049"))["X1"]["angle_deg_Q1"]
    lo, hi = q1[0] + 180.0, q1[1] + 180.0
    led("N-049", lo, "TE wedge lower end")
    led("N-050", hi, "TE wedge upper end")
    led("N-374", hi, "TE wedge upper end (second row)")
    if [round(lo, 1), round(hi, 1)] != [float(v) for v in LC.split_parts(ROWS["N-206"]["원 값"])]:
        stop("TE wedge ends do not round to the fixed document's 48.5–72.1")
    n_in = int(((th >= lo) & (th <= hi)).sum())
    if n_in != B["overall"]["in_wedge"]:
        stop(f"b1: {n_in} drawn units in the wedge, stored {B['overall']['in_wedge']}")
    d.update(wedge=(lo, hi))

    # drawn values with their sources (the 81 unit dots are listed as one row in values.md)
    for c, rid in (("x", "N-247"), ("y", "N-248")):
        drawn("a", f"{c}: share of units linear in TE, all units (black, filled)", jrel,
              f"json:Wongkornchaovalit.overall.{c}.linear_ratio", d["wong"][c], ids=[rid])
        for p in parts[c]:
            drawn("a", f"{c}: participant share ({'filled' if p['at'] else 'open'}, gray)", jrel,
                  f"json:Wongkornchaovalit.by_subject.{p['sub']}.{c}.linear_ratio", p["share"], ids=[p["rid"]])
    drawn("a", "fixed threshold line", mrel, "md:329:0.8 이", wthr, ids=["N-246"])
    drawn("b", "median direction (black line)", brel, "json:overall.theta_median", d["theta_med"], ids=["N-372"])
    drawn("b", "TE wedge lower end (°)", xrel, "expr:{json:X1.angle_deg_Q1[0]}+180", lo, ids=["N-049"])
    drawn("b", "TE wedge upper end (°)", xrel, "expr:{json:X1.angle_deg_Q1[1]}+180", hi, ids=["N-374"])
    drawn("b", "median slopes «+», x (1/s; per ms = /1000)", rrel, "json:obs.bx", obs["bx"], ids=[])
    drawn("b", "median slopes «+», y (1/s; per ms = /1000)", rrel, "json:obs.by", obs["by"], ids=[])
    drawn("b", f"literature arrow: T2_in of gong.{LIT['key']} (ms)", DREL, f"json:gong.{LIT['key']}.T2_in", g["T2_in"], ids=[])
    drawn("b", f"literature arrow: T2_en of gong.{LIT['key']} (ms)", DREL, f"json:gong.{LIT['key']}.T2_en", g["T2_en"], ids=[])
    drawn("b", f"literature arrow: T2_iso of gong.{LIT['key']} (ms, assumed)", DREL, f"json:gong.{LIT['key']}.T2_iso", g["T2_iso"], ids=[])
    return d


def say_above(ax, x, s, ids, dy_pt=1.6, **kw):
    """Text just above the top of the axes, at data x."""
    kw.setdefault("fontsize", PT_TEXT)
    kw.setdefault("color", INK2)
    kw.setdefault("ha", "center")
    t = ax.annotate(s, (x, 1.0), xycoords=("data", "axes fraction"), xytext=(0, dy_pt), textcoords="offset points",
                    va="bottom", annotation_clip=False, **kw)
    REG.append((t, tuple(ids)))
    return t


# ----------------------------------------------------------------------------- drawing helpers
INTERNAL = re.compile(r"\bT-[ABC]\b|\bV4\b|\bC[37]\b|\bL0\b|\bF[12]\b|stage-2|grid L1|\bL1\b|ODI \([ab]\)|"
                      r"Wongkornchaovalit linearity|\b[0-9a-f]{7}\b|\b\d{5}_[A-Z0-9]+\b|osf\.io|\bsub-\d+")
TAGS = []                                 # role-tag Text artists (their frames are checked like texts)


def text_x1_mm(fig, t):
    fig.canvas.draw()
    bb = t.get_window_extent(fig.canvas.get_renderer())
    return bb.x1 / fig.dpi * 25.4


def tag(cv, x, y, text, dashed=False):
    t = cv.text(x, y, text, fontsize=PT_TEXT, color=INK2, ha="left", va="baseline",
                bbox=dict(boxstyle="square,pad=0.28", fc="none", ec=INK2, lw=0.5, ls=TAG_DASHED if dashed else TAG_SOLID))
    TAGS.append(t)
    return t


def head(fig, cv, x, y, title, tags, ids=(), letter=None, below=False):
    """Panel title (7 pt; optional bold letter 8.5 pt) with its role tags beside it (fixed words, solid or dashed frame);
    below=True puts the tags on the next line where there is no room beside the title."""
    x0 = x
    if letter:
        cv.text(x, y, f"({letter})", fontsize=PT_LETTER, fontweight="bold", color=INK, ha="left", va="baseline")
        x0 = x + 5.4
    t = say(cv, x0, y, title, ids, ha="left", va="baseline", color=INK)
    if below:
        xt, yt = x0 + 0.7, y - 4.4
    else:
        xt, yt = text_x1_mm(fig, t) + 2.2, y
    for text, dashed in tags:
        tt = tag(cv, xt, yt, text, dashed)
        xt = text_x1_mm(fig, tt) + 3.0


def mk(ax, x, y, filled, color, marker, ms, z=3, mew=1.0, edge=None):
    ax.plot([x], [y], marker, ms=ms, ls="none", color=color, mfc=color if filled else SURF, mec=edge or color, mew=mew, zorder=z,
            clip_on=False)



THR_LW = 0.6                                                    # threshold lines: solid, INK2
FIG_W, FIG_H = 190.0, 70.0


# ----------------------------------------------------------------------------- figure
def draw(d):
    W_MM, H_MM = FIG_W, FIG_H
    fig = new_fig(W_MM, H_MM)
    cv = canvas(fig)
    lay = {}
    TY = H_MM - 4.2

    # ======== (a) TE-linearity check, overall (black; the fixed reading) and by participant (gray)
    head(fig, cv, 1.5, TY, f"TE-linearity check, {T('N-164')} participants", [("DESCRIPTIVE, READING FIXED", False)],
         ids=("N-164",), letter="a")
    B1X, B1Y, B1W, B1H = 36.0, 24.0, 48.0, 30.0
    ax4 = ax_mm(fig, B1X, B1Y, B1W, B1H)
    ax4.set_yticks([0, 1])
    ax4.set_yticklabels(["x (logit NDI)", "y (free-water coordinate)"])
    ax4.tick_params(axis="y", length=0, pad=3)
    ax4.spines["left"].set_visible(False)
    ax4.grid(True, axis="y", color=GRID, lw=LW_THIN)
    ax4.set_axisbelow(True)
    ax4.set_ylim(1.72, -0.45)
    ax4.axvline(d["wong_thr"], color=INK2, lw=THR_LW, zorder=1)
    COORD = {"x": dict(marker="o", ms=CMS["NDI"], rid="N-247"), "y": dict(marker="s", ms=CMS["FISO"], rid="N-248")}
    PART_DY, PART_MS = 0.40, 2.8
    for k, c in enumerate(("x", "y")):
        v, m = d["wong"][c], COORD[c]
        if not v >= d["wong_thr"]:
            stop(f"multi-TE {c}: overall share {v} below the fixed reading — the figure draws it as met")
        mk(ax4, v, k, True, INK, m["marker"], m["ms"] * 1.1, mew=0.9)
        say(ax4, v + 0.011, k, T(m["rid"]), (m["rid"],), ha="left", va="center")
        for p in d["wong_parts"][c]:
            mk(ax4, p["share"], k + PART_DY, p["at"], GREY, m["marker"], PART_MS, mew=0.8)
    low = min(d["wong_parts"]["x"], key=lambda p: p["share"])
    say(ax4, low["share"] - 0.009, PART_DY, T(low["rid"]), (low["rid"],), ha="right", va="center", color=INK2)
    top_x = max(p["share"] for p in d["wong_parts"]["x"])
    say(ax4, top_x + 0.010, PART_DY, "participants", (), ha="left", va="center", color=INK2)
    ax4.set_xlim(0.70, 1.0)
    ax4.set_xticks([0.7, 0.8, 0.9, 1.0])
    ax4.set_xticklabels(["0.7", "0.8", "0.9", "1"])
    ax4.set_xlabel("share of units linear in TE")
    say_above(ax4, d["wong_thr"], f"fixed threshold {T('N-246')}", ("N-246",), dy_pt=1.6)

    # ======== (b) TE slopes per unit, one scale for both axes
    head(fig, cv, 92.0, TY, "TE slopes and direction", [("EXPLORATORY", True)], letter="b")
    SX0, SY0, SS = 106.0, 11.0, 46.0
    XL, YL = (-0.003, 0.019), (-0.007, 0.015)
    if abs((XL[1] - XL[0]) - (YL[1] - YL[0])) > 1e-12:
        stop("slope panel: the two axes would not share one scale")
    axs = ax_mm(fig, SX0, SY0, SS, SS)
    axs.set_xlim(*XL)
    axs.set_ylim(*YL)
    lo, hi = d["wedge"]
    R = 0.05
    ang = np.radians(np.linspace(lo, hi, 61))
    axs.add_patch(Polygon(np.vstack([[0.0, 0.0], np.column_stack([R * np.cos(ang), R * np.sin(ang)])]), closed=True,
                          facecolor=TE_BLUE, alpha=0.22, edgecolor="none", lw=0, zorder=0.6))
    axs.axhline(0.0, color=AXIS, lw=LW_THIN, zorder=1)
    axs.axvline(0.0, color=AXIS, lw=LW_THIN, zorder=1)
    tm = math.radians(d["theta_med"])
    axs.plot([0.0, R * math.cos(tm)], [0.0, R * math.sin(tm)], color=INK, lw=LW_DATA, zorder=2, solid_capstyle="butt")
    bx = np.array([u["bx"] for u in d["units"]])
    by = np.array([u["by"] for u in d["units"]])
    if not ((XL[0] < bx) & (bx < XL[1]) & (YL[0] < by) & (by < YL[1])).all():
        stop("slope panel: a unit lies outside the axes")
    axs.plot(bx, by, "o", ms=2.6, ls="none", mfc=DOT, mec=SURF, mew=0.3, zorder=3)
    pbx, pby = d["pvec"]
    axs.annotate("", xy=(pbx, pby), xytext=(0.0, 0.0), arrowprops=dict(arrowstyle="-|>", color=TE_BLUE, lw=LW_DATA,
                 shrinkA=0, shrinkB=0, mutation_scale=8), zorder=4)
    mx, my = d["med_slopes"]
    axs.plot([mx], [my], "+", ms=6.5, mew=1.2, color=INK, zorder=5)            # median slopes (b_x, b_y medians)
    axs.set_xticks([0.0, 0.005, 0.01, 0.015])
    axs.set_xticklabels(["0", "0.005", "0.01", "0.015"])
    axs.set_yticks([-0.005, 0.0, 0.005, 0.01, 0.015])
    axs.set_yticklabels(["−0.005", "0", "0.005", "0.01", "0.015"])
    axs.set_xlabel("TE slope of x (per ms)")
    axs.set_ylabel("TE slope of y (per ms)")
    RX = SX0 + SS + 1.6

    def ymm(yv):
        return SY0 + (yv - YL[0]) / (YL[1] - YL[0]) * SS
    say(cv, RX, ymm(0.0135), f"TE wedge {T('N-206')}", ("N-206",), ha="left", va="center", color=TE_BLUE)
    say(cv, RX, ymm(pby), f"literature T2, {T(LIT['dir_id'])}", (LIT["dir_id"],), ha="left", va="center", color=TE_BLUE)
    y_exit = XL[1] * math.tan(tm)
    if not (YL[0] < y_exit < YL[1]):
        stop("slope panel: the median ray does not leave through the right edge")
    say(cv, RX, ymm(y_exit), f"median direction {T('N-372')}", ("N-372",), ha="left", va="center")
    cv.plot([RX + 1.6], [ymm(my) - 4.6], "+", ms=6.5, mew=1.2, color=INK)
    say(cv, RX + 4.0, ymm(my) - 4.6, "median slopes, as multiples", (), ha="left", va="center")
    say(cv, RX + 4.0, ymm(my) - 7.8, f"of the arrow: x {T(LIT['rx'])}, y {T(LIT['ry'])}", (LIT["rx"], LIT["ry"]), ha="left", va="center")

    # one scale for both axes, measured on the drawn figure
    fig.canvas.draw()
    o_ = axs.transData.transform((0.0, 0.0))
    ux = (axs.transData.transform((0.01, 0.0)) - o_)[0]
    uy = (axs.transData.transform((0.0, 0.01)) - o_)[1]
    if abs(ux - uy) > 1e-6 * ux:
        stop(f"slope panel: {ux} px vs {uy} px per 0.01 per ms")
    lay["slope_mm_per_0.01_per_ms"] = round(ux / fig.dpi * 25.4, 4)
    lay["units_drawn"] = len(bx)
    d["layout"] = lay
    return fig


def tick_label_ids(fig, drawn_only):
    ids = set()
    for a in fig.axes:
        for axis in (a.xaxis, a.yaxis):
            ticks = axis._update_ticks() if drawn_only else axis.get_major_ticks() + axis.get_minor_ticks()
            if drawn_only and not a.axison:
                continue
            for tk in ticks:
                ids.update({id(tk.label1), id(tk.label2)})
    return ids


def check_visible_numbers(fig):
    """Every number in a visible non-tick text must be one of the ledger 원고 표기 (manuscript notation) registered for that text."""
    tick = tick_label_ids(fig, drawn_only=False)
    reg = {id(t): ids for t, ids in REG}
    shown = []
    for t in fig.findobj(matplotlib.text.Text):
        s = t.get_text()
        if not s or not t.get_visible() or id(t) in tick:
            continue
        nums = LC.numbers_in(s)
        ids = reg.get(id(t), ())
        allowed = [v for i in ids for v in LC.numbers_in(ROWS[i]["원고 표기"])]
        for v in nums:
            hit = [k for k, a in enumerate(allowed) if math.isclose(a, v, rel_tol=0, abs_tol=1e-12)]
            if not hit:
                stop(f"visible number {v} in {s!r} is not a registered ledger 원고 표기 ({ids})")
            allowed.pop(hit[0])
        for i in ids:
            if LC.norm_text(ROWS[i]["원고 표기"]) not in LC.norm_text(s):
                stop(f"text {s!r} does not contain the 원고 표기 of {i}")
        if nums:
            shown.append((s.replace("\n", " "), ids))
    return shown


def box_poly(bb):
    return np.array([[bb.x0, bb.y0], [bb.x1, bb.y0], [bb.x1, bb.y1], [bb.x0, bb.y1]], float)


def text_poly(t, rnd):
    """Display-space corners of a text's box: its window extent when the rotation is a multiple of 90°, else the unrotated
    box (same alignment) turned about the anchor (rotation_mode 'anchor'), so that 45° labels are tested as drawn."""
    rot = t.get_rotation() % 360.0
    if abs(rot % 90.0) < 1e-9:
        return box_poly(t.get_window_extent(rnd))
    if t.get_rotation_mode() != "anchor":
        stop(f"rotated text {t.get_text()!r} without rotation_mode 'anchor'")
    t.set_rotation(0.0)
    ub = t.get_window_extent(rnd)
    t.set_rotation(rot)
    ax_, ay_ = t.get_transform().transform(t.get_unitless_position())
    c, s = math.cos(math.radians(rot)), math.sin(math.radians(rot))
    return (box_poly(ub) - [ax_, ay_]) @ np.array([[c, s], [-s, c]]) + [ax_, ay_]


def sat_gap(P, Q):
    """Separation of two convex polygons along their edge normals (separating-axis test): negative = overlap depth. For two
    axis-parallel boxes this is max(−ix, −iy) of the v11 test."""
    best = -np.inf
    for poly in (P, Q):
        for k in range(len(poly)):
            e = poly[(k + 1) % len(poly)] - poly[k]
            nn = math.hypot(e[0], e[1])
            if nn == 0:
                continue
            nrm = np.array([-e[1], e[0]]) / nn
            p, q = P @ nrm, Q @ nrm
            best = max(best, q.min() - p.max(), p.min() - q.max())
    return best


def check_layout(fig, d):
    """No two drawn texts (role-tag frames as frames) overlap, no text covers a data mark (tolerance 0.2 mm), every text is
    7 pt or more, and no figure text names an internal test, arm, hash, session or participant. Rotated texts (the 45°
    column names) are tested as rotated boxes; unrotated texts exactly as in v11."""
    fig.canvas.draw()
    rnd = fig.canvas.get_renderer()
    every, drawn_ = tick_label_ids(fig, drawn_only=False), tick_label_ids(fig, drawn_only=True)
    px = fig.dpi / 25.4
    tol = 0.2 * px
    tag_ids = {id(t) for t in TAGS}
    texts = []
    for t in fig.findobj(matplotlib.text.Text):
        s = t.get_text()
        if not s or not t.get_visible() or (id(t) in every and id(t) not in drawn_):
            continue
        if INTERNAL.search(s):
            stop(f"figure text {s!r} names an internal test, arm, hash or session")
        if t.get_fontsize() < PT_TEXT - 1e-6:
            stop(f"figure text {s!r} at {t.get_fontsize()} pt")
        bb = t.get_bbox_patch().get_window_extent(rnd) if id(t) in tag_ids else t.get_window_extent(rnd)
        if bb.width > 0 and bb.height > 0:
            poly = box_poly(bb) if id(t) in tag_ids else text_poly(t, rnd)
            texts.append((s.replace("\n", " "), poly))
    clash = []
    gaps = []
    for i in range(len(texts)):
        for j in range(i + 1, len(texts)):
            g = sat_gap(texts[i][1], texts[j][1])
            if g < -tol:
                clash.append((texts[i][0], texts[j][0]))
            gaps.append((g / px, texts[i][0], texts[j][0]))
    marks = []
    for ln in fig.findobj(matplotlib.lines.Line2D):
        if not ln.get_visible() or ln.get_marker() not in ("o", "s", "^", "+"):
            continue
        if same_color(ln.get_markerfacecolor(), SURF) and same_color(ln.get_markeredgecolor(), SURF):
            continue                                            # surface rings
        xy = ln.get_transform().transform(np.column_stack([np.asarray(ln.get_xdata(), float), np.asarray(ln.get_ydata(), float)]))
        r = ln.get_markersize() / 2.0 * fig.dpi / 72.0
        marks += [(x, y, r) for x, y in xy]
    covered = []
    for s, poly in texts:
        for x, y, r in marks:
            if sat_gap(poly, np.array([[x - r, y - r], [x + r, y - r], [x + r, y + r], [x - r, y + r]])) < -tol:
                covered.append((s, round(x / px, 1), round(y / px, 1)))
    if clash or covered:
        stop(f"layout: overlapping texts {clash[:6]} · texts over marks {covered[:6]}")
    gaps.sort()
    d["layout"].update(n_texts=len(texts), n_marks=len(marks), nearest_text_gap_mm=round(gaps[0][0], 2),
                       nearest_texts=(gaps[0][1], gaps[0][2]))


def esc(x):
    """Markdown table cell: escape the column separator."""
    return str(x).replace("|", "\\|")


# ----------------------------------------------------------------------------- caption
# Panels and symbols only (S45 §7; S41 §3): the observations (0.815 and 0.901, 25.4°, 1.12 and 0.31) are printed in the
# figure and given in the body.
CAPTION = ('![**Fig. 6. Across echo times (multi-TE series, {N-164} participants).** (a) Share of ROI units linear in TE against '
           'the threshold of {N-246} fixed before computing (descriptive); black, all units; gray, each participant; open, below '
           'the threshold. (b) Per-unit TE slopes on equal axes (exploratory); shaded, TE wedge (Section 2.2); arrow, slopes from '
           'literature T2; line, median direction; +, median slopes.](figures/F6_te.pdf){#fig:fig6 .wide width=100%}')
NUMBER_WORDS = sorted(set(LC.WORDS) | {"single", "once", "twice"}, key=len, reverse=True)


def caption():
    def sub(m):
        up, i = m.group(1), m.group(2)
        if i not in ROWS:
            stop(f"caption: {i} is not in the ledger")
        r = ROWS[i]
        v = LC.eval_field(SOURCES, r)
        if v is not None and not LC.same_value(v, r["원 값"]):
            stop(f"caption: {i} does not reproduce from its source")
        if i not in USED:
            stop(f"caption: {i} was not checked against the stored value")
        t = r["원고 표기"]
        return f"{t[:1].upper() + t[1:] if up else t} [{i}]"
    cap_text = re.sub(r"\{(\^?)(N-\d{3})\}", sub, CAPTION)
    spans = []
    for m in LC.TAG.finditer(cap_text):
        end = m.start()
        while end > 0 and cap_text[end - 1] in "  ":
            end -= 1
        st = LC.suffix_start(cap_text, end, ROWS[m.group(1)]["원고 표기"])
        if st is None:
            stop(f"caption: text before {m.group(1)} is not its ledger text")
        spans.append((st, m.end()))
    masked = list(cap_text)
    for a, b in spans:
        masked[a:b] = " " * (b - a)
    masked = "".join(masked)
    for ex in LC.EXCEPTIONS:
        masked = ex.sub(lambda mm: re.sub(r"[^\n]", " ", mm.group(0)), masked)
    hits = [m.group(0) for w in NUMBER_WORDS for m in re.finditer(rf"(?<!\w){re.escape(w)}(?!\w)", masked, re.I)]
    if hits:
        stop(f"caption: number words without a ledger tag: {hits}")
    forb = [m.group(0) for rx, _ in LC.FORBIDDEN for m in rx.finditer(cap_text)]
    if forb:
        stop(f"caption uses words not to use: {forb}")
    if INTERNAL.search(re.sub(r"\[N-\d{3}\]", "", cap_text)):
        stop(f"caption names an internal test: {INTERNAL.search(cap_text).group(0)}")
    body = re.sub(r"\[N-\d{3}\]", "", re.sub(r"^!\[|\]\(figures/.*$", "", cap_text))
    return cap_text, len(re.findall(r"\S+", body))


def ledger_check(path):
    n = log_ledger_check_copies()                              # every copy ledger_check.py reads, logged before it runs
    r = subprocess.run([sys.executable, "-I", "-B", str(LEDGER / "ledger_check.py"), str(path), "--strict-all", "1"],
                       capture_output=True, text=True)
    summary = r.stdout.strip().splitlines()[0] if r.stdout.strip() else r.stderr.strip()
    if r.returncode != 0:
        stop(f"ledger_check on {path.name}: {r.stdout.strip()} {r.stderr.strip()}")
    counts = {k: int(v) for k, v in re.findall(r"(실패|엄격 구역|쓰지 않을 말) (\d+)", summary)}
    if counts.get("실패") != 0 or counts.get("엄격 구역") != 0 or counts.get("쓰지 않을 말") != 0:
        stop(f"ledger_check summary not clean: {summary}")
    return summary, n


# ----------------------------------------------------------------------------- reads.md (log lines of this run)
def log_lines():
    """(line, entry) of this run's entries in drafts/x14/REREAD_LOG.jsonl (matched by RUN_ID)."""
    out = []
    for k, s in enumerate(REREAD_LOG.read_text(encoding="utf-8").splitlines(), 1):
        if RUN_ID in s:
            try:
                e = json.loads(s)
            except json.JSONDecodeError:
                continue
            if e.get("무엇", "").endswith(f"실행 {RUN_ID}"):
                out.append((k, e))
    return out


def write_reads(path):
    ents = log_lines()
    fig_ents = [(k, e) for k, e in ents if "ledger_check.py --strict-all 1 가 읽는" not in e["무엇"]]
    lc_ents = [(k, e) for k, e in ents if "ledger_check.py --strict-all 1 가 읽는" in e["무엇"]]
    for p, info in LOGGED.items():
        if not any(e["다시읽음"] == str(p) and e["sha256"] == info["sha256"] for _, e in fig_ents):
            stop(f"reads.md: {p} has no log entry of this run")
    if len(lc_ents) != len(LC_LOGGED):
        stop(f"reads.md: {len(lc_ents)} ledger_check entries in the log, {len(LC_LOGGED)} logged")
    if not PREVIEW and len(LC_LOGGED) + len(LC_REFUSED) != len({str(x["path"]) for x in LC_LOGGED + LC_REFUSED}):
        stop("reads.md: a ledger_check copy counted twice")
    L = [f"# Fig. 6 (v13) — reads of run {RUN_ID} (written by make_F6.py)", "",
         "Every ledger copy and derived file below was logged with `drafts/x15/ledger/reread_log.py` in "
         "`drafts/x14/REREAD_LOG.jsonl` before it was read (the script reads the bytes after the log entry and stops if their "
         "sha256 differs from the logged one). No image and nothing under `source_data/` was opened.", "",
         f"## Read by the figure script ({len(fig_ents)} entries)", "",
         "| Log line | Time | File | sha256 (12) | Bytes | What |", "|---|---|---|---|---|---|"]
    for k, e in fig_ents:
        L.append(f"| {k} | {e['시각']} | `{os.path.relpath(e['다시읽음'], WORK)}` | `{e['sha256'][:12]}` | {e['바이트']} | "
                 f"{esc(e['무엇'].split(' · ')[-2])} |")
    L += ["", "## Ledger copies compared by ledger_check.py --strict-all 1 on caption.md", "",
          "Not logged: under the read rule of 2026-11-08 (`drafts/x15/v12/READ_POLICY_2026-11-08.md`; direction S46 V1) ledger_check.py comparing "
          "the ledger's copies is not a re-read of results.", ""]
    L += ["", "## Read without a log entry (not derived results)", "",
          "`drafts/x15/ledger/numbers.csv` (ledger rows), `drafts/x15/ledger/ledger_check.py` (imported read-only), "
          f"`figcolors.py`, and this run's own `{NAME}.png` (for the colour-vision sheet)."]
    path.write_text("\n".join(L) + "\n", encoding="utf-8")
    return len(fig_ents), len(lc_ents), (fig_ents[0][0], fig_ents[-1][0]) if fig_ents else None


# ----------------------------------------------------------------------------- values file
def write_values(d, shown, cap_text, n_words, check_line, path):
    L = []
    w = L.append
    pbx, pby = d["pvec"]
    w("# F6 (v13) — drawn values and visible numbers (written by make_F6.py)")
    w("")
    w("Every value below was read from a frozen ledger source copy (`drafts/x15/ledger/sources/<file>/<sha12>`, sha12 = first 12 "
      "hex of the sha256 recorded in the ledger; the script stops if the copy's sha256 differs; every read logged first, see "
      "reads.md) and drawn as is, except the literature slope vector (computed from stored T2, guarded as in v11 and v12: its "
      "direction = the stored dir_deg, and the stored slope ratios times it = the stored median slopes). «no ledger row (drawn "
      "only)» means the value is drawn but not printed and no ledger row cites that field.")
    w("")
    w("## Sources read (frozen ledger copies; sha256 checked)")
    w("")
    w("| Source (ledger 근거 파일) | Frozen copy | sha256 |")
    w("|---|---|---|")
    for (k, sha), p in READ.items():
        w(f"| `{k}` | `{p.relative_to(SOURCES)}` | `{sha}` |")
    w("")
    w("Also read (logged): `decomp/results/remeasure/judgements.json` in the repository, for its sha256 only (equal to the frozen "
      "copy), and the ledger copies that the ledger rows asserted below cite (eval_field). Read without a log entry: "
      "`drafts/x15/ledger/numbers.csv`; `drafts/x15/ledger/ledger_check.py` (imported read-only); `figcolors.py`.")
    w("")
    w("## (a) TE-linearity check")
    w("")
    w(f"All units: x {d['wong']['x']!r} (N-247; {d['wong_n']['x']} of 81), y {d['wong']['y']!r} (N-248; {d['wong_n']['y']} of 81); "
      "participants: " + "; ".join(f"{c} {p['sub']} {p['share']!r} ({p['n_lin']}/27, {p['rid']}, {'filled' if p['at'] else 'open'})"
                                   for c in ("x", "y") for p in d["wong_parts"][c]) +
      f"; fixed threshold {d['wong_thr']} (N-246, line 329 of the fixed document). Participant IDs appear here and in the "
      "reading records only, not in the figure or the caption.")
    w("")
    w("## (b) TE slopes and direction")
    w("")
    w(f"81 unit dots (b1 `by_participant.*.units[*].bx · by`, per ms; {d['in_wedge']} inside the wedge, {d['n_by_neg']} with "
      f"b_y < 0, {d['n_bx_neg']} with b_x < 0 — counts for this file only, not printed); median direction {d['theta_med']!r}° "
      f"(N-372); literature arrow ({pbx!r}, {pby!r}) per ms from `gong.{LIT['key']}` T2 {d['lit_T2']} ms, direction "
      f"{d['lit']!r}° ({LIT['dir_id']}); «+» of the median slopes = rv obs ({d['obs']['bx']}, {d['obs']['by']}) per s = ratios "
      f"{d['ratio'][0]!r} · {d['ratio'][1]!r} ({LIT['rx']} · {LIT['ry']}) times the arrow (the median of the 81 drawn units: "
      f"{d['med_b1']['bx']!r} · {d['med_b1']['by']!r} per s, within the rounding of rv obs); TE wedge {d['wedge'][0]!r}° – "
      f"{d['wedge'][1]!r}° (N-049 · N-050; label N-206); one scale for both axes ({d['layout']['slope_mm_per_0.01_per_ms']} mm "
      "per 0.01 per ms, measured on the drawn figure).")
    w("")
    w("## Every drawn value with its source")
    w("")
    w("| Panel | Element | File | sha12 | Field | Value | Ledger ID |")
    w("|---|---|---|---|---|---|---|")
    for e in DRAWN:
        w(f"| {e['panel']} | {esc(e['element'])} | `{e['file']}` | `{e['sha']}` | `{esc(e['field'])}` | {esc(repr(e['value']))} | "
          f"{', '.join(e['ids']) if e['ids'] else 'no ledger row (drawn only)'} |")
    w(f"| b | 81 unit dots (b_x, b_y), listed in the stored order | `data/remeasure/explore/b1/b1_wong_direction.json` | "
      f"`{ROWS['N-372']['근거 파일 sha256'][:12]}` | `json:by_participant.*.units[*].bx · by` | 81 pairs | "
      "no ledger row (drawn only); count N-011, median direction N-372 |")
    w("")
    w("## Visible numbers (figure text other than tick labels)")
    w("")
    w("| Printed text | Ledger ID → 원고 표기 · raw value |")
    w("|---|---|")
    for s, ids in shown:
        w(f"| {esc(s)} | " + " · ".join(f"{i} → «{esc(ROWS[i]['원고 표기'])}» · {esc(ROWS[i]['원 값'])}" for i in ids) + " |")
    w("")
    w("Tick labels (axis scales) are exempt.")
    w("")
    w("## Ledger rows asserted by the script")
    w("")
    w("| ID | 원고 표기 | raw value (ledger) | value read here | what |")
    w("|---|---|---|---|---|")
    for i, (v, what) in USED.items():
        w(f"| {i} | {esc(ROWS[i]['원고 표기'])} | {esc(ROWS[i]['원 값'])} | {esc(repr(v))} | {esc(what)} |")
    w("")
    w("## Unfavourable results and where they stand")
    w("")
    w("| result | role | place |")
    w("|---|---|---|")
    w(f"| median TE direction {T('N-372')} outside the TE wedge {T('N-206')} (literature {T(LIT['dir_id'])}) | exploratory | "
      "(b): black line and its label, the shaded wedge and the blue arrow |")
    w(f"| y slope {T(LIT['ry'])} of the literature prediction | exploratory | (b): «+» label |")
    w(f"| lowest participant share {T('N-364')} below the reading threshold | descriptive | (a): labelled open gray mark |")
    w("")
    w(f"## Caption ({n_words} words in the figure line without IDs) and its ledger check")
    w("")
    w(cap_text)
    w("")
    w(f"    {check_line}")
    w("")
    w("## Layout check")
    w("")
    lay = d["layout"]
    w(f"{lay['n_texts']} drawn texts (tag frames as frames), none overlapping another and none over any of the {lay['n_marks']} "
      f"data marks (tolerance 0.2 mm); nearest two texts {lay['nearest_text_gap_mm']} mm apart ({lay['nearest_texts'][0]!r} · "
      f"{lay['nearest_texts'][1]!r}).")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


# ----------------------------------------------------------------------------- main
PYCACHES = (REPO / "decomp/manuscript/figures/src/__pycache__", LEDGER / "__pycache__", HERE / "__pycache__")
PYCACHE_BEFORE = {p for p in PYCACHES if p.exists()}        # this folder is shared with the other figure scripts


def main():
    d = load()
    cap_text, n_words = caption()
    fig = draw(d)
    shown = check_visible_numbers(fig)
    check_layout(fig, d)
    print("layout", d["layout"])
    print("visible numbers:")
    for s, ids in shown:
        print("  ", s, "->", ", ".join(f"{i} «{ROWS[i]['원고 표기']}»" for i in ids))
    if PREVIEW:                                                # layout check only: PNG to the given path, nothing else written
        probs = audit(fig, NAME)
        fig.savefig(sys.argv[sys.argv.index("--preview") + 1], dpi=300)
        print("preview written; audit problems:", probs, "; logged reads this run:", len(LOGGED))
        return
    save(fig, NAME)
    NOTES.mkdir(parents=True, exist_ok=True)
    cap_path = NOTES / "caption.md"
    cap_path.write_text(cap_text + "\n", encoding="utf-8")
    check_line, n_lc = ledger_check(cap_path)
    print("ledger check (caption):", check_line, f"(copies it compares: {n_lc}; not logged, read rule of 2026-11-08) "
          f"and not logged: {len(LC_REFUSED)} {[os.path.relpath(x['path'], WORK) for x in LC_REFUSED]})")
    write_values(d, shown, cap_text, n_words, check_line, NOTES / "values.md")
    n_fig, n_lcl, rng = write_reads(NOTES / "reads.md")
    cvd = FC.cvd_sheet(OUT / f"{NAME}.png", OUT / "cvd" / "F6_cvd.png", title="Fig. 6 (v13): original and three colour-vision simulations")
    for d_ in PYCACHES:
        if d_.exists() and d_ not in PYCACHE_BEFORE:
            stop(f"{d_} appeared during the run — a write outside FIG_OUT")
    for f in (OUT / f"{NAME}.pdf", OUT / f"{NAME}.png", cvd, cap_path, NOTES / "values.md", NOTES / "reads.md"):
        print(f"sha256 {f.name} {sha_bytes(f.read_bytes())}")
    print(f"read log {REREAD_LOG.name}: run {RUN_ID}: {n_fig} figure reads (lines {rng[0]}–{rng[1]}), {n_lcl} ledger_check copies")
    print(f"caption words (figure line without IDs): {n_words}")
    print("ledger rows asserted:", len(USED), "·", " ".join(USED))
    print(f"AUDIT PASS {NAME}: {FIG_W:.0f} x {FIG_H:.0f} mm, font {FONT}, guards and layout checks passed")


if __name__ == "__main__":
    main()
