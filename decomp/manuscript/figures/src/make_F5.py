"""Figure 5 of manuscript v13 (Results 4.2): across scanner models and protocols, one comparison per figure (review direction
S45 §7, row 5). Drawn from stored values only; nothing is fitted, estimated or tested here.

(a) ODI against retest — DESCRIPTIVE: the 19 ON-Harmony scanner pairs with common participants, the upper triangle of the
    v12 matrix only. Rows = columns = the 7 scanners ordered by TE (the shorter TE for NOT2ING; ties by name); one name
    column serves the rows, the column names stand above the matrix. Cell = share of ROI units at retest level, sequential
    blue (figcolors.SEQ, 0–1 linear), outlined (white inner line and black border) at or above the fixed 0.90; hatched = no
    common participant; diagonal = the same scanner (neutral). Left of the matrix one track only: the TE of each scanner
    (dots on a TE axis; NOT2ING: its two protocols joined by a segment). Vendor, TR, receive coil and direction scheme are
    not drawn here (they stay in the judge matrix of Fig. 7).
(b) FISO sign and magnitude order — UNBLINDED CHECK: the 8 pairs with a predicted sign; predicted (largest on the primary
    relaxation grid, fixed table) against observed |Δ logit FISO|; circle in figcolors.AGREE where the observed sign agreed
    with the predicted sign, square in figcolors.OPPOSE where it was opposed (shape and colour; the stored sign results, no new
    test). No statistic is printed in the figure (Kendall τb and the ODI line of v12 are in the body).

v13 (2026-11-08): rewritten in place from the v12 script of this folder (sha256 eafc3960fc56803f6adf8ec74829e7e4bf2b160fd933eb631f67578b6504ce84,
copy in the v13 scratch fig5/), which drew v12 Fig. 5 (a) matrix with both triangles and five acquisition tracks, the FISO
scatter with τb, and (b) the multi-TE series. The multi-TE series is now Fig. 6 (make_F6.py, same folder). Kept unchanged:
the readers and guards of the fixed prediction table, the stored pair values, the sign test and τb (recomputed from the drawn
points = stored), the stored «not supported» readings, ODI against retest (share, need, flags, count and names), the TE of
the acquisitions (plan.py ONH + the 28 stage-3 cards; sidecar TE rows; the registered scanner table's TE column), the matrix
cells against the fixed table, the visible-number, layout, artwork and caption checks. Dropped with the elements they
guarded: vendor (access-record Manufacturer summary), TR, receive coil (rv_4a_coil.json, no longer read), direction scheme,
the lower-triangle size rule, ODI against the predicted shift (τb line), the prototype comparison (--proto-diff).
New in v13: every ledger copy and derived file is logged with the ledger's reread_log.py (drafts/x14/REREAD_LOG.jsonl)
before it is read (approved re-read of existing results for the figures, author decision 2026-11-08, direction S45), and reads.md lists the log lines of the run.
The copies ledger_check.py compares for the caption check are not logged (read rule of 2026-11-08, drafts/x15/v12/READ_POLICY_2026-11-08.md).
Writes, into FIG_OUT only: F5_scanners.pdf/.png (audit must pass) and cvd/F5_cvd.png; into FIG_NOTES (default FIG_OUT):
caption.md (checked with ledger_check.py --strict-all 1), values.md, reads.md.
Run (CPU only), from any folder ([--preview <png path>] for a layout check only; nothing else is written then):
  PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=<scratch> CUDA_VISIBLE_DEVICES= DMRI_REPO=<repo> DMRI_WORK=<data> \
    FIG_OUT=<out> FIG_NOTES=<notes> nice -n 10 <repo>/.venv/bin/python -I -B make_F5.py
"""
import sys

sys.dont_write_bytecode = True            # no __pycache__ beside this file, in the repository or in the ledger folder

import ast  # noqa: E402
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
    raise SystemExit(f"STOP F5: {msg}")


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
NAME = "F5_scanners"
PREVIEW = "--preview" in sys.argv
RUN_ID = "F5-" + datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ") + ("-preview" if PREVIEW else "")
WHAT = ("v13 표시 — post hoc display · displayed after the verdict · no test · 판정 불변 #34 · Fig. 5 · 그림용 기존 결과 다시 읽기"
        "(저자 결정 2026-11-08 · 방향 S45)")

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.colors import Normalize, same_color  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from scipy.stats import kendalltau  # noqa: E402

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


# ----------------------------------------------------------------------------- read stored values (v12 load, the ON-Harmony part)
def load():
    d = {}
    # judgements.json (repository; frozen copy = same sha256 as the repository file)
    jrel = "decomp/results/remeasure/judgements.json"
    J = json.loads(frozen(jrel, "N-016"))
    if sha_bytes(logged_read(REPO / jrel, "저장소 judgements.json — 고정 사본과 sha256 대조만")) != ROWS["N-016"]["근거 파일 sha256"]:
        stop("repository judgements.json differs from the frozen copy")
    O = J["ON-Harmony 쌍 검정"]

    # fixed prediction table (MEASUREMENT_FIXED.md §2-1)
    mrel = "decomp/remeasure/MEASUREMENT_FIXED.md"
    lines = frozen(mrel, "N-206").decode("utf-8").splitlines()
    head_ = [k for k, s in enumerate(lines) if s.startswith("**장비 단위 21 쌍 — 부호 예측 · 예측 크기")]
    if len(head_) != 1 or not lines[head_[0] + 2].startswith("| 장비 쌍(21) | 공통 피험자 |"):
        stop("fixed prediction table not found where expected")
    table = []
    for k in range(head_[0] + 4, head_[0] + 4 + 21):
        c = [x.strip() for x in lines[k].split("|")[1:-1]]
        if len(c) != 8:
            stop(f"fixed table line {k + 1}: {lines[k][:80]}")
        cell = c[5].split("·")
        if len(cell) != 2 or not all(re.fullmatch(r"\d\.\d{4}", x.strip()) for x in cell):
            stop(f"fixed table line {k + 1}: L1 · L4 cell {c[5]!r} is not two four-decimal values")
        l1, l4 = (float(x) for x in cell)
        sign = {"+": 1, "−": -1, "—": 0}[c[6]]
        table.append(dict(line=k + 1, pair=c[0], n_common=int(c[1]), L1=l1, L4=l4, sign=sign, why=c[7], cell_L=c[5], cell_sign=c[6]))
    if lines[head_[0] + 4 + 21].startswith("| ") or len({r["pair"] for r in table}) != 21:
        stop("fixed prediction table does not have exactly 21 rows")
    d["table_lines"] = (table[0]["line"], table[-1]["line"])
    if "부호 예측 8 쌍(+ 5 · − 3)" not in lines[91]:          # R1 (line 92): 8 pairs, + 5 · − 3
        stop("R1 line of the fixed document changed")
    n_pos, n_neg = sum(r["sign"] > 0 for r in table), sum(r["sign"] < 0 for r in table)
    led("N-219", [f"+ {n_pos}", f"− {n_neg}"], "fixed table: predicted signs + · −")
    m97 = re.search(r"단위 중 (\d+\.\d+) 이상이면 «재검사 수준»\((\d+) 단위가 모두 있으면 (\d+) 단위 이상\)", lines[96])   # R6 (line 97)
    if not m97:
        stop("threshold line of the fixed document changed")
    d.update(odi_thr_read=float(m97.group(1)), odi_units_read=int(m97.group(2)), odi_need_read=int(m97.group(3)))

    # observed pair values (D_j − D_i for the stored name D_i~D_j)
    pj = {p["pair"]: p for p in O["pairs"]}
    if len(pj) != 21:
        stop("judgements: not 21 pairs")
    sign_rows = []
    for r in table:
        a, b = r["pair"].split("~")
        keys = [k for k in (f"{a}~{b}", f"{b}~{a}") if k in pj]
        if len(keys) != 1:
            stop(f"pair {r['pair']} not found once in judgements")
        p = pj[keys[0]]
        flip = keys[0] != f"{a}~{b}"
        if (r["n_common"] >= 3) != p["judgeable_fiso"] or p["judgeable_fiso"] != p["judgeable_odi"]:
            stop(f"{r['pair']}: judgeable flag disagrees with the fixed table")
        r.update(key=keys[0], flip=flip)
        if r["sign"] == 0:
            continue
        if not p["judgeable_fiso"] or p["fiso"] is None:
            stop(f"{r['pair']}: sign-predicted pair without an observed value")
        obs_ba = -p["fiso"] if flip else p["fiso"]            # observed B − A for the table's «A~B»
        hi, lo = (b, a) if r["sign"] > 0 else (a, b)         # name the row so that the predicted difference is positive
        sign_rows.append(dict(r, label=f"{hi} − {lo}", obs_ba=obs_ba, obs_pos=r["sign"] * obs_ba, abs_obs=abs(p["fiso"]),
                              agree=bool(np.sign(obs_ba) == r["sign"]), fiso_stored=p["fiso"]))
    fsg, ft1, ft4 = O["fiso_sign"], O["fiso_tau_L1"], O["fiso_tau_L4"]
    agree = sum(r["agree"] for r in sign_rows)
    S = float(sum(r["sign"] * np.sign(r["obs_ba"]) for r in sign_rows))
    if agree != fsg["agree"] or len(sign_rows) != fsg["n_pairs"] or S != fsg["S"]:
        stop(f"sign test: recomputed agree {agree}/{len(sign_rows)} S {S} vs stored {fsg}")
    x1 = [r["L1"] for r in sign_rows]
    y = [r["abs_obs"] for r in sign_rows]
    same(kendalltau(x1, y, variant="b").statistic, ft1["tau"], "tau_b L1 from the fixed table and the stored pair values")
    same(kendalltau([r["L4"] for r in sign_rows], y, variant="b").statistic, ft4["tau"], "tau_b L4")
    if ft1["n_pairs"] != len(sign_rows):
        stop("tau: pair count")
    led("N-016", fsg["agree"], "FISO sign: pairs whose observed sign matched")
    led("N-017", fsg["n_pairs"], "FISO sign: pairs with a predicted sign")
    led("N-224", fsg["S"], "FISO sign: S")
    led("N-225", fsg["p"], "FISO sign: one-sided permutation p")
    led("N-220", fsg["n_perm"], "scanner permutations")
    led("N-226", ft1["tau"], "FISO magnitude order: Kendall tau_b (L1), recomputed from the drawn points (not printed)")
    led("N-227", ft1["p"], "FISO tau_b one-sided p (approximate)")
    led("N-228", ft4["tau"], "FISO tau_b with L4")
    led("N-146", O["odi_a"]["n_pairs"], "judgeable scanner pairs")
    for k in ("fiso_sign", "fiso_tau_L1"):                    # v11 load_extra: the stored readings of the two FISO checks
        if "지지되지 않음" not in O[k]["reading"]:
            stop(f"{k}: stored reading «{O[k]['reading']}» is not «not supported»")
    n_judge = sum(1 for r in table if r["n_common"] >= 3)
    no_sign = sorted(r["key"] for r in table if r["n_common"] >= 3 and r["sign"] == 0)
    if n_judge != O["odi_a"]["n_pairs"] or len(no_sign) != n_judge - len(sign_rows):
        stop("judgeable pairs without a predicted sign: count")
    d["no_sign"] = [(k, pj[k]["fiso"]) for k in no_sign]   # not drawn; listed in the values file (caption: no number)
    sign_rows.sort(key=lambda r: (-r["L1"], r["label"]))     # largest predicted shift first
    d.update(sign_rows=sign_rows, table=table, pj=pj)

    # ODI against retest
    ob = O["odi_b"]
    thr = led("N-610", d["odi_thr_read"], "ODI against retest: threshold, fixed document line 97")
    per = ob["per_pair"]
    judge_odi = {k for k, p in pj.items() if p["judgeable_odi"]}
    if set(per) != judge_odi:
        stop("ODI against retest: per-pair keys are not the judgeable pairs")
    odi_rows = []
    for k, e in per.items():
        led("N-179", e["n_units"], "ROI units per pair")
        if e["need"] != math.ceil(thr * e["n_units"] - 1e-9) or (e["n_units"], e["need"]) != (d["odi_units_read"], d["odi_need_read"]):
            stop(f"ODI against retest {k}: need {e['need']} of {e['n_units']}")
        if not LC.same_value(e["need"], LC.split_parts(ROWS["N-237"]["원 값"])[1]):
            stop(f"ODI against retest {k}: need {e['need']} differs from N-237")
        USED["N-237"] = (e["need"], "ODI against retest: units needed of 27 (second part of the raw value; 0.90 is N-610)")
        share = e["n_retest"] / e["n_units"]
        if e["retest_level"] != (e["n_retest"] >= e["need"]) or e["retest_level"] != (share >= thr):
            stop(f"ODI against retest {k}: retest flag disagrees with the share")
        odi_rows.append(dict(pair=k, share=share, n_retest=e["n_retest"], n_units=e["n_units"], need=e["need"], at=e["retest_level"]))
    at = sorted(r["pair"] for r in odi_rows if r["at"])
    if len(at) != ob["n_retest_level"] or at != sorted(ob["retest_pairs"]) or len(odi_rows) != ob["n_pairs"]:
        stop("ODI against retest: count or names of pairs at retest level")
    if at != ["OXF1PRI~OXF2PRI", "OXF2PRI~OXF3TRI"]:
        stop(f"ODI against retest: pairs at retest level are {at}")
    led("N-018", ob["n_retest_level"], "ODI against retest: pairs at retest level")
    led("N-019", ob["n_pairs"], "ODI against retest: pairs described (printed in (a))")
    led("N-360", ob["n_pairs"] - ob["n_retest_level"], "ODI against retest: pairs above retest level")
    odi_rows.sort(key=lambda r: (-r["share"], r["pair"]))
    d.update(odi_rows=odi_rows, odi_thr=thr)
    return d


# ----------------------------------------------------------------------------- acquisitions: TE only (v12 load_acq, TE part)
def device(acq):
    return acq.split("-")[0]


def load_te(d):
    """TE of the eight acquisitions (seven scanners) with the v12 guards, and the matrix cells."""
    A = {}
    # (1) the judge's acquisition table, decomp/s4/plan.py ONH (design · TE range ms · TR range s · coil given)
    prel = "decomp/s4/plan.py"
    lines = frozen(prel, "N-250").decode("utf-8").splitlines()
    k0 = [k for k, s in enumerate(lines) if s.startswith("ONH = {")]
    if len(k0) != 1:
        stop("plan.py: the ONH table is not found once")
    k1 = next((k for k in range(k0[0], len(lines)) if lines[k].strip() == "}"), None)
    if k1 is None:
        stop("plan.py: the ONH table does not close")
    onh = ast.literal_eval("\n".join(re.sub(r"\s#.*$", "", s) for s in lines[k0[0]:k1 + 1]).split("=", 1)[1].strip())
    if len(onh) != 8 or not all(isinstance(v, tuple) and len(v) == 4 and len(v[1]) == 2 for v in onh.values()):
        stop(f"plan.py ONH: unexpected structure {onh}")
    A["plan_lines"] = (k0[0] + 1, k1 + 1)
    led("N-250", [min(v[1][0] for v in onh.values()), max(v[1][1] for v in onh.values())],
        "plan.py ONH: TE range of the eight acquisitions (the judge's acquisition table)")
    row_ok("N-122", "eight acquisitions (seven scanners + the two NOT2ING protocols)")
    if len(onh) != int(float(ROWS["N-122"]["원 값"])):
        stop("plan.py ONH: not eight acquisitions")

    # (2) stage-3 cards: the 28 ON-Harmony acquisition pairs, dTE (ms), B − A for «A~B»
    crel = "decomp/results/stage3/cards/cards.json"
    C = json.loads(frozen(crel, "N-116"))
    cards = C[:28]
    if any(c["set"] != "ON-Harmony" for c in cards) or sum(c["set"] == "ON-Harmony" for c in C) != 28:
        stop("cards.json: the ON-Harmony cards are not the first 28")
    led("N-117", len(cards), "stage-3 cards: ON-Harmony acquisition pairs")
    led("N-116", sum(bool(c["design_differs"]) for c in cards), "stage-3 cards: pairs that differ in direction scheme (file check)")
    cp = {}
    for i, c in enumerate(cards):
        a, b = c["pair"].split("~")
        if a not in onh or b not in onh or frozenset((a, b)) in cp:
            stop(f"cards.json[{i}]: pair {c['pair']}")
        cp[frozenset((a, b))] = dict(i=i, a=a, b=b, dTE=c["dTE"])
    if len(cp) != 28:
        stop("cards.json: not the 28 pairs of eight acquisitions")

    # (3) absolute TE: OXF1PRI's single-valued plan.py entry plus the card differences; every card must reproduce
    ref = "OXF1PRI"
    te_lo, te_hi = onh[ref][1]
    if te_lo != te_hi:
        stop("plan.py: OXF1PRI TE is a range")
    TE, src_te = {ref: te_lo}, {ref: f"plan.py ONH[{ref}] TE {te_lo} ms (anchor)"}
    for acq in onh:
        if acq == ref:
            continue
        c = cp[frozenset((ref, acq))]
        sg = 1.0 if c["a"] == ref else -1.0
        TE[acq] = TE[ref] + sg * c["dTE"]
        src_te[acq] = f"OXF1PRI + cards.json[{c['i']}] ({c['a']}~{c['b']}) dTE"
    for c in cp.values():
        same(TE[c["b"]] - TE[c["a"]], c["dTE"], f"card {c['a']}~{c['b']}: dTE reproduces from the drawn TEs", 1e-9)
    for acq, (_des, ter, _trr, _coil) in onh.items():
        same(TE[acq], (ter[0] + ter[1]) / 2.0, f"{acq}: drawn TE = middle of the plan.py TE range", 1e-9)
    # sidecar records: NOT2ING first-session EchoTime (N-124 · N-126) and the TE range of the NODDI-capable scanners (N-142)
    for acq, i in (("NOT2ING-89", "N-124"), ("NOT2ING-98", "N-126")):
        te_side = float(row_ok(i, f"{acq}: first-session sidecar EchoTime (s), reference for the drawn card TE")) * 1000.0
        if abs(te_side - TE[acq]) >= 0.5:
            stop(f"{acq}: sidecar TE {te_side} ms does not round to the card TE {TE[acq]}")
    rng = [float(v) for v in LC.split_parts(row_ok("N-142", "sidecar TE range of the NODDI-capable scanners (reference)"))]
    if not all(rng[0] - 1e-9 <= v <= rng[1] + 1e-9 for v in TE.values()):
        stop(f"a drawn TE lies outside the sidecar range {rng}")
    devices = sorted({device(a) for a in onh})
    led("N-361", len(devices), "scanners (rows and columns of (a))")

    # (4) the registered scanner table of LOCK_B §198-1 (the copy of N-142): its TE column must contain the drawn values
    lrel = "prereg/LOCK_B_inputs.md"
    lb = frozen(lrel, "N-142").decode("utf-8").splitlines()
    h = [k for k, s in enumerate(lb) if s.startswith("| # | 순위 단위 | 기종 | TE (ms) | TR (s) |")]
    if len(h) != 1:
        stop("LOCK_B: the §198-1 scanner table is not found once")
    lbr = {}
    for k in range(h[0] + 2, len(lb)):
        if not lb[k].startswith("| "):
            break
        c = [x.strip() for x in lb[k].split("|")[1:-1]]
        if not c[0].isdigit():
            continue                                          # the excluded NOT3GEM row («—»)
        lbr[c[1].strip("*")] = dict(line=k + 1, model=c[2].replace("**", ""), te=c[3])
    if sorted(lbr) != devices:
        stop(f"LOCK_B §198-1: scanners {sorted(lbr)}")
    for dev in devices:
        cell = lbr[dev]["te"]
        nv = [float(x) for x in re.findall(r"\d+(?:\.\d+)?", cell)]
        got = sorted(TE[a] for a in onh if device(a) == dev)
        if "–" in cell:
            ok = len(nv) == 2 and all(nv[0] - 1e-6 <= v <= nv[1] + 1e-6 for v in got)
        elif "/" in cell:
            ok = len(nv) == len(got) and all(abs(a_ - b_) <= 1e-6 for a_, b_ in zip(nv, got))
        else:
            ok = len(nv) == 1 and all(abs(nv[0] - v) <= 1e-6 for v in got)
        if not ok:
            stop(f"{dev}: drawn TE {got} vs LOCK_B {cell}")

    dev_acqs = {dev: sorted((a for a in onh if device(a) == dev), key=lambda a: TE[a]) for dev in devices}
    order = sorted(devices, key=lambda dv: (TE[dev_acqs[dv][0]], dv))
    A.update(onh=onh, TE=TE, src_te=src_te, cards=cp, lockb=lbr, dev_acqs=dev_acqs, order=order)
    d["acq"] = A
    for acq in sorted(onh, key=lambda a: (TE[a], a)):
        if acq == ref:
            drawn("a, TE track", f"{acq}: TE (ms), anchor", prel, f"md:{A['plan_lines'][0] + 1} ONH[{ref}][1]", TE[acq], ids=[])
            continue
        c = cp[frozenset((ref, acq))]
        sg = "+" if c["a"] == ref else "−"
        drawn("a, TE track", f"{acq}: TE (ms) = OXF1PRI {TE[ref]} {sg} dTE of card {c['a']}~{c['b']}", crel,
              f"json:[{c['i']}].dTE", TE[acq], ids=[])

    # (5) matrix cells: the 21 device pairs of the fixed table
    cells = {}
    for r in d["table"]:
        a, b = r["pair"].split("~")
        if a not in devices or b not in devices:
            stop(f"fixed table pair {r['pair']}")
        cells[frozenset((a, b))] = dict(pair=r["pair"], key=r["key"], n_common=r["n_common"], sign=r["sign"],
                                        judgeable=r["n_common"] >= 3, line=r["line"])
    if len(cells) != 21:
        stop("matrix: not 21 device pairs")
    for r in d["odi_rows"]:
        cells[frozenset(r["pair"].split("~"))]["odi"] = r
    for r in d["sign_rows"]:
        cells[frozenset(r["pair"].split("~"))]["fiso"] = r
    hatched = sorted((k for k, c in cells.items() if not c["judgeable"]), key=lambda k: sorted(k))
    if [sorted(k) for k in hatched] != [["NOT4GEP", "OXF3TRI"], ["OXF3TRI", "OXF4GEP"]]:
        stop(f"pairs without common participants: {[sorted(k) for k in hatched]}")
    led("N-147", sum(d["pj"][cells[k]["key"]]["n_subjects_fiso"] for k in hatched), "hatched pairs: common participants (sum)")
    for k, c in cells.items():
        if c["judgeable"] != ("odi" in c) or ("fiso" in c) != (c["sign"] != 0 and c["judgeable"]):
            stop(f"matrix cell {c['pair']}: ODI share or FISO mark does not follow the fixed table")
    if sum("fiso" in c for c in cells.values()) != int(float(ROWS["N-017"]["원 값"])):
        stop("FISO points are not the pairs with a predicted sign")
    jrel = "decomp/results/remeasure/judgements.json"
    pj_idx = {p: k for k, p in enumerate(d["pj"])}
    for k, c in sorted(cells.items(), key=lambda kv: kv[1]["pair"]):
        if "odi" in c:
            o = c["odi"]
            drawn("a, matrix", f"{c['key']}: ODI share at retest level = n_retest / n_units = {o['n_retest']}/{o['n_units']}",
                  jrel, f"json:ON-Harmony 쌍 검정.odi_b.per_pair.{c['key']}.n_retest", o["share"], ids=[])
        if "fiso" in c:
            f = c["fiso"]
            drawn("b, scatter y", f"{c['key']}: observed Δ logit FISO (stored; |Δ| = {f['abs_obs']!r} drawn; sign "
                  f"{'agreed' if f['agree'] else 'opposed'} with the predicted {'+' if f['sign'] > 0 else '−'} for {c['pair']})", jrel,
                  f"json:ON-Harmony 쌍 검정.pairs[{pj_idx[c['key']]}].fiso", f["fiso_stored"], ids=[])
            drawn("b, scatter x", f"{c['pair']}: predicted |Δ logit FISO|, largest on the primary relaxation grid (fixed table "
                  f"line {c['line']})", "decomp/remeasure/MEASUREMENT_FIXED.md", f"md:{c['line']}", f["L1"], ids=[])
    d["cells"], d["hatched"] = cells, hatched
    return d


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


def rowaxis(ax, n):
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def mk(ax, x, y, filled, color, marker, ms, z=3, mew=1.0, edge=None):
    ax.plot([x], [y], marker, ms=ms, ls="none", color=color, mfc=color if filled else SURF, mec=edge or color, mew=mew, zorder=z,
            clip_on=False)


FIG_W, FIG_H = 190.0, 96.0
SQ = math.sqrt(math.pi) / 2.0             # square side = SQ × circle diameter gives equal areas
PT_MARK = 5.6                             # (b): circle diameter (pt); square side SQ × it (equal areas)
OFF_BELOW_MM, OFF_MM = 1.9, 1.1           # (b): marks closer than this (mm) hide each other; such a pair is drawn ± OFF_MM (S46 G1)


def odi_cell(ax, x0, y0, w, h, share, outlined, z=1.0):
    """One ODI cell: sequential blue fill (0–1 linear); at or above the fixed threshold a white inner line and a black
    border (two cues besides the colour)."""
    ax.add_patch(Rectangle((x0, y0), w, h, fc=FC.SEQ(Normalize(0.0, 1.0)(share)), ec=SURF, lw=0.8, zorder=z))
    if outlined:
        ax.add_patch(Rectangle((x0 + 0.13 * w, y0 + 0.13 * h), 0.74 * w, 0.74 * h, fc="none", ec=SURF, lw=1.0, zorder=z + 0.2))
        ax.add_patch(Rectangle((x0, y0), w, h, fc="none", ec=INK, lw=1.2, zorder=z + 0.3))


def hatch_cell(ax, x0, y0, w, h, z=1.0):
    ax.add_patch(Rectangle((x0, y0), w, h, fc=SURF, ec=MUTED, lw=0.4, hatch="////", zorder=z))


def te_dots(ax, xs, y, z=3):
    """TE dots of one scanner; two protocols joined by a segment (NOT2ING)."""
    if len(xs) > 1:
        ax.plot([min(xs), max(xs)], [y, y], color=INK2, lw=1.4, zorder=z - 1, solid_capstyle="butt")
    ax.plot(xs, [y] * len(xs), "o", ms=3.4, ls="none", mfc=INK2, mec=SURF, mew=0.5, zorder=z)


# ----------------------------------------------------------------------------- figure
def draw(d):
    W_MM, H_MM = FIG_W, FIG_H
    fig = new_fig(W_MM, H_MM)
    cv = canvas(fig)
    lay = {}
    A = d["acq"]
    order = A["order"]
    n = len(order)
    idx = {dev: k for k, dev in enumerate(order)}
    CELL = 7.6                                                  # matrix cell (mm) = track row pitch
    HB = n * CELL
    TY = H_MM - 4.2
    YT = H_MM - 20.0                                            # top of the row block
    Y0 = YT - HB

    def yrow(k):
        return YT - (k + 0.5) * CELL

    # ======== (a) ODI against retest
    head(fig, cv, 1.5, TY, f"ODI, {T('N-019')} scanner pairs; scanners ordered by TE", [("DESCRIPTIVE", False)],
         ids=("N-019",), letter="a")
    NAME_R = 17.4
    for k, dev in enumerate(order):
        say(cv, NAME_R, yrow(k), dev, (), ha="right", va="center", color=INK)
    # ---- TE track (the only acquisition track here)
    TX, TW = 20.0, 24.0
    say(cv, TX + TW / 2.0, YT + 1.4, "TE (ms)", (), ha="center", va="bottom", color=INK2)
    ax = ax_mm(fig, TX, Y0, TW, HB)
    rowaxis(ax, n)
    lo, hi = 62.0, 102.0
    ax.set_xlim(lo, hi)
    for k, dev in enumerate(order):
        ax.plot([lo, hi], [k, k], color=GRID, lw=LW_THIN, zorder=1, solid_capstyle="butt")
        vals = [A["TE"][a] for a in A["dev_acqs"][dev]]
        if not all(lo < v < hi for v in vals):
            stop(f"TE track: {dev} {vals} outside {lo}–{hi}")
        te_dots(ax, vals, k)
    ax.set_xticks((70, 80, 90, 100))
    ax.set_xticklabels(("70", "80", "90", "100"))
    ax.tick_params(axis="x", length=2, pad=1.5)
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_color(AXIS)
    # ---- the matrix: upper triangle only (row scanner i, column scanner j, i < j)
    MX0 = 48.0
    axm = ax_mm(fig, MX0, Y0, HB, HB)
    axm.set_xlim(-0.5, n - 0.5)
    axm.set_ylim(n - 0.5, -0.5)
    axm.set_xticks([])
    axm.set_yticks([])
    for s in axm.spines.values():
        s.set_visible(False)
    for j, dev in enumerate(order):                             # column names above the matrix (45°, from the column centre)
        say(cv, MX0 + (j + 0.5) * CELL - 1.2, YT + 1.0, dev, (), rotation=45, rotation_mode="anchor", ha="left", va="bottom",
            color=INK)
    for k in range(n):
        axm.add_patch(Rectangle((k - 0.5, k - 0.5), 1, 1, fc=NEUTRAL, ec=SURF, lw=0.8, zorder=0.8))
    n_cells = 0
    for key, c in d["cells"].items():
        a, b = sorted(key, key=lambda dv: idx[dv])
        i, j = idx[a], idx[b]
        if not c["judgeable"]:
            hatch_cell(axm, j - 0.5, i - 0.5, 1, 1)
            continue
        o = c["odi"]
        odi_cell(axm, j - 0.5, i - 0.5, 1, 1, o["share"], o["at"])
        n_cells += 1
    if n_cells != int(float(ROWS["N-019"]["원 값"])):
        stop(f"(a): {n_cells} ODI cells drawn")
    lay["odi_cells"] = n_cells
    lay["hatched_cells"] = len(d["hatched"])

    # ---- keys of (a), two rows under the TE axis and the matrix
    KY1, KY2 = Y0 - 10.0, Y0 - 17.0
    say(cv, 1.5, KY1, "share of units at retest level", (), ha="left", va="center", color=INK)
    cax = ax_mm(fig, 40.0, KY1 - 1.2, 28.0, 2.4)
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0.0, 1.0), cmap=FC.SEQ), cax=cax, orientation="horizontal")
    cax.set_xticks([0.0, 0.5, 1.0])
    cax.set_xticklabels(["0", "0.5", "1"])
    cax.tick_params(axis="x", length=2, pad=1.5, colors=INK2)
    cb.outline.set_linewidth(0.4)
    cb.outline.set_edgecolor(AXIS)
    odi_cell(cv, 73.0, KY1 - 2.4, 4.8, 4.8, 1.0, True, z=2)
    say(cv, 79.6, KY1, f"at or above {T('N-610')}", ("N-610",), ha="left", va="center", color=INK2)
    hatch_cell(cv, 1.5, KY2 - 2.4, 4.8, 4.8, z=2)
    say(cv, 7.8, KY2, "no common participant", (), ha="left", va="center", color=INK2)
    te_dots(cv, [40.0, 45.0], KY2)
    say(cv, 46.6, KY2, "NOT2ING: both protocols", (), ha="left", va="center", color=INK2)

    # ======== (b) FISO sign and magnitude order
    BX = 116.0
    head(fig, cv, BX, TY, f"FISO, {T('N-017')} pairs with a predicted sign", [("UNBLINDED CHECK", False)], ids=("N-017",),
         letter="b")
    sr = d["sign_rows"]
    SX, SW = 131.0, 55.0
    SH = 52.0
    SY = YT - SH
    # key of (b): one line under the title (shape and colour)
    LX, LY = SX + 1.6, YT + 7.4
    mk(cv, LX, LY, True, FC.AGREE, "o", PT_MARK, mew=0.4, edge=SURF)
    say(cv, LX + 2.2, LY, "observed sign agreed with the prediction", (), ha="left", va="center", color=INK2)
    mk(cv, LX, LY - 3.8, True, FC.OPPOSE, "s", PT_MARK * SQ, mew=0.4, edge=SURF)
    say(cv, LX + 2.2, LY - 3.8, "opposed", (), ha="left", va="center", color=INK2)
    ax3 = ax_mm(fig, SX, SY, SW, SH)
    X3, Y3 = (0.0, 0.45), (-0.015, 0.30)
    mmx, mmy = SW / (X3[1] - X3[0]), SH / (Y3[1] - Y3[0])
    dist = sorted((math.hypot((p["L1"] - q["L1"]) * mmx, (p["abs_obs"] - q["abs_obs"]) * mmy), p["key"], q["key"])
                  for i, p in enumerate(sr) for q in sr[i + 1:])
    # Direction S46 G1: marks closer than a mark's width hide each other. That pair is drawn OFF_MM left and right of its
    # predicted value so that all 8 marks show; the stored values, the sign test and τb are not changed (values.md says which).
    OFF = {}
    if dist[0][0] < OFF_BELOW_MM:
        pa, pb = sorted((r for r in sr if r["key"] in dist[0][1:]), key=lambda r: (r["L1"], r["key"]))
        OFF = {pa["key"]: -OFF_MM / mmx, pb["key"]: OFF_MM / mmx}
    xdraw = lambda r: r["L1"] + OFF.get(r["key"], 0.0)
    gaps = [math.hypot((xdraw(p) - xdraw(q)) * mmx, (p["abs_obs"] - q["abs_obs"]) * mmy) for i, p in enumerate(sr) for q in sr[i + 1:]]
    if min(gaps) < OFF_BELOW_MM:
        stop(f"(b): marks still hide each other after the offset ({min(gaps):.2f} mm)")
    for r in sorted(sr, key=lambda r: r["agree"]):            # agreeing marks last (thin white edge)
        m = "o" if r["agree"] else "s"
        mk(ax3, xdraw(r), r["abs_obs"], True, FC.AGREE if r["agree"] else FC.OPPOSE, m, PT_MARK if m == "o" else PT_MARK * SQ,
           z=4 if r["agree"] else 3, mew=0.4, edge=SURF)
    for r in sr:                                               # square outlines on top: a square under a circle stays visible
        if not r["agree"]:
            ax3.plot([xdraw(r)], [r["abs_obs"]], "s", ms=PT_MARK * SQ, ls="none", mfc="none", mec=FC.OPPOSE, mew=0.7, zorder=5,
                     clip_on=False)
    if not all(X3[0] < xdraw(r) < X3[1] and Y3[0] < r["abs_obs"] < Y3[1] for r in sr):
        stop("(b): a point lies outside the axes")
    ax3.set_xlim(*X3)
    ax3.set_ylim(*Y3)
    ax3.set_xticks([0.0, 0.1, 0.2, 0.3, 0.4])
    ax3.set_xticklabels(["0", "0.1", "0.2", "0.3", "0.4"])
    ax3.set_yticks([0.0, 0.1, 0.2, 0.3])
    ax3.set_yticklabels(["0", "0.1", "0.2", "0.3"])
    hgrid(ax3)
    ax3.set_xlabel("predicted |Δ logit FISO|\n(largest on the primary relaxation grid)", linespacing=1.1)
    ax3.set_ylabel("observed |Δ logit FISO|")
    lay["offset_keys"] = sorted(OFF)
    lay["closest_pair_mm"] = (round(dist[0][0], 2), dist[0][1], dist[0][2])
    lay["n_points"] = len(sr)
    lay["row_pitch_mm"] = CELL
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
# Panels and symbols only (S45 §7; S41 §3): the observations (2 of 8 signs, τb, 17 of 19 above retest) are in the body.
CAPTION = ('![**Fig. 5. Across scanner models and protocols.** (a) ODI: share of ROI units at retest level in {N-019} ON-Harmony '
           'scanner pairs (descriptive); outlined, ≥{N-610}; hatched, no common participant; scanners ordered by TE. (b) FISO: '
           'predicted against observed |Δ logit FISO| in the {N-017} pairs with a predicted sign (unblinded check); circle, sign '
           'agreed; square, opposed; overlapping points offset horizontally.](figures/F5_scanners.pdf){#fig:fig5 .wide width=100%}')
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
    L = [f"# Fig. 5 (v13) — reads of run {RUN_ID} (written by make_F5.py)", "",
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
    A = d["acq"]
    w("# F5 (v13) — drawn values and visible numbers (written by make_F5.py)")
    w("")
    w("Every value below was read from a frozen ledger source copy (`drafts/x15/ledger/sources/<file>/<sha12>`, sha12 = first 12 "
      "hex of the sha256 recorded in the ledger; the script stops if the copy's sha256 differs; every read logged first, see "
      "reads.md) and drawn as is, except: the absolute TE of the acquisitions (OXF1PRI's single-valued entry of the judge's "
      "acquisition table plus the card differences, guarded against all 28 cards, the plan.py ranges, the sidecar TE rows and "
      "the registered scanner table); the ODI share (n_retest / n_units); |Δ logit FISO| (absolute value of the stored "
      "difference). «no ledger row (drawn only)» means the value is drawn but not printed and no ledger row cites that field.")
    w("")
    w("## Sources read (frozen ledger copies; sha256 checked)")
    w("")
    w("| Source (ledger 근거 파일) | Frozen copy | sha256 |")
    w("|---|---|---|")
    for (k, sha), p in READ.items():
        w(f"| `{k}` | `{p.relative_to(SOURCES)}` | `{sha}` |")
    w("")
    w("Also read (logged): `decomp/results/remeasure/judgements.json` in the repository, for its sha256 only (equal to the frozen "
      "copy). Read without a log entry: `drafts/x15/ledger/numbers.csv`; `drafts/x15/ledger/ledger_check.py` (imported "
      "read-only); `figcolors.py`.")
    w("")
    w("## (a) TE track")
    w("")
    w("| Row (TE order) | Acquisition | TE ms (drawn) | plan.py TE range | Registered scanner table TE (LOCK_B §198-1) | Source of TE |")
    w("|---|---|---|---|---|---|")
    for k, dev in enumerate(A["order"]):
        for a in A["dev_acqs"][dev]:
            lbm = A["lockb"][dev]
            w(f"| {k + 1} | {a} | {A['TE'][a]!r} | {A['onh'][a][1]} | {esc(lbm['te'])} (line {lbm['line']}) | {esc(A['src_te'][a])} |")
    w("")
    w(f"Order rule: TE of the scanner's shorter-TE protocol, ties by name ({' · '.join(A['order'])}); rows and columns of the "
      "matrix follow it. The TE of OXF4GEP varies within sessions (68.3–71.0 ms) and that of NOT1ACH (70.0–71.2 ms) — LOCK_B "
      "§9 · §198-1; the dot is the card value, the middle of the range. NOT2ING: the two protocols' card values joined by a "
      "segment. Axis 62–102 ms, ticks 70–100 (scale labels, not results).")
    w("")
    w("## (a) matrix — upper triangle (row scanner earlier in TE order than column scanner); 21 device pairs")
    w("")
    w("| Pair (fixed table, A~B) | Stored key | Common participants | Cell (row, column) | ODI share drawn (n_retest/n_units) · outlined |")
    w("|---|---|---|---|---|")
    idx = {dv: k for k, dv in enumerate(A["order"])}
    for key, c in sorted(d["cells"].items(), key=lambda kv: kv[1]["line"]):
        a, b = sorted(key, key=lambda dv: idx[dv])
        pos = f"{a}, {b}"
        if not c["judgeable"]:
            w(f"| {c['pair']} | {c['key']} | {c['n_common']} | {pos} | hatched |")
            continue
        o = c["odi"]
        w(f"| {c['pair']} | {c['key']} | {c['n_common']} | {pos} | {o['share']!r} ({o['n_retest']}/{o['n_units']}) · "
          f"{'yes' if o['at'] else 'no'} |")
    w("")
    w("Diagonal cells (same scanner): neutral fill, no value. Lower triangle: empty (the v12 FISO signs are in (b)).")
    w("")
    w("## (b) FISO — 8 pairs with a predicted sign")
    w("")
    w("Orientation: the fixed table names each pair A~B and predicts the sign of B − A; the stored value of key D_i~D_j is "
      "D_j − D_i, so «observed B − A» is the stored value, negated where the stored key is the reverse of the table name.")
    w("")
    w("| Pair (fixed table) | Stored key | x: predicted \\|Δ\\| (L1, fixed table line) | Predicted sign of B − A | Observed B − A "
      "(stored value) | y: \\|Δ\\| drawn | Agreed · mark |")
    w("|---|---|---|---|---|---|---|")
    for r in d["sign_rows"]:
        w(f"| {r['pair']} | {r['key']} | {r['L1']} (line {r['line']}) | {'+' if r['sign'] > 0 else '−'} | {r['obs_ba']!r} "
          f"(stored {r['fiso_stored']!r}{', reversed key' if r['flip'] else ''}) | {r['abs_obs']!r} | "
          f"{'yes · circle (AGREE)' if r['agree'] else 'no · square (OPPOSE)'} |")
    w("")
    c_ = d["layout"]["closest_pair_mm"]
    w(f"Kendall τb {T('N-226')} (N-226) recomputed from these 8 points equals the stored value (not printed in the figure; in "
      f"the body). The two closest marks are {c_[1]} and {c_[2]}, {c_[0]} mm apart (centres) at their values; they are drawn "
      f"{OFF_MM} mm left and right of their predicted values so that all 8 marks show (direction S46 G1; drawn offset for "
      f"{', '.join(d['layout']['offset_keys'])}; the stored values, the sign test and τb do not change). Judgeable pairs without a predicted sign (not drawn): " +
      "; ".join(f"{k} {v!r}" for k, v in d["no_sign"]) + ".")
    w("")
    w("## Every drawn value with its source")
    w("")
    w("| Panel | Element | File | sha12 | Field | Value | Ledger ID |")
    w("|---|---|---|---|---|---|---|")
    for e in DRAWN:
        w(f"| {e['panel']} | {esc(e['element'])} | `{e['file']}` | `{e['sha']}` | `{esc(e['field'])}` | {esc(repr(e['value']))} | "
          f"{', '.join(e['ids']) if e['ids'] else 'no ledger row (drawn only)'} |")
    w("")
    w("## Visible numbers (figure text other than tick labels)")
    w("")
    w("| Printed text | Ledger ID → 원고 표기 · raw value |")
    w("|---|---|")
    for s, ids in shown:
        w(f"| {esc(s)} | " + " · ".join(f"{i} → «{esc(ROWS[i]['원고 표기'])}» · {esc(ROWS[i]['원 값'])}" for i in ids) + " |")
    w("")
    w("Tick labels (the TE axis, the scatter axes and the colour bar) and the scanner names (row and column labels) are exempt.")
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
    w(f"| FISO sign: {T('N-016')} of {T('N-017')} signs agreed | unblinded check | (b): two circles, six squares |")
    w(f"| FISO magnitude order: τb {T('N-226')} | unblinded check | (b): the scatter itself; the value in the body (4.2, Table 1) |")
    w(f"| ODI above retest level in {T('N-360')} of {T('N-019')} pairs | descriptive | (a): 17 cells without the outline |")
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
      f"{lay['nearest_texts'][1]!r}). TE-track rows and matrix rows share one pitch ({lay['row_pitch_mm']} mm).")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


# ----------------------------------------------------------------------------- main
PYCACHES = (REPO / "decomp/manuscript/figures/src/__pycache__", LEDGER / "__pycache__", HERE / "__pycache__")
PYCACHE_BEFORE = {p for p in PYCACHES if p.exists()}        # this folder is shared with the other figure scripts


def main():
    d = load_te(load())
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
    cvd = FC.cvd_sheet(OUT / f"{NAME}.png", OUT / "cvd" / "F5_cvd.png", title="Fig. 5 (v13): original and three colour-vision simulations")
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
