"""Draws the order of analysis steps and fixing commits (Supplementary Fig. S1 in earlier versions; the list is in docs/fixed_products_and_commits.md), reading each value from its record and stopping on disagreement with the manuscript."""
import datetime as dt
import hashlib
import json
import os
import pathlib
import re
import subprocess
import sys

sys.dont_write_bytecode = True
HERE = pathlib.Path(__file__).resolve().parent
REPO = pathlib.Path(os.environ.get("DMRI_REPO") or "<repo>")
sys.path.insert(0, str(REPO / "decomp" / "manuscript" / "figures" / "src"))
import figstyle as fs  # noqa: E402
import matplotlib.text  # noqa: E402
from matplotlib.font_manager import FontProperties  # noqa: E402
if os.environ.get("FIG_OUT"):          # default: figstyle's own folder (decomp/manuscript/figures/)
    fs.OUTDIR = str(pathlib.Path(os.environ["FIG_OUT"]).resolve())

WORK = pathlib.Path("data")
MS_PATH = REPO / "decomp" / "manuscript" / "manuscript.md"
SUPP_PATH = REPO / "decomp" / "manuscript" / "supplement.md"
HANDOFF = "prereg/artifacts/handoff_a.json"
NAME = "FigS1_commit_chain"
KST = dt.timezone(dt.timedelta(hours=9))
UTC = dt.timezone.utc
FORBIDDEN = (REPO / "prereg" / "artifacts" / "verdict", WORK / "remeasure" / "onharmony",
             WORK / "stage1_work" / "s4" / "v4")
REPORT_S4, REPORT_S4_LINE = "decomp/results/s4/REPORT.md", 138
ONH_LOG = "prereg/artifacts/decomp/onharmony_image_access_log.json"
EXT_LOG = "prereg/artifacts/decomp/external_data_access_log.json"
TONG_JUDGE = WORK / "remeasure" / "tong" / "judge" / "tong_judge.json"
WONG_JUDGE = WORK / "remeasure" / "wong" / "judge" / "wong_judge.json"

# What each commit must be: words its subject carries and files it adds (A) or changes (M).
COMMITS = {
    "5029f4b": (("S4 사양",), (("A", "decomp/SPEC_S4.md"),)),
    "85eb8b0": (("가중치 동결",), (("A", "decomp/results/s4/train_freeze.json"),)),
    "b625516": (("S4 결과를 열기 전", "V4"), (("M", "decomp/SPEC_S4.md"),)),
    "b886d76": (("V4 부록 A", "S4 결과를 열기 전"), (("M", "decomp/SPEC_S4.md"),)),
    "c8825e2": (("V4 코일 조합",), (("M", "decomp/remeasure/judge_v4.py"),)),
    "08a189c": (("⑪ 사전 점검",), (("A", "decomp/remeasure/judge_onh.py"),)),
    "5f95b18": (("고정 문서",), (("A", "decomp/remeasure/MEASUREMENT_FIXED.md"),)),
    "867b56c": (("⑪ 끝", "ON-Harmony 119 / 120", "Tong 36 / 36", "Wong 3 / 3"),
                (("M", ONH_LOG),)),
    "5632626": (("판정 결과 커밋",), (("A", "decomp/results/remeasure/judgements.json"), ("M", ONH_LOG),
                                      ("M", EXT_LOG))),
    "06b7361": (("C7 재실행",), (("M", ONH_LOG),)),
}
H_END = "867b56c"          # closed processing (⑪)
H_JUDG, H_C7 = "5632626", "06b7361"
NO_TEXT_TIME = {H_END, H_C7}   # the texts give no "<hash>, <date> <time> KST" for these two


def stop(msg):
    sys.exit(f"STOP {NAME}: {msg}")


def read_text(p):
    p = pathlib.Path(p).resolve()
    if any(p == f or f in p.parents for f in FORBIDDEN):
        stop(f"refusing to open {p}")
    return p.read_text(encoding="utf-8")


def read_bytes(p):
    p = pathlib.Path(p).resolve()
    if any(p == f or f in p.parents for f in FORBIDDEN):
        stop(f"refusing to open {p}")
    return p.read_bytes()


def git(*args):
    env = dict(os.environ, GIT_OPTIONAL_LOCKS="0", GIT_PAGER="cat")
    r = subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True, env=env)
    if r.returncode:
        stop(f"git {' '.join(args)}: {r.stderr.strip()}")
    return r.stdout


TAG = re.compile(r" ?\[N-\d{3,4}\]")
TEXTS = {p.name: [TAG.sub("", x) for x in read_text(p).split("\n")] for p in (MS_PATH, SUPP_PATH)}
MS = "\n".join("\n".join(v) for v in TEXTS.values())
ALL_LINES = [(f"{n}:{i}", x) for n, v in TEXTS.items() for i, x in enumerate(v, 1)]


def ms_lines_with(pattern):
    rx = re.compile(pattern)
    return [(loc, m, line) for loc, line in ALL_LINES for m in rx.finditer(line)]


def check_ms_times(hits, t, what):
    """Each hit: groups (date or None, HH:MM[:SS]). A missing date must appear in the same line."""
    for loc, m, line in hits:
        d, hm = m.group(1), m.group(2)
        want = t.strftime("%H:%M:%S") if len(hm) == 8 else t.strftime("%H:%M")
        if hm != want:
            stop(f"{what}: {loc} says {hm} KST, record says {want}")
        if d is not None and d != t.strftime("%Y-%m-%d"):
            stop(f"{what}: {loc} says {d}, record says {t:%Y-%m-%d}")
        if d is None and t.strftime("%Y-%m-%d") not in line:
            stop(f"{what}: {loc} gives {hm} without the date {t:%Y-%m-%d}")
    return [loc for loc, _, _ in hits]


# ----------------------------------------------------------------------------- records
def commit(h, subject_words, touches):
    out = git("show", "-s", "--format=%h %ci %s", h).rstrip("\n")
    m = re.fullmatch(r"([0-9a-f]{7,40}) (\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} [+-]\d{4}) (.*)", out)
    if not m or not m.group(1).startswith(h):
        stop(f"git show {h}: unexpected output {out!r}")
    t = dt.datetime.strptime(m.group(2), "%Y-%m-%d %H:%M:%S %z").astimezone(KST)
    subj = m.group(3)
    for w in subject_words:
        if w not in subj:
            stop(f"{h}: subject {subj!r} lacks {w!r}")
    changed = {tuple(x.split("\t", 1)) for x in git("show", "--name-status", "--format=", h).split("\n") if x}
    for st, path in touches:
        if (st, path) not in changed:
            stop(f"{h}: does not {'add' if st == 'A' else 'change'} {path}")
    hits = ms_lines_with(re.escape(h) + r"[,;]?\s*\(?(?:(\d{4}-\d{2}-\d{2})\s+)?(\d{2}:\d{2}(?::\d{2})?)\s+KST")
    lines = check_ms_times(hits, t, h)
    if h not in NO_TEXT_TIME and not lines:
        stop(f"{h}: no '<hash>, <date> <time> KST' found in the manuscript or supplement to check against")
    return dict(h=m.group(1)[:7], t=t, subject=subj, ms=lines, raw=out,
                files=", ".join(f"{st} {p}" for st, p in touches))


def s4_opened():
    if f"`{REPORT_S4}:{REPORT_S4_LINE}`" not in MS:
        stop(f"the manuscript no longer cites {REPORT_S4}:{REPORT_S4_LINE}")
    lines = git("show", f"HEAD:{REPORT_S4}").split("\n")
    rx = re.compile(r"결과를 처음 연 시각 (\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2})Z\((\d{2}:\d{2}:\d{2}) KST\)")
    found = [(i, rx.search(x)) for i, x in enumerate(lines, 1) if rx.search(x)]
    if len(found) != 1 or found[0][0] != REPORT_S4_LINE:
        stop(f"{REPORT_S4}: opening time found at lines {[i for i, _ in found]}, cited :{REPORT_S4_LINE}")
    m = found[0][1]
    t = dt.datetime.fromisoformat(m.group(1)).replace(tzinfo=UTC).astimezone(KST)
    if t.strftime("%H:%M:%S") != m.group(2):
        stop(f"{REPORT_S4}:{REPORT_S4_LINE}: UTC {m.group(1)}Z and KST {m.group(2)} disagree")
    bl = git("blame", "-L", f"{REPORT_S4_LINE},{REPORT_S4_LINE}", "--porcelain", "HEAD", "--", REPORT_S4)
    bh = bl.split()[0]
    bt = dt.datetime.fromtimestamp(int(re.search(r"^committer-time (\d+)$", bl, re.M).group(1)), KST)
    if not bt > t:
        stop(f"{REPORT_S4}:{REPORT_S4_LINE} committed ({bh[:7]} {bt}) before the time it records ({t})")
    hits = ms_lines_with(r"opened\s+(?:at\s+|\()?(?:(\d{4}-\d{2}-\d{2})\s+)?(\d{2}:\d{2}(?::\d{2})?)\s+KST")
    lines_ms = check_ms_times(hits, t, "S4 results opened")
    if not lines_ms:
        stop("no 'opened ... KST' in the manuscript to check against")
    return dict(t=t, line=REPORT_S4_LINE, by=bh[:7], by_t=bt, ms=lines_ms)


def utc_line(text_lines, anchor, utc):
    """1-based line of the '"utc": "<utc>"' that opens the entry holding the anchor line."""
    ia = next((i for i, x in enumerate(text_lines) if anchor in x), None)
    if ia is None:
        stop(f"anchor {anchor!r} not found")
    for i in range(ia, -1, -1):
        if f'"utc": "{utc}"' in text_lines[i]:
            return i + 1
    stop(f"utc {utc} not found above {anchor!r}")


def processing_start():
    text = git("show", f"{H_END}:{ONH_LOG}")
    starts = [e for e in json.loads(text)["열람"] if e.get("세션") == "⑪ 시작"]
    if len(starts) != 1 or not starts[0].get("무엇", "").startswith("⑪ 시작"):
        stop(f"{ONH_LOG}@{H_END}: {len(starts)} entries «⑪ 시작»")
    utc = starts[0]["utc"]
    tu = dt.datetime.fromisoformat(utc).astimezone(UTC)
    z = tu.strftime("%Y-%m-%dT%H:%M:%SZ")
    hits = [loc for loc, x in ALL_LINES if f"entry of {z}" in x]
    if not hits:
        stop(f"the texts do not state the processing-start entry as 'entry of {z}'")
    for loc, m, _ in ms_lines_with(r"entry of (\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z)"):
        if m.group(1) != z:
            stop(f"{loc} gives entry of {m.group(1)}, record {z}")
    line = utc_line(text.split("\n"), '"세션": "⑪ 시작"', utc)
    return dict(t=tu.astimezone(KST), utc=utc, line=line, ms=hits)


def walk(o):
    if isinstance(o, dict):
        if "utc" in o and "파일" in o:
            yield o
        for v in o.values():
            yield from walk(v)
    elif isinstance(o, list):
        for v in o:
            yield from walk(v)


def judgement(result):
    rel = "../" + str(result.relative_to(WORK.parent))
    text = read_text(REPO / EXT_LOG)
    hits = [e for e in walk(json.loads(text)) if any(f.get("file") == rel for f in e["파일"])]
    if len(hits) != 1:
        stop(f"{EXT_LOG}: {len(hits)} entries for {rel}")
    e = hits[0]
    rec = next(f["sha256"] for f in e["파일"] if f["file"] == rel)
    sha = hashlib.sha256(read_bytes(result)).hexdigest()
    if sha != rec:
        stop(f"{result.name}: sha256 {sha[:16]} != access record {rec[:16]}")
    head = json.loads(git("show", f"HEAD:{EXT_LOG}"))
    committed = any(h["utc"] == e["utc"] and any(f.get("file") == rel for f in h["파일"]) for h in walk(head))
    tu = dt.datetime.fromisoformat(e["utc"]).astimezone(UTC)
    line = utc_line(text.split("\n"), f'"file": "{rel}"', e["utc"])
    return dict(t=tu.astimezone(KST), tu=tu, sha=sha, committed=committed, line=line, ms=[])


def onh_entry(fname, commit_h):
    """Recorded time of the HEAD access-log entry whose file ends with fname; commit_h must have added it."""
    def hits(rev):
        return [e for e in json.loads(git("show", f"{rev}:{ONH_LOG}"))["열람"]
                if any(str(f.get("file", "")).endswith("/" + fname) for f in (e.get("파일") or []))]
    head = hits("HEAD")
    if len(head) != 1:
        stop(f"{ONH_LOG}@HEAD: {len(head)} entries for {fname}")
    e = head[0]
    if not any(x["utc"] == e["utc"] for x in hits(commit_h)):
        stop(f"{fname}: entry utc {e['utc']} is not in {ONH_LOG}@{commit_h}")
    if any(x["utc"] == e["utc"] for x in hits(f"{commit_h}^")):
        stop(f"{fname}: entry utc {e['utc']} was already in {commit_h}^ (not added by {commit_h})")
    tu = dt.datetime.fromisoformat(e["utc"]).astimezone(UTC)
    return dict(t=tu.astimezone(KST), utc=e["utc"], what=e.get("무엇", ""), by=commit_h, ms=[])


def c7_record():
    rec = json.loads(git("show", f"{H_JUDG}:{HANDOFF}"))
    k = [x for x in rec if x.startswith("판정(")]
    if len(k) != 1:
        stop(f"{HANDOFF}@{H_JUDG}: {len(k)} «판정» records")
    v4, rerun = rec[k[0]]["V4"], rec[k[0]]["C7 재실행(다음 · 보조 · 커밋 따로)"]
    m = re.search(r"run c7 (\d{2})-(\d{2}) (\d{2}:\d{2}:\d{2}) 시작 · (\d{2})-(\d{2}) (\d{2}:\d{2}:\d{2}) KST 종료 1", v4)
    f = re.search(r"읽기 고정\((\d{4}-\d{2}-\d{2}) (\d{2}:\d{2}:\d{2}) KST", rerun)
    if not m or not f or "v4_c7.json 없음" not in v4:
        stop(f"{HANDOFF}@{H_JUDG}: C7 stop or reading-fixed time not found as expected")
    year = C[H_JUDG]["t"].year
    stopped = dt.datetime.strptime(f"{year}-{m.group(4)}-{m.group(5)} {m.group(6)}",
                                   "%Y-%m-%d %H:%M:%S").replace(tzinfo=KST)
    fixed = dt.datetime.strptime(f"{f.group(1)} {f.group(2)}", "%Y-%m-%d %H:%M:%S").replace(tzinfo=KST)
    s_ms = [loc for loc, x in ALL_LINES if f"({stopped:%Y-%m-%d %H:%M:%S} KST, exit 1)" in x]
    f_ms = [loc for loc, x in ALL_LINES if f"read as fixed on {fixed:%Y-%m-%d %H:%M:%S} KST" in x]
    if not s_ms or not f_ms:
        stop(f"the supplement does not state the C7 stop ({stopped}) and the reading-fixed time ({fixed})")
    return dict(stopped=dict(t=stopped, ms=s_ms), fixed=dict(t=fixed, ms=f_ms))


# ----------------------------------------------------------------------------- read and check
C = {h: commit(h, *v) for h, v in COMMITS.items()}
fixed_doc = git("show", f"{C['5f95b18']['h']}:decomp/remeasure/MEASUREMENT_FIXED.md")
if f"커밋 {C['08a189c']['h']}" not in fixed_doc:
    stop("the fixed document of 5f95b18 does not name the pipeline commit 08a189c")
if "119 of the 120 ON-Harmony sessions were processed, together with all 36 Tong sessions and the data " \
   "of all three Wongkornchaovalit participants" not in MS:
    stop("manuscript's processing counts differ from the subject of 867b56c (119/120 · 36/36 · 3/3)")
OPENED = s4_opened()
PSTART = processing_start()
WONG = judgement(WONG_JUDGE)
TONG = judgement(TONG_JUDGE)
want = f"on {TONG['t']:%Y-%m-%d} at {TONG['t']:%H:%M:%S} KST ({TONG['tu']:%H:%M:%S} UTC)"
TONG["ms"] = [loc for loc, x in ALL_LINES if want in x]
if not TONG["ms"] or TONG["sha"][:16] not in MS:
    stop(f"the texts do not state the Tong judgment as '{want}' with sha256 {TONG['sha'][:16]}…")
for loc, m, _ in ms_lines_with(r"external_data_access_log\.json:(\d+)–(\d+)"):
    if not int(m.group(1)) <= TONG["line"] <= int(m.group(2)):
        stop(f"{loc} cites {EXT_LOG}:{m.group(1)}–{m.group(2)}; Tong entry is at :{TONG['line']}")
V4C3 = onh_entry("v4_c3.json", H_JUDG)
PAIRS = onh_entry("onh_pairs.json", H_JUDG)
C7RERUN = onh_entry("v4_c7.json", H_C7)
C7 = c7_record()
if not (C7["fixed"]["t"] < C[H_JUDG]["t"] < C7RERUN["t"]):
    stop("the C7 reading was not fixed before the judgment commit, or the rerun result predates it")
if not (WONG["committed"] and TONG["committed"]):
    stop("a judgment entry of the external-data access record is not in HEAD")

def history(path):
    out = []
    for x in git("log", "--format=%h %ci", "--", path).split("\n"):
        if x:
            h, t = x.split(" ", 1)
            out.append((h, dt.datetime.strptime(t, "%Y-%m-%d %H:%M:%S %z").astimezone(KST)))
    return out


# "fixed before": no commit changed the specification or the freeze record after the S4 results were
# opened, and none changed the fixed document after the release commit
for path, limit, what in (("decomp/SPEC_S4.md", OPENED["t"], "S4 results opened"),
                          ("decomp/results/s4/train_freeze.json", OPENED["t"], "S4 results opened"),
                          ("decomp/remeasure/MEASUREMENT_FIXED.md", C["5f95b18"]["t"], "release 5f95b18")):
    late = [f"{h} {t:%Y-%m-%d %H:%M:%S}" for h, t in history(path) if t > limit]
    if late:
        stop(f"{path} changed after {what}: {late}")
SPEC_HISTORY = history("decomp/SPEC_S4.md")

v4a, v4b = C["b625516"], C["b886d76"]
if v4a["t"].date() != v4b["t"].date():
    stop("the two V4-rule commits fall on different dates; the label shows one date")
REC = "access record"
ROW_A = [
    dict(name="Judge specification\nfirst commit", src=C["5029f4b"]["h"], t=[C["5029f4b"]["t"]]),
    dict(name="Weights\nfrozen", src=C["85eb8b0"]["h"], t=[C["85eb8b0"]["t"]]),
    dict(name="V4 rules\nfixed", src=f"{v4a['h']} · {v4b['h']}", t=[v4a["t"], v4b["t"]]),
    dict(name="Validation results\nfirst opened", src="judge report", rec=True, t=[OPENED["t"]]),
    dict(name="V4 open points\nfixed in code", src=C["c8825e2"]["h"], t=[C["c8825e2"]["t"]]),
    dict(name="Pipeline\nfixed", src=C["08a189c"]["h"], t=[C["08a189c"]["t"]]),
    dict(name="Real-data rules\nreleased", src=C["5f95b18"]["h"], t=[C["5f95b18"]["t"]]),
]
ROW_B = [
    dict(name="Processing\nstarted", src=REC, rec=True, t=[PSTART["t"]]),
    dict(name="Processing\nended", src=C[H_END]["h"], t=[C[H_END]["t"]]),
    dict(name="Wongkornchaovalit\njudgment", src=REC, rec=True, t=[WONG["t"]]),
    dict(name="Tong\njudgment", src=REC, rec=True, t=[TONG["t"]]),
    dict(name="V4, C3 arm\nresult", src=REC, rec=True, t=[V4C3["t"]]),
]
ROW_C = [
    dict(name="V4, C7 arm stopped\n(no result)", src="handoff record", rec=True, open=True,
         t=[C7["stopped"]["t"]]),
    dict(name="ON-Harmony\npair tests", src=REC, rec=True, t=[PAIRS["t"]]),
    dict(name="Verdicts\ncommitted", src=C[H_JUDG]["h"], t=[C[H_JUDG]["t"]]),
    dict(name="C7 arm rerun\n(auxiliary)", src=REC, rec=True, t=[C7RERUN["t"]]),
]
done = [t for s in ROW_A + ROW_B + ROW_C for t in s["t"]]
if any(b <= a for a, b in zip(done, done[1:])):
    stop("the steps are not in time order: " + ", ".join(f"{t:%m-%d %H:%M:%S}" for t in done))

# ----------------------------------------------------------------------------- draw
W, H = 190.0, 86.0
X0, ROW_W, ROW_H = 8.0, 178.0, 20.0     # row axes in mm; panel letter 7 mm left of X0
Y_MARK, GAP, LINE = 11.0, 2.0, 2.95     # mark height in the row; label gap; 7-pt line pitch (mm)
XM = 13.0                                # first and last mark inset from the row ends
N_SLOT = max(len(ROW_A), len(ROW_B), len(ROW_C))


def lp_metrics_mm(fig):
    """Ascent and descent of "lp" at 7 pt in mm, so every label line sits on a fixed baseline
    (va="top"/"bottom" would move a line that has, or lacks, a descender)."""
    px = 25.4 / fig.dpi
    _, h, d = fig.canvas.get_renderer().get_text_width_height_descent(
        "lp", FontProperties(family=fs.FONT, size=fs.PT_TEXT), ismath=False)
    return (h - d) * px, d * px


def draw_row(fig, y0, letter, header, steps):
    ax = fs.ax_mm(fig, X0, y0, ROW_W, ROW_H)
    ax.set_xlim(0, ROW_W)
    ax.set_ylim(0, ROW_H)
    ax.axis("off")
    asc, desc = lp_metrics_mm(fig)
    xs = [XM + i * (ROW_W - 2 * XM) / (N_SLOT - 1) for i in range(len(steps))]
    for i in range(len(steps) - 1):
        ax.plot(xs[i:i + 2], [Y_MARK] * 2, color=fs.AXIS, lw=0.75, zorder=1, solid_capstyle="butt")
    for x, s in zip(xs, steps):
        ax.plot([x], [Y_MARK], ls="None", marker="o", ms=5.2, mew=0.8, mec=fs.INK,
                mfc=fs.SURF if s.get("open") else fs.INK, zorder=3)
        names = s["name"].split("\n")
        for k, txt in enumerate(names):          # last name line just above the mark
            ax.text(x, Y_MARK + GAP + desc + (len(names) - 1 - k) * LINE, txt, ha="center",
                    va="baseline", fontsize=fs.PT_TEXT, color=fs.INK)
        rows = [
            (s["src"], fs.INK2 if s.get("rec") else fs.INK, "italic" if s.get("rec") else "normal"),
            (" · ".join(sorted({f"{t:%Y-%m-%d}" for t in s["t"]})), fs.INK2, "normal"),
            (" · ".join(f"{t:%H:%M:%S}" for t in s["t"]), fs.INK2, "normal")]
        for k, (txt, col, sty) in enumerate(rows):
            ax.text(x, Y_MARK - GAP - asc - k * LINE, txt, ha="center", va="baseline",
                    fontsize=fs.PT_TEXT, color=col, style=sty)
    fs.panel(ax, letter, header)
    return ax


fig = fs.new_fig(W, H)
draw_row(fig, 58.0, "a", "Judge and real-data rules", ROW_A)
draw_row(fig, 30.0, "b", "Real-data processing and first judgments", ROW_B)
draw_row(fig, 2.0, "c", "C7 arm, pair tests, verdict commit", ROW_C)
key = fs.canvas(fig)
ky = H - 3.6
kx = 128.0
key.plot([kx], [ky], ls="None", marker="o", ms=5.2, mew=0.8, mec=fs.INK, mfc=fs.INK)
key.text(kx + 2.2, ky, "done", ha="left", va="center", fontsize=fs.PT_TEXT, color=fs.INK)
key.plot([kx + 12.5], [ky], ls="None", marker="o", ms=5.2, mew=0.8, mec=fs.INK, mfc=fs.SURF)
key.text(kx + 14.7, ky, "stopped", ha="left", va="center", fontsize=fs.PT_TEXT, color=fs.INK)
fig.canvas.draw()
_r = fig.canvas.get_renderer()
_boxes = [(x.get_text(), x.get_window_extent(_r)) for x in fig.findobj(matplotlib.text.Text) if x.get_text()]
for _i in range(len(_boxes)):
    for _j in range(_i + 1, len(_boxes)):
        _a, _b = _boxes[_i][1], _boxes[_j][1]
        if _a.x0 < _b.x1 - 0.5 and _b.x0 < _a.x1 - 0.5 and _a.y0 < _b.y1 - 0.5 and _b.y0 < _a.y1 - 0.5:
            stop(f"labels overlap: {_boxes[_i][0]!r} and {_boxes[_j][0]!r}")
print(f"labels {len(_boxes)} · no overlap")
key.text(W - 4.0, ky, "times KST (UTC+9)", ha="right", va="center", fontsize=fs.PT_TEXT, color=fs.INK2)

# ----------------------------------------------------------------------------- source table
R = "git -C dmri_cov show -s --format='%h %ci %s'"
table = [
    ("S4 specification first commit", C["5029f4b"], f"{R} 5029f4b · adds decomp/SPEC_S4.md"),
    ("Weights frozen", C["85eb8b0"], f"{R} 85eb8b0 · adds decomp/results/s4/train_freeze.json"),
    ("V4 rules fixed (1)", C["b625516"], f"{R} b625516 · changes decomp/SPEC_S4.md"),
    ("V4 rules fixed (2)", C["b886d76"], f"{R} b886d76 · changes decomp/SPEC_S4.md"),
    ("V4 open points in code", C["c8825e2"], f"{R} c8825e2 · changes decomp/remeasure/judge_v4.py"),
    ("S4 results first opened", dict(h="—", t=OPENED["t"], ms=OPENED["ms"]),
     f"git show HEAD:{REPORT_S4} line {OPENED['line']} (written by {OPENED['by']} "
     f"{OPENED['by_t']:%Y-%m-%d %H:%M:%S})"),
    ("Pipeline fixed", C["08a189c"], f"{R} 08a189c · named in MEASUREMENT_FIXED.md@5f95b18"),
    ("Real-data rules released", C["5f95b18"], f"{R} 5f95b18 · adds decomp/remeasure/MEASUREMENT_FIXED.md"),
    ("Processing started", dict(h="—", t=PSTART["t"], ms=PSTART["ms"]),
     f"git show {H_END}:{ONH_LOG} line {PSTART['line']} (entry «⑪ 시작», utc {PSTART['utc']})"),
    ("Processing ended", C[H_END], f"{R} {H_END} · changes {ONH_LOG}"),
    ("Wongkornchaovalit judgement", dict(h="—", t=WONG["t"], ms=WONG["ms"]),
     f"{EXT_LOG}:{WONG['line']} ({'HEAD' if WONG['committed'] else 'working copy, not in HEAD'}); "
     f"sha256 {WONG['sha'][:16]}… = {WONG_JUDGE.name}"),
    ("Tong judgement", dict(h="—", t=TONG["t"], ms=TONG["ms"]),
     f"{EXT_LOG}:{TONG['line']} ({'HEAD' if TONG['committed'] else 'working copy, not in HEAD'}); "
     f"sha256 {TONG['sha'][:16]}… = {TONG_JUDGE.name}"),
    ("V4, C3 arm result", dict(h="—", t=V4C3["t"], ms=[]),
     f"{ONH_LOG}@HEAD entry utc {V4C3['utc']} (added by {H_JUDG})"),
    ("V4, C7 arm stopped", dict(h="—", t=C7["stopped"]["t"], ms=C7["stopped"]["ms"]),
     f"{HANDOFF}@{H_JUDG} «판정».V4"),
    ("ON-Harmony pair tests", dict(h="—", t=PAIRS["t"], ms=[]),
     f"{ONH_LOG}@HEAD entry utc {PAIRS['utc']} (added by {H_JUDG})"),
    ("Judgments committed", C[H_JUDG], f"{R} {H_JUDG} · adds decomp/results/remeasure/judgements.json"),
    ("C7 reading fixed (caption)", dict(h="—", t=C7["fixed"]["t"], ms=C7["fixed"]["ms"]),
     f"{HANDOFF}@{H_JUDG} «판정».C7 재실행"),
    ("C7 arm rerun (auxiliary)", dict(h="—", t=C7RERUN["t"], ms=[]),
     f"{ONH_LOG}@HEAD entry utc {C7RERUN['utc']} (added by {H_C7})"),
]
print(f"{'step':28s} {'hash':8s} {'KST':19s}  manuscript lines checked · source")
for name, rec, src in table:
    lines = ",".join(str(i) for i in rec["ms"]) or "—"
    print(f"{name:28s} {rec['h']:8s} {rec['t']:%Y-%m-%d %H:%M:%S}  {lines} · {src}")
for h, c in C.items():
    print(f"  subject {c['raw']}")
print("  decomp/SPEC_S4.md commits, all before the opening:",
      " · ".join(f"{h} {t:%m-%d %H:%M:%S}" for h, t in reversed(SPEC_HISTORY)))

fs.save(fig, NAME)
for ext in ("pdf", "png"):
    f = fs.OUTDIR + os.sep + f"{NAME}.{ext}"
    print(f"sha256 {NAME}.{ext} {hashlib.sha256(pathlib.Path(f).read_bytes()).hexdigest()}")
