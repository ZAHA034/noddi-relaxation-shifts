"""Supplementary Fig. S5 (v13; review direction S45 §7, author decision 2026-11-08): where the noise-floor mask lay, and
the noise floor it measured — v12c Fig. 7c moved to the supplement, unchanged.

Left: projections of the stored noise-floor mask of every first Siemens session in template space (rules
drafts/x15/v12/V12B_S43_RULES.md §1.7; drafts/x15/v12/out/c1_projections.npz, sha256 recorded in c1_mask_location.json),
3 scanners × 3 views (axial from above, coronal from behind, sagittal from the left; one scale); gray, the template
brain-mask projection; colour (figcolors.SEQ from 0.25 to 1), the share of the scanner's sessions with a mask voxel on the
line; dashed, the boundary of the projected lateral-ventricle label. Right: one dot per first session at its session
statistic R (rules V12C_S44_RULES.md §2; read_sessions.jsonl, logged re-read; per scanner the count equals the summary's
and the median equals R_med, N-452, N-454, N-456), the scanner's median with its bootstrap interval and effective L
(N-443, N-446, N-449), dashed lines at the Rician floor and the assumed L (N-429, N-430). The functions are those of v12c
make_F7.py (sha256 aac943e9d657…), copied verbatim; only the title line (no panel letter), the labels («S5 …») and main()
are new. No image, mask, transform or template is read; no new computation.
Writes (FIG_OUT only): FigS5_noise_floor_mask.pdf/.png, FS5_caption.md, FS5_values.md, FS5_reads.md, FS5_reads.json,
cvd/FS5_cvd.png.
Run (CPU only):
  PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=<scratch> CUDA_VISIBLE_DEVICES= FIG_OUT=<out> DMRI_REPO=<repo> nice -n 10 \
    <repo>/.venv/bin/python -I -B make_FS5.py
"""
import sys

sys.dont_write_bytecode = True            # no __pycache__ beside the scripts, in the ledger folder or in the repository

import datetime  # noqa: E402
import hashlib  # noqa: E402
import importlib.util  # noqa: E402
import io  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import os  # noqa: E402
import pathlib  # noqa: E402
import subprocess  # noqa: E402

os.environ["CUDA_VISIBLE_DEVICES"] = ""   # CPU only

_HERE = pathlib.Path(__file__).resolve().parent
# v13: the judge-matrix script of v12 (v12 make_F6.py, copied as judge_matrix_v12.py) holds the shared reading, guards,
# caption check and save of the after-the-verdict figures; loaded by path (nothing is drawn or written on import).
_spec = importlib.util.spec_from_file_location("judge_matrix_shared", _HERE / "judge_matrix_v12.py")
B = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(B)

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402
import matplotlib.colors  # noqa: E402
import matplotlib.patches  # noqa: E402
import matplotlib.transforms as mtransforms  # noqa: E402
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator  # noqa: E402

fs, FC = B.fs, B.FC
REPO, WORK = B.REPO, B.WORK
DRAFTS = WORK / "stage1_work" / "drafts"
X14 = DRAFTS / "x14"
LEDGER_DIR = B.LEDGER_DIR
REMEASURE = WORK / "remeasure"
B2K = REMEASURE / "explore" / "b2000"
RULES_REC = DRAFTS / "x15" / "v12" / "RULES_RECORD.jsonl"
REREAD_LOG = X14 / "REREAD_LOG.jsonl"
OUT = pathlib.Path(B.TARGET)
SCRIPT = pathlib.Path(__file__).resolve()
SIE, OTHERS = B.SIEMENS, B.OTHERS
INK, INK2, MUTED, GRID, AXIS, SURF = fs.INK, fs.INK2, fs.MUTED, fs.GRID, fs.AXIS, fs.SURF
PT = fs.PT_TEXT
check, show = B.check, B.show

RUN_ID = "FS5-" + datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
WHAT = "v13 보충 그림 S5 — 소음 바닥 마스크 투영 · 회기 R(판정 뒤 표시) · 새 계산 없음 · 판정 불변 #34"
NAME = "FigS5_noise_floor_mask"
LEDGER_SRC_IDS = ("N-413", "N-116", "N-382", "N-416", "N-428", "N-029", "N-028", "N-423", "N-452", "N-493", "N-459")
L_IDS = {"OXF1PRI": "N-443", "OXF2PRI": "N-446", "OXF3TRI": "N-449"}
RULES_B = DRAFTS / "x15" / "v12" / "V12B_S43_RULES.md"
RULES_B_SHA = "d56e5563ef178cb7b786f614f08a9613d3529e697780a24b79cfda69186484b2"
C1_DIR = pathlib.Path(os.environ.get("C1_OUT") or DRAFTS / "x15" / "v12" / "out")
C1_JSON, C1_NPZ = C1_DIR / "c1_mask_location.json", C1_DIR / "c1_projections.npz"
VIEWS = ("axial", "coronal", "sagittal")
VIEW_GEO = {"axial": (2, 0, +1, 1), "coronal": (1, 0, +1, 2), "sagittal": (0, 1, -1, 2)}
VIEW_HEAD = {"axial": "axial\n(from above)", "coronal": "coronal\n(from behind)", "sagittal": "sagittal\n(from the left)"}
SEQ_LO = 0.25
GREY_BRAIN = "#e4e3df"
CROP_MARGIN = 2
RULES_C = DRAFTS / "x15" / "v12" / "V12C_S44_RULES.md"
RULES_C_SHA = "4b73860865207b2235d447d209bbdeabe743def383f9ed7e261dedaa624a668a"
R_IDS = {"OXF1PRI": "N-452", "OXF2PRI": "N-454", "OXF3TRI": "N-456"}
SR_DOT_MM = 1.1
SR = {}
W = 190
H1 = 0.0                                  # one row only (v12c: the row below panels (a) and (b))
K = 0.12
PJ_X0, PJ_GAP = 15.0, 1.6
PJ_TOP_OFF = 14.0
PJ_BOTTOM_ROOM = 15.0
H2 = None
H = None
HX, HW = 96.0, 90.0
XLO, XHI = 1.0, 6.0
DOT = FC.SEQ_STEPS[4]
C1 = {}
READS7 = []
CASE_FOR_CAPTION = "ㄴ"                   # the caption phrase «where most mask voxels lay» holds for this (c-1) case only


def stop(msg):
    sys.exit(f"STOP FS5: {msg}")


def now():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def sha_bytes(b):
    return hashlib.sha256(b).hexdigest()


def rel_repo(p):
    return os.path.relpath(pathlib.Path(p).resolve(), REPO)


def log_entry(log, pred):
    """(1-based line, entry) of the last entry of a JSONL log that satisfies pred."""
    lines = log.read_bytes().decode("utf-8").split("\n")
    for i in range(len(lines) - 1, -1, -1):
        s = lines[i].strip()
        if not s:
            continue
        try:
            e = json.loads(s)
        except json.JSONDecodeError:
            continue
        if pred(e):
            return i + 1, e
    stop(f"{log.name}: the entry just written was not found")


def reread(path, what, want_sha=None):
    """A derived file: logged with the ledger's reread_log.py (it hashes and appends under its lock), then read once; the
    bytes used are the bytes logged."""
    path = pathlib.Path(path).resolve()
    w = f"{WHAT} · {what} · 실행 {RUN_ID}"
    r = subprocess.run([sys.executable, "-I", "-B", str(LEDGER_DIR / "reread_log.py"), str(path), w], capture_output=True,
                       text=True)
    if r.returncode != 0:
        stop(f"reread_log.py {path}: {r.stdout.strip()} {r.stderr.strip()}")
    h = r.stdout.strip().splitlines()[-1].strip()
    raw = path.read_bytes()
    if sha_bytes(raw) != h:
        stop(f"{path} changed between the logged hash and the read")
    if want_sha and h != want_sha:
        stop(f"{path}: sha256 {h} is not the recorded {want_sha}")
    ln, e = log_entry(REREAD_LOG, lambda x: x.get("다시읽음") == str(path) and x.get("sha256") == h and x.get("무엇") == w)
    READS7.append(dict(kind="derived file, re-read logged by reread_log.py", file=str(path), sha256=h, bytes=len(raw),
                       log=f"drafts/x14/{REREAD_LOG.name}", line=ln, time=e.get("시각"), what=what, run=RUN_ID))
    return raw


def log_originals():
    """The repository files read_values() opens itself (judgements.json, cards.json) and the remeasure originals behind the
    ledger copies it reads, logged as re-read before the run; each must still equal its ledger copy."""
    seen = {}
    for i in LEDGER_SRC_IDS:
        r = B.LEDGER[i]
        rel, sha = r["근거 파일"], r["근거 파일 sha256"]
        if rel in seen:
            continue
        p = (REPO / rel).resolve()
        if not (str(p).startswith(str(REMEASURE.resolve()) + os.sep) or rel.startswith("decomp/results/")):
            continue
        reread(p, "장부 사본의 원본 — 판정 뒤 그림 값(7a · 7b · 7c 표지) · 장부 sha256 대조", want_sha=sha)
        seen[rel] = sha
    return seen


def session_R(v):
    """V12C_S44_RULES.md §2: per scanner, the session statistic R of every first session that entered the scanner summary
    (read_sessions.jsonl, last line per session as the summary, R not null). The re-read is logged; the counts must equal the
    summary's (ledger copy of b2k_result_siemens.json: n_records, n_value, n_sigma_q_drop 0) and the median must equal R_med
    (N-452, N-454, N-456)."""
    if sha_bytes(RULES_C.read_bytes()) != RULES_C_SHA:
        stop("V12C_S44_RULES.md differs from the recorded sha256")
    recs_r = [json.loads(s) for s in RULES_REC.read_text(encoding="utf-8").splitlines() if s.strip()]
    hit = [r for r in recs_r if r.get("sha256") == RULES_C_SHA and str(r.get("파일", "")).endswith("V12C_S44_RULES.md")]
    if len(hit) != 1:
        stop("RULES_RECORD.jsonl does not hold V12C_S44_RULES.md once")
    raw = reread(B2K / "read_sessions.jsonl", "단계 3 회기 기록 — 첫 회기마다 회기 통계 R(그림 7c 오른쪽 띠의 점 · v12c E3 (d))")
    recs = [json.loads(s) for s in raw.decode("utf-8").splitlines() if s.strip()]
    last = {}
    for r in recs:
        last[r["session"]] = r
    BS = v["srcs7"]["BS"]
    out = {}
    for d in SIE:
        names = sorted(k for k, r in last.items() if r.get("device") == d)
        inc = [k for k in names if last[k].get("R") is not None]
        e = BS.get("devices", d)
        if e.get("n_sigma_q_drop", 0) != 0 or len(names) != e["n_records"] or len(inc) != e["n_value"]:
            stop(f"{d}: sessions {len(names)} / with R {len(inc)} differ from the summary's n_records {e['n_records']} / "
                 f"n_value {e['n_value']} (n_sigma_q_drop {e.get('n_sigma_q_drop')})")
        R = np.asarray([float(last[k]["R"]) for k in inc])
        med = float(np.median(R))
        if abs(med - float(v["dev"][d]["R"])) > 1e-12 * max(1.0, abs(med)):
            stop(f"{d}: median of the session R {med!r} differs from R_med {v['dev'][d]['R']!r}")
        check(R_IDS[d], med, f"7c: guard — median of the {len(inc)} session R = R_med")
        out[d] = dict(R=R, sessions=inc, n=len(inc), median=med)
    SR.update(rules_recorded=hit[0]["시각"], sha=sha_bytes(raw), n={d: out[d]["n"] for d in SIE})
    return recs, out


def load_c1(recs):
    """c1_mask_location.json and c1_projections.npz (written by c1_mask_location.py under the S43 rules): rules file and
    record, npz sha256 = the JSON, session set = read_sessions.jsonl (Siemens sessions with a stored R), counts within the
    session numbers. Our own outputs (drafts/x15/v12/out), not remeasure files: hashed, not logged."""
    if sha_bytes(RULES_B.read_bytes()) != RULES_B_SHA:
        stop("V12B_S43_RULES.md differs from the recorded sha256")
    recs_r = [json.loads(s) for s in RULES_REC.read_text(encoding="utf-8").splitlines() if s.strip()]
    hit = [r for r in recs_r if r.get("sha256") == RULES_B_SHA and str(r.get("파일", "")).endswith("V12B_S43_RULES.md")]
    if len(hit) != 1:
        stop("RULES_RECORD.jsonl does not hold V12B_S43_RULES.md once")
    jraw = C1_JSON.read_bytes()
    J = json.loads(jraw)
    if J["rules"]["sha256"] != RULES_B_SHA or J["rules"]["recorded"] != hit[0]["시각"]:
        stop("c1_mask_location.json was not made under the recorded S43 rules")
    nraw = C1_NPZ.read_bytes()
    if sha_bytes(nraw) != J["projections"]["sha256"]:
        stop("c1_projections.npz differs from the sha256 in c1_mask_location.json")
    last = {}
    for r in recs:
        last[r["session"]] = r
    want = {dev: sorted(k for k, r in last.items() if r.get("device") == dev and r.get("R") is not None) for dev in SIE}
    have = {dev: sorted(r["session"] for r in J["sessions"] if r["scanner"] == dev) for dev in SIE}
    if want != have or any(J["session_set"]["n_sessions"][dev] != len(want[dev]) for dev in SIE):
        stop("the sessions of c1_mask_location.json differ from read_sessions.jsonl")
    z = np.load(io.BytesIO(nraw), allow_pickle=False)
    P = {k: z[k] for k in z.files}
    aff = np.asarray(P["template_affine"], float)
    for dev in SIE:
        for view in VIEWS:
            c = P[f"count__{dev}__{view}"]
            if c.min() < 0 or c.max() > len(want[dev]):
                stop(f"{dev} {view}: counts outside 0 … {len(want[dev])}")
    for r in J["sessions"]:                            # the stored warped voxels are the ones the JSON counted
        if len(P[f"warped__{r['session']}"]) != r["n_template_voxels"]:
            stop(f"{r['session']}: stored warped voxels differ from the JSON")
    C1.update(json_sha=sha_bytes(jraw), npz_sha=sha_bytes(nraw), run=J["run"]["id"], case=J["case"], text=J["text"],
              n={dev: len(want[dev]) for dev in SIE}, rules_recorded=hit[0]["시각"], script_sha=J["run"]["script_sha256"])
    for f, raw_, what in ((C1_JSON, jraw, "(c-1) record — session set, case (the caption phrase), npz sha256"),
                          (C1_NPZ, nraw, "(c-2) projections — counts per scanner and view, brain mask, label")):
        READS7.append(dict(kind="own output of c1_mask_location.py (not a remeasure file; sha256 checked, not logged)",
                           file=os.path.relpath(f.resolve(), DRAFTS.parent.parent.parent), sha256=sha_bytes(raw_),
                           bytes=len(raw_), log="—", line="—", time=now(), what=what, run=RUN_ID))
    return dict(P=P, aff=aff, n=C1["n"], J=J)


def world_to_voxel(aff):
    A = np.asarray(aff, float)[:3, :3]
    dom = [int(np.argmax(np.abs(A[:, j]))) for j in range(3)]
    if sorted(dom) != [0, 1, 2]:
        stop("template voxel axes are not aligned with x, y, z")
    return {w: dom.index(w) for w in range(3)}, A


def orient_proj(P2, view, aff):
    """A projection (remaining voxel axes in increasing order) -> display array: rows top to bottom = vertical world axis
    decreasing (A or S at the top); columns left to right = horizontal world axis increasing (axial, coronal: left
    hemisphere on the left, neurological) or decreasing (sagittal: anterior on the left). Transpose and flips only."""
    wv, A = world_to_voxel(aff)
    p, h, sgn, vv = VIEW_GEO[view]
    vp, vh, vvx = wv[p], wv[h], wv[vv]
    rem = [a for a in range(3) if a != vp]
    D = P2 if rem[0] == vvx else P2.T
    if A[vv, vvx] > 0:
        D = D[::-1, :]
    if (A[h, vh] < 0) if sgn > 0 else (A[h, vh] > 0):
        D = D[:, ::-1]
    return D


def outline_loops(M):
    """Closed boundary loops of a 2-D boolean array along pixel edges (pixel (r, c) = [c, c + 1] × [r, r + 1]), collinear
    vertices merged; returns lists of (x, y) vertices in pixel units (y downward)."""
    Pd = np.pad(M, 1)
    nxt = {}
    for r, c in zip(*np.nonzero(M)):
        if not Pd[r, c + 1]:
            nxt.setdefault((c + 1, r), []).append((c, r))
        if not Pd[r + 1, c]:
            nxt.setdefault((c, r), []).append((c, r + 1))
        if not Pd[r + 2, c + 1]:
            nxt.setdefault((c, r + 1), []).append((c + 1, r + 1))
        if not Pd[r + 1, c + 2]:
            nxt.setdefault((c + 1, r + 1), []).append((c + 1, r))
    n_edges = sum(len(x) for x in nxt.values())
    loops, used = [], 0
    while nxt:
        start = next(iter(nxt))
        loop, cur = [start], start
        while True:
            outs = nxt[cur]
            n = outs.pop()
            used += 1
            if not outs:
                del nxt[cur]
            if n == start:
                break
            loop.append(n)
            cur = n
        simp = []
        for i, q in enumerate(loop):
            a, b = loop[i - 1], loop[(i + 1) % len(loop)]
            if (q[0] - a[0]) * (b[1] - q[1]) - (q[1] - a[1]) * (b[0] - q[0]) != 0:
                simp.append(q)
        loops.append(simp + simp[:1])
    if used != n_edges:
        stop("outline tracing lost edges")
    return loops


def c_geometry(c1):
    """Oriented projections, one crop box per view (the same for the three scanners), tile sizes; sets H2 and H."""
    global H2, H
    P, aff = c1["P"], c1["aff"]
    views = {}
    for view in VIEWS:
        brain = orient_proj(P[f"brain__{view}"], view, aff)
        lv = orient_proj(P[f"lv__{view}"], view, aff)
        cnt = {dev: orient_proj(P[f"count__{dev}__{view}"], view, aff) for dev in SIE}
        u = brain | lv
        for dev in SIE:
            u = u | (cnt[dev] > 0)
        rr, cc = np.flatnonzero(u.any(1)), np.flatnonzero(u.any(0))
        r0, r1 = max(0, rr.min() - CROP_MARGIN), min(u.shape[0], rr.max() + 1 + CROP_MARGIN)
        c0, c1_ = max(0, cc.min() - CROP_MARGIN), min(u.shape[1], cc.max() + 1 + CROP_MARGIN)
        crop = (slice(r0, r1), slice(c0, c1_))
        views[view] = dict(brain=brain[crop], lv=lv[crop], cnt={d: cnt[d][crop] for d in SIE}, box=(int(r0), int(r1),
                           int(c0), int(c1_)), full=list(u.shape), nrow=int(r1 - r0), ncol=int(c1_ - c0))
    row_px = max(g["nrow"] for g in views.values())
    row_mm = row_px * K
    block = 3 * row_mm + 2 * PJ_GAP
    H2 = round(PJ_TOP_OFF + block + PJ_BOTTOM_ROOM, 1)
    H = H1 + H2
    return dict(views=views, row_mm=row_mm, block=block, top=H2 - PJ_TOP_OFF)


def paint(view_g, dev, n_dev):
    """RGBA display array: gray where the brain-mask projection is, the share colour where a session's mask is on the line."""
    br, cnt = view_g["brain"], view_g["cnt"][dev]
    img = np.zeros(br.shape + (4,), float)
    img[br] = matplotlib.colors.to_rgba(GREY_BRAIN)
    s = cnt.astype(float) / n_dev
    on = cnt > 0
    img[on] = FC.SEQ(SEQ_LO + (1.0 - SEQ_LO) * s[on])
    return img, s


def dashed_outline(ax, loops, z=4):
    for lp in loops:
        xs, ys = [q[0] for q in lp], [q[1] for q in lp]
        ax.plot(xs, ys, color=SURF, lw=1.1, solid_joinstyle="miter", solid_capstyle="butt", zorder=z, clip_on=False)
        ln, = ax.plot(xs, ys, color=INK, lw=0.55, ls=(0, (1.6, 1.1)), dash_capstyle="butt", zorder=z + 0.1, clip_on=False)


def title_line(fig, cv, x, y, title, tag_text):
    """A single-panel figure: title and role tag on one line, no panel letter."""
    t = cv.text(x, y, title, fontsize=PT, color=INK, ha="left", va="center")
    B.tag(fig, cv, x + B.width_mm(fig, t) + 2.6, y, tag_text, dashed=True)


def draw_c(fig, cv, v, sr, c1, g):
    """Projections of the noise-floor mask in template space (left) and the session-R strips (right); v12c Fig. 7c. Labels «7c …» are relabelled «S5 …»."""
    title_line(fig, cv, 2.0, H2 - 3.2, "The noise-floor mask in template space, and the noise floor", "POST HOC")
    top, row_mm = g["top"], g["row_mm"]
    xs, x = {}, PJ_X0
    for view in VIEWS:
        xs[view] = x
        x += g["views"][view]["ncol"] * K + PJ_GAP
    x_end = x - PJ_GAP
    for view in VIEWS:
        cv.text(xs[view] + g["views"][view]["ncol"] * K / 2, top + 0.8, VIEW_HEAD[view], fontsize=PT, color=INK,
                ha="center", va="bottom", linespacing=1.05)
    painted = {}
    for i, dev in enumerate(SIE):
        rtop = top - i * (row_mm + PJ_GAP)
        cv.text(PJ_X0 - 1.5, rtop - row_mm / 2, dev, fontsize=PT, color=INK, ha="right", va="center")
        for view in VIEWS:
            vg = g["views"][view]
            tw, th = vg["ncol"] * K, vg["nrow"] * K
            ax = fs.ax_mm(fig, xs[view], rtop - (row_mm - th) / 2 - th, tw, th)
            img, s = paint(vg, dev, c1["n"][dev])
            ax.imshow(img, interpolation="nearest", origin="upper", extent=(0, vg["ncol"], vg["nrow"], 0))
            ax.set_xlim(0, vg["ncol"])
            ax.set_ylim(vg["nrow"], 0)
            ax.axis("off")
            loops = outline_loops(vg["lv"])
            dashed_outline(ax, loops)
            on = vg["cnt"][dev] > 0
            painted[(dev, view)] = int(on.sum())
            full = c1["P"][f"count__{dev}__{view}"]
            if int((full > 0).sum()) != int(on.sum()) or int(full.sum()) != int(vg["cnt"][dev].sum()):
                stop(f"{dev} {view}: the crop cut painted pixels")
            B.DRAWN.append(dict(panel="7c projection", element=f"{dev} {view}: share of the {c1['n'][dev]} first sessions with "
                                f"a warped mask voxel on the line ({int(on.sum())} painted pixels; share {float(s[on].min())!r} "
                                f"– {float(s[on].max())!r}; colour FC.SEQ({SEQ_LO} + {1 - SEQ_LO} × share))",
                                file="drafts/x15/v12/out/c1_projections.npz", sha=C1["npz_sha"][:12],
                                field=f"npz:count__{dev}__{view} / {c1['n'][dev]}",
                                value=f"counts 1 – {int(vg['cnt'][dev].max())} of {c1['n'][dev]}", ids=[]))
    for view in VIEWS:
        vg = g["views"][view]
        B.DRAWN.append(dict(panel="7c projection", element=f"{view}: gray, projection of MNI152_T1_1mm_brain_mask "
                            f"({int(vg['brain'].sum())} pixels); crop rows {vg['box'][0]}–{vg['box'][1]}, columns "
                            f"{vg['box'][2]}–{vg['box'][3]} of the oriented {vg['full'][0]} × {vg['full'][1]} projection",
                            file="drafts/x15/v12/out/c1_projections.npz", sha=C1["npz_sha"][:12], field=f"npz:brain__{view}",
                            value=f"{int(vg['brain'].sum())} pixels", ids=[]))
        B.DRAWN.append(dict(panel="7c projection", element=f"{view}: dashed, boundary of the projection of the "
                            f"lateral-ventricle label (Harvard-Oxford values 3, 14; {int(vg['lv'].sum())} pixels), along "
                            "pixel edges", file="drafts/x15/v12/out/c1_projections.npz", sha=C1["npz_sha"][:12],
                            field=f"npz:lv__{view}", value=f"{int(vg['lv'].sum())} pixels", ids=[]))
    # ---- key (words only, rule §1.7)
    yb = top - g["block"]
    y1, y2 = yb - 4.2, yb - 8.6
    sw, sh = 4.0, 2.4
    cv.add_patch(matplotlib.patches.Rectangle((PJ_X0, y1 - sh / 2), sw, sh, fc=GREY_BRAIN, ec="none"))
    t1 = cv.text(PJ_X0 + sw + 1.2, y1, "template brain mask", fontsize=PT, color=INK, ha="left", va="center")
    x2 = PJ_X0 + sw + 1.2 + B.width_mm(fig, t1) + 4.0
    cv.plot([x2, x2 + 6.0], [y1, y1], color=SURF, lw=1.1, solid_capstyle="butt")
    cv.plot([x2, x2 + 6.0], [y1, y1], color=INK, lw=0.55, ls=(0, (1.6, 1.1)), dash_capstyle="butt")
    cv.text(x2 + 7.2, y1, "lateral-ventricle label", fontsize=PT, color=INK, ha="left", va="center")
    gw = 12.0
    kax = fs.ax_mm(fig, PJ_X0, y2 - sh / 2, gw, sh)
    kax.imshow(FC.SEQ(np.linspace(SEQ_LO, 1.0, 64))[None, :, :], aspect="auto", interpolation="nearest",
               extent=(0, 1, 0, 1))
    kax.axis("off")
    cv.text(PJ_X0 + gw + 1.2, y2, "fewer → more sessions with a mask voxel on the line", fontsize=PT, color=INK,
            ha="left", va="center")

    # ---- v12c: session statistic R of every first session (one dot per session), the v11 Siemens marks, the reference lines
    rmin = min(float(sr[d]["R"].min()) for d in SIE)
    rmax = max(float(sr[d]["R"].max()) for d in SIE)
    xlo = min(XLO, float(np.floor(rmin)))
    xhi = max(XHI, float(np.ceil(rmax)))
    row_h = row_mm                                             # the strips are the projection rows
    mm_per_r = HW / (xhi - xlo)
    L4 = show("7c label of the L = 4 line", ["N-429"], check("N-429", v["coil_l"][32]))
    L8 = show("7c label of the L = 8 line", ["N-430"], check("N-430", v["coil_l"][64]))
    labels = {1.0: "Rician", 4.0: f"L = {L4}", 8.0: f"L = {L8}"}
    axes, spread = [], {}
    for i, d in enumerate(SIE):
        rtop = top - i * (row_h + PJ_GAP)
        ax = fs.ax_mm(fig, HX, rtop - row_h, HW, row_h)
        axes.append(ax)
        ax.set_xlim(xlo, xhi)
        ax.set_ylim(0, row_h)                                  # y in mm within the strip
        for L, xr in v["ref"].items():
            ax.axvline(xr, color=AXIS, lw=fs.LW_THIN, ls=(0, (3.0, 2.0)), zorder=1)
        r = sr[d]["R"]
        y_mid, half = 0.40 * row_h, 0.40 * row_h - SR_DOT_MM / 2
        off = fs.swarm(r, SR_DOT_MM / mm_per_r, SR_DOT_MM * 0.92, half)
        spread[d] = float(np.abs(off).max())
        ax.plot(r, y_mid + off, ls="none", marker="o", ms=SR_DOT_MM * 72 / 25.4, mfc=DOT, mec="none", zorder=3,
                clip_on=False)
        e = v["dev"][d]
        ym = 0.86 * row_h
        ax.plot([e["ci"][0], e["ci"][1]], [ym, ym], color=INK, lw=1.0, solid_capstyle="butt", zorder=4)
        B.mark(ax, e["R"], ym, color=INK, marker="o", ms=3.6, filled=True)
        t_L = show(f"7c effective L beside {d}", [L_IDS[d]], check(L_IDS[d], e["L"]))
        ax.text(e["ci"][1] + 0.08, ym, f"L {t_L}", fontsize=PT, color=INK, ha="left", va="center", zorder=6)
        ax.spines["left"].set_visible(False)
        ax.set_yticks([])
        ax.xaxis.set_major_locator(FixedLocator(list(range(int(xlo), int(xhi) + 1))))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.xaxis.set_major_formatter(FuncFormatter(lambda t, _: f"{t:g}"))
        if i < len(SIE) - 1:
            ax.tick_params(axis="x", labelbottom=False)
        ax.text(-0.02, 0.40, d, transform=ax.transAxes, fontsize=PT, color=INK, ha="right", va="center")
        B.DRAWN.append(dict(panel="7c strip", element=f"{d}: session statistic R (median of the per-voxel R in the "
                            f"session's mask) of each of the {sr[d]['n']} first sessions in the scanner summary, one dot per "
                            f"session (vertical offsets only to separate dots); range {float(r.min())!r} – {float(r.max())!r}; "
                            f"median {sr[d]['median']!r} = R_med ({R_IDS[d]}, checked)",
                            file="data/remeasure/explore/b2000/read_sessions.jsonl", sha=SR["sha"][:12],
                            field=f"jsonl: R (device {d}, R not null; last line per session)", value=f"{len(r)} values",
                            ids=[]))
    blend = mtransforms.blended_transform_factory(axes[0].transData, axes[0].transAxes)
    for L, xr in v["ref"].items():
        axes[0].text(xr, 1.0 + 0.04, labels[L], transform=blend, fontsize=PT, color=INK2, ha="center", va="bottom")
    vol = show("7c x-axis label", ["N-508"], check("N-508"))
    bval = check("N-507")                          # «b = 2000 s/mm²»: the figure prints the shell name, the caption the unit
    if not bval.endswith(" s/mm²"):
        stop(f"N-507 is {bval!r}")
    bshown = show("7c x-axis label (shell name; the unit s/mm² is in the caption)", ["N-507"], bval[:-len(" s/mm²")])
    axes[-1].set_xlabel(f"R, mean / SD over the {vol} volumes at {bshown}")
    cv.text(HX + HW / 2, yb - 11.6, "Philips, GE: not determinable", fontsize=PT, color=INK2, ha="center", va="center")
    return dict(xlo=xlo, xhi=xhi, row_h=row_h, mm_per_r=mm_per_r, spread=spread, x_end=x_end, painted=painted,
                views={vw: dict(box=g["views"][vw]["box"], nrow=g["views"][vw]["nrow"], ncol=g["views"][vw]["ncol"])
                       for vw in VIEWS})


def relabel(s):
    """v12c Fig. 7c labels -> Supplementary Fig. S5; the Philips / GE values (read and guarded, not drawn) as guards."""
    if s.startswith("7c"):
        rest = s[2:]
        if rest.lstrip(": ").split(":")[0] in OTHERS:
            return "S5 (guard, not drawn)" + rest
        return "S5" + rest
    return s


def relabel_all():
    for e in B.DRAWN:
        e["panel"] = relabel(e["panel"])
    for e in B.SHOWN:
        e["where"] = relabel(e["where"])
    for i, roles in B.USED.items():
        B.USED[i] = [relabel(r) for r in roles]


def reads_md(head):
    L = [head, "", f"Run `{RUN_ID}` ({now()}). No image, mask, transform or template was read by this script.", "",
         "| # | kind | file | sha256 | bytes | log : line | time | what |", "|---|---|---|---|---|---|---|---|"]
    for i, x in enumerate(READS7, 1):
        f = x["file"] if not str(x["file"]).startswith("/") else rel_repo(x["file"])
        L.append(f"| {i} | {x['kind']} | `{f}` | `{x['sha256']}` | {x['bytes']} | {x['log']}:{x['line']} | {x.get('time')} | "
                 f"{x['what']} |")
    L += ["", "Ledger copies read by judge_matrix_v12.read_values() (drafts/x15/ledger/sources/<file>/<sha12>; checked "
          "against the ledger sha256, not logged — they are the ledger's own copies):", ""]
    for rel, sha in sorted(B.READS.items()):
        L.append(f"- `{rel}` `{sha[:12]}`")
    return "\n".join(L) + "\n"


def merged_srcs(*dicts):
    out, seen = {}, set()
    for d in dicts:
        for k, s in d.items():
            if s.rel not in seen:
                seen.add(s.rel)
                out[f"{len(out)}_{k}"] = s
    return out


CAPTION = ("![**Fig. S5. Where the noise-floor mask lay, and the noise floor it measured.** Displayed after the verdict; "
           "no test. Left: mask of all first Siemens sessions, projected in template space (axial from above, coronal from "
           "behind, sagittal from the left; one scale): darker, more sessions; gray, the template brain mask; dashed "
           "outline, the lateral-ventricle label, where most mask voxels lay (Supplementary Table S10 (g)). Right: "
           "session-median R of every first session (dots) and the scanner's effective L (large dot); dashed lines, the "
           "Rician floor and the assumed L; display settings in Supplementary Table S9 (c).]"
           "(figures/FigS5_noise_floor_mask.pdf){#fig:figs5 .wide width=100%}")
NOTES5 = [
    "v13 (direction S45 §7): v12c Fig. 7c moved to the supplement unchanged (same functions, values and guards); the title "
    "line keeps the v12 title without a panel letter; the main text keeps only the three effective L values.",
    "Caption: the v12c Fig. 7c sentences with the projection views named (axial from above, coronal from behind, sagittal "
    "from the left) and the gray silhouette defined; «where most mask voxels lay» holds for (c-1) case ㄴ of "
    "V12B_S43_RULES.md §1.6 (the script stops if c1_mask_location.json gives another case).",
    "Bootstrap 95 % intervals of the large dots: session bootstrap widths 0.030, 0.057, 0.066 (N-453, N-455, N-457) at 18.0 "
    "mm per unit = 0.54 – 1.19 mm, within the 3.6 pt dots with their 0.8 pt edge; drawn behind the dots.",
    "Projection counts and session R values have no ledger rows (the medians of the session R are N-452, N-454, N-456): "
    "drawn only, never printed. One dot per session; the vertical offsets only separate the dots.",
    "The Rician line is labelled «Rician» (position N-462 = N-458, plan line 1638); the L = 4 and L = 8 lines print N-429 "
    "and N-430.",
]


def main():
    log_originals()
    v = B.read_values()                                   # every data guard of judge_matrix_v12.py; nothing written before
    recs, sr = session_R(v)
    c1 = load_c1(recs)
    if C1["case"] != CASE_FOR_CAPTION:
        stop(f"(c-1) case {C1['case']} is not the case {CASE_FOR_CAPTION} the caption phrase was written for")
    if C1["n"] != SR["n"]:
        stop(f"session numbers of the projections {C1['n']} and of the strip {SR['n']} differ")
    g = c_geometry(c1)
    fig = fs.new_fig(W, H)
    cv = fs.canvas(fig)
    geo = draw_c(fig, cv, v, sr, c1, g)
    for e in B.DRAWN:                                     # Philips / GE: read and guarded, not drawn (as v12c)
        if e["panel"] == "7c" and e["element"].split(":")[0] in OTHERS:
            e["panel"] = "S5 (guard, not drawn)"
    relabel_all()
    extra = ["## Projections and session R", "",
             f"Projections: `drafts/x15/v12/out/c1_projections.npz` (sha256 {C1['npz_sha']}) = the sha256 recorded in "
             f"`c1_mask_location.json` (sha256 {C1['json_sha']}; run {C1['run']}); rules `V12B_S43_RULES.md` §1.7 (recorded "
             f"{C1['rules_recorded']}). Sessions per scanner: " + ", ".join(f"{dev} {C1['n'][dev]}" for dev in SIE) +
             f". (c-1) case {C1['case']} (V {C1['text']['V']}, S {C1['text']['S']}, O {C1['text']['O']}; not printed).", "",
             f"Tiles: one scale, {K} mm on paper per mm of the template; crop boxes: " + "; ".join(
                 f"{vw} rows {b_['box'][0]}–{b_['box'][1]}, columns {b_['box'][2]}–{b_['box'][3]} ({b_['ncol']} × "
                 f"{b_['nrow']} pixels)" for vw, b_ in geo["views"].items()) + ". Painted pixels per tile: " + "; ".join(
                 f"{dev} {vw} {n_}" for (dev, vw), n_ in geo["painted"].items()) + ".", "",
             f"R strips (rules V12C_S44_RULES.md §2, recorded {SR['rules_recorded']}): session statistic R from "
             f"`remeasure/explore/b2000/read_sessions.jsonl` (sha256 {SR['sha']}; logged re-read).", "",
             "| scanner | sessions (dots) | R range | median of the dots | = R_med (ledger) |", "|---|---|---|---|---|"]
    for dev in SIE:
        e = sr[dev]
        extra.append(f"| {dev} | {e['n']} | {float(e['R'].min())!r} – {float(e['R'].max())!r} | {e['median']!r} | "
                     f"{R_IDS[dev]} ({B.LEDGER[R_IDS[dev]]['원고 표기']}) |")
    extra += ["", f"Dot strips: R axis {geo['xlo']:g} – {geo['xhi']:g} ({geo['mm_per_r']:.3f} mm per unit), one dot of "
              f"{SR_DOT_MM} mm per session; strip height {geo['row_h']:.2f} mm (= projection row). Figure 190 × {H:.1f} mm.",
              ""]
    srcs = {k: v["srcs7"][k] for k in ("BS", "BA", "KS", "PL", "CF")}
    audit = B.finish(fig, "S5", NAME, CAPTION, srcs, NOTES5, extra)
    head = ("# v13 Supplementary Fig. S5 — files read (by make_FS5.py)\n\nProjections: c1_projections.npz (sha256 "
            f"{C1['npz_sha']}) and c1_mask_location.json (sha256 {C1['json_sha']}; its own reads are listed in it).")
    (OUT / "FS5_reads.md").write_text(reads_md(head), encoding="utf-8")
    (OUT / "FS5_reads.json").write_text(json.dumps(READS7, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    cvd = FC.cvd_sheet(OUT / f"{NAME}.png", OUT / "cvd" / "FS5_cvd.png")
    print("cvd sheet:", cvd.name, hashlib.sha256(cvd.read_bytes()).hexdigest()[:16])
    print("files re-read (logged):", len(READS7), "· ledger copies read:", len(B.READS))
    print(audit)


if __name__ == "__main__":
    main()
