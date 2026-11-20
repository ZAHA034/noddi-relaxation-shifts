"""make_F4 draws Supplementary Fig. S4 (Figure 4 of earlier versions): the share of units at retest level by ROI unit in the negative control, using guards and ledger checks."""
import sys

sys.dont_write_bytecode = True        # first statement: no __pycache__ anywhere

import csv  # noqa: E402
import datetime  # noqa: E402
import fcntl  # noqa: E402
import gzip  # noqa: E402
import hashlib  # noqa: E402
import importlib.util  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import os  # noqa: E402
import pathlib  # noqa: E402
import re  # noqa: E402
import subprocess  # noqa: E402
from decimal import Decimal, ROUND_HALF_UP  # noqa: E402

os.environ["CUDA_VISIBLE_DEVICES"] = ""   # CPU only


def stop(msg):
    sys.exit(f"STOP F4: {msg}")


# ---------------------------------------------------------------------------------------------------------- paths
if not os.environ.get("DMRI_REPO") or not os.environ.get("FIG_OUT"):
    stop("set DMRI_REPO (repository) and FIG_OUT (output folder)")
for _k in ("REMEASURE_WORK", "REMEASURE_FSLDIR"):
    if os.environ.get(_k):
        stop(f"{_k} is set; the processing paths of rconfig.py must be its defaults")
REPO = pathlib.Path(os.environ["DMRI_REPO"]).resolve()
WORK = pathlib.Path(os.environ.get("DMRI_WORK") or REPO / "data").resolve()
DRAFTS = WORK / "stage1_work" / "drafts"
LEDGER_DIR = DRAFTS / "x15" / "ledger"
SOURCES = LEDGER_DIR / "sources"
X14 = DRAFTS / "x14"
V12 = DRAFTS / "x15" / "v12"
HERE = pathlib.Path(__file__).resolve().parent
OUT = pathlib.Path(os.environ["FIG_OUT"]).resolve()
NAME = "F4_negctrl_roi_prisma"
EXT_LOG = X14 / "EXTERNAL_READ_LOG.jsonl"
REREAD = LEDGER_DIR / "reread_log.py"
RULES = V12 / "V12_DISPLAY_RULES.md"
RULES_SHA = "a4197909b5bb54f2056fdc48f1f2087cb02c73723521d5e0639cb3c07c1cc0f8"
RULES_B = V12 / "V12B_S43_RULES.md"                     # v12b rules (§0, §3), recorded 2026-11-07 03:14:31 KST
RULES_B_SHA = "d56e5563ef178cb7b786f614f08a9613d3529e697780a24b79cfda69186484b2"
RULES_C = V12 / "V12C_S44_RULES.md"                     # v12c rules (§3: the 0.5 – 1.0 scale), recorded 2026-11-07 17:13:29 KST
RULES_C_SHA = "4b73860865207b2235d447d209bbdeabe743def383f9ed7e261dedaa624a668a"
SCALE_LO = 0.5                    # V12C §3: SEQ spans 0.5 – 1.0; below 0.5 the lightest colour (extend min)
RUN_ID = "F4-" + datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
WHAT = "v13g 표시 · 그림 4 (a) 뺌(방향 S53 · 저자 승인 2026-11-09) · 표 S18 값 대조 — post hoc display · displayed after the verdict · no test · 판정 불변 #34"
S18_JSON = V12 / "out" / "F4_S18.json"
S18_MD = pathlib.Path(os.environ.get("F4_S18_MD") or OUT.parent / "fig4" / "table_S18.md")
V12_SHARES = os.environ.get("F4_V12_SHARES")           # the v12 run's F4_shares.json (the shares must equal it)
PYCACHE_WATCH = [REPO / "decomp/remeasure/__pycache__", REPO / "decomp/stage1/__pycache__", LEDGER_DIR / "__pycache__",
                 HERE / "__pycache__"]


def pycache_state():
    """Names, sizes and mtimes in the __pycache__ folders the imports could touch (must not change during a run)."""
    st = {}
    for d in PYCACHE_WATCH:
        st[str(d)] = sorted((p.name, p.stat().st_size, p.stat().st_mtime_ns) for p in d.iterdir()) if d.is_dir() else None
    return st


PYCACHE_BEFORE = pycache_state()

# repository code, imported read-only (the fixed test's own functions)
sys.path.insert(0, str(REPO / "decomp" / "remeasure"))
import rconfig as R  # noqa: E402
import guard as G  # noqa: E402
import roi as ROI  # noqa: E402
import judge_tong as JT  # noqa: E402

if R.REPO.resolve() != REPO or R.WORK.resolve() != (REPO / "data" / "remeasure").resolve() \
        or R.FSLDIR.resolve() != (REPO.parent / "fsl").resolve():
    stop(f"rconfig paths {R.REPO} · {R.WORK} · {R.FSLDIR} are not the defaults next to {REPO}")
SKELETON_FILE = R.FSLDIR / (R.SKELETON + ".nii.gz")
FSL_ATLAS = R.FSLDIR / "data" / "atlases" / "JHU" / "JHU-ICBM-labels-1mm.nii.gz"   # compared with R.ATLAS only, not drawn

# Standard short forms of the JHU ICBM-DTI-81 white-matter labels, keyed by the unit name that roi.label_units gives (the label
# name without « L» / « R»). The script stops unless the keys are exactly the 27 units of the logged label file.
ABBR = {
    "Middle cerebellar peduncle": "MCP",
    "Pontine crossing tract (a part of MCP)": "PCT",
    "Genu of corpus callosum": "GCC",
    "Body of corpus callosum": "BCC",
    "Splenium of corpus callosum": "SCC",
    "Fornix (column and body of fornix)": "FX",
    "Corticospinal tract": "CST",
    "Medial lemniscus": "ML",
    "Inferior cerebellar peduncle": "ICP",
    "Superior cerebellar peduncle": "SCP",
    "Cerebral peduncle": "CP",
    "Anterior limb of internal capsule": "ALIC",
    "Posterior limb of internal capsule": "PLIC",
    "Retrolenticular part of internal capsule": "RLIC",
    "Anterior corona radiata": "ACR",
    "Superior corona radiata": "SCR",
    "Posterior corona radiata": "PCR",
    "Posterior thalamic radiation (include optic radiation)": "PTR",
    "Sagittal stratum (include inferior longitidinal fasciculus and inferior fronto-occipital fasciculus)": "SS",
    "External capsule": "EC",
    "Cingulum (cingulate gyrus)": "CGC",
    "Cingulum (hippocampus)": "CGH",
    "Fornix (cres) / Stria terminalis (can not be resolved with current resolution)": "FX/ST",
    "Superior longitudinal fasciculus": "SLF",
    "Superior fronto-occipital fasciculus (could be a part of anterior internal capsule)": "SFO",
    "Uncinate fasciculus": "UNC",
    "Tapetum": "TAP",
}

_fc_spec = importlib.util.spec_from_file_location("figcolors", HERE / "figcolors.py")
FC = importlib.util.module_from_spec(_fc_spec)
_fc_spec.loader.exec_module(FC)

_spec = importlib.util.spec_from_file_location("ledger_check", LEDGER_DIR / "ledger_check.py")
LC = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(LC)          # read-only use: eval_field · same_value · rounding_ok · numbers_in · FORBIDDEN

LEDGER = {r["ID"]: r for r in csv.DictReader((LEDGER_DIR / "numbers.csv").open(encoding="utf-8"))}

LOGGED = []                           # (log file, line number, file read, sha256, what) for reads.md


# ---------------------------------------------------------------------------------------------------------- read logs
def now():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def utcnow():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def append_log(path, entries):
    """Append JSON lines to an existing log under an exclusive lock on the log itself, in one write; returns line numbers."""
    data = "".join(json.dumps(e, ensure_ascii=False, default=str) + "\n" for e in entries).encode("utf-8")
    fd = os.open(path, os.O_RDWR | os.O_APPEND)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX)
        try:
            os.lseek(fd, 0, os.SEEK_SET)
            n_before = 0
            while True:
                blk = os.read(fd, 1 << 20)
                if not blk:
                    break
                n_before += blk.count(b"\n")
            if os.write(fd, data) != len(data):
                stop(f"short write to {path}")
        finally:
            fcntl.flock(fd, fcntl.LOCK_UN)
    finally:
        os.close(fd)
    return list(range(n_before + 1, n_before + 1 + len(entries)))


def read_template(path, what, basis):
    """A template image: bytes read once, sha256 and size logged (EXTERNAL_READ_LOG, b2k_read.py entry form) BEFORE the
    values are used; the image is then parsed from the same bytes."""
    raw = pathlib.Path(path).read_bytes()
    h = hashlib.sha256(raw).hexdigest()
    entry = dict(시각=now(), utc=utcnow(), 세션="템플릿(참가자 자료 아님)", 파일=os.path.relpath(path, REPO), sha256=h,
                 바이트=len(raw), 무엇=f"{WHAT} · {what} · 실행 {RUN_ID}", 근거=basis, 스크립트=str(pathlib.Path(__file__).resolve()))
    lines = append_log(EXT_LOG, [entry])
    LOGGED.append((str(EXT_LOG.relative_to(DRAFTS)), lines[0], os.path.relpath(path, REPO), h, what))
    import nibabel as nib
    img = nib.Nifti1Image.from_bytes(gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw)
    return img, h, len(raw)


def reread(path, what):
    """A derived table or reference file: logged by the ledger's reread_log.py first (it hashes and appends under its lock),
    then read here; the bytes read must have the sha256 the log recorded."""
    r = subprocess.run([sys.executable, "-I", "-B", str(REREAD), str(path), f"{WHAT} · {what} · 실행 {RUN_ID}"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        stop(f"reread_log.py {path}: {r.stdout.strip()} {r.stderr.strip()}")
    logged = r.stdout.strip().splitlines()[-1]
    raw = pathlib.Path(path).read_bytes()
    h = hashlib.sha256(raw).hexdigest()
    if h != logged:
        stop(f"{path}: sha256 read {h[:12]} differs from the logged {logged[:12]}")
    LOGGED.append(("x14/REREAD_LOG.jsonl", None, str(path), h, what))
    return raw, h


def reread_line_numbers():
    """Fill in the REREAD_LOG line numbers of this run's entries (matched by RUN_ID)."""
    lines = (X14 / "REREAD_LOG.jsonl").read_text(encoding="utf-8").splitlines()
    mine = {}
    for i, ln in enumerate(lines, 1):
        if RUN_ID in ln:
            rec = json.loads(ln)
            mine.setdefault((rec["다시읽음"], rec["sha256"]), []).append(i)
    out = []
    for log, line, f, h, what in LOGGED:
        if log == "x14/REREAD_LOG.jsonl":
            k = (str(pathlib.Path(f).resolve()), h)
            if not mine.get(k):
                stop(f"no REREAD_LOG line of this run for {f}")
            line = mine[k].pop(0)
        out.append((log, line, f, h, what))
    LOGGED[:] = out


# ---------------------------------------------------------------------------------------------------------- template, units, slices
def label_table(xml_raw):
    """48 label names and the 27 units, by the repository's own functions (roi.label_units · unit_order); the same rule is
    applied here to the logged bytes and must agree. Side of each label from its name (… L / … R; none = midline)."""
    text = xml_raw.decode("utf-8", errors="replace")
    names = {int(m.group(1)): m.group(2).strip() for m in re.finditer(r'<label index="(\d+)"[^>]*>([^<]+)</label>', text)}
    mine48 = [re.sub(r"\s+[LR]$", "", names[v]).strip() for v in range(1, ROI.N_LABELS + 1)]
    units48, order = ROI.label_units(R.ATLAS_XML), ROI.unit_order(R.ATLAS_XML)
    if mine48 != units48:
        stop("label names from the logged XML differ from roi.label_units")
    if len(order) != 27:
        stop(f"{len(order)} units, not 27")
    labs = {u: [i for i, nm in enumerate(units48, start=1) if nm == u] for u in order}
    if set(ABBR) != set(order) or len(set(ABBR.values())) != len(order):
        stop(f"abbreviations do not cover the 27 units one to one: {sorted(set(ABBR) ^ set(order))}")
    firsts = [labs[u][0] for u in order]
    if firsts != sorted(firsts) or len(set(firsts)) != len(firsts):
        stop(f"unit_order is not the order of the first label numbers: {firsts}")
    side = {}
    for v in range(1, ROI.N_LABELS + 1):
        m = re.search(r"\s([LR])$", names[v].strip())
        side[v] = m.group(1) if m else "M"
    for u in order:
        s = sorted(side[i] for i in labs[u])
        if s not in (["L", "R"], ["M"]):
            stop(f"unit {u}: labels {labs[u]} sides {s}")
    return order, labs, side, names


def templates(basis):
    """Read the skeleton and the atlas (logged), check one grid, and build the unit volume on the skeleton (>= 2000)."""
    import numpy as np
    xml_raw, xml_sha = reread(R.ATLAS_XML, "JHU 라벨 이름(저장소 참조 파일 · 영상 아님) — 48 라벨 → 27 단위")
    order, labs, side, names = label_table(xml_raw)
    sk_img, sk_sha, sk_n = read_template(SKELETON_FILE, "FMRIB58 FA 골격(문턱 2000) — 슬라이스 선택(v12 기록과 대조) · 그리지 않음", basis)
    at_img, at_sha, at_n = read_template(R.ATLAS, "JHU ICBM 라벨 1 mm — 슬라이스 선택 · 단위 표시", basis)
    skel = np.asarray(sk_img.dataobj)
    atlas = np.asarray(at_img.dataobj)
    if skel.shape != atlas.shape or skel.ndim != 3:
        stop(f"grids differ: skeleton {skel.shape}, atlas {atlas.shape}")
    aff = sk_img.affine
    for img in (sk_img, at_img):
        if not np.allclose(img.affine, aff) or not np.allclose(img.get_sform(), img.get_qform()):
            stop("affines differ (skeleton / atlas, sform / qform)")
    if not (np.allclose(aff[:3, :3], np.diag(np.diag(aff[:3, :3]))) and np.allclose(np.abs(np.diag(aff[:3, :3])), 1.0)):
        stop(f"the template affine is not a 1 mm grid without rotation: {aff.tolist()}")
    if not np.array_equal(atlas, np.round(atlas)) or atlas.min() < 0 or atlas.max() > ROI.N_LABELS:
        stop("atlas values are not labels 0 – 48")
    lab = atlas.astype(int)
    smask = ROI.skeleton_mask(skel)                                  # skel >= 2000, as roi.py (roi_extract.sh:14)
    uidx = np.full(ROI.N_LABELS + 1, -1, int)
    for k, u in enumerate(order):
        for v in labs[u]:
            uidx[v] = k
    uvol = np.where(smask & (lab > 0), uidx[lab], -1)                 # unit index on skeleton ∩ label voxels, else -1
    # drawn voxels: paired units on their left label only, midline units in place
    left_or_mid = np.zeros(ROI.N_LABELS + 1, bool)
    for v in range(1, ROI.N_LABELS + 1):
        left_or_mid[v] = side[v] in ("L", "M")
    drawn = np.where(smask & (lab > 0) & left_or_mid[lab], uvol, -1)   # v12: skeleton voxels (slice choice, as recorded)
    # v12b display: the FULL label regions (not only their skeleton voxels), same hemisphere convention
    drawn_full = np.where((lab > 0) & left_or_mid[lab], uidx[lab], -1)
    right_full = (lab > 0) & ~left_or_mid[lab]                       # right labels of paired units: gray
    if np.any((drawn >= 0) & (drawn_full != drawn)):
        stop("a skeleton voxel drawn in v12 has another unit in the full-label display")
    # the FSL copy of the atlas: sha256 (and, if different, voxels and affine) against the repository copy; never drawn
    fs_img, fs_sha, fs_n = read_template(FSL_ATLAS, "FSL 의 JHU ICBM 라벨 1 mm — 저장소 사본과 sha256 대조만(표시 · 계산에 쓰지 않음)",
                                         basis)
    cmp_ = dict(repo=dict(file=os.path.relpath(R.ATLAS, REPO), sha256=at_sha, bytes=at_n),
                fsl=dict(file=os.path.relpath(FSL_ATLAS, REPO), sha256=fs_sha, bytes=fs_n), same_sha256=fs_sha == at_sha)
    atlas_diff = None
    if fs_sha != at_sha:
        fs = np.asarray(fs_img.dataobj)
        cmp_.update(same_shape=fs.shape == atlas.shape,
                    same_voxels=bool(fs.shape == atlas.shape and np.array_equal(fs, atlas)),
                    same_affine=bool(np.allclose(fs_img.affine, at_img.affine)),
                    dtype=dict(repo=str(at_img.get_data_dtype()), fsl=str(fs_img.get_data_dtype())),
                    note="sha256 differ: the repository copy (the one the ROI extraction used) is drawn")
        if fs.shape == atlas.shape:
            atlas_diff = fs != atlas
            pairs = {}
            for a, b in zip(atlas[atlas_diff].tolist(), fs[atlas_diff].tolist()):
                pairs[(a, b)] = pairs.get((a, b), 0) + 1
            fs_lab = np.round(fs).astype(int) if np.array_equal(fs, np.round(fs)) else None
            cmp_.update(n_voxels_differ=int(atlas_diff.sum()),
                        n_labelled=dict(repo=int((atlas > 0).sum()), fsl=int((fs > 0).sum())),
                        n_differ_on_skeleton=int((atlas_diff & ROI.skeleton_mask(skel)).sum()),
                        fsl_values_are_labels_0_48=bool(fs_lab is not None and fs_lab.min() >= 0 and fs_lab.max() <= ROI.N_LABELS),
                        pairs_repo_to_fsl_top=[dict(repo=a, fsl=b, n=n) for (a, b), n in
                                               sorted(pairs.items(), key=lambda t: -t[1])[:12]],
                        labels_touched=sorted({int(a) for a, _ in pairs} | {int(b) for _, b in pairs if float(b).is_integer()}))
            touched = sorted({a for a, _ in pairs})
            cmp_.update(n_voxels_of_touched_labels=dict(
                repo={str(a): int((atlas == a).sum()) for a in touched},
                fsl_same_value={str(a): int((fs == a).sum()) for a in touched},
                fsl_new_value={str(b): int((fs == b).sum()) for b in sorted({b for _, b in pairs})}))
    # sides by coordinates: every L label at x < 0 mm, every R label at x > 0 mm (centroid and voxel counts across the midline)
    ii, jj, kk = np.nonzero(lab > 0)
    xmm = aff[0, 0] * ii + aff[0, 3]
    sides = {}
    for v in range(1, ROI.N_LABELS + 1):
        sel = lab[ii, jj, kk] == v
        cx = float(xmm[sel].mean())
        wrong = int((xmm[sel] > 0).sum()) if side[v] == "L" else int((xmm[sel] < 0).sum()) if side[v] == "R" else 0
        sides[v] = dict(name=names[v], side=side[v], centroid_x_mm=cx, n_vox=int(sel.sum()), n_other_side=wrong)
        if (side[v] == "L" and not cx < 0) or (side[v] == "R" and not cx > 0):
            stop(f"label {v} «{names[v]}»: centroid x {cx:.1f} mm on the other side")
    return dict(order=order, labs=labs, side=side, names=names, skel=skel, smask=smask, lab=lab, uvol=uvol, drawn=drawn,
                drawn_full=drawn_full, right_full=right_full, atlas_compare=cmp_, atlas_diff=atlas_diff,
                aff=aff, shape=skel.shape, sides=sides,
                files=[(os.path.relpath(R.ATLAS_XML, REPO), xml_sha, len(xml_raw)),
                       (os.path.relpath(SKELETON_FILE, REPO), sk_sha, sk_n), (os.path.relpath(R.ATLAS, REPO), at_sha, at_n),
                       (os.path.relpath(FSL_ATLAS, REPO), fs_sha, fs_n)])


def n_units(sl):
    import numpy as np
    u = np.unique(sl)
    return int((u >= 0).sum())


def pick(counts, mm, target, what):
    """Most units; ties to the nearest of `target` mm; a tie left after that cannot be resolved by the rule (stop)."""
    best = max(counts.values())
    cands = sorted(k for k, n in counts.items() if n == best)
    d = {k: abs(mm[k] - target) for k in cands}
    dmin = min(d.values())
    win = [k for k in cands if math.isclose(d[k], dmin, abs_tol=1e-9)]
    if len(win) != 1:
        stop(f"{what}: slices {win} tie at {best} units and {dmin} mm from {target} mm; the rule gives no further order")
    return win[0], best, cands


def select_slices(T):
    """V12_DISPLAY_RULES §2: axial = the slice holding skeleton voxels of the most of the 27 units; sagittal = the slice with
    x < 0 by the same count; ties to the nearest of z = +20 mm and x = -20 mm. Counted on the units' skeleton voxels (both
    labels of a paired unit); the count on the drawn voxels only (left label of paired units) is recorded beside it."""
    aff, (nx, ny, nz) = T["aff"], T["shape"]
    zmm = {k: float(aff[2, 2] * k + aff[2, 3]) for k in range(nz)}
    xmm = {i: float(aff[0, 0] * i + aff[0, 3]) for i in range(nx)}
    ax = {k: n_units(T["uvol"][:, :, k]) for k in range(nz)}
    ax_drawn = {k: n_units(T["drawn"][:, :, k]) for k in range(nz)}
    left = [i for i in range(nx) if xmm[i] < 0]
    sg = {i: n_units(T["uvol"][i, :, :]) for i in left}
    sg_drawn = {i: n_units(T["drawn"][i, :, :]) for i in left}
    k, k_n, k_c = pick(ax, zmm, 20.0, "axial")
    i, i_n, i_c = pick(sg, xmm, -20.0, "sagittal (x < 0)")
    k2, _, _ = pick(ax_drawn, zmm, 20.0, "axial, drawn voxels only")
    i2, _, _ = pick(sg_drawn, xmm, -20.0, "sagittal, drawn voxels only")
    top = lambda c, mm: [(int(s), mm[s], c[s]) for s in sorted(c, key=lambda s: (-c[s], abs(mm[s] - (20.0 if mm is zmm else -20.0))))[:8]]  # noqa: E731
    return dict(axial=dict(k=int(k), z_mm=zmm[k], n_units=k_n, tied=[int(c) for c in k_c], top=top(ax, zmm),
                           drawn_only=dict(k=int(k2), z_mm=zmm[k2], n_units=ax_drawn[k2], n_units_at_k=ax_drawn[k])),
                sagittal=dict(i=int(i), x_mm=xmm[i], n_units=i_n, tied=[int(c) for c in i_c], top=top(sg, xmm),
                              drawn_only=dict(i=int(i2), x_mm=xmm[i2], n_units=sg_drawn[i2], n_units_at_i=sg_drawn[i])))


def recorded_selection():
    """The slice choice as recorded in fig4/selection.md before any share was computed (the block between the markers)."""
    p = pathlib.Path(os.environ.get("F4_SELECTION_MD") or OUT.parent / "fig4" / "selection.md")
    if not p.is_file():
        stop(f"{p} missing — run --select-slices first and record the choice")
    m = re.search(r"<!-- F4 selection -->\s*```json\s*(\{.*?\})\s*```", p.read_text(encoding="utf-8"), re.S)
    if not m:
        stop(f"{p}: no recorded selection block")
    return json.loads(m.group(1)), p


def rules_unchanged():
    for f, sha in ((RULES, RULES_SHA), (RULES_B, RULES_B_SHA), (RULES_C, RULES_C_SHA)):
        h = hashlib.sha256(f.read_bytes()).hexdigest()
        if h != sha:
            stop(f"{f.name} sha256 {h[:12]} is not the recorded {sha[:12]}")
    rec = [json.loads(ln) for ln in (V12 / "RULES_RECORD.jsonl").read_text(encoding="utf-8").splitlines() if ln.strip()]
    if not any(r.get("sha256") == RULES_B_SHA and r.get("파일", "").endswith("V12B_S43_RULES.md") for r in rec):
        stop("V12B_S43_RULES.md and its sha256 are not in RULES_RECORD.jsonl")
    if not any(r.get("sha256") == RULES_C_SHA and r.get("파일", "").endswith("V12C_S44_RULES.md") for r in rec):
        stop("V12C_S44_RULES.md and its sha256 are not in RULES_RECORD.jsonl")


def compare_v12(S, files):
    """The shares of this run against the v12 run's F4_shares.json (re-read logged): every unit and coordinate, the counts and
    the share, and the pooled counts; the tables (file, sha256) are compared and reported. Any difference in a share: stop."""
    if not V12_SHARES:
        stop("set F4_V12_SHARES to the v12 run's F4_shares.json")
    raw, h = reread(pathlib.Path(V12_SHARES), "v12 실행의 F4_shares.json — 이번 몫과 단위 · 좌표마다 같은지 대조")
    old = json.loads(raw.decode("utf-8"))
    diffs = []
    for co in JT.COORDS:
        if set(old["by_unit"][co]) != set(S["by_unit"][co]):
            stop(f"{co}: ROI units differ from the v12 shares")
        for u, o in old["by_unit"][co].items():
            n = S["by_unit"][co][u]
            if (o["n_retest"], o["n_den"], o["share"]) != (n["n_retest"], n["n_den"], n["share"]):
                diffs.append((co, u, o["n_retest"], o["n_den"], n["n_retest"], n["n_den"]))
        p, q = old["pooled"][co], S["pooled"][co]
        if (p["n_retest"], p["n_den"], p["share"]) != (q["n_retest"], q["n_den"], q["share"]):
            diffs.append((co, "pooled", p["n_retest"], p["n_den"], q["n_retest"], q["n_den"]))
    if diffs:
        stop(f"shares differ from the v12 F4_shares.json: {diffs}")
    old_t = [(t["file"], t["sha256"]) for t in old["tables"]]
    new_t = [(f, sh) for f, sh, _ in files]
    return dict(file=f"F4_shares.json of the v12 run (FIG_OUT of run {old['run']}; scratch, not kept with the drafts)", sha256=h,
                run=old["run"], n_compared=sum(len(S["by_unit"][co]) for co in JT.COORDS),
                equal=True, same_tables=old_t == new_t, n_tables=len(new_t))


# ---------------------------------------------------------------------------------------------------------- (b) shares
TB_IDS = {"x": ("N-324", "N-334"), "y": ("N-326", "N-335"), "ODI": ("N-325", "N-336")}   # pooled share, units at retest level


def load_tong():
    """judge_tong.load(R.WORK, table="주" (main analysis), fname="roi.json") step by step, with every table logged before it is read: sessions without the completion marker DONE.json are left out and listed; c10 = c10r1; the retest set = c10r1 – r3."""
    values, rep, missing, files = {}, {}, [], []
    for sb in R.TONG_SUBJECTS:
        for ses in R.TONG_SESSIONS:
            wd = R.WORK / "tong" / f"{sb}_{ses}"
            p = wd / "roi.json"
            if not (wd / "DONE.json").exists() or not p.exists():
                missing.append(f"{sb}_{ses}")
                continue
            raw, h = reread(p, f"Tong ROI 표 «주»(파생 표) — 고정 크기 시험의 단위를 ROI 단위별로 묶음 · {sb}_{ses}")
            uc = JT.unit_coords(json.loads(raw.decode("utf-8")), "주")
            files.append((os.path.relpath(p, REPO), h, len(raw)))
            if ses in JT.REPEATS:
                rep.setdefault(sb, []).append(uc)
            if ses.endswith("r1"):
                values.setdefault(JT.session_center(ses), {})[sb] = uc
    return values, rep, missing, files


def roi_shares(values, rep, order):
    """Share of units at retest level per ROI unit and coordinate: the loop of judge_tong.tb (same pairs, participants, Δ,
    s, c = t0.975(ν), √2), counted per ROI unit. The pooled counts must equal judge_tong.tb on the same tables and the ledger."""
    import itertools
    units = sorted({u for c in JT.CENTERS for sb in values.get(c, {}) for u in values[c][sb]})
    if set(units) != set(order):
        stop(f"ROI units in the tables differ from the 27 template units: {sorted(set(units) ^ set(order))}")
    s = JT.retest_sd(rep, units)
    repo_tb = JT.tb(values, s)                                         # the fixed test's own function, same inputs
    pairs = list(itertools.combinations(JT.CENTERS, 2))
    subs = sorted({sb for c in JT.CENTERS for sb in values.get(c, {})})
    per = {co: {u: [0, 0] for u in units} for co in JT.COORDS}       # [at retest level, units with Δ and s]
    dfs = {co: sorted({s[co][u][1] for u in units}) for co in JT.COORDS}
    dropped = {co: 0 for co in JT.COORDS}
    for ci, co in enumerate(JT.COORDS):
        for a, b in pairs:
            for sb in subs:
                va, vb = values.get(a, {}).get(sb, {}), values.get(b, {}).get(sb, {})
                for u in set(va) & set(vb):
                    su, nu = s[co].get(u, (None, 0))
                    if su is None:
                        dropped[co] += 1
                        continue
                    c = JT.c_of(nu)
                    per[co][u][1] += 1
                    per[co][u][0] += int(abs(vb[u][ci] - va[u][ci]) <= c * math.sqrt(2) * su)
    computed = now()
    pooled = {}
    for co in JT.COORDS:
        hit = sum(per[co][u][0] for u in units)
        den = sum(per[co][u][1] for u in units)
        r = repo_tb[co]
        if (hit, den, dropped[co]) != (r["n_retest"], r["n_den"], r["units_dropped"]) or hit / den != r["ratio"]:
            stop(f"{co}: per-ROI count pools to {hit}/{den} (dropped {dropped[co]}), judge_tong.tb gives {r}")
        i_share, i_n = TB_IDS[co]
        if not (LC.same_value(hit / den, LEDGER[i_share]["원 값"]) and LC.same_value(hit, LEDGER[i_n]["원 값"])
                and LC.same_value(den, LEDGER["N-214"]["원 값"]) and LC.same_value(r["n_total"], LEDGER["N-214"]["원 값"])):
            stop(f"{co}: pooled {hit}/{den} = {hit / den!r} differs from the ledger {i_share} {LEDGER[i_share]['원 값']} · "
                 f"{i_n} {LEDGER[i_n]['원 값']} · N-214 {LEDGER['N-214']['원 값']} — not drawn")
        if hit / den != float(LEDGER[i_share]["원 값"]):
            stop(f"{co}: pooled share {hit / den!r} is not exactly the ledger raw value {LEDGER[i_share]['원 값']}")
        if dfs[co] != [6]:
            stop(f"{co}: retest degrees of freedom {dfs[co]}, the fixed test used 6")
        pooled[co] = dict(n_retest=hit, n_den=den, share=hit / den, ledger=(i_share, i_n, "N-214"), repo_tb=r)
    n_per = {per[co][u][1] for co in JT.COORDS for u in units}
    if n_per != {len(pairs) * len(subs)}:
        stop(f"units per ROI unit {n_per}, expected {len(pairs)} pairs × {len(subs)} participants")
    out = {co: {u: dict(n_retest=per[co][u][0], n_den=per[co][u][1], share=per[co][u][0] / per[co][u][1],
                        s=s[co][u][0], df=s[co][u][1]) for u in units} for co in JT.COORDS}
    return dict(computed_at=computed, pooled=pooled, by_unit=out, pairs=len(pairs), participants=subs, units=order,
                c=JT.c_of(6))


def shares_json(S, files, missing, sel):
    return dict(what="Fig. 4(b) shares of units at retest level by ROI unit (post hoc display; no test)", run=RUN_ID,
                rules=dict(file=str(RULES.relative_to(DRAFTS)), sha256=RULES_SHA),
                rules_v12b=dict(file=str(RULES_B.relative_to(DRAFTS)), sha256=RULES_B_SHA), selection=sel,
                computed_at=S["computed_at"],
                tables=[dict(file=f, sha256=h, bytes=n) for f, h, n in files], missing_sessions=missing,
                pooled=S["pooled"], c_t0975_df6=S["c"], pairs=S["pairs"], participants=S["participants"],
                by_unit={co: {u: S["by_unit"][co][u] for u in S["units"]} for co in JT.COORDS})


def main_shares():
    """First computation of the shares (no drawing): the slice choice is asserted against selection.md first."""
    rules_unchanged()
    basis = G.require()
    T = templates(basis)
    sel = select_slices(T)
    rec, rec_path = recorded_selection()
    if (sel["axial"]["k"], sel["sagittal"]["i"]) != (rec["axial"]["k"], rec["sagittal"]["i"]):
        stop(f"slice choice {sel['axial']['k']} / {sel['sagittal']['i']} differs from {rec_path}")
    values, rep, missing, files = load_tong()
    if missing:
        stop(f"sessions without completion marker: {missing} (the fixed test listed none)")
    S = roi_shares(values, rep, T["order"])
    reread_line_numbers()
    if pycache_state() != PYCACHE_BEFORE:
        stop("a __pycache__ folder changed during the run")
    OUT.mkdir(parents=True, exist_ok=True)
    js = shares_json(S, files, missing, dict(axial=sel["axial"], sagittal=sel["sagittal"]))
    js["logged"] = LOGGED
    (OUT / "F4_shares.json").write_text(json.dumps(js, ensure_ascii=False, indent=1, default=str) + "\n", encoding="utf-8")
    print("shares computed at", S["computed_at"], "run", RUN_ID)
    for co in JT.COORDS:
        p = S["pooled"][co]
        print(f"pooled {co}: {p['n_retest']}/{p['n_den']} = {p['share']!r} == ledger {p['ledger'][0]} "
              f"{LEDGER[p['ledger'][0]]['원 값']} · {p['ledger'][1]} {LEDGER[p['ledger'][1]]['원 값']}")
    print("pycache unchanged; logged reads:", len(LOGGED))


def main_select():
    rules_unchanged()
    basis = G.require()
    T = templates(basis)
    sel = select_slices(T)
    reread_line_numbers()
    if pycache_state() != PYCACHE_BEFORE:
        stop("a __pycache__ folder changed during the run")
    print(json.dumps(dict(run=RUN_ID, guard=basis, files=T["files"], selection=sel,
                          units=[dict(unit=u, labels=T["labs"][u], sides=[T["side"][v] for v in T["labs"][u]]) for u in T["order"]],
                          label_sides={v: T["sides"][v] for v in sorted(T["sides"])}, logged=LOGGED),
                     ensure_ascii=False, indent=1, default=str))


# ========================================================================================================== style (S42)
import numpy as np  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.colors import Normalize, to_rgba  # noqa: E402
from matplotlib.patches import Polygon, Rectangle  # noqa: E402
from matplotlib.ticker import FixedLocator, FuncFormatter  # noqa: E402

FONT = "Liberation Sans"
PT_TEXT, PT_LETTER = 7.0, 8.5
LW_AXIS, LW_DATA, LW_THIN, LW_MIN = 0.5, 1.0, 0.5, 0.25
MM = 1 / 25.4
W_FIG = 105.0                     # v13g: the former panel (b) alone
SURF = "#ffffff"
INK, INK2, MUTED, GRID, NEUTRAL = FC.INK, FC.INK2, FC.MUTED, FC.GRID, FC.NEUTRAL
GREY = MUTED                      # marks without a verdict (Table 2: gray = descriptive)
SKEL_GREY = GRID                  # (b): skeleton not drawn by unit
CMS = {"x": 4.0, "y": 3.7, "ODI": 4.3}  # mark sizes of v11 (figstyle.CMS NDI · FISO · ODI)


def font_ok():
    r = subprocess.run(["fc-list", ":family=Liberation Sans", "family"], capture_output=True, text=True)
    if "Liberation Sans" not in r.stdout or FONT not in {f.name for f in font_manager.fontManager.ttflist}:
        stop("font Liberation Sans not found by fc-list or matplotlib")
    return r.stdout.strip().splitlines()[0]


plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": [FONT],
    "font.size": PT_TEXT, "axes.labelsize": PT_TEXT, "axes.titlesize": PT_TEXT,
    "xtick.labelsize": PT_TEXT, "ytick.labelsize": PT_TEXT, "legend.fontsize": PT_TEXT,
    "axes.edgecolor": MUTED, "axes.linewidth": LW_AXIS, "axes.labelcolor": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.color": INK2, "ytick.color": INK2, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "xtick.major.width": LW_AXIS, "ytick.major.width": LW_AXIS, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "xtick.minor.size": 0, "ytick.minor.size": 0, "xtick.major.pad": 2, "ytick.major.pad": 2, "axes.labelpad": 3,
    "lines.linewidth": LW_DATA, "patch.linewidth": LW_THIN, "legend.frameon": False,
    "mathtext.default": "regular", "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "lines.solid_capstyle": "round", "axes.unicode_minus": True, "image.interpolation": "nearest",
})


def new_fig(w_mm, h_mm):
    if abs(w_mm - W_FIG) > 1e-6 or h_mm > 240.0:
        stop(f"figure {w_mm} x {h_mm} mm")
    return plt.figure(figsize=(w_mm * MM, h_mm * MM))


def ax_mm(fig, x, y, w, h, **kw):
    W, H = fig.get_size_inches() / MM
    return fig.add_axes([x / W, y / H, w / W, h / H], **kw)


def canvas(fig):
    W, H = fig.get_size_inches() / MM
    c = fig.add_axes([0, 0, 1, 1], zorder=-1)
    c.set_xlim(0, W)
    c.set_ylim(0, H)
    c.axis("off")
    return c


def audit(fig, name):
    """Elsevier artwork rules as figstyle.audit (v11): width, height, text >= 7 pt inside the figure, no sub/superscripts,
    line widths 0.25 – 1.5 pt (lines and line collections)."""
    problems = []
    W, H = fig.get_size_inches() / MM
    if abs(W - W_FIG) > 0.05 or H > 240.05:
        problems.append(f"size {W:.1f} x {H:.1f} mm")
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    fb = fig.bbox
    hidden = set()
    for a in fig.axes:
        for axis in (a.xaxis, a.yaxis):
            ticks = axis.get_major_ticks() + axis.get_minor_ticks()
            drawn = {id(t.label1) for t in axis._update_ticks()} | {id(t.label2) for t in axis._update_ticks()} if a.axison else set()
            hidden.update({id(t.label1) for t in ticks if id(t.label1) not in drawn} | {id(t.label2) for t in ticks if id(t.label2) not in drawn})
    for t in fig.findobj(matplotlib.text.Text):
        s = t.get_text()
        if not s or not t.get_visible() or id(t) in hidden:
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
    for lc in fig.findobj(LineCollection):
        for lw in np.atleast_1d(lc.get_linewidth()):
            if 0 < lw < LW_MIN - 1e-6 or lw > 1.5 + 1e-6:
                problems.append(f"line collection width {lw:.2f} pt")
    if problems:
        print(f"AUDIT FAIL {name}:")
        for p in problems:
            print("  -", p)
    else:
        print(f"AUDIT PASS {name}: {W:.0f} x {H:.0f} mm, font {FONT}")
    return problems


def save(fig, name):
    probs = audit(fig, name)
    if probs:
        stop(f"{name}: {len(probs)} rule problem(s); not saved")
    base = OUT / name
    fig.savefig(str(base) + ".pdf", dpi=600, metadata={"CreationDate": None, "ModDate": None})   # voxel images at 600 dpi, nearest
    fig.savefig(str(base) + ".png", dpi=300)
    plt.close(fig)
    print("wrote", f"{base}.pdf/.png")


# ========================================================================================================== ledger guards (v11)
COORDS = ("x", "y", "ODI")
LABEL = {"x": "x (logit NDI)", "y": "y (free-water coordinate)", "ODI": "ODI"}
MARK = {"x": dict(marker="o", ms=CMS["x"]), "y": dict(marker="s", ms=CMS["y"]), "ODI": dict(marker="^", ms=CMS["ODI"])}
VERDICT = {"재검사 수준": "at retest level", "넘음": "above retest level"}     # stored magnitude-test verdicts (guard only)
T_B_IDS = TB_IDS
MIN_PANEL_GAP = 4.0               # mm between texts of (a) and (b) on the same height (a separator stands between)
TAG_SOLID = "-"
TAG_DASHED = (0, (2.2, 1.4))
TAG_CLEAR = 0.6                   # mm kept free around each role-tag frame
THRESHOLD_ID = "N-330"
INTERNAL = re.compile(r"\bT-[ABC]\b|\bTong T\b|\bV4\b|\bC[37]\b|\bL0\b|\bF2\b|stage-2|grid L1|ODI \([ab]\)|"
                      r"\b[0-9a-f]{7}\b|\b\d{5}_[A-Z0-9]+\b|osf\.io|sub-\d")

CHECKED = []                      # (ID, file, field, raw value, printed text or None, role) for F4_values.md
VISIBLE = []                      # (where, text, IDs) — every figure text that carries a number
LAYOUT = {}                       # results of the layout check, for F4_values.md
NUMBER_WORD_TAGS = []             # (manuscript text, ID) of caption tags that are number words, for F4_values.md
TAGS = []                         # role-tag Text artists (their frames are checked like texts)
MAPS = []                         # (name, axes) of the voxel maps of (b) — checked like texts for overlap


class Source:
    """One ledger copy. `rel` is the name in the ledger column «근거 파일» (basis file); the copy's sha256 must equal the ledger's."""

    def __init__(self, rel, sha):
        self.rel = rel
        self.path = SOURCES / rel.replace("../", "__up__/") / sha[:12]
        self.raw = self.path.read_bytes()
        self.sha = hashlib.sha256(self.raw).hexdigest()
        self.repo_checked = False
        if self.sha != sha:
            stop(f"ledger copy {self.path} has sha256 {self.sha[:12]}, the ledger {sha[:12]}")

    def json(self):
        return json.loads(self.raw)


_SRC = {}


def source(rel, sha):
    if (rel, sha) not in _SRC:
        _SRC[(rel, sha)] = Source(rel, sha)
    return _SRC[(rel, sha)]


def src_of(id_):
    r = LEDGER[id_]
    return source(r["근거 파일"], r["근거 파일 sha256"])


def norm(s):
    return str(s).replace("−", "-").replace(" ", "").replace(" ", "").replace(" ", "").lower()


def same(value, raw):
    try:
        return math.isclose(float(value), float(str(raw).replace("−", "-").replace(",", "")), rel_tol=1e-12, abs_tol=1e-15)
    except ValueError:
        return False


def fmt(v, k, signed=False):
    """Round half up to k decimals, as the ledger's rounding rules do; minus sign U+2212, optional plus sign."""
    d = Decimal(repr(float(v))).quantize(Decimal(1).scaleb(-k), rounding=ROUND_HALF_UP)
    s = f"{abs(d):.{k}f}"
    return ("−" + s) if d < 0 else ("+" + s) if (signed and d > 0) else s


def cap(s):
    return s[:1].upper() + s[1:]


def reproduce(id_):
    """The ledger row reproduces from its own copy (ledger_check.eval_field) and its rounding rule holds."""
    r = LEDGER.get(id_)
    if r is None:
        stop(f"{id_} is not in the ledger")
    try:
        v = LC.eval_field(SOURCES, r)
    except Exception as ex:                                   # noqa: BLE001 — any failure to reproduce is a stop
        stop(f"{id_} does not reproduce from its copy: {type(ex).__name__}: {ex}")
    if v is not None and not LC.same_value(v, r["원 값"]):
        stop(f"{id_}: copy gives {v!r}, ledger raw value {r['원 값']!r}")
    ok, why = LC.rounding_ok(r["원 값"], r["원고 표기"], r["반올림 규칙"])
    if ok is False:
        stop(f"{id_}: manuscript text «{r['원고 표기']}» does not follow the rounding rule ({why})")
    return v


def checked(id_, src, field, value, shown=None, role=""):
    """Stop unless the ledger row cites this file (same sha256) and field, its raw value equals `value`, the row reproduces
    from its copy and, for a printed number, its manuscript text equals `shown`."""
    r = LEDGER.get(id_)
    if r is None:
        stop(f"{id_} is not in the ledger")
    if r["근거 파일"] != src.rel or r["근거 필드"] != field:
        stop(f"{id_} cites {r['근거 파일']} · {r['근거 필드']}, not {src.rel} · {field}")
    if r["근거 파일 sha256"] != src.sha:
        stop(f"{id_}: sha256 of {src.path} is {src.sha[:12]}, ledger has {r['근거 파일 sha256'][:12]}")
    if not same(value, r["원 값"]):
        stop(f"{id_}: value {value!r} read from {src.rel} {field} differs from the ledger raw value {r['원 값']!r}")
    if shown is not None and norm(shown) != norm(r["원고 표기"]):
        stop(f"{id_}: printed «{shown}» differs from the ledger text «{r['원고 표기']}»")
    reproduce(id_)
    CHECKED.append((id_, src.rel, field, value, shown, role))
    return value


def get(doc, field):
    kind, _, path = field.partition(":")
    assert kind == "json", field
    cur = doc
    for m in re.finditer(r"\[(\d+)\]|([^.\[\]]+)", path):
        cur = cur[int(m.group(1))] if m.group(1) is not None else cur[m.group(2)]
    return cur


def md_number(src, field):
    """md:<line>:<fragment> — the fragment must stand on that line; its single number is the value."""
    _, line, frag = field.split(":", 2)
    text = src.raw.decode("utf-8").splitlines()[int(line) - 1]
    if frag not in text:
        stop(f"{src.rel}:{line} does not contain «{frag}»")
    nums = re.findall(r"\d+(?:\.\d+)?", frag)
    if len(nums) != 1:
        stop(f"{src.rel}:{line}: expected one number in «{frag}»")
    return float(nums[0])


def shown_text(where, text, ids):
    """Register a figure text that carries numbers: every number in it (digits or number words) must come from the manuscript
    text of a listed ledger row, and every listed row's manuscript text must stand in it."""
    have = []
    for i in ids:
        if i not in LEDGER:
            stop(f"{where}: {i} is not in the ledger")
        if norm(LEDGER[i]["원고 표기"]) not in norm(text):
            stop(f"{where}: «{text}» does not contain the ledger text «{LEDGER[i]['원고 표기']}» of {i}")
        have += LC.numbers_in(LEDGER[i]["원고 표기"], ms=True)
    for x in LC.numbers_in(text):
        if not any(math.isclose(x, y, rel_tol=1e-12, abs_tol=1e-12) for y in have):
            stop(f"{where}: number {x} in «{text}» has no listed ledger row")
    VISIBLE.append((where, text, tuple(ids)))
    return text


def read_values():
    """(a): the guards of v11 make_F3.read_values for the negative-control panel, unchanged except the panel names in the
    roles. v13: the reads and guards that only panel (c) (two Prismas) needed are gone with the panel."""
    judg = source("decomp/results/remeasure/judgements.json", LEDGER["N-324"]["근거 파일 sha256"])
    fixed = source("decomp/remeasure/MEASUREMENT_FIXED.md", LEDGER["N-330"]["근거 파일 sha256"])
    for s in (judg, fixed):                                   # the copies are the files of the current repository
        live = hashlib.sha256((REPO / s.rel).read_bytes()).hexdigest()
        if live != s.sha:
            stop(f"repository {s.rel} has sha256 {live[:12]}, its ledger copy {s.sha[:12]}")
        s.repo_checked = True
    J = judg.json()
    tong = J["Tong"]["main"]
    ta, tb = tong["T_A"], tong["T_B"]

    # --- direction test: only the undeterminable branch can be drawn
    n_pairs = checked("N-204", judg, "json:Tong.main.T_A.n_pairs", ta["n_pairs"], role="direction test, center pairs (guard only)")
    n_und = checked("N-323", judg, "json:Tong.main.T_A.n_undetermined", ta["n_undetermined"],
                    shown=f"{ta['n_undetermined']} of {ta['n_pairs']}", role="direction test, effects within the margin (note in a)")
    limit = checked("N-332", judg, "json:Tong.main.T_A.undetermined_limit", ta["undetermined_limit"],
                    role="direction test, undeterminable limit (guard only)")
    margin = checked("N-333", judg, "json:Tong.main.T_A.margin", ta["margin"], role="direction test, margin (guard only)")
    if ta["n_determined"] + n_und != n_pairs or limit != n_pairs // 2 + 1 or n_und < limit or n_und != n_pairs:
        stop(f"direction-test counts {ta}")                     # the note says 45 of 45: every pair within the margin
    if ta["judgeable"] or ta["verdict"] != "판정 불가" or ta["p"] is not None or tong["T_C"] != "판정 불가":
        stop("stored direction test / combined reading allow more than the undeterminable branch")
    if not math.isclose(margin, 1.96 * math.sqrt(2) * max(tong["s_xy_unit_median"]), rel_tol=1e-12):
        stop("direction-test margin differs from 1.96 * sqrt(2) * max(s_xy_unit_median)")

    # --- magnitude test: threshold from the fixed document, then shares, counts and verdicts per coordinate
    thr_field = LEDGER["N-330"]["근거 필드"]
    thr_read = md_number(fixed, thr_field)
    thr = checked("N-330", fixed, thr_field, thr_read, shown=fmt(thr_read, 2),
                  role="magnitude-test threshold (line and label in a; outline rule and labels in b)")
    n_tot = checked("N-214", judg, "json:Tong.main.T_B.x.n_total", tb["x"]["n_total"],
                    role="magnitude-test units per coordinate (guard; a and b — y and ODI asserted equal)")
    a = {}
    for k in COORDS:
        b = tb[k]
        if b["n_total"] != n_tot or b["n_den"] != n_tot or b["units_dropped"] != 0 or not b["judgeable"]:
            stop(f"magnitude test {k}: {b}")
        if b["verdict"] not in VERDICT:
            stop(f"magnitude test {k}: verdict {b['verdict']}")
        i_share, i_n = T_B_IDS[k]
        share = checked(i_share, judg, f"json:Tong.main.T_B.{k}.ratio", b["ratio"], shown=fmt(b["ratio"], 3),
                        role=f"magnitude-test share, {k} (mark and label in a; pooled check of b)")
        n_rt = checked(i_n, judg, f"json:Tong.main.T_B.{k}.n_retest", b["n_retest"],
                       role=f"magnitude-test units at retest level, {k} (guard; pooled check of b)")
        if not math.isclose(share, n_rt / b["n_den"], rel_tol=1e-12):
            stop(f"magnitude test {k}: share {share} != {n_rt}/{b['n_den']}")
        at = share >= thr
        if at != (b["verdict"] == "재검사 수준"):
            stop(f"magnitude test {k}: stored verdict {b['verdict']} at share {share} and threshold {thr}")
        a[k] = dict(share=share, n_rt=n_rt, n_tot=b["n_total"], at=at, verdict=b["verdict"])
    if [k for k in COORDS if not a[k]["at"]] != ["y"]:
        stop("the figure shows x and ODI at retest level and y above it")
    if tong["df_by_coord"] != {"x": [6], "y": [6], "ODI": [6]}:
        stop(f"stored retest degrees of freedom {tong['df_by_coord']}")

    for k in ("N-151", "N-015"):                               # the numbers of the (a) title
        reproduce(k)
    srcs = sorted(_SRC.values(), key=lambda s: s.rel)
    return dict(thr=thr, n_pairs=n_pairs, n_und=n_und, margin=margin, n_tot=n_tot, a=a, srcs=srcs, stored_tb=tb)


# ========================================================================================================== (b) display arrays
CROP_MARGIN = 3                   # voxels around the labelled voxels of the slice (layout)
OUTLINE_LW = 0.35                 # pt; outline of units with share < 0.90
OUTLINE_GAP = 1                   # voxels between a unit and its outline (0 = on the unit's own voxel edges; see values list)
STRIP_CW, STRIP_GX, STRIP_GY = 4.4, 0.9, 0.9    # strip cell width, gaps between columns and rows (mm)
STRIP_OUT = 0.2                   # mm between a strip cell and its outline (outside the cell, as in the maps)
STRIP_BOTTOM = 8.0                # mm; lower end of the strip (rows fill from the column headers down to here)
STRIP_LABEL_W = 8.0               # mm kept for the abbreviations left of the cells
STRIP_LABEL_GAP = 1.0             # mm between an abbreviation and its first cell


def view_arrays(T, sel):
    """Display-oriented 2-D arrays (row up, column right) of the chosen slices: unit index on the drawn voxels of the FULL
    label regions (left label of paired units, midline units in place), the right labels of paired units (gray) and, for the
    v12 comparison only, the skeleton; axial in neurological orientation (left on the left, anterior up), sagittal viewed from
    the left (anterior left). Crop: the labelled voxels of the slice and CROP_MARGIN voxels."""
    aff, (nx, ny, nz) = T["aff"], T["shape"]
    if not (aff[0, 0] < 0 < aff[1, 1] and aff[2, 2] > 0):
        stop(f"orientation of the template affine {np.diag(aff).tolist()} is not the FSL MNI one assumed for display")
    k, i = sel["axial"]["k"], sel["sagittal"]["i"]
    out = {}
    # axial: rows = j (posterior -> anterior, up), columns = i reversed (x from -91 to +90 mm: left on the left)
    out["axial"] = dict(drawn=T["drawn_full"][::-1, :, k].T, other=T["right_full"][::-1, :, k].T,
                        skel_drawn=T["drawn"][::-1, :, k].T)
    # sagittal at x < 0: rows = k (inferior -> superior, up), columns = j reversed (anterior on the left)
    out["sagittal"] = dict(drawn=T["drawn_full"][i, ::-1, :].T, other=T["right_full"][i, ::-1, :].T,
                           skel_drawn=T["drawn"][i, ::-1, :].T)
    for v in out.values():
        v["lab"] = (v["drawn"] >= 0) | v["other"]
        r, c = np.nonzero(v["lab"])
        v["box"] = (max(r.min() - CROP_MARGIN, 0), min(r.max() + CROP_MARGIN, v["lab"].shape[0] - 1),
                    max(c.min() - CROP_MARGIN, 0), min(c.max() + CROP_MARGIN, v["lab"].shape[1] - 1))
        if np.any((v["drawn"] >= 0) & v["other"]):
            stop("a voxel is both a coloured unit and a gray right label")
        sk = v["skel_drawn"] >= 0
        if np.any(v["drawn"][sk] != v["skel_drawn"][sk]):
            stop("the full-label display gives a v12 skeleton voxel another unit")
    return out


def units_not_drawn(views, order):
    """ROI units (in unit order) with no drawn voxel in either displayed slice (full labels): only the strip shows them."""
    shown = set()
    for vv in views.values():
        shown |= {int(q) for q in np.unique(vv["drawn"]) if q >= 0}
    return [u for n, u in enumerate(order) if n not in shown]


def rgba(view, share_by_index):
    img = np.zeros(view["drawn"].shape + (4,))
    img[view["other"]] = to_rgba(SKEL_GREY)
    norm_ = Normalize(SCALE_LO, 1.0, clip=True)                       # below SCALE_LO: the lightest colour
    for idx, s in share_by_index.items():
        m = view["drawn"] == idx
        if m.any():
            img[m] = FC.SEQ(norm_(s))
    return img


def outline_segments(mask):
    """Edges between voxels of `mask` and voxels outside it (row/column index coordinates, voxel centres at integers)."""
    pad = np.pad(mask, 1)
    core = pad[1:-1, 1:-1]
    segs = []
    for dr, dc in ((0, -1), (0, 1), (-1, 0), (1, 0)):
        nb = pad[1 + dr:pad.shape[0] - 1 + dr, 1 + dc:pad.shape[1] - 1 + dc]
        for r, c in zip(*np.nonzero(core & ~nb)):
            if dc:
                x = c + 0.5 * dc
                segs.append(((x, r - 0.5), (x, r + 0.5)))
            else:
                y = r + 0.5 * dr
                segs.append(((c - 0.5, y), (c + 0.5, y)))
    return segs


# ========================================================================================================== draw
def head(cv, x, y, letter, title_parts, tag_text, dashed, where):
    """Panel letter (8.5 pt bold), title (7 pt, ink) and, under the title, the role tag (7 pt capitals in a frame; solid for a
    test fixed before computing, dashed for analyses after the results). `title_parts`: text, or (text, ledger IDs)."""
    if letter:
        cv.text(x, y, f"({letter})", fontsize=PT_LETTER, fontweight="bold", color=INK, ha="left", va="baseline")
    text, ids = title_parts if isinstance(title_parts, tuple) else (title_parts, ())
    if ids:
        text = shown_text(where, text, ids)
    x_t = x + 5.4 if letter else x
    cv.text(x_t, y, text, fontsize=PT_TEXT, color=INK, ha="left", va="baseline")
    t = cv.text(x_t + 0.7, y - 3.7, tag_text, fontsize=PT_TEXT, color=INK2, ha="left", va="baseline",
                bbox=dict(boxstyle="square,pad=0.28", fc="none", ec=INK2, lw=0.5, ls=TAG_DASHED if dashed else TAG_SOLID))
    TAGS.append(t)
    return t


H_FIG = 100.0                     # v13: (a) and (b) only; the height is set by the 27-row strip of (b)
SEP = 0.0                         # v13g: no left panel and no separator
HSEP = None                       # v13: no horizontal separator (panel (c) is gone)
OY_A = H_FIG - 55.0               # top of (a) as in v11/v12 (a 55 mm panel under the top edge)
A_ROW = 8.0                       # v13: row height of the magnitude-test rows in (a) (v12: 6.0), the left column being free
A_GAP = 4.0                       # v13: extra space between the magnitude-test axis and the rule above the direction test


def draw(v, views, shares, order):
    W, H = W_FIG, H_FIG
    fig = new_fig(W, H)
    cv = canvas(fig)
    TITLE_Y = H - 4.2

    # ---------------- (b) shares by ROI unit, coloured on the full atlas labels (values from skeleton voxels), and the strip
    BX = 1.5                                               # v13g: W_FIG − BX as before (103.5 mm), so the scale is unchanged
    head(cv, BX, TITLE_Y, None, "By ROI unit, on the full atlas labels", "POST HOC", True, "(b) panel title")
    y_hdr = TITLE_Y - 9.4
    top = y_hdr - 1.6
    # strip geometry (mm): 27 rows x 3 columns at the right edge, abbreviations right-aligned before the cells
    SX1 = W - 1.0
    SX0 = SX1 - 3 * STRIP_CW - 2 * STRIP_GX
    s_pitch = (top - STRIP_BOTTOM) / len(order)
    MX0, MX1, GAPX, GAPY = BX + 4.6, SX0 - STRIP_LABEL_W - 3.0, 2.2, 2.4
    col_w = (MX1 - MX0 - 2 * GAPX) / 3
    hw ={k: (vv["box"][3] - vv["box"][2] + 1, vv["box"][1] - vv["box"][0] + 1) for k, vv in views.items()}   # (cols, rows)
    scale = col_w / max(hw["axial"][0], hw["sagittal"][0])                                                # mm per voxel
    rows_y = {"axial": top - hw["axial"][1] * scale}
    rows_y["sagittal"] = rows_y["axial"] - GAPY - hw["sagittal"][1] * scale
    idx = {u: n for n, u in enumerate(order)}
    below = {}
    for ci, co in enumerate(COORDS):
        x_col = MX0 + ci * (col_w + GAPX)
        cv.text(x_col + col_w / 2, y_hdr, co, ha="center", va="baseline", color=INK)
        share_idx = {idx[u]: shares[co][u] for u in order}
        below[co] = sorted(u for u in order if shares[co][u] < v["thr"])
        for vw in ("axial", "sagittal"):
            vv = views[vw]
            r0, r1, c0, c1 = vv["box"]
            wv, hv = (c1 - c0 + 1) * scale, (r1 - r0 + 1) * scale
            axm = ax_mm(fig, x_col + (col_w - wv) / 2, rows_y[vw], wv, hv)
            axm.imshow(rgba(vv, share_idx), origin="lower", interpolation="nearest", resample=False,
                       extent=(-0.5, vv["drawn"].shape[1] - 0.5, -0.5, vv["drawn"].shape[0] - 0.5))
            segs = []
            for u in below[co]:
                m = vv["drawn"] == idx[u]
                if OUTLINE_GAP:
                    from scipy import ndimage as ndi
                    m = ndi.binary_dilation(m, iterations=OUTLINE_GAP) if m.any() else m
                segs += outline_segments(m)
            if segs:
                axm.add_collection(LineCollection(segs, colors=INK, linewidths=OUTLINE_LW, capstyle="projecting",
                                                  joinstyle="miter", zorder=3))
            axm.set_xlim(c0 - 0.5, c1 + 0.5)
            axm.set_ylim(r0 - 0.5, r1 + 0.5)
            axm.set_aspect("equal")
            axm.axis("off")
            MAPS.append((f"(b) {co} {vw}", axm))
            if ci == 0:
                lab_l, lab_r = ("L", "R") if vw == "axial" else ("A", "P")
                yy = rows_y[vw] + 1.2
                cv.text(x_col + (col_w - wv) / 2 + 0.3, yy, lab_l, ha="left", va="center", color=INK2)
                cv.text(x_col + (col_w + wv) / 2 - 0.3, yy, lab_r, ha="right", va="center", color=INK2)
    for vw in ("axial", "sagittal"):
        hv = (views[vw]["box"][1] - views[vw]["box"][0] + 1) * scale
        cv.text(BX + 2.3, rows_y[vw] + hv / 2, "axial" if vw == "axial" else "sagittal, left", rotation=90,
                ha="center", va="center", color=INK2)
    # strip: every ROI unit, atlas label order, same colour scale, the same outline below the threshold; no numbers
    n_rows = len(order)
    if n_rows != int(float(LEDGER["N-179"]["원 값"])):
        stop(f"the strip has {n_rows} rows, the ledger N-179 {LEDGER['N-179']['원 값']} ROI units")
    reproduce("N-179")
    s_bot = top - n_rows * s_pitch
    PADS = STRIP_OUT + 0.3
    sax = ax_mm(fig, SX0 - PADS, s_bot - PADS, SX1 - SX0 + 2 * PADS, top - s_bot + 2 * PADS)
    sax.set_xlim(SX0 - PADS, SX1 + PADS)
    sax.set_ylim(s_bot - PADS, top + PADS)
    sax.axis("off")
    norm_ = Normalize(SCALE_LO, 1.0, clip=True)                       # same scale as the maps
    strip_rows = []
    for r, u in enumerate(order):
        y_c = top - (r + 0.5) * s_pitch
        hc = s_pitch - STRIP_GY
        for ci, co in enumerate(COORDS):
            x0 = SX0 + ci * (STRIP_CW + STRIP_GX)
            sax.add_patch(Rectangle((x0, y_c - hc / 2), STRIP_CW, hc, facecolor=FC.SEQ(norm_(shares[co][u])), edgecolor="none",
                                    lw=0))
            if shares[co][u] < v["thr"]:
                sax.add_patch(Rectangle((x0 - STRIP_OUT, y_c - hc / 2 - STRIP_OUT), STRIP_CW + 2 * STRIP_OUT,
                                        hc + 2 * STRIP_OUT, facecolor="none", edgecolor=INK, lw=OUTLINE_LW,
                                        joinstyle="miter"))
        cv.text(SX0 - STRIP_LABEL_GAP, y_c, ABBR[u], ha="right", va="center", color=INK)
        strip_rows.append((r + 1, u, ABBR[u]))
    for ci, co in enumerate(COORDS):
        cv.text(SX0 + ci * (STRIP_CW + STRIP_GX) + STRIP_CW / 2, y_hdr, co, ha="center", va="baseline", color=INK)
    cv.text(SX1, y_hdr + 3.6, "all ROI units", ha="right", va="baseline", color=INK2)
    MAPS.append(("(b) strip", sax))
    # colour scale and key (under the maps)
    # V12C §3: the bar spans SCALE_LO – 1.0; at its lower end a triangle in the lightest colour (extend min: shares below
    # SCALE_LO take that colour)
    y_cb = rows_y["sagittal"] - 9.0
    CB_H, CB_TRI = 2.6, 3.0                                # bar height, length of the extend-min triangle (mm)
    cb_x0 = MX0 + CB_TRI
    cb_w = MX1 - cb_x0
    cv.add_patch(Polygon([(MX0, y_cb + CB_H / 2), (cb_x0, y_cb), (cb_x0, y_cb + CB_H)], closed=True,
                         facecolor=FC.SEQ(0.0), edgecolor="none", lw=0))
    cax = ax_mm(fig, cb_x0, y_cb, cb_w, CB_H)
    cax.imshow(np.linspace(SCALE_LO, 1.0, 512)[None, :], cmap=FC.SEQ, vmin=SCALE_LO, vmax=1.0, aspect="auto",
               interpolation="nearest", extent=(SCALE_LO, 1.0, 0.0, 1.0))
    cax.axvline(v["thr"], color=INK, lw=0.8)
    cax.set_xlim(SCALE_LO, 1.0)
    cax.set_yticks([])
    for sp in ("left", "top", "right"):
        cax.spines[sp].set_visible(False)
    cax.xaxis.set_major_locator(FixedLocator([0.5, 0.6, 0.7, 0.8, 1.0]))
    cax.xaxis.set_major_formatter(FuncFormatter(lambda t, _: "1" if t == 1 else f"{t:.2f}"))
    cax.set_xlabel("share of units at retest level, by ROI unit")
    x_t = cb_x0 + (v["thr"] - SCALE_LO) / (1.0 - SCALE_LO) * cb_w
    cv.text(x_t, y_cb + 3.3, shown_text("(b) threshold mark above the colour scale", fmt(v["thr"], 2), [THRESHOLD_ID]),
            ha="center", va="bottom", color=INK)
    # keys under the colour scale: an outlined unit (a coloured patch, its outline outside it), the gray right labels
    KX = MX0
    KW, KH = 3.0, 1.6                                      # key patch, mm (schematic, not to scale)
    KG = 0.45 if OUTLINE_GAP else 0.0
    y_k1 = y_cb - 11.2
    y_k2 = y_k1 - 4.0
    kax = ax_mm(fig, KX, y_k1 - 1.4, 4.4, 2.8)
    kax.set_xlim(0.0, 4.4)
    kax.set_ylim(0.0, 2.8)
    kax.add_patch(Rectangle((0.7, 1.4 - KH / 2), KW, KH, facecolor=FC.SEQ(0.7), edgecolor="none"))
    kax.add_patch(Rectangle((0.7 - KG, 1.4 - KH / 2 - KG), KW + 2 * KG, KH + 2 * KG, facecolor="none", edgecolor=INK,
                            lw=OUTLINE_LW))
    kax.axis("off")
    cv.text(KX + 5.4, y_k1, shown_text("(b) key, outline", f"outlined: below {fmt(v['thr'], 2)} (maps and strip)",
                                       [THRESHOLD_ID]), ha="left", va="center", color=INK)
    kax2 = ax_mm(fig, KX, y_k2 - 1.4, 4.4, 2.8)
    kax2.set_xlim(0.0, 4.4)
    kax2.set_ylim(0.0, 2.8)
    kax2.add_patch(Rectangle((0.7, 1.4 - KH / 2), KW, KH, facecolor=SKEL_GREY, edgecolor="none"))
    kax2.axis("off")
    cv.text(KX + 5.4, y_k2, "gray: right label of paired units", ha="left", va="center", color=INK)
    cv.text(MX0, y_k2 - 4.6, "paired units on their left label; midline units in place", ha="left", va="center", color=INK2)
    cv.text(MX0, y_k2 - 8.0, "labels drawn in full; values from skeleton voxels", ha="left", va="center", color=INK2)
    not_drawn = units_not_drawn(views, order)                # ROI units the two slices cannot show: only the strip shows them
    geo = dict(SEP=SEP, HSEP=HSEP, x_thr=None, share_labels=[], scale=scale, below=below,
               rows_y=rows_y, col_w=col_w, cb=(cb_x0, y_cb, cb_w, SCALE_LO), not_drawn=not_drawn, strip_rows=strip_rows,
               strip=dict(x0=SX0, x1=SX1, top=top, bottom=s_bot, pitch=s_pitch, cell_w=STRIP_CW,
                          cell_h=s_pitch - STRIP_GY, out=STRIP_OUT), maps_x=(MX0, MX1))
    return fig, geo


# ========================================================================================================== checks on the drawn figure
def figure_texts(fig):
    """Visible texts with their boxes in mm; tick labels marked (an axis draws only some of its tick objects)."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    k = 25.4 / fig.dpi
    all_ticks, drawn_ticks = set(), set()
    for a in fig.axes:
        for axis in (a.xaxis, a.yaxis):
            for tk in axis.get_major_ticks() + axis.get_minor_ticks():
                all_ticks.update({id(tk.label1), id(tk.label2)})
            if a.axison:
                for tk in axis._update_ticks():
                    drawn_ticks.update({id(tk.label1), id(tk.label2)})
    out = []
    for t in fig.findobj(matplotlib.text.Text):
        s = t.get_text()
        if not s or not t.get_visible() or (id(t) in all_ticks and id(t) not in drawn_ticks):
            continue
        bb = t.get_window_extent(r)
        out.append(dict(obj=t, text=s, tick=id(t) in all_ticks, x0=bb.x0 * k, x1=bb.x1 * k, y0=bb.y0 * k, y1=bb.y1 * k))
    return out


def check_figure(fig, geo):
    texts = figure_texts(fig)
    r = fig.canvas.get_renderer()
    k = 25.4 / fig.dpi
    for t in TAGS:                                             # tag frames count as texts, with a clearance of TAG_CLEAR mm
        bb = t.get_bbox_patch().get_window_extent(r)
        for e in texts:
            if e["obj"] is t:
                e.update(x0=bb.x0 * k - TAG_CLEAR, x1=bb.x1 * k + TAG_CLEAR, y0=bb.y0 * k - TAG_CLEAR, y1=bb.y1 * k + TAG_CLEAR)
    registered = {t for _, t, _ in VISIBLE}
    for t in texts:
        if not t["tick"] and LC.numbers_in(t["text"]) and t["text"] not in registered:
            stop(f"figure text «{t['text']}» carries a number with no ledger row listed")
        if INTERNAL.search(t["text"]):
            stop(f"figure text «{t['text']}» names an internal test, arm, hash, session or participant")
        forb = [m.group(0) for rx, _ in LC.FORBIDDEN for m in rx.finditer(t["text"])]
        if forb:
            stop(f"figure text «{t['text']}» uses words not to use: {forb}")
    drawn = {t["text"] for t in texts}
    for where, t, _ in VISIBLE:
        if t not in drawn:
            stop(f"{where}: «{t}» was registered but is not in the figure")
    # maps of (b): their drawn boxes take part in the overlap test (no text on a map, no map on a map)
    boxes = []
    for name, a in MAPS:
        bb = a.get_window_extent(r)
        boxes.append(dict(obj=a, text=f"[map {name}]", tick=False, x0=bb.x0 * k, x1=bb.x1 * k, y0=bb.y0 * k, y1=bb.y1 * k))
    items = texts + boxes
    tol = 0.05
    for i, p in enumerate(items):
        for q in items[i + 1:]:
            if p["x0"] < q["x1"] - tol and q["x0"] < p["x1"] - tol and p["y0"] < q["y1"] - tol and q["y0"] < p["y1"] - tol:
                if p["text"] in ("L", "R", "A", "P") or q["text"] in ("L", "R", "A", "P"):
                    continue                                   # orientation letters sit in the empty corner of their map
                stop(f"overlap: «{p['text']}» and «{q['text']}»")
    for t in texts:                                            # share labels in (a) stay clear of the threshold line
        if any(t["obj"] is s for s in geo["share_labels"]):
            if t["x0"] - 0.5 < geo["x_thr"] < t["x1"] + 0.5:
                stop(f"share label «{t['text']}» touches the threshold line")
    left = [t for t in items if t["x1"] <= geo["SEP"]]
    right = [t for t in items if t["x0"] >= geo["SEP"]]
    if len(left) + len(right) != len(items):
        stop("a text crosses the vertical separator: " + ", ".join(t["text"] for t in items if t not in left and t not in right))
    cross = [t["text"] for t in left if geo["HSEP"] is not None and t["y0"] < geo["HSEP"] < t["y1"]]
    if cross:
        stop(f"texts cross the horizontal separator: {cross}")
    gaps = [(q["x0"] - p["x1"], p["text"], q["text"]) for p in left for q in right if p["y0"] < q["y1"] and q["y0"] < p["y1"]]
    gmin = min(gaps) if gaps else (float("inf"), "", "")    # v13g: one column
    if gmin[0] < MIN_PANEL_GAP:
        stop(f"left column – (b) gap {gmin[0]:.1f} mm between «{gmin[1]}» and «{gmin[2]}» is below {MIN_PANEL_GAP} mm")
    small = [(t["text"], t["obj"].get_fontsize()) for t in texts if t["obj"].get_fontsize() < PT_TEXT - 1e-6]
    if small:
        stop(f"text below {PT_TEXT} pt: {small}")
    LAYOUT.update(n_texts=len(texts), n_maps=len(boxes), gap=gmin, sep=geo["SEP"], hsep=geo["HSEP"],
                  scale=geo["scale"], fig_h=H_FIG)
    print(f"layout ok: {len(texts)} texts (tag frames included) and {len(boxes)} maps, no overlap; nearest texts across the "
          f"separator on one height {gmin[0]:.1f} mm apart («{gmin[1]}» · «{gmin[2]}»); map scale {geo['scale']:.3f} mm/voxel")


# ========================================================================================================== caption and value list
# (a) and the title: the S41 §3 caption as decided. (b): the phrases V12_DISPLAY_RULES prescribes, verbatim («by ROI unit;
# post hoc, descriptive; the verdict is in (a)» · «outlined, below 0.90, the threshold of the pooled fixed test; no verdict by
# ROI unit» · «displayed after the verdict; no test»); v12b (V12B_S43_RULES §4, direction S43 §2.2): «drawn on the full atlas
# labels; values from skeleton voxels», the strip with all N-179 ROI units, and the pointer to Supplementary Table S18; the v12
# clause «ROI units outside the slices shown are not drawn» is gone because the strip shows every unit (checked in draw()).
# v13 (direction S45, author decision 2026-11-08): panel (c) and its sentence are gone (the same measured difference is in the
# new Fig. 7b); the rest of the caption is unchanged.
CAPTION = ('![**Fig. 4. Within {N-015} scanner model and protocol, by ROI unit** (post hoc, descriptive; no test). Share of units at '
           'retest level on the full atlas labels of an axial and a sagittal slice; values from skeleton voxels. The strip gives all '
           '{N-179} ROI units (Supplementary Table S18); outlined, below {N-330}, the threshold of the pooled fixed test; the verdict '
           'is in Table 1.](figures/F4_negctrl_roi_prisma.pdf){#fig:fig4 width=100%}')
CAPTION_B_PHRASES = ["(post hoc, descriptive; no test)", "on the full atlas labels", "values from skeleton voxels",
                     "The strip gives all 27 [N-179] ROI units", "Supplementary Table S18",
                     "outlined, below 0.90 [N-330], the threshold of the pooled fixed test", "the verdict is in Table 1"]
OLD_NOT_DRAWN = ("ROI units outside the slices shown are not drawn", "ROI units outside these slices are not drawn")
OLD_C = ("Prisma scanners", "(c)", "σ_d")                      # v13: nothing of panel (c) may stay in the caption
NUMBER_WORDS = sorted(set(LC.WORDS) | {"single", "once", "twice"}, key=len, reverse=True)


def number_words(cap_text):
    """Number words outside tagged ledger text. Returns the untagged hits and the tagged rows whose manuscript text is or
    contains a number word."""
    spans, tagged = [], []
    for m in LC.TAG.finditer(cap_text):
        i = m.group(1)
        end = m.start()
        while end > 0 and cap_text[end - 1] in "  ":
            end -= 1
        st = LC.suffix_start(cap_text, end, LEDGER[i]["원고 표기"])
        if st is None:
            stop(f"caption: text before {i} is not its ledger text")
        spans.append((st, m.end()))
        if any(re.search(rf"(?<!\w){re.escape(w)}(?!\w)", LEDGER[i]["원고 표기"], re.I) for w in NUMBER_WORDS):
            tagged.append((LEDGER[i]["원고 표기"], i))
    masked = list(cap_text)
    for a, b in spans:
        masked[a:b] = " " * (b - a)
    masked = "".join(masked)
    for ex in LC.EXCEPTIONS:
        masked = ex.sub(lambda mm: re.sub(r"[^\n]", " ", mm.group(0)), masked)
    hits = [(m.group(0), cap_text[max(0, m.start() - 30):m.end() + 30]) for w in NUMBER_WORDS
            for m in re.finditer(rf"(?<!\w){re.escape(w)}(?!\w)", masked, re.I)]
    return hits, tagged


def caption():
    def sub(m):
        up, i = m.group(1), m.group(2)
        if i not in LEDGER:
            stop(f"caption: {i} is not in the ledger")
        reproduce(i)
        t = LEDGER[i]["원고 표기"]
        return f"{cap(t) if up else t} [{i}]"
    return re.sub(r"\{(\^?)(N-\d{3})\}", sub, CAPTION)


def ledger_check(path, cap_text):
    """Run the project's ledger check on the caption with every region strict; stop unless it passes cleanly."""
    r = subprocess.run([sys.executable, "-I", "-B", str(LEDGER_DIR / "ledger_check.py"), str(path), "--strict-all", "1"],
                       capture_output=True, text=True)
    summary = r.stdout.strip().splitlines()[0] if r.stdout.strip() else r.stderr.strip()
    if r.returncode != 0:
        stop(f"ledger_check on {path.name}: {r.stdout.strip()} {r.stderr.strip()}")
    counts = {k: int(n) for k, n in re.findall(r"(실패|엄격 구역|쓰지 않을 말) (\d+)", summary)}
    if counts.get("실패") != 0 or counts.get("엄격 구역") != 0 or counts.get("쓰지 않을 말") != 0:
        stop(f"ledger_check summary not clean: {summary}")
    forb = [m.group(0) for rx, _ in LC.FORBIDDEN for m in rx.finditer(cap_text)]
    if forb:
        stop(f"caption uses words not to use: {forb}")
    warns = []
    for i, row in LEDGER.items():
        ok, why = LC.rounding_ok(row["원 값"], row["원고 표기"], row["반올림 규칙"])
        if ok is None and why != "숫자 아닌 값":
            warns.append(i)
    return summary, warns


def values_md(v, S, T, sel, files, cap_text, check_line, warns, cached, cmpv12, s18):
    used_fig = sorted({i for _, _, ids in VISIBLE for i in ids} | {c[0] for c in CHECKED})
    used_cap = [m.group(1) for m in re.finditer(r"\[(N-\d{3})\]", cap_text)]
    L = ["# Figure 4 (v13) — drawn values and visible numbers", "",
         f"Made by `make_F4.py` (drafts/x15/ms/S41_files/src), run {RUN_ID}"
         + (" (shares from the cache of an earlier logged run — preview only)" if cached else "") + f". Ledger: "
         f"`drafts/x15/ledger/numbers.csv` ({len(LEDGER)} rows). Rules: `{RULES.relative_to(DRAFTS)}` (sha256 {RULES_SHA[:12]}…) "
         f"and, for (b) in v12b, `{RULES_B.relative_to(DRAFTS)}` §0 · §3 (sha256 {RULES_B_SHA[:12]}…) and, for its colour scale in "
         f"v12c, `{RULES_C.relative_to(DRAFTS)}` §3 (sha256 {RULES_C_SHA[:12]}…), all asserted on this run. "
         "v13: panels (a) and (b) only (panel (c) and its reads are gone). (a): every row below was asserted by the script as in v11 "
         "(same file and sha256 and field as the ledger row, value equal to the row's raw value, the row reproduced from its own "
         "ledger copy, printed text equal to the row's manuscript text). (b): shares by ROI unit, drawn (not printed); their "
         "pooled counts were asserted equal to the ledger rows of the fixed test and every share equal to the v12 run's; the same "
         "81 shares are written for Supplementary Table S18 (out/F4_S18.json), where they get ledger rows.", "",
         "## Sources read for (a) (ledger copies; nothing under remeasure/ opened for this panel)", "",
         "| file (ledger name) | read from | sha256 |", "|---|---|---|"]
    for s in v["srcs"]:
        where = "ledger copy `sources/" + str(s.path.relative_to(SOURCES)) + "`"
        if s.repo_checked:
            where += "; the repository file has the same sha256"
        L.append(f"| `{s.rel}` | {where} | `{s.sha}` |")
    L += ["", "## (a): values drawn and/or printed", "",
          "| ledger ID | printed in the figure | raw value | file · field | sha12 | role |", "|---|---|---|---|---|---|"]
    for id_, rel, field, val, shown, role in CHECKED:
        sha12 = LEDGER[id_]["근거 파일 sha256"][:12]                # asserted equal to the copy read (lit rows: the ledger's)
        L.append(f"| {id_} | {shown if shown is not None else '— (not printed)'} | {val!r} | `{rel}` `{field}` | `{sha12}` | "
                 f"{role} |")
    L += ["", "## (b): template, atlas, slices and tables", "",
          "| file (relative to the repository) | sha256 | bytes | use |", "|---|---|---|---|"]
    use = {"JHU-labels.xml": "label names: 48 labels -> 27 units (roi.label_units · unit_order); reread_log",
           "FMRIB58_FA-skeleton_1mm.nii.gz": "skeleton >= 2000: slice choice as recorded in v12 (not drawn in v12b); "
                                             "EXTERNAL_READ_LOG after guard.require()",
           "prereg/artifacts/atlas/JHU-ICBM-labels-1mm.nii.gz": "JHU labels of the ROI extraction: slice choice and the full label "
                                                                "regions drawn; EXTERNAL_READ_LOG after guard.require()",
           "JHU-ICBM-labels-1mm.nii.gz": "FSL copy of the JHU labels: sha256 compared with the repository copy only (not drawn); "
                                         "EXTERNAL_READ_LOG after guard.require()"}
    for f, h, n in T["files"]:
        L.append(f"| `{f}` | `{h}` | {n:,} | {use.get(f, use.get(pathlib.Path(f).name))} |")
    for f, h, n in files:
        L.append(f"| `{f}` | `{h}` | {n:,} | ROI table «주» (derived; re-read logged) |")
    ac = T["atlas_compare"]
    L += ["", "Atlas copies: repository `" + ac["repo"]["file"] + f"` sha256 `{ac['repo']['sha256']}` ({ac['repo']['bytes']:,} bytes); "
          "FSL `" + ac["fsl"]["file"] + f"` sha256 `{ac['fsl']['sha256']}` ({ac['fsl']['bytes']:,} bytes) — "
          + ("the same file." if ac["same_sha256"] else
             f"different files (same voxels: {ac['same_voxels']}; same affine: {ac['same_affine']}); the repository copy, the one "
             "the ROI extraction used, is drawn. Details (computed by this run from the logged bytes): `"
             + json.dumps({k: ac.get(k) for k in ("dtype", "n_voxels_differ", "n_labelled", "n_differ_on_skeleton",
                                                  "n_differ_axial_slice", "n_differ_sagittal_slice",
                                                  "fsl_values_are_labels_0_48", "pairs_repo_to_fsl_top", "labels_touched",
                                                  "n_voxels_of_touched_labels")},
                          ensure_ascii=False, default=str) + "`."), "",
          f"Slices (rule V12_DISPLAY_RULES §2, recorded in v12 fig4/selection.md before any share was computed; asserted on this "
          f"run; V12B_S43_RULES §3 keeps them): axial voxel k = {sel['axial']['k']} (z = {sel['axial']['z_mm']:+.0f} mm), "
          f"sagittal voxel i = {sel['sagittal']['i']} (x = {sel['sagittal']['x_mm']:+.0f} mm). Not printed in the figure (no "
          "ledger row).", "",
          f"Shares computed at {S['computed_at']}. c = t0.975(6) = {S['c']!r}; {S['pairs']} center pairs × "
          f"{len(S['participants'])} participants = {S['pairs'] * len(S['participants'])} units per ROI unit and coordinate.", "",
          "### Pooled check (asserted before drawing)", "", "| coordinate | pooled units at retest level / units | pooled share | "
          "ledger | judge_tong.tb on the same tables |", "|---|---|---|---|---|"]
    for co in COORDS:
        p = S["pooled"][co]
        L.append(f"| {co} | {p['n_retest']} / {p['n_den']} | {p['share']!r} | {p['ledger'][0]} {LEDGER[p['ledger'][0]]['원 값']} · "
                 f"{p['ledger'][1]} {LEDGER[p['ledger'][1]]['원 값']} · N-214 {LEDGER['N-214']['원 값']} | "
                 f"{p['repo_tb']['n_retest']} / {p['repo_tb']['n_den']} |")
    L += ["", "### Comparison with the v12 shares (asserted before drawing)", "",
          f"`{cmpv12['file']}` (sha256 `{cmpv12['sha256']}`, v12 run {cmpv12['run']}, re-read logged): all "
          f"{cmpv12['n_compared']} unit × coordinate shares (n_retest, n_den, share) and the pooled counts equal; the 36 tables "
          f"(file, sha256) {'the same' if cmpv12['same_tables'] else 'NOT the same'} as in v12.", "",
          "### Shares by ROI unit — drawn in the maps and the strip; printed only in Table S18 (ledger rows there)", "",
          "Source of every value in this table: the 36 ROI tables listed above (file and sha256 per table), field «주».<unit>."
          "NDI · FWF · ODI -> x · y · ODI (judge_tong.unit_coords); units at retest level / units with Δ and s per ROI unit, "
          "s = the c10 retest SD of the unit (judge_tong.retest_sd, 6 degrees of freedom). Colour = share on SEQ 0 – 1; "
          "* = below 0.90, outlined (maps and strip). Row = strip row (atlas label order).", "",
          "| row | ROI unit | abbr. | first label | x | y | ODI | s x | s y | s ODI | in axial map | in sagittal map |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    k_ax, i_sg = sel["axial"]["k"], sel["sagittal"]["i"]
    ax_units = {int(q) for q in np.unique(T["drawn_full"][:, :, k_ax]) if q >= 0}
    sg_units = {int(q) for q in np.unique(T["drawn_full"][i_sg, :, :]) if q >= 0}
    for n, u in enumerate(T["order"]):
        cells = [f"{S['by_unit'][co][u]['n_retest']}/{S['by_unit'][co][u]['n_den']} = {S['by_unit'][co][u]['share']:.4f}"
                 + (" *" if S["by_unit"][co][u]["share"] < v["thr"] else "") for co in COORDS]
        sds = [f"{S['by_unit'][co][u]['s']:.4f}" for co in COORDS]
        L.append(f"| {n + 1} | {u} | {ABBR[u]} | {T['labs'][u][0]} | " + " | ".join(cells) + " | " + " | ".join(sds)
                 + f" | {'yes' if n in ax_units else '—'} | {'yes' if n in sg_units else '—'} |")
    nd = [u for n, u in enumerate(T["order"]) if n not in ax_units | sg_units]
    if nd != LAYOUT.get("not_drawn"):
        stop(f"units in neither map differ between the drawing ({LAYOUT.get('not_drawn')}) and the value list ({nd})")
    L += ["", f"ROI units with no voxel in either map slice (full labels; the strip shows them): {len(nd)} of {len(T['order'])}"
          + (": " + "; ".join(f"{ABBR[u]} ({u})" for u in nd) if nd else "") + ". The v12 note «ROI units outside these slices "
          "are not drawn» and its caption clause are gone because the strip shows all units (asserted: strip rows = the 27 units "
          "of the label file in atlas order = N-179).", "",
          "| coordinate | ROI units below 0.90 (strip outlines) |", "|---|---|"]
    for co in COORDS:
        L.append(f"| {co} | {len(s18['below'][co])}: {', '.join(s18['below'][co])} |")
    L += ["", "### Table S18 outputs (this run)", "",
          f"- `{os.path.relpath(s18['json'], DRAFTS)}` sha256 `{s18['json_sha256']}` — basis file of the ledger rows (units.<abbr>."
          "<x|y|ODI>.share raw, .text printed; denominator.value). "
          + ("Kept: re-read (logged) and its values (" + ", ".join(S18_VALUE_KEYS) + ") equal this run's, so it was not "
             "rewritten (its «run», «written» and «script» fields are those of the run that wrote it)." if s18["kept"]
             else "Written by this run."),
          f"- `{s18['md']}` sha256 `{s18['md_sha256']}` — pandoc block, {s18['n_tokens']} ID tokens (81 shares and the "
          f"denominator {s18['den']}). " + ("Kept: same text, not rewritten." if s18["md_kept"] else "Written by this run."), "",
          f"Colour scale of (b) (V12C_S44_RULES §3): figcolors.SEQ on {SCALE_LO} – 1.0 in the maps, the strip and the colour bar; "
          f"shares below {SCALE_LO} would take the lightest colour (triangle at the lower end of the bar); the lowest share drawn "
          f"is {min(S['by_unit'][co][u]['share'] for co in COORDS for u in T['order']):.4f}, so no unit is below {SCALE_LO}.", "",
          "## Visible numbers in the figure, one line each", "", "| where | text | ledger ID |", "|---|---|---|"]
    for where, t, ids in VISIBLE:
        L.append(f"| {where} | {t} | {' · '.join(ids)} |")
    L += ["| (a) axis tick labels; (b) colour-scale tick labels 0.50 … 1 | | axis tick labels (exempt) |",
          "| (b) strip | no numbers (abbreviations and colour only) | — |", "",
          "## Numbers in the caption", "", "| text | ledger ID | grade in the ledger |", "|---|---|---|"]
    for i in used_cap:
        L.append(f"| {LEDGER[i]['원고 표기']} | {i} | {LEDGER[i]['근거 등급']} |")
    L += ["| S18 (in «Table S18») | table number, exempt in `ledger_check.py` (EXCEPTIONS) | — |", "",
          "Tagged number words: " + ", ".join(f"«{t}» {i}" for t, i in NUMBER_WORD_TAGS) + ".", "",
          "## Ledger check of the caption", "", "`ledger_check.py F4_caption.md --strict-all 1` (run by `make_F4.py`):", "",
          f"    {check_line}", "",
          f"Of the {n_warn_summary(check_line)} warnings, {len(warns)} are ledger-wide rows whose rounding rule is not "
          f"machine-readable ({', '.join(warns)}; "
          f"{'none' if not set(warns) & set(used_fig + used_cap) else ', '.join(sorted(set(warns) & set(used_fig + used_cap)))} "
          f"of them used in F4) and {n_warn_summary(check_line) - len(warns)} the format warning that a file holding only a "
          "caption gets for having no Highlights (it is no warning in the manuscript copy).", "",
          "## Layout check", "",
          f"Figure {W_FIG:.0f} × {H_FIG:.0f} mm. "
          f"{LAYOUT['n_texts']} figure texts (role-tag frames measured as their frames) and {LAYOUT['n_maps']} drawn boxes (six "
          f"voxel maps and the strip), no two overlapping (orientation letters sit in the empty corner of their map); the nearest "
          f"texts across the vertical separator on one height are «{LAYOUT['gap'][1]}» and «{LAYOUT['gap'][2]}», "
          f"{LAYOUT['gap'][0]:.1f} mm apart (limit {MIN_PANEL_GAP} mm). Separator at x = {LAYOUT['sep']} mm (v13: no horizontal "
          f"separator). (a): magnitude-test rows {A_ROW} mm apart, {A_GAP} mm more above the direction-test rule than v12. Maps at {LAYOUT['scale']:.3f} mm per voxel (axial and sagittal at one scale), "
          "voxels drawn with interpolation «nearest» (each voxel a uniform square; the PDF stores the maps at 600 dpi by nearest "
          "replication), no smoothing. Coloured voxels: the full label region of the unit in the slice (left label of a paired "
          "unit, a midline unit in place); gray: the right labels of paired units. "
          + (f"Outline {OUTLINE_LW} pt along the edge of the unit's drawn voxels dilated by {OUTLINE_GAP} voxel in the slice "
             "plane, as in v12: an annotation only (the coloured voxels are unchanged), drawn outside the unit so that its colour "
             "stays visible." if OUTLINE_GAP else f"Outline {OUTLINE_LW} pt on the unit's own voxel edges.")
          + f" Strip: {len(T['order'])} rows × 3 columns, cells {LAYOUT['strip']['cell_w']:.1f} × "
          f"{LAYOUT['strip']['cell_h']:.2f} mm at a row pitch of {LAYOUT['strip']['pitch']:.2f} mm, x {LAYOUT['strip']['x0']:.1f} – "
          f"{LAYOUT['strip']['x1']:.1f} mm; outline {OUTLINE_LW} pt, {LAYOUT['strip']['out']} mm outside the cell.", ""]
    return "\n".join(L) + "\n"


def n_warn_summary(line):
    m = re.search(r"경고 (\d+)", line)
    return int(m.group(1)) if m else -1


# ========================================================================================================== main
CACHE = os.environ.get("F4_CACHE")    # preview only: a cache folder written by an earlier logged run of this script


def write_cache(T, sel, S, files, cmpv12):
    d = pathlib.Path(CACHE)
    d.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(d / "templates.npz", drawn=T["drawn"], smask=T["smask"], aff=T["aff"], drawn_full=T["drawn_full"],
                        right_full=T["right_full"])
    (d / "meta.json").write_text(json.dumps(dict(run=RUN_ID, rules=RULES_SHA, rules_b=RULES_B_SHA, order=T["order"],
                                                 labs=T["labs"], files=T["files"], sel=sel, tables=files, S=S,
                                                 atlas_compare=T["atlas_compare"], v12=cmpv12),
                                            ensure_ascii=False, default=str), encoding="utf-8")


def read_cache():
    d = pathlib.Path(CACHE)
    meta = json.loads((d / "meta.json").read_text(encoding="utf-8"))
    if meta["rules"] != RULES_SHA or meta.get("rules_b") != RULES_B_SHA:
        stop("cache made under other rules")
    z = np.load(d / "templates.npz")
    T = dict(drawn=z["drawn"], smask=z["smask"], aff=z["aff"], drawn_full=z["drawn_full"], right_full=z["right_full"],
             shape=z["drawn"].shape, order=meta["order"], labs=meta["labs"], files=meta["files"],
             atlas_compare=meta["atlas_compare"])
    S = meta["S"]
    for co in COORDS:                                         # the pooled check again, from the cached counts
        hit = sum(S["by_unit"][co][u]["n_retest"] for u in T["order"])
        den = sum(S["by_unit"][co][u]["n_den"] for u in T["order"])
        if hit / den != float(LEDGER[TB_IDS[co][0]]["원 값"]) or not same(hit, LEDGER[TB_IDS[co][1]]["원 값"]):
            stop(f"cached shares do not pool to the ledger for {co}")
    if not meta["v12"]["equal"]:
        stop("the cached shares were not equal to the v12 shares")
    return T, meta["sel"], S, meta["tables"], meta["run"], meta["v12"]


# ========================================================================================================== Table S18
def half_up_3(n, d):
    """n/d rounded half up to 3 decimals from the exact fraction (no binary rounding)."""
    from fractions import Fraction
    m = math.floor(Fraction(n, d) * 1000 + Fraction(1, 2))
    return f"{m // 1000}.{m % 1000:03d}"


S18_VALUE_KEYS = ("units", "denominator", "pooled", "threshold", "order", "tables", "rounding", "definition", "fields", "what",
                  "grade")


def s18_outputs(S, T, files, cmpv12, thr):
    """drafts/x15/v12/out/F4_S18.json (the basis file of the ledger rows of Table S18) and the pandoc block of the table, with
    «{N:<abbreviation>:<coordinate>}» after every share and «{N:den}» after the denominator, for the ledger IDs."""
    order = T["order"]
    units, dens = {}, set()
    for r, u in enumerate(order, 1):
        ab = ABBR[u]
        ent = dict(row=r, name=u, abbreviation=ab, first_label=T["labs"][u][0], labels=list(T["labs"][u]))
        for co in COORDS:
            b = S["by_unit"][co][u]
            txt = half_up_3(b["n_retest"], b["n_den"])
            if txt != fmt(b["share"], 3) or b["share"] != b["n_retest"] / b["n_den"]:
                stop(f"S18 {ab} {co}: exact rounding {txt} vs float {fmt(b['share'], 3)}")
            ent[co] = dict(n_retest=b["n_retest"], n_den=b["n_den"], share=b["share"], text=txt, below_threshold=b["share"] < thr)
            dens.add(b["n_den"])
        units[ab] = ent
    if len(dens) != 1 or dens != {S["pairs"] * len(S["participants"])}:
        stop(f"S18: denominators {dens}, expected {S['pairs']} pairs × {len(S['participants'])} participants")
    den = dens.pop()
    script = pathlib.Path(__file__).resolve()
    js = dict(what="Supplementary Table S18 — share of units at retest level by ROI unit and coordinate (the values of Fig. 4b)",
              grade="post hoc, descriptive; displayed after the verdict; no test; the verdict is the pooled magnitude test (Fig. 4a)",
              run=RUN_ID, written=now(),
              rules=[dict(file=str(RULES.relative_to(DRAFTS)), sha256=RULES_SHA),
                     dict(file=str(RULES_B.relative_to(DRAFTS)), sha256=RULES_B_SHA, section="§3 (Table S18), §0.4")],
              script=dict(file=str(script.relative_to(DRAFTS)), sha256=hashlib.sha256(script.read_bytes()).hexdigest()),
              definition=("unit = participant × ROI unit × center pair of the fixed magnitude test; at retest level when "
                          "|Δ| <= t0.975(df) · √2 · s, s the c10 retest SD of the ROI unit (c10r1 – c10r3, 6 degrees of freedom); "
                          "share = units at retest level / units with Δ and s, per ROI unit and coordinate; ROI values = means "
                          "over JHU label ∩ FMRIB58 skeleton (>= 2000) voxels (roi.py, table «주»)"),
              threshold=dict(ledger="N-330", value=thr, rule="«*» when share < threshold"),
              rounding="half up to 3 decimals, from the exact fraction n_retest / n_den",
              denominator=dict(value=den, per="ROI unit and coordinate", center_pairs=S["pairs"],
                               participants=len(S["participants"])),
              pooled={co: dict(n_retest=S["pooled"][co]["n_retest"], n_den=S["pooled"][co]["n_den"],
                               share=S["pooled"][co]["share"], ledger=list(S["pooled"][co]["ledger"])) for co in COORDS},
              v12_comparison=cmpv12,
              tables=[dict(file=f, sha256=h, bytes=n) for f, h, n in files],
              order=[ABBR[u] for u in order],
              fields="units.<abbreviation>.<x|y|ODI>.share (raw) · .text (printed) · .n_retest · .n_den; denominator.value",
              units=units)
    S18_JSON.parent.mkdir(parents=True, exist_ok=True)
    kept = False
    if S18_JSON.is_file():
        # v12c: the value file is a ledger basis file; it is re-read (logged) and kept unchanged when its values equal this run's
        raw, sha = reread(S18_JSON, "표 S18 값 파일(장부 근거 파일) — 이번 실행의 값과 같은지 대조 · 같으면 다시 쓰지 않음")
        old = json.loads(raw.decode("utf-8"))
        diff = [k for k in S18_VALUE_KEYS if json.dumps(old.get(k), sort_keys=True) != json.dumps(js[k], sort_keys=True)]
        if diff:
            stop(f"Table S18 value file {S18_JSON.name} (sha256 {sha[:12]}) differs from this run in {diff}; not rewritten")
        kept = True
    else:
        S18_JSON.write_text(json.dumps(js, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
        sha = hashlib.sha256(S18_JSON.read_bytes()).hexdigest()
    rel = os.path.relpath(S18_JSON, REPO)                     # «data/…», the ledger's form for files outside the repository
    L = ["### Table S18 {.unnumbered}", "",
         "Table: **Table S18.** Share of units at retest level by ROI unit and coordinate, as drawn in Fig. 4b (post hoc, "
         "descriptive; displayed after the verdict; no test). A unit is a participant and a center pair of the magnitude test; "
         f"values from skeleton voxels. \\*, below {fmt(thr, 2)} [{THRESHOLD_ID}], the threshold of the pooled fixed test; no "
         "verdict by ROI unit; the verdict is the pooled magnitude test in Fig. 4a. ROI units in atlas label order, with the "
         f"names of the atlas label file. [src: `{rel}` (sha256 {sha[:12]}…) units.*.x · y · ODI.share]{{.src}}", "",
         "| ROI unit | Abbreviation | x | y | ODI |",
         "|----------------------------------------|--------|----------|----------|----------|"]
    for u in order:
        e = units[ABBR[u]]
        cells = [f"{e[co]['text']} {{N:{ABBR[u]}:{co}}}" + ("\\*" if e[co]["below_threshold"] else "") for co in COORDS]
        L.append(f"| {u} | {ABBR[u]} | " + " | ".join(cells) + " |")
    L += ["", f"*Note.* Each share counts the units at retest level among the {den} {{N:den}} units of its ROI unit and "
          "coordinate.", ""]
    md = "\n".join(L)
    if any(p in md for p in OLD_NOT_DRAWN) or [m.group(0) for rx, _ in LC.FORBIDDEN for m in rx.finditer(md)]:
        stop("Table S18: words not to use")
    S18_MD.parent.mkdir(parents=True, exist_ok=True)
    md_kept = False
    if S18_MD.is_file():
        md_raw, _ = reread(S18_MD, "표 S18 블록(v12b 출력) — 이번 실행의 글과 같은지 대조 · 같으면 다시 쓰지 않음")
        md_kept = md_raw.decode("utf-8") == md
    if not md_kept:
        if S18_MD.is_file():
            stop(f"{S18_MD} exists with other content; not overwritten")
        S18_MD.write_text(md, encoding="utf-8")
    n_tok = len(re.findall(r"\{N:[^}]+\}", md))
    if n_tok != 3 * len(order) + 1:
        stop(f"Table S18: {n_tok} ID tokens, expected {3 * len(order) + 1}")
    return dict(json=str(S18_JSON), json_sha256=sha, md=str(S18_MD), md_sha256=hashlib.sha256(md.encode("utf-8")).hexdigest(),
                kept=kept, md_kept=md_kept, n_tokens=n_tok, den=den, below={co: [ABBR[u] for u in order if units[ABBR[u]][co]["below_threshold"]] for co in COORDS})


def reads_md(path):
    L = ["# Fig. 4 (v13) — reads of this run, each logged before it was used", "",
         f"Run {RUN_ID}. Templates and atlas images: `drafts/x14/EXTERNAL_READ_LOG.jsonl` (after guard.require(), entry form of "
         "b2k_read.py); derived tables and reference files: `drafts/x14/REREAD_LOG.jsonl` (through `ledger/reread_log.py`).", "",
         "| log | line | file | sha256 | what |", "|---|---|---|---|---|"]
    for log, line, f, h, what in LOGGED:
        L.append(f"| {log} | {line} | `{f}` | `{h}` | {what} |")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def cvd_sheet():
    """Colour-vision sheet of the figure (direction S53 §5)."""
    s = FC.cvd_sheet(OUT / f"{NAME}.png", OUT / "cvd" / "F4_cvd.png")
    print("cvd sheet:", s.name, hashlib.sha256(s.read_bytes()).hexdigest()[:16])


def main():
    font_line = font_ok()
    rules_unchanged()
    v = read_values()
    cap_text = caption()                                      # caption checks that need no file come before any write
    hits, tagged_words = number_words(cap_text)
    if hits:
        stop(f"caption: number words without a ledger tag: {hits}")
    NUMBER_WORD_TAGS.extend(tagged_words)
    forb = [m.group(0) for rx, _ in LC.FORBIDDEN for m in rx.finditer(cap_text)]
    if forb:
        stop(f"caption uses words not to use: {forb}")
    if INTERNAL.search(re.sub(r"\[N-\d{3}\]", "", cap_text)):
        stop(f"caption names an internal test, arm, hash or session: {INTERNAL.search(cap_text).group(0)}")
    missing_b = [p for p in CAPTION_B_PHRASES if p not in cap_text]
    if missing_b:
        stop(f"caption (b): phrases missing: {missing_b}")
    if any(p in cap_text for p in OLD_NOT_DRAWN):
        stop("caption: the v12 clause on ROI units outside the slices is still there")
    if any(p in cap_text for p in OLD_C):
        stop("caption: a phrase of the removed panel (c) is still there")
    preview = "--preview" in sys.argv
    cached = False
    if preview and CACHE and (pathlib.Path(CACHE) / "meta.json").is_file():
        T, sel, S, files, cache_run, cmpv12 = read_cache()
        cached = True
        print(f"preview from the cache of run {cache_run} (no file under remeasure/ or template read in this run)")
    else:
        basis = G.require()
        T = templates(basis)
        sel = select_slices(T)
        rec, rec_path = recorded_selection()
        if (sel["axial"]["k"], sel["sagittal"]["i"]) != (rec["axial"]["k"], rec["sagittal"]["i"]):
            stop(f"slice choice {sel['axial']['k']} / {sel['sagittal']['i']} differs from {rec_path}")
        values, rep, missing, files = load_tong()
        if missing:
            stop(f"sessions without completion marker: {missing} (the fixed test listed none)")
        S = roi_shares(values, rep, T["order"])               # stops unless the pooled shares equal N-324 · N-326 · N-325
        cmpv12 = compare_v12(S, files)                         # stops unless every share equals the v12 run's
        sel = dict(axial=sel["axial"], sagittal=sel["sagittal"])
        if T["atlas_diff"] is not None:                       # FSL copy vs repository copy, in the two displayed slices
            T["atlas_compare"].update(n_differ_axial_slice=int(T["atlas_diff"][:, :, sel["axial"]["k"]].sum()),
                                      n_differ_sagittal_slice=int(T["atlas_diff"][sel["sagittal"]["i"], :, :].sum()))
        if CACHE:
            write_cache(T, sel, S, files, cmpv12)
    print(f"v12 shares: equal ({cmpv12['n_compared']} unit × coordinate shares and the pooled counts; v12 run {cmpv12['run']}); "
          f"same 36 tables (file, sha256): {cmpv12['same_tables']}")
    ac = T["atlas_compare"]
    print(f"atlas: repository {ac['repo']['sha256'][:12]} · FSL {ac['fsl']['sha256'][:12]} · same sha256: {ac['same_sha256']}"
          + ("" if ac["same_sha256"] else f" · same voxels {ac['same_voxels']} · same affine {ac['same_affine']} · "
             + json.dumps({k: ac.get(k) for k in ("dtype", "n_voxels_differ", "n_labelled", "n_differ_on_skeleton",
                                                  "n_differ_axial_slice", "n_differ_sagittal_slice",
                                                  "fsl_values_are_labels_0_48", "pairs_repo_to_fsl_top", "labels_touched",
                                                  "n_voxels_of_touched_labels")},
                          ensure_ascii=False, default=str)))
    for co in COORDS:                                         # the shares drawn pool to the stored verdict of (a)
        hit = sum(S["by_unit"][co][u]["n_retest"] for u in T["order"])
        den = sum(S["by_unit"][co][u]["n_den"] for u in T["order"])
        if hit != v["a"][co]["n_rt"] or den != v["a"][co]["n_tot"] or hit / den != v["a"][co]["share"]:
            stop(f"(b) {co}: {hit}/{den} does not pool to the drawn share of (a)")
    views = view_arrays(T, sel)
    shares = {co: {u: S["by_unit"][co][u]["share"] for u in T["order"]} for co in COORDS}
    fig, geo = draw(v, views, shares, T["order"])
    check_figure(fig, geo)
    if [u for _, u, _ in geo["strip_rows"]] != list(T["order"]) or len(geo["strip_rows"]) != len(ABBR):
        stop("the strip does not show every ROI unit in atlas label order")
    LAYOUT["not_drawn"] = geo["not_drawn"]
    LAYOUT["strip"] = geo["strip"]
    LAYOUT["maps_x"] = geo["maps_x"]
    nd_below = {co: sum(1 for u in geo["not_drawn"] if shares[co][u] < v["thr"]) for co in COORDS}
    print(f"ROI units with no voxel in either map slice (full labels): {len(geo['not_drawn'])} of {len(T['order'])} "
          f"({', '.join(ABBR[u] for u in geo['not_drawn']) or '—'}; below 0.90 among them: {nd_below}); the strip shows all "
          f"{len(geo['strip_rows'])}")
    if preview:
        fig.savefig(sys.argv[sys.argv.index("--preview") + 1], dpi=300)
        audit(fig, NAME)
        print("below 0.90:", {co: len(geo["below"][co]) for co in COORDS})
        return
    OUT.mkdir(parents=True, exist_ok=True)
    save(fig, NAME)
    (OUT / "F4_caption.md").write_text(cap_text + "\n", encoding="utf-8")
    check_line, warns = ledger_check(OUT / "F4_caption.md", cap_text)
    print("ledger check (caption):", check_line)
    reread_line_numbers()
    s18 = s18_outputs(S, T, files, cmpv12, v["thr"])
    print(f"Table S18: {s18['json']} sha256 {s18['json_sha256']} ({'kept: values equal, not rewritten' if s18['kept'] else 'written'}); "
          f"{s18['md']} sha256 {s18['md_sha256']} ({'kept, same text' if s18['md_kept'] else 'written'}); {s18['n_tokens']} ID "
          f"tokens; denominator {s18['den']}; below 0.90: " + ", ".join(f"{co} {len(s18['below'][co])}" for co in COORDS))
    js = shares_json(S, files, [], sel)
    js["logged"] = LOGGED
    js["v12_comparison"] = cmpv12
    js["atlas_compare"] = T["atlas_compare"]
    (OUT / "F4_shares.json").write_text(json.dumps(js, ensure_ascii=False, indent=1, default=str) + "\n", encoding="utf-8")
    (OUT / "F4_values.md").write_text(values_md(v, S, T, sel, files, cap_text, check_line, warns, cached, cmpv12, s18),
                                      encoding="utf-8")
    reads_md(OUT / "F4_reads.md")
    if pycache_state() != PYCACHE_BEFORE:
        stop("a __pycache__ folder changed during the run")
    cvd_sheet()
    for ext in ("pdf", "png"):
        f = OUT / f"{NAME}.{ext}"
        print(f"sha256 {f.name} {hashlib.sha256(f.read_bytes()).hexdigest()}")
    print("font:", font_line)
    print("logged reads this run:", len(LOGGED), "· below 0.90:", {co: len(geo["below"][co]) for co in COORDS})
    print(f"AUDIT PASS F4_negctrl_roi_prisma: {len(CHECKED)} ledger rows checked, {len(VISIBLE)} visible numbers, pooled "
          + ", ".join(f"{co} {S['pooled'][co]['n_retest']}/{S['pooled'][co]['n_den']}" for co in COORDS)
          + f" = ledger, {cmpv12['n_compared']} shares = v12, strip {len(geo['strip_rows'])} = N-179")


if __name__ == "__main__":
    if "--select-slices" in sys.argv:
        main_select()
    elif "--compute-shares" in sys.argv:
        main_shares()
    else:
        main()
