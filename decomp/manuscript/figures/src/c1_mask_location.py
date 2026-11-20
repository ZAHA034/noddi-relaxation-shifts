"""(c-1) Where the stored noise-floor mask lay in template space, and the projections of panel 7c (option B).

Post hoc, descriptive: defined after the verdict and recorded before computing (drafts/x15/v12/V12B_S43_RULES.md §0 and
§1.1 - §1.7, recorded 2026-11-07 03:14:31 KST, sha256 d56e5563ef17...; RULES_RECORD.jsonl). No new effective L, no test, no
threshold; verdicts unchanged (#34). Nothing is drawn here; make_F7.py draws the projections this script stores.

Steps (each stops the script on any disagreement; the grid checks come before any mask or transform is read and before any
value is computed):
  1. rules file sha256 = the RULES_RECORD.jsonl entry;
  2. records (each logged with the ledger's reread_log.py before it is read): headers.json (sha256 = b2k_read.py HEADERS_SHA),
     read_sessions.jsonl, b2k_result_siemens.json (sha256 = the ledger's rows of the scanner summary);
  3. sessions (§1.1): every session of read_sessions.jsonl (last line per session, as the summary) on OXF1PRI, OXF2PRI,
     OXF3TRI that entered the scanner summary (a stored R); per scanner the number must equal the summary's n_value, no
     session may have been dropped for sigma_q (n_sigma_q_drop 0, otherwise the set cannot be named) and the median of the
     sessions' R must equal the summary's R_med;
  4. grid (§1.2): sens.same_grid true in each session's record; the header of remeasure/onharmony/<session>/dti_FA.nii.gz
     (the input of fa2std_warp) has the raw shape of headers.json and its affine within 1e-3 of the raw affine;
  5. templates (§1.3, §1.4): after guard.require(), FMRIB58_FA_1mm, MNI152_T1_1mm_brain_mask and
     HarvardOxford-sub-maxprob-thr25-1mm (and the label XML, for the names of values 3 and 14) are read and logged in
     drafts/x14/EXTERNAL_READ_LOG.jsonl (entry form of make_F4.py); their shapes and affines must agree (1e-3) and the
     label image must hold the values 3 and 14;
  6. per session: the stored mask «V» (masks/<session>.npz, reread_log.py; «Vf» unused) as a 0/1 uint8 NIfTI with the
     dti_FA header (scratch only); its fa2std_warp.nii.gz (logged in REREAD_LOG.jsonl in the entry form of reread_log.py,
     which refuses .nii.gz) applied with FSL applywarp -r FMRIB58_FA_1mm -w fa2std_warp --interp=nn (cmds.applywarp,
     nearest neighbour); applywarp is given byte-identical scratch copies of the logged warp and reference bytes, so the
     bytes used are the bytes logged; the output must be on the template grid with values 0 and 1 only;
  7. categories (§1.4, priority order, one per template voxel): 1 lateral ventricle (label value 3 or 14), 2 outside the
     brain mask (value 0), 3 inside the brain mask within 5 mm of the nearest outside voxel (distance_transform_edt, 1 mm
     voxels), 4 the rest of the brain mask; counts and shares per session; per scanner median, min, max; pooled medians
     V (category 1), S (category 2 + category 3 per session, then median), O (category 4); case (§1.6): ㄱ if S >= 0.50,
     else ㄴ if V >= 0.50, else ㄷ;
  8. projections (§1.7): per scanner and view (axial along z, coronal along y, sagittal along x; template voxel axes) the
     number of sessions with a warped mask voxel on the line; the brain-mask projection; the projection of category 1.
Writes (C1_OUT, default drafts/x15/v12/out): c1_mask_location.json and c1_projections.npz (its sha256 in the JSON).
Temporary NIfTIs only in C1_TMP (scratch; removed after each session).
Run (CPU only):
  PYTHONDONTWRITEBYTECODE=1 CUDA_VISIBLE_DEVICES= DMRI_REPO=<repo> C1_TMP=<scratch>/tmp nice -n 10 \
    <repo>/.venv/bin/python -I -B c1_mask_location.py
"""
import sys

sys.dont_write_bytecode = True            # no __pycache__ beside the scripts, in the ledger folder or in the repository

import ast  # noqa: E402
import collections  # noqa: E402
import csv  # noqa: E402
import datetime  # noqa: E402
import fcntl  # noqa: E402
import gzip  # noqa: E402
import hashlib  # noqa: E402
import io  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
import pathlib  # noqa: E402
import re  # noqa: E402
import shutil  # noqa: E402
import statistics  # noqa: E402
import subprocess  # noqa: E402
from decimal import Decimal, ROUND_HALF_UP  # noqa: E402

os.environ["CUDA_VISIBLE_DEVICES"] = ""   # CPU only


def stop(msg):
    sys.exit(f"STOP C1: {msg}")


if not os.environ.get("DMRI_REPO") or not os.environ.get("C1_TMP"):
    stop("set DMRI_REPO and C1_TMP (scratch folder for the temporary NIfTIs)")
if os.environ.get("REMEASURE_WORK"):
    stop("REMEASURE_WORK is set; the processing folder must be rconfig's default")
REPO = pathlib.Path(os.environ["DMRI_REPO"]).resolve()
WORK = (REPO / "data").resolve()
DRAFTS = WORK / "stage1_work" / "drafts"
X14 = DRAFTS / "x14"
LEDGER_DIR = DRAFTS / "x15" / "ledger"
V12 = DRAFTS / "x15" / "v12"
OUTDIR = pathlib.Path(os.environ.get("C1_OUT") or V12 / "out").resolve()
TMP = pathlib.Path(os.environ["C1_TMP"]).resolve()
for _p in (REPO, WORK, REPO.parent / "source_data"):
    if str(TMP).startswith(str(_p) + os.sep) or TMP == _p:
        stop(f"C1_TMP {TMP} lies inside {_p}")
TMP.mkdir(parents=True, exist_ok=True)
SCRIPT = pathlib.Path(__file__).resolve()

sys.path.insert(0, str(REPO / "decomp" / "remeasure"))
import rconfig as R  # noqa: E402
import guard as G  # noqa: E402
import numpy as np  # noqa: E402
import nibabel as nib  # noqa: E402
from scipy import ndimage as ndi  # noqa: E402

if R.WORK.resolve() != (WORK / "remeasure").resolve():
    stop(f"rconfig.WORK {R.WORK} is not {WORK / 'remeasure'}")
FSLDIR = R.FSLDIR.resolve()
B2K = R.WORK / "explore" / "b2000"
ONH = R.WORK / "onharmony"
RULES = V12 / "V12B_S43_RULES.md"
RULES_REC = V12 / "RULES_RECORD.jsonl"
RULES_SHA = "d56e5563ef178cb7b786f614f08a9613d3529e697780a24b79cfda69186484b2"
REREAD_LOG = X14 / "REREAD_LOG.jsonl"
REREAD_LOCK = X14 / "reread.lock"          # the lock reread_log.py takes
EXT_LOG = X14 / "EXTERNAL_READ_LOG.jsonl"
SIE = ("OXF1PRI", "OXF2PRI", "OXF3TRI")
RUN_ID = "C1-" + datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
WHAT = "v12b 표시 · (c-1) — post hoc · 판정 불변 #34"
TOL = 1e-3                                 # §1.2 · §1.3: affines equal within 1e-3
BAND_MM = 5.0                              # §1.4 category 3, fixed before computing
LV_VALUES = (3, 14)                        # §1.4 category 1: XML index 2 -> 3 (Left), 13 -> 14 (Right Lateral Ventricle)
LV_NAMES = {2: "Left Lateral Ventricle", 13: "Right Lateral Ventricle"}
T_FA = FSLDIR / "data" / "standard" / "FMRIB58_FA_1mm.nii.gz"
T_BM = FSLDIR / "data" / "standard" / "MNI152_T1_1mm_brain_mask.nii.gz"
T_HO = FSLDIR / "data" / "atlases" / "HarvardOxford" / "HarvardOxford-sub-maxprob-thr25-1mm.nii.gz"
T_HO_XML = FSLDIR / "data" / "atlases" / "HarvardOxford-Subcortical.xml"
VIEWS = ("axial", "coronal", "sagittal")   # projection along z, y, x (§1.7)
VIEW_AXIS = {"axial": 2, "coronal": 1, "sagittal": 0}     # world axis (x 0, y 1, z 2) the view projects along
READS = []


def now():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def sha_bytes(b):
    return hashlib.sha256(b).hexdigest()


def rel(p):
    return os.path.relpath(pathlib.Path(p).resolve(), REPO)


def half_up(x, k=2):
    return str(Decimal(repr(float(x))).quantize(Decimal(1).scaleb(-k), rounding=ROUND_HALF_UP))


# ----------------------------------------------------------------------------- logged reads
def _last_line(log, pred):
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


def reread_tool(path, what, want_sha=None):
    """A derived file: logged by the ledger's reread_log.py (hash, append under its lock), then read once."""
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
    ln, e = _last_line(REREAD_LOG, lambda x: x.get("다시읽음") == str(path) and x.get("sha256") == h and x.get("무엇") == w)
    READS.append(dict(kind="derived file (reread_log.py)", file=rel(path), sha256=h, bytes=len(raw),
                      log=f"drafts/x14/{REREAD_LOG.name}", line=ln, time=e.get("시각"), what=what))
    return raw


def reread_direct(path, what):
    """A derived NIfTI (reread_log.py refuses .nii.gz): bytes read once, then an entry in reread_log.py's own form appended
    to REREAD_LOG.jsonl under reread_log.py's lock, BEFORE any value of the bytes is used."""
    path = pathlib.Path(path).resolve()
    if "source_data" in path.parts:
        stop(f"{path}: not a derived file")
    raw = path.read_bytes()
    h = sha_bytes(raw)
    w = f"{WHAT} · {what} · 실행 {RUN_ID}"
    rec = dict(시각=now(), 다시읽음=str(path), sha256=h, 바이트=len(raw), 무엇=w, 스크립트=str(SCRIPT))
    with open(REREAD_LOCK, "a") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        try:
            cur = REREAD_LOG.read_bytes()
            if cur and not cur.endswith(b"\n"):
                stop(f"{REREAD_LOG.name} does not end with a newline")
            ln = cur.count(b"\n") + 1
            with REREAD_LOG.open("a", encoding="utf-8") as fh:
                fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
                fh.flush()
                os.fsync(fh.fileno())
        finally:
            fcntl.flock(lk, fcntl.LOCK_UN)
    READS.append(dict(kind="derived NIfTI (reread_log.py entry form; the helper refuses .nii.gz)", file=rel(path), sha256=h,
                      bytes=len(raw), log=f"drafts/x14/{REREAD_LOG.name}", line=ln, time=rec["시각"], what=what))
    return raw


def ext_read(path, what, basis):
    """FSL template or atlas file: after guard.require(), bytes read once, logged in EXTERNAL_READ_LOG.jsonl (make_F4.py
    entry form, lock on the log itself) BEFORE any value is used."""
    raw = pathlib.Path(path).read_bytes()
    h = sha_bytes(raw)
    entry = dict(시각=now(), utc=utc(), 세션="템플릿(참가자 자료 아님)", 파일=rel(path), sha256=h, 바이트=len(raw),
                 무엇=f"{WHAT} · {what} · 실행 {RUN_ID}", 근거=basis, 스크립트=str(SCRIPT))
    data = (json.dumps(entry, ensure_ascii=False, default=str) + "\n").encode("utf-8")
    fd = os.open(EXT_LOG, os.O_RDWR | os.O_APPEND)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX)
        try:
            os.lseek(fd, 0, os.SEEK_SET)
            n0, last = 0, b""
            while True:
                blk = os.read(fd, 1 << 20)
                if not blk:
                    break
                n0 += blk.count(b"\n")
                last = blk
            if last and not last.endswith(b"\n"):
                stop(f"{EXT_LOG.name} does not end with a newline")
            if os.write(fd, data) != len(data):
                stop(f"short write to {EXT_LOG}")
        finally:
            fcntl.flock(fd, fcntl.LOCK_UN)
    finally:
        os.close(fd)
    READS.append(dict(kind="FSL template / atlas (EXTERNAL_READ_LOG, after guard.require())", file=rel(path), sha256=h,
                      bytes=len(raw), log=f"drafts/x14/{EXT_LOG.name}", line=n0 + 1, time=entry["시각"], what=what))
    return raw


def nii(raw):
    return nib.Nifti1Image.from_bytes(gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw)


# ----------------------------------------------------------------------------- steps
def rules_check():
    if sha_bytes(RULES.read_bytes()) != RULES_SHA:
        stop("V12B_S43_RULES.md differs from the recorded sha256")
    recs = [json.loads(s) for s in RULES_REC.read_text(encoding="utf-8").splitlines() if s.strip()]
    hit = [r for r in recs if r.get("sha256") == RULES_SHA and str(r.get("파일", "")).endswith("V12B_S43_RULES.md")]
    if len(hit) != 1:
        stop("RULES_RECORD.jsonl does not hold the rules file once")
    return hit[0]["시각"]


def py_literal(path, name):
    tree = ast.parse(path.read_text(encoding="utf-8"))
    hits = [n for n in tree.body if isinstance(n, ast.Assign) and len(n.targets) == 1
            and isinstance(n.targets[0], ast.Name) and n.targets[0].id == name]
    if len(hits) != 1:
        stop(f"{path.name}: {len(hits)} top-level assignments of {name}")
    return ast.literal_eval(hits[0].value)


def summary_sha():
    """The sha256 the ledger records for the scanner summary (rows N-443 ... N-457 cite it)."""
    led = {r["ID"]: r for r in csv.DictReader((LEDGER_DIR / "numbers.csv").open(encoding="utf-8"))}
    shas = {led[i]["근거 파일 sha256"] for i in ("N-443", "N-446", "N-449", "N-452", "N-454", "N-456")}
    files = {led[i]["근거 파일"] for i in ("N-443", "N-446", "N-449", "N-452", "N-454", "N-456")}
    if len(shas) != 1 or files != {"data/remeasure/explore/b2000/b2k_result_siemens.json"}:
        stop(f"ledger rows of the scanner summary: files {files}, sha256 {shas}")
    return shas.pop()


def sessions_set(H_, recs, S_):
    """§1.1: the first sessions on the three Siemens scanners that entered the scanner summary."""
    last, n_lines = {}, collections.Counter()
    for r in recs:
        last[r["session"]] = r                    # as the summary: the last line of a session
        n_lines[r["session"]] += 1
    dup = sorted(k for k, n in n_lines.items() if n > 1)
    out, chk = {}, {}
    for dev in SIE:
        names = sorted(k for k, r in last.items() if r.get("device") == dev)
        for k in names:
            h = H_["sessions"].get(k)
            if h is None or h["device"] != dev or last[k]["unit"] != dev:
                stop(f"{k}: not a {dev} session of headers.json")
        inc = [k for k in names if last[k].get("R") is not None]
        e = S_["devices"][dev]
        n_hdr = sum(1 for h in H_["sessions"].values() if h["device"] == dev)
        med = statistics.median([last[k]["R"] for k in inc]) if inc else None
        c = dict(n_headers=n_hdr, n_read=len(names), n_with_R=len(inc), summary_n_expected=e.get("n_expected"),
                 summary_n_records=e.get("n_records"), summary_n_value=e.get("n_value"),
                 summary_n_sigma_q_drop=e.get("n_sigma_q_drop"), summary_R_med=e.get("R_med"), median_R_of_set=med)
        chk[dev] = c
        if e.get("n_sigma_q_drop", 0) != 0:
            stop(f"{dev}: the summary dropped {e.get('n_sigma_q_drop')} sessions for sigma_q; which ones is not recorded "
                 f"here, so the session set cannot be named ({c})")
        if len(names) != e.get("n_records") or len(inc) != e.get("n_value"):
            stop(f"{dev}: session counts differ from the scanner summary: {c}")
        if med is None or abs(med - e["R_med"]) > 1e-12 * max(1.0, abs(e["R_med"])):
            stop(f"{dev}: the median R of the session set differs from the summary's R_med: {c}")
        out[dev] = inc
    return out, last, chk, dup


def grid_checks(H_, last, sessions):
    """§1.2: sens.same_grid true, and dti_FA (the input of fa2std_warp) on the raw grid of headers.json. Stops at the first
    disagreement, before any mask or transform is read."""
    res, fa_hdr = {}, {}
    for dev in SIE:
        for k in sessions[dev]:
            sens = last[k].get("sens") or {}
            if sens.get("same_grid") is not True:
                stop(f"GRID §1.2: {k} ({dev}) sens.same_grid = {sens.get('same_grid')!r} (sens keys {sorted(sens)}; "
                     f"error {sens.get('error')!r})")
    for dev in SIE:
        for k in sessions[dev]:
            h = H_["sessions"][k]
            raw = reread_direct(ONH / k / "dti_FA.nii.gz", "처리 격자 — dti_FA 헤더(모양 · 아핀; 값은 쓰지 않음 · fa2std_warp 의 입력)")
            img = nii(raw)
            shape_raw = [int(x) for x in h["shape"][:3]]
            aff_raw = np.asarray(h["affine"], float)
            d = float(np.abs(img.affine - aff_raw).max())
            res[k] = dict(device=dev, same_grid=True, raw_shape=shape_raw, fa_shape=[int(x) for x in img.shape],
                          max_abs_affine_diff=d, fa_sha256=sha_bytes(raw))
            if len(img.shape) != 3 or [int(x) for x in img.shape] != shape_raw:
                stop(f"GRID §1.2: {k} ({dev}) dti_FA shape {list(img.shape)} differs from the raw shape {shape_raw}")
            if not d <= TOL:
                stop(f"GRID §1.2: {k} ({dev}) dti_FA affine differs from the raw affine by {d} (> {TOL}); "
                     f"dti_FA {img.affine.tolist()} raw {aff_raw.tolist()}")
            fa_hdr[k] = img.header.copy()
    return res, fa_hdr


def templates(basis):
    """§1.3 · §1.4: the three template images on one grid; the label image holds 3 and 14."""
    xml = ext_read(T_HO_XML, "Harvard-Oxford 피질하 라벨 이름(XML) — 값 3 · 14 = 측뇌실 좌 · 우 대조", basis).decode("utf-8")
    names = {int(m.group(1)): m.group(2).strip() for m in re.finditer(r'<label index="(\d+)"[^>]*>([^<]+)</label>', xml)}
    for idx, nm in LV_NAMES.items():
        if names.get(idx) != nm:
            stop(f"HO XML: index {idx} is {names.get(idx)!r}, not {nm!r}")
    raws = {"FA": ext_read(T_FA, "FMRIB58_FA_1mm — applywarp 기준 영상(-r) · 템플릿 격자 대조", basis),
            "BM": ext_read(T_BM, "MNI152_T1_1mm_brain_mask — 범주 2 · 3(뇌 마스크 바깥 · 5 mm 띠) · 회색 실루엣", basis),
            "HO": ext_read(T_HO, "HarvardOxford-sub-maxprob-thr25-1mm — 범주 1(측뇌실 값 3 · 14) · 점선 윤곽", basis)}
    imgs = {k: nii(v) for k, v in raws.items()}
    shp = {k: [int(x) for x in im.shape] for k, im in imgs.items()}
    aff = {k: im.affine for k, im in imgs.items()}
    diffs = {f"{a}-{b}": float(np.abs(aff[a] - aff[b]).max()) for a, b in (("FA", "BM"), ("FA", "HO"), ("BM", "HO"))}
    info = dict(shapes=shp, affines={k: v.tolist() for k, v in aff.items()}, max_abs_affine_diff=diffs,
                zooms={k: [float(z) for z in im.header.get_zooms()] for k, im in imgs.items()},
                sform_code={k: int(im.header["sform_code"]) for k, im in imgs.items()},
                qform_code={k: int(im.header["qform_code"]) for k, im in imgs.items()},
                sha256={k: sha_bytes(v) for k, v in raws.items()}, ho_xml_names={str(i): names[i] for i in LV_NAMES})
    if len({tuple(s) for s in shp.values()}) != 1:
        stop(f"GRID §1.3: template shapes differ: {shp}")
    if any(not d <= TOL for d in diffs.values()):
        stop(f"GRID §1.3: template affines differ: {diffs}; {info['affines']}")
    HO = np.asanyarray(imgs["HO"].dataobj)
    vals = set(np.unique(HO).tolist())
    info["ho_values"] = sorted(int(v) for v in vals)
    if not all(v in vals for v in LV_VALUES):
        stop(f"§1.4: the label image lacks {[v for v in LV_VALUES if v not in vals]} (values {sorted(vals)})")
    BMd = np.asanyarray(imgs["BM"].dataobj)
    info["brain_mask_values"] = sorted(float(v) for v in np.unique(BMd))
    zooms = [float(z) for z in imgs["BM"].header.get_zooms()[:3]]
    if zooms != [1.0, 1.0, 1.0]:
        stop(f"§1.4: template voxels are {zooms} mm, not 1 mm")
    return imgs, raws, info, HO, BMd


def categories(HO, BMd):
    """§1.4, priority order: 1 lateral ventricle, 2 outside the brain mask, 3 inside within 5 mm of outside, 4 the rest."""
    LV = np.isin(HO, LV_VALUES)
    BM = BMd != 0
    dist = ndi.distance_transform_edt(BM, sampling=(1.0, 1.0, 1.0))
    cat = np.full(BM.shape, 4, np.uint8)
    cat[BM & (dist <= BAND_MM)] = 3
    cat[~BM] = 2
    cat[LV] = 1                                   # highest priority last
    border = int(BM[0].sum() + BM[-1].sum() + BM[:, 0].sum() + BM[:, -1].sum() + BM[:, :, 0].sum() + BM[:, :, -1].sum())
    stats = dict(n_template={str(c): int((cat == c).sum()) for c in (1, 2, 3, 4)}, n_lv_label=int(LV.sum()),
                 n_lv_outside_brain_mask=int((LV & ~BM).sum()), n_brain_mask=int(BM.sum()), brain_mask_voxels_on_grid_faces=border)
    return cat, LV, BM, stats


def world_axes(aff):
    """Voxel axis closest to world x, y, z (asserted one-to-one)."""
    A = np.asarray(aff, float)[:3, :3]
    dom = [int(np.argmax(np.abs(A[:, j]))) for j in range(3)]           # world axis of voxel axis j
    if sorted(dom) != [0, 1, 2]:
        stop(f"template voxel axes are not aligned with x, y, z: {A.tolist()}")
    return {w: dom.index(w) for w in range(3)}                        # world axis -> voxel axis


def fsl_env():
    return dict(os.environ, CUDA_VISIBLE_DEVICES="", OMP_NUM_THREADS="1", FSLDIR=str(FSLDIR), FSLOUTPUTTYPE="NIFTI_GZ",
                PATH=f"{FSLDIR / 'bin'}:{os.environ.get('PATH', '')}")


def warp_session(k, dev, rec, h, fa_hdr, ref_copy, t_shape, t_aff):
    raw = reread_tool(B2K / "masks" / f"{k}.npz", "저장 소음 바닥 마스크 V(원 dwi 격자 복셀 번호) — 템플릿으로 옮겨 범주 셈 · 투영")
    z = np.load(io.BytesIO(raw), allow_pickle=False)
    if "V" not in z.files:
        stop(f"{k}: stored mask without «V»")
    V = np.asarray(z["V"], np.int64)
    shape = [int(x) for x in h["shape"][:3]]
    if V.ndim != 2 or V.shape[1] != 3 or len(V) != rec["n_V"]:
        stop(f"{k}: stored mask {V.shape}, not the stored n_V {rec['n_V']}")
    if (V < 0).any() or (V >= np.asarray(shape)).any() or len(np.unique(V, axis=0)) != len(V):
        stop(f"{k}: stored mask outside the grid or with a voxel twice")
    M = np.zeros(shape, np.uint8)
    M[tuple(V.T)] = 1
    hdr = fa_hdr.copy()
    hdr.set_data_dtype(np.uint8)
    f_in, f_w, f_out = TMP / f"{k}_maskV.nii.gz", TMP / f"{k}_fa2std_warp.nii.gz", TMP / f"{k}_maskV_std.nii.gz"
    nib.save(nib.Nifti1Image(M, None, header=hdr), str(f_in))
    chk = nib.load(str(f_in))
    if np.abs(chk.affine - fa_hdr.get_best_affine()).max() > 0 or int(chk.header["sform_code"]) != int(fa_hdr["sform_code"]) \
            or int(chk.header["qform_code"]) != int(fa_hdr["qform_code"]) or int(np.asanyarray(chk.dataobj).sum()) != len(V):
        stop(f"{k}: the 0/1 mask NIfTI does not carry the dti_FA header or the mask")
    wraw = reread_direct(ONH / k / "fa2std_warp.nii.gz", "회기 변환 fa2std_warp(기존 · ROI 추출과 같은 변환) — applywarp -w "
                                                          "(로그한 바이트 그대로 scratch 사본으로 넘김)")
    f_w.write_bytes(wraw)
    if sha_bytes(f_w.read_bytes()) != sha_bytes(wraw):
        stop(f"{k}: warp copy differs from the logged bytes")
    cmd = ["nice", "-n", "10", str(FSLDIR / "bin" / "applywarp"), "-i", str(f_in), "-r", str(ref_copy), "-w", str(f_w),
           "-o", str(f_out), "--interp=nn"]
    r = subprocess.run(cmd, capture_output=True, text=True, env=fsl_env())
    if r.returncode != 0:
        stop(f"{k}: applywarp failed: {r.stdout.strip()} {r.stderr.strip()}")
    out = nib.load(str(f_out))
    d = np.asanyarray(out.dataobj)
    u = set(np.unique(d).tolist())
    if [int(x) for x in out.shape] != list(t_shape) or float(np.abs(out.affine - t_aff).max()) > TOL:
        stop(f"{k}: applywarp output not on the template grid ({out.shape}, {out.affine.tolist()})")
    if not u <= {0, 1}:
        stop(f"{k}: applywarp output values {sorted(u)[:10]} are not 0/1")
    W = np.argwhere(d == 1).astype(np.int16)
    meta = dict(applywarp=" ".join(["applywarp", "-i", "<mask V as 0/1 NIfTI, dti_FA header>", "-r",
                                    "<copy of FMRIB58_FA_1mm.nii.gz>", "-w", f"<copy of onharmony/{k}/fa2std_warp.nii.gz>",
                                    "-o", "<scratch>", "--interp=nn"]),
                out_dtype=str(out.header.get_data_dtype()), warp_sha256=sha_bytes(wraw), mask_sha256=sha_bytes(raw),
                stderr=r.stderr.strip()[:300])
    for f in (f_in, f_w, f_out):
        f.unlink()
    return V, W, meta


def main():
    t0 = now()
    rules_t = rules_check()
    out_json, out_npz = OUTDIR / "c1_mask_location.json", OUTDIR / "c1_projections.npz"
    if (out_json.exists() or out_npz.exists()) and os.environ.get("C1_OVERWRITE") != "1":
        stop(f"{out_json} or {out_npz} exists (set C1_OVERWRITE=1 to replace)")
    headers_sha = py_literal(DRAFTS / "x15" / "b2000" / "b2k_read.py", "HEADERS_SHA")
    H_ = json.loads(reread_tool(B2K / "headers.json", "단계 1 기록 — 회기 · 원 dwi 격자(모양 · 아핀)", want_sha=headers_sha))
    recs = [json.loads(s) for s in reread_tool(B2K / "read_sessions.jsonl", "단계 3 회기 기록 — 회기 · n_V · R · sens.same_grid")
            .decode("utf-8").splitlines() if s.strip()]
    S_ = json.loads(reread_tool(B2K / "b2k_result_siemens.json", "판정 뒤 분석 No. 9 스캐너 요약 — 회기 수 · R_med(회기 집합 대조)",
                                want_sha=summary_sha()))
    sessions, last, set_chk, dup = sessions_set(H_, recs, S_)
    grid, fa_hdr = grid_checks(H_, last, sessions)
    basis = G.require()
    imgs, traws, tinfo, HO, BMd = templates(basis)
    t_shape = tinfo["shapes"]["FA"]
    t_aff = imgs["FA"].affine
    cat, LV, BM, cstats = categories(HO, BMd)
    wax = world_axes(t_aff)
    ref_copy = TMP / "FMRIB58_FA_1mm.nii.gz"
    ref_copy.write_bytes(traws["FA"])
    if sha_bytes(ref_copy.read_bytes()) != tinfo["sha256"]["FA"]:
        stop("reference copy differs from the logged bytes")
    time_grid = now()
    rows, warped, metas = [], {}, {}
    for dev in SIE:
        for k in sessions[dev]:
            V, W, meta = warp_session(k, dev, last[k], H_["sessions"][k], fa_hdr[k], ref_copy, t_shape, t_aff)
            n = np.bincount(cat[tuple(W.T.astype(np.int64))], minlength=5)[1:5] if len(W) else np.zeros(4, int)
            tot = int(n.sum())
            if tot == 0:
                stop(f"{k}: no template voxel of value 1 after applywarp")
            sh = [float(x) / tot for x in n]
            rows.append(dict(session=k, scanner=dev, n_mask_voxels_native=int(len(V)), n_template_voxels=tot,
                             n_by_category={str(c + 1): int(n[c]) for c in range(4)},
                             share_by_category={str(c + 1): sh[c] for c in range(4)}, S_share=sh[1] + sh[2], **meta))
            warped[k] = W
            print(f"{k} {dev} native {len(V)} template {tot} shares " + " ".join(f"{x:.3f}" for x in sh), flush=True)
    ref_copy.unlink()
    # ---- summaries
    def stat(vals):
        return dict(median=float(statistics.median(vals)), min=float(min(vals)), max=float(max(vals)), n=len(vals))
    keys = ("1", "2", "3", "4")
    per_scanner = {}
    for dev in SIE:
        rs = [r for r in rows if r["scanner"] == dev]
        per_scanner[dev] = {f"cat{c}": stat([r["share_by_category"][c] for r in rs]) for c in keys}
        per_scanner[dev]["S"] = stat([r["S_share"] for r in rs])
    pooled = {f"cat{c}": stat([r["share_by_category"][c] for r in rows]) for c in keys}
    pooled["S"] = stat([r["S_share"] for r in rows])
    Vv, Sv, Ov = pooled["cat1"]["median"], pooled["S"]["median"], pooled["cat4"]["median"]
    case = "ㄱ" if Sv >= 0.50 else ("ㄴ" if Vv >= 0.50 else "ㄷ")
    text = {"per_scanner": {dev: {k_: {s: half_up(v[s]) for s in ("median", "min", "max")} for k_, v in d.items()}
                            for dev, d in per_scanner.items()},
            "pooled": {k_: {s: half_up(v[s]) for s in ("median", "min", "max")} for k_, v in pooled.items()},
            "V": half_up(Vv), "S": half_up(Sv), "O": half_up(Ov)}
    # ---- projections (§1.7)
    proj = {"template_affine": np.asarray(t_aff, float), "template_shape": np.asarray(t_shape, np.int64)}
    for view in VIEWS:
        ax = wax[VIEW_AXIS[view]]
        proj[f"brain__{view}"] = BM.any(axis=ax)
        proj[f"lv__{view}"] = LV.any(axis=ax)
        for dev in SIE:
            cnt = None
            for k in sessions[dev]:
                M = np.zeros(t_shape, bool)
                M[tuple(warped[k].T.astype(np.int64))] = True
                p = M.any(axis=ax).astype(np.int16)
                cnt = p if cnt is None else cnt + p
            proj[f"count__{dev}__{view}"] = cnt
    for k, W in warped.items():
        proj[f"warped__{k}"] = W
    OUTDIR.mkdir(parents=True, exist_ok=True)
    buf = io.BytesIO()
    np.savez_compressed(buf, **proj)
    npz_raw = buf.getvalue()
    out_npz.write_bytes(npz_raw)
    res = dict(
        what="(c-1) location of the stored noise-floor mask in template space; post hoc, descriptive; defined after the "
             "verdict and recorded before computing; no new effective L; verdicts unchanged (#34)",
        rules=dict(file="drafts/x15/v12/V12B_S43_RULES.md", sha256=RULES_SHA, recorded=rules_t, sections="§0, §1.1-§1.7"),
        run=dict(id=RUN_ID, started=t0, grid_checks_passed=time_grid, finished=now(), script=rel(SCRIPT),
                 script_sha256=sha_bytes(SCRIPT.read_bytes()), fsldir=rel(FSLDIR), guard_basis=basis,
                 python=sys.version.split()[0], numpy=np.__version__, nibabel=nib.__version__,
                 scipy=__import__("scipy").__version__),
        definitions=dict(
            sessions="read_sessions.jsonl (last line per session) on OXF1PRI, OXF2PRI, OXF3TRI with a stored R = the sessions "
                     "of b2k_result_siemens.json (n_value; n_sigma_q_drop 0; median R = R_med)",
            mask="masks/<session>.npz «V» (raw dwi grid), written as 0/1 uint8 NIfTI with the dti_FA header",
            transform="applywarp -r FMRIB58_FA_1mm -w <session>/fa2std_warp --interp=nn (cmds.applywarp with nearest neighbour)",
            categories={"1": "lateral ventricle: HarvardOxford-sub-maxprob-thr25-1mm value 3 or 14",
                        "2": "outside MNI152_T1_1mm_brain_mask (value 0), not category 1",
                        "3": "inside the brain mask with distance_transform_edt to the nearest outside voxel <= 5 mm "
                             "(1 mm voxels), not category 1",
                        "4": "the rest of the brain mask"},
            shares="per session: template voxels of value 1 by category / all template voxels of value 1",
            V="median over all sessions of the category-1 share", S="median over all sessions of (category 2 + category 3 share)",
            O="median over all sessions of the category-4 share",
            case_rule="ㄱ if S >= 0.50; else ㄴ if V >= 0.50; else ㄷ (§1.6, fixed before computing)",
            rounding="text fields: half up to two decimals (Decimal of repr)"),
        session_set=dict(checks=set_chk, duplicate_lines_in_read_sessions=dup,
                         n_sessions={dev: len(sessions[dev]) for dev in SIE}),
        grid=dict(sessions=grid, templates=tinfo, tolerance=TOL),
        template_categories=cstats,
        sessions=rows,
        per_scanner=per_scanner, pooled=pooled, V=Vv, S=Sv, O=Ov, case=case, text=text,
        projections=dict(file=out_npz.name, sha256=sha_bytes(npz_raw), views={v: f"projection along world axis "
                                                                                  f"{'xyz'[VIEW_AXIS[v]]} = template voxel "
                                                                                  f"axis {wax[VIEW_AXIS[v]]}" for v in VIEWS},
                         keys=sorted(proj), pixel="count of the scanner's sessions with a warped mask voxel on the line "
                                                  "(share = count / n_sessions)"),
        reads=READS)
    out_json.write_text(json.dumps(res, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    print(f"written {out_json} (sha256 {sha_bytes(out_json.read_bytes())}) and {out_npz} (sha256 {sha_bytes(npz_raw)})")
    print("sessions:", {dev: len(sessions[dev]) for dev in SIE}, "· case", case, "· V", Vv, "· S", Sv, "· O", Ov)
    for dev in SIE:
        print(dev, " · ".join(f"{k_} {v['median']:.3f} [{v['min']:.3f}–{v['max']:.3f}]" for k_, v in per_scanner[dev].items()))
    print("pooled", " · ".join(f"{k_} {v['median']:.3f} [{v['min']:.3f}–{v['max']:.3f}]" for k_, v in pooled.items()))
    print("reads logged:", len(READS))


if __name__ == "__main__":
    main()
