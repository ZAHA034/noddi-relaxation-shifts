"""Figure 1 (v13; direction S45 §3.1, author decision 2026-11-08: Figure 1 keeps panel (a) only; v13d, direction S50 A2: no title
line and no role chip above the maps, column heads and row names stay): one brain at three echo times. Displayed after the
verdict, no test (Section 3.5). Selection: V12_DISPLAY_RULES.md §1 (a), recorded 2026-11-06 20:17:09 KST before any map value
was read (the choices were written to fig1/selection.md before the first read).

Multi-TE data, the participant with the smallest ID: NDI, FISO and ODI maps of the first series TE62, TE92 and TE132; one row
per index, one column per TE, one colour bar per row. One axial slice, MNI z = +20 mm. Gray scale, linear; per index one
window, the 1st and 99th percentiles of the pooled voxels of all displayed maps of that index (now the three maps of the row)
inside the FSL MNI152 brain mask on the slice; values beyond are shown at the ends; voxels drawn with interpolation "nearest".
The windows that also pooled the v12 panel (c) maps are written to F1_values.md for reference only (not drawn).
Width 140 mm (one panel); tiles 32 mm wide; neurological orientation.

Maps: only native-space maps are stored (run.py to_std_and_roi deletes the template-space maps after ROI extraction), so each
map is taken to the template exactly as the ROI extraction did: applywarp -i <map> -r FMRIB58_FA_1mm -w <the same folder's
fa2std_warp> --interp=trilinear (cmds.applywarp). The output is a display resample, kept in the cache folder.
Reads: guard.require() first (decomp/remeasure/guard.py); every map and warp file is logged with time, sha256 and bytes to
drafts/x14/ONH_NEW_READ_LOG.jsonl (ON-Harmony) or drafts/x14/EXTERNAL_READ_LOG.jsonl (Tong, multi-TE), under the lock that
b2k_read.py uses, before applywarp opens it. The read part below (selection, maps, the v12 (b) ROI tables, cache) is the v12
read part unchanged, so that the v12 cache (F1_CACHE) is used as it is: no new map read unless --reread is given. Panels (b)
and (c) of v12/v12b are not drawn any more; their cached values are checked by the cache guard and not used.
Guards: rules file sha256; the selection re-derived from names and metadata equals the recorded one; the slice index from
the template affine; windows recomputed from the three displayed maps of each index; every tile shows its cached slice with
its row's window and one scale on both image axes; every printed number is a colour-bar window value (display setting) or
a data name; texts >= 7 pt Liberation Sans, inside the figure, no overlaps; every PDF font Liberation Sans and embedded.
Writes into FIG_OUT only: F1_maps.pdf/.png, F1_values.md, F1_reads.json; cache in F1_CACHE (default FIG_OUT/F1_cache).
Run (CPU only):
  PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=<scratch> CUDA_VISIBLE_DEVICES= DMRI_REPO=<repo> FIG_OUT=<out> F1_CACHE=<map cache> \
    nice -n 10 <repo>/.venv/bin/python -I -B make_F1.py [--reread]
(a cache from an earlier run of this script is used unless --reread is given; a first run always reads and logs).
"""
import sys

sys.dont_write_bytecode = True            # no __pycache__ in the repository or the ledger folder

import datetime  # noqa: E402
import fcntl  # noqa: E402
import hashlib  # noqa: E402
import importlib.util  # noqa: E402
import itertools  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import os  # noqa: E402
import pathlib  # noqa: E402
import re  # noqa: E402
import subprocess  # noqa: E402
import ast  # noqa: E402

os.environ["CUDA_VISIBLE_DEVICES"] = ""   # CPU only

import numpy as np  # noqa: E402

REPO = pathlib.Path(os.environ.get("DMRI_REPO") or pathlib.Path(__file__).resolve().parents[4])
WORK = pathlib.Path(os.environ.get("DMRI_WORK") or REPO / "data")
DRAFTS = WORK / "stage1_work" / "drafts"
X14 = DRAFTS / "x14"
LEDGER = DRAFTS / "x15" / "ledger"
SOURCES = LEDGER / "sources"
REMEASURE = WORK / "remeasure"
OUT = pathlib.Path(os.environ.get("FIG_OUT") or DRAFTS / "figs_v12")
OUT.mkdir(parents=True, exist_ok=True)
CACHE = pathlib.Path(os.environ.get("F1_CACHE") or OUT / "F1_cache")
NAME = "F1_maps"
RULES = DRAFTS / "x15" / "v12" / "V12_DISPLAY_RULES.md"
RULES_SHA = "a4197909b5bb54f2056fdc48f1f2087cb02c73723521d5e0639cb3c07c1cc0f8"   # RULES_RECORD.jsonl 2026-11-06T20:17:09+0900
LOCK = REMEASURE / "explore" / "b2000" / "read_log.lock"                         # the lock b2k_read.py uses for these logs
WHAT = "v12 표시 — post hoc display · displayed after the verdict · no test · 판정 불변 #34 · Fig. 1({})"

sys.path.insert(0, str(REPO / "decomp" / "remeasure"))
import rconfig as R  # noqa: E402
import guard as G  # noqa: E402
import judge_tong as JT  # noqa: E402
import judge_wong as JW  # noqa: E402

spec = importlib.util.spec_from_file_location("figcolors", pathlib.Path(__file__).with_name("figcolors.py"))
FC = importlib.util.module_from_spec(spec)
spec.loader.exec_module(FC)


def stop(msg):
    print(f"STOP F1: {msg}", flush=True)
    raise SystemExit(1)


def sha_file(p):
    h, n = hashlib.sha256(), 0
    with open(p, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
            n += len(blk)
    return h.hexdigest(), n


def now():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


# ------------------------------------------------------------------ ledger
def ledger_rows():
    import csv
    with open(LEDGER / "numbers.csv", encoding="utf-8") as fh:
        return {r["ID"]: r for r in csv.DictReader(fh)}


ROWS = ledger_rows()


def frozen_copy(rel, row_id):
    """The ledger's frozen copy of a source file, checked against the sha256 the ledger row records."""
    r = ROWS[row_id]
    if r["근거 파일"] != rel:
        stop(f"{row_id}: source file is {r['근거 파일']}, not {rel}")
    want = r["근거 파일 sha256"]
    p = SOURCES / rel.replace("../", "__up__/") / want[:12]
    if not p.exists():
        stop(f"ledger copy missing: {p}")
    got, _ = sha_file(p)
    if got != want:
        stop(f"ledger copy {p} sha256 {got} != ledger {want}")
    return p


# ------------------------------------------------------------------ selection (V12_DISPLAY_RULES §1), recorded before reading
SELECTED = dict(
    a=dict(participant="sub-01", series=("TE62", "TE92", "TE132"),
           folders={"TE62": "wong/sub-01", "TE92": "wong/sub-01/own_TE92", "TE132": "wong/sub-01/own_TE132"}),
    b=dict(participant="sub-01", series=tuple(f"TE{t}" for t in (62, 72, 82, 92, 102, 112, 122, 132))),
    tong=dict(participant="sub-1", centers=("c01", "c02"),
              folders={"c01": "tong/sub-1_ses-c01r1", "c02": "tong/sub-1_ses-c02r1"}),
    onh=dict(participant="sub-03286", picked=("NOT1ACH", "OXF2PRI", "NOT2ING-98"),
             folders={"NOT1ACH": "onharmony/sub-03286_ses-NOT1ACH001", "OXF2PRI": "onharmony/sub-03286_ses-OXF2PRI001",
                      "NOT2ING-98": "onharmony/sub-03286_ses-NOT2ING001"}),
    k=92, z_mm=20.0)


def rederive_selection():
    """The rules applied again to names and metadata only (no image or table value); must equal SELECTED."""
    if sha_file(RULES)[0] != RULES_SHA:
        stop("V12_DISPLAY_RULES.md differs from the recorded sha256")
    if sorted(R.WONG_SUBJECTS)[0] != SELECTED["a"]["participant"]:
        stop("multi-TE participant")
    for s, rel in SELECTED["a"]["folders"].items():
        own = REMEASURE / "wong" / SELECTED["a"]["participant"] / f"own_{s}"
        want = own if own.is_dir() else REMEASURE / "wong" / SELECTED["a"]["participant"]
        if REMEASURE / rel != want or not (want / f"roi_{s}.json").exists():
            stop(f"multi-TE folder of {s}")
    tong = sorted(R.TONG_SUBJECTS, key=lambda s: int(s.split("-")[1]))
    if tong[0] != SELECTED["tong"]["participant"] or tong != sorted(R.TONG_SUBJECTS):
        stop("Tong participant")
    centers = sorted({x[4:7] for x in R.TONG_SESSIONS})[:2]
    if tuple(centers) != SELECTED["tong"]["centers"]:
        stop("Tong centers")
    plan = frozen_copy("decomp/s4/plan.py", "N-250")
    onh_tab = None
    for node in ast.parse(plan.read_text(encoding="utf-8")).body:
        if isinstance(node, ast.Assign) and any(getattr(t, "id", None) == "ONH" for t in node.targets):
            onh_tab = ast.literal_eval(node.value)
    failed = {p.parent.name for p in (REMEASURE / "onharmony").glob("sub-*/FAILED.json")}
    excl = {f"{a}_{b}" for a, b in R.MAIN_EXCLUDE_ONH}
    main = [(a, b) for a, b in R.onh_sessions() if f"{a}_{b}" not in failed | excl]

    def unit(sub, ses):
        d = R.device_of(ses)
        return R.not2ing_protocol(sub) if d == "NOT2ING" else d

    chosen = None
    for sub in sorted({a for a, _ in main}):
        acq = sorted({unit(a, b) for a, b in main if a == sub})
        if len({sum(onh_tab[u][1]) / 2 for u in acq}) >= 3:
            chosen = sub
            break
    if chosen != SELECTED["onh"]["participant"]:
        stop("ON-Harmony participant")
    orders = [sorted(acq, key=lambda u, i=i: (f(onh_tab[u][1]), u)) for i, f in enumerate((min, max, lambda r: sum(r) / 2))]
    if not (orders[0] == orders[1] == orders[2]):
        stop("ON-Harmony TE order depends on the card TE used")
    o = orders[2]
    if (o[0], o[(len(o) - 1) // 2], o[-1]) != SELECTED["onh"]["picked"]:
        stop(f"ON-Harmony acquisitions {o}")
    import nibabel as nib
    A = nib.load(str(R.FSLDIR / (R.FA_REF + ".nii.gz"))).affine
    Am = nib.load(str(R.FSLDIR / "data" / "standard" / "MNI152_T1_1mm_brain_mask.nii.gz")).affine
    if not np.allclose(A, Am):
        stop("template and mask grids differ")
    kk = np.linalg.inv(A) @ np.array([0.0, 0.0, SELECTED["z_mm"], 1.0])
    if abs(kk[2] - SELECTED["k"]) > 1e-9:
        stop("slice index")
    return dict(order=o, card_TE={u: onh_tab[u][1] for u in o}, n_main=len(main), failed=sorted(failed), excluded=sorted(excl),
                affine=A)


# ------------------------------------------------------------------ the maps (a) and (c)
IDX = ("NDI", "FWF", "ODI")                                       # FWF = FISO (fit_f2.py output name)


def map_items():
    it = []
    for s in SELECTED["a"]["series"]:
        d = SELECTED["a"]["folders"][s]
        for f in IDX:
            it.append(dict(key=f"a|{s}|{f}", log="EXT", session=f"{SELECTED['a']['participant']}_{s}", panel="a",
                           native=f"{d}/f2_{s}_{f}.nii.gz", warp=f"{d}/fa2std_warp.nii.gz", index=f))
    for c in SELECTED["tong"]["centers"]:
        d = SELECTED["tong"]["folders"][c]
        for f in IDX[:2]:
            it.append(dict(key=f"tong|{c}|{f}", log="EXT", session=pathlib.Path(d).name, panel="c",
                           native=f"{d}/f2_{f}.nii.gz", warp=f"{d}/fa2std_warp.nii.gz", index=f))
    for u in SELECTED["onh"]["picked"]:
        d = SELECTED["onh"]["folders"][u]
        for f in IDX[:2]:
            it.append(dict(key=f"onh|{u}|{f}", log="ONH", session=pathlib.Path(d).name, panel="c",
                           native=f"{d}/f2_{f}.nii.gz", warp=f"{d}/fa2std_warp.nii.gz", index=f))
    return it


def append_log(path, entries):
    LOCK.parent.mkdir(parents=True, exist_ok=True)
    with open(LOCK, "a") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        try:
            n0 = sum(1 for _ in open(path, encoding="utf-8")) if path.exists() else 0
            with path.open("a", encoding="utf-8") as fh:
                for e in entries:
                    fh.write(json.dumps(e, ensure_ascii=False, default=str) + "\n")
        finally:
            fcntl.flock(lk, fcntl.LOCK_UN)
    return list(range(n0 + 1, n0 + 1 + len(entries)))


def read_maps(basis):
    """Log every file first (under the lock), then applywarp exactly as the ROI extraction did, then keep slice k."""
    import nibabel as nib
    items = map_items()
    logs = {"EXT": X14 / "EXTERNAL_READ_LOG.jsonl", "ONH": X14 / "ONH_NEW_READ_LOG.jsonl"}
    files, entries = {}, {"EXT": [], "ONH": []}
    for it in items:
        for role in ("native", "warp"):
            rel = it[role]
            if rel in files:
                continue
            p = REMEASURE / rel
            hx, n = sha_file(p)
            files[rel] = dict(sha256=hx, bytes=n, log=it["log"])
            what = WHAT.format(it["panel"]) + (f" · {'FISO' if it['index'] == 'FWF' else it['index']} 지도(applywarp 입력)"
                                               if role == "native" else " · 변환 fa2std_warp(applywarp 입력)")
            entries[it["log"]].append(dict(시각=now(), utc=utc(), 세션=it["session"], 파일=os.path.relpath(p, REPO), sha256=hx,
                                           바이트=n, 무엇=what, 근거=basis, 스크립트=str(pathlib.Path(__file__).resolve())))
    lines = {}
    for k_, es in entries.items():
        ln = append_log(logs[k_], es)
        for e, l in zip(es, ln):
            lines[e["파일"]] = (logs[k_].name, l, e["시각"])
    CACHE.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, CUDA_VISIBLE_DEVICES="", OMP_NUM_THREADS="1", FSLDIR=str(R.FSLDIR), FSLOUTPUTTYPE="NIFTI_GZ",
               PATH=f"{R.FSLDIR / 'bin'}:{os.environ.get('PATH', '')}")
    ref = str(R.FSLDIR / R.FA_REF)
    slices, manifest = {}, dict(made=now(), k=SELECTED["k"], ref=os.path.relpath(ref, REPO.parent), items={}, files={})
    for it in items:
        src, warp = REMEASURE / it["native"], REMEASURE / it["warp"]
        for rel in (it["native"], it["warp"]):                    # the bytes that applywarp reads are the bytes logged
            if sha_file(REMEASURE / rel)[0] != files[rel]["sha256"]:
                stop(f"{rel} changed between logging and use")
        out = CACHE / ("std_" + it["key"].replace("|", "_"))
        cmd = [str(R.FSLDIR / "bin" / "applywarp"), "-i", str(src), "-r", ref, "-w", str(warp), "-o", str(out), "--interp=trilinear"]
        r = subprocess.run(["nice", "-n", "19", *cmd], env=env, capture_output=True, text=True)
        if r.returncode != 0:
            stop(f"applywarp {it['key']}: {r.stderr[-500:]}")
        outp = pathlib.Path(str(out) + ".nii.gz")
        img = nib.load(str(outp))
        if img.shape != (182, 218, 182) or not np.allclose(img.affine, nib.load(ref + ".nii.gz").affine):
            stop(f"{it['key']}: output grid")
        sl = np.asarray(img.dataobj[:, :, SELECTED["k"]], dtype=np.float64)
        if not np.isfinite(sl).all():
            stop(f"{it['key']}: non-finite values on the slice")
        slices[it["key"]] = sl
        oh, ob = sha_file(outp)
        manifest["items"][it["key"]] = dict(it, std_sha256=oh, std_bytes=ob, std_dtype=str(img.get_data_dtype()),
                                            command=" ".join(["applywarp"] + [os.path.relpath(c, REPO.parent) if c.startswith("/") else c
                                                                              for c in cmd[1:]]))
        outp.unlink()                                             # the 3-D resample is not kept; its sha256 is
    mask = np.asarray(nib.load(str(R.FSLDIR / "data" / "standard" / "MNI152_T1_1mm_brain_mask.nii.gz")).dataobj[:, :, SELECTED["k"]])
    slices["mask"] = mask.astype(np.float64)
    for rel, v in files.items():
        manifest["files"][rel] = dict(v, log_file=lines[os.path.relpath(REMEASURE / rel, REPO)][0],
                                      log_line=lines[os.path.relpath(REMEASURE / rel, REPO)][1],
                                      logged_at=lines[os.path.relpath(REMEASURE / rel, REPO)][2])
    np.savez_compressed(CACHE / "slices.npz", **{k_.replace("|", "__"): v for k_, v in slices.items()})
    manifest["slices_sha256"] = sha_file(CACHE / "slices.npz")[0]
    manifest["basis"] = basis
    return manifest


# ------------------------------------------------------------------ (b) ROI units at eight TEs
TES = (62, 72, 82, 92, 102, 112, 122, 132)


def b1_logged_sha():
    """sha256 of the ROI tables as logged by the run that stored the per-unit slopes (b1_wong_direction.py)."""
    want = {}
    with open(X14 / "REREAD_LOG.jsonl", encoding="utf-8") as fh:
        for line in fh:
            e = json.loads(line)
            if e.get("무엇") == "B1 Wong TE 궤적 방향(judge_wong.load)" and "/wong/sub-01/" in e.get("다시읽음", ""):
                want[pathlib.Path(e["다시읽음"]).name] = e["sha256"]
    return want


def read_b():
    sub = SELECTED["b"]["participant"]
    if JW.TES != TES or JW.TE_C != float(ROWS["N-244"]["원 값"]):
        stop("judge_wong TEs or centre differ")
    want = b1_logged_sha()
    coords, shas = {}, {}
    for t in TES:
        s = f"TE{t}"
        folder = SELECTED["a"]["folders"].get(s) or (f"wong/{sub}/own_{s}" if (REMEASURE / f"wong/{sub}/own_{s}").is_dir() else f"wong/{sub}")
        p = REMEASURE / folder / f"roi_{s}.json"
        r = subprocess.run([sys.executable, "-I", "-B", str(LEDGER / "reread_log.py"), str(p), WHAT.format("b") + f" · ROI 표 {s}"],
                           capture_output=True, text=True, env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"))
        if r.returncode != 0:
            stop(f"reread_log {p}: {r.stderr[-300:]}")
        b = p.read_bytes()
        h = hashlib.sha256(b).hexdigest()
        if h != r.stdout.strip() or h != want.get(p.name):
            stop(f"{p.name}: sha256 {h} differs from the logged one ({want.get(p.name)})")
        shas[s] = dict(file=os.path.relpath(p, REPO.parent), sha256=h)
        coords[s] = JT.unit_coords(json.loads(b.decode("utf-8")))
    units = sorted(set.intersection(*[set(coords[f"TE{t}"]) for t in TES]))
    if len(units) != int(ROWS["N-179"]["원 값"]):
        stop(f"{len(units)} units, ledger N-179 says {ROWS['N-179']['원 값']}")
    X = np.array([[coords[f"TE{t}"][u][0] for t in TES] for u in units])
    Y = np.array([[coords[f"TE{t}"][u][1] for t in TES] for u in units])
    B = json.loads(frozen_copy("data/remeasure/explore/b1/b1_wong_direction.json", "N-372").read_text(encoding="utf-8"))
    bu = {u["unit"]: u for u in B["by_participant"][sub]["units"]}
    if sorted(bu) != units:
        stop("units differ from the stored per-unit slopes")
    te = np.array(TES, float)
    worst = 0.0
    for i, u in enumerate(units):
        bx, by = JW.fit_unit(te, X[i])["slope"], JW.fit_unit(te, Y[i])["slope"]
        for got, ref in ((bx, bu[u]["bx"]), (by, bu[u]["by"])):
            d = abs(got - ref) / max(abs(ref), 1e-12)
            worst = max(worst, d)
            if d > 1e-9:
                stop(f"unit {u}: slope {got} != stored {ref}")
    return dict(units=units, X=X, Y=Y, shas=shas, worst_rel=worst, n_b1=B["by_participant"][sub]["n"])


# ------------------------------------------------------------------ cache
def load_or_read(reread):
    basis = G.require()                                           # stops (GuardError) unless the guard passes
    man_p = CACHE / "manifest.json"
    if man_p.exists() and not reread:
        man = json.loads(man_p.read_text(encoding="utf-8"))
        if sha_file(CACHE / "slices.npz")[0] != man["slices_sha256"] or sha_file(CACHE / "b.npz")[0] != man["b"]["npz_sha256"]:
            stop("cache files differ from the manifest")
        if man.get("script_read_part") != read_part_sha():
            stop("the read part of this script changed since the cache was made — run with --reread")
        z = np.load(CACHE / "slices.npz")
        slices = {k_.replace("__", "|"): z[k_] for k_ in z.files}
        zb = np.load(CACHE / "b.npz", allow_pickle=False)
        bdat = dict(units=[str(u) for u in zb["units"]], X=zb["X"], Y=zb["Y"])
        return man, slices, bdat, basis
    man = read_maps(basis)
    bdat = read_b()
    np.savez_compressed(CACHE / "b.npz", units=np.array(bdat["units"]), X=bdat["X"], Y=bdat["Y"])
    man["b"] = dict(roi=bdat["shas"], worst_rel_slope_diff=bdat["worst_rel"], n_units=len(bdat["units"]),
                    npz_sha256=sha_file(CACHE / "b.npz")[0])
    man["script_read_part"] = read_part_sha()
    (CACHE / "manifest.json").write_text(json.dumps(man, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    z = np.load(CACHE / "slices.npz")
    slices = {k_.replace("__", "|"): z[k_] for k_ in z.files}
    return man, slices, bdat, basis


def read_part_sha():
    src = pathlib.Path(__file__).read_text(encoding="utf-8")
    a, b = src.index("# ------------------------------------------------------------------ selection"), src.index("# ------------------------------------------------------------------ cache")
    return hashlib.sha256(src[a:b].encode("utf-8")).hexdigest()


# ------------------------------------------------------------------ windows
SHOWN = [f"a|{s}|{f}" for f in IDX for s in SELECTED["a"]["series"]]   # the displayed maps (v13: panel (a) only)


def windows(slices, keys_of=None):
    """V12_DISPLAY_RULES §1 «window»: per index, the 1st and 99th percentiles of the pooled voxels of all displayed maps of that index inside the template brain mask on the slice. keys_of(f) -> the maps pooled (default: the displayed ones)."""
    m = slices["mask"] > 0
    out = {}
    for f in IDX:
        keys = keys_of(f) if keys_of else [k_ for k_ in SHOWN if k_.endswith("|" + f)]
        if not keys_of and len(keys) != len(SELECTED["a"]["series"]):
            stop(f"{f}: {len(keys)} maps pooled, {len(SELECTED['a']['series'])} are displayed")
        vals = np.concatenate([slices[k_][m] for k_ in keys])
        lo, hi = np.percentile(vals, [1, 99])
        out[f] = dict(lo=float(lo), hi=float(hi), n_vox=int(vals.size), n_maps=len(keys), n_mask=int(m.sum()), keys=keys)
    return out


def windows_v12(slices):
    """For reference only (not drawn): the v12 windows, which pooled the panel (c) maps for NDI and FISO."""
    return windows(slices, lambda f: [k_ for k_ in slices if k_ != "mask" and k_.endswith("|" + f)])


# ------------------------------------------------------------------ drawing
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.image import AxesImage  # noqa: E402
from matplotlib.text import Text  # noqa: E402
from matplotlib.transforms import Bbox  # noqa: E402

FONT = "Liberation Sans"
try:
    font_manager.findfont(font_manager.FontProperties(family=FONT), fallback_to_default=False)
except ValueError:
    stop(f"font {FONT} is not installed")
plt.rcParams.update({
    "font.family": FONT, "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "axes.linewidth": 0.6, "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "image.interpolation": "nearest", "image.resample": False,
})
MM = 1 / 25.4
W = 140.0                                                         # one panel: the 140 mm width
MIN_PT = 7.0
INK, INK2 = FC.INK, FC.INK2
TW, GAP, CB_W = 32.0, 0.8, 2.0                                    # tile width, gap between tiles, colour bar width (mm)
NAMES = {"TE 62 ms", "TE 92 ms", "TE 132 ms"}                     # column heads printed in the figure (as the caption and N-827; S46 G3)
LABEL = {"NDI": "NDI", "FWF": "FISO", "ODI": "ODI"}
REG = {}                                                          # id(text) -> role (window | name | LR)
H_FIG = 140.0                                                     # set in draw()


def crop_box(mask):
    """World-mm box of the template brain mask on the slice, plus 3 voxels (display crop only)."""
    disp = mask[::-1, :].T                                        # rows = y (j), columns = world x increasing (neurological)
    rows, cols = np.nonzero(disp > 0)
    x0, x1 = cols.min() - 91 - 3, cols.max() - 91 + 3
    y0, y1 = rows.min() - 126 - 3, rows.max() - 126 + 3
    return (x0 - 0.5, x1 + 0.5), (y0 - 0.5, y1 + 0.5)


def show(ax, sl, lo, hi, box):
    disp = sl[::-1, :].T
    ax.imshow(disp, cmap="gray", vmin=lo, vmax=hi, interpolation="nearest", origin="lower",
              extent=(-91.5, 90.5, -126.5, 91.5))
    ax.set_xlim(*box[0])
    ax.set_ylim(*box[1])
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def ftext(fig, x_mm, y_mm, s, role=None, **kw):
    t = fig.text(x_mm / W, y_mm / H_FIG, s, **kw)
    if role:
        REG[id(t)] = role
    return t


def colorbar(fig, x_mm, y_mm, w_mm, h_mm, lo, hi):
    ax = fig.add_axes([x_mm / W, y_mm / H_FIG, w_mm / W, h_mm / H_FIG])
    ax.imshow(np.linspace(0, 1, 256)[:, None], cmap="gray", vmin=0, vmax=1, aspect="auto", origin="lower",
              interpolation="nearest", extent=(0, 1, 0, 1))
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_linewidth(0.4)
        s.set_color(INK2)
    tx = x_mm + w_mm + 0.8
    ftext(fig, tx, y_mm + h_mm, f"{hi:.3f}", role="window", fontsize=7, ha="left", va="top", color=INK)
    ftext(fig, tx, y_mm, f"{lo:.3f}", role="window", fontsize=7, ha="left", va="bottom", color=INK)
    return ax


def draw(slices, win):
    global H_FIG
    box = crop_box(slices["mask"])
    cw, ch = box[0][1] - box[0][0], box[1][1] - box[1][0]
    TH = TW * ch / cw
    TOP, BOT = 5.0, 1.5                                           # top edge -> top of the tiles (column names); bottom margin
    H_FIG = TOP + 3 * TH + 2 * GAP + BOT
    a_top = H_FIG - TOP
    block = 4.0 + 3 * TW + 2 * GAP + 1.5 + CB_W + 0.8 + 6.4       # row labels .. colour-bar labels
    xa = (W - block) / 2 + 4.0                                    # the block centred in the 140 mm width
    fig = plt.figure(figsize=(W * MM, H_FIG * MM))
    for j, s in enumerate(SELECTED["a"]["series"]):
        ftext(fig, xa + j * (TW + GAP) + TW / 2, a_top + 1.0, f"TE {s[2:]} ms", role="name", fontsize=7, ha="center", va="bottom", color=INK)
    tiles = []
    for i, f in enumerate(IDX):
        y = a_top - (i + 1) * TH - i * GAP
        ftext(fig, xa - 1.2, y + TH / 2, LABEL[f], fontsize=7, ha="right", va="center", color=INK, rotation=90)
        for j, s in enumerate(SELECTED["a"]["series"]):
            ax = fig.add_axes([(xa + j * (TW + GAP)) / W, y / H_FIG, TW / W, TH / H_FIG])
            show(ax, slices[f"a|{s}|{f}"], win[f]["lo"], win[f]["hi"], box)
            tiles.append((ax, f"a|{s}|{f}", f))
        colorbar(fig, xa + 3 * TW + 2 * GAP + 1.5, y + 0.12 * TH, CB_W, 0.76 * TH, win[f]["lo"], win[f]["hi"])
    y_first = a_top - TH
    ftext(fig, xa + 1.0, y_first + 1.0, "L", role="LR", fontsize=7, ha="left", va="bottom", color="white")
    ftext(fig, xa + TW - 1.0, y_first + 1.0, "R", role="LR", fontsize=7, ha="right", va="bottom", color="white")
    return fig, dict(box=box, TW=TW, TH=TH, H=H_FIG, xa=xa, tiles=tiles)


FORBID = re.compile(r"\bV4\b|T-[ABC]\b|\bC[37]\b|\bL0\b|\bF[12]\b|FWF|stage-2|grid L1|sub-|ses-|osf\.io|\b[0-9a-f]{7,40}\b")


def audit(fig, geo, slices, win):
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    fb = fig.bbox
    texts, small, outside, numbers, forb = [], [], [], [], []
    for t in fig.findobj(Text):
        s = t.get_text()
        if not s.strip() or not t.get_visible():
            continue
        if t.get_fontsize() < MIN_PT - 1e-9 or t.get_fontname() != FONT:
            small.append(f"{s} ({t.get_fontsize()} pt, {t.get_fontname()})")
        if FORBID.search(s):
            forb.append(s)
        bb = Text.get_window_extent(t, rend)
        if t.get_bbox_patch() is not None:
            bb = Bbox.union([bb, t.get_bbox_patch().get_window_extent(rend)])
        if bb.x0 < fb.x0 - 0.5 or bb.x1 > fb.x1 + 0.5 or bb.y0 < fb.y0 - 0.5 or bb.y1 > fb.y1 + 0.5:
            outside.append(s)
        role = REG.get(id(t))
        if re.search(r"\d", s):
            if role == "window":
                if not any(s in (f"{w['lo']:.3f}", f"{w['hi']:.3f}") for w in win.values()):
                    numbers.append(s)
            elif role == "name":
                if s not in NAMES:
                    numbers.append(s)
            else:
                numbers.append(s)
        texts.append((t, bb))
    if small:
        stop(f"texts below {MIN_PT} pt or not in {FONT}: {small}")
    if outside:
        stop(f"texts outside the figure: {outside}")
    if numbers:
        stop(f"printed numbers that are neither window values nor data names: {numbers}")
    if forb:
        stop(f"internal names in the figure: {forb}")
    over = [(a.get_text(), b.get_text()) for (a, b1), (b, b2) in itertools.combinations(texts, 2)
            if b1.x0 < b2.x1 - 0.5 and b2.x0 < b1.x1 - 0.5 and b1.y0 < b2.y1 - 0.5 and b2.y0 < b1.y1 - 0.5]
    if over:
        stop(f"overlapping texts: {over}")
    intr = []                                                     # no text (with its frame) reaches into a tile or colour bar
    for t, bb in texts:
        if REG.get(id(t)) == "LR":
            continue
        for ax in fig.axes:
            ab = ax.get_window_extent(rend)
            if bb.x0 < ab.x1 - 0.5 and ab.x0 < bb.x1 - 0.5 and bb.y0 < ab.y1 - 0.5 and ab.y0 < bb.y1 - 0.5:
                intr.append(t.get_text())
    if intr:
        stop(f"texts reaching into the tiles or colour bars: {intr}")
    # tiles: the nine displayed maps, each its cached slice, its row's window, the crop box, one scale on both image axes
    tiles = geo["tiles"]
    if sorted(k for _, k, _ in tiles) != sorted(SHOWN) or len(tiles) != 9:
        stop("tiles differ from the displayed maps")
    for ax, k, f in tiles:
        ims = [im for im in ax.get_children() if isinstance(im, AxesImage)]
        if len(ims) != 1:
            stop(f"{k}: {len(ims)} images")
        im = ims[0]
        if not np.array_equal(np.asarray(im.get_array()), slices[k][::-1, :].T):
            stop(f"{k}: the tile is not the cached slice")
        if im.get_clim() != (win[f]["lo"], win[f]["hi"]) or im.get_interpolation() != "nearest":
            stop(f"{k}: window or interpolation")
        if tuple(ax.get_xlim()) != tuple(geo["box"][0]) or tuple(ax.get_ylim()) != tuple(geo["box"][1]):
            stop(f"{k}: crop")
        ab = ax.get_window_extent(rend)
        px, py = ab.width / (geo["box"][0][1] - geo["box"][0][0]), ab.height / (geo["box"][1][1] - geo["box"][1][0])
        if abs(px - py) > 1e-6 * py:
            stop(f"{k}: image scale differs between axes ({px} vs {py} px per mm)")
    for f in IDX:
        if win[f]["n_maps"] != 3 or sorted(win[f]["keys"]) != sorted(k for _, k, g in tiles if g == f):
            stop(f"{f}: the window does not pool exactly the displayed maps of the row")
    if abs(fig.get_size_inches()[0] / MM - W) > 1e-6:
        stop("figure width")
    return len(texts)


def font_audit(pdf):
    """Every font in the PDF is Liberation Sans and embedded (TrueType, pdf.fonttype 42)."""
    data = pathlib.Path(pdf).read_bytes()
    names = set(re.findall(rb"/BaseFont\s*/([A-Za-z0-9+_.-]+)", data))
    bad = [n for n in names if not re.fullmatch(rb"(?:[A-Z]{6}\+)?LiberationSans(?:-Bold)?", n)]
    if not names or bad:
        stop(f"PDF fonts {sorted(names)}: not all Liberation Sans")
    if data.count(b"/FontFile2") < len(names):
        stop(f"PDF fonts {sorted(names)}: not all embedded")
    return sorted(n.decode() for n in names)


def write_values(man, sel, win, win12, geo, n_texts, fonts):
    L = ["# Fig. 1 (v13, panel (a) only) — values, sources and audit", "",
         f"- Made {now()} by make_F1.py; rules sha256 {RULES_SHA[:12]}; slice k = {SELECTED['k']} (MNI z = +{SELECTED['z_mm']:.0f} mm).",
         f"- Figure {W:.0f} × {geo['H']:.1f} mm; fonts in the PDF {fonts} (embedded).",
         f"- Display crop (world mm): x {tuple(float(v) for v in geo['box'][0])}, y {tuple(float(v) for v in geo['box'][1])}; "
         f"tile {geo['TW']:.1f} × {geo['TH']:.2f} mm; neurological orientation (left hemisphere on the left).",
         f"- Figure texts audited: {n_texts}; all ≥ {MIN_PT} pt {FONT}, inside, no overlap; printed numbers = window values or data names only.",
         f"- Selection re-derived from names and metadata (all v12 selections, as recorded): ON-Harmony main first sessions {sel['n_main']}; "
         "only panel (a) is drawn.",
         "", "## Windows drawn (display settings; printed on the colour bars to three decimals; no ledger row)", "",
         "Pooled: the three displayed maps of the row (TE62, TE92, TE132), voxels inside the MNI152 brain mask on slice k.", "",
         "| index | 1st percentile | 99th percentile | maps pooled | voxels pooled | mask voxels on slice |", "|---|---|---|---|---|---|"]
    for f in IDX:
        w = win[f]
        L.append(f"| {LABEL[f]} | {w['lo']!r} | {w['hi']!r} | {w['n_maps']} | {w['n_vox']} | {w['n_mask']} |")
    L += ["", "## For reference only (not drawn): the v12 windows, which also pooled the v12 panel (c) maps", "",
          "| index | 1st percentile | 99th percentile | maps pooled | shown as |", "|---|---|---|---|---|"]
    for f in IDX:
        w = win12[f]
        L.append(f"| {LABEL[f]} | {w['lo']!r} | {w['hi']!r} | {w['n_maps']} | [{w['lo']:.3f}, {w['hi']:.3f}] |")
    L += ["", "## Maps (drawn; no ledger row — drawn only)", "",
          "| key | native map (sha256) | warp (sha256) | template-space resample sha256 | log line |", "|---|---|---|---|---|"]
    for k_ in SHOWN:
        it = man["items"][k_]
        fn, fw = man["files"][it["native"]], man["files"][it["warp"]]
        L.append(f"| {k_.replace('FWF', 'FISO').replace('|', ' · ')} | `{it['native']}` ({fn['sha256'][:12]}) | `{it['warp']}` ({fw['sha256'][:12]}) | "
                 f"{it['std_sha256'][:12]} | {fn['log_file']}:{fn['log_line']} |")
    L += ["", f"Resample command (each map, as cmds.applywarp): `{man['items'][SHOWN[0]]['command']}`"]
    (OUT / "F1_values.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def main():
    reread = "--reread" in sys.argv
    sel = rederive_selection()
    man, slices, bdat, basis = load_or_read(reread)
    win = windows(slices)
    win12 = windows_v12(slices)
    fig, geo = draw(slices, win)
    n_texts = audit(fig, geo, slices, win)
    meta = {"CreationDate": None, "Producer": None}
    fig.savefig(OUT / f"{NAME}.pdf", metadata=meta, dpi=600)
    fig.savefig(OUT / f"{NAME}.png", dpi=300, metadata={"Software": None})
    fonts = font_audit(OUT / f"{NAME}.pdf")
    write_values(man, sel, win, win12, geo, n_texts, fonts)
    used = {man["items"][k_][r] for k_ in SHOWN for r in ("native", "warp")}
    (OUT / "F1_reads.json").write_text(json.dumps(dict(files={k_: v for k_, v in man["files"].items() if k_ in used}, basis=man["basis"],
                                                       made=man["made"], cache_slices_sha256=man["slices_sha256"]),
                                                  ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    print(f"AUDIT PASS F1_maps — texts {n_texts} · tiles {len(geo['tiles'])} · windows (displayed maps) " +
          " · ".join(f"{LABEL[f]} [{win[f]['lo']:.3f}, {win[f]['hi']:.3f}]" for f in IDX) +
          f" · figure {W:.0f} × {geo['H']:.1f} mm · fonts {fonts} · guard head {basis['head'][:12]}")


if __name__ == "__main__":
    main()
