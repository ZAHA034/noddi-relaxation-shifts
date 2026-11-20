"""b2k step 3 for the b2000 noise floor: reads voxels of first-session b0 and b2000 volumes and computes the ventricle mask and within-direction mean/SD (R)."""
import collections, datetime, fcntl, hashlib, json, os, pathlib, subprocess, sys, time

import numpy as np
from scipy import ndimage as ndi

REPO = pathlib.Path("<REPO>")
sys.path.insert(0, str(REPO / "decomp" / "remeasure"))
import rconfig as R  # noqa: E402
import guard as G  # noqa: E402
import nibabel as nib  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import b2k_core as C  # noqa: E402

DRAFTS = pathlib.Path("<WORK>/stage1_work/drafts")
X14 = DRAFTS / "x14"
HEADERS_SHA = "806180efe1369e0b39d62067e03a388fedc32001c13ebe545cac4bc276679380"   # §27-2 ㄱ
OUT = R.WORK / "explore" / "b2000"
MASKS = OUT / "masks"
TAG = "봉인 해제 — 사후 탐색 · 판정 불변 #34"
LOCK = OUT / "read_log.lock"
HDR = json.loads((OUT / "headers.json").read_text(encoding="utf-8"))
EXPECT = {"NOT1ACH": 20, "NOT2ING-89": 4, "NOT2ING-98": 16, "NOT4GEP": 10, "OXF1PRI": 20, "OXF2PRI": 20, "OXF3TRI": 10, "OXF4GEP": 9}
META = dict(read_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
            core_sha256=hashlib.sha256((HERE / "b2k_core.py").read_bytes()).hexdigest(),
            headers_sha256=hashlib.sha256((OUT / "headers.json").read_bytes()).hexdigest())
L_FINE = np.exp(np.linspace(np.log(0.5), np.log(32.0), 400))
R_FINE = np.array([C.chi_ratio(L) for L in L_FINE])


def now():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def sha(p: pathlib.Path) -> tuple[str, int]:
    h, n = hashlib.sha256(), 0
    with open(p, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
            n += len(blk)
    return h.hexdigest(), n


def append(path: pathlib.Path, entries: list[dict]) -> None:
    LOCK.parent.mkdir(parents=True, exist_ok=True)
    with open(LOCK, "a") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        try:
            with path.open("a", encoding="utf-8") as fh:
                for e in entries:
                    fh.write(json.dumps(e, ensure_ascii=False, default=str) + "\n")
        finally:
            fcntl.flock(lk, fcntl.LOCK_UN)


def diagnostics(b2: np.ndarray, b0: np.ndarray, V: np.ndarray) -> dict:
    """§27-2 descriptive diagnostics — (1) R with zeros included · excluded fraction (2) R of the inner part after one more erosion versus the rim (3) R after dropping outliers (|z| > 5 · MAD) of the per-volume V means."""
    out = {}
    X = b2[V]
    if X.size == 0:
        return out
    sd = X.std(1, ddof=1)
    ok = sd > 0
    if ok.sum():
        out["R_incl0"] = float(np.median(X[ok].mean(1) / sd[ok]))
    out["frac_excluded"] = float(((X <= 0).any(1) | (b0[V] <= 0).any(1)).mean())
    Vin = ndi.binary_erosion(V, structure=ndi.generate_binary_structure(3, 1), iterations=1)
    for k, Vk in (("inner", Vin), ("rim", V & ~Vin)):
        mk = C.floor_metrics(b2, Vk, b0)
        out[f"R_{k}"], out[f"n_{k}"] = mk["R"], mk["n_V"]
    good = ~((X <= 0).any(1) | (b0[V] <= 0).any(1)) & ok
    if good.sum() >= C.MIN_VOX:
        vm = X[good].mean(0)
        med = np.median(vm)
        mad = 1.4826 * np.median(np.abs(vm - med))
        z = (vm - med) / mad if mad > 0 else np.zeros_like(vm)
        flag = np.abs(z) > 5
        out["n_vol_flag"], out["vol_flag"] = int(flag.sum()), np.flatnonzero(flag).tolist()
        if flag.any() and (~flag).sum() >= 3:
            Y = X[good][:, ~flag]
            out["R_vol_excl"] = float(np.median(Y.mean(1) / Y.std(1, ddof=1)))
    return out


def one(name: str) -> dict:
    t0 = time.time()
    h = HDR["sessions"][name]
    sub, ses = name.split("_", 1)
    basis = G.require()
    f_nii = R.ONH_SEALED / sub / ses / "dwi" / f"{name}_dir-AP_dwi.nii.gz"
    f_bval = R.ONH_META / sub / ses / "dwi" / f"{name}_dir-AP_dwi.bval"
    f_json = R.ONH_META / sub / ses / "dwi" / f"{name}_dir-AP_dwi.json"
    entries, bad = [], []
    for kind, f in (("nii", f_nii), ("bval", f_bval), ("json", f_json)):
        hx, n = sha(f)
        what = f"{TAG} · b2000 단계 3 화소(gzip 전체를 풀고 b0 · b2000 부피만 씀)" if kind == "nii" else f"{TAG} · b2000 단계 3"
        if hx != h["sha256"][kind]:
            what += " · sha256 불일치 — 중단"
            bad.append(kind)
        entries.append(dict(시각=now(), utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"), 세션=name,
                            파일=os.path.relpath(f, REPO), sha256=hx, 바이트=n, 무엇=what, 근거=basis, 스크립트=__file__))
    append(X14 / "ONH_NEW_READ_LOG.jsonl", entries)
    if bad:
        raise RuntimeError(f"{name} sha256 가 단계 1 과 다르다: {bad}")
    bval = np.loadtxt(f_bval).ravel()
    i0, i2 = C.shells(bval)
    if not (i0.tolist() == h["i_b0"] and i2.tolist() == h["i_b2k"]):
        raise RuntimeError(f"{name} b0 · b2000 부피 위치가 단계 1 과 다르다")
    img = nib.load(str(f_nii))
    data = np.asanyarray(img.dataobj)                         # Header scale applied (this data has no scale)
    dtype = str(data.dtype)
    b0 = np.asarray(data[..., i0], np.float64)
    b2 = np.asarray(data[..., i2], np.float64)
    del data                                                  # Volumes outside b0 and b2000 are not used
    aff, zooms = img.affine, img.header.get_zooms()[:3]
    B0 = b0.mean(-1)
    V, info = C.ventricle_mask(B0, aff, zooms)
    m = C.floor_metrics(b2, V, b0)
    try:
        extra = dict(diag=diagnostics(b2, b0, V))            # Descriptive diagnostic (§27-2 · not used in the verdict)
    except Exception as ex:
        extra = dict(diag=dict(error=f"{type(ex).__name__}: {ex}"[:300]))
    if m["R"] is not None:
        Lid = float(np.interp(m["R"], R_FINE, L_FINE))
        extra.update(B0_V=float(B0[V].mean()), D2_V=float(b2[V].mean()), L_ideal_population=Lid,
                     sigma_stored_population=float(np.sqrt(m["mean_M2"] / (2 * Lid))))      # Population-curve reference (for display) · S2 · σ_q is output by the summary as a ROI-unit calibration table
    # ── Re-access (R21 «exploratory after the result») — the main metric is kept even on failure ──
    sens = {}
    try:
        wd = R.WORK / "onharmony" / name
        rr = []
        paths = dict(fwf=wd / "f2_FWF.nii.gz", bm=wd / "nodif_brain_mask.nii.gz", fa=wd / "f2_FA.nii.gz", snr=wd / "snr.json")
        for k, p in paths.items():
            hx, n = sha(p)
            rr.append(dict(시각=now(), 다시읽음=str(p), sha256=hx, 바이트=n, 무엇=f"b2000 단계 3 — 민감도 · S2(R21 «결과 뒤 탐색» · §27)", 스크립트=__file__))
        append(X14 / "REREAD_LOG.jsonl", rr)
        fw, bm, fa = nib.load(str(paths["fwf"])), nib.load(str(paths["bm"])), nib.load(str(paths["fa"]))
        same_grid = all(x.shape[:3] == B0.shape and int(np.prod(x.shape)) == int(np.prod(B0.shape)) and np.allclose(x.affine, aff, atol=1e-3)
                        for x in (fw, bm, fa))
        sens["same_grid"] = bool(same_grid)
        sens["tag"] = "R21 결과 뒤 탐색"
        if same_grid:
            FWF = np.asarray(fw.dataobj, np.float64).reshape(B0.shape)
            BM = np.asarray(bm.dataobj).reshape(B0.shape) > 0
            FA = np.asarray(fa.dataobj, np.float64).reshape(B0.shape)
            H = C.head_mask(B0)
            Cb, _, _ = C.central_box(H, C.ras_grid(B0.shape, np.asarray(aff, float)))
            Vf = ndi.binary_erosion((FWF >= 0.95) & BM & Cb, structure=ndi.generate_binary_structure(3, 1), iterations=1)
            mf = C.floor_metrics(b2, Vf, b0)
            W = BM & (FA >= 0.3)
            snr = float(json.loads(paths["snr"].read_bytes())["snr"])
            inter = int((V & Vf).sum())
            sens.update(fwf=dict(R=mf["R"], n_V=mf["n_V"], n_used=mf["n_used"], status=mf["status"], r_q25=mf.get("r_q25"), r_q75=mf.get("r_q75")),
                        dice=(2 * inter / (int(V.sum()) + int(Vf.sum()))) if (V.sum() + Vf.sum()) else None,
                        S2=dict(snr=snr, n_W=int(W.sum()), S_W=float(B0[W].mean()) if W.any() else None,
                                sigma_snr=float(B0[W].mean() / snr) if (W.any() and snr) else None))
            np.savez_compressed(MASKS / f"{name}.npz", V=np.argwhere(V).astype(np.int16), Vf=np.argwhere(Vf).astype(np.int16))
        else:
            np.savez_compressed(MASKS / f"{name}.npz", V=np.argwhere(V).astype(np.int16))

    except Exception as ex:
        sens = dict(error=f"{type(ex).__name__}: {ex}"[:300], tag="R21 결과 뒤 탐색")
        np.savez_compressed(MASKS / f"{name}.npz", V=np.argwhere(V).astype(np.int16))
    return dict(session=name, unit=h["unit"], device=h["device"], dtype=dtype, n_b0=int(len(i0)), n_b2k=int(len(i2)), mask=info,
                **m, **extra, sens=sens, sec=round(time.time() - t0, 1), made=now())


def one_safe(name: str) -> dict:
    """Confines the failure of one session to that session (a guard stop is raised unchanged)."""
    try:
        return one(name)
    except G.GuardError:
        raise
    except Exception as ex:
        h = HDR["sessions"][name]
        return dict(session=name, unit=h["unit"], device=h["device"], R=None, n_V=None, status="실패",
                    error=f"{type(ex).__name__}: {ex}"[:300], made=now())


def main(units: list[str]) -> None:
    from multiprocessing import get_context
    if not units or not set(units) <= set(EXPECT):
        raise SystemExit(f"단위 이름이 §27 과 다르다: {units}")
    if dict(collections.Counter(v["unit"] for v in HDR["sessions"].values())) != EXPECT:
        raise SystemExit("headers.json 세션 수가 §27 과 다르다")
    if META["headers_sha256"] != HEADERS_SHA:
        raise SystemExit("headers.json sha256 가 §27-2 고정값과 다르다(봉인 파일은 열지 않음)")
    dirty = subprocess.run(["git", "--no-optional-locks", "diff", "--quiet", "HEAD", "--", "x15/b2000", "EXPLORATORY_PLAN_2026-11-01.md"],
                           cwd=DRAFTS).returncode
    if dirty:
        raise SystemExit("drafts 의 x15/b2000 또는 계획서가 커밋판과 다르다 — §27-2 고정 커밋 뒤에만 읽는다(봉인 파일은 열지 않음)")
    drafts_head = subprocess.run(["git", "--no-optional-locks", "rev-parse", "HEAD"], cwd=DRAFTS, capture_output=True, text=True).stdout.strip()
    MASKS.mkdir(parents=True, exist_ok=True)
    path = OUT / "read_sessions.jsonl"
    done = set()
    if path.exists():
        for ln in path.read_text(encoding="utf-8").splitlines():
            try:
                r = json.loads(ln)
            except json.JSONDecodeError:
                continue
            if r.get("status") in ("값", "마스크 부족"):
                done.add(r["session"])
    names = [k for k, v in HDR["sessions"].items() if v["unit"] in units and k not in done]
    missing = [str(R.WORK / "onharmony" / k / f) for k in names for f in ("f2_FWF.nii.gz", "nodif_brain_mask.nii.gz", "f2_FA.nii.gz", "snr.json")
               if not (R.WORK / "onharmony" / k / f).is_file()]
    if missing:
        raise SystemExit(f"다시 읽을 파일이 없다(봉인 파일은 열지 않음): {missing[:5]} …")
    append(OUT / "b2k_read_meta.jsonl", [dict(made=now(), units=units, n_sessions=len(names), skipped_done=len(done), drafts_head=drafts_head, **META)])
    nproc = min(int(os.environ.get("B2K_NPROC", "2")), 12)
    t0 = time.time()
    with get_context("fork").Pool(nproc) as pool:
        for n, rec in enumerate(pool.imap_unordered(one_safe, names, chunksize=1), 1):
            append(path, [rec])
            print(f"{n}/{len(names)} {rec['session']} {rec.get('status')} R={rec['R']} nV={rec['n_V']} {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
