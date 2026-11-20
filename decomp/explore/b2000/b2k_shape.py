"""Post hoc descriptive b2000 re-access (b2k): examines the cause of 0 values and the distribution shape, keeping aggregates only and leaving the M′ verdict unchanged."""
import collections, datetime, fcntl, hashlib, json, math, os, pathlib, subprocess, sys, time

import numpy as np
from scipy import ndimage as ndi
from scipy import optimize, special, stats

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import b2k_core as C  # noqa: E402
import b2k_summary as S  # noqa: E402

REPO = pathlib.Path("<REPO>")
DRAFTS = pathlib.Path("<WORK>/stage1_work/drafts")
X14 = DRAFTS / "x14"
OUT = S.OUT
TAG = "봉인 재열람 — 사후 기술 · 판정 불변 #34"
HEADERS_SHA = "806180efe1369e0b39d62067e03a388fedc32001c13ebe545cac4bc276679380"
SKEW_SHA = "2fb3f1668fe96adf9fc5986290c16a857bdaa399a4ecd4986bf3a3d97959de26"
EXPECT = S.EXPECT
SEED = 20261103
N_MIN = 40                      # If a voxel has fewer nonzero values than this, it is dropped from (b)
L_BOUNDS = (0.3, 30.0)
NBOOT_V = 200
HIST_EDGES = np.linspace(0.0, 2.0, 41)
SIZE_BINS = ((1, 1), (2, 9), (10, 99), (100, 10 ** 12))
THIRD, SIXTH, D_ML_MAX = 1.0 / 3.0, 1.0 / 6.0, 0.05
GAP_EXP_MIN = 20.0
DEADLINE_S = 3600
SIEMENS = ("OXF1PRI", "OXF2PRI", "OXF3TRI")


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
    lock = OUT / "read_log.lock"
    with open(lock, "a") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        try:
            with path.open("a", encoding="utf-8") as fh:
                for e in entries:
                    fh.write(json.dumps(e, ensure_ascii=False, default=str) + "\n")
        finally:
            fcntl.flock(lk, fcntl.LOCK_UN)


# ── Formulas of paths 2 and 3 ──────────────────────────────────────────────────────────────
def ml_parts(vals: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(n, Σ log m, Σ m²) of the nonzero values for each voxel."""
    n = np.asarray([len(v) for v in vals], float)
    sl = np.asarray([np.log(v).sum() for v in vals])
    s2 = np.asarray([(v * v).sum() for v in vals])
    return n, sl, s2


def loglik(L: float, n, sl, s2) -> float:
    """Log likelihood of the central chi (2L degrees of freedom) with a common L and a σ_v profile per voxel (s = 2σ² = Σm²/(n L))."""
    ll = n * math.log(2.0) + (2 * L - 1) * sl - n * L - n * L * np.log(s2 / (n * L)) - n * special.gammaln(L)
    return float(ll.sum())


def loglik_trunc(L: float, n, sl, s2, tau: float, iters: int = 30) -> float:
    """[§27-7a] Noncentral chi with the left side cut at τ · σ_v profiled again below the cut (fixed point s = Σm²/(n(L + h(τ²/s))) · h(x) = x^L e^{−x}/(Γ(L) Q(L, x)))."""
    s = s2 / (n * L)
    for _ in range(iters):
        x = tau * tau / s
        Q = special.gammaincc(L, x)
        h = np.exp(L * np.log(np.maximum(x, 1e-300)) - x - special.gammaln(L)) / np.maximum(Q, 1e-300)
        s = s2 / (n * (L + h))
    x = tau * tau / s
    ll = n * math.log(2.0) + (2 * L - 1) * sl - s2 / s - n * L * np.log(s) - n * special.gammaln(L) - n * np.log(np.maximum(special.gammaincc(L, x), 1e-300))
    return float(ll.sum())


def ml_fit(fun) -> tuple[float, str]:
    r = optimize.minimize_scalar(lambda x: -fun(math.exp(x)), bounds=tuple(math.log(b) for b in L_BOUNDS), method="bounded", options=dict(xatol=1e-4))
    L = float(math.exp(r.x))
    flag = "경계" if min(abs(r.x - math.log(b)) for b in L_BOUNDS) < 1e-3 else ""
    return L, flag


def skew_curve() -> dict:
    rng = np.random.default_rng(SEED)
    sk, ku = [], []
    for L in S.PURE_L:
        M = C.ncchi(np.zeros((S.PURE_NVOX, S.N2)), 1.0, L, rng)
        sk.append(float(np.median(stats.skew(M, axis=1, bias=False))))
        ku.append(float(np.median(stats.kurtosis(M, axis=1, bias=False))))
    out = dict(made=now(), L=list(S.PURE_L), skew=sk, kurt=ku, n_vox=S.PURE_NVOX, N2=S.N2, seed=SEED,
               script_sha256=sha(pathlib.Path(__file__))[0], core_sha256=sha(HERE / "b2k_core.py")[0])
    (OUT / "skew_curve_N2_50.json").write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    print(" · ".join(f"L{L:.3g} {s:.3f}" for L, s in zip(S.PURE_L, sk)))
    return out


def load_skew() -> tuple[np.ndarray, np.ndarray]:
    d = json.loads((OUT / "skew_curve_N2_50.json").read_text(encoding="utf-8"))
    return np.asarray(d["L"]), np.asarray(d["skew"])


def l_from_skew(sk: float, Ls: np.ndarray, curve: np.ndarray) -> tuple[float | None, str]:
    """Skewness decreases as L grows (monotone decrease) — invert · outside the range a censoring mark · if the value is not finite, no verdict is possible."""
    if not np.all(np.diff(curve) < 0):
        raise RuntimeError("왜도 곡선이 단조 감소가 아니다")
    if not math.isfinite(sk):
        return None, "판정 불가"
    if sk > curve[0]:
        return float(Ls[0]), "아래 검열"
    if sk < curve[-1]:
        return float(Ls[-1]), "위 검열"
    return float(np.exp(np.interp(sk, curve[::-1], np.log(Ls)[::-1]))), ""


def shape_stats(b2V: np.ndarray, rng: np.random.Generator, Ls_sk, cv_sk, tau: float | None) -> dict:
    """(b) — b2V = the b2000 values inside V, shape (voxels, 50). Only voxels with at least N_MIN nonzero values · tau = the fifth smallest positive value inside V (§27-7a)."""
    vals = [row[row > 0] for row in b2V]
    keep = [i for i, v in enumerate(vals) if len(v) >= N_MIN]
    zc = np.bincount((b2V == 0).sum(1), minlength=b2V.shape[1] + 1)
    out = dict(n_vox_V=int(len(vals)), n_vox_used=len(keep), zero_count_per_voxel_hist=zc.tolist())
    if len(keep) < C.MIN_VOX:
        out.update(status="화소 부족")
        return out
    vk = [vals[i] for i in keep]
    Xk = b2V[keep]
    n, sl, s2 = ml_parts(vk)
    R_same = float(np.median([v.mean() / v.std(ddof=1) for v in vk]))
    Lml, fml = ml_fit(lambda L: loglik(L, n, sl, s2))
    if tau is not None and tau > 0:
        Lmt, fmt = ml_fit(lambda L: loglik_trunc(L, n, sl, s2, tau))
    else:
        Lmt, fmt = Lml, fml
    bs = []
    for _ in range(NBOOT_V):
        idx = rng.integers(0, len(keep), len(keep))
        bs.append(ml_fit(lambda L: loglik(L, n[idx], sl[idx], s2[idx]))[0])
    sk = np.asarray([stats.skew(v, bias=False) for v in vk])
    ku = np.asarray([stats.kurtosis(v, bias=False) for v in vk])
    skm = float(np.nanmedian(sk)) if np.isfinite(sk).any() else float("nan")
    Lsk, fl = l_from_skew(skm, Ls_sk, cv_sk)
    sig = np.sqrt(s2 / (2 * Lml * n))
    z = np.concatenate([v / s for v, s in zip(vk, sig)])
    h_obs, _ = np.histogram(z, bins=HIST_EDGES)
    cdf = special.gammainc(Lml, HIST_EDGES ** 2 / 2)
    h_exp = np.diff(cdf) * len(z)
    e2 = float(h_exp[:2].sum())
    gap = float(h_obs[:2].sum() / e2) if e2 >= GAP_EXP_MIN else None
    hz, hneg, hnan = int((Xk == 0).sum()), int((Xk < 0).sum()), int(np.isnan(Xk).sum())
    out.update(status="값", R_same=R_same, L_ML=Lml, L_ML_flag=fml, L_ML_ci_sampling_only=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
               L_ML_trunc=Lmt, L_ML_trunc_flag=fmt, tau=tau, tau_over_sigma_med=(tau / float(np.median(sig))) if tau else None,
               skew_med=skm, skew_q=[float(np.nanpercentile(sk, 25)), float(np.nanpercentile(sk, 75))], skew_nan=int(np.isnan(sk).sum()),
               kurt_med=float(np.nanmedian(ku)), L_skew=Lsk, L_skew_flag=fl,
               hist_zero=hz, hist_obs=h_obs.tolist(), hist_exp=[float(x) for x in h_exp], n_z=int(len(z)), n_z_ge2=int((z >= 2.0).sum()),
               n_neg_used=hneg, n_nan_used=hnan, hist_identity=bool(Xk.size == len(z) + hz + hneg + hnan),
               gap_low2=gap, gap_exp_low2=e2, gap_obs_low2=int(h_obs[:2].sum()), gap_status="값" if gap is not None else "정의 안 됨(기대 계수 < 20)")
    return out


def zero_stats(b0: np.ndarray, b2: np.ndarray, V: np.ndarray, H: np.ndarray, integer: bool) -> dict:
    """(a) — based on the stored values. Per-volume and per-slice fraction of zeros · negative values · missing slices · runs · integer-valued counts."""
    out = {}
    for nm, A in (("b0", b0), ("b2k", b2)):
        out[f"{nm}_zero_frac_V"] = [float((A[..., k][V] == 0).mean()) for k in range(A.shape[-1])]
        out[f"{nm}_zero_frac_H"] = [float((A[..., k][H] == 0).mean()) for k in range(A.shape[-1])]
        out[f"{nm}_neg_V"] = [int((A[..., k][V] < 0).sum()) for k in range(A.shape[-1])]
        out[f"{nm}_neg_H"] = [int((A[..., k][H] < 0).sum()) for k in range(A.shape[-1])]
    nz = b2.shape[2]
    Hs = H.sum(axis=(0, 1))
    slice_frac = np.zeros((b2.shape[3], nz))
    for k in range(b2.shape[3]):
        z = (b2[..., k] == 0) & H
        slice_frac[k] = np.divide(z.sum(axis=(0, 1)), Hs, out=np.zeros(nz), where=Hs > 0)
    drop = slice_frac >= 0.5
    out["slice_zero_frac_H"] = np.round(slice_frac, 4).tolist()
    out["dropout_slices"] = [[int(k), int(s)] for k, s in np.argwhere(drop)]
    zV_total, zV_drop = 0, 0
    sizes = collections.Counter()
    st = ndi.generate_binary_structure(3, 1)
    for k in range(b2.shape[3]):
        z = (b2[..., k] == 0) & H
        zv = z & V
        cnt = int(zv.sum())
        if cnt == 0:
            continue
        zV_total += cnt
        zV_drop += int(zv[:, :, drop[k]].sum()) if drop[k].any() else 0
        lab, _ = ndi.label(z, structure=st)
        csize = np.bincount(lab.ravel())
        sv = csize[lab[zv]]
        for lo, hi in SIZE_BINS:
            sizes[f"{lo}-{hi if hi < 10 ** 12 else 'inf'}"] += int(((sv >= lo) & (sv <= hi)).sum())
    out.update(zero_V_b2k_total=zV_total, zero_V_on_dropout=zV_drop, zero_V_component_size=dict(sizes))
    v = b2[V].ravel()
    if integer:
        vi = v[(v >= 0) & (v < 128)].astype(np.int64)
        out["int_hist_0_127"] = np.bincount(vi, minlength=128)[:128].tolist()
    pos = np.sort(v[v > 0])
    out["min_positive_V"] = float(pos[0]) if pos.size else None
    out["fifth_positive_V"] = float(pos[4]) if pos.size >= 5 else None
    return out


# ── Reading ──────────────────────────────────────────────────────────────────────
def one(name: str) -> dict:
    sys.path.insert(0, str(REPO / "decomp" / "remeasure"))
    import rconfig as R  # noqa: E402
    import guard as G  # noqa: E402
    import nibabel as nib  # noqa: E402
    t0 = time.time()
    HDR = json.loads((OUT / "headers.json").read_text(encoding="utf-8"))
    h = HDR["sessions"][name]
    sub, ses = name.split("_", 1)
    basis = G.require()
    f_nii = R.ONH_SEALED / sub / ses / "dwi" / f"{name}_dir-AP_dwi.nii.gz"
    f_bval = R.ONH_META / sub / ses / "dwi" / f"{name}_dir-AP_dwi.bval"
    f_json = R.ONH_META / sub / ses / "dwi" / f"{name}_dir-AP_dwi.json"
    entries, bad = [], []
    for kind, f in (("nii", f_nii), ("bval", f_bval), ("json", f_json)):
        hx, n = sha(f)
        what = f"{TAG} · b2000 0 값 · 분포 모양(gzip 전체를 풀고 b0 · b2000 부피만 씀)" if kind == "nii" else f"{TAG} · b2000 0 값 · 분포 모양"
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
    slope, inter = float(img.dataobj.slope), float(img.dataobj.inter)
    integer = bool(np.issubdtype(img.get_data_dtype(), np.integer))
    data = np.asanyarray(img.dataobj)                                  # Scale-applied values (use these for R · shape)
    b0 = np.asarray(data[..., i0], np.float64)
    b2 = np.asarray(data[..., i2], np.float64)
    del data
    scaled = (slope != 1.0) or (inter != 0.0)
    b0s = np.rint((b0 - inter) / slope) if scaled else b0               # Stored values (use these for 0 · negative · integer coefficients)
    b2s = np.rint((b2 - inter) / slope) if scaled else b2
    B0 = b0.mean(-1)
    V, info = C.ventricle_mask(B0, img.affine, img.header.get_zooms()[:3])
    H = C.head_mask(B0)
    rec = dict(session=name, unit=h["unit"], device=h["device"], dtype_disk=str(img.get_data_dtype()), scl_slope=slope, scl_inter=inter,
               scaled_by_nibabel=bool(scaled), n_V=int(V.sum()))
    try:
        saved = np.load(OUT / "masks" / f"{name}.npz")["V"]
        Vs = np.zeros(V.shape, bool)
        Vs[tuple(saved.T)] = True
        same = bool((Vs == V).all())
    except Exception as ex:
        same, rec["mask_error"] = False, f"{type(ex).__name__}: {ex}"[:200]
    rec["mask_same"] = same
    m = C.floor_metrics(b2, V, b0)
    X = b2[V]
    sd = X.std(1, ddof=1)
    ok = sd > 0
    rec.update(R=m["R"], R_status=m["status"], R_incl0=float(np.median(X[ok].mean(1) / sd[ok])) if ok.any() else None)
    del X, sd, ok
    try:
        rec["zero"] = zero_stats(b0s, b2s, V, H, integer)
    except Exception as ex:
        rec["zero"] = dict(error=f"{type(ex).__name__}: {ex}"[:300])
    if same:
        try:
            Ls_sk, cv_sk = load_skew()
            tau = rec["zero"].get("fifth_positive_V") if "error" not in rec["zero"] else None
            tau_scaled = (tau * slope + inter) if (tau is not None and scaled) else tau
            rng = np.random.default_rng([SEED, int(hashlib.sha256(name.encode()).hexdigest()[:8], 16)])
            rec["shape"] = shape_stats(b2[V], rng, Ls_sk, cv_sk, tau_scaled)
        except Exception as ex:
            rec["shape"] = dict(error=f"{type(ex).__name__}: {ex}"[:300])
    else:
        rec["shape"] = dict(status="마스크 불일치")
    rec.update(status="값", sec=round(time.time() - t0, 1), made=now())
    return rec


def one_safe(name: str) -> dict:
    try:
        return one(name)
    except Exception as ex:
        if type(ex).__name__ == "GuardError":
            raise
        return dict(session=name, status="실패", error=f"{type(ex).__name__}: {ex}"[:300], made=now())


def read(units: list[str]) -> None:
    from multiprocessing import TimeoutError, get_context
    HDR = json.loads((OUT / "headers.json").read_text(encoding="utf-8"))
    if not units or not set(units) <= set(EXPECT):
        raise SystemExit(f"단위 이름이 §27 과 다르다: {units}")
    if dict(collections.Counter(v["unit"] for v in HDR["sessions"].values())) != EXPECT:
        raise SystemExit("headers.json 세션 수가 §27 과 다르다")
    if sha(OUT / "headers.json")[0] != HEADERS_SHA:
        raise SystemExit("headers.json sha256 가 고정값과 다르다(봉인 파일은 열지 않음)")
    if sha(OUT / "skew_curve_N2_50.json")[0] != SKEW_SHA:
        raise SystemExit("skew_curve_N2_50.json sha256 가 §27-7a 고정값과 다르다(봉인 파일은 열지 않음)")
    if subprocess.run(["git", "--no-optional-locks", "diff", "--quiet", "HEAD", "--", "x15/b2000", "EXPLORATORY_PLAN_2026-11-01.md"], cwd=DRAFTS).returncode:
        raise SystemExit("drafts 의 x15/b2000 또는 계획서가 커밋판과 다르다 — §27-7a 고정 커밋 뒤에만 연다(봉인 파일은 열지 않음)")
    untracked = subprocess.run(["git", "--no-optional-locks", "ls-files", "--others", "--exclude-standard", "--", "x15/b2000"], cwd=DRAFTS,
                               capture_output=True, text=True).stdout.strip()
    if untracked:
        raise SystemExit(f"x15/b2000 에 추적되지 않은 파일이 있다: {untracked[:200]}")
    head = subprocess.run(["git", "--no-optional-locks", "rev-parse", "HEAD"], cwd=DRAFTS, capture_output=True, text=True).stdout.strip()
    path = OUT / "shape_sessions.jsonl"
    done = set()
    if path.exists():
        for ln in path.read_text(encoding="utf-8").splitlines():
            try:
                r = json.loads(ln)
            except json.JSONDecodeError:
                continue
            if r.get("status") != "실패":
                done.add(r["session"])
    names = [k for k, v in HDR["sessions"].items() if v["unit"] in units and k not in done]
    missing = [k for k in names if not (OUT / "masks" / f"{k}.npz").is_file()]
    if missing:
        raise SystemExit(f"저장 마스크가 없다(봉인 파일은 열지 않음): {missing[:5]}")
    append(OUT / "shape_meta.jsonl", [dict(made=now(), units=units, n_sessions=len(names), skipped_done=len(done), drafts_head=head,
                                           shape_sha256=sha(pathlib.Path(__file__))[0], core_sha256=sha(HERE / "b2k_core.py")[0],
                                           summary_sha256=sha(HERE / "b2k_summary.py")[0], skew_curve_sha256=SKEW_SHA)])
    nproc = min(int(os.environ.get("B2K_NPROC", "2")), 12)
    t0 = time.time()
    pool = get_context("fork").Pool(nproc)
    it = pool.imap_unordered(one_safe, names, chunksize=1)
    k = 0
    try:
        while k < len(names):
            left = DEADLINE_S - (time.time() - t0)
            if left <= 0:
                raise TimeoutError
            rec = it.next(timeout=left)
            k += 1
            if "unit" not in rec:                                      # one_safe failure line — attach the ROI unit
                hh = HDR["sessions"][rec["session"]]
                rec.update(unit=hh["unit"], device=hh["device"])
            append(path, [rec])
            sh = rec.get("shape") or {}
            print(f"{k}/{len(names)} {rec['session']} {rec.get('status')} R={rec.get('R')} L_ML={sh.get('L_ML')} L_skew={sh.get('L_skew')} {time.time() - t0:.0f}s", flush=True)
        pool.close()
    except TimeoutError:
        pool.terminate()
        append(OUT / "shape_meta.jsonl", [dict(made=now(), note="1 시간 상한 도달 — 멈춤", done=k, of=len(names))])
        print(f"1 시간 상한 도달 — {k}/{len(names)} 에서 멈춤", flush=True)
    finally:
        pool.join()


# ── Aggregation ──────────────────────────────────────────────────────────────────
SIGN_CANDIDATES = {
    "(−, +)": "대칭 덧붙은 분산 · 부분 용적 조직 신호 · 0 가까운 값",
    "(≈0, +)": "남은 신호 · 방향 따라 변하는 부분 용적 · 약한 대칭 분산",
    "(+, −)": "σ 섞임 · 이상 부피 · 조직이 낀 부피 · 몇 부피에만 덧붙은 분산 · 결맞은 rank-1 · 낮은 값 잘림",
    "(−, −)": "화소마다 다른 L(열잡음)",
}


def sign_pair(d_ml: float, d_sk: float, tol: float = 0.05) -> str:
    a = "≈0" if abs(d_ml) < tol else ("+" if d_ml > 0 else "−")
    b = "+" if d_sk > 0 else "−"
    key = f"({a}, {b})"
    return key if key in SIGN_CANDIDATES else f"({a}, {b}) · 후보 표에 없음"


def summary() -> dict:
    for name, want in S.FIXED_SHA.items():                         # §27-2 (a) fixed input
        if sha(OUT / name)[0] != want:
            raise RuntimeError(f"{name} sha256 가 §27-2 고정값과 다르다")
    if sha(OUT / "skew_curve_N2_50.json")[0] != SKEW_SHA:
        raise RuntimeError("skew_curve_N2_50.json sha256 가 §27-7a 고정값과 다르다")
    rows, _ = S.read_lines(OUT / "shape_sessions.jsonl")
    recs = {r["session"]: r for r in rows}
    saved = {r["session"]: r for r in S.read_lines(OUT / "read_sessions.jsonl")[0]}
    B = S.load_bank()
    fixed = json.loads((OUT / "kat_summary_siemens.json").read_text(encoding="utf-8"))
    for u in S.NOMINAL_N:
        if json.loads(json.dumps(S.jsonable(S.calib(B, u)))) != fixed["calibration"][u]:
            raise RuntimeError(f"{u} 보정표가 §27-2 고정판과 다르다")
    dev, ell = {}, {k: {} for k in ("ratio", "same", "ML", "skew")}
    for u in S.VENDOR:
        rs = [r for r in recs.values() if r.get("unit") == u]
        if not rs:
            continue
        cal = S.calib(B, u)
        Lp, g = S.g_curve(cal, "없음")
        e = dict(n_expected=EXPECT[u], n_records=len(rs), n_failed=sum(r.get("status") == "실패" for r in rs),
                 mask_same=sum(bool(r.get("mask_same")) for r in rs),
                 R_same_as_step3=sum(1 for r in rs if r.get("R") is not None and abs(r["R"] - (saved.get(r["session"], {}).get("R") or -9)) <= 1e-9),
                 R_incl0_same_as_step3=sum(1 for r in rs if r.get("R_incl0") is not None and
                                           abs(r["R_incl0"] - ((saved.get(r["session"], {}).get("diag") or {}).get("R_incl0") or -9)) <= 1e-9),
                 scaled_sessions=sum(bool(r.get("scaled_by_nibabel")) for r in rs))
        sess = {}
        for r in rs:
            sh = r.get("shape") or {}
            if sh.get("status") != "값" or r.get("R") is None:
                continue
            Lr, fr = S.leff(r["R"], Lp, g)
            Li, fi = S.leff(r["R_incl0"], Lp, g) if r.get("R_incl0") is not None else (None, "없음")
            Ls_, fs = S.leff(sh["R_same"], Lp, g)
            sess[r["session"]] = dict(L_ratio=(Lr, fr), L_ratio_incl0=(Li, fi), L_ratio_same=(Ls_, fs), L_ML=(sh["L_ML"], sh["L_ML_flag"]),
                                      L_ML_trunc=(sh["L_ML_trunc"], sh["L_ML_trunc_flag"]), L_skew=(sh["L_skew"], sh["L_skew_flag"]))
        if sess:
            med = {k: float(np.median([s[k][0] for s in sess.values() if s[k][0] is not None])) for k in
                   ("L_ratio", "L_ratio_incl0", "L_ratio_same", "L_ML", "L_ML_trunc", "L_skew")}
            e["median"] = med
            e["censored"] = {k: sum(1 for s in sess.values() if s[k][1]) for k in ("L_ratio", "L_ratio_same", "L_ML", "L_ML_trunc", "L_skew")}
            fin = all(math.isfinite(med[k]) and med[k] > 0 for k in ("L_ratio_same", "L_ML_trunc", "L_skew"))
            if fin:
                d_ml, d_sk = math.log2(med["L_ML_trunc"] / med["L_ratio_same"]), math.log2(med["L_skew"] / med["L_ratio_same"])
                e.update(delta_ML=d_ml, delta_skew=d_sk, delta_ML_vs_main=math.log2(med["L_ML_trunc"] / med["L_ratio"]),
                         delta_skew_vs_main=math.log2(med["L_skew"] / med["L_ratio"]),
                         delta_session_median=dict(ML=float(np.median([math.log2(s["L_ML_trunc"][0] / s["L_ratio_same"][0]) for s in sess.values()])),
                                                   skew=float(np.median([math.log2(s["L_skew"][0] / s["L_ratio_same"][0]) for s in sess.values() if s["L_skew"][0]]))))
                if u in SIEMENS:
                    agree = abs(d_ml) <= THIRD and abs(d_sk) <= THIRD
                    e["agreement"] = "일치 — 이 감도 안에서 중심 카이 모양과 맞음" if agree else "어긋남 — 중심 카이 모양과 다름(원인은 이 검사로 가르지 않음)"
                    e["sign_pair"] = sign_pair(d_ml, d_sk)
                    e["sign_candidates"] = SIGN_CANDIDATES.get(e["sign_pair"], "후보 표에 없음")
                else:
                    e["agreement"] = "기술 · 재표본 미상(라벨 없음)"
            else:
                e["agreement"] = "판정 불가(장비 중앙이 유한하지 않음)"
            if u in S.NOMINAL_N:
                for k, key in (("ratio", "L_ratio"), ("same", "L_ratio_same"), ("ML", "L_ML_trunc"), ("skew", "L_skew")):
                    ell[k][u] = {S.subj(n): S.ell_bounds(*s[key]) if s[key][1] in ("아래 검열", "위 검열") else
                                 ((-math.inf, math.inf) if s[key][1] in ("경계", "판정 불가") else (math.log2(s[key][0]),) * 2)
                                 for n, s in sess.items() if s[key][0] is not None}
        ok_sh = [r["shape"] for r in rs if (r.get("shape") or {}).get("status") == "값"]
        if ok_sh:
            ho = np.sum([x["hist_obs"] for x in ok_sh], axis=0)
            he = np.sum([x["hist_exp"] for x in ok_sh], axis=0)
            hz = int(sum(x["hist_zero"] for x in ok_sh))
            e["hist"] = dict(edges=HIST_EDGES.tolist(), zero=hz, obs=ho.tolist(), exp=he.tolist(), exp_low2=float(he[:2].sum()),
                             gap_low2=float(ho[:2].sum() / he[:2].sum()) if he[:2].sum() >= GAP_EXP_MIN else None)
            e["skew_med"] = float(np.nanmedian([x["skew_med"] for x in ok_sh]))
            e["kurt_med"] = float(np.nanmedian([x["kurt_med"] for x in ok_sh]))
            e["tau_over_sigma_med"] = float(np.median([x["tau_over_sigma_med"] for x in ok_sh if x.get("tau_over_sigma_med") is not None]))
        z = [r["zero"] for r in rs if "zero" in r and "error" not in r["zero"]]
        if z:
            tot = sum(x["zero_V_b2k_total"] for x in z)
            sizes = collections.Counter()
            for x in z:
                sizes.update(x["zero_V_component_size"])
            contiguous = sizes.get("10-99", 0) + sizes.get("100-inf", 0)
            e["zero"] = dict(
                V_b2k_zero_total=tot, on_dropout_frac=(sum(x["zero_V_on_dropout"] for x in z) / tot) if tot else None,
                component_size=dict(sizes), contiguous_frac=(contiguous / tot) if tot else None, n_dropout_slices=sum(len(x["dropout_slices"]) for x in z),
                sessions_with_dropout=sum(1 for x in z if x["dropout_slices"]),
                b0_zero_frac_V_max=float(max(max(x["b0_zero_frac_V"]) for x in z)),
                b2k_zero_frac_V_med=float(np.median([np.mean(x["b2k_zero_frac_V"]) for x in z])),
                b2k_zero_frac_H_med=float(np.median([np.mean(x["b2k_zero_frac_H"]) for x in z])),
                neg_total=int(sum(sum(x["b2k_neg_H"]) + sum(x["b0_neg_H"]) for x in z)),
                fifth_positive_V_med=float(np.median([x["fifth_positive_V"] for x in z if x.get("fifth_positive_V") is not None])),
                int_hist_0_127=np.sum([x["int_hist_0_127"] for x in z if "int_hist_0_127" in x], axis=0).tolist() if any("int_hist_0_127" in x for x in z) else None)
            # Annotation rule (§27-7a · set in advance)
            if tot == 0:
                e["zero"]["annotation"] = "0 없음"
            elif e["zero"]["contiguous_frac"] is not None and e["zero"]["contiguous_frac"] >= 0.8:
                e["zero"]["annotation"] = "0 은 이어진 영역(가림 · 결손) — 0 을 뺀 R 에 잘림 치우침 없음"
            elif e.get("hist", {}).get("gap_low2") is not None and e["hist"]["gap_low2"] < 0.5 and (e.get("tau_over_sigma_med") or 0) > 0.02:
                e["zero"]["annotation"] = "0 은 낮은 값 잘림과 맞음 — 0 을 뺀 R 은 높아지는 쪽(L 을 크게)으로 치우침"
            else:
                e["zero"]["annotation"] = "원인 미상 — 0 포함 판정 같음(§27-4)"
        dev[u] = e
    out = dict(made=now(), devices=dev, sign_candidates=SIGN_CANDIDATES, script_sha256=sha(pathlib.Path(__file__))[0])
    if all(u in ell["ratio"] for u in ("OXF1PRI", "OXF2PRI")):
        main = json.loads((OUT / "b2k_result_all.json").read_text(encoding="utf-8"))["M_prime"]["contrasts"]["c1"]
        cs = {}
        for k in ("ratio", "same", "ML", "skew"):
            cs[k] = S.contrast(ell[k], "OXF2PRI", "OXF1PRI", 15, np.random.default_rng([SEED, 99]))
            cs[k].pop("subjects", None)
        out["c1"] = cs
        out["c1_ratio_equals_main"] = bool(cs["ratio"]["c"] == main["c"] and cs["ratio"]["ci"] == main["ci"])
        # Pair difference D = c₁^variant − c₁^same voxel — same subjects drawn with replacement (seed [SEED, 99])
        common = sorted(set(ell["same"]["OXF2PRI"]) & set(ell["same"]["OXF1PRI"]) & set(ell["ML"]["OXF2PRI"]) & set(ell["ML"]["OXF1PRI"])
                        & set(ell["skew"]["OXF2PRI"]) & set(ell["skew"]["OXF1PRI"]))
        if len(common) >= 15:
            def arr(k, u):
                return np.asarray([ell[k][u][s][0] for s in common])
            A = {k: (arr(k, "OXF2PRI"), arr(k, "OXF1PRI")) for k in ("same", "ML", "skew")}
            c = {k: float(np.median(a) - np.median(b)) for k, (a, b) in A.items()}
            rng = np.random.default_rng([SEED, 99])
            idx = rng.integers(0, len(common), size=(2000, len(common)))
            cb = {k: np.median(a[idx], axis=1) - np.median(b[idx], axis=1) for k, (a, b) in A.items()}
            D = {k: dict(point=c[k] - c["same"], ci=[float(np.percentile(cb[k] - cb["same"], 2.5)), float(np.percentile(cb[k] - cb["same"], 97.5))])
                 for k in ("ML", "skew")}
            out["D"] = dict(n_common=len(common), D_ML=D["ML"], D_skew=D["skew"])
            both = all(dev[u].get("agreement", "").startswith("일치") for u in ("OXF1PRI", "OXF2PRI"))
            drop = both and abs(D["skew"]["point"]) <= SIXTH and abs(D["ML"]["point"]) <= D_ML_MAX
            out["c1_caveat"] = ("뺀다 — 두 Prisma 일치 · |D_skew| ≤ ⅙ · |D_ML| ≤ 0.05" if drop else
                                "남긴다 — c₁ · c₁^ML · c₁^skew 를 함께 쓴다(방향 없이)")
        else:
            out["c1_caveat"] = "남긴다(쌍 부족)"
    (OUT / "shape_summary.json").write_text(json.dumps(S.jsonable(out), ensure_ascii=False, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    return out


# ── Tests before reading ─────────────────────────────────────────────────────────────
def selftest() -> None:
    rng = np.random.default_rng(5)
    Ls_sk, cv_sk = load_skew()
    print("(나) 순수 잡음 · 화소 σ 0.5 – 2 배 · 화소 600 × 50 · τ = 다섯째 양수")
    for L in (0.5, 1.0, 2.0, 4.0, 8.0):
        sig = rng.uniform(0.5, 2.0, (600, 1)) * 50
        M = C.ncchi(np.zeros((600, 50)), sig, L, rng)
        tau = float(np.sort(M[M > 0])[4])
        s = shape_stats(M, rng, Ls_sk, cv_sk, tau)
        print(f"  L {L}: L_ML {s['L_ML']:.3f} · L_ML,trunc {s['L_ML_trunc']:.3f} · L_skew {s['L_skew']:.3f} · R_same {s['R_same']:.3f} · 틈 {s['gap_low2']} ({s['gap_status']})")
    print("(나) 덧붙은 가우스 분산(σ_add = 0.3σ · 0.6σ) · L 1")
    for a in (0.3, 0.6):
        sig = rng.uniform(0.5, 2.0, (600, 1)) * 50
        M = np.abs(C.ncchi(np.zeros((600, 50)), sig, 1.0, rng) + rng.normal(0, 1, (600, 50)) * a * sig)
        s = shape_stats(M, rng, Ls_sk, cv_sk, float(np.sort(M[M > 0])[4]))
        print(f"  σ_add {a}: R_same {s['R_same']:.3f} · L_ML,trunc {s['L_ML_trunc']:.3f} · L_skew {s['L_skew']:.3f}")
    print("(가) 합성 0 무늬 · L 1 · σ 50 · 정수 저장 · 격자 40³ · 머리 = 전체 · V = 가운데 상자")
    V = np.zeros((40, 40, 40), bool)
    V[12:28, 12:28, 12:28] = True
    H = np.ones_like(V)
    for kind in ("없음", "슬라이스 결손", "낮은 값 잘림(8 DN)", "흩어진 무작위 0(0.5 %)", "이어진 가림"):
        b2 = np.rint(C.ncchi(np.zeros((40, 40, 40, 50)), 50.0, 1.0, rng))
        b0 = np.full((40, 40, 40, 5), 500.0)
        if kind == "슬라이스 결손":
            b2[:, :, 20, 7] = 0
            b2[:, :, 15, 30] = 0
        elif kind.startswith("낮은 값"):
            b2[b2 < 8.0] = 0
        elif kind.startswith("흩어진"):
            b2[rng.random(b2.shape) < 0.005] = 0
        elif kind == "이어진 가림":
            b2[14:22, 14:22, 14:22, 11] = 0
        z = zero_stats(b0, b2, V, H, True)
        s = shape_stats(b2[V], rng, Ls_sk, cv_sk, z["fifth_positive_V"])
        tot = z["zero_V_b2k_total"]
        cont = (z["zero_V_component_size"].get("10-99", 0) + z["zero_V_component_size"].get("100-inf", 0)) / tot if tot else 0
        print(f"  {kind}: V 0 {tot} · 결손 슬라이스 {len(z['dropout_slices'])} · 이어진 몫 {cont:.2f} · 틈 {s.get('gap_low2')} · "
              f"다섯째 양수 {z['fifth_positive_V']:.0f} · 정수 0 – 9 계수 {z['int_hist_0_127'][:10]} · L_ML {s['L_ML']:.2f} · L_ML,trunc {s['L_ML_trunc']:.2f} · L_skew {s['L_skew']:.2f}")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "skewcurve":
        skew_curve()
    elif cmd == "selftest":
        selftest()
    elif cmd == "read":
        read(sys.argv[2:])
    elif cmd == "summary":
        print(json.dumps(S.jsonable(summary()), ensure_ascii=False, indent=1)[:6000])
