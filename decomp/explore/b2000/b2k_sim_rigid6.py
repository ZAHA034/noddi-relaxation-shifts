"""Variant of b2k_sim with rigid6 motion only (eddy-current term 0) for the b2000 noise-floor known-answer test, writing to plan_rigid6/."""
from __future__ import annotations

import datetime
import hashlib
import json
import math
import os
import pathlib
import sys
import time

import numpy as np
from scipy import ndimage as ndi

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import b2k_core as C  # noqa: E402

WORK = pathlib.Path("<WORK>/remeasure")
ONH = WORK / "onharmony"
ORIG = WORK / "explore" / "b2000"
OUT = ORIG / "plan_rigid6"                                         # §27-11 — does not overwrite the original bank · summary
X14 = pathlib.Path("<WORK>/stage1_work/drafts/x14")
L_GRID = (0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0)
N_DRAW = 200
BVAL = 2000.0
LPAR, DISO_T = 1.7e-3, 0.8e-3
VENDOR = {"NOT1ACH": "Philips", "NOT2ING-89": "Philips", "NOT2ING-98": "Philips", "NOT4GEP": "GE", "OXF4GEP": "GE",
          "OXF1PRI": "Siemens", "OXF2PRI": "Siemens", "OXF3TRI": "Siemens"}
RESAMPLE = {"Siemens": ("없음",), "Philips": ("없음", "선형", "3차"), "GE": ("없음", "선형", "3차")}
SUB = np.array([[sx, sy, sz] for sx in (-0.25, 0.25) for sy in (-0.25, 0.25) for sz in (-0.25, 0.25)])   # Subsample in units of voxel size


def now():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


# ── ROI unit information (header table · re-access: eddy_parameters · snr.json) ────────────────────
def unit_info(log: bool = True) -> dict:
    H = json.loads((ORIG / "headers.json").read_text(encoding="utf-8"))
    units, reread = {}, []
    for u in VENDOR:
        ss = {k: v for k, v in H["sessions"].items() if v["unit"] == u}
        shapes = [tuple(v["shape"][:3]) for v in ss.values()]
        shape = max(set(shapes), key=shapes.count)
        zooms = [tuple(v["zooms"][:3]) for v in ss.values()]
        zoom = max(set(zooms), key=zooms.count)
        trt = float(np.median([float(v["side"]["TotalReadoutTime"]) for v in ss.values()]))
        traces, snrs = [], []
        for name, v in ss.items():
            ep = ONH / name / "eddy_out.eddy_parameters"
            sj = ONH / name / "snr.json"
            b1, b2 = ep.read_bytes(), sj.read_bytes()
            reread += [(ep, b1), (sj, b2)]
            P = np.loadtxt(ep)
            assert P.shape == (v["n_vol"], 16), (name, P.shape)
            traces.append(dict(b0=P[v["i_b0"]], b2k=P[v["i_b2k"]]))
            snrs.append(float(json.loads(b2)["snr"]))
        units[u] = dict(shape=shape, zooms=zoom, n_b0=int(np.median([v["n_b0"] for v in ss.values()])), n_b2k=50, trt=trt,
                        traces=traces, snr=(min(snrs), 2.0 * max(snrs)), snr_rec=(min(snrs), max(snrs)), n_sessions=len(ss))
    if log:
        with (X14 / "REREAD_LOG.jsonl").open("a", encoding="utf-8") as fh:
            for p, b in reread:
                fh.write(json.dumps(dict(시각=now(), 다시읽음=str(p), sha256=hashlib.sha256(b).hexdigest(), 바이트=len(b),
                                         무엇="§27-11 계획서판 정답 시험(강체 6 만)의 움직임 · SNR 범위", 스크립트=__file__), ensure_ascii=False) + "\n")
    return units


# ── Model ───────────────────────────────────────────────────────────────────
def fib_dirs(n: int) -> np.ndarray:
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = math.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], 1)


def rot(a: np.ndarray) -> np.ndarray:
    cx, cy, cz = np.cos(a)
    sx, sy, sz = np.sin(a)
    Rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    Ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
    Rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
    return Rz @ Ry @ Rx


def draw(rng: np.random.Generator, unit: dict, vendor: str) -> dict:
    vs = rng.uniform(0.7, 1.5)
    p = dict(brain_axes=np.array([68.0, 82.0, 58.0]) * rng.uniform(0.92, 1.08), brain_c=rng.normal(0, 5.0, 3),
             vent_scale=vs, v3_scale=rng.uniform(0.7, 1.3), kappa=rng.uniform(1.6, 4.0), d_csf=rng.uniform(2.6e-3, 3.2e-3),
             puls=rng.uniform(0.0, 0.3), lperp=rng.uniform(0.3e-3, 0.7e-3), iso=rng.uniform(0.0, 0.5),
             snr=rng.uniform(*unit["snr"]), trace=int(rng.integers(len(unit["traces"]))),
             fib_v0=rng.normal(size=3), fib_ph=rng.uniform(0, 2 * math.pi, 3), sig_ph=rng.uniform(0, 2 * math.pi, 2),
             rs_off=rng.uniform(0.0, 0.5, 3), rs_grad=rng.uniform(-0.2, 0.2, 3))
    p["fib_v0"] /= np.linalg.norm(p["fib_v0"])
    return p


def vent_ellipsoids(p: dict) -> list[tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """(center, half-axes, rotation) — mm relative to the brain center · R(+x) A(+y) S(+z)."""
    s, s3 = p["vent_scale"], p["v3_scale"]
    E = []
    for side in (-1, 1):
        E.append((np.array([side * 7.5, 5.0, 18.0]), np.array([6.0, 28.0, 9.0]) * s, rot(np.array([0, 0, side * math.radians(10)]))))
        E.append((np.array([side * 16.0, -22.0, 10.0]), np.array([7.0, 12.0, 8.0]) * s, np.eye(3)))
    E.append((np.array([0.0, 0.0, -2.0]), np.array([2.0, 12.0, 8.0]) * s3, np.eye(3)))
    return E


def fractions(P: np.ndarray, p: dict) -> tuple[np.ndarray, np.ndarray]:
    """P (…, 3) = model coordinates relative to the brain center → (whether CSF, whether tissue)."""
    inb = ((P / p["brain_axes"]) ** 2).sum(-1) <= 1.0
    inv = np.zeros(P.shape[:-1], bool)
    for c, a, Rm in vent_ellipsoids(p):
        Q = (P - c) @ Rm          # Rotated frame (row vector · the columns of Rm are the axes)
        inv |= ((Q / a) ** 2).sum(-1) <= 1.0
    return inb & inv, inb & ~inv


def fiber(P: np.ndarray, p: dict) -> np.ndarray:
    ph = p["fib_ph"]
    w = np.stack([np.sin(P[..., 0] / 25 + ph[0]), np.sin(P[..., 1] / 30 + ph[1]), np.sin(P[..., 2] / 20 + ph[2])], -1)
    u = p["fib_v0"] + 0.8 * w
    return u / np.linalg.norm(u, axis=-1, keepdims=True)


def ec_field(X: np.ndarray, e: np.ndarray) -> np.ndarray:
    """Eddy-current field (Hz) of the eddy second-order model — columns 7–16: linear 3 · second order 6 · constant 1 (coordinates in mm relative to the field-of-view center)."""
    x, y, z = X[..., 0], X[..., 1], X[..., 2]
    return (e[0] * x + e[1] * y + e[2] * z + e[3] * x * x + e[4] * y * y + e[5] * z * z + e[6] * x * y + e[7] * x * z + e[8] * y * z + e[9])


def signals(X: np.ndarray, zooms, p: dict, unit: dict, kind: str, dirs: np.ndarray | None, rng: np.random.Generator) -> np.ndarray:
    """X (…, 3) voxel centers (mm relative to the field-of-view center) → noise-free signal for each volume (…, N). kind = b0 · b2k."""
    tr = unit["traces"][p["trace"]][kind]
    n = len(tr)
    out = np.empty(X.shape[:-1] + (n,), np.float32)
    eps = np.clip(rng.normal(0, p["puls"], n), -0.9, None)
    res = p["kappa"] * math.exp(-BVAL * p["d_csf"])
    for k in range(n):
        t, a, e = tr[k, :3], tr[k, 3:6], np.zeros(10)          # §27-11 plan version: eddy-current term 0 (only the 3 translations and 3 rotations of eddy_parameters)
        Rm = rot(a)
        csf = np.zeros(X.shape[:-1], np.float32)
        tis = np.zeros(X.shape[:-1], np.float32)
        tsig = None
        for s in SUB:
            Xs = X + s * np.asarray(zooms)
            if kind == "b2k":
                Xs = Xs.copy()
                Xs[..., 1] -= ec_field(X, e) * unit["trt"] * zooms[1]      # Eddy current: shift along the phase-encoding (y) direction
            P = (Xs - t) @ Rm - p["brain_c"]                                # Inverse of a rigid body (row vector · applies R^T)
            c, ti = fractions(P, p)
            csf += c
            tis += ti
            if kind == "b2k" and tsig is None:
                u = fiber(P, p)
                g = dirs[k]
                cos2 = (u @ g) ** 2
                tsig = ((1 - p["iso"]) * np.exp(-BVAL * (p["lperp"] + (LPAR - p["lperp"]) * cos2)) + p["iso"] * math.exp(-BVAL * DISO_T)).astype(np.float32)
        csf /= len(SUB)
        tis /= len(SUB)
        if kind == "b0":
            out[..., k] = tis + csf * p["kappa"] * (1 + eps[k])
        else:
            out[..., k] = tis * tsig + csf * res * (1 + eps[k])
    return out


def resample(vol: np.ndarray, p: dict, zooms, how: str, reg: tuple[slice, slice, slice] | None = None) -> np.ndarray:
    """Resampling of the magnitude image corresponding to the within-scanner geometric adjustment — per axis a sub-voxel shift + a gentle slope (the same displacement for all volumes). If reg is given, only the output of that region is produced (the input is the grid as received)."""
    if how == "없음":
        return vol if reg is None else vol[reg]
    order = 1 if how == "선형" else 3
    shp = vol.shape[:3]
    rng_ = [np.arange(n, dtype=np.float64) for n in shp] if reg is None else [np.arange(s.start, s.stop, dtype=np.float64) for s in reg]
    idx = np.meshgrid(*rng_, indexing="ij")
    coords = []
    for ax in range(3):
        rel = (idx[ax] - (shp[ax] - 1) / 2) * zooms[ax] / 100.0
        coords.append(idx[ax] + p["rs_off"][ax] + p["rs_grad"][ax] * rel)
    out = np.empty(idx[0].shape + (vol.shape[3],), vol.dtype)
    for k in range(vol.shape[3]):
        v = ndi.spline_filter(vol[..., k], order=3) if order == 3 else vol[..., k]
        out[..., k] = ndi.map_coordinates(v, coords, order=order, mode="nearest", prefilter=False)
    return out


def grid(unit: dict) -> tuple[np.ndarray, np.ndarray]:
    shp, z = unit["shape"], np.asarray(unit["zooms"], float)
    aff = np.diag(list(z) + [1.0])
    aff[:3, 3] = -(np.asarray(shp) - 1) / 2 * z
    return C.ras_grid(shp, aff), aff


def region(unit: dict, p: dict, X: np.ndarray) -> tuple[slice, slice, slice]:
    """Central box of the noise-free brain mask + a margin of 8 voxels — b2000 is generated only inside it (the mask V must lie inside it · anything outside is counted)."""
    inb = (((X - p["brain_c"]) / p["brain_axes"]) ** 2).sum(-1) <= 1.0
    Cb, _, _ = C.central_box(inb, X)
    ii = np.argwhere(Cb)
    lo = np.maximum(ii.min(0) - 8, 0)
    hi = np.minimum(ii.max(0) + 9, np.asarray(Cb.shape))
    return tuple(slice(int(a), int(b)) for a, b in zip(lo, hi))


def one_draw(args) -> list[dict]:
    u, i, unit = args
    vendor = VENDOR[u]
    rng = np.random.default_rng([list(VENDOR).index(u), i, 20261103])
    p = draw(rng, unit, vendor)
    X, aff = grid(unit)
    reg = region(unit, p, X)
    dirs = fib_dirs(unit["n_b2k"])
    dirs = dirs[rng.permutation(len(dirs))]
    t0 = time.time()
    S0 = signals(X, unit["zooms"], p, unit, "b0", None, rng)                 # Whole grid
    S2 = signals(X[reg], unit["zooms"], p, unit, "b2k", dirs, rng)           # Region only
    sig0 = 1.0 / p["snr"]
    ph = p["sig_ph"]
    smap = sig0 * (1 + 0.3 * np.sin(X[..., 0] / 40 + ph[0]) * np.cos(X[..., 1] / 50 + ph[1]))
    recs = []
    for L in L_GRID:
        n0 = C.ncchi(S0, smap[..., None], L, rng).astype(np.float32)
        n2 = C.ncchi(S2, smap[reg][..., None], L, rng).astype(np.float32)
        for how in RESAMPLE[vendor]:
            a0 = resample(n0, p, unit["zooms"], how)
            a2 = n2
            if how != "없음":
                full = np.zeros(X.shape[:3] + (n2.shape[3],), np.float32)
                full[reg] = n2
                a2 = resample(full, p, unit["zooms"], how, reg)
            B0 = a0.mean(-1)
            V, info = C.ventricle_mask(B0, aff, unit["zooms"])
            Vr = V[reg]
            out_of_reg = int(V.sum() - Vr.sum())
            m = C.floor_metrics(a2, Vr, a0[reg])
            recs.append(dict(unit=u, draw=i, L=L, resample=how, out_of_region=out_of_reg, n_K=info["n_K"], **m))
    for r in recs:
        r.update(snr=round(p["snr"], 3), kappa=round(p["kappa"], 3), puls=round(p["puls"], 3), vent_scale=round(p["vent_scale"], 3),
                 trace=p["trace"], sec=round(time.time() - t0, 1))
    return recs


def selftest() -> None:
    """Ideal case (motion · partial volume · residual signal 0): whether R is near the theoretical value in a large block of pure noise · whether the mask code finds the model ventricle."""
    rng = np.random.default_rng(1)
    for L in (0.5, 1, 2, 4, 8):
        M = C.ncchi(np.zeros((40, 40, 40, 50)), 1.0, L, rng)
        V = np.ones((40, 40, 40), bool)
        m = C.floor_metrics(M, V)
        print(f"L {L}: 중앙 R {m['R']:.3f} · 모집단 비 {C.chi_ratio(L):.3f} (N2 50 의 화소 비는 표본 SD 때문에 조금 다르다)")
    unit = dict(shape=(104, 104, 72), zooms=(2.0192, 2.0192, 2.0), n_b0=5, n_b2k=50, trt=0.069,
                traces=[dict(b0=np.zeros((5, 16)), b2k=np.zeros((50, 16)))], snr=(15.0, 15.0))
    p = draw(np.random.default_rng(2), unit, "Siemens")
    p.update(puls=0.0, vent_scale=1.0)
    X, aff = grid(unit)
    S0 = signals(X, unit["zooms"], p, unit, "b0", None, rng)
    B0 = C.ncchi(S0, 1.0 / 15, 4, rng).mean(-1)
    V, info = C.ventricle_mask(B0, aff, unit["zooms"])
    csf, _ = fractions(X - p["brain_c"], p)
    print("모형 뇌실 화소", int(csf.sum()), "· 마스크", info, "· 마스크 안 CSF 비율", round(float(csf[V].mean()), 3) if V.any() else None)


def bank(units_sel: list[str]) -> None:
    from multiprocessing import get_context
    U = unit_info()
    tasks = [(u, i, U[u]) for u in units_sel for i in range(N_DRAW)]
    OUT.mkdir(parents=True, exist_ok=True)
    tag = "_".join(units_sel) if len(units_sel) < len(VENDOR) else "all"
    path = OUT / f"kat_bank_{tag}.jsonl"
    t0 = time.time()
    with get_context("fork").Pool(int(os.environ.get("B2K_NPROC", "12"))) as pool, path.open("w", encoding="utf-8") as fh:
        for n, recs in enumerate(pool.imap_unordered(one_draw, tasks, chunksize=1), 1):
            for r in recs:
                fh.write(json.dumps(r, ensure_ascii=False) + "\n")
            if n % 20 == 0:
                fh.flush()
                print(f"{n}/{len(tasks)} {time.time() - t0:.0f}s", flush=True)
    meta = dict(made=now(), script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
                core_sha256=hashlib.sha256((HERE / "b2k_core.py").read_bytes()).hexdigest(), units=units_sel, n_draw=N_DRAW,
                L_grid=L_GRID, unit_params={u: {k: v for k, v in U[u].items() if k != "traces"} for u in units_sel}, bank=str(path))
    (OUT / f"kat_meta_{tag}.json").write_text(json.dumps(meta, ensure_ascii=False, indent=1, default=str) + "\n", encoding="utf-8")


if __name__ == "__main__":
    if sys.argv[1] == "selftest":
        selftest()
    elif sys.argv[1] == "bank":
        bank(sys.argv[2:] or list(VENDOR))
