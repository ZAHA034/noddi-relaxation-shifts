"""Simulates one s4 case, one person's ROI at two sites: interpolated compartment signals, coil noise, AMICO F2 fit, extended coordinates and ROI mean."""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import sys
import time
import zlib

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s4cfg as S  # noqa: E402
import tables as TB  # noqa: E402

ORACLE_SD = 0.005                                           # §7 ii [design]
LOW_SHELL = 1000.0


def pixels(rng, centre, K=S.K_PIX):
    """ROI centers (n, 3: NDI · FISO · ODI) → voxels (n, K, 3) — heterogeneity SD · clipped to the widened range (§3-6 1)."""
    pix = centre[:, None, :] + rng.normal(0.0, 1.0, (len(centre), K, 3)) * S.PIX_SD
    lo = np.array([S.PIX_CLIP["NDI"][0], S.PIX_CLIP["FISO"][0], S.PIX_CLIP["ODI"][0]])
    hi = np.array([S.PIX_CLIP["NDI"][1], S.PIX_CLIP["FISO"][1], S.PIX_CLIP["ODI"][1]])
    return np.clip(pix, lo, hi)


def weights(TE, TR, t2, t1):
    """Compartment weights (n, 3) — the same formula as config.weights (row by row)."""
    return np.exp(-TE[:, None] / t2) * (1.0 - np.exp(-TR[:, None] / t1))


def signal(T, pix, a, gen, dpar_idx, fdir):
    """Noise-free signal (normalized by b0) (n·K, measurements) — the formula of SPEC_v1 §3 · compartment signals by table interpolation."""
    n, K, _ = pix.shape
    nu, f, o = pix[..., 0], pix[..., 1], pix[..., 2]
    ain, aen, aiso = a[:, 0:1], a[:, 1:2], a[:, 2:3]
    nu_app = nu * ain / (nu * ain + (1 - nu) * aen)
    g = np.asarray(gen)[:, None]
    ratio = np.where(g == "T-in", 1 - nu, np.where(g == "T-app", 1 - nu_app, 1 - S.C.NU_NONE))
    di, fd = np.repeat(dpar_idx, K), np.repeat(fdir, K)
    Eic = TB.interp_ic(T, di, o.ravel(), fd)
    Eec = TB.interp_ec(T, di, ratio.ravel(), o.ravel(), fd)
    nu, f = nu.ravel()[:, None], f.ravel()[:, None]
    ain, aen, aiso = (np.repeat(v[:, 0], K)[:, None] for v in (ain, aen, aiso))
    Sg = (1 - f) * (nu * ain * Eic + (1 - nu) * aen * Eec) + f * aiso * T["iso"][None, :]
    S0 = (1 - f) * (nu * ain + (1 - nu) * aen) + f * aiso
    return (Sg / S0).astype(np.float32)


def noisy(rng, Sg, snr, L):
    """Noncentral chi — signal in one channel and noise with σ = 1/SNR in L channels (the same distribution as synth.add_noise: the remaining channels are σ²·χ²_{2(L−1)})."""
    sig = (1.0 / snr)[:, None]
    re = Sg + rng.standard_normal(Sg.shape, dtype=np.float32) * sig
    acc = re ** 2 + (rng.standard_normal(Sg.shape, dtype=np.float32) * sig) ** 2
    acc += (sig ** 2) * rng.chisquare(2 * (L[:, None] - 1), Sg.shape).astype(np.float32)
    return np.sqrt(acc)


def amico_fit(study, subj, vol, bvals, bvecs, nthreads):
    """F2 — AMICO 2.1.1 NODDI (SPEC_v1 §4 · same settings as `decomp/stage1/amico_fit.py`). vol: (X, Y, 1, V) → (ndi, fwf, odi), each (X, Y)."""
    import amico
    import nibabel as nib
    d = os.path.join(study, subj)
    os.makedirs(d, exist_ok=True)
    nib.save(nib.Nifti1Image((vol * 1000.0).astype(np.float32), np.eye(4)), os.path.join(d, "DWI.nii"))
    nib.save(nib.Nifti1Image(np.ones(vol.shape[:3], np.uint8), np.eye(4)), os.path.join(d, "mask.nii"))
    g = bvecs.copy()
    nrm = np.linalg.norm(g, axis=1)
    g[nrm > 0] /= nrm[nrm > 0, None]
    np.savetxt(os.path.join(d, "DWI.scheme"), np.c_[g, bvals], fmt="%.6f", header="VERSION: BVECTOR", comments="")
    amico.setup()
    ae = amico.Evaluation(study, subj)
    ae.set_config("nthreads", nthreads)
    ae.set_config("BLAS_nthreads", 1)
    ae.load_data(dwi_filename="DWI.nii", scheme_filename="DWI.scheme", mask_filename="mask.nii", b0_thr=10)
    ae.set_model("NODDI")
    ae.model.set(dPar=S.DPAR_FIT * 1e-3, dIso=3.0e-3)
    ae.generate_kernels(regenerate=False)
    ae.load_kernels()
    ae.fit()
    ae.save_results()
    o = os.path.join(d, "AMICO", "NODDI")
    out = [nib.load(os.path.join(o, f"fit_{nm}.nii.gz")).get_fdata()[:, :, 0].astype(np.float32) for nm in ("NDI", "FWF", "ODI")]
    shutil.rmtree(d)
    return out


def oracle_obs(rng, P, idx):
    """Oracle observations (n, 6) = (ν′, f′_iso, ODI) of sites A · B + N(0, 0.005²)."""
    th = P["theta"][idx]
    out = []
    for s, cen in (("A", th[:, 3:6]), ("B", th[:, 3:6] + th[:, 0:3])):
        a = weights(P[f"acq_TE_{s}"][idx], P[f"acq_TR_{s}"][idx], P[f"t2{s}"][idx], P[f"t1{s}"][idx])
        nu_app, f_app = S.C.apparent(cen[:, 0], cen[:, 1], (a[:, 0], a[:, 1], a[:, 2]))
        out.append(np.c_[nu_app, f_app, cen[:, 2]])
    o = np.concatenate(out, 1)
    return o + rng.normal(0.0, ORACLE_SD, o.shape)


def chunk(P, idx, rng, Ttab, study, tag, nthreads, oracle_only=False):
    n = len(idx)
    res = dict(idx=idx, oracle=oracle_obs(rng, P, idx).astype(np.float32))
    if oracle_only:
        return res
    import lib  # noqa: E402 — decomp/stage1
    import synth  # noqa: E402
    th = P["theta"][idx]
    dpi = np.array([TB.DPAR_T.index(float(v)) for v in P["dpar"][idx]])
    raws, quality = {}, {}
    for s, cen in (("A", th[:, 3:6]), ("B", th[:, 3:6] + th[:, 0:3])):
        pix = pixels(rng, cen)
        a = weights(P[f"acq_TE_{s}"][idx], P[f"acq_TR_{s}"][idx], P[f"t2{s}"][idx], P[f"t1{s}"][idx])
        design = P[f"acq_design_{s}"][idx]
        snr = P[f"acq_snr_{s}"][idx]
        L = np.array([S.COIL_L[int(c)] for c in P[f"acq_coil_{s}"][idx]])
        raw = np.full((n, S.K_PIX, 7), np.nan, np.float32)
        for dsg in S.DESIGNS:
            m = design == dsg
            if not m.any():
                continue
            ii = np.flatnonzero(m)
            nom = Ttab[(dsg, 1.0)]
            Sg = np.empty((len(ii) * S.K_PIX, len(nom["bvals"])), np.float32)
            for bs in np.unique(P["bscale"][idx][ii]):
                sel = np.flatnonzero(P["bscale"][idx][ii] == bs)
                rows = (sel[:, None] * S.K_PIX + np.arange(S.K_PIX)[None, :]).ravel()
                Sg[rows] = signal(Ttab[(dsg, float(bs))], pix[ii[sel]], a[ii[sel]], P["gen"][idx][ii[sel]], dpi[ii[sel]], P[f"fdir_{s}"][idx][ii[sel]])
            Sn = noisy(rng, Sg, np.repeat(snr[ii], S.K_PIX), np.repeat(L[ii], S.K_PIX))
            vol = Sn.reshape(len(ii), S.K_PIX, -1).transpose(1, 0, 2)[:, :, None, :]          # (K, n_g, 1, V)
            ndi, fwf, odi = amico_fit(study[dsg], f"{tag}_{s}", vol, nom["bvals"], nom["bvecs"], nthreads)
            ext = synth.ext_coords(Sn.astype(float), nom["bvals"], nom["bvecs"], LOW_SHELL).reshape(len(ii), S.K_PIX, 4)
            raw[ii] = np.concatenate([np.stack([ndi.T, fwf.T, odi.T], -1), ext], -1)
        fin = np.isfinite(raw).all(-1)
        quality[f"finite_{s}"] = fin.sum(1)
        mu = np.array([raw[i][fin[i]].mean(0) for i in range(n)])
        Sig = np.array([np.cov(raw[i][fin[i]].T) for i in range(n)])
        roi = S.roi_noise(rng, mu, Sig, P["acq_neff"][idx])
        y, clipped = lib.coords(roi[:, 0], roi[:, 1], roi[:, 2], roi[:, 3:7], "C7")
        raws[s] = y.astype(np.float32)
        quality[f"clipped_{s}"] = clipped
    res["x_c7"] = np.concatenate([raws["A"], raws["B"]], 1)
    res["x_c3"] = np.concatenate([raws["A"][:, :3], raws["B"][:, :3]], 1)
    res.update(quality)
    return res


def run(name, csize=1000, nthreads=16, oracle_only=False, limit=None):
    sys.path.insert(0, str(S.C.RESULTS.parent / "stage1"))
    P = dict(np.load(S.WORK / f"plan_{name}.npz", allow_pickle=False))
    n = len(P["theta"]) if limit is None else min(limit, len(P["theta"]))          # limit: for checking and speed measurement (first n cases only)
    out = S.WORK / f"sim_{name}"
    out.mkdir(parents=True, exist_ok=True)
    log = S.WORK / "logs"
    log.mkdir(exist_ok=True)
    Ttab, study = {}, {}
    if not oracle_only:
        bsc = sorted({float(b) for b in P["bscale"]} | {1.0})
        for dsg in S.DESIGNS:
            for bs in bsc:
                Ttab[(dsg, bs)] = TB.load(dsg, bs)
            study[dsg] = str(S.WORK / "amico" / dsg)
            os.makedirs(study[dsg], exist_ok=True)
    nchunk = (n + csize - 1) // csize
    seed = S.SEEDS.get(name, S.SEEDS[name.split("_")[0]])          # A different seed sequence per dataset (together with the crc32 of the name)
    for c in range(nchunk):
        f = out / f"c{c:04d}.npz"
        if f.exists():
            continue
        t0 = time.time()
        idx = np.arange(c * csize, min(n, (c + 1) * csize))
        rng = np.random.default_rng([seed, zlib.crc32(name.encode()), c])
        r = chunk(P, idx, rng, Ttab, study, f"{name}_c{c:04d}", nthreads, oracle_only)
        np.savez(str(f) + ".tmp.npz", **r)
        os.replace(str(f) + ".tmp.npz", f)
        with open(log / f"sim_{name}.jsonl", "a", encoding="utf-8") as fh:
            fh.write(json.dumps(dict(UTC=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), chunk=c, of=nchunk, n=len(idx),
                                     seconds=round(time.time() - t0, 1), load1=os.getloadavg()[0])) + "\n")
    parts = [dict(np.load(out / f"c{c:04d}.npz")) for c in range(nchunk)]
    sim = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    assert np.array_equal(sim["idx"], np.arange(n))
    keep = ("theta", "cond", "level", "gen", "scen", "acq_kind") + tuple(k for k in ("layer", "cell") if k in P)
    sim.update({k: P[k][:n] for k in keep})
    path = S.WORK / f"sim_{name}.npz"
    np.savez(path, **sim)
    print(json.dumps(dict(file=path.name, n=n, sha256=hashlib.sha256(path.read_bytes()).hexdigest())), flush=True)


if __name__ == "__main__":
    argv = sys.argv[1:]
    lim = int(argv[argv.index("--limit") + 1]) if "--limit" in argv else None
    args = [a for i, a in enumerate(argv) if not a.startswith("--") and not (i and argv[i - 1] == "--limit")]
    run(args[0], int(args[1]) if len(args) > 1 else 1000, int(args[2]) if len(args) > 2 else 16, "--oracle-only" in argv, lim)
