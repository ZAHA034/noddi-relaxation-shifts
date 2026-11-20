"""Synthesizes noise-free signals from explicit tissue and relaxation values and fits them with decomp/stage1/amico_fit.py as a subprocess (d∥ 1.7, d_iso 3.0 µm²/ms, b0 threshold 10), to be called only from envs/reference (AMICO 2.1.1)."""
from __future__ import annotations

import os
import pathlib
from concurrent.futures import ProcessPoolExecutor
import subprocess
import sys

import nibabel as nib
import numpy as np

from ._paths import use

_STAGE1 = use("decomp", "stage1")
import config as C  # noqa: E402
import synth  # noqa: E402

CHUNK = 400_000
AMICO_FIT = _STAGE1 / "amico_fit.py"


class CardProto(synth.Proto):
    """Builds synth.Proto with the b values and directions of the card — only config.protocol, which Proto calls, is changed temporarily."""

    def __init__(self, card: dict):
        orig = synth.C.protocol
        synth.C.protocol = lambda name: (np.asarray(card["bvals"], float).copy(), np.asarray(card["bvecs"], float).copy())
        try:
            super().__init__(card["name"])
        finally:
            synth.C.protocol = orig


def noiseless_explicit(pr, nu, odi, fiso, t2, t1, TE, TR, gen, fdir, dpar=1.7):
    """The same formula as synth.noiseless — returns (S/S0, S0, ν′, f′). If TR or t1 is None, there is no T1 term."""
    a = C.weights(TE, t2) if TR is None or t1 is None else C.weights(TE, t2, TR, t1)
    nu_app, f_app = C.apparent(nu, fiso, a)
    kp = C.odi_to_kappa(odi)
    E_ic = pr.E_ic(dpar, kp, fdir)
    E_ec = pr.E_ec(dpar, C.dperp(gen, dpar, nu, nu_app), kp, fdir)
    S = (1 - fiso) * (nu * a[0] * E_ic + (1 - nu) * a[1] * E_ec) + fiso * a[2] * pr.E_iso()
    S0 = (1 - fiso) * (nu * a[0] + (1 - nu) * a[1]) + fiso * a[2]
    return S / S0, S0, nu_app, f_app


class _CachedProto:
    """A wrapper that remembers E_ec by (d∥, d⊥, κ, fiber direction) — the value is Proto.E_ec itself (the same array)."""

    def __init__(self, pr):
        self.pr, self._ec = pr, {}

    def E_ic(self, *a):
        return self.pr.E_ic(*a)

    def E_iso(self):
        return self.pr.E_iso()

    def E_ec(self, dpar, dperp, kappa, fdir):
        key = (dpar, dperp, kappa, fdir)
        if key not in self._ec:
            self._ec[key] = self.pr.E_ec(dpar, dperp, kappa, fdir)
        return self._ec[key]


_W = {}


def _work(args):
    card, cells = args
    if _W.get("name") != card["name"]:
        _W.update(name=card["name"], pr=_CachedProto(CardProto(card)))
    pr = _W["pr"]
    out = [noiseless_explicit(pr, TE=card["TE"], TR=card.get("TR"), **c) for c in cells]
    return (np.asarray([o[0] for o in out], np.float32), np.asarray([o[2] for o in out]), np.asarray([o[3] for o in out]))


def signals(card: dict, cells: list[dict], n_workers: int = 48):
    """Noise-free signal · ν′ · f′ — computed split over n_workers processes (order kept)."""
    k = max(1, min(n_workers, len(cells) // 200 or 1))
    step = -(-len(cells) // k)
    jobs = [(card, cells[i:i + step]) for i in range(0, len(cells), step)]
    if k == 1:
        parts = [_work(j) for j in jobs]
    else:
        with ProcessPoolExecutor(max_workers=k) as ex:
            parts = list(ex.map(_work, jobs))
    return (np.concatenate([p[0] for p in parts]), np.concatenate([p[1] for p in parts]), np.concatenate([p[2] for p in parts]))


def write_chunk(card: dict, cells: list[dict], outdir, n_workers: int = 48) -> pathlib.Path:
    """Writes the noise-free signal of cells (tissue · relaxation · generator · fiber direction) in the same file form as synth.py. TE · TR come from the card."""
    out = pathlib.Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    sig, nu_app, f_app = signals(card, cells, n_workers)
    if not np.isfinite(sig).all():
        raise ValueError("비유한 신호")
    n = len(sig)
    nib.save(nib.Nifti1Image((sig * 1000.0).reshape(n, 1, 1, -1), np.eye(4)), str(out / "data.nii"))
    nib.save(nib.Nifti1Image(np.ones((n, 1, 1), np.uint8), np.eye(4)), str(out / "nodif_brain_mask.nii"))
    bvals, bvecs = np.asarray(card["bvals"], float), np.asarray(card["bvecs"], float)
    np.savetxt(out / "bvals", bvals[None, :], fmt="%g")
    np.savetxt(out / "bvecs", bvecs.T, fmt="%.6f")
    g = bvecs.copy()
    nrm = np.linalg.norm(g, axis=1)
    g[nrm > 0] /= nrm[nrm > 0, None]
    np.savetxt(out / "DWI.scheme", np.c_[g, bvals], fmt="%.6f", header="VERSION: BVECTOR", comments="")
    np.savez(out / "rows.npz", nu_app=np.asarray(nu_app), f_app=np.asarray(f_app), card=card["name"])
    return out


def fit_chunk(outdir, dpar: float = 1.7, nthreads: int = 32, keep_data: bool = False) -> dict:
    out = pathlib.Path(outdir)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    with open(out / f"f2_{dpar}.log", "w") as log:
        subprocess.run([sys.executable, "-I", "-B", str(AMICO_FIT), str(out), str(dpar), str(nthreads)],
                       check=True, stdout=log, stderr=subprocess.STDOUT, env=env)
    res = np.load(out / f"f2_{dpar}.npz")
    if not keep_data:
        (out / "data.nii").unlink(missing_ok=True)
    return {k: np.asarray(res[k], float) for k in ("ndi", "fwf", "odi")}


def fit_cells(card: dict, cells: list[dict], workdir, nthreads: int = 32, dpar: float = 1.7, n_workers: int = 48) -> dict:
    """Synthesizes and fits cells in chunks of 400,000 rows → {ndi, fwf, odi, nu_app, f_app} (in the order of cells)."""
    work = pathlib.Path(workdir)
    parts = {k: [] for k in ("ndi", "fwf", "odi", "nu_app", "f_app")}
    for k0 in range(0, len(cells), CHUNK):
        d = write_chunk(card, cells[k0:k0 + CHUNK], work / f"c{k0 // CHUNK:04d}", n_workers=n_workers)
        fit = fit_chunk(d, dpar=dpar, nthreads=nthreads)
        meta = np.load(d / "rows.npz")
        for k in ("ndi", "fwf", "odi"):
            parts[k].append(fit[k])
        parts["nu_app"].append(meta["nu_app"])
        parts["f_app"].append(meta["f_app"])
    return {k: np.concatenate(v) for k, v in parts.items()}
