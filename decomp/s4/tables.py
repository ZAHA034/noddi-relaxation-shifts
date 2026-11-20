"""Builds s4 per-design compartment signal tables and interpolates them linearly per voxel by ODI and d⊥/d∥, in envs/reference."""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s4cfg as S  # noqa: E402

DPAR_T = tuple(S.DPAR_TRUE) + (3.4,)
ODI_T = np.round(np.arange(0.04, 0.50 + 1e-9, 0.01), 4)              # 47 points
RATIO_T = np.round(np.linspace(0.0, 1.0, 41), 4)                     # form d⊥/d∥ = 1 − ν · 41 points
NDIR = len(S.C.FIBRE_DIRS)


def build(design: str, bscale: float = 1.0):
    """bscale ≠ 1 is for the injection test (§7 v) — only the synthesis changes the b values (the fit uses the recorded b values)."""
    import synth                                                      # only in envs/reference (AMICO)
    pr = synth.Proto(S.DESIGNS[design]) if bscale == 1.0 else _rescaled(synth, S.DESIGNS[design], bscale)
    t0 = time.time()
    ic = np.zeros((len(DPAR_T), len(ODI_T), NDIR, len(pr.bvals)), np.float32)
    ec = np.zeros((len(DPAR_T), len(RATIO_T), len(ODI_T), NDIR, len(pr.bvals)), np.float32)
    for i, dp in enumerate(DPAR_T):
        for j, o in enumerate(ODI_T):
            k = S.C.odi_to_kappa(float(o))
            for d in range(NDIR):
                ic[i, j, d] = pr.E_ic(dp, k, d)
                for r, rt in enumerate(RATIO_T):
                    ec[i, r, j, d] = pr.E_ec(dp, dp * float(rt), k, d)
    iso = pr.E_iso().astype(np.float32)
    S.WORK.mkdir(parents=True, exist_ok=True)
    path = S.WORK / f"tables_{design}{_tag(bscale)}.npz"
    np.savez(path, ic=ic, ec=ec, iso=iso, dpar=np.array(DPAR_T), odi=ODI_T, ratio=RATIO_T, bvals=pr.bvals, bvecs=pr.bvecs)
    h = hashlib.sha256(path.read_bytes()).hexdigest()
    print(json.dumps(dict(design=design, file=path.name, sha256=h, seconds=round(time.time() - t0, 1), shape_ic=ic.shape, shape_ec=ec.shape)), flush=True)
    return h


def _tag(bscale):
    return "" if bscale == 1.0 else f"_b{bscale:.2f}"


def _rescaled(synth, name, bscale):
    """A Proto with only the b values multiplied by bscale — the scheme is built with the same code as synth.Proto."""
    bvals, bvecs = S.C.protocol(name)
    orig = S.C.protocol
    S.C.protocol = lambda n: (bvals * bscale, bvecs)
    try:
        return synth.Proto(name)
    finally:
        S.C.protocol = orig


def _w(grid, v):
    """Lower index and weight of linear interpolation (outside the grid, clipped to the end)."""
    v = np.clip(v, grid[0], grid[-1])
    i = np.clip(np.searchsorted(grid, v, side="right") - 1, 0, len(grid) - 2)
    return i, (v - grid[i]) / (grid[i + 1] - grid[i])


def interp_ic(T, di, odi, fd):
    """E_ic — di (d∥ index) · odi · fd (direction index) are (n,) → (n, measurements)."""
    j, w = _w(T["odi"], odi)
    A = T["ic"][di, j, fd]
    B = T["ic"][di, j + 1, fd]
    return A * (1 - w)[:, None] + B * w[:, None]


def interp_ec(T, di, ratio, odi, fd):
    """E_ec — ratio = d⊥/d∥ · bilinear interpolation → (n, measurements)."""
    r, wr = _w(T["ratio"], ratio)
    j, wo = _w(T["odi"], odi)
    E = T["ec"]
    v00, v01 = E[di, r, j, fd], E[di, r, j + 1, fd]
    v10, v11 = E[di, r + 1, j, fd], E[di, r + 1, j + 1, fd]
    wr, wo = wr[:, None], wo[:, None]
    return (v00 * (1 - wo) + v01 * wo) * (1 - wr) + (v10 * (1 - wo) + v11 * wo) * wr


def load(design: str, bscale: float = 1.0) -> dict:
    return dict(np.load(S.WORK / f"tables_{design}{_tag(bscale)}.npz"))


if __name__ == "__main__":
    sys.path.insert(0, str(S.C.RESULTS.parent / "stage1"))
    for d in S.DESIGNS:
        for bs in (1.0, 0.95, 1.05):
            build(d, bs)
