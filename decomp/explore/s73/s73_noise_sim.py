"""S73 §5 synthesis, noise and AMICO fit in envs/reference: noiseless signals × 1000, noncentral χ noise, F2 fit, voxel means per unit and replicate."""
from __future__ import annotations

import json
import os
import pathlib
import sys

import nibabel as nib
import numpy as np

REPO = pathlib.Path(os.environ.get("S73_REPO") or pathlib.Path(__file__).resolve().parents[3])
sys.path.insert(0, str(REPO / "decomp"))
from rules_env import cards as K  # noqa: E402
from rules_env import synthfit as SF  # noqa: E402


def noisy(S: np.ndarray, L: int, sigma: float, rng) -> np.ndarray:
    if sigma <= 0:
        return S.copy()
    a = S / np.sqrt(L)
    acc = np.zeros_like(S, dtype=np.float64)
    for _ in range(L):
        acc += (a + rng.normal(0.0, sigma, S.shape)) ** 2 + rng.normal(0.0, sigma, S.shape) ** 2
    return np.sqrt(acc)


def run_job(j: dict, work: pathlib.Path, nthreads: int) -> dict:
    card = K.card(j["card"])
    card["name"] = j["key"]
    cells = [dict(nu=u[0], odi=u[1], fiso=u[2], t2=tuple(j["t2"]), t1=tuple(j["t1"]), gen=j["gen"], fdir=int(fd), dpar=1.7)
             for u in j["units"] for fd in j["fdirs"]]
    sig, _, _ = SF.signals(card, cells, n_workers=16)                    # (unit × fiber direction, volume)
    nU, nF, nV, nR = len(j["units"]), len(j["fdirs"]), int(j["n_per_fdir"]), int(j["n_rep"])
    base = (sig.astype(np.float64) * 1000.0).reshape(nU, nF, 1, -1)
    rng = np.random.default_rng(int(j["seed"]))
    rows = np.empty((nR, nU, nF, nV, base.shape[-1]), np.float32)
    for r in range(nR):
        rows[r] = noisy(np.broadcast_to(base, (nU, nF, nV, base.shape[-1])), int(j["L"]), float(j["sigma"]), rng)
    out = work / j["key"].replace("|", "_").replace(":", "_")
    out.mkdir(parents=True, exist_ok=True)
    n = nR * nU * nF * nV
    nib.save(nib.Nifti1Image(rows.reshape(n, 1, 1, -1), np.eye(4)), str(out / "data.nii"))
    nib.save(nib.Nifti1Image(np.ones((n, 1, 1), np.uint8), np.eye(4)), str(out / "nodif_brain_mask.nii"))
    bvals, bvecs = np.asarray(card["bvals"], float), np.asarray(card["bvecs"], float)
    g = bvecs.copy()
    nrm = np.linalg.norm(g, axis=1)
    g[nrm > 0] /= nrm[nrm > 0, None]
    np.savetxt(out / "DWI.scheme", np.c_[g, bvals], fmt="%.6f", header="VERSION: BVECTOR", comments="")
    fit = SF.fit_chunk(out, dpar=1.7, nthreads=nthreads)
    res = {k: fit[k].reshape(nR, nU, nF * nV).mean(-1) for k in ("ndi", "fwf", "odi")}
    res["n_rows"] = n
    return res


def main(jp, op, work, nthreads=24) -> None:
    spec = json.loads(pathlib.Path(jp).read_text(encoding="utf-8"))
    out = {}
    for j in spec["jobs"]:
        r = run_job(j, pathlib.Path(work), int(nthreads))
        for k, v in r.items():
            out[f"{j['key']}|{k}"] = v
        print(f"끝 {j['key']} · 행 {r['n_rows']:,}", flush=True)
    np.savez(op, **out)


if __name__ == "__main__":
    main(*sys.argv[1:])
