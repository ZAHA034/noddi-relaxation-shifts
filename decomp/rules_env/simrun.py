"""Runs the simulation jobs in envs/reference (AMICO) only: reads the job bundle (json) written by analyses in the main environment as a subprocess call and writes an npz with ndi, fwf, odi, nu_app and f_app per job."""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np

from . import cards as K
from . import synthfit as SF


def run(job_path, out_path, workdir, nthreads: int = 32, n_workers: int = 48) -> None:
    spec = json.loads(pathlib.Path(job_path).read_text(encoding="utf-8"))
    out = {}
    for j in spec["jobs"]:
        card = K.card(j["card"])
        for k in ("TE", "TR"):
            if j.get(k) is not None:
                card[k] = float(j[k])
        card["name"] = j["key"]
        cells = [dict(c, t2=tuple(c["t2"]), t1=(tuple(c["t1"]) if c.get("t1") is not None else None)) for c in j["cells"]]
        res = SF.fit_cells(card, cells, pathlib.Path(workdir) / j["key"].replace("|", "_").replace(":", "_"), nthreads=nthreads, n_workers=n_workers)
        for k, v in res.items():
            out[f"{j['key']}|{k}"] = v
        print(f"끝 {j['key']} · 행 {len(cells):,}", flush=True)
    np.savez(out_path, **out)


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2], sys.argv[3])
