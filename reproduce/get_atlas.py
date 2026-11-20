#!/usr/bin/env python
"""Check a copy of the atlas that step 3 needs and put it in place. Steps 1 and 2 do not need it.

    python reproduce/get_atlas.py --from PATH

PATH is the label image, or a folder that holds it. JHU-labels.xml must lie beside it. A file whose SHA-256 differs is
refused and nothing is written. The 50-label atlas of FSL 6.0.5 and later has a different SHA-256 and a different ROI definition.
"""
from __future__ import annotations

import argparse
import hashlib
import pathlib

FILES = {
    "JHU-ICBM-labels-1mm.nii.gz": "fac584ec75ff2a8631710d3345df96733ed87d9bde3387f5b462f8d22914ed69",
    "JHU-labels.xml": "2d76ce80d1b0a50dccda2698d5eec55c8984a7f1bb438f79111d67a26fc4dc1c",
}
OUT = pathlib.Path(__file__).resolve().parents[1] / "prereg" / "artifacts" / "atlas"


def read(name: str, source: str) -> bytes:
    path = pathlib.Path(source)
    if path.is_file():                                              # the label image itself: the XML lies beside it
        path = path if name.endswith(".nii.gz") else path.parent / name
    else:
        path = path / name
    return path.read_bytes()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from", dest="source", required=True, help="the label image, or a folder that holds it")
    args = ap.parse_args()
    data = {name: read(name, args.source) for name in FILES}
    for name, want in FILES.items():
        got = hashlib.sha256(data[name]).hexdigest()
        if got != want:
            raise SystemExit(f"{name}: SHA-256 {got} differs from the expected {want}; nothing was written")
    OUT.mkdir(parents=True, exist_ok=True)
    for name, blob in data.items():
        (OUT / name).write_bytes(blob)
        print(name, "OK", len(blob), "bytes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
