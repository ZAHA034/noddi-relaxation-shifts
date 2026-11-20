"""Test settings: no stale bytecode, and CPU only.

CPython checks a .pyc file only by the modification time and size of its source, so a file saved twice within one second can run old code.
Bytecode is therefore not written, and the cache path points to an empty folder, so sources are always compiled again.
This must run before any test module is imported.
"""
from __future__ import annotations

import sys
import tempfile

sys.dont_write_bytecode = True
# An empty folder that lives as long as the process
_PYCACHE_DIR = tempfile.TemporaryDirectory(prefix="dmri_cov_pycache_")
sys.pycache_prefix = _PYCACHE_DIR.name

collect_ignore_glob = ["mutation/*"]

# Public repository: data locations. The scheme files are in data/protocols (see reproduce/judge_v4_run.py).
import os  # noqa: E402
import pathlib  # noqa: E402

_REPO = pathlib.Path(__file__).resolve().parents[1]
os.environ["CUDA_VISIBLE_DEVICES"] = ""                              # the judge code stops if a GPU is visible
os.environ.setdefault("DECOMP_WORK", str(_REPO / "data"))
os.environ.setdefault("REMEASURE_WORK", str(_REPO / "data" / "remeasure"))
sys.path.insert(0, str(_REPO / "decomp" / "stage1"))
import config as _config  # noqa: E402

_config.PROTO_DIR = _REPO / "data" / "protocols"

# The atlas is not redistributed (THIRD_PARTY_NOTICES.md). Tests that only need the 27 unit names take them from a
# deposited ROI table; the one test of the label-to-unit mapping is skipped until the atlas is in place (reproduce/get_atlas.py).
import glob  # noqa: E402
import json  # noqa: E402

import pytest  # noqa: E402

sys.path.insert(0, str(_REPO / "decomp" / "remeasure"))
import rconfig as _rconfig  # noqa: E402

ATLAS_MISSING = not _rconfig.ATLAS_XML.exists()
_STUB_DIR = tempfile.TemporaryDirectory(prefix="atlas_labels_")
if ATLAS_MISSING:
    _table = sorted(glob.glob(str(_REPO / "data/remeasure/onharmony/*/roi.json")))[0]
    _names = list(json.load(open(_table, encoding="utf-8"))["\uc8fc"])
    _stub = pathlib.Path(_STUB_DIR.name) / "labels.xml"
    _stub.write_text("".join(f'<label index="{i}">{_names[(i - 1) % len(_names)]}</label>\n' for i in range(1, 49)), encoding="utf-8")
    _rconfig.ATLAS_XML = _stub


def pytest_collection_modifyitems(config, items):
    if ATLAS_MISSING:
        for item in items:
            if item.name.startswith("test_27_ROI_units"):
                item.add_marker(pytest.mark.skip(reason="needs the atlas: python reproduce/get_atlas.py --from <file>"))
