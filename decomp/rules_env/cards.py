"""Defines the acquisition-condition cards (name to TE, TR, bvals, bvecs and scheme) for the ON-Harmony card units, the E5 acquisitions, the Wong multi-TE series and the Tong data."""
from __future__ import annotations

import functools
import itertools
import json
import pathlib

import numpy as np

from ._paths import use

use("decomp", "stage1")
import analyze2 as A2  # noqa: E402
import analyze3 as A3  # noqa: E402 — confirmed TR range of NOT2ING (same as the algebra card · FISO sign rule)
import config as C  # noqa: E402

_CARDS = pathlib.Path(__file__).resolve().parent / "cards"
ONH_UNITS = ("NOT1ACH", "NOT2ING-89", "NOT2ING-98", "NOT4GEP", "OXF1PRI", "OXF2PRI", "OXF3TRI", "OXF4GEP")   # analysis.CARD_UNITS
DESIGN_TO_PROTOCOL = {"S": "P2a", "V": "P2b"}
WONG_TES = (62, 72, 82, 92, 102, 112, 122, 132)


@functools.lru_cache(maxsize=None)
def _protocol(name: str):
    b, v = C.protocol(name)
    return np.asarray(b, float), np.asarray(v, float)


@functools.lru_cache(maxsize=None)
def _json(name: str) -> dict:
    return json.loads((_CARDS / f"{name}.json").read_text(encoding="utf-8"))


def _mid(lohi) -> float:
    return round((lohi[0] + lohi[1]) / 2.0, 6)


def card_meta(name: str) -> dict:
    """Only the names of the TE · TR · direction scheme (no arrays) — the bval · bvec files of P2a · P2b (outside the repository · in the same place as the step-1 code) are not read."""
    kind, _, key = name.partition(":")
    if kind == "ONH":
        if key not in A2.SCANNERS:
            raise KeyError(name)
        design, te, tr = A2.SCANNERS[key]
        tr = A3.NOT2ING_FIXED.get(key, tr)
        return dict(name=name, TE=_mid(te), TR=1000.0 * (tr[0] + tr[1]) / 2.0, scheme=DESIGN_TO_PROTOCOL[design])
    if kind == "E5":
        if key not in C.ACQ_E5:
            raise KeyError(name)
        scheme, te, tr = C.ACQ_E5[key]
        return dict(name=name, TE=float(te), TR=float(tr), scheme=scheme)
    if kind == "WONG" or name == "TONG":
        c = card(name)
        return {k: c[k] for k in ("name", "TE", "TR", "scheme")}
    raise KeyError(name)


def card(name: str) -> dict:
    kind, _, key = name.partition(":")
    if kind in ("ONH", "E5"):
        m = card_meta(name)
        b, v = _protocol(m["scheme"])
        return dict(m, bvals=b, bvecs=v)
    if kind == "WONG":
        te = int(key.removeprefix("TE")) if key.startswith("TE") and key[2:].isdigit() else None
        if te not in WONG_TES:
            raise KeyError(name)
        j = _json("WONG")
        return dict(name=name, TE=float(te), TR=float(j["TR"]), bvals=np.asarray(j["bvals"], float),
                    bvecs=np.asarray(j["bvecs"], float), scheme=j["scheme"])
    if name == "TONG":
        j = _json("TONG")
        return dict(name=name, TE=float(j["TE"]), TR=float(j["TR"]), bvals=np.asarray(j["bvals"], float),
                    bvecs=np.asarray(j["bvecs"], float), scheme=j["scheme"])
    raise KeyError(name)


def pairs(set_name: str) -> list[tuple[str, str]]:
    if set_name == "onh28":
        return [(f"ONH:{a}", f"ONH:{b}") for a, b in itertools.combinations(ONH_UNITS, 2)]
    if set_name == "e5":
        return [("E5:A", f"E5:{q}") for q in ("Q1", "Q2", "Q3", "Q4")]
    if set_name == "wong_end":
        return [("WONG:TE62", "WONG:TE132")]
    raise KeyError(set_name)
