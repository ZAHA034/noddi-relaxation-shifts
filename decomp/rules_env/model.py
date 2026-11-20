"""Implements the algebra of manuscript formulas (2), (3), (6), (7) and Appendix A proposition 4 (weights, apparent fractions, coordinates, translation and inversion), calling the weights, apparent fraction and coordinates of decomp/stage1/config.py as they are."""
from __future__ import annotations

import json
import math
import pathlib

import numpy as np

from ._paths import use

use("decomp", "stage1")
import config as C  # noqa: E402

REQUIRED_KEYS = ("acquisitions", "tissues", "groups")


def weights(TE, t2, TR=None, t1=None) -> np.ndarray:
    """Compartment weights (in, en, iso) — config.weights as they are (if TR is given, t1 must be given too)."""
    return C.weights(TE, t2, TR, t1)


def apparent(nu, fiso, w):
    return C.apparent(nu, fiso, w)


def coords(nu, fiso):
    return C.frac_coords(nu, fiso)


def translation(w) -> tuple[float, float]:
    """The shift of equation (6), (Δx, Δy) = (ln(w_in/w_en), ln(w_iso/w_in))."""
    return math.log(w[0] / w[1]), math.log(w[2] / w[0])


def invert(nu_app, f_app, w) -> tuple[float, float]:
    """Inverse of equation (3) — the true (ν, f) from the apparent (ν′, f′) and the weights."""
    nu = 1.0 / (1.0 + (1.0 - nu_app) / nu_app * w[0] / w[1])
    A = nu * w[0] + (1.0 - nu) * w[1]
    odds = f_app / (1.0 - f_app) * A / w[2]
    return nu, odds / (1.0 + odds)


def _group_weights(site: dict, group: dict) -> np.ndarray:
    t1 = group.get("t1")
    tr = site.get("TR") if t1 is not None else None
    w = weights(site["TE"], group["t2"], tr, t1 if tr is not None else None)
    rho = group.get("rho")
    return w * np.asarray(rho, float) if rho is not None else w


def x_shift(site: dict, group: dict) -> float:
    """Logit NDI shift ln(w_in/w_en) of the group (T2 · T1 · proton density) at a site (TE · TR)."""
    w = _group_weights(site, group)
    return math.log(w[0] / w[1])


def prop4_change(site_a: dict, site_b: dict, ctrl: dict, pat: dict) -> float:
    """Proposition 4: change of the (patient − control) logit NDI difference between sites A − B after the additive calibration — (TE_A − TE_B)·(ΔR₂,p − ΔR₂,c) of equation (7) and the T1 term."""
    return (x_shift(site_a, pat) - x_shift(site_a, ctrl)) - (x_shift(site_b, pat) - x_shift(site_b, ctrl))


def load_input(path) -> dict:
    d = json.loads(pathlib.Path(path).read_text(encoding="utf-8"))
    missing = [k for k in REQUIRED_KEYS if k not in d]
    if missing:
        raise ValueError(f"입력에 없는 키: {', '.join(missing)}")
    return d
