"""Reads the already computed ROI values (roi.json) with the same rules as the existing verdict code and returns per-scale ROI unit values, reusing the loaders and judge functions of decomp/remeasure; it adds only the scale transforms and the per-scale retest SD, and opens no images."""
from __future__ import annotations

import json
import math
import pathlib

import numpy as np

from ._paths import use

use("decomp", "stage1")
use("decomp", "fixed")
use("decomp", "remeasure")
import analysis as AN  # noqa: E402
import fiso_sign as FS  # noqa: E402
import judge_onh as JO  # noqa: E402
import judge_tong as JT  # noqa: E402
import judge_v4 as V4  # noqa: E402
import judge_wong as JW  # noqa: E402
import lib  # noqa: E402
import rconfig as R  # noqa: E402

SCALES = ("raw_ndi", "x", "raw_fiso", "logit_fiso", "y", "odi")


# ── Scale ────────────────────────────────────────────────────────────────────
def _rows(roi_res: dict, table: str):
    for u, row in roi_res[table].items():
        v = [row.get(k) for k in ("NDI", "FWF", "ODI")]
        if any(a is None or not math.isfinite(a) for a in v):
            continue                                          # Same ROI unit set as judge_tong.unit_coords
        yield u, v


def scale_values(roi_res: dict, scale: str, table: str = "주") -> dict[str, float]:
    """roi.json → {unit: value}. x · y are lib.coords as they are (clip 1e-4), and logit_fiso is analysis.logit_fiso."""
    if scale not in SCALES:
        raise ValueError(scale)
    out = {}
    for u, (n, f, o) in _rows(roi_res, table):
        if scale == "raw_ndi":
            out[u] = float(n)
        elif scale == "raw_fiso":
            out[u] = float(f)
        elif scale == "odi":
            out[u] = float(o)
        elif scale == "logit_fiso":
            out[u] = float(AN.logit_fiso(f))
        else:
            xy, _ = lib.coords(np.array([n]), np.array([f]), np.array([o]))
            out[u] = float(xy[0, 0] if scale == "x" else xy[0, 1])
    return out


def clipped(roi_res: dict, table: str = "주") -> set[str]:
    """Units clipped by lib.coords (NDI or FWF outside [1e-4, 1 − 1e-4])."""
    out = set()
    for u, (n, f, o) in _rows(roi_res, table):
        _, cl = lib.coords(np.array([n]), np.array([f]), np.array([o]))
        if cl[0]:
            out.add(u)
    return out


def retest_sd_scalar(rep: dict[str, list[dict[str, float]]], units) -> dict[str, tuple[float | None, int]]:
    """The same formula as judge_tong.retest_sd (variance between repeats for each subject → mean → square root · sum of degrees of freedom), applied to units with a single value."""
    out = {}
    for u in units:
        var, df = [], 0
        for sess in rep.values():
            v = [d[u] for d in sess if u in d]
            if len(v) >= 2:
                var.append(float(np.var(v, ddof=1)))
                df += len(v) - 1
        out[u] = (math.sqrt(sum(var) / len(var)), df) if var else (None, 0)
    return out


# ── Tong ────────────────────────────────────────────────────────────────────
def read_tong(work: pathlib.Path = R.WORK, fname: str = "roi.json") -> dict:
    """Same session choice as judge_tong.load (end marker · loader guard) — returns the original table of roi.json."""
    import guard as G  # noqa: E402
    G.require()
    roi, missing = {}, []
    for sb in R.TONG_SUBJECTS:
        for ses in R.TONG_SESSIONS:
            wd = work / "tong" / f"{sb}_{ses}"
            p = wd / fname
            if not (wd / "DONE.json").exists() or not p.exists():
                missing.append(f"{sb}_{ses}")
                continue
            roi[(sb, ses)] = json.loads(p.read_text(encoding="utf-8"))
    return dict(roi=roi, missing=missing)


def tong_structs(t: dict, scale: str | None = None) -> tuple[dict, dict]:
    """(values, rep) — same form as judge_tong.load. If scale is None, a bundle of (x, y, ODI) (for the verdict code); otherwise the values on that scale."""
    values, rep = {}, {}
    for (sb, ses), r in t["roi"].items():
        uc = JT.unit_coords(r) if scale is None else scale_values(r, scale)
        if ses in JT.REPEATS:
            rep.setdefault(sb, []).append(uc)
        if ses.endswith("r1"):
            values.setdefault(JT.session_center(ses), {})[sb] = uc
    return values, rep


# ── Wong ────────────────────────────────────────────────────────────────────
def read_wong(work: pathlib.Path = R.WORK) -> dict:
    """Same series choice as judge_wong.load (a series registered separately is own_<series>/) — original table of roi.json."""
    import guard as G  # noqa: E402
    G.require()
    roi, own, missing = {}, {}, []
    for sb in R.WONG_SUBJECTS:
        wd = work / "wong" / sb
        if not (wd / "DONE.json").exists():
            missing.append(sb)
            continue
        own[sb] = set(json.loads((wd / "DONE.json").read_text(encoding="utf-8")).get("따로_정합", []))
        for s in R.WONG_SERIES:
            p = (wd / f"own_{s}" if s in own[sb] else wd) / f"roi_{s}.json"
            if p.exists():
                roi.setdefault(sb, {})[s] = json.loads(p.read_text(encoding="utf-8"))
    return dict(roi=roi, own=own, missing=missing)


def wong_series(w: dict, scale: str | None = None) -> dict:
    """{subject: {series: {unit: value}}} — if scale is None, a bundle of (x, y, ODI) (for judge_wong.judge)."""
    return {sb: {s: (JT.unit_coords(r) if scale is None else scale_values(r, scale)) for s, r in ser.items()} for sb, ser in w["roi"].items()}


def s_rt_corrected(ser: dict[str, dict[str, float]]) -> tuple[float | None, int]:
    """(repeat (TE62R2 · TE132R2) − original series) from the unit difference, subtract the mean difference for each repeat pair and take √(Σe² / (2·(N − number of pairs))) — degrees of freedom N − number of pairs."""
    e, k = [], 0
    for r, o in JW.REPEATS.items():
        if r in ser and o in ser:
            us = sorted(set(ser[r]) & set(ser[o]))
            d = np.array([ser[r][u] - ser[o][u] for u in us], float)
            if len(d):
                e.append(d - d.mean())
                k += 1
    if not e:
        return None, 0
    e = np.concatenate(e)
    df = len(e) - k
    return (math.sqrt(float(np.sum(e ** 2)) / (2 * df)), df) if df > 0 else (None, 0)


# ── ON-Harmony ──────────────────────────────────────────────────────────────
def onh_states(work: pathlib.Path = R.WORK) -> dict:
    """Same session state as judge_onh.run · main-analysis sessions (MAIN_EXCLUDE left out) — the values are not read (processing records only)."""
    sref_keys = [(a, b) for a in R.SREF_SUBJECTS for b in R.SREF_SESSIONS]
    states = V4.session_states(work, list(dict.fromkeys(R.onh_sessions() + sref_keys)))
    by_unit, dropped = V4.usable_sessions(work)
    main_by_unit = {u: {s: v for s, v in d.items() if (s, v) not in R.MAIN_EXCLUDE_ONH} for u, d in by_unit.items()}
    return dict(states=states, main_by_unit=main_by_unit, dropped=dropped,
                usable_keys=[k for k, v in states.items() if v == V4.USABLE])


def read_onh(work: pathlib.Path, st: dict) -> dict:
    """roi.json of the usable sessions (judge_onh.read_rois) and whether each unit · pair can be judged · ΔTE · ΔTR of the pairs (fiso_sign.device_rows)."""
    import guard as G  # noqa: E402
    G.require()
    rois = JO.read_rois(work, st["usable_keys"])
    cnt = AN.common_counts(R.onh_sessions())
    rows_by_pair = {}
    for r in FS.device_rows(cnt):
        a, b = r["pair"].split("~")
        i, j = sorted((JO.DEVICES.index(a), JO.DEVICES.index(b)))
        rows_by_pair[f"{JO.DEVICES[i]}~{JO.DEVICES[j]}"] = r
    return dict(rois=rois, units=JO.units_of(rois, "주"), cnt=cnt, rows_by_pair=rows_by_pair)


def onh_device_values(by_unit: dict, rois: dict, scale: str, table: str = "주", not2ing: str = "both") -> dict:
    """Same form as judge_onh.device_values {scanner: {subject: {unit: value}}} — only the scale conversion is scale_values."""
    out = {d: {} for d in JO.DEVICES}
    for unit, subs in by_unit.items():
        if unit == "NOT2ING-89" and not2ing == "98":
            continue
        dev = JO.DT.device_of(unit)
        for sub, ses in subs.items():
            out[dev][sub] = scale_values(rois[f"{sub}_{ses}"], scale, table)
    return out


def onh_pairs(st: dict, o: dict, scale: str) -> dict:
    """Summary of the 21 scanner pairs (analysis.pair_summary · B − A) · pairs that can be judged (judge_onh.keep_pairs) · ΔTE · ΔTR of the pairs."""
    summ = JO.pair_values(onh_device_values(st["main_by_unit"], o["rois"], scale), o["units"])
    keep, moved = JO.keep_pairs(o["cnt"], summ)
    return dict(summ=summ, keep=keep, moved=moved, names=[JO.pair_name(k) for k in range(len(summ))])


def s_onh(states: dict, rois: dict, units, scale: str) -> tuple[dict, int]:
    """analysis.s_ref from the OXF1PRI repeats (the usable sessions among R.SREF_SUBJECTS × SREF_SESSIONS) — judge_onh.s_ref_f2 extended to scales."""
    vals = {}
    for sub in R.SREF_SUBJECTS:
        for ses in R.SREF_SESSIONS:
            key = f"{sub}_{ses}"
            if states.get(key) == V4.USABLE:
                vals.setdefault(sub, []).append(scale_values(rois[key], scale))
    return AN.s_ref(vals, units)
