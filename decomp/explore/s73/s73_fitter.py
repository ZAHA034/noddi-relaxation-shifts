"""S73 §3 (post hoc): regresses the fitter-output shift of x, y and raw NDI on the true tissue state in each cell of S1a and compares its tissue-dependent range with the retest width; stored stage-3 fits only."""
from __future__ import annotations

import hashlib
import json
import math
import os
import pathlib
import sys

import numpy as np

REPO = pathlib.Path(os.environ.get("S73_REPO") or pathlib.Path(__file__).resolve().parents[3])
sys.path[:0] = [str(REPO / "decomp"), str(REPO / "decomp" / "stage1")]
import analyze as A  # noqa: E402
import analyze3 as A3  # noqa: E402
import config as C  # noqa: E402
import lib  # noqa: E402
from rules_env import analyses as RA  # noqa: E402

Q_CLASS = RA.Q_CLASS
PAIRS = ("Q1", "Q2", "Q3", "Q4")


def reg(X: np.ndarray, d: np.ndarray) -> dict:
    """Least squares of d on [1, X] — coefficients (without the intercept) · R² · largest minus smallest fitted value."""
    D = np.c_[np.ones(len(X)), X]
    b, *_ = np.linalg.lstsq(D, d, rcond=None)
    p = D @ b
    ss = float(((d - d.mean()) ** 2).sum())
    r2 = float(1 - ((d - p) ** 2).sum() / ss) if ss > 0 else float("nan")
    return dict(intercept=float(b[0]), coef=[float(v) for v in b[1:]], r2=r2, delta_max=float(p.max() - p.min()))


def load_raw(exp: str = "E6"):
    """Fiber-direction means of fitted raw NDI and apparent ν′, grouped as in load_table (acquisition · generator · tissue · η), and the SHA-256 of the files read."""
    nT, nE = len(C.TISSUE_E5), len(C.ETA_E6)
    ndi = {a: {g: np.zeros((nT, nE)) for g in C.GENERATORS} for a in C.ACQ_E5}
    nua = {a: {g: np.zeros((nT, nE)) for g in C.GENERATORS} for a in C.ACQ_E5}
    cnt = {a: {g: np.zeros((nT, nE)) for g in C.GENERATORS} for a in C.ACQ_E5}
    for prot in ("P2a", "P2b"):
        cols, _ = A.load(exp, prot, A3.FIT)
        for a in C.ACQ_E5:
            for g in C.GENERATORS:
                m = (cols["acq"] == a) & (cols["gen"] == g)
                if m.any():
                    ix = (cols["tissue"][m], cols["eta6"][m])
                    np.add.at(ndi[a][g], ix, cols["ndi"][m].astype(float))
                    np.add.at(nua[a][g], ix, cols["nu_app"][m].astype(float))
                    np.add.at(cnt[a][g], ix, 1)
    for a in C.ACQ_E5:
        for g in C.GENERATORS:
            ndi[a][g] /= cnt[a][g]
            nua[a][g] /= cnt[a][g]
    files = sorted(p for prot in ("P2a", "P2b") for p in (A.FITS / exp / prot).glob("c*/*") if p.name in ("DONE", "rows.npz", "f2_1.7.npz"))
    sha = {str(p.relative_to(A.FITS)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    return ndi, nua, sha


def controls(xt: np.ndarray, Xa: np.ndarray, d_alg: np.ndarray, width: float) -> dict:
    """Definition §3.3 — (i) the algebraic shift (exactly common) · (ii) the same plus 0.10·(x_true − mean)."""
    r1 = reg(Xa, d_alg)
    r2 = reg(Xa, d_alg + 0.10 * (xt - xt.mean()))
    ok1 = all(abs(c) < 1e-10 for c in r1["coef"])
    ok2 = abs(r2["coef"][0] - 0.100) <= 1e-6 and r2["delta_max"] / width > 1
    return dict(i=dict(coef=r1["coef"], ok=ok1), ii=dict(coef_x=r2["coef"][0], ratio=r2["delta_max"] / width, ok=ok2))


def main(out_dir: str) -> dict:
    out = pathlib.Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    from rules_env import empirical as EM
    _, st, o, s = RA.onh_ctx()
    med = lambda sc: float(np.median([v for v in s[sc][0].values() if v is not None]))  # noqa: E731
    sd = dict(x=med("x"), y=med("y"), raw_ndi=med("raw_ndi"))
    width = {k: 1.96 * math.sqrt(2) * v for k, v in sd.items()}
    feat, frac = A3.load_table("E6")
    ndi, nua, sha = load_raw("E6")
    refs = list(C.refs_e5())
    idx = A3.lsets()["L1|cv0.03|rho0.0|iii"][0]
    truth = np.array([C.TISSUE_E5[r] for r in refs], float)          # (NDI, ODI, FISO)
    nu_t, odi_t, f_t = truth[:, 0], truth[:, 1], truth[:, 2]
    xt, yt = (np.asarray(v, float) for v in C.frac_coords(nu_t, f_t))
    Xa = np.c_[xt, yt, odi_t]                                        # explanatory variables of (a)
    Xb = np.c_[nu_t, f_t, odi_t]                                     # explanatory variables of (b)
    res = dict(definition="decomp/plans/S73_POST_HOC_PLAN_2026-11-16.md §3", retest_sd=sd, width=width, n_tissues=len(refs),
               eta=[A3.ETAN[e] for e in idx], generators=list(C.GENERATORS), pair_class={q: Q_CLASS[q] for q in PAIRS},
               input_sha256=sha, n_input_files=len(sha), controls={}, cells={}, by_acquisition={}, summary={})
    # failure controls (definition §3.3) — first, in every cell
    bad = []
    for g in C.GENERATORS:
        for q in PAIRS:
            for e in idx:
                d_alg = frac[q][g][refs, e, 0] - frac["A"][g][refs, e, 0]
                c = controls(xt, Xa, d_alg, width["x"])
                res["controls"][f"{g}|{q}|{A3.ETAN[e]}"] = c
                if not (c["i"]["ok"] and c["ii"]["ok"]):
                    bad.append(f"{g}|{q}|{A3.ETAN[e]}")
    res["controls_ok"] = not bad
    if bad:
        res["stopped"] = f"실패 대조가 기대대로 움직이지 않은 칸 {len(bad)}: {bad[:5]}"
        (out / "s73_fitter.json").write_text(json.dumps(res, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
        raise SystemExit("⛔ " + res["stopped"])
    for g in C.GENERATORS:
        for q in PAIRS:
            for e in idx:
                k = f"{g}|{q}|{A3.ETAN[e]}"
                dx = feat[q][g][refs, e, 0] - feat["A"][g][refs, e, 0]
                dy = feat[q][g][refs, e, 1] - feat["A"][g][refs, e, 1]
                draw = ndi[q][g][refs, e] - ndi["A"][g][refs, e]
                dmod = nua[q][g][refs, e] - nua["A"][g][refs, e]
                rx, ry, rr, rm = reg(Xa, dx), reg(Xa, dy), reg(Xb, draw), reg(Xb, dmod)
                res["cells"][k] = dict(x=dict(rx, ratio=rx["delta_max"] / width["x"]), y=dict(ry, ratio=ry["delta_max"] / width["y"]),
                                       raw_ndi=dict(rr, ratio=rr["delta_max"] / width["raw_ndi"]),
                                       raw_ndi_signal_model=dict(rm, ratio=rm["delta_max"] / width["raw_ndi"]))
    # (c) slopes by acquisition
    for g in C.GENERATORS:
        for e in idx:
            for a in ("A",) + PAIRS:
                xf = feat[a][g][refs, e, 0]
                b, a0 = np.polyfit(xt, xf, 1)
                res["by_acquisition"][f"{g}|{a}|{A3.ETAN[e]}"] = dict(alpha=float(a0), beta=float(b),
                                                                     max_abs_resid=float(np.abs(xf - (a0 + b * xt)).max()))
            for q in PAIRS:
                bq = res["by_acquisition"][f"{g}|{q}|{A3.ETAN[e]}"]["beta"] / res["by_acquisition"][f"{g}|A|{A3.ETAN[e]}"]["beta"]
                res["by_acquisition"][f"{g}|{q}|{A3.ETAN[e]}"]["beta_ratio_to_A"] = float(bq)
    for q in PAIRS:
        ks = [k for k in res["cells"] if k.split("|")[1] == q]
        row = {}
        for co in ("x", "y", "raw_ndi", "raw_ndi_signal_model"):
            r = np.array([res["cells"][k][co]["ratio"] for k in ks])
            row[co] = dict(n_cells=len(r), share_le_1=float(np.mean(r <= 1)), max_ratio=float(r.max()), min_ratio=float(r.min()))
        br = np.array([res["by_acquisition"][k.split("|")[0] + f"|{q}|" + k.split("|", 2)[2]]["beta_ratio_to_A"] for k in ks])
        row["beta_ratio_to_A"] = dict(min=float(br.min()), max=float(br.max()))
        res["summary"][f"{q} ({Q_CLASS[q]})"] = row
    allx = [v["x"]["max_ratio"] for v in res["summary"].values()]
    res["reading"] = dict(x_within_retest_in_every_pair_class=bool(all(v["x"]["share_le_1"] == 1.0 for v in res["summary"].values())),
                          R_x_max=float(max(allx)), R_raw_max=float(max(v["raw_ndi"]["max_ratio"] for v in res["summary"].values())),
                          pair_classes_exceeding_x=[k for k, v in res["summary"].items() if v["x"]["share_le_1"] < 1.0],
                          pair_classes_exceeding_y=[k for k, v in res["summary"].items() if v["y"]["share_le_1"] < 1.0])
    (out / "s73_fitter.json").write_text(json.dumps(res, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    return res


if __name__ == "__main__":
    r = main(sys.argv[1])
    print(json.dumps(dict(controls_ok=r["controls_ok"], reading=r["reading"], summary=r["summary"]), ensure_ascii=False, indent=1))
