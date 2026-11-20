"""Runs the analyses S1, S2, K8, E1 – E4 of the environment's analysis plan (ANALYSIS_PLAN.md) and E5, E6, S3 of its §7 (R2) in the main environment, calling the AMICO fit (simrun in envs/reference) as a subprocess."""
from __future__ import annotations

import contextlib
import io
import json
import math
import os
import pathlib
import subprocess
import sys

import numpy as np

from . import metrics as MT
from . import model as M
from ._paths import repo_root, use

use("decomp", "stage1")
import config as C  # noqa: E402
import lib  # noqa: E402

N_FDIR = 4
FDIRS = C.FDIR_E5                                      # (0, 2, 4, 6) — the same four fiber directions as in stages 2 and 3
GENS = C.GENERATORS
CENTRE_NAME = "L0|iii"
L1_PREFIX = "L1|cv0.03|rho0.0|iii|"


# ── Common ────────────────────────────────────────────────────────────────────
def eta17() -> list[tuple[str, tuple, tuple]]:
    """Center (L0|iii) + 16 points of L1 — the same relaxation setting as step-3 X9."""
    return [e for e in C.ETA_E6 if e[0] == CENTRE_NAME or e[0].startswith(L1_PREFIX)]


def centre_eta() -> tuple[tuple, tuple]:
    _, t2, t1 = next(e for e in C.ETA_E6 if e[0] == CENTRE_NAME)
    return t2, t1


def coords_mean(fit: dict, n_fdir: int = N_FDIR) -> dict:
    """Fit values (fiber direction is the fastest axis) → coordinates averaged over fiber directions — as in step-3 load_table, coordinates are computed row by row and then averaged."""
    ndi, fwf, odi = (np.asarray(fit[k], float) for k in ("ndi", "fwf", "odi"))
    xy, _ = lib.coords(ndi, fwf, odi)
    f = np.clip(fwf, lib.CLIP, 1 - lib.CLIP)
    cols = dict(x=xy[:, 0], y=xy[:, 1], odi=odi, raw_ndi=ndi, raw_fiso=fwf, logit_fiso=np.log(f / (1 - f)))
    return {k: v.reshape(-1, n_fdir).mean(1) for k, v in cols.items()}


def l1_draws(t2c, t1c, n: int, seed: int) -> list:
    """L1 between-person relaxation variation (Supplementary S-M2) — two T2 compartments with CV 0.03 · correlation 0 · cut at Mahalanobis radius 2, common tissue-T1 factor with CV 0.04 (±2 SD)."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        while True:
            z = rng.normal(size=2)
            if math.hypot(*z) <= 2.0:
                break
        t2 = (t2c[0] * (1 + 0.03 * z[0]), t2c[1] * (1 + 0.03 * z[1]), t2c[2])
        if t1c is None:
            out.append((t2, None))
            continue
        w = float(np.clip(rng.normal(), -2.0, 2.0))
        out.append((t2, (t1c[0] * (1 + 0.04 * w), t1c[1] * (1 + 0.04 * w), t1c[2])))
    return out


def run_sim(jobs: list[dict], tag: str, workroot) -> dict:
    """Runs simrun in envs/reference — writes a job-bundle json and reads the result npz."""
    wd = pathlib.Path(workroot) / tag
    wd.mkdir(parents=True, exist_ok=True)
    jp, op = wd / "jobs.json", wd / "fits.npz"
    jp.write_text(json.dumps({"jobs": jobs}), encoding="utf-8")
    src = str(pathlib.Path(__file__).resolve().parents[1])
    code = f"import sys; sys.path.insert(0, {src!r}); from rules_env import simrun; simrun.run({str(jp)!r}, {str(op)!r}, {str(wd)!r})"
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", RULES_ENV_REPO=str(repo_root()))
    env.pop("VIRTUAL_ENV", None)
    with open(wd / "simrun.log", "w") as log:
        subprocess.run(["uv", "run", "--frozen", "python", "-I", "-B", "-c", code], cwd=repo_root() / "envs" / "reference",
                       env=env, check=True, stdout=log, stderr=subprocess.STDOUT)
    return dict(np.load(op))


def cell(nu, odi, fiso, t2, t1, gen, fdir) -> dict:
    return dict(nu=float(nu), odi=float(odi), fiso=float(fiso), t2=list(t2), t1=(list(t1) if t1 is not None else None),
                gen=gen, fdir=int(fdir), dpar=1.7)


# ── S2 ──────────────────────────────────────────────────────────────────────
def _combat(dat, batch, group, mean_only):
    import pandas as pd
    from neuroCombat import neuroCombat
    cov = pd.DataFrame({"batch": [str(b) for b in batch]} | ({"group": list(group)} if group is not None else {}))
    with contextlib.redirect_stdout(io.StringIO()):
        return neuroCombat(dat=dat, covars=cov, batch_col="batch", categorical_cols=(["group"] if group is not None else None),
                           mean_only=mean_only)


def combat_apply(dat, batch, est) -> np.ndarray:
    """The same formula as neuroCombatFromTraining(0.2.12) — that function breaks under numpy 2 on int(1-element array), so it is written out directly. Standardization (sample mean of stand.mean + mod.mean · var.pooled) → subtract and divide by the γ* · δ* of the batch → original scale."""
    if np.abs(np.asarray(est["mod.mean"], float)).max() > 1e-12:     # 2026-11-15 after the code review: this study uses it only for estimates trained without covariates
        raise ValueError("combat_apply: covariate-trained estimates are outside its use here (S2 trains on controls without covariates)")
    levels = [str(b) for b in est["batches"]]
    wh = [levels.index(str(b)) for b in batch]
    vp = np.sqrt(np.asarray(est["var.pooled"], float)).reshape(-1, 1)
    sm = (np.asarray(est["stand.mean"], float)[:, 0] + np.asarray(est["mod.mean"], float).mean(axis=1)).reshape(-1, 1)
    b = (np.asarray(dat, float) - sm) / vp
    b = (b - np.asarray(est["gamma.star"], float)[wh, :].T) / np.sqrt(np.asarray(est["delta.star"], float)[wh, :].T)
    return b * vp + sm


def s2_site_diffs(X: dict, n_sites: int, method: str, noise_sd, rng) -> np.ndarray:
    """X[(site, "c"|"p")] = (person, unit) values → after applying the method, the unit mean D_s of the (patient − control) difference for each site."""
    noise_sd = np.asarray(noise_sd, float)
    V = {k: v + rng.normal(0.0, 1.0, v.shape) * noise_sd[None, :] for k, v in X.items()}
    if method in ("none", "ctrl_offset"):
        if method == "ctrl_offset":
            grand = np.mean([V[(s, "c")].mean(0) for s in range(n_sites)], 0)
            for s in range(n_sites):
                off = V[(s, "c")].mean(0) - grand
                V[(s, "c")], V[(s, "p")] = V[(s, "c")] - off, V[(s, "p")] - off
        return np.array([float(np.mean(V[(s, "p")].mean(0) - V[(s, "c")].mean(0))) for s in range(n_sites)])
    if method in ("combat_loc", "combat_ls"):
        cols, batch, grp = [], [], []
        for s in range(n_sites):
            for g in ("c", "p"):
                cols.append(V[(s, g)])
                batch += [s] * len(V[(s, g)])
                grp += [g] * len(V[(s, g)])
        dat = np.concatenate(cols, 0).T
        H = _combat(dat, batch, grp, mean_only=(method == "combat_loc"))["data"].T
        batch, grp = np.array(batch), np.array(grp)
        return np.array([float(np.mean(H[(batch == s) & (grp == "p")].mean(0) - H[(batch == s) & (grp == "c")].mean(0)))
                         for s in range(n_sites)])
    if method == "combat_train_ctrl_loc":
        dc = np.concatenate([V[(s, "c")] for s in range(n_sites)], 0).T
        bc = np.concatenate([[s] * len(V[(s, "c")]) for s in range(n_sites)])
        tr = _combat(dc, bc, None, mean_only=True)
        dp = np.concatenate([V[(s, "p")] for s in range(n_sites)], 0).T
        bp = np.concatenate([[s] * len(V[(s, "p")]) for s in range(n_sites)])
        hp = combat_apply(dp, bp, tr["estimates"]).T
        hc = tr["data"].T
        return np.array([float(np.mean(hp[bp == s].mean(0) - hc[bc == s].mean(0))) for s in range(n_sites)])
    raise ValueError(method)


# ── Real data auxiliary ─────────────────────────────────────────────────────────────
WONG_TES = (62, 72, 82, 92, 102, 112, 122, 132)


def slope_level_noise(sig: float, tes=WONG_TES) -> tuple[float, float]:
    """Noise variance of the OLS slope σ²/Σ(TE − T̄)² for the repeat SD σ (Σ = 4200 ms² for TE 62–132 ms) and noise variance of the TE-mean level σ²/len(TE)."""
    te = np.asarray(tes, float)
    return sig ** 2 / float(((te - te.mean()) ** 2).sum()), sig ** 2 / len(te)


def te_slopes(ser: dict) -> tuple[dict, dict]:
    """{series: {unit: value}} → TE slope (OLS) and 8-TE mean for units that have all 8 TEs."""
    te = np.array(WONG_TES, float)
    units = sorted(set.intersection(*[set(ser[f"TE{t}"]) for t in WONG_TES]))
    k, xbar = {}, {}
    for u in units:
        v = np.array([ser[f"TE{t}"][u] for t in WONG_TES], float)
        k[u] = float(np.polyfit(te, v, 1)[0])
        xbar[u] = float(v.mean())
    return k, xbar


def grouped_loo_offsets(d: dict) -> dict:
    """d[(bundle, unit)] → the median of the other units of the same bundle (excluding itself)."""
    groups = {}
    for (g, u), v in d.items():
        groups.setdefault(g, {})[u] = v
    out = {}
    for g, uv in groups.items():
        for u in uv:
            others = [v for w, v in uv.items() if w != u]
            out[(g, u)] = float(np.median(others)) if others else float("nan")
    return out


# ── Execution part (order of the fixed analysis document, fixed before computing) ─────────────────────────────────────────
Q_CLASS = {"Q1": "TE only", "Q4": "TR only", "Q3": "scheme only", "Q2": "TE, TR and scheme"}


def _dump(out_dir, name, obj) -> pathlib.Path:
    p = pathlib.Path(out_dir) / f"{name}.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, ensure_ascii=False, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)),
                 encoding="utf-8")
    return p


def onh_ctx():
    from . import empirical as EM
    st = EM.onh_states()
    o = EM.read_onh(EM.R.WORK, st)
    s = {sc: EM.s_onh(st["states"], o["rois"], o["units"], sc) for sc in EM.SCALES}
    return EM, st, o, s


def onh_matched_units(EM, st, o) -> dict:
    """Per-unit median (NDI, FISO, ODI) of the OXF1PRI main-analysis subjects → OXF1PRI card · inverse of equation (3) for the center relaxation → true (ν, f, ODI)."""
    from . import cards as K
    subs = st["main_by_unit"]["OXF1PRI"]
    t2c, t1c = centre_eta()
    m = K.card_meta("ONH:OXF1PRI")
    w = M.weights(m["TE"], t2c, m["TR"], t1c)
    out, dropped = {}, []
    for u in o["units"]:
        rows = [o["rois"][f"{sb}_{ses}"]["주"].get(u, {}) for sb, ses in subs.items()]
        v = [[r.get(k) for k in ("NDI", "FWF", "ODI")] for r in rows]
        v = np.array([x for x in v if all(a is not None and math.isfinite(a) for a in x)], float)
        if not len(v):
            dropped.append(u)
            continue
        n, f, od = np.median(v, 0)
        nu, ff = M.invert(n, f, w)
        if 0 < nu < 1 and 0 < ff < 1:
            out[u] = (float(nu), float(ff), float(od))
        else:
            dropped.append(u)
    return dict(units=out, dropped=dropped, n_subjects=len(subs))


def s1a(out_dir, work_dir) -> dict:
    """Stored step-3 E6 fit results — X9 reproduction (stop) → ORS · proportional bias · ρ · rank reversal at the reference tissue 210."""
    use("decomp", "stage1")
    import analyze3 as A3
    EM, st, o, s = onh_ctx()
    s_x = float(np.median([v for v in s["x"][0].values() if v is not None]))
    feat, _ = A3.load_table("E6")
    refs = list(C.refs_e5())
    idx = A3.lsets()["L1|cv0.03|rho0.0|iii"][0]
    ec = A3.ETAN.index("L0|iii")
    stored = json.loads((repo_root() / "decomp" / "results" / "stage3" / "x9.json").read_text(encoding="utf-8"))["X9"]
    reg = {}
    for g in GENS:
        for q in ("Q1", "Q2", "Q3", "Q4"):
            lo, hi = A3.cone_angles(q, idx)
            ink = []
            for r in refs:
                d = feat[q][g][r, ec, :3] - feat["A"][g][r, ec, :3]
                a = math.degrees(math.atan2(d[1], d[0])) % 180
                ink.append(lo - 1 <= a <= hi + 1 if lo <= hi else True)
            reg[f"{g}|{q}"] = dict(recomputed=float(np.mean(ink)), stored=stored[f"{g}|{q}"]["xy_dir_in_fraction_cone"])
    reg_ok = all(abs(v["recomputed"] - v["stored"]) < 5e-5 for v in reg.values())
    res = dict(x9_regression=reg, x9_regression_ok=reg_ok, s_x_median=s_x, s_x_df=s["x"][1], c=1.96, n_units=len(refs),
               eta=[A3.ETAN[e] for e in idx], settings={})
    if not reg_ok:
        _dump(out_dir, "s1a", res)
        raise SystemExit("⛔ X9 재현 실패 — 멈춤")
    rc = 1.96 * math.sqrt(2) * s_x
    for g in GENS:
        for q in ("Q1", "Q2", "Q3", "Q4"):
            for e in idx:
                xa, xq = feat["A"][g][refs, e, 0], feat[q][g][refs, e, 0]
                d = {str(r): float(v) for r, v in zip(refs, xq - xa)}
                r_ = MT.ors(d, {k: s_x for k in d}, 1.96)
                res["settings"][f"{g}|{q}|{A3.ETAN[e]}"] = dict(ors=r_["share"], retest_level=r_["retest_level"], offset=r_["offset"],
                                                                prop_bias=MT.prop_bias(xa, xq), spearman=MT.spearman(xa, xq),
                                                                rank_flip=MT.rank_flip_share(xa, xq, rc))
    res["primary"] = {f"{g}|{q} ({Q_CLASS[q]})": float(np.mean([res["settings"][f"{g}|{q}|{A3.ETAN[e]}"]["retest_level"] is True for e in idx]))
                      for g in GENS for q in ("Q1", "Q2", "Q3", "Q4")}
    _dump(out_dir, "s1a", res)
    return res


def _wong_matched(EM, w) -> dict:
    t2c, t1c = centre_eta()
    wt = M.weights(62.0, t2c, 5800.0, t1c)
    out, dropped = {}, []
    for sb, ser in sorted(w["roi"].items()):
        for u, (n, f, od) in EM._rows(ser["TE62"], "주"):
            nu, ff = M.invert(n, f, wt)
            if 0 < nu < 1 and 0 < ff < 1:
                out[(sb, u)] = (float(nu), float(ff), float(od))
            else:
                dropped.append(f"{sb}|{u}")
    return dict(units=out, dropped=dropped)


def _cells_units_eta(units: dict, etas) -> list[dict]:
    """Order: unit × η × generator × fiber direction (fastest axis)."""
    return [cell(nu, od, ff, t2, t1, g, fd) for (nu, ff, od) in units.values() for _, t2, t1 in etas for g in GENS for fd in FDIRS]


def s1b(out_dir, work_dir) -> dict:
    """ROI-matched units — τ_fit · β_fit of the 8 TEs of Wong, ORS of the 28 pairs of ON-Harmony."""
    from . import cards as K
    EM, st, o, s = onh_ctx()
    w = EM.read_wong()
    wm, om = _wong_matched(EM, w), onh_matched_units(EM, st, o)
    etas = eta17()
    jobs = [dict(key=f"WONG:TE{t}", card=f"WONG:TE{t}", cells=_cells_units_eta(wm["units"], etas)) for t in WONG_TES]
    jobs += [dict(key=f"ONH:{u}", card=f"ONH:{u}", cells=_cells_units_eta(om["units"], etas)) for u in K.ONH_UNITS]
    fits = run_sim(jobs, "s1b", work_dir)

    def table(key, n_units):
        cm = coords_mean({k: fits[f"{key}|{k}"] for k in ("ndi", "fwf", "odi")})
        return {k: v.reshape(n_units, len(etas), len(GENS)) for k, v in cm.items()}

    wk = list(wm["units"])
    WT = {t: table(f"WONG:TE{t}", len(wk)) for t in WONG_TES}
    res = dict(wong=dict(n_units=len(wk), dropped=wm["dropped"], by={}), onh=dict(n_units=len(om["units"]), dropped=om["dropped"], by={}),
               eta=[e[0] for e in etas])
    te = np.array(WONG_TES, float)
    for co in ("x", "y"):
        for gi, g in enumerate(GENS):
            for ei, (en, _, _) in enumerate(etas):
                V = np.stack([WT[t][co][:, ei, gi] for t in WONG_TES], 1)            # (ROI units, TE)
                k = np.polyfit(te, V.T, 1)[0]
                xbar = V.mean(1)
                res["wong"]["by"][f"{co}|{g}|{en}"] = dict(tau_fit=float(np.std(k, ddof=1)), beta_fit=MT.beta_level(k, xbar, 0.0),
                                                          mean_slope=float(k.mean()))
    ou = list(om["units"])
    OT = {u: table(f"ONH:{u}", len(ou)) for u in K.ONH_UNITS}
    s_x = s["x"][0]
    for a, b in K.pairs("onh28"):
        ua, ub = a.split(":")[1], b.split(":")[1]
        for gi, g in enumerate(GENS):
            for ei, (en, _, _) in enumerate(etas):
                xa, xb = OT[ua]["x"][:, ei, gi], OT[ub]["x"][:, ei, gi]
                d = {u: float(v) for u, v in zip(ou, xb - xa)}
                r_ = MT.ors(d, s_x, 1.96)
                res["onh"]["by"][f"{ua}~{ub}|{g}|{en}"] = dict(ors=r_["share"], retest_level=r_["retest_level"], n=r_["n"],
                                                               prop_bias=MT.prop_bias(xa, xb), spearman=MT.spearman(xa, xb))
    _dump(out_dir, "s1b", res)
    return res


S2_GROUPS = {   # Name: (T2 center, T1 center, true ν difference)
    "fig3_c": ((82.5, 52.0, 1000.0), None, 0.0),
    "fig3_p": ((82.5, 70.0, 1000.0), None, -0.03),
    "cen_c": (None, None, 0.0),                                     # Center relaxation (filled in below)
    "l2en_p": ("l2en", None, -0.03),
    "l2both_p": ("l2both", None, -0.03),
    "frac_p": (None, None, -0.03),
}
S2_SCEN = {"fig3": ("fig3_c", "fig3_p"), "l2en": ("cen_c", "l2en_p"), "l2both": ("cen_c", "l2both_p"), "frac": ("cen_c", "frac_p")}
S2_METHODS = ("none", "ctrl_offset", "combat_loc", "combat_train_ctrl_loc", "combat_ls")
S2_LOCATION = ("none", "ctrl_offset", "combat_loc", "combat_train_ctrl_loc")


def s2_groups() -> dict:
    t2c, t1c = centre_eta()
    g = {}
    for name, (t2, t1, dnu) in S2_GROUPS.items():
        if t2 is None:
            t2, t1 = t2c, t1c
        elif t2 == "l2en":
            t2, t1 = (t2c[0], t2c[1] * 1.0984, t2c[2]), (t1c[0], t1c[1] * 1.0441, t1c[2])
        elif t2 == "l2both":
            t2, t1 = (t2c[0] * 1.0984, t2c[1] * 1.0984, t2c[2]), (t1c[0] * 1.0441, t1c[1] * 1.0441, t1c[2])
        g[name] = dict(t2=tuple(t2), t1=(tuple(t1) if t1 is not None else None), dnu=dnu)
    return g


def _s2_combo(args):
    X, n_sites, method, sd, seed, n_mc = args
    rng = np.random.default_rng(seed)
    return np.array([s2_site_diffs(X, n_sites, method, sd, rng) for _ in range(n_mc)])


def s2_subjects(mode: str, gd: dict, n_subj: int, seed: int) -> list:
    """Population mode (pop): one person with the group-center relaxation (the setting of Proposition 4) · L1 mode (l1): L1 variation for every person (the original design of the document fixed before computing)."""
    if mode == "pop":
        return [(tuple(gd["t2"]), (tuple(gd["t1"]) if gd["t1"] is not None else None))]
    return l1_draws(gd["t2"], gd["t1"], n_subj, seed)


def s2_sample_algebra(subj: dict, sites: dict) -> dict:
    """Algebraic group difference computed with the relaxation of the same people (the relaxation part) — subj[(site, "c"|"p")] = [(T2, T1)] · sites[site] = {"TE","TR"}. D_s = mean over patients of x_shift − mean over controls of x_shift (the true fraction difference cancels between sites, so it is subtracted)."""
    out = {}
    for s_, site in sites.items():
        m = {g: float(np.mean([M.x_shift(site, {"t2": t2, "t1": t1}) for t2, t1 in subj[(s_, g)]])) for g in ("c", "p")}
        out[s_] = m["p"] - m["c"]
    return out


def s2(out_dir, work_dir, n_subj: int = 20, n_mc: int = 200) -> dict:
    """Proposition 4 · fit results — two and eight sites, logit · raw, five methods, Monte Carlo ROI noise. Fix before computing (2026-11-14 · ANALYSIS_PLAN fix 1): the main analysis is the population mode (all simulated participants have the group-center relaxation — the setting of Proposition 4), and the original design (L1 variation) is auxiliary — its comparison target is the algebraic value computed with the relaxation of the same people (the error of drawing people cancels)."""
    from concurrent.futures import ProcessPoolExecutor
    from . import cards as K
    EM, st, o, s = onh_ctx()
    om = onh_matched_units(EM, st, o)
    units = [u for u in om["units"] if s["x"][0].get(u) is not None and s["raw_ndi"][0].get(u) is not None]
    groups = s2_groups()
    sitesets = {"two": [dict(key="S2:TE70", card="E5:A", TE=70.0), dict(key="S2:TE100", card="E5:A", TE=100.0)],
                "onh8": [dict(key=f"ONH:{u}", card=f"ONH:{u}") for u in K.ONH_UNITS]}
    site_te_tr = {"S2:TE70": (70.0, 3600.0), "S2:TE100": (100.0, 3600.0)}
    for u in K.ONH_UNITS:
        m = K.card_meta(f"ONH:{u}")
        site_te_tr[f"ONH:{u}"] = (m["TE"], m["TR"])
    all_sites = [x for ss in sitesets.values() for x in ss]
    subj, jobs = {}, []
    for mode in ("pop", "l1"):
        for si, site in enumerate(all_sites):
            cells = []
            for gi, (gname, gd) in enumerate(groups.items()):
                subj[(mode, site["key"], gname)] = s2_subjects(mode, gd, n_subj, seed=20261115 + 1000 * si + gi)
                for t2, t1 in subj[(mode, site["key"], gname)]:
                    for u in units:
                        nu, ff, od = om["units"][u]
                        for g in GENS:
                            for fd in FDIRS:
                                cells.append(cell(nu + gd["dnu"], od, ff, t2, t1, g, fd))
            jobs.append(dict(site, key=f"{site['key']}#{mode}", cells=cells))
    fits = run_sim(jobs, "s2", work_dir)
    nU, nG = len(units), len(GENS)
    vals = {}
    for mode in ("pop", "l1"):
        ns = 1 if mode == "pop" else n_subj
        for site in all_sites:
            key = f"{site['key']}#{mode}"
            cm = coords_mean({k: fits[f"{key}|{k}"] for k in ("ndi", "fwf", "odi")})
            for sc in ("x", "raw_ndi"):
                arr = cm[sc].reshape(len(groups), ns, nU, nG)
                for gi, gname in enumerate(groups):
                    for g_i, g in enumerate(GENS):
                        v = arr[gi, :, :, g_i]
                        vals[(mode, sc, site["key"], gname, g)] = np.repeat(v, n_subj, 0) if mode == "pop" else v
    combos, keys = [], []
    for mode in ("pop", "l1"):
        for sc in ("x", "raw_ndi"):
            sd = np.array([s[sc][0][u] for u in units])
            for ssn, ss in sitesets.items():
                for scen, (cg, pg) in S2_SCEN.items():
                    for g in GENS:
                        X = {}
                        for i, site in enumerate(ss):
                            X[(i, "c")] = vals[(mode, sc, site["key"], cg, g)]
                            X[(i, "p")] = vals[(mode, sc, site["key"], pg, g)]
                        for mth in S2_METHODS:
                            keys.append((mode, sc, ssn, scen, g, mth))
                            combos.append((X, len(ss), mth, sd, 20261116 + len(keys), n_mc))
    with ProcessPoolExecutor(max_workers=48) as ex:
        reps = list(ex.map(_s2_combo, combos))
    res = dict(n_units=nU, units_dropped=om["dropped"], n_subjects=n_subj, n_mc=n_mc, groups=groups, by={},
               amendment="계산 전 고침 1(2026-11-14): 주 분석 = 모집단 모드 · L1 모드는 같은 사람들의 대수 값과 비교")
    for (mode, sc, ssn, scen, g, mth), D in zip(keys, reps):
        cg, pg = S2_SCEN[scen]
        sites = [x["key"] for x in sitesets[ssn]]
        tes = np.array([site_te_tr[k][0] for k in sites])
        sdict = {k: {"TE": site_te_tr[k][0], "TR": site_te_tr[k][1]} for k in sites}
        sub = {}
        for k in sites:
            sub[(k, "c")] = subj[(mode, k, cg)]
            sub[(k, "p")] = subj[(mode, k, pg)]
        alg_d = s2_sample_algebra(sub, sdict)
        alg = np.array([alg_d[k] for k in sites])
        r = dict(D_mean=D.mean(0))
        if ssn == "two":
            I = D[:, 1] - D[:, 0]
            I_alg = float(alg[1] - alg[0])
            eq7 = M.prop4_change({"TE": 100.0, "TR": 3600.0}, {"TE": 70.0, "TR": 3600.0}, groups[cg], groups[pg])
            r.update(interaction=float(I.mean()), q025=float(np.percentile(I, 2.5)), q975=float(np.percentile(I, 97.5)), eq7=eq7,
                     algebra_same_people=I_alg, fit_minus_algebra=float(I.mean() - I_alg),
                     ratio=(float(I.mean() / eq7) if mode == "pop" and sc == "x" and abs(eq7) > 1e-12 else None))
        else:
            sl = np.array([np.polyfit(tes, d, 1)[0] for d in D])
            res_sd = np.array([np.std(d - np.polyval(np.polyfit(tes, d, 1), tes), ddof=2) for d in D])
            a_fit = np.polyfit(tes, alg, 1)
            r.update(slope=float(sl.mean()), slope_q025=float(np.percentile(sl, 2.5)), slope_q975=float(np.percentile(sl, 97.5)),
                     slope_alg=float(a_fit[0]), resid_sd=float(res_sd.mean()),
                     resid_sd_alg=float(np.std(alg - np.polyval(a_fit, tes), ddof=2)), te=tes, alg=alg)
        res["by"][f"{mode}|{sc}|{ssn}|{scen}|{g}|{mth}"] = r
    rat = [v["ratio"] for k, v in res["by"].items() if k.startswith("pop|x|two|") and k.split("|")[3] in ("fig3", "l2en", "l2both")
           and k.split("|")[5] in S2_LOCATION and v["ratio"] is not None]
    res["primary"] = dict(ratio_min=float(min(rat)), ratio_max=float(max(rat)), n=len(rat), reading=MT.s2_reading(rat))
    _dump(out_dir, "s2", res)
    return res


def k8(out_dir, work_dir) -> dict:
    """FISO sign rule (fixed/fiso_sign) — noise-free fit values (true FISO 0.08 · ODI 0.20) instead of the algebraic Δlogit FISO."""
    from . import cards as K
    use("decomp", "fixed")
    import delta_tau as DT
    import fiso_sign as FS
    EM, st, o, s = onh_ctx()
    grids = FS.grids()
    gl = [(gn, t2, t1) for gn, grid in grids.items() for t2, t1 in grid]
    nus = list(FS.NU)
    jobs = [dict(key=f"ONH:{u}", card=f"ONH:{u}",
                 cells=[cell(nu, 0.20, 0.08, t2, t1, g, fd) for _, t2, t1 in gl for nu in nus for g in GENS for fd in FDIRS])
            for u in K.ONH_UNITS]
    fits = run_sim(jobs, "k8", work_dir)
    LF = {}
    for u in K.ONH_UNITS:
        cm = coords_mean({k: fits[f"ONH:{u}|{k}"] for k in ("ndi", "fwf", "odi")})["logit_fiso"].reshape(len(gl), len(nus), len(GENS))
        m = K.card_meta(f"ONH:{u}")
        LF[(round(m["TE"], 6), round(m["TR"], 6))] = cm
    gidx = {gn: [i for i, (n_, _, _) in enumerate(gl) if n_ == gn] for gn in grids}
    cnt = o["cnt"]
    orig = FS.fiso_sign
    alg = DT.design_signs(cnt)
    out = dict(algebra=alg.tolist(), fitted={})
    for gi_, g in enumerate(GENS):
        def fitted_sign(teA, trA, teB, trB, nu=None, _g=gi_):
            a, b = LF[(round(teA, 6), round(trA, 6))], LF[(round(teB, 6), round(trB, 6))]
            res_ = {}
            for gn, ii in gidx.items():
                v = (b[ii, :, _g] - a[ii, :, _g]).ravel()
                res_[gn] = dict(sign=FS.sign_of(v), max_abs=float(np.abs(v).max()))
            return res_
        FS.fiso_sign = fitted_sign
        try:
            out["fitted"][g] = DT.design_signs(cnt).tolist()
        finally:
            FS.fiso_sign = orig
    fv = EM.JO.device_values(st["main_by_unit"], o["rois"], "FWF", "logit", "주")
    fs_ = EM.JO.pair_values(fv, o["units"])
    keep_f, _ = EM.JO.keep_pairs(cnt, fs_)
    obs = np.array([np.sign(v["value"]) if v["value"] is not None else np.nan for v in fs_])
    pred_pairs = [k for k in keep_f if alg[k] != 0]
    out["observed"] = obs.tolist()
    out["keep"] = keep_f
    out["algebra_agree_observed"] = int(sum(obs[k] == alg[k] for k in pred_pairs))
    out["n_algebra_predicted"] = len(pred_pairs)
    out["regression_ok"] = (out["algebra_agree_observed"], out["n_algebra_predicted"]) == (2, 8)
    for g in GENS:
        f = np.array(out["fitted"][g])
        fp = [k for k in keep_f if f[k] != 0]
        out[f"{g}|n_fitted_predicted"] = len(fp)
        out[f"{g}|fitted_agree_observed"] = int(sum(obs[k] == f[k] for k in fp))
        out[f"{g}|fitted_agree_algebra_on_algebra_pairs"] = int(sum(f[k] == alg[k] for k in pred_pairs))
    out["pairs"] = [EM.JO.pair_name(k) for k in range(len(fs_))]
    _dump(out_dir, "k8", out)
    return out


# ── Real data E1 – E4 ──────────────────────────────────────────────────────────
def _boot_beta(kv, xv, var_x, B=2000, seed=0):
    rng = np.random.default_rng(seed)
    n = len(kv)
    reps = []
    for _ in range(B):
        i = rng.integers(0, n, n)
        reps.append(MT.beta_level(kv[i], xv[i], var_x))
    return float(np.percentile(reps, 2.5)), float(np.percentile(reps, 97.5))


def e1(out_dir, work_dir) -> dict:
    from . import empirical as EM
    w = EM.read_wong()
    fiso = EM.wong_series(w, "raw_fiso")
    s1b_path = pathlib.Path(out_dir) / "s1b.json"
    s1b_ = json.loads(s1b_path.read_text(encoding="utf-8")) if s1b_path.exists() else None
    te = np.array(WONG_TES, float)
    res = dict(by={}, pooled={})
    for co in ("x", "y"):
        ser_all = EM.wong_series(w, co)
        kc, su, cl, vks, dfs = [], [], [], [], []
        for pi, sb in enumerate(sorted(ser_all)):
            ser = ser_all[sb]
            k, xbar = te_slopes(ser)
            us = sorted(k)
            kv, xv = np.array([k[u] for u in us]), np.array([xbar[u] for u in us])
            sig, df = EM.s_rt_corrected(ser)
            var_k, var_x = slope_level_noise(sig)
            tau, beta = MT.tau_hat(kv, var_k), MT.beta_level(kv, xv, var_x)
            pt, lo, hi = MT.te_tolerance_interval(np.full(len(us), sig), kv, var_k, df, clusters=np.arange(len(us)), B=2000, seed=0)
            d = {u: ser["TE132"][u] - ser["TE62"][u] for u in us if u in ser["TE132"] and u in ser["TE62"]}
            obs = MT.ors(d, {u: sig for u in d}, 1.96)
            fl = np.array([np.mean([fiso[sb][f"TE{t}"][u] for t in WONG_TES]) for u in us])
            kr = kv - np.polyval(np.polyfit(fl, kv, 1), fl)
            V = np.array([[ser[f"TE{t}"][u] for t in WONG_TES] for u in us])
            resid = V - np.polyval(np.polyfit(te, V.T, 1), te[:, None]).T
            res["by"][f"{co}|{sb}"] = dict(
                n_units=len(us), s_rt_corrected=sig, df=df, tau=tau, delta_r2_sd_per_s=tau * 1000.0, beta=beta,
                beta_ci95=_boot_beta(kv, xv, var_x), dte_star=pt, dte_ci90=[lo, hi], mean_slope=float(kv.mean()),
                ors_62_132_observed=obs["share"], ors_62_132_predicted=MT._te_share(70.0, np.full(len(us), sig), tau, 1.96),
                fiso_partial=dict(tau=MT.tau_hat(kr, var_k), beta=MT.beta_level(kr, xv, var_x)),
                resid_sigma=float(math.sqrt(np.sum(resid ** 2) / (resid.size - 2 * len(us)))))
            kc += list(kv - kv.mean())
            su += [sig] * len(us)
            cl += [pi] * len(us)
            vks.append((var_k, len(us)))
            dfs.append(df)
        var_k_pool = sum(v * n for v, n in vks) / sum(n for _, n in vks)
        # Version that was run (2026-11-14): finite-n factor 1 + 1/80 with 81 ROI units in the pooled array — differs from n (batch 27) in plan §1 (code review F2). The calculation with the error fixed is e1_fix (§8).
        pt, lo, hi = MT.te_tolerance_interval(np.array(su), np.array(kc), var_k_pool, sum(dfs), clusters=np.array(cl), B=2000, seed=0)
        tau_pool = MT.tau_hat(np.array(kc), var_k_pool)
        res["pooled"][co] = dict(tau=tau_pool, delta_r2_sd_per_s=tau_pool * 1000.0, dte_star=pt, dte_ci90=[lo, hi], df=sum(dfs),
                                 n_units=len(kc))
        if s1b_:
            res["pooled"][co]["s1b_T-in_centre"] = s1b_["wong"]["by"].get(f"{co}|T-in|L0|iii")
    res["primary"] = dict(dte_star_x=res["pooled"]["x"]["dte_star"], ci90=res["pooled"]["x"]["dte_ci90"])
    _dump(out_dir, "e1", res)
    return res


def _onh_scales(EM, st, o, s):
    used = [f"{sb}_{ses}" for subs in st["main_by_unit"].values() for sb, ses in subs.items()]
    excl = set().union(*[EM.clipped(o["rois"][k]) for k in used])
    keep = EM.onh_pairs(st, o, "x")["keep"]
    out = {}
    for sc in EM.SCALES:
        P = EM.onh_pairs(st, o, sc)
        for k in keep:
            bu = {u: v for u, v in P["summ"][k]["by_unit"].items()
                  if v is not None and (sc not in ("x", "y", "logit_fiso") or u not in excl)}
            out[(sc, P["names"][k])] = dict(by_unit=bu, raw=MT.ors(bu, s[sc][0], 1.96, loo=False), ors=MT.ors(bu, s[sc][0], 1.96))
    return out, keep, sorted(excl)


def e2(out_dir, work_dir) -> dict:
    from scipy import stats
    EM, st, o, s = onh_ctx()
    tab, keep, excl = _onh_scales(EM, st, o, s)
    names = [EM.JO.pair_name(k) for k in keep]
    res = dict(pairs=names, excluded_clipped_units=excl, s_df={sc: s[sc][1] for sc in EM.SCALES}, by={})
    for (sc, p), v in tab.items():
        res["by"][f"{sc}|{p}"] = dict(raw=v["raw"]["share"], ors=v["ors"]["share"], n=v["ors"]["n"], retest_level_ors=v["ors"]["retest_level"])
    dte = np.array([abs(o["rows_by_pair"][p]["dTE"]) for p in names])
    ors_x = np.array([tab[("x", p)]["ors"]["share"] for p in names])
    t = stats.kendalltau(dte, ors_x, variant="b")
    res["tau_b_dTE_orsx"] = dict(tau=float(t.statistic), p=float(t.pvalue))
    e1p = pathlib.Path(out_dir) / "e1.json"
    if e1p.exists():
        tau_pool = json.loads(e1p.read_text(encoding="utf-8"))["pooled"]["x"]["tau"]
        sx = np.array([v for v in s["x"][0].values() if v is not None])
        res["predicted_from_e1"] = {p: dict(dTE=float(d), predicted=MT._te_share(float(d), sx, tau_pool, 1.96), observed=float(x))
                                    for p, d, x in zip(names, dte, ors_x)}
    _dump(out_dir, "e2", res)
    return res


def _tong_scales(EM):
    import itertools
    t = EM.read_tong()
    out = {"_clipped_unit_sessions": int(sum(len(EM.clipped(r)) for r in t["roi"].values()))}   # Not dropped under the fixed analysis document; only the count is reported
    for sc in EM.SCALES:
        values, rep = EM.tong_structs(t, sc)
        units = sorted({u for c in values for sb in values[c] for u in values[c][sb]})
        sd = EM.retest_sd_scalar(rep, units)
        d, thr = {}, {}
        for a, b in itertools.combinations(EM.JT.CENTERS, 2):
            for sb in sorted({x for c in (a, b) for x in values.get(c, {})}):
                va, vb = values.get(a, {}).get(sb, {}), values.get(b, {}).get(sb, {})
                for u in set(va) & set(vb):
                    su, df = sd[u]
                    if su is None:
                        continue
                    d[((sb, a, b), u)] = vb[u] - va[u]
                    thr[((sb, a, b), u)] = EM.JT.c_of(df) * math.sqrt(2) * su
        off = grouped_loo_offsets(d)
        out[sc] = dict(n=len(d), raw=float(np.mean([abs(d[k]) <= thr[k] for k in d])),
                       ors=float(np.mean([abs(d[k] - off[k]) <= thr[k] for k in d])), d=d, thr=thr, off=off)
    return out


def e3(out_dir, work_dir) -> dict:
    from . import empirical as EM
    tb = _tong_scales(EM)
    res = dict(by={sc: dict(n=v["n"], raw=v["raw"], ors=v["ors"]) for sc, v in tb.items() if sc in EM.SCALES},
               clipped_unit_sessions=tb["_clipped_unit_sessions"])
    res["regression_ok"] = (round(res["by"]["x"]["raw"], 3), round(res["by"]["odi"]["raw"], 3), round(res["by"]["y"]["raw"], 3)) == (0.931, 0.934, 0.844)
    _dump(out_dir, "e3", res)
    return res


def _wong_pair_scales(EM, w, ta: int, tb: int, scales) -> dict:
    """One TE pair of Wong — for each scale, the difference between units, the value after subtracting the within-participant LOO offset, and the repeat SD after subtracting the mean difference."""
    out = {}
    for sc in scales:
        ser = EM.wong_series(w, sc)
        d, s_ = {}, {}
        for sb, se in ser.items():
            sig, _ = EM.s_rt_corrected(se)
            for u in set(se[f"TE{ta}"]) & set(se[f"TE{tb}"]):
                d[(sb, u)] = se[f"TE{tb}"][u] - se[f"TE{ta}"][u]
                s_[(sb, u)] = sig
        off = grouped_loo_offsets(d)
        out[sc] = dict(dc={k: v - off[k] for k, v in d.items()}, s=s_)
    return out


def e4(out_dir, work_dir) -> dict:
    import itertools
    from . import empirical as EM
    w = EM.read_wong()
    res = dict(wong={}, onh={}, tong={})
    pairs_sc = (("x", "raw_ndi"), ("y", "logit_fiso"))
    for ta, tb in itertools.combinations(WONG_TES, 2):
        P = _wong_pair_scales(EM, w, ta, tb, EM.SCALES)
        row = {sc: MT.ors(P[sc]["dc"], P[sc]["s"], 1.96, loo=False)["share"] for sc in EM.SCALES}
        for a, b in pairs_sc:
            cl = {k: k[0] for k in P[a]["dc"]}
            row[f"{a}-{b}"] = MT.ors_diff(P[a]["dc"], P[a]["s"], P[b]["dc"], P[b]["s"], 1.96, cl, B=2000, seed=0, loo=False)
        res["wong"][f"TE{ta}-TE{tb}"] = row
    sm = {}
    for sb, se in EM.wong_series(w, "raw_ndi").items():
        us = sorted(set(se["TE62"]) & set(se["TE132"]))
        r = MT.scale_models([se["TE62"][u] for u in us], [se["TE132"][u] - se["TE62"][u] for u in us])
        sm[sb] = dict(a=r["a"], c=r["c"], ss_raw=r["ss_raw"], ss_logit=r["ss_logit"], n=len(us))
    res["wong_scale_models_62_132"] = sm
    EMo, st, o, s = onh_ctx()
    tab, keep, excl = _onh_scales(EMo, st, o, s)
    names = [EMo.JO.pair_name(k) for k in keep]
    for p in names:
        row = {sc: tab[(sc, p)]["ors"]["share"] for sc in EM.SCALES}
        for a, b in pairs_sc:
            row[f"{a}-{b}"] = MT.ors_diff(tab[(a, p)]["by_unit"], s[a][0], tab[(b, p)]["by_unit"], s[b][0], 1.96, {}, B=2000, seed=0)
        res["onh"][p] = row
    res["onh_median"] = {sc: float(np.median([res["onh"][p][sc] for p in names])) for sc in EM.SCALES}
    res["onh_excluded_clipped_units"] = excl
    tb_ = _tong_scales(EM)
    res["tong"] = {sc: dict(raw=v["raw"], ors=v["ors"], n=v["n"]) for sc, v in tb_.items() if sc in EM.SCALES}
    res["wong_clipped_units_62_132"] = {sb: sorted(EM.clipped(ser["TE62"]) | EM.clipped(ser["TE132"])) for sb, ser in w["roi"].items()}
    pr = res["wong"]["TE62-TE132"]
    res["primary"] = dict(wong_62_132_x_minus_raw=pr["x-raw_ndi"], wong_62_132_y_minus_logitfiso=pr["y-logit_fiso"], onh_median=res["onh_median"])
    _dump(out_dir, "e4", res)
    return res


# ── R2: analyses defined after seeing E4 (ANALYSIS_PLAN §7 · all exploratory) ─────────────────────────────
def _expit(z):
    return 1.0 / (1.0 + np.exp(-np.asarray(z, float)))


def loo_median(V: np.ndarray) -> np.ndarray:
    """V (..., n) → along the last axis, the median of the rest with each element leaving itself out (the same value as grouped_loo_offsets · vectorized)."""
    V = np.asarray(V, float)
    n = V.shape[-1]
    if n < 2:
        return np.full(V.shape, np.nan)
    order = np.argsort(V, axis=-1, kind="stable")
    S = np.take_along_axis(V, order, -1)
    rank = np.empty_like(order)
    np.put_along_axis(rank, order, np.arange(n) * np.ones(order.shape, int), -1)
    m = n - 1

    def t(j):                                          # The j-th of the rest (sorted) = s_j (j < r) · s_{j+1} (j ≥ r)
        a = np.take(S, j, axis=-1)[..., None]
        b = np.take(S, j + 1, axis=-1)[..., None] if j + 1 < n else a
        return np.where(j < rank, a, b)
    if m % 2:
        return t((m - 1) // 2)
    return 0.5 * (t(m // 2 - 1) + t(m // 2))


def group_share(D: dict, s: dict, c: float = 1.96) -> np.ndarray:
    """D[bundle] = differences of shape (..., n) · s[bundle] = threshold SD → retest ratio after subtracting the within-bundle LOO median (the rule of E4 · the ratio pooled over bundles, for each leading axis)."""
    hit, tot = 0.0, 0
    for g, V in D.items():
        V = np.asarray(V, float)
        hit = hit + (np.abs(V - loo_median(V)) <= c * math.sqrt(2.0) * s[g]).sum(-1)
        tot += V.shape[-1]
    return np.asarray(hit, float) / tot


def e5_predict(par: dict, M_: int = 2000, seed: int = 0, use_beta: bool = True, c: float = 1.96, half: float = 35.0) -> dict:
    """E5 expected ratio — for each participant k_u = k̄ + β(x̄_u − mean) + e_u (e ~ N(0, τ_e²) · τ_e² = τ̂² − β²·(var x̄ − σ²/8)), true x(62) = x̄_u − 35 k_u · x(132) = x̄_u + 35 k_u, measurement noise N(0, σ²) added to x → raw = expit(x). The threshold is the same per-scale s as in E4."""
    rng = np.random.default_rng(seed)
    Dx, Dr, sx, sr = {}, {}, {}, {}
    for sb in sorted(par):
        p = par[sb]
        xb = np.asarray(p["xbar"], float)
        n = len(xb)
        b = float(p["beta"]) if use_beta else 0.0
        te = e5_tau_e(float(p["tau"]), b, xb, float(p["var_x"]))
        k = float(p["kbar"]) + b * (xb - xb.mean()) + rng.normal(0.0, 1.0, (M_, n)) * te
        x62 = xb - half * k + rng.normal(0.0, float(p["sigma"]), (M_, n))
        x132 = xb + half * k + rng.normal(0.0, float(p["sigma"]), (M_, n))
        Dx[sb], Dr[sb] = x132 - x62, _expit(x132) - _expit(x62)
        sx[sb], sr[sb] = float(p["s_x"]), float(p["s_raw"])
    ox, orr = group_share(Dx, sx, c), group_share(Dr, sr, c)
    d = ox - orr
    return dict(ors_x=float(ox.mean()), ors_raw=float(orr.mean()), diff=float(d.mean()),
                diff_q025=float(np.percentile(d, 2.5)), diff_q975=float(np.percentile(d, 97.5)), M=M_, seed=seed, use_beta=use_beta)


def e5_tau_e(tau: float, beta: float, xbar, var_x: float) -> float:
    """§7: slope SD not explained by the level — τ_e² = max(τ̂² − β²·(var x̄ − σ²/8), 0)."""
    vtrue = max(float(np.var(np.asarray(xbar, float), ddof=1)) - var_x, 0.0)
    return math.sqrt(max(tau ** 2 - beta ** 2 * vtrue, 0.0))


def gate_close(a: float, b: float, tol: float = 1e-9) -> bool:
    """Comparison of the stop condition (recomputed value = stored value)."""
    return a is not None and b is not None and abs(float(a) - float(b)) < tol


E6_INPUT_SHA256 = {"fits.npz": "f0a9c634cacc7b3b171adadcb4e2da2dc7f809dd652f9452094a2485567db384",
                   "jobs.json": "f9e75337c643e025aaaa6d2b0eed2a920cbcdb3162c91cb614d9edd90607475c"}   # ANALYSIS_PLAN §7


def e6_inputs_ok(wd, expected: dict = E6_INPUT_SHA256) -> bool:
    """Whether the sha256 of the E6 input (the stored S1b fit values) equals the value in plan §7 (added after code review F6 of 2026-11-15 · agreement confirmed in a post hoc comparison)."""
    import hashlib
    return all(hashlib.sha256((pathlib.Path(wd) / n).read_bytes()).hexdigest() == h for n, h in expected.items())


def e6_cells_match(got: list, want: list) -> bool:
    """Whether the stored job list and the cells of the rebuilt ROI-matched units agree in all cells (including t2 · t1 · dpar)."""
    return len(got) == len(want) and all(a == b for a, b in zip(got, want))


def e5_first_order(kbar: float, beta: float, nu_mid_mean: float) -> dict:
    """First order: relative rate of change of the logit shift R = β/k̄ (per unit of x) · relative rate of change of the raw shift R − (2ν̄ − 1). The scale with the smaller value is the more uniform one."""
    R, T = beta / kbar, 2.0 * nu_mid_mean - 1.0
    return dict(R=R, T=T, logit_rate=R, raw_rate=R - T, raw_more_uniform=bool(abs(R - T) < abs(R)))


def e5_reading(pred_diff: float, lo: float, hi: float) -> str:
    """(a) The predicted difference of the gradient structure lies inside the 95 % interval of E4 · (b) outside it."""
    return "가" if lo <= pred_diff <= hi else "나"


def e6_reading(ors_x: float, ors_raw: float, ratio: float = 0.90) -> str:
    """(a) x passes · raw fails — E4 had the power to tell the scales apart · (b) both pass or both fail · (c) raw passes · x fails."""
    px, pr = ors_x >= ratio, ors_raw >= ratio
    return "가" if (px and not pr) else ("다" if (pr and not px) else "나")


def e6_predict(x62: dict, x132: dict, sigma: dict, s_x: dict, s_raw: dict, M_: int = 2000, seed: int = 0, c: float = 1.96) -> dict:
    """E6 expected ratio — measurement noise N(0, σ²) added to the noise-free fit x (the unit array for each participant) → raw = expit(x) · the rule and threshold of E4."""
    rng = np.random.default_rng(seed)
    Dx, Dr = {}, {}
    for sb in sorted(x62):
        a, b = np.asarray(x62[sb], float), np.asarray(x132[sb], float)
        A = a + rng.normal(0.0, sigma[sb], (M_, len(a)))
        B = b + rng.normal(0.0, sigma[sb], (M_, len(b)))
        Dx[sb], Dr[sb] = B - A, _expit(B) - _expit(A)
    ox, orr = group_share(Dx, s_x, c), group_share(Dr, s_raw, c)
    d = ox - orr
    nf = {sb: dict(sd_x=float(np.std(np.asarray(x132[sb]) - np.asarray(x62[sb]), ddof=1)),
                   sd_raw=float(np.std(_expit(x132[sb]) - _expit(x62[sb]), ddof=1))) for sb in sorted(x62)}
    return dict(ors_x=float(ox.mean()), ors_raw=float(orr.mean()), diff=float(d.mean()), diff_q025=float(np.percentile(d, 2.5)),
                diff_q975=float(np.percentile(d, 97.5)), noise_free_shift_sd=nf, reading=e6_reading(float(ox.mean()), float(orr.mean())))


def s3_algebra(nu_c: float, dnu: float, te_a: float = 100.0, te_b: float = 70.0, t2in: float = 82.5, t2en: float = 52.0) -> dict:
    """S3 — change of the group difference between sites for two groups that share the relaxation (true ν_c · ν_c + Δν): raw (D_A − D_B) · logit (should be 0) · first-order approximation. The same model as Table S15: config.weights (TE term only) · apparent (free water 0)."""
    def app(nu, te):
        return float(C.apparent(nu, 0.0, C.weights(te, np.array([t2in, t2en, 1000.0])))[0])

    def lg(p):
        return math.log(p / (1.0 - p))
    nu_p = nu_c + dnu
    a_c, a_p, b_c, b_p = app(nu_c, te_a), app(nu_p, te_a), app(nu_c, te_b), app(nu_p, te_b)
    raw = (a_p - a_c) - (b_p - b_c)
    logit = (lg(a_p) - lg(a_c)) - (lg(b_p) - lg(b_c))
    m_c, m_p = app(nu_c, 0.5 * (te_a + te_b)), app(nu_p, 0.5 * (te_a + te_b))       # Slope measured at the intermediate TE (midpoint approximation)
    first = (te_a - te_b) * (1.0 / t2en - 1.0 / t2in) * (m_p * (1.0 - m_p) - m_c * (1.0 - m_c))
    return dict(nu_c=nu_c, dnu=dnu, te_a=te_a, te_b=te_b, apparent_c_B=b_c, apparent_p_B=b_p, raw=raw, logit=logit, midpoint=first)


S3_GRID = dict(nu_c=(0.46, 0.60), dnu=(-0.03, -0.10, -0.20), te_a=(80.0, 90.0, 100.0))


def e5(out_dir, work_dir) -> dict:
    """E5 — does the per-unit TE-slope structure of E1 (β̂ · τ̂) explain E4 (the scale comparison of 62 versus 132 ms)?"""
    from . import empirical as EM
    w = EM.read_wong()
    sx_all, sr_all = EM.wong_series(w, "x"), EM.wong_series(w, "raw_ndi")
    P = _wong_pair_scales(EM, w, 62, 132, ("x", "raw_ndi"))
    obs = MT.ors(P["x"]["dc"], P["x"]["s"], 1.96, loo=False)["share"] - MT.ors(P["raw_ndi"]["dc"], P["raw_ndi"]["s"], 1.96, loo=False)["share"]
    e4p = pathlib.Path(out_dir) / "e4.json"
    e4v = json.loads(e4p.read_text(encoding="utf-8"))["primary"]["wong_62_132_x_minus_raw"] if e4p.exists() else None
    reg_ok = e4v is not None and gate_close(obs, e4v["diff"])
    res = dict(observed_x_minus_raw=obs, e4=e4v, regression_ok=reg_ok, by={})
    if not reg_ok:
        _dump(out_dir, "e5", res)
        raise SystemExit("⛔ E4 재현 실패 — 멈춤")
    par = {}
    for sb in sorted(sx_all):
        ser, ser_r = sx_all[sb], sr_all[sb]
        k, xbar = te_slopes(ser)
        us = sorted(k)
        kv, xv = np.array([k[u] for u in us]), np.array([xbar[u] for u in us])
        sig, _ = EM.s_rt_corrected(ser)
        sig_r, _ = EM.s_rt_corrected(ser_r)
        var_k, var_x = slope_level_noise(sig)
        tau, beta = MT.tau_hat(kv, var_k), MT.beta_level(kv, xv, var_x)
        nu_mid = np.array([(ser_r["TE62"][u] + ser_r["TE132"][u]) / 2.0 for u in us])
        par[sb] = dict(xbar=xv, kbar=float(kv.mean()), beta=beta, tau=tau, var_x=var_x, sigma=sig, s_x=sig, s_raw=sig_r)
        res["by"][sb] = dict(n_units=len(us), kbar=float(kv.mean()), beta=beta, tau=tau, nu_mid_mean=float(nu_mid.mean()), s_x=sig, s_raw=sig_r,
                             first_order=e5_first_order(float(kv.mean()), beta, float(nu_mid.mean())))
    res["predicted"] = e5_predict(par, use_beta=True)
    res["predicted_beta0"] = e5_predict(par, use_beta=False)
    res["reading"] = e5_reading(res["predicted"]["diff"], e4v["lo"], e4v["hi"])
    _dump(out_dir, "e5", res)
    return res


def e6(out_dir, work_dir) -> dict:
    """E6 — stored S1b noise-free fit values (no new fit): the scale comparison of 62 versus 132 ms when the ROI-matched units of Wong share the relaxation (one η for all units)."""
    from . import empirical as EM
    w = EM.read_wong()
    wm = _wong_matched(EM, w)
    keys = list(wm["units"])
    etas = eta17()
    wd = pathlib.Path(work_dir) / "s1b"
    jobs = {j["key"]: j for j in json.loads((wd / "jobs.json").read_text(encoding="utf-8"))["jobs"]}
    want = _cells_units_eta(wm["units"], etas)
    if not e6_inputs_ok(wd):                                    # Added after F6 on 2026-11-15 (not in the version that was run · identical in the post hoc comparison)
        raise SystemExit("⛔ E6 입력 sha256 이 계획 §7 과 다름 — 멈춤")
    for t in (62, 132):
        if not e6_cells_match(jobs[f"WONG:TE{t}"]["cells"], want):   # After F6: comparison of all cells (the version that was run had only nu · fiso · odi · gen · fdir)
            raise SystemExit("⛔ 저장된 S1b 일 목록이 ROI 맞춤 단위와 다름 — 멈춤")
    fits = dict(np.load(wd / "fits.npz"))

    def table(key):
        cm = coords_mean({k: fits[f"{key}|{k}"] for k in ("ndi", "fwf", "odi")})
        return cm["x"].reshape(len(keys), len(etas), len(GENS))
    X62, X132 = table("WONG:TE62"), table("WONG:TE132")
    sx_all, sr_all = EM.wong_series(w, "x"), EM.wong_series(w, "raw_ndi")
    sig = {sb: EM.s_rt_corrected(sx_all[sb])[0] for sb in sx_all}
    s_raw = {sb: EM.s_rt_corrected(sr_all[sb])[0] for sb in sr_all}
    idx = {sb: [i for i, (b, _) in enumerate(keys) if b == sb] for sb in sorted({b for b, _ in keys})}
    res = dict(n_units=len(keys), dropped=wm["dropped"], s_x=sig, s_raw=s_raw, by={})
    for gi, g in enumerate(GENS):
        for ei, (en, _, _) in enumerate(etas):
            res["by"][f"{g}|{en}"] = e6_predict({sb: X62[ii, ei, gi] for sb, ii in idx.items()}, {sb: X132[ii, ei, gi] for sb, ii in idx.items()},
                                                sig, sig, s_raw, M_=2000 if (g == "T-in" and en == CENTRE_NAME) else 500, seed=0)
    res["primary"] = res["by"][f"T-in|{CENTRE_NAME}"]
    res["reading_counts"] = {r: sum(v["reading"] == r for v in res["by"].values()) for r in ("가", "나", "다")}
    _dump(out_dir, "e6", res)
    return res


def s3(out_dir, work_dir) -> dict:
    """S3 — the fake site × group term on the raw scale (algebra · no data)."""
    rows = {}
    for nc in S3_GRID["nu_c"]:
        for dn in S3_GRID["dnu"]:
            for ta in S3_GRID["te_a"]:
                rows[f"{nc}|{dn}|{ta:.0f}"] = s3_algebra(nc, dn, te_a=ta)
    res = dict(model="config.weights (TE term only) · apparent (free water 0) · T2,in 82.5 · T2,en 52 ms · TE_B 70 ms", rows=rows,
               max_abs_logit=max(abs(v["logit"]) for v in rows.values()))
    _dump(out_dir, "s3", res)
    return res


E1_FIX_VARIANTS = {   # §8: name → (whether the batch size is used as n, whether it is a median factor, whether it is pooled degrees of freedom)
    "as_run": (False, False, False), "primary": (True, False, False), "sens_df": (True, False, True),
    "sens_median": (True, True, False), "sens_both": (True, True, True)}


def e1_fix_kwargs(spec, n_group, n_groups: int) -> dict:
    use_n, median, pooled_df = spec
    kw = {}
    if use_n:
        kw["n_group"] = n_group
    if median:
        kw["infl_kind"] = "median"
    if pooled_df:
        kw["tau_ddof"] = n_groups
    return kw


def e1_fix(out_dir, work_dir) -> dict:
    """ANALYSIS_PLAN §8 — fix of an implementation error in the pooled ΔTE* of E1 (the finite-n factor is taken per unit from the number of units of that participant) · three sensitivity analyses. Stop condition: the per-participant ΔTE* (already right, since it is a single bundle) must be recomputed with the same function and agree with e1.json to within 1e-9."""
    from . import empirical as EM
    w = EM.read_wong()
    stored = json.loads((pathlib.Path(out_dir) / "e1.json").read_text(encoding="utf-8"))
    res = dict(rule="ANALYSIS_PLAN §8", pooled={}, regression={})
    ok = True
    for co in ("x", "y"):
        ser_all = EM.wong_series(w, co)
        kc, su, cl, ng, vks, dfs = [], [], [], [], [], []
        for pi, sb in enumerate(sorted(ser_all)):
            ser = ser_all[sb]
            k, xbar = te_slopes(ser)
            us = sorted(k)
            kv = np.array([k[u] for u in us])
            sig, df = EM.s_rt_corrected(ser)
            var_k = slope_level_noise(sig)[0]
            pt_sb = MT.te_tolerance(np.full(len(us), sig), MT.tau_hat(kv, var_k))
            ref = stored["by"][f"{co}|{sb}"]["dte_star"]
            res["regression"][f"{co}|{sb}"] = dict(recomputed=pt_sb, stored=ref)
            ok &= abs(pt_sb - ref) < 1e-9
            kc += list(kv - kv.mean()); su += [sig] * len(us); cl += [pi] * len(us); ng += [len(us)] * len(us)
            vks.append((var_k, len(us))); dfs.append(df)
        var_k_pool = sum(v * n for v, n in vks) / sum(n for _, n in vks)
        S, K, CL, NG, P = np.array(su), np.array(kc), np.array(cl), np.array(ng, float), len(vks)
        out = {}
        for name, spec in E1_FIX_VARIANTS.items():
            kw = e1_fix_kwargs(spec, NG, P)
            pt, lo, hi = MT.te_tolerance_interval(S, K, var_k_pool, sum(dfs), clusters=CL, B=2000, seed=0, **kw)
            out[name] = dict(dte_star=pt, dte_ci90=[lo, hi], tau=MT.tau_hat(K, var_k_pool, kw.get("tau_ddof", 1)))
        out["as_run_matches_stored"] = abs(out["as_run"]["dte_star"] - stored["pooled"][co]["dte_star"]) < 1e-9 and \
            all(abs(a - b) < 1e-9 for a, b in zip(out["as_run"]["dte_ci90"], stored["pooled"][co]["dte_ci90"]))
        ok &= out["as_run_matches_stored"]
        res["pooled"][co] = out
    res["regression_ok"] = bool(ok)
    _dump(out_dir, "e1_fix", res)
    if not ok:
        raise SystemExit("⛔ E1 재현 실패 — 멈춤")
    return res

RUNNERS = dict(s1a=s1a, s1b=s1b, s2=s2, k8=k8, e1=e1, e2=e2, e3=e3, e4=e4, e5=e5, e6=e6, s3=s3, e1_fix=e1_fix)

if __name__ == "__main__":
    RUNNERS[sys.argv[1]](pathlib.Path(sys.argv[2]), pathlib.Path(sys.argv[3]))
