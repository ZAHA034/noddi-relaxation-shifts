"""Computes s4 Appendix B Fisher sensitivity: the MDE^lin ratio before and after adding nuisance factors ψ to the linear baseline; not used for verdicts."""
from __future__ import annotations

import hashlib
import json
import os
import resource
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s4cfg as S  # noqa: E402
import simulate as SIM  # noqa: E402
import tables as TB  # noqa: E402

LABEL = "선형 근사 — 1 단계에서 선형 지표가 식별성을 과대평가함이 확인됨 · 판정에 쓰지 않음"
PSI = ("kappa_A", "kappa_B", "g_A", "g_B", "p_in", "p_iso")
FACTORS = dict(G3=("kappa_A", "kappa_B"), G2=("g_A", "g_B"), G1=("p_in", "p_iso"))
WIDTHS = dict(G3=(0.005, 0.01, 0.02, 0.03, 0.05), G2=(0.02, 0.05, 0.10, 0.15, 0.20), G1=(0.02, 0.05, 0.10, 0.20, 0.30))   # Factor SD [design]
KAPPA_H, G_H, P_H = 0.05, 0.05, 0.05        # Difference step — κ: center of S4's b multiplier table (0.95 · 1.05) · g: one-sided 0 → 0.05 · p: center ± 0.05
GM = (0.35, 0.10, 0.50)                      # Gray-matter tissue (NDI · FISO · ODI) [design · not compared with the literature] · relaxation is that person's own, as is
VARIANTS = (("base", {}), ("kappa_A+", dict(bs_A=1.05)), ("kappa_A-", dict(bs_A=0.95)), ("kappa_B+", dict(bs_B=1.05)),
            ("kappa_B-", dict(bs_B=0.95)), ("g_A", dict(g_A=G_H)), ("g_B", dict(g_B=G_H)), ("p_in+", dict(p_in=P_H)),
            ("p_in-", dict(p_in=-P_H)), ("p_iso+", dict(p_iso=P_H)), ("p_iso-", dict(p_iso=-P_H)))
N_OBS, N_OBS_FALLBACK, CAP_CPU_H = 10, 5, 2.0
SEED = 20261107
U_NAMES = ("NDI", "FISO", "ODI", "u_relax")
SUFFIX = dict(maps="|L0", cards="|L0|SNR30|Neff50")


# ── Cell · variant ───────────────────────────────────────────────────────────────────
def fisher_cells(cells_maps, cells_cards):
    """Cell 77 — [(dataset, cell name, case numbers of that cell (ascending))] · maps are «…|L0» · cards are the main cell «…|L0|SNR30|Neff50»."""
    out = []
    for name, cells in (("maps", cells_maps), ("cards", cells_cards)):
        for cl in sorted(set(np.asarray(cells).tolist())):
            if cl.endswith(SUFFIX[name]):
                out.append((name, cl, np.flatnonzero(cells == cl)))
    return out


def configs():
    """Width settings 20 — each factor separately (3 × 5) + all three together (width of the same rank · 5). [(tag, {factor: SD})]."""
    out = [(f"{f}|{i}", {f: w}) for f, ws in WIDTHS.items() for i, w in enumerate(ws)]
    return out + [(f"all|{i}", {f: WIDTHS[f][i] for f in WIDTHS}) for i in range(5)]


def active(cfg):
    """Width setting → ψ column number (PSI order) · prior SD."""
    cols = [PSI.index(p) for f in cfg for p in FACTORS[f]]
    return np.array(cols), np.array([cfg[f] for f in cfg for _ in FACTORS[f]], float)


# ── J_ψ simulation (envs/reference · AMICO) ───────────────────────────────────────
def s0(pix, a):
    """b0 signal (before normalization) (n·K, 1) — the same formula as S0 of simulate.signal."""
    K = pix.shape[1]
    nu, f = pix[..., 0].reshape(-1, 1), pix[..., 1].reshape(-1, 1)
    ain, aen, aiso = (np.repeat(a[:, j], K)[:, None] for j in range(3))
    return (1 - f) * (nu * ain + (1 - nu) * aen) + f * aiso


def gm_mix(Sg, Sgm, s0t, s0g, g):
    """Voxel signal = (1 − g) · tissue + g · gray matter — mixed before b0 normalization and normalized again (Appendix B factor 2)."""
    wt, wg = (1 - g) * s0t, g * s0g
    return ((wt * Sg + wg * Sgm) / (wt + wg)).astype(np.float32)


def sim_variant(P, idx, var, Ttab, study, tag, nthreads, seed):
    """Same order and same random-number consumption as simulate.chunk — only the ψ of the variant var is changed (common random numbers · same seed). var: bs_A · bs_B (b scaling of the synthesis — the fit uses the recorded b) · g_A · g_B (gray-matter partial volume) · p_in · p_iso (e^p in the compartment weights · common to the two sites). Returns: x_c7 (n, 14) · x_c3 (n, 6) · clipped (n,)."""
    import lib  # noqa: E402 — decomp/stage1
    import synth  # noqa: E402
    rng = np.random.default_rng(seed)
    n = len(idx)
    SIM.oracle_obs(rng, P, idx)                                                     # Same random-number consumption as chunk() (the values are not used)
    th = P["theta"][idx]
    dpi = np.array([TB.DPAR_T.index(float(v)) for v in P["dpar"][idx]])
    pd = np.exp([var.get("p_in", 0.0), 0.0, var.get("p_iso", 0.0)])
    raws, clipped = {}, {}
    for s, cen in (("A", th[:, 3:6]), ("B", th[:, 3:6] + th[:, 0:3])):
        pix = SIM.pixels(rng, cen)
        a = SIM.weights(P[f"acq_TE_{s}"][idx], P[f"acq_TR_{s}"][idx], P[f"t2{s}"][idx], P[f"t1{s}"][idx]) * pd[None, :]
        design = P[f"acq_design_{s}"][idx]
        snr = P[f"acq_snr_{s}"][idx]
        L = np.array([S.COIL_L[int(c)] for c in P[f"acq_coil_{s}"][idx]])
        bs, g = float(var.get(f"bs_{s}", 1.0)), float(var.get(f"g_{s}", 0.0))
        raw = np.full((n, S.K_PIX, 7), np.nan, np.float32)
        for dsg in S.DESIGNS:
            m = design == dsg
            if not m.any():
                continue
            ii = np.flatnonzero(m)
            nom, T = Ttab[(dsg, 1.0)], Ttab[(dsg, bs)]
            gen, fd = P["gen"][idx][ii], P[f"fdir_{s}"][idx][ii]
            Sg = SIM.signal(T, pix[ii], a[ii], gen, dpi[ii], fd)
            if g:
                gm = np.broadcast_to(np.asarray(GM, float), pix[ii].shape).copy()
                Sg = gm_mix(Sg, SIM.signal(T, gm, a[ii], gen, dpi[ii], fd), s0(pix[ii], a[ii]), s0(gm, a[ii]), g)
            Sn = SIM.noisy(rng, Sg, np.repeat(snr[ii], S.K_PIX), np.repeat(L[ii], S.K_PIX))
            vol = Sn.reshape(len(ii), S.K_PIX, -1).transpose(1, 0, 2)[:, :, None, :]          # (K, n_g, 1, V)
            ndi, fwf, odi = SIM.amico_fit(study[dsg], f"{tag}_{s}", vol, nom["bvals"], nom["bvecs"], nthreads)
            ext = synth.ext_coords(Sn.astype(float), nom["bvals"], nom["bvecs"], SIM.LOW_SHELL).reshape(len(ii), S.K_PIX, 4)
            raw[ii] = np.concatenate([np.stack([ndi.T, fwf.T, odi.T], -1), ext], -1)
        fin = np.isfinite(raw).all(-1)
        mu = np.array([raw[i][fin[i]].mean(0) for i in range(n)])
        Sig = np.array([np.cov(raw[i][fin[i]].T) for i in range(n)])
        roi = S.roi_noise(rng, mu, Sig, P["acq_neff"][idx])
        y, clipped[s] = lib.coords(roi[:, 0], roi[:, 1], roi[:, 2], roi[:, 3:7], "C7")
        raws[s] = y.astype(np.float32)
    return np.concatenate([raws["A"], raws["B"]], 1), np.c_[raws["A"][:, :3], raws["B"][:, :3]], clipped["A"] | clipped["B"]


def _cpu_s():
    a, b = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
    return a.ru_utime + a.ru_stime + b.ru_utime + b.ru_stime


def _log(obj):
    obj = dict(UTC=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), **obj)
    print(json.dumps(obj, ensure_ascii=False), flush=True)
    with open(S.WORK / "logs" / "fisher_sim.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(obj, ensure_ascii=False) + "\n")


def run_sim(nthreads=16):
    sys.path.insert(0, str(S.C.RESULTS.parent / "stage1"))
    (S.WORK / "logs").mkdir(parents=True, exist_ok=True)
    cpu0, wall0 = _cpu_s(), time.time()
    plans = {nm: dict(np.load(S.WORK / f"plan_{nm}.npz", allow_pickle=False)) for nm in ("maps", "cards")}
    cells = fisher_cells(plans["maps"]["cell"], plans["cards"]["cell"])
    assert len(cells) == 77, len(cells)
    Ttab, study = {}, {}
    for dsg in S.DESIGNS:
        for bs in (1.0, 0.95, 1.05):
            Ttab[(dsg, bs)] = TB.load(dsg, bs)
        study[dsg] = str(S.WORK / "amico" / dsg)
        os.makedirs(study[dsg], exist_ok=True)

    def pick(nobs):
        sel = [(nm, cl, ii[:nobs]) for nm, cl, ii in cells]
        return sel, {k: np.concatenate([plans[nm][k][ii] for nm, _, ii in sel]) for k in plans["maps"]}

    sel1, P1 = pick(1)
    i3 = np.arange(3)                                                               # Self-check: reference variant = simulate.chunk (same seed)
    xa = SIM.chunk(P1, i3, np.random.default_rng([SEED, 9]), Ttab, study, "fisher_chk_a", nthreads)["x_c7"]
    xb = sim_variant(P1, i3, {}, Ttab, study, "fisher_chk_b", nthreads, [SEED, 9])[0]
    chk = float(np.max(np.abs(xa - xb)))
    if not chk <= 1e-5:
        raise SystemExit(f"⛔ 기준 변형이 simulate.chunk 와 다르다(최대 차 {chk}) — 멈춘다")
    c0, w0 = _cpu_s(), time.time()                                                  # Estimate before starting: 1 observation per cell · one variant = 77 cases
    sim_variant(P1, np.arange(len(P1["theta"])), {}, Ttab, study, "fisher_probe", nthreads, [SEED, 1])
    per_cpu, per_wall = (_cpu_s() - c0) / 77, (time.time() - w0) / 77
    spent = (_cpu_s() - cpu0) / 3600
    est = {nobs: spent + per_cpu * 77 * nobs * len(VARIANTS) / 3600 for nobs in (N_OBS, N_OBS_FALLBACK)}
    nobs = N_OBS if est[N_OBS] <= CAP_CPU_H else N_OBS_FALLBACK
    _log(dict(step="estimate", self_check_max_abs_diff=chk, cpu_s_per_case=round(per_cpu, 3), wall_s_per_case=round(per_wall, 3),
              threads=nthreads, estimate_cpu_hours={str(k): round(v, 3) for k, v in est.items()}, cap_cpu_hours=CAP_CPU_H, n_obs=nobs))
    if est[nobs] > CAP_CPU_H:
        raise SystemExit("⛔ 시작 전 추정이 상한 2 CPU-시간을 넘는다(관측 5 로도) — 멈추고 알린다")
    sel, P = pick(nobs)
    allidx = np.arange(len(P["theta"]))
    X7, X3, CL = [], [], []
    for name, var in VARIANTS:
        t0 = time.time()
        x7, x3, cl = sim_variant(P, allidx, var, Ttab, study, f"fisher_{name}", nthreads, [SEED, 0])
        X7.append(x7), X3.append(x3), CL.append(cl)
        _log(dict(step="variant", variant=name, n=len(allidx), seconds=round(time.time() - t0, 1), cpu_hours=round((_cpu_s() - cpu0) / 3600, 3)))
    path = S.WORK / "fisher_sim.npz"
    np.savez(path, x_c7=np.stack(X7), x_c3=np.stack(X3), clipped=np.stack(CL), variants=np.array([v for v, _ in VARIANTS]),
             cell=np.repeat([cl for _, cl, _ in sel], nobs), dataset=np.repeat([nm for nm, _, _ in sel], nobs),
             plan_idx=np.concatenate([ii for _, _, ii in sel]), n_obs=nobs)
    meta = dict(step="done", file=path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest(), n_obs=nobs,
                cpu_hours=round((_cpu_s() - cpu0) / 3600, 3), wall_hours=round((time.time() - wall0) / 3600, 3), threads=nthreads,
                wall_x_threads_hours=round((time.time() - wall0) * nthreads / 3600, 3), self_check_max_abs_diff=chk)
    json.dump(meta | dict(estimate_cpu_hours={str(k): v for k, v in est.items()}), open(S.WORK / "fisher_sim.json", "w"), ensure_ascii=False, indent=1)
    _log(meta)


# ── Ratio (main environment · numpy) ────────────────────────────────────────────────────────
def regress(x, th):
    """§13 12 — linear regression (with intercept) of x (M, d) on θ* (M, 6) → J_θ (d, 6) · residual covariance R (d, d). The same formula as validate.fisher_cell."""
    A = np.c_[np.ones(len(th)), th]
    beta, *_ = np.linalg.lstsq(A, x, rcond=None)
    return beta[1:].T, np.cov((x - A @ beta).T) + 1e-9 * np.eye(x.shape[1])


def post(Jt, R, covT, Jp=None, w=None):
    """I = [J_θ J_ψ]ᵀ R⁻¹ [J_θ J_ψ] + blockdiag(Σ_θ⁻¹, diag(w⁻²)) → (I⁻¹)_θθ (6, 6) · posterior SD of ψ (q,). If Jp is absent, the baseline (= fisher_cell)."""
    J = Jt if Jp is None else np.c_[Jt, Jp]
    q = J.shape[1] - Jt.shape[1]
    P0 = np.zeros((J.shape[1], J.shape[1]))
    P0[:Jt.shape[1], :Jt.shape[1]] = np.linalg.inv(covT)
    if q:
        P0[Jt.shape[1]:, Jt.shape[1]:] = np.diag(1.0 / np.asarray(w, float) ** 2)
    Iinv = np.linalg.inv(J.T @ np.linalg.solve(R, J) + P0)
    return Iinv[:Jt.shape[1], :Jt.shape[1]], np.sqrt(np.diag(Iinv)[Jt.shape[1]:])


def mde_ratio(Jt, R, covT, Jp, w, U):
    """For each direction U (k, 3), MDE^lin(ψ, w) / MDE^lin(baseline) (2.8 cancels) · posterior SD / prior SD of ψ (q,)."""
    S0, _ = post(Jt, R, covT)
    S1, sd = post(Jt, R, covT, Jp, w)
    q = lambda Sg: np.sqrt(np.einsum("ki,ij,kj->k", U, Sg[:3, :3], U))
    return q(S1) / q(S0), sd / np.asarray(w, float)


def jpsi(X, names):
    """x per variant (variant, n, d) → J_ψ (d, 6) — column order PSI · common-random-number difference of the observed mean."""
    v = dict(zip(names, X))
    cols = ((v["kappa_A+"] - v["kappa_A-"]) / (2 * KAPPA_H), (v["kappa_B+"] - v["kappa_B-"]) / (2 * KAPPA_H),
            (v["g_A"] - v["base"]) / G_H, (v["g_B"] - v["base"]) / G_H,
            (v["p_in+"] - v["p_in-"]) / (2 * P_H), (v["p_iso+"] - v["p_iso-"]) / (2 * P_H))
    return np.stack([np.asarray(c, float).mean(0) for c in cols], 1)


def _f(v):
    v = np.asarray(v, float)
    return [None if not np.isfinite(t) else float(t) for t in v]


def _stat(A):
    """Cell × direction → median · minimum · maximum for each direction (cells where it is not defined are dropped — u_relax of matching pairs)."""
    A = np.asarray(A, float)
    out = dict(median=[], min=[], max=[], n_cells=int(len(A)))
    for j in range(A.shape[1]):
        v = A[np.isfinite(A[:, j]), j]
        for k, fn in (("median", np.median), ("min", np.min), ("max", np.max)):
            out[k].append(float(fn(v)) if len(v) else None)
    return out


def run_ratio():
    import validate as V  # noqa: E402 — direction (u_relax) · main prior covariance (same seed as cells_summary) · writing results
    F = dict(np.load(S.WORK / "fisher_sim.npz", allow_pickle=False))
    meta = json.load(open(S.WORK / "fisher_sim.json", encoding="utf-8"))
    names = [str(v) for v in F["variants"]]
    covT = V.cov_prior(np.random.default_rng(4))
    sims = {nm: V.load_set(nm) for nm in ("maps", "cards")}
    cells = fisher_cells(sims["maps"][1]["cell"], sims["cards"][1]["cell"])
    assert len(cells) == 77, len(cells)
    cfgs = configs()
    res, rmin = {}, np.inf
    for nm, cl, ii in cells:
        sim, plan = sims[nm]
        th = sim["theta"][ii]
        ur = np.nanmean(V.directions(th, plan, ii)[:, 3], 0)                       # u_relax — same cell-mean direction as in cells_summary
        ur = ur / np.linalg.norm(ur) if np.isfinite(ur).all() and np.linalg.norm(ur) > 1e-12 else np.full(3, np.nan)
        U = np.r_[np.eye(3), ur[None]]
        m = F["cell"] == cl
        res[cl] = {}
        for arm in ("c3", "c7"):
            Jt, R = regress(V.xobs(sim, arm)[ii], th)
            Jp = jpsi(F[f"x_{arm}"][:, m], names)
            ent = {}
            for tag, cfg in cfgs:
                cols, w = active(cfg)
                r, sdr = mde_ratio(Jt, R, covT, Jp[:, cols], w, U)
                rmin = min(rmin, float(np.nanmin(r)))
                e = dict(ratio=_f(r))
                if "G1" in cfg:
                    e["sd_ratio_p"] = dict(p_in=float(sdr[list(cols).index(PSI.index("p_in"))]),
                                           p_iso=float(sdr[list(cols).index(PSI.index("p_iso"))]))
                ent[tag] = e
            res[cl][arm] = ent
    groups = dict(maps_L0=[cl for nm, cl, _ in cells if nm == "maps"], cards_main=[cl for nm, cl, _ in cells if nm == "cards"],
                  all=[cl for _, cl, _ in cells])
    curves = {}
    for arm in ("c3", "c7"):
        curves[arm] = {}
        for g, cls in groups.items():
            curves[arm][g] = {}
            for f in ("G3", "G2", "G1", "all"):
                pts = []
                for i in range(5):
                    tag = f"{f}|{i}"
                    p = dict(width=dict(cfgs)[tag], ratio=_stat([res[cl][arm][tag]["ratio"] for cl in cls]))
                    if f in ("G1", "all"):
                        p["sd_ratio_p"] = {k: _stat([[res[cl][arm][tag]["sd_ratio_p"][k]] for cl in cls]) for k in ("p_in", "p_iso")}
                    pts.append(p)
                curves[arm][g][f] = pts
    clip = F["clipped"]
    return dict(label=LABEL, spec="decomp/SPEC_S4.md 부록 B(중계 ZZZo-3)", report="비만 — 절대 MDE 는 S4 값만 쓴다",
                directions=list(U_NAMES), psi=list(PSI), factors=FACTORS, widths=WIDTHS, steps=dict(kappa=KAPPA_H, g=G_H, p=P_H),
                gm_tissue=dict(zip(("NDI", "FISO", "ODI"), GM)), n_cells=len(cells), n_obs_per_cell=int(F["n_obs"]),
                variants=names, seed=SEED, sim=dict(file="fisher_sim.npz", sha256=meta["sha256"], cpu_hours=meta["cpu_hours"],
                                                    wall_hours=meta["wall_hours"], threads=meta["threads"],
                                                    self_check_max_abs_diff=meta["self_check_max_abs_diff"]),
                clipped_fraction_by_variant={n: float(c.mean()) for n, c in zip(names, clip)},
                check_ratio_ge_1=dict(min_ratio=rmin, ok=bool(rmin >= 1 - 1e-9)),
                note_sd_ratio_p="p 방향의 사후 SD / 사전 SD — 부록 B 의 명제 확인용. θ* 의 사전이 정보를 가지면 정확한 등가에서도 1 보다 작을 수 있다(ZZZo-1 (4) 와 같은 까닭)",
                curves=curves, cells=res)


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[0] == "sim":
        run_sim(int(a[1]) if len(a) > 1 else 16)
    elif a[0] == "ratio":
        import validate as V  # noqa: E402
        V.jdump("fisher", run_ratio())
    else:
        raise SystemExit(a[0])
