"""Exploratory analysis b: finds the 25° solution family on a noise-free T2_in × r × T2_iso × true d∥ grid fitted by AMICO F2."""
import hashlib, itertools, json, pathlib, shutil, sys, time
import numpy as np
X14 = pathlib.Path("<WORK>/stage1_work/drafts/x14")
sys.path.insert(0, str(X14 / "B1S"))
import b1s_sim as B  # noqa: E402 — Wong direction scheme · pix_signal (verification reference) · tissue · TE · TR · T1
C, SIM, lib = B.C, B.SIM, B.lib
OUT = pathlib.Path("<WORK>/remeasure/explore/rv")
T2IN = (60.0, 70.0, 82.0, 100.0, 125.0, 150.0, 180.0, 210.0, 250.0)
RS = (1.0, 1.3, 1.6, 2.0, 2.4, 2.8, 3.2)
T2ISO = (80.0, 100.0, 119.0, 150.0, 200.0, 300.0, 500.0, 1000.0, 2000.0)
DPARS = (1.7, 2.0, 2.26)
OBS = dict(bx=7.939, by=3.424, bx_ci=(7.293, 8.802), by_ci=(2.303, 4.187))                 # s⁻¹ · B1 81 ROI units · issued before computing
NTHREADS = 12


def signal(pr, nu, f, odi, fd, t2, dpar, te, E=None):
    """Same formula as B1-S pix_signal — takes the T2 triple (in, en, iso) and the true d∥ as parameters."""
    a = C.weights(float(te), t2, B.TR, (B.T1T, B.T1T, B.T1ISO))
    if E is None:
        kp = C.odi_to_kappa(odi)
        E = (pr.E_ic(dpar, kp, fd), np.array([pr.E_ec(dpar, C.dperp("T-in", dpar, n_, None), kp, fd) for n_ in np.atleast_1d(nu)]), pr.E_iso())
    Eic, Eec, Eiso = E
    nu, f = np.atleast_1d(nu)[:, None], np.atleast_1d(f)[:, None]
    Sg = (1 - f) * (nu * a[0] * Eic[None] + (1 - nu) * a[1] * Eec) + f * a[2] * Eiso[None]
    S0 = (1 - f) * (nu * a[0] + (1 - nu) * a[1]) + f * a[2]
    return Sg / S0, E


def main():
    pr = B.synth.Proto("WONG_SUB01_TE92")
    nu0, odi0, f0 = 0.6, 0.20, 0.10
    ref, _ = B.pix_signal(pr, nu0, f0, odi0, 0, 82.0, 82.0 / 1.6, 92.0)
    mine, _ = signal(pr, nu0, f0, odi0, 0, (82.0, 82.0 / 1.6, 1000.0), 1.7, 92.0)
    vs = float(np.max(np.abs(ref - mine)))
    print("V-S′ 최대 차", vs, flush=True)
    if vs > 1e-12:
        sys.exit("⛔ V-S′ — 다시 쓴 식이 pix_signal 과 다르다")
    study = OUT / "amico" / "W"
    if not (study / "kernels").exists():
        shutil.copytree(B.OUT / "amico" / "W" / "kernels", study / "kernels")    # Copy the B1-S kernel (same direction scheme) — do not write into the b1s folder
    conds = list(itertools.product(T2IN, RS, T2ISO, DPARS))
    tis = B.TISSUES
    E = {(ti, dp): signal(pr, nu, f, odi, ti % len(C.FIBRE_DIRS), (82.0, 50.0, 1000.0), dp, 92.0)[1] for ti, (nu, odi, f) in enumerate(tis) for dp in DPARS}
    n = len(conds) * len(tis)
    coords = np.full((len(B.TES), len(conds), len(tis), 3), np.nan)
    t0 = time.time()
    for j, te in enumerate(B.TES):
        Sg = np.empty((n, len(B.BVALS)), np.float32)
        k = 0
        for (t2in, r, t2iso, dp) in conds:
            for ti, (nu, odi, f) in enumerate(tis):
                Sg[k] = signal(pr, nu, f, odi, ti % len(C.FIBRE_DIRS), (t2in, t2in / r, t2iso), dp, te, E[(ti, dp)])[0][0]
                k += 1
        vol = Sg.reshape(len(conds), len(tis), 1, -1)                                   # ≤ 32767 per NIfTI axis — so as two axes (condition, tissue)
        ndi, fwf, odi_ = SIM.amico_fit(str(study), f"rv_b_te{te}", vol, B.BVALS, B.BVECS, NTHREADS)
        y, _ = lib.coords(ndi.ravel(), fwf.ravel(), odi_.ravel(), None, "C3")
        coords[j] = y.reshape(len(conds), len(tis), 3)
        print(json.dumps(dict(TE=te, n=n, seconds=round(time.time() - t0, 1), finite=float(np.isfinite(y).all(1).mean()))), flush=True)
    tes = np.array(B.TES, float)
    tc = tes - tes.mean()
    slope = (tc[:, None, None, None] * (coords - coords.mean(0))).sum(0) / (tc ** 2).sum()      # (condition, tissue, 3) · unit ms⁻¹
    bx, by = 1000 * np.nanmedian(slope[:, :, 0], 1), 1000 * np.nanmedian(slope[:, :, 1], 1)
    match = (bx >= OBS["bx_ci"][0]) & (bx <= OBS["bx_ci"][1]) & (by >= OBS["by_ci"][0]) & (by <= OBS["by_ci"][1])
    rows = [dict(t2in=c[0], r=c[1], t2en=c[0] / c[1], t2iso=c[2], dpar=c[3], bx=float(x), by=float(yv), theta=float(np.degrees(np.arctan2(yv, x))), match=bool(m))
            for c, x, yv, m in zip(conds, bx, by, match)]
    lit = [r for r in rows if r["t2in"] == 82.0 and r["r"] == 1.6 and r["t2iso"] == 1000.0 and r["dpar"] == 1.7][0]
    alg = []
    for t2iso in np.geomspace(80, 2000, 60):                                              # Algebra batch: 1/T2_in = b_y + 1/T2_iso · 1/T2_en = 1/T2_in + b_x
        inv_in = OBS["by"] / 1000 + 1 / t2iso
        alg.append(dict(t2iso=float(t2iso), t2in=float(1 / inv_in), t2en=float(1 / (inv_in + OBS["bx"] / 1000))))
    res = dict(made=time.strftime("%Y-%m-%dT%H:%M:%S%z"), script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(), V_S_prime=vs,
               obs=OBS, grid=dict(t2in=T2IN, r=RS, t2iso=T2ISO, dpar=DPARS), n_conditions=len(conds), n_match=int(match.sum()),
               match_by_dpar={str(dp): int(sum(r["match"] for r in rows if r["dpar"] == dp)) for dp in DPARS},
               check_lit_82_52_1000_d17=dict(theta=lit["theta"], bx=lit["bx"], by=lit["by"], b1s_noisefree_theta=57.4),
               matches=[r for r in rows if r["match"]], algebraic_family=alg, rows=rows)
    (OUT / "rv_b_family.json").write_text(json.dumps(res, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
    for a, dp in zip(ax, DPARS):
        m = [r for r in rows if r["dpar"] == dp and r["match"]]
        a.plot([r["t2iso"] for r in alg], [r["t2in"] for r in alg], color="0.5", lw=1, label="algebraic family (no fit)")
        sc = a.scatter([r["t2iso"] for r in m], [r["t2in"] for r in m], c=[r["r"] for r in m], cmap="viridis", vmin=1.0, vmax=3.2, s=30, label="fit-through match")
        a.set_xscale("log"); a.set_xlabel("T2_iso (ms)"); a.set_title(f"true d∥ {dp} µm²/ms · fit d∥ 1.7")
    ax[0].set_ylabel("T2_in (ms)"); ax[0].legend(fontsize=7, loc="upper left")
    fig.colorbar(sc, ax=ax, label="r = T2_in / T2_en")
    fig.suptitle("b — (T2_in, T2_en, T2_iso) consistent with Wong TE slopes (noise-free, exploratory)")
    fig.savefig(OUT / "rv_b_family.png", dpi=150, bbox_inches="tight")
    print(json.dumps(dict(n_match=res["n_match"], match_by_dpar=res["match_by_dpar"], lit_check=res["check_lit_82_52_1000_d17"]), ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
