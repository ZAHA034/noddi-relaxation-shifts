"""Step 2 of b1s: simulates the Wong processing with AMICO, using the noiseless signal formula, Wong sub-01 TE92 directions and noncentral chi noise."""
import hashlib, json, pathlib, sys, time
import numpy as np
REPO = pathlib.Path("<REPO>")
sys.path.insert(0, str(REPO / "decomp" / "s4"))
import s4cfg as S  # noqa: E402
import simulate as SIM  # noqa: E402
sys.path.insert(0, str(REPO / "decomp" / "stage1"))
import config as C  # noqa: E402
import lib  # noqa: E402
import synth  # noqa: E402
OUT = pathlib.Path("<WORK>/remeasure/explore/b1s")
TES = (62, 72, 82, 92, 102, 112, 122, 132)
TR, T1T, T1ISO, DPAR, K = 5800.0, 1100.0, C.T1_ISO, 1.7, 200
T2SETS = {"lit_82_52": (82.0, 52.0), "wedge_48.5": (95.0, 95.0 / 1.8), "wedge_72.1": (70.0, 70.0 / 1.3), "pos_r2.9": (82.0, 82.0 / 2.9)}
TISSUES = [(nu, odi, f) for nu in (0.45, 0.6, 0.75) for odi in (0.08, 0.20) for f in (0.0, 0.02, 0.05, 0.10)]
CONDS = [("noisefree", 0.0, 4, 1)] + [(f"snr×{s}|L{L}", s, L, 5) for s in (0.5, 1.0, 2.0) for L in (4, 8)]   # L ≥ 2 (χ² degrees of freedom 2(L−1) > 0 in simulate.noisy) · matches the 32 · 64 channels of S4
SEED, NTHREADS = 20261161, 12
sch = np.load(OUT / "wong_sub01_TE92_scheme.npz")
BVALS, BVECS = sch["bvals"], sch["bvecs"]
_orig = C.protocol
C.protocol = lambda name: (BVALS.copy(), BVECS.copy()) if name == "WONG_SUB01_TE92" else _orig(name)
SNR = {int(k): v["snr"] for k, v in json.loads((OUT / "b1s_snr.json").read_text(encoding="utf-8"))["te"].items()}


def pix_signal(pr, nu, f, odi, fd, t2in, t2en, te, E=None):
    """The same formula as synth.noiseless — continuous ν · f for each voxel (one ODI · direction per ROI) · returns (K, V) S/S0."""
    a = C.weights(float(te), (t2in, t2en, 1000.0), TR, (T1T, T1T, T1ISO))
    kp = C.odi_to_kappa(odi)
    if E is None:
        E = (pr.E_ic(DPAR, kp, fd), np.array([pr.E_ec(DPAR, C.dperp("T-in", DPAR, n_, None), kp, fd) for n_ in np.atleast_1d(nu)]), pr.E_iso())
    Eic, Eec, Eiso = E
    nu, f = np.atleast_1d(nu)[:, None], np.atleast_1d(f)[:, None]
    Sg = (1 - f) * (nu * a[0] * Eic[None] + (1 - nu) * a[1] * Eec) + f * a[2] * Eiso[None]
    S0 = (1 - f) * (nu * a[0] + (1 - nu) * a[1]) + f * a[2]
    return Sg / S0, E


def vcheck(pr):
    """V-S — grid tissue (ν 0.6 · ODI 0.20 · FISO 0.10) · T2 set (82, 1.6, 1000) · TE 92 · TR 5800 · T1 (1100, 1.0) · T-in · direction 0: is it the same as synth.noiseless?"""
    ti = C.TISSUE.index((0.60, 0.20, 0.10)); k = C.T2_SETS.index((82.0, 1.6, 1000.0)); t1i = C.T1_SETS.index((1100.0, 1.0))
    ref, _, _, _ = synth.noiseless(pr, dict(tissue=ti, t2set=k, TE=92.0, TR=TR, t1set=t1i, dpar=DPAR, gen="T-in", fdir=0))
    mine, _ = pix_signal(pr, 0.60, 0.10, 0.20, 0, 82.0, 82.0 / 1.6, 92.0)
    return float(np.max(np.abs(np.asarray(ref, float) - mine[0])))


def main():
    pr = synth.Proto("WONG_SUB01_TE92")
    dv = vcheck(pr)
    print("V-S 최대 차", dv, flush=True)
    if dv > 1e-12:
        sys.exit("⛔ V-S — 다시 쓴 식이 synth.noiseless 와 다르다")
    rng = np.random.default_rng(SEED)
    lo, hi = np.array([S.PIX_CLIP["NDI"][0], S.PIX_CLIP["FISO"][0]]), np.array([S.PIX_CLIP["NDI"][1], S.PIX_CLIP["FISO"][1]])
    pix, Es = [], []
    for ti, (nu, odi, f) in enumerate(TISSUES):
        p = np.clip(np.array([nu, f])[None] + rng.normal(0, 1, (K, 2)) * S.PIX_SD[[0, 1]], lo, hi)
        fd = ti % len(C.FIBRE_DIRS)
        _, E = pix_signal(pr, p[:, 0], p[:, 1], odi, fd, 82.0, 52.0, 62.0)
        pix.append((p, odi, fd)); Es.append(E)
    meta, out = [], []
    for c_i, (cname, scale, L, nrep) in enumerate(CONDS):
        for t2n, (t2in, t2en) in T2SETS.items():
            for ti in range(len(TISSUES)):
                for rep in range(nrep):
                    meta.append(dict(cond=cname, t2set=t2n, tissue=TISSUES[ti], rep=rep, scale=scale, L=L))
    n = len(meta)
    res = np.full((n, len(TES), 3), np.nan)
    t0 = time.time()
    for j, te in enumerate(TES):
        Sg = np.empty((n, K, len(BVALS)), np.float32)
        for i, m in enumerate(meta):
            ti = TISSUES.index(m["tissue"]); p, odi, fd = pix[ti]
            t2in, t2en = T2SETS[m["t2set"]]
            Sg[i] = pix_signal(pr, p[:, 0], p[:, 1], odi, fd, t2in, t2en, te, Es[ti])[0]
        flat = Sg.reshape(n * K, -1)
        snr = np.repeat(np.array([SNR[te] * m["scale"] for m in meta]), K)
        Ls = np.repeat(np.array([m["L"] for m in meta]), K)
        noisy = flat.copy()
        on = snr > 0
        if on.any():
            noisy[on] = SIM.noisy(np.random.default_rng([SEED, te]), flat[on], snr[on], Ls[on])
        vol = noisy.reshape(n, K, -1).transpose(1, 0, 2)[:, :, None, :]
        ndi, fwf, odi_ = SIM.amico_fit(str(OUT / "amico" / "W"), f"b1s_te{te}", vol, BVALS, BVECS, NTHREADS)
        maps = np.stack([ndi.T, fwf.T, odi_.T], -1)                       # (n, K, 3)
        b0 = noisy.reshape(n, K, -1)[:, :, BVALS <= 10].mean(-1)
        ok = np.isfinite(maps).all(-1) & (b0 > 0)
        roi = np.array([maps[i][ok[i]].mean(0) for i in range(n)])
        y, _ = lib.coords(roi[:, 0], roi[:, 1], roi[:, 2], None, "C3")
        res[:, j] = y
        print(json.dumps(dict(TE=te, n=n, seconds=round(time.time() - t0, 1), valid_frac=round(float(ok.mean()), 4))), flush=True)
    np.savez(OUT / "b1s_sim.npz", res=res, tes=np.array(TES), meta=np.array([json.dumps(m) for m in meta]))
    (OUT / "b1s_sim_meta.json").write_text(json.dumps(dict(script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(), V_S=dv, n=n, K=K,
                                                          snr=SNR, conds=CONDS, t2sets=T2SETS, tissues=TISSUES, seed=SEED), ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    print("끝", flush=True)


if __name__ == "__main__":
    main()
