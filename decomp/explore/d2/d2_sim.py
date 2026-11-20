"""Simulates d2: the scatter the judge assumes for the same scanner twice (Δ = 0, L0, OXF1PRI acquisition), calling the frozen simulator unchanged."""
import hashlib, json, pathlib, shutil, sys, time, zlib
import numpy as np
REPO = pathlib.Path("<REPO>")
sys.path.insert(0, str(REPO / "decomp" / "s4"))
import s4cfg as S  # noqa: E402
import simulate as SIM  # noqa: E402
import tables as TB  # noqa: E402
sys.path.insert(0, str(S.C.RESULTS.parent / "stage1"))           # Same place as in simulate.run (lib · synth)
import lib  # noqa: E402
OUT = pathlib.Path("<WORK>/remeasure/explore/d2")
SEED_PLAN, SEED_CHUNK, SEED_NEFF = 20261121, 20261122, 20261123
NREP, CSIZE, NTHREADS = 500, 1000, 16
NEFFS = (16.0, 25.0, 50.0, 100.0, 200.0)
SCEN_SNR = ("snr30", "snr20", "snr50", "measured")
OXF1PRI = dict(design="S", TE=92.0, TR=3600.0, coil=32)          # decomp/s4/plan.py ONH["OXF1PRI"] = ("S", (92, 92), (3.6, 3.6), 32)
assert S.DESIGNS["S"] == "P2a" and S.COIL_L[32] == 4


def build_plan():
    obs = json.loads((OUT / "d2_observed.json").read_text(encoding="utf-8"))
    cen, meta = [], []
    for sub, d in obs["subjects"].items():
        for u in d["units"]:
            cen.append(d["centre_ses001"][u]); meta.append((sub, u, d["snr_subject"]))
    cen = np.array(cen, float)
    lo = np.array([S.PIX_CLIP["NDI"][0], S.PIX_CLIP["FISO"][0], S.PIX_CLIP["ODI"][0]])
    hi = np.array([S.PIX_CLIP["NDI"][1], S.PIX_CLIP["FISO"][1], S.PIX_CLIP["ODI"][1]])
    cen_c = np.clip(cen, lo, hi)
    n_clip = int((cen_c != cen).any(1).sum())
    rows = [(si, ci, j) for si in range(len(SCEN_SNR)) for ci in range(len(cen)) for j in range(NREP)]
    n = len(rows)
    si, ci, j = (np.array(v) for v in zip(*rows))
    snr = np.array([30.0, 20.0, 50.0, np.nan])[si]
    snr[si == 3] = np.array([m[2] for m in meta])[ci[si == 3]]
    K9 = [(g, sc) for g in S.GENERATORS for sc in S.SCENARIOS]
    gen = np.array([K9[k % 9][0] for k in j]); scen = np.array([K9[k % 9][1] for k in j])
    rng = np.random.default_rng(SEED_PLAN)
    t2A, t1A = S.sample_relax_A(rng, n, scen)
    lv, p = S.sample_level(rng, n, np.full(n, "L0"))
    t2B, t1B = S.relax_B(rng, lv, p, t2A, t1A)
    theta = np.zeros((n, 6)); theta[:, 3:6] = cen_c[ci]
    P = dict(theta=theta, t2A=t2A, t1A=t1A, t2B=t2B, t1B=t1B, dpar=rng.choice(S.DPAR_TRUE, n),
             fdir_A=rng.integers(0, len(S.C.FIBRE_DIRS), n), fdir_B=rng.integers(0, len(S.C.FIBRE_DIRS), n), bscale=np.ones(n),
             gen=gen, scen=scen, level=lv,
             acq_TE_A=np.full(n, OXF1PRI["TE"]), acq_TE_B=np.full(n, OXF1PRI["TE"]), acq_TR_A=np.full(n, OXF1PRI["TR"]), acq_TR_B=np.full(n, OXF1PRI["TR"]),
             acq_design_A=np.full(n, OXF1PRI["design"]), acq_design_B=np.full(n, OXF1PRI["design"]),
             acq_coil_A=np.full(n, OXF1PRI["coil"]), acq_coil_B=np.full(n, OXF1PRI["coil"]),
             acq_snr_A=snr, acq_snr_B=snr.copy(), acq_neff=np.full(n, 50.0), scen_snr=si, centre_idx=ci, rep=j)
    return P, meta, n_clip


def v2_plan(P):
    ok = dict(delta0=bool(np.all(P["theta"][:, :3] == 0)), L0=bool(np.all(P["level"] == "L0")),
              relax_same=bool(np.array_equal(P["t2A"], P["t2B"]) and np.array_equal(P["t1A"], P["t1B"])),
              TE=bool(np.all(P["acq_TE_A"] == 92.0) and np.all(P["acq_TE_B"] == 92.0)), TR=bool(np.all(P["acq_TR_A"] == 3600.0) and np.all(P["acq_TR_B"] == 3600.0)),
              design=bool(np.all(P["acq_design_A"] == "S") and np.all(P["acq_design_B"] == "S")), coil=bool(np.all(P["acq_coil_A"] == 32) and np.all(P["acq_coil_B"] == 32)),
              bscale=bool(np.all(P["bscale"] == 1.0)), snr_finite=bool(np.isfinite(P["acq_snr_A"]).all()), k9_balanced=bool(len(set(zip(P["gen"], P["scen"]))) == 9),
              dpar_in_table=bool(all(float(v) in TB.DPAR_T for v in np.unique(P["dpar"]))))
    if not all(ok.values()):
        sys.exit(f"⛔ V2 계획 대조 실패: {ok}")
    return ok


def setup():
    study_root = OUT / "amico"
    k_src, k_dst = S.WORK / "amico" / "S" / "kernels", study_root / "S" / "kernels"
    if not k_dst.exists():
        shutil.copytree(k_src, k_dst)                            # Do not write into the S4 working folder — copy the kernel to a new study folder
    Ttab = {("S", 1.0): TB.load("S", 1.0)}
    return Ttab, {"S": str(study_root / "S")}


class Capture:
    """Wraps S.roi_noise — holds on to (m, Σ̂) and calls the original function unchanged with the same arguments and the same random numbers."""
    def __init__(self):
        self.orig, self.got = S.roi_noise, []

    def __enter__(self):
        def wrap(rng, m, Sg, neff, K=S.K_PIX):
            self.got.append((np.array(m, copy=True), np.array(Sg, copy=True)))
            return self.orig(rng, m, Sg, neff, K)
        S.roi_noise = wrap
        return self

    def __exit__(self, *a):
        S.roi_noise = self.orig


def sub(P, idx):
    return {k: v[idx] for k, v in P.items()}


def v1_bitwise(P, Ttab, study):
    idx = np.arange(40)
    Pt = sub(P, np.r_[np.arange(20), len(P["theta"]) // 2 + np.arange(20)])   # 20 from each of the two SNR batches
    r1 = SIM.chunk(Pt, idx, np.random.default_rng([SEED_CHUNK, 1]), Ttab, study, "d2v1a", NTHREADS)
    with Capture() as cap:
        r2 = SIM.chunk(Pt, idx, np.random.default_rng([SEED_CHUNK, 1]), Ttab, study, "d2v1b", NTHREADS)
    res = dict(x_c3=bool(np.array_equal(r1["x_c3"], r2["x_c3"], equal_nan=True)), x_c7=bool(np.array_equal(r1["x_c7"], r2["x_c7"], equal_nan=True)),
               captured=len(cap.got), finite_A=int(np.sum(r1["finite_A"])), finite_B=int(np.sum(r1["finite_B"])))
    if not (res["x_c3"] and res["x_c7"] and res["captured"] == 2):
        sys.exit(f"⛔ V1 비트 대조 실패: {res}")
    return res


def run_chunk(P, c, Ttab, study):
    f = OUT / f"d2_sim_c{c:04d}.npz"
    if f.exists():
        return False
    n = len(P["theta"]); idx = np.arange(c * CSIZE, min(n, (c + 1) * CSIZE))
    Pc = sub(P, idx)
    loc = np.arange(len(idx))
    with Capture() as cap:
        r = SIM.chunk(Pc, loc, np.random.default_rng([SEED_CHUNK, zlib.crc32(b"d2"), c]), Ttab, study, f"d2_c{c:04d}", NTHREADS)
    assert len(cap.got) == 2
    d = np.empty((len(NEFFS), len(idx), 3))
    for ki, ne in enumerate(NEFFS):
        x = []
        for si_, (m, Sg) in enumerate(cap.got):
            roi = S.roi_noise(np.random.default_rng([SEED_NEFF, c, si_]), m, Sg, np.full(len(idx), ne))   # The five N_eff values use the same seed (common random numbers)
            y, _ = lib.coords(roi[:, 0], roi[:, 1], roi[:, 2], roi[:, 3:7], "C7")
            x.append(y[:, :3])
        d[ki] = x[1] - x[0]
    np.savez(str(f) + ".tmp.npz", idx=idx, d=d, x_c3_primary=r["x_c3"], finite_A=r["finite_A"], finite_B=r["finite_B"])
    pathlib.Path(str(f) + ".tmp.npz").replace(f)
    return True


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "verify"
    P, meta, n_clip = build_plan()
    v2 = v2_plan(P)
    Ttab, study = setup()
    v1 = v1_bitwise(P, Ttab, study)
    rec = dict(made=time.strftime("%Y-%m-%dT%H:%M:%S%z"), script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
               simulate_sha256=hashlib.sha256((REPO / "decomp/s4/simulate.py").read_bytes()).hexdigest(),
               s4cfg_sha256=hashlib.sha256((REPO / "decomp/s4/s4cfg.py").read_bytes()).hexdigest(),
               n_cases=len(P["theta"]), n_centres=len(meta), centres_clipped=n_clip, V2=v2, V1=v1, neffs=NEFFS, scen_snr=SCEN_SNR,
               measured_snr={m[0]: m[2] for m in meta}, seeds=dict(plan=SEED_PLAN, chunk=SEED_CHUNK, neff=SEED_NEFF), nrep=NREP)
    (OUT / "d2_sim_verify.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(dict(V1=v1, V2_all=all(v2.values()), centres_clipped=n_clip, n_cases=rec["n_cases"]), ensure_ascii=False), flush=True)
    if mode != "run":
        return
    np.savez(OUT / "d2_plan.npz", **{k: v for k, v in P.items()})
    nchunk = (len(P["theta"]) + CSIZE - 1) // CSIZE
    for c in range(nchunk):
        t0 = time.time()
        if run_chunk(P, c, Ttab, study):
            print(json.dumps(dict(chunk=c, of=nchunk, seconds=round(time.time() - t0, 1))), flush=True)
    print("끝", flush=True)


if __name__ == "__main__":
    main()
