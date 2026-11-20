"""Exploratory analysis c: repeats the D2 comparison with simulation centers matched to the observed apparent NDI and FWF levels by inverting config.apparent."""
import glob, hashlib, json, math, pathlib, shutil, sys, time, zlib
import numpy as np
X14 = pathlib.Path("<WORK>/stage1_work/drafts/x14")
OUT = pathlib.Path("<WORK>/remeasure/explore/rv/c_level")
CHUNKS = list(range(0, 27)) + list(range(81, 108))


def level_centres(cen, C):
    a = C.weights(92.0, (82.0, 82.0 / 1.6, 1000.0), 3600.0, (1100.0, 1100.0, C.T1_ISO))
    a_in, a_en, a_iso = a
    nu_p, f_p = cen[:, 0], cen[:, 1]
    nu = nu_p * a_en / (a_in * (1 - nu_p) + nu_p * a_en)
    A = nu * a_in + (1 - nu) * a_en
    f = f_p * A / (a_iso * (1 - f_p) + f_p * A)
    nu_b, f_b = C.apparent(nu, f, a)
    return np.c_[nu, f, cen[:, 2]], float(np.max(np.abs(np.c_[nu_b, f_b] - cen[:, :2]))), a


def run():
    sys.path.insert(0, str(X14 / "D2"))
    import d2_sim as D2
    S, SIM, lib = D2.S, D2.SIM, D2.lib
    OUT.mkdir(parents=True, exist_ok=True)
    P, meta, n_clip0 = D2.build_plan()
    obs = json.loads((D2.OUT / "d2_observed.json").read_text(encoding="utf-8"))
    cen = np.array([obs["subjects"][sub]["centre_ses001"][u] for sub, u, _ in meta], float)
    newc, back_err, a = level_centres(cen, S.C)
    lo = np.array([S.PIX_CLIP["NDI"][0], S.PIX_CLIP["FISO"][0], S.PIX_CLIP["ODI"][0]])
    hi = np.array([S.PIX_CLIP["NDI"][1], S.PIX_CLIP["FISO"][1], S.PIX_CLIP["ODI"][1]])
    newc_c = np.clip(newc, lo, hi)
    P["theta"][:, 3:6] = newc_c[P["centre_idx"]]
    v2 = D2.v2_plan(P)
    study = OUT / "amico"
    if not (study / "S" / "kernels").exists():
        shutil.copytree(S.WORK / "amico" / "S" / "kernels", study / "S" / "kernels")
    Ttab = {("S", 1.0): D2.TB.load("S", 1.0)}
    st = {"S": str(study / "S")}
    v1 = D2.v1_bitwise(P, Ttab, st)
    rec = dict(made=time.strftime("%Y-%m-%dT%H:%M:%S%z"), script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
               d2_sim_sha256=hashlib.sha256((X14 / "D2" / "d2_sim.py").read_bytes()).hexdigest(), weights_in_en_iso=[float(x) for x in a],
               inverse_check_max_abs=back_err, centres_clipped=int((newc_c != newc).any(1).sum()), centres_clipped_d2=n_clip0,
               centres_level=newc_c.tolist(), centres_obs=cen.tolist(), V1=v1, V2=v2, chunks=CHUNKS)
    (OUT / "c_verify.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(dict(V1=v1, V2_all=all(v2.values()), inverse_check=back_err, clipped=rec["centres_clipped"], clipped_d2=n_clip0), ensure_ascii=False), flush=True)
    for c in CHUNKS:
        f = OUT / f"c_c{c:04d}.npz"
        if f.exists():
            continue
        t0 = time.time()
        n = len(P["theta"]); idx = np.arange(c * D2.CSIZE, min(n, (c + 1) * D2.CSIZE))
        Pc = D2.sub(P, idx)
        with D2.Capture() as cap:
            SIM.chunk(Pc, np.arange(len(idx)), np.random.default_rng([D2.SEED_CHUNK, zlib.crc32(b"d2"), c]), Ttab, st, f"rvc_c{c:04d}", D2.NTHREADS)
        assert len(cap.got) == 2
        d = np.empty((len(D2.NEFFS), len(idx), 4))
        for ki, ne in enumerate(D2.NEFFS):
            x, fw = [], []
            for si_, (m, Sg) in enumerate(cap.got):
                roi = S.roi_noise(np.random.default_rng([D2.SEED_NEFF, c, si_]), m, Sg, np.full(len(idx), ne))
                y, _ = lib.coords(roi[:, 0], roi[:, 1], roi[:, 2], roi[:, 3:7], "C7")
                x.append(y[:, :3]); fw.append(roi[:, 1])
            d[ki, :, :3] = x[1] - x[0]
            d[ki, :, 3] = fw[1] - fw[0]
        lev = cap.got[0][0][:, 1]                                                        # Simulated apparent FWF (ROI mean m of site A · before roi_noise)
        np.savez(str(f) + ".tmp.npz", idx=idx, d=d, fwf_level=lev)
        pathlib.Path(str(f) + ".tmp.npz").replace(f)
        print(json.dumps(dict(chunk=c, seconds=round(time.time() - t0, 1))), flush=True)
    print("끝", flush=True)


def summary():
    import datetime
    O2 = pathlib.Path("<WORK>/remeasure/explore/d2")
    obs = json.loads((O2 / "d2_observed.json").read_text(encoding="utf-8"))
    P = np.load(O2 / "d2_plan.npz", allow_pickle=False)
    ver = json.loads((O2 / "d2_sim_verify.json").read_text(encoding="utf-8"))
    NEFFS, SC = ver["neffs"], ver["scen_snr"]
    n = len(P["theta"])
    D = np.full((len(NEFFS), n, 4), np.nan)
    LEV = np.full(n, np.nan)
    for f in sorted(glob.glob(str(OUT / "c_c*.npz"))):
        z = np.load(f); D[:, z["idx"]] = z["d"]; LEV[z["idx"]] = z["fwf_level"]
    units = [(sub, u) for sub, d in obs["subjects"].items() for u in d["units"]]
    obs_sd = np.array([obs["subjects"][sub]["diff_sd"][obs["subjects"][sub]["units"].index(u)] for sub, u in units])
    # SD of the observed absolute FWF difference — six repeat sessions of roi.json (re-access record)
    W = pathlib.Path("<WORK>/remeasure/onharmony")
    RR = X14 / "REREAD_LOG.jsonl"
    now = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    fw = {}
    for sub in obs["subjects"]:
        vals = []
        for i in range(1, 7):
            p = W / f"{sub}_ses-OXF1PRI00{i}" / "roi.json"
            b = p.read_bytes()
            with RR.open("a", encoding="utf-8") as fh:
                fh.write(json.dumps(dict(시각=now, 다시읽음=str(p), sha256=hashlib.sha256(b).hexdigest(), 바이트=len(b), 무엇="rv c — FWF 절대 차이 SD(관측)", 스크립트=__file__), ensure_ascii=False) + "\n")
            vals.append(json.loads(b)["주"])
        fw[sub] = vals
    obs_fwf = np.array([math.sqrt(2) * np.std([fw[sub][i][u]["FWF"] for i in range(6)], ddof=1) for sub, u in units])
    obs_lvl = np.array([obs["subjects"][sub]["centre_ses001"][u][1] for sub, u in units])
    fornix = np.array(["Fornix" in u for _, u in units])
    ci, si = P["centre_idx"], P["scen_snr"]
    k50 = NEFFS.index(50.0)
    out = dict(made=now, script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(), n_units=len(units), by_scen={})
    for sname in ("snr30", "measured"):
        s = SC.index(sname)
        sim_sd = np.array([[D[k, (si == s) & (ci == c)].std(0, ddof=1) for c in range(len(units))] for k in range(len(NEFFS))])   # (neff, unit, 4)
        lvl = np.array([np.nanmedian(LEV[(si == s) & (ci == c)]) for c in range(len(units))])
        res = {}
        for j, cn in enumerate(("x", "y", "ODI", "FWF_abs")):
            o = obs_sd[:, j] if j < 3 else obs_fwf
            r = o / sim_sd[k50, :, j]
            pooled = lambda k, msk=np.ones(len(units), bool): float(math.sqrt(np.mean(o[msk] ** 2)) / math.sqrt(np.mean(sim_sd[k, msk, j] ** 2)))  # noqa: E731
            res[cn] = dict(ratio_median=float(np.median(r)), ratio_q25=float(np.quantile(r, 0.25)), ratio_q75=float(np.quantile(r, 0.75)),
                           ratio_median_without_fornix=float(np.median(r[~fornix])), pooled_neff50=pooled(k50), pooled_neff50_without_fornix=pooled(k50, ~fornix),
                           contour_pooled={str(int(NEFFS[k])): round(pooled(k), 3) for k in range(len(NEFFS))})
        res["level_check"] = dict(sim_apparent_fwf_median=float(np.median(lvl)), obs_fwf_median=float(np.median(obs_lvl)),
                                  median_abs_diff=float(np.median(np.abs(lvl - obs_lvl))))
        out["by_scen"][sname] = res
    d2 = json.loads((O2 / "d2_repeat.json").read_text(encoding="utf-8"))
    out["d2_original_snr30_neff50"] = {cn: dict(ratio_median=d2["by_coord"][cn]["ratio_median"], pooled=d2["by_coord"][cn]["pooled_ratio_snr30_neff50"],
                                                pooled_without_fornix=d2["by_coord"][cn]["pooled_ratio_without_fornix"]) for cn in ("x", "y", "ODI")}
    yp = out["by_scen"]["snr30"]["y"]["pooled_neff50"]
    out["reading_y"] = ("y 초과는 수준 불일치 · logit 증폭으로 설명됨" if yp <= 1.2 else "수준을 맞춰도 y 초과가 남음" if yp >= 2.0 else "부분")
    (OUT.parent / "rv_c_level.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({s: {k: (round(v["ratio_median"], 3), round(v["pooled_neff50"], 3)) for k, v in r.items() if k != "level_check"} for s, r in out["by_scen"].items()}, ensure_ascii=False))
    print("level", {s: r["level_check"] for s, r in out["by_scen"].items()}, "original", out["d2_original_snr30_neff50"], "reading", out["reading_y"])


if __name__ == "__main__":
    dict(run=run, summary=summary).get(sys.argv[1] if len(sys.argv) > 1 else "", lambda: sys.exit("쓰임: run|summary"))()
