"""Runs d3 v3: continuous r profiles, interval union over six d∥ × T2_iso variants, and a real-noise known-answer test with coverage and reading gates."""
import datetime, hashlib, json, math, pathlib, sys
import numpy as np
from multiprocessing import get_context
from scipy.optimize import brentq, minimize_scalar
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / "B1S"))
import d3_kat as K1  # noqa: E402
import d3_real as DR  # noqa: E402 — inner(D, sig, t2in, r, x0, kw) · form for reading measured data
import b1s_sim as B  # noqa: E402 — Wong direction Proto · pix_signal (V-S 0.0) · SNR(TE)
import s4cfg as S  # noqa: E402
import simulate as SIM  # noqa: E402
import synth  # noqa: E402
import config as C  # noqa: E402
OUT = pathlib.Path("<WORK>/remeasure/explore/d3")
TES = B.TES
VARIANTS = {f"dpar{d}|t2iso{int(t)}": dict(dpar=d, t2iso=t) for d in (1.5, 1.7, 2.0) for t in (1000.0, 2000.0)}
TRUTHS, SCALES, NG, K, NROI, SEED = (1.6, 2.64), (1.0, 2.0), 20, 200, 27, 20261171
CO = np.exp(np.linspace(np.log(0.5), np.log(6.0), 13))
READ = lambda lo, hi: ("T2 비로는 B1의 얕은 각도가 설명되지 않음" if (lo <= 1.6 <= hi and not lo <= 2.6 <= hi) else  # noqa: E731
                       "T2 비가 문헌보다 큼" if (lo <= 2.6 <= hi and not lo <= 1.6 <= hi) else "판정 보류")
LIB = []


def build_lib():
    pr = synth.Proto("WONG_SUB01_TE92")
    rng = np.random.default_rng([SEED, 1])
    lo, hi = np.array([S.PIX_CLIP["NDI"][0], S.PIX_CLIP["FISO"][0]]), np.array([S.PIX_CLIP["NDI"][1], S.PIX_CLIP["FISO"][1]])
    for ti in C.WM:
        nu, odi, f = C.TISSUE[ti]
        for fd in range(len(C.FIBRE_DIRS)):
            p = np.clip(np.array([nu, f])[None] + rng.normal(0, 1, (K, 2)) * S.PIX_SD[[0, 1]], lo, hi)
            _, E = B.pix_signal(pr, p[:, 0], p[:, 1], odi, fd, 82.0, 52.0, 62.0)
            LIB.append(dict(p=p, odi=odi, fd=fd, E=E))
    return B.vcheck(pr)


LCO = [float(v) for v in np.log(CO)]


def prof_cont(Ds, sig, kw):
    """(item ㄱ, the first Korean-lettered item) · fix of §15-2 — F(lr) is a deterministic function of log r. The 13 coarse grid points start in ascending order from the solution of the previous point (the first point (ν, f) = (0.6, 0.1)) · any other lr starts from the solution of the nearest coarse point · an lr already computed uses the remembered value."""
    memo, start = {}, []

    def inner_all(lr, x0s):
        tot, xs = 0.0, []
        for k, D in enumerate(Ds):
            cache = {}

            def g(lt):
                v, x = DR.inner(D, sig, math.exp(lt), math.exp(lr), x0s[k], kw); cache[lt] = (v, x); return v
            m = minimize_scalar(g, bounds=(math.log(20.0), math.log(300.0)), method="bounded", options=dict(xatol=1e-3))
            v, x = cache[m.x] if m.x in cache else DR.inner(D, sig, math.exp(m.x), math.exp(lr), x0s[k], kw)
            tot += v; xs.append(x)
        return tot, xs

    x0 = [np.tile([0.6, 0.1], D.shape[0]) for D in Ds]
    for lr in LCO:
        v, x0 = inner_all(lr, x0)
        memo[lr] = v; start.append(x0)

    def F(lr):
        lr = float(lr)
        if lr not in memo:
            j = int(np.argmin(np.abs(np.array(LCO) - lr)))
            memo[lr] = inner_all(lr, start[j])[0]
        return memo[lr]

    vals = np.array([memo[lr] for lr in LCO])
    i = int(np.argmin(vals))
    m = minimize_scalar(F, bounds=(LCO[max(i - 1, 0)], LCO[min(i + 1, len(LCO) - 1)]), method="bounded", options=dict(xatol=1e-4))
    lrhat = float(m.x) if F(m.x) < vals[i] else LCO[i]
    rmin = F(lrhat)
    h = lambda lr: F(lr) - rmin - 3.84  # noqa: E731
    los = [lr for lr in LCO if lr < lrhat and F(lr) - rmin > 3.84]
    his = [lr for lr in LCO if lr > lrhat and F(lr) - rmin > 3.84]
    lo = math.exp(brentq(h, max(los), lrhat, xtol=1e-4)) if los else None
    hi = math.exp(brentq(h, lrhat, min(his), xtol=1e-4)) if his else None
    return dict(rhat=math.exp(lrhat), lo=lo if lo is not None else float(CO[0]), hi=hi if hi is not None else float(CO[-1]),
                open_lo=lo is None, open_hi=hi is None, rss_min=float(rmin), rss_coarse=[float(v) for v in vals], n_eval=len(memo))


def hull(res):
    lo, hi = min(v["lo"] for v in res.values()), max(v["hi"] for v in res.values())
    return lo, hi


def worker(args):
    gi, r_true, scale, sig = args
    D = make_D(gi, r_true, scale)
    want = "T2 비로는 B1의 얕은 각도가 설명되지 않음" if r_true == 1.6 else "T2 비가 문헌보다 큼"
    try:
        res = {name: prof_cont([D], sig, kw) for name, kw in VARIANTS.items()}
    except Exception as e:  # §15-2 (4): the exception batch is «not covered · not read as matching»
        return gi, dict(r_true=r_true, scale=scale, error=repr(e), covers=False, correct=False, reading="예외")
    lo, hi = hull(res)
    return gi, dict(r_true=r_true, scale=scale, variants=res, hull=[lo, hi], covers=bool(lo <= r_true <= hi), reading=READ(lo, hi), correct=READ(lo, hi) == want)


def make_D(gi, r_true, scale, noisy=True):
    """ROI mean of one bundle (27, 8 TE, 2 shells) — voxel signal → (noise) → ext_coords → mean over finite voxels."""
    rng = np.random.default_rng([SEED, 2, gi])
    ent = rng.integers(0, len(LIB), NROI)
    D = np.empty((NROI, len(TES), 2))
    for j, te in enumerate(TES):
        sg = np.concatenate([B.pix_signal(None, LIB[e]["p"][:, 0], LIB[e]["p"][:, 1], LIB[e]["odi"], LIB[e]["fd"], 82.0, 82.0 / r_true, te, LIB[e]["E"])[0] for e in ent]).astype(np.float32)
        nz = SIM.noisy(np.random.default_rng([SEED, 3, gi, te]), sg, np.full(len(sg), B.SNR[te] * scale), np.full(len(sg), 4)) if noisy else sg
        s = synth.ext_coords(nz.astype(float), B.BVALS, B.BVECS, 700.0)[:, 2:4].reshape(NROI, K, 2)
        ok = np.isfinite(s).all(-1)
        D[:, j] = np.array([s[i][ok[i]].mean(0) for i in range(NROI)])
    return D


def real_fit(sig):
    data = {}
    now = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    for sub in ("sub-01", "sub-02", "sub-03"):
        own = set(json.loads((DR.WONG / sub / "DONE.json").read_text(encoding="utf-8")).get("따로_정합", []))
        tabs = []
        for s in DR.TES:
            p = (DR.WONG / sub / f"own_{s}" if s in own else DR.WONG / sub) / f"roi_{s}.json"
            b = p.read_bytes()
            with DR.LOG.open("a", encoding="utf-8") as fh:
                fh.write(json.dumps(dict(시각=now, 다시읽음=str(p), sha256=hashlib.sha256(b).hexdigest(), 바이트=len(b), 무엇="D3-v3 실측 적합(묶음)", 스크립트=__file__), ensure_ascii=False) + "\n")
            tabs.append(json.loads(b)["주"])
        units = sorted(u for u in tabs[0] if all(u in t and all(t[u].get(k) is not None and math.isfinite(t[u][k]) for k in ("s_low", "s_2000")) for t in tabs))
        data[sub] = np.array([[[t[u]["s_low"], t[u]["s_2000"]] for t in tabs] for u in units], float)
    jobs = [(name, key, Ds) for name, kw in VARIANTS.items() for key, Ds in list({s: [data[s]] for s in data}.items()) + [("overall", list(data.values()))]]
    with get_context("fork").Pool(16) as pool:
        out = pool.map(_real_job, [(name, key, Ds, sig) for name, key, Ds in jobs])
    res = {}
    for name, key, r in out:
        res.setdefault(key, {})[name] = r
    out = {}
    for key, v in res.items():
        n_units = {s: int(data[s].shape[0]) for s in data} if key == "overall" else int(data[key].shape[0])
        bad = [n for n, x in v.items() if "error" in x]
        out[key] = dict(variants=v, n_units=n_units, errors=bad) if bad else dict(variants=v, hull=list(hull(v)), reading=READ(*hull(v)), n_units=n_units)
    return out


def _real_job(a):
    name, key, Ds, sig = a
    try:
        return name, key, prof_cont(Ds, sig, VARIANTS[name])
    except Exception as e:  # §15-2 (5)
        return name, key, dict(error=repr(e))


def main():
    rec = dict(made=datetime.datetime.now().astimezone().isoformat(timespec="seconds"), script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
               plan_sha256_at_definition="5b8085ca8d248507…", variants=list(VARIANTS), truths=TRUTHS, scales=SCALES, n_groups=NG)
    rec["V_S"] = build_lib()
    if rec["V_S"] > 1e-12:
        sys.exit("⛔ V-S")
    sig, nd = K1.reread_sigma()
    rec["sigma_shell"] = dict(b700=float(sig[0]), b2000=float(sig[1]), n_diffs=nd)
    jobs, gi = [], 0
    for r in TRUTHS:
        for sc in SCALES:
            for _ in range(NG):
                jobs.append((gi, r, sc, sig)); gi += 1
    res = {}
    with get_context("fork").Pool(16) as pool, (OUT / "d3_v3_groups.jsonl").open("w", encoding="utf-8") as fh:
        for gi, x in pool.imap_unordered(worker, jobs):
            res[gi] = x
            fh.write(json.dumps(dict(gi=gi, **x), ensure_ascii=False) + "\n"); fh.flush()
    rec["errors"] = sorted(gi for gi, x in res.items() if "error" in x)
    rec["KAT"] = {}
    for r in TRUTHS:
        for sc in SCALES:
            v = [x for x in res.values() if x["r_true"] == r and x["scale"] == sc]
            rec["KAT"][f"r{r}|snr×{sc}"] = dict(n=len(v), coverage_hull=float(np.mean([x["covers"] for x in v])), correct_reading=float(np.mean([x["correct"] for x in v])),
                                                n_error=sum(1 for x in v if "error" in x),
                                                rhat_main_median=float(np.median([x["variants"]["dpar1.7|t2iso1000"]["rhat"] for x in v if "error" not in x])),
                                                hull_median=[float(np.median([x["hull"][0] for x in v if "error" not in x])), float(np.median([x["hull"][1] for x in v if "error" not in x]))],
                                                open_hi_main=sum(1 for x in v if "error" not in x and x["variants"]["dpar1.7|t2iso1000"]["open_hi"]),
                                                readings={k: sum(1 for x in v if x["reading"] == k) for k in {x["reading"] for x in v}})
    gate = all(g["coverage_hull"] >= 0.85 and g["correct_reading"] >= 0.8 for g in rec["KAT"].values())
    rec["gate_for_real_fit"] = dict(passed=bool(gate), rule="두 SNR 배수 모두 · 참값 1.6 · 2.64 각각 합 구간 덮음 ≥ 0.85 · 맞게 읽힌 비율 ≥ 0.8(계획서 §15)")
    (OUT / "d3_v3.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(dict(KAT={k: (v["coverage_hull"], v["correct_reading"], round(v["rhat_main_median"], 3), [round(x, 3) for x in v["hull_median"]]) for k, v in rec["KAT"].items()}, gate=gate), ensure_ascii=False), flush=True)
    if not gate:
        print("넘어가는 조건 미달 — 실측 없음", flush=True); return
    rec["real"] = real_fit(sig)
    (OUT / "d3_v3.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({k: (v["hull"], v["reading"]) for k, v in rec["real"].items()}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
