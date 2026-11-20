"""Runs the d3 known-answer test without measurement fits: spherical-mean formula check, implementation, failure controls, model adequacy and identification, in envs/reference."""
import datetime, hashlib, json, math, os, pathlib, sys
import numpy as np
from multiprocessing import get_context
from scipy.optimize import least_squares
from scipy.special import erf
REPO = pathlib.Path("<REPO>")
sys.path.insert(0, str(REPO / "decomp" / "stage1"))
import config as C  # noqa: E402
OUT = pathlib.Path("<WORK>/remeasure/explore/d3")
WONG = pathlib.Path("<WORK>/remeasure/wong")
LOG = pathlib.Path("<WORK>/stage1_work/drafts/x14/REREAD_LOG.jsonl")
TES = np.array([62.0, 72.0, 82.0, 92.0, 102.0, 112.0, 122.0, 132.0])
BMM = np.array([700.0, 2000.0])                      # s/mm² — Wong shell (rconfig.LOW_SHELL["wong"] = 700)
TR, T2ISO, DPAR_FIT = 5800.0, 1000.0, 1.7
LB, UB = np.array([0.01, 0.001, 20.0, 0.5]), np.array([0.99, 0.8, 300.0, 5.0])
STARTS = (0.8, 1.2, 1.6, 2.2, 3.0)
RGRID = np.exp(np.linspace(np.log(0.6), np.log(4.0), 41))
NREP_C, SEED = 30, 20261131
TOL_A = 1e-3


def smean(b, dpar, dper):
    """Spherical mean (stick · cylinder · free water) — b ms/µm² · d µm²/ms · dper is a scalar or a (TE,) array."""
    b = np.asarray(b, float)
    dper = np.asarray(dper, float)
    stick = np.sqrt(np.pi / (4 * b * dpar)) * erf(np.sqrt(b * dpar))
    dd = (dpar - dper)[..., None]
    cyl = np.exp(-b * dper[..., None]) * np.where(dd > 1e-12, np.sqrt(np.pi / (4 * b * np.maximum(dd, 1e-12))) * erf(np.sqrt(b * np.maximum(dd, 1e-12))), 1.0)
    return stick, cyl, np.exp(-b * C.DISO)


def model_s(th, dpar=DPAR_FIT, t2iso=T2ISO, gen="T-in"):
    """s(TE, b) = log[Σ w_c A_c(b) / Σ w_c] · (8, 2)."""
    nu, f, t2in, r = th
    w_in = (1 - f) * nu * np.exp(-TES / t2in)
    w_en = (1 - f) * (1 - nu) * np.exp(-TES / (t2in / r))
    w_iso = f * np.exp(-TES / t2iso)
    nua = w_in / (w_in + w_en)
    dper = np.array([C.dperp(gen, dpar, nu, x) for x in nua])
    stick, cyl, iso = smean(BMM / 1000.0, dpar, dper)
    S = (w_in[:, None] * stick[None, :] + w_en[:, None] * cyl + w_iso[:, None] * iso[None, :]) / (w_in + w_en + w_iso)[:, None]
    return np.log(S)


def resid(th, data, sig, kw):
    return ((data - model_s(th, **kw)) / sig[None, :]).ravel()


def fit(data, sig, kw):
    best = None
    for r0 in STARTS:
        res = least_squares(resid, np.array([0.6, 0.1, 80.0, r0]), bounds=(LB, UB), args=(data, sig, kw), xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=4000)
        if best is None or res.cost < best.cost:
            best = res
    return best


def profile(data, sig, best, kw):
    rss, prev = [], best.x[:3].copy()
    for r in RGRID:
        f = lambda t3: resid(np.r_[t3, r], data, sig, kw)  # noqa: E731
        cand = [least_squares(f, x0, bounds=(LB[:3], UB[:3]), xtol=1e-10, ftol=1e-10, gtol=1e-10, max_nfev=2000) for x0 in (prev, best.x[:3])]
        res = min(cand, key=lambda z: z.cost)
        rss.append(2 * res.cost)
        prev = res.x
    return np.array(rss)


def ci_of(rss, best):
    m = min(2 * best.cost, rss.min())
    inside = rss - m <= 3.84
    return float(RGRID[inside].min()), float(RGRID[inside].max()), bool(inside[0] or inside[-1])


def truths():
    t2k = [k for k, (t, r, i) in enumerate(C.T2_SETS) if i == 1000.0]
    return [(ti, k) for ti in C.WM for k in t2k]


def gen_A(ti, k):
    nu, odi, f = C.TISSUE[ti]
    t2in, r, _ = C.T2_SETS[k]
    return model_s((nu, f, t2in, r))


def gen_B(pr, ti, k, t1i, dpar, idx):
    import synth  # noqa: E402 — project generator (AMICO Watson synthesis)
    rows = []
    for te in TES:
        cell = dict(tissue=ti, t2set=k, TE=float(te), TR=TR, t1set=t1i, dpar=dpar, gen="T-in", fdir=idx % len(C.FIBRE_DIRS))
        S, _, _, _ = synth.noiseless(pr, cell)
        rows.append(synth.ext_coords(np.asarray(S, float)[None, :], pr.bvals, pr.bvecs, 700.0)[0, 2:4])
    return np.array(rows, float)


def worker_C(args):
    i, clean, r_true, sig, kw = args
    rng = np.random.default_rng([SEED, i])
    out = []
    for _ in range(NREP_C):
        data = clean + rng.normal(0.0, 1.0, clean.shape) * sig[None, :]
        best = fit(data, sig, kw)
        lo, hi, touch = ci_of(profile(data, sig, best, kw), best)
        out.append((float(best.x[3]), lo, hi, touch, lo <= r_true <= hi))
    return i, out


def reread_sigma():
    now = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    d = {700: [], 2000: []}
    for sub in ("sub-01", "sub-02", "sub-03"):
        pairs = ((WONG / sub / "roi_TE62.json", WONG / sub / "own_TE62R2" / "roi_TE62R2.json"),
                 (WONG / sub / "own_TE132" / "roi_TE132.json", WONG / sub / "own_TE132R2" / "roi_TE132R2.json"))
        for a, b in pairs:
            js = []
            for p in (a, b):
                raw = p.read_bytes()
                with LOG.open("a", encoding="utf-8") as fh:
                    fh.write(json.dumps(dict(시각=now, 다시읽음=str(p), sha256=hashlib.sha256(raw).hexdigest(), 바이트=len(raw), 무엇="D3 정답 시험 — 반복 수준 잡음 SD(셸별)", 스크립트=__file__), ensure_ascii=False) + "\n")
                js.append(json.loads(raw)["주"])
            for u in js[0]:
                if u in js[1]:
                    for key, b_ in (("s_low", 700), ("s_2000", 2000)):
                        x, y = js[0][u].get(key), js[1][u].get(key)
                        if x is not None and y is not None and math.isfinite(x) and math.isfinite(y):
                            d[b_].append(y - x)
    sig = np.array([math.sqrt(np.mean(np.square(d[700])) / 2), math.sqrt(np.mean(np.square(d[2000])) / 2)])
    return sig, {700: len(d[700]), 2000: len(d[2000])}


def criterion_A(est, tr):
    bad = [i for i, (e, t) in enumerate(zip(est, tr)) if not (abs(e[3] - t[3]) <= TOL_A and abs(e[0] - t[0]) <= TOL_A and abs(e[1] - t[1]) <= TOL_A)]
    return len(bad) == 0, len(bad)


def main():
    rec = dict(made=datetime.datetime.now().astimezone().isoformat(timespec="seconds"), script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
               design="drafts/x14/DESIGN_D2_D3KAT.md", rgrid=[float(RGRID[0]), float(RGRID[-1]), len(RGRID)], tol_A=TOL_A)
    # V5 — check against the spherical-average formula (dense direction average · 20,000 directions)
    g = C.fibonacci(20000)
    n = np.array([0.0, 0.0, 1.0])
    cos2 = (g @ n) ** 2
    v5 = []
    for b in BMM / 1000.0:
        st, cy, _ = smean(np.array([b]), 1.7, 0.6)
        v5.append(abs(float(np.mean(np.exp(-b * 1.7 * cos2))) - float(st[0])))
        v5.append(abs(float(np.mean(np.exp(-b * 0.6 - b * (1.7 - 0.6) * cos2))) - float(np.ravel(cy)[0])))
    rec["V5_max_abs_diff"] = max(v5)
    if max(v5) > 1e-4:
        sys.exit(f"⛔ V5 실패 {v5}")
    T = truths()
    tr = [(C.TISSUE[ti][0], C.TISSUE[ti][2], C.T2_SETS[k][0], C.T2_SETS[k][1]) for ti, k in T]
    sig1 = np.array([1.0, 1.0])
    # A-1 — implementation
    est = [fit(gen_A(ti, k), sig1, {}).x for ti, k in T]
    okA, nbadA = criterion_A(est, tr)
    rec["A1"] = dict(passed=okA, n_truths=len(T), n_violating=nbadA, max_abs_dr=float(max(abs(e[3] - t[3]) for e, t in zip(est, tr))))
    if not okA:
        sys.exit(f"⛔ A-1 실패 {rec['A1']}")
    # A-2 — failure control (deliberately wrong input)
    perm = np.random.default_rng(SEED).permutation(len(TES))
    assert not np.array_equal(perm, np.arange(len(TES)))
    rec["A2"] = {}
    for name, tf in (("TE 순서 섞음", lambda D: D[perm]), ("셸 바꿈", lambda D: D[:, ::-1])):
        e2 = [fit(tf(gen_A(ti, k)), sig1, {}).x for ti, k in T]
        ok2, nb2 = criterion_A(e2, tr)
        rec["A2"][name] = dict(criterion_passed=ok2, n_violating=nb2, perm=perm.tolist() if "TE" in name else None)
        if ok2:
            sys.exit(f"⛔ A-2 «{name}» 가 실패하지 않았다 — 시험이 틀린 입력을 못 잡는다")
    # B — project generator · noise 0
    import synth  # noqa: E402
    pr = synth.Proto("P1")
    rec["P1"] = dict(n_vol=int(len(pr.bvals)), b_values=sorted({float(b) for b in pr.bvals}))
    t1 = {name: C.T1_SETS.index(v) for name, v in (("1100|1.0", (1100.0, 1.0)), ("1100|1.3", (1100.0, 1.3)))}
    rows = dict(B0=("1100|1.0", 1.7, {}), B1_T1ratio13=("1100|1.3", 1.7, {}), B2_data_dpar20=("1100|1.0", 2.0, {}),
                B3_fit_dpar15=("1100|1.0", 1.7, dict(dpar=1.5)), B4_fit_dpar20=("1100|1.0", 1.7, dict(dpar=2.0)))
    cleanB, rec["B"] = {}, {}
    for name, (t1n, dp, kw) in rows.items():
        key = (t1n, dp)
        if key not in cleanB:
            cleanB[key] = [gen_B(pr, ti, k, t1[t1n], dp, i) for i, (ti, k) in enumerate(T)]
        dr = np.array([fit(D, sig1, kw).x[3] - t[3] for D, t in zip(cleanB[key], tr)])
        rec["B"][name] = dict(median=float(np.median(dr)), q25=float(np.quantile(dr, 0.25)), q75=float(np.quantile(dr, 0.75)), max_abs=float(np.max(np.abs(dr))))
    # C — identification (generator B main line + measured repeat-level noise)
    sig, nd = reread_sigma()
    rec["sigma_shell"] = dict(b700=float(sig[0]), b2000=float(sig[1]), n_diffs=nd)
    jobs = [(i, cleanB[("1100|1.0", 1.7)][i], tr[i][3], sig, {}) for i in range(len(T))]
    with get_context("fork").Pool(8) as pool:
        resC = dict(pool.map(worker_C, jobs))
    allr = [x for i in range(len(T)) for x in resC[i]]
    touch = float(np.mean([x[3] for x in allr]))
    width = [x[2] - x[1] for x in allr]
    cover = float(np.mean([x[4] for x in allr]))
    byr = {}
    for i, t in enumerate(tr):
        byr.setdefault(round(t[3], 2), []).extend(resC[i])
    rec["C"] = dict(n_cases=len(allr), touch_fraction=touch, coverage=cover, width_median=float(np.median(width)),
                    width_q25=float(np.quantile(width, 0.25)), width_q75=float(np.quantile(width, 0.75)),
                    by_r={str(r): dict(touch=float(np.mean([x[3] for x in v])), coverage=float(np.mean([x[4] for x in v])),
                                       rhat_median=float(np.median([x[0] for x in v])), width_median=float(np.median([x[2] - x[1] for x in v]))) for r, v in sorted(byr.items())},
                    rule="«식별 안 됨» = 사전 범위(r ≤ 1.8) 참값에서 95 % 구간이 격자 끝(0.6 · 4.0)에 닿는 경우 ≥ 50 %",
                    identified=bool(touch < 0.5))
    (OUT / "d3_kat.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(dict(V5=rec["V5_max_abs_diff"], A1=rec["A1"], A2={k: v["n_violating"] for k, v in rec["A2"].items()}, B0=rec["B"]["B0"],
                          sigma=rec["sigma_shell"], C={k: rec["C"][k] for k in ("touch_fraction", "coverage", "width_median", "identified")}), ensure_ascii=False))


if __name__ == "__main__":
    main()
