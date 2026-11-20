"""Extends the d3 known-answer test (kat2): true r up to 3.0, fitter bounds r 0.4 – 8, profile 0.5 – 6.0, plus a group test."""
import datetime, hashlib, json, math, pathlib, sys
import numpy as np
from multiprocessing import get_context
from scipy.optimize import least_squares, minimize_scalar
from scipy.sparse import lil_matrix
from scipy.special import erf
REPO = pathlib.Path("<REPO>")
sys.path.insert(0, str(REPO / "decomp" / "stage1"))
import config as C  # noqa: E402
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import d3_kat as K1  # noqa: E402 — formula of the previous version · generator · reading the noise SD (used as is)
OUT = pathlib.Path("<WORK>/remeasure/explore/d3")
TES, BMM = K1.TES, K1.BMM
R_EXT = (1.0, 1.15, 1.3, 1.45, 1.6, 1.8, 2.0, 2.3, 2.64, 3.0)
LB, UB = np.array([0.01, 0.001, 20.0, 0.4]), np.array([0.99, 0.8, 300.0, 8.0])
STARTS = (0.8, 1.2, 1.6, 2.2, 3.0, 4.0)
RGRID = np.exp(np.linspace(np.log(0.5), np.log(6.0), 51))
NREP_C, NGROUP, SEED, TOL_A = 30, 30, 20261151, 1e-3
G_R, G_T2IN = (1.3, 1.6, 2.0, 2.64, 3.0), 82.0
for t2in in C.T2IN_GRID:                                   # Appended only to the running list — config.py is not edited
    for r in R_EXT:
        if (t2in, r, 1000.0) not in C.T2_SETS:
            C.T2_SETS.append((t2in, r, 1000.0))


def fit(data, sig, kw):
    best = None
    for r0 in STARTS:
        res = least_squares(K1.resid, np.array([0.6, 0.1, 80.0, r0]), bounds=(LB, UB), args=(data, sig, kw), xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=4000)
        if best is None or res.cost < best.cost:
            best = res
    return best


def profile(data, sig, best, kw):
    rss, prev = [], best.x[:3].copy()
    for r in RGRID:
        f = lambda t3: K1.resid(np.r_[t3, r], data, sig, kw)  # noqa: E731
        cand = [least_squares(f, x0, bounds=(LB[:3], UB[:3]), xtol=1e-10, ftol=1e-10, gtol=1e-10, max_nfev=2000) for x0 in (prev, best.x[:3])]
        res = min(cand, key=lambda z: z.cost)
        rss.append(2 * res.cost); prev = res.x
    return np.array(rss)


def ci_of(rss, best_cost2):
    m = min(best_cost2, rss.min())
    inside = rss - m <= 3.84
    return float(RGRID[inside].min()), float(RGRID[inside].max()), bool(inside[0] or inside[-1])


# ── Batch: n ROIs share (T2_in, r) · (ν_i, f_i) per ROI · the inner part all at once with a block-diagonal Jacobian ──
def model_vec(nu, f, t2in, r, dpar=1.7, t2iso=1000.0):
    a_in, a_en, a_iso = np.exp(-TES / t2in), np.exp(-TES / (t2in / r)), np.exp(-TES / t2iso)
    w_in = ((1 - f) * nu)[:, None] * a_in[None]; w_en = ((1 - f) * (1 - nu))[:, None] * a_en[None]; w_iso = f[:, None] * a_iso[None]
    b = BMM / 1000.0
    stick = np.sqrt(np.pi / (4 * b * dpar)) * erf(np.sqrt(b * dpar))
    dper = dpar * (1 - nu)
    dd = np.maximum(dpar - dper, 1e-12)[:, None]
    cyl = np.exp(-b[None] * dper[:, None]) * np.sqrt(np.pi / (4 * b[None] * dd)) * erf(np.sqrt(b[None] * dd))
    iso = np.exp(-b * C.DISO)
    S = (w_in[..., None] * stick[None, None] + w_en[..., None] * cyl[:, None, :] + w_iso[..., None] * iso[None, None]) / (w_in + w_en + w_iso)[..., None]
    return np.log(S)


def inner(D, sig, t2in, r, x0):
    n = D.shape[0]
    sp = lil_matrix((n * 16, 2 * n), dtype=int)
    for i in range(n):
        sp[i * 16:(i + 1) * 16, 2 * i:2 * i + 2] = 1
    fun = lambda z: ((D - model_vec(z[0::2], z[1::2], t2in, r)) / sig[None, None, :]).ravel()  # noqa: E731
    lb, ub = np.tile([0.01, 0.001], n), np.tile([0.99, 0.8], n)
    res = least_squares(fun, np.clip(x0, lb + 1e-9, ub - 1e-9), bounds=(lb, ub), jac_sparsity=sp, xtol=1e-10, ftol=1e-10, gtol=1e-10, max_nfev=3000)
    return 2 * res.cost, res.x


def pooled_profile(D, sig):
    n = D.shape[0]
    x0 = np.tile([0.6, 0.1], n)
    out, keep = [], x0
    for r in RGRID:
        cache = {}
        def g(lt):
            v, x = inner(D, sig, math.exp(lt), r, keep)
            cache[lt] = (v, x)
            return v
        m = minimize_scalar(g, bounds=(math.log(20.0), math.log(300.0)), method="bounded", options=dict(xatol=1e-3))
        v, x = cache.get(m.x, (None, None))
        if v is None:
            v, x = inner(D, sig, math.exp(m.x), r, keep)
        out.append(v); keep = x
    return np.array(out)


def worker_C(args):
    i, clean, r_true, sig = args
    rng = np.random.default_rng([SEED, 1, i])
    out = []
    for _ in range(NREP_C):
        data = clean + rng.normal(0.0, 1.0, clean.shape) * sig[None, :]
        best = fit(data, sig, {})
        lo, hi, touch = ci_of(profile(data, sig, best, {}), 2 * best.cost)
        out.append((float(best.x[3]), lo, hi, touch, lo <= r_true <= hi))
    return i, out


def worker_G(args):
    gi, D, r_true, sig = args
    rss = pooled_profile(D, sig)
    lo, hi, touch = ci_of(rss, rss.min())
    read = ("T2 비로는 B1의 얕은 각도가 설명되지 않음" if (lo <= 1.6 <= hi and not lo <= 2.6 <= hi) else
            "T2 비가 문헌보다 큼" if (lo <= 2.6 <= hi and not lo <= 1.6 <= hi) else "판정 보류")
    return gi, dict(r_true=r_true, rhat=float(RGRID[int(np.argmin(rss))]), ci=[lo, hi], touch=touch, covers=bool(lo <= r_true <= hi), reading=read)


def main():
    rec = dict(made=datetime.datetime.now().astimezone().isoformat(timespec="seconds"), script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
               d3_kat_sha256=hashlib.sha256((pathlib.Path(__file__).resolve().parent / "d3_kat.py").read_bytes()).hexdigest(),
               r_truth=R_EXT, rgrid=[float(RGRID[0]), float(RGRID[-1]), len(RGRID)], bounds_r=[float(LB[3]), float(UB[3])], tol_A=TOL_A)
    g = C.fibonacci(20000); cos2 = (g @ np.array([0.0, 0.0, 1.0])) ** 2
    v5 = []
    for b in BMM / 1000.0:
        st, cy, _ = K1.smean(np.array([b]), 1.7, 0.6)
        v5 += [abs(float(np.mean(np.exp(-b * 1.7 * cos2))) - float(st[0])), abs(float(np.mean(np.exp(-b * 0.6 - b * 1.1 * cos2))) - float(np.ravel(cy)[0]))]
    rec["V5_max_abs_diff"] = max(v5)
    if max(v5) > 1e-4:
        sys.exit("⛔ V5")
    T = [(ti, k) for ti in C.WM for k, (t, r, i) in enumerate(C.T2_SETS) if i == 1000.0 and r in R_EXT]
    tr = [(C.TISSUE[ti][0], C.TISSUE[ti][2], C.T2_SETS[k][0], C.T2_SETS[k][1]) for ti, k in T]
    sig1 = np.array([1.0, 1.0])
    est = [fit(K1.gen_A(ti, k), sig1, {}).x for ti, k in T]
    bad = [i for i, (e, t) in enumerate(zip(est, tr)) if not (abs(e[3] - t[3]) <= TOL_A and abs(e[0] - t[0]) <= TOL_A and abs(e[1] - t[1]) <= TOL_A)]
    rec["A1"] = dict(passed=not bad, n_truths=len(T), n_violating=len(bad), max_abs_dr=float(max(abs(e[3] - t[3]) for e, t in zip(est, tr))))
    if bad:
        (OUT / "d3_kat2.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1) + "\n", encoding="utf-8"); sys.exit(f"⛔ A-1 {rec['A1']}")
    perm = np.random.default_rng(SEED).permutation(len(TES)); assert not np.array_equal(perm, np.arange(len(TES)))
    rec["A2"] = {}
    for name, tf in (("TE 순서 섞음", lambda D: D[perm]), ("셸 바꿈", lambda D: D[:, ::-1])):
        e2 = [fit(tf(K1.gen_A(ti, k)), sig1, {}).x for ti, k in T]
        nb = sum(1 for e, t in zip(e2, tr) if not (abs(e[3] - t[3]) <= TOL_A and abs(e[0] - t[0]) <= TOL_A and abs(e[1] - t[1]) <= TOL_A))
        rec["A2"][name] = dict(n_violating=nb, criterion_passed=(nb == 0))
        if nb == 0:
            sys.exit(f"⛔ A-2 «{name}» 실패하지 않음")
    import synth  # noqa: E402
    pr = synth.Proto("P1")
    t1 = {nm: C.T1_SETS.index(v) for nm, v in (("1100|1.0", (1100.0, 1.0)), ("1100|1.3", (1100.0, 1.3)))}
    rows = dict(B0=("1100|1.0", 1.7, {}), B1_T1ratio13=("1100|1.3", 1.7, {}), B2_data_dpar20=("1100|1.0", 2.0, {}),
                B3_fit_dpar15=("1100|1.0", 1.7, dict(dpar=1.5)), B4_fit_dpar20=("1100|1.0", 1.7, dict(dpar=2.0)))
    cleanB, rec["B"] = {}, {}
    for name, (t1n, dp, kw) in rows.items():
        key = (t1n, dp)
        if key not in cleanB:
            cleanB[key] = [K1.gen_B(pr, ti, k, t1[t1n], dp, i) for i, (ti, k) in enumerate(T)]
        dr = np.array([fit(D, sig1, kw).x[3] - t[3] for D, t in zip(cleanB[key], tr)])
        byr = {str(r): float(np.median(dr[[j for j, t in enumerate(tr) if t[3] == r]])) for r in R_EXT}
        rec["B"][name] = dict(median=float(np.median(dr)), q25=float(np.quantile(dr, 0.25)), q75=float(np.quantile(dr, 0.75)), max_abs=float(np.max(np.abs(dr))), median_by_r=byr)
    sig, nd = K1.reread_sigma()
    rec["sigma_shell"] = dict(b700=float(sig[0]), b2000=float(sig[1]), n_diffs=nd)
    jobs = [(i, cleanB[("1100|1.0", 1.7)][i], tr[i][3], sig) for i in range(len(T))]
    with get_context("fork").Pool(8) as pool:
        resC = dict(pool.map(worker_C, jobs))
    byr = {}
    for i, t in enumerate(tr):
        byr.setdefault(t[3], []).extend(resC[i])
    rec["C"] = {str(r): dict(touch=float(np.mean([x[3] for x in v])), coverage=float(np.mean([x[4] for x in v])),
                             rhat_median=float(np.median([x[0] for x in v])), width_median=float(np.median([x[2] - x[1] for x in v]))) for r, v in sorted(byr.items())}
    # Batch known-answer test — generator B · ROI 27 · random tissue · random fiber direction · true value r ∈ G_R · T2_in 82
    rng = np.random.default_rng([SEED, 2])
    gjobs, gi = [], 0
    for r in G_R:
        k = C.T2_SETS.index((G_T2IN, r, 1000.0))
        for _ in range(NGROUP):
            tis = rng.choice(C.WM, 27); fds = rng.integers(0, len(C.FIBRE_DIRS), 27)
            D = np.array([K1.gen_B(pr, int(ti), k, t1["1100|1.0"], 1.7, int(fd)) for ti, fd in zip(tis, fds)])
            D = D + rng.normal(0.0, 1.0, D.shape) * sig[None, None, :]
            gjobs.append((gi, D, r, sig)); gi += 1
    with get_context("fork").Pool(8) as pool:
        resG = dict(pool.map(worker_G, gjobs))
    rec["G"] = {}
    for r in G_R:
        v = [x for x in resG.values() if x["r_true"] == r]
        reads = {}
        for x in v:
            reads[x["reading"]] = reads.get(x["reading"], 0) + 1
        rec["G"][str(r)] = dict(n=len(v), coverage=float(np.mean([x["covers"] for x in v])), touch=float(np.mean([x["touch"] for x in v])),
                                rhat_median=float(np.median([x["rhat"] for x in v])), width_median=float(np.median([x["ci"][1] - x["ci"][0] for x in v])), readings=reads)
    gate = rec["A1"]["passed"] and all(not v["criterion_passed"] for v in rec["A2"].values()) and rec["G"]["1.6"]["coverage"] >= 0.85 and rec["G"]["2.64"]["coverage"] >= 0.85
    rec["gate_for_real_fit"] = dict(passed=bool(gate), rule="A-1 통과 · A-2 둘 다 실패 · 묶음 덮음 ≥ 0.85(참값 1.6 · 2.64) — 계획서 §13(2026-11-02 00:09:25)")
    (OUT / "d3_kat2.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(dict(V5=rec["V5_max_abs_diff"], A1=rec["A1"], A2={k: v["n_violating"] for k, v in rec["A2"].items()}, B0_by_r=rec["B"]["B0"]["median_by_r"],
                          C={r: (round(v["touch"], 3), round(v["coverage"], 3), round(v["width_median"], 3)) for r, v in rec["C"].items()},
                          G={r: (v["coverage"], v["width_median"], v["readings"]) for r, v in rec["G"].items()}, gate=rec["gate_for_real_fit"]), ensure_ascii=False))


if __name__ == "__main__":
    main()
