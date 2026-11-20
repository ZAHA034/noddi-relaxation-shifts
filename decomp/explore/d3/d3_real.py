"""Fits d3 to the Wong measurement only if the d3_kat2 gate is true: participant and whole-group fits sharing T2_in and r, with sensitivity variants."""
import datetime, hashlib, json, math, pathlib, sys
import numpy as np
from multiprocessing import get_context
from scipy.optimize import least_squares, minimize_scalar
from scipy.sparse import lil_matrix
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import d3_kat as K1  # noqa: E402
import d3_kat2 as K2  # noqa: E402
OUT = pathlib.Path("<WORK>/remeasure/explore/d3")
WONG = pathlib.Path("<WORK>/remeasure/wong")
LOG = pathlib.Path("<WORK>/stage1_work/drafts/x14/REREAD_LOG.jsonl")
TES = ("TE62", "TE72", "TE82", "TE92", "TE102", "TE112", "TE122", "TE132")
SENS = {"main": dict(dpar=1.7, t2iso=1000.0), "dpar1.5": dict(dpar=1.5, t2iso=1000.0), "dpar2.0": dict(dpar=2.0, t2iso=1000.0), "t2iso2000": dict(dpar=1.7, t2iso=2000.0)}


def inner(D, sig, t2in, r, x0, kw):
    n = D.shape[0]
    sp = lil_matrix((n * 16, 2 * n), dtype=int)
    for i in range(n):
        sp[i * 16:(i + 1) * 16, 2 * i:2 * i + 2] = 1
    fun = lambda z: ((D - K2.model_vec(z[0::2], z[1::2], t2in, r, **kw)) / sig[None, None, :]).ravel()  # noqa: E731
    lb, ub = np.tile([0.01, 0.001], n), np.tile([0.99, 0.8], n)
    res = least_squares(fun, np.clip(x0, lb + 1e-9, ub - 1e-9), bounds=(lb, ub), jac_sparsity=sp, xtol=1e-10, ftol=1e-10, gtol=1e-10, max_nfev=3000)
    return 2 * res.cost, res.x


def profile(args):
    key, sub, D, sig, kw = args
    keep, out = np.tile([0.6, 0.1], D.shape[0]), []
    for r in K2.RGRID:
        cache = {}
        def g(lt):
            v, x = inner(D, sig, math.exp(lt), r, keep, kw); cache[lt] = (v, x); return v
        m = minimize_scalar(g, bounds=(math.log(20.0), math.log(300.0)), method="bounded", options=dict(xatol=1e-3))
        v, x = cache.get(m.x) or inner(D, sig, math.exp(m.x), r, keep, kw)
        out.append(v); keep = x
    return key, sub, np.array(out)


def ci_read(rss):
    inside = rss - rss.min() <= 3.84
    lo, hi = float(K2.RGRID[inside].min()), float(K2.RGRID[inside].max())
    read = ("T2 비로는 B1의 얕은 각도가 설명되지 않음" if (lo <= 1.6 <= hi and not lo <= 2.6 <= hi) else
            "T2 비가 문헌보다 큼" if (lo <= 2.6 <= hi and not lo <= 1.6 <= hi) else "판정 보류")
    return dict(ci=[lo, hi], touch=bool(inside[0] or inside[-1]), r_grid_argmin=float(K2.RGRID[int(np.argmin(rss))]), reading=read)


def main():
    gate = json.loads((OUT / "d3_kat2.json").read_text(encoding="utf-8")).get("gate_for_real_fit", {})
    if not gate.get("passed"):
        sys.exit(f"⛔ 넘어가는 조건이 참이 아니다 — 실측 없이 멈춤: {gate}")
    now = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    data = {}
    for sub in ("sub-01", "sub-02", "sub-03"):
        own = set(json.loads((WONG / sub / "DONE.json").read_text(encoding="utf-8")).get("따로_정합", []))
        tabs = []
        for s in TES:
            p = (WONG / sub / f"own_{s}" if s in own else WONG / sub) / f"roi_{s}.json"
            b = p.read_bytes()
            with LOG.open("a", encoding="utf-8") as f:
                f.write(json.dumps(dict(시각=now, 다시읽음=str(p), sha256=hashlib.sha256(b).hexdigest(), 바이트=len(b), 무엇="D3 실측 적합(묶음)", 스크립트=__file__), ensure_ascii=False) + "\n")
            tabs.append(json.loads(b)["주"])
        units = sorted(u for u in tabs[0] if all(u in t and all(t[u].get(k) is not None and math.isfinite(t[u][k]) for k in ("s_low", "s_2000")) for t in tabs))
        data[sub] = (units, np.array([[[t[u]["s_low"], t[u]["s_2000"]] for t in tabs] for u in units], float))
    sig, nd = K1.reread_sigma()
    jobs = [(key, sub, data[sub][1], sig, kw) for key, kw in SENS.items() for sub in data]
    with get_context("fork").Pool(8) as pool:
        prof = pool.map(profile, jobs)
    res = dict(made=now, script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(), gate=gate, sigma_shell=dict(b700=float(sig[0]), b2000=float(sig[1]), n_diffs=nd),
               n_units={s: len(data[s][0]) for s in data}, rgrid=[float(K2.RGRID[0]), float(K2.RGRID[-1]), len(K2.RGRID)], by_sens={})
    for key in SENS:
        P = {sub: rss for k, sub, rss in prof if k == key}
        res["by_sens"][key] = dict({sub: ci_read(P[sub]) for sub in P}, overall=ci_read(sum(P.values())))
    (OUT / "d3_real.json").write_text(json.dumps(res, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    for key, v in res["by_sens"].items():
        print(key, {s: (round(x["ci"][0], 2), round(x["ci"][1], 2), x["reading"]) for s, x in v.items()})


if __name__ == "__main__":
    main()
