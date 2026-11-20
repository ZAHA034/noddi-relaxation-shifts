"""Post hoc check2 of d3 v3: asks whether the true value falling far outside the interval at × 2 comes from the inner optimization."""
import hashlib, json, math, pathlib, sys
import numpy as np
from multiprocessing import get_context
from scipy.optimize import minimize_scalar
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / "B1S"))
import d3_v3 as V  # noqa: E402
import d3_real as DR  # noqa: E402
OUT = V.OUT
MAIN = "dpar1.7|t2iso1000"
LT = np.linspace(math.log(20.0), math.log(300.0), 41)
STARTS = ((0.6, 0.1), (0.4, 0.05), (0.8, 0.3))


def exh(D, sig, r, kw):
    grid = np.array([[DR.inner(D, sig, math.exp(lt), r, np.tile(s0, D.shape[0]), kw)[0] for s0 in STARTS] for lt in LT])
    i, j = np.unravel_index(int(np.argmin(grid)), grid.shape)
    x0 = np.tile(STARTS[j], D.shape[0])
    m = minimize_scalar(lambda lt: DR.inner(D, sig, math.exp(lt), r, x0, kw)[0], bounds=(LT[max(i - 1, 0)], LT[min(i + 1, len(LT) - 1)]), method="bounded", options=dict(xatol=1e-3))
    return float(min(grid.min(), m.fun)), float(math.exp(LT[i]))


def job(a):
    gi, r_true, scale, rhat, rss_min, sig = a
    D = V.make_D(gi, r_true, scale)
    kw = V.VARIANTS[MAIN]
    at_true, t2_true = exh(D, sig, r_true, kw)
    at_hat, t2_hat = exh(D, sig, rhat, kw)
    return dict(gi=gi, r_true=r_true, scale=scale, rhat=rhat, rss_min_reported=rss_min, exh_at_true=at_true, exh_at_rhat=at_hat,
                d_exh=at_true - at_hat, exh_rhat_minus_reported=at_hat - rss_min, t2in_grid_true=t2_true, t2in_grid_rhat=t2_hat)


def main():
    V.build_lib()
    s = json.loads((OUT / "d3_v3.json").read_text(encoding="utf-8"))["sigma_shell"]  # Same σ as in the main computation · no re-access
    sig = np.array([s["b700"], s["b2000"]])
    rows = [json.loads(x) for x in (OUT / "d3_v3_groups.jsonl").read_text(encoding="utf-8").splitlines() if x.strip()]
    jobs = [(x["gi"], x["r_true"], x["scale"], x["variants"][MAIN]["rhat"], x["variants"][MAIN]["rss_min"], sig) for x in rows if x["scale"] == 2.0 and "variants" in x]
    with get_context("fork").Pool(16) as pool:
        res = sorted(pool.map(job, jobs), key=lambda d: (d["r_true"], d["gi"]))
    bad = [d["gi"] for d in res if not d["d_exh"] > 3.84]
    rec = dict(rule="ΔRSS_exh = min(r = 참값) − min(r = r̂) > 3.84 이 모든 묶음에서 맞으면 «참값이 구간 밖» 은 안쪽 최적화 탓이 아니다(결과를 본 뒤 · 구현 점검)",
               script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(), sigma=sig.tolist(), n=len(res), failed_gi=bad, passed=not bad,
               d_exh_min={str(rt): float(min(d["d_exh"] for d in res if d["r_true"] == rt)) for rt in V.TRUTHS},
               d_exh_median={str(rt): float(np.median([d["d_exh"] for d in res if d["r_true"] == rt])) for rt in V.TRUTHS},
               exh_rhat_minus_reported_min=float(min(d["exh_rhat_minus_reported"] for d in res)), rows=res)
    (OUT / "d3_v3_check2.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")
    print(json.dumps({k: v for k, v in rec.items() if k != "rows"}, ensure_ascii=False))


if __name__ == "__main__":
    main()
