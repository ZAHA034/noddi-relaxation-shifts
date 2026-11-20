"""cause1b: exploratory supplement repeating the raw-value comparison on all OXF1PRI~OXF2PRI ROI units regardless of the judge, plus check 2(c); not a verdict."""
from __future__ import annotations

import hashlib
import json
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import cause1_raw as C1  # noqa: E402

J = C1.J


def main():
    rng = np.random.default_rng(C1.SEED)
    prev = json.loads((C1.OUT / "cause1_raw.json").read_text(encoding="utf-8"))
    d = dict(np.load(C1.D1 / "d1_c3_units" / f"{C1.PAIR}.npz", allow_pickle=False))
    use191 = ~d["main_excluded"] & ~d["flag"] & d["ok"].all(1)
    D_all = d["X3"][:, 3:6] - d["X3"][:, 0:3]
    subs = d["subjects"]
    sets = {"all": np.ones(len(subs), bool), "all_minus_main_excluded": ~d["main_excluded"], "d1_191": use191}
    # Simulator expectation — same formula as ① · same evaluation points (the point estimate is the same as in ① · the bootstrap uses this script's random numbers)
    sim = np.load(J.S.WORK / "sim_test.npz", allow_pickle=False)
    N = J.S.COND_NAMES
    lay = np.array([l == "match|L0" for l in sim["layer"].tolist()])
    C, TH, XC = sim["cond"].astype(float), sim["theta"], sim["x_c3"].astype(float)
    fin = np.isfinite(XC).all(1)
    sA_s, sB_s = np.exp(C[:, N.index("log_snr_A")]), np.exp(C[:, N.index("log_snr_B")])
    uA, uB = prev["eval_point"]["uA"], prev["eval_point"]["uB"]
    p2 = np.r_[1, 0, 0, 0, uA, uB]
    exp_, boots, filt = {}, {}, {}
    for cn, (ca, cb) in C1.COMBOS.items():
        base = lay & (C[:, N.index("coil_A_64")] == ca) & (C[:, N.index("coil_B_64")] == cb)
        filt[cn] = dict(n_cases=int(base.sum()), n_nonfinite_dropped=int((base & ~fin).sum()),
                        n_clipped_A=int((base & sim["clipped_A"]).sum()), n_clipped_B=int((base & sim["clipped_B"]).sum()),
                        n_clipped_A_or_B=int((base & (sim["clipped_A"] | sim["clipped_B"])).sum()))
        sel = base & fin
        Ds = XC[sel, 3:6] - XC[sel, 0:3]
        X2 = C1.design(TH[sel], 1 / sA_s[sel] ** 2, 1 / sB_s[sel] ** 2, "inv2")
        e, bb = {}, {c: [] for c in C1.COMP}
        for c, j in C1.COMP.items():
            e[c] = float(C1.fit_eval(X2, Ds[:, j], p2)[0])
            assert abs(e[c] - prev["expected"][cn][c]["primary"]) < 1e-9, "① 의 기대와 다르다"
        for _ in range(C1.BOOT):
            ix = rng.integers(0, len(Ds), len(Ds))
            for c, j in C1.COMP.items():
                bb[c].append(float(C1.fit_eval(X2[ix], Ds[ix, j], p2)[0]))
        exp_[cn], boots[cn] = e, bb
    res = dict(script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(), pair=C1.PAIR, n_rows=int(len(subs)),
               n_main_excluded=int(d["main_excluded"].sum()), filter_2c=filt, expected=exp_, sets={})
    for sn, m in sets.items():
        us = np.unique(subs[m])
        out = dict(n_units=int(m.sum()), n_subjects=int(len(us)))
        for c, j in (("x", 0), ("y", 1), ("ODI", 2)):
            v = D_all[m, j]
            bm = []
            for _ in range(C1.BOOT):
                pick = rng.choice(us, len(us), replace=True)
                bm.append(float(np.mean(np.concatenate([v[subs[m] == s] for s in pick]))))
            e64, e32 = exp_["32|64"][c], exp_["32|32"][c]
            br = (np.array(bm) - np.array(boots["32|32"][c])) / (np.array(boots["32|64"][c]) - np.array(boots["32|32"][c]))
            out[c] = dict(mean=float(v.mean()), median=float(np.median(v)), sd=float(v.std(ddof=1)),
                          boot95=[float(np.quantile(bm, 0.025)), float(np.quantile(bm, 0.975))],
                          r=float((v.mean() - e32) / (e64 - e32)), r_boot95=[float(np.quantile(br, 0.025)), float(np.quantile(br, 0.975))])
        res["sets"][sn] = out
    (C1.OUT / "cause1b_allunits.json").write_text(json.dumps(res, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(dict(n_rows=res["n_rows"], n_main_excluded=res["n_main_excluded"], filter_2c=filt), ensure_ascii=False))
    for sn, out in res["sets"].items():
        print(sn, out["n_units"], out["n_subjects"], " | ".join(
            f"{c} {out[c]['mean']:+.4f} [{out[c]['boot95'][0]:+.4f}, {out[c]['boot95'][1]:+.4f}] SD {out[c]['sd']:.3f} r {out[c]['r']:+.2f} [{out[c]['r_boot95'][0]:+.2f}, {out[c]['r_boot95'][1]:+.2f}]"
            for c in ("x", "ODI", "y")))


if __name__ == "__main__":
    main()
