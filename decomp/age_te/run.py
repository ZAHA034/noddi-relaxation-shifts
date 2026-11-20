"""age_te run — runs only after the plan decomp/age_te/ANALYSIS_PLAN.md is fixed by a commit (stops if the working copy of the plan differs from the version committed at HEAD); E3 without --participants, E1 only with it."""
import hashlib
import json
import math
import pathlib
import subprocess
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "remeasure"))
import judge_tong as JT  # noqa: E402
import judge_v4 as V4  # noqa: E402
import judge_wong as JW  # noqa: E402
import rconfig as R  # noqa: E402

from age_te import model as M  # noqa: E402
from age_te.participants import assert_no_ids, parse_participants  # noqa: E402

N_PERM = 10000


def acq_of(sub, ses):
    d = R.device_of(ses)
    return R.not2ing_protocol(sub) if d == "NOT2ING" else d


def onh_rows(work):
    """The same first sessions as at the verdict (main-analysis exclusion 1 · retry failures · collisions · no end mark left out) → per session, per-unit (x, y, ODI) (judge_tong.unit_coords)."""
    sess = [s for s in R.onh_sessions() if s not in R.MAIN_EXCLUDE_ONH]
    states = V4.session_states(work, sess)
    rows = []
    for sub, ses in sess:
        key = f"{sub}_{ses}"
        if states.get(key) != V4.USABLE:
            continue
        per = JT.unit_coords(json.loads((work / "onharmony" / key / "roi.json").read_text(encoding="utf-8")))
        rows.append(dict(sub=sub, acq=acq_of(sub, ses), units={u: tuple(map(float, v)) for u, v in per.items()}))
    return rows


def unit_means(rows):
    """Main analysis: mean over the units that have values in all sessions (common units) — prevents differences arising from different unit sets between sessions (plan Section 4)."""
    common = sorted(set.intersection(*[set(r["units"]) for r in rows]))
    for r in rows:
        a = np.array([r["units"][u] for u in common], float)
        r["x"], r["y"], r["odi"] = (float(v) for v in a.mean(0))
    return common


def covariate_block(rows, z, label, n_perm=N_PERM):
    """E1 — main (x) · negative control (ODI) · descriptive (y) · falsification controls (TE permutation · departure of the acquisition-wise slopes from a line in TE) · three sensitivity analyses (vendor-wise slope terms · within Philips · unit level)."""
    rows = [r for r in rows if r["sub"] in z]
    persons, acqs = [r["sub"] for r in rows], [r["acq"] for r in rows]
    zz = {p: z[p] for p in set(persons)}
    out = dict(label=label, n_sessions=len(rows), n_persons=len(zz))
    for feat in ("x", "odi", "y"):
        v = [r[feat] for r in rows]
        f = M.fit(v, persons, acqs, zz)
        f["mde80"] = M.mde(f["se"], f["df"])
        if feat in ("x", "odi"):
            f["perm_p"] = M.permutation_p(v, persons, acqs, zz, n_perm=n_perm)
            f["te_permutation"] = M.te_permutation_p(v, persons, acqs, zz)
            f["acq_slopes_lack_of_fit"] = M.acq_slopes_lack_of_fit(v, persons, acqs, zz)
        out[feat] = f
    zbar = float(np.mean(list(zz.values())))
    extra = [np.array([(M.VENDOR[a] == vv) * (zz[p] - zbar) for p, a in zip(persons, acqs)], float) for vv in ("Philips", "Siemens")]
    out["x_vendor_terms"] = M.fit([r["x"] for r in rows], persons, acqs, zz, extra)
    ph = [r for r in rows if M.VENDOR[r["acq"]] == "Philips"]
    out["x_within_philips"] = M.fit([r["x"] for r in ph], [r["sub"] for r in ph], [r["acq"] for r in ph], {r["sub"]: zz[r["sub"]] for r in ph})
    out["x_per_unit"] = M.fit_units([r["units"] for r in rows], persons, acqs, zz, coord=0)
    return out


def e3_multite(work):
    """E3 — each person's TE slope in the multi-TE data (common-unit means of x · y · ODI on the 8 TEs, repeat series left out of the regression) · homogeneity between people · comparison of the between-person SD (τ) with the assumed value of the design simulation. Two SEs — residual of the line fit (main) · measurement noise measured from the differences of repeat scans (sensitivity)."""
    data, own, missing = JW.load(work)
    out, het = {}, {}
    for k, (sb, series) in enumerate(sorted(data.items()), 1):
        tes = [t for t in JW.TES if f"TE{t}" in series]
        common = sorted(set.intersection(*[set(series[s]) for s in series]))
        d = dict(n_te=len(tes), n_units=len(common), separately_registered=len(own.get(sb, ())))
        t = np.asarray(tes, float)
        sxx = float(np.sum((t - t.mean()) ** 2))
        for ci, co in enumerate(("x", "y", "odi")):
            mean_at = {s: float(np.mean([series[s][u][ci] for u in common])) for s in series}
            f = M.te_slope(tes, [mean_at[f"TE{t}"] for t in tes])
            diffs = [mean_at[r] - mean_at[o] for r, o in JW.REPEATS.items() if r in mean_at and o in mean_at]
            sm = math.sqrt(sum(x * x for x in diffs) / (2 * len(diffs))) if diffs else None
            f["se_retest"] = sm / math.sqrt(sxx) if sm is not None else None
            f["n_repeat_pairs"] = len(diffs)
            d[co] = f
        out[f"participant_{k}"] = d
    for co in ("x", "y", "odi"):
        sl = [d[co]["slope"] for d in out.values()]
        se = [d[co]["se"] for d in out.values()]
        het[co] = M.slope_heterogeneity(sl, se)
        het[co]["between_sd"] = M.between_sd(sl, se)
        if all(d[co]["se_retest"] for d in out.values()):
            ser = [d[co]["se_retest"] for d in out.values()]
            het[co]["retest_se"] = M.slope_heterogeneity(sl, ser)
            het[co]["retest_se"]["between_sd"] = M.between_sd(sl, ser)
    return dict(by_participant=out, heterogeneity=het, n_missing=len(missing),
                assumed_between_sd_x=M.assumed_slope_sd(), assumed_shift_sd_at_30ms=30.0 * M.assumed_slope_sd(),
                bounding_shift_sd_at_30ms=dict(independent=30.0 * M.assumed_slope_sd(rho=0.0), common=30.0 * M.assumed_slope_sd(rho=1.0)))


def plan_fixed(repo):
    head = subprocess.run(["git", "show", "HEAD:decomp/age_te/ANALYSIS_PLAN.md"], cwd=repo, capture_output=True)
    return head.returncode == 0 and head.stdout == (HERE / "ANALYSIS_PLAN.md").read_bytes()


def main(argv=None):
    a = list(sys.argv[1:] if argv is None else argv)
    ptsv = None
    if "--participants" in a:
        k = a.index("--participants")
        ptsv = pathlib.Path(a[k + 1])
        del a[k:k + 2]
    work, out = pathlib.Path(a[0]), pathlib.Path(a[1])
    if not plan_fixed(HERE.parents[1]):
        raise SystemExit("계획이 HEAD 에 고정되어 있지 않다 — 계산하지 않는다")
    res = dict(plan_sha256=hashlib.sha256((HERE / "ANALYSIS_PLAN.md").read_bytes()).hexdigest())
    ids = []
    if ptsv is not None:
        rows = onh_rows(work)
        common = unit_means(rows)
        age, _female, fmt = parse_participants(ptsv.read_text(encoding="utf-8"))   # sex is not used (E2 dropped · change record 2)
        subs = sorted({r["sub"] for r in rows})
        ids = subs
        agez = {s: age[s] for s in subs if s in age}
        res.update(n_sessions=len(rows), n_persons=len(subs), n_common_units=len(common),
                   sessions_per_acquisition={q: sum(r["acq"] == q for r in rows) for q in M.TE_ACQ},
                   participants_format={k: v for k, v in fmt.items() if not k.startswith("sex")},
                   age_decades=dict(n=len(agez), mean=round(float(np.mean(list(agez.values()))), 3) if agez else None,
                                    sd=round(float(np.std(list(agez.values()), ddof=1)), 3) if len(agez) > 1 else None))
        if len(agez) >= 3:
            res["E1_age"] = covariate_block(rows, agez, "age, per decade")
            g = res["E1_age"]["x"]
            res["E1_age"]["x_t2en_reading"] = dict(point=M.t2en_change(g["gamma"]), ci_low=M.t2en_change(g["ci95"][0]),
                                                   ci_high=M.t2en_change(g["ci95"][1]), mde80_lengthening=M.t2en_change(-g["mde80"]),
                                                   mde80_shortening=M.t2en_change(g["mde80"]))
    res["E3_multite"] = e3_multite(work)
    text = json.dumps(res, ensure_ascii=False, indent=1)
    assert_no_ids(text, ids=ids)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text + "\n", encoding="utf-8")
    print(json.dumps({k: res[k] for k in ("n_sessions", "n_persons", "n_common_units") if k in res}, ensure_ascii=False))


if __name__ == "__main__":
    main()
