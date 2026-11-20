"""Computes the b1 Wong TE trajectory direction per ROI unit (post hoc, descriptive) from roi_TE*.json, stopping if it disagrees with E1 ψ."""
import datetime, hashlib, json, math, pathlib, sys
import numpy as np
REPO = pathlib.Path("<REPO>")
sys.path.insert(0, str(REPO / "decomp" / "remeasure"))
import judge_wong as JW  # noqa: E402
import rconfig as R  # noqa: E402
OUT = pathlib.Path("<WORK>/remeasure/explore/b1")
E1 = pathlib.Path("<WORK>/stage1_work/drafts/manuscript_tex_2026-10-31/exploratory/e1_results.json")
D42 = pathlib.Path("<WORK>/stage1_work/drafts/x14/D4/d4_2_forward.json")
LOG = pathlib.Path("<WORK>/stage1_work/drafts/x14/REREAD_LOG.jsonl")
W_LO, W_HI = 48.5, 72.1
TH_C = (W_LO + W_HI) / 2
SEED, NDRAW = 20261141, 10000
now = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
assert hashlib.sha256(D42.read_bytes()).hexdigest().startswith("dfc2aa513b1fc7bc"), "D4.2 sha 가 다르다"
data, own, missing = JW.load()
assert not missing, missing
with LOG.open("a", encoding="utf-8") as f:
    for sb in R.WONG_SUBJECTS:
        for s in R.WONG_SERIES:
            p = (R.WORK / "wong" / sb / (f"own_{s}" if s in own[sb] else "")) / f"roi_{s}.json"
            b = p.read_bytes()
            f.write(json.dumps(dict(시각=now, 다시읽음=str(p), sha256=hashlib.sha256(b).hexdigest(), 바이트=len(b), 무엇="B1 Wong TE 궤적 방향(judge_wong.load)", 스크립트=__file__), ensure_ascii=False) + "\n")
e1 = json.loads(E1.read_text(encoding="utf-8"))["units"]
te = np.array(JW.TES, float); t = te - JW.TE_C; stt = float((t ** 2).sum())
rng = np.random.default_rng(SEED)
c, s_ = math.cos(math.radians(TH_C)), math.sin(math.radians(TH_C))
res, maxdiff = {}, 0.0
for sb in R.WONG_SUBJECTS:
    ser = data[sb]
    units = sorted(set.intersection(*[set(ser[f"TE{x}"]) for x in JW.TES]))
    rep = []
    for r2, orig in JW.REPEATS.items():
        for u in units:
            if u in ser.get(r2, {}) and u in ser[orig]:
                a, b = np.asarray(ser[orig][u], float), np.asarray(ser[r2][u], float)
                rep.append((-s_ * b[0] + c * b[1]) - (-s_ * a[0] + c * a[1]))
    sd_rep = float(np.sqrt(np.mean(np.square(rep))))
    rows = []
    for u in units:
        X = np.array([ser[f"TE{x}"][u][0] for x in JW.TES], float); Y = np.array([ser[f"TE{x}"][u][1] for x in JW.TES], float)
        bx, by = JW.fit_unit(te, X)["slope"], JW.fit_unit(te, Y)["slope"]
        th = math.degrees(math.atan2(by, bx))
        k1 = [k for k in e1[sb][u] if k.lower().startswith("psi")]
        assert len(k1) == 1, (sb, u, list(e1[sb][u]))
        maxdiff = max(maxdiff, abs(th - float(e1[sb][u][k1[0]])))
        if abs(th - float(e1[sb][u][k1[0]])) > 1e-9:
            sys.exit(f"⛔ θ_obs ≠ E1 ψ: {sb} {u} {th} {e1[sb][u][k1[0]]}")
        rx = X - np.polyval(np.polyfit(t, X, 1), t); ry = Y - np.polyval(np.polyfit(t, Y, 1), t)
        Rr = np.cov(np.vstack([rx, ry]), ddof=2)
        dr = rng.multivariate_normal([bx, by], Rr / stt, NDRAW)
        ths = np.degrees(np.arctan2(dr[:, 1], dr[:, 0]))
        ths = th + ((ths - th + 180.0) % 360.0 - 180.0)
        lo, hi = (float(v) for v in np.percentile(ths, [2.5, 97.5]))
        bperp = -s_ * bx + c * by
        rows.append(dict(unit=u, bx=bx, by=by, theta=th, ci=[lo, hi], in_wedge=W_LO <= th <= W_HI, ci_overlaps_wedge=(lo <= W_HI and hi >= W_LO),
                         bperp_per10ms=10 * bperp, ratio_perp_drift_to_rep=abs(bperp) * 70.0 / sd_rep))
    th_all = np.array([r["theta"] for r in rows])
    res[sb] = dict(n=len(rows), theta_median=float(np.median(th_all)), theta_min=float(th_all.min()), theta_max=float(th_all.max()),
                   in_wedge=int(sum(r["in_wedge"] for r in rows)), ci_overlaps_wedge=int(sum(r["ci_overlaps_wedge"] for r in rows)),
                   bperp_per10ms_median=float(np.median([r["bperp_per10ms"] for r in rows])), sd_rep_p=sd_rep, n_rep_diffs=len(rep),
                   ratio_median=float(np.median([r["ratio_perp_drift_to_rep"] for r in rows])),
                   ratio_min=float(min(r["ratio_perp_drift_to_rep"] for r in rows)), ratio_max=float(max(r["ratio_perp_drift_to_rep"] for r in rows)), units=rows)
allr = [r for sb in res for r in res[sb]["units"]]
overall = dict(n=len(allr), theta_median=float(np.median([r["theta"] for r in allr])), in_wedge=int(sum(r["in_wedge"] for r in allr)),
               ci_overlaps_wedge=int(sum(r["ci_overlaps_wedge"] for r in allr)), ratio_median=float(np.median([r["ratio_perp_drift_to_rep"] for r in allr])),
               max_abs_diff_vs_E1_deg=maxdiff)
out = dict(label="B1 — post hoc after E1 (exploratory · descriptive · no threshold)", made=now, script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
           plan_sha256_at_definition="589bafb09a64e233…", wedge=[W_LO, W_HI], theta_c=TH_C, seed=SEED, ndraw=NDRAW, by_participant=res, overall=overall)
(OUT / "b1_wong_direction.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
L = ["| 참가자 | 단위 | θ_obs 중앙 [범위] | 쐐기 안 | 구간이 쐐기와 겹침 | b⊥ 중앙(10 ms 당) | SD_rep(p) | 비 |b⊥|·70 ms ÷ SD_rep 중앙 [범위] |", "|---|---|---|---|---|---|---|---|"]
for sb, r in res.items():
    L.append(f"| {sb} | {r['n']} | {r['theta_median']:.1f}° [{r['theta_min']:.1f}, {r['theta_max']:.1f}] | {r['in_wedge']}/{r['n']} | {r['ci_overlaps_wedge']}/{r['n']} | {r['bperp_per10ms_median']:+.4f} | {r['sd_rep_p']:.4f} ({r['n_rep_diffs']}) | {r['ratio_median']:.2f} [{r['ratio_min']:.2f}, {r['ratio_max']:.2f}] |")
L.append(f"| 전체 | {overall['n']} | {overall['theta_median']:.1f}° | {overall['in_wedge']}/{overall['n']} | {overall['ci_overlaps_wedge']}/{overall['n']} | — | — | {overall['ratio_median']:.2f} |")
(OUT / "b1_wong_direction.md").write_text("\n".join(L) + "\n", encoding="utf-8")
print("\n".join(L)); print("E1 대조 최대 차", maxdiff)
