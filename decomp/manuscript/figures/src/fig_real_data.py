"""Fig. 3 (S4 cards, C3), Fig. 5 (Tong) and Fig. S2 (Wongkornchaovalit) from stored values — no template marks.
The drawing code is fig_data_templates_v3.draw_fig3 / draw_fig5 / draw_figS1; only the data dict is real.
No new statistic is computed. Each figure stops unless its values agree with the stored judge output:
  Fig. 3  identified counts = NDI 1 (the matched pair) · FISO 0 · ODI 28 (cards rule, REPORT.md:161).
  Fig. 5  tong_judge.json main only: T-A counts, margin = 1.96 · √2 · max(s_xy_unit_median), limit 23,
          D share 0.181, undeterminable branch; T-B share = n_retest / n_den and verdict at 0.90; T-C.
  Fig. S2 per-unit values from the stored ROI files (judge_wong.load); per-unit fits, intercepts and
          curvature by the judge's own functions; per participant s_rt, n_diffs, n_units, linear share
          and the intercept-range quartiles, and the overall linear share = wong_judge.json.
Run: MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1 <pdenv>/python -I -B fig_real_data.py [fig3] [fig5] [figS2]
"""
import hashlib
import json
import math
import pathlib
import sys

sys.dont_write_bytecode = True
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402
import fig_data_templates_v3 as T  # noqa: E402

REPO = pathlib.Path("<repo>")
CARDS = REPO / "decomp" / "results" / "s4" / "cards_c3.json"
TONG_JUDGE = pathlib.Path("data/remeasure/tong/judge/tong_judge.json")
WONG_JUDGE = pathlib.Path("data/remeasure/wong/judge/wong_judge.json")
OUT = HERE.parent


def same(a, b, what, tol=1e-12):
    if a is None or b is None or abs(a - b) > tol:
        sys.exit(f"STOP {what}: {a} != {b}")


def report(name):
    for ext in ("pdf", "png"):
        f = OUT / f"{name}.{ext}"
        print(f"sha256 {f.name} {hashlib.sha256(f.read_bytes()).hexdigest()}")


def fig3():
    cells = json.loads(CARDS.read_text(encoding="utf-8"))["cells"]

    def cell(p, snr, neff):
        for a, b in (p, p[::-1]):
            k = f"{a}~{b}|L0|SNR{snr}|Neff{neff}"
            if k in cells:
                return cells[k]
        raise KeyError(p)

    C, pdep, mref, mcorner = ({c: {} for c in T.COMPONENTS} for _ in range(4))
    for p in T.PAIRS:
        ref = cell(p, 30, 50)
        for c in T.COMPONENTS:
            C[c][p] = ref[c]["C_worst"]
            pdep[c][p] = "사전 의존" in ref[c]["flags"]
            mref[c][p] = ref[c]["MDE_worst"]
            mcorner[c][p] = [cell(p, int(s), int(n))[c]["MDE_worst"] for s, n in T.CORNER_CELLS]
    ids = {c: [p for p in T.PAIRS if C[c][p] >= T.C_ID and not pdep[c][p]] for c in T.COMPONENTS}
    n_id = {c: len(v) for c, v in ids.items()}
    if n_id != {"NDI": 1, "FISO": 0, "ODI": 28} or ids["NDI"] != [T.REF_PAIR]:
        sys.exit(f"STOP fig3: identified {n_id} NDI {ids['NDI']}")
    d = dict(pairs=T.PAIRS, ref_pair=T.REF_PAIR, stratum={p: T.acq_type(p) for p in T.PAIRS}, C=C, prior_dep=pdep,
             ref_cell=T.REF_CELL, corner_cells=T.CORNER_CELLS, mde_ref=mref, mde_corner=mcorner)
    T.draw_fig3(d, template=False)
    print("fig3 ok", n_id, "prior-dependent", {c: sum(pdep[c].values()) for c in T.COMPONENTS})
    report("Fig3_S4_identifiability")


def fig5():
    j = json.loads(TONG_JUDGE.read_text(encoding="utf-8"))["main"]
    ta, tb, tc = j["T_A"], j["T_B"], j["T_C"]
    if ta["n_pairs"] != 45 or ta["n_determined"] + ta["n_undetermined"] != ta["n_pairs"]:
        sys.exit(f"STOP fig5: T-A counts {ta}")
    if ta["undetermined_limit"] != T.TA_UNDET_MIN or ta["undetermined_limit"] != ta["n_pairs"] // 2 + 1:
        sys.exit(f"STOP fig5: limit {ta['undetermined_limit']}")
    same(ta["margin"], 1.96 * math.sqrt(2) * max(j["s_xy_unit_median"]), "fig5 margin")
    same(ta["p_random"], T.D_SHARE, "fig5 D share", 1e-9)
    if ta["judgeable"] or ta["verdict"] != "판정 불가" or ta["p"] is not None or tc != "판정 불가":
        sys.exit("STOP fig5: stored values allow only the undeterminable branch of T-A")
    if ta["n_undetermined"] < ta["undetermined_limit"]:
        sys.exit("STOP fig5: T-A stored as undeterminable with fewer pairs within the margin than the limit")
    tbd = {}
    for k in ("x", "y", "ODI"):
        b = tb[k]
        if not b["judgeable"] or b["n_den"] != b["n_total"] or b["units_dropped"]:
            sys.exit(f"STOP fig5: T-B {k} {b}")
        same(b["ratio"], b["n_retest"] / b["n_den"], f"fig5 T-B {k} share")
        if (b["ratio"] >= 0.90) != (b["verdict"] == "재검사 수준"):
            sys.exit(f"STOP fig5: T-B {k} verdict {b['verdict']} at {b['ratio']}")
        tbd[k] = (b["n_retest"], b["n_den"], b["ratio"])
    d = dict(n_pairs=ta["n_pairs"], n_within=ta["n_undetermined"], margin=ta["margin"],
             limit=ta["undetermined_limit"], d_share=ta["p_random"], tb=tbd, tc=None)
    T.draw_fig5(d, template=False)
    print("fig5 ok", {k: v[:2] for k, v in tbd.items()}, "T-A within margin", ta["n_undetermined"])
    report("Fig5_Tong_negative_control")


def figS2():
    sys.path.insert(0, str(REPO / "decomp" / "remeasure"))
    import judge_wong as JW
    data, own, missing = JW.load()
    if missing:
        sys.exit(f"STOP figS2: missing {missing}")
    stored = json.loads(WONG_JUDGE.read_text(encoding="utf-8"))["main"]
    judge = stored["by_subject"]
    TEs = np.array(JW.TES, float)
    tc = TEs - JW.TE_C
    rep_te = np.array([62.0, 132.0])
    subs = sorted(data)
    out = {}
    for ci, co in enumerate(("x", "y")):
        units, lin_all = [], []
        for si, sb in enumerate(subs):
            ser = data[sb]
            srt, nd = JW.s_rt(ser, ci)
            same(srt, judge[sb][co]["s_rt"], f"figS2 {sb} {co} s_rt")
            if nd != judge[sb][co]["n_diffs"]:
                sys.exit(f"STOP figS2: {sb} {co} n_diffs {nd}")
            us = sorted(set.intersection(*[set(ser[f"TE{t}"]) for t in JW.TES]))
            lin_flags, bratio = [], []
            for u in us:
                v = np.array([ser[f"TE{t}"][u][ci] for t in JW.TES])
                rep = np.array([ser[r][u][ci] for r in ("TE62R2", "TE132R2")])
                f = JW.fit_unit(TEs, v)
                ratio = f["curvature"] / (2 * math.sqrt(2) * srt)
                b = {nm: [t for t in bt if f"TE{t}" in ser] for nm, bt in JW.BUNDLES.items()}
                icpt = tuple(JW.intercept0(b[nm], [ser[f"TE{t}"][u][ci] for t in b[nm]]) for nm in ("낮은 둘", "높은 둘", "전체"))
                units.append(dict(subject=si, v=v, rep=rep, lin=np.array([f["slope"], f["level97"]]), ratio=ratio,
                                  icpt=icpt, single=(float(v.min()), float(v.max()))))
                lin_flags.append(f["curvature"] <= 2 * math.sqrt(2) * srt)
                bratio.append((max(icpt) - min(icpt)) / (v.max() - v.min()))
            if len(us) != judge[sb][co]["n_units"]:
                sys.exit(f"STOP figS2: {sb} {co} n_units {len(us)}")
            same(float(np.mean(lin_flags)), judge[sb][co]["linear_ratio"], f"figS2 {sb} {co} linear share")
            br = judge[sb][co]["bundle_range_over_te_range"]
            for q, val in (("median", np.median(bratio)), ("q25", np.quantile(bratio, 0.25)), ("q75", np.quantile(bratio, 0.75))):
                same(float(val), br[q], f"figS2 {sb} {co} intercept range {q}")
            lin_all += lin_flags
        same(float(np.mean(lin_all)), stored["overall"][co]["linear_ratio"], f"figS2 overall {co} linear share")
        same(float(np.mean([u["ratio"] <= 1 for u in units])), stored["overall"][co]["linear_ratio"],
             f"figS2 overall {co} share at ratio 1 (panel d)")
        out[co] = units
    T.draw_figS1(dict(TE=TEs, TE_rep=rep_te, units=out), template=False)
    print("figS2 ok", {co: len(u) for co, u in out.items()},
          {co: round(stored["overall"][co]["linear_ratio"], 4) for co in ("x", "y")})
    report("FigS2_Wong_multiTE")


if __name__ == "__main__":
    which = sys.argv[1:] or ["fig3", "fig5", "figS2"]
    for w in which:
        {"fig3": fig3, "fig5": fig5, "figS2": figS2}[w]()
