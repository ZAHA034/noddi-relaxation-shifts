"""Supplementary Methods S-M9 and Table S12 (S-M10 and Table S17 until v16z5, 2026-11-17): every result number in the text of S-M9 is recomputed from the stored result files
of the exploratory analyses after the verdict and compared, in context, with the text; Table S12 is regenerated and compared line by line.
Nothing is fitted; one algebraic value (the raw-scale change of fraction-only groups) is recomputed with rules_env.model from the
ROI-matched ON-Harmony states, as it was for the text.

Run from the repository root (or set RULES_ENV_SRC to the folder that holds rules_env):
  PYTHONDONTWRITEBYTECODE=1 RULES_ENV_REPO=<repo> <repo>/.venv/bin/python -I -B decomp/manuscript/checks/check_env_numbers.py \
      [supplement.md] [results folder]
Exit status 0 only if every fragment is found and Table S12 is identical.
"""
import json
import math
import os
import pathlib
import sys

sys.dont_write_bytecode = True
HERE = pathlib.Path(__file__).resolve().parent
ROOT = pathlib.Path(os.environ.get("RULES_ENV_REPO") or HERE.parents[2]).resolve()
SUPP = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "decomp" / "manuscript" / "supplement.md"
RES = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "decomp" / "results" / "rules_env"
sys.path.insert(0, os.environ.get("RULES_ENV_SRC") or str(ROOT / "decomp"))

import numpy as np  # noqa: E402

from rules_env import analyses as AS  # noqa: E402
from rules_env import metrics as MT  # noqa: E402
from rules_env import model as M  # noqa: E402

J = {n: json.loads((RES / f"{n}.json").read_text(encoding="utf-8")) for n in ("s1a", "s1b", "s2", "k8", "e1", "e2", "e3", "e4")}
for n in ("e1_fix", "e5", "e6", "s3"):
    if (RES / f"{n}.json").exists():
        J[n] = json.loads((RES / f"{n}.json").read_text(encoding="utf-8"))
TEXT = SUPP.read_text(encoding="utf-8")
SM10 = TEXT[TEXT.index("### S-M9."):TEXT.index("### S-M10.") if "### S-M10." in TEXT else TEXT.index("## Supplementary tables")]   # 2026-11-17 v16z7b: S-M10 (proofs) comes after S-M9
HAS_R2 = "#### S-M9.5" in TEXT
if HAS_R2 and not all(n in J for n in ("e5", "e6", "s3")):          # 2026-11-15 (F1): the R2 text must not go unchecked
    sys.exit("S-M9.5 is in the supplement but e5/e6/s3 results are missing: the R2 numbers cannot be checked")
if "e1_fix" not in J:
    sys.exit("e1_fix.json is missing: the pooled ΔTE* of S-M9.4 cannot be checked")


def F(v, d, plus=False):
    s = f"{abs(v):.{d}f}"
    neg = v < 0 and float(s) != 0.0
    return ("−" if neg else ("+" if plus else "")) + s


def R(vals, d, plus=False, sep="–"):
    return f"{F(min(vals), d, plus)}{sep}{F(max(vals), d, plus)}"


GENS = ("T-in", "T-app", "T-none")
LOC = ("none", "ctrl_offset", "combat_loc", "combat_train_ctrl_loc")
S2 = J["s2"]["by"]
frags, conds = [], []

# ── fitter output (S1a · S1b)
st = J["s1a"]["settings"].values()
frags += [f"the offset-subtracted share was {F(min(v['ors'] for v in st), 3)} in every pair class",
          f"the common shift ranged from {F(min(v['offset'] for v in st), 3)} to {F(max(v['offset'] for v in st), 3, True)}",
          f"the proportional bias from {F(min(v['prop_bias'] for v in st), 3)} to {F(max(v['prop_bias'] for v in st), 3, True)}",
          f"the lowest rank correlation was {F(min(v['spearman'] for v in st), 3)}"]
conds.append(("no rank reversal (S1a)", max(v["rank_flip"] for v in st) == 0))
ob = J["s1b"]["onh"]["by"]
pairs = {k.split("|")[0] for k in ob}
by_g = {g: [v for k, v in ob.items() if k.split("|")[1] == g] for g in GENS}
conds += [("28 ON-Harmony pairs (S1b)", len(pairs) == 28),
          ("all retest level, T-in and T-app (S1b)", all(v["retest_level"] is True for g in ("T-in", "T-app") for v in by_g[g]))]
frags.append(f"(lowest share {F(min(v['ors'] for v in by_g['T-in']), 3)} and {F(min(v['ors'] for v in by_g['T-app']), 3)}) and in "
             f"{F(np.mean([v['retest_level'] is True for v in by_g['T-none']]), 2)} of setting-pair cells with fixed tortuosity (lowest {F(min(v['ors'] for v in by_g['T-none']), 3)})")
wb = J["s1b"]["wong"]["by"]
frags.append(f"of {F(wb['x|T-in|L0|iii']['tau_fit'], 5)} per ms (true-fraction tortuosity, central relaxation; {R([v['tau_fit'] for k, v in wb.items() if k.startswith('x|T-in|')], 5)} across settings), "
             f"against {F(J['e1']['pooled']['x']['tau'], 4)} per ms observed")

# ── Proposition 4 at fitter output (S2)
def ratios(sc, gens=GENS):
    return [S2[f"pop|x|two|{sc}|{g}|{m}"]["ratio"] for g in gens for m in LOC]


cells = {f"pop|x|two|{sc}|{g}|{m}": S2[f"pop|x|two|{sc}|{g}|{m}"] for sc in ("fig3", "l2en", "l2both") for g in GENS for m in LOC}
out = sorted((v for v in cells.values() if not 0.80 <= v["ratio"] <= 1.25), key=lambda v: -v["ratio"])
frags.append(f"was {R(ratios('fig3'), 2)} times Eq. (7) for the relaxation of Fig. 2, {R(ratios('l2en'), 2)} for diffuse pathology of the extra-neurite compartment "
             f"and {R(ratios('l2both'), 2)} for diffuse pathology of both tissue compartments")
frags.append(f"({F(out[0]['ratio'], 3)}, {F(out[1]['ratio'], 3)} and {F(out[2]['ratio'], 3)}), where Eq. (7) was {F(out[0]['eq7'], 3)} and the fitted change "
             f"{R([v['interaction'] for v in out], 3)}" if len({F(v['interaction'], 3) for v in out}) > 1 else
             f"({F(out[0]['ratio'], 3)}, {F(out[1]['ratio'], 3)} and {F(out[2]['ratio'], 3)}), where Eq. (7) was {F(out[0]['eq7'], 3)} and the fitted change {F(out[0]['interaction'], 3)}")
conds += [("three cells outside, all l2both × T-none", len(out) == 3 and all("|l2both|T-none|" in k for k, v in cells.items() if not 0.80 <= v["ratio"] <= 1.25)),
          (f"all {len(cells)} cells ratio ≥ 1", len(cells) == 36 and all(v["ratio"] >= 1 for v in cells.values()))]
e8 = [S2[f"pop|x|onh8|{sc}|T-in|{m}"] for sc in ("fig3", "l2en", "l2both") for m in LOC]
fr8 = [S2[f"pop|x|onh8|frac|T-in|{m}"] for m in LOC]
frags.append(f"with a slope {R([v['slope'] / v['slope_alg'] for v in e8], 2)} times the algebraic slope (true-fraction tortuosity, the four location terms), and after one group × TE term "
             f"the residual SD across sites ({R([v['resid_sd'] for v in e8], 4)}) equaled that of the scenario with no relaxation difference "
             f"({R([v['resid_sd'] for v in fr8], 4) if len({F(v['resid_sd'], 4) for v in fr8}) > 1 else F(fr8[0]['resid_sd'], 4)})")
fr_raw, fr_x = S2["pop|raw_ndi|two|frac|T-in|none"], S2["pop|x|two|frac|T-in|none"]
EM_, st_, o_, s_ = AS.onh_ctx()
om = AS.onh_matched_units(EM_, st_, o_)
t2c, t1c = AS.centre_eta()


def alg_raw(te, dnu=-0.03):
    w = M.weights(te, t2c, 3600.0, t1c)
    return float(np.mean([M.apparent(nu + dnu, ff, w)[0] - M.apparent(nu, ff, w)[0] for nu, ff, _ in om["units"].values()]))


frags.append(f"changed between TE 70 and 100 ms by {F(fr_raw['interaction'], 4, True)} (2.5–97.5 % of replicates {F(fr_raw['q025'], 4, True)} to {F(fr_raw['q975'], 4, True)}; true-fraction tortuosity), "
             f"against {F(alg_raw(100.0) - alg_raw(70.0), 4, True)} from the algebra; on the logit scale the same change was {F(fr_x['interaction'], 4, True)} "
             f"({F(fr_x['q025'], 4, True)} to {F(fr_x['q975'], 4, True)})")
l1 = [S2[f"l1|x|two|{sc}|{g}|{m}"] for sc in ("fig3", "l2en", "l2both") for g in GENS for m in LOC]
frags.append(f"differed from the algebra of the same participants by {F(min(v['fit_minus_algebra'] for v in l1), 3)} to {F(max(v['fit_minus_algebra'] for v in l1), 3, True)} (scenarios, generators and location terms), "
             f"whereas that algebra differed from Eq. (7) by up to {F(max(abs(v['algebra_same_people'] - v['eq7']) for v in l1), 3)}")

# ── K8
k8 = J["k8"]
frags.append(f"predicted {k8['T-in|n_fitted_predicted']}, {k8['T-app|n_fitted_predicted']} and {k8['T-none|n_fitted_predicted']} pairs with true-fraction, apparent-fraction and fixed tortuosity, of which "
             f"{k8['T-in|fitted_agree_observed']}, {k8['T-app|fitted_agree_observed']} and {k8['T-none|fitted_agree_observed']} agreed with the observed signs, against "
             f"{k8['algebra_agree_observed']} of {k8['n_algebra_predicted']}")

# ── E1
e1 = J["e1"]
bx = [v for k, v in e1["by"].items() if k.startswith("x|")]
px, py = J["e1_fix"]["pooled"]["x"]["primary"], J["e1_fix"]["pooled"]["y"]["primary"]    # corrected pooled values (ANALYSIS_PLAN §8)
px_t, py_t = e1["pooled"]["x"], e1["pooled"]["y"]
s_on = float(J["s1a"]["s_x_median"])
frags.append(f"was {F(px_t['tau'], 4)} per ms (apparent ΔR₂ SD {F(px_t['delta_r2_sd_per_s'], 1)} s⁻¹) and that of y {F(py_t['tau'], 4)} per ms, against repeat SDs of "
             f"{R([v['s_rt_corrected'] for v in bx], 3)} for x; ΔTE* was {F(px['dte_star'], 1)} ms for x (90 % interval {R(px['dte_ci90'], 1)}) and {F(py['dte_star'], 1)} ms for y "
             f"({R(py['dte_ci90'], 1)}), {R([v['dte_star'] for v in bx], 1)} ms for x by participant, and {round(MT.te_tolerance(np.full(27, s_on), px['tau']))} ms with the ON-Harmony retest SD of x")
conds.append(("β intervals exclude 0 (x, all participants)", len(bx) == 3 and all(v["beta_ci95"][0] > 0 for v in bx)))
frags.append(f"(β {R([v['beta'] for v in bx], 4, True, ' to ')} per ms per unit of x; 95 % intervals excluding 0) and was nearly unchanged when the FISO level was subtracted "
             f"({R([v['fiso_partial']['beta'] for v in bx], 4, True, ' to ')}); at TE 62 against 132 ms, the observed share ({R([v['ors_62_132_observed'] for v in bx], 2)}) "
             f"matched the share predicted from τ ({R([v['ors_62_132_predicted'] for v in bx], 2)})")
frags.append(f"Residual SDs of the linear fits of x ({R([v['resid_sigma'] for v in bx], 3)})")
fx = J["e1_fix"]["pooled"]
conds += [("e1_fix regression gates", J["e1_fix"]["regression_ok"] is True and fx["x"]["as_run_matches_stored"] and fx["y"]["as_run_matches_stored"])]
frags.append(f"the computation first run used the 81 pooled units and gave {F(fx['x']['as_run']['dte_star'], 1)} ms ({R(fx['x']['as_run']['dte_ci90'], 1)}) for x and "
             f"{F(fx['y']['as_run']['dte_star'], 1)} ms ({R(fx['y']['as_run']['dte_ci90'], 1)}) for y")
frags.append(f"ΔTE* for x was {F(fx['x']['sens_df']['dte_star'], 1)}, {F(fx['x']['sens_median']['dte_star'], 1)} and {F(fx['x']['sens_both']['dte_star'], 1)} ms")

# ── E2–E4
w = J["e4"]["wong"]["TE62-TE132"]
frags.append(f"the offset-subtracted share was {F(w['raw_ndi'], 3)} for raw NDI and {F(w['x'], 3)} for x (difference {F(w['x-raw_ndi']['diff'], 3)}; "
             f"{F(w['x-raw_ndi']['lo'], 3)} to {F(w['x-raw_ndi']['hi'], 3)}), {F(w['logit_fiso'], 3)} for logit FISO and {F(w['y'], 3)} for y "
             f"({F(w['y-logit_fiso']['diff'], 3, True)}; {F(w['y-logit_fiso']['lo'], 3)} to {F(w['y-logit_fiso']['hi'], 3, True)})")
sm = J["e4"]["wong_scale_models_62_132"]
conds.append(("raw model better in all participants", all(v["ss_raw"] < v["ss_logit"] for v in sm.values())))
frags.append("(residual sums " + ", ".join(f"{F(sm[sb]['ss_raw'], 3)} against {F(sm[sb]['ss_logit'], 3)}" for sb in sorted(sm)[:-1]) +
             f" and {F(sm[sorted(sm)[-1]]['ss_raw'], 3)} against {F(sm[sorted(sm)[-1]]['ss_logit'], 3)})")
oth = [v["x-raw_ndi"] for k, v in J["e4"]["wong"].items() if k != "TE62-TE132"]
nb, na = sum(v["hi"] < 0 for v in oth), sum(v["lo"] > 0 for v in oth)
words = {0: "none", 1: "one"}
frags.append(f"across the other {len(oth)} TE pairs the x-minus-raw difference ranged from {F(min(v['diff'] for v in oth), 3)} to {F(max(v['diff'] for v in oth), 3, True)}, "
             f"with an interval below 0 in {words.get(nb, nb)} pairs and above 0 in {words.get(na, na)}")
md = J["e4"]["onh_median"]
onh = J["e4"]["onh"].values()
nb2, na2 = sum(v["x-raw_ndi"]["hi"] < 0 for v in onh), sum(v["x-raw_ndi"]["lo"] > 0 for v in onh)
tb = J["e2"]["tau_b_dTE_orsx"]
frags.append(f"the median shares were {F(md['raw_ndi'], 3)} for raw NDI, {F(md['x'], 3)} for x, {F(md['raw_fiso'], 3)} for raw FISO, {F(md['logit_fiso'], 3)} for logit FISO, "
             f"{F(md['y'], 3)} for y and {F(md['odi'], 3)} for ODI, with the x-minus-raw interval below 0 in {words.get(nb2, nb2)} pair and above 0 in {words.get(na2, na2)}; "
             f"the share for x was not related to |ΔTE| (τ_b {F(tb['tau'], 2)}, p {F(tb['p'], 2)})")
e3 = J["e3"]["by"]
frags.append("changed the shares from " + ", ".join(f"{F(e3[k]['raw'], 3)} to {F(e3[k]['ors'], 3)} ({nm})" for k, nm in
                                                   (("raw_ndi", "raw NDI"), ("x", "x"), ("raw_fiso", "raw FISO"), ("logit_fiso", "logit FISO"), ("y", "y"))) +
             f" and {F(e3['odi']['raw'], 3)} to {F(e3['odi']['ors'], 3)} (ODI)")
conds.append(("no unit clipped (E3)", J["e3"]["clipped_unit_sessions"] == 0))

# ── S-M9.5 (R2: E5 · E6 · S3) — merged into this script on 2026-11-15 (code review F1); required whenever S-M9.5 is present
if HAS_R2:
    e5, e6, s3 = J["e5"], J["e6"], J["s3"]
    fo = [v["first_order"] for v in e5["by"].values()]
    frags.append(f"To first order, R was {R([f['R'] for f in fo], 2)} against 2ν̄ − 1 of {R([f['T'] for f in fo], 2)}, so the raw shift grew by "
                 f"{R([f['raw_rate'] for f in fo], 2)} of its size per unit of x and the logit shift by {R([f['R'] for f in fo], 2)}")
    conds += [("E4 reproduced (E5)", e5["regression_ok"] is True), ("raw more uniform, all participants (E5)", all(f["raw_more_uniform"] for f in fo)),
              ("E5 reading (가)", e5["reading"] == "가")]
    p, p0 = e5["predicted"], e5["predicted_beta0"]
    frags.append(f"was {F(p['diff'], 3)} ({F(p['diff_q025'], 3)} to {F(p['diff_q975'], 3, True)} over simulated replicates), within the observed interval, against "
                 f"{F(p0['diff'], 3, True)} without the level dependence")
    pr = e6["primary"]
    other = {k: v for k, v in e6["by"].items() if v["reading"] != "가"}
    frags.append(f"the expected share was {F(pr['ors_x'], 3)} for x and {F(pr['ors_raw'], 3)} for raw NDI (true-fraction tortuosity, central relaxation; difference {F(pr['diff'], 3, True)}, "
                 f"{F(pr['diff_q025'], 3)} to {F(pr['diff_q975'], 3, True)} over replicates), {e6['reading_counts']['가']} of {len(e6['by'])} settings with the same reading")
    conds += [("E6 reading (가)", pr["reading"] == "가"), ("E6 other six all T-app, both below 0.90",
              len(other) == 6 and all(k.startswith("T-app|") and v["ors_x"] < 0.90 and v["ors_raw"] < 0.90 for k, v in other.items()))]
    r = lambda nc, dn, ta: s3["rows"][f"{nc}|{dn}|{ta}"]["raw"]
    frags.append(f"of {F(r(0.46, -0.1, 100), 4, True)} for a difference of 0.10 and less than 0.0001 for 0.20 (control NDI 0.46; {F(r(0.6, -0.1, 100), 4, True)} and "
                 f"{F(r(0.6, -0.2, 100), 4, True)} at 0.60)")
    frags.append(f"for a difference of 0.03, as in the Proposition 4 analysis, it was {F(r(0.46, -0.03, 100), 4, True)} and {F(r(0.6, -0.03, 100), 4, True)}")
    conds += [("S3 |raw| < 0.0001 at 0.46, −0.20, TE 100", abs(r(0.46, -0.2, 100)) < 1e-4), ("S3 logit change 0", s3["max_abs_logit"] < 1e-12)]

bad = [f for f in frags if f not in SM10]
badc = [n for n, ok in conds if not ok]
print(f"S-M9 fragments: {len(frags) - len(bad)} of {len(frags)} found; conditions: {len(conds) - len(badc)} of {len(conds)} hold")
for f in bad:
    print("  NOT FOUND:", f)
for n in badc:
    print("  FAILS:", n)

# ── Table S12, regenerated
mk = HERE / "make_S17.py"
if mk.exists():
    ns = {"__name__": "make_S17_check", "RES": RES}
    exec(compile(mk.read_text(encoding="utf-8"), str(mk), "exec"), ns)
    want = ns["table_lines"](J)
    have = TEXT[TEXT.index("### Table S12"):]
    have = have[:have.index("\n### ", 5)] if "\n### " in have[5:] else have
    lines = have.split("\n")
    start = next((i for i, ln in enumerate(lines) if ln == want[0]), None)
    got = lines[start:start + len(want)] if start is not None else []
    after = lines[start + len(want):] if start is not None else []
    extra = [ln for ln in after if ln.startswith("| ")]                     # 2026-11-15 (F5): no table row may follow the generated block
    miss = [w for w, g in zip(want, got + [None] * (len(want) - len(got))) if w != g]
    print(f"Table S12: {len(want) - len(miss)} of {len(want)} generated lines identical in order; extra rows after the block: {len(extra)}")
    for ln in miss[:10]:
        print("  DIFFERS:", (ln or "")[:160])
    for ln in extra[:5]:
        print("  EXTRA:", ln[:160])
    bad += miss + extra
sys.exit(1 if (bad or badc) else 0)
