"""v17k (R49): dependence of the relaxation sensitivity value on the assumed control T2, and the isotropic-T2 path of apparent FISO
(review direction S91, items 6 and 7; exploratory; signal model only; no new data).
    python3 -I -B s93_values.py [--plant | --plant2]   -> s93_values.json in this folder
Input: the first row of Table 3 in s92_values.json (Collorone 2020: control and patient NDI, TE; published values).
(a) Extra-neurite lengthening (%) over the prior of control compartmental T2 (Supplementary Methods S-M2: T2,in 70-95 ms,
    r = T2,in/T2,en 1.0-1.8), i.e. at T2,en 38.9 and 95 ms.
(b) Minimum lengthening of both compartments by one ratio over the same prior grid, and the bound r* = 1 - s*T2,in
    (T2,in 70 and 95 ms) below which no ratio gives d.
(c) The change of the intra- to extra-neurite proton-density ratio that gives the same d, e^d - 1 (independent of TE).
(d) Apparent FISO 0.080 (ledger N-509) when the isotropic T2 of 115 ms (N-1074) is 10-20 % longer, at TE 70-100 ms
    (tissue weights unchanged; equal T1 and proton density; odds multiplied as in Eq. (3)).
Known answers (explicit failures): at the central T2 (82.5 and 52 ms), (a) and (b) equal the Table 3 values of s92 within 1e-9;
(d) with no lengthening returns 0.080. --plant uses d*TE instead of d/TE in (a) and must fail; --plant2 shortens T2 in (d) and must fail.
"""
import json, math, pathlib, sys

HERE = pathlib.Path(__file__).resolve().parent
PLANT, PLANT2 = "--plant" in sys.argv, "--plant2" in sys.argv
logit = lambda p: math.log(p / (1.0 - p))
fails = []
s92 = json.loads((HERE / "s92_values.json").read_text(encoding="utf-8"))["raw"]
NC, NP, TE = s92["t3_0_ndi_con_min"], s92["t3_0_ndi_pat_min"], s92["t3_0_te"]
d = logit(NP) - logit(NC)
s = (d * TE) if PLANT else d / TE                                   # per ms


def ext(t2en):
    return (1.0 / (1.0 / t2en + s) / t2en - 1.0) * 100.0


def both(t2in, t2en):
    q = 1.0 + s / (1.0 / t2en - 1.0 / t2in)
    return None if q <= 0 else (1.0 / q - 1.0) * 100.0


# known answers: the central values of Table 3
for name, mine, ref in (("ext", ext(52.0), s92["t3_0_ext_min"]), ("both", both(82.5, 52.0), s92["t3_0_both_min"])):
    if mine is None or abs(mine - ref) > 1e-9:
        fails.append(f"known answer {name}: {mine} != s92 {ref}")
if fails:
    raise SystemExit("STOP s93: " + "; ".join(fails))
T2IN = [70.0 + 0.25 * i for i in range(101)]                        # 70-95 ms
R = [1.0 + 0.002 * i for i in range(401)]                             # 1.0-1.8
ext_vals = [ext(ti / r) for ti in T2IN for r in R]
both_vals = [b for ti in T2IN for r in R if r > 1.0 for b in [both(ti, ti / r)] if b is not None]
rstar = {ti: 1.0 - s * ti for ti in (70.0, 95.0)}
raw = dict(d=d, s_per_ms=s, te=TE, ndi_con=NC, ndi_pat=NP,
           ext_prior_min=min(ext_vals), ext_prior_max=max(ext_vals), ext_t2en_min=70.0 / 1.8, ext_t2en_max=95.0,
           both_prior_min=min(both_vals), r_nosol_t2in70=rstar[70.0], r_nosol_t2in95=rstar[95.0],
           pd_pct=(math.exp(d) - 1.0) * 100.0, pd_pct_abs=abs((math.exp(d) - 1.0) * 100.0))
if not (raw["ext_prior_min"] < s92["t3_0_ext_min"] < raw["ext_prior_max"]):
    fails.append("central extra-neurite value outside the prior range")
# (d) isotropic-T2 path
F0, T2ISO = 0.080, 115.0


def fiso(te, p):
    t2 = T2ISO * ((1.0 - p) if PLANT2 else (1.0 + p))
    odds = F0 / (1.0 - F0) * math.exp(-te / t2) / math.exp(-te / T2ISO)
    return odds / (1.0 + odds)


if abs(fiso(85.0, 0.0) - F0) > 1e-15:
    fails.append("known answer FISO: not 0.080 without lengthening")
fv = {f"te{te}_p{int(p * 100)}": fiso(te, p) for te in (70.0, 100.0) for p in (0.10, 0.20)}
if not all(v > F0 for v in fv.values()):
    fails.append("apparent FISO did not rise with a longer isotropic T2")
raw.update(fiso_f0=F0, fiso_t2iso=T2ISO, fiso_p_min_pct=10.0, fiso_p_max_pct=20.0, fiso_min=min(fv.values()), fiso_max=max(fv.values()),
           **{f"fiso_{k}": v for k, v in fv.items()})
if fails:
    raise SystemExit("STOP s93: " + "; ".join(fails))
out = dict(source="s93_values.py (v17k · S91 6 · 7; exploratory; signal model only; inputs from s92_values.json)", raw=raw, fails=fails)
if not (PLANT or PLANT2):
    (HERE / "s93_values.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
for k in ("ext_prior_min", "ext_prior_max", "both_prior_min", "r_nosol_t2in70", "r_nosol_t2in95", "pd_pct", "fiso_min", "fiso_max"):
    print(k, round(raw[k], 6))
