"""Builds the d4_2 forward sensitivity table: how apparent values move with r, T2_in, T1 ratio, TE and TR, using repository formulas and no data."""
import datetime, hashlib, json, math, pathlib, sys
import numpy as np
REPO = pathlib.Path("<REPO>")
sys.path.insert(0, str(REPO / "decomp" / "stage1"))
import config as C  # noqa: E402
OUT = pathlib.Path(__file__).resolve().parent
T2ISO = 1000.0                                   # S4 (s4cfg.T2ISO) · manuscript «free-water T2 1000 ms»
R_GRID = [1.0, 1.3, 1.6, 1.8]                    # Grid of the S4 prior U(1.0, 1.8)
REF = {"Gong 82/52 ms(참고)": (82.0, 52.0), "E1 함축 r 2.64(참고 · T2_in 82.5 로 놓음)": (82.5, 82.5 / 2.64)}
T2IN_GRID = [70.0, 82.5, 95.0]                   # S4 prior U(70, 95)
T1C, T1RATIO = [850.0, 975.0, 1100.0], {"i 1.3": 1.3, "ii 1/1.3": 1 / 1.3, "iii 1.0": 1.0}
NU0, F0 = 0.6, 0.10
TE = np.arange(60.0, 132.0 + 1e-9, 2.0)


def t1s(c, ratio):
    return np.array([c * math.sqrt(ratio), c / math.sqrt(ratio), C.T1_ISO])


def xy(nu, f, te, tr, t2in, t2en, t1):
    a = C.weights(np.asarray(te, float)[..., None], np.array([t2in, t2en, T2ISO]), tr, t1)
    nua, fa = C.apparent(nu, f, (a[..., 0], a[..., 1], a[..., 2]))
    return C.frac_coords(nua, fa), (nua, fa), a


def smean(b, dpar, dper):
    """Spherical mean — b (ms/µm²) · d (µm²/ms)."""
    stick = math.sqrt(math.pi / (4 * b * dpar)) * math.erf(math.sqrt(b * dpar))
    dd = dpar - dper
    cyl = math.exp(-b * dper) * (math.sqrt(math.pi / (4 * b * dd)) * math.erf(math.sqrt(b * dd)) if dd > 1e-12 else 1.0)
    return stick, cyl, math.exp(-b * C.DISO)


def s_of_te(bmm, te, tr, t2in, t2en, t1, nu, f, gen, dpar):
    b = bmm / 1000.0
    out = []
    for t in np.atleast_1d(te):
        a = C.weights(np.array(t), np.array([t2in, t2en, T2ISO]), tr, t1)
        w = np.array([(1 - f) * nu * a[0], (1 - f) * (1 - nu) * a[1], f * a[2]])
        nua = w[0] / (w[0] + w[1])
        A = smean(b, dpar, C.dperp(gen, dpar, nu, nua))
        out.append(math.log(float(np.dot(w, A) / w.sum())))
    return np.array(out)


res = {"what": __doc__.splitlines()[0], "made": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
       "script_sha256": hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
       "config_sha256": hashlib.sha256((REPO / "decomp/stage1/config.py").read_bytes()).hexdigest(),
       "fixed": dict(T2ISO=T2ISO, T1_ISO=C.T1_ISO, DISO=C.DISO, nu=NU0, f_iso=F0)}
# (A) TE slope of x′ = logit ν′: analytically −(1/T2_in − 1/T2_en), independent of ν and f, and exactly linear in TE (numerical check)
A = []
for t2in in T2IN_GRID:
    for r in R_GRID:
        (x, _), _, _ = xy(NU0, F0, TE, 7000.0, t2in, t2in / r, t1s(975.0, 1.0))
        sl = np.polyfit(TE, x, 1)[0]
        lin = float(np.max(np.abs(x - np.polyval(np.polyfit(TE, x, 1), TE))))
        A.append(dict(T2_in=t2in, r=r, T2_en=round(t2in / r, 2), dx_per_10ms=round(10 * sl, 4), analytic=round(-10 * (1 / t2in - r / t2in), 4), max_nonlinearity=lin))
for k, (t2in, t2en) in REF.items():
    (x, _), _, _ = xy(NU0, F0, TE, 7000.0, t2in, t2en, t1s(975.0, 1.0))
    A.append(dict(ref=k, T2_in=t2in, T2_en=round(t2en, 2), r=round(t2in / t2en, 3), dx_per_10ms=round(10 * np.polyfit(TE, x, 1)[0], 4)))
res["A_x_TE_slope"] = A
# (B) TE slope of y′ = logit f′ − ln ν′ (at each location) · departs from a straight line (62 – 132 ms) — per tissue and per r
B = []
for nu in (0.5, 0.6, 0.7):
    for f in (0.05, 0.10, 0.20):
        for r in R_GRID:
            (_, y), _, _ = xy(nu, f, TE, 7000.0, 82.5, 82.5 / r, t1s(975.0, 1.0))
            g = np.gradient(y, TE)
            sel = TE >= 62
            fit = np.polyfit(TE[sel], y[sel], 1)
            B.append(dict(nu=nu, f_iso=f, r=r, dy_per_10ms_at={str(int(t)): round(10 * float(g[TE == t][0]), 4) for t in (62, 80, 100, 132)},
                          max_nonlinearity_62_132=round(float(np.max(np.abs(y[sel] - np.polyval(fit, TE[sel])))), 4)))
res["B_y_TE"] = B
# (C) Share of TR — ∂x′/∂TR · ∂y′/∂TR (per 1 s) · TE 80 · T1-ratio scenarios · common T1 c
Cc = []
for name, ratio in T1RATIO.items():
    for c in T1C:
        for tr in (3000.0, 7000.0, 11000.0):
            h = 1.0
            (x1, y1), _, _ = xy(NU0, F0, 80.0, tr - h, 82.5, 82.5 / 1.3, t1s(c, ratio))
            (x2, y2), _, _ = xy(NU0, F0, 80.0, tr + h, 82.5, 82.5 / 1.3, t1s(c, ratio))
            Cc.append(dict(T1_ratio=name, T1_c=c, TR=tr, dx_per_s=round(float((x2 - x1) / (2 * h) * 1000), 5), dy_per_s=round(float((y2 - y1) / (2 * h) * 1000), 5)))
res["C_TR"] = Cc
# (D) TE slope of the direction-averaged s(b, TE) — Wong (b 700 · 2000 · TE 62 – 132) · ON-Harmony (b 1000 · 2000 · TE 60 – 110) · three generators · three d∥ values
D = []
for gen in C.GENERATORS:
    for dpar in (1.5, 1.7, 2.0):
        for r in R_GRID + [None]:
            t2in, t2en = (82.5, 82.5 / r) if r else REF["Gong 82/52 ms(참고)"]
            row = dict(generator=gen, d_par=dpar, r=(r if r else "Gong 82/52"))
            for lab, bs, lo, hi in (("Wong", (700, 2000), 62, 132), ("ON-Harmony", (1000, 2000), 60, 110)):
                te = np.arange(lo, hi + 1e-9, 2.0)
                for bmm in bs:
                    s = s_of_te(bmm, te, 7000.0, t2in, t2en, t1s(975.0, 1.0), NU0, F0, gen, dpar)
                    row[f"{lab}_b{bmm}_ds_per_10ms"] = round(10 * float(np.polyfit(te, s, 1)[0]), 5)
            D.append(row)
res["D_s_TE"] = D
(OUT / "d4_2_forward.json").write_text(json.dumps(res, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
print("A", [(a.get("T2_in"), a.get("r"), a["dx_per_10ms"]) for a in A][:6])
print("A ref", [(a["ref"], a["r"], a["dx_per_10ms"]) for a in A if "ref" in a])
print("A max nonlin", max(a.get("max_nonlinearity", 0) for a in A))
print("B nonlin range", min(b["max_nonlinearity_62_132"] for b in B), max(b["max_nonlinearity_62_132"] for b in B))
print("C |dx/dTR| max", max(abs(c["dx_per_s"]) for c in Cc), "at TR 3000 ratio i:", [c["dx_per_s"] for c in Cc if c["TR"] == 3000.0 and c["T1_ratio"] == "i 1.3"])
print("D T-in 1.7:", [(d["r"], d["Wong_b700_ds_per_10ms"], d["Wong_b2000_ds_per_10ms"]) for d in D if d["generator"] == "T-in" and d["d_par"] == 1.7])
