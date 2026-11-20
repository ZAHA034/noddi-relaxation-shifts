"""The T2 path of the «Free water» paragraph of Section 5.2 (exploratory; existing code, CPU, no new data): how the diffuse tissue T2 lengthening of Section 2.4 moves apparent FISO."""
import hashlib
import json
import math
import pathlib
import platform
import sys
from decimal import Decimal, ROUND_HALF_UP

import numpy as np

sys.dont_write_bytecode = True
REPO = pathlib.Path(sys.argv[1]).resolve()
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(REPO / "decomp" / "stage1"))
import config as C  # noqa: E402

CFG = REPO / "decomp" / "stage1" / "config.py"
S4CFG = REPO / "decomp" / "s4" / "s4cfg.py"
assert "L2_T2, L2_T1 = 0.0984, 0.0441" in S4CFG.read_text(encoding="utf-8")
sha = lambda p: hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()  # noqa: E731

F0 = 0.080                                   # ledger N-509 (same starting point as the calculation of Proposition 2)
NU, T2IN, T2EN, T2ISO = 0.46, 82.5, 52.0, 1000.0   # same inputs as Section 2.4 · the residual table
K = 1.0984                                   # s4cfg L2_T2 = 0.0984 (diffuse · the same multiplier as Section 2.4)
TES = (70.0, 100.0)
CASES = {"both": (K, K), "extra_only": (1.0, K)}   # (T2,in multiplier, T2,en multiplier)


def tissue_weight(te, kin, ken):
    return NU * math.exp(-te / (T2IN * kin)) + (1 - NU) * math.exp(-te / (T2EN * ken))


def f_new(f, rho):
    return f / (f + (1 - f) * rho)


def r3(x):
    return float(Decimal(repr(x)).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP))


res = {}
for case, (kin, ken) in CASES.items():
    row = {}
    for te in TES:
        rho = tissue_weight(te, kin, ken) / tissue_weight(te, 1.0, 1.0)
        f1 = f_new(F0, rho)
        row[f"te{int(te)}"] = dict(TE_ms=te, rho=rho, fiso_app=f1, change=f1 - F0, printed=f"{r3(f1):.3f}")
    grid = [70.0 + 0.5 * k for k in range(61)]
    g = [f_new(F0, tissue_weight(t, kin, ken) / tissue_weight(t, 1.0, 1.0)) - F0 for t in grid]
    row["grid_te70_100_step0.5"] = dict(all_negative=all(v < 0 for v in g),
                                        abs_increases_with_te=all(abs(g[i + 1]) > abs(g[i]) for i in range(60)),
                                        min_change=min(g), max_change=max(g))
    assert row["grid_te70_100_step0.5"]["all_negative"] and row["grid_te70_100_step0.5"]["abs_increases_with_te"], (case, row)  # (3)
    assert abs(min(g) - row["te100"]["change"]) < 1e-15 and abs(max(g) - row["te70"]["change"]) < 1e-15
    res[case] = row


# (1) directly with Eq. (3) of config.py
def solve_true_fiso(a):
    lo, hi = 1e-9, 0.999
    for _ in range(200):
        mid = (lo + hi) / 2
        if C.apparent(NU, mid, a)[1] < F0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


chk = {}
for t1 in (850.0, 1100.0):
    for case, (kin, ken) in CASES.items():
        for te in TES:
            a0 = C.weights(te, (T2IN, T2EN, T2ISO), TR=3600.0, t1=(t1, t1, C.T1_ISO))
            a1 = C.weights(te, (T2IN * kin, T2EN * ken, T2ISO), TR=3600.0, t1=(t1, t1, C.T1_ISO))
            ft = solve_true_fiso(a0)
            assert abs(C.apparent(NU, ft, a0)[1] - F0) < 1e-12
            direct = float(C.apparent(NU, ft, a1)[1])
            closed = res[case][f"te{int(te)}"]["fiso_app"]
            assert abs(direct - closed) < 1e-12, (t1, case, te, direct, closed)
            chk[f"T1{int(t1)}|{case}|TE{int(te)}"] = dict(true_fiso=ft, direct=direct, closed=closed, diff=direct - closed)

# (2) the T1 path (Proposition 2 · diffuse) on the same path — the same value as the closed form λ of s75
t1_path = {}
for t1 in (850.0, 1100.0):
    lam = (1 - math.exp(-3600.0 / (1.0441 * t1))) / (1 - math.exp(-3600.0 / t1))
    a0 = C.weights(0.0, (T2IN, T2EN, T2ISO), TR=3600.0, t1=(t1, t1, C.T1_ISO))
    a1 = C.weights(0.0, (T2IN, T2EN, T2ISO), TR=3600.0, t1=(t1 * 1.0441, t1 * 1.0441, C.T1_ISO))
    ft = solve_true_fiso(a0)
    direct = float(C.apparent(NU, ft, a1)[1])
    assert abs(direct - f_new(F0, lam)) < 1e-12, (t1, direct, f_new(F0, lam))
    t1_path[f"T1{int(t1)}"] = dict(fiso_app=direct, change=direct - F0)

out = dict(
    what=("Size and sign of the tissue-T2 path on apparent FISO for Section 5.2 (2026-11-17; exploratory: computed after closely related values had "
          "been seen; existing code; no data read): a diffuse lengthening of tissue T2 by 9.84 % (both tissue compartments, or the extra-neurite "
          "compartment only; assignment assumed), starting from an apparent FISO of 0.080 as for Proposition 2, at TE 70-100 ms; tissue compartments share T1"),
    model="decomp/stage1/config.py weights and apparent (Eq. (3)); closed form f' = f / (f + (1 - f) rho), rho = ratio of tissue weights",
    inputs=dict(f0_apparent=F0, nu=NU, T2_in_ms=T2IN, T2_en_ms=T2EN, T2_iso_ms=T2ISO, factor=K, TE_ms=list(TES), TR_ms_for_check=3600.0,
                T1_tissue_ms_for_check=[850.0, 1100.0], T1_iso_ms_for_check=C.T1_ISO),
    results=res, check_direct_eq3=chk, t1_path_same_code=t1_path,
    sources={"config.py": sha(CFG), "s4cfg.py": sha(S4CFG)},
    python=platform.python_version(), numpy=np.__version__,
)
out["printed"] = dict(     # the text says «lowers it by …» — the decrease of apparent FISO (positive)
    t2_factor_pct=round((K - 1) * 100, 10),
    both_te70_drop=-res["both"]["te70"]["change"], both_te100_drop=-res["both"]["te100"]["change"],
    en_te70_drop=-res["extra_only"]["te70"]["change"], en_te100_drop=-res["extra_only"]["te100"]["change"],
)
(HERE / "s76_values.json").write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
for case, row in res.items():
    print(case, " · ".join(f"TE {int(row[k]['TE_ms'])}: {row[k]['fiso_app']:.6f} ({row[k]['change']:+.6f}; {row[k]['printed']})" for k in ("te70", "te100")))
print("T1 path (same code):", {k: round(v["fiso_app"], 6) for k, v in t1_path.items()})
print("direct check max |diff|:", max(abs(v["diff"]) for v in chk.values()))
