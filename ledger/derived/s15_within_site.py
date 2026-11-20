"""Within-site column of the supplementary residual table: patient − control difference of apparent NDI at the same TE in five scenarios (Eqs. (3) · (4); derived, no data), by the same code path and inputs as residual_s15.py."""
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
sys.path.insert(0, str(REPO / "decomp" / "stage1"))
import config as C  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
DERIVED = HERE                                              # x15/ledger/derived (where this file is) — read only
S15_JSON, S15_PY = DERIVED / "residual_s15.json", DERIVED / "residual_s15.py"
S74_JSON, S74_PY = DERIVED / "s74_clinical.json", DERIVED / "s74_clinical.py"
CFG = REPO / "decomp" / "stage1" / "config.py"

NU, T2IN_C, T2EN_C, T2ISO = 0.46, 82.5, 52.0, 1000.0       # same as residual_s15.py
L2 = 1.0984                                                 # diffuse multiplier of ROWS in residual_s15.py
TES = (70.0, 80.0, 90.0, 100.0)
ROWS = [  # same order · same formula as ROWS in residual_s15.py (key, patient T2,in, patient T2,en)
    ("diffuse_en", T2IN_C, T2EN_C * 1.0984),
    ("diffuse_both", T2IN_C * 1.0984, T2EN_C * 1.0984),
    ("fig3", T2IN_C, 70.0),
    ("focal_small", T2IN_C * 1.02, T2EN_C * 2.1),
    ("focal_large", T2IN_C * 0.77, T2EN_C * 3.4),
]


def nu_app(te, t2in, t2en):                                 # same formula as nu_app in residual_s15.py
    return float(C.apparent(NU, 0.0, C.weights(te, np.array([t2in, t2en, T2ISO])))[0])


def lg(p):
    return math.log(p / (1 - p))


def r3(x):                                                  # same rounding as residual_s15.py
    return float(Decimal(repr(x)).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP))


def pr3(x):                                                 # printed form in the manuscript (Unicode minus)
    return f"{r3(x):.3f}".replace("-", "−")


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


s15 = json.loads(S15_JSON.read_text(encoding="utf-8"))
s74 = json.loads(S74_JSON.read_text(encoding="utf-8"))

# (1) are the inputs the same as the record of the table
i15 = s15["inputs"]
assert (i15["nu"], i15["T2_in_controls_ms"], i15["T2_en_controls_ms"], i15["TE_B_ms"], i15["TE_A_ms"]) == (NU, T2IN_C, T2EN_C, 70.0, [80.0, 90.0, 100.0]), i15
assert list(s15["rows"]) == [k for k, _, _ in ROWS], list(s15["rows"])
for key, t2in_p, t2en_p in ROWS:
    r = s15["rows"][key]
    assert (r["T2_in_patients_ms"], r["T2_en_patients_ms"]) == (t2in_p, t2en_p), (key, r["T2_in_patients_ms"], r["T2_en_patients_ms"])
assert sha(CFG) == s15["config_py_sha256"], "config.py 가 표 S15 를 낸 뒤 바뀜"

out = {"what": "Supplementary Table S15, column 'within one site' (direction §5d): patient − control apparent NDI at one TE, both groups at that TE "
               "(illustration of Section 2.4; Eqs. (3)-(4); derived; no data)",
       "model": "decomp/stage1/config.py weights (TE term only, no TR term) and apparent (free water 0), the code path of residual_s15.py (Table S15) and Fig. 3",
       "inputs": {"nu_both_groups": NU, "T2_in_controls_ms": T2IN_C, "T2_en_controls_ms": T2EN_C, "T2_iso_ms_unused_fiso0": T2ISO,
                  "TE_ms": list(TES), "l2_factor_code": L2, "l2_pct_code": (L2 - 1) * 100,
                  "brandhofe_ratio_113.03/102.9": 113.03 / 102.9, "brandhofe_pct": (113.03 / 102.9 - 1) * 100,
                  "brandhofe_pct_printed": f"{(113.03 / 102.9 - 1) * 100:.2f}", "l2_pct_code_printed": f"{(L2 - 1) * 100:.2f}"},
       "read_only_sources": {"config_py": {"path": "decomp/stage1/config.py (repository)", "sha256": sha(CFG)},
                             "residual_s15_json": {"path": "x15/ledger/derived/residual_s15.json", "sha256": sha(S15_JSON)},
                             "residual_s15_py": {"path": "x15/ledger/derived/residual_s15.py", "sha256": sha(S15_PY)},
                             "s74_clinical_json": {"path": "x15/ledger/derived/s74_clinical.json", "sha256": sha(S74_JSON)},
                             "s74_clinical_py": {"path": "x15/ledger/derived/s74_clinical.py", "sha256": sha(S74_PY)}},
       "runtime": {"python": platform.python_version(), "numpy": np.__version__},
       "rows": {}}

dR_c = 1 / T2EN_C - 1 / T2IN_C
for key, t2in_p, t2en_p in ROWS:
    dR_p = 1 / t2en_p - 1 / t2in_p
    row = {"label": s15["rows"][key]["label"], "source": s15["rows"][key]["source"], "T2_in_patients_ms": t2in_p, "T2_en_patients_ms": t2en_p, "within": {}}
    for te in TES:
        c, p = nu_app(te, T2IN_C, T2EN_C), nu_app(te, t2in_p, t2en_p)
        raw, logit, logit_eq4 = p - c, lg(p) - lg(c), te * (dR_p - dR_c)
        assert abs(logit - logit_eq4) < 1e-9, (key, te, logit, logit_eq4)                       # (3)
        row["within"][f"te{int(te)}"] = {"TE_ms": te, "control_app_ndi": c, "patient_app_ndi": p, "raw": raw, "raw_printed": pr3(raw),
                                         "logit": logit, "logit_eq4": logit_eq4, "logit_printed": pr3(logit)}
    w = row["within"]
    for te_a in (80, 90, 100):                                                                   # (2)
        ref = s15["rows"][key][f"dte{te_a - 70}"]
        d_raw, d_log = w[f"te{te_a}"]["raw"] - w["te70"]["raw"], w[f"te{te_a}"]["logit"] - w["te70"]["logit"]
        assert abs(d_raw - ref["raw"]) < 1e-12, (key, te_a, d_raw, ref["raw"])
        assert abs(d_log - ref["logit"]) < 1e-9, (key, te_a, d_log, ref["logit"])
    grid = [70.0 + 0.5 * k for k in range(61)]                                                   # TE 70–100 ms, step 0.5 ms (check of the range sentence)
    g = [nu_app(te, t2in_p, t2en_p) - nu_app(te, T2IN_C, T2EN_C) for te in grid]
    row["grid_te70_100_step0.5"] = {"abs_min": min(abs(v) for v in g), "abs_max": max(abs(v) for v in g),
                                     "te_at_abs_min": grid[min(range(61), key=lambda k: abs(g[k]))], "te_at_abs_max": grid[max(range(61), key=lambda k: abs(g[k]))],
                                     "abs_increases_with_te": all(abs(g[k + 1]) > abs(g[k]) for k in range(60)), "all_negative": all(v < 0 for v in g)}
    out["rows"][key] = row

# (4) Fig. 3 row = 0.086 · 0.120 of Section 2.4
f3 = out["rows"]["fig3"]["within"]
chk4 = {"te70": f3["te70"]["raw"], "ledger_N649": 0.086, "te100": f3["te100"]["raw"], "ledger_N650": 0.120,
        "equals_residual_s15_fig3_check": abs(-f3["te70"]["raw"] - s15["fig3_check"]["gap_TE70"]) < 1e-12 and abs(-f3["te100"]["raw"] - s15["fig3_check"]["gap_TE100"]) < 1e-12}
chk4["ok"] = r3(-f3["te70"]["raw"]) == 0.086 and r3(-f3["te100"]["raw"]) == 0.120 and chk4["equals_residual_s15_fig3_check"]
out["fig3_check"] = chk4
if not chk4["ok"]:
    sys.exit(f"그림 3 대조 실패: {chk4}")

# (5) check against s74_clinical.json — 113.03/102.9 on the same code path
K = 113.03 / 102.9
i74 = s74["inputs"]
assert (i74["nu"], i74["t2_in"], i74["t2_en"], i74["ratio"]) == (NU, T2IN_C, T2EN_C, K), i74
cmp74 = {}
for te in (70, 100):
    c = nu_app(float(te), T2IN_C, T2EN_C)
    both_k, extra_k = nu_app(float(te), T2IN_C * K, T2EN_C * K) - c, nu_app(float(te), T2IN_C, T2EN_C * K) - c
    ref = s74["by_te"][str(te)]
    same_path = abs(c - ref["control"]) < 1e-12 and abs(both_k - ref["both"]) < 1e-12 and abs(extra_k - ref["extra_only"]) < 1e-12
    assert same_path, (te, c, both_k, extra_k, ref)
    t_both, t_extra = out["rows"]["diffuse_both"]["within"][f"te{te}"]["raw"], out["rows"]["diffuse_en"]["within"][f"te{te}"]["raw"]
    cmp74[str(te)] = {"both": {"table_path_l2_1.0984": t_both, "s74_clinical_113.03/102.9": ref["both"], "difference": t_both - ref["both"],
                               "printed_equal": pr3(t_both) == pr3(ref["both"])},
                      "extra_only": {"table_path_l2_1.0984": t_extra, "s74_clinical_113.03/102.9": ref["extra_only"], "difference": t_extra - ref["extra_only"],
                                     "printed_equal": pr3(t_extra) == pr3(ref["extra_only"])},
                      "s74_reproduced_by_this_code_path_with_113.03/102.9": same_path}
out["s74_clinical_comparison"] = cmp74


# range of the sentence (size · four TEs · must equal the 0.5 ms grid)
def rng(keys):
    v = [abs(out["rows"][k]["within"][f"te{int(t)}"]["raw"]) for k in keys for t in TES]
    gmin = min(out["rows"][k]["grid_te70_100_step0.5"]["abs_min"] for k in keys)
    gmax = max(out["rows"][k]["grid_te70_100_step0.5"]["abs_max"] for k in keys)
    assert abs(min(v) - gmin) < 1e-15 and abs(max(v) - gmax) < 1e-15                              # the ends of the range are TE 70 · 100 ms
    return {"abs_min": min(v), "abs_max": max(v), "abs_min_printed": pr3(min(v)), "abs_max_printed": pr3(max(v))}


out["sentence_2_4"] = {"both": rng(["diffuse_both"]), "extra_only": rng(["diffuse_en"])}
out["sentence_5_2"] = {"diffuse": rng(["diffuse_both", "diffuse_en"])}

(HERE / "s15_within_site.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
print("row            TE 70 · 80 · 90 · 100 ms: raw (logit)")
for k, v in out["rows"].items():
    print(f"{k:13s}", " · ".join(f"{v['within'][f'te{int(t)}']['raw']:+.6f} ({v['within'][f'te{int(t)}']['logit']:+.4f})" for t in TES),
          "| printed", " · ".join(v["within"][f"te{int(t)}"]["raw_printed"] for t in TES), "| |d| monotone", v["grid_te70_100_step0.5"]["abs_increases_with_te"])
print("fig3 check", chk4["ok"], "· s74 (5) path reproduced", all(x["s74_reproduced_by_this_code_path_with_113.03/102.9"] for x in cmp74.values()))
for te, x in cmp74.items():
    print(f"TE {te}: both {x['both']['table_path_l2_1.0984']:.6f} vs s74 {x['both']['s74_clinical_113.03/102.9']:.6f} (diff {x['both']['difference']:+.2e}) · "
          f"extra {x['extra_only']['table_path_l2_1.0984']:.6f} vs s74 {x['extra_only']['s74_clinical_113.03/102.9']:.6f} (diff {x['extra_only']['difference']:+.2e})")
print("2.4", out["sentence_2_4"], "\n5.2", out["sentence_5_2"])
print("pct", out["inputs"]["brandhofe_pct"], out["inputs"]["l2_pct_code"])
