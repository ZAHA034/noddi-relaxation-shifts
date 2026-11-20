"""Supplementary table of the size of the residual of Proposition 4 (derived, no data): the logit residual of Eq. (7) and the raw-scale change in the signal model of decomp/stage1/config.py, with the row of Fig. 2 as failure control."""
import hashlib
import json
import math
import pathlib
import sys
from decimal import Decimal, ROUND_HALF_UP

import numpy as np

sys.dont_write_bytecode = True
REPO = pathlib.Path(sys.argv[1])
sys.path.insert(0, str(REPO / "decomp" / "stage1"))
import config as C  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
NU, T2IN_C, T2EN_C, T2ISO = 0.46, 82.5, 52.0, 1000.0
TE_B, TE_A = 70.0, (80.0, 90.0, 100.0)
ROWS = [  # key, name in the table, patient T2,in, patient T2,en, source
    ("diffuse_en", "Diffuse, extra-neurite only (T2,en + 9.84 %)", T2IN_C, T2EN_C * 1.0984, "S-M2 level L2"),
    ("diffuse_both", "Diffuse, both compartments (T2,in and T2,en + 9.84 %)", T2IN_C * 1.0984, T2EN_C * 1.0984, "S-M2 level L2"),
    ("fig3", "Illustration of Section 2.4 and Fig. 3 (T2,en 70 ms)", T2IN_C, 70.0, "Section 2.4; Fig. 3"),
    ("focal_small", "Focal lesion, small end (T2,en × 2.1, T2,in × 1.02)", T2IN_C * 1.02, T2EN_C * 2.1, "S-M2 level L3"),
    ("focal_large", "Focal lesion, large end (T2,en × 3.4, T2,in × 0.77)", T2IN_C * 0.77, T2EN_C * 3.4, "S-M2 level L3"),
]


def nu_app(te, t2in, t2en):
    return float(C.apparent(NU, 0.0, C.weights(te, np.array([t2in, t2en, T2ISO])))[0])


def lg(p):
    return math.log(p / (1 - p))


def r3(x):
    return float(Decimal(repr(x)).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP))


out = {"what": "Supplementary Table S15 (v16a): change in the patient-control difference between two sites that differ in TE (illustration of Proposition 4; derived; no data)",
       "model": "decomp/stage1/config.py weights (TE term only, no TR term) and apparent (free water 0), as for Fig. 3 (decomp/manuscript/recompute/a_values.py)",
       "inputs": {"nu": NU, "T2_in_controls_ms": T2IN_C, "T2_en_controls_ms": T2EN_C, "TE_B_ms": TE_B, "TE_A_ms": list(TE_A)},
       "config_py_sha256": hashlib.sha256((REPO / "decomp" / "stage1" / "config.py").read_bytes()).hexdigest(),
       "rows": {}}
dR_c = 1 / T2EN_C - 1 / T2IN_C
for key, label, t2in_p, t2en_p, src in ROWS:
    dR_p = 1 / t2en_p - 1 / t2in_p
    row = {"label": label, "source": src, "T2_in_patients_ms": t2in_p, "T2_en_patients_ms": t2en_p}
    for te_a in TE_A:
        dte = te_a - TE_B
        logit = dte * (dR_p - dR_c)                                               # Eq. (7)
        logit_model = (lg(nu_app(te_a, t2in_p, t2en_p)) - lg(nu_app(te_a, T2IN_C, T2EN_C))) - \
                      (lg(nu_app(TE_B, t2in_p, t2en_p)) - lg(nu_app(TE_B, T2IN_C, T2EN_C)))
        assert abs(logit - logit_model) < 1e-9, (key, dte, logit, logit_model)  # Eq. (7) is exact in the signal model
        raw = (nu_app(te_a, t2in_p, t2en_p) - nu_app(te_a, T2IN_C, T2EN_C)) - \
              (nu_app(TE_B, t2in_p, t2en_p) - nu_app(TE_B, T2IN_C, T2EN_C))
        row[f"dte{int(dte)}"] = {"dTE_ms": dte, "logit": logit, "raw": raw}
    out["rows"][key] = row

# failure control: the Fig. 3 row = the values of Section 2.4 of the manuscript
f = out["rows"]["fig3"]
gap70 = -(nu_app(70.0, T2IN_C, 70.0) - nu_app(70.0, T2IN_C, T2EN_C))
gap100 = -(nu_app(100.0, T2IN_C, 70.0) - nu_app(100.0, T2IN_C, T2EN_C))
chk = {"logit_dte30": f["dte30"]["logit"], "ledger_N106": -0.148, "gap_TE70": gap70, "ledger_N649": 0.086, "gap_TE100": gap100, "ledger_N650": 0.120,
       "raw_dte30": f["dte30"]["raw"], "raw_dte30_equals_gap70_minus_gap100": abs(f["dte30"]["raw"] - (gap70 - gap100)) < 1e-12}
chk["ok"] = (r3(chk["logit_dte30"]) == -0.148 and r3(gap70) == 0.086 and r3(gap100) == 0.120 and chk["raw_dte30_equals_gap70_minus_gap100"])
out["fig3_check"] = chk
if not chk["ok"]:
    sys.exit(f"그림 3 대조 실패: {chk}")
# range of the sentence of Section 5.2 (ΔTE 30 ms · size)
d30 = {k: abs(v["dte30"]["raw"]) for k, v in out["rows"].items()}
out["sentence_5_2"] = {"dTE_ms": 30.0,
                       "diffuse_min": min(d30["diffuse_en"], d30["diffuse_both"]), "diffuse_max": max(d30["diffuse_en"], d30["diffuse_both"]),
                       "focal_min": min(d30["focal_small"], d30["focal_large"]), "focal_max": max(d30["focal_small"], d30["focal_large"])}
(HERE / "residual_s15.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
for k, v in out["rows"].items():
    print(f"{k:13s}", " · ".join(f"ΔTE {v[d]['dTE_ms']:.0f}: {v[d]['logit']:+.4f} / {v[d]['raw']:+.4f}" for d in ("dte10", "dte20", "dte30")))
print("그림 3 대조", chk["ok"], f"logit {chk['logit_dte30']:.5f} · gap70 {gap70:.5f} · gap100 {gap100:.5f} · raw30 {chk['raw_dte30']:.5f}")
print("5.2", {k: round(v, 5) for k, v in out["sentence_5_2"].items()})
