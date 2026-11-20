"""Values printed in the manuscript from the results of the multi-site design simulation (decomp/design_sim); no new computation, per-cell metrics copied and rounded."""
import hashlib
import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
OUT = pathlib.Path(sys.argv[1])
S = json.loads((OUT / "summary.json").read_text(encoding="utf-8"))
C = json.loads((OUT / "checks.json").read_text(encoding="utf-8"))
if not C["all_ok"]:
    raise SystemExit("STOP: 검사 K-0 · K-a – K-e 가 모두 맞지 않음")
cells = {}
for p in sorted((OUT / "cells").glob("cell_*.json")):
    d = json.loads(p.read_text(encoding="utf-8"))
    c = d["cell"]
    cells[(c["D"], c["S"], float(c["V"]), c["H"], c["N"])] = d["metrics"]


def m(D, S_, V=0.0, H="H0", N="bal"):
    return cells[(D, S_, float(V), H, N)]


R = S["readings"]
raw = dict(
    m1_cov_e=R["R-a S2e"]["M1_coverage"], m1_cov_b=R["R-a S2b"]["M1_coverage"],
    m1_bias_e=R["R-a S2e"]["M1_bias"], m1_bias_b=R["R-a S2b"]["M1_bias"],
    m1c_cov_e=m("D1", "S2e")["M1c"]["coverage"], m2_cov_e=m("D1", "S2e")["M2"]["coverage"],
    m3_cov_e=R["R-b e"]["M3_te0_coverage_S2"], m3_cov_b=R["R-b b"]["M3_te0_coverage_S2"],
    m3_pow_e=R["R-b e"]["M3_te0_power_S1"], m3_pow_b=R["R-b b"]["M3_te0_power_S1"],
    m3_se_d1=R["R-c S1e"]["D1_se"], m3_se_d2=R["R-c S1e"]["D2_se"],
    m3_pow_d1_e=R["R-c S1e"]["D1_power"], m3_pow_d2_e=R["R-c S1e"]["D2_power"],
    m3_pow_d1_b=R["R-c S1b"]["D1_power"], m3_pow_d2_b=R["R-c S1b"]["D2_power"],
    rd_slope_reject=R["R-d"]["M3_slope_reject"],
    m5b_cov_e=m("D1", "S2e")["M5b"]["coverage"], m5a_cov_e=m("D1", "S2e")["M5a"]["coverage"],
    m1_cov_s0=m("D1", "S0")["M1"]["coverage"], m2_cov_d0_s0=m("D0", "S0")["M2"]["coverage"],
    m1_cov_e_v10=m("D1", "S2e", V=0.10)["M1"]["coverage"], m1_cov_e_unbal=m("D1", "S2e", N="unbal")["M1"]["coverage"],
)
printed = {k: (round(v, 2) if "cov" in k or "pow" in k or "reject" in k else round(v, 3)) for k, v in raw.items()}
sentences = {k: v.get("sentence") for k, v in R.items() if isinstance(v, dict) and "sentence" in v}
rec = dict(source=str(OUT.name), summary_sha256=hashlib.sha256((OUT / "summary.json").read_bytes()).hexdigest(),
           checks_all_ok=C["all_ok"], raw=raw, printed=printed, sentences=sentences)
(HERE / "s82_values.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")
print(json.dumps(rec["printed"], ensure_ascii=False))
print(json.dumps(sentences, ensure_ascii=False, indent=1))
