"""Compares the original b2000 result (with eddy-current term) and the plan_rigid6 result (rigid 6 only), descriptive only, without changing the verdict."""
import json, pathlib, subprocess, sys
ORIG = pathlib.Path("<WORK>/remeasure/explore/b2000")
NEW = ORIG / "plan_rigid6"
RR = pathlib.Path(__file__).resolve().parents[1] / "ledger" / "reread_log.py"
subprocess.run([sys.executable, "-I", "-B", str(RR), str(ORIG / "b2k_result.json"), "§27-11 계획서판 비교 — 원판 결과"], check=True)
A = json.loads((ORIG / "b2k_result.json").read_text(encoding="utf-8"))
B = json.loads((NEW / "b2k_result.json").read_text(encoding="utf-8"))
out = {"devices": {}, "contrasts": {}, "decision": {}}
for u in ("OXF1PRI", "OXF2PRI", "OXF3TRI"):
    a, b = A["devices"][u], B["devices"][u]
    out["devices"][u] = {k: {"orig": a.get(k), "plan": b.get(k)} for k in ("R_med", "n_value", "n_sigma_q_drop")}
    out["devices"][u]["L_eff"] = {how: {"orig": a["L_eff"][how], "plan": b["L_eff"].get(how)} for how in a["L_eff"]}
    out["devices"][u]["separable"] = {"orig": a.get("separable"), "plan": b.get("separable")}
    out["devices"][u]["spread"] = {"orig": a.get("spread"), "plan": b.get("spread")}
for k in A["M_prime"]["contrasts"]:
    out["contrasts"][k] = {"orig": A["M_prime"]["contrasts"][k], "plan": B["M_prime"]["contrasts"].get(k)}
out["decision"] = {"orig": A["M_prime"]["decision"], "plan": B["M_prime"]["decision"]}
out["literal_27"] = {"orig": A.get("literal_27", {}).get("decision"), "plan": B.get("literal_27", {}).get("decision")}
(NEW / "compare_orig_vs_plan.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
print(json.dumps(out, ensure_ascii=False, indent=1)[:6000])
