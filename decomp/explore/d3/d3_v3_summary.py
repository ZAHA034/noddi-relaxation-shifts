"""Builds d3 v3 tables from stored JSON results into explore/d3/d3_v3.md, without hand-copied numbers."""
import json, pathlib
import numpy as np
OUT = pathlib.Path("<WORK>/remeasure/explore/d3")
R = json.loads((OUT / "d3_v3.json").read_text(encoding="utf-8"))
CK = json.loads((OUT / "d3_v3_check.json").read_text(encoding="utf-8"))
f3 = lambda x: f"{x:.3f}"  # noqa: E731
L = [f"# D3-v3 요약(스크립트 산출 · {R['made']} 계산 · 스크립트 sha256 {R['script_sha256'][:16]}…)", "",
     f"- 정의: 계획서 §15(sha256 {R['plan_sha256_at_definition']}) · 변형 {len(R['variants'])} · 참값 {R['truths']} · SNR 배수 {R['scales']} · 묶음 {R['n_groups']} × 4.",
     f"- V-S {R['V_S']} · 잔차 σ b700 {R['sigma_shell']['b700']:.4f} · b2000 {R['sigma_shell']['b2000']:.4f}(차이 {R['sigma_shell']['n_diffs']}).",
     f"- 구현 점검 V-N({CK['rule_VN']}): " + " · ".join(f"참값 {x['r_true']} → r̂ {x['rhat']:.4f}" for x in CK["VN"]) + f" · V-D({CK['rule_VD']}): {CK['VD']['equal']} · 둘 다 {'통과' if CK['passed'] else '실패'}",
     f"- 예외 묶음(§15-2 (4) · «덮지 않음 · 맞게 읽히지 않음»으로 셈): {len(R['errors'])}" + (f" — gi {R['errors']}" if R["errors"] else ""), "",
     "| 참값 r · SNR 배수 | 합 구간 덮음 | 맞게 읽힌 비율 | r̂ 중앙(주 변형) | 합 구간 중앙 [lo, hi] | 주 변형 위 열림 | 예외 | 읽기 |", "|---|---|---|---|---|---|---|---|"]
for k, v in R["KAT"].items():
    L.append(f"| {k.replace('|', ' · ')} | {v['coverage_hull']:.2f} | {v['correct_reading']:.2f} | {f3(v['rhat_main_median'])} | [{f3(v['hull_median'][0])}, {f3(v['hull_median'][1])}] | {v['open_hi_main']}/{v['n'] - v['n_error']} | {v['n_error']} | " + " · ".join(f"{a} {b}" for a, b in sorted(v["readings"].items())) + " |")
g = R["gate_for_real_fit"]
L += ["", f"- 넘어가는 조건({g['rule']}): **{'충족 → 실측 적합' if g['passed'] else '미달 → 실측 적합 안 함'}**"]
if "real" in R:
    L += ["", "| 묶음 | 합 구간 [lo, hi] | 읽기 | 변형별 r̂ (lo, hi) |", "|---|---|---|---|"]
    for k, v in R["real"].items():
        if v.get("errors"):
            L.append(f"| {k} | 예외 {v['errors']} | — | — |"); continue
        L.append(f"| {k} | [{f3(v['hull'][0])}, {f3(v['hull'][1])}] | {v['reading']} | " + " · ".join(f"{n.replace('|', ' ')} {f3(x['rhat'])} ({f3(x['lo'])}, {f3(x['hi'])}){'↑열림' if x['open_hi'] else ''}{'↓열림' if x['open_lo'] else ''}" for n, x in v["variants"].items()) + " |")
(OUT / "d3_v3.md").write_text("\n".join(L) + "\n", encoding="utf-8")
print("\n".join(L))
