"""Builds exploratory, descriptive d5 tables from d5_*.json and vr_gate.json into explore/d5/d5.md, without hand-copied numbers; not a verdict."""
import json, pathlib
import numpy as np
O = pathlib.Path("<WORK>/remeasure/explore/d5")
J = lambda n: json.loads((O / n).read_text(encoding="utf-8")) if (O / n).exists() else None  # noqa: E731
g, d1, d2, d3, d4 = J("vr_gate.json"), J("d5_1_device.json"), J("d5_2_transfer.json"), J("d5_3_outliers.json"), J("d5_4_link.json")
L = ["# D5 요약(탐색 · 기술 · 판정 아님 · 계획서 §17 · §17-1 · 스크립트 산출)", ""]
L += [f"- V-R(지도 재추출이 저장된 roi.json «주» 와 같음): {'통과' if g['passed'] else '실패'} · 세션 {g['n']} · 최대 차 {g['max_abs_diff']}", ""]
if d1:
    L += ["## ① 기기 판별(피험자 하나 빼기 · 수축 LDA · 균형 정확도 · 기기 7)", "", "| 특징 묶음 | 세션 | 특징 수 | 균형 정확도 | 우연 기준 95 % | 우연보다 높음 |", "|---|---|---|---|---|---|"]
    for nm, r in d1["variants"]["main"].items():
        c = r.get("chance_95")
        L.append(f"| {nm} | {r['n_sessions']} | {r['n_features']} | {r['balanced_accuracy']:.3f} | {'—' if c is None else f'{c:.3f}'} | {r.get('above_chance')} |")
    n2 = d1["not2ing_89_vs_98_T"]
    L += ["", f"- 부: NOT2ING-89 대 98(같은 스캐너 · TE 만 다름 · T) 균형 정확도 {n2['balanced_accuracy']:.3f} · 우연 기준 95 % {n2['chance_95']:.3f} · n {n2['n']}(89 = {n2['n_89']})",
          "- 민감도(주 칸 제외 세션 넣음): " + " · ".join(f"{k} {v['balanced_accuracy']:.3f}" for k, v in d1["variants"]["with_main_excluded"].items()), ""]
    for nm in ("T", "E"):
        r = d1["variants"]["main"][nm]
        L += [f"혼동 행렬({nm} · 행 = 참 · 열 = 예측 · 순서 {' · '.join(r['labels'])}):", "", "```"] + [" ".join(f"{x:3d}" for x in row) for row in r["confusion"]] + ["```", ""]
if d2:
    L += ["## ② 기기 사이 대응 예측(R = RMS(오차 ÷ OXF1PRI 재검사 차이 SD) · 피험자 하나 빼기 · M1 = ROI 마다 이동)", "",
          "| 쌍 | 피험자 | SNR A · B 중앙 | 기대 R | x: M0 → M1 [95 %] 읽기 | y: M0 → M1 [95 %] 읽기 | ODI: M0 → M1 [95 %] 읽기 | M1 R_between x · y · ODI |", "|---|---|---|---|---|---|---|---|"]
    for p, v in d2["pairs"].items():
        cells = []
        for c in ("x", "y", "ODI"):
            m0, m1 = v["coords"][c]["M0"], v["coords"][c]["M1"]
            cells.append(f"{m0['R']:.2f} → {m1['R']:.2f} [{m1['ci95'][0]:.2f}, {m1['ci95'][1]:.2f}] {m1['reading']}")
        rb = " · ".join(f"{v['coords'][c]['M1']['R_between']:.2f}" for c in ("x", "y", "ODI"))
        L.append(f"| {p} | {v['n_subjects']} | {v['snr_A_median']:.1f} · {v['snr_B_median']:.1f} | {v['expected_R_M1']:.2f} | " + " | ".join(cells) + f" | {rb} |")
    rd = [v["coords"][c]["M1"]["reading"] for v in d2["pairs"].values() for c in ("x", "y", "ODI")]
    L += ["", "- 읽기 수(쌍 × 좌표 " + str(len(rd)) + "): " + " · ".join(f"{k} {rd.count(k)}" for k in sorted(set(rd))),
          "- 중앙 R(23 쌍): " + " · ".join(f"{c} M0 {np.median([v['coords'][c]['M0']['R'] for v in d2['pairs'].values()]):.2f} · M1 {np.median([v['coords'][c]['M1']['R'] for v in d2['pairs'].values()]):.2f} · "
                                         f"M1g {np.median([v['coords'][c]['M1g']['R'] for v in d2['pairs'].values()]):.2f} · M2 {np.median([v['coords'][c]['M2']['R'] for v in d2['pairs'].values() if v['coords'][c]['M2']]):.2f}" for c in ("x", "y", "ODI")),
          f"- 재검사 차이 SD(OXF1PRI 반복 · ROI 중앙): " + " · ".join(f"{k} {v:.4f}" for k, v in d2["sigma_d_median"].items()),
          "- 한계: 재검사 척도는 OXF1PRI 한 기기(SNR 15 – 20)의 것이다 — SNR 이 낮은 기기의 쌍은 대응이 완벽해도 R 이 클 수 있다.", ""]
if d3:
    L += ["## ③ 장비 안 특이점 · 전형 분류", ""]
    for nm, v in d3["session"].items():
        L.append(f"- 세션 {nm}(p {v['p']}): 보정 문턱 {v['cutoff']:.1f}(χ² 0.99 = {v['chi2_99']:.1f}) · 특이 {v['n_flag']}(정규 기대 {v['expected_false_flags']:.1f}) · 안정 특이 {v['n_stable']} "
                 + (f"— {', '.join(s['session'] + f' (d² {s['d2']:.1f} · 비율 {s['freq']:.2f})' for s in v['stable_sessions'])}" if v["stable_sessions"] else ""))
        for k, x in v["main_excluded_check"].items():
            L.append(f"  - 점검 · 주 칸 제외 세션 {k}: d² {x['d2']:.1f} · 특이 {x['flag']}")
    u = d3.get("unit")
    if u:
        L += [f"- 단위(세션 × ROI · {u['n_rows']} 행 · p 6): 보정 문턱 {u['cutoff']:.1f} · 특이 {u['n_flag']}(정규 기대 {u['expected_false_flags']:.0f}) · 안정 특이 {u['n_stable']}",
              "  - **해석 주의(결과를 본 뒤 · A):** 정규 기대의 수십 배다. 특징 분포가 정규와 멀다 — 유효 비율이 정확히 1.0 인 행이 많아(MAD 0 → 대체 척도) 1 보다 조금만 작아도 크게 벗어나고, FWF 90 백분위에 0 인 행이 있다. "
              "그래서 단위 «안정 특이» 는 «드문 이상» 이 아니라 «분포 꼬리 · 경계값» 을 많이 담는다. 미리 정한 절차의 결과로 그대로 두고, ④ 에서 U-특이를 쓸 때 이 단서를 붙인다."]
    t = d3.get("typology", {})
    if t:
        L += [f"- 전형 분류(가) 장비 사이: GMM k {t['between']['k']} · 기기와 ARI {t['between']['ari_vs_device']:.3f}",
              f"- 전형 분류(나) 장비 안: GMM k {t['within']['k']} · 크기 {t['within']['sizes']} · 부트스트랩 ARI 중앙 {t['within']['stability_ari_median']} · 기기와 ARI {t['within']['ari_vs_device']:.3f}", ""]
if d4:
    L += ["## ④ D1 과 잇기", "", "| 표지 · 결과 | 칸 | 칸 사이 중앙 차 [95 %] | 단위 가중 차 [95 %] |", "|---|---|---|---|"]
    for k, v in d4["aggregate"].items():
        f = lambda x, ci: "—" if x is None else f"{x:+.3f} [{ci[0]:+.3f}, {ci[1]:+.3f}]"  # noqa: E731
        L.append(f"| {k} | {v['n_cells']} | {f(v['median_diff'], v['median_ci95'])} | {f(v['weighted_diff'], v['weighted_ci95'])} |")
    L += ["", "- 쌍 단위 Spearman(② M1 R 대 D1 ODI 칸): " + json.dumps(d4["pair_level_spearman_R_M1_vs_D1_ODI"], ensure_ascii=False)]
(O / "d5.md").write_text("\n".join(L) + "\n", encoding="utf-8")
print("\n".join(L))
