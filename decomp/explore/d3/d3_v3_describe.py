"""Counts boundary solutions and fit mismatch sizes from stored d3 v3 outputs, descriptive and defined after seeing the results, writing JSON and Markdown."""
import json, pathlib
import numpy as np
OUT = pathlib.Path("<WORK>/remeasure/explore/d3")
MAIN = "dpar1.7|t2iso1000"
N_DATA = 27 * 8 * 2  # ROI 27 · TE 8 (b1s_sim.TES) · two shells (s_low · s_2000)
LO, HI, TOL = 0.5, 6.0, 1e-3
rows = [json.loads(x) for x in (OUT / "d3_v3_groups.jsonl").read_text(encoding="utf-8").splitlines() if x.strip()]
CK = json.loads((OUT / "d3_v3_check.json").read_text(encoding="utf-8"))
rec = dict(note="결과를 본 뒤 · 계산 전에 정하지 않은 기술 · 넘어가는 조건과 읽기를 바꾸지 않음", main=MAIN, n_data=N_DATA,
           noise_free_VN={str(x["r_true"]): dict(rhat=x["rhat"], rss_min=x["rss_min"]) for x in CK["VN"]}, cells={})
for rt in (1.6, 2.64):
    for sc in (1.0, 2.0):
        g = [x for x in rows if x["r_true"] == rt and x["scale"] == sc and "variants" in x]
        m = np.array([x["variants"][MAIN]["rhat"] for x in g])
        allv = np.array([v["rhat"] for x in g for v in x["variants"].values()])
        inner = m[(m > LO + TOL) & (m < HI - TOL)]
        rss = np.array([x["variants"][MAIN]["rss_min"] for x in g])
        # Difference between the value at the coarse point (chain start) and the refit value — number of variant fits with min(coarse-point RSS) − RSS_min > 3.84
        jump = [min(v["rss_coarse"]) - v["rss_min"] for x in g for v in x["variants"].values()]
        rec["cells"][f"r{rt}|snr×{sc}"] = dict(
            coarse_above_min_gt_3_84=f"{sum(j > 3.84 for j in jump)}/{len(jump)}", coarse_above_min_max=float(max(jump)),
            n=len(g),
            main_at_upper=int((m >= HI - TOL).sum()), main_at_lower=int((m <= LO + TOL).sum()), main_interior=int(inner.size),
            main_interior_median=float(np.median(inner)) if inner.size else None,
            all_variants_at_bound=f"{int(((allv >= HI - TOL) | (allv <= LO + TOL)).sum())}/{allv.size}",
            hull_reaches_lower=int(sum(x["hull"][0] <= LO + TOL for x in g)), hull_reaches_upper=int(sum(x["hull"][1] >= HI - TOL for x in g)),
            main_rss_min_median=float(np.median(rss)), main_rss_per_datum_median=float(np.median(rss) / N_DATA))
(OUT / "d3_v3_describe.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")
f3 = lambda x: "—" if x is None else f"{x:.3f}"  # noqa: E731
L = ["# D3-v3 기술(결과를 본 뒤 · 계산 전에 정하지 않은 기술 · 넘어가는 조건과 읽기를 바꾸지 않음)", "",
     f"- 주 변형 {MAIN} · 묶음마다 자료 {N_DATA} 개(ROI 27 · TE 8 · 셸 둘) · RSS 는 잔차 ÷ σ 의 제곱합.",
     "- 잡음 0(V-N) 의 RSS_min: " + " · ".join(f"참값 {k} → {v['rss_min']:.3f}(r̂ {v['rhat']:.4f})" for k, v in rec["noise_free_VN"].items()), "",
     "| 참값 r · SNR 배수 | 주 변형 r̂: 위 끝 6.0 · 아래 끝 0.5 · 안쪽(중앙) | 여섯 변형 가운데 끝에 닿은 r̂ | 합 구간이 0.5 에 닿음 · 6.0 에 닿음 | 주 변형 RSS_min 중앙(자료 하나당) | 거친 점 최솟값 − RSS_min > 3.84(최대) |",
     "|---|---|---|---|---|---|"]
for k, v in rec["cells"].items():
    L.append(f"| {k.replace('|', ' · ')} | {v['main_at_upper']} · {v['main_at_lower']} · {v['main_interior']}({f3(v['main_interior_median'])}) | {v['all_variants_at_bound']} | "
             f"{v['hull_reaches_lower']}/{v['n']} · {v['hull_reaches_upper']}/{v['n']} | {v['main_rss_min_median']:.0f}({v['main_rss_per_datum_median']:.1f}) | "
             f"{v['coarse_above_min_gt_3_84']}({v['coarse_above_min_max']:.1f}) |")
(OUT / "d3_v3_describe.md").write_text("\n".join(L) + "\n", encoding="utf-8")
print("\n".join(L))
