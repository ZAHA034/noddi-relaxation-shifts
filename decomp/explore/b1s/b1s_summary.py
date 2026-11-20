"""Summarizes b1s results (descriptive only, not a verdict): angles and ODI TE slopes, after validation gates on synthetic lines and noise-free angles."""
import json, math, pathlib, sys
import numpy as np
REPO = pathlib.Path("<REPO>")
sys.path.insert(0, str(REPO / "decomp" / "remeasure"))
import judge_wong as JW  # noqa: E402
O = pathlib.Path("<WORK>/remeasure/explore/b1s")
z = np.load(O / "b1s_sim.npz", allow_pickle=False)
res, tes, meta = z["res"], z["tes"].astype(float), [json.loads(m) for m in z["meta"]]
assert np.array_equal(tes, np.array(JW.TES, float))
ALG = {k: math.degrees(math.atan((1 - t2in / 1000.0) / (t2in / t2en - 1))) for k, (t2in, t2en) in
       json.loads((O / "b1s_sim_meta.json").read_text(encoding="utf-8"))["t2sets"].items()}
theta = lambda X, Y: math.degrees(math.atan2(JW.fit_unit(tes, Y)["slope"], JW.fit_unit(tes, X)["slope"]))  # noqa: E731
chk = {}
for th in (25.0, 58.0, 72.0):
    X = 0.3 + math.cos(math.radians(th)) * 0.004 * (tes - 97); Y = -1.2 + math.sin(math.radians(th)) * 0.004 * (tes - 97)
    chk[str(th)] = abs(theta(X, Y) - th)
if max(chk.values()) > 1e-9:
    sys.exit(f"⛔ 검증 (1) B1 코드 {chk}")
rows = []
for i, m in enumerate(meta):
    X, Y, Q = res[i, :, 0], res[i, :, 1], res[i, :, 2]
    if not np.isfinite(res[i]).all():
        continue
    rows.append(dict(m, theta=theta(X, Y), odi_slope_10ms=10 * JW.fit_unit(tes, Q)["slope"]))
out = dict(alg_theta=ALG, b1_code_check=chk, groups={})
for cond in dict.fromkeys(r["cond"] for r in rows):
    for t2 in ALG:
        v = [r for r in rows if r["cond"] == cond and r["t2set"] == t2]
        th = np.array([r["theta"] for r in v]); od = np.array([r["odi_slope_10ms"] for r in v])
        byf = {str(f): float(np.median([r["theta"] for r in v if r["tissue"][2] == f])) for f in (0.0, 0.02, 0.05, 0.10)}
        out["groups"][f"{cond}|{t2}"] = dict(n=len(v), theta_median=float(np.median(th)), theta_q25=float(np.quantile(th, 0.25)), theta_q75=float(np.quantile(th, 0.75)),
                                            theta_median_by_fiso=byf, odi_slope_median=float(np.median(od)))
nf = {t2: out["groups"][f"noisefree|{t2}"]["theta_median"] for t2 in ALG}
out["check_noisefree"] = {t2: dict(theta=nf[t2], algebra=ALG[t2], diff=nf[t2] - ALG[t2], within_3deg=abs(nf[t2] - ALG[t2]) <= 3.0) for t2 in ALG}
if not all(v["within_3deg"] for v in out["check_noisefree"].values()):
    (O / "b1s_summary.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    sys.exit(f"⛔ 검증 (2) 잡음 0 의 각도가 대수 ±3° 밖 — 읽기 없이 멈춤: {json.dumps(out['check_noisefree'], ensure_ascii=False)}")
out["reading"] = {L: ("(나)로 설명할 수 있음" if out["groups"][f"snr×1.0|L{L}|lit_82_52"]["theta_median"] <= 38.0 else "(나)로는 설명할 수 없음") for L in (4, 8)}
(O / "b1s_summary.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
L = ["| 조건 | " + " | ".join(f"{t2} (대수 {ALG[t2]:.1f}°)" for t2 in ALG) + " |", "|---|" + "---|" * len(ALG)]
for cond in dict.fromkeys(r["cond"] for r in rows):
    L.append(f"| {cond} | " + " | ".join(f"{out['groups'][f'{cond}|{t2}']['theta_median']:.1f}° [{out['groups'][f'{cond}|{t2}']['theta_q25']:.1f}, {out['groups'][f'{cond}|{t2}']['theta_q75']:.1f}]" for t2 in ALG) + " |")
(O / "b1s_summary.md").write_text("\n".join(L) + "\n", encoding="utf-8")
print("\n".join(L)); print("읽기", out["reading"]); print("ODI 기울기(lit · ×1 · L4)", out["groups"]["snr×1.0|L4|lit_82_52"]["odi_slope_median"], "· FISO 별 θ(lit · ×1 · L4)", out["groups"]["snr×1.0|L4|lit_82_52"]["theta_median_by_fiso"])
