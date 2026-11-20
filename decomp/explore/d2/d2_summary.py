"""Summarizes d2 descriptively: ratio of observed to simulated difference SDs per ROI unit and coordinate, with a contour over SNR and N_eff; not a verdict."""
import glob, json, math, pathlib
import numpy as np
O = pathlib.Path("<WORK>/remeasure/explore/d2")
obs = json.loads((O / "d2_observed.json").read_text(encoding="utf-8"))
P = np.load(O / "d2_plan.npz", allow_pickle=False)
ver = json.loads((O / "d2_sim_verify.json").read_text(encoding="utf-8"))
NEFFS, SC = ver["neffs"], ver["scen_snr"]
files = sorted(glob.glob(str(O / "d2_sim_c*.npz")))
n = len(P["theta"])
D = np.full((len(NEFFS), n, 3), np.nan)
for f in files:
    z = np.load(f); D[:, z["idx"]] = z["d"]
assert np.isfinite(D).all(), "덩어리가 모자라다"
units = [(sub, u) for sub, d in obs["subjects"].items() for u in d["units"]]
obs_sd = np.array([obs["subjects"][sub]["diff_sd"][obs["subjects"][sub]["units"].index(u)] for sub, u in units])   # (54, 3)
fornix = np.array(["Fornix" in u for _, u in units])
ci, si = P["centre_idx"], P["scen_snr"]
sim_sd = np.full((len(SC), len(NEFFS), len(units), 3), np.nan)
for s in range(len(SC)):
    for c in range(len(units)):
        m = (si == s) & (ci == c)
        sim_sd[s, :, c] = D[:, m].std(1, ddof=1)
COORD = ("x", "y", "ODI")
k50, s30 = NEFFS.index(50.0), SC.index("snr30")
out = dict(made=__import__("datetime").datetime.now().astimezone().isoformat(timespec="seconds"), n_units=len(units), coords=COORD, by_coord={})
for j, cn in enumerate(COORD):
    r = obs_sd[:, j] / sim_sd[s30, k50, :, j]
    rnf = r[~fornix]
    pooled = lambda s, k, msk=np.ones(len(units), bool): float(math.sqrt(np.mean(obs_sd[msk, j] ** 2)) / math.sqrt(np.mean(sim_sd[s, k, msk, j] ** 2)))  # noqa: E731
    contour = {SC[s]: {str(int(NEFFS[k])): round(pooled(s, k), 3) for k in range(len(NEFFS))} for s in range(len(SC))}
    cross = {}
    for s in range(len(SC)):
        v = np.array([pooled(s, k) for k in range(len(NEFFS))])     # The ratio grows as N_eff grows (the simulated SD becomes smaller)
        lg = np.log(NEFFS)
        idx = np.flatnonzero((v[:-1] - 1) * (v[1:] - 1) <= 0)
        cross[SC[s]] = (float(np.exp(lg[idx[0]] + (1 - v[idx[0]]) * (lg[idx[0] + 1] - lg[idx[0]]) / (v[idx[0] + 1] - v[idx[0]]))) if len(idx) else
                        ("모든 N_eff 에서 > 1" if v.min() > 1 else "모든 N_eff 에서 < 1"))
    out["by_coord"][cn] = dict(ratio_median=float(np.median(r)), ratio_q25=float(np.quantile(r, 0.25)), ratio_q75=float(np.quantile(r, 0.75)),
                               ratio_median_without_fornix=float(np.median(rnf)), pooled_ratio_snr30_neff50=pooled(s30, k50),
                               pooled_ratio_without_fornix=pooled(s30, k50, ~fornix), contour_pooled_ratio=contour, neff_where_ratio_1=cross)
(O / "d2_repeat.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
L = ["| 좌표 | 비 중앙 [사분위] (SNR 30 · N_eff 50) | 뇌활 뺌 | 합동 비 | 비 = 1 인 N_eff: SNR 30 · 20 · 50 · 실측 |", "|---|---|---|---|---|"]
for cn, v in out["by_coord"].items():
    cr = " · ".join(str(round(x, 1)) if isinstance(x, float) else x for x in v["neff_where_ratio_1"].values())
    L.append(f"| {cn} | {v['ratio_median']:.2f} [{v['ratio_q25']:.2f}, {v['ratio_q75']:.2f}] | {v['ratio_median_without_fornix']:.2f} | {v['pooled_ratio_snr30_neff50']:.2f} | {cr} |")
(O / "d2_repeat.md").write_text("\n".join(L) + "\n", encoding="utf-8")
print("\n".join(L))
for cn, v in out["by_coord"].items():
    print(cn, json.dumps(v["contour_pooled_ratio"], ensure_ascii=False))
