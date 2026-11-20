"""cause2: exploratory descriptive trends of OXF1PRI~OXF2PRI ROI units against log b2000 SNR and position; not a verdict."""
from __future__ import annotations

import datetime
import hashlib
import json
import pathlib
import re
import sys
import xml.etree.ElementTree as ET

import numpy as np

REPO = pathlib.Path("<REPO>")
DR = pathlib.Path("<WORK>/stage1_work/drafts")
sys.path.insert(0, str(REPO / "decomp" / "remeasure"))
import judge_v4 as J  # noqa: E402
import roi as ROI  # noqa: E402

PAIR = "OXF1PRI~OXF2PRI"
OUT = J.R.WORK / "explore" / "cause"
D1 = J.R.WORK / "explore" / "d1"
V4 = J.R.WORK / "onharmony" / J.V4_DIR
REREAD = DR / "x14" / "REREAD_LOG.jsonl"
SHA_ATTR = "5184137854abc69710904889503eaf8dbb111ca9736e46af13b49067ea65b90d"
BOOT, SEED = 2000, 20261103
NOW = datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def reread(path, what, expect=None) -> bytes:
    b = pathlib.Path(path).read_bytes()
    s = hashlib.sha256(b).hexdigest()
    with REREAD.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(dict(시각=datetime.datetime.now().astimezone().isoformat(timespec="seconds"), 다시읽음=str(path), sha256=s, 바이트=len(b),
                                 무엇=what, 스크립트=__file__), ensure_ascii=False) + "\n")
    if expect is not None and s != expect:
        raise RuntimeError(f"⛔ {path} sha256 {s[:16]}… ≠ {expect[:16]}…")
    return b


def rank(v):
    o = np.argsort(np.argsort(v, kind="mergesort"), kind="mergesort").astype(float)
    # Equal values get the mean rank
    _, inv, cnt = np.unique(v, return_inverse=True, return_counts=True)
    s = np.zeros(len(cnt))
    np.add.at(s, inv, o)
    return (s / cnt)[inv]


def spearman(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return float("nan")
    return float(np.corrcoef(rank(a[ok]), rank(b[ok]))[0, 1])


def ols_std(P, o):
    Z = (P - P.mean(0)) / P.std(0, ddof=1)
    X = np.c_[np.ones(len(o)), Z]
    beta, *_ = np.linalg.lstsq(X, o, rcond=None)
    pred = X @ beta
    r2 = 1 - ((o - pred) ** 2).sum() / ((o - o.mean()) ** 2).sum()
    return beta[1:] / o.std(ddof=1), float(r2)


def label_xyz(xml_path):
    root = ET.fromstring(reread(xml_path, "원인 ② — JHU 라벨 대표 화소 좌표(저장소 참조 파일 · 영상 아님)"))
    out = {}
    for lab in root.iter("label"):
        i = int(lab.get("index"))
        if i > 0:
            out[i] = np.array([90 - float(lab.get("x")), float(lab.get("y")) - 126, float(lab.get("z")) - 72])
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    # Position
    xyz = label_xyz(J.R.ATLAS_XML)
    names = ROI.label_units()
    units = ROI.unit_order()
    ucoord = {u: np.mean([xyz[i] for i, nm in enumerate(names, start=1) if nm == u], 0) for u in units}
    center = np.mean([ucoord[u] for u in units], 0)
    # ROI units · results
    d = dict(np.load(D1 / "d1_c3_units" / f"{PAIR}.npz", allow_pickle=False))
    use = ~d["main_excluded"] & ~d["flag"] & d["ok"].all(1)
    c1 = dict(np.load(OUT / "cause1_units.npz", allow_pickle=False))
    subs, rois = d["subjects"][use], d["rois"][use]
    assert (subs == c1["subjects"]).all() and (rois == c1["rois"]).all(), "① 단위와 D1 단위가 다르다"
    with np.errstate(invalid="ignore", divide="ignore"):
        Z = d["mean"][use][:, 0] / d["sd"][use][:, 0]                         # (n, 9, 3) · combination 0 = 32|64 (the only candidate of this cell)
    zN, zO = np.nanmedian(Z[:, :, 0], 1), np.nanmedian(Z[:, :, 2], 1)
    resid = c1["D"] - c1["exp_32_64"]
    outcomes = {"z_NDI": zN, "z_ODI": zO, "resid_x": resid[:, 0], "resid_ODI": resid[:, 2]}
    # b2000 SNR — roi.json «주» (main analysis) s_2000 (two sessions)
    at = json.loads(reread(V4 / "attr_table.json", "원인 ② — 속성 표(두 세션 이름)", SHA_ATTR))["pairs"][PAIR]
    root = J.R.WORK / "onharmony"
    s2 = {}
    for sub, x in at.items():
        for k, ses in enumerate(x["sessions"]):
            r = json.loads(reread(root / f"{sub}_{ses}" / "roi.json", "원인 ② — roi.json «주» s_2000(판정 때 읽은 것 · OXF1PRI~OXF2PRI 만)"))["주"]
            s2[(sub, k)] = {u: (r.get(u) or {}).get("s_2000") for u in units}
    lsA = np.array([np.log(c1["snrA"][i]) + (s2[(s, 0)][u] if s2[(s, 0)][u] is not None else np.nan) for i, (s, u) in enumerate(zip(subs, rois))])
    lsB = np.array([np.log(c1["snrB"][i]) + (s2[(s, 1)][u] if s2[(s, 1)][u] is not None else np.nan) for i, (s, u) in enumerate(zip(subs, rois))])
    preds = {"log_b2000_snr": 0.5 * (lsA + lsB), "z_mm": np.array([ucoord[u][2] for u in rois]),
             "dist_mm": np.array([np.linalg.norm(ucoord[u] - center) for u in rois])}
    P = np.c_[preds["log_b2000_snr"], preds["z_mm"], preds["dist_mm"]]
    ok = np.isfinite(P).all(1)
    usub = np.unique(subs)
    res = dict(made=NOW, script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(), pair=PAIR, n_units=int(len(subs)), n_used=int(ok.sum()),
               n_subjects=int(len(usub)), center_mm=[float(v) for v in center],
               predictor_spearman={f"{a}~{b}": spearman(preds[a][ok], preds[b][ok]) for a, b in (("log_b2000_snr", "z_mm"), ("log_b2000_snr", "dist_mm"), ("z_mm", "dist_mm"))},
               predictor_summary={k: dict(median=float(np.nanmedian(v)), q25=float(np.nanpercentile(v, 25)), q75=float(np.nanpercentile(v, 75))) for k, v in preds.items()},
               outcomes={})
    idx = {s: np.flatnonzero(ok & (subs == s)) for s in usub}
    for on, y in outcomes.items():
        sgn = float(np.sign(np.nanmean(y[ok])))
        o = y * sgn
        okk = ok & np.isfinite(o)
        rho = {k: spearman(v[okk], o[okk]) for k, v in preds.items()}
        beta, r2 = ols_std(P[okk], o[okk])
        bs = {k: [] for k in preds}
        bb, br2 = [], []
        for _ in range(BOOT):
            ix = np.concatenate([idx[s] for s in rng.choice(usub, len(usub), replace=True)])
            ix = ix[np.isfinite(o[ix])]
            for k, v in preds.items():
                bs[k].append(spearman(v[ix], o[ix]))
            b_, r_ = ols_std(P[ix], o[ix])
            bb.append(b_)
            br2.append(r_)
        bb = np.array(bb)
        ci = {k: [float(np.nanquantile(v, 0.025)), float(np.nanquantile(v, 0.975))] for k, v in bs.items()}
        read = {
            "(가) 신호 약할수록 치우침 큼": "예측과 맞음" if ci["log_b2000_snr"][1] < 0 else ("반대 방향" if ci["log_b2000_snr"][0] > 0 else "가를 수 없음"),
            "(나) 위아래 위치를 따라감": "예측과 맞음" if (ci["z_mm"][0] > 0 or ci["z_mm"][1] < 0) else "가를 수 없음",
            "(다) 깊이(중심 거리)를 따라감": "예측과 맞음" if (ci["dist_mm"][0] > 0 or ci["dist_mm"][1] < 0) else "가를 수 없음",
        }
        res["outcomes"][on] = dict(sign_aligned=sgn, mean=float(np.nanmean(y[okk])), n=int(okk.sum()), spearman=rho, spearman_boot95=ci,
                                   ols_std_beta=dict(zip(preds, [float(x) for x in beta])),
                                   ols_std_beta_boot95={k: [float(np.quantile(bb[:, i], 0.025)), float(np.quantile(bb[:, i], 0.975))] for i, k in enumerate(preds)},
                                   ols_r2=r2, ols_r2_boot95=[float(np.quantile(br2, 0.025)), float(np.quantile(br2, 0.975))], reading=read)
    (OUT / "cause2_trend.json").write_text(json.dumps(res, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    L = [f"# 원인 ② 신호 세기 · 위치 경향 — {PAIR}(탐색 · 기술 · 판정 아님 · {NOW})", "",
         f"- 단위 {res['n_used']}/{res['n_units']} · 피험자 {res['n_subjects']}(묶음 · 독립 아님) · 중심(mm) {[round(v, 1) for v in center]}",
         "- 예측자끼리 ρ: " + " · ".join(f"{k} {v:+.2f}" for k, v in res["predictor_spearman"].items()), "",
         "| 결과(치우침 쪽 정렬) | 평균 | ρ log b2000 SNR [95 %] | ρ z(mm) [95 %] | ρ 중심 거리 [95 %] | 함께 넣은 β(SNR · z · 거리) | R² | 읽기 (가) · (나) · (다) |",
         "|---|---|---|---|---|---|---|---|"]
    for on, v in res["outcomes"].items():
        s, c, b = v["spearman"], v["spearman_boot95"], v["ols_std_beta"]
        L.append(f"| {on} | {v['mean']:+.3f} | {s['log_b2000_snr']:+.2f} [{c['log_b2000_snr'][0]:+.2f}, {c['log_b2000_snr'][1]:+.2f}] | "
                 f"{s['z_mm']:+.2f} [{c['z_mm'][0]:+.2f}, {c['z_mm'][1]:+.2f}] | {s['dist_mm']:+.2f} [{c['dist_mm'][0]:+.2f}, {c['dist_mm'][1]:+.2f}] | "
                 f"{b['log_b2000_snr']:+.2f} · {b['z_mm']:+.2f} · {b['dist_mm']:+.2f} | {v['ols_r2']:.2f} | " + " · ".join(v["reading"].values()) + " |")
    (OUT / "cause2_trend.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
