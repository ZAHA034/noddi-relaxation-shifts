"""Post hoc descriptive cause2 step: Spearman correlation of the measurement D (OXF2PRI − OXF1PRI) with signal-strength and position predictors, within and across subjects."""
from __future__ import annotations

import hashlib
import json
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import cause2_trend as C2  # noqa: E402

J, ROI = C2.J, C2.ROI


def main():
    rng = np.random.default_rng(C2.SEED)
    xyz = C2.label_xyz(J.R.ATLAS_XML)
    names, units = ROI.label_units(), ROI.unit_order()
    ucoord = {u: np.mean([xyz[i] for i, nm in enumerate(names, start=1) if nm == u], 0) for u in units}
    center = np.mean([ucoord[u] for u in units], 0)
    c1 = dict(np.load(C2.OUT / "cause1_units.npz", allow_pickle=False))
    subs, rois = c1["subjects"], c1["rois"]
    at = json.loads(C2.reread(C2.V4 / "attr_table.json", "② 뒤 기술 — 속성 표(두 세션 이름)", C2.SHA_ATTR))["pairs"][C2.PAIR]
    root = J.R.WORK / "onharmony"
    s2 = {}
    for sub, x in at.items():
        for k, ses in enumerate(x["sessions"]):
            r = json.loads(C2.reread(root / f"{sub}_{ses}" / "roi.json", "② 뒤 기술 — roi.json «주» s_2000(OXF1PRI~OXF2PRI 만)"))["주"]
            s2[(sub, k)] = {u: (r.get(u) or {}).get("s_2000") for u in units}
    g = lambda s, k, u: s2[(s, k)][u] if s2[(s, k)][u] is not None else np.nan  # noqa: E731
    lsA = np.array([np.log(c1["snrA"][i]) + g(s, 0, u) for i, (s, u) in enumerate(zip(subs, rois))])
    lsB = np.array([np.log(c1["snrB"][i]) + g(s, 1, u) for i, (s, u) in enumerate(zip(subs, rois))])
    preds = {"log_b2000_snr": 0.5 * (lsA + lsB), "z_mm": np.array([ucoord[u][2] for u in rois]),
             "dist_mm": np.array([np.linalg.norm(ucoord[u] - center) for u in rois])}
    usub = np.unique(subs)
    idx = {s: np.flatnonzero(subs == s) for s in usub}
    out = dict(script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(), label="결과를 본 뒤 · 기술 · 판정 아님", outcomes={})
    D = c1["D"]
    for on, y in (("D_x", D[:, 0]), ("D_ODI", D[:, 2])):
        ydm = y - np.array([y[idx[s]].mean() for s in subs])               # Within subject
        for tag, yy in ((on, y), (on + "_within_subject", ydm)):
            rho = {k: C2.spearman(v, yy) for k, v in preds.items()}
            bs = {k: [] for k in preds}
            for _ in range(C2.BOOT):
                ix = np.concatenate([idx[s] for s in rng.choice(usub, len(usub), replace=True)])
                for k, v in preds.items():
                    bs[k].append(C2.spearman(v[ix], yy[ix]))
            out["outcomes"][tag] = dict(mean=float(np.mean(yy)), sd=float(np.std(yy, ddof=1)), spearman=rho,
                                        spearman_boot95={k: [float(np.nanquantile(v, 0.025)), float(np.nanquantile(v, 0.975))] for k, v in bs.items()})
    (C2.OUT / "cause2_posthoc_D.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    for k, v in out["outcomes"].items():
        s, c = v["spearman"], v["spearman_boot95"]
        print(f"{k:24s} 평균 {v['mean']:+.4f} SD {v['sd']:.4f} | " + " · ".join(f"{p} {s[p]:+.2f} [{c[p][0]:+.2f}, {c[p][1]:+.2f}]" for p in preds))


if __name__ == "__main__":
    main()
