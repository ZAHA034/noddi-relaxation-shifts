"""cause1: exploratory raw-value decomposition of OXF1PRI~OXF2PRI coordinate differences against simulator expectations from frozen test data; not a verdict."""
from __future__ import annotations

import datetime
import hashlib
import json
import pathlib
import sys

import numpy as np

REPO = pathlib.Path("<REPO>")
DR = pathlib.Path("<WORK>/stage1_work/drafts")
sys.path.insert(0, str(REPO / "decomp" / "remeasure"))
import judge_v4 as J  # noqa: E402

PAIR = "OXF1PRI~OXF2PRI"
OUT = J.R.WORK / "explore" / "cause"
D1 = J.R.WORK / "explore" / "d1"
V4 = J.R.WORK / "onharmony" / J.V4_DIR
REREAD = DR / "x14" / "REREAD_LOG.jsonl"
SHA = dict(attr_table="5184137854abc69710904889503eaf8dbb111ca9736e46af13b49067ea65b90d",
           snr_table="1714e6df9e2221df4abc59d58dc40f85c81d5e10c4406a4bcb2e36caddbd8cec")
COMBOS = {"32|64": (0, 1), "32|32": (0, 0)}                        # (coil_A_64, coil_B_64)
COMP = {"x": 0, "y": 1, "ODI": 2}
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


def design(theta, uA, uB, kind):
    f = (lambda u: u) if kind == "inv2" else (lambda u: np.sqrt(u))   # inv2 = 1/SNR² · inv1 = 1/SNR
    return np.c_[np.ones(len(theta)), theta[:, 0:3], f(uA), f(uB)]


def fit_eval(X, Y, x0):
    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    return x0 @ beta, beta


def pair_snr(at, snr):
    """(SNR_A, SNR_B) for each subject — the two sessions of the attribute table."""
    out = {}
    for sub, x in at["pairs"][PAIR].items():
        sa, sb = x["sessions"]
        out[sub] = (snr[f"{sub}_{sa}"]["snr"], snr[f"{sub}_{sb}"]["snr"])
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    at = json.loads(reread(V4 / "attr_table.json", "원인 ① — 속성 표(쌍 · 피험자의 두 세션)", SHA["attr_table"]))
    snr = json.loads(reread(V4 / "snr_table.json", "원인 ① — 세션 SNR(기대값을 낼 SNR · SNR 방향 확인)", SHA["snr_table"]))["sessions"]
    ps = pair_snr(at, snr)
    # Measurement
    d = dict(np.load(D1 / "d1_c3_units" / f"{PAIR}.npz", allow_pickle=False))
    use = ~d["main_excluded"] & ~d["flag"] & d["ok"].all(1)
    X3 = d["X3"][use]
    subs = d["subjects"][use]
    Dr = X3[:, 3:6] - X3[:, 0:3]
    sA = np.array([ps[s][0] for s in subs], float)
    sB = np.array([ps[s][1] for s in subs], float)
    usub = np.unique(subs)
    real = {}
    for c, j in COMP.items():
        v = Dr[:, j]
        bm = []
        for _ in range(BOOT):
            pick = rng.choice(usub, len(usub), replace=True)
            bm.append(float(np.mean(np.concatenate([v[subs == s] for s in pick]))))
        real[c] = dict(mean=float(v.mean()), median=float(np.median(v)), sd=float(v.std(ddof=1)), q25=float(np.percentile(v, 25)), q75=float(np.percentile(v, 75)),
                       boot95=[float(np.quantile(bm, 0.025)), float(np.quantile(bm, 0.975))], boot=bm)
    ex = lambda x: 1 / (1 + np.exp(-x))  # noqa: E731
    dndi_app = ex(X3[:, 3]) - ex(X3[:, 0])
    real["NDI_apparent"] = dict(mean=float(dndi_app.mean()), median=float(np.median(dndi_app)), sd=float(dndi_app.std(ddof=1)))
    # Simulator
    sim = np.load(J.S.WORK / "sim_test.npz", allow_pickle=False)
    N = J.S.COND_NAMES
    lay = np.array([l == "match|L0" for l in sim["layer"].tolist()])
    C, TH, XC = sim["cond"].astype(float), sim["theta"], sim["x_c3"].astype(float)
    fin = np.isfinite(XC).all(1)
    snrA_s, snrB_s = np.exp(C[:, N.index("log_snr_A")]), np.exp(C[:, N.index("log_snr_B")])
    uA_real, uB_real = float(np.median(1 / np.array([ps[s][0] for s in usub]) ** 2)), float(np.median(1 / np.array([ps[s][1] for s in usub]) ** 2))
    uA_cl = float(np.median(1 / np.maximum(np.array([ps[s][0] for s in usub]), 15.0) ** 2))
    uB_cl = float(np.median(1 / np.maximum(np.array([ps[s][1] for s in usub]), 15.0) ** 2))
    exp_, boots, unit_exp = {}, {}, {}
    for cn, (ca, cb) in COMBOS.items():
        sel = lay & fin & (C[:, N.index("coil_A_64")] == ca) & (C[:, N.index("coil_B_64")] == cb)
        Ds = XC[sel, 3:6] - XC[sel, 0:3]
        th, uA, uB = TH[sel], 1 / snrA_s[sel] ** 2, 1 / snrB_s[sel] ** 2
        X2, X1 = design(th, uA, uB, "inv2"), design(th, uA, uB, "inv1")
        p2 = np.r_[1, 0, 0, 0, uA_real, uB_real]
        p1 = np.r_[1, 0, 0, 0, np.sqrt(uA_real), np.sqrt(uB_real)]
        pcl = np.r_[1, 0, 0, 0, uA_cl, uB_cl]
        win = (snrA_s[sel] <= 30) & (snrB_s[sel] <= 30)
        Xw = np.c_[np.ones(win.sum()), th[win, 0:3]]
        e = {"n_cases": int(sel.sum()), "n_window_le30": int(win.sum())}
        bb = {c: [] for c in COMP}
        for c, j in COMP.items():
            m2, beta2 = fit_eval(X2, Ds[:, j], p2)
            e[c] = dict(primary=float(m2), s1_inv_snr=float(fit_eval(X1, Ds[:, j], p1)[0]), s2_window_intercept=float(fit_eval(Xw, Ds[win, j], np.r_[1, 0, 0, 0])[0]),
                        s3_clamped15=float(fit_eval(X2, Ds[:, j], pcl)[0]), beta=[float(x) for x in beta2])
        for _ in range(BOOT):
            ix = rng.integers(0, len(Ds), len(Ds))
            for c, j in COMP.items():
                bb[c].append(float(fit_eval(X2[ix], Ds[ix, j], p2)[0]))
        for c in COMP:
            e[c]["boot95"] = [float(np.quantile(bb[c], 0.025)), float(np.quantile(bb[c], 0.975))]
        exp_[cn], boots[cn] = e, bb
        # Expectation per ROI unit (SNR of the two sessions of the subject of that unit)
        ue = {}
        for c, j in COMP.items():
            beta = np.array(e[c]["beta"])
            ue[c] = (np.c_[np.ones(len(sA)), np.zeros((len(sA), 3)), 1 / sA ** 2, 1 / sB ** 2] @ beta)
        unit_exp[cn] = ue
    # r and reading
    rr = {}
    for c in ("x", "ODI", "y"):
        e64, e32 = exp_["32|64"][c]["primary"], exp_["32|32"][c]["primary"]
        dl = e64 - e32
        bdl = np.array(boots["32|64"][c]) - np.array(boots["32|32"][c])
        br = (np.array(real[c]["boot"]) - np.array(boots["32|32"][c])) / bdl
        dl_ci = [float(np.quantile(bdl, 0.025)), float(np.quantile(bdl, 0.975))]
        defined = not (dl_ci[0] <= 0 <= dl_ci[1])
        rr[c] = dict(delta_L=float(dl), delta_L_boot95=dl_ci, r=float((real[c]["mean"] - e32) / dl) if dl != 0 else None, r_defined=bool(defined),
                     r_boot95=[float(np.quantile(br, 0.025)), float(np.quantile(br, 0.975))],
                     dist_to_32_32_in_real_sd=float((real[c]["mean"] - e32) / real[c]["sd"]), dist_to_32_64_in_real_sd=float((real[c]["mean"] - e64) / real[c]["sd"]))
    rx, ro = rr["x"]["r"], rr["ODI"]["r"]
    if not (rr["x"]["r_defined"] and rr["ODI"]["r_defined"]):
        reading = "가를 수 없음(기대 차 Δ_L 의 구간이 0 을 포함하는 성분이 있음)"
    elif not (-0.5 <= rx <= 1.5) or not (-0.5 <= ro <= 1.5) or abs(rx - ro) > 1:
        reading = "(나) 쪽 — 두 기대 모두와 멀거나 NDI · ODI 의 r 이 1 넘게 다름"
    elif -0.5 <= rx <= 0.5 and -0.5 <= ro <= 0.5:
        reading = "(가) 쪽 — 두 성분 모두 32|32 기대에 가깝고 32|64 기대와 멂"
    elif 0.5 <= rx <= 1.5 and 0.5 <= ro <= 1.5:
        reading = "시뮬레이터 가정(L 8 대 4)과 맞음 — 두 성분 모두 32|64 기대 쪽"
    else:
        reading = "가를 수 없음"
    snr_dir = dict(OXF1PRI_median=float(np.median([ps[s][0] for s in usub])), OXF2PRI_median=float(np.median([ps[s][1] for s in usub])),
                   n_subjects=int(len(usub)), dSNR_B_minus_A_median=float(np.median([ps[s][1] - ps[s][0] for s in usub])))
    for c in COMP:
        real[c].pop("boot")
    res = dict(made=NOW, script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(), pair=PAIR, n_units=int(use.sum()), n_subjects=int(len(usub)),
               eval_point=dict(uA=uA_real, uB=uB_real, snrA_equiv=uA_real ** -0.5, snrB_equiv=uB_real ** -0.5), real=real, expected=exp_, r=rr, reading=reading,
               snr_direction=snr_dir, note="탐색 · 판정 아님 · 판단기 사후를 쓰지 않음 · 시험 자료는 순방향 값(x_c3)만")
    (OUT / "cause1_raw.json").write_text(json.dumps(res, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    np.savez(OUT / "cause1_units.npz", subjects=subs, rois=d["rois"][use], D=Dr, snrA=sA, snrB=sB,
             exp_32_64=np.c_[unit_exp["32|64"]["x"], unit_exp["32|64"]["y"], unit_exp["32|64"]["ODI"]],
             exp_32_32=np.c_[unit_exp["32|32"]["x"], unit_exp["32|32"]["y"], unit_exp["32|32"]["ODI"]])
    f = lambda v: f"{v:+.3f}"  # noqa: E731
    L = [f"# 원인 ① 원값 수준 분해 — {PAIR}(탐색 · 판정 아님 · {NOW})", "",
         f"- 단위 {res['n_units']} · 피험자 {res['n_subjects']} · 평가 SNR(1/SNR² 중앙의 환산) A {res['eval_point']['snrA_equiv']:.1f} · B {res['eval_point']['snrB_equiv']:.1f}",
         f"- SNR 방향: OXF1PRI 중앙 {snr_dir['OXF1PRI_median']:.2f} · OXF2PRI 중앙 {snr_dir['OXF2PRI_median']:.2f} · ΔSNR(B − A) 중앙 {snr_dir['dSNR_B_minus_A_median']:+.2f}",
         f"- 읽기(계산 전 규칙): **{reading}**", "",
         "| 성분 | 실측 평균 [피험자 묶음 95 %] | 실측 SD | 기대 32\\|64 [95 %] | 기대 32\\|32 [95 %] | Δ_L | r [95 %] | 민감도 32\\|64 · 32\\|32 (s1 · s2 · s3) |", "|---|---|---|---|---|---|---|---|"]
    for c in ("x", "ODI", "y"):
        R_, E6, E3 = real[c], exp_["32|64"][c], exp_["32|32"][c]
        L.append(f"| {c} | {f(R_['mean'])} [{f(R_['boot95'][0])}, {f(R_['boot95'][1])}] | {R_['sd']:.3f} | {f(E6['primary'])} [{f(E6['boot95'][0])}, {f(E6['boot95'][1])}] | "
                 f"{f(E3['primary'])} [{f(E3['boot95'][0])}, {f(E3['boot95'][1])}] | {f(rr[c]['delta_L'])} | "
                 f"{rr[c]['r']:+.2f} [{rr[c]['r_boot95'][0]:+.2f}, {rr[c]['r_boot95'][1]:+.2f}]{'' if rr[c]['r_defined'] else ' (정의 안 됨)'} | "
                 f"{f(E6['s1_inv_snr'])} · {f(E6['s2_window_intercept'])} · {f(E6['s3_clamped15'])} / {f(E3['s1_inv_snr'])} · {f(E3['s2_window_intercept'])} · {f(E3['s3_clamped15'])} |")
    L.append(f"\n- 겉보기 NDI 차이(참고): 평균 {f(real['NDI_apparent']['mean'])} · 중앙 {f(real['NDI_apparent']['median'])} · SD {real['NDI_apparent']['sd']:.3f}")
    (OUT / "cause1_raw.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
