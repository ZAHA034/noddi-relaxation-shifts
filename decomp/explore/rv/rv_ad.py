"""Exploratory analyses a and d: location-matched expected misspecification exceedance rate, and a repeat of the D5 between-scanner prediction using SNR-scaled retest; not a verdict."""
from __future__ import annotations

import datetime
import hashlib
import json
import math
import pathlib
import sys

import numpy as np

REPO = pathlib.Path("<REPO>")
X14 = pathlib.Path("<WORK>/stage1_work/drafts/x14")
sys.path.insert(0, str(REPO / "decomp" / "remeasure"))
import judge_v4 as J  # noqa: E402

OUT = J.R.WORK / "explore" / "rv"
V4 = J.R.WORK / "onharmony" / J.V4_DIR
REREAD = X14 / "REREAD_LOG.jsonl"
NOW = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
SHA = dict(attr_table="5184137854abc69710904889503eaf8dbb111ca9736e46af13b49067ea65b90d",
           snr_table="1714e6df9e2221df4abc59d58dc40f85c81d5e10c4406a4bcb2e36caddbd8cec",
           v4_c3="f6003589c01ad0bd693a6c809cb78ec6e69c70907be55703e22774d47f240a2a")
COORD = ("x", "y", "ODI")


def _blocked(*a, **k):
    raise RuntimeError("⛔ V4 결과 파일 · 열람 기록을 쓰지 않는다")


J._write_once = J._write_json = J._log = _blocked


def reread(path, what, expect=None) -> bytes:
    b = pathlib.Path(path).read_bytes()
    s = hashlib.sha256(b).hexdigest()
    with REREAD.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(dict(시각=NOW, 다시읽음=str(path), sha256=s, 바이트=len(b), 무엇=what, 스크립트=__file__), ensure_ascii=False) + "\n")
    if expect and s != expect:
        raise RuntimeError(f"⛔ {path} sha256 이 열람 기록과 다르다")
    return b


def save(name, obj):
    OUT.mkdir(parents=True, exist_ok=True)
    obj = dict(made=NOW, script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(), **obj)
    (OUT / name).write_text(json.dumps(obj, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")


# ── a ──────────────────────────────────────────────────────────────────────────
def mode_a():
    at = json.loads(reread(V4 / "attr_table.json", "rv a — 속성 표", SHA["attr_table"]))
    v4 = json.loads(reread(V4 / "v4_c3.json", "rv a — 저장된 칸 초과율", SHA["v4_c3"]))
    root = J.R.WORK / "onharmony"
    for p in sorted({root / f"{sub}_{s}" / "roi.json" for d in at["pairs"].values() for sub, x in d.items() for s in x["sessions"]}):
        reread(p, "rv a — roi.json «주»(collect)")
    rows = J.collect(J.R.WORK, at["pairs"])
    sim = np.load(J.S.WORK / "sim_test.npz", allow_pickle=False)
    tc = np.load(J.S.WORK / "v4" / "test_cover_c3.npz", allow_pickle=False)
    res = dict(rule="a — 실측 단위 X3 에 가장 가까운 시험 경우(같은 획득 유형 층 · 층 SD 로 표준화한 6 차원 유클리드) k 개의 오지정 비율 → 칸 평균 = 기대 · "
                    "관측 = 저장된 exceed_rate · 비 = 기대 ÷ 관측 · «관측의 몇 %» = 칸 비의 중앙값(주) · 합산 비(계획서 §19)", ks=[25, 50, 100], cells={})
    for key, c4 in v4["main"]["cells"].items():
        if c4.get("verdict") == "자료 없음":
            continue
        pair = key.split("|")[0]
        t = J.acq_type(pair)
        idx = np.flatnonzero(sim["layer"] == f"{t}|L0")
        Xt = sim["x_c3"][idx].astype(float)
        ft = tc[f"{t}_flag"].astype(float)
        assert len(ft) == len(Xt)
        sd = Xt.std(0)
        r = [x for x in rows[pair] if not x["main_excluded"]]
        Xr = np.array([np.r_[x["x"][0:3], x["x"][7:10]] for x in r], float)
        D = np.sqrt((((Xr[:, None, :] - Xt[None, :, :]) / sd) ** 2).sum(-1))
        order = np.argsort(D, 1)
        exp = {k: float(ft[order[:, :k]].mean(1).mean()) for k in res["ks"]}
        res["cells"][key] = dict(acq_type=t, n_units=len(r), observed=c4["exceed_rate"], expected=exp,
                                 ratio={k: (v / c4["exceed_rate"] if c4["exceed_rate"] else None) for k, v in exp.items()},
                                 nn_dist_median=float(np.median(D[np.arange(len(D)), order[:, 49]])), test_flag_rate_layer=float(ft.mean()))
    C = res["cells"]
    for k in res["ks"]:
        rat = [v["ratio"][k] for v in C.values() if v["ratio"][k] is not None]
        res[f"summary_k{k}"] = dict(ratio_median=float(np.median(rat)), ratio_min=float(min(rat)), ratio_max=float(max(rat)),
                                    pooled_ratio=float(sum(v["expected"][k] * v["n_units"] for v in C.values()) / sum(v["observed"] * v["n_units"] for v in C.values())),
                                    expected_range=[float(min(v["expected"][k] for v in C.values())), float(max(v["expected"][k] for v in C.values()))])
    res["observed_range"] = [float(min(v["observed"] for v in C.values())), float(max(v["observed"] for v in C.values()))]
    save("rv_a_location.json", res)
    print(json.dumps({k: v for k, v in res.items() if k.startswith("summary") or k == "observed_range"}, ensure_ascii=False))
    for k, v in C.items():
        print(k, v["n_units"], "obs %.3f exp50 %.3f ratio %.2f" % (v["observed"], v["expected"][50], v["ratio"][50]))


# ── d ──────────────────────────────────────────────────────────────────────────
def mode_d():
    sys.path.insert(0, str(X14 / "D5"))
    import importlib.util
    spec = importlib.util.spec_from_file_location("d5a", X14 / "D5" / "d5_analysis.py")
    A = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(A)
    E, T, _, units = A.load_feats()
    at = json.loads(reread(V4 / "attr_table.json", "rv d — 쌍 × 피험자의 두 세션", SHA["attr_table"]))
    snr = json.loads(reread(V4 / "snr_table.json", "rv d — 세션 SNR", SHA["snr_table"]))["sessions"]
    d2 = json.loads(reread(J.R.WORK / "explore" / "d2" / "d2_repeat.json", "rv d — 시뮬레이터 SNR 척도(D2 등고선)"))
    alpha_c = {c: math.log(d2["by_coord"][c]["contour_pooled_ratio"]["snr30"]["50"] / d2["by_coord"][c]["contour_pooled_ratio"]["snr20"]["50"]) / math.log(1.5) for c in COORD}
    ss, dfs = np.zeros((27, 3)), np.zeros((27, 3))
    for s in ("sub-13192", "sub-16974"):
        V = np.array([T[f"{s}_ses-OXF1PRI00{i}"] for i in range(1, 7)])
        ss += np.nansum((V - np.nanmean(V, 0)) ** 2, 0)
        dfs += np.maximum(np.isfinite(V).sum(0) - 1, 0)
    sw = np.sqrt(ss / np.where(dfs > 0, dfs, np.nan))
    jens = np.sqrt((dfs - 2) / dfs)
    snr_ref = float(np.median([snr[f"{s}_ses-OXF1PRI00{i}"]["snr"] for s in ("sub-13192", "sub-16974") for i in range(1, 7)]))
    rng = np.random.default_rng(20261115)
    res = dict(rule="d — D5 ② 의 LOSO 오차를 그대로 · R = RMS(e/σ_d,s,r) × √((ν−2)/ν) · σ_d,s,r = σ_w,r √((SNR_ref/SNR_A)^(2α) + (SNR_ref/SNR_B)^(2α)) · "
                    "α = 1(주) · α_c(시뮬레이터 · D2 등고선) · 읽기는 ② 와 같음 · «결과를 본 뒤»(계획서 §19)", snr_ref=snr_ref, alpha_c=alpha_c,
               df_median=float(np.median(dfs)), pairs={})
    for pair, d in sorted(at["pairs"].items()):
        subs = sorted(s for s, v in d.items() if not v["main_excluded"])
        if len(subs) < 4:
            continue
        A_ = np.array([T[f"{s}_{d[s]['sessions'][0]}"] for s in subs])
        B_ = np.array([T[f"{s}_{d[s]['sessions'][1]}"] for s in subs])
        sA = np.array([snr[f"{s}_{d[s]['sessions'][0]}"]["snr"] for s in subs])
        sB = np.array([snr[f"{s}_{d[s]['sessions'][1]}"]["snr"] for s in subs])
        n = len(subs)
        exp = math.sqrt(n / (n - 1))
        pr = dict(n_subjects=n, expected_R_M1=exp, coords={})
        for j, c in enumerate(COORD):
            a, b = A_[:, :, j], B_[:, :, j]
            out = {}
            for kind in ("M0", "M1"):
                e = b - A.loso_predict(a, b, kind)
                for lab, al in (("alpha1", 1.0), ("alpha_c", alpha_c[c])):
                    scale = np.sqrt((snr_ref / sA) ** (2 * al) + (snr_ref / sB) ** (2 * al))[:, None]
                    z = e / (sw[None, :, j] * scale) * jens[None, :, j]
                    Rv = float(np.sqrt(np.nanmean(z ** 2)))
                    boots = [float(np.sqrt(np.nanmean(z[rng.integers(0, n, n)] ** 2))) for _ in range(2000)]
                    lo, hi = float(np.quantile(boots, 0.025)), float(np.quantile(boots, 0.975))
                    o = dict(R=Rv, ci95=[lo, hi])
                    if kind == "M1":
                        o["reading"] = ("재검사 수준으로 대응됨" if hi <= 1.5 * exp else "재검사 수준 밖" if lo > 1.5 * exp else "가를 수 없음")
                    out[f"{kind}|{lab}"] = o
            pr["coords"][c] = out
        res["pairs"][pair] = pr
    cnt = {}
    for lab in ("alpha1", "alpha_c"):
        for c in COORD:
            rd = [v["coords"][c][f"M1|{lab}"]["reading"] for v in res["pairs"].values()]
            cnt[f"{lab}|{c}"] = {k: rd.count(k) for k in sorted(set(rd))}
    res["reading_counts"] = cnt
    save("rv_d_retest.json", res)
    print(json.dumps(dict(snr_ref=snr_ref, alpha_c=alpha_c, reading_counts=cnt), ensure_ascii=False))


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    dict(a=mode_a, d=mode_d).get(cmd, lambda: sys.exit("쓰임: a|d"))()
