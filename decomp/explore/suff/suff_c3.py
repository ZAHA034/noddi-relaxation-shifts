"""Exploratory sufficiency test on c3: redraws OXF1PRI~OXF2PRI along the verdict path under three coil conditions (32|64, 32|32, 64|64); CPU only, not a verdict."""
from __future__ import annotations

import datetime
import hashlib
import json
import math
import os
import pathlib
import sys
import time

import numpy as np

REPO = pathlib.Path("<REPO>")
DR = pathlib.Path("<WORK>/stage1_work/drafts")
sys.path.insert(0, str(REPO / "decomp" / "remeasure"))
sys.path.insert(0, str(DR / "x14" / "D1"))
import judge_v4 as J  # noqa: E402
from d1_c3 import Rec  # noqa: E402  — draw_fn wrapping (no change in computation or random-number consumption)

PAIR = "OXF1PRI~OXF2PRI"
COMBOS = ((32, 64), (32, 32), (64, 64))
COMPS = {"NDI": 0, "ODI": 2}
THR = {"NDI": 0.80, "ODI": 0.75}                                   # Symmetric threshold of that cell (e · d1_rulethr.json threshold_true_combo_ours)
TAIL = 0.05
OUT = J.R.WORK / "explore" / "suff"
D1 = J.R.WORK / "explore" / "d1"
RV = J.R.WORK / "explore" / "rv"
V4 = J.R.WORK / "onharmony" / J.V4_DIR
REREAD = DR / "x14" / "REREAD_LOG.jsonl"
SHA = dict(attr_table="5184137854abc69710904889503eaf8dbb111ca9736e46af13b49067ea65b90d",
           v4_c3="f6003589c01ad0bd693a6c809cb78ec6e69c70907be55703e22774d47f240a2a")
NOW = datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def _blocked(*a, **k):
    raise RuntimeError("⛔ 충분성 시험은 V4 결과 파일 · 열람 정본을 쓰지 않는다")


J._write_once = J._write_json = J._log = _blocked


def reread(path, what, expect=None) -> bytes:
    b = pathlib.Path(path).read_bytes()
    s = hashlib.sha256(b).hexdigest()
    with REREAD.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(dict(시각=datetime.datetime.now().astimezone().isoformat(timespec="seconds"), 다시읽음=str(path), sha256=s, 바이트=len(b),
                                 무엇=what, 스크립트=__file__), ensure_ascii=False) + "\n")
    if expect is not None and s != expect:
        raise RuntimeError(f"⛔ {path} sha256 {s[:16]}… ≠ 열람 기록 {expect[:16]}…")
    return b


def status(msg: str):
    (OUT / "STATUS_SUFF").write_text(f"{datetime.datetime.now():%F %T} {msg}\n", encoding="utf-8")


def cov90(mu, s):
    F = lambda x: 0.5 * (1 + math.erf(x / math.sqrt(2)))  # noqa: E731
    q = 1.6448536269514722
    return F((q - mu) / s) - F((-q - mu) / s)


def stats(P, Z, C):
    """P · Z (n, 9) · C (n, 9) bool → for each k → k median."""
    k = lambda v: float(np.nanmedian(v))  # noqa: E731
    return dict(n=int(len(P)), coverage=k(C.mean(0)), z_mean=k(Z.mean(0)), z_sd=k(Z.std(0, ddof=1)), pit_mean=k(P.mean(0)),
                tail=k(((P < TAIL) | (P > 1 - TAIL)).mean(0)))


def expected():
    """Expectation before computing (the review-session estimate, using the match-stratum values) — z under the same L assumption ≈ D1 observed z (32|64) + the per-side shift of the test data · if the z SD is unchanged, normal-approximation coverage."""
    side = json.loads((RV / "rv_e2.json").read_text(encoding="utf-8"))["coil_sensitivity_4"]["side"]
    obs = json.loads((D1 / "d1_post_sym.json").read_text(encoding="utf-8"))["sym"]
    shift = {"32|32": "side B|64->32", "64|64": "side A|32->64"}           # 32|32 = B set smaller than the actual (64) · 64|64 = A set larger than the actual (32)
    out = {}
    for comp in COMPS:
        r = obs[f"{PAIR}|{comp}"]["real"]
        for cn, sk in shift.items():
            m = side[f"{comp}|{sk}"]["match"]
            mu = r["z_mean"] + m["z_mean"]
            out[f"{cn}|{comp}"] = dict(observed_32_64_z=r["z_mean"], sim_shift_match=m["z_mean"], sim_shift_key=f"{comp}|{sk}|match", expected_z=mu,
                                       z_sd_kept=r["z_sd"], expected_coverage_normal=cov90(mu, r["z_sd"]))
    return out


def reading(res):
    """§24-1: support = all four have coverage ≥ threshold and |z mean| < 1 · weakening = for any component, in both combinations (coverage < 0.5 or |z mean| > 2) · anything else = judgment withheld."""
    four = [res[f"{cn}|{c}"] for cn in ("32|32", "64|64") for c in COMPS]
    if all(v["coverage"] >= THR[v["comp"]] and abs(v["z_mean"]) < 1 for v in four):
        return "지지"
    for c in COMPS:
        if all(res[f"{cn}|{c}"]["coverage"] < 0.5 or abs(res[f"{cn}|{c}"]["z_mean"]) > 2 for cn in ("32|32", "64|64")):
            return "약화"
    return "판단 보류"


def main():
    assert os.environ.get("CUDA_VISIBLE_DEVICES", None) == "", "CUDA_VISIBLE_DEVICES= 로 띄운다(CPU 만)"
    assert os.environ.get("S4_THREADS") == "16", "S4_THREADS=16(D1 · 원 실행과 같게)"
    for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        assert v not in os.environ, f"{v} 는 정하지 않는다(D1 과 같게)"
    OUT.mkdir(parents=True, exist_ok=True)
    status("시작")
    exp = expected()
    at = json.loads(reread(V4 / "attr_table.json", "충분성 — 속성 표(판정 때 읽은 것)", SHA["attr_table"]))
    v4 = json.loads(reread(V4 / "v4_c3.json", "충분성 — 쌍 시드(재현 문)", SHA["v4_c3"]))
    sub_at = {PAIR: at["pairs"][PAIR]}
    root = J.R.WORK / "onharmony"
    for f in sorted({root / f"{sub}_{s}" / "roi.json" for sub, x in sub_at[PAIR].items() for s in x["sessions"]}):
        reread(f, "충분성 — roi.json «주»(collect · 판정 때 읽은 것 · OXF1PRI~OXF2PRI 만)")
    draw_fn, logq_fn, K = J.real_fns("c3")
    miss = json.loads((J.RESULTS / "miss.json").read_text(encoding="utf-8"))
    rows = J.collect(J.R.WORK, sub_at)[PAIR]
    pi = J.card_pairs().index(PAIR)
    seed = v4["seeds_torch"][f"SNR30|Neff50|{PAIR}"]
    assert seed == 6000 + pi
    X7 = np.array([r["x"] for r in rows])
    X3 = np.c_[X7[:, 0:3], X7[:, 7:10]]
    te, tr = np.array([r["te"] for r in rows], float), np.array([r["tr"] for r in rows], float)
    conds = [J.cond_rows(PAIR, te, tr, cc, 30.0, 50.0) for cc in COMBOS]
    rec = Rec(draw_fn, len(rows), len(COMBOS), len(K))
    t0 = time.time()
    status(f"사후 표본 · 조합 {len(COMBOS)} × k {len(K)} · 단위 {len(rows)}")
    J.seed_torch(seed)
    cov, ok, _, okbc = J.cover_by_k(rec, X3, conds, K)
    assert rec.c == len(COMBOS) * len(K)
    # Statement of reproduction — first combination = D1 stored value
    d = dict(np.load(D1 / "d1_c3_units" / f"{PAIR}.npz", allow_pickle=False))
    same = lambda a, b: bool(np.array_equal(a, b, equal_nan=True))  # noqa: E731
    gate = dict(pit=same(rec.pit[:, 0], d["pit"][:, 0]), mean=same(rec.mean[:, 0], d["mean"][:, 0]), sd=same(rec.sd[:, 0], d["sd"][:, 0]),
                okc=same(rec.ok[:, 0], d["okc"][:, 0]), covc=same(rec.cov[:, 0], d["covc"][:, 0]), n_units=len(rows) == len(d["flag"]))
    gate["passed"] = all(gate.values())
    if not gate["passed"]:
        (OUT / "suff_gate.json").write_text(json.dumps(dict(gate=gate, seconds=round(time.time() - t0, 1)), ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
        status("멈춤 — 재현 문 실패(값 없음)")
        sys.exit(2)
    flag3, _, fbc = J.miss_flags(logq_fn, X3, conds, K, miss["layers"][f"{J.acq_type(PAIR)}|L0"]["threshold"])
    use = ~d["main_excluded"] & ~d["flag"] & d["ok"].all(1)                 # same ROI unit as D1
    res = {}
    for j, cc in enumerate(COMBOS):
        cn = f"{cc[0]}|{cc[1]}"
        uj = use & rec.ok[:, j].all(1)
        for comp, ci in COMPS.items():
            with np.errstate(invalid="ignore", divide="ignore"):
                Z = rec.mean[:, j, :, ci] / rec.sd[:, j, :, ci]
                Z0 = rec.mean[:, 0, :, ci] / rec.sd[:, 0, :, ci]
            s = stats(rec.pit[uj][:, j, :, ci], Z[uj], rec.cov[uj][:, j, :, ci])
            s.update(comp=comp, threshold=THR[comp], n_flag_this_combo=int(fbc[use, j].sum()),
                     dz_vs_32_64=float(np.nanmedian((Z[uj] - Z0[uj]).mean(0))))
            if cn != "32|64":
                s["expected"] = exp[f"{cn}|{comp}"]
            res[f"{cn}|{comp}"] = s
    rd = reading(res)
    out = dict(made=NOW, script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(), pair=PAIR, combos=[f"{a}|{b}" for a, b in COMBOS],
               seed=seed, seconds=round(time.time() - t0, 1), gate=gate, n_rows=len(rows), n_use_d1=int(use.sum()),
               rule="§24-1 — 지지 = 같은 L 두 조합 × NDI · ODI 넷 모두 덮음 ≥ 대칭 문턱(NDI 0.80 · ODI 0.75)이고 |z 평균| < 1 · "
                    "약화 = 어느 성분이든 두 조합 모두에서 덮음 < 0.5 또는 |z 평균| > 2 · 그 밖 = 판단 보류",
               reading=rd, cells=res, note="탐색 · 판정 아님 · 판단기 점 추정은 내지 않음(z · PIT 는 보정 진단)")
    (OUT / "suff_c3.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    L = [f"# 충분성 시험 — {PAIR} 같은 L 가정(탐색 · 판정 아님 · {NOW})", "", f"- 재현 문(32|64 = D1 저장값): {'통과' if gate['passed'] else '실패'} · 단위 {len(rows)} · D1 쓰는 단위 {int(use.sum())} · {out['seconds']} 초",
         f"- 읽기(계산 전 규칙): **{rd}**", "",
         "| 조합 · 성분 | 단위 | 덮음 | 문턱 | z 평균 | z SD | PIT 평균 | 양 끝 | Δz(− 32\\|64) | 기대 z · 기대 이동 · 기대 덮음 |", "|---|---|---|---|---|---|---|---|---|---|"]
    for k, v in res.items():
        e = v.get("expected")
        es = f"{e['expected_z']:+.2f} · {e['sim_shift_match']:+.2f} · {e['expected_coverage_normal']:.2f}" if e else "—"
        L.append(f"| {k.replace('|', ' ', 1).replace('|', ' · ')} | {v['n']} | {v['coverage']:.3f} | {v['threshold']} | {v['z_mean']:+.2f} | {v['z_sd']:.2f} | "
                 f"{v['pit_mean']:.2f} | {v['tail']:.2f} | {v['dz_vs_32_64']:+.2f} | {es} |")
    (OUT / "suff_c3.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    status(f"끝 · 읽기 {rd}")
    print("\n".join(L), flush=True)


if __name__ == "__main__":
    main()
