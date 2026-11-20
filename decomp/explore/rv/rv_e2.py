"""Exploratory e2 tables: per-cell measurement coverage and thresholds, symmetric-cell thresholds from true-combination test cases, and coil-assumption sensitivity; not a verdict."""
import datetime, hashlib, json, pathlib, sys
import numpy as np
REPO = pathlib.Path("<REPO>")
sys.path.insert(0, str(REPO / "decomp" / "remeasure"))
import judge_v4 as J  # noqa: E402
X14 = pathlib.Path("<WORK>/stage1_work/drafts/x14")
D1 = J.R.WORK / "explore" / "d1"
OUT = J.R.WORK / "explore" / "rv"
COMBOS4 = ((32, 32), (32, 64), (64, 32), (64, 64))
TYPES = ("match", "te", "tetr", "tetr_design")
SYM = {"OXF1PRI~OXF2PRI": (32, 64), "OXF1PRI~OXF3TRI": (32, 32), "OXF2PRI~OXF3TRI": (64, 32)}
NOW = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
SHA_V4 = "f6003589c01ad0bd693a6c809cb78ec6e69c70907be55703e22774d47f240a2a"


def kmed_cov(cov):           # cov (m, 9) → mean for each k → median
    return float(np.median(cov.mean(0))) if len(cov) else None


def main():
    b = (J.R.WORK / "onharmony" / "v4" / "v4_c3.json").read_bytes()
    assert hashlib.sha256(b).hexdigest() == SHA_V4
    with (X14 / "REREAD_LOG.jsonl").open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(dict(시각=NOW, 다시읽음=str(J.R.WORK / "onharmony" / "v4" / "v4_c3.json"), sha256=SHA_V4, 바이트=len(b), 무엇="rv e2 — 칸 문턱 · 덮음 · 묶음 수", 스크립트=__file__), ensure_ascii=False) + "\n")
    v4 = json.loads(b)
    seed = v4["seed_numpy"]
    rt = json.loads((D1 / "d1_rulethr.json").read_text(encoding="utf-8"))
    ref = {t: dict(np.load(D1 / f"d1_ref_{t}.npz", allow_pickle=False)) for t in TYPES}
    res = dict(made=NOW, script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(), table_1b={}, symmetric_2c={}, coil_sensitivity_4={})
    # 1b
    for k, v in rt["cells"].items():
        nc = len(v["candidates"])
        pair = k.split("|")[0]
        lab = "보정됨(후보 1)" if nc == 1 else ("약함(문턱 0.10 – 0.20)" if nc == 2 else "정보 없음(문턱 0)")
        if pair in SYM and nc > 1:
            lab += " · 4a 뒤 대칭 칸(2c)"
        res["table_1b"][k] = dict(candidates=nc, coverage_real=v["coverage_real"], threshold_rule=v["threshold_rule_ours"], threshold_stored=v["threshold_stored"], label=lab)
    # 2c
    for k, c4 in v4["main"]["cells"].items():
        pair = k.split("|")[0]
        if pair not in SYM or c4.get("verdict") == "자료 없음":
            continue
        comp = k.split("|")[1]
        ci, t = J.COMPS.index(comp), J.acq_type(pair)
        rf = ref[t]
        true = np.array([tuple(x) == SYM[pair] for x in rf["true_combo"].tolist()])
        j = COMBOS4.index(SYM[pair])
        cov, ok, flag = rf["cov"][:, j, :, ci], rf["ok"][:, j], rf["score"][:, j] > float(rf["threshold"])
        n_cl = int(c4["n_subjects"])
        g = np.random.default_rng([seed, J.ARMS.index("c3"), J.S.ACQ_TYPES.index(t), ci + 1, n_cl, 0])
        use = true & ok.all(1) & ~flag
        res["symmetric_2c"][k] = dict(acq_type=t, true_combo=f"{SYM[pair][0]}|{SYM[pair][1]}", n_cases=int(true.sum()), n_used=int(use.sum()),
                                      threshold_true_combo_only=J.threshold_from_test(cov[true], ok[true], flag[true], n_cl, g),
                                      threshold_stored=c4["threshold"], test_coverage_kmed=kmed_cov(cov[use]), coverage_real_all_candidates=c4["coverage"],
                                      note="실측 쪽 «그 조합의 덮음률» 은 D1 real 의 조합별 덮음으로(오늘 밤)")
    # 4
    def stats(rf, sel, ja, ci):
        ok = rf["ok"][:, ja].all(1) & sel
        if not ok.any():
            return None
        cov = rf["cov"][ok, ja, :, ci]
        z = (rf["mean"][ok, ja, :, ci] - rf["theta"][ok][:, None, ci]) / rf["sd"][ok, ja, :, ci]
        pit = rf["pit"][ok, ja, :, ci]
        return dict(n=int(ok.sum()), coverage=float(np.median(cov.mean(0))), z_mean=float(np.median(z.mean(0))), pit_mean=float(np.median(pit.mean(0))))
    full, side = {}, {}
    for comp, ci in (("NDI", 0), ("ODI", 2)):
        for t in TYPES:
            rf = ref[t]
            tc = [tuple(x) for x in rf["true_combo"].tolist()]
            for tt in COMBOS4:
                sel = np.array([x == tt for x in tc])
                for aa in COMBOS4:
                    full[f"{comp}|{t}|true {tt[0]}|{tt[1]}|assumed {aa[0]}|{aa[1]}"] = stats(rf, sel, COMBOS4.index(aa), ci)
        for s_i, s_name in ((1, "B"), (0, "A")):
            for tr_c in (32, 64):
                for as_c in (32, 64):
                    agg = []
                    for t in TYPES + ("pooled",):
                        rows = []
                        for tt in COMBOS4:
                            if tt[s_i] != tr_c:
                                continue
                            aa = list(tt); aa[s_i] = as_c; aa = tuple(aa)
                            for tl in (TYPES if t == "pooled" else (t,)):
                                v = full[f"{comp}|{tl}|true {tt[0]}|{tt[1]}|assumed {aa[0]}|{aa[1]}"]
                                if v:
                                    rows.append(v)
                        if rows:
                            w = np.array([r["n"] for r in rows], float)
                            agg.append((t, dict(n=int(w.sum()), coverage=float(np.average([r["coverage"] for r in rows], weights=w)),
                                                z_mean=float(np.average([r["z_mean"] for r in rows], weights=w)), pit_mean=float(np.average([r["pit_mean"] for r in rows], weights=w)))))
                    side[f"{comp}|side {s_name}|{tr_c}->{as_c}"] = dict(agg)
    res["coil_sensitivity_4"] = dict(rule="§22-4 · 아홉 k 모두 찬 경우 · k 중앙 · 쪽별 = 한쪽만 틀림(다른 쪽은 참과 같음) · pooled = 층 넷의 경우 수 가중 평균", side=side, full=full)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "rv_e2.json").write_text(json.dumps(res, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    f = lambda x: "—" if x is None else f"{x:.3f}"  # noqa: E731
    L = ["# e 회신 처리(계획서 §22 · 스크립트 산출 · 기술 · 판정 아님)", "", "## 1b 칸별 — 실측 덮음 · 규칙 맞춘 문턱 · 원래 문턱", "",
         "| 칸 | 후보 | 실측 덮음 | 규칙 맞춘 문턱 | 원래 문턱 | 표지 |", "|---|---|---|---|---|---|"]
    for k, v in res["table_1b"].items():
        L.append(f"| {k.replace('|', ' · ')} | {v['candidates']} | {f(v['coverage_real'])} | {f(v['threshold_rule'])} | {f(v['threshold_stored'])} | {v['label']} |")
    L += ["", "## 2c 대칭 칸 — 참 조합 경우만의 문턱", "", "| 칸 | 참 조합 | 시험 경우(씀) | 그 조합만의 문턱 | 원래 문턱 | 시험 덮음 | 실측 덮음(후보 모두 규칙) |", "|---|---|---|---|---|---|---|"]
    for k, v in res["symmetric_2c"].items():
        L.append(f"| {k.replace('|', ' · ')} | {v['true_combo']} | {v['n_cases']}({v['n_used']}) | {f(v['threshold_true_combo_only'])} | {f(v['threshold_stored'])} | {f(v['test_coverage_kmed'])} | {f(v['coverage_real_all_candidates'])} |")
    L += ["", "## 4 코일 가정 민감도(쪽별 · 층 넷 합침) — 덮음 · z 평균 · PIT 평균", "", "| 성분 · 쪽 · 참 → 가정 | n | 덮음 | z 평균 | PIT 평균 |", "|---|---|---|---|---|"]
    for k, v in side.items():
        p = v.get("pooled")
        if p:
            L.append(f"| {k.replace('|', ' · ')} | {p['n']} | {p['coverage']:.3f} | {p['z_mean']:+.3f} | {p['pit_mean']:.3f} |")
    L += ["", "match 층(OXF1PRI~OXF2PRI 의 층)만:", "", "| 성분 · 쪽 · 참 → 가정 | n | 덮음 | z 평균 | PIT 평균 |", "|---|---|---|---|---|"]
    for k, v in side.items():
        p = v.get("match")
        if p:
            L.append(f"| {k.replace('|', ' · ')} | {p['n']} | {p['coverage']:.3f} | {p['z_mean']:+.3f} | {p['pit_mean']:.3f} |")
    (OUT / "rv_e2.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
