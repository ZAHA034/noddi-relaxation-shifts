"""Reads stored d1 outputs to describe the four symmetric cells and the H7 selection within ROI units; exploratory, descriptive, not a verdict."""
import datetime, hashlib, json, math, pathlib
import numpy as np
D1 = pathlib.Path("<WORK>/remeasure/explore/d1")
UNITS = D1 / "d1_c3_units"
COMBOS4 = ((32, 32), (32, 64), (64, 32), (64, 64))
COMPS = ("NDI", "FWF", "ODI")
TAIL = 0.05
SYM = {"OXF1PRI~OXF2PRI|NDI": (32, 64), "OXF1PRI~OXF2PRI|ODI": (32, 64), "OXF1PRI~OXF3TRI|ODI": (32, 32), "OXF2PRI~OXF3TRI|ODI": (64, 32)}
THR = json.loads((D1 / "d1_rulethr.json").read_text(encoding="utf-8"))["cells"]
NOW = datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def kmed(x):
    return float(np.nanmedian(x))


def stats(P, Z, C=None):
    """P · Z (n, 9) · C (n, 9) bool → value for each k → k median."""
    out = dict(n=int(len(P)), pit_mean=kmed(P.mean(0)), tail=kmed(((P < TAIL) | (P > 1 - TAIL)).mean(0)), z_mean=kmed(Z.mean(0)),
               z_sd=kmed(Z.std(0, ddof=1)), absz_median=kmed(np.median(np.abs(Z), 0)))
    if C is not None:
        out["coverage"] = kmed(C.mean(0))
    return out


def cov90(mu, s):
    q = 1.6448536269514722
    F = lambda x: 0.5 * (1 + math.erf(x / math.sqrt(2)))  # noqa: E731
    return F((q - mu) / s) - F((-q - mu) / s)


def main():
    d1 = json.loads((D1 / "d1_c3.json").read_text(encoding="utf-8"))
    res = dict(made=NOW, script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
               inputs=dict(d1_c3_json_sha256=hashlib.sha256((D1 / "d1_c3.json").read_bytes()).hexdigest()), sym={}, h7={}, norm={})
    refs = {}
    for cell in d1["cells"]:
        pair, comp = cell.split("|")
        ci = COMPS.index(comp)
        d = dict(np.load(UNITS / f"{pair}.npz", allow_pickle=False))
        use = ~d["main_excluded"] & ~d["flag"] & d["ok"].all(1)                     # same ROI unit as used in the D1 summary
        names = d["combos"].tolist()
        with np.errstate(invalid="ignore", divide="ignore"):
            Zc = d["mean"][..., ci] / d["sd"][..., ci]                                 # Measured true value 0
        per = {}
        for j, cn in enumerate(names):
            per[cn] = stats(d["pit"][use][:, j, :, ci], Zc[use][:, j, :], d["covc"][use][:, j, :, ci])
        res["h7"][cell] = dict(combos=per, has_L4_only="32|32" in names)
        if cell in SYM:
            tc = f"{SYM[cell][0]}|{SYM[cell][1]}"
            t = d1["cells"][cell]["acq_type"]
            if t not in refs:
                refs[t] = dict(np.load(D1 / f"d1_ref_{t}.npz", allow_pickle=False))
            rf = refs[t]
            j4 = COMBOS4.index(SYM[cell])
            sel = np.array([tuple(x) == SYM[cell] for x in rf["true_combo"].tolist()])
            ok = sel & rf["ok"][:, j4].all(1) & ~(rf["score"][:, j4] > float(rf["threshold"]))
            with np.errstate(invalid="ignore", divide="ignore"):
                Zr = ((rf["mean"][:, j4] - rf["theta"][:, None, :]) / rf["sd"][:, j4])[..., ci]
            ref = stats(rf["pit"][ok][:, j4, :, ci], Zr[ok], rf["cov"][ok][:, j4, :, ci])
            real = per[tc]
            res["sym"][cell] = dict(true_combo=tc, real=real, ref=ref, threshold_true_combo_ours=THR[cell]["threshold_true_combo_ours"],
                                    threshold_stored=THR[cell]["threshold_stored"])
            mu, s, s0 = real["z_mean"], real["z_sd"], ref["z_sd"]
            res["norm"][cell] = dict(rule="[결과를 본 뒤] z ~ N(μ, s²) 의 90 % 덮음 · μ · s 는 실측 참 조합 · 기준은 시험 경우 z SD",
                                     both=cov90(mu, s), bias_only=cov90(mu, s0), spread_only=cov90(0.0, s), reference=cov90(ref["z_mean"], s0),
                                     observed=real["coverage"])
    out = D1 / "d1_post_sym.json"
    out.write_text(json.dumps(res, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    f = lambda x: f"{x:.2f}"  # noqa: E731
    print("## 대칭 칸 — 참 조합(실측 · 기준)")
    for c, v in res["sym"].items():
        r, b, nm = v["real"], v["ref"], res["norm"][c]
        print(f"{c} {v['true_combo']} n {r['n']} | 덮음 {f(r['coverage'])} 문턱(참 조합) {v['threshold_true_combo_ours']} 저장 {v['threshold_stored']} | "
              f"PIT {f(r['pit_mean'])}({f(b['pit_mean'])}) 양끝 {f(r['tail'])}({f(b['tail'])}) z {r['z_mean']:+.2f}({b['z_mean']:+.2f}) zSD {f(r['z_sd'])}({f(b['z_sd'])}) "
              f"| 정규: 둘 {f(nm['both'])} μ만 {f(nm['bias_only'])} s만 {f(nm['spread_only'])} 기준 {f(nm['reference'])} 관측 {f(nm['observed'])} · 기준 n {b['n']}")
    print("## H7 — 칸 안 조합마다(덮음 · PIT 평균 · |z| 중앙)")
    for c, v in res["h7"].items():
        print(c, " · ".join(f"{cn}: {f(s['coverage'])} {f(s['pit_mean'])} {f(s['absz_median'])}" for cn, s in v["combos"].items()))
    print(out, hashlib.sha256(out.read_bytes()).hexdigest())


if __name__ == "__main__":
    main()
