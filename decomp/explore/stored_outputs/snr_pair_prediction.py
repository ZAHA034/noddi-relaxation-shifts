"""Computes per-cell ranges of the noise hypothesis predicted ratio 30 ÷ SNR_pair, plus ΔSNR and SNR ratio, from stored attribute and SNR tables only."""
import datetime, hashlib, json, math, pathlib, statistics as st
V = pathlib.Path("<WORK>/remeasure/onharmony/v4")
X = pathlib.Path(__file__).resolve().parent
now = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
def load(name, purpose):
    p = V / name; b = p.read_bytes(); s = hashlib.sha256(b).hexdigest()
    with (X / "REREAD_LOG.jsonl").open("a", encoding="utf-8") as f:
        f.write(json.dumps(dict(시각=now, 다시읽음=str(p), sha256=s, 바이트=len(b), 무엇=purpose, 스크립트=__file__), ensure_ascii=False) + "\n")
    return json.loads(b), s
attr, sa = load("attr_table.json", "1-3 예측 비 · 1-2 ΔSNR — 쌍 × 피험자의 두 세션")
snr, ss = load("snr_table.json", "1-3 예측 비 · 1-2 ΔSNR — 세션 SNR")
v4, sv = load("v4_c3.json", "1-3 — 자료 있는 칸 목록")
S = snr["sessions"]; lo, hi = snr["train_range"]
cells = [k for k, c in v4["main"]["cells"].items() if c.get("verdict") != "자료 없음"]
rows = []
for k in sorted(cells):
    pair = k.split("|")[0]
    pr, d, rt, out_any = [], [], [], 0
    for sub, at in sorted(attr["pairs"].get(pair, {}).items()):
        if at["main_excluded"]:
            continue
        sA, sB = (S[f"{sub}_{ses}"]["snr"] for ses in at["sessions"])
        sp = math.sqrt(2.0 / (1.0 / sA ** 2 + 1.0 / sB ** 2))
        pr.append(30.0 / sp); d.append(sB - sA); rt.append(sB / sA)
        out_any += int(not (lo <= sA <= hi) or not (lo <= sB <= hi))
    rows.append(dict(cell=k, n_subjects=len(pr), n_subjects_out_of_train=out_any,
                     pred_ratio=dict(min=round(min(pr), 3), median=round(st.median(pr), 3), max=round(max(pr), 3)),
                     dSNR=dict(min=round(min(d), 2), median=round(st.median(d), 2), max=round(max(d), 2)),
                     SNR_ratio_B_over_A=dict(min=round(min(rt), 3), median=round(st.median(rt), 3), max=round(max(rt), 3))))
allp = [r["pred_ratio"] for r in rows]
out = dict(what=__doc__.splitlines()[0], made=now, inputs=dict(attr_table=sa, snr_table=ss, v4_c3=sv), rows=rows,
           overall=dict(pred_ratio_min=min(p["min"] for p in allp), pred_ratio_max=max(p["max"] for p in allp),
                        cell_median_range=[min(p["median"] for p in allp), max(p["median"] for p in allp)]))
(X / "snr_pair_prediction.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
lab = lambda s: s.replace("|", " · ")
L = ["| 칸 | 피험자 | 범위 밖 피험자 | 예측 비 30 ÷ SNR_pair 최소 · 중앙 · 최대 | ΔSNR(B − A) 최소 · 중앙 · 최대 | SNR 비 B/A 중앙 |", "|---|---|---|---|---|---|"]
for r in rows:
    p, dd, q = r["pred_ratio"], r["dSNR"], r["SNR_ratio_B_over_A"]
    L.append(f"| {lab(r['cell'])} | {r['n_subjects']} | {r['n_subjects_out_of_train']} | {p['min']:.2f} · {p['median']:.2f} · {p['max']:.2f} | {dd['min']:+.1f} · {dd['median']:+.1f} · {dd['max']:+.1f} | {q['median']:.2f} |")
(X / "snr_pair_prediction.md").write_text("\n".join(L) + "\n", encoding="utf-8")
print(json.dumps(out["overall"], ensure_ascii=False)); print("\n".join(L[:6]))
