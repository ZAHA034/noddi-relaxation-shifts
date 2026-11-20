"""Describes the coil-only contrast between OXF1PRI~X and OXF2PRI~X for six scanners X, from stored outputs only, with no tests."""
import datetime, hashlib, json, pathlib, statistics as st
V = pathlib.Path("<WORK>/remeasure/onharmony/v4")
X = pathlib.Path(__file__).resolve().parent
now = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
def load(name):
    p = V / name; b = p.read_bytes(); s = hashlib.sha256(b).hexdigest()
    with (X / "REREAD_LOG.jsonl").open("a", encoding="utf-8") as f:
        f.write(json.dumps(dict(시각=now, 다시읽음=str(p), sha256=s, 바이트=len(b), 무엇="4b 코일만 다른 대조(저장된 출력)", 스크립트=__file__), ensure_ascii=False) + "\n")
    return json.loads(b), s
v4, s1 = load("v4_c3.json"); attr, s2 = load("attr_table.json"); snr, s3 = load("snr_table.json")
cells, A, S = v4["main"]["cells"], attr["pairs"], snr["sessions"]
FIG = {"NOT1ACH": (161, 64), "NOT2ING-98": (197, 49), "OXF3TRI": (162, 142), "NOT4GEP": (104, 47), "OXF4GEP": (87, 40), "NOT2ING-89": (47, 10)}
rows = []
for x, (f1, f2) in FIG.items():
    r = dict(X=x)
    for dev, fig in (("OXF1PRI", f1), ("OXF2PRI", f2)):
        pair = f"{dev}~{x}"; c = cells[f"{pair}|ODI"]
        subs = {s for s, a in A[pair].items() if not a["main_excluded"]}
        oxf_ses = [A[pair][s]["sessions"][0] for s in sorted(subs)]
        snrs = [S[f"{s}_{ses}"]["snr"] for s, ses in zip(sorted(subs), oxf_ses)]
        r[dev] = dict(cell=f"{pair}|ODI", n_subjects=c["n_subjects"], n_units=c["n_units"], misspec=c["n_excluded_misspec"], k_unfilled=c["n_units_dropped_k"],
                      used=c["n_units_used"], identity=c["n_units"] == c["n_units_used"] + c["n_excluded_misspec"] + c["n_units_dropped_k"],
                      coverage_median=round(c["coverage"], 3), k_range=[round(c["coverage_min"], 3), round(c["coverage_max"], 3)],
                      oxf_snr_median=round(st.median(snrs), 2), subjects=sorted(subs), fig_used=fig, fig_matches=(c["n_units_used"] == fig))
    s1_, s2_ = set(r["OXF1PRI"]["subjects"]), set(r["OXF2PRI"]["subjects"])
    r["subjects_overlap"] = dict(n_both=len(s1_ & s2_), n_only_OXF1PRI=len(s1_ - s2_), n_only_OXF2PRI=len(s2_ - s1_))
    rows.append(r)
out = dict(made=now, inputs=dict(v4_c3=s1, attr_table=s2, snr_table=s3), rows=rows)
(X / "coil_contrast.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
L = ["| X | 피험자(같은 사람 · 1 만 · 2 만) | OXF1PRI~X 단위 · 오지정 · k · 쓴 | OXF2PRI~X 단위 · 오지정 · k · 쓴 | 덮음 중앙 [아홉 k 범위] 1 / 2 | OXF 세션 SNR 중앙 1 / 2 | 그림 숫자 같음 |", "|---|---|---|---|---|---|---|"]
for r in rows:
    a, b, o = r["OXF1PRI"], r["OXF2PRI"], r["subjects_overlap"]
    L.append(f"| {r['X']} | {o['n_both']} · {o['n_only_OXF1PRI']} · {o['n_only_OXF2PRI']} | {a['n_units']} · {a['misspec']} · {a['k_unfilled']} · {a['used']} | {b['n_units']} · {b['misspec']} · {b['k_unfilled']} · {b['used']} | "
             f"{a['coverage_median']:.3f} [{a['k_range'][0]:.3f}–{a['k_range'][1]:.3f}] / {b['coverage_median']:.3f} [{b['k_range'][0]:.3f}–{b['k_range'][1]:.3f}] | {a['oxf_snr_median']:.2f} / {b['oxf_snr_median']:.2f} | {a['fig_matches'] and b['fig_matches']} |")
(X / "coil_contrast.md").write_text("\n".join(L) + "\n", encoding="utf-8")
print("\n".join(L)); print("identity all:", all(r[d]["identity"] for r in rows for d in ("OXF1PRI", "OXF2PRI")))
