"""Builds the per-cell sensitivity table from the stored V4 output v4_c3.json, with no new computation or image reading; post hoc exploratory, not a verdict."""
import hashlib, json, pathlib, datetime, collections
V4 = pathlib.Path("<WORK>/remeasure/onharmony/v4/v4_c3.json")
X = pathlib.Path(__file__).resolve().parent
LOGF = X / "REREAD_LOG.jsonl"
EXPECT = "f6003589"                      # Access record · first 8 characters written on the receipt
b = V4.read_bytes(); sha = hashlib.sha256(b).hexdigest()
assert sha.startswith(EXPECT), sha
now = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
with LOGF.open("a", encoding="utf-8") as f:
    f.write(json.dumps(dict(시각=now, 다시읽음=str(V4), sha256=sha, 바이트=len(b), 무엇=b"\xe2\x91\xa1".decode() + " 민감도 표(저장된 출력)", 스크립트=__file__), ensure_ascii=False) + "\n")
d = json.loads(b)
cells = d["main"]["cells"]
COND = ["SNR20|Neff25", "SNR20|Neff100", "SNR50|Neff25", "SNR50|Neff100"]
rows, cnt = [], {c: collections.Counter() for c in COND}
for k, v in cells.items():
    if v.get("verdict") == "자료 없음":
        continue
    sp = v.get("sensitivity_passed") or {}
    r = dict(cell=k, n_units=v["n_units"], n_used=v["n_units_used"], n_excl_misspec=v["n_excluded_misspec"], n_drop_k=v["n_units_dropped_k"],
             n_subjects=v["n_subjects"], few=v["too_few_subjects"], coverage=v["coverage"], threshold=v["threshold"], verdict=v["verdict"],
             **{c: sp.get(c, "없음") for c in COND})
    rows.append(r)
    for c in COND:
        cnt[c][{True: "통과", False: "실패", None: "쓴 단위 0"}.get(sp.get(c), "없음")] += 1
lab = lambda s: s.replace("|", " · ")
fmt = lambda x: "—" if x is None else ("통과" if x is True else "실패" if x is False else x)
L = [f"# V4 민감도 표 — 저장된 출력만(사후 탐색 · 판정 아님) · {now}", "",
     f"- 입력: `v4_c3.json` sha256 {sha[:16]}… · {len(b):,} B(다시 읽음 — `drafts/x14/REREAD_LOG.jsonl`) · 새 계산 0 · 새 영상 읽기 0.",
     "- 주 칸 = SNR 30 · N_eff 50. 민감도 조건의 문턱은 주 칸 것 그대로(judge_v4 A-1 2 A 읽기).",
     "- **민감도 조건마다의 쓴 단위 · 오지정 제외 · 덮음은 저장돼 있지 않다** — 판정 코드가 통과 표시만 남겼다. «쓴 단위 0» = 표시가 None(그 조건에서 오지정 · k 미충족으로 모두 빠짐 — 피험자 수는 주 칸과 같음).", "",
     "## 조건별 셈(자료 있는 칸 24)", "", "| 조건 | 통과 | 실패 | 쓴 단위 0 |", "|---|---|---|---|"]
L += [f"| {lab(c)} | {cnt[c]['통과']} | {cnt[c]['실패']} | {cnt[c]['쓴 단위 0']} |" for c in COND]
L += ["", "## 칸마다(주 칸 수 + 민감도 표시)", "",
      "| 칸 | 단위 | 쓴 단위 | 오지정 제외 | k 미충족 제외 | 피험자 | 덮음 | 문턱 | 판정 | " + " | ".join(lab(c) for c in COND) + " |",
      "|---|" + "---|" * (8 + len(COND))]
for r in sorted(rows, key=lambda r: r["cell"]):
    cov = "—" if r["coverage"] is None else f"{r['coverage']:.3f}"
    L.append(f"| {lab(r['cell'])} | {r['n_units']} | {r['n_used']} | {r['n_excl_misspec']} | {r['n_drop_k']} | {r['n_subjects']} | {cov} | {r['threshold']:.2f} | {r['verdict']} | "
             + " | ".join("쓴 단위 0" if r[c] is None else fmt(r[c]) for c in COND) + " |")
(X / "v4_sensitivity_table_stored.md").write_text("\n".join(L) + "\n", encoding="utf-8")
(X / "v4_sensitivity_table_stored.json").write_text(json.dumps(dict(input_sha256=sha, made=now, rows=rows, counts={c: dict(cnt[c]) for c in COND}), ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
print(json.dumps({c: dict(cnt[c]) for c in COND}, ensure_ascii=False))
