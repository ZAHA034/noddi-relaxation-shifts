"""Counts whether the 73 of 109 first sessions outside the training range (15 – 60) lie below 15 or above 60, from stored snr_table.json."""
import datetime, hashlib, json, pathlib
P = pathlib.Path("<WORK>/remeasure/onharmony/v4/snr_table.json")
X = pathlib.Path(__file__).resolve().parent
b = P.read_bytes(); sha = hashlib.sha256(b).hexdigest()
now = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
with (X / "REREAD_LOG.jsonl").open("a", encoding="utf-8") as f:
    f.write(json.dumps(dict(시각=now, 다시읽음=str(P), sha256=sha, 바이트=len(b), 무엇="2-7 학습 범위 밖 세션의 위 · 아래 수", 스크립트=__file__), ensure_ascii=False) + "\n")
d = json.loads(b); lo, hi = d["train_range"]; s = d["sessions"]
rep = {f"sub-{a}_ses-OXF1PRI00{i}" for a in ("13192", "16974") for i in range(2, 7)}
first = {k: v for k, v in s.items() if k not in rep}
out = [v["snr"] for v in first.values() if not v["in_train"]]
below = sorted(x for x in out if x < lo); above = sorted(x for x in out if x > hi)
res = dict(snr_table_sha256=sha, made=now, n_sessions=len(s), n_first=len(first), n_out=len(out), n_below_15=len(below), n_above_60=len(above),
           below_range=[round(min(below), 2), round(max(below), 2)] if below else None, above_range=[round(min(above), 2), round(max(above), 2)] if above else None,
           first_in_train_range=[round(min(v["snr"] for v in first.values() if v["in_train"]), 2), round(max(v["snr"] for v in first.values() if v["in_train"]), 2)])
(X / "snr_range_counts.json").write_text(json.dumps(res, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
print(json.dumps(res, ensure_ascii=False))
