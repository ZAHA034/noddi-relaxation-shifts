"""Step 4a: extracts only ReceiveCoilName, CoilCombinationMethod and ParallelReductionFactorInPlane from dwi sidecars of the 109 V4 first sessions; exploratory, descriptive."""
import collections, datetime, hashlib, json, os, pathlib, re, sys
REPO = pathlib.Path("<REPO>")
sys.path.insert(0, str(REPO / "decomp" / "remeasure"))
import rconfig as R  # noqa: E402
import guard as G  # noqa: E402
X14 = pathlib.Path("<WORK>/stage1_work/drafts/x14")
OUT = R.WORK / "explore" / "rv"
KEYS = ("ReceiveCoilName", "CoilCombinationMethod", "ParallelReductionFactorInPlane")
CH = (8, 12, 16, 20, 32, 48, 64)
NOW = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
SHA_ATTR = "5184137854abc69710904889503eaf8dbb111ca9736e46af13b49067ea65b90d"


def channels(name):
    if not isinstance(name, str):
        return None
    for m in re.findall(r"\d+", name):
        if int(m) in CH:
            return int(m)
    return None


def main():
    basis = G.require()
    ap = R.WORK / "onharmony" / "v4" / "attr_table.json"
    b = ap.read_bytes()
    assert hashlib.sha256(b).hexdigest() == SHA_ATTR
    with (X14 / "REREAD_LOG.jsonl").open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(dict(시각=NOW, 다시읽음=str(ap), sha256=SHA_ATTR, 바이트=len(b), 무엇="4a — 첫 세션 목록", 스크립트=__file__), ensure_ascii=False) + "\n")
    at = json.loads(b)
    sess = sorted({(sub, s) for d in at["pairs"].values() for sub, x in d.items() for s in x["sessions"]})
    rows, logs = {}, []
    for sub, ses in sess:
        p = R.ONH_META / sub / ses / "dwi" / f"{sub}_{ses}_dir-AP_dwi.json"
        raw = p.read_bytes()
        meta = json.loads(raw)
        got = {k: meta.get(k, "칸 없음") for k in KEYS}
        del meta
        logs.append(dict(시각=NOW, utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"), 세션=f"{sub}_{ses}",
                         파일=os.path.relpath(p, REPO), sha256=hashlib.sha256(raw).hexdigest(), 바이트=len(raw), 무엇="4a 사이드카 세 칸(ReceiveCoilName · CoilCombinationMethod · ParallelReductionFactorInPlane)",
                         근거=basis, 스크립트=__file__))
        dev = R.device_of(ses)
        rows[f"{sub}_{ses}"] = dict(device=dev, unit=(R.not2ing_protocol(sub) if dev == "NOT2ING" else dev), **got,
                                    channels=channels(got["ReceiveCoilName"]))
    with (X14 / "ONH_NEW_READ_LOG.jsonl").open("a", encoding="utf-8") as fh:
        for e in logs:
            fh.write(json.dumps(e, ensure_ascii=False) + "\n")
    summ = {}
    for dev in sorted({v["device"] for v in rows.values()}):
        vs = [v for v in rows.values() if v["device"] == dev]
        summ[dev] = {k: dict(collections.Counter(str(v[k]) for v in vs)) for k in KEYS + ("channels",)}
        summ[dev]["n_sessions"] = len(vs)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "rv_4a_coil.json").write_text(json.dumps(dict(made=NOW, script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
                                                         n_files=len(rows), keys=KEYS, by_device=summ, sessions=rows), ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    for dev, s in summ.items():
        print(dev, s["n_sessions"], json.dumps({k: s[k] for k in KEYS + ("channels",)}, ensure_ascii=False))
    print("읽은 파일", len(logs))


if __name__ == "__main__":
    main()
