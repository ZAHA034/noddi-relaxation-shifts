"""Computes d2 observed six-session SDs, session-pair difference SDs, centers and subject SNR from the 12 OXF1PRI repeat sessions, per subject, ROI and coordinate."""
import datetime, hashlib, json, math, pathlib, sys
import numpy as np
REPO = pathlib.Path("<REPO>")
sys.path.insert(0, str(REPO / "decomp" / "remeasure"))
import judge_v4 as JV  # noqa: E402 — session_x (same coordinates as the judge) · torch is not called at the module head
W = pathlib.Path("<WORK>/remeasure")
OUT = W / "explore" / "d2"
LOG = pathlib.Path("<WORK>/stage1_work/drafts/x14/REREAD_LOG.jsonl")
SUBS, SES = ("sub-13192", "sub-16974"), [f"ses-OXF1PRI00{i}" for i in range(1, 7)]
now = datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def reread(p, what):
    b = p.read_bytes()
    s = hashlib.sha256(b).hexdigest()
    with LOG.open("a", encoding="utf-8") as f:
        f.write(json.dumps(dict(시각=now, 다시읽음=str(p), sha256=s, 바이트=len(b), 수정시각=datetime.datetime.fromtimestamp(p.stat().st_mtime).astimezone().isoformat(timespec="seconds"),
                                무엇=what, 스크립트=__file__), ensure_ascii=False) + "\n")
    return b, s


snr_b, snr_sha = reread(W / "onharmony" / "v4" / "snr_table.json", "D2 — 피험자 SNR")
SNR = json.loads(snr_b)["sessions"]
out = dict(what=__doc__.splitlines()[0], made=now, script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
           snr_table_sha256=snr_sha, subjects={})
for sub in SUBS:
    x3, roi_rows, files = [], [], {}
    for ses in SES:
        p = W / "onharmony" / f"{sub}_{ses}" / "roi.json"
        b, s = reread(p, "D2 — OXF1PRI 반복 흩어짐(관측)")
        files[ses] = s
        roi = json.loads(b)
        units = list(roi["주"])
        x7 = JV.session_x(roi, units)                      # Judge function as is (C7) — only ROI units that have all seven maps
        x3.append({u: v[:3].tolist() for u, v in x7.items()})
        roi_rows.append(roi["주"])
    common = sorted(set.intersection(*[set(d) for d in x3]))
    A = np.array([[d[u] for u in common] for d in x3])      # (6 sessions, ROI units, 3)
    sd = A.std(0, ddof=1)
    snrs = [SNR[f"{sub}_{ses}"]["snr"] for ses in SES]
    out["subjects"][sub] = dict(
        units=common, n_units=len(common), roi_sha256=files,
        centre_ses001={u: [roi_rows[0][u]["NDI"], roi_rows[0][u]["FWF"], roi_rows[0][u]["ODI"]] for u in common},
        x3_by_session={ses: [d[u] for u in common] for ses, d in zip(SES, x3)},
        sd=sd.tolist(), diff_sd=(math.sqrt(2) * sd).tolist(),
        snr_sessions=snrs, snr_subject=math.sqrt(len(snrs) / sum(1 / s ** 2 for s in snrs)))
    print(sub, "단위", len(common), "SNR_subj", round(out["subjects"][sub]["snr_subject"], 3), "차이 SD 중앙(x,y,ODI)", np.round(np.median(math.sqrt(2) * sd, 0), 4))
(OUT / "d2_observed.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
print("→", OUT / "d2_observed.json")
