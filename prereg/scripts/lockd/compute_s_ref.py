"""Computes and freezes s_ref, the reference retest SD s of the stratum 1 interval, from OXF1PRI repeats of 2 people × 6 sessions."""
from __future__ import annotations

import hashlib
import json
import pathlib
import re
import statistics
import sys

from scipy.stats import chi2

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import _harness as H  # noqa: E402
from dmri_cov.utils import ext  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = ROOT / "prereg" / "artifacts" / "unseal" / "s_ref.frozen.json"
REF = "OXF1PRI"


def main() -> int:
    roi_dir = ext.locate("data", "roi")
    files = sorted(roi_dir.glob(f"*_ses-{REF}*.roi.tsv"))
    units, order = H._label_units(), H.roi_names()
    by_subj: dict[str, list[list[list[float]]]] = {}
    manifest = []
    for f in files:
        sub = re.match(r"(sub-[^_]+)_ses-", f.name).group(1)
        rows = [[float(x) for x in l.split()] for l in f.read_text().splitlines() if l.strip()]
        assert len(rows) == 48, f.name
        by_subj.setdefault(sub, []).append(rows)
        manifest.append({"file": f.name, "sha256": hashlib.sha256(f.read_bytes()).hexdigest()})
    df = sum(len(v) - 1 for v in by_subj.values())
    assert len(by_subj) == 2 and df == 10, (len(by_subj), df)
    ci = [(df / chi2.ppf(0.975, df)) ** 0.5, (df / chi2.ppf(0.025, df)) ** 0.5]
    s_ref: dict[str, dict[str, float]] = {}
    for ax in H.AXES:
        col = H._ROI_TSV_COL[ax]
        per_unit: dict[str, list[float]] = {u: [] for u in order}
        for sess in by_subj.values():
            vals = {u: [] for u in order}
            for rows in sess:
                acc: dict[str, list[float]] = {}
                for lab in range(48):
                    acc.setdefault(units[lab], []).append(rows[lab][col])
                for u in order:
                    vals[u].append(statistics.mean(acc[u]))
            for u in order:
                per_unit[u].append(statistics.variance(vals[u]))
        s_ref[ax] = {u: statistics.mean(v) ** 0.5 for u, v in per_unit.items()}
    doc = {
        "문서": "s_ref — 층 1 구간의 참조 재검사 SD(LOCK_B §143-2 · §110). ⛔ roi_s.frozen 의 판정선 s_ROI 와 다른 양",
        "입력": {"참조": REF, "피험자": sorted(by_subj), "세션 수": {k: len(v) for k, v in by_subj.items()},
               "파일": manifest, "스캐너 간 값": "쓰지 않음(입력 디렉터리에 없다, P9-9)"},
        "자유도": df,
        "s 의 95% 구간 배수(자유도만)": {"하한": round(ci[0], 4), "상한": round(ci[1], 4)},
        "⛔ 얇다": "df 10 — roi_s(df 45) 보다 얇다. 이 불확실성은 G5 폭 분모 z·s_ref 와 층 1 v 의 뺄셈 두 곳에 동시에 들어간다(저자 «승인 넷» 2 ①)",
        "단위": "원래 단위(로짓 아님) · 축 = Frame.AXES · ODI 는 판정 대상이 아니다(기술)",
        "ROI별 s_ref": s_ref,
    }
    OUT.write_text(json.dumps(doc, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    for ax in H.AXES:
        v = list(s_ref[ax].values())
        print(f"{ax}: 27 ROI s_ref 범위 [{min(v):.5f}, {max(v):.5f}]")
    print(f"df {df} · 95% 배수 [{ci[0]:.3f}, {ci[1]:.3f}] → {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
