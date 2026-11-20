"""Compares Wong TR with the published Methods by reading only pixdim[4] and time units from the NIfTI headers of the sub-01 8 TE series."""
import datetime, gzip, json, os, pathlib, struct
W = pathlib.Path("<WORK>/external/wongkornchaovalit2025")
REG = pathlib.Path("<REPO>/prereg/artifacts/decomp/wongkornchaovalit2025_sha256.txt")
X = pathlib.Path(__file__).resolve().parent
METHODS_TR_MS = 5800.0
reg = {}
for ln in REG.read_text(encoding="utf-8").splitlines():
    s, b, p = ln.split(None, 2); reg[p] = (s, int(b))
now = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
rows = []
for te in (62, 72, 82, 92, 102, 112, 122, 132):
    rel = f"MTE-dMRI/sub-01/dwi/sub-01_acq-TE{te}_dwi.nii.gz"
    p = W / rel
    s, nbytes = reg[rel]
    assert p.stat().st_size == nbytes, (rel, p.stat().st_size, nbytes)
    with gzip.open(p, "rb") as g:
        h = g.read(348)
    assert len(h) == 348
    e = "<" if struct.unpack("<i", h[:4])[0] == 348 else ">"
    assert struct.unpack(e + "i", h[:4])[0] == 348 and h[344:347] in (b"n+1", b"ni1"), "NIfTI-1 머리가 아님"
    dim = struct.unpack(e + "8h", h[40:56]); pixdim = struct.unpack(e + "8f", h[76:108]); units = h[123]
    tunit = {8: "s", 16: "ms", 24: "us", 0: "미상"}.get(units & 0x38, f"코드 {units & 0x38}")
    tr_ms = pixdim[4] * {"s": 1000.0, "ms": 1.0, "us": 0.001}.get(tunit, float("nan"))
    rows.append(dict(series=f"TE{te}", file=rel, registered_sha256=s, bytes=nbytes, header_bytes_read=348, dim4=dim[4], pixdim4=pixdim[4], time_unit=tunit, TR_ms_from_header=tr_ms))
    with (X / "EXTERNAL_READ_LOG.jsonl").open("a", encoding="utf-8") as f:
        f.write(json.dumps(dict(시각=now, 자료="Wongkornchaovalit 2025", 읽음=str(p), 읽은_것="NIfTI-1 머리 348 B(pixdim[4] · 시간 단위 · dim[4]) — 화소 0",
                                sha256_등록값=s, 바이트=nbytes, 무엇="Wong TR 대조(D3 전 · 중계 2026-11-01 결정)", 스크립트=__file__), ensure_ascii=False) + "\n")
info = [r["pixdim4"] for r in rows]
uninformative = all(v in (0.0, 1.0) for v in info)
same = len({round(r["TR_ms_from_header"], 3) for r in rows}) == 1
verdict = ("머리 정보 없음(0 · 1) — Methods 만 씀" if uninformative else
           ("머리 TR 이 계열마다 같고 Methods 와 같음" if same and abs(rows[0]["TR_ms_from_header"] - METHODS_TR_MS) < 1.0 else "⛔ 다름 — 멈추고 보고"))
out = dict(made=now, methods=dict(TR_ms=METHODS_TR_MS, delta_ms="15.2/25.2", source="decomp/stage1/analyze3.py:339–340 (Wongkornchaovalit et al. Sci Data 2025 · PMC12246157)"),
           header_rows=rows, header_uninformative=uninformative, header_all_same=same, result=verdict)
(X / "wong_tr_header.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
for r in rows:
    print(r["series"], "dim4", r["dim4"], "pixdim4", round(r["pixdim4"], 4), r["time_unit"], "→ TR", round(r["TR_ms_from_header"], 1), "ms")
print("결과:", verdict)
