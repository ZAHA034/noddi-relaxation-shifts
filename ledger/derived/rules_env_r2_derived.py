"""Derived ledger values from the second set of rules_env outputs (s3.json, e6.json): one maximum and one interval, no new computation; logs each result json before reading it."""
import hashlib
import json
import pathlib
import subprocess
import sys

sys.dont_write_bytecode = True
HERE = pathlib.Path(__file__).resolve().parent
LEDGER = HERE.parent
OUT = LEDGER.parent / "env" / "out"


def logged(name):
    p = OUT / f"{name}.json"
    r = subprocess.run([sys.executable, "-I", "-B", str(LEDGER / "reread_log.py"), str(p), "장부 도출 값(rules_env_r2_derived.py) — R2 결과 합치기 · 새 계산 없음"],
                       capture_output=True, text=True)
    h = r.stdout.strip().splitlines()[-1] if r.returncode == 0 and r.stdout.strip() else None
    raw = p.read_bytes()
    if h != hashlib.sha256(raw).hexdigest():
        raise SystemExit(f"STOP: {p} 기록 실패 또는 바뀜")
    return json.loads(raw), h


s3, h = logged("s3")
e6, h6 = logged("e6")
rows30 = {k: v for k, v in s3["rows"].items() if v["te_a"] == 100.0}
kmax = max(rows30, key=lambda k: rows30[k]["raw"])
res = dict(sources={"s3.json": h, "e6.json": h6}, s3=dict(max_raw_dte30=rows30[kmax]["raw"], max_raw_dte30_cell=kmax, n_cells_dte30=len(rows30)),
           e6=dict(diff_interval=[e6["primary"]["diff_q025"], e6["primary"]["diff_q975"]]))
(HERE / "rules_env_r2_derived.json").write_text(json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")
print(json.dumps(res, ensure_ascii=False))
