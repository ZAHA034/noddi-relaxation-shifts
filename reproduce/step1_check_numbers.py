#!/usr/bin/env python
"""Step 1 of the reproduction: check the numbers of the main text against their source files.

ledger/numbers.csv lists every number of the main text; ledger/selectors.csv gives the file and field it comes from.
This script reads each source file, recomputes the field and compares it with the recorded value.
It also compares the SHA-256 of each source with the value recorded when the number was written,
and of each locked file with MANIFEST.csv. It needs nothing outside this repository.

Rows without a machine check are literature values and definitions. Working records are represented
by the lines or words the ledger cites, so their rows have no SHA-256 to compare.
"""
import csv
import hashlib
import pathlib
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "reproduce"))
import ledger_fields as LF  # noqa: E402

rows = list(csv.DictReader(open(REPO / "ledger" / "selectors.csv", encoding="utf-8")))
count = dict(checked=0, sha256_compared=0, no_machine_check=0, not_public=0)
bad = []
for r in rows:
    src, fld, sha = r["source_file"], r["source_field"], r["source_sha256"]
    if src == "not public":
        count["not_public"] += 1
        continue
    if not src or fld.startswith(("lit:", "none:")):
        count["no_machine_check"] += 1
        continue
    if sha:
        count["sha256_compared"] += 1
        if hashlib.sha256((REPO / src).read_bytes()).hexdigest() != sha:
            bad.append((r["id"], "source file differs from the recorded SHA-256"))
            continue
    # the field readers of ledger_fields.py use the Korean column names: source file, source field, SHA-256 of the source file, exact value
    row = {"근거 파일": src, "근거 필드": fld, "근거 파일 sha256": "", "원 값": r["raw_value"]}
    try:
        value = LF.eval_field(REPO, row)
        ok = LF.same_value(value, r["raw_value"])
    except Exception as e:
        ok, value = False, f"{type(e).__name__}: {e}"
    count["checked"] += 1
    if not ok:
        bad.append((r["id"], f"recorded {r['raw_value'][:60]!r}, source gives {str(value)[:60]!r}"))
locked = 0
for r in csv.DictReader(open(REPO / "MANIFEST.csv", encoding="utf-8")):
    if r["status"] == "identical":
        locked += 1
        if hashlib.sha256((REPO / r["path"]).read_bytes()).hexdigest() != r["sha256_here"]:
            bad.append((r["path"], "file differs from MANIFEST.csv"))
print(f"rows {len(rows)}:", count, f"· files identical to the development repository: {locked}")
for b in bad[:30]:
    print("  FAILED", *b)
print("RESULT", "all checked values agree" if not bad else f"{len(bad)} failures")
raise SystemExit(1 if bad else 0)
