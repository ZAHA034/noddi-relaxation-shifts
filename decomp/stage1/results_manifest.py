"""Lists sha256 of committed results under decomp/results/** so stage1 analyses cannot overwrite them and tests can compare them before and after runs."""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C  # noqa: E402

ROOT = C.RESULTS.parents[1]


def tracked() -> list[str]:
    out = subprocess.run(["git", "ls-files", "--", "decomp/results"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.split()
    me = C.COMMITTED.relative_to(ROOT).as_posix()
    return sorted(p for p in out if p != me)


def current() -> dict[str, str]:
    return {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in tracked()}


def main() -> int:
    cur = current()
    if "--write" in sys.argv:
        doc = {"무엇": "커밋된 decomp 결과 파일의 sha256 — config.guard_write 가 덮어쓰기를 막고 tests/test_decomp_outputs.py 가 대조한다(중계 ZZZ-3 나)",
               "files": cur}
        C.COMMITTED.write_text(json.dumps(doc, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
        print("썼다:", C.COMMITTED.relative_to(ROOT), "·", len(cur), "파일")
        return 0
    old = json.loads(C.COMMITTED.read_text(encoding="utf-8"))["files"]
    bad = sorted(set(old) ^ set(cur)) + sorted(p for p in set(old) & set(cur) if old[p] != cur[p])
    print("일치" if not bad else f"⛔ 다름 {bad}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
