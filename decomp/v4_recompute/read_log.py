"""Write-ahead read log: the line that names a file is on disk before the file is read.

If the log cannot be written, ReadLogError stops the run and nothing is read. Each file gets two lines: one before it is
read, with its size, and one after, with its SHA-256.
"""
from __future__ import annotations

import datetime
import fcntl
import hashlib
import json
import os
import pathlib


class ReadLogError(RuntimeError):
    pass


class ReadLog:
    def __init__(self, path, what: str):
        self.path = pathlib.Path(path)
        self.what = what
        self.seen: set[str] = set()
        self._append([])                                            # stops here if the log cannot be written

    def _append(self, records: list[dict]) -> None:
        try:
            with open(self.path, "a", encoding="utf-8") as fh:
                fcntl.flock(fh, fcntl.LOCK_EX)
                fh.write("".join(json.dumps(r, ensure_ascii=False) + "\n" for r in records))
                fh.flush()
                os.fsync(fh.fileno())
        except OSError as e:
            raise ReadLogError(f"The read log {self.path} cannot be written, so nothing is read: {e}") from e

    def _base(self, p: pathlib.Path) -> dict:
        return {"시각": datetime.datetime.now().astimezone().isoformat(timespec="seconds"), "다시읽음": str(p), "무엇": self.what, "스크립트": "read_log.py"}

    def before_read(self, paths) -> None:
        """Write a line for every new path first, in one write. Only then hash each file and write its second line."""
        new = [pathlib.Path(p) for p in dict.fromkeys(str(x) for x in paths) if p not in self.seen]
        self._append([{**self._base(p), "단계": "before read", "바이트": p.stat().st_size} for p in new])
        self.seen |= {str(p) for p in new}
        for p in new:
            self._append([{**self._base(p), "단계": "after read", "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}])
