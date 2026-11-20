#!/usr/bin/env python
"""Make the local Git history that the loader guard reads, for a copy unpacked from the Zenodo archive.

decomp/remeasure/guard.py checks with Git that the analysis code and the fixed document are committed and unchanged,
and that the fixed document was committed after the code. An unpacked archive has no Git history, so this script commits
the files in that order. The commits exist only on your machine; their dates are not evidence of when the rules were fixed.
"""
from __future__ import annotations

import importlib.util
import pathlib
import subprocess
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]


def git(*args: str) -> None:
    subprocess.run(["git", "-c", "user.name=local", "-c", "user.email=local@example.invalid", *args], cwd=REPO, check=True)


def main() -> int:
    if len(sys.argv) > 1:
        print(__doc__)
        print("Usage: python3 reproduce/make_local_git.py   (no options; run once in the unpacked folder)")
        return 0 if sys.argv[1] in ("-h", "--help") else 2
    if (REPO / ".git").exists():
        raise SystemExit("This folder already has a Git history; nothing to do.")
    spec = importlib.util.spec_from_file_location("guard", REPO / "decomp/remeasure/guard.py")
    guard = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(guard)
    git("init", "-q")
    git("add", "--", *guard.PIPELINE)
    git("commit", "-q", "-m", "Add analysis code and environment files")
    git("add", "--", guard.FIXED_DOC)
    git("commit", "-q", "-m", "Add the fixed analysis document")
    git("add", "-A")
    git("commit", "-q", "-m", "Add specifications, results, tables and reproduction scripts")
    print("Guard:", guard.check(REPO))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
