"""Finds the repository root: two levels above this file in the repository layout (decomp/rules_env), or the environment variable RULES_ENV_REPO in the development layout."""
import os
import pathlib
import sys

_HERE = pathlib.Path(__file__).resolve().parent


def repo_root() -> pathlib.Path:
    cand = _HERE.parents[1]
    if (cand / "decomp" / "stage1" / "config.py").is_file():
        return cand
    env = os.environ.get("RULES_ENV_REPO", "")
    if env and (pathlib.Path(env) / "decomp" / "stage1" / "config.py").is_file():
        return pathlib.Path(env)
    raise RuntimeError("저장소를 찾지 못함 — decomp/rules_env 배치이거나 RULES_ENV_REPO 를 지정")


def use(*parts: str) -> pathlib.Path:
    """Puts a directory inside the repository at the very front of the import path (e.g. use("decomp", "stage1"))."""
    p = repo_root().joinpath(*parts)
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
    return p
