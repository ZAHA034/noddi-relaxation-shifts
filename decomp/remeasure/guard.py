"""새 실측(⑨-b) — 로더 가드와 열람 기록(R10_SEAL.md:21 · :24 · 중계 ZZZz-2 · 저자 결정 ZZZz-0 · ZZZZc-0). 해시를 박지 않는다 — git 의 상태로 본다.

영상을 열기 전에 `require()` 가 넷을 확인한다. 하나라도 어긋나면 GuardError 로 멈춘다.
  1. 고정 문서(FIXED_DOC)가 HEAD 에 추적되는가
  2. 고정 문서의 작업 사본이 HEAD 의 커밋판과 같은가
  3. 파이프라인 파일(PIPELINE)이 모두 HEAD 에 있고 작업 사본이 커밋판과 같은가 [A 덧붙임 — 커밋하지 않은 코드로 돌지 않게]
  4. 고정 문서를 마지막으로 바꾼 커밋이 파이프라인을 마지막으로 바꾼 커밋보다 뒤인가(같은 커밋은 뒤가 아니다) — ⑩ 뒤 코드를 바꾸면 고정 문서를 다시 커밋해야(이탈 기록) 열린다
열람 기록은 `log_access()` 로 정본 파일에 이어 적는다 — 세션마다 시각 · 근거 커밋 · 파일 · sha256(R10_SEAL.md:21).
"""
from __future__ import annotations

import datetime
import fcntl
import hashlib
import json
import os
import pathlib
import subprocess

HERE = pathlib.Path(__file__).resolve().parent
FIXED_DOC = "decomp/remeasure/MEASUREMENT_FIXED.md"
PIPELINE = ("decomp/remeasure/rconfig.py", "decomp/remeasure/guard.py", "decomp/remeasure/cmds.py",
            "decomp/remeasure/roi.py", "decomp/remeasure/fit_f2.py", "decomp/remeasure/run.py",
            "decomp/remeasure/analysis.py", "decomp/remeasure/sessions_onharmony.tsv",
            "decomp/fixed/fiso_sign.py", "decomp/fixed/delta_tau.py",                  # 예측(설계) — 중계 ZZZZf-2 ①
            "decomp/stage1/lib.py", "decomp/stage1/synth.py",                          # S4 와 같은 좌표 · 자름 · 확장 좌표
            "decomp/stage1/config.py", "decomp/stage1/analyze.py", "decomp/stage1/analyze2.py",
            "decomp/stage1/analyze3.py",                                                # 예측 값을 정하는 의존(카드 대수 · 장비 기록)
            "decomp/remeasure/judge_tong.py", "decomp/remeasure/judge_wong.py",
            "decomp/remeasure/judge_v4.py", "decomp/remeasure/judge_onh.py",            # 판정 코드(중계 ZZZZh-3 · ON-Harmony 쌍 검정 — ⑪ 사전 점검)
            "decomp/s4/validate.py", "decomp/s4/s4cfg.py", "decomp/s4/plan.py",         # V4 가 부르는 동결 판단기 코드
            "decomp/results/s4/train_freeze.json", "decomp/results/s4/sim_hashes.json",
            "decomp/results/s4/miss.json", "decomp/results/s4/cards_c3.json",
            "decomp/results/s4/cards_c7.json",                                          # 동결 기록 · 오지정 문턱 · 칸 분류(가중치 sha256 은 judge_v4.check_frozen 이 대조)
            "uv.lock", "pyproject.toml", "envs/reference/uv.lock",
            "envs/reference/pyproject.toml")                                            # 환경 고정(권고 ZZZZh-7 4 — 정확한 판 · 해시)
GUARD_EXIT = 77                       # 가드로 멈춘 하위 프로세스의 종료 코드 — 처리 실패와 가른다(⑪ 중 «영구 제외» 로 잘못 세지 않게)


class GuardError(RuntimeError):
    pass


def _git(root: pathlib.Path, *a: str, check: bool = True) -> str:
    r = subprocess.run(["git", *a], cwd=root, capture_output=True, text=True)
    if check and r.returncode != 0:
        raise GuardError(f"git {' '.join(a)} 실패: {r.stderr.strip()}")
    return r.stdout.strip()


def _same_as_head(root: pathlib.Path, rel: str) -> bool:
    head = subprocess.run(["git", "show", f"HEAD:{rel}"], cwd=root, capture_output=True)
    return head.returncode == 0 and (root / rel).is_file() and (root / rel).read_bytes() == head.stdout


def check(root: pathlib.Path = HERE.parents[1], fixed_doc: str = FIXED_DOC, pipeline: tuple[str, ...] = PIPELINE) -> dict:
    root = pathlib.Path(root)
    if not _git(root, "ls-tree", "HEAD", "--", fixed_doc):
        raise GuardError(f"⛔ 고정 문서 {fixed_doc} 가 HEAD 에 없다 — ⑩ 전에는 영상을 열지 않는다(#51 · R10)")
    if not _same_as_head(root, fixed_doc):
        raise GuardError(f"⛔ 고정 문서 {fixed_doc} 의 작업 사본이 커밋판과 다르다")
    for rel in pipeline:
        if not _git(root, "ls-tree", "HEAD", "--", rel):
            raise GuardError(f"⛔ 파이프라인 파일 {rel} 이 HEAD 에 없다")
        if not _same_as_head(root, rel):
            raise GuardError(f"⛔ 파이프라인 파일 {rel} 의 작업 사본이 커밋판과 다르다 — 커밋하지 않은 코드로 돌지 않는다")
    fixed_c = _git(root, "log", "-1", "--format=%H", "--", fixed_doc)
    pipe_c = _git(root, "log", "-1", "--format=%H", "--", *pipeline)
    after = fixed_c != pipe_c and subprocess.run(["git", "merge-base", "--is-ancestor", pipe_c, fixed_c], cwd=root).returncode == 0
    if not after:
        raise GuardError(f"⛔ 고정 문서 커밋({fixed_c[:7]})이 파이프라인 커밋({pipe_c[:7]})보다 뒤가 아니다 — 코드를 바꿨으면 이탈을 적고 고정 문서를 다시 커밋한다")
    return dict(head=_git(root, "rev-parse", "HEAD"), fixed_commit=fixed_c, pipeline_commit=pipe_c)


def require(root: pathlib.Path = HERE.parents[1]) -> dict:
    """영상을 열기 전에 부른다 — 통과하면 근거(커밋 셋)를 돌려준다."""
    return check(root)


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def log_access(log_path: pathlib.Path, key: tuple[str, ...], entry: dict, lock_path: pathlib.Path) -> None:
    """정본 열람 기록(JSON)의 key 경로 끝 목록에 entry 를 잠금 안에서 덧붙인다 — 동시 작업(동시 3)이 겹쳐 쓰지 않게."""
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with open(lock_path, "a") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        try:
            d = json.loads(log_path.read_text(encoding="utf-8"))
            node = d
            for k in key[:-1]:
                node = node.setdefault(k, {})
            node.setdefault(key[-1], []).append(entry)
            tmp = log_path.with_suffix(".tmp")
            tmp.write_text(json.dumps(d, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
            tmp.replace(log_path)
        finally:
            fcntl.flock(lk, fcntl.LOCK_UN)


def access_entry(what: str, files: list[pathlib.Path], base: pathlib.Path, basis: dict, session: str) -> dict:
    """열람 항목 하나 — 시각(UTC) · 세션 · 근거(#51 · 고정 문서 커밋) · 파일(base = 저장소 맨 위 기준 상대 경로 «../…») · sha256."""
    return {"utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"), "세션": session, "무엇": what,
            "근거": f"#51 · 고정 문서 {basis['fixed_commit'][:12]} · 파이프라인 {basis['pipeline_commit'][:12]}",
            "파일": [{"file": os.path.relpath(f, base), "sha256": sha256(pathlib.Path(f))} for f in files]}
