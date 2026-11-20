"""새 실측(⑨-a) — 세션 처리 구동(ON-Harmony · Tong · Wongkornchaovalit). CPU 만(저자 결정 ZZZZd-1 (B)).

실자료 경로(dmri_data · dmri_work/external · 기본 WORK)를 만지면 먼저 guard.require()(⑩ 뒤에만 통과) · 열 때마다 정본 열람 기록에 이어 적는다.
가짜 자료(그 밖의 경로)는 가드 · 기록 없이 돈다 — 정답 시험과 연기 시험용.

    uv run --frozen python decomp/remeasure/run.py onharmony <sub> <ses>        # 한 세션(하위 세션을 띄우고 시도를 적는다 · 일괄과 같은 규칙)
    uv run --frozen python decomp/remeasure/run.py onharmony-batch [동시 수]    # 첫 세션 110 + 반복 10 · 동시 3 · 상한 벽시계 5 일
    uv run --frozen python decomp/remeasure/run.py tong-extract <sub>          # 분할 zip → 세션마다 dwi(제공자 전처리본)만
    uv run --frozen python decomp/remeasure/run.py tong <sub> <ses>
    uv run --frozen python decomp/remeasure/run.py wong <sub>

종료 코드 — 처리 실패로 세지 않는 멈춤 셋(A 제안 a12): 77 가드 멈춤 · 75 여유 디스크 부족 · 76 잠김(같은 작업 폴더나 같은 세션을 다른 프로세스가 돌리는 중).
"""
from __future__ import annotations

import contextlib
import datetime
import fcntl
import json
import os
import pathlib
import shutil
import subprocess
import sys
import time
import traceback

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cmds as K  # noqa: E402
import guard as G  # noqa: E402
import rconfig as R  # noqa: E402

REAL_ROOTS = (R.ROOT / "dmri_data", R.ROOT / "dmri_work" / "external", R.ROOT / "dmri_work" / "remeasure")
LOCK = R.ROOT / "dmri_work" / "remeasure_locks" / "access_log.lock"     # 저장소 밖 — 동시 작업의 기록 잠금
ENV = dict(os.environ, CUDA_VISIBLE_DEVICES="", OMP_NUM_THREADS=str(R.EDDY_NTHR), FSLDIR=str(R.FSLDIR), FSLOUTPUTTYPE="NIFTI_GZ",
           PATH=f"{R.FSLDIR / 'bin'}:{os.environ.get('PATH', '')}")                # bet 이 bet2 를 PATH 에서 찾는다
STOP_EXIT = 75                                                               # 여유 디스크 부족 — 멈춤(옛 경로 timing_run.sh 의 «중단» 갈래 · 재점검 RUN-1)
BUSY_EXIT = 76                                                               # 잠김 — 다른 프로세스가 같은 작업 폴더 · 같은 세션을 돌리는 중(재점검 RUN-4)
NOT_COUNTED = (0, G.GUARD_EXIT, STOP_EXIT, BUSY_EXIT)                        # 처리 실패로 세지 않는 종료 코드
WORKER = "REMEASURE_WORKER"                                                  # 일괄 · 한 세션 명령이 띄운 하위 세션의 표지 — 시도는 띄운 쪽이 적는다


class StopRun(RuntimeError):
    """처리 실패가 아닌 멈춤 — 여유 디스크 부족(종료 75). 세지 않고 멈춘다 · 디스크를 비운 뒤 같은 명령으로 다시."""


class Busy(RuntimeError):
    """잠금이 잡혀 있다(종료 76) — 처리 실패가 아니다 · 그것이 끝난 뒤 같은 명령으로 다시."""


_HELD: list[int] = []                                                        # 세션 잠금 fd — 단계의 하위 프로세스(eddy 등)에 넘긴다


@contextlib.contextmanager
def flock_or_busy(path: pathlib.Path, what: str, hand_down: bool = False):
    """잠금(fcntl.flock · 기다리지 않음 — 잡혀 있으면 Busy). 잠금 파일은 세션 폴더 밖(정리 단계가 지우지 않음).
    hand_down = 단계의 하위 프로세스에 fd 를 넘긴다 — 파이썬이 먼저 죽어도 손자(eddy)가 끝날 때까지 잠금이 남는다."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fh = open(path, "a")
    try:
        fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        fh.close()
        raise Busy(f"⛔ {what} — 잠금 {path} 이 잡혀 있다 · 그것이 끝난 뒤 같은 명령으로 다시(처리 실패가 아니다)") from None
    if hand_down:
        _HELD.append(fh.fileno())
    try:
        yield
    finally:
        if hand_down:
            _HELD.remove(fh.fileno())
        fh.close()


def _now() -> float:
    return time.time()


def _sleep(s: float) -> None:
    time.sleep(s)


def _utc() -> str:
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


# ── 가드 · 기록 ──────────────────────────────────────────────────────────────
def is_real(p: pathlib.Path) -> bool:
    p = pathlib.Path(os.path.abspath(p))
    return any(p == r or r in p.parents for r in REAL_ROOTS)


def guard_if_real(paths) -> dict | None:
    """실자료 경로가 하나라도 있으면 가드를 통과해야 한다 — 통과하면 근거, 가짜 자료면 None."""
    return G.require() if any(is_real(p) for p in paths) else None


def log(kind: str, what: str, files, basis: dict | None, session: str) -> None:
    if basis is None:
        return
    path, key = {"onharmony": (R.ACCESS_LOG_ONH, ("열람",)),
                 "tong": (R.ACCESS_LOG_EXT, ("Tong 2020", "열람(새 분석)")),
                 "wong": (R.ACCESS_LOG_EXT, ("Wongkornchaovalit 2025", "열람(새 분석)"))}[kind]
    G.log_access(path, key, G.access_entry(what, [pathlib.Path(f) for f in files], R.REPO, basis, session), LOCK)


# ── 단계 실행 ────────────────────────────────────────────────────────────────
def free_gb(p: pathlib.Path) -> float:
    return shutil.disk_usage(p).free / 2 ** 30


def low_disk(p: pathlib.Path) -> bool:
    """p(없으면 있는 가장 가까운 위 폴더)의 여유 디스크가 문턱(옛 경로와 같은 20 GB) 아래인가."""
    p = pathlib.Path(os.path.abspath(p))
    while not p.exists() and p != p.parent:
        p = p.parent
    return free_gb(p) < R.MINFREE_GB


def step(wd: pathlib.Path, name: str, argv: list[str]) -> None:
    """한 단계 — nice · 로그 · 시간 · 여유 디스크 확인(옛 경로와 같은 20 GB · 모자라면 처리 실패가 아니라 멈춤 — StopRun · 재점검 RUN-1)."""
    if free_gb(wd) < R.MINFREE_GB:
        raise StopRun(f"⛔ 여유 디스크 {free_gb(wd):.0f} GB < {R.MINFREE_GB} GB — 멈춘다(처리 실패가 아니다 · 디스크를 비운 뒤 같은 명령으로 다시)")
    t0 = time.time()
    with open(wd / f"log_{name}.txt", "w") as fh:
        r = subprocess.run(["nice", "-n", str(R.NICE), *argv], cwd=wd, stdout=fh, stderr=subprocess.STDOUT, env=ENV, pass_fds=tuple(_HELD))
    with open(wd / "timings.tsv", "a") as fh:
        fh.write(f"{name}\t{time.time() - t0:.1f}\t{r.returncode}\n")
    if r.returncode == G.GUARD_EXIT:
        raise G.GuardError(f"⛔ [{name}] 가드로 멈춤 — 처리 실패가 아니다(가드 상태를 되돌린 뒤 같은 명령으로 다시) · {wd / f'log_{name}.txt'}")
    if r.returncode != 0:
        if free_gb(wd) < R.MINFREE_GB:                                          # 단계 도중 디스크가 찼다 — 세지 않는다
            raise StopRun(f"⛔ [{name}] 종료 {r.returncode} · 여유 디스크 {free_gb(wd):.0f} GB < {R.MINFREE_GB} GB — 디스크 부족으로 멈춘다(처리 실패가 아니다)")
        raise RuntimeError(f"⛔ [{name}] 실패(종료 {r.returncode}) — {wd / f'log_{name}.txt'}")


def f2(wd: pathlib.Path, dwi: str, bvals: str, bvecs: str, mask: str, low: float, prefix: str) -> None:
    """F2 · 확장 좌표(envs/reference)."""
    step(wd, f"f2_{prefix}", ["uv", "run", "--frozen", "--no-sync", "--project", str(R.REPO / "envs" / "reference"),
                              "python", str(HERE / "fit_f2.py"), str(wd), dwi, bvals, bvecs, mask, str(low), prefix, str(R.EDDY_NTHR)])


def to_std_and_roi(wd: pathlib.Path, prefix: str, warp: pathlib.Path, mask: str, out: str) -> dict:
    """원 공간 지도 일곱 + 뇌 마스크(nn) + 적합 유효 지도 둘(trilinear) → 표준 공간 → 골격 ROI 27(주 · 민감도) → <out>.json · 표준 공간 지도는 지운다."""
    import roi  # noqa: E402
    std = {}
    for f in R.FEATURES:
        step(wd, f"warp_{prefix}_{f}", K.applywarp(f"{prefix}_{f}", f"std_{prefix}_{f}", str(warp)))
        std[f] = wd / f"std_{prefix}_{f}.nii.gz"
    step(wd, f"warp_{prefix}_mask", K.applywarp(mask, f"std_{prefix}_mask", str(warp), interp="nn"))
    for v in ("extvalid", "f2valid"):
        step(wd, f"warp_{prefix}_{v}", K.applywarp(f"{prefix}_{v}", f"std_{prefix}_{v}", str(warp)))
    aux = [wd / f"std_{prefix}_{v}.nii.gz" for v in ("mask", "extvalid", "f2valid")]
    res = roi.extract(std, *aux)
    (wd / f"{out}.json").write_text(json.dumps(res, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    for p in list(std.values()) + aux:
        p.unlink(missing_ok=True)
    return res


def done(wd: pathlib.Path, info: dict) -> None:
    info["끝_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    (wd / "DONE.json").write_text(json.dumps(info, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")


# ── ON-Harmony ───────────────────────────────────────────────────────────────
def onh_inputs(sub: str, ses: str, sealed: pathlib.Path = R.ONH_SEALED, meta: pathlib.Path = R.ONH_META) -> dict:
    s, m = sealed / sub / ses / "dwi", meta / sub / ses / "dwi"
    b = f"{sub}_{ses}_dir"
    return dict(ap=s / f"{b}-AP_dwi.nii.gz", pa=s / f"{b}-PA_dwi.nii.gz", bval=m / f"{b}-AP_dwi.bval", bvec=m / f"{b}-AP_dwi.bvec",
                json_ap=m / f"{b}-AP_dwi.json", json_pa=m / f"{b}-PA_dwi.json")


KEEP_ONH = ("acqparams.txt", "index.txt", "bvals", "bvecs", "nodif_brain_mask.nii.gz", "dti_FA.nii.gz", "fa2std.mat",
            "fa2std_warp.nii.gz", "roi.json", "snr.json", "timings.tsv", "DONE.json")


def prune_onh(wd: pathlib.Path) -> None:
    """세션 중간물 정리 — 원 공간 F2 · 확장 지도 일곱 · 마스크 · FA · 변환 · eddy 품질 글(영상 아님) · 로그만 남긴다(prune_session.sh 의 보존 원칙)."""
    for p in wd.iterdir():
        keep = (p.name in KEEP_ONH or p.name.startswith("f2_") or p.name.startswith("log_")
                or (p.name.startswith("eddy_out.eddy_") and not p.name.endswith(".nii.gz")))
        if not keep:
            shutil.rmtree(p) if p.is_dir() and not p.is_symlink() else p.unlink()


def onh_session(sub: str, ses: str, work: pathlib.Path = R.WORK, inputs: dict | None = None) -> dict:
    """한 세션 — 세션 잠금 안에서(같은 세션을 두 프로세스가 함께 돌리지 않는다 · 잠금은 단계의 하위 프로세스에도 넘긴다 · 재점검 RUN-4)."""
    with flock_or_busy(work / "onharmony" / "locks" / f"{sub}_{ses}.lock", f"{sub}_{ses} 를 다른 프로세스가 돌리는 중", hand_down=True):
        return _onh_session(sub, ses, work, inputs)


def _onh_session(sub: str, ses: str, work: pathlib.Path, inputs: dict | None) -> dict:
    inp = inputs or onh_inputs(sub, ses)
    wd = work / "onharmony" / f"{sub}_{ses}"
    if (wd / "DONE.json").exists():
        return json.loads((wd / "DONE.json").read_text(encoding="utf-8"))
    if (wd / "FAILED.json").exists():                                     # 한 세션 명령으로 «다시도 실패» 를 우회하지 못하게(재점검 C-1)
        raise RuntimeError(f"⛔ {sub}_{ses} 는 다시도 실패로 빠진 세션이다(FAILED.json) — 돌리지 않는다(중계 ZZZZe-2 ⑤)")
    basis = guard_if_real([inp["ap"], inp["pa"]])
    wd.mkdir(parents=True, exist_ok=True)
    log("onharmony", "새 분석 전처리 입력(dwi AP · PA)", [inp["ap"], inp["pa"]], basis, f"{sub}_{ses}")
    for k, nm in (("ap", "AP.nii.gz"), ("pa", "PA.nii.gz")):
        (wd / nm).unlink(missing_ok=True)
        os.symlink(os.path.abspath(inp[k]), wd / nm)
    shutil.copyfile(inp["bval"], wd / "bvals")
    shutil.copyfile(inp["bvec"], wd / "bvecs")
    jap, jpa = (json.loads(pathlib.Path(inp[k]).read_text()) for k in ("json_ap", "json_pa"))
    bvals = (wd / "bvals").read_text().split()
    (wd / "acqparams.txt").write_text(K.acqparams(jap["PhaseEncodingDirection"], jpa["PhaseEncodingDirection"], jap["TotalReadoutTime"]))
    (wd / "index.txt").write_text(K.index_txt(len(bvals)))
    for i, argv in enumerate(K.prep(K.first_b0(bvals))):
        step(wd, f"01_prep{i}", argv)
    step(wd, "02_topup", K.topup())
    for i, argv in enumerate(K.mask()):
        step(wd, f"03_mask{i}", argv)
    step(wd, "04_eddy", K.eddy())
    step(wd, "05_dtifit", K.dtifit("eddy_out.nii.gz", "eddy_out.eddy_rotated_bvecs", "bvals"))
    snr = session_snr(wd)
    (wd / "snr.json").write_text(json.dumps(snr, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    for i, argv in enumerate(K.register()):
        step(wd, f"06_reg{i}", argv)
    f2(wd, "eddy_out.nii.gz", "bvals", "eddy_out.eddy_rotated_bvecs", "nodif_brain_mask.nii.gz", R.LOW_SHELL["onharmony"], "f2")
    to_std_and_roi(wd, "f2", wd / "fa2std_warp", "nodif_brain_mask", "roi")
    prune_onh(wd)
    info = dict(자료="ON-Harmony", 세션=f"{sub}_{ses}", 근거=basis, eddy=R.EDDY_EXE, 스레드=R.EDDY_NTHR, SNR=snr)
    done(wd, info)
    return info


def session_snr(wd: pathlib.Path, data: str = "eddy_out.nii.gz", bvals: str = "bvals") -> dict:
    """세션 SNR(중계 ZZZZf-6 1 · analysis.snr) — eddy 보정 영상의 b0 부피(b ≤ 50) · 뇌 마스크 · dtifit FA."""
    import nibabel as nib
    import analysis as A  # noqa: E402
    bv = np.loadtxt(wd / bvals).ravel()
    img = nib.load(str(wd / data))
    b0 = np.stack([np.asarray(img.dataobj[..., int(i)], dtype=np.float32) for i in np.flatnonzero(bv <= R.B0_MAX_PREP)], -1)
    mask = np.asarray(nib.load(str(wd / "nodif_brain_mask.nii.gz")).dataobj) > 0
    fa = np.asarray(nib.load(str(wd / "dti_FA.nii.gz")).dataobj, dtype=np.float32)
    return A.snr(b0, mask, fa)


def after_exit(rc: int, tries: int) -> str:
    """하위 세션이 끝난 뒤 할 일 — «끝» · 멈춤 둘(«가드 멈춤» 77 · «디스크 멈춤» 75 — 세지 않고 일괄을 멈춘다) ·
    «잠김»(76 — 세지 않고 이번 일괄에서 다시 띄우지 않는다) · «다시»(처리 실패 — 한 번만) · «실패»(다시도 실패 — 빼고 보고)."""
    if rc == 0:
        return "끝"
    if rc == G.GUARD_EXIT:
        return "가드 멈춤"
    if rc == STOP_EXIT:
        return "디스크 멈춤"
    if rc == BUSY_EXIT:
        return "잠김"
    return "다시" if tries < 2 else "실패"


def prior_failures(work: pathlib.Path) -> dict[tuple[str, str], list[int]]:
    """앞 일괄 · 한 세션 명령의 처리 실패(재점검 C-1 · RUN-3) — `onharmony/attempts.jsonl`(세션 폴더 밖 · 정리 단계가 지우지 않음)의
    종료 코드 가운데 NOT_COUNTED(0 끝 · 77 가드 · 75 디스크 · 76 잠김)가 아닌 것만 센다. 끝 줄 없이 끊긴 시도(재부팅 등)는 세지 않는다."""
    out: dict[tuple[str, str], list[int]] = {}
    p = work / "onharmony" / "attempts.jsonl"
    if p.exists():
        for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
            try:
                r = json.loads(line)
                key, rc = (r["sub"], r["ses"]), int(r["rc"])
            except (ValueError, KeyError, TypeError):                  # 디스크가 차서 반쪽만 쓰인 줄 — 끊긴 시도처럼 세지 않는다(RUN5-2)
                continue
            if rc not in NOT_COUNTED:
                out.setdefault(key, []).append(rc)
    return out


def parent_rc(work: pathlib.Path, key: tuple[str, str], rc: int) -> int:
    """부모 쪽 판정(RUN5-1) — 하위 세션의 보고(출력 쓰기)가 실패해도 처리 결과로 판정한다:
    끝 표지(DONE.json)가 있으면 종료 코드와 무관하게 끝(0) · 처리 실패에 여유 디스크 부족이 겹치면 멈춤(75)."""
    if rc in NOT_COUNTED:
        return rc
    if (work / "onharmony" / f"{key[0]}_{key[1]}" / "DONE.json").exists():
        return 0
    return STOP_EXIT if low_disk(work / "onharmony") else rc


def record_attempt(work: pathlib.Path, key: tuple[str, str], rc: int, by: str, raw: int | None = None) -> None:
    """하위 세션의 끝 줄 — 띄운 쪽(일괄 · 한 세션 명령)만 적는다(하위 세션은 적지 않는다 — 두 번 세지 않게).
    rc = 부모 쪽 판정(parent_rc) · raw = 하위 세션의 실제 종료 코드(다를 때만). 앞 줄이 반쪽이면 줄을 바꾼 뒤 적는다."""
    p = work / "onharmony" / "attempts.jsonl"
    lead = ""
    if p.exists() and p.stat().st_size:
        with open(p, "rb") as fh:
            fh.seek(-1, os.SEEK_END)
            lead = "" if fh.read(1) == b"\n" else "\n"
    rec = dict(utc=_utc(), sub=key[0], ses=key[1], rc=rc, 띄운_쪽=by, **({} if raw is None or raw == rc else {"실제_종료_코드": raw}))
    with open(p, "a") as fh:
        fh.write(lead + json.dumps(rec, ensure_ascii=False) + "\n")


def mark_failed(work: pathlib.Path, key: tuple[str, str], rcs: list[int]) -> None:
    """«다시도 실패» 표지 — 세션 폴더의 FAILED.json. 다시 띄운 일괄도 · 한 세션 명령도 이 세션을 돌리지 않는다."""
    wd = work / "onharmony" / f"{key[0]}_{key[1]}"
    wd.mkdir(parents=True, exist_ok=True)
    (wd / "FAILED.json").write_text(json.dumps(dict(세션=f"{key[0]}_{key[1]}", 처리_실패=len(rcs), 종료_코드=rcs, utc=_utc(),
                                                    규칙="다시도 실패 — 빼고 보고한다 · 다시 돌리지 않는다(중계 ZZZZe-2 ⑤)"),
                                               ensure_ascii=False, indent=1) + "\n", encoding="utf-8")


def settle_failed(work: pathlib.Path, prev: dict) -> list[str]:
    """처리 실패가 둘인데 표지가 없는 세션(표지를 쓰기 전에 끊긴 몫)에 FAILED.json 을 쓴다. 끝 표지(DONE.json)가 함께 있으면 쓰지 않고
    «충돌» 로 돌려준다 — 잠금 뒤에는 생기지 않아야 한다 · 보고하고 주 분석에서 뺀다(재점검 RUN-4)."""
    root, conflict = work / "onharmony", []
    for key, rcs in prev.items():
        if len(rcs) < 2:
            continue
        wd = root / f"{key[0]}_{key[1]}"
        if (wd / "DONE.json").exists():
            conflict.append(f"{key[0]}_{key[1]}")
        elif not (wd / "FAILED.json").exists():
            mark_failed(work, key, rcs)
    return sorted(conflict)


def _launch(sub: str, ses: str) -> subprocess.Popen:
    """하위 세션 — 띄운 쪽과 같은 프로세스 묶음(Ctrl-C 가 함께 닿는다) · WORKER 표지(시도는 띄운 쪽이 적는다)."""
    return subprocess.Popen([sys.executable, str(HERE / "run.py"), "onharmony", sub, ses], env=dict(ENV, **{WORKER: "1"}))


def run_lock(work: pathlib.Path):
    """작업 폴더 잠금 — 일괄과 한 세션 명령이 함께 잡는다(둘이 같은 작업 폴더에서 함께 돌지 않는다 · 재점검 RUN-3 · RUN-4)."""
    return flock_or_busy(work / "onharmony" / "locks" / "run.lock", "다른 일괄 또는 한 세션 명령이 이 작업 폴더에서 도는 중")


def single_onh(sub: str, ses: str, work: pathlib.Path = R.WORK) -> int:
    """한 세션 명령(사람이 직접 · 재점검 RUN-3) — 일괄과 같은 규칙: 작업 폴더 잠금 안에서 · 다시도 실패한 세션은 돌리지 않는다 ·
    하위 세션으로 돌리고 그 종료 코드를 attempts.jsonl 에 적는다 · 처리 실패가 둘이 되면 FAILED.json. 하위 세션의 종료 코드를 돌려준다."""
    key, root = (sub, ses), work / "onharmony"
    with run_lock(work):
        settle_failed(work, prior_failures(work))
        wd = root / f"{sub}_{ses}"
        if (wd / "DONE.json").exists():
            return 0
        if (wd / "FAILED.json").exists():
            raise RuntimeError(f"⛔ {sub}_{ses} 는 다시도 실패로 빠진 세션이다(FAILED.json) — 돌리지 않는다(중계 ZZZZe-2 ⑤)")
        if low_disk(root):
            raise StopRun(f"⛔ 여유 디스크 < {R.MINFREE_GB} GB — 띄우지 않는다(처리 실패가 아니다 · 디스크를 비운 뒤 같은 명령으로 다시)")
        n = len(prior_failures(work).get(key, [])) + 1
        _blog(root / "batch_log.jsonl", "한 세션 명령 · 시작", sub, ses, n)
        raw = _launch(sub, ses).wait()
        rc = parent_rc(work, key, raw)
        record_attempt(work, key, rc, "한 세션 명령", raw)
        _blog(root / "batch_log.jsonl", f"한 세션 명령 · 끝 · 종료 {rc}", sub, ses, n)
        rcs = prior_failures(work).get(key, [])
        if rc not in NOT_COUNTED and len(rcs) >= 2:
            mark_failed(work, key, rcs)
        return rc


def batch_onh(n_parallel: int = R.N_PARALLEL, work: pathlib.Path = R.WORK, cap_s: float = R.CAP_WALL_S) -> dict:
    """첫 세션 110 + s_ref_F2 용 OXF1PRI 반복 10(중계 ZZZZf-3) — 동시 n · 끝난 세션은 건너뜀 · 처리 실패는 같은 명령으로 한 번만 다시 ·
    상한을 넘으면 새 세션을 띄우지 않는다. 작업 폴더 잠금 안에서 돈다 · 처리 실패 수는 일괄을 다시 띄워도 이어 센다(attempts.jsonl) ·
    앞에서 다시도 실패한 세션(FAILED.json)은 띄우지 않는다 · 가드 · 디스크 멈춤은 세지 않고 일괄을 멈춘다 · 잠김은 세지 않는다 ·
    요약은 앞 요약을 읽어 이어 쓴다(재점검 C-1 · RUN-1 – RUN-4 · A 제안 a12)."""
    G.require()
    with run_lock(work):
        return _batch_onh(n_parallel, work, cap_s)


def _batch_onh(n_parallel: int, work: pathlib.Path, cap_s: float) -> dict:
    root = work / "onharmony"
    blog = root / "batch_log.jsonl"
    every = R.onh_sessions() + R.sref_extra()
    prev = prior_failures(work)
    conflict = settle_failed(work, prev)
    todo = [(a, b) for a, b in every if not any((root / f"{a}_{b}" / f).exists() for f in ("DONE.json", "FAILED.json"))]
    start, t0 = _utc(), _now()
    running, failed, busy, capped = {}, [], [], []
    stopped, why, first = {"가드": [], "디스크": []}, {"가드": None, "디스크": None}, None    # 멈춤 — 종류마다 세션 · 까닭 · 먼저 멈춘 종류
    tries = {k: len(prev.get(k, [])) for k in todo}                          # 앞 처리 실패부터 이어 센다
    while todo or running:
        while todo and len(running) < n_parallel and _now() - t0 < cap_s and first is None:
            try:
                G.require()                                                   # 띄우기 전마다 — 가드가 닫혔으면 새 세션을 띄우지 않는다
            except G.GuardError as e:
                first, why["가드"] = "가드", str(e)
                _blog(blog, f"⛔ 가드 — 새 세션을 띄우지 않는다: {e}", "", "", 0)
                break
            if low_disk(root):                                                 # 띄우기 전마다 — 디스크가 모자라면 띄우지 않는다(RUN-1)
                first, why["디스크"] = "디스크", f"여유 디스크 < {R.MINFREE_GB} GB(띄우기 전 확인)"
                _blog(blog, f"⛔ 여유 디스크 < {R.MINFREE_GB} GB — 새 세션을 띄우지 않는다", "", "", 0)
                break
            sub, ses = todo.pop(0)
            tries[(sub, ses)] += 1
            running[(sub, ses)] = _launch(sub, ses)
            _blog(blog, "시작", sub, ses, tries[(sub, ses)])
        if _now() - t0 >= cap_s and todo and first is None:                   # 멈춘 뒤 남은 세션은 멈춤의 «남은_세션» 에 둔다(RUN5-5)
            _blog(blog, f"⛔ 상한 {cap_s / 3600:.0f} 시간 — 띄우지 않는다: " + " · ".join(f"{a}_{b}" for a, b in todo), "", "", 0)
            capped.extend(todo)                                               # 상한 뒤에 돌아온 첫 실패도 더한다(RUN-2)
            todo = []
        for key, pr in list(running.items()):
            raw = pr.poll()
            if raw is None:
                continue
            del running[key]
            rc = parent_rc(work, key, raw)                                     # 끝 표지 · 디스크로 판정(RUN5-1)
            record_attempt(work, key, rc, "일괄", raw)                         # 먼저 적는다 — 일괄이 끊겨도 다음 일괄이 이어 센다
            _blog(blog, f"끝 · 종료 {raw}" + (f" → 판정 {rc}" if rc != raw else ""), *key, tries[key])
            act = after_exit(rc, tries[key])
            if act == "다시":
                todo.append(key)
            elif act == "실패":
                failed.append(key)
                mark_failed(work, key, prior_failures(work)[key])
            elif act == "잠김":
                busy.append(key)
                tries[key] -= 1                                               # 세지 않는다 · 이번 일괄에서 다시 띄우지 않는다
            elif act in ("가드 멈춤", "디스크 멈춤"):
                kind = act.split()[0]
                stopped[kind].append(key)
                tries[key] -= 1                                               # 세지 않는다(attempts.jsonl 에서도)
                why[kind] = why[kind] or f"{key[0]}_{key[1]} 가 {act}(종료 {rc})"
                first = first or kind
        if first is not None and not running:
            break
        _sleep(30)
    names = lambda ks: [f"{a}_{b}" for a, b in ks]  # noqa: E731
    rec = {k: (dict(까닭=why[k], 세션=names(stopped[k]), 남은_세션=names(todo)) if why[k] else None) for k in why}
    sp = root / "batch_summary.json"
    try:
        old = json.loads(sp.read_text(encoding="utf-8")) if sp.exists() else {}
    except ValueError:                                                         # 쓰다 끊긴 요약 — 버리지 않고 옆에 두고 새로 시작(RUN5-2)
        os.replace(sp, sp.with_name(f"batch_summary.broken_{_utc().replace(':', '')}.json"))
        old = {}
    this = dict(시작_utc=start, 걸린_시간_s=round(_now() - t0, 1), 이번에_다시도_실패=names(failed), 먼저_멈춘_까닭=first,
                가드로_멈춤=rec["가드"], 디스크로_멈춤=rec["디스크"], 잠김=names(busy), 상한으로_남김=names(capped))
    summary = dict(실패_다시도_실패=[f"{a}_{b}" for a, b in every if (root / f"{a}_{b}" / "FAILED.json").exists()],   # 누적 — 앞 일괄 것 포함
                   충돌=conflict,
                   규칙="처리 실패 세션은 같은 명령으로 한 번만 다시 돌린다 · 다시도 실패하면 빼고 보고한다(중계 ZZZZe-2 ⑤) · "
                        "처리 실패 수는 일괄을 다시 띄워도 이어 센다 · 가드(77) · 디스크(75) 멈춤과 잠김(76)은 세지 않는다 · "
                        "충돌(끝 표지와 실패 둘)은 보고하고 주 분석에서 뺀다(a12)",
                   주_분석_제외={f"{a}_{b}": w for (a, b), w in R.MAIN_EXCLUDE_ONH.items()},
                   s_ref_F2_반복=[f"{a}_{b}" for a, b in R.sref_extra()],
                   먼저_멈춘_까닭=first, 가드로_멈춤=rec["가드"], 디스크로_멈춤=rec["디스크"], 잠김=this["잠김"], 걸린_시간_s=this["걸린_시간_s"],
                   실행=old.get("실행", []) + [this])                         # 맨 위 멈춤 · 잠김 · 시간은 이번 일괄 것(종료 코드 판단용)
    tmp = sp.with_suffix(".tmp")
    tmp.write_text(json.dumps(summary, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    os.replace(tmp, sp)                                                        # 다 쓴 뒤 바꾼다 — 잘린 요약이 남지 않게(RUN5-2)
    return summary


def _blog(path: pathlib.Path, what: str, sub: str, ses: str, n: int) -> None:
    with open(path, "a") as fh:
        fh.write(json.dumps(dict(utc=_utc(), 무엇=what, sub=sub, ses=ses, 시도=n), ensure_ascii=False) + "\n")


# ── Tong(제공자 전처리본 · topup · eddy 없음 — ZZZZd-0) ────────────────────────
def tong_extract(sub: str, work: pathlib.Path = R.WORK) -> list[pathlib.Path]:
    """분할 zip 을 잠시 합쳐(합친 파일은 바로 지움) 세션마다 dwi 폴더만 꺼낸다."""
    zips = [R.TONG_ZIPS / f"{sub}.z01", R.TONG_ZIPS / f"{sub}.zip"]
    basis = guard_if_real(zips)
    out = work / "tong" / "raw"
    out.mkdir(parents=True, exist_ok=True)
    merged = work / "tong" / f"{sub}_merged.zip"
    if low_disk(work):
        raise StopRun(f"⛔ 여유 디스크 < {R.MINFREE_GB} GB — zip 을 합치지 않는다(처리 실패가 아니다)")
    try:
        subprocess.run(["zip", "-q", "-s", "0", str(zips[1]), "--out", str(merged)], check=True)
        subprocess.run(["unzip", "-q", "-o", str(merged), f"{sub}/ses-*/dwi/*", "-d", str(out)], check=True)
    except subprocess.CalledProcessError as e:
        if low_disk(work):                                                     # 합치거나 꺼내다 디스크가 찼다 — 세지 않는다(RUN5-3)
            raise StopRun(f"⛔ {e.cmd[0]} 종료 {e.returncode} · 여유 디스크 < {R.MINFREE_GB} GB — 멈춘다(처리 실패가 아니다)") from e
        raise
    finally:
        merged.unlink(missing_ok=True)
    files = sorted((out / sub).glob("ses-*/dwi/*"))
    log("tong", "새 분석 — 제공자 전처리본 dwi 꺼냄(합친 zip 은 지움)", zips + files, basis, sub)
    return files


def select_vols(bvals: np.ndarray, bmax: float) -> np.ndarray:
    return np.flatnonzero(bvals <= bmax)


def write_subset(dwi: pathlib.Path, bvals: np.ndarray, bvecs: np.ndarray, idx: np.ndarray, out: pathlib.Path, name: str) -> None:
    import nibabel as nib
    img = nib.load(str(dwi))
    data = np.asarray(img.dataobj)[..., idx]
    nib.save(nib.Nifti1Image(data, img.affine, img.header), str(out / f"{name}.nii.gz"))
    np.savetxt(out / f"{name}.bval", bvals[idx][None, :], fmt="%g")
    np.savetxt(out / f"{name}.bvec", bvecs[:, idx], fmt="%.6f")


def tong_session(sub: str, ses: str, work: pathlib.Path = R.WORK, dwi_dir: pathlib.Path | None = None) -> dict:
    d = dwi_dir or (work / "tong" / "raw" / sub / ses / "dwi")
    wd = work / "tong" / f"{sub}_{ses}"
    if (wd / "DONE.json").exists():
        return json.loads((wd / "DONE.json").read_text(encoding="utf-8"))
    dwi = d / f"{sub}_{ses}_dwi.nii.gz"
    basis = guard_if_real([dwi])
    wd.mkdir(parents=True, exist_ok=True)
    log("tong", "새 분석 처리 입력(dwi)", [dwi, d / f"{sub}_{ses}_dwi.bval", d / f"{sub}_{ses}_dwi.bvec"], basis, f"{sub}_{ses}")
    bvals = np.loadtxt(d / f"{sub}_{ses}_dwi.bval").ravel()
    bvecs = np.loadtxt(d / f"{sub}_{ses}_dwi.bvec")
    write_subset(dwi, bvals, bvecs, select_vols(bvals, R.TONG_MAIN_BMAX), wd, "main")              # 주 = b0 · 1000 · 2000(D-w3)
    write_subset(dwi, bvals, bvecs, select_vols(bvals, R.B0_MAX_PREP), wd, "b0")
    for i, argv in enumerate(K.mask("b0")):
        step(wd, f"03_mask{i}", argv)
    step(wd, "05_dtifit", K.dtifit("main.nii.gz", "main.bvec", "main.bval"))
    for i, argv in enumerate(K.register()):
        step(wd, f"06_reg{i}", argv)
    f2(wd, "main.nii.gz", "main.bval", "main.bvec", "nodif_brain_mask.nii.gz", R.LOW_SHELL["tong"], "f2")
    to_std_and_roi(wd, "f2", wd / "fa2std_warp", "nodif_brain_mask", "roi")
    f2(wd, str(dwi), str(d / f"{sub}_{ses}_dwi.bval"), str(d / f"{sub}_{ses}_dwi.bvec"), "nodif_brain_mask.nii.gz", R.LOW_SHELL["tong"], "f2all")
    to_std_and_roi(wd, "f2all", wd / "fa2std_warp", "nodif_brain_mask", "roi_allshell")               # 민감도 = 세 셸 모두(D-w3)
    for nm in ("main.nii.gz", "b0.nii.gz"):
        (wd / nm).unlink(missing_ok=True)
    info = dict(자료="Tong 2020", 세션=f"{sub}_{ses}", 근거=basis, 처리="제공자 전처리본 · topup · eddy 없음(ZZZZd-0)")
    done(wd, info)
    return info


# ── Wongkornchaovalit(제공자 전처리본 · 기준 계열 정합 하나 — 검토자 ZZZZd-3) ─────
def align_r(ref: np.ndarray, other: np.ndarray, mask: np.ndarray) -> float:
    """TE 계열 정렬 점검 지표 — 마스크 안 평균 b0 정규화 신호 지도의 화소 상관(⑩ 뒤 계산 · 문턱 rconfig.WONG_ALIGN_R_MIN)."""
    return float(np.corrcoef(ref[mask], other[mask])[0, 1])


def mean_dw_ratio(data: np.ndarray, bvals: np.ndarray) -> np.ndarray:
    b0 = data[..., bvals <= R.B0_MAX_PREP].mean(-1)
    dw = data[..., bvals > R.B0_MAX_PREP].mean(-1)
    return np.where(b0 > 0, dw / np.maximum(b0, 1e-6), 0.0)


def own_registration(align: dict[str, float], r_min: float = R.WONG_ALIGN_R_MIN) -> list[str]:
    """정렬 r 이 문턱에 못 미치는 계열(NaN 포함) — 그 계열의 FA 로 따로 정합한다 · 보고 · 민감도 = 그 계열 뺌(중계 ZZZZe-2)."""
    return [s for s, r in align.items() if not r >= r_min]


def _prep_reg(wd: pathlib.Path, base: str) -> None:
    """계열 하나로 뇌 마스크(b0 평균 bet) · dtifit FA · 정합(옛 옵션) — 기준 계열과 따로 정합하는 계열이 같은 단계를 쓴다."""
    bv = np.loadtxt(f"{base}.bval").ravel()
    write_subset(pathlib.Path(f"{base}.nii.gz"), bv, np.loadtxt(f"{base}.bvec"), select_vols(bv, R.B0_MAX_PREP), wd, "b0")
    for i, argv in enumerate(K.mask("b0")):
        step(wd, f"03_mask{i}", argv)
    step(wd, "05_dtifit", K.dtifit(f"{base}.nii.gz", f"{base}.bvec", f"{base}.bval"))
    for i, argv in enumerate(K.register()):
        step(wd, f"06_reg{i}", argv)
    (wd / "b0.nii.gz").unlink(missing_ok=True)


def wong_subject(sub: str, work: pathlib.Path = R.WORK, root: pathlib.Path = R.WONG) -> dict:
    """기준 계열(첫 TE62)의 마스크 · 정합을 열 계열에 쓴다. 정렬 r < 문턱인 계열은 그 계열로 따로 마스크 · 정합(own_<계열>/)."""
    import nibabel as nib
    d = root / sub / "dwi"
    wd = work / "wong" / sub
    if (wd / "DONE.json").exists():
        return json.loads((wd / "DONE.json").read_text(encoding="utf-8"))
    ser = {s: d / f"{sub}_acq-{s}_dwi" for s in R.WONG_SERIES}
    basis = guard_if_real([p.with_suffix(".nii.gz") for p in ser.values()])
    wd.mkdir(parents=True, exist_ok=True)
    log("wong", "새 분석 처리 입력(dwi 10 계열)", [pathlib.Path(f"{p}.nii.gz") for p in ser.values()], basis, sub)
    ref = ser[R.WONG_REF]
    _prep_reg(wd, str(ref))
    msk = np.asarray(nib.load(str(wd / "nodif_brain_mask.nii.gz")).dataobj) > 0
    ref_map = mean_dw_ratio(np.asarray(nib.load(f"{ref}.nii.gz").dataobj, dtype=np.float32), np.loadtxt(f"{ref}.bval").ravel())
    align = {s: align_r(ref_map, mean_dw_ratio(np.asarray(nib.load(f"{p}.nii.gz").dataobj, dtype=np.float32),
                                               np.loadtxt(f"{p}.bval").ravel()), msk) for s, p in ser.items()}
    own = own_registration(align)
    for s, p in ser.items():
        sd = wd / f"own_{s}" if s in own else wd
        if s in own:
            sd.mkdir(exist_ok=True)
            _prep_reg(sd, str(p))
        f2(sd, f"{p}.nii.gz", f"{p}.bval", f"{p}.bvec", "nodif_brain_mask.nii.gz", R.LOW_SHELL["wong"], f"f2_{s}")
        to_std_and_roi(sd, f"f2_{s}", sd / "fa2std_warp", "nodif_brain_mask", f"roi_{s}")
    info = dict(자료="Wongkornchaovalit 2025", 피험자=sub, 근거=basis, 정렬_r=align, 정렬_문턱=R.WONG_ALIGN_R_MIN,
                따로_정합=own, 따로_정합_규칙="문턱 미달 계열은 그 계열의 FA 로 따로 정합 · 보고 · 민감도 = 그 계열 뺌(중계 ZZZZe-2)")
    done(wd, info)
    return info


def main(argv: list[str]) -> None:
    try:
        _main(argv)
    except G.GuardError as e:                                         # 멈춤 셋은 처리 실패와 다른 코드로 끝낸다(batch_onh after_exit)
        print(e, flush=True)
        sys.exit(G.GUARD_EXIT)
    except StopRun as e:
        print(e, flush=True)
        sys.exit(STOP_EXIT)
    except Busy as e:
        print(e, flush=True)
        sys.exit(BUSY_EXIT)
    except Exception:
        if low_disk(R.WORK):                                           # 여유 디스크가 모자랄 때 난 오류는 명령 · 종류와 무관하게 멈춤(RUN-1 · RUN5-3)
            with contextlib.suppress(Exception):                       # 보고 쓰기마저 실패해도 75 로 끝낸다
                traceback.print_exc()
                print(f"⛔ 여유 디스크 < {R.MINFREE_GB} GB — 처리 실패가 아니라 멈춤", flush=True)
            sys.exit(STOP_EXIT)
        raise


def _main(argv: list[str]) -> None:
    cmd = argv[0]
    if cmd == "onharmony":
        if os.environ.get(WORKER) == "1":                              # 일괄 · 한 세션 명령이 띄운 하위 세션
            print(json.dumps(onh_session(argv[1], argv[2]), ensure_ascii=False))
        else:                                                           # 사람이 직접 — 하위 세션을 띄우고 시도를 적는다(RUN-3)
            rc = single_onh(argv[1], argv[2])
            if rc:
                sys.exit(rc)
    elif cmd == "onharmony-batch":
        s = batch_onh(int(argv[1]) if len(argv) > 1 else R.N_PARALLEL)
        print(json.dumps(s, ensure_ascii=False))
        if s["먼저_멈춘_까닭"] == "가드":                               # 멈춘 일괄도 멈춤 코드로 끝난다(77 · 75 · 76)
            raise G.GuardError(f"⛔ 일괄이 가드로 멈춤 — {s['가드로_멈춤']['까닭']}")
        if s["먼저_멈춘_까닭"] == "디스크":
            raise StopRun(f"⛔ 일괄이 디스크 부족으로 멈춤 — {s['디스크로_멈춤']['까닭']}")
        if s["잠김"]:
            raise Busy(f"⛔ 잠긴 세션 {len(s['잠김'])} — 다른 프로세스가 끝난 뒤 같은 명령으로 다시")
    elif cmd == "tong-extract":
        print(len(tong_extract(argv[1])), "파일")
    elif cmd == "tong":
        print(json.dumps(tong_session(argv[1], argv[2]), ensure_ascii=False))
    elif cmd == "wong":
        print(json.dumps(wong_subject(argv[1]), ensure_ascii=False))
    else:
        raise SystemExit(f"모르는 명령 {cmd}")


if __name__ == "__main__":
    main(sys.argv[1:])
