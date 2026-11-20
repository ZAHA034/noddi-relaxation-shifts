"""f1_yield runs F1 (cuDIMOT) on GPU 1 in stage1, yielding at once when another process appears and resuming later, logging per chunk."""
from __future__ import annotations

import datetime
import json
import math
import os
import pathlib
import shutil
import signal
import subprocess
import sys
import time

import nibabel as nib
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parents[1]
ROOT = REPO.parent
sys.path.insert(0, str(HERE))
import config as C  # noqa: E402

WORK = C.WORK
F1DIR = WORK / "f1"
LEDGER = F1DIR / "ledger.jsonl"
STATUS = F1DIR / "status.json"
GPU = os.environ.get("F1_GPU", "1")
YIELD_TO_OTHERS = os.environ.get("F1_YIELD", "1") != "0"
IGNORE_PIDS = {149499}
POLL = 30.0                      # Seconds
ORDER = ("E1", "E1-TR", "E3", "E2a", "E2c", "E2b", "E5", "E6", "E6n")   # stage 1 first
SIDE = 100                       # Re-stacking: (SIDE, SIDE, Z, V)


def now():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def log(event, **kw):
    F1DIR.mkdir(parents=True, exist_ok=True)
    rec = {"UTC": now(), "event": event, **kw}
    with open(LEDGER, "a", encoding="utf-8") as f:
        f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    print(json.dumps(rec, ensure_ascii=False), flush=True)


def write_status(**kw):
    STATUS.write_text(json.dumps({"UTC": now(), **kw}, ensure_ascii=False, indent=1), encoding="utf-8")


def gpu_uuid():
    out = subprocess.run(["nvidia-smi", "--query-gpu=index,uuid", "--format=csv,noheader"], capture_output=True, text=True, check=True).stdout
    for line in out.strip().splitlines():
        i, u = [x.strip() for x in line.split(",")]
        if i == GPU:
            return u
    raise RuntimeError(f"GPU {GPU} 없음")


def foreign(uuid, our_pgid=None):
    """Among the compute processes on the GPU, the ones that are ours (our process group) · PIDs not on the ignore list. In the mode that does not yield, always an empty list."""
    if not YIELD_TO_OTHERS:
        return []
    out = subprocess.run(["nvidia-smi", "--query-compute-apps=gpu_uuid,pid", "--format=csv,noheader"],
                         capture_output=True, text=True).stdout
    pids = []
    for line in out.strip().splitlines():
        if not line.strip():
            continue
        u, p = [x.strip() for x in line.split(",")]
        if u != uuid:
            continue
        p = int(p)
        if p in IGNORE_PIDS:
            continue
        if our_pgid is not None:
            try:
                if os.getpgid(p) == our_pgid:
                    continue
            except ProcessLookupError:
                continue
        pids.append(p)
    return pids


def jobs():
    out = []
    for exp in ORDER:
        for prot in C.PROTOCOLS:
            d = WORK / "fits" / exp / prot
            if not d.is_dir():
                continue
            for c in sorted(d.glob("c[0-9][0-9][0-9][0-9]")):
                if (c / "DONE").exists():
                    out.append((exp, prot, int(c.name[1:]), c))
    return out


def env_gpu():
    e = dict(os.environ)
    e["FSLDIR"] = str(ROOT / "fsl")
    e["PATH"] = f"{ROOT / 'fsl' / 'bin'}:{e['PATH']}"
    e["FSLOUTPUTTYPE"] = "NIFTI_GZ"
    e["CUDIMOT"] = str(ROOT / "tools" / "cudimot" / "v0.1.0-fsl6.0.7.23-sm89-watson01")
    e["CUDA_VISIBLE_DEVICES"] = GPU
    return e


def prepare(exp, prot, cid, cdir, tmp):
    """Signal resynthesis · control · restacking. Returns: number of rows n."""
    if tmp.exists():
        shutil.rmtree(tmp)
    raw = tmp / "raw"
    subprocess.run(["uv", "run", "--frozen", "python", str(HERE / "synth.py"), exp, prot, str(cid), str(raw)],
                   cwd=REPO / "envs" / "reference", check=True, capture_output=True)
    new, old = np.load(raw / "rows.npz"), np.load(cdir / "rows.npz")
    if not (np.array_equal(new["cell"], old["cell"]) and np.array_equal(new["ext"], old["ext"])):
        raise RuntimeError("재합성 대조 실패 — rows.npz 의 cell · ext 가 다르다")
    sig = np.asarray(nib.load(raw / "data.nii").dataobj, dtype=np.float32)
    n, V = sig.shape[0], sig.shape[-1]
    sig = sig.reshape(n, V)
    Z = math.ceil(n / (SIDE * SIDE))
    pad = SIDE * SIDE * Z - n
    vol = np.concatenate([sig, np.zeros((pad, V), np.float32)]).reshape(SIDE, SIDE, Z, V)
    mask = np.r_[np.ones(n, np.uint8), np.zeros(pad, np.uint8)].reshape(SIDE, SIDE, Z)
    subj = tmp / "subj"
    subj.mkdir(parents=True)
    nib.save(nib.Nifti1Image(vol, np.eye(4)), subj / "data.nii.gz")
    nib.save(nib.Nifti1Image(mask, np.eye(4)), subj / "nodif_brain_mask.nii.gz")
    for f in ("bvals", "bvecs"):
        shutil.copy(raw / f, subj / f)
    shutil.rmtree(raw)
    return n


def harvest(tmp, cdir, n):
    d = tmp / "subj.NODDI_Watson"
    g = lambda nm: nib.load(d / f"{nm}.nii.gz").get_fdata().reshape(-1)[:n].astype(np.float32)
    np.savez(cdir / "f1.npz", ndi=g("mean_fintra"), fwf=g("mean_fiso"), odi=g("OD"))
    (cdir / "F1_DONE").write_text(now() + "\n")


def kill_group(pgid):
    for sig, wait in ((signal.SIGTERM, 10), (signal.SIGKILL, 5)):
        try:
            os.killpg(pgid, sig)
        except ProcessLookupError:
            return
        t = time.time()
        while time.time() - t < wait:
            try:
                os.killpg(pgid, 0)
            except ProcessLookupError:
                return
            time.sleep(0.5)


def wait_free(uuid, why):
    """Until other people's processes · YIELD are gone (confirmed twice in a row). False if STOP."""
    clean = 0
    logged = False
    while clean < 2:
        if (F1DIR / "STOP").exists():
            return False
        f = foreign(uuid)
        y = (F1DIR / "YIELD").exists()
        if f or y:
            clean = 0
            if not logged:
                log("대기", why=why, foreign=f, yield_file=y)
                logged = True
            write_status(state="대기", foreign=f, yield_file=y)
        else:
            clean += 1
        time.sleep(POLL)
    if logged:
        log("비었음")
    return True


def run():
    F1DIR.mkdir(parents=True, exist_ok=True)
    uuid = gpu_uuid()
    todo = [j for j in jobs() if not (j[3] / "F1_DONE").exists()]
    total = len(jobs())
    log("시작", pid=os.getpid(), gpu=GPU, yield_to_others=YIELD_TO_OTHERS, remaining=len(todo), total=total, ignore=sorted(IGNORE_PIDS))
    tmp = WORK / "f1" / "tmp"
    i = 0
    while i < len(todo):
        exp, prot, cid, cdir = todo[i]
        name = f"{exp}/{prot}/c{cid:04d}"
        if not wait_free(uuid, "시작 전"):
            log("멈춤(STOP)", chunk=name)
            return
        try:
            n = prepare(exp, prot, cid, cdir, tmp)
        except Exception as e:                                   # noqa: BLE001
            log("오류", chunk=name, stage="준비", error=str(e)[:300])
            i += 1
            continue
        log("덩어리 시작", chunk=name, rows=n, done=total - len(todo) + i)
        t0 = time.time()
        with open(tmp / "cascade.log", "w") as lf:
            p = subprocess.Popen(["bash", str(REPO / "prereg" / "scripts" / "cudimot_cascade.sh"), str(tmp / "subj")],
                                 env=env_gpu(), stdout=lf, stderr=subprocess.STDOUT, start_new_session=True)
        pgid = p.pid
        write_status(state="실행", chunk=name, pgid=pgid, rows=n, done=total - len(todo) + i, total=total)
        reason = None
        while p.poll() is None:
            time.sleep(POLL)
            f = foreign(uuid, pgid)
            if f:
                reason = f"남의 프로세스 {f}"
            elif (F1DIR / "YIELD").exists():
                reason = "YIELD 파일"
            elif (F1DIR / "STOP").exists():
                reason = "STOP 파일"
            if reason:
                kill_group(pgid)
                p.wait()
                log("중단", chunk=name, reason=reason, elapsed_s=round(time.time() - t0))
                break
        if reason:
            if reason == "STOP 파일":
                log("멈춤(STOP)", chunk=name)
                write_status(state="멈춤", chunk=name)
                return
            continue                                             # Redo the same chunk from the start
        if p.returncode != 0:
            log("오류", chunk=name, stage="cuDIMOT", returncode=p.returncode, log_tail=(tmp / "cascade.log").read_text()[-400:])
            i += 1
            continue
        harvest(tmp, cdir, n)
        shutil.rmtree(tmp)
        log("덩어리 끝", chunk=name, elapsed_s=round(time.time() - t0))
        i += 1
    log("전부 끝", total=total)
    write_status(state="끝", total=total)


def status():
    js = jobs()
    done = sum((c / "F1_DONE").exists() for *_, c in js)
    print(f"F1 덩어리 {done}/{len(js)} 끝")
    if STATUS.exists():
        print(STATUS.read_text(encoding="utf-8"))
    if LEDGER.exists():
        print("".join(LEDGER.read_text(encoding="utf-8").splitlines(True)[-5:]))


if __name__ == "__main__":
    {"run": run, "status": status}[sys.argv[1]]()
