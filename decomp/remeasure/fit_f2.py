"""새 실측(⑨-a) — F2(AMICO 2.1.1 NODDI · S4 와 같은 설정 · decomp/s4/simulate.py amico_fit) 적합과 확장 좌표 지도
(MD · FA · s_low · s_2000 — decomp/stage1/synth.py ext_coords 와 같은 식). `envs/reference` 에서 돈다.

    python fit_f2.py <작업 폴더> <dwi> <bvals> <bvecs> <마스크> <낮은 셸 b> <출력 접두> [스레드]
출력(원 공간): <작업 폴더>/<접두>_{NDI,FWF,ODI,MD,FA,s_low,s_2000,extvalid,f2valid}.nii.gz
  (extvalid = 확장 좌표가 정의된 화소 1 · f2valid = F2 적합이 뜻 있는 화소 1 — 둘이 «적합 유효 지도»(중계 ZZZZe-2 ②))
실자료 경로면 guard.require() 를 먼저 부른다(run.py 가 이미 불렀어도 다시 — 이 스크립트만 따로 돌 때를 막는다).
"""
from __future__ import annotations

import os
import pathlib
import shutil
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "stage1"))
import rconfig as R  # noqa: E402

CHUNK = 20000                                                      # ext_coords 는 화소 × 행 × 7 × 7 배열을 만든다 — 메모리 때문에 나눈다


def ext_maps(sig: np.ndarray, bvals: np.ndarray, bvecs: np.ndarray, low: float, ext_fn) -> tuple[np.ndarray, np.ndarray]:
    """화소 신호(n, V) → (확장 좌표 (n, 4), 유효 (n,)) — ext_fn(= synth.ext_coords) 을 CHUNK 화소씩.
    0 이하 값이 하나라도 있는 화소(마스크 가장자리 등 — 시뮬레이터 신호에는 없다)는 로그 적합이 풀리지 않아 빼고(유효 0 · 값 0),
    그래도 남는 특이 화소는 덩어리를 반씩 나눠 찾아 뺀다."""
    out = np.zeros((len(sig), 4), np.float32)
    ok = np.all(sig > 0, axis=1)
    idx = np.flatnonzero(ok)
    for a in range(0, len(idx), CHUNK):
        sel = idx[a:a + CHUNK]
        val = _safe(ext_fn, sig[sel].astype(float), bvals, bvecs, low)
        fin = np.isfinite(val).all(1)
        out[sel[fin]] = val[fin]
        ok[sel[~fin]] = False
    return out, ok


def _safe(fn, s: np.ndarray, bvals, bvecs, low) -> np.ndarray:
    try:
        return np.asarray(fn(s, bvals, bvecs, low), float)
    except np.linalg.LinAlgError:
        if len(s) == 1:
            return np.full((1, 4), np.nan)
        h = len(s) // 2
        return np.vstack([_safe(fn, s[:h], bvals, bvecs, low), _safe(fn, s[h:], bvals, bvecs, low)])


def f2_valid(msk: np.ndarray, b0mean: np.ndarray, maps: list[np.ndarray]) -> np.ndarray:
    """F2 적합 유효 지도(중계 ZZZZe-2 ②) — 마스크 안 · 평균 b0 > 0 · 세 지도 모두 유한.
    AMICO 는 평균 b0 ≤ 0 인 화소의 신호를 0 으로 만든 뒤 적합한다(load_data 의 b0 정규화 · b0_min_signal 0) — 그 값은 뜻이 없다."""
    ok = msk & (b0mean > 0)
    for m in maps:
        ok &= np.isfinite(m)
    return ok


def scheme(bvals: np.ndarray, bvecs: np.ndarray) -> np.ndarray:
    """AMICO 방식 표(gx gy gz b) — 방향은 단위 벡터로(0 벡터는 그대로) · synth.main 과 같다."""
    g = bvecs.copy()
    nrm = np.linalg.norm(g, axis=1)
    g[nrm > 0] /= nrm[nrm > 0, None]
    return np.c_[g, bvals]


def main():
    import nibabel as nib
    wd, dwi, bval_f, bvec_f, mask_f, low, prefix = sys.argv[1:8]
    nth = int(sys.argv[8]) if len(sys.argv) > 8 else 8
    wd, low = pathlib.Path(wd), float(low)
    import run  # noqa: E402 — 실자료면 가드
    run.guard_if_real([pathlib.Path(dwi)])
    import amico
    import synth  # noqa: E402 — decomp/stage1(amico 필요)
    bvals = np.loadtxt(bval_f).ravel()
    bvecs = np.loadtxt(bvec_f)
    bvecs = bvecs.T if bvecs.shape[0] == 3 else bvecs
    img = nib.load(dwi)
    msk = np.asarray(nib.load(mask_f).dataobj) > 0
    study = wd / f"amico_{prefix}"
    subj = study / "subj"
    subj.mkdir(parents=True, exist_ok=True)
    for src, dst in ((dwi, "DWI.nii.gz"), (mask_f, "mask.nii.gz")):
        if not (subj / dst).exists():
            os.symlink(os.path.abspath(src), subj / dst)
    np.savetxt(subj / "DWI.scheme", scheme(bvals, bvecs), fmt="%.6f", header="VERSION: BVECTOR", comments="")
    amico.setup()
    ae = amico.Evaluation(str(study), "subj")
    ae.set_config("nthreads", nth)
    ae.set_config("BLAS_nthreads", 1)
    ae.load_data(dwi_filename="DWI.nii.gz", scheme_filename="DWI.scheme", mask_filename="mask.nii.gz", b0_thr=R.AMICO_B0_THR)
    ae.set_model("NODDI")
    ae.model.set(dPar=R.DPAR, dIso=R.DISO)
    ae.generate_kernels(regenerate=True)                                       # 세션마다 방식표(회전된 bvec)가 달라 핵을 새로 만든다
    ae.load_kernels()
    ae.fit()
    ae.save_results()
    o = subj / "AMICO" / "NODDI"
    for nm in ("NDI", "FWF", "ODI"):
        shutil.copyfile(o / f"fit_{nm}.nii.gz", wd / f"{prefix}_{nm}.nii.gz")
    data = np.asarray(img.dataobj, dtype=np.float32)
    f2m = [np.asarray(nib.load(str(wd / f"{prefix}_{nm}.nii.gz")).dataobj, dtype=np.float32) for nm in ("NDI", "FWF", "ODI")]
    okf2 = f2_valid(msk, data[..., bvals <= R.AMICO_B0_THR].mean(-1), f2m)               # AMICO 가 b0 로 보는 부피(b0_thr)
    nib.save(nib.Nifti1Image(okf2.astype(np.float32), img.affine), str(wd / f"{prefix}_f2valid.nii.gz"))
    ext, valid = ext_maps(data[msk], bvals, bvecs, low, synth.ext_coords)
    for j, nm in enumerate(("MD", "FA", "s_low", "s_2000", "extvalid")):
        vol = np.zeros(msk.shape, np.float32)
        vol[msk] = ext[:, j] if j < 4 else valid
        nib.save(nib.Nifti1Image(vol, img.affine), str(wd / f"{prefix}_{nm}.nii.gz"))      # float32 그대로(원 머리의 정수형을 물려받지 않게)
    shutil.rmtree(study)
    print(f"F2 · 확장 좌표 끝: {wd.name} {prefix}", flush=True)


if __name__ == "__main__":
    import guard as G  # noqa: E402
    try:
        main()
    except G.GuardError as e:                                       # 가드 멈춤은 처리 실패가 아니다 — 정한 코드로 끝낸다(run.step 이 다시 올린다)
        print(e, flush=True)
        sys.exit(G.GUARD_EXIT)
