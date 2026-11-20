"""새 실측(⑨-a) — FSL 명령 조립(순수 함수 · 실행하지 않는다). 옛 등록 경로(prereg/scripts/lockd/timing_run.sh:57–64)와 같은 옵션이고,
다른 것은 eddy 실행 파일(eddy_cuda → eddy_cpu · --nthr — 저자 결정 ZZZZd-1 (B))과 acqparams 의 위상 방향을 사이드카에서 읽는 것뿐이다.
"""
from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import rconfig as R  # noqa: E402


def fsl(exe: str) -> str:
    return str(R.FSLDIR / "bin" / exe)


def first_b0(bvals) -> int:
    """AP 의 첫 b0 색인 — b ≤ 50(timing_run.sh:49)."""
    return [i for i, b in enumerate(bvals) if float(b) <= R.B0_MAX_PREP][0]


def acqparams(pe_ap: str, pe_pa: str, trt: float) -> str:
    """두 줄(AP · PA) — 위상 방향은 사이드카 PhaseEncodingDirection · 읽기 시간은 AP 의 TotalReadoutTime(timing_run.sh:47 · :53)."""
    return f"{R.PE_VEC[pe_ap]} {trt}\n{R.PE_VEC[pe_pa]} {trt}\n"


def index_txt(n_vol: int) -> str:
    return " ".join(["1"] * n_vol) + "\n"


def prep(idx0: int) -> list[list[str]]:
    """b0 한 쌍(AP 첫 b0 + PA 첫 볼륨) — timing_run.sh:50 – 52."""
    return [[fsl("fslroi"), "AP.nii.gz", "b0_AP", str(idx0), "1"],
            [fsl("fslroi"), "PA.nii.gz", "b0_PA", "0", "1"],
            [fsl("fslmerge"), "-t", "b0s", "b0_AP", "b0_PA"]]


def topup() -> list[str]:
    return [fsl("topup"), "--imain=b0s", "--datain=acqparams.txt", f"--config={R.FSLDIR / R.TOPUP_CONFIG}",
            "--out=topup_res", "--iout=hifi_b0", "--fout=field"]


def mask(src: str = "hifi_b0") -> list[list[str]]:
    """평균 b0 → bet(f 0.25 · -m) — timing_run.sh:59."""
    return [[fsl("fslmaths"), src, "-Tmean", f"{src}_mean"], [fsl("bet"), f"{src}_mean", "nodif_brain", "-m", "-f", str(R.BET_F)]]


def eddy(nthr: int = R.EDDY_NTHR) -> list[str]:
    """timing_run.sh:60 과 같은 옵션(기본값 · 반복 5) — 실행 파일만 eddy_cpu 와 --nthr."""
    return [fsl(R.EDDY_EXE), "--imain=AP.nii.gz", "--mask=nodif_brain_mask", "--acqp=acqparams.txt", "--index=index.txt",
            "--bvecs=bvecs", "--bvals=bvals", "--topup=topup_res", "--out=eddy_out", "--verbose", f"--nthr={nthr}"]


def dtifit(data: str, bvecs: str, bvals: str, mask_img: str = "nodif_brain_mask") -> list[str]:
    """timing_run.sh:61."""
    return [fsl("dtifit"), "-k", data, "-o", "dti", "-m", mask_img, "-r", bvecs, "-b", bvals]


def register() -> list[list[str]]:
    """FA → FMRIB58_FA_1mm(flirt 12 자유도 + fnirt) — timing_run.sh:64 그대로(옛 warp 는 쓰지 않는다 · 검토자 ZZZZd-3)."""
    ref = str(R.FSLDIR / R.FA_REF)
    return [[fsl("flirt"), "-in", "dti_FA", "-ref", ref, "-omat", "fa2std.mat", "-out", "fa2std_lin", "-dof", str(R.FLIRT_DOF)],
            [fsl("fnirt"), "--in=dti_FA", "--aff=fa2std.mat", f"--config={R.FNIRT_CONFIG}", "--cout=fa2std_warp", "--iout=fa2std_nl"]]


def applywarp(src: str, out: str, warp: str = "fa2std_warp", interp: str = "trilinear") -> list[str]:
    """지도 → 표준 공간(roi_extract.sh:26 — trilinear · 마스크는 nn)."""
    return [fsl("applywarp"), "-i", src, "-r", str(R.FSLDIR / R.FA_REF), "-w", warp, "-o", out, f"--interp={interp}"]
