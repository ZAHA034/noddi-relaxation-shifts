"""새 실측(⑨-a) — 처리 파이프라인의 값 · 경로 · 세션 목록(중계 ZZZz-2 · ZZZZd · 고정 문서 §4). 영상은 ⑩(고정 문서 커밋) 뒤에만 연다 — `guard.require()`.

처리는 모두 CPU 다(저자 결정 ZZZZd-1 (B)) · FSL 6.0.7.23(저장소의 부모/fsl) · F2 = AMICO 2.1.1(`envs/reference` · S4 와 같은 설정).
옛 등록 경로와 같은 값은 출처 줄을 단다(`prereg/scripts/lockd/timing_run.sh` · `roi_extract.sh` · `measure_high_odi.py`).
경로는 저장소의 부모를 기준으로 만든다 — 절대경로를 적지 않는다(tests/test_no_absolute_paths.py).
"""
from __future__ import annotations

import os
import pathlib

REPO = pathlib.Path(__file__).resolve().parents[2]
ROOT = REPO.parent                                                       # 저장소의 부모(옛 스크립트와 같은 배치)
FSLDIR = pathlib.Path(os.environ.get("REMEASURE_FSLDIR", "") or ROOT / "fsl")
WORK = pathlib.Path(os.environ.get("REMEASURE_WORK", "") or ROOT / "dmri_work" / "remeasure")
ONH_SEALED = ROOT / "dmri_data" / "onharmony_sealed" / "ds004712"
ONH_META = ROOT / "dmri_data" / "onharmony_meta" / "ds004712"
TONG_ZIPS = ROOT / "dmri_work" / "external" / "tong2020_figshare8851955v6"
WONG = ROOT / "dmri_work" / "external" / "wongkornchaovalit2025" / "MTE-dMRI"
ATLAS = REPO / "prereg" / "artifacts" / "atlas" / "JHU-ICBM-labels-1mm.nii.gz"
ATLAS_XML = REPO / "prereg" / "artifacts" / "atlas" / "JHU-labels.xml"
SESSIONS_ONH = pathlib.Path(__file__).resolve().parent / "sessions_onharmony.tsv"
ACCESS_LOG_ONH = REPO / "prereg" / "artifacts" / "decomp" / "onharmony_image_access_log.json"     # 저자 결정 ZZZZc-0
ACCESS_LOG_EXT = REPO / "prereg" / "artifacts" / "decomp" / "external_data_access_log.json"       # 저자 결정 ZZZz-0

# ── 옛 등록 경로와 같은 값 ───────────────────────────────────────────────────
B0_MAX_PREP = 50.0                    # 첫 b0 = b ≤ 50(timing_run.sh:49)
TOPUP_CONFIG = "etc/flirtsch/b02b0.cnf"   # FSLDIR 기준(timing_run.sh:58)
BET_F = 0.25                          # timing_run.sh:59
FA_REF = "data/standard/FMRIB58_FA_1mm"                                 # timing_run.sh:64 · roi_extract.sh:26
FNIRT_CONFIG = "FA_2_FMRIB58_1mm"     # timing_run.sh:64
FLIRT_DOF = 12                        # timing_run.sh:64
SKELETON = "data/standard/FMRIB58_FA-skeleton_1mm"                      # roi_extract.sh:9
SKELETON_THR = 2000                   # roi_extract.sh:14
PE_VEC = {"j-": "0 -1 0", "j": "0 1 0", "i": "1 0 0", "i-": "-1 0 0"}    # acqparams 줄(사이드카 PhaseEncodingDirection) — 옛 경로는 AP «0 −1 0» · PA «0 1 0» 고정(timing_run.sh:53)

# ── CPU(저자 결정 ZZZZd-1 (B)) ───────────────────────────────────────────────
EDDY_EXE = "eddy_cpu"                 # 옵션은 옛 등록 경로와 같다(timing_run.sh:60 — 기본값 · 반복 5) + 스레드
EDDY_NTHR = 8
N_PARALLEL = 3
NICE = 19
CAP_WALL_S = 5 * 24 * 3600            # 벽시계 5 일 — 넘으면 새 세션을 띄우지 않고 멈춘다
MINFREE_GB = 20                       # 옛 경로와 같다(timing_run.sh:20)

# ── F2(S4 와 같은 설정 · decomp/s4/simulate.py amico_fit) ─────────────────────
DPAR = 1.7e-3
DISO = 3.0e-3
AMICO_B0_THR = 10
FEATURES = ("NDI", "FWF", "ODI", "MD", "FA", "s_low", "s_2000")         # C3 = 앞 셋 · C7 = 일곱(lib.coords)

# ── 자료별 ────────────────────────────────────────────────────────────────────
LOW_SHELL = {"onharmony": 1000.0, "tong": 1000.0, "wong": 700.0}        # 확장 좌표의 낮은 셸(synth.ext_coords)
TONG_SUBJECTS = ("sub-1", "sub-2", "sub-3")
TONG_SESSIONS = tuple(f"ses-c{i:02d}r1" for i in range(1, 11)) + ("ses-c10r2", "ses-c10r3")
TONG_MAIN_BMAX = 2050.0               # 주 = b0 · 1000 · 2000(D-w3) · 민감도 = 세 셸 모두
WONG_SUBJECTS = ("sub-01", "sub-02", "sub-03")
WONG_SERIES = ("TE62", "TE72", "TE82", "TE92", "TE102", "TE112", "TE122", "TE132", "TE62R2", "TE132R2")
WONG_REF = "TE62"                     # 정합 기준 계열 = 첫 TE62 촬영(검토자 ZZZZd-3)
WONG_ALIGN_R_MIN = 0.90               # TE 계열 정렬 점검 문턱(검토자 받음 ZZZZe-2) — 미달 계열은 그 계열의 FA 로 따로 정합 · 보고 · 민감도 = 그 계열 뺌


# 주 분석에서 빼는 ON-Harmony 세션 [A 제안 a3 — 고정 문서 v5 · 중계 ZZZZe-2 확인 뒤] — 처리는 110 세션 모두 같은 방식(사이드카대로) · 주 분석(V4 · 쌍 검정) 뺌 · 민감도 = 넣음
MAIN_EXCLUDE_ONH = {
    ("sub-16794", "ses-NOT4GEP001"): ("README:122 — 역방향 위상 b0 를 찍지 않았고 제공자는 Synb0-DisCo 합성 b0(AP 방향 · 왜곡 없음)로 topup 을 돌렸다 · "
                                      "받은 dir-PA 사이드카는 실측 pepolar 계열처럼 보인다(SeriesNumber 7 · AcquisitionTime 15:06:38 · ImageType ORIGINAL · "
                                      "PhaseEncodingDirection j) — 둘이 어긋나고, 영상을 열기 전에는 dir-PA 의 정체를 확인할 수 없다"),
}


# NOT2ING 의 두 프로토콜(중계 ZZZZf-2 ①) — 첫 세션 사이드카 EchoTime(메타데이터 · 2026-10-26T17:14:09Z 열람 기록): 0.089075 s(TE 89 · 가속 2)는 이 넷 ·
# 나머지 16 명은 0.097967 s(16793 은 0.098356 · TE 98 · 가속 1.5) — 두 프로토콜을 모두 찍은 피험자는 없다(15320 의 반복 002 – 006 도 TE 98)
NOT2ING_TE89 = ("sub-13192", "sub-13305", "sub-14229", "sub-14230")
PAIR_MIN_SUBJECTS = 3                 # 장비 쌍의 공통 피험자가 이보다 적으면 «판정 불가 쌍» — 미리 뺀다(중계 ZZZZf-2 ③)

# s_ref_F2(중계 ZZZZf-3 · [검토자 제안 — ⑩ 허가로 정하 확정]) — 옛 s_ref(5262f43 · compute_s_ref.py)와 같은 입력: OXF1PRI 반복 2 명 × 6 세션.
# 첫 세션 110 에 없는 10 세션(002 – 006)을 같은 파이프라인으로 더 처리한다.
SREF_SUBJECTS = ("sub-13192", "sub-16974")
SREF_SESSIONS = tuple(f"ses-OXF1PRI00{i}" for i in range(1, 7))

# 세션 SNR(중계 ZZZZf-6 1 · [검토자 제안 · S4 결과 뒤 추가] · 방법은 A 제안 · 검토자 확인) — S4 학습 범위 SNR U_log[15, 60](decomp/SPEC_S4.md:76)
SNR_TRAIN = (15.0, 60.0)
SNR_FA_MIN = 0.3                      # 백질 화소(뇌 마스크 ∩ dtifit FA ≥ 0.3)에서 낸다 [A 제안]


def onh_sessions(path: pathlib.Path = SESSIONS_ONH) -> list[tuple[str, str]]:
    """첫 세션 110(첫 세션 규칙 judge_P3.py:26 · NOT3GEM 없음) — (sub, ses)."""
    rows = [ln.split("\t") for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    return [(a, b) for a, b in rows]


def device_of(ses: str) -> str:
    return ses[len("ses-"):-3]


def sref_extra(path: pathlib.Path = SESSIONS_ONH) -> list[tuple[str, str]]:
    """s_ref_F2 를 위해 더 처리할 OXF1PRI 반복 세션 — 첫 세션 110 에 없는 것(중계 ZZZZf-3)."""
    have = set(onh_sessions(path))
    return [(a, b) for a in SREF_SUBJECTS for b in SREF_SESSIONS if (a, b) not in have]


def not2ing_protocol(sub: str) -> str:
    """그 피험자의 NOT2ING 첫 세션 프로토콜(카드 단위 이름) — NOT2ING-89 · NOT2ING-98."""
    return "NOT2ING-89" if sub in NOT2ING_TE89 else "NOT2ING-98"
