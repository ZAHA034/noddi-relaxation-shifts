"""S4 — 사양 `decomp/SPEC_S4.md` §3 의 값과 표본 추출(사전 · 조건 · 성가신 변수). numpy 만 — 메인 환경과 `envs/reference` 양쪽에서 불러온다.

이 파일의 값은 SPEC_S4 를 옮긴 것이다. 다르면 SPEC_S4 가 정본이다.
"""
from __future__ import annotations

import math
import os
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "stage1"))
import config as C  # noqa: E402

WORK = C.WORK / "s4"                                   # 저장소 밖 작업 폴더(시뮬레이션 자료 · 표 · 가중치)
RESULTS = C.RESULTS / "s4"                             # 저장소 안 결과(해시 · 표 · 카드 · 지도)

# ── θ* = (ΔNDI, ΔFISO, ΔODI, NDI_A, FISO_A, ODI_A) — SPEC_S4 §3-1 · §3-2 ─────────────────
THETA_NAMES = ("dNDI", "dFISO", "dODI", "NDI_A", "FISO_A", "ODI_A")
BOX = {"NDI": (0.30, 0.80), "FISO": (0.01, 0.30), "ODI": (0.08, 0.40)}        # 조직_A · 조직_B 범위
D_MAIN = np.array([0.10, 0.05, 0.10])                   # 주 보고 사전 반폭(ΔNDI · ΔFISO · ΔODI)
D_TRAIN = 2.0 * D_MAIN                                  # 학습 사전 반폭(두 배)
PIX_SD = np.array([0.03, 0.02, 0.03])                   # ROI 안 이질성 SD(NDI · FISO · ODI) — §3-6
PIX_CLIP = {"NDI": (0.25, 0.85), "FISO": (0.0, 0.35), "ODI": (0.04, 0.50)}   # 화소 자름 — 표의 범위

# ── 이완 — SPEC_S4 §3-3 ───────────────────────────────────────────────────────
T2IN_RANGE, R_RANGE, T2ISO = (70.0, 95.0), (1.0, 1.8), 1000.0
T1C_RANGE, T1ISO = (850.0, 1100.0), C.T1_ISO
SCEN = {"i": 1.3, "ii": 1 / 1.3, "iii": 1.0}              # T1_in/T1_en
L_LEVELS = ("L0", "L1", "L2", "L3")
L_WEIGHTS = (0.25, 0.35, 0.20, 0.20)
L1_ST2 = (0.0, 0.06)
L1_RHO = (0.0, 0.8)
L1_ST1 = 0.04
L2_T2, L2_T1 = 0.0984, 0.0441
L3_EN, L3_IN, L3_T1 = (2.1, 3.4), (0.77, 1.02), (1.17, 1.56)

# ── 생성기 가정 — §3-4 ────────────────────────────────────────────────────────
GENERATORS = C.GENERATORS                               # T-in · T-app · T-none
SCENARIOS = ("i", "ii", "iii")
DPAR_TRUE = C.DPAR_TRUE                                 # 1.7 · 2.0 · 2.26 · 2.5 · 3.0
DPAR_FIT = 1.7

# ── 획득 — §3-5 ───────────────────────────────────────────────────────────────
TE_RANGE, TR_RANGE = (60.0, 110.0), (3000.0, 11000.0)   # ms
ACQ_MIX = {"match": 0.30, "te": 0.15, "tr": 0.15, "free": 0.40}
DESIGNS = {"S": "P2a", "V": "P2b"}                      # 설계 → 기록된 bval · bvec(protocols/)
COILS = (32, 64)
COIL_L = {32: 4, 64: 8}                                 # 비중심 카이 L [설계 — §3-5]
SNR_RANGE = (15.0, 60.0)                                # 로그 균일
NEFF_RANGE = (16.0, 200.0)                              # 로그 균일 · N_eff ≤ K
K_PIX = 256

# ── 층(§7) ─────────────────────────────────────────────────────────────────────
ACQ_TYPES = ("match", "te", "tr", "tetr", "tetr_design")

SEEDS = dict(train=20261101, test=20261102, oracle=20261103, cards=20261104, maps=20261105, inject=20261106)


def feasible(theta: np.ndarray) -> np.ndarray:
    """θ* 행마다 조직_A · 조직_B 가 BOX 안인가."""
    A = theta[:, 3:6]
    B = A + theta[:, 0:3]
    lo = np.array([BOX["NDI"][0], BOX["FISO"][0], BOX["ODI"][0]])
    hi = np.array([BOX["NDI"][1], BOX["FISO"][1], BOX["ODI"][1]])
    return np.all((A >= lo) & (A <= hi) & (B >= lo) & (B <= hi), axis=1)


def sample_theta(rng: np.random.Generator, n: int, half: np.ndarray = D_TRAIN) -> np.ndarray:
    """절단 균일 사전 — A 는 BOX 균일 · Δ 는 ±half 균일 · B 가 BOX 밖이면 다시 뽑는다(§3-2)."""
    out = np.empty((0, 6))
    lo = np.array([BOX["NDI"][0], BOX["FISO"][0], BOX["ODI"][0]])
    hi = np.array([BOX["NDI"][1], BOX["FISO"][1], BOX["ODI"][1]])
    while len(out) < n:
        m = max(2 * (n - len(out)), 1024)
        A = rng.uniform(lo, hi, (m, 3))
        D = rng.uniform(-half, half, (m, 3))
        t = np.c_[D, A]
        out = np.vstack([out, t[feasible(t)]])
    return out[:n]


def in_main(theta: np.ndarray, scale: float = 1.0) -> np.ndarray:
    """주 사전(반폭 D_MAIN × scale)의 지지 안인가 — 중요도 가중은 이 지시함수다(§3-2 · [대수])."""
    return feasible(theta) & np.all(np.abs(theta[:, 0:3]) <= D_MAIN * scale + 1e-12, axis=1)


def prior_var(rng: np.random.Generator, U: np.ndarray, scale: float = 1.0, n: int = 1_000_000) -> np.ndarray:
    """Var_prior(uᵀΔ) — 주 사전(절단 포함 · 반폭 × scale)에서 몬테카를로. U: (k, 3) 방향 → (k,)."""
    t = sample_theta(rng, n, D_MAIN * scale)
    return np.var(t[:, 0:3] @ np.asarray(U).T, axis=0)


def sample_relax_A(rng: np.random.Generator, n: int, scen: np.ndarray):
    """사이트 A 의 사람: (T2_in, T2_en, T2_iso) · (T1_in, T1_en, T1_iso) — §3-3."""
    t2in = rng.uniform(*T2IN_RANGE, n)
    r = rng.uniform(*R_RANGE, n)
    c = rng.uniform(*T1C_RANGE, n)
    ratio = np.array([SCEN[s] for s in scen])
    t2 = np.c_[t2in, t2in / r, np.full(n, T2ISO)]
    t1 = np.c_[c * np.sqrt(ratio), c / np.sqrt(ratio), np.full(n, T1ISO)]
    return t2, t1


def _trunc_bivariate(rng, n, rho):
    out = np.empty((n, 2))
    for i in range(n):
        while True:
            z = rng.standard_normal(2)
            L = np.array([[1.0, 0.0], [rho[i], math.sqrt(1 - rho[i] ** 2)]])
            if z @ z <= 4.0:                                  # 마할라노비스 ≤ 2(표준 정규 공간에서)
                out[i] = L @ z
                break
    return out


def sample_level(rng: np.random.Generator, n: int, levels=None):
    """B 의 이완 수준과 그 매개변수(조건 a 에 들어간다) — §3-3."""
    lv = levels if levels is not None else rng.choice(L_LEVELS, n, p=L_WEIGHTS)
    p = {k: np.zeros(n) for k in ("s_t2", "rho", "s_t1", "d_t2", "d_t1", "alloc_en", "m_en", "m_in", "m_t1")}
    p["m_en"][:] = p["m_in"][:] = p["m_t1"][:] = 1.0
    for i, l in enumerate(lv):
        if l == "L1":
            p["s_t2"][i] = rng.uniform(*L1_ST2)
            p["rho"][i] = rng.choice(L1_RHO)
            p["s_t1"][i] = L1_ST1
        elif l == "L2":
            p["d_t2"][i], p["d_t1"][i] = L2_T2, L2_T1
            p["alloc_en"][i] = float(rng.random() < 0.5)          # (b) 세포 밖만 / (a) 두 구획
        elif l == "L3":
            p["m_en"][i] = rng.uniform(*L3_EN)
            p["m_in"][i] = rng.uniform(*L3_IN)
            p["m_t1"][i] = rng.uniform(*L3_T1)
    return np.asarray(lv), p


def relax_B(rng: np.random.Generator, lv, p, t2A, t1A):
    """B 의 이완 = A 의 이완 ⊕ 수준(무작위 z 는 성가신 변수)."""
    n = len(lv)
    t2B, t1B = t2A.copy(), t1A.copy()
    z = _trunc_bivariate(rng, n, p["rho"])
    zt1 = np.clip(rng.standard_normal(n), -2, 2)
    for i, l in enumerate(lv):
        if l == "L1":
            t2B[i, 0] *= 1 + p["s_t2"][i] * z[i, 0]
            t2B[i, 1] *= 1 + p["s_t2"][i] * z[i, 1]
            t1B[i, 0:2] *= 1 + p["s_t1"][i] * zt1[i]
        elif l == "L2":
            if p["alloc_en"][i]:
                t2B[i, 1] *= 1 + p["d_t2"][i]
                t1B[i, 1] *= 1 + p["d_t1"][i]
            else:
                t2B[i, 0:2] *= 1 + p["d_t2"][i]
                t1B[i, 0:2] *= 1 + p["d_t1"][i]
        elif l == "L3":
            t2B[i, 1] *= p["m_en"][i]
            t2B[i, 0] *= p["m_in"][i]
            t1B[i, 0:2] *= p["m_t1"][i]
    return t2B, t1B


def sample_acq(rng: np.random.Generator, n: int, kinds=None):
    """두 사이트의 획득 — §3-5. kinds 를 주면 그 층(ACQ_TYPES)으로 뽑는다(시험 자료)."""
    if kinds is None:
        kinds = rng.choice(list(ACQ_MIX), n, p=list(ACQ_MIX.values()))
    TEa, TRa = rng.uniform(*TE_RANGE, n), rng.uniform(*TR_RANGE, n)
    TEb, TRb = TEa.copy(), TRa.copy()
    dA = rng.choice(list(DESIGNS), n)
    dB = dA.copy()
    for i, k in enumerate(kinds):
        if k in ("te", "tetr", "tetr_design", "free"):
            TEb[i] = rng.uniform(*TE_RANGE)
        if k in ("tr", "tetr", "tetr_design", "free"):
            TRb[i] = rng.uniform(*TR_RANGE)
        if k == "free":
            dB[i] = rng.choice(list(DESIGNS))
        if k == "tetr_design":
            dB[i] = "V" if dA[i] == "S" else "S"
    coilA, coilB = rng.choice(COILS, n), rng.choice(COILS, n)
    snrA = np.exp(rng.uniform(*np.log(SNR_RANGE), n))
    snrB = np.exp(rng.uniform(*np.log(SNR_RANGE), n))
    neff = np.exp(rng.uniform(*np.log(NEFF_RANGE), n))
    return dict(kind=np.asarray(kinds), TE_A=TEa, TE_B=TEb, TR_A=TRa, TR_B=TRb, design_A=dA, design_B=dB,
                coil_A=coilA, coil_B=coilB, snr_A=snrA, snr_B=snrB, neff=np.minimum(neff, K_PIX))


def design_b0(design: str):
    """설계 파일이 정하는 b0 의 b 값 · 수(조건에 그대로)."""
    bvals, _ = C.protocol(DESIGNS[design])
    b0 = bvals <= 50
    return float(np.round(bvals[b0].max())), int(b0.sum())


def cond_vector(acq: dict, lv, p, gen, scen) -> np.ndarray:
    """조건 a(§3-1) — 숫자 열. 순서는 COND_NAMES."""
    n = len(lv)
    b0 = {d: design_b0(d) for d in DESIGNS}
    cols = [acq["TE_A"], acq["TE_B"], acq["TR_A"] / 1000, acq["TR_B"] / 1000,
            (acq["design_A"] == "V").astype(float), (acq["design_B"] == "V").astype(float),
            np.array([b0[d][0] for d in acq["design_A"]]), np.array([b0[d][0] for d in acq["design_B"]]),
            np.array([b0[d][1] for d in acq["design_A"]], float), np.array([b0[d][1] for d in acq["design_B"]], float),
            (acq["coil_A"] == 64).astype(float), (acq["coil_B"] == 64).astype(float),
            np.log(acq["snr_A"]), np.log(acq["snr_B"]), np.log(acq["neff"])]
    cols += [(np.asarray(lv) == l).astype(float) for l in L_LEVELS]
    cols += [p[k] for k in ("s_t2", "rho", "s_t1", "d_t2", "d_t1", "alloc_en", "m_en", "m_in", "m_t1")]
    cols += [(np.asarray(gen) == g).astype(float) for g in GENERATORS]
    cols += [(np.asarray(scen) == s).astype(float) for s in SCENARIOS]
    out = np.c_[tuple(cols)].astype(np.float32)
    assert out.shape == (n, len(COND_NAMES))
    return out


COND_NAMES = (("TE_A", "TE_B", "TR_A_s", "TR_B_s", "design_A_V", "design_B_V", "b0val_A", "b0val_B", "nb0_A", "nb0_B",
               "coil_A_64", "coil_B_64", "log_snr_A", "log_snr_B", "log_neff")
              + tuple(f"L_{l}" for l in L_LEVELS)
              + ("s_t2", "rho", "s_t1", "d_t2", "d_t1", "alloc_en", "m_en", "m_in", "m_t1")
              + tuple(f"gen_{g}" for g in GENERATORS) + tuple(f"scen_{s}" for s in SCENARIOS))


def roi_noise(rng: np.random.Generator, m: np.ndarray, S: np.ndarray, neff: np.ndarray, K: int = K_PIX) -> np.ndarray:
    """x = m + ε · ε ~ N(0, Σ̂ (1/N_eff − 1/K)) — §3-6 5. m: (n, d) · S: (n, d, d) · neff: (n,)."""
    f = np.clip(1.0 / neff - 1.0 / K, 0.0, None)
    L = np.linalg.cholesky(S * f[:, None, None] + 1e-12 * np.eye(S.shape[1])[None])
    return m + np.einsum("nij,nj->ni", L, rng.standard_normal(m.shape))
