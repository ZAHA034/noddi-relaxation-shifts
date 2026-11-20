"""decomp 1 단계 — 격자 · 프로토콜 · 생성기 축 · 실험 행 목록 (사양 `decomp/SPEC_v1.md` §5–§7).

numpy 만 쓴다 — 메인 환경(분석 · 시험)과 `envs/reference`(합성 · AMICO) 양쪽에서 불러온다.
이 파일의 값은 SPEC_v1 의 값을 옮긴 것이다. 다르면 SPEC_v1 이 정본이다.
"""
from __future__ import annotations

import itertools
import math
import os
import pathlib

import numpy as np

WORK = pathlib.Path(os.environ.get("DECOMP_WORK", "") or pathlib.Path(__file__).resolve().parents[3] / "dmri_work" / "decomp_stage1")   # 저장소 밖(부모 폴더) · 점검용 덮어쓰기
PROTO_DIR = pathlib.Path(__file__).resolve().parents[3] / "dmri_work" / "decomp_stage1" / "protocols"

# ── 결과 갈래 · 출력 경로 (중계 ZZZ-3 나 · 2026-10-24) ──────────────────────
# 갈래 = 적합기. 정본(F2)은 "" · F1 판정(SPEC_v3 §0 예외)은 "_F1" — 결과 디렉터리 results/stage{n}{TAG} 와 작업 폴더의 방향 저장소 이름에 붙는다.
# 출력은 이 파일 위치(저장소)에서 정한다 — 저장소 밖 사본에서 돌려 결과가 사본에 남은 3 단계 소실(10-23)을 되풀이하지 않는다.
TAG = os.environ.get("DECOMP_TAG", "")
TAG_FITTER = {"": "F2", "_F1": "F1"}
RESULTS = pathlib.Path(__file__).resolve().parents[1] / "results"
COMMITTED = RESULTS / "committed_sha256.json"          # 커밋된 결과 파일 → sha256(results_manifest.py 가 쓴다)


def out_dir(stage: int) -> pathlib.Path:
    return RESULTS / f"stage{stage}{TAG}"


def store_path() -> pathlib.Path:
    """1 단계 방향 저장소(작업 폴더) — 갈래마다 따로. 정본은 옛 이름 그대로(directions_store.npz)."""
    return WORK / f"directions_store{TAG}.npz"


def guard_write(path) -> pathlib.Path:
    """⛔ 커밋된 결과 파일은 덮어쓰지 않는다 — COMMITTED 에 있는 경로면 멈춘다. 쓸 경로를 그대로 돌려준다."""
    import json
    p = pathlib.Path(path).resolve()
    done = json.loads(COMMITTED.read_text(encoding="utf-8"))["files"] if COMMITTED.exists() else {}
    rel = p.relative_to(RESULTS.parents[1]).as_posix() if p.is_relative_to(RESULTS.parents[1]) else str(p)
    if rel in done:
        raise SystemExit(f"⛔ 커밋된 결과를 덮어쓰려 했다: {rel} — 갈래(DECOMP_TAG)를 확인하라")
    return p

# ── 격자 (SPEC_v1 §6 = v0.2 §6) ────────────────────────────────────────────
TE_GRID = np.arange(60.0, 135.0 + 1e-9, 5.0)                 # 16 점
T2IN_GRID = (70.0, 82.0, 95.0)
R_GRID = (1.0, 1.15, 1.3, 1.45, 1.6, 1.8)
T2ISO_GRID = (1000.0, 2000.0)
T2_SETS = [(t, r, i) for t in T2IN_GRID for r in R_GRID for i in T2ISO_GRID]   # 36
T2_CENTRE = (82.0, 1.6, 1000.0)
T2_E2A = [T2_CENTRE, (82.0, 1.0, 1000.0), (82.0, 1.3, 1000.0), (82.0, 1.8, 1000.0)]
TR_GRID = (3600.0, 5800.0, 10620.0)                           # ms
T1_SETS = [(tt, ratio) for tt in (850.0, 1100.0) for ratio in (1 / 1.3, 1.0, 1.3)]   # (T1_tissue, T1_in/T1_en) 6
T1_ISO = 4300.0
NDI_GRID = (0.30, 0.45, 0.60, 0.75)
ODI_GRID = (0.08, 0.20, 0.40)
FISO_GRID = (0.02, 0.10, 0.25)
TISSUE = [(v, o, f) for v in NDI_GRID for o in ODI_GRID for f in FISO_GRID]   # 36 (NDI, ODI, FISO)
WM = [i for i, (v, o, f) in enumerate(TISSUE) if v >= 0.45 and o <= 0.20 and f <= 0.10]   # 12
DPAR_TRUE = (1.7, 2.0, 2.26, 2.5, 3.0)                        # µm²/ms
DPAR_FIT_F2 = (1.7, 2.0, 2.26, 2.5)
DISO = 3.0
SNR_N1 = (15.0, 20.0, 30.0, 50.0)
SNR_N2 = (20.0, 30.0, 50.0)
TE_REF = (70.0, 100.0)
R_REP = 200
GENERATORS = ("T-in", "T-app", "T-none")
NU_NONE = 0.60                                                # T-none 의 d⊥ = d∥(1 − 0.60)
PROTOCOLS = ("P1", "P2a", "P2b")
N2_REF = dict(f_in=0.5, T2_in=90.0, T2_en=60.0, f_iso=0.04, odi=0.20, TE=98.0, dpar=1.7)
CHUNK = int(os.environ.get("DECOMP_CHUNK", "400000"))                # 점검용 덮어쓰기(본 실행은 400,000)


def t2_triplet(t2set):
    """(T2_in, r, T2_iso) → (T2_in, T2_en, T2_iso)."""
    t2in, r, t2iso = t2set
    return t2in, t2in / r, t2iso


def t1_triplet(t1set):
    """(T1_tissue, T1_in/T1_en) → (T1_in, T1_en, T1_iso). 기하 평균이 T1_tissue."""
    tt, ratio = t1set
    return tt * math.sqrt(ratio), tt / math.sqrt(ratio), T1_ISO


def weights(TE, t2, TR=None, t1=None):
    """구획 가중 a_c = e^{−TE/T2_c} · (1 − e^{−TR/T1_c}) — (in, en, iso)."""
    a = np.exp(-TE / np.asarray(t2, float))
    if TR is not None:
        a = a * (1.0 - np.exp(-TR / np.asarray(t1, float)))
    return a


def apparent(nu, fiso, a):
    """겉보기 신호 분율 (ν′, f′_iso) — SPEC_v1 §1-4."""
    a_in, a_en, a_iso = a
    A = nu * a_in + (1 - nu) * a_en
    return nu * a_in / A, fiso * a_iso / (fiso * a_iso + (1 - fiso) * A)


def frac_coords(nu, fiso):
    """신호 분율의 (x, y) — x = logit ν, y = logit f_iso − ln ν."""
    return np.log(nu / (1 - nu)), np.log(fiso / (1 - fiso)) - np.log(nu)


def dperp(gen, dpar, nu, nu_app):
    if gen == "T-in":
        return dpar * (1 - nu)
    if gen == "T-app":
        return dpar * (1 - nu_app)
    if gen == "T-none":
        return dpar * (1 - NU_NONE)
    raise ValueError(gen)


def odi_to_kappa(odi):
    return 1.0 / math.tan(odi * math.pi / 2.0)


def fibonacci(n):
    k = np.arange(n) + 0.5
    z = 1 - 2 * k / n
    r = np.sqrt(1 - z * z)
    ph = math.pi * (1 + 5 ** 0.5) * k
    return np.c_[r * np.cos(ph), r * np.sin(ph), z]


FIBRE_DIRS = fibonacci(16)[fibonacci(16)[:, 2] >= 0][:8]          # 반구 8 방향(결정적)


def repulsion_dirs(n=30, iters=2000, step=0.01):
    """정전기 반발(대척점 대칭 쿨롱) — 피보나치 시작 · 경사 하강 · 결정적."""
    p = fibonacci(n)
    for _ in range(iters):
        g = np.zeros_like(p)
        for s in (1.0, -1.0):
            d = p[:, None, :] - s * p[None, :, :]
            r2 = (d ** 2).sum(-1)
            if s == 1.0:
                np.fill_diagonal(r2, np.inf)
            g += (d / r2[..., None] ** 1.5).sum(1)
        g -= (g * p).sum(1, keepdims=True) * p          # 접평면 성분만
        p = p + step * g / np.abs(g).max()
        p /= np.linalg.norm(p, axis=1, keepdims=True)
    return p


def protocol(name):
    """(bvals[V], bvecs[V,3]) — P1 은 생성, P2 는 메타 미러에서 복사한 파일(로그: prereg/artifacts/decomp/onharmony_access_log.json)."""
    if name == "P1":
        d = repulsion_dirs()
        bvecs = np.vstack([np.zeros((4, 3)), d, d])
        bvals = np.r_[np.zeros(4), np.full(30, 700.0), np.full(30, 2000.0)]
        return bvals, bvecs
    bvals = np.loadtxt(PROTO_DIR / f"{name}_bval").ravel()
    bvecs = np.loadtxt(PROTO_DIR / f"{name}_bvec")
    bvecs = bvecs.T if bvecs.shape[0] == 3 else bvecs
    return bvals, bvecs


def low_shell(name):
    return 700.0 if name == "P1" else 1000.0


# ── 실험 행 (SPEC_v1 §7) ───────────────────────────────────────────────────
# 행 = 합성 신호 하나. 열: exp, prot, gen, TE, TR(없으면 nan), t2(3), t1(3 · 없으면 nan),
#      tissue, fdir, dpar_true, noise('N0'|'N1'|'N2'|'L4'|'L8'), snr, rep
FIELDS = ("prot", "gen", "TE", "TR", "t2set", "t1set", "tissue", "fdir", "dpar", "noise", "snr", "rep")


def cells(exp):
    """(행 생성기) — 반복(rep)을 펼치기 전의 칸 목록. 각 칸은 dict."""
    nt, nd = range(len(TISSUE)), range(len(FIBRE_DIRS))
    if exp == "E1":
        for p, g, te, t2, ti, fd in itertools.product(PROTOCOLS, GENERATORS, TE_GRID, range(36), nt, nd):
            yield dict(prot=p, gen=g, TE=te, TR=np.nan, t2set=t2, t1set=-1, tissue=ti, fdir=fd, dpar=1.7, noise="N0", snr=np.inf, nrep=1)
    elif exp == "E1-TR":
        c = T2_SETS.index(T2_CENTRE)
        for p, g, tr, t1, ti, fd in itertools.product(PROTOCOLS, GENERATORS, TR_GRID, range(6), nt, nd):
            yield dict(prot=p, gen=g, TE=92.0, TR=tr, t2set=c, t1set=t1, tissue=ti, fdir=fd, dpar=1.7, noise="N0", snr=np.inf, nrep=1)
    elif exp == "E3":
        c = T2_SETS.index(T2_CENTRE)
        for p, g, dp, te, ti, fd in itertools.product(PROTOCOLS, GENERATORS, DPAR_TRUE, TE_GRID, nt, nd):
            yield dict(prot=p, gen=g, TE=te, TR=np.nan, t2set=c, t1set=-1, tissue=ti, fdir=fd, dpar=dp, noise="N0", snr=np.inf, nrep=1)
    elif exp == "E2a":
        for p, g, s, t2, te, ti, fd in itertools.product(PROTOCOLS, GENERATORS, SNR_N1, T2_E2A, TE_REF, nt, nd):
            yield dict(prot=p, gen=g, TE=te, TR=np.nan, t2set=T2_SETS.index(t2), t1set=-1, tissue=ti, fdir=fd,
                       dpar=1.7, noise="N1", snr=s, nrep=R_REP)
    elif exp == "E2c":
        c = T2_SETS.index(T2_CENTRE)
        for p, L, te, ti, fd in itertools.product(("P2a", "P2b"), ("L4", "L8"), TE_REF, nt, nd):
            yield dict(prot=p, gen="T-in", TE=te, TR=np.nan, t2set=c, t1set=-1, tissue=ti, fdir=fd, dpar=1.7, noise=L, snr=30.0, nrep=R_REP)
    elif exp == "E5":
        for acq, (p, te, tr) in ACQ_E5.items():
            for g, t2, t1, ti, fd in itertools.product(GENERATORS, ETA_T2, range(6), range(len(TISSUE_E5)), FDIR_E5):
                yield dict(prot=p, gen=g, TE=te, TR=tr, t2set=T2_SETS.index(t2), t1set=t1, tissue=ti, fdir=fd, dpar=1.7,
                           noise="N0", snr=np.inf, nrep=1, dense=True, acq=acq)
    elif exp == "E6":
        for acq, (p, te, tr) in ACQ_E5.items():
            for g, e, ti, fd in itertools.product(GENERATORS, range(len(ETA_E6)), range(len(TISSUE_E5)), FDIR_E5):
                yield dict(prot=p, gen=g, TE=te, TR=tr, t2set=-1, t1set=-1, eta6=e, tissue=ti, fdir=fd, dpar=1.7,
                           noise="N0", snr=np.inf, nrep=1, dense=True, acq=acq)
    elif exp == "E6n":
        etas = [n for n, _, _ in ETA_E6]
        for g, (ra, rb, ea, eb, case) in itertools.product(GENERATORS, e6n_cases()):
            for acq, (p, te, tr) in ACQ_E5.items():
                for side, ti, en in (("A", ra, ea), ("B", rb, eb)):
                    if side == "A" and acq != "A":
                        continue
                    yield dict(prot=p, gen=g, TE=te, TR=tr, t2set=-1, t1set=-1, eta6=etas.index(en), tissue=ti, fdir=ra % 4 * 2,
                               dpar=1.7, noise="N1", snr=30.0, nrep=10, dense=True, acq=acq, side=side, case=case, ref=ra)
    elif exp == "E2b":
        c = T2_SETS.index(T2_CENTRE)
        for p, s, te, ti, fd in itertools.product(PROTOCOLS, SNR_N2, TE_GRID, nt, nd):
            yield dict(prot=p, gen="T-in", TE=te, TR=np.nan, t2set=c, t1set=-1, tissue=ti, fdir=fd, dpar=1.7, noise="N2", snr=s, nrep=R_REP)
    else:
        raise ValueError(exp)


def t2_of(idx):
    """칸의 t2set 색인 → (T2_in, r, T2_iso). E2a 의 네 집합도 모두 T2_SETS 안에 있다."""
    return T2_SETS[idx]


def n_rows(exp):
    return sum(c["nrep"] for c in cells(exp))


def chunk_rows(exp, prot):
    """(칸 목록, 덩어리별 (칸 시작, 칸 끝)) — 반복은 칸 안에 붙어 있다. 덩어리 = CHUNK 행 이하."""
    cl = [c for c in cells(exp) if c["prot"] == prot]
    bounds, start, acc = [], 0, 0
    for i, c in enumerate(cl):
        if acc + c["nrep"] > CHUNK and acc > 0:
            bounds.append((start, i))
            start, acc = i, 0
        acc += c["nrep"]
    bounds.append((start, len(cl)))
    return cl, bounds


# ── 2 단계(SPEC_v2 §11) — 조밀 격자 E5 ─────────────────────────────────────
N_DENSE = 26                                                  # SPEC_v2 §4: 9,720·N² ≤ 6,928,272 인 가장 큰 N
NDI_DENSE = tuple(float(v) for v in np.linspace(0.30, 0.80, N_DENSE))
FISO_DENSE = tuple(float(v) for v in np.linspace(0.01, 0.30, N_DENSE))
TISSUE_E5 = [(v, o, f) for v in NDI_DENSE for o in ODI_GRID for f in FISO_DENSE]       # 2,028 (NDI, ODI, FISO)
ETA_T2 = [(t, r, 1000.0) for t in T2IN_GRID for r in (1.3, 1.6, 1.8)]                  # 9
FDIR_E5 = (0, 2, 4, 6)
ACQ_E5 = {"A": ("P2a", 92.0, 3600.0), "Q1": ("P2a", 68.0, 3600.0), "Q2": ("P2b", 70.6, 10215.0),
          "Q3": ("P2b", 92.0, 3600.0), "Q4": ("P2a", 92.0, 10215.0)}                    # LOCK_B EN :7748–7754


# ── 3 단계(SPEC_v3 §2 · §9) — 사람 사이 이완 변이 축 · E6 · E6n ──────────────
T1_TISSUE_C = 1100.0                                          # SPEC_v3 §2 중심 T1_tissue(2 단계 기준 η 와 같게)
SCEN = {"i": 1.3, "ii": 1 / 1.3, "iii": 1.0}                  # T1_in/T1_en 시나리오(SPEC_v3 §1-3)
T2C = (82.0, 82.0 / 1.6, 1000.0)                              # 중심 (T2_in, T2_en, T2_iso)


def _t1(ratio, scale_t=1.0, scale_iso=1.0, en_only=False):
    tin, ten = T1_TISSUE_C * math.sqrt(ratio), T1_TISSUE_C / math.sqrt(ratio)
    if en_only:
        return (tin, ten * scale_t, T1_ISO * scale_iso)
    return (tin * scale_t, ten * scale_t, T1_ISO * scale_iso)


def _l1_points(cv, rho, scen="iii", sens=False):
    """SPEC_v3 §4 — 표준 정규 공간의 반지름 {1, 2} × 각 8(0, 45, …, 315°) = 16 점 → ρ 의 콜레스키 → (T2_in, T2_en).
    점 k 의 조직 T1(두 구획 같이) 배율 1 + 0.04·{−2, 0, +2}[k % 3]. T1_iso 는 고정 — sens(L1s)에서만 배율 1 + 0.08·{0, +2, −2}[k % 3]."""
    L = np.linalg.cholesky(np.array([[1.0, rho], [rho, 1.0]]) + 1e-12 * np.eye(2))
    pts = []
    for k, (rad, ang) in enumerate([(r_, a_) for r_ in (1.0, 2.0) for a_ in range(0, 360, 45)]):
        u = np.array([rad * math.cos(math.radians(ang)), rad * math.sin(math.radians(ang))])
        z = L @ u
        st = (-2, 0, 2)[k % 3]; si = (0, 2, -2)[k % 3]
        pts.append(((T2C[0] * (1 + cv * z[0]), T2C[1] * (1 + cv * z[1]), T2C[2]), _t1(SCEN[scen], 1 + 0.04 * st, (1 + 0.08 * si) if sens else 1.0)))
    return pts


def eta_e6():
    """η 199 개 — (이름, (T2_in, T2_en, T2_iso), (T1_in, T1_en, T1_iso)). SPEC_v3 §9."""
    out = [(f"L0|{s}", T2C, _t1(SCEN[s])) for s in ("i", "ii", "iii")]
    for cv in (0.015, 0.03, 0.04, 0.06):
        for rho in (0.0, 0.8):
            out += [(f"L1|cv{cv}|rho{rho}|iii|{k}", t2, t1) for k, (t2, t1) in enumerate(_l1_points(cv, rho, "iii"))]
    for s in ("i", "ii"):
        out += [(f"L1|cv0.03|rho0.0|{s}|{k}", t2, t1) for k, (t2, t1) in enumerate(_l1_points(0.03, 0.0, s))]
    out += [(f"L1s|cv0.03|rho0.0|iii|{k}", t2, t1) for k, (t2, t1) in enumerate(_l1_points(0.03, 0.0, "iii", sens=True))]
    out.append(("L2|a", (T2C[0] * 1.0984, T2C[1] * 1.0984, T2C[2]), _t1(1.0, 1.0441)))
    out.append(("L2|b", (T2C[0], T2C[1] * 1.0984, T2C[2]), _t1(1.0, 1.0441, en_only=True)))
    for me in (2.1, 2.9, 3.4):                                # SPEC_v3 §2 L3 — Lampinen 2019 · 2020 에서 [유도]
        for mi in (0.77, 0.87, 1.02):
            for m1 in (1.17, 1.56):
                out.append((f"L3|en{me}|in{mi}|t1{m1}", (T2C[0] * mi, T2C[1] * me, T2C[2]), _t1(1.0, m1)))
    return out


ETA_E6 = eta_e6()                                             # 199
N_REFS_NOTE = "기준 조직 = TISSUE_E5 중 NDI ∈ [0.45, 0.75] · FISO ∈ [0.02, 0.10] · ODI ∈ {0.08, 0.20}"


def refs_e5():
    return [i for i, (v, o, f) in enumerate(TISSUE_E5) if 0.45 <= v <= 0.75 and 0.02 <= f <= 0.10 and o in (0.08, 0.20)]


def e6n_cases():
    """SPEC_v3 §5 — 기준 조직마다 경우 6: L0 Δ0 · L0 + δ(NDI ± 2 칸 · FISO ± 2 칸) · L1(CV 3 % · ρ 0 · iii, 점 = 기준 색인 % 16)."""
    nN, nF = len(NDI_DENSE), len(FISO_DENSE)
    out = []
    for r in refs_e5():
        iv, io, jf = r // (3 * nF), (r // nF) % 3, r % nF
        out.append((r, r, "L0|iii", "L0|iii", "L0d0"))
        for dv, df, nm in ((2, 0, "dNDI+2"), (-2, 0, "dNDI-2"), (0, 2, "dFISO+2"), (0, -2, "dFISO-2")):
            v2, f2 = iv + dv, jf + df
            if 0 <= v2 < nN and 0 <= f2 < nF:
                out.append((r, (v2 * 3 + io) * nF + f2, "L0|iii", "L0|iii", nm))
        out.append((r, r, "L0|iii", f"L1|cv0.03|rho0.0|iii|{r % 16}", "L1d0"))
    return out


EXPERIMENTS = ("E1", "E1-TR", "E3", "E2a", "E2c", "E2b")      # SPEC_v1 §7 의 순서(E4 는 재표집)
F2_DPARS = {"E3": DPAR_FIT_F2}                                # 그 밖은 1.7
