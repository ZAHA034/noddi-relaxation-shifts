"""decomp 1 단계 분석 — SPEC_v1 §8–§11 · §16 · §17. 메인 환경에서 돈다.

사용:  python analyze.py directions   → 방향 · SG1–SG6 · Π 지도 · §16 → results/stage1/directions.json(+ npz)
       python analyze.py e4 [nproc]   → E4(SG7) → results/stage1/e4.json
       python analyze.py report       → results/stage1/summary.json · REPORT.md
구현 선택(결과 전 고정 — SPEC_v1 에 없는 세부): 이 파일 머리의 IMPL 목록.
"""
from __future__ import annotations

import glob
import json
import math
import os
import pathlib
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C  # noqa: E402
import lib  # noqa: E402

REPO = pathlib.Path(__file__).resolve().parents[2]
OUT = C.out_dir(1)                     # 정본 results/stage1 · F1 판 results/stage1_F1(중계 ZZZ-3 나)
FITS = C.WORK / "fits"
FITTERS = tuple(os.environ.get("DECOMP_FITTERS", "F2").split(","))   # 기본 F2 만(CPU) — F1 은 GPU 허락 뒤 "F1,F2"
if C.TAG not in C.TAG_FITTER or FITTERS != (C.TAG_FITTER[C.TAG],):
    # 갈래 하나 = 적합기 하나(중계 ZZZ-3 나) — 둘을 한 번에 돌리면 E4 시드가 밀려 정본 F2 값이 바뀐다
    raise SystemExit(f"⛔ 갈래 DECOMP_TAG={C.TAG!r} 와 적합기 {FITTERS} 가 맞지 않는다 — {C.TAG_FITTER}")
STAGES = ("C3", "C5", "C7")

IMPL = [
    "u_TE 처럼 TE 전 구간에 걸친 양은 Σ 를 TE_ref 100(격자 평균 97.5 에 가장 가까움)에서 고른다.",
    "SG1 은 섬유 방향별 칸에서, 할선과 그 칸의 u_TE 사이 각을 칸의 상태 Σ 로 잰다.",
    "SG2 의 상태 간 각은 WM 12 상태 Σ 의 평균(같은 프로토콜 · 중심 T2 · TE_ref 100)으로, 프로토콜 간 각은 C3 에서 세 프로토콜 Σ 평균으로(C5 · C7 은 P2a · P2b 만 — 셸 구성이 같은 것끼리).",
    "SG4 의 Σ 는 쌍 두 프로토콜 Σ 의 평균(중심 T2 · TE_ref 100). 기본 적합기 판정은 F1 · F2 모두 내고 F2 는 적합 d∥ 1.7.",
    "u_SNR(s) = (잡음 칸의 반복 · 섬유 방향 평균) − (E1 의 같은 TE · T2 · 조직 무잡음 적합의 섬유 방향 평균).",
    "범주 대표: N-SNR 은 SNR 네 수준(중심 T2 · TE_ref 100), N-M×P 는 참 d∥ 2.0 · 2.26 · 2.5(쌍 A: (P2a, 100) 대 (P2b, 70)).",
    "Π 지도의 쌍 유형: T1 TE 만 {N-TE, N-SNR} · T2 설계만 {N-P, N-M×P, N-SNR} · T3 TE + 설계 {N-TE, N-SNR, N-P, N-M×P} · T4 다섯 전부 · T5 SNR 만 {N-SNR}.",
    "E4 의 U 원 좌표 벡터: SG2 통과면 WM 상태 평균(중심 T2), 실패면 쌍 평균 좌표에 가장 가까운 WM 상태(백색화 거리) — (적합기, 생성기 U) 별 SG2(C3) 결과로 정한다.",
    "E4 의 참 W_d · 참 Δ: 피험자 5,000 명 표본에서. 대립의 스파이크 방향 = 스캐너 B 프로토콜 · WM 평균 t_NDI(자료 생성기).",
]


# ── 적재 ─────────────────────────────────────────────────────────────────
def load(exp, prot, fitter, dpar_fit=1.7):
    """(행 열 dict, 좌표 C7 (N, 7) float32, 실패 표시 dict)."""
    rows, ys = [], []
    for d in sorted(glob.glob(str(FITS / exp / prot / "c*"))):
        if not os.path.exists(f"{d}/DONE"):
            raise SystemExit(f"⛔ 끝나지 않은 덩어리 {d}")
        r = dict(np.load(f"{d}/rows.npz", allow_pickle=False))
        f = np.load(f"{d}/f1.npz" if fitter == "F1" else f"{d}/f2_{dpar_fit}.npz")
        y, _ = lib.coords(f["ndi"].astype(float), f["fwf"].astype(float), f["odi"].astype(float), r["ext"].astype(float), "C7")
        r["ndi"], r["fwf"], r["odi"] = f["ndi"], f["fwf"], f["odi"]
        rows.append(r)
        ys.append(y.astype(np.float32))
    cols = {k: np.concatenate([r[k] for r in rows]) for k in rows[0] if np.ndim(rows[0][k]) >= 1}
    return cols, np.concatenate(ys)


def grid_index(cols, keys, sizes):
    """행 → 격자 색인 튜플."""
    return tuple(np.asarray(cols[k]) for k in keys)


def to_array(exp, fitter, dpar_fit=1.7):
    """무잡음 묶음을 격자 배열로. E1: (gen, prot, TE, t2, tissue, fdir, 7) · E1-TR: (gen, prot, TR, t1, tissue, fdir, 7) ·
    E3: (gen, prot, dpar, TE, tissue, fdir, 7). 실패 표시 배열도 같은 꼴."""
    if exp == "E1":
        shape = (3, 3, 16, 36, 36, 8)
    elif exp == "E1-TR":
        shape = (3, 3, 3, 6, 36, 8)
    else:
        shape = (3, 3, 5, 16, 36, 8)
    Y = np.full(shape + (7,), np.nan, np.float32)
    FL = {k: np.zeros(shape, bool) for k in ("sat", "rec")}
    APP = np.full(shape + (2,), np.nan, np.float32)
    for pi, prot in enumerate(C.PROTOCOLS):
        cols, y = load(exp, prot, fitter, dpar_fit)
        g = np.array([C.GENERATORS.index(s) for s in cols["gen"]])
        if exp == "E1":
            a = np.searchsorted(C.TE_GRID, cols["TE"]); b = cols["t2set"]
        elif exp == "E1-TR":
            a = np.searchsorted(C.TR_GRID, cols["TR"]); b = cols["t1set"]
        else:
            a = np.array([C.DPAR_TRUE.index(round(float(v), 2)) for v in cols["dpar"]]); b = np.searchsorted(C.TE_GRID, cols["TE"])
        idx = (g, np.full_like(g, pi), a, b, cols["tissue"], cols["fdir"])
        Y[idx] = y
        ndi, fwf, odi = cols["ndi"], cols["fwf"], cols["odi"]
        FL["sat"][idx] = (ndi >= 0.98) | (ndi <= 0.02) | (fwf >= 0.98) | (fwf <= 2e-4) | (odi <= ODI_MIN)   # §183-2 1
        FL["rec"][idx] = (np.abs(ndi - cols["nu_app"]) > 0.1) | (np.abs(fwf - cols["f_app"]) > 0.1)
        APP[idx] = np.c_[cols["nu_app"], cols["f_app"]]
    return Y, FL, APP


def e2a_array(fitter, gen):
    """E2a: (prot, snr4, t2e4, teref2, tissue, fdir, rep, 7) float32 — 한 생성기."""
    out = np.full((3, 4, 4, 2, 36, 8, C.R_REP, 7), np.nan, np.float32)
    t2e = [C.T2_SETS.index(t) for t in C.T2_E2A]
    for pi, prot in enumerate(C.PROTOCOLS):
        cols, y = load("E2a", prot, fitter)
        m = cols["gen"] == gen
        s = np.array([C.SNR_N1.index(v) for v in cols["snr"][m]])
        t = np.array([t2e.index(v) for v in cols["t2set"][m]])
        te = np.array([C.TE_REF.index(v) for v in cols["TE"][m]])
        out[pi, s, t, te, cols["tissue"][m], cols["fdir"][m], cols["rep"][m]] = y[m]
    return out


# ── 방향 계산 ─────────────────────────────────────────────────────────────
ODI_MIN = 2 / math.pi * math.atan(1 / (64 - 2e-4))            # κ ≥ 64 − 2e-4 ⇔ ODI ≤ 이 값
TE90 = int(np.where(C.TE_GRID == 90.0)[0][0])
TE65 = int(np.where(C.TE_GRID == 65.0)[0][0])
CENTRE = C.T2_SETS.index(C.T2_CENTRE)
T2E = [C.T2_SETS.index(t) for t in C.T2_E2A]
R_OF = [C.T2_SETS[i][1] for i in range(36)]


def nearest_t2e(t2set):
    """§8 Σ 선택: r 로만 · 동률이면 큰 쪽."""
    r = R_OF[t2set]
    rs = [t[1] for t in C.T2_E2A]
    d = [abs(r - x) for x in rs]
    best = min(d)
    cand = [i for i, x in enumerate(d) if abs(x - best) < 1e-12]
    return max(cand, key=lambda i: rs[i])


def nearest_teref(te):
    d = [abs(te - x) for x in C.TE_REF]
    return 1 if abs(d[0] - d[1]) < 1e-12 else int(np.argmin(d))


def sigma_table(E2):
    """E2a SNR 30 의 공분산: (prot, t2e, teref, tissue, 7, 7) — 반복 공분산을 섬유 방향 평균."""
    s30 = C.SNR_N1.index(30.0)
    X = E2[:, s30].astype(float)                                   # (prot, t2e, teref, tissue, fdir, rep, 7)
    Xc = X - X.mean(-2, keepdims=True)
    S = np.einsum("...ri,...rj->...ij", Xc, Xc) / (X.shape[-2] - 1)
    return S.mean(-3)


def sub(S, p):
    return S[..., :p, :p]


def directions(fitter, gen, Y1, YTR, Y3, E2, Sig):
    gi = C.GENERATORS.index(gen)
    out = {}
    m1 = np.nanmean(Y1[gi], axis=-2).astype(float)                  # (prot, TE, t2, tissue, 7)
    gTE = np.stack([[[[lib.ls_slope(C.TE_GRID, m1[p, :, t, s]) for s in range(36)] for t in range(36)]] for p in range(3)])[:, 0]
    out["gTE"] = gTE                                                # (prot, t2, tissue, 7) — 원 좌표 · /ms
    mTR = np.nanmean(YTR[gi], axis=-2).astype(float)               # (prot, TR3, t1 6, tissue, 7)
    trs = np.array(C.TR_GRID)
    out["gTR"] = np.stack([[[lib.ls_slope(trs, mTR[p, :, a, s]) for s in range(36)] for a in range(6)] for p in range(3)])
    out["trsec"] = np.stack([mTR[:, 1] - mTR[:, 0], mTR[:, 2] - mTR[:, 1]], 0)      # (2, prot, t1, tissue, 7)
    m3 = np.nanmean(Y3[gi], axis=-2).astype(float)                  # (prot, dpar5, TE, tissue, 7)
    tA, tB = list(C.TE_GRID).index(100.0), list(C.TE_GRID).index(70.0)
    d17 = C.DPAR_TRUE.index(1.7)
    mxp = {}
    for nm, (ta, tb) in {"A": (tA, tB), "B": (tA, tA)}.items():
        base = m3[2, d17, tb] - m3[1, d17, ta]
        mxp[nm] = np.stack([(m3[2, C.DPAR_TRUE.index(d), tb] - m3[1, C.DPAR_TRUE.index(d), ta]) - base for d in (2.0, 2.26, 2.5)], 1)  # (tissue, 3, 7)
    out["gMxP"] = mxp
    out["gP"] = m1[2, TE90] - m1[1, TE90]                           # (t2, tissue, 7)
    # u_SNR(s): (prot, snr, t2e, teref, tissue, 7)
    noisy = np.nanmean(E2, axis=(-3, -2)).astype(float)             # (prot, snr, t2e, teref, tissue, 7)
    clean = np.stack([[m1[p, list(C.TE_GRID).index(te), T2E[t]] for te in C.TE_REF] for p in range(3) for t in range(4)]).reshape(3, 4, 2, 36, 7)
    out["gSNR"] = noisy - clean[:, None]
    # 조직 방향 t_j (prot, tissue, 3, 7) — TE 90 · 중심 T2
    tj = np.zeros((3, 36, 3, 7))
    for s, (v, o, f) in enumerate(C.TISSUE):
        for j, (grid, val) in enumerate(((C.NDI_GRID, v), (C.ODI_GRID, o), (C.FISO_GRID, f))):
            i = grid.index(val)
            i2 = i + 1 if i + 1 < len(grid) else i - 1
            tv = [v, o, f]; tv[j] = grid[i2]
            s2 = C.TISSUE.index(tuple(tv))
            tj[:, s, j] = (m1[:, TE90, CENTRE, s2] - m1[:, TE90, CENTRE, s]) / (grid[i2] - grid[i])
    out["tj"] = tj
    out["m1"] = m1
    return out


def sig_for(Sig, prot, t2set, tissue, te=100.0, p=7):
    return sub(Sig[prot, nearest_t2e(t2set), nearest_teref(te), tissue], p)


def sg_all(fitter, gen, D, Sig, Y1):
    """SG1–SG6 · Π 지도 · §16 — 좌표 단계별."""
    gi = C.GENERATORS.index(gen)
    res = {}
    for st in STAGES:
        p = lib.STAGE_P[st]
        R = {}
        # SG1 — 섬유 방향별 칸 · 할선과 칸의 u_TE (백색화 뒤 보통 각)
        a1 = []
        for pr in range(3):
            for t in range(36):
                for s in range(36):
                    W12 = lib.inv_sqrt(sig_for(Sig, pr, t, s, 100.0, p))
                    y = Y1[gi, pr, :, t, s, :, :p].astype(float)             # (TE, fdir, p)
                    if not np.isfinite(y).all():
                        continue
                    yw = np.einsum("ij,tfj->tfi", W12, y)
                    tc = C.TE_GRID - C.TE_GRID.mean()
                    u = np.einsum("t,tfi->fi", tc, yw - yw.mean(0))           # (fdir, p) — 기울기 방향
                    sec = np.diff(yw, axis=0)                                  # (15, fdir, p)
                    c = np.abs(np.einsum("sfi,fi->sf", sec, u)) / (np.linalg.norm(sec, axis=-1) * np.linalg.norm(u, axis=-1) + 1e-300)
                    a1 += list(np.degrees(np.arccos(np.clip(c, 0, 1))).max(0))
        a1 = np.array(a1)
        R["SG1"] = dict(value=float((a1 <= 5).mean()), pass_=bool((a1 <= 5).mean() >= 0.95), n=len(a1),
                        alt={"3°": float((a1 <= 3).mean()), "10°": float((a1 <= 10).mean())}, p95=float(np.percentile(a1, 95)))
        # SG2
        angs, meds = [], {}
        for pr in range(3):
            Sp = np.mean([sig_for(Sig, pr, CENTRE, s, 100.0, p) for s in C.WM], 0)
            W12 = lib.inv_sqrt(Sp)
            V = np.stack([W12 @ D["gTE"][pr, CENTRE, s, :p] for s in C.WM])
            V /= np.linalg.norm(V, axis=1, keepdims=True)
            V *= np.sign(V @ V[0])[:, None]
            m = V.mean(0); m /= np.linalg.norm(m)
            angs += [math.degrees(math.acos(min(1, abs(v @ m)))) for v in V]
            meds[pr] = (np.linalg.inv(W12) @ m)
        prots = (0, 1, 2) if st == "C3" else (1, 2)
        Spool = np.mean([np.mean([sig_for(Sig, pr, CENTRE, s, 100.0, p) for s in C.WM], 0) for pr in prots], 0)
        pa = max(lib.wangle(meds[a], meds[b], Spool) for a in prots for b in prots if a < b)
        q95 = float(np.percentile(angs, 95))
        R["SG2"] = dict(state_p95=q95, protocol_max=pa, pass_=bool(q95 <= 5 and pa <= 5),
                        alt={"3°": bool(q95 <= 3 and pa <= 3), "10°": bool(q95 <= 10 and pa <= 10)})
        # SG3
        a3 = []
        for pr in range(3):
            for r in (1.3, 1.6, 1.8):
                t2 = C.T2_SETS.index((82.0, r, 1000.0)); te_i = [C.T2_E2A.index(x) for x in C.T2_E2A if x[1] == r][0]
                for k, te in enumerate(C.TE_REF):
                    for s in range(36):
                        S = sub(Sig[pr, te_i, k, s], p)
                        for sn in range(4):
                            a3.append(lib.wangle(D["gTE"][pr, t2, s, :p], D["gSNR"][pr, sn, te_i, k, s, :p], S))
        a3 = np.array(a3)
        R["SG3"] = dict(value=float((a3 >= 30).mean()), pass_=bool((a3 >= 30).mean() >= 0.9), n=len(a3), median=float(np.median(a3)),
                        alt={"20°": float((a3 >= 20).mean()), "45°": float((a3 >= 45).mean())})
        # 범주 대표 · SG4 · Π
        a4, PI = [], {}
        for s in range(36):
            for pr in (1, 2):
                S = sig_for(Sig, pr, CENTRE, s, 100.0, p)
                W12 = lib.inv_sqrt(S)
                g = {"N-TE": D["gTE"][pr, CENTRE, s, :p][:, None],
                     "N-SNR": D["gSNR"][pr, :, 0, 1, s, :p].T,
                     "N-P": D["gP"][CENTRE, s, :p][:, None],
                     "N-MxP": D["gMxP"]["A"][s, :, :p].T,
                     "N-TR": D["gTR"][pr, :, s, :p].T}
                rep = {k: lib.group_basis([v], W12) for k, v in g.items()}
                T = W12 @ D["tj"][pr, s, :, :p].T
                types = {"T1": ("N-TE", "N-SNR"), "T2": ("N-P", "N-MxP", "N-SNR"), "T3": ("N-TE", "N-SNR", "N-P", "N-MxP"),
                         "T4": ("N-TE", "N-TR", "N-SNR", "N-P", "N-MxP"), "T5": ("N-SNR",)}
                for tn, cats in types.items():
                    Ub = np.concatenate([rep[c] for c in cats], 1)
                    Ub = Ub[:, np.linalg.norm(Ub, axis=0) > 0]
                    k, Q, P, sv = lib.rank_projector(Ub)
                    Pi, Pij, pa_ = lib.preservation(P, T) if p - k >= 1 else (0.0, np.zeros(3), np.zeros(0))
                    if p - k < 3:
                        Pi = 0.0
                    PI.setdefault(tn, []).append(dict(prot=C.PROTOCOLS[pr], tissue=s, wm=s in C.WM, k=k, Pi=Pi, Pij=[float(x) for x in Pij],
                                                      sv=[float(x) for x in sv]))
                if pr == 1:
                    Sm = (S + sig_for(Sig, 2, CENTRE, s, 100.0, p)) / 2
                    W12m = lib.inv_sqrt(Sm)
                    um = lib.group_basis([g["N-MxP"]], W12m)[:, :1]
                    rest = np.concatenate([lib.group_basis([g[c]], W12m)[:, :1] for c in ("N-TE", "N-SNR", "N-P")], 1)
                    a4.append(lib.span_angle(um[:, 0], rest))
        a4 = np.array(a4)
        R["SG4"] = dict(value=float((a4 >= 30).mean()), pass_=bool((a4 >= 30).mean() >= 0.9), n=len(a4),
                        alt={"20°": float((a4 >= 20).mean()), "45°": float((a4 >= 45).mean())})
        # SG5 — r 곡선(중앙 1.6 대비)
        sg5 = {}
        for r in C.R_GRID:
            t2r, t2c = C.T2_SETS.index((82.0, r, 1000.0)), CENTRE
            sg5[str(r)] = float(np.median([lib.wangle(D["gTE"][pr, t2r, s, :p], D["gTE"][pr, t2c, s, :p], sig_for(Sig, pr, CENTRE, s, 100.0, p))
                                           for pr in (1, 2) for s in C.WM]))
        R["SG5"] = sg5
        # SG6 — 선택 규칙 U = T3 · WM · P2a/P2b
        t3 = [d for d in PI["T3"] if d["wm"]]
        frac = float(np.mean([d["Pi"] >= lib.PI_STAR for d in t3]))
        R["SG6"] = dict(value=frac, pass_=bool(frac >= 0.9), alt={"0.3": float(np.mean([d["Pi"] >= 0.3 for d in t3])),
                                                                  "0.7": float(np.mean([d["Pi"] >= 0.7 for d in t3]))})
        R["Pi_map"] = {tn: dict(frac_below=float(np.mean([d["Pi"] < lib.PI_STAR for d in v if d["wm"]])),
                                median_Pi=float(np.median([d["Pi"] for d in v if d["wm"]])),
                                median_Pij=[float(np.median([d["Pij"][j] for d in v if d["wm"]])) for j in range(3)],
                                k=[int(x) for x in np.unique([d["k"] for d in v])]) for tn, v in PI.items()}
        # §16 (i)(ii) — T4 · 전 조직 · 중심 T2
        R["S16"] = dict(fiso_below_0p1=float(np.mean([d["Pij"][2] < 0.1 for d in PI["T4"]])),
                        ndi_median=float(np.median([d["Pij"][0] for d in PI["T4"]])),
                        ndi_median_T3=float(np.median([d["Pij"][0] for d in PI["T3"]])))
        res[st] = R
    return res


def run_directions():
    OUT.mkdir(parents=True, exist_ok=True)
    allres, store = {}, {}
    for F in FITTERS:
        Y1, FL1, APP1 = to_array("E1", F)
        YTR, _, _ = to_array("E1-TR", F)
        Y3, FL3, _ = to_array("E3", F)
        for gen in C.GENERATORS:
            E2 = e2a_array(F, gen)
            Sig = sigma_table(E2)
            D = directions(F, gen, Y1, YTR, Y3, E2, Sig)
            r = sg_all(F, gen, D, Sig, Y1)
            gi = C.GENERATORS.index(gen)
            # 실패 표시 · §16 (iii)
            x_fit = Y1[gi, ..., 0].astype(float)
            x_app = np.log(APP1[gi, ..., 0] / (1 - APP1[gi, ..., 0]))
            shift = np.abs(Y1[gi, :, TE65, ..., 0] - Y1[gi, :, TE90, ..., 0])
            r["flags"] = dict(saturation=float(FL1["sat"][gi].mean()), recovery=float(FL1["rec"][gi].mean()),
                              shift_outlier=float((shift > 0.5).mean()), E3_saturation=float(FL3["sat"][gi].mean()),
                              E3_recovery=float(FL3["rec"][gi].mean()))
            r["S16_iii_median_abs_x_minus_logit_app"] = float(np.nanmedian(np.abs(x_fit - x_app)))
            allres[f"{F}|{gen}"] = r
            store[f"{F}|{gen}"] = {k: v for k, v in D.items() if k in ("gTE", "gTR", "gP", "gSNR", "tj", "m1")} | {"gMxP_A": D["gMxP"]["A"], "gMxP_B": D["gMxP"]["B"], "Sig": Sig}
            print(F, gen, {st: {k: r[st][k]["pass_"] for k in ("SG1", "SG2", "SG3", "SG4", "SG6")} for st in STAGES}, flush=True)
    json.dump(dict(impl=IMPL, results=allres), open(C.guard_write(OUT / "directions.json"), "w"), ensure_ascii=False, indent=1, default=float)
    np.savez_compressed(C.store_path(), **{f"{k}::{kk}": vv for k, v in store.items() for kk, vv in v.items()})


# ── E4 ───────────────────────────────────────────────────────────────────
PAIRS = {"N1": ((1, 100.0, 30.0), (1, 70.0, 30.0)), "N2": ((1, 100.0, 30.0), (2, 100.0, 20.0)),
         "N3": ((1, 100.0, 30.0), (2, 70.0, 20.0)), "N4": ((1, 100.0, 30.0), (1, 100.0, 50.0))}
LAMBDAS = (0.0, 5.0, 20.0, 80.0)
M_SETS = 1000


def _state_idx(st):
    pr, te, snr = st
    return pr, C.SNR_N1.index(snr), C.TE_REF.index(te)


def e4_groups(store, F, genU, pair, p, W12_true, Dmean):
    """쌍의 범주 열 묶음(원 좌표) — §9-7 포함 규칙 + C1."""
    S = store[f"{F}|{genU}"]
    (pa, tea, _), (pb, teb, _) = PAIRS[pair]
    dte = teb - tea
    groups = []
    wm = C.WM
    # 작동점 — IMPL 8
    if S.get("sg2_pass", True):
        pick = lambda arr: arr[wm].mean(0)
    else:
        m1 = S["m1"][pa, list(C.TE_GRID).index(tea), CENTRE][:, :p]
        dist = [np.linalg.norm(W12_true @ (m1[s] - Dmean)) for s in wm]
        sbest = wm[int(np.argmin(dist))]
        pick = lambda arr: arr[sbest]
    if dte != 0:
        g = pick(S["gTE"][pa, CENTRE][:, :p]) * dte
        mk = max(np.linalg.norm(W12_true @ (S["gTE"][pa, C.T2_SETS.index((t, r, 1000.0)), s, :p] * dte))
                 for t in (70.0, 82.0, 95.0) for r in (1.3, 1.6, 1.8) for s in wm)
        if mk >= 1:
            groups.append(g[:, None])
    groups.append(pick(np.moveaxis(S["gSNR"][pa, :, 0, 1], 0, 1)[..., :p]).T)          # N-SNR 언제나
    if pa != pb:
        gP = pick(S["gP"][CENTRE][:, :p])
        if max(np.linalg.norm(W12_true @ S["gP"][CENTRE, s, :p]) for s in wm) >= 1:
            groups.append(gP[:, None])
        mxp = S["gMxP_A" if dte != 0 else "gMxP_B"][..., :p]
        if max(np.linalg.norm(W12_true @ mxp[s, j]) for s in wm for j in range(3)) >= 1:
            groups.append(pick(mxp).T)
    return groups


_STV = None


def _init_e4():
    global _STV
    S = dict(np.load(C.store_path()))
    stv = {}
    for k, v in S.items():
        a, b = k.split("::")
        stv.setdefault(a, {})[b] = v
    sg = json.load(open(OUT / "directions.json"))["results"]
    for k in stv:
        stv[k]["sg2_pass"] = sg[k]["C3"]["SG2"]["pass_"]
    _STV = stv


def e4_config(args):
    F, gD, gU, pair, n, st, lam, seed = args
    p = lib.STAGE_P[st]
    E2 = np.load(C.WORK / f"e2a_{F}_{gD}.npy", mmap_mode="r")
    stv = _STV
    rng = np.random.default_rng(seed)
    A, Bs = (_state_idx(x) for x in PAIRS[pair])
    t2e = range(4)

    def draw(N):
        s = rng.choice(C.WM, N); f = rng.integers(0, 8, N); t = rng.integers(0, 4, N)
        ra = rng.integers(0, C.R_REP, N); rb = (ra + 1 + rng.integers(0, C.R_REP - 1, N)) % C.R_REP
        ya = np.asarray(E2[A[0], A[1], t, A[2], s, f, ra, :p], float)
        yb = np.asarray(E2[Bs[0], Bs[1], t, Bs[2], s, f, rb, :p], float)
        return yb - ya

    big = draw(5000)
    Wt = np.cov(big.T); W12t = lib.inv_sqrt(Wt); Dt = big.mean(0)
    groups = e4_groups(stv, F, gU, pair, p, W12t, Dt)
    Pt, kt = lib.proj_from_basis(lib.group_basis(groups, W12t))
    spike = np.zeros(p)
    if lam > 0:
        t = stv[f"{F}|{gD}"]["tj"][Bs[0]][C.WM, 0, :p].mean(0)
        c = math.sqrt(lam / n) / max(np.linalg.norm(Pt @ (W12t @ t)), 1e-12)
        spike = c * t
    psi_true = float(np.linalg.norm(Pt @ (W12t @ (Dt + spike))) ** 2)
    den_true = float(np.linalg.norm(W12t @ (Dt + spike)) ** 2)
    rec = []
    for m in range(M_SETS):
        Dm = draw(n) + spike
        r = lib.index_stats(Dm, None, rng, groups=groups)
        if r.get("exhausted"):
            rec.append(dict(exhausted=True)); continue
        rec.append(dict(rej=r["p_value"] < 0.05, cov=r["Psi_ci"][0] <= psi_true <= r["Psi_ci"][1] if r["Psi_ci"][1] > 0 else psi_true <= 1e-12 or r["Psi_ci"][1] >= psi_true,
                        psi=r["Psi"], I=r["I"], rep=r["reported"]))
    ok = [x for x in rec if not x.get("exhausted")]
    exh = 1 - len(ok) / len(rec)
    if not ok:
        return dict(cfg=args[:7], exhausted=True, exhausted_frac=1.0, k_true=int(kt))
    Is = [x["I"] for x in ok if x["I"] is not None]
    return dict(cfg=args[:7], k_true=int(kt), psi_true=psi_true, I_true=psi_true / den_true if den_true > 0 else None,
                rej=float(np.mean([x["rej"] for x in ok])), cov=float(np.mean([x["cov"] for x in ok])),
                psi_bias=float(np.mean([x["psi"] for x in ok]) - psi_true),
                I_bias=(float(np.mean(Is) - psi_true / den_true) if Is and den_true > 0 else None),
                reported=float(np.mean([x["rep"] for x in ok])), n_sets=len(ok), exhausted_frac=exh)


def run_e4(nproc=16):
    for F in FITTERS:
        for gD in C.GENERATORS:
            path = C.WORK / f"e2a_{F}_{gD}.npy"
            if not path.exists():
                np.save(path, e2a_array(F, gD))
    cfgs, seed = [], 0
    for F in FITTERS:
        for gD in C.GENERATORS:
            for gU in C.GENERATORS:
                for pair in PAIRS:
                    for n in (10, 20):
                        for st in STAGES:
                            for lam in LAMBDAS:
                                seed += 1
                                cfgs.append((F, gD, gU, pair, n, st, lam, 20261020 + seed))
    with ProcessPoolExecutor(nproc, initializer=_init_e4) as ex:
        res = list(ex.map(e4_config, cfgs, chunksize=1))
    json.dump(res, open(C.guard_write(OUT / "e4.json"), "w"), ensure_ascii=False, indent=1, default=float)
    print("E4 끝", len(res))


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "directions":
        run_directions()
    elif cmd == "e4":
        run_e4(int(sys.argv[2]) if len(sys.argv) > 2 else 16)
    else:
        raise SystemExit(cmd)
