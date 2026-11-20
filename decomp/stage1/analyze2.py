"""decomp 2 단계 분석 — SPEC_v2 §11 의 X1–X6 · §13 형식. 메인 환경 · CPU 만 · F2 만.

사용:  python analyze2.py x1 | x2 | x3 [nproc] | x4 | x6 | all
결과: decomp/results/stage2/*.json (REPORT.md 는 summary 에서 옮긴다)
구현 선택(결과 전 고정 — SPEC_v2 에 없는 세부): IMPL2.
"""
from __future__ import annotations

import glob
import json
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze as A  # noqa: E402
import config as C  # noqa: E402
import lib  # noqa: E402

OUT2 = C.out_dir(2)                    # 정본 results/stage2 · F1 판 results/stage2_F1(중계 ZZZ-3 나)
#: 갈래의 적합기(중계 ZZZ-3 나) — 정본 F2 · F1 판은 DECOMP_TAG=_F1 DECOMP_FITTERS=F1. 정본 출력은 바이트 그대로가 되게 문자열을 갈래 F2 일 때 바꾸지 않는다
FIT = A.FITTERS[0]
FEATS = {"frac": (0, 2), "S": (3, 7), "C3": (0, 3), "C7": (0, 7)}   # 열 범위(C7 배열 기준). frac 은 겉보기 분율의 (x, y)
PAIRS2 = ("Q0", "Q1", "Q2", "Q3", "Q4")
IMPL2 = [
    "E5 특징은 섬유 방향 4 개의 평균(무잡음) — 조직 · 이완 · 획득마다 하나.",
    "Σ(SPEC_v2 §4): 기준 조직에 가장 가까운 1 단계 거친 격자 상태(NDI · ODI · FISO 유클리드) · 중심 T2 · TE_ref 는 획득 TE 에 가까운 쪽(92 → 100 · 68 · 70.6 → 70) · 획득의 프로토콜 · 같은 생성기. frac 특징은 C3 의 앞 2×2 블록, S 는 C7 의 뒤 4×4 블록.",
    "H2(같은 η 판): η 54 마다 기준 θ 의 특징을 그 η 에서 내고 같은 η 의 후보와 비교 — 방향은 고유 분율 (x, y) 의 선(mod 180°), 해상도 = θ 에서 본 θ′ 방향이 θ′ 를 NDI · FISO 로 한 칸 옮길 때 바뀌는 각의 최댓값.",
    "Q0(A 대 A)는 정의상 초과 0 — 검산으로만 보고.",
    "X6 의 g_TE · g_TR 은 두 프로토콜(P2a · P2b) 중 큰 쪽, 장비 TE · TR 이 범위면 범위 끝끼리의 최대 차.",
]


# ── E5 적재 ──────────────────────────────────────────────────────────────
def load_e5():
    """→ feat[acq][gen]: (tissue 2028, eta 54, 7) C7 좌표(F2) · frac[acq][gen]: (tissue, eta, 2) 겉보기 분율 (x, y)."""
    nT, nE = len(C.TISSUE_E5), len(C.ETA_T2) * 6
    feat = {a: {g: np.zeros((nT, nE, 7)) for g in C.GENERATORS} for a in C.ACQ_E5}
    frac = {a: {g: np.zeros((nT, nE, 2)) for g in C.GENERATORS} for a in C.ACQ_E5}
    cnt = {a: {g: np.zeros((nT, nE)) for g in C.GENERATORS} for a in C.ACQ_E5}
    for prot in ("P2a", "P2b"):
        cols, y = A.load("E5", prot, FIT)
        t2i = np.array([[C.T2_SETS.index(t) for t in C.ETA_T2].index(v) for v in cols["t2set"]])
        eta = t2i * 6 + cols["t1set"]
        nu, fa = cols["nu_app"], cols["f_app"]
        fx = np.c_[np.log(nu / (1 - nu)), np.log(fa / (1 - fa)) - np.log(nu)]
        for a in C.ACQ_E5:
            for g in C.GENERATORS:
                m = (cols["acq"] == a) & (cols["gen"] == g)
                if not m.any():
                    continue
                np.add.at(feat[a][g], (cols["tissue"][m], eta[m]), y[m].astype(float))
                np.add.at(frac[a][g], (cols["tissue"][m], eta[m]), fx[m])
                np.add.at(cnt[a][g], (cols["tissue"][m], eta[m]), 1)
    for a in C.ACQ_E5:
        for g in C.GENERATORS:
            feat[a][g] /= cnt[a][g][..., None]
            frac[a][g] /= cnt[a][g][..., None]
    return feat, frac


def store1():
    raw = dict(np.load(C.store_path()))
    st = {}
    for k, v in raw.items():
        a, b = k.split("::")
        st.setdefault(a, {})[b] = v
    return st


def coarse_state(nu, odi, fiso):
    return int(np.argmin([(nu - v) ** 2 + (odi - o) ** 2 + (fiso - f) ** 2 for v, o, f in C.TISSUE]))


def acq_sigma(Sig, acq, state):
    prot, te, _ = C.ACQ_E5[acq]
    pr = C.PROTOCOLS.index(prot)
    k = int(np.argmin([abs(te - x) for x in C.TE_REF]))
    return Sig[pr, 0, k, state]                                       # (7, 7) — 중심 T2(E2a 집합 0)


def refs():
    out = []
    for i, (v, o, f) in enumerate(C.TISSUE_E5):
        if 0.45 <= v <= 0.75 and 0.02 <= f <= 0.10 and o in (0.08, 0.20):
            out.append(i)
    return out


ETA_C = [C.T2_SETS.index(t) for t in C.ETA_T2].index(C.T2_SETS.index((82.0, 1.6, 1000.0))) * 6 + C.T1_SETS.index((1100.0, 1.0))
DSTEP = {"NDI": C.NDI_DENSE[1] - C.NDI_DENSE[0], "FISO": C.FISO_DENSE[1] - C.FISO_DENSE[0]}


# ── X1 ───────────────────────────────────────────────────────────────────
def x1():
    res = {}
    T = np.array(C.TISSUE_E5)
    eta = [(C.t2_triplet(t2), C.t1_triplet(t1)) for t2 in C.ETA_T2 for t1 in C.T1_SETS]
    errs, ang = [], {"Q1": [], "Q4": [], "Q2": []}
    dets = []
    A0 = C.ACQ_E5["A"]
    for t2, t1 in eta:
        aA = C.weights(A0[1], t2, A0[2], t1)
        xA, yA = C.frac_coords(*C.apparent(T[:, 0], T[:, 2], aA))
        vecs = {}
        for q in ("Q1", "Q2", "Q3", "Q4"):
            _, te, tr = C.ACQ_E5[q]
            aB = C.weights(te, t2, tr, t1)
            xB, yB = C.frac_coords(*C.apparent(T[:, 0], T[:, 2], aB))
            pred = np.array([math.log((aB[0] / aB[1]) / (aA[0] / aA[1])), math.log((aB[2] / aB[0]) / (aA[2] / aA[0]))])
            errs.append(max(np.abs(xB - xA - pred[0]).max(), np.abs(yB - yA - pred[1]).max()))
            vecs[q] = pred
            if q in ang:
                ang[q].append(math.degrees(math.atan2(pred[1], pred[0])))
        vte, vtr = vecs["Q1"], vecs["Q4"]
        dets.append(abs(vte[0] * vtr[1] - vte[1] * vtr[0]) / (np.linalg.norm(vte) * np.linalg.norm(vtr)))
    for q in ang:
        res[f"angle_deg_{q}"] = [min(ang[q]), max(ang[q])]
    res["max_translation_error"] = float(max(errs))
    res["sin_angle_TE_TR_min_max"] = [float(min(dets)), float(max(dets))]
    return res


# ── X2 · X5 ──────────────────────────────────────────────────────────────
def x2(feat, frac, st):
    out = {}
    T = np.array(C.TISSUE_E5)
    xt, yt = C.frac_coords(T[:, 0], T[:, 2])                         # 고유 분율의 (x, y) — H2 방향
    for g in C.GENERATORS:
        Sig = st[f"{FIT}|{g}"]["Sig"]
        for fname, (lo, hi) in FEATS.items():
            p = hi - lo
            thr = stats.chi2.ppf(0.95, p)
            for q in PAIRS2:
                bq = "A" if q == "Q0" else q
                rows = []
                for r in refs():
                    v, o, f = C.TISSUE_E5[r]
                    s = coarse_state(v, o, f)
                    if fname == "frac":
                        SA, SB = acq_sigma(Sig, "A", s)[:2, :2], acq_sigma(Sig, bq, s)[:2, :2]
                        fA, FB, FA_ = frac["A"][g][r, ETA_C], frac[bq][g].reshape(-1, 2), frac["A"][g].reshape(-1, 2)
                    else:
                        SA, SB = acq_sigma(Sig, "A", s)[lo:hi, lo:hi], acq_sigma(Sig, bq, s)[lo:hi, lo:hi]
                        fA = feat["A"][g][r, ETA_C, lo:hi]
                        FB, FA_ = feat[bq][g][..., lo:hi].reshape(-1, p), feat["A"][g][..., lo:hi].reshape(-1, p)
                    Si = np.linalg.inv(SA + SB)
                    def conf(F):
                        d = F - fA
                        return np.einsum("ni,ij,nj->n", d, Si, d) <= thr
                    cAB = conf(FB).reshape(len(C.TISSUE_E5), -1).any(1)
                    cAA = conf(FA_).reshape(len(C.TISSUE_E5), -1).any(1)
                    rec = {}
                    for j, ax in ((0, "NDI"), (2, "FISO")):
                        mab = np.abs(T[cAB, j] - T[r, j]).max() if cAB.any() else 0.0
                        maa = np.abs(T[cAA, j] - T[r, j]).max() if cAA.any() else 0.0
                        rec[ax] = bool(mab - maa >= 2 * DSTEP[ax] - 1e-12)
                        edge = (C.NDI_DENSE[0], C.NDI_DENSE[-1]) if j == 0 else (C.FISO_DENSE[0], C.FISO_DENSE[-1])
                        rec[ax + "_edge"] = bool(cAB.any() and np.isin(T[cAB, j], edge).any())
                    rec["empty"] = not cAB.any()
                    rows.append(rec)
                out[f"{g}|{fname}|{q}"] = dict(n_refs=len(rows), NDI_unident=float(np.mean([x["NDI"] for x in rows])),
                                                FISO_unident=float(np.mean([x["FISO"] for x in rows])),
                                                NDI_edge=float(np.mean([x["NDI_edge"] for x in rows])), FISO_edge=float(np.mean([x["FISO_edge"] for x in rows])),
                                                empty=float(np.mean([x["empty"] for x in rows])))
    return out, h2_same_eta(feat, frac, st)


def _ang_res(T, xt, yt, r, i):
    a = math.degrees(math.atan2(yt[i] - yt[r], xt[i] - xt[r])) % 180
    res = 0.0
    for dv, df in ((DSTEP["NDI"], 0), (0, DSTEP["FISO"])):
        v2 = min(max(T[i, 0] + dv, 1e-4), 1 - 1e-4); f2 = min(max(T[i, 2] + df, 1e-4), 1 - 1e-4)
        x2_, y2_ = C.frac_coords(v2, f2)
        a2 = math.degrees(math.atan2(y2_ - yt[r], x2_ - xt[r])) % 180
        res = max(res, min(abs(a2 - a), 180 - abs(a2 - a)))
    return a, res


def h2_same_eta(feat, frac, st):
    """SPEC_v2 §10 H2 — 같은 η 판: η 마다 C_AB^η · C_AA^η, 합집합 초과 θ′ 의 방향이 쐐기 [48.5°, 72.1°] ± 해상도 안인가."""
    T = np.array(C.TISSUE_E5)
    xt, yt = C.frac_coords(T[:, 0], T[:, 2])
    out = {}
    for g in C.GENERATORS:
        Sig = st[f"{FIT}|{g}"]["Sig"]
        for fname in ("frac", "C3", "C7"):
            lo, hi = FEATS[fname]
            p = hi - lo
            thr = stats.chi2.ppf(0.95, p)
            inside, nref_ex = [], 0
            for r in refs():
                v, o, f = C.TISSUE_E5[r]
                s = coarse_state(v, o, f)
                if fname == "frac":
                    SA, SB = acq_sigma(Sig, "A", s)[:2, :2], acq_sigma(Sig, "Q1", s)[:2, :2]
                    FA_all, FB_all = frac["A"][g], frac["Q1"][g]
                else:
                    SA, SB = acq_sigma(Sig, "A", s)[lo:hi, lo:hi], acq_sigma(Sig, "Q1", s)[lo:hi, lo:hi]
                    FA_all, FB_all = feat["A"][g][..., lo:hi], feat["Q1"][g][..., lo:hi]
                Si = np.linalg.inv(SA + SB)
                uab = np.zeros(len(T), bool); uaa = np.zeros(len(T), bool)
                for e in range(FA_all.shape[1]):
                    fA = FA_all[r, e]
                    dB = FB_all[:, e] - fA; dA = FA_all[:, e] - fA
                    uab |= np.einsum("ni,ij,nj->n", dB, Si, dB) <= thr
                    uaa |= np.einsum("ni,ij,nj->n", dA, Si, dA) <= thr
                ex = np.where(uab & ~uaa)[0]
                nref_ex += len(ex) > 0
                for i in ex:
                    if i == r:
                        continue
                    a, res = _ang_res(T, xt, yt, r, i)
                    inside.append((48.5 - res) <= a <= (72.1 + res))
            out[f"{g}|{fname}"] = dict(n_excess=len(inside), refs_with_excess=nref_ex, frac_inside=float(np.mean(inside)) if inside else None,
                                       H2_pass=(float(np.mean(inside)) >= 0.9) if inside else None)
    return out


# ── X3 ───────────────────────────────────────────────────────────────────
def _x3_cfg(args):
    F, g, pair, n, st_, lam, seed = args
    A.M_SETS = 1000
    r = _x3_one(F, g, pair, n, st_, lam, seed)
    return r


def _x3_one(F, gD, pair, n, st_, lam, seed):
    p = lib.STAGE_P[st_]
    E2 = np.load(C.WORK / f"e2a_{F}_{gD}.npy", mmap_mode="r")
    stv = A._STV
    rng = np.random.default_rng(seed)
    Aa, Bs = (A._state_idx(x) for x in A.PAIRS[pair])

    def draw(N):
        s = rng.choice(C.WM, N); f = rng.integers(0, 8, N); t = rng.integers(0, 4, N)
        ra = rng.integers(0, C.R_REP, N); rb = (ra + 1 + rng.integers(0, C.R_REP - 1, N)) % C.R_REP
        return np.asarray(E2[Bs[0], Bs[1], t, Bs[2], s, f, rb, :p], float) - np.asarray(E2[Aa[0], Aa[1], t, Aa[2], s, f, ra, :p], float)

    big = draw(5000)
    Wt = np.cov(big.T); W12t = lib.inv_sqrt(Wt); Dt = big.mean(0)
    groups = A.e4_groups(stv, F, gD, pair, p, W12t, Dt)
    Pt, _ = lib.proj_from_basis(lib.group_basis(groups, W12t))
    spike = np.zeros(p)
    if lam > 0:
        t = stv[f"{F}|{gD}"]["tj"][Bs[0]][C.WM, 0, :p].mean(0)
        spike = math.sqrt(lam / n) / max(np.linalg.norm(Pt @ (W12t @ t)), 1e-12) * t
    psi_true = float(np.linalg.norm(Pt @ (W12t @ (Dt + spike))) ** 2)
    cov_old, cov_new = [], []
    for m in range(1000):
        r = lib.index_stats(draw(n) + spike, None, rng, groups=groups, basic_ci=True)
        if r.get("exhausted"):
            continue
        lo, hi = r["Psi_ci"]; cov_old.append(lo <= psi_true <= hi if hi > 0 else psi_true <= 1e-12)
        lo2, hi2 = r["Psi_ci_basic"]; cov_new.append(lo2 <= psi_true <= hi2)
    return dict(cfg=[F, gD, pair, n, st_, lam], psi_true=psi_true, cov_v1=float(np.mean(cov_old)), cov_basic=float(np.mean(cov_new)), n_sets=len(cov_new))


def x3(nproc):
    cfgs, seed = [], 0
    for g in C.GENERATORS:
        for pair in ("N2", "N4"):
            for n in (10, 20):
                for st_ in A.STAGES:
                    for lam in A.LAMBDAS:
                        seed += 1
                        cfgs.append((FIT, g, pair, n, st_, lam, 20261021 + seed))
    with ProcessPoolExecutor(nproc, initializer=A._init_e4) as ex:
        res = list(ex.map(_x3_cfg, cfgs, chunksize=1))
    ok = [r["cov_basic"] >= 0.88 for r in res]
    return dict(cells=res, H_psi_pass=bool(all(ok)), frac_cells_ok=float(np.mean(ok)))


# ── X4 ───────────────────────────────────────────────────────────────────
def x4(feat, st):
    out = {}
    W2 = [[C.T2_SETS.index(t) for t in C.ETA_T2].index(C.T2_SETS.index((82.0, r, 1000.0))) * 6 + C.T1_SETS.index((1100.0, 1.0)) for r in (1.3, 1.6, 1.8)]
    widths = {"W1": [ETA_C], "W2": W2, "W3": list(range(54))}
    for g in C.GENERATORS:
        Sig = st[f"{FIT}|{g}"]["Sig"]
        for fname in ("C3", "C7"):
            lo, hi = FEATS[fname]
            for q in ("Q1", "Q2", "Q4"):
                for wn, ws in widths.items():
                    vals = []
                    for r in refs():
                        v, o, f = C.TISSUE_E5[r]
                        s = coarse_state(v, o, f)
                        S = acq_sigma(Sig, "A", s)[lo:hi, lo:hi] + acq_sigma(Sig, q, s)[lo:hi, lo:hi]
                        W12 = lib.inv_sqrt(S)
                        dhat = feat[q][g][r, :, lo:hi] - feat["A"][g][r, :, lo:hi]
                        vals += [float(np.linalg.norm(W12 @ (dhat[e] - dhat[ETA_C]))) for e in ws]
                    out[f"{g}|{fname}|{q}|{wn}"] = dict(median=float(np.median(vals)), p90=float(np.percentile(vals, 90)), max=float(np.max(vals)))
    return out


# ── X6 ───────────────────────────────────────────────────────────────────
SCANNERS = {  # LOCK_B EN :7748–7754 (TE ms · TR s) · 등급 :7055–7064 · b0 :7044–7049
    "OXF1PRI": ("S", (92.0, 92.0), (3.60, 3.60)), "OXF2PRI": ("S", (92.0, 92.0), (3.60, 3.60)),
    "OXF3TRI": ("S", (96.4, 96.4), (3.60, 3.60)), "NOT2ING-89": ("V", (89.0, 89.0), (4.00, 4.38)),
    "NOT2ING-98": ("V", (98.0, 98.0), (4.00, 4.38)), "NOT4GEP": ("V", (68.0, 68.0), (4.90, 4.90)),
    "OXF4GEP": ("V", (68.3, 71.0), (8.00, 8.00)), "NOT1ACH": ("V", (70.0, 71.2), (9.81, 10.62)),
}


def x6(st):
    out = []
    s50 = C.SNR_N1.index(50.0)
    sig50 = {}
    for g in C.GENERATORS:
        E2 = np.load(C.WORK / f"e2a_{FIT}_{g}.npy", mmap_mode="r")
        X = np.asarray(E2[:, s50, 0], float)                         # (prot, teref, tissue, fdir, rep, 7) — 중심 T2
        Xc = X - X.mean(-2, keepdims=True)
        sig50[g] = (np.einsum("...ri,...rj->...ij", Xc, Xc) / (X.shape[-2] - 1)).mean(-3)   # (prot, teref, tissue, 7, 7)
    names = list(SCANNERS)
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = SCANNERS[names[i]], SCANNERS[names[j]]
            dte = max(abs(a[1][0] - b[1][1]), abs(a[1][1] - b[1][0]))
            dtr = 1000 * max(abs(a[2][0] - b[2][1]), abs(a[2][1] - b[2][0]))
            rec = dict(pair=f"{names[i]}~{names[j]}", dTE_max=dte, dTR_max_ms=dtr, design_differs=a[0] != b[0])
            for stg in ("C3", "C7"):
                p = lib.STAGE_P[stg]
                mte = mtr = 0.0
                for g in C.GENERATORS:
                    S = st[f"{FIT}|{g}"]
                    for pr in (1, 2):
                        for s in C.WM:
                            W12 = lib.inv_sqrt(2 * sig50[g][pr, 1, s, :p, :p])
                            for t in C.ETA_T2:
                                mte = max(mte, np.linalg.norm(W12 @ (S["gTE"][pr, C.T2_SETS.index(t), s, :p] * dte)))
                            for k in range(6):
                                mtr = max(mtr, np.linalg.norm(W12 @ (S["gTR"][pr, k, s, :p] * dtr)))
                rec[f"m_TE_{stg}"] = float(mte); rec[f"m_TR_{stg}"] = float(mtr)
                rec[f"class_{stg}"] = "이완 불일치" if (mte >= 1 or mtr >= 1) else "이완 일치"
            out.append(rec)
    return out


def main():
    OUT2.mkdir(parents=True, exist_ok=True)
    cmd = sys.argv[1]
    st = store1()
    res = {}
    if cmd in ("x1", "all"):
        res["X1"] = x1()
    if cmd in ("x2", "x4", "all"):
        feat, frac = load_e5()
        if cmd in ("x2", "all"):
            res["X2"], res["X5_H2"] = x2(feat, frac, st)
        if cmd in ("x4", "all"):
            res["X4"] = x4(feat, st)
    if cmd in ("x3", "all"):
        res["X3"] = x3(int(sys.argv[2]) if len(sys.argv) > 2 else 16)
    if cmd in ("x6", "all"):
        res["X6"] = x6(st)
    res["impl"] = IMPL2 if FIT == "F2" else IMPL2 + [f"적합기 {FIT} 판(SPEC_v3 §0 예외 · 중계 ZZZ-5) — 같은 코드 · 같은 기준 · 적합기만 바꿨다"]
    path = C.guard_write(OUT2 / f"{cmd}.json")
    json.dump(res, open(path, "w"), ensure_ascii=False, indent=1, default=float)
    print("씀", path)


if __name__ == "__main__":
    main()
