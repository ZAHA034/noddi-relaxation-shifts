"""decomp 3 단계 분석 — SPEC_v3 §3–§9 (X7–X13). 메인 환경 · CPU 만 · F2 만.

사용: python analyze3.py x13 | x10 | x7 | x8 | x9 | x11 | x12 | all
결과: decomp/results/stage3/<이름>.json · cards/
구현 선택(결과 전 고정 — SPEC_v3 에 없는 세부): IMPL3.
"""
from __future__ import annotations

import itertools
import json
import math
import os
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze as A  # noqa: E402
import analyze2 as A2  # noqa: E402
import config as C  # noqa: E402
import lib  # noqa: E402

OUT3 = C.out_dir(3)                    # 정본 results/stage3 · F1 판 results/stage3_F1(중계 ZZZ-3 나)
FIT = A2.FIT
FEATS = A2.FEATS                                    # frac · S · C3 · C7
ETAN = [n for n, _, _ in C.ETA_E6]
IMPL3 = [
    "E6 특징은 섬유 방향 4 개 평균(무잡음). E6n 관측은 그 경우의 섬유 방향 하나(기준 색인 % 4 의 짝수 색인) · Rician SNR 30 · F2 적합.",
    "frac 특징의 관측은 겉보기 분율 (x, y) + N(0, C3 Σ 의 2×2 블록) — 결정적 시드 20261021.",
    "X8 의 분석가 가정 η 집합 = L1(CV 3 % · ρ 0 · iii) 16 점 + 중심 — 모든 경우에 같다(건강 성인 가정). 같은 사람(L0) 이라는 정보는 쓰지 않는다.",
    "X8 구간: S_A · S_B 를 따로 구해 Δ = [min S_B − max S_A, max S_B − min S_A]. 빈 집합이면 «구간 없음» 으로 세고 포함으로 치지 않는다.",
    "X9 의 «K 안» 판정: 적합 좌표 (x, y) 이동 방향(선 · mod 180°)이 분율 수준 원뿔 각 범위(η = L1 전체 위의 대수 이동 방향) 안인가.",
    "X12 의 K 생성 공간: 쌍의 특징 이동 벡터 Δ_q(θ, η)(기준 조직 · L1 η) 를 Σ_A 로 백색화해 쌓은 행렬의 앞 두 왼쪽 특이벡터.",
    "카드의 원뿔 각은 L4(사전 격자 54) 와 L1 두 판을 낸다. 외부 자료 카드는 대수 칸만(적합기 수준 «미평가»).",
]


def lsets():
    """L 수준 이름 → η 색인 목록(SPEC_v3 §4). 기준 η 는 각 집합의 시나리오 중심."""
    s = {}
    for sc in ("i", "ii", "iii"):
        s[f"L0|{sc}"] = ([ETAN.index(f"L0|{sc}")], ETAN.index(f"L0|{sc}"))
    for cv in (0.015, 0.03, 0.04, 0.06):
        for rho in (0.0, 0.8):
            key = f"L1|cv{cv}|rho{rho}|iii"
            s[key] = ([ETAN.index("L0|iii")] + [i for i, n in enumerate(ETAN) if n.startswith(key + "|")], ETAN.index("L0|iii"))
    for sc in ("i", "ii"):
        key = f"L1|cv0.03|rho0.0|{sc}"
        s[key] = ([ETAN.index(f"L0|{sc}")] + [i for i, n in enumerate(ETAN) if n.startswith(key + "|")], ETAN.index(f"L0|{sc}"))
    s["L1s|cv0.03|rho0.0|iii"] = ([ETAN.index("L0|iii")] + [i for i, n in enumerate(ETAN) if n.startswith("L1s|")], ETAN.index("L0|iii"))
    for n in ETAN:
        if n.startswith("L2|") or n.startswith("L3|"):
            s[n] = ([ETAN.index(n)], ETAN.index("L0|iii"))
    return s


def load_table(exp="E6"):
    nT, nE = len(C.TISSUE_E5), len(C.ETA_E6)
    feat = {a: {g: np.zeros((nT, nE, 7)) for g in C.GENERATORS} for a in C.ACQ_E5}
    frac = {a: {g: np.zeros((nT, nE, 2)) for g in C.GENERATORS} for a in C.ACQ_E5}
    cnt = {a: {g: np.zeros((nT, nE)) for g in C.GENERATORS} for a in C.ACQ_E5}
    for prot in ("P2a", "P2b"):
        cols, y = A.load(exp, prot, FIT)
        nu, fa = cols["nu_app"], cols["f_app"]
        fx = np.c_[np.log(nu / (1 - nu)), np.log(fa / (1 - fa)) - np.log(nu)]
        for a in C.ACQ_E5:
            for g in C.GENERATORS:
                m = (cols["acq"] == a) & (cols["gen"] == g)
                if m.any():
                    np.add.at(feat[a][g], (cols["tissue"][m], cols["eta6"][m]), y[m].astype(float))
                    np.add.at(frac[a][g], (cols["tissue"][m], cols["eta6"][m]), fx[m])
                    np.add.at(cnt[a][g], (cols["tissue"][m], cols["eta6"][m]), 1)
    for a in C.ACQ_E5:
        for g in C.GENERATORS:
            feat[a][g] /= cnt[a][g][..., None]
            frac[a][g] /= cnt[a][g][..., None]
    return feat, frac


def sig_block(Sig, acq, state, fname):
    S = A2.acq_sigma(Sig, acq, state)
    if fname == "frac":
        return S[:2, :2]
    lo, hi = FEATS[fname]
    return S[lo:hi, lo:hi]


def tab(feat, frac, acq, g, fname):
    if fname == "frac":
        return frac[acq][g]
    lo, hi = FEATS[fname]
    return feat[acq][g][..., lo:hi]


STATE = None


def ref_states():
    return {r: A2.coarse_state(*C.TISSUE_E5[r]) for r in C.refs_e5()}


# ── X13 · 대수 재현 ──────────────────────────────────────────────────────
def x13():
    rng = np.random.default_rng(20261021)
    out = {}
    TE, nu = 90.0, 0.6
    rows = []
    for cv in (0.015, 0.03, 0.04):
        for rho in (0.0, 1.0):
            z = rng.standard_normal((1_000_000, 2))
            z2 = rho * z[:, 0] + math.sqrt(max(0.0, 1 - rho ** 2)) * z[:, 1]
            t2in, t2en = 82.0 * (1 + cv * z[:, 0]), 52.0 * (1 + cv * z2)
            ai, ae = np.exp(-TE / t2in), np.exp(-TE / t2en)
            nup = nu * ai / (nu * ai + (1 - nu) * ae)
            rows.append(dict(cv=cv, rho=rho, ndi_app_mean=float(nup.mean()), ndi_app_sd=float(nup.std())))
    out["가_NDI_SD"] = rows
    z = rng.standard_normal((1_000_000, 2))
    r = (82.0 * (1 + 0.03 * z[:, 0])) / (52.0 * (1 + 0.03 * z[:, 1]))
    out["나_r_pm2SD"] = [float(r.mean() - 2 * r.std()), float(r.mean() + 2 * r.std()), float(np.percentile(r, 2.275)), float(np.percentile(r, 97.725))]
    les = {"건강": (82, 52), "병변 평균": (73, 142), "가장 약한": (84, 110), "가장 강한": (63, 173)}
    ang = {}
    for k, (ti, te_) in les.items():
        dx, dy = 1 / te_ - 1 / ti, 1 / ti - 1 / 2000
        ang[k] = dict(r=ti / te_, angle_deg=math.degrees(math.atan2(dy, dx)), ndi_rises_with_TE=bool(dx > 0))
    out["ㄷ_병변"] = ang
    out["T2iso_근거"] = dict(dy_ratio_1000_vs_2000=(1 / 1000 - 1 / 2000) / (1 / 82 - 1 / 1000),
                            between_subject_frac=(1 / 1730 - 1 / 2160) / (1 / 82 - 1 / 2030))
    return out


# ── X10 · 쌍 분류 격자 ────────────────────────────────────────────────────
NOT2ING_FIXED = {"NOT2ING-89": (4.004, 4.031), "NOT2ING-98": (4.031, 4.378)}


def x10(st):
    out = {}
    sig = {}
    for g in C.GENERATORS:
        E2 = np.load(C.WORK / f"e2a_{FIT}_{g}.npy", mmap_mode="r")
        for si, snr in enumerate(C.SNR_N1):
            X = np.asarray(E2[:, si, 0], float)
            Xc = X - X.mean(-2, keepdims=True)
            sig[(g, snr)] = (np.einsum("...ri,...rj->...ij", Xc, Xc) / (X.shape[-2] - 1)).mean(-3)
    for version in ("v2_범위", "확정"):
        sc = dict(A2.SCANNERS)
        if version == "확정":
            for k, tr in NOT2ING_FIXED.items():
                sc[k] = (sc[k][0], sc[k][1], tr)
        names = list(sc)
        recs = []
        for i, j in itertools.combinations(range(len(names)), 2):
            a, b = sc[names[i]], sc[names[j]]
            dte = max(abs(a[1][0] - b[1][1]), abs(a[1][1] - b[1][0]))
            dtr = 1000 * max(abs(a[2][0] - b[2][1]), abs(a[2][1] - b[2][0]))
            rec = dict(pair=f"{names[i]}~{names[j]}", dTE_max=dte, dTR_max_ms=dtr)
            for stg in ("C3", "C7"):
                p = lib.STAGE_P[stg]
                for snr in C.SNR_N1:
                    mte = mtr = 0.0
                    for g in C.GENERATORS:
                        S = st[f"{FIT}|{g}"]
                        for pr in (1, 2):
                            for s in C.WM:
                                W12 = lib.inv_sqrt(2 * sig[(g, snr)][pr, 1, s, :p, :p])
                                for t in C.ETA_T2:
                                    mte = max(mte, np.linalg.norm(W12 @ (S["gTE"][pr, C.T2_SETS.index(t), s, :p] * dte)))
                                for k in range(6):
                                    mtr = max(mtr, np.linalg.norm(W12 @ (S["gTR"][pr, k, s, :p] * dtr)))
                    for thr in (0.5, 1.0, 2.0):
                        rec[f"{stg}|SNR{int(snr)}|{thr}σ"] = "일치" if (mte < thr and mtr < thr) else "불일치"
                    rec[f"{stg}|SNR{int(snr)}|m"] = [float(mte), float(mtr)]
            recs.append(rec)
        out[version] = recs
    # 장비 기준 21 쌍 요약: NOT2ING 은 두 프로토콜 중 하나라도 «일치» 면 그 격자 칸에서 «일치 가능»
    return out


# ── X7 · 관측 동치(L 축) ─────────────────────────────────────────────────
def x7(feat, frac, st):
    T = np.array(C.TISSUE_E5)
    rs = ref_states()
    out = {}
    for lname, (idx, ec) in lsets().items():
        for g in C.GENERATORS:
            Sig = st[f"{FIT}|{g}"]["Sig"]
            for fname in FEATS:
                lo, hi = (0, 2) if fname == "frac" else FEATS[fname]
                p = hi - lo
                thr = stats.chi2.ppf(0.95, p)
                TA = tab(feat, frac, "A", g, fname)
                for q in A2.PAIRS2:
                    bq = "A" if q == "Q0" else q
                    TB = tab(feat, frac, bq, g, fname)
                    rows = []
                    for r, s in rs.items():
                        Si = np.linalg.inv(sig_block(Sig, "A", s, fname) + sig_block(Sig, bq, s, fname))
                        fA = TA[r, ec]
                        def conf(Tb):
                            d = Tb[:, idx] - fA
                            return (np.einsum("tei,ij,tej->te", d, Si, d) <= thr).any(1)
                        cAB, cAA = conf(TB), conf(TA)
                        rec = {}
                        for j, ax in ((0, "NDI"), (2, "FISO")):
                            mab = np.abs(T[cAB, j] - T[r, j]).max() if cAB.any() else 0.0
                            maa = np.abs(T[cAA, j] - T[r, j]).max() if cAA.any() else 0.0
                            rec[ax] = bool(mab - maa >= 2 * A2.DSTEP[ax] - 1e-12)
                            edge = (C.NDI_DENSE[0], C.NDI_DENSE[-1]) if j == 0 else (C.FISO_DENSE[0], C.FISO_DENSE[-1])
                            rec[ax + "_base_edge"] = bool(cAA.any() and np.isin(T[cAA, j], edge).any())
                        rec["empty"] = not cAB.any()
                        rows.append(rec)
                    o = dict(n=len(rows))
                    for ax in ("NDI", "FISO"):
                        be = float(np.mean([x[ax + "_base_edge"] for x in rows]))
                        o[ax + "_base_edge"] = be
                        o[ax + "_saturated"] = be >= 0.5
                        o[ax + "_unident"] = float(np.mean([x[ax] for x in rows]))
                    o["empty"] = float(np.mean([x["empty"] for x in rows]))
                    out[f"{lname}|{g}|{fname}|{q}"] = o
        print("X7", lname, flush=True)
    return out


# ── X8 · 부분 식별 구간 포함률 ─────────────────────────────────────────────
def load_obs():
    """E6n → obs[(gen, case_i, side, acq)] = (반복 10, 7) 특징 · 경우 목록(C.e6n_cases 순서)."""
    cases = C.e6n_cases()
    cidx = {(c[0], c[4]): i for i, c in enumerate(cases)}
    acc = {}
    for prot in ("P2a", "P2b"):
        cols, y = A.load("E6n", prot, FIT)
        for k in range(len(y)):
            key = (str(cols["gen"][k]), cidx[(int(cols["ref"][k]), str(cols["case"][k]))], str(cols["side"][k]), str(cols["acq"][k]))
            acc.setdefault(key, []).append(y[k].astype(float))
    return {k: np.array(v) for k, v in acc.items()}, cases


def x8(feat, frac, st, obs, cases):
    T = np.array(C.TISSUE_E5)
    rng = np.random.default_rng(20261021)
    rs = ref_states()
    idx = lsets()["L1|cv0.03|rho0.0|iii"][0]
    out = {}
    for gt in C.GENERATORS:                       # 관측의 생성기
        for ga in C.GENERATORS:                   # 분석가가 가정한 생성기(표)
            Sig = st[f"{FIT}|{ga}"]["Sig"]
            for fname in ("frac", "C3", "C7"):
                lo, hi = (0, 2) if fname == "frac" else FEATS[fname]
                p = hi - lo
                thr = stats.chi2.ppf(0.95, p)
                for q in A2.PAIRS2:
                    bq = "A" if q == "Q0" else q
                    rec = {"NDI": [], "FISO": [], "w_NDI": [], "w_FISO": [], "edge_NDI": [], "edge_FISO": [], "empty": 0}
                    for ci, (ra, rb, ea, eb, case) in enumerate(cases):
                        s = rs[ra]
                        SA, SB = sig_block(Sig, "A", s, fname), sig_block(Sig, bq, s, fname)
                        if fname == "frac":
                            ia, ib = ETAN.index(ea), ETAN.index(eb)
                            oa = frac["A"][gt][ra, ia][None] + rng.multivariate_normal(np.zeros(2), SA, 10)
                            ob = frac[bq][gt][rb, ib][None] + rng.multivariate_normal(np.zeros(2), SB, 10)
                        else:
                            oa = obs[(gt, ci, "A", "A")][:, lo:hi]
                            ob = obs[(gt, ci, "B", bq)][:, lo:hi]
                        TA, TB = tab(feat, frac, "A", ga, fname)[:, idx], tab(feat, frac, bq, ga, fname)[:, idx]
                        SAi, SBi = np.linalg.inv(SA), np.linalg.inv(SB)
                        for k in range(10):
                            dA = TA - oa[k]; dB = TB - ob[k]
                            inA = (np.einsum("tei,ij,tej->te", dA, SAi, dA) <= thr).any(1)
                            inB = (np.einsum("tei,ij,tej->te", dB, SBi, dB) <= thr).any(1)
                            if not inA.any() or not inB.any():
                                rec["empty"] += 1
                                rec["NDI"].append(False); rec["FISO"].append(False)
                                continue
                            for j, ax in ((0, "NDI"), (2, "FISO")):
                                lo_, hi_ = T[inB, j].min() - T[inA, j].max(), T[inB, j].max() - T[inA, j].min()
                                tru = T[rb, j] - T[ra, j]
                                rec[ax].append(lo_ - 1e-12 <= tru <= hi_ + 1e-12)
                                rec["w_" + ax].append(hi_ - lo_)
                                grid = C.NDI_DENSE if j == 0 else C.FISO_DENSE
                                rec["edge_" + ax].append(bool(np.isin(T[inA, j], (grid[0], grid[-1])).any() or np.isin(T[inB, j], (grid[0], grid[-1])).any()))
                    n = len(rec["NDI"])
                    out[f"{gt}|{ga}|{fname}|{q}"] = dict(n=n, cov_NDI=float(np.mean(rec["NDI"])), cov_FISO=float(np.mean(rec["FISO"])),
                                                         width_NDI_median=float(np.median(rec["w_NDI"])) if rec["w_NDI"] else None,
                                                         width_FISO_median=float(np.median(rec["w_FISO"])) if rec["w_FISO"] else None,
                                                         edge_NDI=float(np.mean(rec["edge_NDI"])) if rec["edge_NDI"] else None,
                                                         edge_FISO=float(np.mean(rec["edge_FISO"])) if rec["edge_FISO"] else None,
                                                         empty_frac=rec["empty"] / max(n, 1))
            print("X8", gt, ga, flush=True)
    matched = [v for k, v in out.items() if k.split("|")[0] == k.split("|")[1] and k.split("|")[2] in ("C3", "C7")]   # frac 은 검산(SPEC_v3 §5)
    hpi = all(v["cov_NDI"] >= 0.9 and v["cov_FISO"] >= 0.9 for v in matched)
    return dict(cells=out, HPI_pass=bool(hpi),
                frac_cells_ok=float(np.mean([v["cov_NDI"] >= 0.9 and v["cov_FISO"] >= 0.9 for v in matched])),
                width_ok=dict(NDI=float(np.mean([(v["width_NDI_median"] or 9) < 0.25 for v in matched])),
                              FISO=float(np.mean([(v["width_FISO_median"] or 9) < 0.145 for v in matched]))))


# ── X9 · 카드 «가능» 잔존 ─────────────────────────────────────────────────
def cone_angles(qacq, eta_idx):
    """분율 수준 이동 방향(선 · mod 180°)의 범위 — A 대 q, η 목록 위."""
    _, teA, trA = C.ACQ_E5["A"]
    _, te, tr = C.ACQ_E5[qacq]
    angs = []
    for e in eta_idx:
        _, t2, t1 = C.ETA_E6[e]
        aA, aB = C.weights(teA, t2, trA, t1), C.weights(te, t2, tr, t1)
        dx, dy = math.log((aB[0] / aB[1]) / (aA[0] / aA[1])), math.log((aB[2] / aB[0]) / (aA[2] / aA[0]))
        angs.append(math.degrees(math.atan2(dy, dx)) % 180)
    return min(angs), max(angs)


def x9(feat, st):
    rs = ref_states()
    out = {}
    idx = lsets()["L1|cv0.03|rho0.0|iii"][0]
    ec = ETAN.index("L0|iii")
    for g in C.GENERATORS:
        Sig = st[f"{FIT}|{g}"]["Sig"]
        for q in ("Q1", "Q2", "Q3", "Q4"):
            lo_a, hi_a = cone_angles(q, idx)
            odi_ok, inK = [], []
            for r, s in rs.items():
                S = A2.acq_sigma(Sig, "A", s)[:3, :3] + A2.acq_sigma(Sig, q, s)[:3, :3]
                d = feat[q][g][r, ec, :3] - feat["A"][g][r, ec, :3]
                odi_ok.append(abs(d[2]) <= 2 * math.sqrt(S[2, 2]))
                a = math.degrees(math.atan2(d[1], d[0])) % 180
                inK.append(lo_a - 1 <= a <= hi_a + 1 if (lo_a <= hi_a) else True)
            out[f"{g}|{q}"] = dict(ODI_within_2sd=float(np.mean(odi_ok)), xy_dir_in_fraction_cone=float(np.mean(inK)), cone_deg=[lo_a, hi_a])
    return out


# ── X11 · 카드 ────────────────────────────────────────────────────────────
EXT = {
    "Tong2020": {"centres": [f"C{i:02d}" for i in range(1, 11)], "acq": {f"C{i:02d}": ("Prisma-S", 71.0, 5400.0, "same") for i in range(1, 11)},
                 "src": "Tong et al. Sci Data 2020 (PMC7253426): 10 × Prisma 64ch · TR/TE 5.4 s/71 ms · b 1000/2000/3000 × 30"},
    "MUSHAC": {"acq": {"GE-ST": ("GE", 89.0, None, "41.4/26.0"), "Prisma-ST": ("Siemens", 89.0, 7200.0, "38.3/19.5"),
                       "Prisma-SA": ("Siemens", 80.0, 4500.0, "41.8/28.5"), "Connectom-ST": ("Siemens", 89.0, 7200.0, "31.1/8.5"),
                       "Connectom-SA": ("Siemens", 68.0, 5400.0, "불명")},
               "src": "Tax et al. NeuroImage 2019 Table 2 (PMC6556555) — GE TR 은 심박 동기(값 없음)"},
    "Wongkornchaovalit2025": {"acq": {f"TE{te}": ("Siemens?", float(te), 5800.0, "15.2/25.2") for te in (62, 72, 82, 92, 102, 112, 122, 132)},
                              "src": "Wongkornchaovalit et al. Sci Data 2025 (PMC12246157) — TR 5800 · δ/Δ 15.2/25.2 · 반복 62R2 · 132R2"},
}


def card_algebra(teA, trA, teB, trB):
    """분율 수준 카드 칸: L4(사전 격자 54) · L1(중심 ± 2 SD) 위의 이동 방향 범위 · 부호 · 원뿔 폭."""
    res = {}
    grids = {"L4": [(C.t2_triplet(t2), C.t1_triplet(t1)) for t2 in C.ETA_T2 for t1 in C.T1_SETS],
             "L1": [(t2, t1) for _, t2, t1 in C.ETA_E6 if _.startswith("L1|cv0.03|rho0.0|")]}
    for gname, grid in grids.items():
        vecs = []
        for t2, t1 in grid:
            aA = C.weights(teA, t2, trA, t1) if trA else C.weights(teA, t2)
            aB = C.weights(teB, t2, trB, t1) if trB else C.weights(teB, t2)
            vecs.append((math.log((aB[0] / aB[1]) / (aA[0] / aA[1])), math.log((aB[2] / aB[0]) / (aA[2] / aA[0]))))
        v = np.array(vecs)
        ang = np.degrees(np.arctan2(v[:, 1], v[:, 0]))
        res[gname] = dict(angle_deg=[float(ang.min()), float(ang.max())], cone_width_deg=float(ang.max() - ang.min()),
                          dy_sign=("+" if (v[:, 1] > 0).all() else "−" if (v[:, 1] < 0).all() else "±"),
                          dx_sign=("+" if (v[:, 0] > 0).all() else "−" if (v[:, 0] < 0).all() else "±"),
                          max_shift=float(np.abs(v).max()))
    return res


def x11(x10res, x7res, x9res):
    cards = []
    sc = dict(A2.SCANNERS)
    for k, tr in NOT2ING_FIXED.items():
        sc[k] = (sc[k][0], sc[k][1], tr)
    grid = {r["pair"]: r for r in x10res["확정"]}
    for a, b in itertools.combinations(list(sc), 2):
        A_, B_ = sc[a], sc[b]
        teA, teB = sum(A_[1]) / 2, sum(B_[1]) / 2
        trA, trB = 1000 * sum(A_[2]) / 2, 1000 * sum(B_[2]) / 2
        g = grid[f"{a}~{b}"]
        cards.append(dict(set="ON-Harmony", pair=f"{a}~{b}", dTE=teB - teA, dTR_ms=trB - trA, design_differs=A_[0] != B_[0],
                          algebra=card_algebra(teA, trA, teB, trB),
                          classification_grid={k: v for k, v in g.items() if "σ" in k},
                          n_match_cells=sum(v == "일치" for k, v in g.items() if "σ" in k),
                          unmodelled=["확산 시간 δ/Δ: 불명(사이드카에 키 없음 — LOCK_B EN :7074)", "기울기 비선형 · 재구성 · 소프트웨어(모형 밖)"],
                          trust={"대수": "A1–A4 · 사전 범위 · L 수준 조건부", "생성기": "NDI 는 T-in 에서만", "적합기": "F2 시뮬레이션(X7 · X9 — P2 프로토콜) · F1 대기" if FIT == "F2" else f"{FIT} 시뮬레이션(X7 · X9 — P2 프로토콜) · SPEC_v3 §0 예외"}))
    for setname, d in EXT.items():
        acqs = d["acq"]
        for a, b in itertools.combinations(list(acqs), 2):
            (_, teA, trA, ddA), (_, teB, trB, ddB) = acqs[a], acqs[b]
            cards.append(dict(set=setname, pair=f"{a}~{b}", dTE=teB - teA, dTR_ms=(trB - trA) if (trA and trB) else None,
                              algebra=card_algebra(teA, trA, teB, trB) if (trA and trB) else card_algebra(teA, None, teB, None),
                              unmodelled=(["확산 시간 δ/Δ 다름" if ddA != ddB else "확산 시간 같음"] +
                                          (["TR 심박 동기 — TR 차이 미상(TR 칸은 TR 항 없이 계산)"] if not (trA and trB) else []) +
                                          ["기울기 비선형 · 재구성(모형 밖)"]),
                              trust={"대수": "A1–A4 · 사전 범위 조건부", "적합기": "미평가(이 프로토콜을 합성하지 않음)"}, source=d["src"]))
    return cards


# ── X12 · ComBat 최소형 (SPEC_v3 §7) ───────────────────────────────────────
def x12(feat, st):
    import pandas as pd
    from neuroCombat import neuroCombat
    rng = np.random.default_rng(20261022)
    rs = list(ref_states().items())
    idx = lsets()["L1|cv0.03|rho0.0|iii"][0]
    nF = len(C.FISO_DENSE)
    out = {}

    def shift(r, kind):
        iv, io, jf = r // (3 * nF), (r // nF) % 3, r % nF
        if kind == "F":
            jf = min(jf + 2, nF - 1)
        else:
            iv = max(iv - 2, 0)
        return (iv * 3 + io) * nF + jf

    def label(m):
        return "반전" if m < 0 else "삭제" if m <= 0.1 else "부분" if m < 0.9 else "보존" if m <= 1.1 else "과장"

    for g in C.GENERATORS:
        Sig = st[f"{FIT}|{g}"]["Sig"]
        for q in ("Q1", "Q2", "Q3"):
            for kind in ("F", "N"):
                for design, covs in (("balanced", (True, False)), ("confounded", (False,))):
                    for use_cov in covs:
                        rat, rat_raw, errs = [], [], 0
                        for rep in range(200):
                            n = 20
                            pick = rng.choice(len(rs), 2 * n)
                            etas = rng.choice(idx, 2 * n)
                            site = np.r_[np.zeros(n), np.ones(n)].astype(int)
                            grp = (np.tile(np.r_[np.zeros(10), np.ones(10)], 2) if design == "balanced" else site.copy()).astype(int)
                            X, Xt = [], []
                            for k in range(2 * n):
                                r, s = rs[pick[k]]
                                th = shift(r, kind) if grp[k] else r
                                acq = "A" if site[k] == 0 else q
                                X.append(feat[acq][g][th, etas[k]] + rng.multivariate_normal(np.zeros(7), A2.acq_sigma(Sig, acq, s)))
                                Xt.append(feat["A"][g][th, etas[k]])
                            X, Xt = np.array(X), np.array(Xt)
                            W12 = lib.inv_sqrt(A2.acq_sigma(Sig, "A", rs[0][1]))
                            gw = W12 @ (Xt[grp == 1].mean(0) - Xt[grp == 0].mean(0))
                            covars = pd.DataFrame({"batch": site, "grp": grp})
                            try:
                                H = neuroCombat(dat=X.T, covars=covars, batch_col="batch",
                                                categorical_cols=(["grp"] if use_cov else None))["data"].T
                            except Exception:                                         # noqa: BLE001
                                errs += 1
                                continue
                            for Y, lst in ((H, rat), (X, rat_raw)):
                                gh = W12 @ (Y[grp == 1].mean(0) - Y[grp == 0].mean(0))
                                lst.append(float(gh @ gw / (gw @ gw)))
                        m, mr = (float(np.median(rat)) if rat else None), (float(np.median(rat_raw)) if rat_raw else None)
                        out[f"{g}|{q}|{kind}|{design}|{'cov' if use_cov else 'nocov'}"] = dict(
                            ratio_median=m, label=label(m) if m is not None else None, ratio_raw_median=mr,
                            label_raw=label(mr) if mr is not None else None, n=len(rat), errors=errs)
            print("X12", g, q, flush=True)
    return out


def main():
    OUT3.mkdir(parents=True, exist_ok=True)
    cmd = sys.argv[1]
    st = A2.store1()
    res = {"impl": IMPL3 if FIT == "F2" else IMPL3 + [f"적합기 {FIT} 판(SPEC_v3 §0 예외 · 중계 ZZZ-5) — 같은 코드 · 같은 기준 · 적합기만 바꿨다"]}
    if cmd in ("x13", "all"):
        res["X13"] = x13()
    if cmd in ("x10", "x11", "all"):
        res["X10"] = x10(st)
    if cmd in ("x7", "x8", "x9", "x11", "x12", "all"):
        feat, frac = load_table("E6")
    if cmd in ("x7", "x11", "all"):
        res["X7"] = x7(feat, frac, st) if cmd != "x11" else json.load(open(OUT3 / "x7.json"))["X7"]
    if cmd in ("x8", "all"):
        obs, cases = load_obs()
        res["X8"] = x8(feat, frac, st, obs, cases)
    if cmd in ("x9", "x11", "all"):
        res["X9"] = x9(feat, st)
    if cmd in ("x11", "all"):
        cards = x11(res["X10"], res.get("X7"), res["X9"])
        (OUT3 / "cards").mkdir(exist_ok=True)
        json.dump(cards, open(C.guard_write(OUT3 / "cards" / "cards.json"), "w"), ensure_ascii=False, indent=1, default=float)
        res["X11"] = dict(n_cards=len(cards))
    if cmd in ("x12", "all"):
        res["X12"] = x12(feat, st)
    json.dump(res, open(C.guard_write(OUT3 / f"{cmd}.json"), "w"), ensure_ascii=False, indent=1, default=float)
    print("씀", OUT3 / f"{cmd}.json")


if __name__ == "__main__":
    main()
