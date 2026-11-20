"""S4 — 검증 · 카드 · 지도(SPEC_S4 §5 – §9 · §13). 메인 환경(sbi · torch) · CPU 만.

    uv run --frozen python decomp/s4/validate.py <단계> ...(main() 의 설명)
    결과 → decomp/results/s4/<단계>.json(저장소 · config.guard_write) · 큰 중간값 → WORK/s4/val/
⛔ ML 점추정값을 결과로 쓰지 않는다 — 결과는 수축 c · 사후 SD · MDE · 포함률 · 분류뿐이다(§0).
판정 · 분류는 생성기 가정 k(비틀림 셋 × T1 시나리오 셋 = 아홉)마다의 조건부 사후로 한다 — 아홉의 균등 혼합은 참고로만 낸다(중계 ZZZg-1 · §13 1).
"""
from __future__ import annotations

import json
import math
import os
import sys
import time
import warnings

import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s4cfg as S  # noqa: E402

torch.set_num_threads(int(os.environ.get("S4_THREADS", "16")))
warnings.filterwarnings("ignore", message=".*(outside the prior support|max_sampling_batch_size|triangular_solve).*")     # 학습 지지 밖 표본은 draw() 가 버린다

GEN_SCEN = [(g, s) for g in S.GENERATORS for s in S.SCENARIOS]
K_NAMES = [f"{g}|{s}" for g, s in GEN_SCEN]
CI = (0.50, 0.80, 0.90, 0.95)
CENTRE = dict(t2=(82.0, 82.0 / 1.6, S.T2ISO), t1=(1100.0, 1100.0, S.T1ISO))      # 이완 중심 · 시나리오 iii(§13 2)
VAL = S.pathlib.Path(os.environ["S4_VAL"]) if os.environ.get("S4_VAL") else S.WORK / "val"       # 처음 2 % 시간 재기는 별도 폴더(§10-1)
OUT = S.pathlib.Path(os.environ["S4_OUT"]) if os.environ.get("S4_OUT") else S.RESULTS     # 점검용 덮어쓰기(저장소 밖)
ARMS = ("c3", "c7")
N_ENS = int(os.environ.get("S4_ENS", "5"))
N_KEEP = dict(test=1000, test_k=250, cards=2000, maps=1000, inject=1000)                  # 경우 · k 마다 주 지지 표본 수(§13 14)
MCMC = dict(chains=4, warmup=500, keep=2500, rhat=1.01, ess=500, max_excluded=0.20)      # 참조 사후(§7 ii · ZZZg-2)


# ── 공통 ───────────────────────────────────────────────────────────────────────
def load_posts(model, n):
    return [torch.load(S.WORK / "models" / f"{model}_{i}.obj.pt", weights_only=False) for i in range(n)]


def with_gs(cond, g, s):
    c = cond.copy()
    for gg in S.GENERATORS:
        c[:, S.COND_NAMES.index(f"gen_{gg}")] = float(gg == g)
    for ss in S.SCENARIOS:
        c[:, S.COND_NAMES.index(f"scen_{ss}")] = float(ss == s)
    return c


def cov_prior(rng, scale=1.0, n=1_000_000):
    t = S.sample_theta(rng, n, S.D_MAIN * scale)
    return np.cov(t.T)


def prior_covs(rng, sens=True):
    c = dict(main=cov_prior(rng)[0:3, 0:3])
    if sens:
        c.update(half=cov_prior(rng, 0.5)[0:3, 0:3], two=cov_prior(rng, 2.0)[0:3, 0:3])
    return c


def apparent_shift(theta_A, acq, i0, relax=None):
    """u_relax(§13 2 · 11) — (n, 3). relax: (t2 (n,3), t1 (n,3)) — 없으면 이완 중심."""
    n = len(theta_A)
    t2 = np.tile(CENTRE["t2"], (n, 1)) if relax is None else relax[0]
    t1 = np.tile(CENTRE["t1"], (n, 1)) if relax is None else relax[1]
    out = []
    for s in ("A", "B"):
        a = np.exp(-acq[f"acq_TE_{s}"][i0][:, None] / t2) * (1 - np.exp(-acq[f"acq_TR_{s}"][i0][:, None] / t1))
        out.append(np.c_[S.C.apparent(theta_A[:, 0], theta_A[:, 1], (a[:, 0], a[:, 1], a[:, 2]))])
    d = np.c_[out[1] - out[0], np.zeros(n)]
    nrm = np.linalg.norm(d, axis=1)
    u = np.where(nrm[:, None] > 1e-6, d / np.maximum(nrm, 1e-12)[:, None], np.nan)
    return u, nrm


def u_ref_of(sim, plan):
    """P3 의 기준 방향(§13 9) — 시험 자료 전체의 «TR 만 · L0» 경우의 u_relax 평균(대수 — 사후와 무관)."""
    m = np.flatnonzero(sim["layer"] == "tr|L0")
    if len(m) == 0:
        return np.full(3, np.nan)
    u, _ = apparent_shift(sim["theta"][m, 3:6], plan, m)
    v = np.nanmean(u, 0)
    return v / np.linalg.norm(v)


def directions(theta, acq, i0, u_ref=None):
    """방향(n, 5, 3): e_NDI · e_FISO · e_ODI · u_relax · u_⊥ — u_ref 를 주면 여섯째로 붙인다(시험 자료 · P3)."""
    n = len(theta)
    E = np.tile(np.eye(3), (n, 1, 1))
    ur, _ = apparent_shift(theta[:, 3:6], acq, i0)
    up = np.c_[-ur[:, 1], ur[:, 0], np.zeros(n)]
    D = [E, ur[:, None], up[:, None]] + ([np.tile(u_ref, (n, 1, 1))] if u_ref is not None else [])
    return np.concatenate(D, 1)


U_NAMES = ("NDI", "FISO", "ODI", "u_relax", "u_perp")
U_NAMES_TEST = U_NAMES + ("u_ref",)


def draw(posts, X, cond, n_keep, gs=None, cap=20):
    """한 생성기 가정 k 의 조건부 사후(앙상블 모형들의 균등 혼합) 표본 — 주 사전 지지 안 n_keep 개(§13 1).
    gs=None 이면 cond 의 k 그대로(시험 경우의 참 k) · gs=(비틀림, 시나리오) 면 그 k 로 바꿔 넣는다.
    모형마다 같은 수를 뽑고 학습 지지 밖은 버린다 · 경우마다 뽑는 수는 최대 cap × n_keep(못 이르면 «재가중 불가»).
    반환: 표본 (b, n_keep, 6) · 모형 꼬리표 (b, n_keep) · 학습 지지 표본 (b, n_keep, 6) · 성공 (b,) · 주 지지 비율 (b,)."""
    b, M = len(X), len(posts)
    c = cond if gs is None else with_gs(cond, *gs)
    keep, keep2 = [[] for _ in range(b)], [[] for _ in range(b)]
    n_raw, n_tr, n_in = np.zeros(b), np.zeros(b), np.zeros(b)
    todo = np.arange(b)
    rng = np.random.default_rng(12345)
    while len(todo):
        acc = np.where(n_raw[todo] > 0, n_in[todo] / np.maximum(n_raw[todo], 1), 0.5)
        want = float(np.quantile((n_keep - n_in[todo]) / np.maximum(acc, 0.05), 0.75))
        k = max(1, math.ceil(min(want, cap * n_keep - n_raw[todo].min()) / M))
        xc = torch.as_tensor(np.c_[X[todo], c[todo]], dtype=torch.float32)
        for m, post in enumerate(posts):
            with torch.no_grad():
                smp = post.sample_batched((k,), x=xc, show_progress_bars=False, reject_outside_prior=False).numpy().transpose(1, 0, 2)
            for j, i in enumerate(todo):
                t = smp[j][S.in_main(smp[j], 2.0)]                   # 학습 지지(반폭 D_TRAIN = 2·D_MAIN · 조직 범위) 안만
                inm = S.in_main(t)
                keep2[i].append(np.c_[t, np.full(len(t), m)])
                keep[i].append(np.c_[t[inm], np.full(int(inm.sum()), m)])
                n_tr[i] += len(t)
                n_in[i] += inm.sum()
        n_raw[todo] += k * M
        todo = np.flatnonzero((n_in < n_keep) & (n_raw < cap * n_keep))
    T = np.full((b, n_keep, 6), np.nan, np.float32)
    T2 = np.full((b, n_keep, 6), np.nan, np.float32)
    tag = np.full((b, n_keep), -1, np.int16)
    ok = n_in >= n_keep
    for i in range(b):
        A, A2 = np.concatenate(keep[i]), np.concatenate(keep2[i])
        if ok[i]:
            s = A[rng.choice(len(A), n_keep, replace=False)]
            T[i], tag[i] = s[:, :6], s[:, 6]
        if len(A2):
            s2 = A2[rng.choice(len(A2), min(len(A2), n_keep), replace=False)]
            T2[i, :len(s2)] = s2[:, :6]
    return T, tag, T2, ok, n_in / np.maximum(n_tr, 1)


def cstat(D, U, pv, w=None):
    """c · sd — D (m, 3) Δ 표본 · U (n, 3) 방향(nan 이면 정의 안 됨) · pv (n,) 사전 분산 · w 표본 가중(합 1 · 참고 혼합)."""
    Uv = np.nan_to_num(U)
    P = D @ Uv.T
    var = P.var(0) if w is None else (w[:, None] * (P - w @ P) ** 2).sum(0)
    bad = ~np.isfinite(U).all(1)
    return np.where(bad, np.nan, 1 - var / pv), np.where(bad, np.nan, np.sqrt(var))


def stats_one(T, T2, ok, U, covs, sens=False):
    """한 관측 · 한 k — 방향마다 c · sd(§5) · sens 면 사전 민감도 ×0.5(추가 재가중 · 남은 표본 200 미만이면 nan) · ×2(재가중하지 않은 사후)(§7 iii)."""
    Uv = np.nan_to_num(U)
    pv = {k: np.einsum("ni,ij,nj->n", Uv, v, Uv) for k, v in covs.items()}
    nan = np.full(len(U), np.nan)
    r = dict(c=nan, sd=nan) | (dict(c_half=nan, n_half=0, c_two=nan) if sens else {})
    if not ok:
        return r
    r["c"], r["sd"] = cstat(T[:, 0:3], U, pv["main"])
    if sens:
        h = S.in_main(T, 0.5)
        r["n_half"] = int(h.sum())
        r["c_half"] = cstat(T[h, 0:3], U, pv["half"])[0] if h.sum() >= 200 else nan
        f2 = np.isfinite(T2).all(1)
        r["c_two"] = cstat(T2[f2, 0:3], U, pv["two"])[0] if f2.sum() >= 200 else nan
    return r


def mix_pool(Ts, accs):
    """참고용 균등 혼합(§13 1) — k 마다의 주 지지 표본 · 주 지지 비율 → 표본 · 가중(합 1).
    가중 ∝ 주 지지 비율 / 표본 수 = 학습 사전 사후 아홉의 균등 혼합을 주 사전으로 재가중한 것(옛 §13 1 의 정의)."""
    T = np.concatenate(Ts)
    w = np.concatenate([np.full(len(t), a / len(t)) for t, a in zip(Ts, accs)])
    return T, w / w.sum()


def wquantile(x, w, q):
    """가중 분위 — x (m,) · w (m,) 합 1 · q 목록."""
    o = np.argsort(x)
    cw = np.cumsum(w[o]) - 0.5 * w[o]
    return np.interp(q, cw, x[o])


def classify(c, hi=0.75, lo=0.25):
    return np.where(np.isnan(c), "재가중 불가", np.where(c >= hi, "식별", np.where(c < lo, "비식별", "부분")))


def jdump(name, obj):
    import config as C  # noqa: E402
    OUT.mkdir(parents=True, exist_ok=True)
    path = C.guard_write(OUT / f"{name}.json")
    json.dump(obj, open(path, "w"), ensure_ascii=False, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o))
    print("씀", path, flush=True)


def load_set(name):
    sim = dict(np.load(S.WORK / f"sim_{name}.npz", allow_pickle=False))
    plan = dict(np.load(S.WORK / f"plan_{name}.npz", allow_pickle=False))
    return sim, plan


def xobs(sim, arm):
    return sim["x_c3"] if arm == "c3" else sim["x_c7"]


def _nanmed(v):
    v = np.asarray(v, float)
    return float(np.nanmedian(v)) if np.isfinite(v).any() else None


# ── (i) · (iv) · P1–P6 — 시험 자료 ─────────────────────────────────────────────
def test_part(arm, part=0, parts=1, probe=None, batch=250):
    """시험 자료의 한 몫(색인 part::parts) — WORK/s4/val 에(§7 i · iv · §8 · §13 1 · 5 · 6 · 14).
    참 k 의 조건부 사후(표본 1,000): 포함 · 순위 · 모형별 포함 · TARP 표본 · c · sd
    k 아홉 각각의 조건부 사후(표본 250): c · sd(P1–P6) · 참고: 아홉의 균등 혼합의 c · 포함."""
    torch.manual_seed(1000 + part)                                     # 사후 표본의 torch 난수(몫마다 고정)
    sim, plan = load_set("test")
    n = len(sim["theta"]) if probe is None else probe
    covs = prior_covs(np.random.default_rng(1), sens=False)
    posts = load_posts(arm, N_ENS)
    idx = np.arange(n)[part::parts]
    X, cond, th = xobs(sim, arm)[idx], sim["cond"][idx], sim["theta"][idx]
    U = directions(th, plan, idx, u_ref_of(sim, plan))
    m_, nu, nk = len(idx), U.shape[1], len(GEN_SCEN)
    qs = [q for pp in CI for q in ((1 - pp) / 2, (1 + pp) / 2)]
    cover = np.zeros((m_, 3, len(CI)), bool)
    cover_m = np.zeros((N_ENS, m_, 3, len(CI)), bool)
    cover_mix = np.zeros((m_, 3, len(CI)), bool)
    ranks = np.full((m_, 3), np.nan)
    okk, ok_mix = np.zeros(m_, bool), np.zeros(m_, bool)
    samp = np.full((m_, N_KEEP["test"], 3), np.nan, np.float32)
    c_true, sd_true, c_mix = (np.full((m_, nu), np.nan) for _ in range(3))
    c_k, sd_k = np.full((m_, nk, nu), np.nan), np.full((m_, nk, nu), np.nan)
    acc_k = np.full((m_, nk), np.nan)
    inside = lambda q, i: (th[i, 0:3] >= q[0::2]) & (th[i, 0:3] <= q[1::2])         # q (2·len(CI), 3) → (len(CI), 3)
    t0 = time.time()
    for b0 in range(0, m_, batch):
        ii = np.arange(b0, min(m_, b0 + batch))
        T, tag, T2, ok, _ = draw(posts, X[ii], cond[ii], N_KEEP["test"])
        okk[ii] = ok
        for j, i in enumerate(ii):
            if not ok[j]:
                continue
            D = T[j, :, 0:3]
            samp[i] = D
            cover[i] = inside(np.quantile(D, qs, axis=0), i).T
            for m in range(N_ENS):
                Dm = D[tag[j] == m]
                if len(Dm) >= 20:
                    cover_m[m, i] = inside(np.quantile(Dm, qs, axis=0), i).T
            ranks[i] = (D < th[i, 0:3]).mean(0)
            r = stats_one(T[j], T2[j], True, U[i], covs)
            c_true[i], sd_true[i] = r["c"], r["sd"]
        per = [draw(posts, X[ii], cond[ii], N_KEEP["test_k"], gs) for gs in GEN_SCEN]
        for kk, (Tk, _, T2k, okk_k, acc) in enumerate(per):
            acc_k[ii, kk] = np.where(okk_k, acc, np.nan)
            for j, i in enumerate(ii):
                r = stats_one(Tk[j], T2k[j], okk_k[j], U[i], covs)
                c_k[i, kk], sd_k[i, kk] = r["c"], r["sd"]
        for j, i in enumerate(ii):                                     # 참고: 아홉의 균등 혼합(판정에 쓰지 않는다)
            if all(p[3][j] for p in per):
                P, w = mix_pool([p[0][j] for p in per], [p[4][j] for p in per])
                Uv = np.nan_to_num(U[i])
                c_mix[i] = cstat(P[:, 0:3], U[i], np.einsum("ni,ij,nj->n", Uv, covs["main"], Uv), w)[0]
                cover_mix[i] = inside(np.array([wquantile(P[:, d], w, qs) for d in range(3)]).T, i).T
                ok_mix[i] = True
    VAL.mkdir(parents=True, exist_ok=True)
    np.savez(VAL / f"test_{arm}_p{part}of{parts}.npz", idx=idx, cover=cover, cover_m=cover_m, cover_mix=cover_mix, ranks=ranks, okk=okk,
             ok_mix=ok_mix, samp=samp, U=U, c_true=c_true, sd_true=sd_true, c_mix=c_mix, c_k=c_k, sd_k=sd_k, acc_k=acc_k,
             seconds=time.time() - t0)
    print(json.dumps(dict(arm=arm, part=part, parts=parts, n=m_, seconds=round(time.time() - t0, 1))), flush=True)


def _merge(prefix, parts):
    P = [dict(np.load(VAL / f"{prefix}_p{i}of{parts}.npz")) for i in range(parts)]
    idx = np.concatenate([p["idx"] for p in P])
    o = np.argsort(idx)
    out = {}
    for k in (k for k in P[0] if k != "seconds"):
        if P[0][k].ndim == 0:
            continue
        ax = 1 if k == "cover_m" else 0
        out[k] = np.take(np.concatenate([p[k] for p in P], axis=ax), o, axis=ax)
    out["seconds_sum"] = float(sum(float(p["seconds"]) for p in P))
    return out


def test_summary(arm, parts, probe=None):
    """(i) · (iv) 는 참 k 의 조건부 사후 · P1–P6 은 k 아홉의 층 중앙값의 중앙값(아홉 값 모두 적는다) · 균등 혼합은 참고(ZZZg-1)."""
    sim, plan = load_set("test")
    M = _merge(f"test_{arm}", parts)
    n = len(M["idx"])
    lay_all = sim["layer"][M["idx"]]
    okk, cover, cover_m, ranks, samp, U = M["okk"], M["cover"], M["cover_m"], M["ranks"], M["samp"], M["U"]
    th = sim["theta"][M["idx"]]
    c_k, sd_k = M["c_k"], M["sd_k"]
    res = dict(compute_seconds_sum=M["seconds_sum"], n=n, k_names=K_NAMES, directions=U_NAMES_TEST, layers={})
    for L in sorted(set(lay_all)):
        m = (lay_all == L) & okk
        nl = int(m.sum())
        cov = cover[m].mean(0)
        thr = np.array([pp - 3 * math.sqrt(pp * (1 - pp) / max(nl, 1)) for pp in CI])          # 허용폭은 실제로 남은 경우 수 n(ZZZg-5)
        ml = lay_all == L
        mx = ml & M["ok_mix"]
        ent = dict(n=nl, n_reweight_fail=int((ml & ~okk).sum()), coverage=cov, threshold=thr, fail=bool((cov < thr[None, :]).any()),
                   conservative=(cov > np.array(CI)[None, :] + 0.05), coverage_members=cover_m[:, m].mean(1),
                   sbc_rank_hist=[np.histogram(ranks[m, k], bins=20, range=(0, 1))[0] for k in range(3)],
                   C_median_true_k={u: _nanmed(M["c_true"][m, a]) for a, u in enumerate(U_NAMES_TEST)},
                   sd_median_true_k={u: _nanmed(M["sd_true"][m, a]) for a, u in enumerate(U_NAMES_TEST)},
                   C_median_by_k={u: [_nanmed(c_k[ml, kk, a]) for kk in range(len(K_NAMES))] for a, u in enumerate(U_NAMES_TEST)},
                   sd_median_by_k={u: [_nanmed(sd_k[ml, kk, a]) for kk in range(len(K_NAMES))] for a, u in enumerate(U_NAMES_TEST)},
                   n_by_k=[int(np.isfinite(c_k[ml, kk, 0]).sum()) for kk in range(len(K_NAMES))],
                   reference_uniform_mixture=dict(note="참고 — 판정에 쓰지 않는다(ZZZg-1)", n=int(mx.sum()),
                                                  coverage=cover_mix_mean(M["cover_mix"], mx),
                                                  C_median={u: _nanmed(M["c_mix"][mx, a]) for a, u in enumerate(U_NAMES_TEST)}))
        try:
            from sbi.diagnostics.tarp import _run_tarp, get_tarp_references
            good = m & np.isfinite(samp).all((1, 2))
            tt = torch.as_tensor(th[good, 0:3], dtype=torch.float32)
            ecp, _ = _run_tarp(torch.as_tensor(samp[good].transpose(1, 0, 2)), tt, get_tarp_references(tt))
            ent["tarp_ecp"] = ecp.numpy()
        except Exception as e:                                         # 진단 — 통과 기준이 아니다
            ent["tarp_error"] = repr(e)
        res["layers"][L] = ent
    res["n_layers_fail"] = int(sum(v["fail"] for v in res["layers"].values()))
    res["S4_fail_i"] = res["n_layers_fail"] >= 10
    ui = {u: a for a, u in enumerate(U_NAMES_TEST)}
    by_k = lambda m, u, arr=c_k: [_nanmed(arr[m, kk, ui[u]]) for kk in range(len(K_NAMES))]          # 층 중앙값 아홉
    med9 = lambda v: _nanmed([x if x is not None else np.nan for x in v])
    lay = lambda L: (lay_all == L)
    P = {}
    p1 = {u: by_k(lay("match|L0"), u) for u in ("NDI", "FISO", "ODI")}
    P["P1"] = dict(by_k=p1, median_of_k={u: med9(v) for u, v in p1.items()},
                   holds=all((med9(v) if med9(v) is not None else -1) >= 0.75 for v in p1.values()))
    p2 = {f"{k}|{l}": by_k(lay(f"{k}|{l}"), "ODI") for k in S.ACQ_TYPES for l in ("L0", "L1", "L2")}
    P["P2"] = dict(by_k=p2, median_of_k={L: med9(v) for L, v in p2.items()},
                   holds=all((med9(v) if med9(v) is not None else -1) >= 0.75 for v in p2.values()))
    q3 = dict(C_u_relax_TR=by_k(lay("tr|L0"), "u_relax"), C_u_ref_match=by_k(lay("match|L0"), "u_ref"),
              C_ODI_TR=by_k(lay("tr|L0"), "ODI"), C_ODI_match=by_k(lay("match|L0"), "ODI"))
    v3 = {k: med9(v) for k, v in q3.items()}
    P["P3"] = dict(by_k=q3, median_of_k=v3, holds=bool(None not in v3.values() and v3["C_u_relax_TR"] <= v3["C_u_ref_match"] - 0.2
                                                       and abs(v3["C_ODI_TR"] - v3["C_ODI_match"]) <= 0.1))
    q4 = dict(C_NDI_TE=by_k(lay("te|L0"), "NDI"), C_NDI_match=by_k(lay("match|L0"), "NDI"))
    v4 = {k: med9(v) for k, v in q4.items()}
    P["P4"] = dict(by_k=q4, median_of_k=v4, holds=bool(None not in v4.values() and v4["C_NDI_TE"] < v4["C_NDI_match"]))
    condv = sim["cond"][M["idx"]]
    m1 = lay("match|L1") & (condv[:, S.COND_NAMES.index("rho")] == 0)
    m0 = lay("match|L0")
    q5 = dict(sd_L1_rho0=by_k(m1, "NDI", sd_k), sd_L0=by_k(m0, "NDI", sd_k))
    s1, s0 = med9(q5["sd_L1_rho0"]), med9(q5["sd_L0"])
    dsig = math.sqrt(max(s1 ** 2 - s0 ** 2, 0)) if None not in (s1, s0) else None
    srms = float(np.sqrt(np.mean(condv[m1, S.COND_NAMES.index("s_t2")] ** 2))) if m1.any() else None
    expct = 0.397 * srms if srms is not None else None
    P["P5"] = dict(by_k=q5, median_of_k=dict(sd_L1_rho0=s1, sd_L0=s0), dsigma=dsig, expected=expct,
                   holds=bool(None not in (dsig, expct) and s1 > s0 and 0.5 * expct <= dsig <= 1.5 * expct))
    ii = np.flatnonzero(np.isin(lay_all, [f"{k}|L3" for k in ("te", "tetr", "tetr_design")]))
    frac = None
    if len(ii):
        cc = condv[ii]
        men, minn, mt1 = (cc[:, S.COND_NAMES.index(k)] for k in ("m_en", "m_in", "m_t1"))
        t2l = np.tile(CENTRE["t2"], (len(ii), 1)) * np.c_[minn, men, np.ones(len(ii))]
        t1l = np.tile(CENTRE["t1"], (len(ii), 1)) * np.c_[mt1, mt1, np.ones(len(ii))]
        gi = M["idx"][ii]
        u0, _ = apparent_shift(th[ii, 3:6], plan, gi)
        u3, _ = apparent_shift(th[ii, 3:6], plan, gi, (t2l, t1l))
        ok6 = np.isfinite(u0[:, 0]) & np.isfinite(u3[:, 0])
        frac = float(np.mean(np.sign(u0[ok6, 0]) != np.sign(u3[ok6, 0]))) if ok6.any() else None
    P["P6"] = dict(frac_sign_flip=frac, holds=bool(frac is not None and frac > 0.5), note="사후를 쓰지 않는다(대수) — k 와 무관")
    res["P"] = P
    res["S4_fail_P"] = not (P["P1"]["holds"] and P["P2"]["holds"])
    return res


def cover_mix_mean(cover_mix, m):
    return cover_mix[m].mean(0) if m.any() else None


# ── (ii) 신탁 · 참조 사후 ──────────────────────────────────────────────────────
REF_DIM = 12             # w(6: 사이트 A · B 의 무잡음 ν′ · f′_iso · ODI) · A 의 이완(T2_in · r · c) · z(2) · z′(1) — z · z′ 는 L1 에서만 우도에 든다
REF_WIDTH = np.array([0.01] * 6 + [5.0, 0.2, 50.0, 1.0, 1.0, 1.0])                        # 슬라이스 첫 폭(버림 구간에서 조정)
REF_PARAMS = ("s_t2", "rho", "s_t1", "d_t2", "d_t1", "alloc_en", "m_en", "m_in", "m_t1")


def ref_consts(sim, plan, idx):
    """참조 사후의 경우별 상수(행 = 경우) — 관측 · 획득 · T1 시나리오 비 · 이완 수준과 그 매개변수(조건 a 에서)."""
    cond = sim["cond"][idx]
    col = lambda k: cond[:, S.COND_NAMES.index(k)].astype(float)
    scen = np.array(S.SCENARIOS)[np.argmax(np.c_[tuple(col(f"scen_{s}") for s in S.SCENARIOS)], 1)]
    lv = np.array(S.L_LEVELS)[np.argmax(np.c_[tuple(col(f"L_{l}") for l in S.L_LEVELS)], 1)]
    K = dict(obs=sim["oracle"][idx].astype(float), ratio=np.array([S.SCEN[s] for s in scen]), lv=lv)
    for s in ("A", "B"):
        K[f"TE_{s}"], K[f"TR_{s}"] = plan[f"acq_TE_{s}"][idx].astype(float), plan[f"acq_TR_{s}"][idx].astype(float)
    K.update({k: col(k) for k in REF_PARAMS})
    return K


def _relax(state, K):
    """A · B 의 이완 — s4cfg.sample_relax_A · relax_B 와 같은 식 · 무작위 z · z′ 는 상태에서."""
    n = len(state)
    t2in, rr, c = state[:, 6], state[:, 7], state[:, 8]
    t2A = np.c_[t2in, t2in / rr, np.full(n, S.T2ISO)]
    sq = np.sqrt(K["ratio"])
    t1A = np.c_[c * sq, c / sq, np.full(n, S.T1ISO)]
    t2B, t1B = t2A.copy(), t1A.copy()
    lv = K["lv"]
    L1, L2, L3 = lv == "L1", lv == "L2", lv == "L3"
    z0, z1 = state[:, 9], K["rho"] * state[:, 9] + np.sqrt(1 - K["rho"] ** 2) * state[:, 10]
    zt1 = np.clip(state[:, 11], -2.0, 2.0)
    t2B[:, 0] *= np.where(L1, 1 + K["s_t2"] * z0, 1.0)
    t2B[:, 1] *= np.where(L1, 1 + K["s_t2"] * z1, 1.0)
    t1B[:, 0:2] *= np.where(L1, 1 + K["s_t1"] * zt1, 1.0)[:, None]
    en = L2 & (K["alloc_en"] == 1)
    both = L2 & (K["alloc_en"] != 1)
    t2B[:, 1] *= np.where(en, 1 + K["d_t2"], 1.0)
    t1B[:, 1] *= np.where(en, 1 + K["d_t1"], 1.0)
    t2B[:, 0:2] *= np.where(both, 1 + K["d_t2"], 1.0)[:, None]
    t1B[:, 0:2] *= np.where(both, 1 + K["d_t1"], 1.0)[:, None]
    t2B[:, 1] *= np.where(L3, K["m_en"], 1.0)
    t2B[:, 0] *= np.where(L3, K["m_in"], 1.0)
    t1B[:, 0:2] *= np.where(L3, K["m_t1"], 1.0)[:, None]
    return (t2A, t1A), (t2B, t1B)


def ref_theta(state, K, fixed=None):
    """상태 → θ* · log|∂θ*/∂w| · 정의역 안(§13 7). w 에서 겉보기 분율을 거꾸로 풀어 참 분율을 얻는다:
    ν = ν′a_en / (ν′a_en + (1 − ν′)a_in) · f = f′A / (f′A + (1 − f′)a_iso) · A = νa_in + (1 − ν)a_en.
    fixed=((t2, t1) A, (t2, t1) B) 이면 이완을 그 값으로 고정(정답이 알려진 시험)."""
    w = state[:, 0:6]
    rA, rB = _relax(state, K) if fixed is None else fixed
    ok = (state[:, 6] >= S.T2IN_RANGE[0]) & (state[:, 6] <= S.T2IN_RANGE[1]) & (state[:, 7] >= S.R_RANGE[0]) & (state[:, 7] <= S.R_RANGE[1])
    ok &= (state[:, 8] >= S.T1C_RANGE[0]) & (state[:, 8] <= S.T1C_RANGE[1])           # 이완 균일 사전의 지지(고정할 때도 — 평평한 무한 방향을 막는다)
    out, logj = [], np.zeros(len(state))
    with np.errstate(all="ignore"):
        for s, (t2, t1) in (("A", rA), ("B", rB)):
            ws = w[:, 0:3] if s == "A" else w[:, 3:6]
            a = np.exp(-K[f"TE_{s}"][:, None] / t2) * (1 - np.exp(-K[f"TR_{s}"][:, None] / t1))
            a_in, a_en, a_iso = a[:, 0], a[:, 1], a[:, 2]
            nup, fp = ws[:, 0], ws[:, 1]
            ok &= (nup > 0) & (nup < 1) & (fp > 0) & (fp < 1)
            d1 = nup * a_en + (1 - nup) * a_in
            nu = nup * a_en / d1
            A = nu * a_in + (1 - nu) * a_en
            d2 = fp * A + (1 - fp) * a_iso
            out.append(np.c_[nu, fp * A / d2, ws[:, 2]])
            logj += np.log(a_in * a_en) - 2 * np.log(d1) + np.log(A * a_iso) - 2 * np.log(d2)
    th = np.c_[out[1] - out[0], out[0]]
    return th, logj, ok & np.isfinite(th).all(1) & np.isfinite(logj)


def ref_logpost(state, K, sd, fixed=None):
    """참조 사후의 로그 밀도(상수 뺌) — 정확한 가우스 우도 N(obs; w, sd²) + θ* 의 주 사전(절단 균일 → 지지 지시) + 이완 균일 · z · z′ 정규
    + 변수 바꿈의 야코비안. z 는 |z|² ≤ 4 로 절단 · z′ 는 우도에서 ±2 로 자른다(s4cfg.relax_B 와 같다)."""
    th, logj, ok = ref_theta(state, K, fixed)
    z, zp = state[:, 9:11], state[:, 11]
    ok &= (z ** 2).sum(1) <= 4.0
    ok[ok] &= S.in_main(th[ok])
    lp = -0.5 * (((state[:, 0:6] - K["obs"]) / sd) ** 2).sum(1) + logj - 0.5 * (z ** 2).sum(1) - 0.5 * zp ** 2
    return np.where(ok & np.isfinite(lp), lp, -np.inf)


def ref_init(rng, K, fixed=None):
    """사슬 시작점(행 = 사슬) — 사전에서 흩어 뽑는다(ZZZg-2): θ* · 이완 · z · z′ 를 사전에서 뽑고 w = 그 점의 무잡음 겉보기 값."""
    n = len(K["obs"])
    th = S.sample_theta(rng, n, S.D_MAIN)
    z = rng.standard_normal((4 * n + 64, 2))
    z = z[(z ** 2).sum(1) <= 4.0][:n]
    st = np.c_[np.zeros((n, 6)), rng.uniform(*S.T2IN_RANGE, n), rng.uniform(*S.R_RANGE, n), rng.uniform(*S.T1C_RANGE, n), z, rng.standard_normal(n)]
    rA, rB = _relax(st, K) if fixed is None else fixed
    for s, (t2, t1), cen in (("A", rA, th[:, 3:6]), ("B", rB, th[:, 3:6] + th[:, 0:3])):
        a = np.exp(-K[f"TE_{s}"][:, None] / t2) * (1 - np.exp(-K[f"TR_{s}"][:, None] / t1))
        nu_app, f_app = S.C.apparent(cen[:, 0], cen[:, 1], (a[:, 0], a[:, 1], a[:, 2]))
        st[:, (0 if s == "A" else 3):(3 if s == "A" else 6)] = np.c_[nu_app, f_app, cen[:, 2]]
    return st


def ref_mcmc(K, seed, sd=None, fixed=None, chains=MCMC["chains"], warmup=MCMC["warmup"], keep=MCMC["keep"]):
    """참조 사후(§7 ii · §13 7 · ZZZg-2) — sbi 0.27.0 의 슬라이스 표집기(SliceSamplerVectorized — MCMCPosterior 의 slice_np_vectorized 와 같은 것)에
    사용자 퍼텐셜(ref_logpost)을 준다. 경우 묶음의 사슬(경우 × chains)을 한 번에 돌린다 · 버림 warmup(폭 조정 구간) · 보존 keep.
    반환: θ* (경우, chains, keep, 6)."""
    from sbi.samplers.mcmc.slice_numpy import SliceSamplerVectorized
    from simulate import ORACLE_SD
    sd = ORACLE_SD if sd is None else sd
    n = len(K["obs"])
    KC = {k: np.repeat(v, chains, axis=0) for k, v in K.items()}                    # 사슬 = 경우마다 chains 개 · 경우 순서대로
    fx = None if fixed is None else tuple(tuple(np.repeat(a, chains, axis=0) for a in r) for r in fixed)
    x0 = ref_init(np.random.default_rng(seed), KC, fx)
    np.random.seed(seed % (2 ** 32))                                                  # sbi 표집기는 np.random 을 쓴다
    smp = SliceSamplerVectorized(lambda s: ref_logpost(s, KC, sd, fx), x0, num_chains=n * chains, tuning=warmup, verbose=False,
                                 init_width=REF_WIDTH).run(keep)                       # (n·chains, keep, 12)
    KK = {k: np.repeat(v, keep, axis=0) for k, v in KC.items()}
    fk = None if fx is None else tuple(tuple(np.repeat(a, keep, axis=0) for a in r) for r in fx)
    th, _, _ = ref_theta(smp.reshape(-1, REF_DIM), KK, fk)
    return th.reshape(n, chains, keep, 6)


def _ranknorm(x):
    from scipy.special import ndtri
    from scipy.stats import rankdata
    r = rankdata(x.ravel(), method="average").reshape(x.shape)
    return ndtri((r - 0.375) / (x.size + 0.25))


def _split(x):
    h = x.shape[1] // 2
    return np.concatenate([x[:, :h], x[:, x.shape[1] - h:]], 0)


def _rhat_raw(x):
    m, n = x.shape
    W = x.var(1, ddof=1).mean()
    B = n * x.mean(1).var(ddof=1)
    return float(np.sqrt(((n - 1) / n * W + B / n) / W))


def _ess_raw(x):
    """ESS — 사슬 사이 자기공분산 평균 · Geyer 초기 단조열(Stan · ArviZ 의 식)."""
    m, n = x.shape
    xc = x - x.mean(1, keepdims=True)
    f = np.fft.rfft(xc, n=2 * n, axis=1)
    acov = np.fft.irfft(f * np.conj(f), axis=1)[:, :n] / n
    chain_var = acov[:, 0] * n / (n - 1)
    mean_var = chain_var.mean()
    var_plus = mean_var * (n - 1) / n + (x.mean(1).var(ddof=1) if m > 1 else 0.0)
    rho = np.zeros(n)
    mac = acov.mean(0)
    rho[0] = 1.0
    even, odd = 1.0, 1 - (mean_var - mac[1]) / var_plus
    rho[1] = odd
    s = 1
    while s < n - 4 and even + odd > 0:
        even = 1 - (mean_var - mac[s + 1]) / var_plus
        odd = 1 - (mean_var - mac[s + 2]) / var_plus
        if even + odd >= 0:
            rho[s + 1], rho[s + 2] = even, odd
        s += 2
    max_s = s
    if even > 0:
        rho[max_s + 1] = even
    for t in range(1, max_s - 2, 2):
        if rho[t + 1] + rho[t + 2] > rho[t - 1] + rho[t]:
            rho[t + 1] = rho[t + 2] = (rho[t - 1] + rho[t]) / 2
    tau = -1 + 2 * rho[:max_s].sum() + rho[max_s + 1]
    tau = max(tau, 1 / np.log10(m * n))
    return float(m * n / tau)


def rhat_ess(ch):
    """ch (사슬, 표본) — (R̂, bulk ESS) · R̂ = 순위 정규화 분할 R̂ 의 벌크 · 접힘 중 큰 것 · bulk ESS = 순위 정규화 분할 사슬의 ESS(§13 7)."""
    z = _ranknorm(_split(ch))
    zf = _ranknorm(_split(np.abs(ch - np.median(ch))))
    return max(_rhat_raw(z), _rhat_raw(zf)), _ess_raw(z)


def ref_posterior_is(rng, cond_row, plan, i, obs, n=1_000_000, sd=None):
    """옛 참조 사후(사전 표본 + 정확한 가우스 우도의 중요도 추출) — ZZZg-2 로 MCMC 로 바꿨다. ESS 참고 확인과 시험의 대조용으로만 남긴다."""
    from simulate import ORACLE_SD
    sd = ORACLE_SD if sd is None else sd
    th = S.sample_theta(rng, n, S.D_MAIN)
    scen = S.SCENARIOS[int(np.argmax([cond_row[S.COND_NAMES.index(f"scen_{s}")] for s in S.SCENARIOS]))]
    lv = S.L_LEVELS[int(np.argmax([cond_row[S.COND_NAMES.index(f"L_{l}")] for l in S.L_LEVELS]))]
    t2A, t1A = S.sample_relax_A(rng, n, np.full(n, scen))
    p = {k: np.full(n, cond_row[S.COND_NAMES.index(k)]) for k in REF_PARAMS}
    t2B, t1B = S.relax_B(rng, np.full(n, lv), p, t2A, t1A) if lv != "L1" else _relax_B_fast(rng, p, t2A, t1A)
    mu = []
    for s, cen, t2, t1 in (("A", th[:, 3:6], t2A, t1A), ("B", th[:, 3:6] + th[:, 0:3], t2B, t1B)):
        a = np.exp(-plan[f"acq_TE_{s}"][i] / t2) * (1 - np.exp(-plan[f"acq_TR_{s}"][i] / t1))
        nu_app, f_app = S.C.apparent(cen[:, 0], cen[:, 1], (a[:, 0], a[:, 1], a[:, 2]))
        mu.append(np.c_[nu_app, f_app, cen[:, 2]])
    mu = np.concatenate(mu, 1)
    ll = -0.5 * (((mu - obs[None, :]) / sd) ** 2).sum(1)
    w = np.exp(ll - ll.max())
    ess = w.sum() ** 2 / (w ** 2).sum()
    return th, w / w.sum(), ess


def _relax_B_fast(rng, p, t2A, t1A):
    """L1 의 무작위 z · z′ — relax_B 와 같은 분포를 벡터로(옛 참조 사후의 10⁶ 표본)."""
    n = len(t2A)
    z = rng.standard_normal((n * 3, 2))
    z = z[(z ** 2).sum(1) <= 4.0][:n]
    while len(z) < n:
        zz = rng.standard_normal((n, 2))
        z = np.vstack([z, zz[(zz ** 2).sum(1) <= 4.0]])[:n]
    rho = p["rho"]
    zc = np.c_[z[:, 0], rho * z[:, 0] + np.sqrt(1 - rho ** 2) * z[:, 1]]
    zt1 = np.clip(rng.standard_normal(n), -2, 2)
    t2B, t1B = t2A.copy(), t1A.copy()
    t2B[:, 0] *= 1 + p["s_t2"] * zc[:, 0]
    t2B[:, 1] *= 1 + p["s_t2"] * zc[:, 1]
    t1B[:, 0:2] *= (1 + p["s_t1"] * zt1)[:, None]
    return t2B, t1B


def oracle_part(part=0, parts=1, probe=None, batch=25):
    """신탁 참조 경우의 한 몫 — 참조 사후(MCMC) · 수렴 진단 · 신탁 NPE(참 k) · C2ST → WORK/s4/val(§7 ii · §13 7)."""
    from sbi.utils.metrics import c2st
    torch.manual_seed(3000 + part)
    sim, plan = load_set("oracle_ref")
    n = len(sim["theta"]) if probe is None else probe
    idx = np.arange(n)[part::parts]
    post = load_posts("oracle", 1)[0]
    rng = np.random.default_rng(2 + part)
    covD = cov_prior(np.random.default_rng(2))[0:3, 0:3]
    rows = []
    t0 = time.time()
    for b0 in range(0, len(idx), batch):
        ii = idx[b0:b0 + batch]
        TH = ref_mcmc(ref_consts(sim, plan, ii), seed=20261103 + int(ii[0]))
        for j, i in enumerate(ii):
            diag = [rhat_ess(TH[j, :, :, d]) for d in range(6)]
            row = dict(i=int(i), kind=str(sim["acq_kind"][i]), rhat_max=max(r for r, _ in diag), ess_bulk_min=min(e for _, e in diag))
            row["converged"] = bool(row["rhat_max"] <= MCMC["rhat"] and row["ess_bulk_min"] >= MCMC["ess"])
            if row["converged"]:
                R = TH[j].reshape(-1, 6)
                U = directions(sim["theta"][i:i + 1], plan, np.array([i]))[0]
                Uv = np.nan_to_num(U)
                pv = np.einsum("ni,ij,nj->n", Uv, covD, Uv)
                row["C_ref"] = cstat(R[:, 0:3], U, pv)[0]
                xs = torch.as_tensor(np.c_[sim["oracle"][i][None], sim["cond"][i][None]], dtype=torch.float32)
                got = []
                for _ in range(40):
                    with torch.no_grad():
                        smp = post.sample((2000,), x=xs, show_progress_bars=False).numpy()
                    got.append(smp[S.in_main(smp)])
                    if sum(len(a) for a in got) >= 1000:
                        break
                G = np.concatenate(got)[:1000]
                if len(G) >= 1000:
                    row["C_npe"] = cstat(G[:, 0:3], U, pv)[0]
                    Rs = R[rng.choice(len(R), 1000, replace=False)]
                    row["c2st"] = float(c2st(torch.as_tensor(G, dtype=torch.float32), torch.as_tensor(Rs, dtype=torch.float32)))
            rows.append(row)
    VAL.mkdir(parents=True, exist_ok=True)
    json.dump(dict(rows=rows, seconds=time.time() - t0), open(VAL / f"oracle_p{part}of{parts}.json", "w"),
              default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o))
    print(json.dumps(dict(step="oracle", part=part, parts=parts, n=len(idx), seconds=round(time.time() - t0, 1))), flush=True)


def oracle_summary(parts):
    """(ii) 판정 — 수렴하지 못한 경우를 빼고 그 수를 적는다 · 20 % 넘게 빠지면 «판정 불가» → «기술적 불가능»(ZZZg-2)."""
    P = [json.load(open(VAL / f"oracle_p{i}of{parts}.json")) for i in range(parts)]
    rows = sorted((r for p in P for r in p["rows"]), key=lambda r: r["i"])
    n = len(rows)
    n_ex = sum(1 for r in rows if not r["converged"])
    res = dict(compute_seconds_sum=sum(p["seconds"] for p in P), n=n, n_excluded_not_converged=n_ex, frac_excluded=n_ex / max(n, 1),
               mcmc=MCMC, rows=rows)
    if n_ex > MCMC["max_excluded"] * n:
        res.update(undecidable_ii=True, technical_impossible=True, **{"pass": None}, S4_fail_ii=None)
        return res
    sel = [r for r in rows if "C_npe" in r]
    da = [abs(r["C_npe"][3] - r["C_ref"][3]) for r in sel if r["kind"] in ("tr", "tetr") and np.isfinite(r["C_npe"][3]) and np.isfinite(r["C_ref"][3])]
    res.update(undecidable_ii=False, technical_impossible=False, n_npe_short=sum(1 for r in rows if r["converged"] and "C_npe" not in r),
               a_median_abs_diff_C_u_relax=float(np.median(da)) if da else None,
               b_median_C_ODI=float(np.median([r["C_npe"][2] for r in sel])) if sel else None,
               c_median_c2st=float(np.median([r["c2st"] for r in sel])) if sel else None)
    res["pass"] = bool(res["a_median_abs_diff_C_u_relax"] is not None and res["a_median_abs_diff_C_u_relax"] <= 0.10
                       and (res["b_median_C_ODI"] or 0) >= 0.9 and (res["c_median_c2st"] or 1) <= 0.65)
    res["S4_fail_ii"] = not res["pass"]
    return res


# ── (v) 오지정 점수 · 주입 ─────────────────────────────────────────────────────
def mix_score(logq):
    """오지정 점수 = −log[(1/9) Σ_k q(x | a, k)] — logq (n, 9)(ZZZg-4)."""
    m = logq.max(1, keepdims=True)
    return -(m[:, 0] + np.log(np.exp(logq - m).mean(1)))


def q_scores(sim):
    """실자료에서 k 를 모르므로 아홉 k 의 균등 혼합 밀도로 낸다(ZZZg-4 · §7 v) — 문턱 · 오경보 · 주입 모두 이 점수."""
    est = torch.load(S.WORK / "models" / "q_0.obj.pt", weights_only=False)
    x = torch.as_tensor(sim["x_c3"], dtype=torch.float32)
    L = np.zeros((len(x), len(GEN_SCEN)))
    for kk, (g, s) in enumerate(GEN_SCEN):
        a = torch.as_tensor(with_gs(sim["cond"], g, s), dtype=torch.float32)
        for b0 in range(0, len(x), 4096):
            with torch.no_grad():
                L[b0:b0 + 4096, kk] = est.log_prob(x[b0:b0 + 4096][None], condition=a[b0:b0 + 4096]).reshape(-1).numpy()
    return mix_score(L)


def run_miss():
    sim, _ = load_set("test")
    sc = q_scores(sim)
    res = dict(score="−log[(1/9) Σ_k q(x | a, k)]", layers={})
    rng = np.random.default_rng(3)
    for L in sorted(set(sim["layer"])):
        idx = np.flatnonzero(sim["layer"] == L)
        idx = rng.permutation(idx)
        h = len(idx) // 2
        thr = float(np.quantile(sc[idx[:h]], 0.99))
        res["layers"][L] = dict(threshold=thr, false_alarm=float(np.mean(sc[idx[h:]] > thr)))
    res["global_threshold"] = float(np.quantile(sc, 0.99))
    res["false_alarm_max"] = max(v["false_alarm"] for v in res["layers"].values())
    res["false_alarm_ok"] = res["false_alarm_max"] <= 0.02
    return res, sc


def layer_of(kind, dA, dB, lv):
    k = {"match": "match", "te": "te", "tr": "tr", "free": "tetr" if dA == dB else "tetr_design"}[kind]
    return f"{k}|{lv}"


def run_inject(miss, probe=None):
    """주입 시험(기술만) — 탐지율 · 탐지되지 않은 경우의 90 % 포함률(참 k 의 조건부 사후 — (i) 과 같다)."""
    torch.manual_seed(4000)
    sim, plan = load_set("inject")
    n = len(sim["theta"])
    plan = {k: v[:n] for k, v in plan.items()}
    sc = q_scores(sim)
    out = {}
    lay = np.array([layer_of(k, a, b, l) for k, a, b, l in zip(plan["acq_kind"], plan["acq_design_A"], plan["acq_design_B"], plan["level"])])
    glob = miss["global_threshold"]
    thr = np.array([miss["layers"].get(L, {}).get("threshold", glob) for L in lay])          # 층이 없으면 전체 99 백분위
    flag = sc > thr
    posts = load_posts("c3", N_ENS)
    for name in np.unique(sim["cell"]):
        m = sim["cell"] == name
        ii = np.flatnonzero(m & ~flag)
        cov = []
        for b0 in range(0, len(ii), 250):
            jj = ii[b0:b0 + 250]
            T, _, _, ok, _ = draw(posts, sim["x_c3"][jj], sim["cond"][jj], N_KEEP["inject"])
            for j, i in enumerate(jj):
                if ok[j]:
                    lo, hi = np.quantile(T[j, :, 0:3], [0.05, 0.95], axis=0)
                    cov.append((sim["theta"][i, 0:3] >= lo) & (sim["theta"][i, 0:3] <= hi))
        out[str(name)] = dict(n=int(m.sum()), detected=float(flag[m].mean()), n_undetected=int(len(ii)),
                              coverage90_undetected=np.mean(cov, 0) if cov else None)
    return out


# ── 카드 · 지도 · 선형 기준선 ─────────────────────────────────────────────────
def fisher_cell(x, th, covT):
    """§13 12 — x (M, d) 를 θ* (M, 6) 에 선형 회귀 → 사후 공분산."""
    A = np.c_[np.ones(len(th)), th]
    beta, *_ = np.linalg.lstsq(A, x, rcond=None)
    J = beta[1:].T
    R = np.cov((x - A @ beta).T) + 1e-9 * np.eye(x.shape[1])
    F = J.T @ np.linalg.solve(R, J)
    return np.linalg.inv(F + np.linalg.inv(covT))


def cells_part(name, arm, part=0, parts=1, probe=None, batch=100):
    """카드 · 지도 자료의 한 몫 — 관측마다 k 아홉 각각의 c · sd(카드는 사전 민감도 ×0.5 · ×2 도) · 참고 균등 혼합의 c → WORK/s4/val(§5 · §7 iii · §9 · ZZZg-1 · 3)."""
    torch.manual_seed(2000 + 100 * ("cards", "maps").index(name) + part)
    sim, plan = load_set(name)
    sens = name == "cards"                                              # (iii) 은 카드 자료에서만(ZZZg-3)
    n = len(sim["theta"]) if probe is None else probe
    covs = prior_covs(np.random.default_rng(4), sens)
    posts = load_posts(arm, N_ENS)
    idx = np.arange(n)[part::parts]
    X, cond, th = xobs(sim, arm)[idx], sim["cond"][idx], sim["theta"][idx]
    U = directions(th, plan, idx)
    m_, nk, nu = len(idx), len(GEN_SCEN), U.shape[1]
    out = {k: np.full((m_, nk, nu), np.nan) for k in (("c", "sd", "c_half", "c_two") if sens else ("c", "sd"))}
    n_half = np.zeros((m_, nk), int)
    c_mix = np.full((m_, nu), np.nan)
    t0 = time.time()
    for b0 in range(0, m_, batch):
        ii = np.arange(b0, min(m_, b0 + batch))
        per = [draw(posts, X[ii], cond[ii], N_KEEP[name], gs) for gs in GEN_SCEN]
        for kk, (Tk, _, T2k, okk_k, _) in enumerate(per):
            for j, i in enumerate(ii):
                r = stats_one(Tk[j], T2k[j], okk_k[j], U[i], covs, sens)
                for k in out:
                    out[k][i, kk] = r[k]
                if sens:
                    n_half[i, kk] = r["n_half"]
        for j, i in enumerate(ii):                                     # 참고: 아홉의 균등 혼합
            if all(p[3][j] for p in per):
                P, w = mix_pool([p[0][j] for p in per], [p[4][j] for p in per])
                Uv = np.nan_to_num(U[i])
                c_mix[i] = cstat(P[:, 0:3], U[i], np.einsum("ni,ij,nj->n", Uv, covs["main"], Uv), w)[0]
    VAL.mkdir(parents=True, exist_ok=True)
    np.savez(VAL / f"{name}_{arm}_p{part}of{parts}.npz", idx=idx, U=U, n_half=n_half, c_mix=c_mix, seconds=time.time() - t0,
             **{f"k_{k}": v for k, v in out.items()})
    print(json.dumps(dict(name=name, arm=arm, part=part, parts=parts, n=m_, seconds=round(time.time() - t0, 1))), flush=True)


def worst_of_k(Ck, MDEk):
    """k 아홉의 칸 값 → 분류용 값(가장 낮은 C · 가장 큰 MDE) · 중앙값 · «가정 의존»(k 마다의 분류가 갈림) — ZZZg-1."""
    C = np.array([np.nan if v is None else v for v in Ck], float)
    E = np.array([np.nan if v is None else v for v in MDEk], float)
    if not np.isfinite(C).all():
        return None
    cls_k = [str(x) for x in classify(C)]
    return dict(C_worst=float(C.min()), MDE_worst=float(np.nanmax(E)), C_median=float(np.median(C)), MDE_median=float(np.nanmedian(E)),
                cls=str(classify(np.array([C.min()]))[0]), cls_sens=str(classify(np.array([C.min()]), 0.5, 0.1)[0]),
                cls_by_k=cls_k, assumption_dependent=len(set(cls_k)) > 1)


def cells_summary(name, arm, parts):
    """칸마다 · 방향마다: k 별 C · σ · MDE → 분류는 아홉 중 가장 나쁜 것 · 중앙값을 함께 · 분류가 갈리면 «가정 의존» ·
    카드는 ×0.5 · ×2 에서 (가장 나쁜 k 의) 분류가 바뀌면 «사전 의존»(§7 iii) · C^lin 은 k 를 주변화한 선형 기준선이라 같은 가정의 균등 혼합과만 비교한다(§6 · §13 12 · ZZZi v3-2)."""
    sim, plan = load_set(name)
    M = _merge(f"{name}_{arm}", parts)
    idx = M["idx"]
    sens = name == "cards"
    U = M["U"]
    X, th = xobs(sim, arm)[idx], sim["theta"][idx]
    covT = cov_prior(np.random.default_rng(4))
    covD = covT[0:3, 0:3]
    cells = {}
    for cl in sorted(set(sim["cell"][idx])):
        m = sim["cell"][idx] == cl
        Sig = fisher_cell(X[m], th[m], covT)
        ent = {}
        for a, u in enumerate(U_NAMES):
            if not np.isfinite(M["k_c"][m, :, a]).any():
                continue
            Ck = [_nanmed(M["k_c"][m, kk, a]) for kk in range(len(K_NAMES))]
            sdk = [_nanmed(M["k_sd"][m, kk, a]) for kk in range(len(K_NAMES))]
            MDEk = [None if s is None else 2.80 * s for s in sdk]
            wk = worst_of_k(Ck, MDEk)
            uu = np.nanmean(U[m, a], 0)
            uu = uu / np.linalg.norm(uu)
            e = dict(C_by_k=Ck, sd_by_k=sdk, MDE_by_k=MDEk, n=int(np.isfinite(M["k_c"][m, 0, a]).sum()),
                     C_lin=float(1 - (uu @ Sig[0:3, 0:3] @ uu) / (uu @ covD @ uu)),
                     reference_uniform_mixture_C=_nanmed(M["c_mix"][m, a]))
            if e["reference_uniform_mixture_C"] is not None:                       # 같은 가정(k 정보 없음)의 차 — §6 · ZZZi v3-2
                e["C_mix_minus_C_lin"] = e["reference_uniform_mixture_C"] - e["C_lin"]
            flags = []
            if wk is None:
                e["cls"] = "재가중 불가"
            else:
                e.update(wk)
                if wk["assumption_dependent"]:
                    flags.append("가정 의존")
            if sens and wk is not None:
                for pr in ("half", "two"):
                    Cp = [_nanmed(M[f"k_c_{pr}"][m, kk, a]) for kk in range(len(K_NAMES))]
                    e[f"C_{pr}_by_k"] = Cp
                    if any(v is None for v in Cp):
                        e[f"cls_{pr}"] = "재가중 불가"
                        flags.append(f"사전 민감도 ×{'0.5' if pr == 'half' else '2'} 재가중 불가")
                    else:
                        e[f"cls_{pr}"] = str(classify(np.array([min(Cp)]))[0])
                        if e[f"cls_{pr}"] != wk["cls"]:
                            flags.append("사전 의존")
                e["n_half_median"] = float(np.median(M["n_half"][m]))
            e["flags"] = sorted(set(flags))
            ent[u] = e
        cells[cl] = ent
    return dict(compute_seconds_sum=M["seconds_sum"], n=len(idx), k_names=K_NAMES, cells=cells)


def main():
    """  test | cards | maps  <arm> --part i --parts n        (몫 계산)
         test-summary | cards-summary | maps-summary <arm> --parts n
         oracle --part i --parts n · oracle-summary --parts n
         miss · inject                                       (작다 — 한 번에)
         공통: --probe N(앞 N 경우만 · 점검)"""
    if torch.cuda.is_available():
        raise SystemExit("⛔ CUDA 가 보인다 — CPU 만(PLAN §9-1 #49). CUDA_VISIBLE_DEVICES= 로 다시 돌려라")
    a = sys.argv[1:]
    opt = lambda k, d=None: int(a[a.index(k) + 1]) if k in a else d
    probe, part, parts = opt("--probe"), opt("--part", 0), opt("--parts", 1)
    tag = "" if probe is None else f"_probe{probe}"
    VAL.mkdir(parents=True, exist_ok=True)
    step = a[0]
    if step == "test":
        test_part(a[1], part, parts, probe)
    elif step in ("cards", "maps"):
        cells_part(step, a[1], part, parts, probe)
    elif step == "test-summary":
        jdump(f"test_{a[1]}" + tag, test_summary(a[1], parts, probe))
    elif step in ("cards-summary", "maps-summary"):
        nm = step.split("-")[0]
        jdump(f"{nm}_{a[1]}" + tag, cells_summary(nm, a[1], parts))
    elif step == "oracle":
        oracle_part(part, parts, probe)
    elif step == "oracle-summary":
        jdump("oracle" + tag, oracle_summary(parts))
    elif step in ("miss", "inject"):
        miss, _ = run_miss()
        jdump("miss" + tag, miss)
        if step == "inject":
            jdump("inject" + tag, run_inject(miss, probe))
    else:
        raise SystemExit(step)


if __name__ == "__main__":
    main()
