"""S4 — 자료별 경우 목록(θ* · 조건 · 성가신 변수) — SPEC_S4 §3 · §7 · §9. numpy 만 · 결정적(자료마다 시드 고정 — §3-8).

자료: train(학습 · 학습 사전) · test(20 층 × 1,000 · 주 사전) · oracle_train(신탁 NPE 학습) · oracle_ref(획득 유형 5 × 100 · 참조 사후)
      · cards(ON-Harmony 28 쌍 · 주 칸 + 민감도 칸 · 관측 M = 100) · maps(ΔTE × ΔTR 7 × 7 × L 4 · M = 100) · inject(학습 범위 밖 넷 × 1,000).
    python plan.py <자료> [n]       # WORK/s4/plan_<자료>.npz
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s4cfg as S  # noqa: E402

ONH = {  # SPEC_v2 §3 · decomp/stage1/analyze2.py SCANNERS(LOCK_B EN :7748–7754) — 설계 · TE 범위(ms) · TR 범위(s) · 코일(LOCK_B EN :2904–2905 — 나머지는 미상)
    "OXF1PRI": ("S", (92.0, 92.0), (3.60, 3.60), 32), "OXF2PRI": ("S", (92.0, 92.0), (3.60, 3.60), 64),
    "OXF3TRI": ("S", (96.4, 96.4), (3.60, 3.60), None), "NOT2ING-89": ("V", (89.0, 89.0), (4.00, 4.38), None),
    "NOT2ING-98": ("V", (98.0, 98.0), (4.00, 4.38), None), "NOT4GEP": ("V", (68.0, 68.0), (4.90, 4.90), None),
    "OXF4GEP": ("V", (68.3, 71.0), (8.00, 8.00), None), "NOT1ACH": ("V", (70.0, 71.2), (9.81, 10.62), None),
}
MAP_DTE = np.linspace(-40.0, 40.0, 7)
MAP_DTR = np.linspace(-7000.0, 7000.0, 7)
MAP_REF = dict(TE=80.0, TR=6000.0, design="V", snr=30.0, neff=50.0)


def _finish(rng, theta, acq, lv, p, gen, scen, extra=None):
    n = len(theta)
    t2A, t1A = S.sample_relax_A(rng, n, scen)
    t2B, t1B = S.relax_B(rng, lv, p, t2A, t1A)
    dpar = rng.choice(S.DPAR_TRUE, n)
    out = dict(theta=theta, cond=S.cond_vector(acq, lv, p, gen, scen), level=lv, gen=np.asarray(gen), scen=np.asarray(scen),
               t2A=t2A, t1A=t1A, t2B=t2B, t1B=t1B, dpar=dpar, fdir_A=rng.integers(0, len(S.C.FIBRE_DIRS), n),
               fdir_B=rng.integers(0, len(S.C.FIBRE_DIRS), n), bscale=np.ones(n), **{f"acq_{k}": v for k, v in acq.items()})
    out.update(extra or {})
    return out


def _draw(rng, n, theta_half, kinds=None, levels=None):
    theta = S.sample_theta(rng, n, theta_half)
    acq = S.sample_acq(rng, n, kinds)
    lv, p = S.sample_level(rng, n, levels)
    gen = rng.choice(S.GENERATORS, n)
    scen = rng.choice(S.SCENARIOS, n)
    return theta, acq, lv, p, gen, scen


def train(n):
    rng = np.random.default_rng(S.SEEDS["train"])
    return _finish(rng, *_draw(rng, n, S.D_TRAIN))


def test(n_per=1000):
    rng = np.random.default_rng(S.SEEDS["test"])
    kinds = np.repeat(S.ACQ_TYPES, len(S.L_LEVELS) * n_per)
    levels = np.tile(np.repeat(S.L_LEVELS, n_per), len(S.ACQ_TYPES))
    return _finish(rng, *_draw(rng, len(kinds), S.D_MAIN, kinds, levels), extra=dict(layer=np.char.add(np.char.add(kinds, "|"), levels)))


def oracle_train(n):
    rng = np.random.default_rng(S.SEEDS["oracle"])
    return _finish(rng, *_draw(rng, n, S.D_TRAIN))


def oracle_ref(n_per=100):
    rng = np.random.default_rng(S.SEEDS["oracle"] + 1)
    kinds = np.repeat(S.ACQ_TYPES, n_per)
    return _finish(rng, *_draw(rng, len(kinds), S.D_MAIN, kinds))


def _pair_acq(rng, ua, ub, m, snr, neff):
    da, tea, tra, ca = ONH[ua]
    db, teb, trb, cb = ONH[ub]
    coil = lambda c: np.full(m, c) if c else rng.choice(S.COILS, m)
    return dict(kind=np.full(m, "onh"), TE_A=rng.uniform(*tea, m), TE_B=rng.uniform(*teb, m),
                TR_A=rng.uniform(*tra, m) * 1000, TR_B=rng.uniform(*trb, m) * 1000, design_A=np.full(m, da), design_B=np.full(m, db),
                coil_A=coil(ca), coil_B=coil(cb), snr_A=np.full(m, snr), snr_B=np.full(m, snr), neff=np.full(m, float(neff)))


def cards(M=100):
    """ON-Harmony 28 쌍 × (주 칸 L0–L3 · SNR 30 · N_eff 50) + (민감도 칸 L0 · SNR {20, 50} × N_eff {25, 100}) — §9."""
    rng = np.random.default_rng(S.SEEDS["cards"])
    units = list(ONH)
    pairs = [(a, b) for i, a in enumerate(units) for b in units[i + 1:]]
    cells = [(l, 30.0, 50.0) for l in S.L_LEVELS] + [("L0", s, ne) for s in (20.0, 50.0) for ne in (25.0, 100.0)]
    parts, tags = [], []
    for (a, b) in pairs:
        for (l, snr, ne) in cells:
            acq = _pair_acq(rng, a, b, M, snr, ne)
            lv, p = S.sample_level(rng, M, np.full(M, l))
            gen, scen = rng.choice(S.GENERATORS, M), rng.choice(S.SCENARIOS, M)
            parts.append(_finish(rng, S.sample_theta(rng, M, S.D_MAIN), acq, lv, p, gen, scen))
            tags += [f"{a}~{b}|{l}|SNR{int(snr)}|Neff{int(ne)}"] * M
    out = {k: np.concatenate([q[k] for q in parts]) for k in parts[0]}
    out["cell"] = np.array(tags)
    return out


def maps(M=100):
    """ΔTE × ΔTR 7 × 7 × L 4 — 기준 획득 TE 80 · TR 6 s · 설계 V · SNR 30 · N_eff 50(§9)."""
    rng = np.random.default_rng(S.SEEDS["maps"])
    parts, tags = [], []
    for dte in MAP_DTE:
        for dtr in MAP_DTR:
            for l in S.L_LEVELS:
                r = MAP_REF
                acq = dict(kind=np.full(M, "map"), TE_A=np.full(M, r["TE"]), TE_B=np.full(M, r["TE"] + dte),
                           TR_A=np.full(M, r["TR"]), TR_B=np.full(M, r["TR"] + dtr), design_A=np.full(M, r["design"]),
                           design_B=np.full(M, r["design"]), coil_A=rng.choice(S.COILS, M), coil_B=rng.choice(S.COILS, M),
                           snr_A=np.full(M, r["snr"]), snr_B=np.full(M, r["snr"]), neff=np.full(M, r["neff"]))
                lv, p = S.sample_level(rng, M, np.full(M, l))
                gen, scen = rng.choice(S.GENERATORS, M), rng.choice(S.SCENARIOS, M)
                parts.append(_finish(rng, S.sample_theta(rng, M, S.D_MAIN), acq, lv, p, gen, scen))
                tags += [f"dTE{dte:+.0f}|dTR{dtr / 1000:+.2f}|{l}"] * M
    out = {k: np.concatenate([q[k] for q in parts]) for k in parts[0]}
    out["cell"] = np.array(tags)
    return out


def inject(n=1000):
    """학습 범위 밖 넷(§7 v): SNR 10 · 조직 d∥ 3.4 · b 값 배율 ±5 %(합성만) · N_eff 400."""
    rng = np.random.default_rng(S.SEEDS["inject"])
    parts, tags = [], []
    for name in ("snr10", "dpar3.4", "bscale", "neff400"):
        th, acq, lv, p, gen, scen = _draw(rng, n, S.D_MAIN)
        if name == "snr10":
            acq["snr_A"][:] = acq["snr_B"][:] = 10.0
        if name == "neff400":
            acq["neff"][:] = 400.0
        q = _finish(rng, th, acq, lv, p, gen, scen)
        if name == "dpar3.4":
            q["dpar"][:] = 3.4
        if name == "bscale":
            q["bscale"][:] = rng.choice((0.95, 1.05), n)
        parts.append(q)
        tags += [name] * n
    out = {k: np.concatenate([q[k] for q in parts]) for k in parts[0]}
    out["cell"] = np.array(tags)
    return out


MAKERS = dict(train=train, test=test, oracle_train=oracle_train, oracle_ref=oracle_ref, cards=cards, maps=maps, inject=inject)


def make(name, n=None):
    f = MAKERS[name]
    d = f(n) if n is not None else f()
    S.WORK.mkdir(parents=True, exist_ok=True)
    path = S.WORK / f"plan_{name}.npz"
    np.savez(path, **d)
    return path, len(d["theta"])


if __name__ == "__main__":
    nm = sys.argv[1]
    p, n = make(nm, int(sys.argv[2]) if len(sys.argv) > 2 else None)
    print(p.name, n)
