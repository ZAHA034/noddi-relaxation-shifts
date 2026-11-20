"""design_sim.checks — the checks that must be able to fail, K-0, K-a – K-e (ANALYSIS_PLAN.md Section 7), run before the main run; each must pass on the clean code and fail on its planted fault, and both must hold for exit status 0."""
from __future__ import annotations

import copy
import json
import os
import pathlib
import sys
import time

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_k, "1")
import numpy as np  # noqa: E402

if __package__ in (None, ""):
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from design_sim import analyses as A  # noqa: E402
from design_sim import model as M  # noqa: E402

S75_BETA = {"b": -0.0006369086121453225, "e": -0.0017227855902291478}     # ledger/derived/s75_values.json power.beta_per_ms
S75_POWER = {"S2b": 0.13619266744099404, "S2e": 0.633805758090803}         # power.by_te.mid … 0.08
K = len(M.ACQ)


def cell(D="D0", S="S0", V=0.0, H="H0", N="bal", idx=0):
    return dict(idx=idx, D=D, S=S, V=V, H=H, N=N)


def k0(factor=M.LENGTHEN):
    dev = max(abs(M.delta(a, factor) - S75_BETA[a]) for a in ("b", "e"))
    return dict(value=dev, ok=dev < 1e-12, criterion="|Δ − s75 기울기| < 1e−12 (b, e)")


def ka(te_term=True):
    worst = 0.0
    for D in M.DESIGNS:
        for S in M.SCENARIOS[1:]:
            for N in M.NS:
                c = cell(D, S, N=N)
                d = M.generate(c, np.random.default_rng(0), noise=False, te_term=te_term)
                want = M.theta_of(S) + d["te_s"] * M.delta(M.SCEN[S][1])
                for s in range(K):
                    m = d["site"] == s
                    got = d["y"][m & (d["g"] == 1)].mean(0) - d["y"][m & (d["g"] == 0)].mean(0)
                    worst = max(worst, float(np.abs(got - want[s]).max()))
    return dict(value=worst, ok=worst < 1e-10, criterion="사이트별 집단 차 = θ + TE_s·Δ, 최대 차 < 1e−10")


def _reps(c, fn, reps=M.REPS, seed=1, **gen):
    rng = np.random.default_rng(np.random.SeedSequence([M.SEED, 9000 + seed]))
    out = []
    for _ in range(reps):
        d = M.generate(c, rng, **gen)
        out.append(fn(d))
    return np.array(out)


def kb(se_scale=1.0):
    def f(d):
        e = A.m2(d["y"], d["g"], d["site"], K)
        e = dict(e, se=e["se"] * se_scale)
        return A.covered(e, d["theta"]).mean()
    cv = _reps(cell("D0", "S0"), f, seed=1)
    v = float(cv.mean())
    return dict(value=v, mcse=float(cv.std(ddof=1) / np.sqrt(len(cv))), ok=0.929 <= v <= 0.971, criterion="M2 포함률 0.929–0.971")


def kc(factor=M.LENGTHEN):
    want = M.TE_C * M.delta("e")                                   # expected value under lengthening
    er = _reps(cell("D0", "S2e"), lambda d: (A.m2(d["y"], d["g"], d["site"], K)["est"] - d["theta"]).mean(), seed=2, factor=factor)
    b, se = float(er.mean()), float(er.std(ddof=1) / np.sqrt(len(er)))
    return dict(value=b, mcse=se, expected=want, ok=abs(b - want) <= 3 * se, criterion="M2 편향 = TE_c·Δ_e ± 3 MC SE")


def kd(scen):
    def f(d):
        m3 = A.m3(d["y"], d["g"], d["site"], K, d["te"], M.TE_C)
        return A.rejected(m3["slope"]).mean()
    rj = _reps(cell("D1", scen), f, seed=3, eps=False, gamma_dev=False, vendor_mean=False)
    v, se = float(rj.mean()), float(rj.std(ddof=1) / np.sqrt(len(rj)))
    return v, se


def ke(eb=True):
    from design_sim.analyses import RA
    d = M.generate(cell("D1", "S0"), np.random.default_rng(5))
    tr = RA._combat(d["y"].T, d["site"], None, mean_only=False)
    est = copy.deepcopy(tr["estimates"])
    if not eb:
        est["gamma.star"], est["delta.star"] = est["gamma_hat"], est["delta_hat"]     # names of neuroCombat 0.2.12 (underscore)
    got = RA.combat_apply(d["y"].T, d["site"], est)
    dev = float(np.abs(got - tr["data"]).max())
    return dict(value=dev, ok=dev < 1e-10, criterion="combat_apply(모든 참가자로 학습) = neuroCombat 출력, 최대 차 < 1e−10")


def main():
    out = pathlib.Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    rows = []

    def add(name, clean, planted, fault):
        rows.append(dict(check=name, clean=clean, planted_fault=fault, planted=planted,
                         passed=bool(clean["ok"]), planted_failed=not bool(planted["ok"])))
        print(f"{name}: 깨끗한 코드 {'통과' if clean['ok'] else '실패'} ({clean['value']:.6g}) · 심은 결함 «{fault}» "
              f"{'실패(맞음)' if not planted['ok'] else '통과(틀림)'} ({planted['value']:.6g})", flush=True)

    add("K-0 Δ 상수", k0(), k0(1.10), "연장 배수 1.0984 → 1.10")
    add("K-a 잡음 없는 대수", ka(), ka(te_term=False), "생성기에서 TE 항을 뺌")
    add("K-b 영가설 보정(D0 · S0 · M2)", kb(), kb(0.5), "SE 를 절반으로")
    clean_c, short_c = kc(), kc(1.0 / M.LENGTHEN)
    short_c["sign_flipped"] = (np.sign(short_c["value"]) == -np.sign(clean_c["value"]))
    add("K-c 해석적 편향(D0 · S2e · M2)", clean_c, short_c, "T2 를 같은 비율로 줄임(부호가 바뀌어야)")
    for scen in ("S2b", "S2e"):
        v, se = kd(scen)
        clean = dict(value=v, mcse=se, expected=S75_POWER[scen], ok=abs(v - S75_POWER[scen]) <= 3 * se,
                     criterion=f"M3 group × TE 기각률 = {S75_POWER[scen]:.3f} ± 3 MC SE")
        v0, se0 = kd("S0")
        planted = dict(value=v0, mcse=se0, expected=S75_POWER[scen], ok=abs(v0 - S75_POWER[scen]) <= 3 * se0)
        add(f"K-d S-M9.7 검정력({scen})", clean, planted, "이완 차이를 0 으로(S0)")
    add("K-e ComBat 구현", ke(), ke(eb=False), "EB 축소를 끔(γ̂, δ̂)")
    ok = all(r["passed"] and r["planted_failed"] for r in rows) and bool(rows[3]["planted"]["sign_flipped"])
    rec = dict(when=time.strftime("%Y-%m-%d %H:%M:%S"), seconds=round(time.time() - t0, 1), all_ok=ok, checks=rows)
    (out / "checks.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1, default=float), encoding="utf-8")
    print("모두 맞음" if ok else "어긋남 있음 — 멈춤")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
