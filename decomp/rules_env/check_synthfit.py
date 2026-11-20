"""Checks rules_env.synthfit in envs/reference (AMICO), outside pytest: signal equality with CardProto and Proto and with synth.noiseless, a refit of 2,000 stored E5 rows to within ±1e-5, and a speed measurement with the Task 7 extrapolation."""
import pathlib
import sys
import time

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
for _cand in (HERE.parent, HERE.parent / "decomp"):
    if (_cand / "rules_env").is_dir():
        sys.path.insert(0, str(_cand))
        break
from rules_env import cards as K  # noqa: E402
from rules_env import synthfit as SF  # noqa: E402

C, synth = SF.C, SF.synth
WORK, E5DIR = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
T7_ROWS = 4_200_000                       # Task 7 estimate (S1(b) ≈ 1.0 M · S2 ≈ 2.6 M · K8 ≈ 0.5 M)
fails = []


def check(name, ok, msg=""):
    print(("통과 " if ok else "실패 ") + name + (" — " + msg if msg else ""), flush=True)
    if not ok:
        fails.append(name)


r = dict(np.load(E5DIR / "rows.npz"))             # Solved only once (npz re-reads the whole array on every access)
stored = dict(np.load(E5DIR / "f2_1.7.npz"))
idx = np.linspace(0, len(r["cell"]) - 1, 2000).astype(int)


def cell_of(i):
    tissue = C.TISSUE_E5[int(r["tissue"][i])]
    return dict(nu=tissue[0], odi=tissue[1], fiso=tissue[2], t2=C.t2_triplet(C.t2_of(int(r["t2set"][i]))),
                t1=C.t1_triplet(C.T1_SETS[int(r["t1set"][i])]), gen=str(r["gen"][i]), fdir=int(r["fdir"][i]), dpar=float(r["dpar"][i]))


def legacy_cell(i):
    return dict(tissue=int(r["tissue"][i]), dense=True, t2set=int(r["t2set"][i]), t1set=int(r["t1set"][i]), TE=float(r["TE"][i]),
                TR=float(r["TR"][i]), gen=str(r["gen"][i]), fdir=int(r["fdir"][i]), dpar=float(r["dpar"][i]))


card = K.card("E5:A")
pr_card, pr_old = SF.CardProto(card), synth.Proto("P2a")
a_ok = all(np.array_equal(synth.noiseless(pr_card, legacy_cell(i))[0], synth.noiseless(pr_old, legacy_cell(i))[0]) for i in idx[:5])
check("(a) CardProto = Proto(P2a)", a_ok)

gens_seen, b_max = set(), 0.0
for j, i in enumerate(idx[::40]):
    g = C.GENERATORS[j % 3]                                  # The stored chunk has only one generator, so the generator is swapped while running
    c, lc = dict(cell_of(i), gen=g), dict(legacy_cell(i), gen=g)
    new = SF.noiseless_explicit(pr_card, TE=card["TE"], TR=card["TR"], **c)
    old = synth.noiseless(pr_old, lc)
    b_max = max(b_max, float(np.max(np.abs(new[0] - old[0]))), abs(new[2] - old[2]), abs(new[3] - old[3]))
    gens_seen.add(g)
check("(b2) 병렬 · 캐시 신호 = 순차 신호", np.array_equal(SF.signals(card, [cell_of(i) for i in idx[:400]], n_workers=4)[0],
      np.asarray([SF.noiseless_explicit(pr_card, TE=card["TE"], TR=card["TR"], **cell_of(i))[0] for i in idx[:400]], np.float32)))
check("(b) noiseless_explicit = synth.noiseless", b_max == 0.0, f"최대 차 {b_max:.3g} · 생성기 {sorted(gens_seen)}")

t0 = time.time()
res = SF.fit_cells(card, [cell_of(i) for i in idx], WORK / "check_c", nthreads=32)
dt = time.time() - t0
dmax = max(float(np.max(np.abs(res[k] - stored[k][idx]))) for k in ("ndi", "fwf", "odi"))
check("(c) 저장된 E5 적합값 재현", dmax <= 1e-5, f"최대 차 {dmax:.3g} · 행 2,000 · {dt:.1f} 초")

t0 = time.time()
n = 20000
cells = [cell_of(i) for i in np.linspace(0, len(r["cell"]) - 1, n).astype(int)]
res = SF.fit_cells(card, cells, WORK / "check_d", nthreads=32)
dt = time.time() - t0
est_h = T7_ROWS / (n / dt) / 3600.0
check("(d) 속도 · 외삽", est_h <= 6.0, f"{n / dt:,.0f} 행/초(합성 + 적합) · Task 7 어림 {T7_ROWS:,} 행 → {est_h:.2f} 시간")

print(f"{5 - len(fails)} checks passed" if not fails else f"실패: {fails}")
sys.exit(1 if fails else 0)
