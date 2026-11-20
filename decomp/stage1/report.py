"""Writes the decomp stage1 report (summary.json, REPORT.md) from directions.json, e4.json and the direction store, per SPEC_v1 §17."""
from __future__ import annotations

import datetime
import itertools
import json
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze as A  # noqa: E402
import config as C  # noqa: E402
import lib  # noqa: E402

CATS = ("N-TE", "N-TR", "N-SNR", "N-P", "N-MxP")


def cat_vectors(S, pr, s, p):
    return {"N-TE": S["gTE"][pr, A.CENTRE, s, :p][:, None], "N-TR": S["gTR"][pr, :, s, :p].T,
            "N-SNR": S["gSNR"][pr, :, 0, 1, s, :p].T, "N-P": S["gP"][A.CENTRE, s, :p][:, None],
            "N-MxP": S["gMxP_A"][s, :, :p].T}


def direction_tables(store):
    """§17-3: angles between the representative vectors of the categories (pairwise · span) · ratio of singular values within a category — WM 12 states · median of P2a/P2b."""
    out = {}
    for key, S in store.items():
        for st in A.STAGES:
            p = lib.STAGE_P[st]
            pair, span, within = {}, {}, {}
            for pr in (1, 2):
                for s in C.WM:
                    W12 = lib.inv_sqrt(A.sig_for(S["Sig"], pr, A.CENTRE, s, 100.0, p))
                    g = cat_vectors(S, pr, s, p)
                    rep = {}
                    for c, v in g.items():
                        u, ratio = lib.representative(list(v.T), W12)
                        rep[c] = u
                        within.setdefault(c, []).append(ratio)
                    for a, b in itertools.combinations(CATS, 2):
                        ang = np.degrees(np.arccos(min(1.0, abs(float(rep[a] @ rep[b])))))
                        pair.setdefault(f"{a}~{b}", []).append(ang)
                    for c in CATS:
                        others = np.stack([rep[o] for o in CATS if o != c], 1)
                        span.setdefault(c, []).append(lib.span_angle(rep[c], others))
            out[f"{key}|{st}"] = dict(pair_deg={k: float(np.median(v)) for k, v in pair.items()},
                                     span_deg={k: float(np.median(v)) for k, v in span.items()},
                                     within_sv_ratio={k: float(np.median(v)) for k, v in within.items()})
    return out


def generator_rotation(store, F):
    """§17-5: rotation of u_TE between generator levels (whitening · Σ of T-in) — WM · center T2 · median of P2a/P2b."""
    out = {}
    for st in A.STAGES:
        p = lib.STAGE_P[st]
        for a, b in itertools.combinations(C.GENERATORS, 2):
            Sa, Sb = store[f"{F}|{a}"], store[f"{F}|{b}"]
            angs = [lib.wangle(Sa["gTE"][pr, A.CENTRE, s, :p], Sb["gTE"][pr, A.CENTRE, s, :p],
                               A.sig_for(Sa["Sig"], pr, A.CENTRE, s, 100.0, p)) for pr in (1, 2) for s in C.WM]
            out[f"{st}|{a}~{b}"] = float(np.median(angs))
    return out


def e4_summary(e4):
    rows = []
    for r in e4:
        F, gD, gU, pair, n, st, lam = r["cfg"]
        rows.append(dict(F=F, gen_data=gD, gen_U=gU, pair=pair, n=n, stage=st, lam=lam, **{k: v for k, v in r.items() if k != "cfg"}))
    return rows


def sg7(rows):
    """SG7: type I error rate of the null (λ = 0) with the same generator (gen_data = gen_U) ∈ [2.5 %, 7.5 %] · coverage of the Ψ interval ≥ 90 % — as a proportion of cells."""
    out = {}
    for F in A.FITTERS:
        for g in C.GENERATORS:
            for st in A.STAGES:
                cells = [r for r in rows if r["F"] == F and r["gen_data"] == g and r["gen_U"] == g and r["stage"] == st
                         and r["lam"] == 0 and not r.get("exhausted")]
                ok = [(2.5 <= 100 * r["rej"] <= 7.5) and r["cov"] >= 0.9 for r in cells]
                out[f"{F}|{g}|{st}"] = dict(n_cells=len(cells), frac_ok=float(np.mean(ok)) if ok else None,
                                            pass_=bool(ok) and all(ok),
                                            size=[round(100 * r["rej"], 1) for r in cells], cov=[round(r["cov"], 3) for r in cells])
    return out


def git(*a):
    return subprocess.run(["git", *a], cwd=A.REPO, capture_output=True, text=True).stdout.strip()


def main():
    d = json.load(open(A.OUT / "directions.json"))
    e4 = json.load(open(A.OUT / "e4.json"))
    raw = dict(np.load(C.store_path()))
    store = {}
    for k, v in raw.items():
        a, b = k.split("::")
        store.setdefault(a, {})[b] = v
    rows = e4_summary(e4)
    W = C.WORK / "logs"
    run_log = (W / "run.log").read_text().strip().splitlines()
    loads = [ln.split()[1:4] for ln in (W / "loadavg.log").read_text().splitlines() if ln.strip()]
    if A.FITTERS != ("F2",):
        # F1 version (exception to SPEC_v3 §0 · relay ZZZ-5) — the wording of the canonical (F2) summary is left as is, and only the F1 version is written separately (relay ZZZ-3 (b))
        f1 = [json.loads(ln) for ln in (C.WORK / "f1" / "ledger.jsonl").read_text(encoding="utf-8").splitlines() if ln.strip()]
        summary = dict(spec_commit="a862a41", code_commit=git("rev-parse", "--short=7", "HEAD"), fitters=list(A.FITTERS),
                       run=dict(fit_ledger="data/f1/ledger.jsonl", fit_first=f1[0]["UTC"], fit_last=f1[-1]["UTC"],
                                fit_chunks=sum(e["event"] == "덩어리 끝" for e in f1), rows={e: C.n_rows(e) for e in C.EXPERIMENTS},
                                gpu="F1 적합만 GPU(2026-10-21 사용자 지시 · §9-1 #48) — 분석은 CPU"),
                       impl=d["impl"], sg=d["results"], directions=direction_tables(store),
                       generator_rotation={F: generator_rotation(store, F) for F in A.FITTERS},
                       e4=rows, sg7=sg7(rows),
                       deviations=[dict(id="D-1 (닫힘 · F1 판)", what="F1(cuDIMOT) 판정 — SPEC_v3 §0 의 예외 «F1 재실행 — 같은 코드 · 같은 기준»(저자 결정 D3 · §9-1 #49)",
                                        original="F1 · F2 두 적합기의 판정(SPEC_v1 §7 · §11)", changed="F2 판(2b198a6)과 따로 낸 F1 판 — 갈래 DECOMP_TAG=_F1")],
                       generated=datetime.datetime.now(datetime.timezone.utc).isoformat())
        json.dump(summary, open(C.guard_write(A.OUT / "summary.json"), "w"), ensure_ascii=False, indent=1, default=float)
        print("summary.json 씀(F1 판)")
        return
    summary = dict(
        spec_commit="a862a41", code_commit="cd028ef", fitters=list(A.FITTERS),
        run=dict(log=run_log, loadavg_1min_range=[min(float(x[0]) for x in loads), max(float(x[0]) for x in loads)],
                 rows={e: C.n_rows(e) for e in C.EXPERIMENTS}, gpu="none (CPU 만 — 2026-10-20 사용자 지시)"),
        impl=d["impl"], sg=d["results"], directions=direction_tables(store),
        generator_rotation={F: generator_rotation(store, F) for F in A.FITTERS},
        e4=rows, sg7=sg7(rows),
        deviations=[dict(id="D-1", what="F1(cuDIMOT) 미실행 — 2026-10-20 사용자 지시(GPU 1 금지 · CPU 만). SPEC_v1 §7 · §11 은 F1 · F2 둘 다",
                         original="F1 · F2 두 적합기의 판정", changed="F2(AMICO) 판정만 — F1 은 GPU 허락 뒤(같은 코드 · run_stage1.sh 의 gpu 인자 · analyze 의 DECOMP_FITTERS=F1,F2)")],
        generated=datetime.datetime.now(datetime.timezone.utc).isoformat())
    json.dump(summary, open(C.guard_write(A.OUT / "summary.json"), "w"), ensure_ascii=False, indent=1, default=float)
    print("summary.json 씀 — REPORT.md 는 summary.json 에서 사람이 읽는 판으로 옮긴다")


if __name__ == "__main__":
    main()
