"""Exploratory d1 on c3: redraws main-cell posterior samples along the verdict path and computes bias and overconfidence per ROI unit, coil combination and k."""
from __future__ import annotations

import datetime
import hashlib
import json
import os
import pathlib
import sys
import time
import warnings

import numpy as np

REPO = pathlib.Path("<REPO>")
X14 = pathlib.Path("<WORK>/stage1_work/drafts/x14")
sys.path.insert(0, str(REPO / "decomp" / "remeasure"))
import judge_v4 as J  # noqa: E402

OUT = J.R.WORK / "explore" / "d1"
UNITS = OUT / "d1_c3_units"
V4 = J.R.WORK / "onharmony" / J.V4_DIR
REREAD = X14 / "REREAD_LOG.jsonl"
SHA = dict(attr_table="5184137854abc69710904889503eaf8dbb111ca9736e46af13b49067ea65b90d",      # Access record 2026-10-30T09:11:12Z
           snr_table="1714e6df9e2221df4abc59d58dc40f85c81d5e10c4406a4bcb2e36caddbd8cec",       # Access record 2026-10-30T09:11:11Z
           v4_c3="f6003589c01ad0bd693a6c809cb78ec6e69c70907be55703e22774d47f240a2a")           # Access record 2026-10-31T10:28:09Z
MAIN = (30.0, 50.0)
SENS = ((20.0, 25.0), (20.0, 100.0), (50.0, 25.0), (50.0, 100.0))
COMBOS4 = ((32, 32), (32, 64), (64, 32), (64, 64))
COIL_A, COIL_B = J.S.COND_NAMES.index("coil_A_64"), J.S.COND_NAMES.index("coil_B_64")
TAIL = 0.05
BOOT, BOOT_SEED = 2000, 20261102
COIL_ONLY_X = ("NOT1ACH", "NOT2ING-98", "OXF3TRI", "NOT4GEP", "OXF4GEP", "NOT2ING-89")       # §2-8 4b · six pairs
NOW = datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def _blocked(*a, **k):
    raise RuntimeError("⛔ D1 은 V4 결과 파일 · 열람 기록을 쓰지 않는다")


J._write_once = J._write_json = J._log = _blocked


def reread(path, what, expect=None) -> bytes:
    b = pathlib.Path(path).read_bytes()
    s = hashlib.sha256(b).hexdigest()
    with REREAD.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(dict(시각=NOW, 다시읽음=str(path), sha256=s, 바이트=len(b), 무엇=what, 스크립트=__file__), ensure_ascii=False) + "\n")
    if expect is not None and s != expect:
        raise RuntimeError(f"⛔ {path} sha256 {s[:16]}… ≠ 열람 기록 {expect[:16]}…")
    return b


def tables(purpose: str):
    at = json.loads(reread(V4 / "attr_table.json", f"D1 {purpose} — 속성 표(판정 때 읽은 것)", SHA["attr_table"]))
    v4 = json.loads(reread(V4 / "v4_c3.json", f"D1 {purpose} — 저장된 칸 요약 · 쌍 시드(재현 문)", SHA["v4_c3"]))
    return at, v4


def log_roi(at: dict, purpose: str) -> int:
    root = J.R.WORK / "onharmony"
    files = sorted({root / f"{sub}_{s}" / "roi.json" for d in at["pairs"].values() for sub, x in d.items() for s in x["sessions"]})
    for p in files:
        reread(p, f"D1 {purpose} — roi.json «주»(collect · 판정 때 읽은 것)")
    return len(files)


def identified() -> list[str]:
    ident = J.identified_cells(json.loads((J.RESULTS / "cards_c3.json").read_text(encoding="utf-8")))
    return sorted(f"{p}|{c}" for c, ps in ident.items() for p in ps)


def combos_of(pair: str) -> list[tuple[int, int]]:
    a, b = pair.split("~")
    return [(ca, cb) for ca in J.coil_options(a) for cb in J.coil_options(b)]


def env_check():
    assert os.environ.get("CUDA_VISIBLE_DEVICES", None) == "", "CUDA_VISIBLE_DEVICES= 로 띄운다(CPU 만)"
    assert os.environ.get("S4_THREADS") == "16", "S4_THREADS=16(원 실행과 같게)"


class Rec:
    """Wraps draw_fn — returns the received (T, ok, acc) unchanged (0 change in computation or random-number consumption) and keeps only a summary of the samples. Order of calls = (outer over combination j, inner over k) of cover_by_k → j, kk = divmod(number of calls, number of k)."""

    def __init__(self, draw_fn, n: int, m: int, nk: int, truth=None):
        self.f, self.nk, self.c = draw_fn, nk, 0
        self.truth = np.zeros((n, 3)) if truth is None else np.asarray(truth, float)
        self.pit = np.full((n, m, nk, 3), np.nan, np.float64)
        self.mean = np.full((n, m, nk, 3), np.nan, np.float64)
        self.sd = np.full((n, m, nk, 3), np.nan, np.float64)
        self.cov = np.zeros((n, m, nk, 3), bool)
        self.ok = np.zeros((n, m, nk), bool)

    def __call__(self, X, cond, n_keep, gs):
        T, ok, acc = self.f(X, cond, n_keep, gs)
        j, kk = divmod(self.c, self.nk)
        self.c += 1
        A = np.asarray(T[:, :, 0:3], dtype=np.float64)
        fin = np.isfinite(A)
        nf = fin.sum(1)
        below = ((A < self.truth[:, None, :]) & fin).sum(1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            self.pit[:, j, kk] = np.where(nf > 0, below / np.maximum(nf, 1), np.nan)
            self.mean[:, j, kk] = np.nanmean(A, 1)
            self.sd[:, j, kk] = np.nanstd(A, 1, ddof=1)
            lo, hi = np.nanquantile(T[:, :, 0:3], 0.05, axis=1), np.nanquantile(T[:, :, 0:3], 0.95, axis=1)
        self.cov[:, j, kk] = (lo <= self.truth) & (self.truth <= hi)
        self.ok[:, j, kk] = np.asarray(ok, bool)
        return T, ok, acc


def compare(u: dict, c4: dict, ci: int) -> dict:
    """Reproduction stop condition — with the verdict «main» variant (subjects excluded from the main cell left out), cell_verdict · coil_counts are recomputed with the same code and compared with the stored values."""
    if c4.get("verdict") == "자료 없음":
        n = int(len(u["cov"]))
        return dict(no_data=True, match=n == 0, mismatch_units=n, mismatch_coil=0)
    us = J.select(u, ~u["main_excluded"])
    r = J.cell_verdict(us["cov"][:, :, ci], us["ok"], us["flag"], us["subjects"], None, np.random.default_rng(0), None, n_boot=0)
    cc = J.coil_counts(us)
    n1, n0 = r["n_units_used"], c4["n_units_used"]
    k1 = [round(x * n1) if x is not None else 0 for x in r["coverage_by_k"]]
    k0 = [round(x * n0) if x is not None else 0 for x in c4["coverage_by_k"]]
    mis = (sum(abs(p - q) for p, q in zip(k1, k0)) + abs(n1 - n0) + abs(r["n_excluded_misspec"] - c4["n_excluded_misspec"])
           + abs(r["n_units_dropped_k"] - c4["n_units_dropped_k"]))
    c0 = c4.get("by_coil_combo", {})
    mc = sum(abs(cc.get(k, {}).get(f, 0) - c0.get(k, {}).get(f, 0)) for k in set(cc) | set(c0) for f in ("n_misspec", "n_k_unfilled"))
    exact = (r["coverage_by_k"] == c4["coverage_by_k"] and n1 == n0 and r["n_excluded_misspec"] == c4["n_excluded_misspec"]
             and r["n_units_dropped_k"] == c4["n_units_dropped_k"] and cc == c0)
    return dict(no_data=False, match=bool(exact), mismatch_units=int(mis), mismatch_coil=int(mc))


# ── real ───────────────────────────────────────────────────────────────────────
def mode_real():
    env_check()
    for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        assert v not in os.environ, f"{v} 는 정하지 않는다(원 실행기와 같게)"
    UNITS.mkdir(parents=True, exist_ok=True)
    at, v4 = tables("real")
    n_roi = log_roi(at, "real")
    draw_fn, logq_fn, K = J.real_fns("c3")
    miss = json.loads((J.RESULTS / "miss.json").read_text(encoding="utf-8"))
    cells = identified()
    assert len(cells) == 29 and set(cells) == set(v4["main"]["cells"]), "식별 칸 29 가 저장값과 다르다"
    rows = J.collect(J.R.WORK, at["pairs"])
    gate = dict(rule="칸 29(자료 있는 24 + 없는 5) 모두 k 아홉 덮음 · n_units_used · n_excluded_misspec · n_units_dropped_k · by_coil_combo 가 "
                     "v4_c3.json 과 정확히 같음 · 한 쌍이라도 어긋나면 그 자리에서 멈춤(계획서 §2-2 · §16)",
                started=NOW, roi_files_logged=n_roi, pairs={}, cells={}, passed=None)
    gp = OUT / "d1_gate.json"
    t0 = time.time()
    for pi, pair in enumerate(J.card_pairs()):
        seed = v4["seeds_torch"][f"SNR30|Neff50|{pair}"]
        assert seed == 6000 + pi, f"시드가 판정 규칙과 다르다 {pair}"
        m, nk = len(combos_of(pair)), len(K)
        rec = Rec(draw_fn, len(rows[pair]), m, nk)
        ts = time.time()
        J.seed_torch(seed)
        u = J.unit_results(pair, rows[pair], "c3", rec, logq_fn, K, miss["layers"][f"{J.acq_type(pair)}|L0"]["threshold"], *MAIN)
        if len(rows[pair]):
            assert rec.c == m * nk
            assert (rec.cov.all(1) == u["cov"]).all() and (rec.ok.all(1) == u["ok"]).all(), "감싼 요약과 덮음이 다르다"
            np.savez(UNITS / f"{pair}.npz", pit=rec.pit, mean=rec.mean, sd=rec.sd, okc=rec.ok, cov=u["cov"], ok=u["ok"], flag=u["flag"],
                     subjects=u["subjects"], rois=np.array([r["unit"] for r in rows[pair]]), main_excluded=u["main_excluded"], noise=u["noise"],
                     combos=np.array(u["combos"]), flag_by_combo=u["flag_by_combo"], ok_by_combo=u["ok_by_combo"], seed=seed,
                     covc=rec.cov, X3=np.array([np.r_[r["x"][0:3], r["x"][7:10]] for r in rows[pair]]))   # §16-2: coverage and position per combination
        mis = 0
        for key in [c for c in cells if c.split("|")[0] == pair]:
            gate["cells"][key] = compare(u, v4["main"]["cells"][key], J.COMPS.index(key.split("|")[1]))
            mis += gate["cells"][key]["mismatch_units"] + gate["cells"][key]["mismatch_coil"] + (0 if gate["cells"][key]["match"] else 1)
        gate["pairs"][pair] = dict(n_units=int(len(rows[pair])), seconds=round(time.time() - ts, 1), seed=seed)
        (OUT / "STATUS_REAL").write_text(f"{datetime.datetime.now():%F %T} 쌍 {pi + 1}/{len(J.card_pairs())} {pair} · 어긋남 {mis}\n", encoding="utf-8")
        if mis:
            gate.update(passed=False, stopped_at=pair, seconds=round(time.time() - t0, 1),
                        mismatch_units_total=sum(c["mismatch_units"] for c in gate["cells"].values()),
                        mismatch_coil_total=sum(c["mismatch_coil"] for c in gate["cells"].values()),
                        cells_not_checked=[c for c in cells if c not in gate["cells"]])
            gp.write_text(json.dumps(gate, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
            print("재현 문 실패", pair, flush=True)
            sys.exit(2)
    gate.update(passed=all(c["match"] for c in gate["cells"].values()) and len(gate["cells"]) == 29, seconds=round(time.time() - t0, 1),
                mismatch_units_total=0, mismatch_coil_total=0)
    gp.write_text(json.dumps(gate, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    print("재현 문", "통과" if gate["passed"] else "실패", flush=True)
    sys.exit(0 if gate["passed"] else 2)


# ── ref ────────────────────────────────────────────────────────────────────────
def needed_types() -> list[str]:
    return [t for t in J.S.ACQ_TYPES if t in {J.acq_type(p) for p in J.card_pairs() if p not in J.NO_DATA}]


def mode_ref():
    env_check()
    OUT.mkdir(parents=True, exist_ok=True)
    draw_fn, logq_fn, K = J.real_fns("c3")
    V = J._validate()
    sim, _ = V.load_set("test")
    miss = json.loads((J.RESULTS / "miss.json").read_text(encoding="utf-8"))
    for ti, t in enumerate(J.S.ACQ_TYPES):
        f = OUT / f"d1_ref_{t}.npz"
        mk = OUT / f"d1_ref_{t}.inprogress"
        if t not in needed_types() or f.exists():
            continue
        if mk.exists():                                                       # §16-2: another process is computing — wait until it finishes
            while not f.exists():
                time.sleep(60)
            continue
        mk.write_text(str(os.getpid()), encoding="utf-8")
        ts = time.time()
        idx = np.flatnonzero(sim["layer"] == f"{t}|L0")
        X, cond, X3 = V.xobs(sim, "c3")[idx], sim["cond"][idx], sim["x_c3"][idx]
        th = sim["theta"][idx, 0:3].astype(np.float64)
        thr = miss["layers"][f"{t}|L0"]["threshold"]
        rec = Rec(draw_fn, len(idx), len(COMBOS4), len(K), truth=th)
        score, seeds = np.empty((len(idx), len(COMBOS4))), {}
        for cj, (ca, cb) in enumerate(COMBOS4):
            cc = cond.copy()
            cc[:, COIL_A], cc[:, COIL_B] = float(ca == 64), float(cb == 64)
            seeds[f"{ca}|{cb}"] = J.seed_torch(7000 + 10 * ti + cj)
            rec.c = cj * len(K)
            for gs in K:
                rec(X, cc, J.N_KEEP, gs)
            score[:, cj] = J.miss_flags(logq_fn, X3, [cc], K, thr)[1]
            (OUT / "STATUS_REF").write_text(f"{datetime.datetime.now():%F %T} {t} 조합 {ca}|{cb} 끝\n", encoding="utf-8")
        true = np.c_[np.where(cond[:, COIL_A] > 0.5, 64, 32), np.where(cond[:, COIL_B] > 0.5, 64, 32)]
        np.savez(f, pit=rec.pit, mean=rec.mean, sd=rec.sd, ok=rec.ok, cov=rec.cov, score=score, threshold=thr, true_combo=true,
                 theta=th, idx=idx, seeds=json.dumps(seeds), seconds=round(time.time() - ts, 1))
        mk.unlink(missing_ok=True)
    print("기준 분포 끝", flush=True)


# ── sens ───────────────────────────────────────────────────────────────────────
def mode_sens():
    env_check()
    at, v4 = tables("sens")
    log_roi(at, "sens")
    _, logq_fn, K = J.real_fns("c3")
    miss = json.loads((J.RESULTS / "miss.json").read_text(encoding="utf-8"))
    cells = identified()
    rows = J.collect(J.R.WORK, at["pairs"])
    out = dict(rule="§9 — 민감도 칸마다 동결 q 밀도로 오지정 점수만(같은 문턱 · 코일 조합 규칙 같음 · main 변형). k 미충족 수는 저장된 통과 표시가 None 이고 "
                    "피험자 ≥ 3 인 칸에서만 «단위 수 − 오지정 수» · 그 밖은 «알 수 없음»(계획서 §2-9)", main_check={}, cells={})
    for snr, neff in (MAIN,) + SENS:
        lab = f"SNR{int(snr)}|Neff{int(neff)}"
        for pair in J.card_pairs():
            keys = [c for c in cells if c.split("|")[0] == pair]
            r = rows[pair]
            if not r:
                for key in keys:
                    out["cells"].setdefault(key, {})[lab] = dict(no_data=True)
                continue
            X7 = np.array([x["x"] for x in r])
            X3 = np.c_[X7[:, 0:3], X7[:, 7:10]]
            te, tr = np.array([x["te"] for x in r], float), np.array([x["tr"] for x in r], float)
            conds = [J.cond_rows(pair, te, tr, cc, snr, neff) for cc in combos_of(pair)]
            flag = J.miss_flags(logq_fn, X3, conds, K, miss["layers"][f"{J.acq_type(pair)}|L0"]["threshold"])[0]
            keep = ~np.array([x["main_excluded"] for x in r])
            subs = np.array([x["sub"] for x in r])[keep]
            n_units, n_mis, n_sub = int(keep.sum()), int(flag[keep].sum()), int(len(np.unique(subs)))
            for key in keys:
                c4 = v4["main"]["cells"][key]
                if (snr, neff) == MAIN:
                    out["main_check"][key] = dict(n_excluded_misspec=n_mis, stored=c4["n_excluded_misspec"], same=n_mis == c4["n_excluded_misspec"])
                    continue
                stored = c4.get("sensitivity_passed", {}).get(lab, "저장값 없음")
                if stored is None and n_sub >= J.R.PAIR_MIN_SUBJECTS:
                    ku = n_units - n_mis
                elif stored is None:
                    ku = "알 수 없음(피험자 < 3 — 통과 표시가 쓴 단위 수를 알려 주지 않음)"
                else:
                    ku = "알 수 없음(표본 필요)"
                out["cells"].setdefault(key, {})[lab] = dict(n_units=n_units, n_subjects=n_sub, n_excluded_misspec=n_mis, stored_passed=stored, n_k_unfilled=ku)
    out["main_check_all_same"] = all(v["same"] for v in out["main_check"].values())
    (OUT / "d1_sens_flags.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    print("민감도 오지정 표시 끝 · 주 칸 오지정 수 재현", out["main_check_all_same"], flush=True)


# ── summary ────────────────────────────────────────────────────────────────────
def pick_rule(pit, mean, sd, truth):
    """If there are several coil combinations, picks the combination whose PIT is farthest from 0.5 (the earlier combination if equal) and also computes z from that combination → two (n, k, 3) · the chosen combination."""
    d = np.where(np.isfinite(pit), np.abs(pit - 0.5), -1.0)
    j = np.argmax(d, axis=1)
    P = np.take_along_axis(pit, j[:, None], 1)[:, 0]
    with np.errstate(invalid="ignore", divide="ignore"):
        Zall = (mean - np.asarray(truth, float)[:, None, None, :]) / sd
    return P, np.take_along_axis(Zall, j[:, None], 1)[:, 0], j


def kstats(P, Z, use) -> dict:
    """One component P · Z (n, k) · use (n,) → for each k the PIT mean · proportion at both ends · z mean · z SD · PIT SD (each (k,))."""
    Pu, Zu = P[use], Z[use]
    if len(Pu) == 0:
        nan = np.full(P.shape[1], np.nan)
        return dict(pit_mean=nan, tail=nan, z_mean=nan, z_sd=nan, pit_sd=nan)
    return dict(pit_mean=Pu.mean(0), tail=((Pu < TAIL) | (Pu > 1 - TAIL)).mean(0), z_mean=Zu.mean(0),
                z_sd=Zu.std(0, ddof=1) if len(Zu) > 1 else np.full(P.shape[1], np.nan), pit_sd=Pu.std(0, ddof=1) if len(Pu) > 1 else np.full(P.shape[1], np.nan))


def decomp(V, rois, use):
    """(a) Between subjects within the same ROI = pooled SD of what remains after subtracting the ROI mean (degrees of freedom = number of units − number of ROIs) · (b) SD of the ROI means · (b) comparison value = (a) ÷ √(mean number of units per ROI) → one value for each (k,)."""
    Vu, ru = V[use], rois[use]
    groups = [Vu[ru == r] for r in np.unique(ru)]
    res2 = sum(((g - g.mean(0)) ** 2).sum(0) for g in groups)
    df = sum(len(g) - 1 for g in groups)
    nan = np.full(V.shape[1], np.nan)
    a = np.sqrt(res2 / df) if df > 0 else nan
    b = np.array([g.mean(0) for g in groups]).std(0, ddof=1) if len(groups) > 1 else nan
    return a, b, a / np.sqrt(np.mean([len(g) for g in groups])) if groups else nan


def med(v) -> dict | None:
    v = np.asarray(v, float)
    if not np.isfinite(v).any():
        return None
    return dict(median=float(np.nanmedian(v)), min=float(np.nanmin(v)), max=float(np.nanmax(v)))


def spearman(x, y):
    from scipy.stats import spearmanr
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4 or np.ptp(x[m]) == 0 or np.ptp(y[m]) == 0:
        return dict(rho=None, n=int(m.sum()))
    return dict(rho=float(spearmanr(x[m], y[m]).statistic), n=int(m.sum()))


def mode_summary():
    gate = json.loads((OUT / "d1_gate.json").read_text(encoding="utf-8"))
    if not gate.get("passed"):
        raise SystemExit("⛔ 재현 문을 통과하지 않았다 — D1 값을 내지 않는다(계획서 §2-2)")
    at, v4 = tables("summary")
    snr = json.loads(reread(V4 / "snr_table.json", "D1 summary — 세션 SNR(§2-7 예측 비 · §2-8 ΔSNR)", SHA["snr_table"]))
    S_, (lo_t, hi_t) = snr["sessions"], snr["train_range"]
    cards = {c["pair"]: c for c in json.loads((REPO / "decomp/results/stage3/cards/cards.json").read_text(encoding="utf-8")) if c.get("set") == "ON-Harmony"}
    miss = json.loads((J.RESULTS / "miss.json").read_text(encoding="utf-8"))
    cells = [c for c in identified() if v4["main"]["cells"][c].get("verdict") != "자료 없음"]
    rng = np.random.default_rng(BOOT_SEED)
    res = dict(made=NOW, script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(), gate=dict(passed=True, seconds=gate.get("seconds")),
               note="탐색 분석 · 판정 아님 · 판단기의 단위별 사후 평균 · 물리 단위 구간은 내지 않음", tail=TAIL, cells={})
    ref_cache = {}
    for key in cells:
        pair, comp = key.split("|")
        ci = J.COMPS.index(comp)
        d = dict(np.load(UNITS / f"{pair}.npz", allow_pickle=False))
        n = len(d["flag"])
        P, Z, _ = pick_rule(d["pit"], d["mean"], d["sd"], np.zeros((n, 3)))
        P, Z = P[:, :, ci], Z[:, :, ci]
        okall = d["ok"].all(1)
        use = ~d["main_excluded"] & ~d["flag"] & okall
        use_in = ~d["main_excluded"] & okall                                         # Sensitivity: misspecified ROI units included
        subs, rois = d["subjects"], d["rois"]
        st, st_in = kstats(P, Z, use), kstats(P, Z, use_in)
        aZ, bZ, eZ = decomp(Z, rois, use)
        aP, bP, eP = decomp(P, rois, use)
        cov_k = d["cov"][use][:, :, ci].mean(0) if use.any() else np.full(9, np.nan)
        # Criterion (§2-4) — same acquisition type · true combination inside the candidates · same misspecification · sample rule
        t = J.acq_type(pair)
        if t not in ref_cache:
            ref_cache[t] = dict(np.load(OUT / f"d1_ref_{t}.npz", allow_pickle=False))
        rf = ref_cache[t]
        cand = combos_of(pair)
        cj = [COMBOS4.index(c) for c in cand]
        incand = np.array([tuple(x) in cand for x in rf["true_combo"].tolist()])
        okc = rf["ok"][:, cj].all(axis=(1, 2))
        flg = (rf["score"][:, cj] > float(rf["threshold"])).all(1)
        ruse = incand & okc & ~flg
        Pr, Zr, _ = pick_rule(rf["pit"][:, cj], rf["mean"][:, cj], rf["sd"][:, cj], rf["theta"])
        jt = np.array([COMBOS4.index(tuple(x)) for x in rf["true_combo"].tolist()])
        Pt = rf["pit"][np.arange(len(jt)), jt][:, :, ci]
        with np.errstate(invalid="ignore", divide="ignore"):
            Zt = ((rf["mean"][np.arange(len(jt)), jt] - rf["theta"][:, None, :]) / rf["sd"][np.arange(len(jt)), jt])[:, :, ci]
        sr, stt = kstats(Pr[:, :, ci], Zr[:, :, ci], ruse), kstats(Pt, Zt, ruse)
        # Noise prediction (§2-7) · ΔSNR (§2-8)
        sess = {s: at["pairs"][pair][s]["sessions"] for s in np.unique(subs)}
        snrA = np.array([S_[f"{s}_{sess[s][0]}"]["snr"] for s in subs], float)
        snrB = np.array([S_[f"{s}_{sess[s][1]}"]["snr"] for s in subs], float)
        sp = np.sqrt(2.0 / (1.0 / snrA ** 2 + 1.0 / snrB ** 2))
        pred = float(np.median(30.0 / sp[use])) if use.any() else None
        with np.errstate(invalid="ignore", divide="ignore"):
            obs_ratio_k = aZ / sr["z_sd"]
        outr = (snrA < lo_t) | (snrA > hi_t) | (snrB < lo_t) | (snrB > hi_t)
        st_out, st_inr = kstats(P, Z, use & outr), kstats(P, Z, use & ~outr)
        aZo = decomp(Z, rois, use & outr)[0]
        # Within-subject ΔSNR relation (descriptive) — mean z of the ROI units used per subject (for each k) → median over k
        us_ = np.unique(subs[use])
        zsub = np.array([np.median(Z[use & (subs == s)].mean(0)) for s in us_]) if len(us_) else np.array([])
        psub = np.array([np.median(P[use & (subs == s)].mean(0)) for s in us_]) if len(us_) else np.array([])
        dsub = np.array([float((snrB - snrA)[subs == s][0]) for s in us_]) if len(us_) else np.array([])
        # Subject-batch bootstrap (§2-3 · descriptive)
        bs = {k: [] for k in ("pit_mean", "tail", "z_mean", "z_sd")}
        su = np.unique(subs[use])
        for _ in range(BOOT if len(su) else 0):
            pick = rng.choice(su, len(su), replace=True)
            ix = np.concatenate([np.flatnonzero(use & (subs == s)) for s in pick])
            b_st = kstats(P[ix], Z[ix], np.ones(len(ix), bool))
            for k in bs:
                bs[k].append(float(np.nanmedian(b_st[k])))
        cd = cards.get(pair, {})
        a_, b_ = pair.split("~")
        # §16-2 — coverage per combination (same ROI unit used · one combination at a time) · position (mean of the two sessions of y) and |z| · rank correlation within the cells at both ends
        combo_cov = {}
        if "covc" in d:
            for j, cn in enumerate(d["combos"].tolist()):
                combo_cov[cn] = med(d["covc"][use][:, j, :, ci].mean(0)) if use.any() else None
        loc = d["X3"][:, [1, 4]].mean(1) if "X3" in d else None
        absz, tailu = np.median(np.abs(Z), 1), ((P < TAIL) | (P > 1 - TAIL)).mean(1)
        loc_rel = dict(rho_absz=spearman(loc[use], absz[use]), rho_tail=spearman(loc[use], tailu[use])) if loc is not None else None
        res["cells"][key] = dict(combo_coverage=combo_cov, location_y_relation=loc_rel,
            acq_type=t, n_units_used=int(use.sum()), n_subjects=int(len(su)), n_rois=int(len(np.unique(rois[use]))), coverage=med(cov_k),
            pit_mean=med(st["pit_mean"]), tail=med(st["tail"]), z_mean=med(st["z_mean"]), z_sd=med(st["z_sd"]), pit_sd=med(st["pit_sd"]),
            boot95={k: ([float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))] if v else None) for k, v in bs.items()},
            with_misspec_units=dict(n=int(use_in.sum()), pit_mean=med(st_in["pit_mean"]), tail=med(st_in["tail"]), z_mean=med(st_in["z_mean"]), z_sd=med(st_in["z_sd"])),
            a_z=med(aZ), b_z=med(bZ), b_z_expected_if_no_roi_diff=med(eZ), a_pit=med(aP), b_pit=med(bP), b_pit_expected_if_no_roi_diff=med(eP),
            ref=dict(n_cases=int(ruse.sum()), candidates=[f"{x}|{y}" for x, y in cand],
                     rule=dict(pit_mean=med(sr["pit_mean"]), tail=med(sr["tail"]), z_mean=med(sr["z_mean"]), z_sd=med(sr["z_sd"]), pit_sd=med(sr["pit_sd"])),
                     true_combo_only=dict(pit_mean=med(stt["pit_mean"]), tail=med(stt["tail"]), z_mean=med(stt["z_mean"]), z_sd=med(stt["z_sd"])),
                     tail_added_by_rule=med(sr["tail"] - stt["tail"])),
            noise=dict(obs_ratio_a_over_ref=med(obs_ratio_k), pred_ratio_30_over_snr_pair=pred,
                       out_of_train=dict(n=int((use & outr).sum()), pit_mean=med(st_out["pit_mean"]), tail=med(st_out["tail"]), z_mean=med(st_out["z_mean"]),
                                         z_sd=med(st_out["z_sd"]), a_z=med(aZo)),
                       in_train=dict(n=int((use & ~outr).sum()), pit_mean=med(st_inr["pit_mean"]), tail=med(st_inr["tail"]), z_mean=med(st_inr["z_mean"]), z_sd=med(st_inr["z_sd"]))),
            attrs=dict(dTE=cd.get("dTE"), abs_dTE=None if cd.get("dTE") is None else abs(cd["dTE"]), dTR_ms=cd.get("dTR_ms"), design_differs=cd.get("design_differs"),
                       coil=f"{J.P.ONH[a_][3] or '미상'}|{J.P.ONH[b_][3] or '미상'}", dSNR_median=float(np.median((snrB - snrA)[use])) if use.any() else None,
                       snr_ratio_median=float(np.median((snrB / snrA)[use])) if use.any() else None),
            within_cell_dSNR=dict(n_subjects=int(len(us_)), rho_z=spearman(dsub, zsub), rho_pit=spearman(dsub, psub)))
    C = res["cells"]
    res["no_data_cells"] = [c for c in identified() if c not in C]
    att = lambda f: [C[k]["attrs"][f] for k in C]  # noqa: E731
    rel = {}
    for f in ("pit_mean", "z_mean", "a_z"):
        y = [np.nan if C[k][f] is None else C[k][f]["median"] for k in C]
        rel[f] = dict(abs_dTE=spearman(att("abs_dTE"), y), dTR_ms=spearman(att("dTR_ms"), y), dSNR_median=spearman(att("dSNR_median"), y),
                      by_design={str(g): med([v for v, dd in zip(y, att("design_differs")) if dd == g]) for g in (False, True)},
                      by_coil={g: med([v for v, cc in zip(y, att("coil")) if cc == g]) for g in sorted(set(att("coil")))})
    res["relations_descriptive"] = rel
    # §16-2 — signed ΔTE (B − A) and bias · additive device model (cell z mean ≈ intercept + d_B − d_A · card unit 8 · leave-one-cell-out prediction R²)
    res["signed_dTE_relation"] = {f: spearman(att("dTE"), [np.nan if C[k][f] is None else C[k][f]["median"] for k in C]) for f in ("z_mean", "pit_mean")}
    odi = [k for k in C if k.endswith("|ODI") and C[k]["z_mean"] is not None]
    units8 = sorted({u for k in odi for u in k.split("|")[0].split("~")})
    Xd = np.zeros((len(odi), len(units8) + 1))
    for i, k in enumerate(odi):
        a_, b_ = k.split("|")[0].split("~")
        Xd[i, 0], Xd[i, 1 + units8.index(b_)], Xd[i, 1 + units8.index(a_)] = 1.0, 1.0, -1.0
    yd = np.array([C[k]["z_mean"]["median"] for k in odi])
    fit = lambda X, y: np.linalg.lstsq(X, y, rcond=None)[0]  # noqa: E731
    pred_in = Xd @ fit(Xd, yd)
    pred_loo = np.array([Xd[i] @ fit(np.delete(Xd, i, 0), np.delete(yd, i)) for i in range(len(yd))])
    ss = float(((yd - yd.mean()) ** 2).sum())
    res["additive_device_model_z_mean_ODI"] = dict(n_cells=len(odi), units=units8, r2_in=float(1 - ((yd - pred_in) ** 2).sum() / ss) if ss > 0 else None,
                                                  r2_loo=float(1 - ((yd - pred_loo) ** 2).sum() / ss) if ss > 0 else None,
                                                  effects=dict(zip(["intercept"] + units8, fit(Xd, yd).tolist())))
    res["coil_only_contrast"] = {x: {side: (f"{side}~{x}|ODI", C.get(f"{side}~{x}|ODI") and {f: C[f"{side}~{x}|ODI"][f] for f in ("n_units_used", "pit_mean", "tail", "z_mean", "z_sd", "a_z", "b_z")})
                                     for side in ("OXF1PRI", "OXF2PRI")} for x in COIL_ONLY_X}
    (OUT / "d1_c3.json").write_text(json.dumps(res, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    write_md(res)
    figures(res)
    print("D1 요약 끝", flush=True)


def write_md(res: dict):
    f2 = lambda x: "—" if x is None else f"{x['median']:.2f} [{x['min']:.2f}, {x['max']:.2f}]"  # noqa: E731
    f1 = lambda x: "—" if x is None else f"{x['median']:.2f}"  # noqa: E731
    C = res["cells"]
    L = [f"# D1 요약(탐색 분석 · 판정 아님 · {res['made']} · 스크립트 sha256 {res['script_sha256'][:16]}…)", "",
         "- 재현 문 통과(칸 29 · 허용 오차 없음) · 값은 k 아홉의 중앙값 [최소, 최대] · PIT 양 끝 = PIT < 0.05 또는 > 0.95.",
         "- 판단기의 단위별 사후 평균 · 물리 단위 구간은 내지 않는다.",
         f"- 자료 없는 칸(계산 없음): {' · '.join(c.replace('|', ' · ') for c in res['no_data_cells'])}", "",
         "## 1. 칸 요약(주 분석 · 오지정 뺌)", "",
         "| 칸 | 단위 · 피험자 | 덮음 | PIT 평균 | 양 끝 | z 평균 | z SD | (a) z | (b) z (차이 없을 때 기대) |", "|---|---|---|---|---|---|---|---|---|"]
    for k, v in C.items():
        L.append(f"| {k.replace('|', ' · ')} | {v['n_units_used']} · {v['n_subjects']} | {f1(v['coverage'])} | {f2(v['pit_mean'])} | {f2(v['tail'])} | {f2(v['z_mean'])} | "
                 f"{f2(v['z_sd'])} | {f1(v['a_z'])} | {f1(v['b_z'])} ({f1(v['b_z_expected_if_no_roi_diff'])}) |")
    L += ["", "## 2. 기준 분포(시험 자료 · 같은 획득 유형 · 참 조합이 후보 안)", "",
          "| 칸 | 경우 | 후보 조합 | PIT 평균 규칙 · 참 조합 | 양 끝 규칙 · 참 조합(규칙이 더한 몫) | z 평균 규칙 | z SD 규칙 |", "|---|---|---|---|---|---|---|"]
    for k, v in C.items():
        r = v["ref"]
        L.append(f"| {k.replace('|', ' · ')} | {r['n_cases']} | {' '.join(r['candidates'])} | {f1(r['rule']['pit_mean'])} · {f1(r['true_combo_only']['pit_mean'])} | "
                 f"{f1(r['rule']['tail'])} · {f1(r['true_combo_only']['tail'])} ({f1(r['tail_added_by_rule'])}) | {f1(r['rule']['z_mean'])} | {f1(r['rule']['z_sd'])} |")
    L += ["", "## 3. 잡음 가설의 예측(계산 전에 적은 식 · 판정 아님)", "",
          "| 칸 | 관측 비 (a) ÷ 기준 z SD | 예측 비 30 ÷ SNR_pair 중앙 | 범위 밖 단위: n · PIT 평균 · 양 끝 · z SD | 범위 안 단위: n · PIT 평균 · 양 끝 · z SD |", "|---|---|---|---|---|"]
    for k, v in C.items():
        o, i = v["noise"]["out_of_train"], v["noise"]["in_train"]
        pr = v["noise"]["pred_ratio_30_over_snr_pair"]
        L.append(f"| {k.replace('|', ' · ')} | {f1(v['noise']['obs_ratio_a_over_ref'])} | {'—' if pr is None else f'{pr:.2f}'} | "
                 f"{o['n']} · {f1(o['pit_mean'])} · {f1(o['tail'])} · {f1(o['z_sd'])} | {i['n']} · {f1(i['pit_mean'])} · {f1(i['tail'])} · {f1(i['z_sd'])} |")
    L += ["", "## 4. 관계 기술(순위 상관 · p 값 · 문턱 없음 · 칸끼리 독립 아님)", ""]
    for f, r in res["relations_descriptive"].items():
        rr = lambda x: "—" if x["rho"] is None else f"{x['rho']:+.2f} (n {x['n']})"  # noqa: E731
        L.append(f"- {f}: |ΔTE| {rr(r['abs_dTE'])} · ΔTR {rr(r['dTR_ms'])} · ΔSNR 중앙 {rr(r['dSNR_median'])} · 설계 같음/다름 {f1(r['by_design']['False'])} / {f1(r['by_design']['True'])} · "
                 "코일 " + " · ".join(f"{g} {f1(x)}" for g, x in r["by_coil"].items()))
    L += ["", "### 코일만 다른 대조(OXF1PRI~X · OXF2PRI~X · ODI)", "", "| X | 쪽 | 단위 | PIT 평균 | 양 끝 | z 평균 | z SD | (a) z | (b) z |", "|---|---|---|---|---|---|---|---|---|"]
    for x, sides in res["coil_only_contrast"].items():
        for side, (name, v) in sides.items():
            if not v:
                L.append(f"| {x} | {side} | 자료 없음 | | | | | | |")
                continue
            L.append(f"| {x} | {side} | {v['n_units_used']} | {f1(v['pit_mean'])} | {f1(v['tail'])} | {f1(v['z_mean'])} | {f1(v['z_sd'])} | {f1(v['a_z'])} | {f1(v['b_z'])} |")
    L += ["", "### 칸 안 ΔSNR 과 z · PIT(피험자 단위 · 순위 상관 · 기술)", "", "| 칸 | 피험자 | ρ(ΔSNR, z) | ρ(ΔSNR, PIT) | ΔSNR 중앙 |", "|---|---|---|---|---|"]
    for k, v in C.items():
        w = v["within_cell_dSNR"]
        rr = lambda x: "—" if x["rho"] is None else f"{x['rho']:+.2f}"  # noqa: E731
        ds = v["attrs"]["dSNR_median"]
        L.append(f"| {k.replace('|', ' · ')} | {w['n_subjects']} | {rr(w['rho_z'])} | {rr(w['rho_pit'])} | {'—' if ds is None else f'{ds:+.1f}'} |")
    (OUT / "d1_c3.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def figures(res: dict):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    C = res["cells"]
    y = np.array([np.nan if C[k]["z_mean"] is None else C[k]["z_mean"]["median"] for k in C])
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.6))
    for a, (f, lab) in zip(ax, (("abs_dTE", "|ΔTE| (ms)"), ("dTR_ms", "ΔTR A→B (ms)"), ("dSNR_median", "ΔSNR B−A (median)"))):
        a.scatter([C[k]["attrs"][f] for k in C], y, s=18)
        a.axhline(0, lw=0.6, color="0.5")
        a.set_xlabel(lab)
        a.set_ylabel("z mean (k median)")
    fig.suptitle("D1 — cell z mean vs acquisition differences (exploratory · descriptive)")
    fig.tight_layout()
    fig.savefig(OUT / "d1_fig_relations.png", dpi=150)
    plt.close(fig)
    obs = np.array([np.nan if C[k]["noise"]["obs_ratio_a_over_ref"] is None else C[k]["noise"]["obs_ratio_a_over_ref"]["median"] for k in C])
    pred = np.array([np.nan if C[k]["noise"]["pred_ratio_30_over_snr_pair"] is None else C[k]["noise"]["pred_ratio_30_over_snr_pair"] for k in C])
    fig, a = plt.subplots(figsize=(4.4, 4.2))
    a.scatter(pred, obs, s=18)
    lim = [0, np.nanmax(np.r_[obs, pred, 1.0]) * 1.1]
    a.plot(lim, lim, lw=0.8, color="0.4")
    a.set_xlabel("predicted ratio 30 / SNR_pair (median)")
    a.set_ylabel("observed ratio (a) / reference z SD")
    a.set_title("D1 — noise hypothesis (exploratory)")
    fig.tight_layout()
    fig.savefig(OUT / "d1_fig_noise.png", dpi=150)
    plt.close(fig)


def mode_rulethr():
    """e (review-session reply of 2026-11-02 · plan §16-2) — for each cell, applies the measurement rule to test data (the same acquisition type · the true combination is among the candidates) and recomputes the threshold (5th percentile) and the 95th percentile of the expected exceedance rate. Placed beside the values obtained on the same version (our ref samples) with the "one true combination" rule · the random numbers use the same keys as judge_v4.drawn_thresholds. The verdict is not changed (descriptive only)."""
    _, v4 = tables("rulethr")
    seed = v4["seed_numpy"]
    out = dict(rule="e — 시험 자료에 실측 코일 규칙(후보 모두 덮어야 덮음 · 후보 모두 문턱 넘어야 오지정 · 후보 모두에서 아홉 k 표본)을 적용한 문턱 · 기대 초과율 95 백분위 · "
                    "같은 표본의 참 조합 하나 규칙과 나란히 · 저장된 판정 문턱 · 실측 덮음과 나란히 · 판정 바꾸지 않음(계획서 §16-2)", seed_numpy=seed, cells={})
    for key, c4 in v4["main"]["cells"].items():
        if c4.get("verdict") == "자료 없음":
            continue
        pair, comp = key.split("|")
        ci, t = J.COMPS.index(comp), J.acq_type(pair)
        rf = dict(np.load(OUT / f"d1_ref_{t}.npz", allow_pickle=False))
        cand = combos_of(pair)
        cj = [COMBOS4.index(c) for c in cand]
        true = [tuple(x) for x in rf["true_combo"].tolist()]
        incand = np.array([x in cand for x in true])
        thr_m = float(rf["threshold"])
        jt = np.array([COMBOS4.index(x) for x in true])
        ar = np.arange(len(jt))
        cov_true, ok_true, flag_true = rf["cov"][ar, jt][:, :, ci], rf["ok"][ar, jt], rf["score"][ar, jt] > thr_m
        cov_rule = rf["cov"][:, cj][:, :, :, ci].all(1)
        ok_rule = rf["ok"][:, cj].all(1)
        flag_rule = (rf["score"][:, cj] > thr_m).all(1)
        n_cl = int(c4["n_subjects"])
        a_ = J.ARMS.index("c3")
        g = lambda kind, ci_: np.random.default_rng([seed, a_, J.S.ACQ_TYPES.index(t), ci_ + 1, n_cl, int(kind == "p95")])  # noqa: E731
        m = incand
        res = dict(acq_type=t, candidates=[f"{x}|{y}" for x, y in cand], n_cases=int(m.sum()), n_clusters=n_cl,
                   threshold_stored=c4["threshold"], threshold_fallback_stored=c4["threshold_fallback"],
                   threshold_true_combo_ours=J.threshold_from_test(cov_true[m], ok_true[m], flag_true[m], n_cl, g("thr", ci)),
                   threshold_rule_ours=J.threshold_from_test(cov_rule[m], ok_rule[m], flag_rule[m], n_cl, g("thr", ci)),
                   p95_stored=c4["exceed_p95"], p95_true_combo_ours=J.exceed_p95(flag_true[m], n_cl, g("p95", -1)), p95_rule_ours=J.exceed_p95(flag_rule[m], n_cl, g("p95", -1)),
                   coverage_real=c4["coverage"], exceed_rate_real=c4["exceed_rate"],
                   test_coverage_true_combo_kmed=float(np.median(cov_true[m & ~flag_true & ok_true.all(1)].mean(0))) if (m & ~flag_true & ok_true.all(1)).any() else None,
                   test_coverage_rule_kmed=float(np.median(cov_rule[m & ~flag_rule & ok_rule.all(1)].mean(0))) if (m & ~flag_rule & ok_rule.all(1)).any() else None,
                   test_flag_rate_true=float(flag_true[m].mean()), test_flag_rate_rule=float(flag_rule[m].mean()))
        out["cells"][key] = res
    (OUT / "d1_rulethr.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    for k, v in out["cells"].items():
        print(k, len(v["candidates"]), "thr stored %.3f true %.3f rule %.3f | real cov %.3f | test cov true %.3f rule %.3f" % (
            v["threshold_stored"], v["threshold_true_combo_ours"] or float("nan"), v["threshold_rule_ours"] or float("nan"), v["coverage_real"] or float("nan"),
            v["test_coverage_true_combo_kmed"] or float("nan"), v["test_coverage_rule_kmed"] or float("nan")))


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    dict(real=mode_real, ref=mode_ref, sens=mode_sens, summary=mode_summary, rulethr=mode_rulethr).get(cmd, lambda: sys.exit("쓰임: real|ref|sens|summary|rulethr"))()
