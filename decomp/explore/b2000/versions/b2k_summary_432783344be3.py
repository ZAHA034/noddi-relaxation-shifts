"""Earlier version of `decomp/explore/b2000/b2k_summary.py`: the version that wrote the recorded output. Summarizes the b2000 noise-floor known-answer test (calibration table, power) and produces the measurement verdict from b2k results."""
from __future__ import annotations

import collections, datetime, hashlib, json, math, pathlib, sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import b2k_core as C  # noqa: E402

OUT = pathlib.Path("<WORK>/remeasure/explore/b2000")
BANK_FILES = ("kat_bank_siemens.jsonl", "kat_bank_all.jsonl", "kat_bank_PG.jsonl")   # Names fixed (§27-2) — the earlier file first · Siemens is a fixed copy
VENDOR = {"NOT1ACH": "Philips", "NOT2ING-89": "Philips", "NOT2ING-98": "Philips", "NOT4GEP": "GE", "OXF4GEP": "GE",
          "OXF1PRI": "Siemens", "OXF2PRI": "Siemens", "OXF3TRI": "Siemens"}
EXPECT = {"NOT1ACH": 20, "NOT2ING-89": 4, "NOT2ING-98": 16, "NOT4GEP": 10, "OXF1PRI": 20, "OXF2PRI": 20, "OXF3TRI": 10, "OXF4GEP": 9}
L_GRID = (0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0)
EQUAL_L = (1.0, 2.0, 4.0, 8.0)
PAIRS = (("OXF2PRI", "OXF1PRI"), ("OXF2PRI", "OXF3TRI"), ("OXF1PRI", "OXF3TRI"))
NOMINAL_N = {"OXF1PRI": 20, "OXF2PRI": 20, "OXF3TRI": 10}
NBOOT, SEED = 2000, 20261103
# §27-2 M′
PURE_L = tuple(float(x) for x in np.exp(np.linspace(math.log(0.5), math.log(64.0), 41)))
PURE_NVOX, N2 = 100_000, 50
CONTRASTS = (("c1", "OXF2PRI", "OXF1PRI", 15), ("c2", "OXF2PRI", "OXF3TRI", 8), ("c3", "OXF1PRI", "OXF3TRI", 8))
HALF, THIRD = 0.5, 1.0 / 3.0
SIGMA_Q_FLAG, SIGMA_Q_DROP = 10.0, 3.0
SPREAD_X = 3.0


def now():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def sha(p: pathlib.Path) -> str:
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()


def jsonable(v):
    """inf · nan as strings (standard JSON · blocked with allow_nan=False)."""
    if isinstance(v, dict):
        return {str(k): jsonable(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [jsonable(x) for x in v]
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        v = float(v)
        if math.isnan(v):
            return "nan"
        if math.isinf(v):
            return "inf" if v > 0 else "-inf"
        return v
    if isinstance(v, np.bool_):
        return bool(v)
    return v


def read_lines(p: pathlib.Path) -> tuple[list[dict], int]:
    """Decodes only complete lines (a last line cut off while being appended is dropped and counted). If a line that is not the last is broken, it stops."""
    raw = p.read_bytes()
    raw = raw[: raw.rfind(b"\n") + 1]
    rows, bad = [], 0
    lines = raw.decode("utf-8").splitlines()
    for i, ln in enumerate(lines):
        try:
            rows.append(json.loads(ln))
        except json.JSONDecodeError:
            if i != len(lines) - 1:
                raise RuntimeError(f"{p}: {i + 1} 번째 줄이 깨졌다(끝 줄이 아님)")
            bad += 1
    return rows, bad


def load_bank(units: tuple[str, ...] | None = None) -> dict:
    """(unit, resample, L) → list of R (cases with too little mask are counted). For the same (unit, simulation, L, resample), only the first is kept."""
    seen, out, fail, info = set(), collections.defaultdict(list), collections.Counter(), {}
    for name in BANK_FILES:
        f = OUT / name
        if not f.exists():
            continue
        rows, bad = read_lines(f)
        info[name] = dict(sha256=sha(f), rows=len(rows), cut_last=bad)
        for r in rows:
            if units and r["unit"] not in units:
                continue
            key = (r["unit"], r["draw"], r["L"], r["resample"])
            if key in seen:
                continue
            seen.add(key)
            k = (r["unit"], r["resample"], float(r["L"]))
            if r["R"] is None:
                fail[k] += 1
            else:
                out[k].append(float(r["R"]))
    return dict(R=out, fail=fail, files=info)


def calib(B: dict, unit: str) -> dict:
    res = {}
    for how in sorted({k[1] for k in B["R"] if k[0] == unit}):
        rows = {}
        for L in L_GRID:
            a = np.asarray(B["R"].get((unit, how, L), []))
            n_fail = B["fail"].get((unit, how, L), 0)
            rows[L] = (dict(n=int(len(a)), n_fail=int(n_fail), med=float(np.median(a)), p025=float(np.percentile(a, 2.5)),
                            p975=float(np.percentile(a, 97.5)), q25=float(np.percentile(a, 25)), q75=float(np.percentile(a, 75)))
                       if len(a) else dict(n=0, n_fail=int(n_fail)))
        res[how] = rows
    return res


def separable(cal: dict) -> dict:
    """Separated = 2.5th percentile of R at L 8 > 97.5th percentile of R at L 4 (for Philips · GE, with all resampling choices pooled)."""
    c4 = [c[4.0]["p975"] for c in cal.values() if c.get(4.0, {}).get("n")]
    c8 = [c[8.0]["p025"] for c in cal.values() if c.get(8.0, {}).get("n")]
    if not c4 or not c8:
        return dict(separable=None)
    return dict(L4_p975_max=max(c4), L8_p025_min=min(c8), separable=bool(min(c8) > max(c4)))


def boot_med(rng, a: np.ndarray, n: int, reps: int) -> np.ndarray:
    return np.median(rng.choice(a, size=(reps, n), replace=True), axis=1)


def delta_min(B: dict, A: str, Bu: str, nA: int, nB: int, per_L: bool = False) -> dict:
    """Δ_min of the letter rule of §27 (R scale · gathers L 1 · 2 · 4 · 8 · exactly as defined). If per_L, it also gives the value for each L (for the report)."""
    rng = np.random.default_rng(SEED)
    d, byL = [], {}
    for L in EQUAL_L:
        a, b = np.asarray(B["R"][(A, "없음", L)]), np.asarray(B["R"][(Bu, "없음", L)])
        x = boot_med(rng, a, nA, NBOOT) - boot_med(rng, b, nB, NBOOT)
        d.append(x)
        byL[L] = float(np.percentile(np.abs(x), 97.5))
    d = np.concatenate(d)
    dm = float(np.percentile(np.abs(d), 97.5))
    out = dict(nA=nA, nB=nB, delta_min=dm)
    if per_L:
        out["delta_min_by_L"] = byL
    if A == "OXF2PRI":
        alt = boot_med(rng, np.asarray(B["R"][(A, "없음", 8.0)]), nA, NBOOT) - boot_med(rng, np.asarray(B["R"][(Bu, "없음", 4.0)]), nB, NBOOT)
        out.update(power=float((alt >= dm).mean()), alt_median=float(np.median(alt)), alt_p025=float(np.percentile(alt, 2.5)))
    return out


def kat(units: list[str]) -> dict:
    B = load_bank()
    cal = {u: calib(B, u) for u in units}
    sep = {u: separable(cal[u]) for u in units}
    vend = collections.defaultdict(list)
    for u in units:
        vend[VENDOR[u]].append(u)
    vsum = {v: dict(units=us, any_separable=any(bool(sep[u]["separable"]) for u in us)) for v, us in vend.items()}
    dm = {}
    if all(u in units for u in NOMINAL_N):
        for A, Bu in PAIRS:
            dm[f"{A}-{Bu}"] = delta_min(B, A, Bu, NOMINAL_N[A], NOMINAL_N[Bu], per_L=True)
    out = dict(made=now(), files=B["files"], script_sha256=sha(pathlib.Path(__file__)), core_sha256=sha(HERE / "b2k_core.py"),
               calibration=cal, separable=sep, vendor=vsum, delta_min_nominal=dm)
    tag = "siemens" if set(units) == set(NOMINAL_N) else ("all" if len(units) == len(VENDOR) else "_".join(units))
    (OUT / f"kat_summary_{tag}.json").write_text(json.dumps(jsonable(out), ensure_ascii=False, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    for u in units:
        for how, rows in cal[u].items():
            s = " · ".join(f"L{L:g} {r['med']:.3f} [{r['p025']:.3f}, {r['p975']:.3f}] n{r['n']}" + (f" 부족{r['n_fail']}" if r["n_fail"] else "")
                           for L, r in rows.items() if r["n"])
            print(f"{u} {how}: {s}")
        print(f"  가른다(L4 대 L8): {sep[u]}")
    for k, v in dm.items():
        print("Δ_min", k, v)
    return out


# ── §27-2 M′: adjustment g_u(L) = P(L) × b_u(L) ─────────────────────────────────────────
def pure() -> dict:
    """Median voxel ratio P(L) of pure noise N2 50 — 41 points log-spaced over L 0.5 – 64 · 10⁵ voxels per point · seed fixed."""
    rng = np.random.default_rng(SEED)
    P = []
    for L in PURE_L:
        M = C.ncchi(np.zeros((PURE_NVOX, N2)), 1.0, L, rng)
        P.append(float(np.median(M.mean(1) / M.std(1, ddof=1))))
    out = dict(made=now(), L=list(PURE_L), P=P, n_vox=PURE_NVOX, N2=N2, seed=SEED, script_sha256=sha(pathlib.Path(__file__)),
               core_sha256=sha(HERE / "b2k_core.py"))
    (OUT / "pure_curve_N2_50.json").write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    print(" · ".join(f"L{L:.3g} {p:.4f}" for L, p in zip(PURE_L, P)))
    return out


def load_pure() -> tuple[np.ndarray, np.ndarray]:
    d = json.loads((OUT / "pure_curve_N2_50.json").read_text(encoding="utf-8"))
    return np.asarray(d["L"]), np.asarray(d["P"])


def g_curve(cal_u: dict, how: str = "없음") -> tuple[np.ndarray, np.ndarray]:
    """Computes the unit calibration curve g_u on the L grid of the pure curve (0.5 – 64). b_u = bank median / P(grid L) · interpolation in log L · for L > 12, b_u(12)."""
    Lp, P = load_pure()
    rows = cal_u[how]
    Lg = np.asarray([L for L in L_GRID if rows[L]["n"]])
    Pg = np.exp(np.interp(np.log(Lg), np.log(Lp), np.log(P)))
    bg = np.asarray([rows[L]["med"] for L in Lg]) / Pg
    b = np.interp(np.log(Lp), np.log(Lg), bg)               # Outside the grid: end value held (no L < 0.5 · for L > 12, b_u(12))
    g = P * b
    if not np.all(np.diff(g) > 0):
        raise RuntimeError("보정 곡선이 단조 증가가 아니다")
    return Lp, g


def leff(R: float, Lp: np.ndarray, g: np.ndarray) -> tuple[float, str]:
    """L_eff = g⁻¹(R) · values outside the range are clamped to 0.5 · 64 (censored) and marked."""
    if R < g[0]:
        return float(Lp[0]), "아래 검열"
    if R > g[-1]:
        return float(Lp[-1]), "위 검열"
    return float(np.exp(np.interp(R, g, np.log(Lp)))), ""


def subj(name: str) -> str:
    return name.split("_")[0]


def contrast(ell: dict, A: str, Bu: str, nmin: int, rng: np.random.Generator) -> dict:
    """c = med_s ℓ(A, s) − med_s ℓ(B, s) — over the common subjects s that have values on both scanners. Interval = 2,000 resamples with replacement at the subject level."""
    common = sorted(set(ell[A]) & set(ell[Bu]))
    out = dict(n_common=len(common), subjects=common)
    if len(common) < nmin:
        out.update(c=None, ci=None, status="쌍 부족")
        return out
    a = np.asarray([ell[A][s] for s in common])
    b = np.asarray([ell[Bu][s] for s in common])
    c = float(np.median(a) - np.median(b))
    idx = rng.integers(0, len(common), size=(NBOOT, len(common)))
    bs = np.median(a[idx], axis=1) - np.median(b[idx], axis=1)
    d = a - b
    k = np.sort(d)
    lo, hi = sign_ci_ranks(len(d))
    out.update(c=c, ci=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))], status="값",
               paired_median=float(np.median(d)), paired_sign_ci=[float(k[lo - 1]), float(k[hi - 1])], paired_sign_ranks=[lo, hi])
    return out


def sign_ci_ranks(n: int) -> tuple[int, int]:
    """Order-statistic interval of the sign test (the narrowest symmetric [k, n − k + 1] with coverage ≥ 0.95)."""
    from scipy.stats import binom
    best = (1, n)
    for k in range(1, n // 2 + 1):
        cover = 1 - 2 * binom.cdf(k - 1, n, 0.5)
        if cover >= 0.95:
            best = (k, n - k + 1)
    return best


def decide(cs: dict) -> dict:
    """§27-2 items ㄹ · ㅁ (the fourth and fifth Korean-lettered items) — gate G (c₃ interval ⊂ [−½, ½]) · closeness (both sides) · size (predicted ±⅓)."""
    c1, c2, c3 = cs["c1"], cs["c2"], cs["c3"]
    if any(x["ci"] is None for x in (c1, c2, c3)):
        return dict(gate=None, side="가를 수 없음(쌍 부족)", size=None)
    G = c3["ci"][0] >= -HALF and c3["ci"][1] <= HALF
    if not G:
        return dict(gate=False, side="32 채널 대조 불일치 — 채널 수 효과로 읽지 않음", size=None)
    if c1["ci"][0] >= HALF and c2["ci"][0] >= HALF:
        side, cstar = "시뮬레이터 가정 쪽", 1.0
    elif all(-HALF < x["ci"][0] and x["ci"][1] < HALF for x in (c1, c2)):
        side, cstar = "(가) 쪽", 0.0
    else:
        return dict(gate=True, side="가를 수 없음", size=None)
    ok = all(cstar - THIRD <= x["ci"][0] and x["ci"][1] <= cstar + THIRD for x in (c1, c2))
    return dict(gate=True, side=side, size="크기 맞음" if ok else "크기 다름", c_star=cstar)


def classify_literal(d: dict, dm: dict) -> str:
    """Exactly the letter rule of §27 (supplementary · no interpretation — §27-2)."""
    d1, d2, d3 = d["OXF2PRI-OXF1PRI"], d["OXF2PRI-OXF3TRI"], d["OXF1PRI-OXF3TRI"]
    m1, m2, m3 = dm["OXF2PRI-OXF1PRI"], dm["OXF2PRI-OXF3TRI"], dm["OXF1PRI-OXF3TRI"]
    if d1 >= m1 and d2 >= m2:
        return "시뮬레이터 가정 쪽(L 8 대 4)"
    if abs(d1) < m1 and abs(d2) < m2 and abs(d3) < m3:
        return "(가) 쪽 — 같은 L 과 맞음"
    return "가를 수 없음"


def result() -> dict:
    fixed = json.loads((OUT / "kat_summary_siemens.json").read_text(encoding="utf-8"))
    B = load_bank(tuple(NOMINAL_N))
    for u in NOMINAL_N:                                      # Whether it is the same bank as the fixed calibration table
        now_cal = json.loads(json.dumps(jsonable(calib(B, u))))
        assert now_cal == fixed["calibration"][u], f"{u} 보정표가 §27-2 고정판과 다르다"
    recs, cut = {}, 0
    rows, cut = read_lines(OUT / "read_sessions.jsonl")
    for r in rows:
        recs[r["session"]] = r                               # The same session is the last line
    status = collections.Counter((r["unit"], r.get("status") or ("값" if r.get("R") is not None else "?")) for r in recs.values())
    dev, ell = {}, {}
    for ui, u in enumerate(VENDOR):
        rs = {k: r for k, r in recs.items() if r["unit"] == u}
        if not rs:
            continue
        cal = calib(load_bank((u,)), u) if u not in NOMINAL_N else calib(B, u)
        sep = separable(cal) if cal else dict(separable=None)
        e = dict(n_expected=EXPECT[u], n_records=len(rs), n_value=sum(r.get("R") is not None for r in rs.values()),
                 n_short=sum(r.get("status") == "마스크 부족" for r in rs.values()), n_failed=sum(r.get("status") == "실패" for r in rs.values()),
                 separable=sep.get("separable"))
        vals = {k: r["R"] for k, r in rs.items() if r.get("R") is not None}
        rng = np.random.default_rng([SEED, ui])
        if vals:
            v = np.asarray(list(vals.values()))
            bm = boot_med(rng, v, len(v), NBOOT)
            e.update(R_med=float(np.median(v)), R_ci=[float(np.percentile(bm, 2.5)), float(np.percentile(bm, 97.5))],
                     R_min=float(v.min()), R_max=float(v.max()), R_iqr=float(np.subtract(*np.percentile(v, [75, 25]))))
            e["L_eff"], e["sigma_q"] = {}, {}
            for how in cal:
                Lp, g = g_curve(cal, how)
                ls = {k: leff(x, Lp, g) for k, x in vals.items()}
                e["L_eff"][how] = dict(device=leff(e["R_med"], Lp, g), ci=[leff(e["R_ci"][0], Lp, g), leff(e["R_ci"][1], Lp, g)],
                                       censored=sum(1 for x in ls.values() if x[1]))
                sq = {k: math.sqrt(rs[k]["mean_M2"] / (2 * ls[k][0])) for k in vals}
                e["sigma_q"][how] = dict(min=float(min(sq.values())), median=float(np.median(list(sq.values()))),
                                         n_flag=sum(x < SIGMA_Q_FLAG for x in sq.values()), n_drop=sum(x < SIGMA_Q_DROP for x in sq.values()))
                if how == "없음" and u in NOMINAL_N:
                    ell[u] = {subj(k): math.log2(ls[k][0]) for k in vals if sq[k] >= SIGMA_Q_DROP}
                    Lnear = min(L_GRID, key=lambda L: abs(math.log(L) - math.log(e["L_eff"][how]["device"][0])))
                    sim_iqr = cal[how][Lnear]["q75"] - cal[how][Lnear]["q25"]
                    e["spread"] = dict(L_near=Lnear, real_iqr=e["R_iqr"], sim_iqr=sim_iqr, ratio=e["R_iqr"] / sim_iqr if sim_iqr > 0 else None,
                                       flag=bool(sim_iqr > 0 and e["R_iqr"] > SPREAD_X * sim_iqr))
            fw = [r["sens"]["fwf"]["R"] for r in rs.values() if ((r.get("sens") or {}).get("fwf") or {}).get("R") is not None]
            if fw:
                e["sens_fwf"] = dict(R_med=float(np.median(fw)), n=len(fw))
            if sep.get("separable") and u in NOMINAL_N:
                Lh = e["L_eff"]["없음"]["device"][0]
                s2 = [math.sqrt(r["mean_M2"] / (2 * Lh)) / r["sens"]["S2"]["sigma_snr"] for r in rs.values()
                      if r.get("R") is not None and ((r.get("sens") or {}).get("S2") or {}).get("sigma_snr")]
                if s2:
                    e["S2"] = dict(L_hat=Lh, median=float(np.median(s2)), n=len(s2), tag="R21 결과 뒤 탐색")
        dev[u] = e
    out = dict(made=now(), status=[[a, b, n] for (a, b), n in sorted(status.items())], read_cut_last=cut, devices=dev,
               bank_files=B["files"], script_sha256=sha(pathlib.Path(__file__)), core_sha256=sha(HERE / "b2k_core.py"))
    if all(u in ell for u in NOMINAL_N):
        rng = np.random.default_rng([SEED, 99])
        cs = {k: contrast(ell, A, Bu, nmin, rng) for k, A, Bu, nmin in CONTRASTS}
        out["M_prime"] = dict(contrasts=cs, decision=decide(cs))
    if all(u in dev and dev[u].get("R_med") is not None for u in NOMINAL_N):
        n = {u: dev[u]["n_value"] for u in NOMINAL_N}
        dm = {f"{A}-{Bu}": delta_min(B, A, Bu, n[A], n[Bu], per_L=True) for A, Bu in PAIRS}
        d = {f"{A}-{Bu}": dev[A]["R_med"] - dev[Bu]["R_med"] for A, Bu in PAIRS}
        out["literal_27"] = dict(note="보충 · 해석하지 않음(§27-2)", delta_min=dm, d=d,
                                 decision=classify_literal(d, {k: v["delta_min"] for k, v in dm.items()}))
    (OUT / "b2k_result.json").write_text(json.dumps(jsonable(out), ensure_ascii=False, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    return out


if __name__ == "__main__":
    if sys.argv[1] == "kat":
        kat(sys.argv[2:] or list(VENDOR))
    elif sys.argv[1] == "pure":
        pure()
    elif sys.argv[1] == "result":
        print(json.dumps(jsonable(result()), ensure_ascii=False, indent=1)[:8000])
