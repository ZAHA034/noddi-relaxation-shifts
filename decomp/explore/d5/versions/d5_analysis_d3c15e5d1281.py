"""Earlier version of `decomp/explore/d5/d5_analysis.py`: the version that wrote the recorded output. Exploratory descriptive d5 analysis with modes for scanner discrimination, between-scanner prediction, within-scanner outliers, linkage to D1 and a Markdown summary."""
from __future__ import annotations

import datetime
import hashlib
import json
import pathlib
import sys
import warnings

import numpy as np

REPO = pathlib.Path("<REPO>")
X14 = pathlib.Path("<WORK>/stage1_work/drafts/x14")
for p in ("decomp/remeasure", "decomp/stage1", "decomp/s4"):
    sys.path.insert(0, str(REPO / p))
import rconfig as R  # noqa: E402

OUT = R.WORK / "explore" / "d5"
V4 = R.WORK / "onharmony" / "v4"
D1 = R.WORK / "explore" / "d1"
REREAD = X14 / "REREAD_LOG.jsonl"
NOW = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
SHA = dict(attr_table="5184137854abc69710904889503eaf8dbb111ca9736e46af13b49067ea65b90d",
           snr_table="1714e6df9e2221df4abc59d58dc40f85c81d5e10c4406a4bcb2e36caddbd8cec")
MAIN_EXCL = {f"{a}_{b}" for a, b in R.MAIN_EXCLUDE_ONH}
COORD = ("x", "y", "ODI")
SREF = [(s, f"ses-OXF1PRI00{i}") for s in ("sub-13192", "sub-16974") for i in range(1, 7)]


def sha(p) -> str:
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()


def reread(path, what, expect=None) -> bytes:
    b = pathlib.Path(path).read_bytes()
    s = hashlib.sha256(b).hexdigest()
    with REREAD.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(dict(시각=NOW, 다시읽음=str(path), sha256=s, 바이트=len(b), 무엇=what, 스크립트=__file__), ensure_ascii=False) + "\n")
    if expect and s != expect:
        raise RuntimeError(f"⛔ {path} sha256 이 열람 기록과 다르다")
    return b


def load_feats():
    E = json.loads((OUT / "features_eddy.json").read_text(encoding="utf-8"))["sessions"]
    r = np.load(OUT / "features_roi.npz", allow_pickle=False)
    T = {k: x for k, x in zip(r["keys"].tolist(), r["X"])}
    m = np.load(OUT / "features_map.npz", allow_pickle=False) if (OUT / "features_map.npz").exists() else None
    M = None if m is None else {k: (P, g) for k, P, g in zip(m["keys"].tolist(), m["P"], m["G"])}
    return E, T, M, r["units"].tolist()


def first_sessions(E: dict) -> list[str]:
    first = {f"{a}_{b}" for a, b in R.onh_sessions()}
    return sorted(k for k in E if k in first)


def save(name, obj):
    obj = dict(made=NOW, script_sha256=sha(__file__), **obj)
    (OUT / name).write_text(json.dumps(obj, ensure_ascii=False, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)) + "\n", encoding="utf-8")


# ── ① Device discrimination ────────────────────────────────────────────────────────────────
def loso_bacc(X, y, g):
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
    from sklearn.metrics import balanced_accuracy_score
    from sklearn.model_selection import LeaveOneGroupOut, cross_val_predict
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    pipe = make_pipeline(StandardScaler(), LDA(solver="lsqr", shrinkage="auto"))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pred = cross_val_predict(pipe, X, y, groups=g, cv=LeaveOneGroupOut())
    return float(balanced_accuracy_score(y, pred)), pred


def perm_within(y, g, rng):
    yp = y.copy()
    for s in np.unique(g):
        i = np.flatnonzero(g == s)
        yp[i] = y[rng.permutation(i)]
    return yp


def device_sets(keys, E, T, M):
    Ex = np.array([[E[k][f] for f in ("E1", "E2", "E3", "E4", "E5", "E6")] for k in keys])
    Tx = np.array([T[k].ravel() for k in keys])
    Tx = Tx[:, np.isfinite(Tx).all(0)]
    sets = {"E": Ex, "T": Tx, "E+T": np.c_[Ex, Tx]}
    if M is not None:
        sets["P̄+G"] = np.array([np.r_[np.nanmedian(M[k][0][:, 0:5], 0), M[k][1]] for k in keys])
    return sets


def mode_device():
    from sklearn.metrics import confusion_matrix
    E, T, M, _ = load_feats()
    vg = OUT / "vr_gate.json"
    if M is not None and not (vg.exists() and json.loads(vg.read_text(encoding="utf-8"))["passed"]):
        M = None                                                              # If V-R does not pass, P̄ + G is not used (plan §17)
    res = dict(rule="LOSO · 표준화 → 수축 LDA · 균형 정확도 · 우연 기준 = 피험자 안 라벨 섞기 1,000 번의 95 백분위(계획서 §17 ①)", variants={})
    for label, excl in (("main", MAIN_EXCL), ("with_main_excluded", set())):
        keys = [k for k in first_sessions(E) if k not in excl]
        y = np.array([E[k]["device"] for k in keys])
        g = np.array([E[k]["sub"] for k in keys])
        sets = device_sets(keys, E, T, M)
        rng = np.random.default_rng(20261103)
        perms = [perm_within(y, g, rng) for _ in range(1000 if label == "main" else 0)]
        out = {}
        for nm, X in sets.items():
            b, pred = loso_bacc(X, y, g)
            labs = sorted(set(y))
            null = [loso_bacc(X, yp, g)[0] for yp in perms]
            out[nm] = dict(n_sessions=len(keys), n_features=int(X.shape[1]), balanced_accuracy=b, chance_95=(float(np.quantile(null, 0.95)) if null else None),
                           null_median=(float(np.median(null)) if null else None), above_chance=(bool(b > np.quantile(null, 0.95)) if null else None),
                           labels=labs, confusion=confusion_matrix(y, pred, labels=labs).tolist())
        # T_c — subject centering (subtract the subject's device mean)
        Tx = sets["T"].copy()
        for s in np.unique(g):
            Tx[g == s] -= Tx[g == s].mean(0)
        b, pred = loso_bacc(Tx, y, g)
        null = [loso_bacc(Tx, yp, g)[0] for yp in perms]
        out["T_c"] = dict(n_sessions=len(keys), n_features=int(Tx.shape[1]), balanced_accuracy=b, chance_95=(float(np.quantile(null, 0.95)) if null else None),
                          above_chance=(bool(b > np.quantile(null, 0.95)) if null else None))
        res["variants"][label] = out
    # Sub: NOT2ING-89 vs 98 (same scanner · T · one session per subject → LOO)
    keys = [k for k in first_sessions(E) if E[k]["device"] == "NOT2ING" and k not in MAIN_EXCL]
    y = np.array([E[k]["unit"] for k in keys])
    g = np.array([E[k]["sub"] for k in keys])
    X = np.array([T[k].ravel() for k in keys])
    X = X[:, np.isfinite(X).all(0)]
    b, _ = loso_bacc(X, y, g)
    rng = np.random.default_rng(20261103)
    null = [loso_bacc(X, rng.permutation(y), g)[0] for _ in range(1000)]
    res["not2ing_89_vs_98_T"] = dict(n=len(keys), n_89=int((y == "NOT2ING-89").sum()), balanced_accuracy=b, chance_95=float(np.quantile(null, 0.95)),
                                     above_chance=bool(b > np.quantile(null, 0.95)))
    save("d5_1_device.json", res)
    print(json.dumps({k: {n: round(v["balanced_accuracy"], 3) for n, v in d.items()} for k, d in res["variants"].items()}, ensure_ascii=False), res["not2ing_89_vs_98_T"])


# ── ② Prediction of correspondence between devices ──────────────────────────────────────────
def sigma_d(T) -> np.ndarray:
    """OXF1PRI repeated 12 sessions → pooled within-subject variance for each ROI · coordinate → σ_d = √2 σ_w (27, 3)."""
    ss, dfs = np.zeros((27, 3)), np.zeros((27, 3))
    for s in ("sub-13192", "sub-16974"):
        V = np.array([T[f"{s}_{ses}"] for ss_, ses in SREF if ss_ == s])        # (6, 27, 3)
        fin = np.isfinite(V)
        cnt = fin.sum(0)
        mu = np.nanmean(V, 0)
        ss += np.nansum((V - mu) ** 2, 0)
        dfs += np.maximum(cnt - 1, 0)
    return np.sqrt(2.0 * ss / np.where(dfs > 0, dfs, np.nan))


def loso_predict(A, B, kind):
    n = len(A)
    P = np.full_like(A, np.nan)
    for i in range(n):
        tr = np.arange(n) != i
        D = B[tr] - A[tr]
        if kind == "M0":
            P[i] = A[i]
        elif kind == "M1":
            P[i] = A[i] + np.nanmean(D, 0)
        elif kind == "M1g":
            P[i] = A[i] + np.nanmean(D)
        elif kind == "M2":
            for r in range(A.shape[1]):
                m = np.isfinite(A[tr, r]) & np.isfinite(B[tr, r])
                if m.sum() >= 6 and np.ptp(A[tr, r][m]) > 0:
                    beta, alpha = np.polyfit(A[tr, r][m], B[tr, r][m], 1)
                    P[i, r] = alpha + beta * A[i, r]
    return P


def mode_transfer():
    E, T, _, units = load_feats()
    at = json.loads(reread(V4 / "attr_table.json", "D5 ② — 쌍 × 피험자의 두 세션", SHA["attr_table"]))
    snr = json.loads(reread(V4 / "snr_table.json", "D5 ② — 세션 SNR(나란히 적음)", SHA["snr_table"]))["sessions"]
    sd = sigma_d(T)
    rng = np.random.default_rng(20261104)
    res = dict(rule="LOSO · M0 그대로 · M1 ROI 이동(주) · M1g 쌍 이동 · M2 ROI 일차(훈련 ≥ 6) · R = RMS(e/σ_d) · 기대 √(n/(n−1)) · 피험자 부트스트랩 2,000 · "
                    "읽기: M1 R 95 % 위 ≤ 1.5 × 기대 → 재검사 수준으로 대응됨 · 아래 > 1.5 × 기대 → 재검사 수준 밖 · 그 밖 가를 수 없음(계획서 §17 ②)",
               sigma_d_median={c: float(np.nanmedian(sd[:, j])) for j, c in enumerate(COORD)}, pairs={})
    for pair, d in sorted(at["pairs"].items()):
        subs = sorted(s for s, v in d.items() if not v["main_excluded"])
        if len(subs) < 4:
            continue
        A = np.array([T[f"{s}_{d[s]['sessions'][0]}"] for s in subs])           # (n, 27, 3)
        B = np.array([T[f"{s}_{d[s]['sessions'][1]}"] for s in subs])
        n = len(subs)
        exp = float(np.sqrt(n / (n - 1)))
        pr = dict(n_subjects=n, expected_R_M1=exp,
                  snr_A_median=float(np.median([snr[f"{s}_{d[s]['sessions'][0]}"]["snr"] for s in subs])),
                  snr_B_median=float(np.median([snr[f"{s}_{d[s]['sessions'][1]}"]["snr"] for s in subs])), coords={})
        for j, c in enumerate(COORD):
            a, b = A[:, :, j], B[:, :, j]
            sb = np.sqrt(np.nanmean(np.nanvar(b, 0, ddof=1)))
            cc = {}
            for kind in ("M0", "M1", "M1g", "M2"):
                if kind == "M2" and n - 1 < 6:
                    cc[kind] = None
                    continue
                e = b - loso_predict(a, b, kind)
                z = e / sd[None, :, j]
                Rv = float(np.sqrt(np.nanmean(z ** 2)))
                boots = []
                for _ in range(2000):
                    i = rng.integers(0, n, n)
                    boots.append(float(np.sqrt(np.nanmean(z[i] ** 2))))
                lo, hi = float(np.quantile(boots, 0.025)), float(np.quantile(boots, 0.975))
                cc[kind] = dict(R=Rv, ci95=[lo, hi], R_between=float(np.sqrt(np.nanmean(e ** 2)) / sb), n_units=int(np.isfinite(z).sum()))
                if kind == "M1":
                    cc[kind]["reading"] = ("재검사 수준으로 대응됨" if hi <= 1.5 * exp else "재검사 수준 밖" if lo > 1.5 * exp else "가를 수 없음")
            pr["coords"][c] = cc
        res["pairs"][pair] = pr
    save("d5_2_transfer.json", res)
    for p, v in res["pairs"].items():
        print(p, v["n_subjects"], {c: (round(x["M0"]["R"], 2), round(x["M1"]["R"], 2), x["M1"]["reading"]) for c, x in v["coords"].items()})


# ── ③ Within-scanner outliers · typical classification ───────────────────────────────────────────────
def robust_z(V, groups):
    Z = np.full_like(V, np.nan, dtype=float)
    for gname in np.unique(groups):
        m = groups == gname
        med = np.nanmedian(V[m], 0)
        mad = 1.4826 * np.nanmedian(np.abs(V[m] - med), 0)
        alt = 1.2533 * np.nanmean(np.abs(V[m] - med), 0)
        scale = np.where(mad > 0, mad, np.where(alt > 0, alt, np.nan))
        Z[m] = np.where(np.isfinite(scale), (V[m] - med) / scale, 0.0)
    return Z


def sub_z(V, groups, inb):
    """From the rows of inb (the subjects drawn), computes the median · MAD for each device and standardizes all rows (if MAD is 0, mean absolute deviation × 1.2533 · if that is also 0, 0)."""
    Z = np.full_like(V, np.nan, dtype=float)
    for gname in np.unique(groups):
        m = groups == gname
        mb = m & inb
        if mb.sum() < 3:
            continue
        med = np.nanmedian(V[mb], 0)
        mad = 1.4826 * np.nanmedian(np.abs(V[mb] - med), 0)
        alt = 1.2533 * np.nanmean(np.abs(V[mb] - med), 0)
        scale = np.where(mad > 0, mad, np.where(alt > 0, alt, np.nan))
        Z[m] = np.where(np.isfinite(scale), (V[m] - med) / scale, 0.0)
    return Z


def d2_full(V, groups, seed):
    from sklearn.covariance import MinCovDet
    Z = robust_z(V, groups)
    return MinCovDet(random_state=seed).fit(Z).mahalanobis(Z)


def d2_sub(V, groups, subj, rng, seed, frac=0.8):
    from sklearn.covariance import MinCovDet
    us = np.unique(subj)
    pick = rng.choice(us, int(round(frac * len(us))), replace=False)
    inb = np.isin(subj, pick)
    Z = sub_z(V, groups, inb)
    ok = np.isfinite(Z).all(1)
    d2 = np.full(len(V), np.nan)
    d2[ok] = MinCovDet(random_state=seed).fit(Z[inb & ok]).mahalanobis(Z[ok])
    return d2


def null_cutoff(groups, subj, p, q, reps, seed, sub=False):
    """q-quantile of the distances obtained by running the same procedure (within-device robust standardization → MCD · if a subsample, that version) on normal null data (same rows · devices · subject structure); χ²_p has thin tails because of small samples and per-device MAD, and so flags too many outliers (synthetic test · plan §17-1)."""
    rng = np.random.default_rng(seed)
    vals = []
    for r in range(reps):
        V = rng.standard_normal((len(groups), p))
        vals.append(d2_sub(V, groups, subj, rng, seed + r) if sub else d2_full(V, groups, seed + r))
    v = np.concatenate(vals)
    return float(np.quantile(v[np.isfinite(v)], q))


def mcd_flags(V, groups, cut, seed=20261105):
    d2 = d2_full(V, groups, seed)
    return d2, d2 > cut


def stability(V, groups, subj, reps, seed, cut_sub):
    """For each subject subsample (80 % without replacement), recomputes the median · MAD · MCD and judges all rows with the calibrated threshold → outlier fraction."""
    rng = np.random.default_rng(seed)
    hit, seen = np.zeros(len(V)), np.zeros(len(V))
    for b in range(reps):
        try:
            d2 = d2_sub(V, groups, subj, rng, seed + b)
        except Exception:
            continue
        ok = np.isfinite(d2)
        hit[ok] += d2[ok] > cut_sub
        seen[ok] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        return hit / seen, seen


def gmm_k(Z, kmax, seed=20261107):
    from sklearn.mixture import GaussianMixture
    fits = [GaussianMixture(k, covariance_type="diag", n_init=10, random_state=seed).fit(Z) for k in range(1, kmax + 1)]
    bic = [f.bic(Z) for f in fits]
    k = int(np.argmin(bic)) + 1
    return k, fits[k - 1], bic


def mode_outliers():
    from sklearn.metrics import adjusted_rand_score as ARI
    E, T, M, units = load_feats()
    gate = json.loads((OUT / "vr_gate.json").read_text(encoding="utf-8"))
    keys = sorted(k for k in E if M is not None and k in M)
    dev = np.array([E[k]["device"] for k in keys])
    subj = np.array([E[k]["sub"] for k in keys])
    logit = lambda p: np.log(p / (1 - p))  # noqa: E731
    Q = np.array([[np.log10(E[k]["E1"]), np.log10(E[k]["E2"]), logit((E[k]["E3_k"] + 0.5) / (E[k]["E3_N"] + 1)), np.log10(E[k]["E4"])] for k in keys])
    res = dict(rule="기기 안 강건 표준화 → 합동 MCD → χ²_p 0.99 · 피험자 부트스트랩 안정 ≥ 0.8(계획서 §17 ③)", n_sessions=len(keys), session={}, unit={}, typology={})
    fams = {"Q": Q}
    if gate["passed"]:
        fams["S"] = np.array([np.r_[np.nanmedian(M[k][0][:, 0:5], 0), M[k][1]] for k in keys])
    else:
        res["S_skipped"] = "V-R 실패 — P · G 를 쓰지 않음"
    flags = {}
    for nm, V in fams.items():
        cut = null_cutoff(dev, subj, V.shape[1], 0.99, 200, 20261111)
        cut_sub = null_cutoff(dev, subj, V.shape[1], 0.99, 200, 20261112, sub=True)
        d2, fl = mcd_flags(V, dev, cut)
        freq, seen = stability(V, dev, subj, 500, 20261106, cut_sub)
        stable = (freq >= 0.8) & fl
        flags[nm] = stable
        res["session"][nm] = dict(p=int(V.shape[1]), cutoff=cut, cutoff_sub=cut_sub, chi2_99=float(__import__("scipy.stats", fromlist=["chi2"]).chi2.ppf(0.99, V.shape[1])),
                                  expected_false_flags=0.01 * len(V), n_flag=int(fl.sum()), n_stable=int(stable.sum()),
                                  stable_sessions=[dict(session=k, device=dv, d2=float(x), freq=float(f)) for k, dv, x, f, s in zip(keys, dev, d2, freq, stable) if s],
                                  flagged_unstable=[dict(session=k, d2=float(x), freq=float(f)) for k, x, f, fl_, s in zip(keys, d2, freq, fl, stable) if fl_ and not s],
                                  main_excluded_check={k: dict(d2=float(x), flag=bool(f), freq=float(q)) for k, x, f, q in zip(keys, d2, fl, freq) if k in MAIN_EXCL},
                                  by_device={dv: int(stable[dev == dv].sum()) for dv in sorted(set(dev))}, d2_all={k: float(x) for k, x in zip(keys, d2)})
    res["session_stable_any"] = sorted(k for k, a in zip(keys, np.any([flags[n] for n in flags], 0)) if a)
    if gate["passed"]:
        rows = [(k, i) for k in keys for i in range(len(units)) if np.isfinite(M[k][0][i]).all()]
        V = np.array([M[k][0][i] for k, i in rows])
        gr = np.array([f"{E[k]['device']}|{units[i]}" for k, i in rows])
        su = np.array([E[k]["sub"] for k, _ in rows])
        cut = null_cutoff(gr, su, 6, 0.99, 50, 20261113)
        cut_sub = null_cutoff(gr, su, 6, 0.99, 50, 20261114, sub=True)
        d2, fl = mcd_flags(V, gr, cut)
        freq, _ = stability(V, gr, su, 200, 20261109, cut_sub)
        stable = (freq >= 0.8) & fl
        res["unit"] = dict(p=6, n_rows=len(rows), cutoff=cut, cutoff_sub=cut_sub, expected_false_flags=0.01 * len(rows), n_flag=int(fl.sum()), n_stable=int(stable.sum()),
                           stable=[dict(session=k, roi=units[i], d2=float(x), freq=float(f)) for (k, i), x, f, s in zip(rows, d2, freq, stable) if s],
                           by_device={dv: int(sum(s for (k, _), s in zip(rows, stable) if E[k]["device"] == dv)) for dv in sorted(set(dev))},
                           by_roi={u: int(sum(s for (_, i), s in zip(rows, stable) if units[i] == u)) for u in units})
        # Typical classification
        Eall = np.array([[E[k][f] for f in ("E1", "E2", "E3", "E4", "E5", "E6")] for k in keys])
        Xa = np.c_[Eall, fams["S"]]
        Xa = (Xa - Xa.mean(0)) / Xa.std(0)
        k1, f1, bic1 = gmm_k(Xa, 8)
        lab1 = f1.predict(Xa)
        Xb = np.c_[robust_z(Q, dev), robust_z(fams["S"], dev)]
        k2, f2, bic2 = gmm_k(Xb, 5)
        lab2 = f2.predict(Xb)
        rng = np.random.default_rng(20261110)
        from sklearn.mixture import GaussianMixture
        aris = []
        for _ in range(200 if k2 > 1 else 0):
            i = rng.integers(0, len(Xb), len(Xb))
            g = GaussianMixture(k2, covariance_type="diag", n_init=10, random_state=20261107).fit(Xb[i])
            aris.append(ARI(lab2, g.predict(Xb)))
        res["typology"] = dict(between=dict(k=k1, bic=bic1, ari_vs_device=float(ARI(dev, lab1)),
                                            table={str(c): {dv: int(((lab1 == c) & (dev == dv)).sum()) for dv in sorted(set(dev))} for c in range(k1)}),
                               within=dict(k=k2, bic=bic2, stability_ari_median=(float(np.median(aris)) if aris else None),
                                           sizes=[int((lab2 == c).sum()) for c in range(k2)], ari_vs_device=float(ARI(dev, lab2))))
    save("d5_3_outliers.json", res)
    print(json.dumps({n: (v["n_flag"], v["n_stable"]) for n, v in res["session"].items()}), res.get("unit", {}).get("n_stable"),
          {k: (v["k"], round(v["ari_vs_device"], 3)) for k, v in res["typology"].items()})


# ── ④ Connecting to D1 ───────────────────────────────────────────────────────────────
def mode_link():
    sys.path.insert(0, str(X14 / "D1"))
    import importlib.util
    spec = importlib.util.spec_from_file_location("d1c3", X14 / "D1" / "d1_c3.py")
    D = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(D)
    gate = json.loads((D1 / "d1_gate.json").read_text(encoding="utf-8"))
    assert gate.get("passed"), "D1 재현 문을 통과하지 않았다 — 잇지 않는다"
    s1 = json.loads((D1 / "d1_c3.json").read_text(encoding="utf-8"))
    o = json.loads((OUT / "d5_3_outliers.json").read_text(encoding="utf-8"))
    tr = json.loads((OUT / "d5_2_transfer.json").read_text(encoding="utf-8"))
    at = json.loads(reread(V4 / "attr_table.json", "D5 ④ — 쌍 × 피험자의 두 세션", SHA["attr_table"]))
    sflag = set(o["session_stable_any"])
    uflag = {(x["session"], x["roi"]) for x in o.get("unit", {}).get("stable", [])}
    rng = np.random.default_rng(20261108)
    cells = {}
    for key in s1["cells"]:
        pair, comp = key.split("|")
        ci = ("NDI", "FISO", "ODI").index(comp)
        d = dict(np.load(D1 / "d1_c3_units" / f"{pair}.npz", allow_pickle=False))
        n = len(d["flag"])
        P, Z, _ = D.pick_rule(d["pit"], d["mean"], d["sd"], np.zeros((n, 3)))
        P, Z = P[:, :, ci], Z[:, :, ci]
        subs, rois = d["subjects"], d["rois"]
        sA = np.array([f"{s}_{at['pairs'][pair][s]['sessions'][0]}" for s in subs])
        sB = np.array([f"{s}_{at['pairs'][pair][s]['sessions'][1]}" for s in subs])
        fS = np.array([a in sflag or b in sflag for a, b in zip(sA, sB)])
        fU = np.array([(a, r) in uflag or (b, r) in uflag for a, b, r in zip(sA, sB, rois)])
        base = ~d["main_excluded"]
        use = base & ~d["flag"] & d["ok"].all(1)
        out_m = d["flag"].astype(float)
        out_z = np.median(np.abs(Z), 1)
        out_t = ((P < 0.05) | (P > 0.95)).mean(1)

        def diff(sel_units, flag, y):
            f, t = sel_units & flag, sel_units & ~flag
            return (float(y[f].mean() - y[t].mean()) if f.any() and t.any() else None), int(f.sum()), int(t.sum())
        cc = {}
        for fname, flag in (("S", fS), ("U", fU)):
            r = {}
            for oname, y, sel in (("misspec", out_m, base), ("abs_z", out_z, use), ("tail", out_t, use)):
                dv, nf, nt = diff(sel, flag, y)
                r[oname] = dict(diff=dv, n_flag=nf, n_typical=nt, typical_mean=(float(y[sel & ~flag].mean()) if (sel & ~flag).any() else None),
                                all_mean=(float(y[sel].mean()) if sel.any() else None))
            cc[fname] = r
        cells[key] = dict(data=dict(fS=fS, fU=fU, base=base, use=use, m=out_m, z=out_z, t=out_t, subs=subs), res=cc)
    # Aggregation · within-cell subject bootstrap
    agg = {}
    for fname in ("S", "U"):
        for oname, yk, sk in (("misspec", "m", "base"), ("abs_z", "z", "use"), ("tail", "t", "use")):
            def stat(idx_by_cell):
                ds, ws = [], []
                for key, c in cells.items():
                    D_ = c["data"]
                    ix = idx_by_cell[key]
                    sel, fl, y = D_[sk][ix], D_[f"f{fname}"][ix], D_[yk][ix]
                    f, t = sel & fl, sel & ~fl
                    if f.any() and t.any():
                        ds.append(y[f].mean() - y[t].mean())
                        ws.append(sel.sum())
                return (float(np.median(ds)) if ds else None), (float(np.average(ds, weights=ws)) if ds else None), len(ds)
            ident = {k: np.arange(len(c["data"]["subs"])) for k, c in cells.items()}
            med0, wm0, nc = stat(ident)
            bm, bw = [], []
            for _ in range(2000):
                idx = {}
                for k, c in cells.items():
                    su = c["data"]["subs"]
                    us = np.unique(su)
                    pick = rng.choice(us, len(us), replace=True)
                    idx[k] = np.concatenate([np.flatnonzero(su == s) for s in pick])
                m_, w_, _ = stat(idx)
                if m_ is not None:
                    bm.append(m_)
                    bw.append(w_)
            ci = lambda v: [float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))] if v else None  # noqa: E731
            agg[f"{fname}|{oname}"] = dict(n_cells=nc, median_diff=med0, median_ci95=ci(bm), weighted_diff=wm0, weighted_ci95=ci(bw))
    # Pair unit — R of ② and the D1 cell summary
    from scipy.stats import spearmanr
    rel = {}
    for c in COORD:
        xs, ys = [], {f: [] for f in ("pit_mean", "tail", "z_sd")}
        for key, v in s1["cells"].items():
            pair, comp = key.split("|")
            if comp != "ODI" or pair not in tr["pairs"]:
                continue
            xs.append(tr["pairs"][pair]["coords"][c]["M1"]["R"])
            for f in ys:
                ys[f].append(np.nan if v[f] is None else v[f]["median"])
        rel[c] = {f: (float(spearmanr(xs, y, nan_policy="omit").statistic) if len(xs) >= 4 else None) for f, y in ys.items()}
    save("d5_4_link.json", dict(rule="칸마다 (특이 − 전형) → 칸 사이 중앙값 · 단위 수 가중 평균 · 칸 안 피험자 부트스트랩 2,000 · 읽기는 계획서 §17 ④",
                                cells={k: c["res"] for k, c in cells.items()}, aggregate=agg, pair_level_spearman_R_M1_vs_D1_ODI=rel))
    print(json.dumps(agg, ensure_ascii=False))


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    dict(device=mode_device, transfer=mode_transfer, outliers=mode_outliers, link=mode_link).get(cmd, lambda: sys.exit("쓰임: device|transfer|outliers|link"))()
