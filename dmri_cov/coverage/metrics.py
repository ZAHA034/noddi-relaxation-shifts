"""Compares arms by always reporting coverage, sharpness and vacuity together, because coverage alone misses intervals that are calibrated but uninformative."""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from typing import Mapping

import numpy as np

from .arms import ArmData, Correction, restrict, widen_values
from .bootstrap import band_is_informative, cluster_ci, m_curve_ci
from .sample_size import z_for_alpha

# Recorded in §107: if the set is at least this fraction of the prior support, it is counted as "공허" ("empty").
VACUITY_THRESHOLD = 0.9

__all__ = ["EvalSet", "ArmResult", "evaluate", "evaluate_m_curve", "VACUITY_THRESHOLD"]


@dataclass(frozen=True)
class EvalSet:
    """What is used in the evaluation. All must be the same length at the **case level**. truth apparent value realized in the reference acquisition (§103 — not the "true value") center posterior center (before calibration) s_post posterior SD subject subject id — the exchangeable unit (§82) feature **which ROI/feature the case is** (0-based). Required when the calibration is per feature. If absent, only a scalar calibration is allowed, and an array calibration is **rejected** — this prevents wrong values from silently coming out through broadcasting (§145). prior_measure prior support measure of each case (vacuity denominator, §107) ref_repeat_measure total measure of the **reference->reference repeat** (sharpness denominator, §110)"""

    truth: np.ndarray
    center: np.ndarray
    s_post: np.ndarray
    subject: np.ndarray
    prior_measure: np.ndarray
    ref_repeat_measure: float
    feature: np.ndarray | None = None

    def __post_init__(self):
        n = np.asarray(self.truth).shape[0]
        for nm in ("center", "s_post", "subject", "prior_measure"):
            if np.asarray(getattr(self, nm)).shape[0] != n:
                raise ValueError(f"EvalSet.{nm} 의 길이가 truth({n}) 와 다르다")
        if not np.all(np.asarray(self.s_post) > 0):
            raise ValueError("s_post 에 0 이하가 있다")
        if not np.all(np.asarray(self.prior_measure) > 0):
            raise ValueError("prior_measure 에 0 이하가 있다 — 공허성의 분모다")
        if self.feature is not None:
            f = np.asarray(self.feature)
            if f.shape[0] != n:
                raise ValueError(f"EvalSet.feature 의 길이가 truth({n}) 와 다르다")
            if f.min() < 0:
                raise ValueError("EvalSet.feature 는 0-기반 비음수 색인이어야 한다")
        if float(self.ref_repeat_measure) <= 0.0:
            raise ValueError(
                f"예리도 분모(참조->참조 반복 측도)가 {self.ref_repeat_measure} 다 — "
                "§110 이 분모를 참조 반복으로 정의했으므로 양수여야 한다"
            )


@dataclass(frozen=True)
class ArmResult:
    """The result of one arm. **It has all three** (§104)."""

    name: str
    coverage: float
    coverage_ci: tuple[float, float]
    sharpness: float
    sharpness_ci: tuple[float, float]
    vacuity: float
    vacuity_ci: tuple[float, float]
    k_median: float
    c_hat_median: float
    n_eval_subjects: int
    n_fit_subjects_min: int
    n_folds: int
    scheme: str

    def __str__(self) -> str:
        return (f"{self.name:22s} cov {self.coverage:.4f} "
                f"[{self.coverage_ci[0]:.4f},{self.coverage_ci[1]:.4f}]  "
                f"sharp {self.sharpness:.4f} "
                f"[{self.sharpness_ci[0]:.4f},{self.sharpness_ci[1]:.4f}]  "
                f"vac {self.vacuity:.4f} "
                f"[{self.vacuity_ci[0]:.4f},{self.vacuity_ci[1]:.4f}]  k {self.k_median:.4f}")


def _subset_paired(data: ArmData, keep: np.ndarray) -> ArmData:
    """Cuts the traveling data **to a set of subjects** — when only subjects outside the fold are kept. **Both** `traveling` and `paired` are cut. If only one were cut, the arms would see different subjects depending on the arm, and the matching of information content in §123-1 would break."""
    kw = {}
    if data.paired is not None:
        diff, subj = data.paired
        m = np.isin(np.asarray(subj), keep)
        kw["paired"] = (np.asarray(diff)[m], np.asarray(subj)[m])
    if data.traveling is not None:
        vals, site, subj = data.traveling
        m = np.isin(np.asarray(subj), keep)
        kw["traveling"] = (np.asarray(vals)[:, m], np.asarray(site)[m], np.asarray(subj)[m])
    return replace(data, **kw) if kw else data


def _broadcast_to_cases(x, ev: EvalSet, idx: np.ndarray, what: str):
    """Moves a **per-feature** quantity produced by the arm to the case axis. If it is an array and `ev.feature` is absent, it is **rejected** — leaving it to numpy broadcasting would silently produce wrong pairings whenever the lengths happen to match (§145)."""
    a = np.asarray(x, dtype=float)
    if a.ndim == 0:
        return np.full(idx.size, float(a))
    if ev.feature is None:
        raise ValueError(
            f"팔이 특징별 {what}({a.shape})를 냈는데 EvalSet.feature 가 없다 — "
            "어느 사례가 어느 ROI 인지 알 수 없다. `EvalSet(..., feature=...)` 를 채워라(§145)."
        )
    f = np.asarray(ev.feature)[idx]
    if f.max() >= a.size:
        raise ValueError(
            f"EvalSet.feature 의 최댓값 {int(f.max())} 이 팔의 {what} 길이 {a.size} 이상이다"
        )
    return a[f]


def _pointwise(ev: EvalSet, corr: Correction, idx: np.ndarray, *, widening: str, alpha: float):
    """Per case (covered?, width, vacuous?). ⛔ **Amendment of 2026-09-21 (§146 F1·F2)**: the earlier version squeezed the case-axis `s_post` into the feature axis with `_feature_s` and then spread it back to the case axis. In that round trip the axes became misaligned, and **the width was wrong by up to 4.96 times** (scalar v + per-case s). Now **c_hat·v are moved to the case axis first** and k is computed **for each case** — the axis conversion goes in only one direction."""
    s_cases = np.asarray(ev.s_post, dtype=float)[idx]
    c_cases = _broadcast_to_cases(corr.c_hat, ev, idx, "c_hat")
    v_cases = _broadcast_to_cases(corr.v, ev, idx, "v")
    k = widen_values(v_cases, s_cases, df=corr.df, t_scope=corr.t_scope,
                     widening=widening, alpha=alpha)
    half = z_for_alpha(alpha) * s_cases * np.asarray(k, dtype=float)
    center = np.asarray(ev.center, dtype=float)[idx] + c_cases
    truth = np.asarray(ev.truth, dtype=float)[idx]
    bad = ~np.isfinite(center) | ~np.isfinite(half)
    if bad.any():                       # F6 — NaN is not reported as "coverage 0"
        raise ValueError(
            f"보정 또는 폭에 유한하지 않은 값이 {int(bad.sum())}/{idx.size} 개 있다 "
            f"(c_hat 유한 {np.isfinite(c_cases).all()}, v 유한 {np.isfinite(v_cases).all()}). "
            "NaN 을 덮지 못한 것으로 세면 '커버리지 0' 이라는 날조된 파국이 된다(§146 F6)."
        )
    lo, hi = center - half, center + half
    cov = (lo <= truth) & (truth <= hi)
    wid = hi - lo
    vac = wid >= VACUITY_THRESHOLD * np.asarray(ev.prior_measure, dtype=float)[idx]
    return cov, wid, vac


def _arm_widening(arm, default: str) -> str:
    """If the arm has declared `widening`, **it is used** (F7). §143-3 fixed `widening="t"` as the first-order arm, but the earlier version forced the kwarg of `evaluate` on all arms, so **the setting of the arm was dead**."""
    w = getattr(arm, "widening", None)
    return w if w in ("plugin", "t") else default


def _check_layer1_s(arm, ev: EvalSet, name: str) -> None:
    """Asserts that the `s_ref` of the stratum-1 arm and the `s_post` of the evaluation set are **the same quantity** (F3·F12). The `z·s·k = t·s_d·√(1+1/n)` identity of §143-2 holds only when `s_post == s_ref`. In the code there is no path along which `s_ref` reaches `interval`, so it is **blocked here**."""
    if getattr(arm, "components", None) != "layer1":
        return
    sr = np.asarray(arm.s_ref, dtype=float)
    sp = np.asarray(ev.s_post, dtype=float)
    if sr.ndim == 0:
        ok = np.allclose(sp, float(sr))
    elif ev.feature is None:
        ok = sp.size == sr.size and np.allclose(sp, sr)
    else:
        ok = np.allclose(sp, sr[np.asarray(ev.feature)])
    if not ok:
        raise ValueError(
            f"팔 '{name}' 은 components='layer1' 인데 EvalSet.s_post 가 그 s_ref 와 다르다. "
            "§143-2 의 층 1 정의는 **s = 참조 재검사 SD** 를 전제하고, 그래야 "
            "z·s·k = t·s_d·√(1+1/n) 가 성립한다(§146 F3)."
        )


def evaluate(arms: Mapping[str, object], data: ArmData, ev: EvalSet, *,
             alpha: float = 0.10, n_folds: int | None = None,
             n_boot: int = 2000, seed: int = 20260921, widening: str = "plugin",
             ci: float = 0.95) -> dict[str, ArmResult]:
    """Evaluates the arms on the same data and the same folds. n_folds=None -> **LOSO** (default — uses the fit sample to the maximum). n_folds=K -> K-fold. ⚠ The coverage obtained here is **conditional** (§123-4). For the marginal guarantee, use `evaluate_m_curve`."""
    if not arms:
        raise ValueError("팔이 비어 있다 — 비교는 처음부터 형태에 들어 있어야 한다(§122-3)")
    if not 0.0 < alpha < 1.0:
        raise ValueError(f"alpha 는 (0,1): {alpha}")

    rng = np.random.default_rng(seed)
    subject = np.asarray(ev.subject)
    subs = np.unique(subject)
    n_sub = subs.size
    if n_sub < 2:
        raise ValueError(f"평가 피험자가 {n_sub} 명이다 — 클러스터 부트스트랩 불가(§82)")

    if n_folds is None:
        folds = [np.array([s]) for s in subs]
        scheme = "LOSO"
    else:
        k = min(int(n_folds), n_sub)
        order = rng.permutation(n_sub)
        folds = [subs[order[i::k]] for i in range(k)]
        scheme = f"{k}-fold"

    lo_q, hi_q = (1.0 - ci) / 2.0, 1.0 - (1.0 - ci) / 2.0
    z = z_for_alpha(alpha)
    out: dict[str, ArmResult] = {}

    for name, arm in arms.items():
        _check_layer1_s(arm, ev, name)
        cov_a = np.empty(0, dtype=bool)
        wid_a = np.empty(0, dtype=float)
        vac_a = np.empty(0, dtype=bool)
        sub_a = np.empty(0, dtype=subject.dtype)
        ks, chs, fit_ns = [], [], []

        for held in folds:
            fit_subs = np.setdiff1d(subs, held)
            if fit_subs.size < 2:
                raise ValueError(
                    f"겹 하나를 빼면 적합 피험자가 {fit_subs.size} 명이다(§82). "
                    "피험자를 늘리거나 겹을 줄여라."
                )
            view = restrict(_subset_paired(data, fit_subs), arm.requires, name)
            corr = arm.fit(view)
            idx = np.flatnonzero(np.isin(subject, held))
            if idx.size == 0:
                continue
            c, w, v = _pointwise(ev, corr, idx, widening=_arm_widening(arm, widening),
                                 alpha=alpha)
            cov_a = np.concatenate([cov_a, c])
            wid_a = np.concatenate([wid_a, w])
            vac_a = np.concatenate([vac_a, v])
            sub_a = np.concatenate([sub_a, subject[idx]])
            ks.append(float(w.mean() / (2.0 * z * np.asarray(ev.s_post)[idx].mean())))
            chs.append(corr.c_hat)
            fit_ns.append(int(corr.n_units))

        point = (float(cov_a.mean()),
                 float(np.median(wid_a) / ev.ref_repeat_measure),
                 float(vac_a.mean()))

        # Bootstrap: resamples **only the evaluation subjects** (on one side only). The bootstrap module blocks this by machine.
        kw = dict(n_boot=n_boot, seed=seed, ci=ci, side="eval", resampled_sides=("eval",))
        b_cov = cluster_ci(cov_a.astype(float), sub_a, statistic=np.mean, **kw)
        b_shp = cluster_ci(wid_a, sub_a,
                           statistic=lambda x: float(np.median(x) / ev.ref_repeat_measure), **kw)
        b_vac = cluster_ci(vac_a.astype(float), sub_a, statistic=np.mean, **kw)

        out[name] = ArmResult(
            name=name,
            coverage=point[0], coverage_ci=(b_cov["lo"], b_cov["hi"]),
            sharpness=point[1], sharpness_ci=(b_shp["lo"], b_shp["hi"]),
            vacuity=point[2], vacuity_ci=(b_vac["lo"], b_vac["hi"]),
            k_median=float(np.nanmedian(ks)) if ks else float("nan"),
            c_hat_median=float(np.median(chs)) if chs else float("nan"),
            n_eval_subjects=int(n_sub),
            n_fit_subjects_min=int(min(fit_ns)) if fit_ns else 0,
            n_folds=len(folds), scheme=scheme,
        )
    return out


def evaluate_m_curve(arms: Mapping[str, object], data: ArmData, ev: EvalSet, *,
                     m_values=(4, 8, 12, 16), n_repeats: int = 50, alpha: float = 0.10,
                     seed: int = 20260921, widening: str = "plugin",
                     predicted_effect: float | None = None) -> dict:
    """The test fixed in §123-4 — **the marginal coverage curve over the cohort size m**. In each repeat, **m of the traveling subjects are drawn and the calibration is fitted**, and it is evaluated on subjects **that do not overlap with that subsample**. Unlike LOSO, the calibration error η is **realized anew in each repeat**, so the marginal guarantee is measured. `predicted_effect` — the reporting format fixed in §139-4. If the band is wider than this value, it is marked `power_ok=False`. This is built so that the result **carries along the fact that passing does not mean "the formula was right"** (§89 prior block). ⚠ **No default value is set.** This value differs from comparison to comparison — the plug-in deficit of §139-4 is 0.0067 (df=8, decomposed form), but in the stratum-1 definition of §143-2 it is **about 0.15**. A fixed default would create a misunderstanding in the other direction. If it is not declared, `power_ok=None` and `power_note` carries **"예측 결손 미선언" (predicted deficit not declared)** — it is left blank rather than hidden. Returns: {arm name: {m: {"coverage", "band", "k", "n_eval", "power_ok", "power_note"}}}"""
    rng = np.random.default_rng(seed)
    for _n, _a in arms.items():
        _check_layer1_s(_a, ev, _n)
    subject = np.asarray(ev.subject)
    subs = np.unique(subject)
    out: dict[str, dict] = {name: {} for name in arms}
    for m in m_values:
        if m + 1 > subs.size:
            raise ValueError(f"m={m} 에 평가 피험자까지 두려면 {m+1}명 이상 필요하다(있는 수 {subs.size})")
        per_arm = {name: {"cov": [], "k": [], "n": []} for name in arms}
        for _ in range(n_repeats):
            perm = rng.permutation(subs)
            fit_subs, eval_subs = perm[:m], perm[m:]
            idx = np.flatnonzero(np.isin(subject, eval_subs))
            if idx.size == 0:
                continue
            for name, arm in arms.items():
                view = restrict(_subset_paired(data, fit_subs), arm.requires, name)
                corr = arm.fit(view)
                c, w, _v = _pointwise(ev, corr, idx,
                                      widening=_arm_widening(arm, widening), alpha=alpha)
                per_arm[name]["cov"].append(float(c.mean()))
                per_arm[name]["k"].append(
                    float(w.mean() / (2.0 * z_for_alpha(alpha) * np.asarray(ev.s_post)[idx].mean()))
                )
                per_arm[name]["n"].append(int(np.unique(subject[idx]).size))
        for name in arms:
            d = per_arm[name]
            if not d["cov"]:
                out[name][m] = {"coverage": float("nan"), "band": None, "k": float("nan"),
                                "n_eval": 0, "power_ok": None, "power_note": "반복 0회"}
                continue
            band = m_curve_ci(d["cov"])          # Repeat quantiles — not resampling (§144)
            if predicted_effect is None:
                ok, note = None, ("예측 결손 미선언 — 밴드가 의미 있는 폭인지 판정할 수 없다. "
                                  "비교마다 다르므로 호출자가 선언해야 한다(§144).")
            else:
                ok, note = band_is_informative(band, predicted_effect)
            out[name][m] = {"coverage": band["point"], "band": (band["lo"], band["hi"]),
                            "k": float(np.mean(d["k"])), "n_eval": int(np.mean(d["n"])),
                            "power_ok": ok, "power_note": note}
    return out
