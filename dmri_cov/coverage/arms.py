"""Defines the comparison arms and mechanically enforces what each arm may access from the target."""

from __future__ import annotations

import math
from dataclasses import dataclass, fields
from typing import Protocol, runtime_checkable

import numpy as np

from .sample_size import _k_none, k_corrected, z_for_alpha

# The full set of fields the harness can provide. An arm declares only some of them with `requires`.
#:
# ⚠ The previous version was a **hand-written list**, and when a field was added to `ArmData`, `restrict`
# rejected as an "알 수 없는 필드" ("unknown field"). The **structure in which the list and the dataclass drift apart** was itself the accident.
# Now **derived** from `ArmData` — there is no need to match two places.
# (This is also why the definition moved down under `ArmData`.)


@dataclass(frozen=True)
class ArmData:
    """All data held by the harness. **It is not handed to an arm whole** — it passes through `restrict`. paired : (diff, subject_id) — paired differences of the traveling subjects and the subject ids. By the rule of §82 the independent unit is the **subject**, not the voxel."""

    source_cal: np.ndarray | None = None
    target_scores: np.ndarray | None = None
    target_covariates: np.ndarray | None = None
    site_labels: np.ndarray | None = None
    # (values, site, subject) — the **raw measurements** of traveling subjects. values is (feature, sample).
    # Because our arm and the ComBat arm require **this same field**, the amount of information is matched mechanically (§123-2).
    traveling: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None
    # (diff, subject) — convenience field derived from `traveling`. Used only for the simple test.
    paired: tuple[np.ndarray, np.ndarray] | None = None
    rho_ball: float | None = None
    # §186-4 **per-session** acquisition coordinates. `{site: [(TE_ms, TR_s), ...]}` or representative values `{site: (TE, TR)}`.
    # Requires only the physics arm — even **without** traveling data, this alone yields the adjustment amount (ladder ①).
    acquisitions: dict | None = None
    # **Prior sample** of ΔR2 for the physics arm to use. No point estimate is hard-coded (§12/§16).
    dr2_draws: np.ndarray | None = None


ALL_FIELDS = frozenset(f.name for f in fields(ArmData))


class _RestrictedView:
    """A view that dies at once if an undeclared field is touched. A machine, not a convention (§122-3)."""

    __slots__ = ("_data", "_allowed", "_arm")

    def __init__(self, data: ArmData, allowed: frozenset[str], arm: str):
        object.__setattr__(self, "_data", data)
        object.__setattr__(self, "_allowed", allowed)
        object.__setattr__(self, "_arm", arm)

    def __getattr__(self, name: str):
        allowed = object.__getattribute__(self, "_allowed")
        arm = object.__getattribute__(self, "_arm")
        if name not in allowed:
            known = {f.name for f in fields(ArmData)}
            if name in known:
                raise AttributeError(
                    f"팔 '{arm}' 은 '{name}' 을 requires 에 선언하지 않았다 — 접근 금지(§123-1). "
                    f"선언된 것: {sorted(allowed)}"
                )
            raise AttributeError(f"'{name}' 은 ArmData 의 필드가 아니다")
        return getattr(object.__getattribute__(self, "_data"), name)

    def __setattr__(self, name, value):
        raise AttributeError("팔은 데이터를 수정할 수 없다")


def restrict(data: ArmData, requires: frozenset[str], arm_name: str) -> _RestrictedView:
    bad = requires - ALL_FIELDS
    if bad:
        raise ValueError(f"팔 '{arm_name}' 이 알 수 없는 필드를 요구한다: {sorted(bad)}")
    return _RestrictedView(data, requires, arm_name)


@dataclass(frozen=True)
class Correction:
    """The estimated calibration shift and the **predictive variance for a new subject**. ⚠ **Amendment of 2026-09-21 (§141)**: the `v` of the earlier version was only `Var(ĉ)`. That is the **confidence interval of the mean shift**, whereas what we need is the **prediction interval of the shift of a new subject**. A new subject has its own shift deviation σ²_het, and **that term does not disappear however much n is increased.** v = σ²_het + Var(ĉ) ~~~~~~ ~~~~~~~ floor shrinks as 1/n"""

    # Point adjustment amount. **A scalar or a `(feature,)` array** (§141-2) — the type annotation
    # `float` is a notational convenience; in reality it may differ for each ROI (§146 F22).
    c_hat: float
    # **Prediction variance** = sigma2_het + v_chat. A scalar or a `(feature,)` array.
    v: float
    n_units: int       # Number of independent units = number of subjects (§82)
    df: int            # Degrees of freedom of the v estimate = n_units − 1
    sigma2_het: float = 0.0  # Between-subject shift heterogeneity. **A floor independent of n**
    v_chat: float = 0.0      # Estimation error of ĉ itself (= σ²_between/n)
    # To **where** in the width the t quantile is applied.
    # "v_only" — s_post is a known quantity (NPE posterior SD) and only v is estimated (decomposed form).
    # "all"    — the **whole** prediction variance is estimated with n−1 degrees of freedom (stratum 1). Then the interval
    # is exactly the same as the textbook paired-difference prediction interval `d̄ ± t_{n−1}·s_d·√(1+1/n)`.
    t_scope: str = "v_only"
    # Value **before** clipping at 0 (stratum 1 `Var(d)(1+1/n) − s_ref²`). `None` on a path that does not clip.
    # ⛔ Author «승인 넷» («four approvals») 2 ② — the same rule as «negative not truncated · raw value reported» of the variance decomposition: the value before truncation is output together
    # Counts the number of times truncation applied (the number of elements with `v_raw < 0`).
    v_raw: float | None = None


def paired_from_traveling(values: np.ndarray, site: np.ndarray, subject: np.ndarray,
                          ref: str, target: str) -> tuple[np.ndarray, np.ndarray]:
    """Builds the **within-subject** paired differences from the traveling measurements. **Sign convention: `diff = ref − target`.** That is, **adding** `diff` to a value measured at the target gives the reference coordinate — consistent with `interval()` doing `center + c_hat`. (Amendment of 2026-09-21: the first version had the sign reversed, `target − ref`.) Only subjects scanned on both scanners are used. diff is two-dimensional, `(feature, n_paired_subject)` — **the calibration differs for each ROI.** Returns: `(diff, subject, within_var, bias_f, inv_n_sess)` — a 5-tuple (§146 F21). within_var **measurement noise** measured in repeated sessions (None if absent) bias_f de-biasing factor of `full`: Var(per-subject value) = σ²_het + within_var·bias_f inv_n_sess session-mean factor of `layer1`, `1/n_ref + 1/n_tgt` (2.0 if each has 1 session)"""
    values = np.asarray(values, dtype=float)
    site = np.asarray(site, dtype=str)
    subject = np.asarray(subject)
    if values.ndim != 2 or values.shape[1] != site.size or site.size != subject.size:
        raise ValueError(
            f"values 는 (feature, sample) 이고 site·subject 와 표본 축이 같아야 한다: "
            f"{values.shape}, {site.shape}, {subject.shape}"
        )
    diffs, subs, within, inv_n = [], [], [], []
    for u in np.unique(subject):
        m = subject == u
        a = values[:, m & (site == ref)]
        b = values[:, m & (site == target)]
        if a.shape[1] == 0 or b.shape[1] == 0:
            continue  # A subject imaged on one side only is not a traveling pair
        diffs.append(a.mean(axis=1) - b.mean(axis=1))   # ref − target, sessions averaged
        subs.append(u)
        inv_n.append(1.0 / a.shape[1] + 1.0 / b.shape[1])   # How much the session mean reduced the noise
        # If there are repeat sessions, the **measurement noise** is measured separately here (for σ²_het de-biasing, §141)
        for arr in (a, b):
            if arr.shape[1] > 1:
                within.append(arr.var(axis=1, ddof=1))
    if not diffs:
        raise ValueError(f"'{ref}' 와 '{target}' 을 모두 찍은 피험자가 없다")
    wv = np.mean(within, axis=0) if within else None
    iv = float(np.mean(inv_n))
    # (diff, subject, within_var, bias_f, inv_n_sess)
    # bias_f — de-biasing factor of `full`: Var(per-subject value) = σ²_het + within_var·bias_f
    # inv_n_sess — session-averaging factor of `layer1` (2.0 if each has 1 session)
    # On the traveling path the two are **the same** (there is only one value per subject and that value is the session mean).
    return np.stack(diffs, axis=1), np.asarray(subs), wv, iv, iv


@runtime_checkable
class Arm(Protocol):
    name: str
    requires: frozenset[str]

    def fit(self, view) -> Correction: ...


class NoCorrectionArm:
    """Arm 0 — does nothing. The baseline that shows the size of the problem (§122-2)."""

    name = "no-correction"
    requires = frozenset()

    def fit(self, view) -> Correction:  # noqa: ARG002
        return Correction(c_hat=0.0, v=0.0, n_units=0, df=0)


class TravelingArm:
    """Arm 4 — **our prescription**. Estimates the calibration shift directly from paired traveling differences and **propagates its variance into the set**. With `propagate_variance=False`, only the point calibration is done and no width is bought — that place is **the same place as the point calibration of ComBat**, and the difference between the two is exactly the contribution of C4 (§123-2)."""

    name = "traveling"
    # Must be **exactly the same set** as the ComBat arm — the amount of information is matched by machine (§123-2).
    # (`paired` is a convenience field derived deterministically from `traveling`, so it adds no information.)
    requires = frozenset({"source_cal", "traveling", "paired"})

    def __init__(self, *, propagate_variance: bool = True, widening: str = "t",
                 components: str = "layer1", ref: str = "OXF1PRI", target: str | None = None,
                 s_ref: float | np.ndarray | None = None):
        """components — the variance components to put into the width. layer1 : **the default and simplest, safest form of stratum 1** (§143-2). `v = Var(d)·(1 + 1/n) − s_ref²`. Here `s` is the **reference retest SD** (the **same quantity** as the denominator of the §110 sharpness) and `d` is the per-subject paired difference. The error of a new subject j is `d_j − ĉ`, and `d_j` has the **same distribution** as the `d_i` of the training subjects, so the predictive variance is exactly `Var(d)(1+1/n)` — **heterogeneity and the noise of both scanners are already in it.** There is no need to estimate the target noise separately, no need to assume that the two scanner noises are equal, and a scanner with no repeat-scan subjects (NOT3GEM) is no problem. `s_ref` must be given. full : σ²_het + Var(ĉ) — decomposed form. Used in designs (simulations) where σ²_het can be separated out using repeated sessions. In stratum 1, `layer1` is safer. chat_only : only Var(ĉ) — **the confidence interval of the mean**. The (wrong) behavior before 2026-09-21, used only for the missing-term demonstration ablation of §141. none : 0 — the same place as the plug-in of ComBat (§123-2). It is the interval of a user who believes "once calibrated, only the reference retest noise is left". Why the default of widening is `"t"` — if the variance is **estimated** with `n−1` degrees of freedom and z is used, the interval is biased downward in small cohorts (§143-3). t₀.₉₅(df=3)=2.353 vs z=1.645."""
        if widening not in ("plugin", "t"):
            raise ValueError(f"widening 은 'plugin' 또는 't': {widening}")
        if components not in ("layer1", "full", "chat_only", "none"):
            raise ValueError(f"components 는 layer1/full/chat_only/none: {components}")
        if components == "layer1":
            if s_ref is None:
                raise ValueError(
                    "components='layer1' 은 참조 재검사 SD `s_ref` 를 요구한다 — "
                    "§110 예리도 분모와 **같은 양**이어야 한다(§143-2)."
                )
            sr = np.asarray(s_ref, dtype=float)
            if np.any(~np.isfinite(sr)) or np.any(sr <= 0.0):
                raise ValueError(
                    f"s_ref 는 유한한 양수여야 한다: {s_ref}. 0 이나 음수를 주면 "
                    "s_ref² 를 못 빼서 v 가 s_post² 위에 **이중계산**된다(§146 F17)."
                )
        self.s_ref = s_ref
        self.propagate_variance = propagate_variance
        self.widening = widening
        self.components = "none" if not propagate_variance else components
        self.ref = ref
        self.target = target

    def _pairs(self, view):
        """If `traveling` is present, pairs are made from it; if not, `paired` is used."""
        allowed = object.__getattribute__(view, "_allowed")
        if "traveling" in allowed and view.traveling is not None:
            vals, site, subj = view.traveling
            tgt = self.target
            if tgt is None:
                others = [c for c in np.unique(np.asarray(site, dtype=str)) if c != self.ref]
                if len(others) != 1:
                    raise ValueError(
                        f"target 을 지정해야 한다 — 순회 데이터에 스캐너가 "
                        f"{sorted(others)} 로 여럿이다."
                    )
                tgt = others[0]
            return paired_from_traveling(vals, site, subj, self.ref, tgt)
        if "paired" not in allowed or view.paired is None:
            raise ValueError("traveling 도 paired 도 없다 — 순회 보정을 적합할 수 없다")
        # Convenience field `paired` = (diff, subject) one-dimensional. Folded per subject into (1, n_subj).
        #
        # ⚠ `paired` is **already a difference**, so it carries no session structure. Several values per subject
        # Since it is a voxel/ROI and not a session, counting it as a session would make `inv_n` of layer1
        # the adjustment fires wrongly (§146). So with **`inv_n = 2.0` (1 session each)**
        # fixed, and returns `within_var=None` — the adjustment term becomes exactly 0.
        # If a session structure is needed, `traveling` must be used.
        d, sb = view.paired
        d = np.asarray(d, dtype=float); sb = np.asarray(sb)
        uu = np.unique(sb)
        per = np.array([[d[sb == u].mean() for u in uu]])        # (1, n_subj)
        cnt = np.array([int((sb == u).sum()) for u in uu], dtype=float)
        if np.all(cnt <= 1):
            return per, uu, None, 0.0, 2.0
        ss = sum(float(((d[sb == u] - per[0, j]) ** 2).sum()) for j, u in enumerate(uu) if cnt[j] > 1)
        dof = float(cnt[cnt > 1].sum() - (cnt > 1).sum())
        wv = np.array([ss / dof]) if dof > 0 else None
        # bias_f = 1/m̄ (averaged over the number of observations per subject), inv_n_sess = 2.0 (no session information)
        return per, uu, wv, float(np.mean(1.0 / cnt)), 2.0

    def fit(self, view) -> Correction:
        diff, subj, within_var, bias_f, inv_n = self._pairs(view)
        diff = np.atleast_2d(np.asarray(diff, dtype=float))   # (feature, n_subject)
        subj = np.asarray(subj)
        n = subj.size
        if diff.shape[1] != n:
            raise ValueError(f"짝 차분 {diff.shape} 와 피험자 {n} 명이 안 맞는다")
        if n < 2:
            raise ValueError(
                f"독립 단위(피험자)가 {n} 명이다 — 보정량의 분산을 추정할 수 없다(§82). "
                "voxel/ROI 수로 n 을 세면 안 된다."
            )
        c_hat = diff.mean(axis=1)
        s2_between = diff.var(axis=1, ddof=1)   # = σ²_het + σ²_noise·E[1/n_a + 1/n_b]

        # σ²_het de-biasing (§141): repeated sessions identify the measurement noise separately.
        if within_var is None or bias_f == 0.0:
            sigma2_het = s2_between          # Conservative — if the noise cannot be subtracted, it is left on the high side
        else:
            sigma2_het = np.maximum(0.0, s2_between - np.asarray(within_var) * bias_f)
        v_chat = s2_between / n

        v_raw = None
        if self.components == "none":
            v = np.zeros_like(c_hat)
        elif self.components == "chat_only":
            v = v_chat
        elif self.components == "layer1":
            # Prediction variance = Var(d_new) + Var(ĉ).  d_new is **1 session each**, so Var(d_new)
            # = σ²_het + 2σ²_noise, and if the fitted subject averaged the sessions
            # s2_between = σ²_het + σ²_noise·E[1/n_a+1/n_b], so it is **shrunken** (§146 F5).
            sr2 = np.asarray(self.s_ref, dtype=float) ** 2
            v_pred = s2_between * (1.0 + 1.0 / n)
            if abs(inv_n - 2.0) > 1e-9:          # The fit side has repeat sessions
                if within_var is None:
                    raise ValueError(
                        f"적합 피험자가 스캐너당 평균 {2/inv_n:.2f} 세션인데 측정 잡음을 "
                        "따로 잴 수 없다(반복 세션 없음). layer1 의 유도는 각 1세션 전제다 — "
                        "세션을 하나만 쓰거나 components='full' 을 써라(§146 F5)."
                    )
                v_pred = v_pred + np.asarray(within_var, dtype=float) * (2.0 - inv_n)
            v_raw = v_pred - sr2
            v = np.maximum(0.0, v_raw)
        else:
            v = sigma2_het + v_chat
        sq = lambda a: float(a[0]) if a.size == 1 else a
        return Correction(c_hat=sq(c_hat), v=sq(np.asarray(v, dtype=float)),
                          n_units=int(n), df=int(n - 1),
                          sigma2_het=sq(np.asarray(sigma2_het, dtype=float)),
                          v_chat=sq(v_chat),
                          t_scope="all" if self.components == "layer1" else "v_only",
                          v_raw=None if v_raw is None else sq(np.asarray(v_raw, dtype=float)))


def widen_values(v, s_post, *, df: int = 0, t_scope: str = "v_only",
                 widening: str = "plugin", alpha: float = 0.10):
    """Width multiplier k — **`v` and `s_post` must be on the same axis** (usually the case axis). ⛔ **Amendment of 2026-09-21 (§146)**: the earlier version mixed `v` on the feature axis with `s` on the case axis, and `metrics._feature_s` changed the axis in the middle; as a result, with scalar `v` + per-case `s`, **the width was wrong by up to 4.96 times** (without exception). Now the axis conversion is finished **by the caller first** and this function computes only element by element. plugin : `k = √(1 + v/s²)` (§106 closed form) t : reflects the estimation error of v through the degrees of freedom. `t_scope="all"` — **all** of the predictive variance is estimated (stratum 1) → t on the whole width. **Even with `v = 0`, `t/z` remains** (§146 F4). `t_scope="v_only"` — `s_post` is a known quantity, so only on the v component."""
    if widening not in ("plugin", "t"):
        raise ValueError(f"widening 은 'plugin' 또는 't': {widening}")
    s_post = np.asarray(s_post, dtype=float)
    v = np.asarray(v, dtype=float)
    if np.any(s_post <= 0.0):
        raise ValueError(f"s_post 는 양수여야 한다: {s_post}")
    if np.any(v < 0.0):
        raise ValueError(f"v 는 음수일 수 없다: {v}")
    if v.ndim and s_post.ndim and v.shape != s_post.shape:
        raise ValueError(
            f"v {v.shape} 와 s_post {s_post.shape} 의 축이 다르다 — 축 변환은 호출자가 "
            "먼저 끝내야 한다(§146)."
        )
    ratio = np.sqrt(1.0 + v / s_post**2)
    if widening == "plugin":
        k = ratio
    else:
        if df < 1:
            raise ValueError(f"자유도가 {df} 라 t 보정을 쓸 수 없다")
        from scipy.stats import t as _t
        z = z_for_alpha(alpha)
        tq = float(_t.ppf(1.0 - alpha / 2.0, df))
        if t_scope == "all":
            k = (tq / z) * ratio          # t/z remain even when v=0
        else:
            k = np.sqrt(s_post**2 + (tq / z) ** 2 * v) / s_post
    return float(k) if np.ndim(k) == 0 else k


def widen(corr: Correction, s_post, *, widening: str = "plugin", alpha: float = 0.10):
    """Convenience wrapper of `Correction`. The **axes of `corr.v` and `s_post` must be the same**."""
    return widen_values(corr.v, s_post, df=corr.df, t_scope=corr.t_scope,
                        widening=widening, alpha=alpha)


def interval(
    center: float, s_post: float, corr: Correction, *, widening: str = "plugin", alpha: float = 0.10
) -> tuple[float, float]:
    """Prediction interval after applying the calibration and widening. Moves the center by `c_hat` and gives the half-width as `z·s·k`."""
    k = widen(corr, s_post, widening=widening, alpha=alpha)
    half = z_for_alpha(alpha) * np.asarray(s_post) * np.asarray(k)
    c = np.asarray(center) + np.asarray(corr.c_hat)
    return (c - half, c + half)


def k_without_correction(bias: float, s_post: float, alpha: float = 0.10) -> float:
    """Width multiplier needed to keep nominal coverage **when no calibration is done** (`k_none` of §106). Compared with the k of `interval`, it shows at once whether calibration helps or harms — if `widen(...) < k_without_correction(...)`, it helps (the harm threshold of §106)."""
    b = abs(float(bias))
    if b == 0.0:
        return 1.0
    return _k_none(b / s_post, alpha)


# **Per-axis transfer coefficient `G_axis`** — the share of the physical-formula shift that actually appears in the estimator output.
#:
# §224 showed it structurally and §224-2 confirmed it by measurement: after b0 normalization, the acquisition parameters (TE·TR)
# Because it enters the signal **only through the two scalars `(ν_app, f_app)`**, the fitter cannot distinguish the T2 path (TE) and
# **no means to distinguish** the T1 path (TR). → **one coefficient per axis, not per term.**
#:
# Measurement (four d∥ × two axes = eight pairs, agreeing within 0.2~3.1%):
# axis   d∥    TE path      T1 path      **G = mean of the two paths**
#:   NDI  1.7   0.8967       0.9126       0.9047
#:   NDI  2.0   0.7463       0.7697       0.7580
# NDI  2.26  0.7040       0.7027       **0.7034**   ← Howard tissue value
#:   NDI  2.5   0.6439       0.6531       0.6485
#:   FISO 1.7   0.9777       0.9469       0.9623
#:   FISO 2.0   0.8791       0.8600       0.8696
#:   FISO 2.26  0.7995       0.8013       **0.8004**
#:   FISO 2.5   0.7350       0.7581       0.7466
#:
# ⛔ **A value measured under no noise** (§209-3). The isotropic path moved only +2~4% at SNR 25, but
# **The tissue path cannot be measured at SNR 25** (the physical term is only 0.0198 logit, so SE 0.68). Of the noise condition
# `G` is **not known** — when citing it, always write «value measured under no noise».
# ⛔ **The «one statement» structure itself rests on a noise-free premise** (§224-1): with noise, TE also changes the absolute signal strength (SNR)
# changes it, and that path is not explained by the apparent fraction. **Under noise, an SNR path is added.**
G_TABLE: dict[str, dict[float, float]] = {
    "ndi":  {1.7: 0.9047, 2.0: 0.7580, 2.26: 0.7034, 2.5: 0.6485},
    "fiso": {1.7: 0.9623, 2.0: 0.8696, 2.26: 0.8004, 2.5: 0.7466},
}

# Discrepancy between the two path measurements |TE path − T1 path| — recorded as the **residual uncertainty** of `G` (§224-2).
G_SPREAD: dict[str, dict[float, float]] = {
    "ndi":  {1.7: 0.0159, 2.0: 0.0234, 2.26: 0.0013, 2.5: 0.0092},
    "fiso": {1.7: 0.0308, 2.0: 0.0191, 2.26: 0.0018, 2.5: 0.0231},
}

# **Recorded support interval** of tissue d∥ (§209-2). 1.7 is the fitter's fixed value, not the support of the tissue prior.
DPAR_SUPPORT = (2.0, 2.5)


def gate(axis: str, tissue_dpar: float) -> float:
    """`G_axis(operating point)` — not a constant but a **function of the tissue d∥** (§224 (a-2)). Outside the grid **it does not interpolate** — that would create information that does not exist (the rule of §208-3). It accepts only values on the grid; anything outside raises ValueError."""
    if axis not in G_TABLE:
        raise ValueError(f"축은 {sorted(G_TABLE)} 중 하나여야 한다: {axis!r}")
    tab = G_TABLE[axis]
    key = min(tab, key=lambda k: abs(k - tissue_dpar))
    if abs(key - tissue_dpar) > 1e-9:
        raise ValueError(
            f"`G` 는 d∥ {sorted(tab)} 에서만 측정됐다. {tissue_dpar} 은 보간 대상이 아니다 — "
            "그 d∥ 에서 직접 재거나, 양 끝값을 병기해 포락선으로 보고하라(§212-3)."
        )
    return tab[key]


def _rho(tr: float, t1: float) -> float:
    """`ρ = 1 − exp(−TR/T1)` — the T1 saturation factor."""
    return -float(np.expm1(-tr / t1))


class PhysicsArm:
    """Arm 3 — **ladder (1): 0 traveling subjects.** Predicts the calibration shift from the acquisition coordinates (TE·TR) alone. **The formula differs for each axis** (§193-1 · §224-3). The axis must be specified. **NDI axis** (§118, purely logit-additive):: logit ν_app = logit ν + TE·ΔR2 + ln(ρ_in/ρ_en) **FISO axis** (§73-1 · §224-3, **not logit-additive**):: logit f_app = logit f + ln(w_iso / A), A = ν·w_in + (1−ν)·w_en = logit f + [−TE/T2_iso − ln Ã(TE)] + ln(ρ_iso/ρ_tissue) ⛔ **The TE term of the FISO axis is not the NDI formula** — sign and size both differ and it **depends on `ν`**. From TE 92→68 the NDI term is −0.2066 while the FISO term is −0.2855 to −0.3575 (ν 0.85 to 0.40). At `ν = 0.7335` it becomes **−0.3006**, reproducing the value fixed in §118-1. **Transfer**: the **sum of the two terms is multiplied once by the per-axis coefficient `G_axis`** (§224-2). ~~A different coefficient for each term~~ was discarded. ⚠ **This arm does not do two things.** 1. **It does not pin a ΔR2 point estimate** (§12/§16). It receives prior samples through `dr2_draws` and puts the spread as it is into the predictive variance. 2. **It does not propagate the transfer coefficient by draws** — `G` is a scalar, so its variance contribution is **0** (§208-3). The width of `G` is reported as **the union (envelope) of the intervals computed at the two end values** (§212-3)."""

    name = "physics"
    requires = frozenset({"acquisitions", "dr2_draws", "site_labels"})

    def __init__(
        self,
        reference: str,
        target: str,
        *,
        axis: str = "ndi",
        tissue_dpar: float = 2.26,
        t1_in: float | None = None,
        t1_en: float | None = None,
        t1_iso: float | None = None,
        t1_tissue: float | None = None,
        nu: float = 0.575,
        t2_in: float = 82.0,
        t2_iso: float = 1000.0,
        gate_value: float | None = None,
        propagate_prior: bool = True,
    ) -> None:
        if axis not in G_TABLE:
            raise ValueError(f"축은 {sorted(G_TABLE)} 중 하나여야 한다: {axis!r}")
        if (t1_in is None) != (t1_en is None):
            raise ValueError("t1_in 과 t1_en 은 **둘 다** 주거나 둘 다 None 이어야 한다(§137)")
        if (t1_iso is None) != (t1_tissue is None):
            raise ValueError("t1_iso 와 t1_tissue 는 **둘 다** 주거나 둘 다 None 이어야 한다(§202-2)")
        if axis == "ndi" and t1_iso is not None:
            raise ValueError("NDI 축에 등방 T1 항을 넣지 않는다 — 다른 축이다(§202-2·§236)")
        if axis == "fiso" and t1_in is not None:
            raise ValueError("FISO 축에 조직 T1 항을 넣지 않는다 — 다른 축이다(§202-2·§236)")
        if not 0.0 < nu < 1.0:
            raise ValueError(f"ν 는 (0, 1) 이어야 한다: {nu}")
        g = gate(axis, tissue_dpar) if gate_value is None else float(gate_value)
        if not 0.0 <= g <= 1.0:
            raise ValueError(
                f"전달 계수는 [0, 1] 이어야 한다: {g}. "
                "⛔ **자르지 않고 여기서 멈춘다** — 조용히 clamp 하면 틀린 답이 경고 없이 나온다. "
                "`G_TABLE` 은 **무잡음 전제**의 값이고 그 범위는 [0.65, 0.97] 이다(§224-2). "
                "⛔ **잡음 조건에서는 `G > 1` 이 될 수 있다** — 실측 d∥ 2.0 · SNR 25 에서 "
                "**정본은 n=550,000 의 §250 표다** — d∥ 2.0 에서 m_tissue = 0.9529, "
                "95% 구간 [0.8453, **1.0605**] 로 **상한이 1 을 넘는다**. "
                "(§247-2 의 1.0565 는 n=100,000 값이고 §250 이 대체한다.) "
                "잡음 값을 쓰려면 **잡음 수준을 팔의 입력으로 받는 구조 변경**이 필요하다(QUEUE)."
            )
        self.reference, self.target, self.axis = reference, target, axis
        self.tissue_dpar, self.gate = tissue_dpar, g
        self.t1_in, self.t1_en = t1_in, t1_en
        self.t1_iso, self.t1_tissue = t1_iso, t1_tissue
        self.nu, self.t2_in, self.t2_iso = nu, t2_in, t2_iso
        self.propagate_prior = propagate_prior

    @staticmethod
    def _acq(acqs: dict, site: str) -> tuple[float, float]:
        """Accepts both `{site: (TE,TR)}` and `{site: [(TE,TR), ...]}`. The latter is the **session mean**."""
        v = acqs.get(site)
        if v is None:
            raise KeyError(f"획득 좌표에 '{site}' 가 없다. 있는 것: {sorted(acqs)}")
        arr = np.atleast_2d(np.asarray(v, dtype=float))
        if arr.shape[-1] != 2:
            raise ValueError(f"'{site}' 의 좌표는 (TE_ms, TR_s) 여야 한다: {v}")
        return float(arr[:, 0].mean()), float(arr[:, 1].mean())

    def _te_term(self, te_r: float, te_t: float, draws: np.ndarray) -> np.ndarray:
        """Target → reference TE term. **The formula differs for each axis** (§224-3)."""
        if self.axis == "ndi":
            return (te_r - te_t) * draws
        t2_en = 1.0 / (1.0 / self.t2_in + draws)          # ΔR2 differs with each draw
        ln_a = lambda te: np.log(
            self.nu * np.exp(-te / self.t2_in) + (1.0 - self.nu) * np.exp(-te / t2_en)
        )
        return (-te_r / self.t2_iso - ln_a(te_r)) - (-te_t / self.t2_iso - ln_a(te_t))

    def _tr_term(self, tr_r: float, tr_t: float) -> float:
        """Target → reference TR term. Uses a different T1 pair for each axis. If not given, 0 (§137)."""
        pair = (self.t1_in, self.t1_en) if self.axis == "ndi" else (self.t1_iso, self.t1_tissue)
        if pair[0] is None:
            return 0.0
        a, b = pair
        return (np.log(_rho(tr_r, a)) - np.log(_rho(tr_r, b))) - (
            np.log(_rho(tr_t, a)) - np.log(_rho(tr_t, b))
        )

    def shift_draws(self, view) -> np.ndarray:
        """**Deterministic samples** of the logit shift (the fixed ΔR2 grid). `fit` and (c) the original-unit variance use **the same samples**."""
        te_r, tr_r = self._acq(view.acquisitions, self.reference)
        te_t, tr_t = self._acq(view.acquisitions, self.target)
        draws = np.asarray(view.dr2_draws, dtype=float).ravel()
        if draws.size == 0:
            raise ValueError("dr2_draws 가 비었다 — 물리 팔은 사전 표본을 요구한다(§12/§16)")
        # **One statement**: add the two terms first, then multiply by the per-axis coefficient **only once** (§224-2).
        return self.gate * (self._te_term(te_r, te_t, draws) + self._tr_term(tr_r, tr_t))

    def fit(self, view) -> Correction:
        draws = np.asarray(view.dr2_draws, dtype=float).ravel()
        shift = self.shift_draws(view)

        c_hat = float(shift.mean())
        # **Zero-point uncertainty** = prior spread. It does not shrink with n — there are 0 traveling subjects, so there is no sample.
        sigma2_prior = float(shift.var(ddof=1)) if (self.propagate_prior and shift.size > 1) else 0.0
        return Correction(
            c_hat=c_hat,
            v=sigma2_prior,
            n_units=0,                  # §82 — 0 traveling subjects
            df=max(draws.size - 1, 1),
            sigma2_het=sigma2_prior,    # Floor that does not shrink with aggregation
            v_chat=0.0,
            t_scope="none",             # There is no sample, so no t multiplier is attached
        )

def native_var_from_draws(value_target, shift_draws) -> np.ndarray:
    """(c) **Variance in original units** over the physical-arm interval — the samples are passed through without any approximation. `Var_draws( expit( logit(v_j) + shift ) )` — the center is the calibrated value `expit(logit(v_j) + ĉ)`, and the same deterministic samples are passed through the inverse logit so that the variance is measured **directly**. ⛔ This is **the enforcement of the wording of §1-4-a**, not a new record («the comparison quantity is taken in original units after the inverse logit»). The delta method `(p(1−p))²·v` was an approximation of this, and at the calibrated center it was wrong **in both directions**, by NDI −5.56 to +4.93 % (author «approval four» 1 · `LADDER_ARMS` §2-ter)."""
    v = np.clip(np.asarray(value_target, dtype=float), 1e-6, 1 - 1e-6)
    x = np.log(v / (1 - v))[..., None] + np.asarray(shift_draws, dtype=float)
    return (1.0 / (1.0 + np.exp(-x))).var(axis=-1, ddof=1)


def loso_scanner(sites, fit_eval):
    """**LOSO-scanner** — fit with one scanner left out and evaluate on the scanner left out. Fixed in §169-5: **the main uncertainty is subject resampling**, and generalization at the scanner level is examined **only as a sensitivity analysis with this jackknife**. **A bootstrap over 7 scanners is unstable.** `fit_eval(train_sites, held_out) -> float` is called once per scanner. The return value is ``(values, mean, jackknife_se)``, and the SE is the **jackknife formula** ``√((n−1)/n · Σ(x_i − x̄)²)`` — not a bootstrap SD. When n is small, this one is the honest choice."""
    uniq = list(dict.fromkeys(sites))
    if len(uniq) < 3:
        raise ValueError(f"LOSO-scanner 는 스캐너 ≥3 을 요구한다: {uniq}")
    vals = np.array([fit_eval([s for s in uniq if s != out], out) for out in uniq], dtype=float)
    if not np.all(np.isfinite(vals)):
        bad = [uniq[i] for i in np.flatnonzero(~np.isfinite(vals))]
        raise ValueError(f"LOSO-scanner 결과가 유한하지 않다: {bad}")
    n = vals.size
    return vals, float(vals.mean()), float(np.sqrt((n - 1) / n * np.sum((vals - vals.mean()) ** 2)))
