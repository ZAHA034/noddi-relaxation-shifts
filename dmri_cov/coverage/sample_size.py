"""Computes the traveling-subject calibration sample size and harm threshold: calibration shifts the center but adds variance, so small samples can do worse."""
from __future__ import annotations

import math
from dataclasses import dataclass

__all__ = [
    "z_for_alpha",
    "correction_variance",
    "correction_variance_3level",
    "min_subjects",
    "k_corrected",
    "harm_threshold",
    "is_beneficial",
    "effective_units_3level",
    "SampleSizeResult",
]


def z_for_alpha(alpha: float = 0.10) -> float:
    """z of the two-sided (1−alpha) interval. alpha=0.10 → 1.6448536269514722 (nominal 0.90, §107)."""
    if not 0.0 < alpha < 1.0:
        raise ValueError(f"alpha 는 (0,1) 이어야 한다: {alpha}")
    # Φ⁻¹(1 − alpha/2) — without scipy (the coverage stratum uses only numpy/torch)
    return _ndtri(1.0 - alpha / 2.0)


def _ndtri(p: float) -> float:
    """Standard normal quantile. Acklam rational approximation + one Halley step (|err| < 1e-15)."""
    if not 0.0 < p < 1.0:
        raise ValueError(f"p 는 (0,1) 이어야 한다: {p}")
    a = (-3.969683028665376e01, 2.209460984245205e02, -2.759285104469687e02,
         1.383577518672690e02, -3.066479806614716e01, 2.506628277459239e00)
    b = (-5.447609879822406e01, 1.615858368580409e02, -1.556989798598866e02,
         6.680131188771972e01, -1.328068155288572e01)
    c = (-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e00,
         -2.549732539343734e00, 4.374664141464968e00, 2.938163982698783e00)
    d = (7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e00,
         3.754408661907416e00)
    plow, phigh = 0.02425, 1.0 - 0.02425
    if p < plow:
        q = math.sqrt(-2.0 * math.log(p))
        x = (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) / \
            ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0)
    elif p > phigh:
        q = math.sqrt(-2.0 * math.log(1.0 - p))
        x = -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) / \
             ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0)
    else:
        q, r = p - 0.5, (p - 0.5) ** 2
        x = (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q / \
            (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1.0)
    # Halley update
    e = 0.5 * math.erfc(-x / math.sqrt(2.0)) - p
    u = e * math.sqrt(2.0 * math.pi) * math.exp(x * x / 2.0)
    return x - u / (1.0 + x * u / 2.0)


def correction_variance(sd_paired: float, rho: float, n_subj: int) -> float:
    """Variance of the calibration shift ĉ, v = SD²_paired · ρ / n_subj (§71-bis-2). With ρ = 0, one should not take v = 0 but go to the **independent-sampling floor** SD²/n_total; here, however, only the hierarchical floor is treated — the term that cannot be lowered by adding voxels when ρ > 0. For ρ = 0, the same formula is used with the total number of observations in the place of `n_subj` (the rule of §82: n is the **number of independent units**)."""
    if sd_paired < 0.0:
        raise ValueError(f"sd_paired 는 음수일 수 없다: {sd_paired}")
    if not 0.0 <= rho <= 1.0:
        raise ValueError(f"rho 는 [0,1] 이어야 한다: {rho}")
    if n_subj < 1:
        raise ValueError(f"n_subj 는 1 이상이어야 한다: {n_subj}")
    return sd_paired * sd_paired * rho / float(n_subj)


def correction_variance_3level(
    sd_paired: float, rho_batch: float, rho_subj: float, n_batch: int, n_subj_per_batch: int
) -> float:
    """v at 3 levels (batch > subject > voxel) (LOCK_D §3-0b). Variance components: batch σ²_B = SD²·ρ_batch, subject σ²_S = SD²·(ρ_subj − ρ_batch), and the residual is taken to vanish when averaged over voxels (as measured in §71-bis-2). v = σ²_B / n_batch + σ²_S / (n_batch · n_subj_per_batch) With ρ_batch = 0 it reduces exactly to the 2-level formula (n_subj = n_batch · n_subj_per_batch)."""
    # The voxel residual term σ²_E/(n_batch·n_subj·n_vox) is omitted — §71-bis-2 is voxel 16-fold
    # Since it was measured that increasing it hardly lowers v, this is the n_vox → ∞ limit. For finite voxels
    # This formula **underestimates** v, so it underestimates n_min (not conservative). Needs verification.
    if not 0.0 <= rho_batch <= rho_subj <= 1.0:
        raise ValueError(f"0 ≤ rho_batch ≤ rho_subj ≤ 1 이어야 한다: {rho_batch}, {rho_subj}")
    if n_batch < 1 or n_subj_per_batch < 1:
        raise ValueError("n_batch·n_subj_per_batch 는 1 이상이어야 한다")
    s2 = sd_paired * sd_paired
    var_b = s2 * rho_batch
    var_s = s2 * (rho_subj - rho_batch)
    return var_b / n_batch + var_s / (n_batch * n_subj_per_batch)


def effective_units_3level(rho_batch: float, rho_subj: float, n_batch: int, n_subj_per_batch: int) -> float:
    """**Effective number of independent units** of the 3-level design — the n that gives the same v when put into the 2-level formula (§82). When the batch component dominates (ρ_batch → ρ_subj), the effective units **collapse to n_batch**. ON-Harmony has 2–6 batches per scanner (LOCK_D §3-0), so this collapse can overturn C4."""
    v = correction_variance_3level(1.0, rho_batch, rho_subj, n_batch, n_subj_per_batch)
    if v <= 0.0:
        return float(n_batch * n_subj_per_batch)
    return rho_subj / v


@dataclass(frozen=True)
class SampleSizeResult:
    """Return value of min_subjects. It carries both sides of the inequality (§89: prevents confusing the measured quantity with the claimed quantity)."""
    n_min: int              # Minimum number of subjects at which the adjustment pays off (rounded up to an integer)
    n_min_exact: float      # z²·ρ·SD²/b² (before rounding up)
    v_at_n_min: float       # Adjusted variance at that n
    bias: float             # |b|
    z: float
    beneficial_at_n_min: bool   # Self-verification — does it actually pay off at n_min
    form: str = "analytic"      # "exact" | "analytic" | "degenerate"
    beta: float = float("nan")  # b/s (exact form only)
    k_none: float = float("nan")

    def __str__(self) -> str:
        return (f"n_min={self.n_min} (exact {self.n_min_exact:.3f}, {self.form}) "
                f"|b|={self.bias:.6g} v={self.v_at_n_min:.6g} "
                f"{'✅' if self.beneficial_at_n_min else '❌'}")


def _k_none(beta: float, alpha: float = 0.10) -> float:
    """Width multiplier that gives nominal coverage without calibration. Because r = beta is fixed, a root search (§106)."""
    z = z_for_alpha(alpha)
    target = 1.0 - alpha

    def cov(k: float) -> float:
        return _Phi(z * k - beta) - _Phi(-z * k - beta)

    if cov(1.0) >= target:
        return 1.0
    lo, hi = 1.0, 2.0
    while cov(hi) < target and hi < 1e6:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if cov(mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _Phi(x: float) -> float:
    return 0.5 * math.erfc(-x / math.sqrt(2.0))


def k_corrected(variance: float, s_post: float) -> float:
    """**Closed form** of the width multiplier needed after calibration: k_corr = √(1 + v/s²) (§106). It is **exact** under the assumption that the residual after calibration is zero-mean Gaussian. If the posterior is not Gaussian it is an approximation, and since the posterior of ODI is bimodal (§41: 22.8 times at the fifth quantile), it is used only as an approximation (§107)."""
    if s_post <= 0.0:
        raise ValueError(f"s_post 는 양수여야 한다: {s_post}")
    return math.sqrt(1.0 + max(variance, 0.0) / (s_post * s_post))


def min_subjects(bias: float, sd_paired: float, rho: float, s_post: float | None = None,
                 alpha: float = 0.10) -> SampleSizeResult:
    """**Minimum number of traveling subjects** for which calibration beats no calibration (§106). Exact form (s_post given): n > SD²ρ / ( s²·(k_none²−1) ), β = b/s Analytic form (s_post=None): n ≈ ρ·SD²/b² (condition β ≲ 0.4) Parameters ---------- bias : the bias b that the calibration will cancel. Sign does not matter. sd_paired : SD of the **paired difference** (§82 — not the SD of a single measurement). rho : within-subject correlation ρ_trav of the paired difference. To be measured in ON-Harmony (LOCK_D §3-4). s_post : posterior SD. **If None, the analytic form** — the practical advantage of this tool is that using the rule does not need the precision of one's own estimator (see the module docstring). alpha : 1 − nominal coverage. Default 0.10 (§107, fixed in dMRI_PLAN.md:310 and confirmed). Notes ----- If b = 0, no n makes calibration worthwhile → n_min = ∞. This is the limit of Richter et al. 2022 "harmonization in the absence of a scanner effect may be harmful", and this inequality is the **quantitative version** of that sentence (§106 preceding batch)."""
    b = abs(float(bias))
    z = z_for_alpha(alpha)
    if b == 0.0:
        return SampleSizeResult(n_min=2**31 - 1, n_min_exact=float("inf"), v_at_n_min=0.0,
                                bias=0.0, z=z, beneficial_at_n_min=False, form="degenerate",
                                beta=0.0, k_none=float("inf"))
    if s_post is None:
        exact = rho * sd_paired * sd_paired / (b * b)
        form, beta, kn = "analytic", float("nan"), float("nan")
    else:
        beta = b / s_post
        kn = _k_none(beta, alpha)
        denom = s_post * s_post * (kn * kn - 1.0)
        if denom <= 0.0:
            exact = float("inf")
        else:
            exact = sd_paired * sd_paired * rho / denom
        form = "exact"
    n = max(1, int(math.floor(exact)) + 1) if math.isfinite(exact) else 2**31 - 1
    v = correction_variance(sd_paired, rho, n) if math.isfinite(exact) else 0.0
    ben = (k_corrected(v, s_post) < kn) if (s_post is not None and math.isfinite(exact)) else \
          (v < b * b)
    return SampleSizeResult(n_min=n, n_min_exact=exact, v_at_n_min=v, bias=b, z=z,
                            beneficial_at_n_min=ben, form=form, beta=beta, k_none=kn)


def harm_threshold(bias: float, sd_paired: float, rho: float, s_post: float | None = None,
                   alpha: float = 0.10) -> int:
    """**Harm threshold** — at or below this number of traveling subjects, calibration is worse than no calibration. It is the other side of the `min_subjects` inequality. It is an inequality, not an empirical observation (§106)."""
    r = min_subjects(bias, sd_paired, rho, s_post, alpha)
    return max(0, r.n_min - 1)


def is_beneficial(bias: float, variance: float, s_post: float, alpha: float = 0.10) -> bool:
    """Does calibration help: k_corr < k_none (§106 exact form). ⛔ `z·√v < |b|` of the version before 2026-09-21 was **wrong** — it treated the residual after calibration as a fixed bias and was over-conservative by a factor of z²=2.706. It was caught in the numerical check."""
    return k_corrected(variance, s_post) < _k_none(abs(float(bias)) / s_post, alpha)
