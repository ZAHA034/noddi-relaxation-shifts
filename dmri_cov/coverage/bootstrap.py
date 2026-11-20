"""Cluster bootstrap that resamples one side only, because calibration uncertainty is already in the interval through k and must not count twice."""

from __future__ import annotations

from typing import Callable, Iterable, Sequence

import numpy as np

__all__ = ["cluster_ci", "m_curve_ci", "BootstrapResult"]


class BootstrapResult(dict):
    """`{"point": float, "lo": float, "hi": float, "se": float, "n_units": int, "method": str}`."""

    def __repr__(self) -> str:  # pragma: no cover - for display
        return (f"{self['point']:.4f} [{self['lo']:.4f}, {self['hi']:.4f}] "
                f"(se {self['se']:.4f}, {self['n_units']} units, {self['method']})")


def cluster_ci(values: np.ndarray, groups: np.ndarray, *,
               statistic: Callable[[np.ndarray], float] = np.mean,
               n_boot: int = 2000, seed: int = 20260921, ci: float = 0.95,
               side: str = "eval", resampled_sides: Sequence[str] = ("eval",)) -> BootstrapResult:
    """Subject cluster bootstrap. values : value per case (covered 0/1, width, vacuous …) groups : subject id per case — the **exchangeable unit** (§82) side : the side on which this statistic is computed. Must be contained in `resampled_sides`. resampled_sides : the sides on which resampling is allowed. **The length must be 1** (fixed in §143)."""
    sides = tuple(resampled_sides)
    if len(sides) != 1:
        raise ValueError(
            f"재표집은 한쪽에만 허용된다 — resampled_sides={sides}. "
            "적합 쪽까지 재표집하면 보정 불확실성을 두 번 센다(k 에 이미 들어 있다)."
        )
    if side not in sides:
        raise ValueError(f"side='{side}' 는 재표집 대상이 아니다(허용: {sides[0]})")

    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups)
    if values.shape[0] != groups.shape[0]:
        raise ValueError(f"values {values.shape} 와 groups {groups.shape} 의 길이가 다르다")
    keys = np.unique(groups)
    if keys.size < 2:
        raise ValueError(f"클러스터가 {keys.size} 개다 — 부트스트랩이 성립하지 않는다(§82)")

    idx_of = {k: np.flatnonzero(groups == k) for k in keys}
    rng = np.random.default_rng(seed)
    draws = np.empty(n_boot, dtype=float)
    for b in range(n_boot):
        pick = rng.integers(0, keys.size, keys.size)
        sel = np.concatenate([idx_of[keys[p]] for p in pick])
        draws[b] = float(statistic(values[sel]))
    lo_q, hi_q = (1.0 - ci) / 2.0, 1.0 - (1.0 - ci) / 2.0
    return BootstrapResult(point=float(statistic(values)),
                           lo=float(np.quantile(draws, lo_q)),
                           hi=float(np.quantile(draws, hi_q)),
                           se=float(draws.std(ddof=1)),
                           n_units=int(keys.size),
                           method=f"cluster bootstrap on '{side}' only")


def m_curve_ci(per_repeat: Iterable[float], *, ci: float = 0.95) -> BootstrapResult:
    """m-curve band — it is **a quantile over repeats**, not resampling. The m-curve of §123-4 **draws the cohort anew in every repeat**, so it is itself a marginalization. Applying a bootstrap on top would count the same spread twice. The quantiles across repeats are used as they are."""
    x = np.asarray(list(per_repeat), dtype=float)
    if x.size < 2:
        raise ValueError(f"반복이 {x.size} 회다 — 밴드를 만들 수 없다")
    lo_q, hi_q = (1.0 - ci) / 2.0, 1.0 - (1.0 - ci) / 2.0
    return BootstrapResult(point=float(x.mean()),
                           lo=float(np.quantile(x, lo_q)),
                           hi=float(np.quantile(x, hi_q)),
                           se=float(x.std(ddof=1) / np.sqrt(x.size)),
                           n_units=int(x.size),
                           method="m-curve repeat quantiles (no resampling)")


def band_is_informative(band: BootstrapResult, predicted_effect: float) -> tuple[bool, str]:
    """The reporting format fixed in §139-4 — it is meaningful only if **the band is narrower than the predicted deficit**. If the band is wider, passing does not mean "the formula was right". That fact is **returned together with the verdict**, so that it cannot be left out of the report."""
    width = band["hi"] - band["lo"]
    if width > abs(predicted_effect):
        return False, (f"검정력 부족 — 밴드 폭 {width:.4f} > 예측 결손 "
                       f"{abs(predicted_effect):.4f}. 통과해도 식이 맞았다는 뜻이 아니다.")
    return True, f"밴드 폭 {width:.4f} ≤ 예측 결손 {abs(predicted_effect):.4f}"
