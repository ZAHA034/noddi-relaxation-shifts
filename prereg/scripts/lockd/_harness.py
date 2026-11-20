"""LOCK_D verdict harness: runs dummy data end to end before unsealing, freezing the contracts on input shape, missing cells and related rules."""
from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import pathlib
import sys
from typing import Any, Literal

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from dmri_cov.coverage.threshold import OPERATING_POINT, decision_line, dkk  # noqa: E402

# ─────────────────────────── Shape (§211 · §223-3) ───────────────────────────

# ON-Harmony: 20 subjects · 7 NODDI-capable scanners · 165 subject-sessions · 155 multi-shell
N_SUBJECT, N_SCANNER, N_SESSION = 20, 7, 165
# TBSS ROI unit — 27 from the 48-label JHU (§212-1 · §223-3)
N_ROI = 27
SCANNERS = ["OXF1PRI", "OXF2PRI", "OXF3TRI", "NOT1ACH", "NOT2ING", "NOT4GEP", "OXF4GEP"]
AXES = ("ndi", "odi", "fiso")

# ⛔ **The two cohorts are disjoint** (§211) — n is 20 / 10 / **0** depending on the pair.
# n=0 pairs: OXF3TRI↔NOT4GEP · OXF3TRI↔OXF4GEP. **Dropped from the verdict as «not measurable».**
EMPTY_PAIRS = {frozenset(("OXF3TRI", "NOT4GEP")), frozenset(("OXF3TRI", "OXF4GEP"))}

# ⛔ **Structural 0 — one group of a split target** (2026-10-02 · review-session reply «application of the same rule» — same handling as n=0 pairs).
# NOT2ING has two TE groups (@lo 89 · @hi 98 — LOCK_B_inputs §150 · `estimand.split_group`), and the four @lo of the §150 table
# (13192 · 13305 · 14229 · 14230) are **all cohort A**. The response of P2 is **cohort B fixed** (§248-1 ②), so
# The «@lo-only» version of the sensitivity analysis in §198-2 is **not defined**. It is not counted as 0 — **it is dropped by name and printed in the output.**
# Key = (verdict, split group, cohort) · value = reason. On real data `check_empty_splits(frame)` confirms this 0 again.
# ⛔ **Undefined resample** (2026-10-02 · relay #2 · `JUDGE_PATHS` LD-26) — in the cohort-stratified subject resample of ledger D
# If NOT2ING@lo (4 cohort A subjects) becomes 0 (about 0.6¹⁰ ≈ 0.006), γ · R are undefined. That resample is **dropped by name and
# the count is printed in the output** — it is not passed on as 0 or NaN (same handling as n = 0 pairs · split groups). The verdict code (`judge_D`) enforces it.
UNDEFINED_RESAMPLES: dict[tuple[str, str], str] = {
    ("원장 D", "NOT2ING@lo"): "코호트 층화 재표본에서 @lo 넷(전부 코호트 A)이 하나도 뽑히지 않으면 그 타깃의 γ 가 없다 — 재표본째 뺀다",
}

EMPTY_SPLITS: dict[tuple[str, str, str], str] = {
    ("P2", "NOT2ING@lo", "B"): "§150 표 @lo 넷 전부 코호트 A — §198-2 «@lo-only» 판 정의 안 됨",
}
# Cohort B = visitors to NOT4GEP · OXF4GEP (§211-1). A subject with such a row is regarded as B.
COHORT_B_MARKER = "NOT4GEP"


def structural_zero_lines() -> list[str]:
    """Structural zeros **by name** — the n=0 pairs (§211) and the separated group (§150). Printed in every output («they are not silently counted as 0»)."""
    order = {s: i for i, s in enumerate(SCANNERS)}
    pairs = [f"{a}↔{b} (n=0, §211)" for a, b in
             sorted((tuple(sorted(p, key=order.get)) for p in EMPTY_PAIRS), key=lambda t: (order[t[0]], order[t[1]]))]
    splits = [f"{pred} {grp} (코호트 {coh} n=0 — {why})" for (pred, grp, coh), why in sorted(EMPTY_SPLITS.items())]
    return pairs + splits


def check_empty_splits(frame: "Frame") -> None:
    """Whether the structural zeros in the real data are really 0 — if not, **it stops** (a version is not silently made)."""
    if frame.acq is None:
        raise ValueError("⛔ acq 가 없어 NOT2ING 군을 가를 수 없다 — 구조적 0 을 확인할 수 없다")
    from dmri_cov.coverage.estimand import split_group
    cohort_b = set(frame.subject[frame.scanner == COHORT_B_MARKER].tolist())
    for (pred, grp, coh), why in EMPTY_SPLITS.items():
        name = grp.split("@")[0]
        hit = [s for s, sc, te in zip(frame.subject.tolist(), frame.scanner.tolist(), frame.acq[:, 0].tolist())
               if sc == name and s in cohort_b and split_group(name, float(te)) == grp]
        if hit:
            raise ValueError(f"⛔ 구조적 0 이 깨졌다: {pred} {grp} — 코호트 {coh} 에 {len(hit)} 행({why}). 판을 만들지 않고 멈춘다")

MissingPolicy = Literal["drop_pair", "drop_roi", "propagate_nan"]

# Names of the 27 ROI units — from the 48-label XML, the left-right tags are stripped and folded (same rule as §212-1).
_ATLAS_XML = ROOT_ATLAS = pathlib.Path(__file__).resolve().parents[3] / "prereg/artifacts/atlas/JHU-labels.xml"


def roi_names() -> list[str]:
    """The 27 ROI names of the second axis of `value`. If there is no atlas, it falls back to `ROI00…`."""
    import re
    if not _ATLAS_XML.is_file():
        return [f"ROI{i:02d}" for i in range(N_ROI)]
    labels = {int(m.group(1)): m.group(2).strip() for m in re.finditer(
        r'<label index="(\d+)"[^>]*>([^<]+)</label>', _ATLAS_XML.read_text(errors="replace"))}
    base, seen = [], set()
    for k in sorted(labels):
        if k == 0:
            continue
        nm = re.sub(r"\s+[LR]$", "", labels[k]).strip()
        if nm not in seen:
            seen.add(nm); base.append(nm)
    return base if len(base) == N_ROI else [f"ROI{i:02d}" for i in range(N_ROI)]


@dataclasses.dataclass(frozen=True)
class Frame:
    """The **only shape** of a verdict input. Both real data and dummies must have this form."""

    subject: np.ndarray        # (n_row,)  subject ID
    scanner: np.ndarray        # (n_row,)  scanner name
    session: np.ndarray        # (n_row,)  session ID
    value: np.ndarray          # (n_row, N_ROI, 3)  ROI × (ndi, odi, fiso)
    source: str                # "dummy" | actual file path
    sha256: str | None         # For real data, the file hash — kept in the R10 log
    # ⛔ **[contract revision 2026-09-28 · §260 revision record · P9-8-75/76/77]**
    # Before, there was only the raw IDP `(n_row, N_ROI, 3)`. So there was **no material to build the physical prediction `ĉ`**
    # was not in the contract**, and `judge_B17` filled the place with random numbers, `judge_B7` with the observation, `judge_B19p` with the overall mean
    # — **none of the three was a recorded quantity, and all directions were «relaxing»**.
    # ⇒ **This revision is an «enforcement», not a recorded «change».** The recorded quantity was the physical predicted shift from the start.
    # Source: **TE·TR = sidecar JSON (R10 exception 1, acquisition metadata — prior access allowed)**.
    # is neither a verdict value nor a new measurement, so it **is not entered in the reading record** (READ_LOG §1).
    acq: np.ndarray | None = None      # (n_row, 2) = (TE_ms, TR_s) · per session

    def __post_init__(self) -> None:
        n = self.subject.size
        for name, a, shape in (("scanner", self.scanner, (n,)), ("session", self.session, (n,)),
                               ("value", self.value, (n, N_ROI, 3))):
            if a.shape != shape:
                raise ValueError(f"{name} 의 형상이 {a.shape} 다 — {shape} 여야 한다")
        bad = set(np.unique(self.scanner)) - set(SCANNERS)
        if bad:
            raise ValueError(f"등록되지 않은 스캐너: {sorted(bad)}. 등록본: {SCANNERS}")
        if self.acq is not None and self.acq.shape != (n, 2):
            raise ValueError(f"acq 의 형상이 {self.acq.shape} 다 — {(n, 2)} 여야 한다 "
                             "((TE_ms, TR_s) 를 **행마다**)")

    @property
    def n_row(self) -> int:
        return int(self.subject.size)

    def axis(self, name: str) -> np.ndarray:
        """(n_row, N_ROI) — takes out one axis only."""
        return self.value[:, :, AXES.index(name)]

    def acquisitions(self) -> dict[str, list[tuple[float, float]]]:
        """Scanner → list of **per-session** (TE, TR). ⛔ **No session representative value is made** — `en/OSF:189` «per-session sidecar TE, **not a scanner representative value**». `PhysicsArm._acq` receives the list and uses the session mean. ★ **On the NDI axis, that mean is identically equal to "the mean of the per-session ĉ"** (linear in TE, P9-8-78) — so there is no choice here. ⛔ On the FISO axis it is not equal."""
        if self.acq is None:
            raise ValueError(
                "⛔ 이 Frame 에 **획득 좌표가 없다**. 물리 예측 `ĉ` 를 만들 수 없다.\n"
                "   B7 · B17 · B19′ 는 **물리 예측 이동**을 요구한다(P9-8-75).\n"
                "   → 실데이터 로더가 사이드카 TE·TR 을 함께 실어야 한다(R10 예외 1).")
        out: dict[str, list[tuple[float, float]]] = {}
        for i, sc in enumerate(self.scanner):
            te, tr = float(self.acq[i, 0]), float(self.acq[i, 1])
            if np.isfinite(te) and np.isfinite(tr):
                out.setdefault(str(sc), []).append((te, tr))
        return out


# ⛔ **Deterministic** sample of the ΔR2 prior — a quantile grid. **No random numbers are used.**
# `JUDGE_PATHS` §1-2 **LD-9** has already recorded: «Because the NDI TE term is **linear** in the draw,
# `ĉ = G·ΔTE·E[ΔR2]`, **computed analytically** with the truncated mean. **No seed needed**».
# ⇒ there are no random numbers on the verdict input path — the requirement of `check_prediction_path.py` is satisfied **by design**.
# The number of grid points is **recorded** here (if changed, update it together with `EXPECTED_COUNTS`).
DR2_GRID_N = 2001

# ⛔ **T1 of the FISO-axis TR term — recorded value** (2026-10-02 · relay #2 · before stage 30).
# OSF EN :347 rule 3 «**It must be included when ΔR2 is small, and on the FISO axis**» · :201 «**The FISO-axis negative
# control includes the TR term**» · the value is the table of LOCK_B_inputs_EN §209 «`T1_iso` 4.0 s · `T1_tissue` 1.0 s» (the value used by the §224 calculation).
# The previous version did not pass the T1 pair, so `PhysicsArm._tr_term` gave **0** — the predicted shift of the negative control (TE matched · TR 1.81 s difference) was
# Sticking to 0 at −0.0032, it matched «if TE is the mechanism, it is near 0» automatically (favorable). This version: +0.0409. ⇒ Direction: favorable → eliminated.
# The NDI axis is not changed — the recorded calculation of ledger D is «TR term = 0» (OSF EN :213 · `JUDGE_PATHS` LD-8).
FISO_T1 = {"t1_iso": 4.0, "t1_tissue": 1.0}

# ⛔ **FISO ν envelope — enforcement of the recorded reporting convention** (2026-10-02 · relay #4 ■4 · before stage 30 · after `_harness.py` was frozen — re-frozen).
# Recorded (attached `LOCK_B_inputs_EN@8c80a1c` §249-1): «when writing a FISO shift, **write ν with it.**» · «The primary report is
# **ν = 0.575 · d∥ 2.26**, **with the envelope at the support interval's ends written alongside**». Support ends = ν «[0.40, 0.85]» (item 2 of the same section)
# × d∥ `arms.DPAR_SUPPORT` (2.0, 2.5) (§209-2). **Not entered into the verdict value — report only.**
# The direction of the gap was written **before** fixing it (`CODE_TERMS_SWEEP` 2-bis · a0aea74): **favorable** — without the envelope the FISO prediction looks more precise than it is.
FISO_NU_MAIN = 0.575
FISO_NU_SUPPORT = (0.40, 0.85)


def dr2_draws() -> np.ndarray:
    """**Deterministic** samples from the fixed ΔR2 prior (ppf grid). The mean agrees with the fixed truncated mean to within 1e-6."""
    from dmri_cov.eval.t2_oracle import DR2_PRIOR
    q = (np.arange(1, DR2_GRID_N + 1) - 0.5) / DR2_GRID_N
    return np.asarray([DR2_PRIOR.ppf(float(x)) for x in q], dtype=float)


def _physics_arm(reference: str, target: str, axis: str, tissue_dpar: float, nu: float = FISO_NU_MAIN):
    from dmri_cov.coverage.arms import PhysicsArm
    t1 = FISO_T1 if axis == "fiso" else {}
    return PhysicsArm(reference=reference, target=target, axis=axis, tissue_dpar=tissue_dpar, nu=nu, **t1)


def physics_by_subject(frame: "Frame", reference: str, target: str, *,
                       axis: str = "ndi", tissue_dpar: float = 2.26,
                       nu: float = FISO_NU_MAIN) -> dict[str, tuple[float, np.ndarray]]:
    """`(ĉ, logit-shift samples)` for each common subject — computed from the sidecar TE · TR of that subject's **first session**. ⛔ [2026-10-02 · relay #3 (4) · before step 30] Fixed: OSF EN :189 «the physics transfer uses the **per-session sidecar TE**, not a scanner representative value» · `JUDGE_PATHS` LD-10 · session choice is §1-5-a (first session). The earlier version used one session-mean TE · TR for each scanner — the per-subject TE agreement of the negative control (fixed «per-subject |ΔTE| 0.04 to 1.54») was lost and the @lo · @hi of NOT2ING were mixed into a mean (`CODE_TERMS_SWEEP` 2-bis — before the fix, the direction was written as "unknown · the more the physics is right, the worse for it")."""
    import types
    if frame.acq is None:
        raise ValueError("⛔ 획득 좌표가 없다 — 세션별 TE · TR 로 ĉ 를 낼 수 없다(R10 예외 1 · 사이드카)")
    cells = first_session_cells(frame)
    acq = {(str(s), str(c)): (float(a[0]), float(a[1])) for s, c, a in zip(cells.subject, cells.scanner, cells.acq)}
    subs = sorted({s for s, c in acq if c == reference} & {s for s, c in acq if c == target})
    draws = dr2_draws()
    arm = _physics_arm(reference, target, axis, tissue_dpar, nu)
    out: dict[str, tuple[float, np.ndarray]] = {}
    for s in subs:
        view = types.SimpleNamespace(acquisitions={reference: acq[(s, reference)], target: acq[(s, target)]},
                                     dr2_draws=draws, site_labels=np.asarray(SCANNERS))
        d = arm.shift_draws(view)
        out[s] = (float(d.mean()), d)
    return out


def _physics_representative(frame: "Frame", reference: str, target: str, axis: str, tissue_dpar: float,
                            nu: float = FISO_NU_MAIN) -> np.ndarray:
    """Only for pairs with **0** common subjects (`EMPTY_PAIRS`) — session-mean TE · TR for each scanner («Representative values are for aggregation and reporting only» — OSF EN :189). Used only to **report** the point calibration shift (G3 — for a pair with 0 subjects, coverage and width are undefined)."""
    import types
    view = types.SimpleNamespace(acquisitions=frame.acquisitions(), dr2_draws=dr2_draws(), site_labels=np.asarray(SCANNERS))
    return _physics_arm(reference, target, axis, tissue_dpar, nu).shift_draws(view)


def physics_chat(frame: "Frame", reference: str, target: str, *,
                 axis: str = "ndi", tissue_dpar: float = 2.26, nu: float = FISO_NU_MAIN) -> float:
    """The fixed **physical transfer** `ĉ = G_axis(d∥) · (TE term + TR term)` (§224-2) — the value of one pair = **the mean of the per-common-subject ĉ**. ⛔ **This is the "predicted shift" that B7 · B17 · B19′ must use.** It is neither the observation (`va − vb`), nor a random number (`N(observation, 0.01)`), nor the overall mean (`nanmean(ndi)`) — none of the three was a fixed quantity, and **all went in the direction of "easing"** (P9-8-75 · P9-8-76). ⛔ `tissue_dpar` is accepted **only at the four points on the fixed support** — `arms.gate()` **rejects** anything outside the grid with `ValueError` (as the no-interpolation clause says). Both ends are given together as the §212-3 envelope. ⛔ [10-02] It is the mean of the per-subject ĉ (first-session TE · TR) (`physics_by_subject`) — the representative value is "for aggregation and reporting only" (OSF EN :189). Only pairs with 0 common subjects use the scanner representative value (`_physics_representative` — the point calibration shift for reporting)."""
    by = physics_by_subject(frame, reference, target, axis=axis, tissue_dpar=tissue_dpar, nu=nu)
    if by:
        return float(np.mean([c for c, _ in by.values()]))
    return float(_physics_representative(frame, reference, target, axis, tissue_dpar, nu).mean())


def physics_draws(frame: "Frame", reference: str, target: str, *,
                  axis: str = "ndi", tissue_dpar: float = 2.26, nu: float = FISO_NU_MAIN) -> np.ndarray:
    """Samples of the logit shift under the **same rule** as `physics_chat` — the element-wise mean of the common-subject samples (a representative value for a pair with 0 subjects). For the interval of each subject, the samples of `physics_by_subject` are used (`ladder30`)."""
    by = physics_by_subject(frame, reference, target, axis=axis, tissue_dpar=tissue_dpar, nu=nu)
    if by:
        return np.mean([d for _, d in by.values()], axis=0)
    return _physics_representative(frame, reference, target, axis, tissue_dpar, nu)


def physics_nu_envelope(frame: "Frame", reference: str, target: str, *, tissue_dpar: float = 2.26) -> dict:
    """What is **carried along** when the FISO shift is used — the main value (ν 0.575 · d∥ 2.26), the `ĉ` at the four corners of the support (ν × d∥), and their envelope (minimum · maximum). It is produced by the same rule: prior mean of the TE · TR · ΔR2 of the first session of each common subject (`physics_chat`). The envelope in the fixed text (§249-1 «−0.0347 to −0.1032») is a value at the fixed operating point measured with the **point** ΔR2, so it is not the same quantity as the value on this path (prior mean · actual pairs) (`CODE_TERMS_SWEEP` 2-quater · 2-quinquies)."""
    from dmri_cov.coverage.arms import DPAR_SUPPORT
    corners = {(nu, dp): physics_chat(frame, reference, target, axis="fiso", tissue_dpar=dp, nu=nu)
               for nu in FISO_NU_SUPPORT for dp in DPAR_SUPPORT}
    v = list(corners.values())
    return {"ν": FISO_NU_MAIN, "d∥": tissue_dpar,
            "ĉ": physics_chat(frame, reference, target, axis="fiso", tissue_dpar=tissue_dpar, nu=FISO_NU_MAIN),
            "포락선": [min(v), max(v)],
            "모서리": {f"ν {nu} · d∥ {dp}": c for (nu, dp), c in corners.items()},
            "규약": "LOCK_B_inputs_EN@8c80a1c §249-1 «write ν with it» — 지지 끝 ν [0.40, 0.85] × d∥ [2.0, 2.5] · 보고 전용"}


def s_ref_from_frame(frame: "Frame", reference: str = "OXF1PRI") -> np.ndarray:
    """`s_ref` for the recovery · synthetic frame — **the same rule** as `compute_s_ref.py` (square root of the mean per-subject session variance of the reference repeats). ⛔ Not used on real data — the real-data `s_ref` is read only from the frozen file (`s_ref.frozen.json`, commit `5262f43`). Returns: `(3, N_ROI)` — in `AXES` order."""
    m = frame.scanner == reference
    subs = [s for s in np.unique(frame.subject[m]) if int((m & (frame.subject == s)).sum()) > 1]
    if not subs:
        raise ValueError("참조 반복 세션이 없다 — s_ref 를 만들 수 없다")
    var = np.mean([frame.value[m & (frame.subject == s)].var(axis=0, ddof=1) for s in subs], axis=0)
    return np.sqrt(var).T


def s_ref_frozen() -> np.ndarray:
    """Frozen real-data `s_ref` — `(3, N_ROI)`, ROI order = `roi_names()`."""
    doc = json.loads((pathlib.Path(__file__).resolve().parents[3] / "prereg" / "artifacts" / "unseal"
                      / "s_ref.frozen.json").read_text(encoding="utf-8"))
    order = roi_names()
    return np.array([[doc["ROI별 s_ref"][ax][u] for u in order] for ax in AXES])


def dummy_frame(seed: int = 20260924, missing_frac: float = 0.03) -> Frame:
    """A synthetic input that **has only the same shape**. ⛔ The values have no meaning — only the wiring is checked. Missing values are put in on purpose (`missing_frac`) — **if the missing-value path does not run, it is caught here.** The structure in which the two cohorts are disjoint (§211) is also reproduced."""
    rng = np.random.default_rng(seed)
    coh_a, coh_b = SCANNERS[:4], SCANNERS[3:]          # No. 3 overlaps, making 6 pairs with n=20
    rows_s, rows_c = [], []
    for subj in range(N_SUBJECT):
        pool = coh_a if subj < N_SUBJECT // 2 else coh_b
        for sc in pool:
            rows_s.append(subj)
            rows_c.append(sc)
    subject = np.array(rows_s)
    scanner = np.array(rows_c)
    session = np.arange(subject.size)
    value = rng.uniform(0.05, 0.85, size=(subject.size, N_ROI, 3))
    mask = rng.random(value.shape) < missing_frac
    value[mask] = np.nan
    # ⛔ The acquisition coordinates also synthesize **only the shape**. The values have no meaning — only the wiring is checked.
    # ⛔ Do not write the recorded TE here: if a dummy imitates a real number, «judged with a dummy» can no longer be
    # told apart by eye (P9-8-75 was exactly that case).
    te = 60.0 + 40.0 * rng.random(subject.size)        # Synthetic TE (ms)
    tr = 3.0 + 0.8 * rng.random(subject.size)          # Synthetic TR (s)
    acq = np.column_stack([te, tr])
    return Frame(subject=subject, scanner=scanner, session=session, value=value,
                 source="dummy", sha256=None, acq=acq)


# ⛔ **Recorded traveling structure** (`LADDER_ARMS.md` §1-4-b table · §1-5-a) — the retrieval frame and the 30 runner **assert** it.
# The 4 common scanners have all 20 subjects; cohort A (10 subjects) has OXF3TRI, cohort B (10 subjects) has OXF4GEP·NOT4GEP.
# ⇒ visited cells 110 · pairs {20 subjects: 6, 10 subjects: 13, 0 subjects: 2} · 9 repeat subjects × 6 sessions ⇒ 155 sessions.
# ⛔ `dummy_frame` does **not reproduce** this structure ({10: 12, 0: 9}, H-16) — it is for wiring only.
COMMON4 = ("OXF1PRI", "OXF2PRI", "NOT1ACH", "NOT2ING")
COHORT_EXTRA = {"A": ("OXF3TRI",), "B": ("OXF4GEP", "NOT4GEP")}
REGISTERED_CELLS, REGISTERED_SESSIONS = 110, 155
REGISTERED_PAIR_STRATA = {20: 6, 10: 13, 0: 2}
# Scanner layout of the 9 repeat subjects (measured in §141-5 — NOT3GEM was dropped from the record).
REPEAT_SCANNERS = ("OXF1PRI", "OXF1PRI", "OXF2PRI", "OXF2PRI", "NOT1ACH", "NOT2ING",
                   "OXF3TRI", "OXF4GEP", "NOT4GEP")


def first_session_cells(frame: "Frame") -> "Frame":
    """§1-5-a — for each visit cell (subject, scanner), **the one with the smallest session number**. The choice does not depend on the values."""
    keep: dict[tuple, int] = {}
    for i in range(frame.n_row):
        key = (str(frame.subject[i]), str(frame.scanner[i]))
        if key not in keep or frame.session[i] < frame.session[keep[key]]:
            keep[key] = i
    idx = np.array(sorted(keep.values()))
    return Frame(subject=frame.subject[idx], scanner=frame.scanner[idx],
                 session=frame.session[idx], value=frame.value[idx], source=frame.source,
                 sha256=frame.sha256, acq=None if frame.acq is None else frame.acq[idx])


def pair_strata(frame: "Frame") -> dict[int, int]:
    """For each scanner pair in the visit cells, counts the **number of shared subjects** and gives {number of people: number of pairs} (all C(7,2) pairs)."""
    cells = first_session_cells(frame)
    who = {sc: set(map(str, cells.subject[cells.scanner == sc])) for sc in SCANNERS}
    out: dict[int, int] = {}
    for a in range(len(SCANNERS)):
        for b in range(a + 1, len(SCANNERS)):
            n = len(who[SCANNERS[a]] & who[SCANNERS[b]])
            out[n] = out.get(n, 0) + 1
    return dict(sorted(out.items(), reverse=True))


def recovery_frame(seed: int = 20260929) -> tuple["Frame", dict]:
    """**Recovery frame** — a synthetic input that keeps the fixed structure and **knows the true values** (for testing the step-30 runner). ⛔ The values are synthetic. TE is synthetic too (it does not imitate the fixed TE — this prevents dummy and real data from becoming impossible to tell apart by eye, P9-8-75). There is one use: **whether an arm recovers the known true value**. True structure (logit): `logit v = base_ROI + subject + scanner effect + subject×scanner heterogeneity + session noise`, and the scanner effect is **exactly the fixed physical transfer** — `e_OXF1PRI = 0`, `e_t = −ĉ(OXF1PRI → t)`. ⇒ The expected value of the true paired shift `logit v_i − logit v_j` of the pair (i, j) is `ĉ(i → j)` (linear TE term)."""
    rng = np.random.default_rng(seed)
    subs = [f"sub-{k:02d}" for k in range(20)]
    plan = []                                            # (subject, scanner, session)
    for k, s in enumerate(subs):
        for sc in COMMON4 + COHORT_EXTRA["A" if k < 10 else "B"]:
            plan.append((s, sc, 1))
    # Repeats: for each scanner, chosen deterministically from the subjects who visited that scanner
    used: set[str] = set()
    for sc in REPEAT_SCANNERS:
        pool = [s for s, c, _ in plan if c == sc and s not in used]
        s = pool[0]
        used.add(s)
        plan += [(s, sc, n) for n in range(2, 7)]
    subject = np.array([p[0] for p in plan]); scanner = np.array([p[1] for p in plan])
    session = np.array([p[2] for p in plan])
    te_sc = {sc: float(t) for sc, t in zip(SCANNERS, 64.0 + 35.0 * rng.random(len(SCANNERS)))}
    tr_sc = {sc: 3.6 for sc in SCANNERS}                 # The TR term is 0 on the verdict path
    acq = np.array([[te_sc[c], tr_sc[c]] for c in scanner])
    probe = Frame(subject=subject, scanner=scanner, session=session,
                  value=np.full((subject.size, N_ROI, 3), 0.5), source="dummy", sha256=None, acq=acq)
    lvl = {"ndi": np.linspace(0.40, 0.75, N_ROI), "odi": np.linspace(0.15, 0.45, N_ROI),
           "fiso": np.linspace(0.03, 0.12, N_ROI)}
    eff = {ax: {sc: (0.0 if sc == "OXF1PRI" else
                     -physics_chat(probe, "OXF1PRI", sc, axis=ax) if ax != "odi" else 0.0)
                for sc in SCANNERS} for ax in AXES}
    subj_eff = {s: rng.normal(0.0, 0.15) for s in subs}
    het = {(s, sc): rng.normal(0.0, 0.02, N_ROI) for s in subs for sc in SCANNERS}
    value = np.empty((subject.size, N_ROI, 3))
    for i in range(subject.size):
        for a, ax in enumerate(AXES):
            u = (np.log(lvl[ax] / (1 - lvl[ax])) + subj_eff[subject[i]] + eff[ax][scanner[i]]
                 + het[(subject[i], scanner[i])] + rng.normal(0.0, 0.03, N_ROI))
            value[i, :, a] = 1.0 / (1.0 + np.exp(-u))
    f = Frame(subject=subject, scanner=scanner, session=session, value=value,
              source="dummy", sha256=None, acq=acq)
    return f, {"scanner_effect_logit": eff, "te": te_sc}


# Column order of the ROI table (`roi_extract.sh`) — `fintra · fiso · OD` → `Frame.AXES` = (ndi, odi, fiso)
_ROI_TSV_COL = {"ndi": 0, "odi": 2, "fiso": 1}
# Recorded exclusion — single shell, so NODDI is impossible (LOCK_B §9 · `npe.EXCLUDE_SCANNERS`)
_EXCLUDED_SCANNERS = ("NOT3GEM",)


def _label_units() -> list[str]:
    """Labels 1..48 → names of 27 units (left/right tags dropped). **The same rule** as `roi_names()` · `compute_roi_s.py`."""
    import re
    labels = {int(m.group(1)): m.group(2).strip() for m in re.finditer(
        r'<label index="(\d+)"[^>]*>([^<]+)</label>', _ATLAS_XML.read_text(errors="replace"))}
    return [re.sub(r"\s+[LR]$", "", labels[v]).strip() for v in range(1, 49)]


# Log that leaves a line each time a real-data `Frame` is read — evidence of «when it was opened» (R25: the state comes from the output).
FRAME_READ_LOG = pathlib.Path(__file__).resolve().parents[3] / "prereg" / "artifacts" / "unseal" / "frame_reads.jsonl"


def _require_open() -> None:
    """Checks from the **committed gate ledger** whether step 30 has been opened — the call arguments are not trusted (R25)."""
    import unseal_gate as gate
    try:
        gate.require("30_scanner_cmp")
    except gate.GateError as e:
        raise PermissionError(f"실데이터는 **봉인 해제 후**(30 단계가 열린 뒤)에만 읽는다 — 게이트가 막았다(R10 · R25):\n{e}") from e


def _log_read(entry: dict) -> None:
    FRAME_READ_LOG.parent.mkdir(parents=True, exist_ok=True)
    with FRAME_READ_LOG.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(entry, ensure_ascii=False) + "\n")


def load_frame(path: str, *, utc_iso: str | None = None, meta_root: str | None = None,
               manifest: list | None = None) -> Frame:
    """Real-data `Frame` — **called only at step 30** (R10). The implementation is before the seal; the call is after the seal is released. Input: `path` — output directory of `roi_extract.sh` (`<sub>_<ses>.roi.tsv`, 48 labels × `fintra·fiso·OD`) `meta_root` — root of the ON-Harmony acquisition metadata (default: `onharmony_meta` of `dmri_cov.utils.ext`). Reads `EchoTime` (s → ms) · `RepetitionTime` (s) from the `dir-AP` dwi sidecar of the session (R10 exception 1). Rules: · Scanner = first 7 characters of the session name (`ses-OXF1PRI001` → `OXF1PRI`). `NOT3GEM` is a **fixed exclusion**, so it is not loaded and is left in `manifest` as «제외» (excluded). · 48 labels → 27 units by a **simple left-right average** (the same rule as `compute_roi_s.py`, order from `roi_names()`). · **All sessions** are loaded — the first session per cell is chosen by `first_session_cells` (§1-5-a). · NOT2ING does not split the scanner name — sessions with a different TE (@lo/@hi) are told apart by the TE of `acq`. · If the sidecar is missing, or `EchoTime`·`RepetitionTime` is missing, it **fails** (nothing is silently dropped). ⛔ **The control is not an argument but an output** (R25) — if the committed gate ledger has not opened step 30, nothing is read. ⛔ On every read, `FRAME_READ_LOG` gets a **"start" line before the read** and a **"done" line after the read** (digest · row count). The time `utc_iso` is given by the caller (the script does not read the clock — a gate convention). If absent, nothing is read. `(relative path, sha256)` is appended to `manifest` — the caller copies it into the access log (OSF §8-1 item 2)."""
    import hashlib
    import json
    import re

    _require_open()
    if not utc_iso:
        raise ValueError("`utc_iso`(열람 시각, UTC)를 밖에서 줘야 한다 — 로그에 남는다")
    _log_read({"event": "시작", "utc": utc_iso, "source": str(path)})
    from dmri_cov.utils import ext
    root = pathlib.Path(path)
    meta = pathlib.Path(meta_root) if meta_root else ext.component_root("onharmony_meta")
    units = _label_units()
    order = roi_names()
    if len(order) != N_ROI or set(order) != set(units):
        raise ValueError("아틀라스 단위 27 개를 만들 수 없다 — JHU-labels.xml 을 확인하라")
    log = manifest if manifest is not None else []

    def _sha(f: pathlib.Path) -> str:
        return hashlib.sha256(f.read_bytes()).hexdigest()

    subj, scan, sess, vals, acq = [], [], [], [], []
    files = sorted(root.glob("*.roi.tsv"))
    if not files:
        raise FileNotFoundError(f"ROI 표가 없다: {root}")
    for f in files:
        m = re.fullmatch(r"(sub-[^_]+)_ses-([A-Z0-9]+)\.roi\.tsv", f.name)
        if not m:
            raise ValueError(f"파일 이름 형식이 아니다: {f.name}")
        sub, ses = m.group(1), m.group(2)
        sc = ses[:7]
        if sc in _EXCLUDED_SCANNERS:
            log.append((f.name, "제외(등록된 제외 스캐너)"))
            continue
        rows = [[float(x) for x in ln.split()] for ln in f.read_text().splitlines() if ln.strip()]
        if len(rows) != 48 or any(len(r) != 3 for r in rows):
            raise ValueError(f"{f.name}: 48 행 × 3 열이어야 한다")
        v = np.empty((N_ROI, 3))
        for a, ax in enumerate(AXES):
            acc: dict[str, list[float]] = {}
            for lab in range(48):
                acc.setdefault(units[lab], []).append(rows[lab][_ROI_TSV_COL[ax]])
            v[:, a] = [float(np.mean(acc[u])) for u in order]
        js = meta / sub / f"ses-{ses}" / "dwi" / f"{sub}_ses-{ses}_dir-AP_dwi.json"
        if not js.is_file():
            raise FileNotFoundError(f"사이드카가 없다: {js}")
        d = json.loads(js.read_text())
        if d.get("EchoTime") is None or d.get("RepetitionTime") is None:
            raise ValueError(f"{js.name}: EchoTime·RepetitionTime 이 없다")
        log += [(f.name, _sha(f)), (str(js.relative_to(meta)), _sha(js))]
        subj.append(sub); scan.append(sc); sess.append(ses); vals.append(v)
        acq.append((float(d["EchoTime"]) * 1000.0, float(d["RepetitionTime"])))
    digest = hashlib.sha256("\n".join(f"{a} {b}" for a, b in log).encode()).hexdigest()
    _log_read({"event": "완료", "utc": utc_iso, "source": str(path), "rows": len(subj),
               "files": len(log), "sha256": digest})
    return Frame(subject=np.array(subj), scanner=np.array(scan), session=np.array(sess),
                 value=np.stack(vals), source=str(root), sha256=digest, acq=np.array(acq))


# ─────────────────────────── Output · verdict line ───────────────────────────

@dataclasses.dataclass
class Verdict:
    """The **only output format** of a verdict result. The sentence template is fixed **now**."""

    prediction: str
    criterion: str                     # Recorded acceptance criterion — **verbatim**
    statistic: dict[str, Any]          # Computed values
    # ⛔ [revision 2026-09-28 · §260 «if the script has to be changed, leave that fact and the reason in the revision record»]
    # Before, there were four — `부합 (consistent) · 빗나감 (missed) · 판정 불가 (not judgeable) · 추정값(판정 없음) (estimate (no verdict))`.
    # `LEDGER_AUDIT` §1 fixes the status vocabulary at **seven** (`부합 (consistent) · 조건부 부합 (conditionally consistent) · 빗나감 (missed) ·
    # not judgeable · not run · discarded · replaced`), and `RESULTS_TABLES` Table 3 as **live verdicts**
    # Actually uses «조건부 부합» («conditionally consistent») (P4 · B7 · B12 · B13). The recorded text governs, so **the code was the outdated side**.
    # ⇒ Accepts the seven. A narrow `Literal` does not create a silence in the record.
    # ⛔ `추정값(판정 없음)` (estimate (no verdict)) is **an eighth word that is not among the seven** (the code made it). It is not deleted now —
    # `estimate_B19p` uses it, and it is **not a verdict state but an indication that «no verdict was given»**.
    # The relation to the seven is left to a **human decision** (batch `rule_scope`, family of `AUDIT_SCOPE` §5-2 #8).
    verdict: Literal["부합", "조건부 부합", "빗나감", "판정 불가",
                     "미실행", "폐기", "교체됨", "추정값(판정 없음)"]
    n_used: int
    excluded: dict[str, Any]           # What was dropped and why
    axis: str | None = None
    note: str = ""
    # ⛔ **True only for verdicts that actually use the decision line.** If a decision line is printed where none is used
    # A reader misreads it as «this verdict was compared with 1.545%» — this goes against the intent of §230-3.
    uses_decision_line: bool = False

    def threshold_block(self) -> dict[str, Any] | None:
        """Axis · stratum · operating point · verdict line **together on one line** (§230-3 · §228-3). Only the verdicts in use."""
        if not self.uses_decision_line or self.axis is None or self.axis not in OPERATING_POINT:
            return None
        b0, s, p = OPERATING_POINT[self.axis]
        # ⛔ [2026-10-04 · re-frozen — relay #22 · #23 ③] The old version had «층» («stratum»): "ROI" — this operating point (`threshold.OPERATING_POINT`: β₀ 0.828 · `s` 0.0237) is
        # it is the recorded value of the **voxel stratum** (§193 · §229-3). The head of the verdict line is a place read by people, so a voxel result could be attributed to the ROI stratum (the axis of §170 · §188).
        return {"축": self.axis, "층": "voxel", "β₀": b0, "s": s, "p": p,
                "판정선_%": None if p is None else round(decision_line(self.axis) * 100, 4)}

    def sentence(self) -> str:
        t = self.threshold_block()
        head = f"[{self.prediction}] **{self.verdict}**"
        thr = (f" (판정선 {t['판정선_%']}% · {t['축']} 축 · {t['층']} 층 · β₀ {t['β₀']} · "
               f"s {t['s']} · p {t['p']})") if t is not None else (
              f" (축 {self.axis} · ⛔ 판정선을 쓰지 않는 판정)" if self.axis else "")
        ex = ""
        if self.excluded:
            n_ex = sum(len(v) if isinstance(v, list) else 0 for v in self.excluded.values())
            pol = self.excluded.get("정책", "—")
            ex = f" 제외 {n_ex} 건(정책 {pol})."
        return f"{head}{thr} — n={self.n_used}.{ex} {self.note}".rstrip()

    def to_json(self) -> str:
        d = dataclasses.asdict(self)
        d["threshold"] = self.threshold_block()
        d["sentence"] = self.sentence()
        return json.dumps(d, ensure_ascii=False, indent=2, default=float)


# ⛔ **Output contract** — one verdict line, delimiter, one JSON block. Do not change the order or the delimiter.
JSON_DELIM = "----- BEGIN VERDICT JSON -----"
JSON_END = "----- END VERDICT JSON -----"


def emit(v: Verdict, frame: Frame) -> None:
    """**The output format is set only here** (Q-P contract). It is not changed after the lock is released."""
    print(v.sentence())
    print()
    print(JSON_DELIM)
    print(v.to_json())
    print(JSON_END)
    if frame.source == "dummy":
        print("\n⛔ **더미 입력이다. 이 값은 판정이 아니다** — 배선만 확인한 것이다(Q-P).")


def cli(desc: str, extra=None) -> argparse.Namespace:
    """Common arguments. `extra(ap)` adds **arguments specific to each verdict** (`--numerator` of B7). ⛔ [revised 2026-09-29 · revision record of §260] The point is **not to give the added arguments a default** — the numerator branch of B7 («predicted»/«observed») is an **item a person confirms** (P9-8-75 (4) #3), so the code must not choose it. The common contract (`--dummy`/`--input`/`--seed`/`--exclude-roi`) is **unchanged**."""
    ap = argparse.ArgumentParser(description=desc)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--dummy", action="store_true", help="합성 입력으로 배선만 확인(봉인 유지)")
    g.add_argument("--input", type=str, help="실데이터 경로 — **봉인 해제 후에만**")
    ap.add_argument("--seed", type=int, default=20260924)
    # §6-1 Δ_ROI exclusion sensitivity (COMBINE_RULE §7): **both the applied and the not-applied versions are reported.**
    # ⛔ The default is the not-applied version — the exclusion is given explicitly «after it is decided to use it».
    ap.add_argument("--exclude-roi", type=str, default="",
                    help="§6-1 로 제외할 ROI 인덱스(쉼표) — 적용본 팔. 비우면 미적용본")
    if extra is not None:
        extra(ap)
    return ap.parse_args()


def get_frame(a: argparse.Namespace) -> Frame:
    return dummy_frame(a.seed) if a.dummy else load_frame(a.input)


def sha256_of(path: str) -> str:
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


__all__ = ["Frame", "Verdict", "JSON_DELIM", "JSON_END", "dummy_frame", "load_frame", "emit", "cli", "get_frame",
           "N_ROI", "SCANNERS", "COMMON4", "REGISTERED_PAIR_STRATA", "first_session_cells", "pair_strata", "recovery_frame", "AXES", "EMPTY_PAIRS", "EMPTY_SPLITS", "structural_zero_lines", "check_empty_splits", "dkk", "decision_line", "OPERATING_POINT", "FISO_T1", "UNDEFINED_RESAMPLES", "physics_by_subject",
           "FISO_NU_MAIN", "FISO_NU_SUPPORT", "physics_nu_envelope"]
