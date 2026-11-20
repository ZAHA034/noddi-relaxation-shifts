# Design simulation addendum M6 — a second TE within each site, fixed before computing (English translation)

> English translation of `decomp/design_sim/PLAN_M6.md` (SHA-256 2ed0593c332bd1b258563f677c29cc4791287fca52021dbefd563926ea996897), checked by script for numbers and identifiers. The Korean file is the record; where the two differ, the Korean file governs.

Written 2026-11-18. In the grades of the manuscript this is **post hoc** (defined after the verdict and fixed before computing). All inputs are assumed values, and the results are used in the manuscript together with «illustration». The original plan (`ANALYSIS_PLAN.md`) and its results (288 cells) are not changed. The main run is made only after this document is committed to the repository. If anything must change before the run, it is written in «Change record» at the end and committed again.

Basis: Section 3 C3 of review direction S89 (proposed by the review session, approved by the author — pending items proceed as the analysis session recommends). The aim is to put numbers on the recommendation of manuscript 5.4 («diffusion-weighted images at a second TE»). The current prescription of the manuscript, extrapolation between sites (M3), has an SE about 7 times that of the pooled difference, and a power of 0.12 at the both-compartment size.

## 1. Design

- Data generation uses `model.generate` as it is. Cells, seed rule, people, sites, ROIs, relaxation and noise are all the same.
- An observation at a second TE is added to it (`m6.py`).
  - Same person, same site, TE₂ = TE_s + ΔTE, ΔTE = 30 ms.
  - The true fraction, the between-person fraction variation b, the site · ROI terms γ and δ, and the site × group term η are the same as in the first observation. Same scanner · same session is assumed.
  - Only the relaxation term changes to TE₂.
  - The measurement noise is drawn anew. Its SD is exp(ΔTE / 70 ms) times that of the first observation (assumed). This is the share by which the signal falls at the longer TE; 70 ms is an approximate tissue T2.
- The time taken to acquire the second observation is not put in the model. It is only recorded.

## 2. Analyses

**M6a — every participant is scanned at the second TE.**
1. For each person i and ROI r, a value extrapolated to TE → 0: y₀ = y₁ − TE_s · (y₂ − y₁) / ΔTE.
2. M2 (OLS with a fixed site term) is applied to this value to give the group coefficient.

In the signal model this value removes exactly each person's relaxation term (including the variation between people), leaving only the T1 · proton-density term. In this simulation that term is common to the groups.

**M6b — only half of the participants of each site · group are scanned at the second TE.** The first half in numbering order within each site · group.
1. In that half, the group difference in TE slope is estimated for each site: the group difference of (y₂ − y₁)/ΔTE.
2. From the site-wise group difference computed from everyone's y₁, the product of TE_s and that slope difference is subtracted.
3. The sites are combined with inverse-variance weights.

The SE adds, for each site, the variances of the two differences (including the covariance of the two estimates within the site).

**Comparison.** M2 and M3 (te0) of the same replicates are kept alongside.

**Target cells (16).** D1 × eight scenarios (S0, S1b, S1e, S2b, S2e, S2i, S3b, S3e) × V (0, 0.10) × H0 × balanced design. Replicates are 1,000 per cell, and the seeds are the SeedSequence of the original cell numbers (the same data). The noise of the second observation is drawn with SeedSequence([20261118, cell number]).

## 3. Reporting

For each method (M6a, M6b, M2, M3):
- bias and its Monte Carlo SE
- coverage of 95 % intervals
- rejection rate in cells with a fraction difference (power) and in cells without one (false positives)
- mean SE, and its ratio to the mean SE of M2 in the same cell

## 4. Readings fixed before computing

**R-e (validity of M6a).** If the coverage of M6a is 0.93–0.97 in all cells where only relaxation differs (S2b, S2e, S2i; V 0 and 0.10), write:
> «with a second TE in every participant, intervals for the fraction-only difference were valid in this design»

If any is outside, write those cells and values.

**R-f (precision).** Report the ratio of the mean SE of M6a · M6b to the mean SE of M2 in the same cell, in S1b and S1e, together with the power, alongside the ratio of the extrapolation M3 (about 7). Descriptive, with no threshold.

**R-g (M6b).** If the coverage of M6b in the same cells is 0.93–0.97, write «half the participants suffice for valid intervals in this design». If outside, write the values.

**R-h (variation of relaxation between people).** Write the coverage and power of M6a in the V = 0.10 cells alongside the V = 0 cells (descriptive).

## 5. Place in the manuscript

- A new sub-paragraph and Table S13 (f) at the end of Supplement S-M11.
- In the main text, one or two numbers in the recommendation sentence of 5.4.
- Only one row is added to Figure 4 and its height is not increased. If that is not possible, the figure stays as it is.
- The grade is post hoc.

## 6. Tests (`tests/test_design_sim_m6.py`, 18 checks)

- **Same data.** The first observation equals the y of `model.generate` to floating point (three cells).
- **Tests with known answers.**
  - Without noise and site terms, M6a · M6b recover θ in six scenarios.
  - M6a is exact even when relaxation differs between people.
- **Failure controls.**
  - Using 0 instead of TE_s in the extrapolation formula is wrong in the relaxation cells.
  - ΔTE = 0 stops.
- **Shared common terms.** With measurement noise set to 0, y2 − y1 is for each person a constant that does not depend on the ROI. In cells where all relaxation is the same, it is one constant for all people.
- **Other.**
  - The noise multiplier of the second observation exp(30/70), the half mask, and the 16 target cells are checked.
  - The SE agrees with the spread over 150 replicates (ratio 0.8–1.25).
- **Mutation tests (before the commit, recorded outside the repository).** 8 faults were planted, and 6 were caught.
  - Caught: the sign of the extrapolation, the slope coefficient, the noise multiplier, and dropping each of b · γ · η from the second observation.
  - Not caught: the covariance term of the M6b variance formula (sign flipped or dropped). This term changes the SE of these cells by only about 4 %, so the SE-ratio test over 150 replicates cannot separate it.
  - So the validity of the intervals of M6b is judged only by the coverage of the main run (R-g).

## Change record

(none)
