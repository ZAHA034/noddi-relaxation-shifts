# Age × scanner TE — analysis plan, fixed before computing (English translation)

> English translation of `decomp/age_te/ANALYSIS_PLAN.md` (SHA-256 8b19e79e7dc56dee46842cb860d15fa1dc0221e64ddf98116a36e90973a1b2c5), checked by script for numbers and identifiers. The Korean file is the record; where the two differ, the Korean file governs.

Written 2026-11-18. In the grades of the manuscript this is **post hoc** (defined after the verdict and fixed before computing). No verdict is reopened. The main run is made only after this document is committed to the repository (`run.py` compares the working copy of this file with the version committed at HEAD and stops if they differ). If anything must change before the run, the place and the reason are written in «Change record» at the end and it is committed again. Nothing is changed after the run. All results are reported, favourable or not.

Basis: proposed by the analysis session (2026-11-18), approved by the author (2026-11-18 — including the first reading of the participant information file, only when the rule of Section 6a says «run»), reply of the review session (2026-11-18, five conditions — change records 1–5).

## 0. Concerns this plan takes up (review session, 2026-11-18)

The review session did not recommend an analysis that measures «a TE slope for each person» in ON-Harmony. There were three reasons.

1. The two protocols of NOT2ING (TE 89 · 98 ms) were scanned in different people. So there is no comparison in which only TE differs within one scanner, and TE always changes together with the scanner.
2. Shifts that differ between people are of a size similar to person × scanner effects arising from head position or registration.
3. The result is likely to end as «none».

In a later reply the review session judged age × TE to be different: it is a question that directly measures the term of the continuous-covariate corollary. It attached conditions, and this plan takes them up as follows.

- E1 does not measure a slope for each person. It measures, pooled over all people, whether age changes the TE slope. Person × scanner noise enters the residual.
- The entanglement of TE and scanner is addressed directly with two falsification controls and one negative control (Section 5).
  - Falsification control 1 is a TE permutation. It uses age × acquisition deviations unrelated to TE as an empirical null.
  - Falsification control 2 is the departure of the acquisition-wise slopes from a line in TE.
  - The negative control is ODI.
- Whether to run is decided before opening the participant file (Section 6a). The criterion is the age-related change of white-matter T2 in verified literature. A test whose result is almost certainly «none» is not run.
- If «none» results, the minimum detectable effect (MDE) is reported with it. That size is translated into a change of extra-neurite T2 by a set rule (Section 6).
- The multi-TE data (E3) are run as well. According to the original report these data have a fixed diffusion time and repeated scans. They are also compared with the assumed value of the manuscript (the between-person relaxation variation of the design simulation).

## 1. Question

The corollary (continuous covariate) of Section 2.4 of the manuscript says:

> slope of apparent logit NDI on a covariate z = slope of the true fraction + TE · d(R2,en − R2,in)/dz.

The manuscript showed the size of the second term only with assumed values. ON-Harmony scanned the same people on several scanners with different TEs. So if person fixed effects remove the true fraction (and the main effect of age), «does the TE slope change with age» can be measured in real data.

- **E1 (primary).** Age × scanner TE. Run only when the rule of Section 6a says «run».
- **E3 (secondary).** The TE slopes of the three people of the multi-TE data. In the signal model this slope is that person's R2,en − R2,in. Whether the three slopes differ from each other, and whether the between-person SD agrees with the assumed value of the manuscript, is examined descriptively. It uses no participant information, so it is run regardless of Section 6a.
- E2 (sex × scanner TE) was dropped (change record 2).

## 2. Data

- **ON-Harmony (E1).**
  - The same first sessions as at the verdict are used. From the 110 sessions of `rconfig.onh_sessions` the 1 exclusion of the main analysis is removed, and only sessions that `judge_v4.session_states` marks «usable» are kept.
  - The acquisition is the scanner. NOT2ING is split in two by protocol (8 card units).
  - TE is the value in `model.TE_ACQ`. Same as Supplement S-M9.7; a TE given as a range uses the middle value.
- **Coordinates.** From `judge_tong.unit_coords` (= `lib.coords`, the same clipping 1e-4). x = logit NDI, y = logit FISO − ln NDI, and ODI.
- **Main value.** For each session, the mean over the units that have values in all sessions (common units). The number of common units is reported.
- **Participant information.** One participant-information tsv of the copy of the dataset metadata is read (E1 only). Only the age column is used (Section 3). This file is not put in the repository (Section 8).
- **Multi-TE data (E3).**
  - Read with `judge_wong.load` (same as at the verdict).
  - Three people × 8 TEs (62–132 ms).
  - The repeat series (TE62R2 · TE132R2) are left out of the regression and used only to measure measurement noise.
  - The series registered separately are included as at the verdict.
- Images (NIfTI) are not read. ROI values are read again from results already computed. No GPU is used.

## 3. Rules for reading participant information (`participants.py`)

- **Age column.** The first column whose name begins with «age».
  - A number is used as it is.
  - An interval «a-b» uses its middle value.
  - Anything else (n/a and the like) drops that person from E1.
  - The unit is 10 years.
- The number of interval values, of missing values, and of people without information are reported.
- If there is no age column or fewer than 3 people can be used, E1 is not run and that fact is reported.
- The sex column remains in the parsing function, but `run.py` does not use it (change record 2).

## 4. Model (E1)

For person i and acquisition a, the value v_ia is modelled as

> v_ia = α_i + μ_a + γ · (TE_a − TE_C) · (z_i − z̄) + ε_ia

- TE_C is the mean of the eight acquisition TEs (84.45625 ms).
- α_i is the person level. It absorbs both the true fraction and the main effect of age.
- μ_a is the acquisition level. It absorbs the main effects of scanner and TE.
- The fit is OLS.

Inference:

- **Interval.** The 95 % interval of γ from the t distribution.
- **Person permutation p.** Age shuffled between people (10,000 times, seed 20261118). p = (1 + #{|γ*| ≥ |γ̂|}) / (10,000 + 1).
- **Minimum detectable effect.** MDE80 = (t_0.975 + t_0.80) · SE. Two-sided α 0.05, power 0.80.
- **Interpretation.** In the signal model, γ of x is d(R2,en − R2,in)/dz (unit: ms⁻¹ per 10 years of age). For interpretation this value is translated in one way. Read as «intra-neurite T2 unchanged and only extra-neurite T2 (control 52 ms, Section 2.4) changing», T2,en′ = 1/(1/T2,en + γ). This is a conversion under an assumption, not a measurement.

## 5. Controls · sensitivity · limits

**Negative control: ODI.** The same model is fitted to ODI.
- In the signal model, compartmental relaxation does not move ODI. The two compartments share the orientation distribution, and compartmental T2 is taken as independent of fibre direction.
- So if age × TE also appears in ODI, that is a sign that the interaction comes from a cause other than relaxation. R-E1a of Section 6 is tied to this result.

**Falsification control 1: TE permutation (empirical null).**
- The data are left as they are, and only the TE values attached to acquisitions are shuffled between acquisitions. γ is computed again with all distinct arrangements. With eight acquisitions this is 20,160 arrangements because OXF1PRI · OXF2PRI share the same TE.
- If age × acquisition deviations unrelated to TE are large, large |γ| arise in any arrangement. So whether the |γ| of the actual arrangement stands out in this distribution is examined. p = #{|γ*| ≥ |γ̂|} / number of arrangements, including the actual arrangement.
- Arrangements that swap acquisitions with nearly the same TE (68 · 69.65 · 70.6 ms, 96.4 · 98 ms) give similar γ. So the smallest possible p is larger than 1/number of arrangements (roughly 24/20,160).
- Run for both x and ODI.

**Falsification control 2: departure of the acquisition-wise slopes from a line in TE.**
- A model with a separate age slope for each acquisition (saturated model) is compared by F with a model linear in TE. The degrees of freedom are the number of acquisitions − 2.
- This value describes the size of age × acquisition deviations not explained by TE.

**Descriptive: y.** Reported descriptively only.

**Three sensitivity analyses.**
- *Vendor-wise covariate-slope terms.* Philips × (z − z̄) and Siemens × (z − z̄) are added (GE reference). In this model γ is fixed only by TE differences within vendors. Those differences come mainly from Philips (70.6–98 ms).
- *Philips-only subset.* Only NOT1ACH · NOT2ING-89 · NOT2ING-98 are used.
- *Unit-level model.* Observations are session × unit, with fixed effects person × unit and acquisition × unit. The SE ignores the correlation between units of the same session, so it is used descriptively only.

**Limits.**
- In ON-Harmony, acquisition conditions change together with TE: diffusion time (not recorded in the sidecars), the b-value scheme, the coil, and the reconstruction. The interaction of these conditions with age cannot be separated from TE in this design.
- The TE permutation only examines «whether the TE order stands out against other orders of the acquisitions». So other conditions that change in the same order as TE still cannot be separated.
- According to the original report, E3 is a TE series within one scanner with fixed diffusion time. So it partly fills this limit.

## 6a. Rule on whether to run (E1 · fixed before seeing the literature and before opening the participant file)

**Literature.** Only studies that meet all of the following are used.
- Studies that measured white-matter T2 of healthy adults with multi-echo spin echo or equivalent quantitative T2.
- Studies that report the relationship with age as numbers.
- Studies in Q1 journals or at top-tier conferences. Archives (preprints) are not used.

**Literature value.** L_max is taken from the literature.
- L_max is the change of white-matter (whole or major tracts) T2 over the age range of ON-Harmony (roughly 20–60 years). The unit is ms/10 years, magnitude only.
- The largest value among the included studies is used.
- A value per year is multiplied by 10.
- If there is only a curve, the difference between ages 20 and 50 is divided by 3.
- The search and reading are recorded outside the repository: search terms, what was opened, what was blocked, and whether it stopped midway.

**Generous conversion (`model.go_rule`).** The literature value is moved toward what E1 can catch most easily.
- Single-exponential T2 is taken as the signal-weighted mean of compartmental T2.
- All of its change is taken to be in the extra-neurite compartment (weight 0.5). Then the extra-neurite T2 change = 2 · L_max.
- γ for this change is set to the smaller of the two directions.

**Power.** Computed with an optimistic SE: SE0 = 2.36 × 10⁻⁴ ms⁻¹/10 years, using σ 0.031 and the balanced approximation of Section 7.

**Decision.**
- If the power is 0.5 or more, «run». This means L_max is about 0.64 ms/10 years or more.
- Otherwise «do not run». Then it is written as «future work» in the manuscript, together with the literature value and the minimum detectable effect of this design.
- If there is no usable literature value at all, also «do not run».

**When «run».** The participant file is opened only after confirming again with the author.

This rule is generous on purpose. If «do not run» results, it means that even under the most favourable reading the probability of detection is below one half.

## 6. Readings fixed before computing

**R-E1a.** First check whether the 95 % interval of x excludes 0 and the person permutation p < 0.05. If so, the reading is set by the two controls.
- If TE permutation p < 0.05 and the interval of ODI covers 0, write: «in these data the dependence of apparent logit NDI on acquisition TE changes with age (post hoc)».
- If TE permutation p is 0.05 or more, lower it to: «an age × acquisition effect is present but not distinguished from age × acquisition deviations unrelated to TE».
- If the interval of ODI excludes 0, lower it to: «an age × acquisition effect not explained by relaxation». If both conditions fail, write both.

In any case, if the sign changes or the interval covers 0 in the model with vendor terms, add «relies on the vendor contrast».

**R-E1b.** Otherwise write: «not detected. The minimum detectable effect at 80 % power is |γ| = MDE80 (in extra-neurite T2 +… / −… ms per 10 years)». Do not write «no effect».

**R-E3.** Read by the χ² p of homogeneity of the x slopes (residual SE of the line fits).
- If p < 0.05, write «the TE slopes (R2,en − R2,in) of the three people differ (descriptive, 3 people)».
- Otherwise write «not seen to differ (3 people, low power)».

In either case, also write:
- the three slopes, the two SEs (residual · repeat), and the estimate τ of the between-person SD.
- the ratio to the assumed value of the design simulation. The assumed value is 0.00121 ms⁻¹, 0.036 logit at 30 ms. It is the same value as the code (`decomp/design_sim/model.py`: lognormal multipliers with coefficient of variation 0.10 on the T2 of the two compartments, within-person correlation 0.8) (change record 4). If the ratio is below 0.5 or above 2, that fact is written.
- the two end values of the correlation: 0.069 if each compartment varies separately, and 0.022 if the two compartments vary together (30 ms, logit).
- homogeneity of ODI (control).

**Adjustment and place.**
- No adjustment for multiple comparisons is made. It is reported stating that E1 is primary and E3 secondary.
- Whatever the results, they are placed in the supplement (a new subsection and a row of Table S5). The main text has only one or two sentences (change record 5). The grade is post hoc.

## 7. Power estimate (before computing, only from values already public)

**Assumptions.**
- σ is set to the retest SD of x for OXF1PRI, 0.031 logit (manuscript 4.1). The actual residual adds person × scanner effects, so 0.06 is also considered.
- Approximated as a balanced design: 110 observations, SD of the eight TEs 11.95 ms, SD of age 1.05 (10 years, summary of the data paper 31.5 ± 10.5 years).

**Estimate.** SE ≈ σ / √(110 · 11.95² · 1.05²).

| σ (logit) | SE (ms⁻¹/10 years) | MDE80 (ms⁻¹/10 years) | Translated to extra-neurite T2 (per 10 years) |
|---|---|---|---|
| 0.031 | 2.36 × 10⁻⁴ (= SE0 of Section 6a) | 6.7 × 10⁻⁴ | about +1.9 / −1.7 ms |
| 0.06 | — | 1.3 × 10⁻³ | about +3.7 / −3.3 ms |

Person fixed effects remove the between-person information, so the actual SE is larger than this. The actual MDE is reported with the SE of the main run.

## 8. Personal information

- The participant-information tsv and per-person ages are not put in the repository, the public package, or the manuscript. The DENY list of the repository (`prereg/scripts/verify_seal_tracked.py`) blocks this file name.
- The output is only a model summary: coefficients, intervals, p, number of people, mean and SD of age.
- `run.py` stops without writing if the output text contains a participant number (`participants.assert_no_ids`).
- The people of E3 are written only as participant_1–3.
- The output file is `decomp/results/age_te/summary.json`. It is committed and added to `decomp/results/committed_sha256.json`.

## 9. Order of execution

1. Commit this plan, the code (`model.py` · `participants.py` · `run.py`) and the tests (`tests/test_age_te.py`).
2. Search the literature and reach the decision of Section 6a. The decision is recorded outside the repository, and is added as `decomp/age_te/E1_GO_DECISION.md` when the results are committed.
3. Run E3. The rereading of the multi-TE ROI values is logged first. The run is made without `--participants`.
4. Run E1 only when the decision is «run» and the author has confirmed again. In this order:
   1. The participant-information tsv is entered in the log of new data reads with its sha256. The ON-Harmony ROI values are entered in the log of rereads.
   2. Run with `--participants`.
5. Commit the results and report all of them.

## 10. Tests (`tests/test_age_te.py`, synthetic data · no real data read)

**Tests with known answers.**
- Without noise, γ is recovered exactly.
- γ stays the same even when the main effects of age and TE are large.
- With noise, the interval coverage is 0.91–0.98, and the SE is 0.85–1.15 times the Monte Carlo spread.
- Under the null, the rejection rate of the t interval is 0.02–0.09 and that of the person permutation is 0.15 or less.

**Failure controls.**
- Fitting with shuffled TE gives a wrong γ.
- With only a vendor × age term, the main model is fooled and the sensitivity model is not.
- Acquisition-wise age slopes unrelated to TE fool the person permutation but not the TE permutation. The departure F catches them.

**Other tests.**
- An interaction linear in TE stands out in the TE permutation (p < 0.05). The departure F is then not small.
- The threshold of the rule on whether to run is checked (0.62 is «do not run», 0.66 is «run»).
- The E3 tools are tested: the reciprocal relation of the extra-neurite T2 conversion, the TE slope, homogeneity, the assumed SD (0 when the variation is 0), τ, MDE.
- The parsing rules for participant information are checked.
- The output guard is checked: it stops when a participant number enters.

**Mutation tests (before the commit, recorded outside the repository).** 16 faults were planted one at a time, and the tests caught all 16.
- 11 faults of the first version (R38p): dropping the acquisition fixed effects, the sign of the conversion, the power of the MDE, the person permutation p, the homogeneity weights, the unit of age, the SE element, the fixed effects of the unit model, the centre of the interval, the output guard, the sex mapping.
- 5 faults of this version: the TE permutation not actually changing TE, dropping the acquisition-wise slopes in the departure model, the weight of the rule on whether to run, a one-sided tail of the power, the comparison of the TE permutation p.
- A slope with TE not centred is a mutation giving the same value (equivalent) and was not counted.

## Change record

All changed on 2026-11-18, before the run. The first version is commit R38p. The basis is the reply of the review session (five conditions) and the author's message.

1. **Falsification controls added.** A TE permutation (using age × acquisition deviations unrelated to TE as an empirical null) and the departure of the acquisition-wise slopes from a line in TE were added. ODI was changed from «weak control [assumed]» to a negative control tied to R-E1a. Reason: to address directly the problem that TE changes together with the scanner.
2. **E2 (sex × scanner TE) dropped.** Reason: the sex groups are small and the power is practically nil. Including it would also require changing the sentence on sex in manuscript 5.3. The review session said «10 versus 10». The summary of the data paper first used by this plan says «7 female». Either way the conclusion is the same.
3. **Rule on whether to run (Section 6a) added.** If the age-related change of white-matter T2 in the literature is hard to catch with this design, E1 is not run. The approval to read the participant file was also tied to this rule. Reason: a test that is almost certainly «none» only increases the number of post hoc analyses (Supplementary Table S5 already has 56 rows).
4. **The reference value of E3 (0.036) reconciled.** Supplement S-M11 and Table S13 (e) of the manuscript only say «CV 0.10 for each compartmental T2». Read as written, this gives 0.069 (separately) or 0.022 (together) at 30 ms. The code of the design simulation uses a within-person correlation of 0.8 between the two compartments (`CORR_V = 0.8`, also written in Section 2 of that plan). So 0.036 is the same value as the code. The correlation 0.8 is written into the manuscript description (next version of the manuscript). The reference of E3 is 0.036, and the two end values are written with it.
5. **Place of the results set.** Whatever the result, it goes to the supplement, and the main text has only one or two sentences. Reason: so that the compression of the manuscript and the S88 fixes are not delayed by this experiment.
