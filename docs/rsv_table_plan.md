# Table applying the relaxation sensitivity value to published clinical results — plan, fixed before seeing the values (English translation)

> English translation of `decomp/rsv_table/PLAN.md` (SHA-256 3183edda2531fdb90b907eb68effa2f5ae3420d6f9f37a276497fe11c86addd3), checked by script for numbers and identifiers. The Korean file is the record; where the two differ, the Korean file governs.

Written 2026-11-18. In the grades of the manuscript this is **post hoc**, and in character an **illustration** (control compartmental T2 is assumed). No new image data are read. Only published numbers are used. The values of the candidates are read only after this document is committed to the repository. What was done before the commit was an access check of whether texts open and whether the fields exist. Its record is in a work note outside the repository (r41/C2_ACCESS_2026-11-18.md).

Basis: Section 3 C2 of review direction S89 and its addendum (proposed by the review session, approved by the author — pending items proceed as the analysis session recommends).

## 1. What is computed

The mean white-matter NDI of the control and patient groups of each study (or their difference) and the TE are used.

- **d.** d = logit(NDI_patients) − logit(NDI_controls). It is the logit difference of group means, so it can differ from the mean of person-level logits. This limit is written in the table notes.
- **Relaxation sensitivity value.** d / TE. The unit is s⁻¹, and it is the group difference of R2,en − R2,in (manuscript 2.4).
- **T2 conversion (two assignments).** Control T2,in 82.5 ms and T2,en 52 ms (manuscript 2.4) are used.
  - (a) extra-neurite only: the % lengthening of T2,en that gives that difference
  - (b) both compartments in the same proportion: the same % lengthening
  - If d < 0 (apparent NDI decrease) it is written as a lengthening, if d > 0 as a shortening.

## 2. Studies

**(i) Clinical results cited in the introduction.** Following the access check, Collorone 2020 (MS, brain NAWM NDI) is included. The other three are not included, and the reasons are written in the table notes.
- Cox 2016: no raw-scale slope and no TE.
- Hall 2026: the values are only in figures.
- Høgestøl 2025: not NODDI.

**(ii) PubMed query (this sentence as it is, with the date of the query recorded).**

```
("neurite orientation dispersion and density imaging"[Title/Abstract] OR NODDI[Title/Abstract]) AND ("white matter"[Title/Abstract]) AND (patients[Title/Abstract]) AND ("healthy controls"[Title/Abstract] OR "healthy control"[Title/Abstract] OR "control subjects"[Title/Abstract])
```

- The results are screened in PubMed's descending order of publication date.
- Stop when 6 studies meet all inclusion criteria, with at most 2 per disease category.
- Every screened record is given «included» or «reason for exclusion». At the stage of judging from abstracts and bibliographic records only, no result values are copied.

**Inclusion criteria.**
1. An original article in a journal of SJR 2025 quartile Q1 (looked up in a copy of the public export).
2. Compares human patients with healthy controls.
3. Reports NODDI NDI of white matter (normal-appearing white matter, whole white matter or major tracts) as a mean for each group, or the mean difference as a number (in a table or the text; excluded if only in figures).
4. Reports the TE of the diffusion images.
5. The full text can be read (PMC, open on the publisher's site, author-accepted manuscript).
6. Archives (preprints) are excluded.

## 3. Reported T2 change of the disease (comparison column)

- For each disease in the table, search with the query: `"<disease name>" AND "white matter" AND ("T2 relaxometry" OR "T2 relaxation time" OR "quantitative T2")`.
- Use one Q1 original article that reports the patient − control difference of white-matter single-exponential (bulk) T2 in %. The value must be read from the full text. Brandhofe 2021 (SVD), already in the manuscript, is used as it is.
- If none is found or it cannot be read, write «not reported». The row is still kept.
- The bulk T2 change is written alongside (a) · (b) without compartmental assignment. The limits of the assignment are written in the notes (the same words as manuscript 2.4).

## 4. Readings fixed before computing

- If the reported bulk T2 lengthening (%) is at least the % of (b), write: «could be explained by relaxation alone (both compartments)». The same rule is written separately for (a).
- If the reported lengthening is smaller than both (a) · (b), write: «would need a larger relaxation difference than reported».
- If there is no reported value, write: «no reported T2 change to compare».
- «artifact», «spurious» and «validation» are not used.
- Manuscript 5.3 gets one sentence: «this table is not a validation and does not test Proposition 4».

## 5. Recording of values

- The extracted values (group mean NDI, their difference, TE, group sizes) are put in a ledger outside the repository together with the table number or page of the source.
- The computation is done by a derived script. That script has a known-answer test (manuscript 2.4: 0.02 NDI · TE 90 ms → about 5 %) and a failure control.
- For each row the table gives the index, region, group sizes, TE, TR (if reported), d, sensitivity value (s⁻¹), % lengthening of (a) · (b), the reported T2 change and its source, and the reading.

## 6. Place in the manuscript

The numeric paragraph of main-text 5.2 (example of the change at a TE difference of 30 ms) is replaced by this table. About 0.4 page of main text. If the table is long, the main text keeps up to three rows and the rest goes to the supplement.

## Change record

1. **2026-11-18, before reading the values.** Screening (in the same way as the access check, values not copied; record outside the repository r41/C2_SCREEN_2026-11-18.md) produced cases the plan had not settled, and these are settled before seeing the values.
   - **(a) Region.**
     - If a study reports the NDI of whole white matter or whole normal-appearing white matter, that is used.
     - If not, all tracts reported in the tables of the main text are used, and the row gives the range (minimum–maximum) of the sensitivity value and the number of tracts.
     - A single tract is not picked.
   - **(b) Adjusted coefficient.** If only an adjusted group coefficient (a difference in raw-scale NDI) is given instead of group means, it is used as the mean difference. The control NDI level must then be in the main text. If not, the row is dropped as «NDI level not readable».
   - **(c) Continued screening.** If a row drops at the stage of extracting values (NDI level or TE not readable), the next records in screening order (from number 94) are screened by the same criteria until the included studies number 6. If the records run out, only what remains is kept.
   - **(d) Full texts that do not open.** Records whose full text did not open on this server («blocked», 18 records) are taken as not meeting criterion 5 (the full text can be read). Their number is written in the table notes.
   - **(e) When the reported T2 change is given only as group means.** The % is computed as (patients − controls)/controls × 100. If there are several regions, the range is written.
   - **(f) Multiple sclerosis rows.** Both (i) Collorone 2020 and the rows of (ii) are kept (the cap of 2 per disease applies only to (ii)).
2. **2026-11-18, after extracting the values · before computing.** The plan had not settled the case of a study that compared several patient groups with the same control group.
   - Studies concerned: Sone 2018 (L-TLE · R-TLE and TLE-HS · MRI−/PET+ — the same 33 people split in two ways), Chen 2025 (NMOSD · MOGAD).
   - Rule: on the same principle as change record 1 (a) (no picking), all patient groups printed alongside the control group in the tables of the main text are used.
   - Patient groups with different diseases (Chen) get a separate row for each disease.
   - Several patient groups with the same disease (Sone) are put in one row with the range over all patient groups × tracts, and the fact that the same patients were split in overlapping ways is written in the notes.
   - This rule was written after seeing some values of the two studies (checking against the original, the forceps rows of Sone Tables 6–7). It is a rule that does not select, so it does not depend on the values.
   - The original does not state whether the means of Takahashi 2026 Table 3 are before or after ComBat. This is written in the notes.
3. **2026-11-19, after computing (after seeing the values).** Two things are added. The readings fixed before computing (Section 4) are not changed and are reported as they are.
   - **(a) Three rows in the main text.** Section 6 only said «up to three rows in the main text». The main-text table carries the row of (i) (Collorone 2020) and the first two rows of (ii) in order of inclusion (Wang 2026, Takahashi 2026). It is an order rule, so it does not depend on the values. All rows are in the supplementary table.
   - **(b) Bulk-conversion column.** The reading of Section 4 compares the reported bulk T2 change directly with the change of compartmental T2 ((a) · (b)). If both compartments are lengthened in the same proportion, bulk T2 also lengthens in nearly the same proportion, so the comparison with (b) is like with like. If only the extra-neurite compartment is lengthened, bulk T2 lengthens less (note of manuscript Table S11: 3–4 % for a lengthening of 9.84 %). So the comparison with (a) makes the needed difference look larger.
   - To show this, the bulk T2 change that the lengthening of (a) gives in the signal model is written with it. The signal model is control ν 0.46, T2,in 82.5 ms, T2,en 52 ms, and the TE ranges of the single-exponential fits are 10–100 · 10–120 · 70–100 ms (the same as the note of Table S11). This value is computed by a derived script of the ledger outside the repository. The known-answer test checks equality with the bulk values of s84 within 1e−12, and the failure control is a copy that lengthens the intra-neurite compartment.
   - This column was added after seeing the values and is not a reading rule. The main text gives the reading fixed before computing and this column side by side.
   - The bulk T2 search of §3 was completed on this day for the following diseases (query form as in Section 3, sorted by relevance). Galactosemia («galactosemia») 0 records, MOG antibody-associated disease («MOG antibody-associated disease») 0 records. Both diseases are «not reported».
