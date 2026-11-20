# Earlier registered study: table of registered predictions and recorded outcomes (English translation)

- **Study.** The earlier study of the same ON-Harmony data, registered on the Open Science Framework as osf.io/5juwv (registered 2026-09-26). It is not the study of this article, and the article tests none of its predictions.
- **Source.** Translation of `docs/earlier_study_outcomes_ko.md`, an excerpt of `prereg/RESULTS_TABLES.md` of the development repository (lines 65–140, commit `0a471e1`, 2026-10-13; SHA-256 of the whole source file `d2c0fe4ba7b83a89fb768ffefa904fb0e81bfee8a5c9ad54d7e2178928569bab`). The Korean excerpt remains the record.
- **What it is.** All 29 registered entries of that study with the outcome that study recorded for each, copied as recorded and not judged again. The judgment outputs named in the table are not included, and they were not opened for this excerpt.
- **Translation.** Made under the same rule as the other translations in `docs/`: the meaning is not changed. The fixed verdict words are translated as: 부합 consistent · 조건부 부합 conditionally consistent · 빗나감 missed · 판정 불가 could not be judged · 미실행 not run · 폐기 discarded · 교체됨 replaced. Section numbers (§), relay names and file names refer to the records of that study.
- **Note by this repository.** The paragraph «Revision» below counts column B as 9 run, 17 not run and 3 discarded or replaced. The rows of the table as last changed give 15 run, 11 not run and 3 discarded or replaced. The footnote on B7 gives its real-data arm as «not run», whereas the B7 row gives «could not be judged». Both are reproduced as recorded.

---

## Table 3 — Main verdict table (**all registered predictions · fixed at 29 rows**)

> ⛔ The first version had 15 rows. The ledger actually has **29 rows** (including the struck-through `~~P5~~` and `~~B19~~`), so the first version violated rule 1.
> The statuses and the times of confirmation were transferred by machine from the **ledger audit** (`prereg/artifacts/audit/ledger_status_audit_20260927.txt`),
> not written by hand (R22). There are seven status words (LEDGER_AUDIT §1).

| Prediction | Test unit | Verdict | Time of withdrawal or confirmation | A. **Rule** fixed before unsealing | B. Was it **run**? | Note |
|---|---|---|---|---|---|---|
| P1 | voxel | **Could not be judged** — signature 3 · deviation D-6 (the registered index ψ is not defined on the real-data judgment path · kept in the denominator) | 2026-10-01 signature 3 · D-6 | ✅ | ❌ not run | No real-data computation — the simulation reading is given only as an estimate (no verdict · D-6 addendum) |
| P2 | rank | **Could not be judged** — transcribed from 40_judge | 40_judge 2026-10-04T19:30:36Z | ✅ | ✅ run |  |
| P3 | interval | **Could not be judged** — transcribed from 40_judge | 40_judge 2026-10-04T19:30:36Z | ✅ | ✅ run |  |
| P4 | voxel | **Conditionally consistent** (physics formula only) — when passed through the estimator, **exceeds the decision line 2.88–4.49 times** | 2026-09-22 `71dbd85` | ✅ | ✅ run | ⚠ The name and the number point in opposite directions → footnote |
| ~~P5~~ | scanner pair | **Discarded** | 2026-09-21 `238d356` | ✅ | — discarded |  |
| P6 | ROI | **Not run** | no verdict | ✅ | ❌ not run |  |
| P6′ | rank | **Missed** — transcribed from 40_judge | 40_judge 2026-10-04T19:30:36Z | ✅ | ✅ run |  |
| B1 | voxel | **Not run** | no verdict · `238d356` | ✅ | ❌ not run |  |
| B2 | voxel | **Not run** | no verdict · `2eb789f` | ✅ | ❌ not run |  |
| B3 | voxel | **Not run** | no verdict · `2eb789f` | ✅ | ❌ not run |  |
| B4 | voxel | **Not run** | no verdict · `2eb789f` | ✅ | ❌ not run | ⛔ Element-by-element comparison closed — **no attribution** (tissue sample and SNR level do not match) → footnote |
| B5 | voxel | ~~consistent~~ → **Not run** | withdrawn 2026-09-27 | ✅ | ❌ not run | ⛔ **Same reason** as B4 (not the registered tissue sample) → footnote |
| B6 | voxel | **Not run** | no verdict | ✅ | ❌ not run |  |
| B7 | scanner pair | **Conditionally consistent** (literature + simulation arm) · **Could not be judged** (real-data arm · transcribed from 40_judge) | 40_judge 2026-10-04T19:30:36Z | ✅ (literature and simulation arms) | ✅ run | ⚠ Two arms → footnote |
| B8 | voxel | **Not run** | no verdict | ✅ | ❌ not run |  |
| B9 | voxel | **Not run** | no verdict · 2026-09-22 `7aa89dd` | ✅ | ❌ not run |  |
| B10 | voxel | **Consistent** (condition: tissue d∥ = fitted d∥) | 2026-09-22 `83e8fa8` | ✅ | ✅ run | ⚠ Rerun under the registered condition = **B10′(i) missed** → footnote |
| B10′ | voxel | ~~composite~~ → **Missed** — (i) missed · (ii) consistent · (iii) missed · (iv) not run (four parts of equal standing · §6-3 · 2026-10-01 · ROW_TABLE) | 2026-09-22 `d420567` | ✅ | ✅ run | ⚠ Composite and conditional → footnote |
| B11 | interval | **Could not be judged** | 2026-09-22 `0fdeaab` | ✅ | ✅ run |  |
| B12 | voxel | **Conditionally consistent** | 2026-09-22 `0fdeaab` | ✅ | ✅ run |  |
| B13 | voxel | **Conditionally consistent** | 2026-09-22 `9af3f6b` | ✅ | ✅ run |  |
| B14 | voxel | **Missed** | 2026-09-22 `213e08f` | ✅ | ✅ run |  |
| B15 | voxel | **Replaced** | — | ✅ | — replaced |  |
| B15′ | voxel | **Not run** | no verdict | ✅ | ❌ not run |  |
| B16 | voxel | **Consistent** | 2026-09-22 `68dd097` | ✅ | ✅ run |  |
| B17(β) | scanner pair | **Missed** — transcribed from 40_judge | 40_judge 2026-10-04T19:30:36Z | ✅ | ✅ run |  |
| B18 | voxel | **Consistent** | 2026-09-22 `ff03776` | ✅ | ✅ run |  |
| ~~B19~~ | scanner pair | **Discarded** | 2026-09-22 `6b5b406` | ✅ | — discarded |  |
| B19′ | participant | (verdict cell left empty) «estimate (no verdict)» — transcribed from 40_judge | 40_judge 2026-10-04T19:30:36Z | ✅ | ✅ run |  |

### ⛔ Why the column was split in two (RV42-03)

Previously there was **one column, «before unsealing»**. That one cell carried **two meanings**: «was the verdict fixed before unsealing?» and
«was the verdict **rule** fixed before unsealing?». Depending on which is read, the same table becomes
**«all 29 rows ✅»** or **«9 ✅ / 20 ❌»**. ⛔ **A column that carries two meanings is read in the favorable way.**

⛔ **Revision**: after splitting, an actual count shows that column B has **3 values**: **✅ run 9 · ❌ not run 17 · — discarded or replaced 3**
(total 29). «20 ❌» was a number obtained only by merging the 3 discarded or replaced entries into ❌. Splitting exposed
even that merge; this is the effect of splitting the column.

> **[Registered] Which meaning is right is not decided. **The column is split in two.****
>
> | Column | What it says | How it is filled |
> |---|---|---|
> | **A. Rule fixed before unsealing** | Whether the verdict **rule** (criterion, threshold, procedure) was registered before the first real-data value | ✅/❌ + commit of fixing |
> | **B. Run** | Whether the judgment **ran** | ✅ run / ❌ not run / — discarded or replaced |

- Reference time: the first real-data value is `10_roi_s` at **2026-09-27T09:33:35Z**. A ✅ in column A means fixed **before that time**.
- ⛔ One exception is declared: the whole-brain mean of 1 repeat pair in P9-8-22 (not a judgment input, `LEDGER_AUDIT` §3).
- Column B is **not interpreted**: it follows mechanically from the word in the verdict cell («not run» → ❌).

★ **This is not only more honest but also stronger.** What we can actually claim is
**«the rules were all fixed before unsealing, and only some were run»**, and that is enough.
The ambiguity disappears **through the structure of the table, not through interpretation**.

**Rows that need a footnote**
- **B7**: the verdict splits: (a) the **literature + simulation arm** = conditionally consistent (2026-09-22) · (b) the **real-data arm** (`judge_B7.py` + `30_scanner_cmp`) = **not run**. The ledger and the snapshot do not distinguish the two arms → written as two lines in the table.
- **B10′**: the verdicts of the four parts (i–iv) differ and do not fold into one label; they are given by part.
- **B11**: ledger «not run» · snapshot «judged (§188)» · condition table «awaiting verdict»: a **3-fold mismatch** → the difference between the interpretation rule (if ambiguous, not run) and the audit verdict (could not be judged) is written in the difference table.
- **B15 / B15′**: the live one is **B15′**. The Korean and English ledgers differ (EN marks it superseded, KR keeps the original text).
- **P4**: the verdict name (conditionally consistent) and the number (exceeds the decision line 2.88–4.49 times) read in opposite directions → both are written in one cell.
- **B4**: an element-by-element comparison (audit §2-C) shows that the values of §188 are **not observations of this prediction**: registered «the same 3000 tissues» against an actual **fixed grid of 120**, and registered b0 SNR **{15,25,40}** against an actual **single value of 25**. The §188 values are given only as description, **without attribution**.
- **B5**: it came from the **same run and the same commit** (`414e104`) as B4, and the denominator of its verdict is the quantity of B4. The two mismatches that blocked the attribution of B4 **apply equally to B5**; leaving the two with different grades is an **open item** (audit §2-C-1, human decision).
- **P1 · P3**: ~~not assigned to any step of the unsealing plan (no judgment script).~~ ⛔ Outdated — revised on 2026-10-13 (relay TT-6 · non-reading batch 58): P3 was judged at step 40 by the judgment script `prereg/scripts/lockd/judge_P3.py` (the P3 row of Table 3 above = transcription of step 40). P1 was disposed of without a judgment script, by signature 3 · deviation D-6 (the P1 row of Table 3 above).

**Fixed verdict words** — consistent / conditionally consistent / missed / could not be judged / not run / discarded / replaced. No other words are used.
**Fixed robustness words** — robust (holds over the whole grid) / assumption-dependent (the corner where it flips is stated) / not applicable.
