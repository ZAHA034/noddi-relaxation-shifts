# Computation record of the judge

This file holds the record that earlier versions of the article's supplementary material printed in Supplementary Methods S-M5.7. The commits and the record files named here are in the private development repository, whose history is available to editors and reviewers on request.

All judge computation ran on CPU within a budget of 24 h of wall-clock time set before training. The validation share was enlarged by the unused training share, a change made before validation (commit c307e67, 2026-10-25 19:37:06 KST, before validation).

**One deviation.** The first validation run was terminated by a host shutdown after 1.68 h, before any stratum had been written. It was restarted from the beginning with identical code, weights and seeds (commit 6518128; record `prereg/artifacts/decomp/s4_val_restart_record.json`). Timing-probe outputs from both runs were identical in all eight arrays other than the timing fields, and are retained but not used (`decomp/results/s4/probe/PROBE_LABEL.json`, commit 8bb2632).
