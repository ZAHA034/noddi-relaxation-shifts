# Fixed products, commits and the order of the analysis steps

This file holds records that earlier versions of the article's supplementary material printed (Supplementary Table S2 and Supplementary Fig. S1 of those versions). Commit hashes refer to the private development repository, whose history is available to editors and reviewers on request. Times are KST (UTC+9); commit times are author dates as in `git log`. Paths are those of the development repository.

## (a) Commits

| Item | Commit | Time (KST) |
|---|---|---|
| Fixed document (`decomp/remeasure/MEASUREMENT_FIXED.md`) | 5f95b18f2d6ee5c8f0ad849f6b9177c7681239df | 2026-10-28 01:56:00 |
| Pipeline (processing and analysis code) | 08a189c65b4bc96fd2a65e4c4e7f0a7a7655595b | 2026-10-28 00:00:25 |
| Frozen judge weights | 85eb8b0 | 2026-10-25 19:13:31 |

Full commit hashes are from the first entry of the image access record (2026-10-27T17:35:55Z), written before the first image was read. That entry names no file, so it has no row in `records/access_records.csv`; the same two hashes are in all 158 processing records `data/remeasure/*/*/DONE.json` (field `근거`, "basis", keys `fixed_commit` and `pipeline_commit`). The pipeline time is stated in `decomp/remeasure/MEASUREMENT_FIXED.md` (line 430).

- Commits made after the real-data results (author dates, KST, as in `git log`): 5632626 (2026-11-01 15:20:30; judgment results, results summary and handoff record); d1e4576 (2026-11-01 15:48:42; card summary and recomputations, with input hashes); 73e1780 (2026-11-01 16:17:02; marker check); 9bf26b2 (2026-11-01 16:44:40; a snapshot of the manuscript taken after the real-data results were computed and before they were written into the text, with build files and figure scripts); 3e41a98 (2026-11-01 17:12:19; application scripts, compression cards and the SHA-256 list of the working records).
- Added after the real-data results, from the stored V4 output without recomputation (`data/remeasure/onharmony/v4/v4_c3.json` (sha256 f6003589…)), 2026-11-01 23:24 KST: the attribute block of the session pairs and the units behind it (printed in earlier versions of the supplementary material; the values are in that file), the unit exclusions and the NDI marker, and, in Supplementary Table S3, the units before → used and excluded (misspecified / k unfilled).

## (b) The 32 guarded files, with the SHA-256 of each file at the pipeline commit (`git show 08a189c:<path> | sha256sum`), as listed in the fixed document; any change after release is a deviation, and the loader guard (`decomp/remeasure/guard.py`, list PIPELINE) does not open until the fixed document is committed again.

| Group | File | SHA-256 |
|---|---|---|
| Processing and prediction code | `decomp/remeasure/rconfig.py` | d35b0a650d1d7a18dc5079e4a1c1442b26309633302482dab06c134756845cfb |
| Processing and prediction code | `decomp/remeasure/guard.py` | 2868199e32401dc280fc6360d1f3965f7f2d6649673ee1ff4946999ca698c784 |
| Processing and prediction code | `decomp/remeasure/cmds.py` | 3e09ca7e81b506344c5ede0be6110a2db9c9828aafc74e1566b770fc8ee3e926 |
| Processing and prediction code | `decomp/remeasure/roi.py` | 7f47cafd8bc574273646cb604861fb71e6c48056d40c201e55c7cd101530e72c |
| Processing and prediction code | `decomp/remeasure/fit_f2.py` | bdcf79a991e145b71f884363bfbf8e99de41712c4af6d4ec9a59212376717c70 |
| Processing and prediction code | `decomp/remeasure/run.py` | 5e8c232fdd3cf864a7962a9c4b8e21580c9a6b84b28444cbbe35ea6f30d22c8e |
| Processing and prediction code | `decomp/remeasure/analysis.py` | faa0e20e84c6c497d5c48394624dfdd10cd7bf0bb67c41f1a84b8db82d0425a3 |
| Processing and prediction code | `decomp/remeasure/sessions_onharmony.tsv` | 9810dbcbda3726efc4b65ca97b0c4d3b538073d01702a66e2d0364b8db0d76ca |
| Processing and prediction code | `decomp/fixed/fiso_sign.py` | b7864ac5ffe16dc346f7c0d8b0d65dffa0bfdeadac159a937bbea75ed5b3bc85 |
| Processing and prediction code | `decomp/fixed/delta_tau.py` | fc2fa38b8d07862b23c968dcdff03e3e4ca3290301aaf5914a27c37b279a5d52 |
| Processing and prediction code | `decomp/stage1/lib.py` | 6c72ad9ce5a2bb2dca5ff23625cf27a42b3b89493cc46fc592f1facd8d8d00c1 |
| Processing and prediction code | `decomp/stage1/synth.py` | 825994813ce8f7ebe15c6f52f60a5d239f37e9b0b044efe1b820ccc95114c347 |
| Processing and prediction code | `decomp/stage1/config.py` | acde78213382e0dd88903425045b1d25db0a84a8c67ea8e0831971a2fc40f8d7 |
| Processing and prediction code | `decomp/stage1/analyze.py` | f1c6f7a51434edb560b55f49c6842cdbd14f07483cea20fef6638dc4be306104 |
| Processing and prediction code | `decomp/stage1/analyze2.py` | 52d3f0d6a316947c7d823d76dbbe7fb1b40195ac05c793734cf5c0298eaebe91 |
| Processing and prediction code | `decomp/stage1/analyze3.py` | c34649196960583087896def6595c42639c4c45e982f1b61c1ec023af7e635f5 |
| Judgment programs | `decomp/remeasure/judge_tong.py` | ac2c541b437b5c33aca6e5ba62d69f6ba86aefa2214f3382c7f2991c552b888b |
| Judgment programs | `decomp/remeasure/judge_wong.py` | b3e43d17449ea77221117e956f5fc812d2d1b62edb17bb1f4e71dba2b5e45cf6 |
| Judgment programs | `decomp/remeasure/judge_v4.py` | f751a5fc02805b7574da1c723a01ae5bfca4974270e15edcc6ecf65f7c7cb018 |
| Judgment programs | `decomp/remeasure/judge_onh.py` | 48da3c736e14bec108aebaee5ffe5b921743a89c7f64834ef5c897a40686cdd7 |
| Judge code called by V4 | `decomp/s4/validate.py` | 808096649695172ffc53e5f94fa736ad8c252b1542931e616cd1baef2b89739d |
| Judge code called by V4 | `decomp/s4/s4cfg.py` | c31892986ed7b04b8d45353633e1be670b0d079d935ab884385791eeda76c75d |
| Judge code called by V4 | `decomp/s4/plan.py` | 4a8c6ec590cd1dbf28e5d5ea0ba4b2dc4377f99682dc09c81cd1bdcb5259c314 |
| Frozen judge records | `decomp/results/s4/train_freeze.json` (weights) | 97a6a5705904640578c73423c614d9456cfc7aad018a4a931a6e31d30962dcd8 |
| Frozen judge records | `decomp/results/s4/sim_hashes.json` (SHA-256 record of the test data) | 73b2ddd237bcbadc34210cb7acfaf18da2c4fadcab66c07f6f9633de5121d628 |
| Frozen judge records | `decomp/results/s4/miss.json` (misspecification thresholds) | f87cf7f4d4b71598fe31220f7819f1e038f7638974508603cba3e421dcfd5bff |
| Frozen judge records | `decomp/results/s4/cards_c3.json` (cell classification, C3) | f9059614ffa8aec87e8eb41938e0b4df739ad94b538793a24f5eb30d1a9aaea5 |
| Frozen judge records | `decomp/results/s4/cards_c7.json` (cell classification, C7) | 28f04e8085c41da4837c73b5cad0e8f59151f6d754e215d295c063a6fc1178df |
| Environment locks | `uv.lock` | 2c9081d8c8b600912f60fc6b7b9e6d3f2fb00e97f1205a8857fc9ea66add6f16 |
| Environment locks | `pyproject.toml` | ff6bba08e9a24cd88df14fde5f0643d8e61d99338009644569d0ed35470e8373 |
| Environment locks | `envs/reference/uv.lock` | bb0739127d6ba65878b7788cb7519d0a83bdb9a08645fbda6b19c51b4f8ad1a5 |
| Environment locks | `envs/reference/pyproject.toml` | 33209ec644a6b91ab2c02e02a8eac6e8132bdc3c7db5788b7bd47ad5eb2b6b90 |

## (c) Frozen judge products (weights frozen at 85eb8b0; all in the Zenodo deposit), checked by `judge_v4.check_frozen` against these values before use. The same values are recorded in `decomp/results/s4/train_freeze.json` (models, sha256_obj) and `decomp/results/s4/sim_hashes.json` (files).

| File | SHA-256 | Bytes |
|---|---|---|
| `data/s4/models/c3_0.obj.pt` | afbdf0f2337760e04604cf1e527a5f71dc3a0e25553481e0aec225e84844ce47 | 2,766,409 |
| `data/s4/models/c3_1.obj.pt` | 17cc88959f1663a3f81949b1040713d4b3831a1a235e2923cf2d0a2cf9fa7340 | 2,766,409 |
| `data/s4/models/c3_2.obj.pt` | e85d74ee0903ac0fbd548bfa56a655630a1f436dea067be97fe8cc357c38adcd | 2,766,409 |
| `data/s4/models/c3_3.obj.pt` | db16b1158f87e78dafedc087dd6b9f9191152c3f4c55c46d14f59f7e409c79a0 | 2,766,409 |
| `data/s4/models/c3_4.obj.pt` | f18824037397e2b4cadc2081ec17014097759ff24ca4d262a3dbbb975d273956 | 2,766,409 |
| `data/s4/models/c7_0.obj.pt` | 76c2715e8097c7aaad76c3838719e2a2503bb5b21c396a4d62d2182a4c0612b7 | 2,770,633 |
| `data/s4/models/c7_1.obj.pt` | d45971c99f055ac890dc723c892cd472f4c5428a24e628cf572f05110b5dfa4a | 2,770,633 |
| `data/s4/models/c7_2.obj.pt` | 9d19218bac2df0d2d758b267fd10791727f6f9895580a98f3ef132845b19ddc2 | 2,770,633 |
| `data/s4/models/c7_3.obj.pt` | a6846ad049a88706d08208a55ebaedc6b6946533945ac89882b98ac2eb767058 | 2,770,633 |
| `data/s4/models/c7_4.obj.pt` | f41550389c5727d46f48a735906b089b33cad74119f7e56247677928348cd029 | 2,770,633 |
| `data/s4/models/q_0.obj.pt` | 7d589063699feac046c79e40aac88998f9358aff27918f3dbb13c8f9cbf9b775 | 1,883,929 |
| `data/s4/sim_test.npz` | c4ce012aa8066597e206faf124bc1fec146237687e57170d1a3488d24c899625 | 9,163,672 |
| `data/s4/plan_test.npz` | f6e014478f532536229b0d77052c05f110895ce7ce90bcfb115a0a208666fd02 | 10,726,394 |

## (d) FSL 6.0.7.23: the 12 executables used, with the full SHA-256 written into the access record at the start of processing (entry of 2026-10-27T17:35:55Z) and the sizes given in the fixed document.
| Executable | SHA-256 | Bytes |
|---|---|---|
| eddy_cpu | a1858f454d77e7869b829d9e84a3917fc782e27657d63511616e3c865a5ed577 | 3,924,848 |
| topup | 0d8ce02a76b99a9d882723e0b6a619b9aa9ab27c786e24910f1c8421903a888a | 109,952 |
| bet | 16b852dc6c5b256fc059c4dd46ecad16550742d80b13a1fac91ed8982ab275b1 | 13,831 |
| dtifit | 5daab302fe5cd058e2818ac1c7e49614c1065224dc84d7cd685fd072125f7eb4 | 485,520 |
| flirt | 36f248e35fda52012190f08606f5f53b5f3f5afd31f3ba837ef8bf5b23a4cb52 | 532,760 |
| fnirt | d08c80d7d0c1acd36b8fc918b86b3d9abf103c4416127e76b8a199139a6bcc79 | 1,061,016 |
| applywarp | f6e0198e501af495e604357d636739c3b42a9fe1d8630e8141e19665e8507a62 | 367,232 |
| fslmaths | 6a70bc3c1d3c12fef8595fb33fc8b2097baebee18cdec7c5e077991c640dc3b8 | 1,111,160 |
| fslstats | 1af365b8a0b2c394bf7139bfbd3260e4c080e9621922f1be49b9f2c1d87b8187 | 233,048 |
| fslroi | e8915edaf9b2c025b248e8ef1813546e4372101dc6414a2226137aa32b405016 | 91,936 |
| fslmerge | 2f4b3f29a91c9936ce8bd17bcc2fcb88b7b8f65404833a21f3c411fc4ba7f249 | 119,528 |
| imcp | 419a31f1b7737b40f0430576bd35210b5d5438f7b27a0f6dcf01d74d3b4db56a | 250 |

*Notes.* All 12 matched the 64-digit values recorded outside the repository before processing; the fixed document's abbreviated SHA-256 of topup ended in "3888a" where that of the executable ends in "a888a" (deviation (i), Supplementary Table S1 (a)).
The version file `fslversion` (SHA-256 c22541346e51f53fab9f84a1e2a2427c59b2615420775c71a8c1f701bae1383e) was recorded in the same entry.

## (e) Order of the analysis steps

Steps in time order. Source: the commit or, where no commit marks the step, the record that holds its time. Earlier versions drew this list as Supplementary Fig. S1 with `decomp/manuscript/figures/src/fig_supp_S1.py`, which reads each time from Git and the records and stops if a time stated in the manuscript or the supplementary material differs from the record.

| Step | Commit or record | Time (KST) |
|---|---|---|
| **(a) Judge and real-data rules** | | |
| Judge specification, first commit | 5029f4b | 2026-10-25 03:48:19 |
| Weights frozen | 85eb8b0 | 2026-10-25 19:13:31 |
| V4 rules fixed | b625516 · b886d76 | 2026-10-26 00:17:46 · 01:00:53 |
| Validation results first opened | judge report (development repository) | 2026-10-26 09:04:55 |
| V4 open points fixed in code | c8825e2 | 2026-10-27 19:11:18 |
| Pipeline fixed | 08a189c | 2026-10-28 00:00:25 |
| Real-data rules released | 5f95b18 | 2026-10-28 01:56:00 |
| **(b) Real-data processing and first judgments** | | |
| Processing started | access record | 2026-10-28 02:35:55 |
| Processing ended | 867b56c | 2026-10-30 18:07:02 |
| Wongkornchaovalit judgment | access record | 2026-10-30 18:11:10 |
| Tong judgment | access record | 2026-10-31 12:56:44 |
| V4, C3 arm result | access record | 2026-10-31 19:28:09 |
| **(c) C7 arm, pair tests, verdict commit** | | |
| V4, C7 arm stopped (no result) | handoff record (development repository) | 2026-11-01 08:22:55 |
| ON-Harmony pair tests | access record | 2026-11-01 14:07:59 |
| Verdicts committed | 5632626 | 2026-11-01 15:20:30 |
| C7 arm rerun (auxiliary) | access record | 2026-11-03 00:12:50 |

(a) The judge specification, the frozen weights and the core V4 rules (Appendix A of the specification) were committed before the judge's validation results were first opened; points the specification left open were fixed in the judge code afterwards (c8825e2; Supplementary Methods S-M7.2). Later amendments of the specification before the opening are not listed. The real-data pipeline and rules were committed before processing began. (b) Processing and the first judgments; each judgment is listed at the time its result entered the access record. (c) The C7 arm of V4 stopped on a numerical error without a result (Supplementary Table S1 (a), deviation (iv)); the ON-Harmony pair tests were run and the verdicts committed (5632626). The C7 arm was then rerun as an auxiliary report, with its reading fixed at 2026-11-01 14:53:03 KST, before the verdict commit; it enters no verdict (Supplementary Table S12). Post hoc and exploratory analyses are not listed; their definition and first-computation times are in Supplementary Table S6.
