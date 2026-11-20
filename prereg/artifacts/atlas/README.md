# Atlas for step 3

The ROI definition uses the 48-label JHU ICBM-DTI-81 white-matter atlas, from FSL `data_atlases`, tag `FSL6-0-1`, folder `JHU/`. It is not redistributed here: it stays under its owners' terms, which this deposit cannot grant (`THIRD_PARTY_NOTICES.md`). Steps 1 and 2 do not need it. Without it, one test of the label-to-unit mapping is skipped.

A copy with these SHA-256 values is the right one, and this command checks a copy and puts it in this folder:

```
python3 reproduce/get_atlas.py --from <file>
```

| File | Bytes | SHA-256 |
|---|---:|---|
| `JHU-ICBM-labels-1mm.nii.gz` | 47,474 | `fac584ec75ff2a8631710d3345df96733ed87d9bde3387f5b462f8d22914ed69` |
| `JHU-labels.xml` | 4,594 | `2d76ce80d1b0a50dccda2698d5eec55c8984a7f1bb438f79111d67a26fc4dc1c` |

The 50-label atlas of FSL 6.0.5 and later has a different SHA-256 and a different ROI definition, so it is refused.
