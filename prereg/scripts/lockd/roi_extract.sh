#!/usr/bin/env bash
# Material of 10_roi_s — skeleton projection + extraction of **27 ROI units**. ⛔ Runs inside the gate.
# ⛔ [2026-10-05 · just before 30 · JUDGE_PATHS §35] The input of stage 30 (the remaining 101 sessions) is also built in the **same stage** — only the session list is an argument (default repeat54.tsv).
# Canonical 48-label atlas → 27 ROI units by left-right mean (§212-1). Also outputs vps (number of skeleton voxels).
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)  # Parent of the repository (absolute path dropped, stage 1)
export FSLDIR=$ROOT/fsl; . $FSLDIR/etc/fslconf/fsl.sh; export PATH=$FSLDIR/bin:$PATH
ATLAS=$ROOT/dmri_cov/prereg/artifacts/atlas/JHU-ICBM-labels-1mm.nii.gz
SKEL=$FSLDIR/data/standard/FMRIB58_FA-skeleton_1mm.nii.gz
OUT=$ROOT/data/roi; mkdir -p "$OUT"

# The skeleton mask and the per-ROI vps do not depend on the session, so they are made **only once**
if [ ! -f "$OUT/vps.txt" ]; then
  fslmaths "$SKEL" -thr 2000 -bin "$OUT/skel_mask"
  fslmaths "$ATLAS" -mas "$OUT/skel_mask" "$OUT/atlas_skel"
  fslstats -K "$ATLAS" "$OUT/skel_mask" -V | awk '{print $1}' > "$OUT/vps.txt"
fi

n=0
while IFS=$'\t' read -r sub ses; do
  W=$ROOT/data/timing/${sub}_${ses}; S=$OUT/${sub}_${ses}
  [ -f "$S.roi.tsv" ] && { n=$((n+1)); continue; }
  mkdir -p "$S.tmp"
  for m in cudimot.NODDI_Watson/mean_fintra cudimot.NODDI_Watson/mean_fiso cudimot.NODDI_Watson/OD; do
    b=$(basename $m)
    applywarp -i "$W/$m" -r "$FSLDIR/data/standard/FMRIB58_FA_1mm" -w "$W/fa2std_warp" -o "$S.tmp/std_$b" --interp=trilinear
    fslmaths "$S.tmp/std_$b" -mas "$OUT/skel_mask" "$S.tmp/skel_$b"
    fslstats -K "$ATLAS" "$S.tmp/skel_$b" -M > "$S.tmp/roi_$b.txt"
  done
  paste "$S.tmp/roi_mean_fintra.txt" "$S.tmp/roi_mean_fiso.txt" "$S.tmp/roi_OD.txt" > "$S.roi.tsv"
  rm -rf "$S.tmp"
  n=$((n+1)); printf '\r  %2d/54 %s' "$n" "${sub}_${ses}"
done < "${1:-$ROOT/data/repeat54.tsv}"
echo; echo "완료 — $(ls $OUT/*.roi.tsv | wc -l) 세션 · 48 라벨 × 3 축"
