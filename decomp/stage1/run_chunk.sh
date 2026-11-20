#!/bin/bash
# decomp stage 1 — one chunk: synthesis → F1 cuDIMOT → F2 AMICO (per fitted d∥) → keep only the fitted values and delete the signal (SPEC_v1 §4 · §7).
# Usage: run_chunk.sh <exp> <prot> <chunk_id> <gpu|none>  — none: F1 skipped (CPU only · user instruction of 2026-10-20)
set -euo pipefail
exp=$1; prot=$2; cid=$3; gpu=$4
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$HERE/../.." && pwd); ROOT=$(cd "$REPO/.." && pwd)
export FSLDIR=$ROOT/fsl; export PATH=$FSLDIR/bin:$PATH; export FSLOUTPUTTYPE=NIFTI_GZ
export CUDIMOT=$ROOT/tools/cudimot/v0.1.0-fsl6.0.7.23-sm89-watson01
WORK=${DECOMP_WORK:-$(cd "$(dirname "$0")/../../.." && pwd)/data}
out=$WORK/fits/$exp/$prot/c$(printf %04d $cid); done=$out/DONE
[ -f "$done" ] && { echo "건너뜀 $out"; exit 0; }
rm -rf "$out"; mkdir -p "$out"
subj=$out/subj
(cd "$REPO/envs/reference" && uv run --frozen python "$HERE/synth.py" "$exp" "$prot" "$cid" "$subj")
mv "$subj/rows.npz" "$out/rows.npz"
# F1 — cuDIMOT cascade (as is) · skipped if gpu=none (CUDA only)
if [ "$gpu" != "none" ]; then
CUDA_VISIBLE_DEVICES=$gpu bash "$REPO/prereg/scripts/cudimot_cascade.sh" "$subj" > "$out/f1.log" 2>&1
(cd "$REPO" && uv run --frozen python - "$subj.NODDI_Watson" "$out/f1.npz" <<'PY'
import sys, numpy as np, nibabel as nib
d, o = sys.argv[1], sys.argv[2]
g = lambda nm: nib.load(f"{d}/{nm}.nii.gz").get_fdata().reshape(-1).astype(np.float32)
np.savez(o, ndi=g("mean_fintra"), fwf=g("mean_fiso"), odi=g("OD"))
PY
)
fi
# F2 — AMICO
dpars="1.7"; [ "$exp" = "E3" ] && dpars="1.7 2.0 2.26 2.5"
for dp in $dpars; do
  (cd "$REPO/envs/reference" && uv run --frozen python "$HERE/amico_fit.py" "$subj" "$dp" 16) > "$out/f2_$dp.log" 2>&1
  mv "$subj/f2_$dp.npz" "$out/"
done
rm -rf "$subj" "$subj.NODDI_Watson"
date -u +%FT%TZ > "$done"
echo "끝 $out"
