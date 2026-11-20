#!/usr/bin/env bash
# 1 subject end-to-end **timing** — produces the number that answers «yes/no» to whether 165 sessions run within the deadline.
#
# What is measured: wall clock per stage · separation of GPU wait and computation · **disk per session** · output size.
# ⛔ Only the registered path is used (LAYER2_PLAN): topup → eddy → dtifit → cuDIMOT watson01 → FA alignment → ROI 27.
# ⛔ Only sessions of the same subject and same scanner are processed — a second-scanner value is not made (P9-9 order guarantee).
set -euo pipefail

SUB=${1:?sub-XXXXX}; SES=${2:?ses-XXXXXXX}
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)  # Parent of the repository (absolute path dropped, stage 1)
SEALED=$ROOT/source_data/onharmony_sealed/ds004712/$SUB/$SES/dwi
META=$ROOT/source_data/onharmony_meta/ds004712/$SUB/$SES/dwi
WORK=$ROOT/data/timing/${SUB}_${SES}
export FSLDIR=$ROOT/fsl
export CUDIMOT=$ROOT/tools/cudimot/v0.1.0-fsl6.0.7.23-sm89-watson01
export CUDA_VISIBLE_DEVICES=${GPU:-0}
. $FSLDIR/etc/fslconf/fsl.sh
export PATH=$FSLDIR/bin:$PATH
TIMES=$WORK/timings.tsv
MINFREE_GB=20

mkdir -p "$WORK"; cd "$WORK"
: > "$TIMES"

free_gb() { df -BG --output=avail "$WORK" | tail -1 | tr -dc '0-9'; }
gpu_used() { nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i "${CUDA_VISIBLE_DEVICES}" 2>/dev/null | head -1; }

step() {  # step <name> <command...>
  local name=$1; shift
  local free; free=$(free_gb)
  if [ "$free" -lt "$MINFREE_GB" ]; then
    echo "⛔ 여유 디스크 ${free} GB < ${MINFREE_GB} GB — 중단" | tee -a "$TIMES"; exit 3
  fi
  local g0 t0; g0=$(gpu_used); t0=$(date +%s.%N)
  echo "── [$name] 시작 $(date -u +%H:%M:%SZ) · 여유 ${free} GB · GPU 사용 ${g0} MiB"
  "$@" > "$WORK/log_${name}.txt" 2>&1 || { echo "⛔ [$name] 실패 — log_${name}.txt 참조"; tail -5 "$WORK/log_${name}.txt"; exit 4; }
  local t1 du_b g1; t1=$(date +%s.%N); du_b=$(du -sb "$WORK" | cut -f1); g1=$(gpu_used)
  printf '%s\t%.1f\t%s\t%s\t%s\n' "$name" "$(echo "$t1-$t0" | bc)" "$du_b" "$g0" "$g1" >> "$TIMES"
  printf '   [%s] %.1f 초 · 작업디렉터리 %.2f GiB · GPU %s→%s MiB\n' "$name" "$(echo "$t1-$t0" | bc)" "$(echo "scale=2; $du_b/1073741824" | bc)" "$g0" "$g1"
}

prep() {
  ln -sf "$SEALED/${SUB}_${SES}_dir-AP_dwi.nii.gz" AP.nii.gz
  ln -sf "$SEALED/${SUB}_${SES}_dir-PA_dwi.nii.gz" PA.nii.gz
  cp "$META/${SUB}_${SES}_dir-AP_dwi.bval" bvals
  cp "$META/${SUB}_${SES}_dir-AP_dwi.bvec" bvecs
  local trt; trt=$(python3 -c "import json;print(json.load(open('$META/${SUB}_${SES}_dir-AP_dwi.json'))['TotalReadoutTime'])")
  # Take the first b0 of AP (b≤50) and the first volume of PA and concatenate them
  local idx0; idx0=$(python3 -c "print([i for i,b in enumerate(open('bvals').read().split()) if float(b)<=50][0])")
  fslroi AP.nii.gz b0_AP "$idx0" 1
  fslroi PA.nii.gz b0_PA 0 1
  fslmerge -t b0s b0_AP b0_PA
  printf '0 -1 0 %s\n0 1 0 %s\n' "$trt" "$trt" > acqparams.txt
  python3 -c "print(' '.join(['1']*len(open('bvals').read().split())))" > index.txt
}

step 01_prep            prep
step 02_topup           topup --imain=b0s --datain=acqparams.txt --config="$FSLDIR/etc/flirtsch/b02b0.cnf" --out=topup_res --iout=hifi_b0 --fout=field
step 03_mask            bash -c 'fslmaths hifi_b0 -Tmean hifi_b0_mean && bet hifi_b0_mean nodif_brain -m -f 0.25'
step 04_eddy            eddy_cuda --imain=AP.nii.gz --mask=nodif_brain_mask --acqp=acqparams.txt --index=index.txt --bvecs=bvecs --bvals=bvals --topup=topup_res --out=eddy_out --verbose
step 05_dtifit          dtifit -k eddy_out.nii.gz -o dti -m nodif_brain_mask -r eddy_out.eddy_rotated_bvecs -b bvals
step 06_cudimot_prep    bash -c 'mkdir -p cudimot && imcp eddy_out cudimot/data && imcp nodif_brain_mask cudimot/nodif_brain_mask && cp eddy_out.eddy_rotated_bvecs cudimot/bvecs && cp bvals cudimot/bvals'
step 07_cudimot         bash "$ROOT/dmri_cov/prereg/scripts/cudimot_cascade.sh" "$WORK/cudimot"
step 08_reg_fa          bash -c "flirt -in dti_FA -ref \$FSLDIR/data/standard/FMRIB58_FA_1mm -omat fa2std.mat -out fa2std_lin -dof 12 && fnirt --in=dti_FA --aff=fa2std.mat --config=FA_2_FMRIB58_1mm --cout=fa2std_warp --iout=fa2std_nl"

echo; echo "=== 요약 ==="; cat "$TIMES"
