#!/bin/bash
# §170 appendix — runs the cuDIMOT NODDI_Watson cascade **sequentially without fsl_sub**.
#
# Pipeline_NODDI_Watson.sh queues the 4 stages with fsl_sub, and in an environment without a cluster the
# fsl_sub rejects an empty `-j ""`, which breaks the passing of options. Here the **same option string**
# Built as it is, then split → fit → merge are called directly. The computation path is identical.
#
# Usage: cudimot_cascade.sh <subjdir>     (needs $subjdir/{data,nodif_brain_mask,bvals,bvecs})
set -e
subjdir=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
model=NODDI_Watson; out=${subjdir}.${model}; bin=${CUDIMOT}/bin
step1=GridSeach; step2=FitFractions  # ← keeps the spelling of the pipeline original (typos included)

rm -rf "$out"
mkdir -p "$out"/{diff_parts,logs,Dtifit,$step1/diff_parts,$step1/logs,$step2/diff_parts,$step2/logs}
${FSLDIR}/bin/imcp "$subjdir"/nodif_brain_mask "$out"

# CFP = bvecs, bvals
cp "$subjdir"/bvecs "$subjdir"/bvals "$out"
printf '%s\n%s\n' "$out/bvecs" "$out/bvals" > "$out"/CFP

# FixP = S0 (mean of the b ≤ 50 volumes)
mkdir -p "$out"/temporal; pos=0
for i in $(cat "$subjdir"/bvals); do
  if [ "${i%%.*}" -le 50 ]; then ${FSLDIR}/bin/fslroi "$subjdir"/data "$out"/temporal/volume_$pos $pos 1; fi
  pos=$((pos+1))
done
${FSLDIR}/bin/fslmerge -t "$out"/temporal/S0s "$out"/temporal/volume*
${FSLDIR}/bin/fslmaths "$out"/temporal/S0s -Tmean "$out"/S0
rm -rf "$out"/temporal
echo "$out/S0" > "$out"/FixP

opts="--bi=1000 --nj=1250 --se=25 --data=$subjdir/data --maskfile=$out/nodif_brain_mask --forcedir --CFP=$out/CFP"

run_stage () {  # $1=name  $2=working directory  $3...=extra options
  local nm=$1 wd=$2; shift 2
  local o="$opts --outputdir=$wd --partsdir=$wd/diff_parts --FixP=$out/FixP $*"
  echo "── $nm"
  $bin/split_parts_$model       $o --idPart=0 --nParts=1 --logdir=$wd/logs/pre
  $bin/$model                   $o --idPart=0 --nParts=1 --logdir=$wd/logs/fit
  $bin/merge_parts_$model       $o --idPart=0 --nParts=1 --logdir=$wd/logs/post
}

# 0) Dtifit → initial th·ph values
$bin/Run_dtifit.sh "$subjdir" "$out" "$bin" > "$out"/logs/dtifit.log 2>&1
DTI=$out/Dtifit

# 1) GridSearch (fiso·fintra·kappa), th·ph fixed at the DTI values
printf '\n\n\n%s\n%s\n' "$DTI/dtifit_V1_th.nii.gz" "$DTI/dtifit_V1_ph.nii.gz" > "$out/$step1"/InitializationParameters
{ echo "search[0]=(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)"
  echo "search[1]=(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)"
  echo "search[2]=(1,2,3,4,5,6,7,8)"; } > "$out/$step1"/GridSearch
run_stage "GridSearch" "$out/$step1" --gridSearch="$out/$step1/GridSearch" --no_LevMar \
          --init_params="$out/$step1/InitializationParameters" --fixed=3,4

# 2) LM — fractions only
printf '%s\n%s\n%s\n%s\n%s\n' "$out/$step1/Param_0_samples" "$out/$step1/Param_1_samples" \
       "$out/$step1/Param_2_samples" "$DTI/dtifit_V1_th.nii.gz" "$DTI/dtifit_V1_ph.nii.gz" \
       > "$out/$step2"/InitializationParameters
run_stage "LM (분율)" "$out/$step2" --init_params="$out/$step2/InitializationParameters" --fixed=2,3,4

# 3) LM — all parameters
printf '%s\n%s\n%s\n%s\n%s\n' "$out/$step2/Param_0_samples" "$out/$step2/Param_1_samples" \
       "$out/$step1/Param_2_samples" "$DTI/dtifit_V1_th.nii.gz" "$DTI/dtifit_V1_ph.nii.gz" \
       > "$out"/InitializationParameters
run_stage "LM (전 모수)" "$out" --init_params="$out/InitializationParameters"

# 4) finish — OD·dyads
$bin/${model}_finish.sh "$out"
echo "완료: $out"
