#!/usr/bin/env bash
# §27-11 plan-version known-answer test — when the bank is finished: rename → input copy (re-access record) → kat → fill the fixed hash slots → result. CPU only · 0 image reads.
set -euo pipefail
export CUDA_VISIBLE_DEVICES= PYTHONDONTWRITEBYTECODE=1
PY=<REPO>/.venv/bin/python
H=<WORK>/stage1_work/drafts/x15/b2000
RR=<WORK>/stage1_work/drafts/x15/ledger/reread_log.py
ORIG=<WORK>/remeasure/explore/b2000; OUT=$ORIG/plan_rigid6
while pgrep -f "^<REPO>/.venv/bin/python -I -B b2k_sim_rigid6.py bank" > /dev/null; do sleep 20; done
grep -q "^600/600" "$OUT/bank.log" || { echo "은행이 600/600 으로 끝나지 않음"; tail -3 "$OUT/bank.log"; exit 1; }
mv "$OUT/kat_bank_OXF1PRI_OXF2PRI_OXF3TRI.jsonl" "$OUT/kat_bank_siemens.jsonl"
for f in read_sessions.jsonl pure_curve_N2_50.json kat_bank_PG.jsonl kat_summary_PG.json; do
  $PY -I -B "$RR" "$ORIG/$f" "§27-11 계획서판 정답 시험 — plan_rigid6/ 에 사본(요약 사본의 입력 · 필립스 · GE 둘은 고정 해시 검사 때문 · 값은 원판과 같음)"
  cp -p "$ORIG/$f" "$OUT/$f"
done
cd "$H"
nice -n 19 $PY -I -B b2k_summary_rigid6.py kat OXF1PRI OXF2PRI OXF3TRI > "$OUT/kat.log" 2>&1
B=$(sha256sum "$OUT/kat_bank_siemens.jsonl" | cut -d' ' -f1); K=$(sha256sum "$OUT/kat_summary_siemens.json" | cut -d' ' -f1)
sed -i "s/@@BANK@@/$B/; s/@@KSUM@@/$K/" b2k_summary_rigid6.py
grep -q "@@" b2k_summary_rigid6.py && { echo "고정 해시 자리가 남음"; exit 1; }
nice -n 19 $PY -I -B b2k_summary_rigid6.py result > "$OUT/result.log" 2>&1
echo "끝 $(date '+%F %T') bank $B kat $K"
