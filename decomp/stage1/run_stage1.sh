#!/bin/bash
# decomp stage 1 run — runs the chunks in the order of SPEC_v1 §7 (E1 → E1-TR → E3 → E2a → E2c → E2b). Can be resumed (DONE marker).
# GPU: used only when the user specifies a card (no automatic selection — user instruction of 2026-10-20: no GPU 1 · CPU first). When specified, pmon (1 s) is recorded (R26). loadavg is always recorded.
# Usage: run_stage1.sh [number of parallel chunks — default 2] [gpu number | none — default none: skip F1, F2 only]
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
WORK=$(cd "$(dirname "$0")/../../.." && pwd)/data; mkdir -p "$WORK/logs"
NPAR=${1:-2}
GPU=${2:-none}
echo "$(date -u +%FT%TZ) GPU $GPU — $(nvidia-smi --query-gpu=index,memory.used,utilization.gpu --format=csv,noheader | tr '\n' ';')" | tee -a "$WORK/logs/run.log"
PMON=
if [ "$GPU" != "none" ]; then nvidia-smi pmon -i "$GPU" -d 1 -o T >> "$WORK/logs/pmon_gpu$GPU.log" 2>&1 & PMON=$!; fi
( while true; do echo "$(date -u +%FT%TZ) $(cat /proc/loadavg)"; sleep 30; done ) >> "$WORK/logs/loadavg.log" & LOADP=$!
trap 'kill $PMON $LOADP 2>/dev/null || true' EXIT
python3 - "$HERE" > "$WORK/logs/jobs.txt" <<'PY'
import sys; sys.path.insert(0, sys.argv[1]); import config as C
for e in C.EXPERIMENTS:
    for p in C.PROTOCOLS:
        if any(c["prot"] == p for c in C.cells(e)):
            for i in range(len(C.chunk_rows(e, p)[1])):
                print(e, p, i)
PY
echo "$(date -u +%FT%TZ) 덩어리 $(wc -l < "$WORK/logs/jobs.txt") 개 · 병렬 $NPAR" | tee -a "$WORK/logs/run.log"
xargs -P "$NPAR" -L 1 bash -c 'bash "'"$HERE"'/run_chunk.sh" "$0" "$1" "$2" '"$GPU"' >> "'"$WORK"'/logs/chunks.log" 2>&1 || echo "⛔ 실패 $0 $1 $2" >> "'"$WORK"'/logs/run.log"' < "$WORK/logs/jobs.txt"
echo "$(date -u +%FT%TZ) 끝" | tee -a "$WORK/logs/run.log"
