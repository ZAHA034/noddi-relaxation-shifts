#!/bin/bash
# decomp stage 2 run — E5 of SPEC_v2 §11 (dense grid · F2 only · CPU only · no GPU). Can be resumed (DONE marker).
# Usage: run_stage2.sh [number of parallel chunks — default 2]
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
WORK=$(cd "$(dirname "$0")/../../.." && pwd)/data; mkdir -p "$WORK/logs"
NPAR=${1:-2}
echo "$(date -u +%FT%TZ) 2 단계 E5 · GPU none" | tee -a "$WORK/logs/run2.log"
( while true; do echo "$(date -u +%FT%TZ) $(cat /proc/loadavg)"; sleep 30; done ) >> "$WORK/logs/loadavg2.log" & LOADP=$!
trap 'kill $LOADP 2>/dev/null || true' EXIT
python3 - "$HERE" > "$WORK/logs/jobs2.txt" <<'PY'
import sys; sys.path.insert(0, sys.argv[1]); import config as C
for p in ("P2a", "P2b"):
    for i in range(len(C.chunk_rows("E5", p)[1])):
        print("E5", p, i)
PY
echo "$(date -u +%FT%TZ) 덩어리 $(wc -l < "$WORK/logs/jobs2.txt") 개 · 병렬 $NPAR" | tee -a "$WORK/logs/run2.log"
xargs -P "$NPAR" -L 1 bash -c 'bash "'"$HERE"'/run_chunk.sh" "$0" "$1" "$2" none >> "'"$WORK"'/logs/chunks2.log" 2>&1 || echo "⛔ 실패 $0 $1 $2" >> "'"$WORK"'/logs/run2.log"' < "$WORK/logs/jobs2.txt"
echo "$(date -u +%FT%TZ) 끝" | tee -a "$WORK/logs/run2.log"
