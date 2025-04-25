#!/usr/bin/env bash
set -Eeuo pipefail

# 1) catch any error and report line/command
trap 'echo "[ERROR] $0 failed on line $LINENO: \"$BASH_COMMAND\"" >&2' ERR  # :contentReference[oaicite:0]{index=0}

# 2) prepare log directory
LOG_DIR=/tmp/entrylog
mkdir -p "$LOG_DIR"

# 3) redirect stdout/stderr to both console and log files
TS=$(date '+%Y%m%d-%H%M%S')
exec > >(tee -a "$LOG_DIR/app_$TS.log") 2> >(tee -a "$LOG_DIR/app-error_$TS.log" >&2)

echo "[INFO] Starting at $(date)"

Rscript /usr/local/bin/run_deseq_for_spikein.R "$@"

# 4) run each script in turn
# for SCRIPT in /usr/local/bin/*; do
#   echo "[INFO] Running $SCRIPT"
#   case "$SCRIPT" in
#     *.R)   Rscript "$SCRIPT"    ;;
#     *)     echo "[WARN] Skipping $SCRIPT" ;;
#   esac
# done

echo "[INFO] All done at $(date)"