#!/usr/bin/env bash
set -Eeuo pipefail
trap 'echo "[ERROR] $0:$LINENO: \"$BASH_COMMAND\" failed" >&2' ERR

LOG_DIR=/tmp/log
mkdir -p "$LOG_DIR"

TS=$(date '+%Y%m%d-%H%M%S')
LOG_OUT="$LOG_DIR/run-$TS.log"
LOG_ERR="$LOG_DIR/run-$TS-error.log"

# tee stdout → console + $LOG_OUT; tee stderr → console + $LOG_ERR
exec > >(tee -a "$LOG_OUT") 2> >(tee -a "$LOG_ERR" >&2)

# finally, run whatever was passed in (your R script + its args)
exec "$@"
