#!/bin/bash
set -euo pipefail

declare -A TOOL_ENVS=(
  [bowtie2]=bowtie2
  [trim_galore]=trim_galore
  [samtools]=samtools
  [file]=utils
  [unzip]=utils
)

MISSING=()

for TOOL in "${!TOOL_ENVS[@]}"; do
  ENV_NAME="${TOOL_ENVS[$TOOL]}"
  if ! mamba run -n "$ENV_NAME" "$TOOL" --help >/dev/null 2>&1; then
    MISSING+=("$TOOL ($ENV_NAME)")
  fi

done

if [ ${#MISSING[@]} -eq 0 ]; then
  echo "All required tools are installed and available in their mamba environments."
  exit 0
else
  echo "Missing tools: ${MISSING[*]}" >&2
  exit 1
fi
