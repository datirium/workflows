#!/bin/bash
set -euo pipefail

declare -A TOOL_ENVS=(
  [bowtie2]=bowtie2
  [trim_galore]=trim_galore
  [samtools]=samtools
  [preseq]=preseq
  [macs3]=macs3
  [bedtools]=bedtools
  [bedGraphToBigWig]=ucsc
  [fastx_quality_stats]=fastx
  [Rscript]=r_base
  [file]=utils
  [unzip]=utils
)

MISSING=()
declare -A ENV_SIZES=()

for TOOL in "${!TOOL_ENVS[@]}"; do
  ENV_NAME="${TOOL_ENVS[$TOOL]}"
  if ! mamba run -n "$ENV_NAME" which "$TOOL" >/dev/null 2>&1; then
    MISSING+=("$TOOL ($ENV_NAME)")
  fi
  if [[ -z "${ENV_SIZES[$ENV_NAME]:-}" ]]; then
    if [ -d "/opt/conda/envs/$ENV_NAME" ]; then
      SIZE=$(du -sh "/opt/conda/envs/$ENV_NAME" 2>/dev/null | awk '{print $1}')
      ENV_SIZES[$ENV_NAME]="$SIZE"
    fi
  fi
done

echo -e "\nMamba environment sizes:" >&2
for ENV in "${!ENV_SIZES[@]}"; do
  echo "  $ENV: ${ENV_SIZES[$ENV]}" >&2
done

if [ ${#MISSING[@]} -eq 0 ]; then
  echo "All required tools are installed and available in their mamba environments."
  exit 0
else
  echo "Missing tools: ${MISSING[*]}" >&2
  exit 1
fi
