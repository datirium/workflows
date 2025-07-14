#!/bin/bash
set -euo pipefail

declare -A TOOL_ENVS=(
  [bowtie2]=aligner
  [STAR]=aligner
  [trim_galore]=trim_galore
  [samtools]=samtools
  [preseq]=preseq
  [macs3]=macs3
  [bedtools]=bedtools
  [bedGraphToBigWig]=ucsc
  [fastx_quality_stats]=fastx
  [file]=utils
  [unzip]=utils
  [make]=compiler
  [qmake]=compiler
  [python3]=python
  [p_pandas]=python
  [p_yaml]=python
  [Rscript]=r_base
  [r_modules]=r_base
  [r_future]=r_base
  [r_argparse]=r_base
  [r_knitr]=r_base
  [r_forcats]=r_base
  [r_UpSetR]=r_base
  [r_ggupset]=r_base
  [r_ChIPseeker]=r_base
  [r_txdbmaker]=r_base
  [iaintersect]=none
  [extract_fastq.sh]=none
  [test_extract_fastq.sh]=none
  [collect_stats_dna.sh]=none
)

MISSING=()
declare -A ENV_SIZES=()

for TOOL in "${!TOOL_ENVS[@]}"; do
  ENV_NAME="${TOOL_ENVS[$TOOL]}"
  if [[ "$ENV_NAME" == "none" ]]; then
    if ! command -v "$TOOL" >/dev/null 2>&1; then
      MISSING+=("$TOOL (global)")
    fi
  elif [[ "$TOOL" == r_* ]]; then
    PKG="${TOOL#r_}"
    if ! mamba run -n "$ENV_NAME" Rscript -e "if (!requireNamespace('$PKG', quietly=TRUE)) quit(status=1)" >/dev/null 2>&1; then
      MISSING+=("R package $PKG ($ENV_NAME)")
    fi
  elif [[ "$TOOL" == p_* ]]; then
    PKG="${TOOL#p_}"
    if ! mamba run -n "$ENV_NAME" python3 -c "import $PKG" >/dev/null 2>&1; then
      MISSING+=("Python module $PKG ($ENV_NAME)")
    fi
  else
    if ! mamba run -n "$ENV_NAME" which "$TOOL" >/dev/null 2>&1; then
      MISSING+=("$TOOL ($ENV_NAME)")
    fi
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
  echo "All required tools and R packages are installed and available in their mamba environments."
  exit 0
else
  echo "Missing tools or R packages: ${MISSING[*]}" >&2
  exit 1
fi
