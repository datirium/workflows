#!/bin/bash
set -euo pipefail

REQUIRED_TOOLS=(
  bowtie2
  fastqc
  cutadapt
  trim_galore
  perl
  unzip
  bzip2
  file
)

MISSING=()

for TOOL in "${REQUIRED_TOOLS[@]}"; do
  if ! command -v "$TOOL" >/dev/null 2>&1; then
    MISSING+=("$TOOL")
  fi
done

if [ ${#MISSING[@]} -eq 0 ]; then
  echo "All required tools are installed and available in PATH."
  exit 0
else
  echo "Missing tools: ${MISSING[*]}" >&2
  exit 1
fi
