#!/bin/bash
set -euo pipefail

usage() {
  echo "Usage: $0 [-o output.fastq] <fastq1> [fastq2 ...]"
  echo "  -o, --output   Output FASTQ file (default: combined.fastq)"
  exit 1
}

OUT="combined.fastq"
FILES=()

# Argument parsing
while [[ $# -gt 0 ]]; do
  case "$1" in
    -o|--output)
      OUT="$2"; shift 2;;
    -h|--help)
      usage;;
    --)
      shift; break;;
    -*)
      echo "Unknown option: $1"; usage;;
    *)
      FILES+=("$1"); shift;;
  esac
done

if [[ ${#FILES[@]} -eq 0 ]]; then
  usage
fi

> "$OUT"  # Truncate output

for FILE in "${FILES[@]}"; do
  echo "Extracting: $FILE"
  if [[ ! -f "$FILE" ]]; then
    echo "Error: input file does not exist: $FILE" >&2
    exit 1
  fi
  case "$(file -b --mime-type "$FILE")" in
    application/gzip)    gunzip -c "$FILE" >> "$OUT" ;;
    application/x-bzip2) bzcat "$FILE" >> "$OUT" ;;
    application/zip)     unzip -p "$FILE" >> "$OUT" ;;
    text/plain)          cat "$FILE" >> "$OUT" ;;
    *)                   echo "Unknown file type: $FILE" >&2; exit 1 ;;
  esac
done