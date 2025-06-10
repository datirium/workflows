#!/bin/bash
# Test script for extract_fastq.sh
# Usage: ./test_extract_fastq.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/data/dm3_chr4/fastq/chr4_100_chrM_25_mapped_reads"
TMP_DIR="$(mktemp -d)"

# Test 1: Single uncompressed FASTQ (upstream)
cp "$DATA_DIR/upstream.fastq" "$TMP_DIR/"
/usr/local/bin/extract_fastq.sh -o "$TMP_DIR/out1.fastq" "$TMP_DIR/upstream.fastq"
diff "$TMP_DIR/upstream.fastq" "$TMP_DIR/out1.fastq"
echo "Test 1 passed: single uncompressed FASTQ (upstream)"

# Test 2: Single gzipped FASTQ (downstream)
cp "$DATA_DIR/downstream.fastq.gz" "$TMP_DIR/"
/usr/local/bin/extract_fastq.sh -o "$TMP_DIR/out2.fastq" "$TMP_DIR/downstream.fastq.gz"
gunzip -c "$TMP_DIR/downstream.fastq.gz" | diff - "$TMP_DIR/out2.fastq"
echo "Test 2 passed: single gzipped FASTQ (downstream)"

# Test 3: Single bzip2 FASTQ (upstream)
cp "$DATA_DIR/upstream.fastq.bz2" "$TMP_DIR/"
/usr/local/bin/extract_fastq.sh -o "$TMP_DIR/out3.fastq" "$TMP_DIR/upstream.fastq.bz2"
bzcat "$TMP_DIR/upstream.fastq.bz2" | diff - "$TMP_DIR/out3.fastq"
echo "Test 3 passed: single bzip2 FASTQ (upstream)"

# Test 4: Multiple files (mixed: uncompressed, gzipped, bzip2)
cp "$DATA_DIR/upstream.fastq" "$TMP_DIR/"
cp "$DATA_DIR/downstream.fastq.gz" "$TMP_DIR/"
cp "$DATA_DIR/upstream.fastq.bz2" "$TMP_DIR/"
/usr/local/bin/extract_fastq.sh -o "$TMP_DIR/out4.fastq" "$TMP_DIR/upstream.fastq" "$TMP_DIR/downstream.fastq.gz" "$TMP_DIR/upstream.fastq.bz2"
cat "$TMP_DIR/upstream.fastq" <(gunzip -c "$TMP_DIR/downstream.fastq.gz") <(bzcat "$TMP_DIR/upstream.fastq.bz2") | diff - "$TMP_DIR/out4.fastq"
echo "Test 4 passed: mixed input (uncompressed, gzipped, bzip2)"

# Test 5: No input files (should fail)
if /usr/local/bin/extract_fastq.sh -o "$TMP_DIR/out6.fastq"; then
  echo "Test 5 failed: script should have failed with no input files" >&2
  exit 1
else
  echo "Test 5 passed: script correctly failed with no input files"
fi

# Test 7: Nonexistent file (should fail)
if /usr/local/bin/extract_fastq.sh -o "$TMP_DIR/out7.fastq" "$TMP_DIR/nonexistent.fastq"; then
  echo "Test 6 failed: script should have failed with nonexistent input file" >&2
  exit 1
else
  echo "Test 6 passed: script correctly failed with nonexistent input file"
fi

# Clean up
echo "All tests passed. Cleaning up."
rm -rf "$TMP_DIR"
