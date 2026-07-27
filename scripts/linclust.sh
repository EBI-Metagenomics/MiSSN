#!/bin/bash

MIN_COVERAGE=0.8 # -c, default 0.8
COVERAGE_MODE=1  # --cov-mode, 1: coverage of target

if [ $# -lt 4 ]; then
  echo "Usage: $0 <fasta_file> <min_seq_id> <min_coverage> <out_dir>"
  exit 1
fi

FASTA_FILE=$1
MIN_SEQ_ID=$2
MIN_COVERAGE=$3
OUT_DIR="${4%/}/linclust_results"

if [ ! -f "$FASTA_FILE" ]; then
  echo "Error: Input FASTA file '$FASTA_FILE' not found."
  exit 1
fi

if [ -z "$MIN_SEQ_ID" ] || ! [[ "$MIN_SEQ_ID" =~ ^(0(\.[0-9]+)?|1(\.0+)?)$ ]]; then
  echo "Error: Given min-seq-id ($MIN_SEQ_ID) is not a valid number between 0.0 and 1.0."
  exit 1
fi

FILENAME=$(basename "$FASTA_FILE")
BASENAME=$(echo "$FILENAME" | sed -E 's/\.(fasta|fa|faa)(\.gz)?$//')

mkdir -p "$OUT_DIR"
TMP_DIR="tmp"
mkdir -p "$TMP_DIR"

PREFIX="$OUT_DIR/${BASENAME}"

echo "[$(date +%H:%M:%S)] Running MMSeqs easy-linclust with $MIN_SEQ_ID minimum sequence identity, $MIN_COVERAGE coverage and coverage mode $COVERAGE_MODE..."

mmseqs easy-linclust "$FASTA_FILE" "$PREFIX" "$TMP_DIR" \
  -c "$MIN_COVERAGE" \
  --cov-mode "$COVERAGE_MODE" \
  --min-seq-id "$MIN_SEQ_ID" \
  -v 0

rm -r "$TMP_DIR"

echo "[$(date +%H:%M:%S)] Results successfully generated in '$OUT_DIR'."
