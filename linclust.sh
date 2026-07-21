#!/bin/bash

COVERAGE=0.8    # -c default 0.8
COVERAGE_MODE=1 # --cov-mode, 1: coverage of target

FASTA_FILE=$1
MIN_SEQ_ID=$2

if [ $# -lt 2 ]; then
  echo usage: $0 fasta_file min_seq_id
  exit 1
fi

if [ ! -f "$FASTA_FILE" ]; then
  echo "Input FASTA file '$FASTA_FILE' not found."
  exit 1
fi

if [ -z "$MIN_SEQ_ID" ] || ! [[ "$MIN_SEQ_ID" =~ ^(0(\.[0-9]+)?|1(\.0+)?)$ ]]; then
  echo "Given min-seq-id ($MIN_SEQ_ID) is not a valid number between 0.0 and 1.0."
  exit 1
fi

FILENAME=$(basename "$FASTA_FILE")
BASENAME=$(echo "$FILENAME" | sed -E 's/\.(fasta|fa|faa)(\.gz)?$//')
SAFE_SEQ_ID=$(LC_NUMERIC=C printf "%.2f" "$MIN_SEQ_ID" | sed 's/\.//')

OUT_DIR="linclust_results"
mkdir -p "$OUT_DIR"
TMP_DIR="tmp"
mkdir -p "$TMP_DIR"

PREFIX="$OUT_DIR/${BASENAME}_${SAFE_SEQ_ID}"

echo "Running MMSeqs easy-linclust with $MIN_SEQ_ID minimum sequence identity, $COVERAGE coverage and coverage mode $COVERAGE_MODE..."

mmseqs easy-linclust "$FASTA_FILE" "$PREFIX" "$TMP_DIR" \
  -c "$COVERAGE" \
  --cov-mode "$COVERAGE_MODE" \
  --min-seq-id "$MIN_SEQ_ID"

rm -r "$TMP_DIR"

echo "Results successfully generated in $OUT_DIR"
