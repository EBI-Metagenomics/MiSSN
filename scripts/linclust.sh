#!/bin/bash

COVERAGE_MODE=1 # --cov-mode, 1: coverage of target

while getopts "f:i:c:o:" option; do
  case $option in
  f) FASTA_FILE=$OPTARG ;;
  i) MIN_SEQ_ID=$OPTARG ;;
  c) MIN_COVERAGE=$OPTARG ;;
  o) OUT_DIR=$OPTARG ;;
  *)
    echo "Usage: $0 -f <fasta_file> -i <min_seq_id> -c <min_coverage> -o <out_dir>" >&2
    exit 1
    ;;
  esac
done

if [ -z "$FASTA_FILE" ] || [ -z "$MIN_SEQ_ID" ] || [ -z "$MIN_COVERAGE" ] || [ -z "$OUT_DIR" ]; then
  echo "Usage: $0 -f <fasta_file> -i <min_seq_id> -c <min_coverage> -o <out_dir>" >&2
  exit 1
fi

OUT_DIR="${OUT_DIR%/}/linclust_results"

if [ ! -f "$FASTA_FILE" ]; then
  echo "Error: Input FASTA file '$FASTA_FILE' not found." >&2
  exit 1
fi

if ! [[ "$MIN_SEQ_ID" =~ ^(0(\.[0-9]+)?|1(\.0+)?)$ ]]; then
  echo "Error: Given min-seq-id ($MIN_SEQ_ID) is not a valid number between 0.0 and 1.0." >&2
  exit 1
fi

if ! [[ "$MIN_COVERAGE" =~ ^(0(\.[0-9]+)?|1(\.0+)?)$ ]]; then
  echo "Error: Given min-coverage ($MIN_COVERAGE) is not a valid number between 0.0 and 1.0." >&2
  exit 1
fi

FILENAME=$(basename "$FASTA_FILE")
BASENAME=$(echo "$FILENAME" | sed -E 's/\.(fasta|fa|faa)(\.gz)?$//')

mkdir -p "$OUT_DIR"
TMP_DIR="tmp"
mkdir -p "$TMP_DIR"

PREFIX="$OUT_DIR/${BASENAME}"

echo "[$(date +%H:%M:%S)] Running MMseqs2 easy-linclust with $MIN_SEQ_ID minimum sequence identity, $MIN_COVERAGE coverage and coverage mode $COVERAGE_MODE..."

mmseqs easy-linclust "$FASTA_FILE" "$PREFIX" "$TMP_DIR" \
  -c "$MIN_COVERAGE" \
  --cov-mode "$COVERAGE_MODE" \
  --min-seq-id "$MIN_SEQ_ID" \
  -v 0

rm -r "$TMP_DIR"

echo "[$(date +%H:%M:%S)] Results successfully generated in '$OUT_DIR'."
