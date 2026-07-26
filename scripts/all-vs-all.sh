#!/bin/bash

if [ $# -lt 2 ]; then
  echo "Usage: $0 <fasta_dir> <out_dir>"
  exit 1
fi

FASTA_DIR=$1
OUT_DIR=$2

if [ ! -d "$FASTA_DIR" ]; then
  echo "Error: Folder containing fasta files not found: $FASTA_DIR"
  exit 1
fi

mkdir -p "$OUT_DIR"

run_diamond() {
  FASTA_FILE=$1
  CLUSTER_NAME=$(basename "$FASTA_FILE" .fasta)

  diamond makedb \
    --in "$FASTA_FILE" \
    --db "$OUT_DIR/${CLUSTER_NAME}" \
    --quiet

  diamond blastp \
    --query "$FASTA_FILE" \
    --db "$OUT_DIR/${CLUSTER_NAME}.dmnd" \
    --out "$OUT_DIR/${CLUSTER_NAME}.tsv" \
    --outfmt 6 qseqid sseqid pident length qlen slen positive evalue bitscore \
    -k 0 \
    --no-self-hits \
    --fast \
    --quiet

  rm "$OUT_DIR/${CLUSTER_NAME}.dmnd"
}

export -f run_diamond

echo "[$(date +%H:%M:%S)] Calculating DIAMOND all-vs-all alignments..."

for file in "$FASTA_DIR"/*.{fasta,fa,faa}; do
  [ -e "$file" ] || continue

  run_diamond "$file"
done

echo "[$(date +%H:%M:%S)] Alignments completed."
