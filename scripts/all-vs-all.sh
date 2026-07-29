#!/bin/bash

while getopts "f:o:" option; do
  case $option in
  f) FASTA_DIR=$OPTARG ;;
  o) OUT_DIR=$OPTARG ;;
  *)
    echo "Usage: $0 -f <fasta_dir> -o <out_dir>" >&2
    exit 1
    ;;
  esac
done

if [ -z "$FASTA_DIR" ] || [ -z "$OUT_DIR" ]; then
  echo "Usage: $0 -f <fasta_dir> -o <out_dir>" >&2
  exit 1
fi

if [ ! -d "$FASTA_DIR" ]; then
  echo "Error: Folder '$FASTA_DIR' containing fasta files not found." >&2
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
