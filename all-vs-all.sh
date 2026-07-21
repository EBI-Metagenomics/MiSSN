#!/bin/bash

export FASTA_DIR="separated_clusters_test"
export OUT_DIR="cluster_alignments"

mkdir -p "$OUT_DIR"

run_diamond() {
  FASTA_FILE=$1
  CLUSTER_NAME=$(basename "$FASTA_FILE" .fasta)

  ./diamond makedb \
    --in "$FASTA_FILE" \
    --db "$OUT_DIR/${CLUSTER_NAME}" \
    --quiet

  ./diamond blastp \
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

echo "Calculating DIAMOND alignments..."

for file in "$FASTA_DIR"/*.fasta; do
  run_diamond "$file"
done

echo "Alignments completed."
