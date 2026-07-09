#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=2
#SBATCH --mem=300G
#SBATCH --mail-user=vdmszblcs@yahoo.com

export FASTA_DIR="separated_clusters"
export OUT_DIR="cluster_alignments"

mkdir -p "$OUT_DIR"

run_diamond() {
  FASTA_FILE=$1
  CLUSTER_NAME=$(basename "$FASTA_FILE" .fasta)

  # Force 1 thread per job so we don't oversubscribe the 20 cores
  ./diamond makedb \
    --in "$FASTA_FILE" \
    --db "$OUT_DIR/${CLUSTER_NAME}" \
    --threads 1 \
    --quiet

  # --evalue test, --ultra-sensitive
  ./diamond blastp \
    --query "$FASTA_FILE" \
    --db "$OUT_DIR/${CLUSTER_NAME}.dmnd" \
    --out "$OUT_DIR/${CLUSTER_NAME}.tsv" \
    --outfmt 6 \
    -k 0 \
    --fast \
    -p 1 \
    --quiet
  rm "$OUT_DIR/${CLUSTER_NAME}.dmnd"
}

export -f run_diamond

echo "Calculating DIAMOND alignments..."

find "$FASTA_DIR" -maxdepth 1 -name "*.fasta" | parallel --bar run_diamond {}

echo "Alignments completed."
