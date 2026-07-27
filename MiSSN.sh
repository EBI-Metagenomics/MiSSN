#!/bin/bash
# MiSSN: A tool for generating MGnify interactive Sequence Similarity Networks

while getopts "f:m:i:c:s:o:" option; do
  case $option in
  f) FASTA_FILE=$OPTARG ;;
  m) METADATA_FILE=$OPTARG ;;
  i) MIN_SEQ_ID=$OPTARG ;;
  c) MIN_COVERAGE=$OPTARG ;;
  s) MIN_CLUSTER_SIZE=$OPTARG ;;
  o) OUT_DIR=$OPTARG ;;
  *)
    echo "Usage: $0 -f <fasta_file> -m <metadata_file> -i <min_seq_id> -c <min_coverage> -s <min_cluster_size> [-o <out_dir>]" >&2
    exit 1
    ;;
  esac
done

if [ -z "$FASTA_FILE" ] || [ -z "$METADATA_FILE" ] || [ -z "$MIN_SEQ_ID" ] || [ -z "$MIN_COVERAGE" ] || [ -z "$MIN_CLUSTER_SIZE" ]; then
  echo "Usage: $0 -f <fasta_file> -m <metadata_file> -i <min_seq_id> -c <min_coverage> -s <min_cluster_size> [-o <out_dir>]" >&2
  exit 1
fi

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

if ! [[ "$MIN_CLUSTER_SIZE" =~ ^[1-9][0-9]*$ ]]; then
  echo "Error: min_cluster_size ($MIN_CLUSTER_SIZE) must be a positive integer." >&2
  exit 1
fi

FILENAME=$(basename "$FASTA_FILE")
BASENAME=$(echo "$FILENAME" | sed -E 's/\.(fasta|fa|faa)(\.gz)?$//')

# Set default output directory if -o was not provided
if [ -z "$OUT_DIR" ]; then
  OUT_DIR="$BASENAME"
fi

# If output directory already exists delete and recreate it
if [ -d "$OUT_DIR" ]; then
  read -p "Output directory '$OUT_DIR' already exists. Do you want to overwrite it? (y/N): " -r
  if [[ $REPLY =~ ^[Yy]$ ]]; then
    rm -rf "$OUT_DIR"
  else
    echo "Aborted by user." >&2
    exit 1
  fi
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

"$SCRIPT_DIR/scripts/linclust.sh" -f "$FASTA_FILE" -i "$MIN_SEQ_ID" -c "$MIN_COVERAGE" -o "$OUT_DIR"

TSV_FILE="$OUT_DIR/linclust_results/${BASENAME}_cluster.tsv"
"$SCRIPT_DIR/scripts/separate_clusters.py" -f "$FASTA_FILE" -t "$TSV_FILE" -o "$OUT_DIR" -s "$MIN_CLUSTER_SIZE"

SEPARATED_CLUSTERS_DIR="$OUT_DIR/separated_clusters"
CLUSTER_ALIGNMENTS_DIR="$OUT_DIR/cluster_alignments"
NETWORKS_DIR="$OUT_DIR/networks"

"$SCRIPT_DIR/scripts/all-vs-all.sh" -f "$SEPARATED_CLUSTERS_DIR" -o "$CLUSTER_ALIGNMENTS_DIR"

"$SCRIPT_DIR/scripts/build_ssn.py" -e "$CLUSTER_ALIGNMENTS_DIR" -m "$METADATA_FILE" -i "$MIN_SEQ_ID" -c "$MIN_COVERAGE" -o "$NETWORKS_DIR"
