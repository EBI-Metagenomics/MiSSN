#!/bin/bash
# MiSSN: A workflow for generating MGnify interactive Sequence Similarity Networks

usage() {
  echo "Usage: $0 -f <fasta_file> -b <biome_metadata_file> -p <pfam_metadata_file> -i <min_seq_id> -c <min_coverage> -s <min_cluster_size> [-o <out_dir>] [-g <color_group_parts>]" >&2
  echo "Try '$0 -h' for more information." >&2
  exit 1
}

print_help() {
  cat <<EOF
MiSSN: A workflow for generating MGnify interactive Sequence Similarity Networks

Usage: $0 -f <fasta_file> -b <biome_metadata_file> -p <pfam_metadata_file> -i <min_seq_id> -c <min_coverage> -s <min_cluster_size> [-o <out_dir>] [-g <color_group_parts>]

Required arguments:
  -f    Path to the input FASTA file (can be .gz).
  -b    Path to the GOLD biome classifications Parquet file.
  -p    Path to the Pfam accessions Parquet file.
  -i    Minimum sequence identity threshold to keep an edge (0.0 to 1.0).
  -c    Minimum alignment coverage threshold to keep an edge (0.0 to 1.0).
  -s    Minimum number of members required to extract a cluster.

Optional arguments:
  -o    Output directory. (Default: dynamically named based on the FASTA file).
  -g    Number of leading lineage parts for color grouping. (Default: 4).
  -h    Show this help message and exit.
EOF

  exit 0
}

while getopts "f:b:p:i:c:s:o:g:h" option; do
  case $option in
  f) FASTA_FILE=$OPTARG ;;
  b) BIOME_METADATA_FILE=$OPTARG ;;
  p) PFAM_METADATA_FILE=$OPTARG ;;
  i) MIN_SEQ_ID=$OPTARG ;;
  c) MIN_COVERAGE=$OPTARG ;;
  s) MIN_CLUSTER_SIZE=$OPTARG ;;
  o) OUT_DIR=$OPTARG ;;
  g) COLOR_GROUP_PARTS=$OPTARG ;;
  h) print_help ;;
  *) usage ;;
  esac
done

if [ -z "$FASTA_FILE" ] || [ -z "$BIOME_METADATA_FILE" ] || [ -z "$PFAM_METADATA_FILE" ] || [ -z "$MIN_SEQ_ID" ] || [ -z "$MIN_COVERAGE" ] || [ -z "$MIN_CLUSTER_SIZE" ]; then
  usage
fi

if [ ! -f "$FASTA_FILE" ]; then
  echo "Error: Input FASTA file '$FASTA_FILE' not found." >&2
  exit 1
fi

if [ ! -f "$BIOME_METADATA_FILE" ]; then
  echo "Error: Biome metadata file '$BIOME_METADATA_FILE' not found." >&2
  exit 1
fi

if [ ! -f "$PFAM_METADATA_FILE" ]; then
  echo "Error: Pfam metadata file '$PFAM_METADATA_FILE' not found." >&2
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

if [ -z "$COLOR_GROUP_PARTS" ]; then
  COLOR_GROUP_PARTS=4
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

"$SCRIPT_DIR/scripts/linclust.sh" \
  -f "$FASTA_FILE" \
  -i "$MIN_SEQ_ID" \
  -c "$MIN_COVERAGE" \
  -o "$OUT_DIR"

TSV_FILE="$OUT_DIR/linclust_results/${BASENAME}_cluster.tsv"
"$SCRIPT_DIR/scripts/separate_clusters.py" \
  -f "$FASTA_FILE" \
  -t "$TSV_FILE" \
  -o "$OUT_DIR" \
  -s "$MIN_CLUSTER_SIZE"

SEPARATED_CLUSTERS_DIR="$OUT_DIR/separated_clusters"
CLUSTER_ALIGNMENTS_DIR="$OUT_DIR/cluster_alignments"
NETWORKS_DIR="$OUT_DIR/networks"

"$SCRIPT_DIR/scripts/all-vs-all.sh" \
  -f "$SEPARATED_CLUSTERS_DIR" \
  -o "$CLUSTER_ALIGNMENTS_DIR"

"$SCRIPT_DIR/scripts/build_ssn.py" \
  -e "$CLUSTER_ALIGNMENTS_DIR" \
  -b "$BIOME_METADATA_FILE" \
  -p "$PFAM_METADATA_FILE" \
  -i "$MIN_SEQ_ID" \
  -c "$MIN_COVERAGE" \
  -o "$NETWORKS_DIR" \
  -g "$COLOR_GROUP_PARTS"
