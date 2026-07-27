# MiSSN: A tool for generating MGnify interactive Sequence Similarity Networks

if [ $# -lt 5 ]; then
  echo "Usage: $0 <fasta_file> <metadata_file> <min_seq_id> <min_coverage> <min_cluster_size> [out_dir]"
  exit 1
fi

FASTA_FILE=$1
METADATA_FILE=$2
MIN_SEQ_ID=$3
MIN_COVERAGE=$4
MIN_CLUSTER_SIZE=$5

if [ ! -f "$FASTA_FILE" ]; then
  echo "Input FASTA file '$FASTA_FILE' not found."
  exit 1
fi

if [ -z "$MIN_SEQ_ID" ] || ! [[ "$MIN_SEQ_ID" =~ ^(0(\.[0-9]+)?|1(\.0+)?)$ ]]; then
  echo "Given min-seq-id ($MIN_SEQ_ID) is not a valid number between 0.0 and 1.0."
  exit 1
fi

if [ -z "$MIN_COVERAGE" ] || ! [[ "$MIN_COVERAGE" =~ ^(0(\.[0-9]+)?|1(\.0+)?)$ ]]; then
  echo "Given min-coverage ($MIN_COVERAGE) is not a valid number between 0.0 and 1.0."
  exit 1
fi

if ! [[ "$MIN_CLUSTER_SIZE" =~ ^[1-9][0-9]*$ ]]; then
  echo "Error: min_cluster_size ($MIN_CLUSTER_SIZE) must be a positive integer." >&2
  exit 1
fi

FILENAME=$(basename "$FASTA_FILE")
BASENAME=$(echo "$FILENAME" | sed -E 's/\.(fasta|fa|faa)(\.gz)?$//')

if [ -n "$6" ]; then
  OUT_DIR="$6"
else
  OUT_DIR="$BASENAME"
fi

# If output directory already exists delete and recreate it
if [ -d "$OUT_DIR" ]; then
  read -p "Output directory '$OUT_DIR' already exists. Do you want to overwrite it? (y/N): " -r
  if [[ $REPLY =~ ^[Yy]$ ]]; then
    rm -rf "$OUT_DIR"
  else
    echo "Aborted by user."
    exit
  fi
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

"$SCRIPT_DIR/scripts/linclust.sh" "$FASTA_FILE" "$MIN_SEQ_ID" "$MIN_COVERAGE" "$OUT_DIR"

TSV_FILE="$OUT_DIR/linclust_results/${BASENAME}_cluster.tsv"
python3 "$SCRIPT_DIR/scripts/separate_clusters.py" "$FASTA_FILE" "$TSV_FILE" "$OUT_DIR" --min-size "$MIN_CLUSTER_SIZE"

"$SCRIPT_DIR/scripts/all-vs-all.sh" "$OUT_DIR/separated_clusters" "$OUT_DIR/cluster_alignments"

python3 "$SCRIPT_DIR/scripts/build_ssn.py" "$OUT_DIR/cluster_alignments" "$METADATA_FILE" "$OUT_DIR/networks" --min-seq-id "$MIN_SEQ_ID" --min-coverage "$MIN_COVERAGE"
