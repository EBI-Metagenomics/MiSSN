#!/bin/bash

FASTA_FILE=$1
TSV_FILE=$2

if [ $# -lt 2 ];
then
        echo usage: $0 fasta_file cluster_tsv_file
        exit 1
fi

if [ ! -f "$FASTA_FILE" ]; then
        echo "Input FASTA file '$FASTA_FILE' not found."
        exit 1
fi

if [ ! -f "$TSV_FILE" ]; then
        echo "Input cluster TSV file '$TSV_FILE' not found."
        exit 1
fi

OUT_DIR="separated_clusters"
MIN_MEMBERS=5

mkdir -p "$OUT_DIR"

CORES=$(nproc 2>/dev/null || echo 20)

echo "Finding unique clusters with at least $MIN_MEMBERS members..."
awk -v min="$MIN_MEMBERS" '
	{ count[$1]++}
	END { for (cluster in count) if (count[cluster] >= min) print cluster }
' "$TSV_FILE" > "$OUT_DIR/target_clusters.txt"

CLUSTER_COUNT=$(wc -l "$OUT_DIR/target_clusters.txt")
echo "Found $CLUSTER_COUNT clusters..."

awk -v out="$OUT_DIR" '
	NR==FNR { valid[$1]; next } # Only build dictionary
	$1 in valid {
		file = out "/" $1 ".lst"
		print $2 >> file
		close(file)
	}
' "$OUT_DIR/target_clusters.txt" "$TSV_FILE"

echo "Creating FASTA index using seqkit..."
seqkit faidx "$FASTA_FILE"

export FASTA_FILE
export OUT_DIR

cat "$OUT_DIR/target_clusters.txt" | xargs -I {} -P 2 bash -c '
	seqkit faidx "$FASTA_FILE" -l "$OUT_DIR/{}.lst" > "$OUT_DIR/{}.fasta"
'

echo "Clusters successfully separated in $OUT_DIR"
