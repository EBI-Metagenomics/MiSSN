#!/bin/bash

FASTA_FILE=$1
TSV_FILE=$2

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

export FASTA_FILE
export OUT_DIR
cat "$OUT_DIR/target_clusters.txt" | xargs -I {} -P "$CORES" bash -c '
	seqtk subseq "$FASTA_FILE" "$OUT_DIR/{}.lst" > "$OUT_DIR/{}.fasta"
'

echo "Clusters separated."
