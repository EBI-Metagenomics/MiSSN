#!/bin/sh

COVERAGE=0.8 # -c default 0.8
COVERAGE_MODE=1 # --cov-mode, 1: coverage of target 

if [ $# -lt 2 ];
then
        echo usage: $0 fasta_file min_seq_id
	exit 1
fi

# TODO: param checks here

FASTA_FILE=$1
MIN_SEQ_ID=$2

OUTPUT_DIR="linclust_results"
mkdir -p "$OUTPUT_DIR"

echo "Running MMSeqs Linclust with $MIN_SEQ_ID minimum sequence identity..."

SAFE_SEQ_ID=$(echo "$MIN_SEQ_ID" | sed 's/\.//')

mmseqs easy-linclust "$FASTA_FILE" "$OUTPUT_DIR/${SAFE_SEQ_ID}_seq_id" tmp \
	-c "$COVERAGE" \
	--cov-mode "$COVERAGE_MODE" \
	--min-seq-id "$MIN_SEQ_ID"

echo "Results successfully generated in $OUTPUT_DIR" 
