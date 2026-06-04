#!/bin/sh

COVERAGE=0.8 # -c default 0.8
COVERAGE_MODE=1 # --cov-mode, 1: coverage of target 
DB_DIR="dataset"
TMP_DIR="tmp"

if [ $# -lt 2 ];
then
        echo usage: $0 fasta_file min_seq_id
	exit 1
fi

# TODO: param checks here

FASTA_FILE=$1
FILENAME=$(basename "$FASTA_FILE")
BASENAME="${FILENAME%.*}"

MIN_SEQ_ID=$2

OUT_DIR="linclust_results"
mkdir -p "$OUT_DIR"

echo "Running MMSeqs Linclust with $MIN_SEQ_ID minimum sequence identity, $COVERAGE coverage and coverage mode $COVERAGE_MODE..."

SAFE_SEQ_ID=$(LC_NUMERIC=C printf "%.2f" "$MIN_SEQ_ID" | sed 's/\.//')

# 1. Create reusable sequence database
mmseqs createdb "$FASTA_FILE" "$DB_DIR/${BASENAME}_seqDB"

# 2. Run Linclust
mkdir -p "$TMP_DIR"
mmseqs linclust "$DB_DIR/${BASENAME}_seqDB" "$DB_DIR/${BASENAME}_cluDB" "$TMP_DIR" \
	-c "$COVERAGE" \
	--cov-mode "$COVERAGE_MODE" \
	--min-seq-id "$MIN_SEQ_ID"

# 3. Generate the cluster TSV
mmseqs createtsv "$DB_DIR/${BASENAME}_seqDB" "$DB_DIR/${BASENAME}_seqDB" \
	"$DB_DIR/${BASENAME}_cluDB" "$OUT_DIR/${BASENAME}_${SAFE_SEQ_ID}_cluster.tsv"

# 4. Optional?: Extract representative FASTA
mmseqs result2repseq "$DB_DIR/${BASENAME}_seqDB" "$DB_DIR/${BASENAME}_cluDB" \
	"$DB_DIR/${BASENAME}_repDB"
mmseqs result2flat "$DB_DIR/${BASENAME}_seqDB" "$DB_DIR/${BASENAME}_seqDB" \
	"$DB_DIR/${BASENAME}_repDB" "$OUT_DIR/${BASENAME}_${SAFE_SEQ_ID}_rep_seq.fasta"

# Extract clustered FASTA
mmseqs result2flat "$DB_DIR/${BASENAME}_seqDB" "$DB_DIR/${BASENAME}_seqDB" \
	"$DB_DIR/${BASENAME}_cluDB" "$OUT_DIR/${BASENAME}_${SAFE_SEQ_ID}_all_seqs.fasta"

echo "Results successfully generated in $OUT_DIR" 
