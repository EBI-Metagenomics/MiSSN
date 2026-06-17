#!/bin/bash

COVERAGE=0.8 # -c default 0.8
COVERAGE_MODE=1 # --cov-mode, 1: coverage of target 

FASTA_FILE=$1
MIN_SEQ_ID=$2

if [ $# -lt 2 ];
then
        echo usage: $0 fasta_file min_seq_id
	exit 1
fi

if [ ! -f "$FASTA_FILE" ]; then
	echo "Input FASTA file '$FASTA_FILE' not found."
	exit 1
fi

if [ -z "$MIN_SEQ_ID" ] || ! [[ "$MIN_SEQ_ID" =~ ^(0(\.[0-9]+)?|1(\.0+)?)$ ]]; then
	echo "Given min-seq-id ($MIN_SEQ_ID) is not a valid number between 0.0 and 1.0."
	exit 1
fi

FILENAME=$(basename "$FASTA_FILE")
BASENAME="${FILENAME%.*}"
SAFE_SEQ_ID=$(LC_NUMERIC=C printf "%.2f" "$MIN_SEQ_ID" | sed 's/\.//')

OUT_DIR="linclust_results"
mkdir -p "$OUT_DIR"
TMP_DIR="tmp"
mkdir -p "$TMP_DIR"
DB_DIR="dataset"
mkdir -p "$DB_DIR"

SEQ_DB="$DB_DIR/${BASENAME}_seqDB"
CLU_DB="$DB_DIR/${BASENAME}_${SAFE_SEQ_ID}_cluDB"
REP_DB="$DB_DIR/${BASENAME}_${SAFE_SEQ_ID}_repDB"

TSV_FILE="$OUT_DIR/${BASENAME}_${SAFE_SEQ_ID}_cluster.tsv"
REP_FASTA="$OUT_DIR/${BASENAME}_${SAFE_SEQ_ID}_rep_seq.fasta"
ALL_FASTA="$OUT_DIR/${BASENAME}_${SAFE_SEQ_ID}_all_seqs.fasta"

echo "Running MMSeqs Linclust with $MIN_SEQ_ID minimum sequence identity, $COVERAGE coverage and coverage mode $COVERAGE_MODE..."

# 1. Create reusable sequence database
if [[ -f "$SEQ_DB" ]]; then
	echo "Sequence DB already exists... skipping step"
else
	mmseqs createdb "$FASTA_FILE" "$SEQ_DB"
fi

# 2. Run Linclust
if [[ -f "$CLU_DB" ]]; then
	echo "Clustering DB for $MIN_SEQ_ID already exists... skipping step"
else 
	mmseqs linclust "$SEQ_DB" "$CLU_DB" "$TMP_DIR" \
	-c "$COVERAGE" \
	--cov-mode "$COVERAGE_MODE" \
	--min-seq-id "$MIN_SEQ_ID"
fi

# 3. Generate the cluster TSV
if [[ -f "$TSV_FILE" ]]; then
	echo "Cluster TSV already exists... skipping step"
else
	mmseqs createtsv "$SEQ_DB" "$SEQ_DB" "$CLU_DB" "$TSV_FILE"
fi

# 4. Extract representative FASTA
if [[ -f "$REP_FASTA" ]]; then
	echo "Representative FASTA already exists... skipping step"
else
	mmseqs result2repseq "$SEQ_DB" "$CLU_DB" "$REP_DB"
	mmseqs result2flat "$SEQ_DB" "$SEQ_DB" "$REP_DB" "$REP_FASTA"
fi

# 5. Extract clustered FASTA
if [[ -f "$ALL_FASTA" ]]; then
	echo "Clustered all-seqs FASTA already exists... skipping step"
else
	mmseqs result2flat "$SEQ_DB" "$SEQ_DB" "$CLU_DB" "$ALL_FASTA"
fi

echo "Results successfully generated in $OUT_DIR" 
