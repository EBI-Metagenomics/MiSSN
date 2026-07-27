#!/usr/bin/env python3

import sys
import os
import time
import gzip
import argparse


def open_text(filepath):
    """Opens a file in text mode, automatically decompressing if it's a .gz file."""
    return (
        gzip.open(filepath, "rt") if filepath.endswith(".gz") else open(filepath, "r")
    )


def open_binary(filepath):
    """Opens a file in binary mode, required for fast line counting."""
    return (
        gzip.open(filepath, "rb") if filepath.endswith(".gz") else open(filepath, "rb")
    )


def count_lines(filepath):
    """Counts lines in a file quickly using binary mode."""
    with open_binary(filepath) as f:
        return sum(1 for _ in f)


def flush_sequence(out_dir, cluster_id, valid_clusters, seq_data):
    """Write a full sequence to disk in one I/O operation."""
    if cluster_id and cluster_id in valid_clusters:
        out_path = os.path.join(out_dir, f"{cluster_id}.fasta")
        # Open in append mode ('a'), write the whole sequence chunk, and close
        with open(out_path, "a") as out_f:
            out_f.write("".join(seq_data))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-f",
        "--fasta-file",
        required=True,
        help="Path to the input FASTA file (can be .gz)",
    )
    parser.add_argument(
        "-t", "--tsv-file", required=True, help="Path to the cluster mapping TSV file"
    )
    parser.add_argument("-o", "--out-dir", required=True, help="Output directory")
    parser.add_argument(
        "-s",
        "--min-size",
        type=int,
        required=True,
        help="Minimum number of members required to extract a cluster",
    )

    args = parser.parse_args()

    out_dir = os.path.join(args.out_dir, "separated_clusters")
    os.makedirs(out_dir, exist_ok=True)

    print(f"[{time.strftime('%H:%M:%S')}] Loading cluster mapping into memory...")
    seq_to_cluster = {}
    cluster_counts = {}

    # Read TSV and map sequences to clusters
    with open_text(args.tsv_file) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                cluster_id, seq_id = parts[0], parts[1]
                seq_to_cluster[seq_id] = cluster_id
                cluster_counts[cluster_id] = cluster_counts.get(cluster_id, 0) + 1

    # Filter for clusters with at least min_size members
    valid_clusters = {
        cid for cid, count in cluster_counts.items() if count >= args.min_size
    }
    print(
        f"\tFound {len(valid_clusters):,} valid clusters with >={args.min_size} members."
    )

    print(f"[{time.strftime('%H:%M:%S')}] Calculating total lines in FASTA file...")
    total_lines = count_lines(args.fasta_file)
    print(f"\tTotal lines to process: {total_lines:,}.")
    print(
        f"[{time.strftime('%H:%M:%S')}] Extracting sequences into individual files in '{out_dir}'..."
    )

    current_cluster = None
    current_seq_data = []
    processed_lines = 0

    # Stream the FASTA file line-by-line
    with open_text(args.fasta_file) as f:
        for line in f:
            processed_lines += 1

            # Print a live progress update every 100,000 lines
            if processed_lines % 100000 == 0:
                percent = (processed_lines / total_lines) * 100
                sys.stderr.write(
                    f"\r\tProgress: {processed_lines:,} / {total_lines:,} lines ({percent:.1f}%) completed..."
                )
                sys.stderr.flush()

            if line.startswith(">"):
                # A new sequence has started. Save the previous one to the hard drive.
                flush_sequence(
                    out_dir, current_cluster, valid_clusters, current_seq_data
                )

                # Extract the sequence ID (ignoring the '>' and any descriptions)
                seq_id = line.strip().split()[0][1:]
                current_cluster = seq_to_cluster.get(seq_id)
                current_seq_data = [line]  # Start storing the new sequence in memory
            else:
                # If it's part of a valid cluster, keep storing the amino acid lines
                if current_cluster:
                    current_seq_data.append(line)

    # Flush the very last sequence in the file
    flush_sequence(out_dir, current_cluster, valid_clusters, current_seq_data)

    sys.stderr.write(
        f"\r\tProgress: {total_lines:,} / {total_lines:,} lines (100.0%) completed.\n"
    )

    print(
        f"[{time.strftime('%H:%M:%S')}] Generated {len(valid_clusters):,} FASTA files in '{out_dir}'."
    )
