#!/usr/bin/env python3

import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import defaultdict
import gzip
from pathlib import Path
import time


def check_positive_int(value):
    """Ensures the input is a valid non-negative integer."""
    ivalue = int(value)
    if ivalue < 0:
        raise argparse.ArgumentTypeError(f"Minimum size cannot be negative: {value}")
    return ivalue


def open_text(filepath):
    """Opens a file in text mode, automatically decompressing if it's a .gz file."""
    return (
        gzip.open(filepath, "rt") if filepath.suffix == ".gz" else open(filepath, "r")
    )


def load_cluster_mapping(tsv_file, min_size):
    """Reads the cluster mapping TSV and returns sequence-to-cluster mapping and valid clusters after filtering."""
    print(f"[{time.strftime('%H:%M:%S')}] Loading cluster mapping into memory...")
    seq_to_cluster = {}
    cluster_counts = {}

    # Read TSV and map sequences to clusters
    # For example: MGYP007654976183 MGYP007331563105 (cluster_id, protein_id)
    with open_text(tsv_file) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                cluster_id, protein_id = parts[0], parts[1]
                seq_to_cluster[protein_id] = cluster_id
                cluster_counts[cluster_id] = cluster_counts.get(cluster_id, 0) + 1

    # Filter for clusters with at least min_size members
    valid_clusters = {cid for cid, count in cluster_counts.items() if count >= min_size}
    print(f"\tFound {len(valid_clusters):,} valid clusters with >={min_size} members.")

    return seq_to_cluster, valid_clusters


def group_sequences(fasta_file, seq_to_cluster, valid_clusters):
    """Streams the FASTA file and groups valid sequences by cluster in memory."""
    print(f"[{time.strftime('%H:%M:%S')}] Grouping valid sequences in memory...")

    # Group all sequence strings by their cluster ID
    cluster_data = defaultdict(list)

    # Stream the FASTA file
    with open_text(fasta_file) as f:
        for title, seq in SimpleFastaParser(f):  # Biopython's fast parser
            protein_id = title.split()[0]

            cluster_id = seq_to_cluster.get(protein_id)
            if cluster_id in valid_clusters:
                # Reconstruct the FASTA format and save it to memory
                fasta_string = f">{title}\n{seq}\n"
                cluster_data[cluster_id].append(fasta_string)

    return cluster_data


def write_cluster_files(cluster_data, out_dir):
    """Writes grouped cluster sequences into individual FASTA files."""
    print(
        f"[{time.strftime('%H:%M:%S')}] Writing {len(cluster_data):,} cluster files to disk..."
    )

    for cluster_id, seq_lines in cluster_data.items():
        out_path = out_dir / f"{cluster_id}.fasta"
        with open(out_path, "w") as out_f:
            out_f.write("".join(seq_lines))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-f",
        "--fasta-file",
        type=Path,
        required=True,
        help="Path to the input FASTA file (can be .gz)",
    )
    parser.add_argument(
        "-t",
        "--tsv-file",
        type=Path,
        required=True,
        help="Path to the cluster mapping TSV file",
    )
    parser.add_argument(
        "-o", "--out-dir", type=Path, required=True, help="Output directory"
    )
    parser.add_argument(
        "-s",
        "--min-size",
        type=check_positive_int,
        required=True,
        help="Minimum number of members required to extract a cluster",
    )

    args = parser.parse_args()

    out_dir = args.out_dir / "separated_clusters"
    out_dir.mkdir(parents=True, exist_ok=True)

    seq_to_cluster, valid_clusters = load_cluster_mapping(args.tsv_file, args.min_size)
    cluster_data = group_sequences(args.fasta_file, seq_to_cluster, valid_clusters)
    write_cluster_files(cluster_data, out_dir)

    print(f"[{time.strftime('%H:%M:%S')}] Generated FASTA files in '{out_dir}'.")
