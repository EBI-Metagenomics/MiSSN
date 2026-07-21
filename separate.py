import sys
import os
import time
import re
import gzip

if len(sys.argv) < 3:
    print("Usage: python3 extract.py <fasta_file> <cluster_tsv>")
    sys.exit(1)

fasta_file = sys.argv[1]
tsv_file = sys.argv[2]
min_members = 5

base_name = os.path.basename(fasta_file)
clean_name = re.sub(r"\.(fasta|fa|faa)(\.gz)?$", "", base_name)
out_dir = os.path.join("separated_clusters", clean_name)

# Create output directory
os.makedirs(out_dir, exist_ok=True)


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


print("Loading cluster mapping into memory...")
seq_to_cluster = {}
cluster_counts = {}

# 1. Read TSV and map sequences to clusters
with open_text(tsv_file) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            cluster_id, seq_id = parts[0], parts[1]
            seq_to_cluster[seq_id] = cluster_id
            cluster_counts[cluster_id] = cluster_counts.get(cluster_id, 0) + 1

# 2. Filter for clusters with at least 5 members
valid_clusters = {cid for cid, count in cluster_counts.items() if count >= min_members}
print(f"Found {len(valid_clusters):,} valid clusters with >={min_members} members.")

print("Calculating total lines in FASTA file for progress tracking...")


# Fast line counting
def count_lines(filepath):
    with open_binary(filepath) as f:
        return sum(1 for _ in f)


total_lines = count_lines(fasta_file)
print(f"Total lines to process: {total_lines:,}")

print(f"Extracting sequences into individual files in ./{out_dir}/...")

current_cluster = None
current_seq_data = []


# Helper function to write a full sequence to disk in one I/O operation
def flush_sequence():
    if current_cluster and current_cluster in valid_clusters:
        out_path = os.path.join(out_dir, f"{current_cluster}.fasta")
        # Open in append mode ('a'), write the whole sequence chunk, and close
        with open(out_path, "a") as out_f:
            out_f.write("".join(current_seq_data))


processed_lines = 0
start_time = time.time()

# 3. Stream the FASTA file line-by-line
with open_text(fasta_file) as f:
    for line in f:
        processed_lines += 1

        # Print a live progress update every 100,000 lines
        if processed_lines % 100000 == 0:
            percent = (processed_lines / total_lines) * 100
            sys.stderr.write(
                f"\rProgress: {processed_lines:,} / {total_lines:,} lines ({percent:.1f}%) completed..."
            )
            sys.stderr.flush()

        if line.startswith(">"):
            # A new sequence has started. Save the previous one to the hard drive.
            flush_sequence()

            # Extract the sequence ID (ignoring the '>' and any descriptions)
            seq_id = line.strip().split()[0][1:]
            current_cluster = seq_to_cluster.get(seq_id)
            current_seq_data = [line]  # Start storing the new sequence in memory
        else:
            # If it's part of a valid cluster, keep storing the amino acid lines
            if current_cluster:
                current_seq_data.append(line)

# Flush the very last sequence in the file
flush_sequence()

elapsed = time.time() - start_time
sys.stderr.write(
    f"\rProgress: {total_lines:,} / {total_lines:,} lines (100.0%) completed...\n"
)
print(f"Done! Generated {len(valid_clusters):,} files in {elapsed:.2f} seconds.")
