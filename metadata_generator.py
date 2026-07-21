import random
import csv
import argparse
import sys
import os
import tarfile
import gzip
import io

BIOMES = [
    "human gut",
    "aquatic",
    "plants",
    "terrestrial",
    "wastewater",
    "food production",
]

OUT_DIR = "dataset"


def generate_mock_csv(fasta_stream, csv_path):
    with open(csv_path, "w", newline="", encoding="utf-8") as output_csv:
        writer = csv.writer(output_csv)
        count = 0

        for line in fasta_stream:
            line = line.strip()

            # We encountered a protein
            if line.startswith(">"):
                protein_id = line[1:].split()[0]
                # Hardcoded max amount of biomes
                num_biomes = 2
                assigned_biomes = random.sample(BIOMES, num_biomes)

                # Create row array
                row = [protein_id] + assigned_biomes
                # Write array as csv row
                writer.writerow(row)
                count += 1

    print(f"Generated mock data for {count} proteins saved to '{csv_path}'.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_fasta",
        help="Input FASTA file",
    )

    args = parser.parse_args()
    input_path = args.input_fasta

    if args.input_fasta == "-" and sys.stdin.isatty():
        parser.print_help()
        sys.exit(1)

    base_name = os.path.basename(input_path)
    # Strip extensions
    for ext in [".tar.gz", ".tgz", ".fasta.gz", ".fa.gz", ".fasta", ".fa", ".gz"]:
        if base_name.endswith(ext):
            base_name = base_name[: -len(ext)]
            break

    csv_path = os.path.join(os.getcwd(), OUT_DIR, f"{base_name}_nodes.csv")

    if input_path.endswith((".tar.gz", ".tgz")):
        with tarfile.open(input_path, "r:gz") as tar:
            fasta_member = next(m for m in tar.getmembers() if m.isfile())

            with tar.extractfile(fasta_member) as f:
                text_stream = io.TextIOWrapper(f, encoding="utf-8")
                generate_mock_csv(text_stream, csv_path)

    elif input_path.endswith(".gz"):
        with gzip.open(input_path, "rt", encoding="utf-8") as f:
            generate_mock_csv(f, csv_path)

    else:
        with open(input_path, "r", encoding="utf-8") as f:
            generate_mock_csv(f, csv_path)
