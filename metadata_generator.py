import random
import csv
import argparse
import sys
import os

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
    with open(csv_path, "w", newline="") as output_csv:
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
        nargs="?",
        type=argparse.FileType("r"),
        default=sys.stdin,
        help="Input FASTA file",
    )

    args = parser.parse_args()

    if args.input_fasta.name == "<stdin>" and sys.stdin.isatty():
        parser.print_help()
        sys.exit(1)

    base_name = args.input_fasta.name.replace(".fasta", "")
    csv_path = os.path.join(os.getcwd(), f"{base_name}_nodes.csv")
    print(csv_path)
    generate_mock_csv(args.input_fasta, csv_path)
