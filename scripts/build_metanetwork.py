#!/usr/bin/env python3
"""Concatenate per-cluster SSNs into node and edge metanetwork Parquets."""

import argparse
from pathlib import Path

import duckdb


def build_metanetwork(out_dir):
    """Merge out_dir/networks/*/*_{nodes,edges}.parquet into out_dir."""
    out_dir = Path(out_dir)
    networks = out_dir / "networks"
    if not networks.is_dir():
        raise FileNotFoundError(f"Networks directory not found: {networks}")

    input_files = {}
    for network_type in ("nodes", "edges"):
        input_files[network_type] = sorted(str(path) for path in networks.glob(f"*/*_{network_type}.parquet"))
        if not input_files[network_type]:
            raise FileNotFoundError(f"No cluster {network_type} Parquet files found in {networks}")

    with duckdb.connect() as conn:
        for network_type, paths in input_files.items():
            output = out_dir / f"metanetwork_{network_type}.parquet"
            conn.execute(
                "COPY (SELECT * FROM read_parquet($paths, hive_partitioning=false)) "
                "TO $output (FORMAT PARQUET)",
                {"paths": paths, "output": str(output)},
            )
            print(f"Merged {len(paths)} {network_type} files into '{output}'.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-o", "--out-dir", type=Path, required=True,
        help="Pipeline output directory containing the networks subdirectory",
    )
    args = parser.parse_args()
    try:
        build_metanetwork(args.out_dir)
    except (OSError, duckdb.Error) as error:
        parser.exit(1, f"Error: {error}\n")
