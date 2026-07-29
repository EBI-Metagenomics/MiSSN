#!/usr/bin/env python3

import argparse
import duckdb
from pathlib import Path
import time


def process_networks(edge_list_dir, metadata_file, min_seq_id, min_coverage, out_dir):
    """Process TSV edge lists using DuckDB to generate node and edge parquet files."""

    if min_seq_id < 0.0 or min_seq_id > 1.0:
        raise ValueError(
            f"Given min_seq_id ({min_seq_id}) is not a valid number between 0.0 and 1.0."
        )
    if min_coverage < 0.0 or min_coverage > 1.0:
        raise ValueError(
            f"Given min_coverage ({min_coverage}) is not a valid number between 0.0 and 1.0."
        )

    out_dir.mkdir(parents=True, exist_ok=True)

    # Initialize DuckDB
    conn = duckdb.connect()

    print(f"[{time.strftime('%H:%M:%S')}] Loading metadata into DuckDB...")
    conn.execute(f"""
    CREATE TABLE global_nodes AS
    SELECT
        mgyp AS id,
        * EXCLUDE(mgyp)
    FROM read_parquet('{metadata_file.as_posix()}')
    """)

    tsv_files = list(edge_list_dir.glob("*.tsv")) + list(edge_list_dir.glob("*.tsv.gz"))
    print(f"\tFound {len(tsv_files)} edge lists to process.")

    for tsv_path in tsv_files:
        base_name = tsv_path.name.replace(".tsv.gz", "").replace(".tsv", "")

        network_dir = out_dir / base_name
        network_dir.mkdir(parents=True, exist_ok=True)

        # Process edges (filter, unify, deduplicate)
        conn.execute(f"""
        CREATE OR REPLACE TEMP TABLE unique_edges AS
        SELECT
            LEAST(qseqid, sseqid) AS source,
            GREATEST(qseqid, sseqid) AS target,
            MAX(pident) AS sequence_identity
        FROM read_csv(
            '{tsv_path.as_posix()}',
            header=false,
            delim='\\t',
            auto_detect=true,
            names=[
                'qseqid',
                'sseqid',
                'pident',
                'length',
                'qlen',
                'slen',
                'positive',
                'evalue',
                'bitscore'
            ]
        )
        WHERE qseqid != sseqid
        AND pident >= {min_seq_id} * 100
        AND positive >= {min_coverage} * GREATEST(qlen, slen)
        GROUP BY source, target;
        """)

        # Extract active nodes
        conn.execute("""
        CREATE OR REPLACE TEMP TABLE active_nodes AS
        SELECT g.*
        FROM global_nodes g
        SEMI JOIN (
            SELECT source AS id FROM unique_edges
            UNION
            SELECT target AS id FROM unique_edges
        ) a ON g.id = a.id
        """)

        # TODO: add color column to the parquet?

        # Export Cosmograph edges and nodes
        edge_file = network_dir / f"{base_name}_edges.parquet"
        node_file = network_dir / f"{base_name}_nodes.parquet"

        conn.execute(f"COPY unique_edges TO '{edge_file.as_posix()}' (FORMAT PARQUET)")

        conn.execute(f"COPY active_nodes TO '{node_file.as_posix()}' (FORMAT PARQUET)")

        # Cleanup
        conn.execute("DROP TABLE unique_edges")
        conn.execute("DROP TABLE active_nodes")

    print(
        f"[{time.strftime('%H:%M:%S')}] Successfully exported all files into '{out_dir}'."
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-e",
        "--edge-list-dir",
        type=Path,
        required=True,
        help="Directory containing edge lists",
    )
    parser.add_argument(
        "-m",
        "--metadata-file",
        type=Path,
        required=True,
        help="Path to the global metadata file",
    )
    parser.add_argument(
        "-i",
        "--min-seq-id",
        type=float,
        required=True,
        help="Minimum sequence identity threshold to keep an edge (0.0 to 1.0)",
    )
    parser.add_argument(
        "-c",
        "--min-coverage",
        type=float,
        required=True,
        help="Minimum alignment coverage threshold to keep an edge (0.0 to 1.0)",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        type=Path,
        required=True,
        help="Output directory",
    )

    args = parser.parse_args()

    process_networks(
        edge_list_dir=args.edge_list_dir,
        metadata_file=args.metadata_file,
        min_seq_id=args.min_seq_id,
        min_coverage=args.min_coverage,
        out_dir=args.out_dir,
    )
