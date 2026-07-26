import argparse
import duckdb
import os
import glob
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

    os.makedirs(out_dir, exist_ok=True)

    # Initialize DuckDB
    conn = duckdb.connect()

    print(f"[{time.strftime('%H:%M:%S')}] Loading metadata into DuckDB...")
    conn.execute(f"""
    CREATE TABLE global_nodes AS
    SELECT
        mgyp AS id,
        * EXCLUDE(mgyp)
    FROM read_parquet('{metadata_file}')
    """)

    tsv_files = glob.glob(os.path.join(edge_list_dir, "*.tsv")) + glob.glob(
        os.path.join(edge_list_dir, "*.tsv.gz")
    )
    print(f"\tFound {len(tsv_files)} edge lists to process.")

    for tsv_path in tsv_files:
        base_name = (
            os.path.basename(tsv_path).replace(".tsv.gz", "").replace(".tsv", "")
        )

        network_dir = os.path.join(out_dir, base_name)
        os.makedirs(network_dir, exist_ok=True)

        # Process edges (filter, unify, deduplicate)
        conn.execute(f"""
        CREATE OR REPLACE TEMP TABLE unique_edges AS
        SELECT
            LEAST(qseqid, sseqid) AS source,
            GREATEST(qseqid, sseqid) AS target,
            MAX(pident) AS sequence_identity
        FROM read_csv(
            '{tsv_path}',
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
        SELECT DISTINCT source AS id FROM unique_edges
        UNION
        SELECT DISTINCT target AS id FROM unique_edges
        """)

        # TODO: add color column to the parquet?

        # Export Cosmograph edges and nodes
        conn.execute(
            f"COPY unique_edges TO '{os.path.join(network_dir, f'{base_name}_edges.parquet')}' (FORMAT PARQUET)"
        )

        node_query = (
            "SELECT g.* FROM global_nodes g SEMI JOIN active_nodes a ON g.id = a.id"
        )
        conn.execute(
            f"COPY ({node_query}) TO '{os.path.join(network_dir, f'{base_name}_nodes.parquet')}' (FORMAT PARQUET)"
        )

        # Cleanup
        conn.execute("DROP TABLE unique_edges")
        conn.execute("DROP TABLE active_nodes")

    print(
        f"[{time.strftime('%H:%M:%S')}] Successfully exported all files into '{out_dir}'."
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "edge_list_dir",
        type=str,
        help="Directory containing edge lists",
    )
    parser.add_argument(
        "metadata_file",
        type=str,
        help="Path to the global metadata file",
    )
    parser.add_argument(
        "out_dir",
        type=str,
        help="Output directory",
    )
    parser.add_argument(
        "--min-seq-id",
        type=float,
        required=True,
        help="Minimum sequence identity threshold to keep an edge (0.0 to 1.0)",
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        required=True,
        help="Minimum alignment coverage threshold to keep an edge (0.0 to 1.0)",
    )

    args = parser.parse_args()

    process_networks(
        edge_list_dir=args.edge_list_dir,
        metadata_file=args.metadata_file,
        min_seq_id=args.min_seq_id,
        min_coverage=args.min_coverage,
        out_dir=args.out_dir,
    )
