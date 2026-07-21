import argparse
import duckdb
import os
import glob
import time

edge_list_dir = "cluster_alignments"
output_dir = "networks"

parser = argparse.ArgumentParser()
parser.add_argument(
    "--metadata-file", type=str, required=True, help="Path to the global metadata file"
)
parser.add_argument(
    "--min_seq_id",
    type=float,
    required=True,
    help="Minimum sequence identity threshold to keep an edge",
)
parser.add_argument(
    "--min_coverage",
    type=float,
    required=True,
    help="Minimum alignment coverage threshold to keep an edge",
)

args = parser.parse_args()

metadata_file = args.metadata_file
min_seq_id = args.min_seq_id
min_coverage = args.min_coverage

os.makedirs(output_dir, exist_ok=True)

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
print(f"[{time.strftime('%H:%M:%S')}] Found {len(tsv_files)} edge lists to process.\n")

for tsv_path in tsv_files:
    base_name = os.path.basename(tsv_path).replace(".tsv.gz", "").replace(".tsv", "")

    network_dir = os.path.join(output_dir, base_name)
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
    AND pident >= {min_seq_id}
    AND positive >= ({min_coverage} / 100.0) * GREATEST(qlen, slen)
    GROUP BY source, target;
    """)

    result = conn.execute("SELECT COUNT(*) FROM unique_edges").fetchone()
    edge_count = result[0] if result else 0
    print(f"[{time.strftime('%H:%M:%S')}] Extracted {edge_count:,} unique edges.")

    # Extract active nodes
    conn.execute("""
    CREATE OR REPLACE TMP TABLE active_nodes AS
    SELECT DISTINCT source AS id FROM unique_edges
    UNION
    SELECT DISTINCT target AS id FROM unique_edges
    """)

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
    f"[{time.strftime('%H:%M:%S')}] Successfully exported all files into '{output_dir}'"
)
