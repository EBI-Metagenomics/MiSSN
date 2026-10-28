#!/usr/bin/env python3

import argparse
import colorsys
import re
import time
from pathlib import Path

import duckdb
import pandas as pd


def clean_lineage(raw):
    raw = raw.strip()
    idx = raw.find("root:")
    if idx > 0:
        return raw[idx:]
    if raw == "root" or raw.startswith("root:"):
        return raw
    return raw


class TreeNode:
    __slots__ = ("children", "depth", "hue", "hue_span", "is_lineage", "name")

    def __init__(self, name, depth):
        self.name = name
        self.children = {}
        self.hue = 0.0
        self.hue_span = 0.0
        self.depth = depth
        self.is_lineage = False  # True if some input lineage ends exactly here


def build_tree(lineages):
    root = TreeNode("root", 0)
    for lin in lineages:
        parts = [p.strip() for p in lin.split(":") if p.strip() != ""]
        node = root
        for depth, part in enumerate(parts[1:], start=1):  # skip literal "root"
            node = node.children.setdefault(part, TreeNode(part, depth))
        node.is_lineage = True
    return root


def collect_group_nodes(node, path_parts, max_depth, out):
    """Collect every (path_tuple, node) with 1 <= depth <= max_depth."""
    if 1 <= node.depth <= max_depth:
        out.append((tuple(path_parts), node))
    if node.depth < max_depth:
        for name, child in node.children.items():
            collect_group_nodes(child, path_parts + [name], max_depth, out)


def assign_flat_group_hues(group_nodes):
    """
    Give every collected group node (spanning all depths from 1 up to the
    group boundary, flattened into one list) its own well-separated hue
    from a single palette that spans the whole wheel. Sorting by full
    path keeps a parent immediately next to its own children in the
    sequence, so they land on neighboring hues, while unrelated branches
    end up far apart. Returns {path_tuple: (band_start, band_end)} for
    seeding the narrower nudge below the group boundary.
    """
    group_nodes.sort(key=lambda pn: pn[0])
    n = len(group_nodes)
    bands = {}
    if n == 0:
        return bands
    span = 1.0 / n
    gap = span * 0.04 if n > 1 else 0.0
    for i, (path, node) in enumerate(group_nodes):
        start = i * span
        end = start + span - gap
        node.hue = (start + end) / 2.0 % 1.0
        node.hue_span = end - start
        bands[path] = (start, end)
    return bands


def assign_nested_hues(node, hue_start, hue_end):
    """
    Below the group boundary: subdivide the parent's own band narrowly so
    children stay close to their group's hue (same color family), just
    nudged apart enough to tell siblings apart. Lightness (applied
    separately, by absolute depth) carries most of the "going deeper"
    variability from here on.
    """
    node.hue = (hue_start + hue_end) / 2.0 % 1.0
    node.hue_span = hue_end - hue_start

    children = sorted(node.children.values(), key=lambda n: n.name)
    n = len(children)
    if n == 0:
        return

    parent_span = hue_end - hue_start
    margin = parent_span * 0.35
    usable_start = hue_start + margin / 2
    usable_span = max(parent_span - margin, 1e-6)
    child_span = usable_span / n
    for i, child in enumerate(children):
        s = usable_start + i * child_span
        e = s + child_span
        assign_nested_hues(child, s, e)


def assign_hues(root, group_parts):
    """
    group_parts: how many lineage parts (including "root") define a color
    group, e.g. 3 -> every "root:X:Y" prefix in the data (X = 2-part,
    Y = 3-part) is a group, and they all share one flat, well-separated
    palette. Depths beyond the group boundary nudge off their nearest
    group's hue instead of getting a brand new one.
    """
    group_max_depth = group_parts - 1  # tree depth of the deepest group level
    group_nodes = []
    collect_group_nodes(root, [], group_max_depth, group_nodes)
    bands = assign_flat_group_hues(group_nodes)

    def recurse(node, path_parts):
        if node.depth == group_max_depth:
            start, end = bands.get(
                tuple(path_parts), (node.hue - 0.005, node.hue + 0.005)
            )
            for child in node.children.values():
                assign_nested_hues(child, start, end)
        else:
            for name, child in node.children.items():
                recurse(child, path_parts + [name])

    recurse(root, [])


# Lightness/saturation ramp per depth level (index 0 = the top-level
# category itself, i.e. depth 1 in the tree). Depth beyond the table
# length is clamped to the last entry.
LIGHTNESS_BY_LEVEL = [0.38, 0.48, 0.58, 0.68, 0.76, 0.83, 0.88]
SATURATION_BY_LEVEL = [0.75, 0.72, 0.68, 0.64, 0.60, 0.56, 0.52]


def depth_style(depth):
    """depth=0 is 'root' itself, depth=1 is the main/top-level category."""
    idx = max(depth - 1, 0)
    idx = min(idx, len(LIGHTNESS_BY_LEVEL) - 1)
    return SATURATION_BY_LEVEL[idx], LIGHTNESS_BY_LEVEL[idx]


def hsl_to_rgb255(hue, saturation, lightness):
    # colorsys uses HLS ordering
    r, g, b = colorsys.hls_to_rgb(hue, lightness, saturation)
    return tuple(round(c * 255) for c in (r, g, b))


def collect_colors(node, prefix_parts, out):
    if node.depth > 0:
        saturation, lightness = depth_style(node.depth)
        rgb = hsl_to_rgb255(node.hue, saturation, lightness)
        lineage = "root:" + ":".join(prefix_parts) if prefix_parts else "root"
        out[lineage] = rgb
    for name, child in sorted(node.children.items()):
        collect_colors(child, prefix_parts + [name], out)


def build_color_table(lineages, group_parts=3):
    if group_parts < 2:
        raise ValueError("group_parts must be >= 2")
    root = build_tree(lineages)
    assign_hues(root, group_parts)
    colors = {}
    collect_colors(root, [], colors)
    # also color the bare "root" lineage if it was present in the input
    if "root" in lineages:
        colors.setdefault("root", (200, 200, 200))
    ordered = {lin: colors.get(lin, (200, 200, 200)) for lin in lineages}
    return ordered


def normalize_color(value):
    """Validate a hex or RGB colour and return the canonical RGB string."""
    value = value.strip()
    if re.fullmatch(r"#[0-9a-fA-F]{6}", value):
        channels = tuple(int(value[i:i + 2], 16) for i in (1, 3, 5))
    else:
        match = re.fullmatch(r"rgb\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)", value)
        if not match:
            raise ValueError(f"Invalid color {value!r}: expected #RRGGBB or rgb(r,g,b)")
        channels = tuple(int(channel) for channel in match.groups())
        if any(channel > 255 for channel in channels):
            raise ValueError(f"Invalid color {value!r}: RGB channels must be 0–255")
    return "rgb({},{},{})".format(*channels)


def setup_accession_colors(conn, color_metadata_file):
    """Load and validate overrides once; enforce one colour per accession."""
    conn.execute("CREATE TEMP TABLE accession_colors (id VARCHAR, color VARCHAR)")
    if color_metadata_file is None:
        return
    rows = conn.execute(
        "SELECT * FROM read_parquet(?) LIMIT 0", [str(color_metadata_file)]
    )
    columns = {column[0]: column[1] for column in rows.description}
    for name in ("accession", "color"):
        if name not in columns or str(columns[name]) != "VARCHAR":
            raise ValueError(f"Color metadata requires a string column named {name!r}")
    # Fetch only the two required columns, even if the input has extra fields.
    records = conn.execute(
        "SELECT accession, color FROM read_parquet(?)", [str(color_metadata_file)]
    ).fetchall()
    colors = {}
    for accession, color in records:
        if accession is None or not accession.strip() or color is None or not color.strip():
            raise ValueError("Color metadata accession and color must not be null or blank")
        normalized = normalize_color(color)
        if accession in colors and colors[accession] != normalized:
            raise ValueError(f"Conflicting colors for accession {accession!r}")
        colors[accession] = normalized
    if colors:
        conn.executemany("INSERT INTO accession_colors VALUES (?, ?)", list(colors.items()))


def setup_database(biome_metadata_file=None, pfam_metadata_file=None, color_metadata_file=None):
    """Aggregate optional annotations and load accession colour overrides."""
    conn = duckdb.connect()
    try:
        for table, column, path in (
            ("temp_biome", "biome", biome_metadata_file),
            ("temp_pfam", "pfam_accession", pfam_metadata_file),
        ):
            conn.execute(f"CREATE TEMP TABLE {table} (id VARCHAR, {column} VARCHAR)")
            if path is not None:
                print(f"[{time.strftime('%H:%M:%S')}] Loading {column} metadata into DuckDB...")
                conn.execute(f"""
                    INSERT INTO {table}
                    SELECT mgyp, string_agg(DISTINCT {column}, ';' ORDER BY {column})
                    FROM read_parquet(?) GROUP BY mgyp
                """, [str(path)])

        conn.execute("""
            CREATE TABLE global_nodes AS
            SELECT
                COALESCE(b.id, p.id) AS id,
                b.biome,
                p.pfam_accession
            FROM temp_biome b
            FULL OUTER JOIN temp_pfam p ON b.id = p.id
        """)
        setup_accession_colors(conn, color_metadata_file)
    except Exception:
        conn.close()
        raise

    return conn


def setup_biome_colors(conn, group_parts):
    """Computes biome colors and registers the mapping in DuckDB."""
    print(f"[{time.strftime('%H:%M:%S')}] Computing biome colors...")
    lineage_rows = conn.execute(
        "SELECT DISTINCT biome FROM global_nodes WHERE biome IS NOT NULL"
    ).fetchall()

    formatted_colors = []
    single_raw_lineages = []
    cleaned_lineages = []

    for row in lineage_rows:
        raw_lin = row[0]
        if not raw_lin:
            continue

        if ";" in raw_lin:
            # It is a combined biome: color it black
            formatted_colors.append((raw_lin, "rgb(0,0,0)"))
        else:
            # It is a single biome: save it for coloring
            single_raw_lineages.append(raw_lin)
            cleaned_lineages.append(clean_lineage(raw_lin))

    color_map = build_color_table(cleaned_lineages, group_parts=group_parts)
    for raw_lin, cleaned_lin in zip(single_raw_lineages, cleaned_lineages):
        r, g, b = color_map.get(cleaned_lin, (200, 200, 200))
        formatted_colors.append((raw_lin, f"rgb({r},{g},{b})"))

    color_df = pd.DataFrame(formatted_colors, columns=pd.Index(["biome", "color"]))
    color_df = color_df.astype("string")

    conn.register("color_lookup", color_df)
    conn.execute(
        "CREATE TEMP TABLE biome_color_map AS SELECT biome, color FROM color_lookup"
    )


def process_single_network(conn, tsv_path, out_dir, min_seq_id, min_coverage):
    """
    Processes an edge list, requiring min_seq_id and min_coverage on both
    sequences per HSP, and extracts nodes and edges to Parquet files.
    DIAMOND's qcovhsp and scovhsp fields are percentages (0 to 100).
    """
    base_name = tsv_path.name.replace(".tsv.gz", "").replace(".tsv", "")

    # Load raw edges to capture all sequences
    try:
        conn.execute(f"""
        CREATE OR REPLACE TEMP TABLE raw_edges AS
        SELECT *
        FROM read_csv(
            '{tsv_path.as_posix()}',
            header=false,
            delim='\\t',
            auto_detect=true,
            types={{'qseqid': 'VARCHAR', 'sseqid': 'VARCHAR'}},
            names=[
                'qseqid',
                'sseqid',
                'pident',
                'length',
                'qlen',
                'slen',
                'positive',
                'evalue',
                'bitscore',
                'qcovhsp',
                'scovhsp'
            ]
        )
        """)
    except duckdb.InvalidInputException:
        print(f"\tSkipping empty or invalid file: {tsv_path.name}")
        return

    network_dir = out_dir / base_name
    network_dir.mkdir(parents=True, exist_ok=True)

    conn.execute(f"""
    CREATE OR REPLACE TEMP TABLE unique_edges AS
    SELECT
        LEAST(qseqid, sseqid) AS source,
        GREATEST(qseqid, sseqid) AS target,
        MAX(pident) AS sequence_identity
    FROM raw_edges
    WHERE qseqid != sseqid
    AND pident >= {min_seq_id} * 100
    AND LEAST(qcovhsp, scovhsp) >= {min_coverage} * 100
    GROUP BY source, target;
    """)

    # Extract active nodes
    conn.execute("""
    CREATE OR REPLACE TEMP TABLE active_nodes AS
    SELECT
        a.id,
        g.biome,
        g.pfam_accession,
        COALESCE(o.color, c.color, 'rgb(200, 200, 200)') AS color
    FROM (
        SELECT CAST(qseqid AS VARCHAR) AS id FROM raw_edges
        UNION
        SELECT CAST(sseqid AS VARCHAR) AS id FROM raw_edges
    ) a
    LEFT JOIN global_nodes g ON g.id = a.id
    LEFT JOIN biome_color_map c ON g.biome = c.biome
    LEFT JOIN accession_colors o ON a.id = o.id
    """)

    # Export Cosmograph edges and nodes
    edge_file = network_dir / f"{base_name}_edges.parquet"
    node_file = network_dir / f"{base_name}_nodes.parquet"

    conn.execute(f"COPY unique_edges TO '{edge_file.as_posix()}' (FORMAT PARQUET)")

    conn.execute(f"COPY active_nodes TO '{node_file.as_posix()}' (FORMAT PARQUET)")

    # Cleanup
    conn.execute("DROP TABLE raw_edges")
    conn.execute("DROP TABLE unique_edges")
    conn.execute("DROP TABLE active_nodes")


def process_networks(
    edge_list_dir,
    biome_metadata_file,
    pfam_metadata_file,
    min_seq_id,
    min_coverage,
    group_parts,
    out_dir,
    color_metadata_file=None,
):
    """Process TSV edge lists using DuckDB to generate node and edge parquet files."""
    out_dir.mkdir(parents=True, exist_ok=True)

    conn = setup_database(biome_metadata_file, pfam_metadata_file, color_metadata_file)
    try:
        setup_biome_colors(conn, group_parts)
        tsv_files = list(edge_list_dir.glob("*.tsv")) + list(edge_list_dir.glob("*.tsv.gz"))
        print(f"\tFound {len(tsv_files)} edge lists to process.")
        for tsv_path in tsv_files:
            process_single_network(conn, tsv_path, out_dir, min_seq_id, min_coverage)
    finally:
        conn.close()
    print(
        f"[{time.strftime('%H:%M:%S')}] Successfully exported all files into '{out_dir}'."
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # Helper function to validate float ranges between 0.0 and 1.0
    def valid_fraction(arg_name):
        def validator(x):
            try:
                val = float(x)
            except ValueError:
                raise argparse.ArgumentTypeError(f"'{x}' is not a valid float.")
            if not (0.0 <= val <= 1.0):
                raise argparse.ArgumentTypeError(
                    f"{arg_name} must be between 0.0 and 1.0 (got {val})."
                )
            return val

        return validator

    parser.add_argument(
        "-e",
        "--edge-list-dir",
        type=Path,
        required=True,
        help="Directory containing edge lists",
    )
    parser.add_argument(
        "-b",
        "--biome-metadata-file",
        type=Path,
        help="Optional Parquet file containing the GOLD biome classifications",
    )
    parser.add_argument(
        "-p",
        "--pfam-metadata-file",
        type=Path,
        help="Optional Parquet file containing Pfam accessions",
    )
    parser.add_argument(
        "-C",
        "--color-metadata-file",
        type=Path,
        help="Optional Parquet with accession and color columns; overrides biome colors",
    )
    parser.add_argument(
        "-i",
        "--min-seq-id",
        type=valid_fraction("min-seq-id"),
        required=True,
        help="Minimum sequence identity threshold to keep an edge (0.0 to 1.0)",
    )
    parser.add_argument(
        "-c",
        "--min-coverage",
        type=valid_fraction("min-coverage"),
        required=True,
        help="Minimum per-HSP coverage of both sequences to keep an edge (0.0 to 1.0)",
    )
    parser.add_argument(
        "-g",
        "--color-group-parts",
        type=lambda x: (
            int(x)
            if int(x) >= 2
            else (_ for _ in ()).throw(argparse.ArgumentTypeError("Must be >= 2"))
        ),
        default=4,
        help="Number of leading lineage parts for color grouping",
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
        biome_metadata_file=args.biome_metadata_file,
        pfam_metadata_file=args.pfam_metadata_file,
        color_metadata_file=args.color_metadata_file,
        min_seq_id=args.min_seq_id,
        min_coverage=args.min_coverage,
        group_parts=args.color_group_parts,
        out_dir=args.out_dir,
    )
