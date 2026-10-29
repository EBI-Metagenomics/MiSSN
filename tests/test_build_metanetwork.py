import subprocess
import sys
from collections import Counter
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from scripts.build_metanetwork import build_metanetwork

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/build_metanetwork.py"
NODE_SCHEMA = pa.schema([
    ("id", pa.string()), ("biome", pa.string()),
    ("pfam_accession", pa.string()), ("color", pa.string()),
])
EDGE_SCHEMA = pa.schema([
    ("source", pa.string()), ("target", pa.string()), ("sequence_identity", pa.float64()),
])
NODES = [
    {"id": "A", "biome": None, "pfam_accession": None, "color": "rgb(200, 200, 200)"},
    {"id": "B", "biome": "root:Environmental", "pfam_accession": "PF001",
     "color": "rgb(255,0,0)"},
]
EDGES = [{"source": "A", "target": "B", "sequence_identity": 90.0}]


def write_cluster(out_dir, cluster, nodes=NODES, edges=EDGES):
    folder = out_dir / "networks" / cluster
    folder.mkdir(parents=True)
    for network_type, rows, schema in (("nodes", nodes, NODE_SCHEMA), ("edges", edges, EDGE_SCHEMA)):
        pq.write_table(pa.Table.from_pylist(rows, schema=schema), folder / f"cluster_{network_type}.parquet")


def row_counts(table):
    return Counter(tuple(row.values()) for row in table.to_pylist())


def test_merge_preserves_contents_and_sources(tmp_path):
    out = tmp_path / "output with spaces and 'quotes'"
    write_cluster(out, "group=one")
    write_cluster(out, "group=two")
    # Neither unrelated Parquets nor root-level files are cluster inputs.
    (out / "networks/group=one/unrelated.parquet").write_text("ignored")
    (out / "networks/unused_nodes.parquet").write_text("ignored")
    (out / "unrelated_edges.parquet").write_text("ignored")
    sources = {path: path.read_bytes() for path in (out / "networks").glob("*/*_*.parquet")}

    for _ in range(2):
        build_metanetwork(out)
        for network_type, rows, schema in (("nodes", NODES, NODE_SCHEMA), ("edges", EDGES, EDGE_SCHEMA)):
            output = out / f"metanetwork_{network_type}.parquet"
            table = pq.read_table(output)
            assert table.schema == schema
            assert row_counts(table) == row_counts(pa.Table.from_pylist(rows * 2, schema=schema))
            assert not (out / "networks" / output.name).exists()
    assert all(path.read_bytes() == content for path, content in sources.items())


@pytest.mark.parametrize("empty_edges", [False, True])
def test_single_cluster(tmp_path, empty_edges):
    edges = [] if empty_edges else EDGES
    write_cluster(tmp_path, "cluster", edges=edges)
    build_metanetwork(tmp_path)
    assert pq.read_table(tmp_path / "metanetwork_nodes.parquet").to_pylist() == NODES
    table = pq.read_table(tmp_path / "metanetwork_edges.parquet")
    assert table.schema == EDGE_SCHEMA
    assert table.to_pylist() == edges


def test_empty_nodes_and_edges(tmp_path):
    write_cluster(tmp_path, "cluster", nodes=[], edges=[])
    build_metanetwork(tmp_path)
    for network_type, schema in (("nodes", NODE_SCHEMA), ("edges", EDGE_SCHEMA)):
        table = pq.read_table(tmp_path / f"metanetwork_{network_type}.parquet")
        assert table.schema == schema
        assert table.num_rows == 0


@pytest.mark.parametrize("missing", ["output", "networks", "nodes", "edges", "both"])
def test_missing_inputs(tmp_path, missing):
    out = tmp_path / "output"
    if missing != "output":
        out.mkdir()
    if missing in {"nodes", "edges", "both"}:
        (out / "networks").mkdir()
    if missing in {"nodes", "edges"}:
        cluster = out / "networks/cluster"
        cluster.mkdir()
        network_type, schema = ("edges", EDGE_SCHEMA) if missing == "nodes" else ("nodes", NODE_SCHEMA)
        pq.write_table(pa.Table.from_pylist([], schema=schema), cluster / f"cluster_{network_type}.parquet")
    with pytest.raises(FileNotFoundError, match="[Nn]etworks|No cluster"):
        build_metanetwork(out)
    assert not (out / "metanetwork_nodes.parquet").exists()
    assert not (out / "metanetwork_edges.parquet").exists()


@pytest.mark.parametrize("network_type", ["nodes", "edges"])
@pytest.mark.parametrize("invalid", ["malformed", "incompatible_schema"])
def test_invalid_inputs(tmp_path, network_type, invalid):
    write_cluster(tmp_path, "first")
    write_cluster(tmp_path, "second")
    path = tmp_path / "networks/second" / f"cluster_{network_type}.parquet"
    if invalid == "malformed":
        path.write_text("not parquet")
    else:
        pq.write_table(pa.table({"wrong_column": ["wrong"]}), path)
    with pytest.raises(duckdb.Error):
        build_metanetwork(tmp_path)


@pytest.mark.parametrize("flag", ["-o", "--out-dir"])
def test_cli(tmp_path, flag):
    out = tmp_path / "output with spaces"
    write_cluster(out, "one")
    result = subprocess.run(
        [sys.executable, str(SCRIPT), flag, str(out)],
        check=True, capture_output=True, text=True,
    )
    assert "Merged 1 nodes files" in result.stdout
    assert pq.read_table(out / "metanetwork_nodes.parquet").to_pylist() == NODES
    assert pq.read_table(out / "metanetwork_edges.parquet").to_pylist() == EDGES


def test_cli_missing_input(tmp_path):
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "-o", str(tmp_path)],
        check=False, capture_output=True, text=True,
    )
    assert result.returncode == 1
    assert "Networks directory not found" in result.stderr


def test_cli_requires_output_directory():
    result = subprocess.run(
        [sys.executable, str(SCRIPT)], check=False, capture_output=True, text=True,
    )
    assert result.returncode == 2
    assert "--out-dir" in result.stderr
