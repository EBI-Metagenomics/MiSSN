import itertools
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from scripts.build_ssn import (
    build_tree,
    process_networks,
    process_single_network,
    setup_biome_colors,
    setup_database,
)


def test_build_tree():
    """Test the hierarchical tree construction from lineage strings."""
    lineages = [
        "root:Host-associated:Human:Digestive system:Oral:Saliva",
        "root:Host-associated:Birds:Digestive system:Digestive tube:Cecum",
        "root:Environmental:Aquatic:Marine",
        "root"
    ]
    root = build_tree(lineages)
    
    assert root.name == "root"
    assert root.depth == 0
    assert root.is_lineage is True  # True because "root" is explicitly in the input list
    
    # Check top-level children
    assert "Host-associated" in root.children
    assert "Environmental" in root.children
    
    # Check the Environmental branch
    env_node = root.children["Environmental"]
    assert env_node.depth == 1
    assert env_node.is_lineage is False  # Passed through, but not an endpoint
    
    # Check the Host-associated branch
    host_node = root.children["Host-associated"]
    assert "Human" in host_node.children
    assert "Birds" in host_node.children
    
    human_saliva_node = (
        host_node.children["Human"]
        .children["Digestive system"]
        .children["Oral"]
        .children["Saliva"]
    )
    assert human_saliva_node.depth == 5
    assert human_saliva_node.is_lineage is True
    
    
@pytest.fixture
def mock_metadata(tmp_path):
    """Fixture to generate temporary biome and Pfam Parquet files."""
    metadata_dir = tmp_path / "metadata"
    metadata_dir.mkdir()
    
    biome_metadata_file = metadata_dir / "biome.parquet"
    # Note first one has multiple biomes to test string aggregation
    pd.DataFrame({
        "mgyp": ["MGYP001342377832", "MGYP001342377832", "MGYP001208449037", "MGYP006224946903"],
        "biome": ["root:Environmental:Aquatic:Marine", "root:Environmental:Aquatic:Marine:Oceanic", "root:Environmental:Aquatic:Marine:Oceanic", "root:Environmental:Aquatic:Marine"]
    }).to_parquet(biome_metadata_file)
    
    pfam_metadata_file = metadata_dir / "pfam.parquet"
    # First has one Pfam, second has multiple, third is missing (should be null).
    pd.DataFrame({
        "mgyp": ["MGYP001342377832", "MGYP001208449037", "MGYP001208449037", "MGYP001208449037"],
        "pfam_accession": ["PF04724", "PF00293", "PF28549", "PF04724"]
    }).to_parquet(pfam_metadata_file)
    
    return biome_metadata_file, pfam_metadata_file


def test_setup_database(mock_metadata):
    """Test loading Parquets, string aggregation, and LEFT JOINs in DuckDB."""
    biome_metadata_file, pfam_metadata_file = mock_metadata
    conn = setup_database(biome_metadata_file, pfam_metadata_file)
    
    df = conn.execute("SELECT * FROM global_nodes ORDER BY id").df()
    
    # Verify the expected columns are present
    assert set(df.columns) == {"id", "biome", "pfam_accession"}
    
    # Check aggregation and joining (3 unique IDs)
    assert len(df) == 3
    
    # Check first ID (multiple biomes, single Pfam)
    mgyp1_row = df[df["id"] == "MGYP001342377832"].iloc[0]
    assert mgyp1_row["biome"] == "root:Environmental:Aquatic:Marine;root:Environmental:Aquatic:Marine:Oceanic"
    assert mgyp1_row["pfam_accession"] == "PF04724"
    
    # Check second ID (multiple Pfams sorted alphabetically)
    mgyp2_row = df[df["id"] == "MGYP001208449037"].iloc[0]
    assert mgyp2_row["pfam_accession"] == "PF00293;PF04724;PF28549"
    
    # Check third ID (missing Pfam)
    mgyp3_row = df[df["id"] == "MGYP006224946903"].iloc[0]
    assert pd.isna(mgyp3_row["pfam_accession"])
    
    conn.close()


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA_COMBINATIONS = list(itertools.product((False, True), repeat=3))


@pytest.fixture
def optional_metadata_case(tmp_path):
    """Include partial annotations, duplicate rows, and unused metadata IDs."""
    edge_dir = tmp_path / "edge lists"
    edge_dir.mkdir()
    rows = []
    for target, identity in (("B", 90), ("C", 90), ("D", 90), ("E", 10)):
        rows.append(f"A\t{target}\t{identity}\t100\t100\t100\t90\t1e-20\t100\t100\t100\n")
    (edge_dir / "cluster.tsv").write_text("".join(rows))
    # Verify connection-level setup can be reused across multiple networks.
    (edge_dir / "second.tsv").write_text(rows[0])
    biome = tmp_path / "biome annotations.parquet"
    pd.DataFrame({
        "mgyp": ["A", "A", "A", "B", "UNUSED"],
        "biome": ["root:Environmental", "root:Host-associated", "root:Environmental",
                  "root:Environmental", "root:Environmental"],
    }).to_parquet(biome)
    pfam = tmp_path / "pfam annotations.parquet"
    pd.DataFrame({
        "mgyp": ["A", "A", "A", "C", "UNUSED"],
        "pfam_accession": ["PF002", "PF001", "PF001", "PF003", "PF004"],
    }).to_parquet(pfam)
    colors = tmp_path / "color annotations.parquet"
    pd.DataFrame({
        "accession": ["A", "A", "B", "C", "UNUSED"],
        "color": ["#ff0000", "rgb(255, 0, 0)", "rgb(0, 128, 255)", "#00FF00", "#ffffff"],
    }).to_parquet(colors)
    return edge_dir, biome, pfam, colors


@pytest.mark.parametrize("use_biome,use_pfam,use_colors", METADATA_COMBINATIONS)
@pytest.mark.parametrize("via_cli", [False, True], ids=["python", "cli"])
def test_optional_metadata_combinations(
    tmp_path, optional_metadata_case, use_biome, use_pfam, use_colors, via_cli
):
    edge_dir, biome, pfam, colors = optional_metadata_case
    out = tmp_path / "network output"
    paths = [biome if use_biome else None, pfam if use_pfam else None,
             colors if use_colors else None]
    if via_cli:
        command = [sys.executable, str(REPO_ROOT / "scripts/build_ssn.py"),
                   "-e", str(edge_dir), "-i", "0.4", "-c", "0.8", "-o", str(out)]
        color_flag = "-C" if use_biome else "--color-metadata-file"
        for flag, path in zip(("-b", "-p", color_flag), paths):
            if path is not None:
                command.extend([flag, str(path)])
        subprocess.run(command, check=True, capture_output=True, text=True)
    else:
        process_networks(edge_dir, paths[0], paths[1], 0.4, 0.8, 4, out, paths[2])

    nodes_path = out / "cluster/cluster_nodes.parquet"
    nodes = pd.read_parquet(nodes_path).set_index("id")
    assert set(nodes.index) == {"A", "B", "C", "D", "E"}
    assert list(nodes.columns) == ["biome", "pfam_accession", "color"]
    schema = pq.read_schema(nodes_path)
    for name in ("id", "biome", "pfam_accession", "color"):
        assert pa.types.is_string(schema.field(name).type)
    if use_biome:
        assert nodes.loc["A", "biome"] == "root:Environmental;root:Host-associated"
        assert nodes.loc["B", "biome"] == "root:Environmental"
    else:
        assert nodes.biome.isna().all()
    assert nodes.loc[["C", "D", "E"], "biome"].isna().all()
    if use_pfam:
        assert nodes.loc["A", "pfam_accession"] == "PF001;PF002"
        assert nodes.loc["C", "pfam_accession"] == "PF003"
    else:
        assert nodes.pfam_accession.isna().all()
    assert nodes.loc[["B", "D", "E"], "pfam_accession"].isna().all()
    if use_colors:
        assert nodes.loc["A", "color"] == "rgb(255,0,0)"
        assert nodes.loc["B", "color"] == "rgb(0,128,255)"
        assert nodes.loc["C", "color"] == "rgb(0,255,0)"
    else:
        assert nodes.loc["A", "color"] == ("rgb(0,0,0)" if use_biome else "rgb(200, 200, 200)")
        if use_biome:
            assert nodes.loc["B", "color"] not in {"rgb(0,0,0)", "rgb(200, 200, 200)"}
        else:
            assert nodes.loc["B", "color"] == "rgb(200, 200, 200)"
        assert nodes.loc["C", "color"] == "rgb(200, 200, 200)"
    assert set(nodes.loc[["D", "E"], "color"]) == {"rgb(200, 200, 200)"}
    edges = pd.read_parquet(out / "cluster/cluster_edges.parquet")
    assert set(zip(edges.source, edges.target, edges.sequence_identity)) == {
        ("A", "B", 90), ("A", "C", 90), ("A", "D", 90)
    }
    second = pd.read_parquet(out / "second/second_nodes.parquet").set_index("id")
    pd.testing.assert_frame_equal(second.sort_index(), nodes.loc[["A", "B"]].sort_index())


@pytest.mark.parametrize("use_biome,use_pfam,use_colors", METADATA_COMBINATIONS)
def test_empty_metadata(tmp_path, optional_metadata_case, use_biome, use_pfam, use_colors):
    edge_dir, biome, pfam, colors = optional_metadata_case
    for path, columns in ((biome, ["mgyp", "biome"]),
                          (pfam, ["mgyp", "pfam_accession"]),
                          (colors, ["accession", "color"])):
        pd.DataFrame({name: pd.Series(dtype="string") for name in columns}).to_parquet(path)
    out = tmp_path / "empty metadata output"
    process_networks(edge_dir, biome if use_biome else None, pfam if use_pfam else None,
                     0.4, 0.8, 4, out, colors if use_colors else None)
    nodes = pd.read_parquet(out / "cluster/cluster_nodes.parquet")
    assert set(nodes.id) == {"A", "B", "C", "D", "E"}
    assert nodes[["biome", "pfam_accession"]].isna().all().all()
    assert set(nodes.color) == {"rgb(200, 200, 200)"}


@pytest.mark.parametrize("data,message", [
    ({"id": ["A"], "color": ["#ffffff"]}, "accession"),
    ({"accession": ["A"]}, "color"),
    ({"accession": [1], "color": ["#ffffff"]}, "string column"),
    ({"accession": ["A"], "color": [1]}, "string column"),
    ({"accession": [None], "color": ["#ffffff"]}, "accession"),
    ({"accession": ["A", "B"], "color": [None, "#ffffff"]}, "null or blank"),
    ({"accession": [" "], "color": ["#ffffff"]}, "null or blank"),
    ({"accession": ["A"], "color": [" "]}, "null or blank"),
    ({"accession": ["A"], "color": ["red"]}, "Invalid color"),
    ({"accession": ["A"], "color": ["#fff"]}, "Invalid color"),
    ({"accession": ["A"], "color": ["rgb(256,0,0)"]}, "Invalid color"),
    ({"accession": ["A"], "color": ["rgb(-1,0,0)"]}, "Invalid color"),
    ({"accession": ["A"], "color": ["rgb(1.5,0,0)"]}, "Invalid color"),
    ({"accession": ["A", "A"], "color": ["#ffffff", "#000000"]}, "Conflicting colors"),
])
def test_invalid_color_metadata(tmp_path, data, message):
    path = tmp_path / "invalid.parquet"
    pd.DataFrame(data).to_parquet(path)
    with pytest.raises(ValueError, match=message):
        setup_database(None, None, path)


def test_partial_color_override_keeps_biome_fallback(tmp_path, optional_metadata_case):
    edge_dir, biome, pfam, colors = optional_metadata_case
    pd.DataFrame({"accession": ["C"], "color": ["#010203"]}).to_parquet(colors)
    out = tmp_path / "partial"
    process_networks(edge_dir, biome, pfam, 0.4, 0.8, 4, out, colors)
    nodes = pd.read_parquet(out / "cluster/cluster_nodes.parquet").set_index("id")
    assert nodes.loc["A", "color"] == "rgb(0,0,0)"
    assert nodes.loc["B", "color"] not in {"rgb(0,0,0)", "rgb(200, 200, 200)"}
    assert nodes.loc["C", "color"] == "rgb(1,2,3)"


def test_numeric_accessions_match_exactly(tmp_path):
    edge_dir = tmp_path / "edges"
    edge_dir.mkdir()
    (edge_dir / "numeric.tsv").write_text(
        "001\t002\t90\t100\t100\t100\t90\t1e-20\t100\t100\t100\n"
    )
    colors = tmp_path / "colors.parquet"
    pd.DataFrame({"accession": ["001", "1"], "color": ["#000000", "#ffffff"]}).to_parquet(colors)
    out = tmp_path / "out"
    process_networks(edge_dir, None, None, 0.4, 0.8, 4, out, colors)
    nodes = pd.read_parquet(out / "numeric/numeric_nodes.parquet").set_index("id")
    assert set(nodes.index) == {"001", "002"}
    assert nodes.loc["001", "color"] == "rgb(0,0,0)"
    assert nodes.loc["002", "color"] == "rgb(200, 200, 200)"
    edges = pd.read_parquet(out / "numeric/numeric_edges.parquet")
    assert list(zip(edges.source, edges.target)) == [("001", "002")]


@pytest.mark.parametrize("column", ["accession", "color"])
def test_color_metadata_rejects_typed_nulls(tmp_path, column):
    data = pd.DataFrame({"accession": ["A"], "color": ["#ffffff"]}, dtype="string")
    data[column] = pd.Series([pd.NA], dtype="string")
    path = tmp_path / "nulls.parquet"
    data.to_parquet(path)
    with pytest.raises(ValueError, match="null or blank"):
        setup_database(None, None, path)


@pytest.mark.parametrize("flag", ["-b", "-p", "-C"])
@pytest.mark.parametrize("malformed", [False, True], ids=["missing", "malformed"])
def test_cli_rejects_bad_metadata(tmp_path, flag, malformed):
    path = tmp_path / "bad.parquet"
    if malformed:
        path.write_text("not a parquet")
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "scripts/build_ssn.py"), "-e", str(tmp_path),
         "-i", "0.4", "-c", "0.8", "-o", str(tmp_path / "out"), flag, str(path)],
        check=False, capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert "parquet" in result.stderr.lower()


@pytest.fixture
def shell_pipeline(tmp_path):
    """Copy only the tracked entry point and replace expensive stages with stubs."""
    root = tmp_path / "pipeline with spaces"
    scripts = root / "scripts"
    scripts.mkdir(parents=True)
    shutil.copyfile(REPO_ROOT / "MiSSN.sh", root / "MiSSN.sh")
    for name in ("linclust.sh", "separate_clusters.py", "all-vs-all.sh"):
        stub = scripts / name
        stub.write_text('#!/bin/bash\nexit 0\n')
        stub.chmod(0o755)
    log = tmp_path / "arguments.json"
    builder = scripts / "build_ssn.py"
    builder.write_text(
        f"#!{sys.executable}\nimport json, os, sys\n"
        "with open(os.environ['MISSN_TEST_ARG_LOG'], 'w') as handle:\n"
        "    json.dump(sys.argv[1:], handle)\n"
        "sys.exit(int(os.environ.get('MISSN_TEST_BUILD_EXIT', '0')))\n"
    )
    builder.chmod(0o755)
    merger = scripts / "build_metanetwork.py"
    merger.write_text(
        "import json, os, sys\n"
        "assert os.path.exists(os.environ['MISSN_TEST_ARG_LOG'])\n"
        "with open(os.environ['MISSN_TEST_MERGE_LOG'], 'w') as handle:\n"
        "    json.dump(sys.argv[1:], handle)\n"
        "sys.exit(int(os.environ.get('MISSN_TEST_MERGE_EXIT', '0')))\n"
    )
    fasta = root / "input sequences.fasta"
    fasta.write_text(">A\nAAAA\n")
    command = ["bash", str(root / "MiSSN.sh"), "-f", str(fasta), "-i", "0.4",
               "-c", "0.8", "-s", "2", "-o", str(root / "output directory")]
    env = {**os.environ, "MISSN_TEST_ARG_LOG": str(log),
           "MISSN_TEST_MERGE_LOG": str(tmp_path / "merge_arguments.json")}
    return command, env, log


@pytest.mark.parametrize("use_biome,use_pfam,use_colors", METADATA_COMBINATIONS)
def test_shell_optional_metadata(shell_pipeline, tmp_path, use_biome, use_pfam, use_colors):
    command, env, log = shell_pipeline
    expected = {}
    for flag, supplied in zip(("-b", "-p", "-C"), (use_biome, use_pfam, use_colors)):
        if supplied:
            path = tmp_path / f"metadata {flag}.parquet"
            path.touch()
            command.extend([flag, str(path)])
            expected[flag] = str(path)
    subprocess.run(command, env=env, check=True, capture_output=True, text=True)
    args = json.loads(log.read_text())
    forwarded = dict(zip(args[::2], args[1::2]))
    for flag in ("-b", "-p", "-C"):
        if flag in expected:
            assert forwarded[flag] == expected[flag]
        else:
            assert flag not in forwarded
    assert forwarded["-g"] == "4"
    assert forwarded["-i"] == "0.4"
    assert forwarded["-c"] == "0.8"
    assert forwarded["-o"].endswith("output directory/networks")
    assert not Path(env["MISSN_TEST_MERGE_LOG"]).exists()


@pytest.mark.parametrize("position", [2, 6, 12], ids=["before", "between", "after"])
def test_shell_metanetwork(shell_pipeline, position):
    command, env, log = shell_pipeline
    command.insert(position, "--metanetwork")
    subprocess.run(command, env=env, check=True, capture_output=True, text=True)
    build_args = json.loads(log.read_text())
    build_output = dict(zip(build_args[::2], build_args[1::2]))["-o"]
    merge_args = json.loads(Path(env["MISSN_TEST_MERGE_LOG"]).read_text())
    assert merge_args == ["-o", str(Path(build_output).parent)]


@pytest.mark.parametrize("merge_enabled", [False, True])
def test_shell_build_failure(shell_pipeline, merge_enabled):
    command, env, log = shell_pipeline
    env["MISSN_TEST_BUILD_EXIT"] = "7"
    if merge_enabled:
        command.append("--metanetwork")
    result = subprocess.run(command, env=env, check=False, capture_output=True, text=True)
    assert result.returncode == 7
    assert log.exists()
    assert not Path(env["MISSN_TEST_MERGE_LOG"]).exists()


def test_shell_merge_failure(shell_pipeline):
    command, env, _ = shell_pipeline
    env["MISSN_TEST_MERGE_EXIT"] = "9"
    result = subprocess.run(command + ["--metanetwork"], env=env, check=False,
                            capture_output=True, text=True)
    assert result.returncode == 9
    assert Path(env["MISSN_TEST_MERGE_LOG"]).exists()


@pytest.mark.parametrize("flag", ["--unknown", "--metanetwork=true"])
def test_shell_rejects_unknown_long_options(shell_pipeline, flag):
    command, env, log = shell_pipeline
    result = subprocess.run(command + [flag], env=env, check=False, capture_output=True, text=True)
    assert result.returncode != 0
    assert not log.exists()


@pytest.mark.parametrize("flag", ["-b", "-p", "-C"])
@pytest.mark.parametrize("path", ["missing file.parquet", ""])
def test_shell_rejects_bad_metadata(shell_pipeline, flag, path):
    command, env, log = shell_pipeline
    result = subprocess.run(command + [flag, path], check=False, env=env, capture_output=True, text=True)
    assert result.returncode != 0
    assert not log.exists()


def test_setup_biome_colors(mock_metadata):
    """Test the assignment of colors within DuckDB, ensuring concatenated biomes get black."""
    biome_metadata_file, pfam_metadata_file = mock_metadata
    conn = setup_database(biome_metadata_file, pfam_metadata_file)
    setup_biome_colors(conn, group_parts=2)
    
    df = conn.execute("SELECT * FROM color_lookup").df()
    
    # Find all rows where the biome contains a semicolon (combined biomes)
    combined_rows = df[df["biome"].str.contains(";")]
    assert len(combined_rows) > 0
    
    # Check that every combined biome is explicitly assigned black
    assert (combined_rows["color"] == "rgb(0,0,0)").all()
    
    # Find rows without a semicolon (single biomes)
    single_rows = df[~df["biome"].str.contains(";")]
    assert len(single_rows) > 0
    
    # Check that every single biome has a computed rgb value
    assert single_rows["color"].str.startswith("rgb(").all()
    assert (single_rows["color"] != "rgb(0,0,0)").all()
    
    conn.close()


def test_process_single_network(tmp_path, mock_metadata):
    """Test TSV parsing, sequence/coverage filtering, and graph output."""
    biome_metadata_file, pfam_metadata_file = mock_metadata
    conn = setup_database(biome_metadata_file, pfam_metadata_file)
    setup_biome_colors(conn, group_parts=3)
    
    # Create a mock DIAMOND output TSV
    tsv_file = tmp_path / "cluster.tsv"
    tsv_content = (
        # 1. Valid edge: both coverages >= 80%, despite fewer positive matches.
        "MGYP001208449037\tMGYP001342377832\t87.4\t250\t300\t296\t230\t6.73e-186\t512.0\t83.3\t84.5\n"
                
        # 2. Fails min_seq_id (30% < 40%), despite sufficient coverage.
        "MGYP001208449037\tMGYP006224946903\t30.0\t150\t180\t168\t100\t7.27e-66\t201.0\t83.3\t89.3\n"
        
        # 3. Fails min_coverage (both below 80%, but above 0.8).
        "MGYP001208449037\tMGYP006224946903\t90.0\t100\t427\t168\t95\t1e-5\t200.0\t23.4\t59.5\n"
        
        # 4. Self loop (should be removed by the SQL logic)
        "MGYP001342377832\tMGYP001342377832\t100.0\t296\t296\t296\t296\t0.0\t500.0\t100\t100\n"
    )
    tsv_file.write_text(tsv_content)
    
    out_dir = tmp_path / "networks"
    process_single_network(
        conn=conn,
        tsv_path=tsv_file,
        out_dir=out_dir,
        min_seq_id=0.4,
        min_coverage=0.8,
    )
    
    cluster_dir = out_dir / "cluster"
    assert cluster_dir.exists()
    
    edges_file = cluster_dir / "cluster_edges.parquet"
    nodes_file = cluster_dir / "cluster_nodes.parquet"
    assert edges_file.exists()
    assert nodes_file.exists()

    edges_df = pd.read_parquet(edges_file)
    assert len(edges_df) == 1
    
    assert (edges_df["source"] != edges_df["target"]).all()
    
    assert set(edges_df.iloc[0][["source", "target"]]) == {"MGYP001342377832", "MGYP001208449037"}

    nodes_df = pd.read_parquet(nodes_file)
    assert len(nodes_df) == 3
    assert set(nodes_df["id"]) == {"MGYP001342377832", "MGYP001208449037", "MGYP006224946903"}
    
    conn.close()
    

def test_process_single_network_accepts_paths_with_apostrophes(tmp_path):
    """DuckDB file paths should remain valid when they contain apostrophes."""
    input_dir = tmp_path / "input's"
    input_dir.mkdir()
    tsv_file = input_dir / "cluster's.tsv"
    tsv_file.write_text(
        "A\tB\t90\t100\t100\t100\t90\t1e-20\t100\t100\t100\n"
    )

    conn = setup_database()
    setup_biome_colors(conn, group_parts=4)
    output_dir = input_dir / "output's"

    process_single_network(conn, tsv_file, output_dir, 0.4, 0.8)

    assert (output_dir / "cluster's" / "cluster's_edges.parquet").exists()
    assert (output_dir / "cluster's" / "cluster's_nodes.parquet").exists()
    conn.close()


def test_integration_process_networks(tmp_path, mock_metadata):
    """End-to-end integration test of the main processing loop."""
    biome_metadata_file, pfam_metadata_file = mock_metadata
    
    # Setup edge directory
    edge_dir = tmp_path / "edges"
    edge_dir.mkdir()
        
    # cluster_A valid edge: pident >= 40%, both coverages >= 80%.
    (edge_dir / "cluster_A.tsv").write_text("MGYP001208449037\tMGYP001342377832\t87.4\t250\t300\t296\t230\t6.73e-186\t512.0\t83.3\t84.5\n")
    
    # cluster_B valid edge: pident >= 40%, both coverages >= 80%.
    (edge_dir / "cluster_B.tsv").write_text("MGYP001208449037\tMGYP006224946903\t85.0\t150\t180\t168\t140\t1e-80\t250.0\t83.3\t89.3\n")
    
    out_dir = tmp_path / "networks"
    
    process_networks(
        edge_list_dir=edge_dir,
        biome_metadata_file=biome_metadata_file,
        pfam_metadata_file=pfam_metadata_file,
        min_seq_id=0.4,
        min_coverage=0.8,
        group_parts=3,
        out_dir=out_dir,
    )
    
    assert (out_dir / "cluster_A" / "cluster_A_nodes.parquet").exists()
    assert (out_dir / "cluster_A" / "cluster_A_edges.parquet").exists()
    assert (out_dir / "cluster_B" / "cluster_B_nodes.parquet").exists()
    assert (out_dir / "cluster_B" / "cluster_B_edges.parquet").exists()

    for cluster in ("cluster_A", "cluster_B"):
        assert len(pd.read_parquet(out_dir / cluster / f"{cluster}_edges.parquet")) == 1


@pytest.mark.parametrize(
    "qcovhsp, scovhsp, expected_edges",
    [
        (80, 80, 1),  # Inclusive threshold; positive matches alone would fail.
        (79.9, 100, 0),  # Subject coverage cannot rescue low query coverage.
        (100, 79.9, 0),  # Query coverage cannot rescue low subject coverage.
        (1, 1, 0),  # Percentages must be compared with 80, not 0.8.
    ],
)
def test_coverage_threshold(tmp_path, mock_metadata, qcovhsp, scovhsp, expected_edges):
    biome_metadata_file, pfam_metadata_file = mock_metadata
    conn = setup_database(biome_metadata_file, pfam_metadata_file)
    setup_biome_colors(conn, group_parts=3)
    tsv_file = tmp_path / "coverage.tsv"
    tsv_file.write_text(
        "MGYP001208449037\tMGYP001342377832\t50\t1000\t1000\t1000\t500\t1e-80\t250"
        f"\t{qcovhsp}\t{scovhsp}\n"
    )

    process_single_network(conn, tsv_file, tmp_path / "networks", 0.4, 0.8)

    edges = pd.read_parquet(tmp_path / "networks" / "coverage" / "coverage_edges.parquet")
    assert len(edges) == expected_edges
    conn.close()
