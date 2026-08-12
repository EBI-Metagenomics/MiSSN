import pytest
import pandas as pd

from scripts.build_ssn import (
    build_tree,
    setup_database,
    setup_biome_colors,
    process_single_network,
    process_networks,
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
    # First has one Pfam, second has multiple, third is missing (should default to None)
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
    assert mgyp3_row["pfam_accession"] == "None"
    
    conn.close()


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
        # 1. Valid edge (pident 87.4% >= 40%, coverage: 350 / 427 = 82.0% >= 80%)
        "MGYP001208449037\tMGYP001342377832\t87.4\t286\t427\t296\t350\t6.73e-186\t512.0\n"
                
        # 2. Fails min_seq_id (pident 30.0% < 40%, coverage: 350 / 427 = 82.0% >= 80%)
        "MGYP001208449037\tMGYP006224946903\t30.0\t167\t427\t168\t350\t7.27e-66\t201.0\n"    
        
        # 3. Fails min_coverage (pident 90.0% >= 40%, coverage: 50 / 427 = 11.7% < 80%)
        "MGYP001208449037\tMGYP006224946903\t90.0\t100\t427\t168\t50\t1e-5\t200.0\n"
        
        # 4. Self loop (should be removed by the SQL logic)
        "MGYP001342377832\tMGYP001342377832\t100.0\t296\t296\t296\t296\t0.0\t500.0\n"
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
    

def test_integration_process_networks(tmp_path, mock_metadata):
    """End-to-end integration test of the main processing loop."""
    biome_metadata_file, pfam_metadata_file = mock_metadata
    
    # Setup edge directory
    edge_dir = tmp_path / "edges"
    edge_dir.mkdir()
        
    # cluster_A valid edge: (pident 87.4 >= 40%, coverage: 350 / 427 = 82.0% >= 80%)
    (edge_dir / "cluster_A.tsv").write_text("MGYP001208449037\tMGYP001342377832\t87.4\t286\t427\t296\t350\t6.73e-186\t512.0\n")
    
    # cluster_B valid edge: (pident 85.0 >= 40%, coverage: 350 / 427 = 82.0% >= 80%)
    (edge_dir / "cluster_B.tsv").write_text("MGYP001208449037\tMGYP006224946903\t85.0\t150\t427\t168\t350\t1e-80\t250.0\n")
    
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