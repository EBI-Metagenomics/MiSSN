import pandas as pd
import networkx as nx
import numpy as np
import os
import glob
import gc

edge_list_dir = "cluster_alignments"
# TODO: param for this
nodes_file = "dataset/gsoc_2026_test_set_nodes.csv"
output_dir = "networks"

os.makedirs(output_dir, exist_ok=True)

nodes_df = pd.read_csv(nodes_file, header=None, index_col=0, dtype=str).fillna("")

tsv_files = glob.glob(os.path.join(edge_list_dir, "*.tsv"))
print(f"Found {len(tsv_files)} edge list tsv files to process.")

for tsv_path in tsv_files:
    base_name = os.path.basename(tsv_path).replace(".tsv", "")

    # The input is in TSV
    edges_df = pd.read_csv(tsv_path, sep="\t", header=None)

    edges_df.rename(
        columns={0: "source", 1: "destination", 2: "sequence_identity"}, inplace=True
    )

    # Remove self edges
    # TODO: these are already removed by --no-self-hits
    edges_df = edges_df[edges_df["source"] != edges_df["destination"]]

    # Unify undirected edges
    edges_df[["node_u", "node_v"]] = np.sort(
        edges_df[["source", "destination"]], axis=1
    )

    # Keep only largest edges
    edges_df = edges_df.sort_values("sequence_identity", ascending=False)
    edges_df = edges_df.drop_duplicates(subset=["node_u", "node_v"], keep="first")

    # TODO: keep only edges with condition met

    G = nx.from_pandas_edgelist(
        edges_df,
        source="node_u",
        target="node_v",
        edge_attr="sequence_identity",
        create_using=nx.Graph(),
    )

    cluster_nodes = list(G.nodes())

    local_nodes_dict = nodes_df.loc[cluster_nodes].to_dict("index")

    nx.set_node_attributes(G, local_nodes_dict)

    network_dir = os.path.join(output_dir, base_name)
    os.makedirs(network_dir, exist_ok=True)

    # Save graphml
    nx.write_graphml(G, os.path.join(network_dir, f"{base_name}.graphml"))

    # Save csv for Cosmograph
    cosmograph_edges = edges_df[["node_u", "node_v", "sequence_identity"]].copy()
    cosmograph_edges.rename(
        columns={"node_u": "source", "node_v": "target"}, inplace=True
    )

    cosmograph_edges.to_csv(
        os.path.join(network_dir, f"{base_name}_edges.csv"), index=False
    )
    cosmograph_edges.to_parquet(
        os.path.join(network_dir, f"{base_name}_edges.parquet"), index=False
    )

    cluster_nodes = list(G.nodes())
    cosmograph_nodes = nodes_df.loc[cluster_nodes].copy().reset_index()

    cosmograph_nodes.rename(columns={cosmograph_nodes.columns[0]: "id"}, inplace=True)
    cosmograph_nodes.columns = cosmograph_nodes.columns.astype(str)

    cosmograph_nodes.to_csv(
        os.path.join(network_dir, f"{base_name}_nodes.csv"), index=False
    )
    cosmograph_nodes.to_parquet(
        os.path.join(network_dir, f"{base_name}_nodes.parquet"), index=False
    )

    del edges_df
    del G
    del local_nodes_dict
    del cosmograph_edges
    del cosmograph_nodes
    del cluster_nodes
    gc.collect()

print(f"Successfully export GraphML, CSV and Parquet files into {output_dir}")
