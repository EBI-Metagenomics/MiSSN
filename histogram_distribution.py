import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

CLUSTER_COL = 0 # In MMSeqs2 Linclust the cluster.tsv contains the cluster_col in the first column
OUT_DIR = 'plots'

def plot_cluster_distribution(file_path):
    base_name = os.path.basename(file_path)
    prefix = base_name.split('_')[0]
    os.makedirs(OUT_DIR, exist_ok=True)
    out_filename = f"{prefix}_seq_id_distribution_histogram1.png"
    out_path = os.path.join(OUT_DIR, out_filename)
    
    df = pd.read_csv(file_path, sep='\t', header=None) # TSV
        
    cluster_sizes = df.groupby(CLUSTER_COL).size()

    non_singletons = cluster_sizes[cluster_sizes > 1]
    print(f"Number of clusters: {cluster_sizes.size}")
    print(f"Smallest non-singleton: {non_singletons.min()}")
    print(f"Largest non-singleton:  {non_singletons.max()}")
    print(f"Mean size:              {cluster_sizes.mean():.2f}")
    print(f"Median size:            {cluster_sizes.median():.2f}")

    plt.style.use('seaborn-v0_8-dark-palette')
    fig, ax = plt.subplots(figsize=(10, 6))

    bins = np.logspace(np.log10(1), np.log10(cluster_sizes.max()), 30)

    # Plot histogram of sizes
    ax.hist(cluster_sizes, bins=bins, histtype='step', linewidth=2,
             color='skyblue', align='left')

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel('Cluster size (log scale)')
    ax.set_ylabel('Frequency (log scale)')
    
    ax.grid(axis='y', linestyle='--', alpha=0.7)

    plt.title(base_name)

    plt.tight_layout()
    plt.savefig(out_path)
    print(f"Successfully saved plot to: {out_path}")
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a histogram from an MMSeqs2 cluster TSV file.")
    parser.add_argument("input_file", help="Path to the input cluster file (e.g., linclust_results/040_seq_id_cluster.tsv)")
    
    args = parser.parse_args()
    
    plot_cluster_distribution(args.input_file)