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
    out_filename = f"{prefix}_seq_id_distribution_boxplot_min_5.png"
    out_path = os.path.join(OUT_DIR, out_filename)
    
    df = pd.read_csv(file_path, sep='\t', header=None) # TSV

    all_cluster_sizes = df.groupby(CLUSTER_COL).size()

    cluster_sizes = all_cluster_sizes[all_cluster_sizes >= 5]

    print(f"Total number of clusters (before filtering): {all_cluster_sizes.size}")
    print(f"Number of clusters (>= 5 members): {cluster_sizes.size}")

    print(f"Smallest size:  {cluster_sizes.min()}")
    print(f"Largest size:   {cluster_sizes.max()}")
    print(f"Mean size:      {cluster_sizes.mean():.2f}")
    print(f"Median size:    {cluster_sizes.median():.2f}")

    plt.style.use('seaborn-v0_8-dark-palette')
    fig, ax = plt.subplots(figsize=(10, 4))

    box_style = dict(facecolor='#4C72B0', color='#4C72B0', linewidth=1.5) 
    median_style = dict(color='#C44E52', linewidth=2)      
    whisker_style = dict(color='#4C72B0', linewidth=1.5)
    capprops_style = dict(color='#4C72B0', linewidth=1.5)

    # Plot horizontal box plot of sizes
    ax.boxplot(cluster_sizes, 
               vert=False, 
               patch_artist=True, 
               widths=0.4,
               showfliers=False,
               boxprops=box_style, 
               medianprops=median_style, 
               whiskerprops=whisker_style,
               capprops=capprops_style,)

    # To see the median line on the left
    ax.set_xlim(left=0.5)
    
    ax.set_xlabel('Cluster size (>= 5 members)')
    
    ax.set_yticks([]) 
    ax.set_ylabel('')
    
    ax.grid(axis='x', linestyle='--', alpha=0.5)

    plt.title(base_name)

    plt.tight_layout()
    plt.savefig(out_path)
    print(f"Successfully saved plot to: {out_path}")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a box plot from an MMSeqs2 cluster TSV file.")
    parser.add_argument("input_file", help="Path to the input cluster file (e.g., linclust_results/040_seq_id_cluster.tsv)")
    
    args = parser.parse_args()
    
    plot_cluster_distribution(args.input_file)
