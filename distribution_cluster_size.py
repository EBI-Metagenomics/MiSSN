import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

CLUSTER_COL=0

def plot_cluster_distribution(file_path):
    df = pd.read_csv(file_path, sep='\t', header=None) # TSV
        
    cluster_sizes = df.groupby(CLUSTER_COL).size()

    non_singletons = cluster_sizes[cluster_sizes > 1]
    print(f"Number of clusters: {cluster_sizes.size}")
    print(f"Smallest non-singleton: {non_singletons.min()}")
    print(f"Largest non-singleton:  {non_singletons.max()}")
    print(f"Mean size:              {cluster_sizes.mean():.2f}")
    print(f"Median size:            {cluster_sizes.median():.2f}")

    fig, ax = plt.subplots(figsize=(10, 6))

    bins = np.logspace(np.log10(1), np.log10(cluster_sizes.max()), 30)

    # Plot histogram of sizes
    ax.hist(cluster_sizes, bins=bins,
             color='skyblue', edgecolor='black', align='left')

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel('Cluster size (log scale)')
    ax.set_ylabel('Frequency (log scale)')
    
    ax.grid(axis='y', linestyle='--', alpha=0.7)

    plt.savefig('cluster_distribution.png')
    
plot_cluster_distribution('cluster_results/cluster_results_options_cluster.tsv')
