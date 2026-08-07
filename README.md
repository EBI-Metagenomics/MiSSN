# MiSSN: A workflow for generating MGnify interactive Sequence Similarity Networks.

The latest release of the [MGnify Proteins Database](https://www.ebi.ac.uk/metagenomics/proteins/) contains over 1.6 billion non-redundant protein records including relevant metagenomics metadata. The visualisation of sets of protein sequences using SSNs is a common approach for extracting novel insights about protein-protein relationships, including functional, structural, and evolutionary hypotheses.

To achieve this, MiSSN automates the network generation process. It uses pre-clustering to reduce sequence redundancy before computing alignments. The network data is enriched with [*GOLD Ecosystem Classification paths*](https://gold.jgi.doe.gov/ecosystem_classification) and [*Pfam accessions*](https://www.ebi.ac.uk/interpro/entry/pfam/), resulting in a context-rich tool for the exploration of MGnify proteins.

# Workflow

```mermaid
---
config:
  layout: dagre
---
flowchart TD
    %% Node Definitions
    A["dataset.fasta"]
    B["MMseqs2 linclust\nlinclust.sh"]
    C["dataset_biome_annotations.parquet"]
    D["linclust_results/\ndataset_cluster.tsv"]
    E["separated_clusters/\nFASTA per cluster"]
    F["Cluster size distribution notebook\nvisualise_size_distribution.ipynb"]
    G["separate_clusters.py"]
    H["DIAMOND blastp\nall-vs-all.sh"]
    I["cluster_alignments/\nTSV edge list\nper cluster"]
    K["build_ssn.py"]
    M["networks/\nedges.parquet\nnodes.parquet\nper cluster"]
    N["SSN visualisation notebook\nvisualise_ssn.ipynb"]
    O["dataset_pfam_annotations.parquet"]

    %% Shapes
    G@{ shape: rect}
    I@{ shape: rect}
    N@{ shape: rect}
    O@{ shape: rect}

    %% Connections / Flow
    A --> B
    B --> D
    D --> F
    D --> G
    G --> E
    E --> H
    H --> I
    I --> K
    C --> K
    O --> K
    K --> M
    M --> N

    %% Styles
    style A fill:#ecfeff,stroke:#22d3ee,stroke-width:2px,color:#0f172a
    style B fill:#f0fdf4,stroke:#4ade80,stroke-width:2px,color:#0f172a
    style C fill:#ecfeff,stroke:#22d3ee,stroke-width:2px,color:#0f172a
    style D fill:#fff7ed,stroke:#fb923c,stroke-width:2px,color:#0f172a
    style E fill:#fff7ed,stroke:#fb923c,stroke-width:2px,color:#0f172a
    style F fill:#f5f3ff,stroke:#a78bfa,stroke-width:2px,color:#0f172a
    style G fill:#f0fdf4,stroke:#4ade80,stroke-width:2px,color:#0f172a
    style H fill:#f0fdf4,stroke:#4ade80,stroke-width:2px,color:#0f172a
    style I fill:#fff7ed,stroke:#fb923c,stroke-width:2px,color:#0f172a
    style K fill:#f0fdf4,stroke:#4ade80,stroke-width:2px,color:#0f172a
    style M fill:#fff7ed,stroke:#fb923c,stroke-width:2px,color:#0f172a
    style N fill:#f5f3ff,stroke:#a78bfa,stroke-width:2px,color:#0f172a
    style O fill:#ecfeff,stroke:#22d3ee,stroke-width:2px,color:#0f172a
```
---

The pipeline consists of four major steps:

1. **Pre-clustering:** Utilizing MMseqs2 (`easy-linclust`) to perform linear-time clustering on the input dataset using **minimum sequence identity** (`-i`) and **minimum coverage** (`-c`) thresholds, reducing sequence redundancy before computing alignments.
2. **Separation:** Filtering the clustering results by a **minimum cluster size** (`-s`) and splitting the retained groups into individual FASTA files.
3. **Calculating all-vs-all alignments**: Running DIAMOND `blastp` on each separated cluster to compute all-vs-all pairwise sequence alignments, outputting the edge lists for the networks.
4. **Building the SSNs**: Filtering the edge lists by the same **minimum sequence identity** and **minimum coverage** thresholds, then enriching the network nodes with metadata—specifically *GOLD Ecosystem Classification paths* and *Pfam accessions*—to generate the final network edge and node Parquet files. This step calculates node colors using the **color group parts** (`-g`) parameter, which controls how many levels deep into the biome hierarchy base colors are assigned (with deeper sub-biomes inheriting shades of their parent).

> [!NOTE]
> **Why we filter twice:** MMseqs2 `easy-linclust` relies on a fast, heuristic algorithm to achieve linear-time scaling. While highly efficient for initial dataset reduction, its sequence identity and coverage boundaries are only approximate. Applying the exact same thresholds again in to the precise DIAMOND alignments guarantees that the final network edges strictly enforce the defined parameters.

The generated SSNs can be visualised using Jupyter notebooks found under the `/notebooks` directory. We provide the following interactive environments to tune your parameters and help you analyze your data:

1. **Cluster size distribution notebook** (`visualise_size_distribution.ipynb`): This notebook is designed to analyze the cluster size distribution, helping you accurately set or adjust the `<min_cluster_size>` parameter for the main pipeline.

2. **SSN visualisation notebook** (`visualise_ssn.ipynb`): An interactive environment to load the generated `.parquet` network files, explore the sequence similarity networks visually, and interactively filter nodes by their biome and Pfam annotations.

# Requirements

This project manages dependencies using `conda`. The environment includes Python libraries for data processing and visualization, alongside compiled bioinformatics tools.

## Prerequisites

You will need to have [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) (or [Mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)) installed on your system.

## First time setup

You can easily create the required environment using the provided `environment.yml` file. Run the following commands in your terminal:

```bash
# Create the conda environment
conda env create -f environment.yml

# Activate the environment
conda activate missn
```

## Key dependencies

- Python 3.13
- Bioinformatics tools: MMseqs2, DIAMOND, Biopython
- Data processing: DuckDB, Pandas, PyArrow, NumPy
- Visualisation & UI: Cosmograph, Matplotlib, Jupyter Widgets

# Usage

Clone the repository and make sure to use an environment that fulfilles the dependencies.

```bash
git clone https://github.com/vid-szabi/MiSSN.git
cd MiSSN
```

You can run the full pipeline using the following command:

```bash
./MiSSN.sh -f <fasta_file> -m <metadata_file> -i <min_seq_id> -c <min_coverage> -s <min_cluster_size> [-o <out_dir>] [-g <color_group_parts>]
```

## Inputs

### Required arguments

- `-f <fasta_file>`: Path to the input protein sequence file in FASTA format (`.fasta`, `.fa`, `.faa`), supporting gzip-compressed inputs (`.gz`).

- `-m <metadata_file>`: Path to the metadata file containing biome annotations, provided in Apache Parquet (`.parquet`) format.

- `-i <min_seq_id>`: Minimum sequence identity threshold required to keep an edge between nodes in the network, specified as a decimal float from `0.0` to `1.0`.

- `-c <min_coverage>`: Minimum alignment coverage threshold required to keep an edge between nodes in the network, specified as a decimal float from `0.0` to `1.0`.

- `-s <min_cluster_size>`: Minimum number of connected members required within a component to extract and output a distinct cluster.

### Optional arguments

- `-o <out_dir>`: Custom path for the output directory where network files and reports will be saved. Defaults to a dynamically generated name based on the input FASTA filename.

- `-g <color_group_parts>`: Number of leading levels of the biome lineage to use for grouping and color-coding nodes within the network visualization. (Default: `4`).

- `-h`: Display the help message and exit the program. 

> [!IMPORTANT]
> The input dataset must use the `.fasta`, `.fa`, or `.faa` extension. Gzip-compressed versions (e.g., `.fasta.gz`) are also supported and will be decompressed on the fly. The metadata files containing the biome and Pfam annotations must be provided in `.parquet` format.

## Outputs

The filenames are prefixed with the sequence ID of the cluster's representative member.

- `<cluster_rep>_nodes.parquet`: Contains the node attributes for the specific cluster, including the sequence IDs, biome lineage, colors (based on the biomes), and Pfam annotations.

| id | biome | color | pfam |
| :--- | :--- | :--- | :--- |
| MGYP010172319147 | root:Environmental:Aquatic:Freshwater | rgb(75,221,153) | TBA |
| MGYP000588563658 | root:Environmental:Aquatic:Estuary:Sediment | rgb(121,226,168)| TBA |

- `<cluster_rep>_edge.parquet`: Contains the network edges between the sequences in that cluster, including sequence identity.

| source | target | sequence identity |
| :--- | :--- | :--- |
| MGYP010172319147 | MGYP011046998155 | 77.2 |
| MGYP000588563658 | MGYP010172319147 | 78.1 |

> [!TIP]
> These output files are formatted ready for visualization. You can load them directly into the provided `visualise_ssn.ipynb` notebook, or import them into the Cosmograph web app or Cytoscape.

# Resources

Example sequence data and metadata annotations are available in the `/dataset` directory allowing you to evaluate the pipeline immediately after installation. All provided data is derived from the MGnify Protein Database release 2026_07.

# Credits

MiSSN was originally written as a [Google Summer of Code 2026 project](https://summerofcode.withgoogle.com/programs/2026/projects/SzccYhba) by [Szabolcs Vidám](https://github.com/vid-szabi) under the mentoring of [Christian Atallah](https://github.com/chrisAta) and [Ekaterina Sakharova](https://github.com/KateSakharova).

# Citation

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

# Licensing

All of the source code making up MiSSN in this repository is licensed under the terms of the [Apache 2.0 license](LICENSE). The example datasets provided are licensed under the terms of the [CC0 1.0 Universal (CC0 1.0) license](https://creativecommons.org/publicdomain/zero/1.0/).
