# MiSSN: A workflow for generating MGnify interactive Sequence Similarity Networks.

The [latest release](https://ebi-metagenomics.github.io/blog/2026/07/17/MGnify-Proteins-Release/) of the [MGnify Proteins Database](https://www.ebi.ac.uk/metagenomics/proteins/) contains over 5.7 billion non-redundant protein sequences, with over 1.6 billion cluster representatives including relevant metagenomics metadata. The visualisation of sets of protein sequences using Sequence Similarity Networks (SSNs) is a common approach for extracting novel insights about protein-protein relationships, including functional, structural, and evolutionary hypotheses.

To achieve this, MiSSN automates the network generation process. It uses pre-clustering to reduce sequence redundancy before computing all-vs-all pairwise alignments. The network can optionally be annotated with [*GOLD biome classifications*](https://gold.jgi.doe.gov/ecosystem_classification) and [*Pfam accessions*](https://www.ebi.ac.uk/interpro/entry/pfam/), resulting in a context-rich tool for the exploration of MGnify Proteins sequence clusters.

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
    C["dataset_biome_annotations.parquet (optional)"]
    D["linclust_results/\ndataset_cluster.tsv"]
    E["separated_clusters/\nFASTA per cluster"]
    F["Cluster size distribution notebook\nvisualise_size_distribution.ipynb"]
    G["separate_clusters.py"]
    H["DIAMOND blastp\nall-vs-all.sh"]
    I["cluster_alignments/\nTSV pairwise alignments\nper cluster"]
    K["build_ssn.py"]
    M["networks/\nedges.parquet\nnodes.parquet\nper cluster"]
    N["SSN visualisation notebook\nvisualise_ssn.ipynb"]
    O["dataset_pfam_annotations.parquet (optional)"]
    P["dataset_colors.parquet (optional)"]

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
    C -.-> K
    O -.-> K
    P -.-> K
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
4. **Building the SSNs**: Filtering the edge lists by the same **minimum sequence identity** and **minimum coverage** thresholds, then optionally enriching the network nodes with *GOLD biome classifications* and *Pfam accessions* to generate the final network edge and node Parquet files. This step calculates node colors using the **color group parts** (`-g`) parameter, which controls how many levels deep into the biome hierarchy base colors are assigned (with deeper sub-biomes inheriting shades of their parent). An optional accession-specific colour file (`-C`) overrides these colours for matching nodes. Nodes without a colour override or biome colour are grey.

> [!NOTE]
> **Why we filter twice:** MMseqs2 `easy-linclust` relies on a fast, heuristic algorithm to achieve linear-time scaling. While highly efficient for initial dataset reduction, its sequence identity and coverage boundaries are only approximate. Moreover, it only aligns sequences against their cluster representative, meaning it does not compute the all-vs-all pairwise identities required to build a complete graph. Running DIAMOND on the separated clusters calculates those missing all-vs-all connections, and applying the same thresholds again guarantees that the final network edges strictly enforce your defined parameters.

The generated SSNs can be visualised using Jupyter notebooks found under the `/notebooks` directory. We provide the following interactive environments to tune your parameters and help you analyze your data:

- **Cluster size distribution notebook** (`visualise_size_distribution.ipynb`): This notebook is designed to analyze the cluster size distribution, helping you accurately set or adjust the `<min_cluster_size>` parameter for the main pipeline.

- **SSN visualisation notebook** (`visualise_ssn.ipynb`): An interactive environment to load the generated `.parquet` network files, explore the sequence similarity networks visually, and interactively filter nodes by their biome and Pfam annotations. The high-performance network rendering is handled by Cosmograph.

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
./MiSSN.sh -f <fasta_file> [-b <biome_metadata_file>] [-p <pfam_metadata_file>] [-C <color_metadata_file>] -i <min_seq_id> -c <min_coverage> -s <min_cluster_size> [-o <out_dir>] [-g <color_group_parts>]
```

**Example**

This command runs the pipeline on the provided FASTA and metadata files, establishing network edges between sequences that share at least 40% sequence identity and 80% alignment coverage, and only extracting clusters that contain 5 or more members. The results will be generated under the `example_dataset/` directory.

```bash
./MiSSN.sh \
  -f example_dataset.fasta.gz \
  -b example_dataset_biome_annotations.parquet \
  -p example_dataset_pfam_annotations.parquet \
  -i 0.4 \
  -c 0.8 \
  -s 5
```

**Without annotation files**

```bash
./MiSSN.sh -f example_dataset.fasta.gz -i 0.4 -c 0.8 -s 5
```

**With accession-specific colours**

```bash
./MiSSN.sh -f example_dataset.fasta.gz -C example_colors.parquet -i 0.4 -c 0.8 -s 5
```

The builder can also run directly on existing alignment files; its metadata options are optional too:

```bash
python scripts/build_ssn.py -e cluster_alignments -i 0.4 -c 0.8 -o networks \
  --color-metadata-file example_colors.parquet
```

## Inputs

### Required arguments

- `-f <fasta_file>`: Path to the input protein sequence file in FASTA format (`.fasta`, `.fa`, `.faa`), supporting gzip-compressed inputs (`.gz`).

- `-i <min_seq_id>`: Minimum sequence identity threshold required to keep an edge between nodes in the network, specified as a decimal float from `0.0` to `1.0`.

- `-c <min_coverage>`: Minimum per-HSP alignment coverage required for both proteins to keep an edge between nodes in the network, specified as a decimal float from `0.0` to `1.0`. For example, `0.8` requires DIAMOND's `qcovhsp` and `scovhsp` to both be at least `80` percent. Regenerate older nine-column alignment TSVs with `scripts/all-vs-all.sh` to include these coverage fields.

- `-s <min_cluster_size>`: Minimum number of connected members required within a component to extract and output a distinct cluster.

### Optional arguments

Biome, Pfam, and colour metadata are independently optional; supply any combination, or omit all three.

- `-b <biome_metadata_file>`: Path to the metadata file containing biome annotations, provided in Apache Parquet (`.parquet`) format.

  **Example**
  | mgyp | biome |
  | :--- | :---  |
  | MGYP007687434724 | root:Host-associated:Birds:Digestive system |
  | MGYP003341907928 | root:Environmental:Aquatic:Marine:Oceanic |

- `-p <pfam_metadata_file>`: Path to the metadata file containing Pfam accessions, provided in Apache Parquet (`.parquet`) format.

  **Example**
  | mgyp | pfam_accession |
  | :--- | :---  |
  | MGYP000084156255 | PF22721 |
  | MGYP000117404121 | PF00756 |

- `-C <color_metadata_file>`: Parquet file containing string columns `accession` and `color`. Accessions must exactly match the sequence IDs in the FASTA/alignment files (the output node `id`). Colours may be `#RRGGBB` or `rgb(r,g,b)` with integer channels from 0–255; output colours are normalised to `rgb(r,g,b)`.

  | accession | color |
  | :--- | :--- |
  | MGYP007520866722 | #FF0000 |
  | MGYP001245521661 | rgb(0,128,255) |

  Matching accession colours take precedence over biome colours, including black for multiple biomes. Nodes without a matching override retain their biome colour, or grey if unavailable. Identical duplicate mappings are collapsed after colour normalisation; conflicting mappings, null/blank values, invalid colours, or missing/non-string required columns cause an error. Empty files with the correct string-column schema are accepted. Unmatched accessions do not create extra nodes.

- `-o <out_dir>`: Custom path for the output directory where network files and reports will be saved. Defaults to a dynamically generated name based on the input FASTA filename.

- `-g <color_group_parts>`: Number of leading levels of the biome lineage to use for grouping and color-coding nodes within the network visualization. (Default: `4`).

  **Example**
  
  If set to `4`, a biome of `root:Engineered:Solid waste:Composting:Bioreactor` is is grouped under `root:Engineered:Solid waste:Composting` to determine its base color.

- `-h`: Display the help message and exit the program. 

> [!IMPORTANT]
> The input dataset must use the `.fasta`, `.fa`, or `.faa` extension. Gzip-compressed versions (e.g., `.fasta.gz`) are also supported and will be decompressed on the fly. When supplied, biome, Pfam, and colour metadata must be valid `.parquet` files. Explicitly supplied missing or malformed files cause an error.

## Outputs

The filenames are prefixed with the sequence ID of the cluster's representative member.

- `<cluster_rep>_nodes.parquet`: Contains the node attributes for the specific cluster, including the sequence IDs, biome lineage, colours (accession overrides, biome colours, or grey), and Pfam annotations.

  **Example**
  | id | biome | pfam_accession | color |
  | :--- | :--- | :--- | :--- |
  | MGYP007520866722 | root:Environmental:Terrestrial:Soil:Contaminated | PF13604 | rgb(121,216,226) |
  | MGYP001245521661 | root:Engineered:Wastewater:Activated Sludge | PF22721 | rgb(118,221,75) |

> [!NOTE]
> If a protein sequence has multiple biomes or Pfam accessions, the values are concatenated into a single string separated by semicolons, like `A;B`. Missing biome and Pfam annotations are Parquet nulls, both when a file is omitted and when a sequence has no matching annotation. Missing Pfam values previously used the string `None`; they now use nulls. Nodes with multiple biomes are coloured black unless an accession colour override is provided. Nodes without a biome or override use `rgb(200, 200, 200)`.

Nodes are the distinct query and subject IDs present in each raw alignment file, including those whose edges fail filtering. Annotation and colour files do not determine node membership. FASTA sequences absent from every alignment row are not included.

- `<cluster_rep>_edges.parquet`: Contains the network edges between the sequences in that cluster, including sequence identity.

  **Example**
  | source | target | sequence identity |
  | :--- | :--- | :--- |
  | MGYP010172319147 | MGYP011046998155 | 77.2 |
  | MGYP000588563658 | MGYP010172319147 | 78.1 |

> [!TIP]
> These output files can be imported into the Cosmograph web app or Cytoscape. Updates to the provided `visualise_ssn.ipynb` notebook are deferred: its current annotation filters exclude null values and may hide unannotated nodes, and its biome filter swatches assume one colour per biome.

# Resources

Example sequence data and metadata annotations are available in the `/dataset` directory allowing you to evaluate the pipeline immediately after installation. All provided data is derived from the MGnify Protein Database release 2026_07.

# Credits

MiSSN was originally written as a [Google Summer of Code 2026 project](https://summerofcode.withgoogle.com/programs/2026/projects/SzccYhba) by [Szabolcs Vidám](https://github.com/vid-szabi) under the mentoring of [Christian Atallah](https://github.com/chrisAta) and [Ekaterina Sakharova](https://github.com/KateSakharova).

# Citation

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

# Licensing

All of the source code making up MiSSN in this repository is licensed under the terms of the [Apache 2.0 license](LICENSE). The example datasets provided are licensed under the terms of the [CC0 1.0 Universal (CC0 1.0) license](https://creativecommons.org/publicdomain/zero/1.0/).
