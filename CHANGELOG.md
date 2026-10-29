# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

## [1.1.0] - 2026-10-29

### Added

- Metanetwork generation through the pipeline's `--metanetwork` option. This option 
  combines all node parquet files into one, and all edge parquet files into one. Use
  this option if feasible to get an overview of the sequence similarity patterns across
  your entire input dataset.
- Accession-specific node colour overrides through `-C` / `--color-metadata-file`,
  accepting `#RRGGBB` and `rgb(r,g,b)` values and taking precedence over biome colours.
  This option is especially useful when an _a priori_ separation of the nodes exists,
  e.g. annotated vs un-annotated sequences, and the user wants to visualise the network
  with this separation in mind.

### Changed

- Biome and Pfam metadata files are now optional, allowing networks to
  be generated with either, both, or neither annotation file.
- Missing biome and Pfam annotations are now represented as Parquet 
  `NULL` values in generated networks.
- Alignment coverage filtering now directly uses DIAMOND's query and subject coverage 
  (`qcovhsp` and `scovhsp`) values for comparing against minimum coverage threshold.
  Specifically, the smallest value between the two is chosen and then compared to the
  set coverage threshold. This change replaces the calculation based on positive matches 
  and sequence lengths.

### Fixed

- Retain nodes present in alignment files even when they have no biome
  annotations. Nodes without a biome colour or an accession-specific override
  default to grey.
- Preserve numeric-looking sequence identifiers as strings, including leading
  zeroes, so annotation and colour mappings match the original identifiers.

## [1.0.0] - 2026-08-17

### Added

- Initial release of MiSSN, a workflow for generating sequence
  similarity networks (SSNs) from sets of protein sequences. MiSSN pre-clusters
  sequences with MMseqs2, calculates within-cluster all-vs-all alignments with
  DIAMOND/blastp, and builds network files in Parquet format, usable by [Cosmograph](https://cosmograph.app/).
- Network annotations with GOLD biome classifications and Pfam accessions,
  plus Jupyter notebooks for exploring networks and cluster size distributions.

[Unreleased]: https://github.com/EBI-Metagenomics/MiSSN/compare/v1.1.0...HEAD
[1.1.0]: https://github.com/EBI-Metagenomics/MiSSN/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/EBI-Metagenomics/MiSSN/releases/tag/v1.0.0
