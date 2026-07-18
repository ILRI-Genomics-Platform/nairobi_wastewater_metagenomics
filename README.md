# Wastewater Genomic Surveillance: Socioeconomic Stratification Analysis
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21429478.svg)](https://doi.org/10.5281/zenodo.21429478)


Analysis code accompanying the manuscript investigating Bacterial abundance and diversity, antimicrobial resistance
(AMR) gene composition in wastewater across low-, middle-, and high-income catchment areas in Nairobi, Kenya, including a composite central
wastewater treatment plant (WWTP) site.

## Overview

This repository contains scripts for:

- Quality control and taxonomic/AMR classification of shotgun metagenomic
  sequencing reads (fastp, Bowtie2, Centrifuge, ResFinder)
- Bacterial species abundance, Diversity and PCoA Ordination clustering
- AMR abundance heatmap visualization with hierarchical clustering
- Gene overlap analysis (Venn diagrams) across socioeconomic categories
- Beta-diversity ordination (PCoA, NMDS) with Bray-Curtis and Robust Aitchison
  distances
- PERMANOVA and pairwise PERMANOVA testing of community composition
  differences across socioeconomic categories
- Consolidation and summary of statistical outputs across all analyses

## Repository structure

```
.
├── scripts/
├── docs/
├── LICENSE
└── README.md
```

## Requirements

- R >= 4.4
- R packages: `dplyr`, `tidyr`, `tibble`, `ggplot2`, `stringr`, `vegan`,
  `ggrepel`, `pheatmap`, `ggvenn`, `patchwork`, `ggplotify`, `DESeq2`, `ashr`,
  `pairwiseAdonis` (GitHub: `pmartinezarbizu/pairwiseAdonis`), `purrr`, `readr`

See `docs/packages.md` for full package version notes and rationale for methodological choices.

## Usage

Scripts are numbered in the order they are intended to be run. Each script
sources shared functions from `scripts/functions.R` and, where relevant,
`scripts/ordination_permanova_functions.R`. Update the `setwd()` and file path
variables at the top of each script to match your local environment before
running.

```r
source("scripts/01.sample_data.R")
```

## Data availability

Raw sequencing data availability and accession numbers are described in the manuscript.   
Intermediate abundance matrices and metadata required to reproduce the statistical analyses are available upon request. Also the custom Centrifuge database is also available on Zenodo: [DOI:10.5281/zenodo.21391161](https://doi.org/10.5281/zenodo.21391161)

## Citation

If you use this code, please cite:

> [Author list]. [Manuscript title]. [Journal, year]. DOI: [add once available]

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE)
file for details. The MIT license permits reuse, modification, and
redistribution of this code (including for commercial purposes) provided the
original copyright notice is retained; it comes with no warranty. If your
institution or funder requires a different license for published research
code (e.g. GPL-3.0 or Apache-2.0), swap the `LICENSE` file accordingly before
publishing the repository — this is a one-file change and does not require
editing anything else here.
