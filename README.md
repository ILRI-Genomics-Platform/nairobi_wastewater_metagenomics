# AMR Wastewater Surveillance: Socioeconomic Stratification Analysis

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
- DESeq2 differential abundance testing of individual AMR/drug classes
- Consolidation and summary of statistical outputs across all analyses

## Repository structure

```
.
├── scripts/
│   ├── functions.R                                  # shared helper functions (package install/load, etc.)
│   ├── ordination_permanova_functions.R              # PCoA/NMDS/envfit/PERMANOVA helper functions
│   ├── 06_ARG_ResFinder_Analysis_dimentionalityReduction.R
│   ├── 06_ARG_ResFinder_Analysis_perDrugClass.R
│   ├── 07_ARG_Ordination_PERMANOVA_AllMetrics.R      # main ordination + PERMANOVA driver
│   ├── 08_ARG_DESeq2_DifferentialAbundance.R         # differential abundance testing
│   ├── 09_parse_PERMANOVA_pairwiseAdonis_outputs.R   # consolidates all statistical outputs
│   └── amr_heatmap_venn.R                            # heatmap + Venn diagram figure generation
├── results_arg/
│   └── r-analysis/
│       ├── data/       # input abundance matrices, metadata, consolidated summary tables
│       └── plots/      # generated figures and per-model statistical output files
├── docs/
│   └── ANALYSIS_METHODS_AND_INTERPRETATION.md        # methods rationale, citations, interpretation guide
├── LICENSE
└── README.md
```

## Requirements

- R >= 4.4
- R packages: `dplyr`, `tidyr`, `tibble`, `ggplot2`, `stringr`, `vegan`,
  `ggrepel`, `pheatmap`, `ggvenn`, `patchwork`, `ggplotify`, `DESeq2`, `ashr`,
  `pairwiseAdonis` (GitHub: `pmartinezarbizu/pairwiseAdonis`), `purrr`, `readr`

Bioconductor packages (`DESeq2`) require R >= 4.4 with a matching Bioconductor
release (3.20). See `docs/ANALYSIS_METHODS_AND_INTERPRETATION.md` for full
package version notes and rationale for methodological choices.

## Usage

Scripts are numbered in the order they are intended to be run. Each script
sources shared functions from `scripts/functions.R` and, where relevant,
`scripts/ordination_permanova_functions.R`. Update the `setwd()` and file path
variables at the top of each script to match your local environment before
running.

```r
source("scripts/07_ARG_Ordination_PERMANOVA_AllMetrics.R")
source("scripts/08_ARG_DESeq2_DifferentialAbundance.R")
source("scripts/09_parse_PERMANOVA_pairwiseAdonis_outputs.R")
```

## Data availability

Raw sequencing data availability and accession numbers are described in the
manuscript. Intermediate abundance matrices and metadata required to
reproduce the statistical analyses in this repository are provided in
`results_arg/r-analysis/data/` [update with accession/Zenodo DOI once
deposited]. The custom Centrifuge database is also available on Zenodo:
[update with accession/Zenodo DOI once deposited]

## Citation

If you use this code, please cite:

> [Author list]. [Manuscript title]. [Journal, year]. DOI: [add once available]

See `docs/ANALYSIS_METHODS_AND_INTERPRETATION.md` for citations of the
statistical methods and software packages used.

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE)
file for details. The MIT license permits reuse, modification, and
redistribution of this code (including for commercial purposes) provided the
original copyright notice is retained; it comes with no warranty. If your
institution or funder requires a different license for published research
code (e.g. GPL-3.0 or Apache-2.0), swap the `LICENSE` file accordingly before
publishing the repository — this is a one-file change and does not require
editing anything else here.

## Contact

[Your name / lab] — [email or ORCID]
