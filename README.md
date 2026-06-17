# Pasture Rhizosphere and Root Metagenomics Study

## Overview

This repository contains the bioinformatics workflows, processed metagenomic data, metadata, statistical analyses, and figure-generation scripts used to investigate microbial taxonomic composition and functional potential in pasture rhizosphere and root compartments across multiple field locations in New Zealand.

The repository accompanies a research manuscript and has been structured to support full computational reproducibility, GitHub publication, and Zenodo archiving.

---

# Repository Structure

```text
Github/
├── bioinformatics/
├── location/
├── taxonomy/
├── function/
├── manuscript_figures/
├── docs/
└── README.md
```

---

# Quick Start

## Requirements

- R >= 4.4.0
- Kraken2
- Bracken
- HUMAnN

### Install Required R Packages

```r
install.packages(c(
  "tidyverse",
  "vegan",
  "ggplot2",
  "ape",
  "pheatmap",
  "patchwork",
  "cowplot",
  "ggrepel",
  "sf",
  "rnaturalearth",
  "rnaturalearthdata",
  "ggspatial",
  "ggsci",
  "digest",
  "here"
))
```

---

# Running the Workflows

Open R or RStudio and set the repository root:

```r
setwd("path/to/Github")
```

## 1. Generate Figure 1

```r
source("location/run_location_figure1_pipeline.R")
```

Generates:

- Figure 1A (Study-site map)
- Merged Figure 1 manuscript panel

Outputs:

```text
manuscript_figures/Main_Figures/
```

---

## 2. Generate Taxonomy Figures

```r
source("taxonomy/scripts/run_taxonomy_pipeline.R")
```

Generates:

- Figure 1B (Taxonomy PCoA)
- Figure 2 (Taxonomic response)
- Supplementary taxonomy figures

Outputs:

```text
taxonomy/results/
manuscript_figures/Main_Figures/
manuscript_figures/Supplementary_Taxonomy_Figures/
```

---

## 3. Generate Functional Potential Figures

```r
source("function/scripts/run_function_pipeline.R")
```

Generates:

- Figure 1C (Pathway PCoA)
- Figure 3 (Functional response)
- Figure 4 (Pathway consistency)
- Figure 5 (Stratified pathway contributors)
- Supplementary Figures S3–S5

Outputs:

```text
function/results/
manuscript_figures/Main_Figures/
manuscript_figures/Supplementary_Function_Figures/
```

---

# Recommended Workflow Order

For complete manuscript reproduction:

```r
source("location/run_location_figure1_pipeline.R")

source("taxonomy/scripts/run_taxonomy_pipeline.R")

source("function/scripts/run_function_pipeline.R")
```

Final publication-ready figures will be available in:

```text
manuscript_figures/Main_Figures/
```

---

# Figure Reproducibility

| Figure | Workflow |
|----------|----------|
| Figure 1A | location |
| Figure 1B | taxonomy |
| Figure 1C | function |
| Figure 2 | taxonomy |
| Figure 3 | function |
| Figure 4 | function |
| Figure 5 | function |
| Supplementary Taxonomy Figures | taxonomy |
| Supplementary Function Figures | function |

---

# Reproducibility

- Repository-relative paths
- Fixed random seed (`set.seed(1)`)
- Publication-quality PDF, PNG and SVG outputs
- Workflow documentation provided in `docs/`

---

# Raw Sequencing Data

BioProject: PRJNA1475531

SRA accession(s): PRJNA1475531

---

# Code Availability

Permanent Zenodo archive:

**DOI:** TO BE ADDED

---

# Contact

**Sereyboth Soth**

Lincoln University, New Zealand

GitHub: https://github.com/Sereyboth94
