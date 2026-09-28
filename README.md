# Pasture Rhizosphere and Root Metagenomics Study

## Overview

This repository contains the bioinformatics workflows, processed metagenomic data, metadata, plant biomass measurements, statistical analyses, and figure-generation scripts used to investigate *Trichoderma* inoculation in pasture rhizosphere and root-associated microbiomes at five New Zealand locations. The Kowhai irrigated and rainfed fields are analysed separately, giving six field environments.

The repository accompanies a research manuscript and has been structured to support full computational reproducibility, GitHub publication, and Zenodo archiving.

---

# Repository Structure

```text
Github/
├── bioinformatics/
├── Biomass/
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
  "here",
  "readxl",
  "lme4",
  "lmerTest",
  "emmeans",
  "car"
))
```

---

# Running the Workflows

Open R or RStudio and set the repository root:

```r
setwd("path/to/Github")
```

## 1. Generate Taxonomy Figures

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

## 2. Generate Functional Potential Figures

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

## 3. Analyse Biomass and Assemble Figure 1

```r
source("Biomass/run_location_figure1_pipeline.R")
```

Generates:

- Figure 1A (Plant biomass)
- Merged Figure 1A–C manuscript figure
- Biomass analysis tables and diagnostic outputs

`Biomass/Plant_Biomass.xlsx` contains plot-level aboveground biomass measurements. Eight blocks were assessed in each field environment, with a control and an inoculated plot per block. The analysis reports treatment means and standard deviations, paired comparisons, and percentage change relative to the control mean. To run only the biomass analysis, use `source("Biomass/Plant_Biomass_ANOVA.R")`.

Outputs:

```text
Biomass/Plant_Biomass_ANOVA_results/
location/results/
manuscript_figures/Main_Figures/
```

---

# Recommended Workflow Order

For complete manuscript reproduction:

```r
source("taxonomy/scripts/run_taxonomy_pipeline.R")
source("function/scripts/run_function_pipeline.R")
source("Biomass/run_location_figure1_pipeline.R")
```

Final publication-ready figures will be available in:

```text
manuscript_figures/Main_Figures/
```

---

# Figure Reproducibility

| Figure | Workflow |
|----------|----------|
| Figure 1A | Biomass |
| Figure 1B | taxonomy |
| Figure 1C | function |
| Figure 1A–C assembly | Biomass |
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

SRA runs: see BioProject PRJNA1475531

---

# Code Availability

Permanent Zenodo archive:

**DOI:** TO BE ADDED

---

# Contact

**Sereyboth Soth**

Lincoln University, New Zealand

GitHub: https://github.com/Sereyboth94
