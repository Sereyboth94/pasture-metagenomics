\# Taxonomy Workflow



\## Overview



This workflow performs taxonomic profiling analyses of wheat root and rhizosphere metagenomes using Bracken abundance tables derived from Kraken2 classifications.



The workflow generates manuscript figures, supplementary figures, statistical analyses, and publication-ready outputs.



\---



\# Input Data



Root compartment:



```text

taxonomy/data/root/

├── Metadata.txt

├── Bracken\_Phylum.tsv

├── Bracken\_Family.tsv

└── Bracken\_Genus.tsv

```



Rhizosphere compartment:



```text

taxonomy/data/rhizosphere/

├── Metadata.txt

├── Bracken\_Phylum.tsv

├── Bracken\_Family.tsv

└── Bracken\_Genus.tsv

```



\---



\# Run Workflow



From the repository root:



```bash

Rscript taxonomy/scripts/run\_taxonomy\_pipeline.R

```



Or from R:



```r

source("taxonomy/scripts/run\_taxonomy\_pipeline.R")

```



\---



\# Analyses Performed



\## Community Structure



\* Relative abundance calculation

\* Square-root transformation

\* Bray–Curtis dissimilarity



\## Ordination



\* Principal Coordinates Analysis (PCoA)

\* Community clustering assessment



\## Statistical Testing



\* PERMANOVA

\* Site effects

\* Inoculation effects

\* Compartment effects



\## Taxonomic Visualization



\* Phylum-level composition

\* Family-level composition

\* Genus-level composition

\* Rare taxa summarization



\---



\# Manuscript Figures



Generated automatically:



\## Figure 1B



Taxonomic PCoA



Output:



```text

manuscript\_figures/Main\_Figures/

```



\## Figure 2



Taxonomic response across locations and treatments.



Output:



```text

manuscript\_figures/Main\_Figures/

```



\---



\# Supplementary Figures



Saved to:



```text

manuscript\_figures/Supplementary\_Taxonomy\_Figures/

```



Includes:



\* Supplementary ordinations

\* Supplementary abundance plots

\* Additional statistical summaries



\---



\# Outputs



Main results:



```text

taxonomy/results/

```



Publication-ready figures:



```text

manuscript\_figures/Main\_Figures/

```



Supplementary figures:



```text

manuscript\_figures/Supplementary\_Taxonomy\_Figures/

```



\---



\# Reproducibility



\* Repository-relative paths

\* Fixed random seed (`set.seed(1)`)

\* Fully standalone workflow

\* Publication-ready PDF, PNG and SVG outputs



