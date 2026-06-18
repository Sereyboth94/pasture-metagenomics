\# Functional Potential Workflow



\## Overview



This workflow performs functional profiling of pasture root and rhizosphere metagenomes using HUMAnN pathway abundance tables.



The workflow generates manuscript figures, supplementary figures, pathway consistency analyses, differential abundance analyses, and contributor analyses.



\---



\# Input Data



Root compartment:



```text

function/data/root/

├── Root\_metadata.tsv.gz

├── Root\_pathabundance\_UNSTRATIFIED.tsv.gz

├── Root\_pathabundance\_UNSTRATIFIED\_CPM.tsv.gz

└── Root\_pathabundance\_STRATIFIED.tsv.gz

```



Rhizosphere compartment:



```text

function/data/rhizosphere/

├── Rhizo\_metadata.tsv

├── Rhizosphere\_pathabundance\_CPM\_unstratified\_Rprefix.tsv.gz

├── Rhizosphere\_pathabundance\_RAW\_unstratified\_Rprefix.tsv.gz

└── Rhizosphere\_pathabundance\_CPM\_stratified\_Rprefix.tsv.gz

```



\---



\# Run Workflow



From repository root:



```bash

Rscript function/scripts/run\_function\_pipeline.R

```



Or from R:



```r

source("function/scripts/run\_function\_pipeline.R")

```



\---



\# Analyses Performed



\## Functional Community Structure



\* Relative abundance profiling

\* CPM normalization

\* Bray–Curtis dissimilarity



\## Ordination



\* Principal Coordinates Analysis (PCoA)



\## Statistical Testing



\* PERMANOVA

\* ALDEx2 differential abundance analysis



\## Pathway Analyses



\* Differential pathway abundance

\* Pathway consistency across locations

\* Inoculation response analyses

\* Stratified pathway contribution analyses



\---



\# Manuscript Figures



Generated automatically.



\## Figure 1C



Functional pathway PCoA.



\## Figure 3



Functional response to inoculation.



\## Figure 4



Pathway consistency across locations.



\## Figure 5



Stratified pathway contributors.



All manuscript figures are copied to:



```text

manuscript\_figures/Main\_Figures/

```



\---



\# Supplementary Figures



Saved to:



```text

manuscript\_figures/Supplementary\_Function\_Figures/

```



Includes:



\* Figure S3

\* Figure S4

\* Figure S5



\---



\# Outputs



Analysis outputs:



```text

function/results/

```



Publication-ready figures:



```text

manuscript\_figures/Main\_Figures/

```



Supplementary figures:



```text

manuscript\_figures/Supplementary\_Function\_Figures/

```



\---



\# Reproducibility



\* Repository-relative paths

\* Fixed random seed (`set.seed(1)`)

\* Standalone workflow

\* Publication-ready PDF, PNG and SVG outputs

\* Figure merge steps performed automatically



