# Metagenomics Analysis Workflow

## Overview

This repository contains the bioinformatics and statistical workflows used to analyse root and rhizosphere metagenomes from New Zealand pasture systems.

## Bioinformatics Pipeline

### 1. Quality Control

Raw paired-end Illumina NovaSeq reads were processed using:
- fastp v0.23.4
- FastQC v0.12.1
- MultiQC v1.24.1

Quality filtering included:
- Adapter trimming
- Minimum quality score: Q20
- Minimum read length: 50 bp

### 2. Host Removal

Reads were aligned against:
- Lolium perenne (GCF_019359855.2)
- Trifolium repens (GCA_030408175.1)
- Homo sapiens (GRCh38)

Bowtie2 v2.5.4 was used with --very-sensitive.

### 3. Taxonomic Profiling

- Kraken2 v2.1.6
- PlusPF database
- Bracken v2.7

Taxonomic abundances generated at:
- Phylum
- Family
- Genus
- Species

### 4. Functional Profiling

- MetaPhlAn4
- HUMAnN4
- UniRef90
- MetaCyc

Outputs:
- Gene family abundance
- Reaction abundance
- Pathway abundance

### 5. Statistical Analyses

Taxonomy:
- Relative abundance
- Alpha diversity
- Bray-Curtis dissimilarity
- PCoA
- PERMANOVA
- ALDEx2

Function:
- Pathway abundance
- Bray-Curtis dissimilarity
- PCoA
- PERMANOVA
- ALDEx2
- Heatmaps
- Stratified pathway analyses

## Computational Environment

Platform: REANNZ HPC
Scheduler: SLURM
Partition: Milan
Resources: up to 32 CPUs and 180 GB RAM
