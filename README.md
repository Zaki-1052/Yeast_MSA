# Yeast Genomic Analysis Project

## Overview

This repository contains a comprehensive bioinformatics pipeline and analysis of genetic mutations in *Saccharomyces cerevisiae* (baker's yeast) under four experimental treatments (WT, WTA, STC, and CAS). The project involves whole-genome sequencing, variant calling, and extensive downstream analyses to characterize mutation patterns and their biological implications. A reference-guided multiple sequence alignment was performed where the whole genomic sequence of the yeast cells were aligned with the reference genome from [here](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000766475.2/)

## Key Findings

1. **Treatment-Specific Mutation Signatures**: Each treatment produces a unique "genetic signature" with characteristic mutation patterns
2. **Clustered Mutations**: 71-84% of mutations occur in clusters rather than being randomly distributed
3. **Hotspot Regions**: Specific genomic regions show higher susceptibility to mutations across treatments
4. **Mutation Type Patterns**: 
   - CAS treatment produces exclusively transition mutations (100%)
   - WT and WTA treatments show complementary mutation patterns (G→A vs. A→G)
   - STC treatment displays the most diverse mutation spectrum

## Directory Structure

### Data

- `data/`: Contains raw and processed sequencing data (note: raw data files not included in repository)

### Reference

- `reference/`: Reference genome files
  - `yeast_w303.fasta`: Yeast W303 strain reference genome
  - Associated index files (`.amb`, `.ann`, `.bwt`, `.fai`, `.pac`, `.sa`)
  - `sequence_report.tsv`: Reference sequence metadata

### Scripts

- `scripts/`: Analysis pipeline scripts
  - `preprocessing.sh`: Quality control and adapter removal
  - `alignment.sh`: Maps reads to reference genome
  - `variants.sh`: Calls genetic variants
  - `comparisons.sh`: Compares variants across treatments
  - `consistency.sh`: Analyzes variant reproducibility across replicates
  - `hotspot_analysis.sh`: Identifies mutation hotspots
  - `position_analysis.sh`: Analyzes mutation positions
  - `sequence_context_analysis.R`: Examines sequence context around mutations
  - `*_visualization.R`: Various visualization scripts

### Results

The `results/` directory contains all analysis outputs organized by analysis type:

#### Quality Control

- `results/fastqc/`: Quality reports for raw sequencing data
- `results/fastp/`: Trimming and filtering reports
- `results/stats/`: Alignment statistics and quality metrics
  - `alignment_summary.txt`: Coverage and mapping statistics

#### Variant Calling

- `results/vcf/`: Variant call format files
  - `individual/`: Variants for each sample
  - `filtered/`: Quality-filtered variants
  - `comparison/`: Variant comparisons between treatments
  - `consistency/`: Replicate consistency analysis
  - `direct_comparison/`: Control vs treatment comparisons
  - `merged/`: Combined variant sets

#### Variant Annotation

- `results/annotation/`: Functional annotations of variants
  - `*_variants.tsv`: Annotated variant tables
  - `*_scaffold_counts.txt`: Variant distribution by scaffold

#### Functional Analysis

- `results/functional/`: Functional impact analysis
  - Treatment-specific variant details
  - Cross-treatment variant comparisons
  - Overlap analysis between treatments

#### Sequence Context Analysis

- `results/sequence_context/`: Analysis of sequence context around mutations
  - `*_contexts.txt`: Extracted sequence contexts
  - `*_logo.png`: Sequence motif logos
  - `*_frequencies.tsv`: Position-specific nucleotide frequencies
  - `*_motifs.tsv`: Enriched sequence motifs
  - `sequence_context_report.html`: Comprehensive sequence context analysis

#### Hotspot Analysis

- `results/analysis/hotspots/`: Mutation clustering analysis
  - `*_density.tsv`: Mutation density calculations
  - `*_top_hotspots.txt`: Identified mutation hotspots
  - `hotspot_summary.txt`: Summary of hotspot findings
  - `plots/`: Visualizations of hotspot distributions
  - `position_analysis/`: Detailed positional analysis

#### Position Clustering

- `results/analysis/position_clustering/`: Analysis of mutation proximity
  - `*_distances.txt`: Distance calculations between mutations
  - `*_position_density.pdf`: Density plot visualizations
  - `*_stats.txt`: Statistical summaries of positions
  - `summary.txt`: Overall clustering findings

#### Mutation Signatures

- `results/analysis/signatures/`: Mutation pattern analysis
  - `*_specific_variants.tsv`: Treatment-specific variant collections
  - `mutation_patterns.txt`: Identified mutation signatures

#### Visualizations

- `results/visualization/`: High-level data visualizations
  - `variant_distribution.png`: Mutation distribution plots
  - `scaffold_distribution.png`: Variant location across genome
  - `variant_sharing_heatmap.png`: Similarity between treatments
  - `indel_percentage.png`: Insertion/deletion statistics

- `results/visualizations/`: Advanced multi-dimensional visualizations
  - `signatures/`: Mutation spectrum visualizations
  - `clustering/`: Mutation clustering visualizations
  - `integrative/`: Combined data visualizations
  - `interpretation/`: Biological mechanism visualizations

#### Summary Reports

- `results/yeast_analysis_summary.Rmd`: R Markdown summary report
- `results/yeast_analysis_summary.html`: HTML rendered report
- `results/summary_document.R`: R script for report generation

## Workflow

The analysis pipeline follows these main steps:

1. **Quality Control and Preprocessing**
   - FastQC for raw read quality assessment
   - fastp for adapter removal and quality trimming

2. **Alignment**
   - BWA-MEM for mapping reads to reference genome
   - Samtools for BAM processing and quality filtering

3. **Variant Calling**
   - Variant identification using samtools/bcftools
   - Filtering for high-confidence variants
   - Comparison against control samples

4. **Variant Analysis**
   - Mutation type classification
   - Sequence context extraction
   - Positional clustering analysis
   - Hotspot identification

5. **Visualization and Interpretation**
   - Mutation spectrum plots
   - Position clustering visualization
   - Sequence motif analysis
   - Hotspot characterization
   - Treatment comparisons and relationships

## Key Visualizations

- **Mutation Spectra**: `results/visualizations/signatures/mutation_spectrum.png`
- **Transition vs Transversion**: `results/visualizations/signatures/transition_transversion.png`
- **Mutation Clustering**: `results/visualizations/clustering/position_distribution.png`
- **Hotspot Positions**: `results/analysis/hotspots/plots/hotspot_positions.pdf`
- **Sequence Motifs**: `results/sequence_context/motif_comparison.png`
- **Treatment Relationships**: `results/visualizations/integrative/treatment_network.png`

## Summary of Treatment Effects

### CAS Treatment
- 100% transition mutations
- Strong C→T bias (40% of mutations)
- High clustering (83.6%)
- Likely affects cytosine deamination repair pathways

### WT Treatment
- Strong G→A bias (45.5%)
- High clustering (83.1%)
- Associated with oxidative damage to guanine bases

### WTA Treatment
- Strong A→G bias (55.6%)
- Complementary to WT's mutation pattern
- Similar transition percentage to WT (77.8%)

### STC Treatment
- Most diverse mutation profile
- Lowest clustering (70.6%)
- Highest transversion rate (40%)
- Suggests disruption of multiple DNA repair pathways

## Comprehensive Summary Report

For a complete analysis with interpretations and biological implications, see:
- `results/yeast_analysis_summary.html`

## Usage

The scripts directory contains all necessary code to reproduce the analysis. Most scripts are designed to be run sequentially:

```bash
# Example workflow
bash scripts/preprocessing.sh
bash scripts/alignment.sh
bash scripts/variants.sh
bash scripts/comparisons.sh
bash scripts/hotspot_analysis.sh
Rscript scripts/sequence_context_analysis.R
Rscript scripts/visualization.R
Rscript scripts/summary_document.R
```

## Requirements

- BWA
- Samtools/BCFtools
- FastQC
- fastp
- R with packages:
  - dplyr
  - ggplot2
  - knitr
  - kableExtra
  - seqLogo
