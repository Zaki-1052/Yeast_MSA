# Yeast Genomic Analysis Project

See [project organization documentation](project_organization.md) and [additional documentation](documentation.md) for more details.

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

The project follows a modular organization:

### Scripts

Scripts are organized by their function:

- `scripts/pipeline/`: Sequential scripts for the complete analysis workflow
- `scripts/analysis/`: Modules for specific types of analyses
- `scripts/visualization/`: Scripts for generating visualizations
- `scripts/utils/`: Common utility functions
- `scripts/reporting/`: Report generation scripts

### Results

Results are organized to mirror the analysis workflow:

- `results/preprocessing/`: Quality control outputs
- `results/alignment/`: Genome alignment results
- `results/variants/`: Variant calling results
- `results/annotation/`: Functional annotations
- `results/analysis/`: Analysis results by category
  - `mutation_spectrum/`: Mutation type analysis
  - `genomic_context/`: Sequence context analysis
  - `population_structure/`: Population relationships
  - `regional_enrichment/`: Hotspot identification
  - `scaffold_distribution/`: Genome-wide distribution
  - `mutational_signatures/`: Signature analysis
  - `statistical_patterns/`: Statistical analyses
- `results/figures/`: Publication-ready visualizations
- `results/reports/`: Analysis reports

See the [project organization document](project_organization.md) for a complete description of the directory structure.

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

## Usage

The scripts directory contains all necessary code to reproduce the analysis. The main pipeline can be run with:

```bash
# Run the entire pipeline
bash scripts/pipeline/07_run_full_pipeline.sh

# Or run individual analysis modules
python scripts/analysis/mutation_spectrum/spectrum_analysis.py
```

## Requirements

- BWA
- Samtools/BCFtools
- FastQC
- fastp
- Python 3.6+ with:
  - pandas
  - numpy
  - matplotlib
  - seaborn
  - biopython
- R with packages:
  - dplyr
  - ggplot2
  - pheatmap
  - RColorBrewer