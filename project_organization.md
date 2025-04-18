# Yeast MSA Project Organization

This document describes the organization of the Yeast MSA (Multiple Sequence Alignment) project, which analyzes genetic mutations in *Saccharomyces cerevisiae* under different experimental treatments.

## Overview

The project is organized to facilitate:
1. Reproducibility of the entire analysis pipeline
2. Clear separation of scripts and their corresponding outputs
3. Logical organization of different analysis types
4. Easy navigation between related components

## Directory Structure

### Scripts Directory (`/scripts/`)

The scripts directory is organized hierarchically, separating pipeline components from analysis modules:

#### Pipeline (`/scripts/pipeline/`)

Sequential scripts for the complete analysis workflow:

| Script                | Purpose                                                |
|-----------------------|--------------------------------------------------------|
| `01_preprocessing.sh` | Quality control and adapter trimming of raw reads      |
| `02_alignment.sh`     | Alignment of reads to the reference genome             |
| `03_variant_calling.sh` | Identifying genetic variants from aligned reads      |
| `04_filtering.sh`     | Quality filtering of variants                          |
| `05_annotation.sh`    | Functional annotation of filtered variants             |
| `06_comparisons.sh`   | Comparing variants across treatments                   |
| `07_run_full_pipeline.sh` | Master script to run the entire pipeline           |

#### Analysis Modules (`/scripts/analysis/`)

Each subdirectory contains scripts for a specific type of analysis:

1. **Mutation Spectrum** (`/scripts/analysis/mutation_spectrum/`)
   - Analysis of mutation types (transitions, transversions)
   - Mutation classification by nucleotide change

2. **Genomic Context** (`/scripts/analysis/genomic_context/`)
   - Extraction and analysis of sequence context around mutations
   - Motif identification and visualization

3. **Population Structure** (`/scripts/analysis/population_structure/`)
   - Analysis of relationships between samples
   - Clustering and dimensionality reduction

4. **Regional Enrichment** (`/scripts/analysis/regional_enrichment/`)
   - Identification of mutation hotspots
   - Statistical testing for enriched regions

5. **Scaffold Distribution** (`/scripts/analysis/scaffold_distribution/`)
   - Genome-wide distribution of mutations
   - Identification of enriched scaffolds

6. **Mutational Signatures** (`/scripts/analysis/mutational_signatures/`)
   - Extraction of characteristic mutation patterns
   - Signature comparison across treatments

7. **Statistical Patterns** (`/scripts/analysis/statistical_patterns/`)
   - Correlation and regression analysis
   - Multivariate statistical methods

#### Visualization (`/scripts/visualization/`)

Scripts for generating visualizations:
- Core plotting functions
- Advanced visualizations
- Publication-ready figures

#### Utilities (`/scripts/utils/`)

Common functions used across multiple scripts:
- Data processing
- File handling
- Statistical testing
- Sequence manipulation

#### Reporting (`/scripts/reporting/`)

Scripts for generating reports:
- R Markdown templates
- Report generation functions

### Results Directory (`/results/`)

The results directory mirrors the analysis workflow, with clear separation of intermediate and final outputs:

#### Preprocessing (`/results/preprocessing/`)

Quality control outputs:
- FastQC reports
- Trimming statistics
- Quality metrics

#### Alignment (`/results/alignment/`)

Genome alignment results:
- BAM files
- Alignment metrics
- Coverage statistics

#### Variants (`/results/variants/`)

Variant calling results:
- Raw and filtered variants
- Merged variant sets
- Cross-sample comparisons

#### Annotation (`/results/annotation/`)

Functional annotation:
- Annotated VCF files
- Gene-level summaries
- Functional impact predictions

#### Analysis (`/results/analysis/`)

Analysis results organized by category:

1. **Mutation Spectrum** (`/results/analysis/mutation_spectrum/`)
   - Mutation type distributions
   - Transition/transversion ratios
   - Treatment-specific spectra

2. **Genomic Context** (`/results/analysis/genomic_context/`)
   - Sequence logos
   - Context matrices
   - Context analysis by mutation type and treatment

3. **Population Structure** (`/results/analysis/population_structure/`)
   - Clustering results
   - PCA plots
   - Sample similarity metrics
   - Variant sharing patterns

4. **Regional Enrichment** (`/results/analysis/regional_enrichment/`)
   - Hotspot identification
   - Enrichment statistics
   - Treatment-specific and shared regions

5. **Scaffold Distribution** (`/results/analysis/scaffold_distribution/`)
   - Density maps
   - Enriched scaffolds
   - Correlation analysis

6. **Mutational Signatures** (`/results/analysis/mutational_signatures/`)
   - Signature profiles
   - Context enrichment
   - Signature similarities

7. **Statistical Patterns** (`/results/analysis/statistical_patterns/`)
   - Correlation results
   - Regression models
   - Clustering and PCA results

#### Figures (`/results/figures/`)

Publication-ready visualizations:
- Main manuscript figures
- Supplementary materials
- Presentation slides

#### Reports (`/results/reports/`)

Analysis reports:
- Project summaries
- Technical documentation
- Treatment-specific reports
- Integrated analyses

## Workflow

The recommended workflow is:

1. Run the pipeline scripts sequentially or use the master script:
   ```bash
   bash scripts/pipeline/07_run_full_pipeline.sh
   ```

2. Run specific analysis modules as needed:
   ```bash
   python scripts/analysis/mutation_spectrum/spectrum_analysis.py
   Rscript scripts/analysis/genomic_context/motif_analysis.R
   ```

3. Generate visualizations:
   ```bash
   Rscript scripts/visualization/figure_generation.R
   ```

4. Produce reports:
   ```bash
   Rscript scripts/reporting/report_generation.R
   ```

## Script-Output Relationships

Each script is designed to output to a specific results directory:

| Script Pattern | Output Directory |
|----------------|------------------|
| `scripts/pipeline/*.sh` | `results/preprocessing/`, `results/alignment/`, `results/variants/` |
| `scripts/analysis/mutation_spectrum/*` | `results/analysis/mutation_spectrum/` |
| `scripts/analysis/genomic_context/*` | `results/analysis/genomic_context/` |
| (and so on for each analysis type) | |

This organization ensures that the relationship between scripts and their outputs is clear and consistent.

## Dependencies

Scripts indicate their dependencies at the top of each file. Generally, the project requires:

- Bash shell environment
- Python 3.6+ with BioPython, NumPy, Pandas
- R 4.0+ with tidyverse, ggplot2, and specific bioinformatics packages
- Bioinformatics tools: BWA, Samtools, BCFtools, FastQC, fastp

## Configuration

Global configuration parameters are stored in `scripts/utils/config.py` and can be modified to adjust analysis parameters.

## Versioning

Each major analysis has a version identifier in its output directory to track changes over time.