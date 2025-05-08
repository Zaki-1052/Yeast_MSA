# Yeast MSA (Multiple Sequence Alignment) Project

## Table of Contents
1. [Project Overview](#project-overview)
2. [Experimental Design](#experimental-design)
3. [Project Structure](#project-structure)
4. [Pipeline Overview](#pipeline-overview)
5. [Analysis Modules](#analysis-modules)
6. [Key Findings](#key-findings)
7. [Future Work: Annotation with snpEff](#future-work-annotation-with-snpeff)
8. [Previous Annotation Analysis](#previous-annotation-analysis)
9. [Usage](#usage)
10. [References](#references)

## Project Overview

This Yeast reference-guided multiple sequence analysis project investigates how yeast (S. cerevisiae, W303 strain) adapts to different environmental stresses through genetic mutations. The project employs comparative genomic analysis to study variations between different yeast strains subjected to different treatments (temperature, low oxygen) and gene modifications. By analyzing the mutational patterns, genomic contexts, and regional enrichments, this project aims to understand the molecular mechanisms underlying adaptation to environmental stress.

The primary goals of this project include:
1. Identifying and characterizing mutations in different treatment conditions
2. Analyzing mutation spectrum patterns and sequence contexts
3. Detecting regional enrichment of mutations across the genome
4. Comparing population structures between different treatment groups
5. Understanding the influence of gene modifications on adaptation mechanisms
6. Preparing for functional annotation of variants using snpEff

## Experimental Design

### Treatment Conditions

The project analyzes four main treatment conditions:

1. **WT-37**: Temperature-adapted wild type yeast
   - Adaptation: Temperature
   - Gene Modification: None

2. **WTA**: Low oxygen-adapted wild type yeast
   - Adaptation: Low Oxygen
   - Gene Modification: None

3. **STC**: STC gene-modified strain with low oxygen adaptation
   - Adaptation: Low Oxygen
   - Gene Modification: STC gene

4. **CAS**: CAS gene-modified strain with temperature adaptation
   - Adaptation: Temperature
   - Gene Modification: CAS gene

Each treatment has control samples (WT-CTRL, STC-CTRL, CAS-CTRL) and three biological replicates (e.g., WT-37-55-1, WT-37-55-2, WT-37-55-3).

### Experimental Groups

The experimental design allows for multiple comparison groups:

- **Treatment vs Control**: Comparing each treatment with its respective control
- **Adaptation Type**: Comparing temperature adaptation vs low oxygen adaptation
- **Gene Modification**: Comparing gene-modified vs non-modified strains
- **Treatment Pairs**: Direct comparisons between specific treatment pairs

## Project Structure

The project is organized into the following main directories:

```
Yeast_MSA/
├── analysis/                   # Analysis results
│   ├── genomic_context_results/    # Sequence context analysis
│   ├── mutation_spectrum_results/  # Mutation spectrum analysis
│   ├── mutational_signatures_results/ # Mutational signatures
│   ├── population_structure_results/ # Population structure analysis
│   ├── regional_enrichment_results/  # Regional mutation hotspots
│   ├── scaffold_distribution_results/ # Variant distribution by scaffold
│   └── statistical_pattern_results/   # Statistical analyses
├── annotation/                # Files related to annotation
│   └── vcf_preparation_report.txt   # Report on VCF preparation
├── mutation_spectrum_analysis/ # Mutation data for spectrum analysis
├── old_version/              # Previous versions (for reference)
├── reference/                # Reference genome and annotations
│   ├── chromosome_mapping.txt       # Chromosome naming maps
│   ├── w303_chromosomal.fasta      # Main reference genome
│   └── w303_annotations/           # Genome annotations
├── results/                  # Pipeline processing results
│   ├── merged/                     # Merged variant data
│   ├── preprocessing/              # Preprocessing results
│   └── stats/                      # Alignment statistics
├── scripts/                  # Analysis and pipeline scripts
│   ├── analysis/                   # Analysis scripts for each module
│   ├── pipeline/                   # Pipeline execution scripts
│   ├── utils/                      # Utility scripts
│   └── vcf/                        # VCF processing scripts
└── vcf/                      # VCF files
    ├── filtered/                   # Filtered VCF files
    ├── individual/                 # Individual sample VCFs
    ├── initial/                    # Initial VCF files
    └── merged/                     # Merged VCF files
```

## Pipeline Overview

The analysis pipeline follows these sequential steps:

1. **Setup** (01_setup.sh): Creates the directory structure for results

2. **Filtering & Normalization** (02_filter_normalize.sh): 
   - Processes VCF files with quality filters (QUAL≥20, DP≥10)
   - Normalizes variant representation

3. **Contig Fixing** (03_fix_contigs.sh): 
   - Ensures consistency in contig definitions
   - Standardizes chromosome naming

4. **Merging VCFs** (04_merge_vcfs.sh): 
   - Combines VCF files from different samples
   - Creates merged datasets for analysis

5. **Sample Testing** (05_test_samples.sh): 
   - Validates sample selection
   - Performs quality control checks

6. **Group Comparisons** (06_compare_groups.sh): 
   - Compares treatment vs control samples
   - Identifies treatment-specific variants

7. **Cross-treatment Comparison** (07_cross_treatment.sh): 
   - Analyzes differences between treatments
   - Identifies shared and unique variants

8. **Unique Variant Identification** (08_identify_unique.sh): 
   - Finds treatment-specific variants
   - Identifies adaptation-specific mutations

9. **Direct Comparisons** (09_direct_comparison.sh): 
   - Performs pairwise VCF comparisons
   - Calculates similarity metrics

10. **Consistency Analysis** (10_consistency_analysis.sh): 
    - Evaluates variation within replicates
    - Measures reproducibility

11. **Summary Report Generation** (11_summary_report.sh): 
    - Creates final report
    - Summarizes key statistics and findings

The complete pipeline can be executed using the `12_run_pipeline.sh` script.

## Analysis Modules

The project includes several analytical modules, each focused on a specific aspect of mutation analysis:

### 1. Mutation Spectrum Analysis

- **Script**: `mutation_spectrum_analysis.py`
- **Purpose**: Examines patterns of single nucleotide mutations
- **Outputs**:
  - Mutation spectrum plots for each treatment
  - Comparative mutation spectrum across treatments
  - Transition/transversion ratios
  - Statistical tests for mutation pattern differences

### 2. Genomic Context Analysis

- **Script**: `genomic_context_analysis.py`
- **Purpose**: Studies the sequence context around mutations
- **Outputs**:
  - Sequence logos for different mutation types
  - GC content analysis in mutation regions
  - Homopolymer and repeat analysis
  - Mutation context heatmaps

### 3. Mutational Signatures Analysis

- **Script**: `mutational_signature_analysis.py`
- **Purpose**: Identifies characteristic patterns of mutations
- **Outputs**:
  - Mutation signature plots for each treatment
  - Enriched mutation contexts
  - Signature similarity analysis
  - Adaptation-specific signatures

### 4. Population Structure Analysis

- **Script**: `population_spectrum_analysis.py`
- **Purpose**: Examines genetic relationships between samples
- **Outputs**:
  - Principal Component Analysis (PCA) plots
  - Multi-dimensional scaling (MDS) plots
  - Dendrogram clustering
  - Variant sharing analysis

### 5. Regional Enrichment Analysis

- **Script**: `regional_enrichment_analysis.py`
- **Purpose**: Identifies genomic regions with mutation hotspots
- **Outputs**:
  - Enriched region tables for each treatment
  - Adaptation-specific enriched regions
  - Shared enrichment regions
  - Enrichment heatmaps

### 6. Scaffold Distribution Analysis

- **Script**: `scaffold_distribution_analysis.py`
- **Purpose**: Analyzes mutation distribution across scaffolds
- **Outputs**:
  - Variant density plots
  - Treatment correlation heatmaps
  - Bubble charts for variant distribution
  - Scaffold enrichment statistics

### 7. Statistical Pattern Analysis

- **Script**: `statistical_pattern_analysis.py`
- **Purpose**: Performs statistical tests on mutation patterns
- **Outputs**:
  - Correlation analyses
  - Regression models
  - PCA biplots
  - Clustering analyses

## Key Findings

The analyses have yielded several important findings:

### 1. Mutation Spectrum 

- Different treatments show distinct mutation spectra, with transition/transversion (Ti/Tv) ratios ranging from 0.22 (WT-37) to 0.58 (CAS)
- Temperature-adapted strains (WT-37 and CAS) show a preference for C>A mutations
- Low oxygen-adapted strains (WTA and STC) show a preference for C>G mutations
- Gene modifications influence mutation spectra, with gene-modified strains showing higher Ti/Tv ratios

### 2. Genomic Context

- GC content around mutations varies slightly between treatments (0.38-0.40)
- A high proportion of mutations occur near homopolymer regions (~92%)
- Temperature adaptation shows enrichment for specific trinucleotide contexts (TCT with C>A mutations)
- Gene-modified strains show different sequence context preferences compared to non-modified strains

### 3. Population Structure

- Samples cluster primarily by adaptation type in PCA and MDS analyses
- Temperature and low oxygen adaptations show distinct genetic profiles
- Within treatments, replicates show high similarity (Jaccard indices 0.60-0.73)
- Gene-modified strains show increased variant counts compared to non-modified strains

### 4. Regional Enrichment

- Several genomic regions show significant enrichment for mutations under specific treatments
- Regions on scaffolds CM007977.1 and CM007967.1 show particularly high fold enrichments
- Some enriched regions are shared across multiple treatments, suggesting common mutational hotspots
- Adaptation-specific regions have been identified for both temperature and low oxygen conditions

### 5. Scaffold Distribution

- Mutations are non-randomly distributed across the genome
- Small scaffolds (CM007980.1, LYZE01000020.1) show high variant density
- Treatment correlation analysis shows highest similarity between WT-37 and STC (ρ = 0.53)
- Temperature-adapted strains show higher global variant density (0.0032 variants/kb) than low oxygen-adapted strains (0.0025 variants/kb)

### 6. Treatment vs Control

- All treatments show significantly fewer variants compared to controls (p < 1e-70)
- Gene-modified strains show more variants than their non-modified counterparts (CAS: 38 vs WT-37: 27; STC: 27 vs WTA: 20)
- Replicates show high consistency, with 70-72% of variants present in all three replicates

## Future Work: Annotation with snpEff

The project is currently preparing for functional annotation of variants using snpEff. The preparation phase has been completed, as documented in `annotation/vcf_preparation_report.txt`. The next steps include:

1. **Building a Custom Database**:
   - Creating a custom snpEff database for the W303 yeast strain
   - Incorporating gene and feature annotations

2. **Functional Annotation**:
   - Annotating variants with predicted functional impacts
   - Categorizing variants (missense, nonsense, synonymous, etc.)
   - Prioritizing variants based on functional impact

3. **Pathway Analysis**:
   - Identifying enriched biological pathways affected by mutations
   - Linking mutations to adaptive responses
   - Investigating gene-environment interactions

4. **Comparative Functional Analysis**:
   - Comparing functional impacts between different treatments
   - Analyzing adaptation-specific functional changes
   - Investigating gene modification effects on functional outcomes

## Previous Annotation Analysis

The `old_version/` directory contains earlier approaches to annotation and analysis:

- **Old Reference Genome**: The previous version used alternative reference sequences
- **Old Analysis Scripts**: Earlier versions of analysis scripts with different methodologies
- **Previous Conclusions**: Initial findings and interpretations, which have been refined in the current version
- **Debugging Logs**: Records of troubleshooting and optimization steps
- **Alternative Approaches**: Different analytical approaches that were explored

The current version represents a significant improvement in terms of:
- More sophisticated analytical methods
- Better integration between different analysis modules
- Improved visualization and interpretation
- More comprehensive biological context

## Usage

To run the complete analysis pipeline:

```bash
bash scripts/pipeline/12_run_pipeline.sh
```

To run individual analysis modules:

```bash
python3 scripts/utils/run_analysis.py
```

The `run_analysis.py` script will execute all analysis modules in the correct order, ensuring dependencies are met.

## References

- Yeast W303 strain reference genome and annotations
- Statistical methods for variant analysis
- Mutation signature analysis methodologies
- Sequence context analysis techniques

---

*Note: This README documentation provides a comprehensive overview of the Yeast MSA project, including its purpose, methodology, results, and future directions. For detailed results from each analysis module, please refer to the corresponding files in the analysis/ directory.*