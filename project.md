# Yeast MSA (Multiple Sequence Alignment) Project

## Table of Contents
1. [Project Overview](#project-overview)
2. [Experimental Design](#experimental-design)
3. [Project Structure](#project-structure)
4. [Pipeline Overview](#pipeline-overview)
5. [Analysis Modules](#analysis-modules)
   - [Core Analysis](#core-analysis)
   - [Gene Analysis](#gene-analysis)
   - [Functional Impact Analysis](#functional-impact-analysis)
   - [Network Analysis](#network-analysis)
   - [Sterol Profile Analysis](#sterol-profile-analysis)
6. [Key Findings](#key-findings)
7. [Visualizations and Reports](#visualizations-and-reports)
8. [Usage](#usage)
   - [Script Execution Order](#script-execution-order)
   - [Analysis Output Examination Guide](#analysis-output-examination-guide)
9. [Dependencies](#dependencies)
10. [References](#references)

## Project Overview

This Yeast reference-guided multiple sequence analysis project investigates how yeast (S. cerevisiae, W303 strain) adapts to different environmental stresses through genetic mutations. The project employs comparative genomic analysis to study variations between different yeast strains subjected to different treatments (temperature, low oxygen) and gene modifications. By analyzing the mutational patterns, genomic contexts, and regional enrichments, this project aims to understand the molecular mechanisms underlying adaptation to environmental stress with particular focus on the ergosterol pathway's role in membrane adaptation.

The primary goals of this project include:
1. Identifying and characterizing mutations in different treatment conditions
2. Analyzing mutation spectrum patterns and sequence contexts
3. Detecting regional enrichment of mutations across the genome
4. Comparing population structures between different treatment groups
5. Understanding the influence of gene modifications on adaptation mechanisms
6. Analyzing functional impacts of variants in the ergosterol pathway genes
7. Investigating network relationships between core pathway genes and affected neighboring genes
8. Characterizing purifying selection and adaptive mechanisms at the gene level
9. Correlating genomic findings with biochemical adaptations in sterol profiles
10. Understanding the hierarchical conservation architecture in essential pathways

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
├── analysis/                        # Analysis results
│   ├── gene_analysis/                  # Gene-specific analysis
│   ├── general_gene_analysis/          # General gene analysis
│   │   ├── gene_genomic_context_results/
│   │   ├── gene_mutation_spectrum_results/
│   │   ├── gene_population_structure_results/
│   │   ├── gene_regional_enrichment_results/
│   │   ├── gene_statistical_pattern_results/
│   │   └── gene_treatment_control_analysis/
│   ├── genomic_context_results/        # Sequence context analysis
│   ├── genes_of_interest/              # Analysis of specific genes
│   ├── mutation_spectrum_results/      # Mutation spectrum analysis
│   ├── mutational_signatures_results/  # Mutational signatures
│   ├── population_structure_results/   # Population structure analysis
│   ├── regional_enrichment_results/    # Regional mutation hotspots
│   ├── scaffold_distribution_results/  # Variant distribution by scaffold
│   └── statistical_pattern_results/    # Statistical analyses
├── filtered_variants/               # Filtered variant data
├── mutation_spectrum_analysis/      # Mutation data for spectrum analysis
├── reference/                       # Reference genome and annotations
│   ├── chromosome_mapping.tsv         # Chromosome naming maps
│   ├── gene_mapping.tsv               # Basic gene annotation mapping
│   ├── gene_mapping_full.tsv          # Comprehensive gene annotation
│   ├── genes_of_interest_mapping.tsv  # Ergosterol pathway genes
│   ├── w303_chromosomal.fasta         # Main reference genome
│   └── w303_annotations/              # Genome annotations
├── results/                         # Pipeline processing results
│   ├── functional_impact/             # Functional impact analysis
│   ├── gene_variants/                 # Gene-specific variant analysis
│   │   ├── gene_specific/               # Variants by gene
│   │   └── treatment_specific/          # Variants by treatment
│   ├── gene_variants_expanded/        # Expanded gene variant analysis
│   ├── gene_verification/             # Gene coordinate verification
│   ├── network_analysis/              # Gene interaction networks
│   │   └── erg_subnetworks/             # Ergosterol gene subnetworks
│   ├── scaffold_variants/             # Scaffold-specific variant analysis
│   ├── stats/                         # Alignment statistics
│   ├── sterol_analysis/               # Sterol profile analysis
│   │   ├── basic_stats/                 # Basic sterol statistics
│   │   ├── comparative/                 # Treatment comparisons
│   │   ├── correlation/                 # Genomic-sterol correlations
│   │   ├── pathway/                     # Pathway analysis results
│   │   └── visualizations/              # Sterol profile visualizations
│   └── treatment_analysis/            # Treatment comparison results
├── scripts/                         # Analysis and pipeline scripts
│   ├── analysis/                      # Core analysis scripts
│   ├── annotation/                    # Annotation scripts
│   │   ├── new/                         # Updated annotation pipeline
│   │   └── old/                         # Legacy annotation scripts
│   ├── functional_impact/             # Functional impact scripts
│   ├── gene_analysis/                 # Gene-specific analysis
│   ├── general_gene_analysis/         # General gene analysis
│   ├── pipeline/                      # Pipeline execution scripts
│   ├── sterols/                       # Sterol analysis scripts
│   │   ├── preprocess_sterols.py        # Sterol data preprocessing
│   │   ├── sterol_analysis.py           # Core sterol analysis
│   │   ├── sterol_pathway.py            # Pathway analysis 
│   │   ├── sterol_integration.py        # Integration with genomic data
│   │   └── run_sterol_analysis.sh       # Complete sterol analysis pipeline
│   ├── utils/                         # Utility scripts
│   ├── variants/                      # Variant extraction scripts
│   └── vcf/                           # VCF processing scripts
└── vcf/                             # VCF files
    ├── annotated/                     # Annotated VCF files
    ├── filtered/                      # Filtered VCF files
    ├── individual/                    # Individual sample VCFs
    ├── initial/                       # Initial VCF files
    └── renamed/                       # Renamed VCF files with standardized chromosome IDs
```

## Pipeline Overview

The complete analysis pipeline consists of two major phases: initial processing (preprocessing, alignment, variant calling) and variant analysis. Both phases are essential to the project's workflow.

### Phase 1: Initial Data Processing

Before running the main pipeline, several preprocessing steps are performed to prepare raw sequencing data for analysis:

1. **Quality Control** (run_fastqc.sh):
   - Assesses raw read quality with FastQC
   - Generates quality reports for all samples
   - Optionally aggregates results with MultiQC

2. **Read Preprocessing** (preprocessing.sh):
   - Trims and filters raw reads using fastp
   - Removes adapters automatically
   - Applies quality filters (Q20)
   - Removes reads shorter than 50bp after trimming

3. **Alignment** (alignment.sh):
   - Aligns reads to the W303 reference genome using BWA-MEM
   - Filters and sorts aligned reads
   - Marks and removes duplicate reads
   - Adds read group information
   - Generates alignment statistics

4. **Variant Calling** (variants.sh):
   - Calls variants using bcftools mpileup and call
   - Uses haploid model appropriate for yeast
   - Normalizes indel representations
   - Filters variants by quality (QUAL≥20) and depth (DP≥10)
   - Generates per-sample variant statistics

5. **Initial Comparison** (vcf_comparison.sh):
   - Performs initial comparisons between treatments and controls
   - Identifies treatment-specific variants
   - Evaluates consistency within replicates
   - Generates preliminary summary reports

6. **Direct Comparison** (direct_comparison.sh):
   - Implements a more reliable approach to variant comparison
   - Compares treatment samples directly against controls
   - Identifies high-confidence variants present in multiple replicates
   - Creates cross-treatment comparison matrices
   - Identifies treatment-unique variants

7. **Consistency Analysis** (consistency.sh):
   - Analyzes the consistency of variants within treatment replicates
   - Categorizes variants by presence in 1, 2, or all 3 replicates
   - Provides statistical summaries of variant reproducibility
   - Generates consistency reports for each treatment group

8. **Position Analysis** (position_analysis.sh):
   - Analyzes the distribution of variant positions
   - Identifies potential clustering of mutations
   - Calculates distances between variants
   - Generates statistics on variant spacing

### Phase 2: Main Analysis Pipeline

After the initial processing, the main pipeline is executed to perform comprehensive analyses:

1. **Setup** (01_setup.sh):  
   - Creates the directory structure for results

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

The complete main pipeline can be executed using the `12_run_pipeline.sh` script.

## Analysis Modules

### Core Analysis

The project includes several core analytical modules, each focused on a specific aspect of mutation analysis:

#### 1. Mutation Spectrum Analysis

- **Script**: `mutation_spectrum_analysis.py`
- **Purpose**: Examines patterns of single nucleotide mutations
- **Outputs**:
  - Mutation spectrum plots for each treatment
  - Comparative mutation spectrum across treatments
  - Transition/transversion ratios
  - Statistical tests for mutation pattern differences

#### 2. Genomic Context Analysis

- **Script**: `genomic_context_analysis.py`
- **Purpose**: Studies the sequence context around mutations
- **Outputs**:
  - Sequence logos for different mutation types
  - GC content analysis in mutation regions
  - Homopolymer and repeat analysis
  - Mutation context heatmaps

#### 3. Mutational Signatures Analysis

- **Script**: `mutational_signature_analysis.py`
- **Purpose**: Identifies characteristic patterns of mutations
- **Outputs**:
  - Mutation signature plots for each treatment
  - Enriched mutation contexts
  - Signature similarity analysis
  - Adaptation-specific signatures

#### 4. Population Structure Analysis

- **Script**: `population_spectrum_analysis.py`
- **Purpose**: Examines genetic relationships between samples
- **Outputs**:
  - Principal Component Analysis (PCA) plots
  - Multi-dimensional scaling (MDS) plots
  - Dendrogram clustering
  - Variant sharing analysis

#### 5. Regional Enrichment Analysis

- **Script**: `regional_enrichment_analysis.py`
- **Purpose**: Identifies genomic regions with mutation hotspots
- **Outputs**:
  - Enriched region tables for each treatment
  - Adaptation-specific enriched regions
  - Shared enrichment regions
  - Enrichment heatmaps

#### 6. Scaffold Distribution Analysis

- **Script**: `scaffold_distribution_analysis.py`
- **Purpose**: Analyzes mutation distribution across scaffolds
- **Outputs**:
  - Variant density plots
  - Treatment correlation heatmaps
  - Bubble charts for variant distribution
  - Scaffold enrichment statistics

#### 7. Statistical Pattern Analysis

- **Script**: `statistical_pattern_analysis.py`
- **Purpose**: Performs statistical tests on mutation patterns
- **Outputs**:
  - Correlation analyses
  - Regression models
  - PCA biplots
  - Clustering analyses

### Gene Analysis

A major enhancement to the project is the addition of gene-specific analysis modules:

#### 1. Gene-Specific Analysis

- **Script Directory**: `scripts/gene_analysis/`
- **Run Script**: `run_gene_analysis.sh`
- **Purpose**: Analyzes mutations in specific genes of interest (especially ergosterol pathway genes)
- **Outputs**:
  - Gene-specific mutation profiles
  - Mutation distribution in ergosterol pathway genes
  - Gene-focused sequence context analysis
  - Gene variant statistics

#### 2. General Gene Analysis

- **Script Directory**: `scripts/general_gene_analysis/`
- **Run Script**: `run_general_gene_analysis.sh`
- **Purpose**: Provides a comprehensive analysis of mutations across all genes
- **Key Components**:
  - `generate_gene_mapping_full.py`: Creates comprehensive gene mapping
  - Gene-specific versions of core analysis modules
- **Outputs**:
  - Gene-centric mutation spectrum results
  - Gene-specific genomic context analysis
  - Gene population structure results
  - Gene regional enrichment results
  - Gene statistical pattern analysis
  - Gene treatment control analysis

### Functional Impact Analysis

A significant addition to the project is the functional impact analysis module:

#### 1. High Impact Variant Analysis

- **Script**: `analyze_high_impact_variants.py`
- **Run Script**: `run_high_impact_analysis.sh`
- **Purpose**: Identifies and analyzes variants with significant functional impacts
- **Features**:
  - Analyzes HIGH and MODERATE impact variants
  - Focuses on ergosterol pathway genes
  - Predicts functional consequences of amino acid changes
  - Maps variants to protein domains
- **Outputs**:
  - High impact variant tables
  - Protein domain impact maps
  - Functional consequence predictions
  - Treatment-specific high impact profiles

#### 2. Treatment-Specific Pattern Analysis

- **Script**: `analyze_treatment_specific_patterns.py`
- **Run Script**: `run_treatment_analysis.sh`
- **Purpose**: Identifies treatment-specific functional impact patterns
- **Outputs**:
  - Treatment-specific functional impact profiles
  - Comparative treatment analyses
  - Statistical significance of treatment-specific patterns

#### 3. Regulatory Element Analysis

- **Scripts**: 
  - `analyze_promoter_elements.py`: Analyzes upstream variants and promoter regions
  - `analyze_tfbs.py`: Identifies impacts on transcription factor binding sites
- **Purpose**: Characterizes variants affecting gene regulation
- **Features**:
  - Maps variants relative to transcription start sites
  - Analyzes distribution patterns of regulatory variants
  - Identifies affected transcription factor binding motifs
  - Predicts regulatory impact based on conservation and position
- **Outputs**:
  - Promoter variant distribution visualizations (`results/regulatory_analysis/promoters/`)
  - TSS distance distribution plots
  - Position-specific enrichment heatmaps
  - Stress response element (SRE) overlap analysis
  - Transcription factor binding site analysis (`results/regulatory_analysis/tfbs/`)

### Network Analysis

The network analysis module is a new addition that maps relationships between genes:

#### 1. Extended Ergosterol Network Analysis

- **Script**: `build_extended_erg_network.py`
- **Run Script**: `run_extended_network_analysis.sh`
- **Purpose**: Constructs and analyzes gene interaction networks
- **Features**:
  - Maps ergosterol pathway genes
  - Identifies affected neighboring genes
  - Calculates gene-gene relationships
  - Visualizes network topology
- **Outputs**:
  - Extended network visualization
  - Treatment-specific network maps
  - Gene subnetwork visualizations
  - Network statistics (centrality, connectivity)
  - Network analysis report

#### 2. Variant Proximity Impact Analysis

- **Script**: `variant_proximity_impact_summary.py`
- **Purpose**: Analyzes the relationship between variant impact and distance to ergosterol genes
- **Features**:
  - Measures minimum distances between variants and ERG genes
  - Groups variants by treatment, impact, and position relative to ERG genes
  - Correlates variant impact with genomic distance
  - Generates visualization of distance-impact relationships
- **Outputs**:
  - Variant proximity impact summary table
  - Variant count heatmap by ERG gene and treatment
  - Distance distribution boxplots
  - Interactive HTML report with distance-impact visualizations

### OSH Gene Family Analysis

The OSH gene family analysis module examines the relationship between sterol synthesis and transport:

#### 1. OSH Gene Analysis

- **Script Directory**: `scripts/osh_analysis/`
- **Run Script**: `run_osh_analysis.sh`
- **Purpose**: Investigates the OSH (OxySterol binding Homology) gene family and its relationship to ergosterol genes
- **Components**:
  - `analyze_osh_genes.py`: Maps OSH family genes in the reference genome
  - `osh_variants.py`: Analyzes variants in and around OSH genes in all treatment conditions
  - `osh_erg_distance.py`: Calculates genomic distances between OSH and ergosterol pathway genes
  - `count_variant_locations.py`: Analyzes the distribution of variants relative to OSH genes
- **Outputs**:
  - OSH gene summary table
  - OSH variant analysis results
  - OSH-ERG distance calculations
  - Visualizations of OSH-ERG relationships
  - OSH gene analysis report (OSH_Results.md)

### Sterol Profile Analysis

The sterol profile analysis module connects genomic findings with biochemical data on yeast membrane composition:

#### 1. Sterol Preprocessing and Basic Statistics

- **Script**: `preprocess_sterols.py`
- **Run Script**: `run_sterol_analysis.sh`
- **Purpose**: Processes and characterizes sterol profile data
- **Features**:
  - Processes raw sterol concentration measurements
  - Adds metadata (treatment, adaptation type, gene modification status)
  - Calculates relative abundances and basic statistics
  - Identifies sterol composition by treatment
- **Outputs**:
  - Standardized sterol data
  - Sterol distributions by treatment and adaptation
  - Concentration statistics and visualization
  - Sterol diversity metrics

#### 2. Comparative Sterol Analysis

- **Script**: `sterol_analysis.py`
- **Purpose**: Compares sterol profiles across conditions
- **Features**:
  - Statistical testing for treatment differences
  - Fold change calculations
  - Adaptation-specific sterols identification
  - Gene modification effects analysis
- **Outputs**:
  - Statistical test results and p-values
  - Fold change heatmaps and visualizations
  - Treatment-specific sterol profiles
  - Unique sterol identification by condition

#### 3. Sterol Pathway Analysis

- **Script**: `sterol_pathway.py`
- **Purpose**: Connects sterol profiles to biochemical pathway steps
- **Features**:
  - Maps sterol intermediates to pathway enzymes
  - Calculates substrate/product ratios
  - Identifies altered pathway branches
  - Analyzes pathway flux patterns
- **Outputs**:
  - Pathway mapping visualizations
  - Adaptation-specific pathway diagrams
  - Flux analysis by treatment
  - Enzyme activity inference

#### 4. Genomic-Sterol Integration

- **Script**: `sterol_integration.py`
- **Purpose**: Integrates sterol profiles with genomic conservation patterns
- **Features**:
  - Correlates sterol changes with variant patterns
  - Connects satellite genes to sterol production
  - Tests hierarchical conservation models
  - Builds integrated adaptation models
- **Outputs**:
  - Correlation analyses
  - Satellite gene-sterol connection maps
  - Integrated visualizations
  - Comprehensive adaptation model

## Combined Analysis Results

The project generates a comprehensive summary of findings in the `combined_analysis_results.txt` file, which integrates results from all analysis modules. This file provides:

1. **Variant Analysis Summary**:
   - Treatment-specific variant counts and high-confidence variants
   - Variant consistency statistics within replicates (70-72% in all 3 replicates)
   - Comparison between different variant analysis methods
   - Cross-treatment comparisons showing shared and unique variants

2. **Genomic Context Results**:
   - GC content analysis (0.38-0.40 across treatments)
   - Homopolymer region analysis (~92% of variants near homopolymers)
   - Dinucleotide repeat analysis
   - Treatment-specific context preferences

3. **Mutation Spectrum Results**:
   - Transition/transversion ratios by treatment
   - Most common mutation types
   - Statistical comparisons of mutation patterns
   - Adaptation-specific mutation preferences

4. **Population Structure Results**:
   - Variant sharing statistics by treatment and adaptation
   - Genetic similarity measurements
   - Clustering and PCA analysis interpretation

5. **Mutational Signatures Results**:
   - Treatment signature similarities
   - Adaptation-specific signature analysis
   - Enriched trinucleotide contexts

6. **Scaffold Distribution Analysis**:
   - Variant density analysis by scaffold
   - Treatment correlations
   - Gene modification effects

7. **Regional Enrichment Results**:
   - Enriched regions identification
   - Adaptation-specific enrichment
   - Shared enrichment regions

8. **Statistical Pattern Results**:
   - Regression analyses
   - Treatment vs control statistical tests
   - Fold changes and significance metrics

This combined analysis file serves as a comprehensive reference for all key results from the project.
