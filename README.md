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
6. [Key Findings](#key-findings)
7. [Additional Analyses](#additional-analyses)
8. [Usage](#usage)
9. [References](#references)

## Project Overview

This Yeast reference-guided multiple sequence analysis project investigates how yeast (S. cerevisiae, W303 strain) adapts to different environmental stresses through genetic mutations. The project employs comparative genomic analysis to study variations between different yeast strains subjected to different treatments (temperature, low oxygen) and gene modifications. By analyzing the mutational patterns, genomic contexts, and regional enrichments, this project aims to understand the molecular mechanisms underlying adaptation to environmental stress.

The primary goals of this project include:
1. Identifying and characterizing mutations in different treatment conditions
2. Analyzing mutation spectrum patterns and sequence contexts
3. Detecting regional enrichment of mutations across the genome
4. Comparing population structures between different treatment groups
5. Understanding the influence of gene modifications on adaptation mechanisms
6. Analyzing functional impacts of variants in the ergosterol pathway genes
7. Investigating network relationships between core pathway genes and affected neighboring genes
8. Characterizing purifying selection and adaptive mechanisms at the gene level

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
│   ├── chromosome_mapping.tsv            # Chromosome naming maps
│   ├── gene_mapping.tsv                  # Basic gene annotation mapping
│   ├── gene_mapping_full.tsv             # Comprehensive gene annotation
│   ├── genes_of_interest_mapping.tsv     # Ergosterol pathway genes
│   ├── w303_chromosomal.fasta           # Main reference genome
│   └── w303_annotations/                # Genome annotations
├── results/                         # Pipeline processing results
│   ├── functional_impact/               # Functional impact analysis
│   ├── gene_variants/                   # Gene-specific variant analysis
│   │   ├── gene_specific/                  # Variants by gene
│   │   └── treatment_specific/             # Variants by treatment
│   ├── gene_variants_expanded/          # Expanded gene variant analysis
│   ├── gene_verification/               # Gene coordinate verification
│   ├── network_analysis/                # Gene interaction networks
│   │   └── erg_subnetworks/                # Ergosterol gene subnetworks
│   ├── scaffold_variants/               # Scaffold-specific variant analysis
│   ├── stats/                           # Alignment statistics
│   └── treatment_analysis/              # Treatment comparison results
├── scripts/                         # Analysis and pipeline scripts
│   ├── analysis/                        # Core analysis scripts
│   ├── annotation/                      # Annotation scripts
│   │   ├── new/                            # Updated annotation pipeline
│   │   └── old/                            # Legacy annotation scripts
│   ├── functional_impact/               # Functional impact scripts
│   ├── gene_analysis/                   # Gene-specific analysis
│   ├── general_gene_analysis/           # General gene analysis
│   ├── pipeline/                        # Pipeline execution scripts
│   ├── utils/                           # Utility scripts
│   ├── variants/                        # Variant extraction scripts
│   └── vcf/                             # VCF processing scripts
└── vcf/                             # VCF files
    ├── annotated/                       # Annotated VCF files
    ├── filtered/                        # Filtered VCF files
    ├── individual/                      # Individual sample VCFs
    ├── initial/                         # Initial VCF files
    └── renamed/                         # Renamed VCF files with standardized chromosome IDs
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

### 7. Gene Analysis

- Strong evidence of purifying selection in ergosterol pathway genes
- No HIGH or MODERATE impact variants found in core ergosterol pathway genes
- Gene modifications influence the distribution of variants around pathway genes
- Non-genic regions show higher mutation rates than genic regions

### 8. Network Analysis

- Extended ergosterol network reveals adaptation through changes in neighboring genes
- Gene-gene interactions suggest coordinated adaptation response
- Treatment-specific network patterns identified
- HIGH impact variants consistently found at specific distances from pathway genes

## Additional Analyses

### Functional Impact Analysis

The functional impact analysis provides insights into how mutations affect protein function:

1. **Protein Domain Impact**: 
   - Maps variants to functional domains in proteins
   - Predicts functional consequences based on amino acid properties
   - Prioritizes variants based on domain criticality

2. **Treatment-Specific Functional Patterns**:
   - Identifies treatment-specific functional impacts
   - Compares functional impacts across treatments
   - Links functional changes to adaptation mechanisms

### SNP Effect Analysis

SNP effect analysis has been performed using snpEff:

1. **Variant Annotation**:
   - Variants annotated with predicted functional effects
   - Classification of variants by impact severity (HIGH, MODERATE, LOW, MODIFIER)
   - Gene and transcript-level annotations

2. **Effect Distribution**:
   - Most variants (>90%) classified as MODIFIER impact
   - LOW impact variants account for ~5% of total
   - MODERATE impact variants at ~3%
   - HIGH impact variants are rare (<1%)

### Ergosterol Pathway Analysis

Specific analysis of the ergosterol pathway provides insights into this essential cellular process:

1. **Pathway Conservation**:
   - Strong evidence of purifying selection in pathway genes
   - No direct HIGH impact mutations in core pathway genes
   - Adaptation occurs through changes in regulatory regions and neighboring genes

2. **Extended Network Effects**:
   - Identifies genes functionally connected to the ergosterol pathway
   - Maps how these connected genes change in response to different stressors
   - Shows how gene modifications (STC, CAS) influence the extended network

## Usage

### Initial Data Processing

Before running the main analysis pipeline, the raw sequencing data must be processed:

```bash
# 1. Quality control check on raw reads
bash scripts/utils/run_fastqc.sh

# 2. Preprocess and clean raw reads
bash scripts/utils/preprocessing.sh

# 3. Align cleaned reads to reference genome
bash scripts/utils/alignment.sh

# 4. Call variants from aligned reads
bash scripts/utils/variants.sh

# 5. Perform initial variant comparison
bash scripts/utils/vcf_comparison.sh

# 6. Perform direct variant comparison 
bash scripts/utils/direct_comparison.sh

# 7. Analyze variant consistency within replicates
bash scripts/utils/consistency.sh

# 8. Analyze variant positions
bash scripts/utils/position_analysis.sh
```

### Running the Complete Pipeline

After the initial processing is complete, run the full analysis pipeline:

```bash
bash scripts/pipeline/12_run_pipeline.sh
```

### Running Gene Analysis

To run the gene-specific analysis:

```bash
bash scripts/gene_analysis/run_gene_analysis.sh
```

### Running General Gene Analysis

To run the comprehensive gene analysis:

```bash
bash scripts/general_gene_analysis/run_general_gene_analysis.sh
```

### Running Functional Impact Analysis

To analyze high impact variants:

```bash
bash scripts/functional_impact/run_high_impact_analysis.sh
```

### Running Network Analysis

To perform network analysis:

```bash
bash scripts/functional_impact/run_extended_network_analysis.sh
```

### Individual Analysis Modules

To run individual analysis modules:

```bash
python3 scripts/utils/run_analysis.py --module [module_name]
```

The `run_analysis.py` script will execute specific analysis modules, ensuring dependencies are met.

### Data Flow

The full data flow in the project follows this sequence:

1. Raw FASTQ files → FastQC → Quality reports
2. Raw FASTQ files → fastp → Cleaned FASTQ files
3. Cleaned FASTQ files → BWA-MEM → SAM → BAM files
4. BAM files → bcftools → Individual VCF files
5. Individual VCF files → Filtering → Normalized & filtered VCF files
6. Filtered VCF files → Merging → Combined VCF for analysis
7. VCF files → Annotation → Annotated VCF files
8. VCF files → Various analysis modules → Results and visualizations

## References

- Yeast W303 strain reference genome and annotations
- Statistical methods for variant analysis
- Mutation signature analysis methodologies
- Sequence context analysis techniques
- Network analysis approaches for genomic data
- Functional impact prediction models
- Ergosterol pathway biochemistry and regulation

---

*Note: This documentation provides a comprehensive overview of the Yeast MSA project, including its purpose, methodology, results, and future directions. For detailed results from each analysis module, please refer to the corresponding files in the analysis/ directory.*