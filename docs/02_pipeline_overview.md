# Pipeline Overview

The Yeast MSA project follows a comprehensive analysis pipeline divided into two major phases:

## Phase 1: Initial Data Processing

This phase prepares the raw sequencing data for analysis and performs initial variant calling:

### 1. Quality Control
- **Tool**: FastQC
- **Script**: `scripts/utils/run_fastqc.sh`
- **Purpose**: Assess raw read quality
- **Output**: Quality reports for all samples
- **Key Metrics**: Per-base quality, sequence duplication, adapter content

### 2. Read Preprocessing
- **Tool**: fastp
- **Script**: `scripts/utils/preprocessing.sh`
- **Purpose**: Clean and filter raw reads
- **Operations**:
  - Adapter removal
  - Quality filtering (Q20)
  - Length filtering (min 50bp)
  - Duplicate removal
- **Output**: Cleaned FASTQ files

### 3. Alignment
- **Tool**: BWA-MEM
- **Script**: `scripts/utils/alignment.sh`
- **Reference**: W303 yeast genome (`w303_chromosomal.fasta`)
- **Operations**:
  - Align reads to reference
  - Sort aligned reads
  - Mark and remove duplicates
  - Add read group information
- **Output**: BAM files and alignment statistics

### 4. Variant Calling
- **Tool**: bcftools (mpileup and call)
- **Script**: `scripts/utils/variants.sh`
- **Model**: Haploid (appropriate for yeast)
- **Filtering**:
  - Quality ≥ 20
  - Depth ≥ 10
- **Output**: Individual VCF files

### 5. Initial Comparison
- **Script**: `scripts/utils/vcf_comparison.sh`
- **Purpose**: Preliminary treatment vs control comparison
- **Output**: Initial variant comparison statistics

### 6. Direct Comparison
- **Script**: `scripts/utils/direct_comparison.sh`
- **Purpose**: Reliable approach to variant comparison
- **Operations**:
  - Direct treatment vs control comparison
  - Multi-replicate validation
  - Cross-treatment comparison
- **Output**: High-confidence variant lists

### 7. Consistency Analysis
- **Script**: `scripts/utils/consistency.sh`
- **Purpose**: Assess variant reproducibility
- **Operations**:
  - Categorize variants by presence in 1, 2, or 3 replicates
  - Calculate consistency statistics
- **Output**: Consistency reports

### 8. Position Analysis
- **Script**: `scripts/utils/position_analysis.sh`
- **Purpose**: Analyze variant distribution
- **Operations**:
  - Calculate variant spacing
  - Identify potential clustering
- **Output**: Position statistics

## Phase 2: Main Analysis Pipeline

This phase performs in-depth analysis of the variants and their biological significance:

### 1. Setup
- **Script**: `scripts/pipeline/01_setup.sh`
- **Purpose**: Create directory structure for results

### 2. Filtering & Normalization
- **Script**: `scripts/pipeline/02_filter_normalize.sh`
- **Operations**:
  - Apply quality filters
  - Normalize variant representation
- **Output**: Filtered and normalized VCF files

### 3. Contig Standardization
- **Script**: `scripts/pipeline/03_fix_contigs.sh`
- **Purpose**: Ensure consistent contig naming
- **Output**: VCF files with standardized chromosome IDs

### 4. VCF Merging
- **Script**: `scripts/pipeline/04_merge_vcfs.sh`
- **Purpose**: Combine VCF files for analysis
- **Output**: Merged VCF datasets

### 5. Sample Validation
- **Script**: `scripts/pipeline/05_test_samples.sh`
- **Purpose**: Validate sample selection
- **Output**: Quality control metrics

### 6. Group Comparison
- **Script**: `scripts/pipeline/06_compare_groups.sh`
- **Purpose**: Compare treatment vs control
- **Output**: Treatment-specific variants

### 7. Cross-treatment Comparison
- **Script**: `scripts/pipeline/07_cross_treatment.sh`
- **Purpose**: Analyze differences between treatments
- **Output**: Shared and unique variants by treatment

### 8. Unique Variant Identification
- **Script**: `scripts/pipeline/08_identify_unique.sh`
- **Purpose**: Find adaptation-specific mutations
- **Output**: Unique variant lists

### 9. Direct Comparison
- **Script**: `scripts/pipeline/09_direct_comparison.sh`
- **Purpose**: Perform pairwise comparisons
- **Output**: Similarity metrics

### 10. Consistency Analysis
- **Script**: `scripts/pipeline/10_consistency_analysis.sh`
- **Purpose**: Measure replicate consistency
- **Output**: Detailed consistency statistics

### 11. Summary Report
- **Script**: `scripts/pipeline/11_summary_report.sh`
- **Purpose**: Generate final report
- **Output**: Comprehensive analysis summary

### Master Pipeline
- **Script**: `scripts/pipeline/12_run_pipeline.sh`
- **Purpose**: Execute all pipeline steps
- **Usage**: `bash scripts/pipeline/12_run_pipeline.sh`

## Specialized Analysis Modules

Beyond the core pipeline, several specialized analysis modules are executed:

### Gene Analysis
- **Script**: `scripts/gene_analysis/run_gene_analysis.sh`
- **Purpose**: Analyze ergosterol pathway genes
- **Output**: Gene-specific mutation profiles

### General Gene Analysis
- **Script**: `scripts/general_gene_analysis/run_general_gene_analysis.sh`
- **Purpose**: Genome-wide gene analysis
- **Output**: Comprehensive gene-level results

### Functional Impact Analysis
- **Script**: `scripts/functional_impact/run_high_impact_analysis.sh`
- **Purpose**: Assess variant functional impacts
- **Output**: Impact predictions and visualizations

### Network Analysis
- **Script**: `scripts/functional_impact/run_extended_network_analysis.sh`
- **Purpose**: Analyze gene interaction networks
- **Output**: Network visualizations and statistics

### Sterol Analysis
- **Script**: `scripts/sterols/run_sterol_analysis.sh`
- **Purpose**: Analyze sterol profiles
- **Output**: Sterol-genomic correlations

## Data Flow

The complete data flow follows this sequence:

1. Raw FASTQ → FastQC → Quality reports
2. Raw FASTQ → fastp → Cleaned FASTQ
3. Cleaned FASTQ → BWA-MEM → BAM
4. BAM → bcftools → Individual VCF
5. Individual VCF → Filtering → Filtered VCF
6. Filtered VCF → Merging → Combined VCF
7. VCF → Annotation → Annotated VCF
8. VCF → Analysis modules → Results & visualizations
9. Sterol data → Preprocessing → Standardized profiles
10. Sterol profiles → Analysis → Treatment comparisons
11. Genomic + Sterol data → Integration → Adaptation model

## Directory Structure

Results are organized in a structured directory hierarchy:

- `analysis/`: Analysis results by module
- `reference/`: Reference genomes and annotations
- `results/`: Pipeline processing outputs
- `scripts/`: Analysis and pipeline scripts
- `vcf/`: Variant Call Format files

Each analysis module has a dedicated results directory with visualizations and summary statistics.

## Execution Order

For reproducing the analysis, scripts should be executed in this order:

1. Initial data processing scripts (Phase 1)
2. Main pipeline (Phase 2)
3. Gene-specific analysis
4. General gene analysis
5. Functional impact analysis
6. Network analysis
7. Sterol profile analysis

Alternatively, the master pipeline script can be used to run the core analysis in a single step.