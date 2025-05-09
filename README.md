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
9. [References](#references)

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

## Key Findings

The analyses have yielded several important findings:

### SNP Effect Analysis

Comprehensive SNP effect analysis was performed using snpEff to annotate and predict the functional effects of variants:

- **Variant Effect Distribution**:
  - MODIFIER impact: 90.2% (primarily intergenic and intronic variants)
  - LOW impact: 5.8% (synonymous variants and other minor effects)
  - MODERATE impact: 3.4% (missense variants primarily)
  - HIGH impact: 0.6% (stop gained, frameshift, splice site disruptors)

- **Gene Region Distribution**:
  - Intergenic: 62.7% 
  - Exonic: 18.3%
  - Intronic: 12.5%
  - Upstream/Downstream: 6.5%

- **Treatment-Specific Patterns**:
  - Temperature-adapted strains show higher proportion of MODERATE impact variants
  - Gene-modified strains show more HIGH impact variants in non-pathway genes
  - Core ergosterol pathway genes show almost exclusively MODIFIER variants

- **Codon Changes**:
  - Most common amino acid changes: Ala→Val, Ser→Asn, Gly→Asp
  - Temperature adaptation shows enrichment for hydrophobic substitutions
  - Low oxygen adaptation shows enrichment for charged residue substitutions

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

### 9. Sterol Profile Analysis

- Temperature adaptation shows significantly higher ergosterol levels (10.25 mean concentration) compared to low oxygen adaptation (2.73 mean concentration, p=0.0109)
- Gene-modified strains have 1.5× more diverse sterol profiles than non-modified strains
- Treatment-specific sterol signatures: Temperature adaptation has 7 unique sterols; low oxygen adaptation has unique Tetrahymanol marker
- Sterol compositions validate the hierarchical conservation pattern found in genomic analysis
- "Satellite genes" at consistent distances (7-50kb) from ergosterol pathway genes correlate with specific sterol production
- Evidence of adaptation through altered sterol profiles despite perfect conservation of pathway genes
- Integration reveals a four-layered conservation architecture from absolute conservation (core) to adaptive flexibility (satellite zone)

## Visualizations and Reports

The project generates a rich variety of reports and visualizations that integrate findings across multiple analysis modules. These outputs provide different views of the data, from raw statistics to interactive visualizations.

### Analysis Reports

Text-based summary reports provide detailed analysis of findings:

1. **Core Analysis Summaries**:
   - `combined_analysis_results.txt`: Master summary of all findings across modules
   - `analysis/mutation_spectrum_results/statistical_test_results.txt`: Statistical analysis of mutation patterns
   - `analysis/mutational_signatures_results/mutational_signatures_summary.txt`: Signature analysis findings
   - `analysis/population_structure_results/population_structure_summary.txt`: Population analysis results
   - `analysis/regional_enrichment_results/regional_enrichment_summary.txt`: Enrichment findings
   - `analysis/statistical_pattern_results/statistical_analysis_summary.txt`: Statistical analysis report

2. **Pipeline Processing Reports**:
   - `results/vcf/analysis_report.txt`: Summary of variant statistics across treatments
   - `results/scaffold_variants/scaffold_variant_summary.txt`: Analysis of variants by scaffold
   - `results/gene_variants/gene_summary.tsv`: Summary of gene-specific variant patterns
   - `results/treatment_analysis/analysis_summary.txt`: Treatment comparison summary
   - `vcf/annotated/annotation_summary.txt`: Summary of variant annotations

3. **Sterol Analysis Reports**:
   - `results/sterol_analysis/basic_stats/sterol_statistics.txt`: Basic sterol composition statistics
   - `results/sterol_analysis/comparative/comparative_analysis_summary.txt`: Treatment comparisons
   - `results/sterol_analysis/correlation/statistical_correlation_results.txt`: Genomic-sterol correlations
   - `results/sterol_analysis/pathway/pathway_analysis_summary.txt`: Sterol pathway flux analysis

### Interactive HTML Reports

Interactive HTML dashboards provide rich visualizations and explorable data:

1. **Ergosterol Variant Analysis Dashboard**:
   - **Path**: `results/reports/ergosterol_variant_analysis.html`
   - **Generator**: `scripts/utils/generate_ergosterol_variant_report.py`
   - **Contents**:
     - Interactive variant distribution visualizations
     - Distance-based analysis from pathway genes
     - Treatment-specific variant patterns
     - Purifying selection evidence
     - Biological significance interpretation
   - **Features**:
     - Interactive image gallery with zoom capability
     - Tabbed visualization sections
     - Dark/light mode toggle
     - Mobile-responsive design

2. **Functional Impact Analysis Dashboard**:
   - **Path**: `results/reports/functional_impact.html`
   - **Generator**: `scripts/utils/generate_functional_impact_report.py`
   - **Contents**:
     - Protein domain impact visualizations
     - Conservation patterns analysis
     - Satellite gene architecture diagrams
     - Adaptation mechanisms models
     - Hierarchical conservation zone visualizations

3. **Sterol Profile Analysis Dashboard**:
   - **Path**: `results/reports/sterols.html`
   - **Generator**: `scripts/sterols/generate_html_report.py`
   - **Contents**:
     - Interactive sterol profile visualizations
     - Treatment comparison heatmaps
     - Pathway visualizations with flux indicators
     - Genomic-sterol integration diagrams
     - Adaptation model visualizations

4. **Variant Analysis Dashboard**:
   - **Path**: `results/reports/variant_analysis.html`
   - **Generator**: `scripts/variants/generate_variant_report.py`
   - **Contents**:
     - Sample comparison visualizations
     - Interactive filtering of variants
     - Annotation statistics and distribution
     - Genome browser-like variant visualization
     - Treatment-specific variant patterns

### Markdown Reports

Detailed narrative analysis with embedded visualizations:

1. **Network Analysis Report**:
   - **Path**: `results/network_analysis/network_analysis_report.md`
   - **Contents**:
     - Extended ergosterol network analysis
     - Subnetwork interactions
     - Treatment-specific network patterns
     - Pathway distance models
     - Network statistics and interpretations

2. **Integrated Findings Report**:
   - **Path**: `results/sterol_analysis/correlation/integrated_findings_report.md`
   - **Contents**:
     - Integration of sterol and genomic findings
     - Hierarchical conservation model explanation
     - Satellite gene-sterol connections
     - Four-zone conservation architecture
     - Comprehensive adaptation model

3. **High Impact Variants Report**:
   - **Path**: `results/functional_impact/high_impact/high_impact_variants_report.md`
   - **Contents**:
     - Analysis of high impact variants
     - Functional domain mappings
     - Evolutionary conservation patterns
     - Treatment-specific functional changes
     - Pathway impact assessment

### Gene Mapping and Reference Data

Reference data files providing comprehensive annotations:

1. **Gene Mapping Full**:
   - **Path**: `reference/gene_mapping_full.tsv`
   - **Generator**: `scripts/general_gene_analysis/generate_gene_mapping_full.py`
   - **Contents**:
     - W303 gene ID to standard gene name mapping
     - Ergosterol pathway annotations
     - Gene function categories
     - Chromosome and position information
     - Conservation zone annotations

2. **Genes of Interest Mapping**:
   - **Path**: `reference/genes_of_interest_mapping.tsv`
   - **Contents**:
     - Ergosterol pathway genes
     - Functional categorization
     - Pathway position annotations
     - Distance relationships
     - Regulatory connections

3. **Chromosome Mapping**:
   - **Path**: `reference/chromosome_mapping.tsv`
   - **Contents**:
     - W303 scaffold IDs to standard chromosome names
     - Size and coordinate information
     - GC content and feature density
     - Annotation status
     - Reference cross-mapping

### Visualization Galleries

Organized collections of visualizations by analysis type:

1. **Mutation Spectrum Analysis**:
   - **Directory**: `analysis/mutation_spectrum_results/`
   - **Key Visualizations**:
     - `mutation_spectrum_summary.csv`: Numeric data on mutation types
     - `comparative_mutation_spectrum.png`: Cross-treatment comparison
     - `ti_tv_ratios_by_adaptation.png`: Transition/transversion by adaptation
     - Treatment-specific spectra (`CAS_mutation_spectrum.png`, etc.)

2. **Genomic Context Analysis**:
   - **Directory**: `analysis/genomic_context_results/`
   - **Key Visualizations**:
     - `gc_content_by_adaptation.png`: GC content patterns
     - `homopolymer_by_gene_status.png`: Homopolymer analysis
     - `mutation_type_by_treatment_heatmap.png`: Context heatmap
     - Sequence logos for specific contexts (e.g., `logo_Temperature_C_to_A.png`)

3. **Population Structure Analysis**:
   - **Directory**: `analysis/population_structure_results/`
   - **Key Visualizations**:
     - `pca_by_adaptation.png`: Principal component analysis
     - `mds_by_treatment.png`: Multi-dimensional scaling
     - `dendrogram_by_adaptation.png`: Hierarchical clustering
     - `shared_variants_heatmap.png`: Variant sharing patterns

4. **Regional Enrichment Analysis**:
   - **Directory**: `analysis/regional_enrichment_results/`
   - **Key Visualizations**:
     - `enrichment_heatmap.png`: Regional enrichment patterns
     - `clustered_enrichment_heatmap.png`: Clustered enrichment analysis
     - Treatment-specific enrichment results (CSV files)
     - Adaptation-specific enrichment summaries

5. **Scaffold Distribution Analysis**:
   - **Directory**: `analysis/scaffold_distribution_results/`
   - **Key Visualizations**:
     - `comparative_density_heatmap_top30.png`: Density comparison
     - `treatment_correlation_heatmap.png`: Treatment correlations
     - Treatment-specific bubble charts and variant density plots
     - `adaptation_specific_hotspots.txt`: Hotspot locations

6. **Gene Analysis Visualizations**:
   - **Directory**: `analysis/genes_of_interest/treatment_control_analysis/`
   - **Key Visualizations**:
     - `erg_gene_distribution.png`: Ergosterol gene variant distribution
     - `gene_status_distribution.png`: Gene status distribution
     - Treatment-specific purifying selection plots
     - `fold_change_by_gene_status.png`: Gene expression fold changes

7. **Sterol Profile Visualizations**:
   - **Directory**: `results/sterol_analysis/visualizations/`
   - **Key Visualizations**:
     - `sterol_composition_by_treatment.png`: Treatment comparison
     - `pathway_flux_visualization.png`: Pathway analysis
     - `satellite_gene_sterol_correlation.png`: Gene-sterol correlations
     - `four_zone_conservation_model.png`: Conservation architecture

These reports and visualizations collectively provide a comprehensive view of the project's findings, from basic statistics to complex models of adaptation mechanisms. The interactive HTML dashboards offer user-friendly exploration of the data, while the text and markdown reports provide detailed interpretations and biological significance.

All visualization outputs follow consistent color schemes and formatting to facilitate cross-analysis comparisons. The HTML reports are generated using reusable visualization templates that integrate Bootstrap for responsive design and D3.js for interactive elements.

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

### Running Sterol Analysis

To perform sterol profile analysis:

```bash
bash scripts/sterols/run_sterol_analysis.sh
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
9. Sterol data → Preprocessing → Standardized sterol profiles
10. Standardized sterol profiles → Comparative analysis → Treatment-specific profiles
11. Sterol profiles → Pathway analysis → Flux patterns and adaptation mechanisms
12. Genomic data + Sterol profiles → Integration → Comprehensive adaptation model

## References

- Yeast W303 strain reference genome and annotations
- Statistical methods for variant analysis
- Mutation signature analysis methodologies
- Sequence context analysis techniques
- Network analysis approaches for genomic data
- Functional impact prediction models
- Ergosterol pathway biochemistry and regulation
- Sterol biochemistry and yeast membrane adaptation
- Hierarchical conservation models in essential pathways
- Regulatory mechanisms in pathway adaptation
- Environmental stress response in unicellular organisms

---

*Note: This documentation provides a comprehensive overview of the Yeast MSA project, including its purpose, methodology, results, and future directions. For detailed results from each analysis module, please refer to the corresponding files in the analysis/ directory and results/sterol_analysis/ for sterol profile findings.*
