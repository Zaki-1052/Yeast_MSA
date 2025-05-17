# Analysis Modules

The Yeast MSA project employs several specialized analysis modules, each focusing on a specific aspect of the genomic data. This document describes each module's purpose, methodology, and outputs.

## Core Analysis Modules

### 1. Mutation Spectrum Analysis

- **Script**: `scripts/analysis/mutation_spectrum_analysis.py`
- **Purpose**: Examines patterns of single nucleotide mutations
- **Methodology**:
  - Classifies mutations by nucleotide change (e.g., C>A, G>T)
  - Calculates transition/transversion ratios
  - Performs statistical tests for pattern differences
  - Compares spectra across treatments and adaptations
- **Key Outputs**:
  - Treatment-specific mutation spectrum plots
  - Comparative spectrum visualization
  - Ti/Tv ratio analysis
  - Statistical significance tests
- **Results Location**: `analysis/mutation_spectrum_results/`

### 2. Genomic Context Analysis

- **Script**: `scripts/analysis/genomic_context_analysis.py`
- **Purpose**: Studies the sequence context around mutations
- **Methodology**:
  - Extracts sequences flanking mutation sites
  - Analyzes GC content in mutation regions
  - Identifies homopolymer and repeat sequences
  - Creates sequence logos for mutation contexts
- **Key Outputs**:
  - Sequence logos for mutation types
  - GC content analysis visualizations
  - Homopolymer association plots
  - Context heatmaps
- **Results Location**: `analysis/genomic_context_results/`

### 3. Mutational Signatures Analysis

- **Script**: `scripts/analysis/mutational_signature_analysis.py`
- **Purpose**: Identifies characteristic patterns of mutations
- **Methodology**:
  - Analyzes trinucleotide context of mutations
  - Identifies enriched contexts
  - Compares signatures across treatments
  - Identifies adaptation-specific signatures
- **Key Outputs**:
  - Signature plots by treatment
  - Enriched context visualization
  - Signature similarity heatmaps
  - Adaptation-specific signature profiles
- **Results Location**: `analysis/mutational_signatures_results/`

### 4. Population Structure Analysis

- **Script**: `scripts/analysis/population_spectrum_analysis.py`
- **Purpose**: Examines genetic relationships between samples
- **Methodology**:
  - Performs Principal Component Analysis (PCA)
  - Applies Multi-dimensional Scaling (MDS)
  - Constructs hierarchical clustering dendrograms
  - Calculates variant sharing statistics
- **Key Outputs**:
  - PCA and MDS plots
  - Dendrograms by adaptation and treatment
  - Variant sharing heatmaps
  - Genetic similarity metrics
- **Results Location**: `analysis/population_structure_results/`

### 5. Regional Enrichment Analysis

- **Script**: `scripts/analysis/regional_enrichment_analysis.py`
- **Purpose**: Identifies genomic regions with mutation hotspots
- **Methodology**:
  - Divides genome into windows
  - Calculates variant density per window
  - Identifies significantly enriched regions
  - Compares enrichment across treatments
- **Key Outputs**:
  - Enriched region tables for treatments
  - Enrichment heatmaps
  - Statistical significance tests
  - Shared enrichment analysis
- **Results Location**: `analysis/regional_enrichment_results/`

### 6. Scaffold Distribution Analysis

- **Script**: `scripts/analysis/scaffold_distribution_analysis.py`
- **Purpose**: Analyzes mutation distribution across scaffolds
- **Methodology**:
  - Maps variants to genomic scaffolds
  - Calculates variant density by scaffold
  - Compares distributions across treatments
  - Identifies scaffold-specific patterns
- **Key Outputs**:
  - Variant density plots
  - Correlation heatmaps
  - Bubble charts for visualization
  - Scaffold-specific statistics
- **Results Location**: `analysis/scaffold_distribution_results/`

### 7. Statistical Pattern Analysis

- **Script**: `scripts/analysis/statistical_pattern_analysis.py`
- **Purpose**: Performs statistical tests on mutation patterns
- **Methodology**:
  - Correlation analysis between treatments
  - Regression modeling
  - PCA for feature importance
  - Significance testing
- **Key Outputs**:
  - Correlation matrices
  - Regression plots
  - PCA biplots
  - Statistical test results
- **Results Location**: `analysis/statistical_pattern_results/`

## Gene Analysis Modules

### 1. Gene-Specific Analysis

- **Script Directory**: `scripts/gene_analysis/`
- **Run Script**: `run_gene_analysis.sh`
- **Purpose**: Analyzes mutations in ergosterol pathway genes
- **Methodology**:
  - Maps variants to ergosterol genes
  - Calculates gene-specific mutation rates
  - Analyzes variant effects in pathway genes
  - Compares mutation patterns across treatments
- **Key Outputs**:
  - Gene mutation profiles
  - Pathway gene distribution plots
  - Gene-specific context analysis
  - Pathway impact assessment
- **Results Location**: `analysis/genes_of_interest/`

### 2. General Gene Analysis

- **Script Directory**: `scripts/general_gene_analysis/`
- **Run Script**: `run_general_gene_analysis.sh`
- **Purpose**: Genome-wide gene-level analysis
- **Methodology**:
  - Comprehensive gene mapping
  - Gene-level mutation spectrum analysis
  - Gene-specific genomic context analysis
  - Gene-focused population structure
- **Key Outputs**:
  - Gene-centric mutation spectra
  - Gene-specific context analysis
  - Gene population structure results
  - Gene-level enrichment analysis
- **Results Location**: `analysis/general_gene_analysis/`

## Functional Impact Analysis

### 1. High Impact Variant Analysis

- **Script**: `scripts/functional_impact/analyze_high_impact_variants.py`
- **Run Script**: `run_high_impact_analysis.sh`
- **Purpose**: Identifies and analyzes functionally significant variants
- **Methodology**:
  - Identifies HIGH and MODERATE impact variants
  - Maps variants to protein domains
  - Predicts functional consequences
  - Analyzes treatment-specific patterns
- **Key Outputs**:
  - High impact variant tables
  - Protein domain impact maps
  - Functional consequence predictions
  - Treatment-specific impact profiles
- **Results Location**: `results/functional_impact/high_impact/`

### 2. Treatment-Specific Pattern Analysis

- **Script**: `scripts/functional_impact/analyze_treatment_specific_patterns.py`
- **Run Script**: `run_treatment_analysis.sh`
- **Purpose**: Identifies treatment-specific functional impacts
- **Methodology**:
  - Compares impact profiles across treatments
  - Identifies treatment-specific patterns
  - Statistical testing of pattern significance
  - Functional enrichment analysis
- **Key Outputs**:
  - Treatment-specific impact reports
  - Comparative treatment visualizations
  - Statistical significance tests
  - Functional pattern heatmaps
- **Results Location**: `results/functional_impact/treatment_patterns/`

### 3. Regulatory Element Analysis

- **Scripts**:
  - `scripts/functional_impact/analyze_promoter_elements.py`
  - `scripts/functional_impact/analyze_tfbs.py`
- **Purpose**: Characterizes variants affecting gene regulation
- **Methodology**:
  - Maps variants to promoter regions
  - Identifies affected transcription factor binding sites
  - Analyzes distribution patterns
  - Predicts regulatory impact
- **Key Outputs**:
  - Promoter variant distribution plots
  - TSS distance analysis
  - TFBS impact assessment
  - Regulatory pattern visualization
- **Results Location**: `results/functional_impact/regulatory/`

## Network Analysis

### Extended Ergosterol Network Analysis

- **Script**: `scripts/functional_impact/build_extended_erg_network.py`
- **Run Script**: `run_extended_network_analysis.sh`
- **Purpose**: Constructs and analyzes gene interaction networks
- **Methodology**:
  - Maps ergosterol pathway genes
  - Identifies affected neighboring genes
  - Calculates gene-gene relationships
  - Visualizes network topology
- **Key Outputs**:
  - Extended network visualization
  - Treatment-specific network maps
  - Gene subnetwork visualizations
  - Network statistics
- **Results Location**: `results/network_analysis/`

## Sterol Profile Analysis

### 1. Sterol Preprocessing

- **Script**: `scripts/sterols/preprocess_sterols.py`
- **Purpose**: Processes sterol concentration measurements
- **Methodology**:
  - Processes raw sterol data
  - Adds treatment metadata
  - Calculates relative abundances
  - Performs basic statistics
- **Key Outputs**:
  - Standardized sterol data
  - Basic sterol statistics
  - Concentration visualizations
  - Diversity metrics
- **Results Location**: `results/sterol_analysis/basic_stats/`

### 2. Comparative Sterol Analysis

- **Script**: `scripts/sterols/sterol_analysis.py`
- **Purpose**: Compares sterol profiles across conditions
- **Methodology**:
  - Statistical testing between treatments
  - Fold change calculations
  - Adaptation-specific sterol identification
  - Gene modification effects analysis
- **Key Outputs**:
  - Statistical test results
  - Fold change heatmaps
  - Treatment-specific profiles
  - Unique sterol identification
- **Results Location**: `results/sterol_analysis/comparative/`

### 3. Sterol Pathway Analysis

- **Script**: `scripts/sterols/sterol_pathway.py`
- **Purpose**: Connects sterol profiles to pathway steps
- **Methodology**:
  - Maps sterols to pathway enzymes
  - Calculates substrate/product ratios
  - Identifies altered pathway branches
  - Analyzes pathway flux
- **Key Outputs**:
  - Pathway mapping visualizations
  - Adaptation-specific pathway diagrams
  - Flux analysis by treatment
  - Enzyme activity inference
- **Results Location**: `results/sterol_analysis/pathway/`

### 4. Genomic-Sterol Integration

- **Script**: `scripts/sterols/sterol_integration.py`
- **Purpose**: Integrates sterol profiles with genomic patterns
- **Methodology**:
  - Correlates sterol changes with variants
  - Connects satellite genes to sterol production
  - Tests hierarchical conservation models
  - Builds integrated adaptation models
- **Key Outputs**:
  - Correlation analyses
  - Satellite gene-sterol maps
  - Integrated visualizations
  - Conservation architecture model
- **Results Location**: `results/sterol_analysis/correlation/`

## Integration and Reporting

### 1. Combined Analysis Results

- **Script**: `scripts/utils/combine_results.py`
- **Purpose**: Integrates findings across modules
- **Methodology**:
  - Collects key outputs from all modules
  - Summarizes significant findings
  - Creates cross-module correlation analysis
  - Generates comprehensive reports
- **Key Outputs**:
  - `combined_analysis_results.txt`
  - Cross-module correlation analyses
  - Integrated hypothesis validation
  - Meta-analysis of findings
- **Results Location**: `results/reports/`

### 2. Interactive HTML Reports

- **Scripts**:
  - `scripts/utils/generate_ergosterol_variant_report.py`
  - `scripts/utils/generate_functional_impact_report.py`
  - `scripts/sterols/generate_html_report.py`
  - `scripts/variants/variant_proximity_impact_summary.py`
  - `scripts/variants/generate_filtered_variants_visualizations.py`
- **Purpose**: Create interactive visualizations
- **Methodology**:
  - Integrates data into interactive format
  - Creates navigable dashboards
  - Implements visualization interactivity
  - Structures findings for exploration
- **Key Outputs**:
  - Ergosterol variant analysis dashboard
  - Functional impact analysis dashboard
  - Sterol profile analysis dashboard
  - Variant analysis dashboard
  - Variant proximity impact dashboard
  - Filtered variants visualization report
- **Results Location**: `results/reports/`

Each of these modules contributes to a comprehensive understanding of yeast adaptation mechanisms, from raw variant identification to biological interpretation and integration with biochemical data.