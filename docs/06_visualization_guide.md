# Visualization Guide

This document provides an overview of the various visualizations and reports available in the Yeast MSA project, explaining how to interpret the key visualizations and where to find them.

## Interactive HTML Reports

The project includes several interactive HTML dashboards that provide rich visualizations and explorable data:

### 1. Ergosterol Variant Analysis Dashboard

- **Path**: `results/reports/ergosterol_variant_analysis.html`
- **Generator**: `scripts/utils/generate_ergosterol_variant_report.py`
- **Contents**:
  - Interactive variant distribution visualizations
  - Distance-based analysis from pathway genes
  - Treatment-specific variant patterns
  - Purifying selection evidence
  - Biological significance interpretation
- **How to Use**:
  - Navigate through tabs to explore different analyses
  - Use the zoom functionality for detailed inspection
  - Toggle between dark/light mode for viewing preference
  - View the image gallery of key visualizations

### 2. Functional Impact Analysis Dashboard

- **Path**: `results/reports/functional_impact.html`
- **Generator**: `scripts/utils/generate_functional_impact_report.py`
- **Contents**:
  - Protein domain impact visualizations
  - Conservation patterns analysis
  - Satellite gene architecture diagrams
  - Adaptation mechanisms models
  - Hierarchical conservation zone visualizations
- **How to Use**:
  - Use the navigation sidebar to explore different sections
  - Hover over visualizations for detailed information
  - Click on genes of interest for detailed impact information
  - Use the search functionality to find specific features

### 3. Sterol Profile Analysis Dashboard

- **Path**: `results/reports/sterols.html`
- **Generator**: `scripts/sterols/generate_html_report.py`
- **Contents**:
  - Interactive sterol profile visualizations
  - Treatment comparison heatmaps
  - Pathway visualizations with flux indicators
  - Genomic-sterol integration diagrams
  - Adaptation model visualizations
- **How to Use**:
  - Select treatments of interest from the dropdown menus
  - Toggle between different visualization types
  - Click on pathway steps to see detailed information
  - Use the comparison tool to directly compare treatments

### 4. Variant Analysis Dashboard

- **Path**: `results/reports/variant_analysis.html`
- **Generator**: `scripts/variants/generate_variant_report.py`
- **Contents**:
  - Sample comparison visualizations
  - Interactive filtering of variants
  - Annotation statistics and distribution
  - Genome browser-like variant visualization
  - Treatment-specific variant patterns
- **How to Use**:
  - Filter variants by impact, type, or location
  - Navigate through the genome using the browser interface
  - Compare samples side-by-side
  - Export filtered variant lists

## Key Visualization Galleries

Organized collections of visualizations by analysis type are available in the analysis directory:

### 1. Mutation Spectrum Analysis

- **Directory**: `analysis/mutation_spectrum_results/`
- **Key Visualizations**:
  - `comparative_mutation_spectrum.png`: Shows mutation type distributions across all treatments in a single plot, enabling direct comparison
  - `ti_tv_ratios_by_adaptation.png`: Displays transition/transversion ratios across adaptation types
  - `mutation_by_adaptation.png`: Shows mutation types categorized by adaptation (temperature vs. low oxygen)
  - Treatment-specific spectra (e.g., `CAS_mutation_spectrum.png`): Individual mutation profiles for each treatment

### 2. Genomic Context Analysis

- **Directory**: `analysis/genomic_context_results/`
- **Key Visualizations**:
  - `gc_content_by_adaptation.png`: Shows GC content distribution by adaptation type
  - `homopolymer_by_gene_status.png`: Illustrates the relationship between homopolymer regions and gene status
  - `mutation_type_by_treatment_heatmap.png`: Heatmap of mutation types across treatments
  - Sequence logos (e.g., `logo_Temperature_C_to_A.png`): Visual representations of sequence context for specific mutation types

### 3. Population Structure Analysis

- **Directory**: `analysis/population_structure_results/`
- **Key Visualizations**:
  - `pca_by_adaptation.png`: Principal component analysis plot showing sample clustering by adaptation type
  - `mds_by_treatment.png`: Multi-dimensional scaling plot showing genetic relationships between treatments
  - `dendrogram_by_adaptation.png`: Hierarchical clustering dendrogram organized by adaptation type
  - `shared_variants_heatmap.png`: Heatmap showing variant sharing patterns between samples

### 4. Regional Enrichment Analysis

- **Directory**: `analysis/regional_enrichment_results/`
- **Key Visualizations**:
  - `enrichment_heatmap.png`: Heatmap showing regional enrichment patterns across the genome
  - `clustered_enrichment_heatmap.png`: Clustered heatmap to identify similar enrichment patterns
  - `enrichment_distribution.png`: Distribution plot of enrichment scores across the genome
  - Treatment-specific enrichment visualizations and CSV files

### 5. Scaffold Distribution Analysis

- **Directory**: `analysis/scaffold_distribution_results/`
- **Key Visualizations**:
  - `comparative_density_heatmap_top30.png`: Heatmap comparing variant density across top 30 scaffolds
  - `treatment_correlation_heatmap.png`: Correlation heatmap between treatments
  - `adaptation_correlation_heatmap.png`: Correlation heatmap between adaptation types
  - Treatment-specific bubble charts (e.g., `CAS_bubble_chart.png`): Visualize variant distribution by scaffold size

### 6. Gene Analysis Visualizations

- **Directory**: `analysis/genes_of_interest/treatment_control_analysis/`
- **Key Visualizations**:
  - `erg_gene_distribution.png`: Distribution of variants across ergosterol genes
  - `gene_status_distribution.png`: Variant distribution by gene status
  - `fold_change_by_gene_status.png`: Fold changes in variant frequency by gene status
  - Treatment-specific purifying selection plots (e.g., `CAS_purifying_selection.png`)

### 7. Sterol Profile Visualizations

- **Directory**: `results/sterol_analysis/visualizations/`
- **Key Visualizations**:
  - `sterol_composition_by_treatment.png`: Sterol composition across treatments
  - `pathway_flux_visualization.png`: Visualization of flux through ergosterol pathway
  - `satellite_gene_sterol_correlation.png`: Correlation between satellite genes and sterol profiles
  - `four_zone_conservation_model.png`: Visualization of the hierarchical conservation model

## Specialized Visualization Types

### Sequence Logos

Sequence logos (in `analysis/genomic_context_results/` and `analysis/mutational_signatures_results/`) visually represent sequence context around mutations:

- **Height of letters**: Represents the frequency of each nucleotide at that position
- **Color coding**: A=green, C=blue, G=yellow, T=red
- **Position numbering**: Center (0) is the mutation site, with flanking nucleotides numbered relative to it
- **Example**: `logo_Temperature_C_to_A.png` shows sequence context for C>A mutations in temperature adaptation

### Heatmaps

Various heatmaps are used throughout the project to visualize complex relationships:

- **Enrichment heatmaps**: Color intensity represents degree of enrichment
- **Correlation heatmaps**: Colors represent correlation strength (typically blue=negative, red=positive)
- **Clustered heatmaps**: Include dendrograms showing hierarchical relationships
- **Treatment heatmaps**: Compare patterns across treatment conditions

### Network Visualizations

Found in `results/network_analysis/`:

- **Node size**: Represents the number of variants or connectivity
- **Node color**: Indicates gene function or pathway membership
- **Edge thickness**: Represents strength of relationship
- **Edge color**: Indicates type of interaction
- **Network layouts**: Typically force-directed to emphasize clusters

### Principal Component Analysis (PCA)

PCA plots in `analysis/population_structure_results/` and other directories:

- **Axes**: Represent principal components explaining most variation
- **Point position**: Represents genetic similarity (closer points = more similar)
- **Point shape/color**: Typically indicates treatment or adaptation type
- **Variance explained**: Often shown as percentages on axes

## How to Navigate Results

### For a High-Level Overview

1. Start with the interactive HTML reports in `results/reports/`
2. Review the `combined_analysis_results.txt` in the same directory
3. Examine the key visualization galleries in the analysis directory

### For Detailed Gene Analysis

1. Check `analysis/genes_of_interest/` for ergosterol pathway gene analysis
2. Examine `results/gene_variants/` for detailed variant information by gene
3. Look at `results/network_analysis/` for gene interaction networks

### For Sterol-Genomic Integration

1. Start with the sterol HTML report in `results/reports/sterols.html`
2. Look at `results/sterol_analysis/correlation/` for genomic-sterol correlations
3. Check `results/sterol_analysis/pathway/` for pathway analysis

### For Treatment Comparisons

1. Examine `analysis/population_structure_results/` for genetic relationships
2. Look at `analysis/treatment_control_analysis/` for direct treatment comparisons
3. Check treatment-specific visualizations throughout the analysis directories

## Tips for Interpretation

1. **Color consistency**: Throughout the project, consistent color schemes are used:
   - Treatment colors: WT-37=blue, WTA=green, STC=orange, CAS=red
   - Adaptation colors: Temperature=red, Low Oxygen=blue
   - Gene modification: Modified=orange, Non-modified=green

2. **Naming conventions**: File names typically include the treatment or condition they represent:
   - Treatment-specific files start with the treatment name (e.g., `CAS_mutation_spectrum.png`)
   - Adaptation-specific files include "adaptation" (e.g., `adaptation_correlation_heatmap.png`)
   - Comparison files start with "comparative" or include "vs"

3. **Statistical significance**: Look for accompanying text files with statistical test results:
   - `statistical_test_results.txt` in analysis directories
   - p-values are typically indicated in visualizations when significant
   - Asterisks often denote significance levels (* p<0.05, ** p<0.01, *** p<0.001)

4. **Interactive features**: The HTML reports include various interactive features:
   - Hover for additional information
   - Click for detailed views
   - Toggle between different visualization modes
   - Filter data to focus on specific features

This visualization guide should help navigate the rich set of visualizations available in the Yeast MSA project, enabling a comprehensive understanding of the findings.