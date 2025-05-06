# Yeast MSA Gene Analysis Project Documentation

## Project Overview

This project involves modifying a series of Python scripts in the `scripts/gene_analysis` directory to incorporate gene-specific functionality. The original scripts in `scripts/analysis` perform various analyses on mutation data from yeast adaptation experiments, focusing on scaffold-level (genomic regions) analysis. The goal is to create enhanced versions of these scripts that map variants to specific genes, allowing for more biologically relevant analyses, with particular focus on detecting purifying selection in ergosterol pathway genes.

## Project Goals

1. Update all scripts in `scripts/gene_analysis` to incorporate gene-specific functionality
2. Implement gene mapping using genomic coordinates from reference annotations
3. Add visualization and analysis of variant patterns within genes
4. Create gene-specific summary reports for each type of analysis
5. Focus particularly on ergosterol pathway genes (genes of interest), which are under purifying selection
6. Detect and analyze both enrichment (positive selection) and depletion (purifying selection) of variants in genes
7. Statistically compare ergosterol pathway genes with other genes to assess selection pressure
8. Ensure consistent data structures and nomenclature across all scripts

## Repository Structure

- `/scripts/analysis/` - Original scaffold-level analysis scripts
- `/scripts/gene_analysis/` - New gene-level analysis scripts being modified
- `/reference/` - Contains reference genome and gene annotation files
  - `gene_mapping.tsv` - Maps gene IDs to genomic coordinates
  - `genes_of_interest_mapping.tsv` - List of ergosterol pathway genes
- `/analysis/` - Output directory for results
- `/vcf/` - Contains variant call format files with mutation data
- `/results/` - Contains intermediate and processed data

## Key Files

- `/reference/gene_mapping.tsv` - Essential file mapping W303 gene IDs to genomic coordinates
- `/reference/genes_of_interest_mapping.tsv` - Lists genes in the ergosterol pathway
- `/scripts/gene_analysis/*.py` - Scripts to be updated for gene-level analysis

## Data Structure

### Gene Data Format

The gene mapping file (`gene_mapping.tsv`) contains the following fields:
- `w303_gene_id` - Unique identifier for the gene in W303 yeast strain
- `sc_gene_id` - Standard S. cerevisiae gene ID
- `erg_name` - Ergosterol pathway gene name (if applicable)
- `locus_tag` - Locus tag identifier
- `w303_scaffold` - Scaffold/chromosome identifier
- `start` - Gene start position
- `end` - Gene end position
- `strand` - Strand direction (+ or -)
- `product` - Gene product description

### Variant Data Format

Variants in VCF files are identified by:
- Scaffold/chromosome
- Position
- Reference allele
- Alternate allele
- Format: `scaffold_position_ref_alt`

## Implementation Strategy

### Common Data Structures

1. `GENE_DATA` - Dictionary mapping gene IDs to their details
2. `SCAFFOLD_GENES` - Dictionary mapping scaffolds to lists of genes
3. `GENES_OF_INTEREST` - Set of gene IDs involved in the ergosterol pathway

### Common Functions

1. `load_gene_mapping()` - Loads gene data from reference files
2. `map_variants_to_genes()` - Maps variants to genes based on coordinates
3. `create_gene_specific_visualizations()` - Creates gene-specific plots
4. `create_gene_specific_summary()` - Generates a gene-specific report

### Script Modification Process

For each script in `scripts/gene_analysis/`:

1. Add gene data structures and common functions
2. Modify the input parsing to extract gene information
3. Update analysis functions to consider gene context
4. Add gene-specific visualization functions
5. Add gene-specific summary report generation
6. Update the main function to include gene-specific workflow
7. Ensure correct output paths for gene-specific results

## Progress Tracker

### Completed Scripts

- [x] `mutation_spectrum_analysis.py` - Added gene mapping, visualization, and reporting
- [x] `genomic_context_analysis.py` - Added gene mapping and context analysis
- [x] `mutational_signature_analysis.py` - Added gene-specific signature analysis
- [x] `population_spectrum_analysis.py` - Added gene mapping, visualization, and reporting
- [x] `regional_enrichment_analysis.py` - Added gene-specific enrichment/depletion analysis with focus on purifying selection
- [x] `scaffold_distribution_analysis.py` - Added gene distribution analysis with purifying selection detection
- [x] `statistical_pattern_analysis.py` - Added gene-level statistical analysis with purifying selection detection

### Remaining Scripts

- [x] `TC_Visualization.py` - Added gene-specific visualization
- [x] `variation.py` - Added gene-specific variation analysis with purifying selection detection

## Changes Made to Scripts

### Changes to `population_spectrum_analysis.py`

The following changes were made to add gene-specific functionality:

1. **Enhanced module docstring** to explain gene-level analysis
2. **Fixed code issues**:
   - Removed duplicate `ADAPTATION_COLORS` dictionary
3. **Added gene data structures**:
   - `GENE_DATA` dictionary for gene information
   - `SCAFFOLD_GENES` dictionary for scaffold-to-genes mapping
   - `GENES_OF_INTEREST` set for ergosterol pathway genes
4. **Added new output directory**:
   - `GENE_OUTPUT_DIR` for gene-specific analysis results
5. **Added new functions**:
   - `load_gene_mapping()` - Loads gene data from reference files
   - `map_variants_to_genes()` - Maps variants to genes based on coordinates
   - `create_gene_specific_visualizations()` - Creates gene-specific charts
   - `create_gene_specific_summary()` - Generates a gene-specific report
6. **Modified main function**:
   - Added gene mapping data loading
   - Added variant-to-gene mapping step
   - Added gene-specific visualization and reporting
7. **Added new gene-specific outputs**:
   - Distribution of variants by gene type
   - Variants in genes vs. non-genic regions by treatment
   - Distribution of variants across gene types by treatment
   - Gene-specific population structure summary report

### Changes to `regional_enrichment_analysis.py`

Extensive modifications were made to support gene-specific analysis with a focus on purifying selection:

1. **Enhanced docstring and imports**:
   - Added imports for statistical analysis (scipy.stats.ttest_ind)
   - Improved documentation to emphasize purifying selection analysis
2. **Added gene data structures**:
   - `GENE_DATA` dictionary for gene information
   - `SCAFFOLD_GENES` dictionary for scaffold-to-genes mapping
   - `GENES_OF_INTEREST` set for ergosterol pathway genes
   - `GENE_OUTPUT_DIR` for gene-specific analysis results
3. **Added new functions**:
   - `load_gene_mapping()` - Loads gene data from reference files
   - `map_variants_to_genes()` - Maps variants to genes based on coordinates
   - `analyze_gene_specific_enrichment()` - Analyzes both enrichment and depletion patterns
   - `plot_gene_enrichment()` - Creates visualizations showing both enrichment and depletion
   - `create_gene_specific_reports()` - Generates detailed gene-specific reports
   - `create_gene_specific_summary()` - Creates comprehensive gene-level summary
4. **Implemented purifying selection detection**:
   - Added log2 fold change calculation for better visualization of depletion
   - Modified p-value calculations to detect statistically significant depletion
   - Added statistical comparison between ergosterol and non-ergosterol genes
   - Created sophisticated visualizations highlighting genes under purifying selection
5. **Enhanced visualization capabilities**:
   - Added boxplots and violin plots to compare gene groups
   - Created bar plots for top enriched and depleted genes
   - Added statistical annotations to visualizations
   - Created specific visualizations for ergosterol pathway genes
6. **Comprehensive reporting**:
   - Added detailed stats on ergosterol pathway gene purifying selection
   - Added statistical comparison with non-ergosterol genes
   - Created specific files for genes under purifying selection
   - Enhanced summary reports with interpretation of selection patterns
7. **Modified main function**:
   - Added gene mapping and analysis workflow
   - Integrated gene-specific analysis with regional enrichment analysis
   - Added comprehensive statistical testing and reporting

### Changes to `scaffold_distribution_analysis.py`

Significant enhancements were made to add gene-level analysis to the scaffold distribution script:

1. **Added gene-related imports and data structures**:
   - Added statistical testing imports (scipy.stats.ttest_ind, mannwhitneyu)
   - Added new gene-related color schemes and visualization parameters
   - Added gene status constants (ERG, Non-ERG, No Gene)
2. **Implemented gene mapping functionality**:
   - Added `load_gene_mapping()` to parse gene annotation files
   - Added `build_gene_coordinates()` to create scaffold-to-genes mapping
   - Added `map_variants_to_genes()` to assign variants to genes
3. **Added gene-specific analysis functions**:
   - `count_variants_by_gene_status()` - Counts variants in ERG, Non-ERG, and intergenic regions
   - `count_variants_per_gene()` - Counts variants within each specific gene
   - `calculate_gene_density()` - Calculates variant density per kilobase of gene length
   - `identify_gene_patterns()` - Identifies statistically enriched or depleted genes
4. **Implemented purifying selection detection**:
   - Added bidirectional p-value calculations to detect both enrichment and depletion
   - Added log2 fold change calculation for visualizing selection direction
   - Implemented separate analysis for enriched vs. depleted genes
   - Added specific detection for purifying selection in ergosterol genes
5. **Enhanced visualization capabilities**:
   - Added `plot_gene_enrichment()` - Creates volcano plots showing enrichment and depletion
   - Added `plot_gene_status_distribution()` - Shows variant distribution by gene status
   - Added statistical visualizations comparing ERG vs. Non-ERG genes
   - Created specific visualizations for ergosterol pathway genes
6. **Comprehensive reporting**:
   - Enhanced `create_summary_table()` with gene-specific statistics
   - Updated `create_summary_report()` with detailed gene analysis
   - Added specific sections on purifying selection in ergosterol pathway
   - Added statistical comparisons between gene categories
7. **Updated main workflow**:
   - Integrated gene mapping into the main analysis pipeline
   - Added gene-specific analysis steps in parallel with scaffold analysis
   - Added comprehensive gene-specific visualization generation
   - Enhanced summary statistics with gene-level metrics

### Changes to `statistical_pattern_analysis.py`

Comprehensive enhancements were made to add gene-level statistical pattern analysis capabilities:

1. **Added gene-related imports and data structures**:
   - Added statistical testing imports (scipy.stats.ttest_ind, mannwhitneyu, poisson)
   - Added statistical correction tools (statsmodels.stats.multitest)
   - Added new gene-related color schemes and visualization parameters
   - Added gene status constants (ERG, Non-ERG, No Gene)
   - Added output directory for gene-specific results

2. **Implemented gene mapping functionality**:
   - Added `load_gene_mapping()` to parse gene annotation files
   - Added `map_variants_to_genes()` to assign variants to genes by coordinates
   - Added gene-specific metadata to variants (gene_id, gene_name, gene_type)

3. **Enhanced data integration**:
   - Updated `integrate_data()` to include gene mapping
   - Added gene-specific statistics calculation (variant counts, densities)
   - Added separate gene-specific output directory (GENE_OUTPUT_DIR)
   - Added statistical comparisons between gene categories

4. **Added purifying selection analysis**:
   - Added `analyze_gene_patterns()` for comprehensive gene-specific statistical analysis
   - Implemented bidirectional statistical testing to detect both enrichment and depletion
   - Added log2 fold change calculation for visualizing selection patterns
   - Implemented specific detection of purifying selection in ergosterol genes
   - Added multiple testing correction with False Discovery Rate (FDR) control

5. **Enhanced statistical capabilities**:
   - Added statistical comparison between ERG and non-ERG genes
   - Added both parametric (t-test) and non-parametric (Mann-Whitney U) testing
   - Added genome-wide mutation rate calculation for expected variant counts
   - Added poisson distribution modeling for statistical significance

6. **Advanced gene-specific visualizations**:
   - Added `create_gene_specific_visualizations()` for comprehensive gene-specific visualization
   - Implemented visualizations of variant distribution by gene status
   - Added treatment-specific gene status visualizations
   - Added log2 fold change visualization for enrichment/depletion
   - Created visualizations of top enriched and depleted genes with gene status coloring

7. **Comprehensive reporting**:
   - Added `create_gene_specific_report()` for detailed gene-specific analysis reports
   - Added interpretation of purifying selection evidence
   - Added statistical summaries comparing ERG vs. non-ERG genes
   - Added listing of top enriched and depleted genes with significance values
   - Included comprehensive biological interpretation of gene-level patterns

### Changes to `TC_Visualization.py`

Comprehensive enhancements were made to add gene-specific visualization functionality:

1. **Enhanced module docstring and imports**:
   - Updated docstring to explain gene-level visualization functionality
   - Added imports for file handling and data structures (os, defaultdict, logging)
   - Added logging configuration to track script execution

2. **Added gene data structures and constants**:
   - Added `GENE_DATA` dictionary for gene information
   - Added `SCAFFOLD_GENES` dictionary for scaffold-to-genes mapping
   - Added `GENES_OF_INTEREST` set for ergosterol pathway genes
   - Added `GENE_COLORS` dictionary for consistent gene status visualization
   - Added output directory for gene-specific visualizations (GENE_OUTPUT_DIR)

3. **Implemented gene mapping functionality**:
   - Added `load_gene_mapping()` to parse gene annotation files
   - Added `map_variants_to_genes()` to assign variants to genes by coordinates
   - Added variant extraction from different file formats with ID parsing

4. **Added gene-specific visualization functions**:
   - Added `create_gene_specific_visualizations()` for comprehensive gene-level visualizations
   - Created variant distribution visualization by gene status (ERG, Non-ERG, No Gene)
   - Added proportion visualization of gene status variants across treatments
   - Added ergosterol pathway gene-specific variant distribution
   - Created fold change visualization by gene status with bubble charts

5. **Added comprehensive reporting**:
   - Added `create_gene_specific_report()` for detailed gene-specific analysis reporting
   - Implemented gene status statistics by treatment
   - Added ergosterol gene-specific variant statistics
   - Added purifying selection analysis with expected vs. observed proportions
   - Created fold difference calculation and interpretation for selection pressure

6. **Enhanced data loading and integration**:
   - Added `load_gene_level_variant_data()` to load and process variant data
   - Implemented flexible file path checking to find available data
   - Added fallback options when gene mapping data is missing
   - Added treatment extraction from sample IDs for better categorization

7. **Updated main workflow**:
   - Modified main function to incorporate gene-specific analysis
   - Added standard visualization generation with the original function
   - Added gene-specific visualization when data is available
   - Implemented error handling and user feedback for missing data

### Changes to `variation.py`

Extensive modifications were made to add gene-specific analysis capabilities with a focus on purifying selection:

1. **Enhanced module docstring and imports**:
   - Updated docstring to explain gene-level analysis functionality
   - Added imports for statistical analysis (ttest_ind, mannwhitneyu, poisson)
   - Added imports for visualization (matplotlib, seaborn)
   - Added logging configuration for better tracking

2. **Added gene data structures and constants**:
   - Added `GENE_DATA` dictionary for gene information
   - Added `SCAFFOLD_GENES` dictionary for scaffold-to-genes mapping
   - Added `GENES_OF_INTEREST` set for ergosterol pathway genes
   - Added `GENE_COLORS` dictionary for consistent gene status visualization
   - Added output directory for gene-specific results (GENE_OUTPUT_DIR)

3. **Implemented gene mapping functionality**:
   - Added `load_gene_mapping()` to parse gene annotation files
   - Added `extract_variants_from_vcf()` to extract and map variants to genes
   - Added comprehensive validation and error handling

4. **Added gene-specific analysis functions**:
   - Added `analyze_gene_specific_patterns()` for detailed gene-level analysis
   - Implemented variant counting by gene status (ERG, Non-ERG, No Gene)
   - Added fold change calculation for each gene type
   - Implemented statistical testing for enrichment/depletion using Poisson distribution
   - Added special detection of purifying selection (observed < expected)
   - Added multiple testing correction with FDR control

5. **Enhanced statistical capabilities**:
   - Added gene-specific statistical comparison between ERG and non-ERG genes
   - Implemented both parametric (t-test) and non-parametric (Mann-Whitney U) tests
   - Added log2 fold change calculation for visualizing selection direction
   - Added differential fold change calculation between treatment and control

6. **Added comprehensive visualization functions**:
   - Added `create_gene_specific_visualizations()` for detailed gene-level visualizations
   - Implemented gene status distribution bar charts
   - Added stacked proportion charts of gene status variants
   - Created ergosterol pathway gene distribution visualizations
   - Added fold change comparison by gene category
   - Created volcano plots showing gene enrichment/depletion patterns

7. **Comprehensive reporting**:
   - Added `create_gene_specific_report()` for detailed gene-specific analysis reporting
   - Created treatment summary with gene status breakdown
   - Added purifying selection interpretation by treatment
   - Added statistical comparison of ERG vs. non-ERG genes
   - Created tables of significantly enriched/depleted genes
   - Added biological interpretation of selection patterns

## Technical Challenges and Solutions

### Challenge 1: Mapping Variants to Genes

**Problem**: Variants are identified by scaffold and position, which needs to be mapped to gene coordinates.

**Solution**: 
- Created the `map_variants_to_genes()` function that:
  - Extracts scaffold and position from variant IDs
  - Checks if position falls within any gene's start-end range
  - Adds gene metadata to variant data
  - Identifies if variant is in a gene of interest

### Challenge 2: Detecting Purifying Selection

**Problem**: Ergosterol pathway genes are under purifying selection, resulting in fewer variants than expected. This requires detection of both enrichment and depletion patterns.

**Solution**:
- Modified enrichment analysis to detect both positive and negative selection:
  - Values > 1 in fold_enrichment indicate positive selection
  - Values < 1 indicate purifying selection
  - Added log2 fold change calculation for better visualization of depletion
  - Implemented statistical comparison between ergosterol and non-ergosterol genes
  - Created visualizations that highlight both enrichment and depletion

### Challenge 2: Handling Duplicate Code

**Problem**: The gene_analysis scripts had duplicate adaptation color definitions.

**Solution**:
- Removed duplicate `ADAPTATION_COLORS` dictionary
- Ensured consistent use of global variables

### Challenge 3: Making Analysis Gene-Aware

**Problem**: Original analysis functions only considered scaffold-level patterns.

**Solution**:
- Added 'in_gene', 'gene_id', and 'gene_type' fields to variant data
- Created gene-specific groupings in analysis functions
- Added gene-specific visualization functions
- Implemented statistical tests to compare gene groups
- Created visualizations that show both enrichment and depletion patterns within genes

## Strategy for Remaining Scripts

### For `variation.py`:

1. Add gene mapping functionality
2. Add gene-specific variation analysis
3. Create gene-level variation summaries
4. Implement detection of purifying selection
5. Create gene-specific visualization and reporting

## Purifying Selection Analysis for Ergosterol Pathway Genes

A key focus of the gene-specific analysis is detecting and analyzing purifying selection in ergosterol pathway genes. This required several specialized approaches implemented in `regional_enrichment_analysis.py`:

### Statistical Methods for Detecting Purifying Selection

1. **Enrichment/Depletion Calculation**:
   - For each gene, we calculate expected variants based on gene length and genome-wide mutation rate
   - Fold enrichment = observed/expected variants
   - Values < 1 indicate purifying selection (fewer mutations than expected)
   - Log2 transformation of fold enrichment to better visualize both enrichment and depletion

2. **Statistical Testing**:
   - For enrichment (variant_count > expected): P-value calculated as 1 - Poisson CDF(variant_count - 1, expected)
   - For depletion (variant_count < expected): P-value calculated as Poisson CDF(variant_count, expected)
   - Multiple testing correction using False Discovery Rate (FDR) method
   - Comparison between ergosterol and non-ergosterol genes using t-tests

### Visualizations for Purifying Selection

1. **Gene Category Comparison**:
   - Box plots comparing log2 fold change between ergosterol and non-ergosterol genes
   - Violin plots showing distribution of selection pressure by gene category
   - Statistical annotations showing significance of differences

2. **Specific Gene Visualizations**:
   - Horizontal bar charts of ergosterol genes showing both enrichment and depletion
   - Color-coded visualizations distinguishing enriched from depleted genes
   - Clear highlighting of genes under strongest purifying selection

3. **Comparative Analysis**:
   - Side-by-side visualization of enriched vs. depleted genes
   - Statistical comparison with non-ergosterol genes
   - Treatment-specific patterns of selection pressure

### Comprehensive Reporting

1. **Gene-Specific Reports**:
   - Detailed analysis of ergosterol pathway genes with selection statistics
   - Lists of genes under significant purifying selection
   - Statistical comparison with background genomic patterns

2. **Interpretive Summaries**:
   - Automated interpretation of selection patterns
   - Assessment of heterogeneity within the ergosterol pathway
   - Detailed statistics on the strength of purifying selection
   - Comparison across treatments and adaptations

## Tips for Implementation

1. **Start with common functions**: Reuse the `load_gene_mapping()` and `map_variants_to_genes()` functions across all scripts
2. **Maintain consistency**: Use the same data structures and naming conventions across all scripts
3. **Test incrementally**: Test each script after modification before moving to the next one
4. **Focus on genes of interest**: Ensure special analysis for ergosterol pathway genes
5. **Create meaningful visualizations**: Make gene-specific visualizations informative and focused
6. **Generate comprehensive reports**: Each script should produce a detailed gene-specific summary

## Next Steps

1. Test all scripts together on a sample dataset
2. Create comprehensive documentation for all functions and scripts
3. Create a final report summarizing all gene-specific analysis findings
4. Generate integrated visualizations showing gene-level patterns across all analyses