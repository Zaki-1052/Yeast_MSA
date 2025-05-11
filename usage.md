
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

### Script Execution Order

Below is the recommended order for executing scripts in this project, organized by module:

#### 1. Initial Data Processing (Phase 1)

These scripts handle raw sequencing data preparation and basic variant calling:

```
scripts/utils/run_fastqc.sh                  # Quality control for raw reads
scripts/utils/preprocessing.sh               # Clean and trim raw reads
scripts/utils/alignment.sh                   # Align reads to reference genome
scripts/utils/variants.sh                    # Call variants
scripts/utils/vcf_comparison.sh              # Initial variant comparison
scripts/utils/direct_comparison.sh           # Direct sample comparison
scripts/utils/consistency.sh                 # Check variant consistency
scripts/utils/position_analysis.sh           # Analyze variant positions
```

#### 2. Main Pipeline (Phase 2)

The complete analysis pipeline, which can be run all at once using the master script:

```
scripts/pipeline/12_run_pipeline.sh          # Master script that runs all pipeline steps
```

Or step by step:

```
scripts/pipeline/01_setup.sh                 # Setup directory structure
scripts/pipeline/02_filter_normalize.sh      # Filter and normalize VCFs
scripts/pipeline/03_fix_contigs.sh           # Fix contig naming
scripts/pipeline/04_merge_vcfs.sh            # Merge sample VCFs
scripts/pipeline/05_test_samples.sh          # Validate sample selection
scripts/pipeline/06_compare_groups.sh        # Compare treatment vs control
scripts/pipeline/07_cross_treatment.sh       # Compare between treatments
scripts/pipeline/08_identify_unique.sh       # Find unique variants
scripts/pipeline/09_direct_comparison.sh     # Pairwise comparisons
scripts/pipeline/10_consistency_analysis.sh  # Check consistency
scripts/pipeline/11_summary_report.sh        # Generate summary report
```

#### 3. Core Analysis Modules

These modules perform specific analyses and can be run individually or through the pipeline:

```
scripts/analysis/mutation_spectrum_analysis.py     # Analyze mutation patterns
scripts/analysis/genomic_context_analysis.py       # Study sequence context around mutations
scripts/analysis/mutational_signature_analysis.py  # Identify mutation signatures
scripts/analysis/population_spectrum_analysis.py   # Examine genetic relationships
scripts/analysis/regional_enrichment_analysis.py   # Find mutation hotspots
scripts/analysis/scaffold_distribution_analysis.py # Analyze mutation distribution
scripts/analysis/statistical_pattern_analysis.py   # Statistical testing
```

#### 4. Gene-Specific Analysis

These scripts focus on analyzing mutations in genes of interest:

```
scripts/gene_analysis/run_gene_analysis.sh        # Master script for gene analysis
```

#### 5. General Gene Analysis

These scripts provide comprehensive analysis across all genes:

```
scripts/general_gene_analysis/generate_gene_mapping_full.py    # Create gene mapping
scripts/general_gene_analysis/run_general_gene_analysis.sh     # Master script for general gene analysis
```

#### 6. Functional Impact Analysis

Analyzes functional impacts of variants:

```
scripts/functional_impact/run_high_impact_analysis.sh         # Analyze high impact variants
scripts/functional_impact/run_treatment_analysis.sh           # Treatment-specific patterns
scripts/functional_impact/run_promoter_analysis.sh            # Promoter region analysis
scripts/functional_impact/run_tfbs_analysis.sh                # TFBS analysis
```

#### 7. Network Analysis

Analyzes gene interaction networks:

```
scripts/functional_impact/run_extended_network_analysis.sh    # Build and analyze gene networks
```

#### 8. Sterol Profile Analysis

Analyzes sterol biochemical data:

```
scripts/sterols/run_sterol_analysis.sh                       # Complete sterol analysis
```

#### 9. Variant Extraction and Reporting

Extract and report on variants:

```
scripts/variants/run_extract_variants.sh                     # Extract variants
scripts/variants/verify_gene_coordinates.sh                  # Verify gene coordinates
scripts/utils/generate_ergosterol_variant_report.py          # Generate variant report
scripts/utils/generate_functional_impact_report.py           # Generate impact report
```

For simplicity, the master scripts for each major module can be run in this order:

1. `scripts/pipeline/12_run_pipeline.sh` (Main pipeline)
2. `scripts/gene_analysis/run_gene_analysis.sh` (Gene-specific analysis)
3. `scripts/general_gene_analysis/run_general_gene_analysis.sh` (General gene analysis)
4. `scripts/functional_impact/run_high_impact_analysis.sh` (Functional impact analysis)
5. `scripts/functional_impact/run_extended_network_analysis.sh` (Network analysis)
6. `scripts/sterols/run_sterol_analysis.sh` (Sterol analysis)
7. `scripts/variants/run_extract_variants.sh` (Variant extraction and reporting)

### Analysis Output Examination Guide

After running the pipeline, you'll want to examine the outputs in a logical order to understand the analysis results. Below is a recommended examination order for the key outputs, organized to progressively build understanding from raw data to integrated findings.

#### 1. Initial Quality Control & Processing Outputs

Start by reviewing the quality of the raw data and preprocessing:

```
results/preprocessing/fastqc/multiqc_report.html     # Overview of sequencing quality
results/stats/alignment_summary.txt                  # Summary of alignment statistics
results/vcf/filtered/*.stats.txt                     # Statistics for filtered VCF files
```

These outputs help verify data quality and provide confidence in downstream analyses.

#### 2. Variant Processing Results

Next, examine the main variant processing outputs to understand what variants were found:

```
results/merged/summary/analysis_report.txt          # Overview of variant analysis
results/merged/summary/treatment_comparison.tsv     # Comparison between treatments
results/gene_variants/all_gene_variants.tsv         # Variants mapped to genes
results/scaffold_variants/scaffold_variant_summary.txt # Variants by scaffold
```

These files show the basic variant identification results and their organization by genomic location.

#### 3. Mutation Spectrum Analysis

Start your biological analysis by looking at mutation patterns:

```
analysis/mutation_spectrum_results/                 # Directory with mutation spectrum outputs
analysis/mutation_spectrum_results/mutation_spectrum_summary.csv   # Summary statistics
analysis/mutation_spectrum_results/comparative_mutation_spectrum.png  # Visual comparison
analysis/mutation_spectrum_results/statistical_test_results.txt   # Statistical analysis
```

The mutation spectrum results show base substitution patterns that characterize each treatment condition.

#### 4. Genomic Context Analysis

Explore the sequence context around mutations:

```
analysis/genomic_context_results/genomic_context_summary.txt     # Summary of context findings
analysis/genomic_context_results/gc_content_by_adaptation.png    # GC content patterns
analysis/genomic_context_results/homopolymer_by_treatment.png    # Homopolymer analysis
analysis/genomic_context_results/mutation_type_by_treatment_heatmap.png  # Context heatmap
```

These outputs reveal sequence preferences and patterns around mutation sites.

#### 5. Mutational Signatures Analysis

Examine characteristic mutation patterns:

```
analysis/mutational_signatures_results/mutational_signatures_summary.txt  # Summary of signatures
analysis/mutational_signatures_results/signature_similarity_heatmap.png   # Signature comparisons
analysis/mutational_signatures_results/*_signature.png                    # Treatment-specific signatures
```

Mutation signatures provide insights into the underlying mutational processes at work.

#### 6. Population Structure Analysis

Look at how samples relate to each other genetically:

```
analysis/population_structure_results/population_structure_summary.txt  # Summary of findings
analysis/population_structure_results/pca_by_adaptation.png            # PCA visualization
analysis/population_structure_results/dendrogram_by_adaptation.png     # Hierarchical clustering
analysis/population_structure_results/shared_variants_heatmap.png      # Variant sharing patterns
```

Population structure analysis reveals genetic relationships between samples and adaptation types.

#### 7. Regional Enrichment Analysis

Identify hotspots of genetic variation:

```
analysis/regional_enrichment_results/regional_enrichment_summary.txt  # Summary of enriched regions
analysis/regional_enrichment_results/enrichment_heatmap.png          # Enrichment visualization
analysis/regional_enrichment_results/*_enriched_regions.csv          # Treatment-specific enrichments
```

These outputs show genomic regions with significant enrichment for mutations.

#### 8. Scaffold Distribution Analysis

Review how variants are distributed across the genome:

```
analysis/scaffold_distribution_results/scaffold_distribution_summary.txt  # Distribution summary
analysis/scaffold_distribution_results/comparative_density_heatmap_top30.png  # Density comparison
analysis/scaffold_distribution_results/treatment_correlation_heatmap.png  # Treatment correlations
```

This analysis shows mutation distribution patterns across scaffolds and treatments.

#### 9. Gene-Specific Analysis

Explore mutations in genes of interest:

```
analysis/gene_mutation_spectrum_results/gene_analysis_report.txt    # Gene mutation report
analysis/genes_of_interest/treatment_control_analysis/gene_specific_report.md  # Focused gene analysis
analysis/general_gene_analysis/gene_treatment_control_analysis/    # Treatment control comparisons
```

Gene-specific analyses focus on mutations in key genes, particularly those in the ergosterol pathway.

#### 10. Functional Impact Analysis

Review the predicted impacts of variants:

```
results/functional_impact/high_impact/high_impact_variants_report.md  # High impact variant analysis
results/functional_impact/variants_by_distance/variants_by_distance_report.md  # Distance analysis
results/functional_impact/key_genomic_regions/key_regions_report.md   # Key region analysis
```

These reports explain the predicted functional consequences of the observed genetic variations.

#### 11. Network Analysis

Explore gene interaction networks:

```
results/network_analysis/network_analysis_report.md        # Network analysis findings
results/network_analysis/extended_erg_network.png         # Complete network visualization
results/network_analysis/erg_subnetworks/                 # Individual pathway gene networks
```

Network analysis reveals relationships between ergosterol pathway genes and affected genes.

#### 12. Sterol Profile Analysis

Examine biochemical adaptations in sterol profiles:

```
results/sterol_analysis/basic_stats/summary_statistics.txt  # Basic sterol statistics
results/sterol_analysis/comparative/comparative_analysis_summary.txt  # Treatment comparisons
results/sterol_analysis/pathway/pathway_mapping.txt   # Sterol pathway analysis
results/sterol_analysis/visualizations/              # Sterol visualizations directory
```

Sterol profile analysis connects genomic findings to biochemical adaptations.

#### 13. Integrated Analysis

Finally, review the integrated findings:

```
results/sterol_analysis/correlation/integrated_findings_report.md  # Genomic-sterol integration
results/reports/combined_analysis_results.txt                     # Complete analysis summary
results/reports/ergosterol_variant_analysis.html                  # Interactive ergosterol report
results/reports/functional_impact.html                            # Interactive functional impact report
results/reports/sterols.html                                      # Interactive sterol report
```

These comprehensive reports integrate findings across all analysis modules, providing a complete picture of adaptation mechanisms.

Each of these outputs corresponds to specific analysis scripts discussed in the Script Execution Order section. The outputs are designed to build upon each other, progressively revealing more complex insights about the yeast adaptation process.

## Dependencies

### Python Packages

The analysis scripts require the following Python packages to be installed:

```bash
pip install numpy pandas matplotlib seaborn scipy biopython scikit-learn statsmodels networkx plotly logomaker
```

#### Core Dependencies
- **numpy**: Numerical computing and array operations
- **pandas**: Data manipulation and analysis
- **matplotlib**: Plotting and visualization
- **seaborn**: Statistical data visualization
- **scipy**: Scientific computing (statistical tests, clustering, etc.)

#### Bioinformatics Dependencies
- **biopython**: Biological sequence analysis and file parsing
- **logomaker**: Sequence logo generation for motif analysis

#### Machine Learning & Statistical Analysis
- **scikit-learn**: Machine learning algorithms for clustering and dimension reduction
- **statsmodels**: Statistical models and hypothesis testing

#### Visualization & Network Analysis
- **networkx**: Network creation and analysis
- **plotly**: Interactive visualizations

### System Dependencies

The pipeline also requires several command-line tools:

- **FastQC**: Quality control for sequencing data
- **BWA**: Burrows-Wheeler Aligner for read mapping
- **bcftools**: VCF file manipulation and variant calling
- **samtools**: SAM/BAM file manipulation
- **snpEff**: Variant annotation and effect prediction

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
