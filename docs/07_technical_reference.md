# Technical Reference

This document provides technical details on the tools, methods, and parameters used in the Yeast MSA project.

## Reference Genome

- **Source**: W303 yeast strain reference genome
- **Assembly**: GCA_000766475.2
- **Location**: `reference/w303_chromosomal.fasta`
- **Annotation**: Gene annotations from WGAP (Wolfe lab yeast annotation pipeline)
- **Chromosomes**: 16 core chromosomes plus mitochondrial genome
- **Total size**: ~12 Mbp

## Core Bioinformatics Tools

### Quality Control and Preprocessing

#### FastQC
- **Version**: FastQC v0.11.9
- **Usage**: Quality assessment of raw reads
- **Parameters**: Default settings
- **Script**: `scripts/utils/run_fastqc.sh`

#### fastp
- **Version**: fastp v0.20.1
- **Usage**: Read trimming and filtering
- **Parameters**:
  - Quality filter: `-q 20` (Phred quality ≥ 20)
  - Length filter: `-l 50` (Minimum read length 50bp)
  - Adapter trimming: `--detect_adapter_for_pe` (Auto-detect adapters)
- **Script**: `scripts/utils/preprocessing.sh`

### Alignment and Variant Calling

#### BWA-MEM
- **Version**: BWA v0.7.17
- **Usage**: Align reads to reference genome
- **Parameters**:
  - Algorithm: `mem` (BWA-MEM algorithm)
  - Threads: `-t 8` (8 threads)
  - Mark shorter split hits as secondary: `-M`
- **Script**: `scripts/utils/alignment.sh`

#### samtools
- **Version**: samtools v1.13
- **Usage**: SAM/BAM file manipulation
- **Key Operations**:
  - Sort: `sort -n` (Name sort for preprocessing)
  - Sort: `sort -o` (Coordinate sort for final BAM)
  - Index: `index` (Create BAI index)
  - Stats: `flagstat` (Alignment statistics)
  - Coverage: `depth` (Coverage statistics)
- **Script**: `scripts/utils/alignment.sh`

#### bcftools
- **Version**: bcftools v1.13
- **Usage**: Variant calling and VCF manipulation
- **Variant Calling Parameters**:
  - mpileup: `-f reference.fasta -Ou` (Generate genotype likelihoods)
  - call: `-mv -Ov` (Multi-allelic caller with variant-only output)
  - ploidy: `--ploidy 1` (Haploid model for yeast)
- **Filtering Parameters**:
  - Quality: `QUAL>=20` (Minimum quality score 20)
  - Depth: `DP>=10` (Minimum read depth 10)
- **Scripts**: 
  - `scripts/utils/variants.sh`
  - `scripts/pipeline/02_filter_normalize.sh`

#### snpEff
- **Version**: snpEff v5.0
- **Usage**: Variant annotation and effect prediction
- **Database**: Custom built from WGAP annotations
- **Parameters**:
  - Genome: `Saccharomyces_cerevisiae_W303`
  - Annotation: `-classic` (Classic annotation style)
  - Stats: `-stats` (Generate summary statistics)
- **Annotation Categories**:
  - Impact: HIGH, MODERATE, LOW, MODIFIER
  - Effect: missense_variant, synonymous_variant, stop_gained, etc.
  - Region: exon, intron, intergenic, upstream, downstream, etc.
- **Script**: Located in annotation pipeline

## Analysis Methods

### Mutation Spectrum Analysis

- **Method**: Count and categorize single nucleotide substitutions
- **Categories**: A>C, A>G, A>T, C>A, C>G, C>T, G>A, G>C, G>T, T>A, T>C, T>G
- **Normalization**: Counts normalized to total mutations per sample
- **Statistical Tests**: Chi-square test for independence between treatments
- **Script**: `scripts/analysis/mutation_spectrum_analysis.py`

### Genomic Context Analysis

- **Method**: Extract and analyze sequence context around mutation sites
- **Window Size**: ±5bp from mutation site
- **Features Analyzed**:
  - GC content in flanking regions
  - Homopolymer presence (runs of 3+ identical nucleotides)
  - Dinucleotide repeat patterns
  - Trinucleotide context
- **Visualization**: Sequence logos created with Logomaker Python package
- **Script**: `scripts/analysis/genomic_context_analysis.py`

### Mutational Signatures Analysis

- **Method**: Analyze trinucleotide context of mutations
- **Context Definition**: One base upstream + mutated base + one base downstream
- **Categories**: 96 possible trinucleotide contexts (4×6×4)
- **Normalization**: Normalized to mutation count and adjusted for trinucleotide frequency
- **Similarity Metric**: Cosine similarity between signature vectors
- **Script**: `scripts/analysis/mutational_signature_analysis.py`

### Population Structure Analysis

- **Method**: Compare genetic relationships between samples
- **Techniques**:
  - Principal Component Analysis (PCA)
  - Multi-dimensional Scaling (MDS)
  - Hierarchical clustering (UPGMA method)
  - Jaccard similarity index
- **Input**: Binary variant presence/absence matrix
- **Software**: scikit-learn for PCA and clustering
- **Script**: `scripts/analysis/population_spectrum_analysis.py`

### Regional Enrichment Analysis

- **Method**: Identify genomic regions with mutation hotspots
- **Window Size**: 5kb sliding windows with 1kb step
- **Enrichment Score**: Observed/expected variant count ratio
- **Significance**: Binomial test with multiple testing correction
- **Threshold**: FDR-corrected p < 0.05 and fold enrichment > 2
- **Script**: `scripts/analysis/regional_enrichment_analysis.py`

### Scaffold Distribution Analysis

- **Method**: Analyze mutation distribution across genome scaffolds
- **Normalization**: Variants per kilobase of scaffold
- **Correlation**: Spearman rank correlation between treatments
- **Visualization**: Bubble charts, heatmaps, density plots
- **Script**: `scripts/analysis/scaffold_distribution_analysis.py`

### Gene Analysis

- **Method**: Map variants to genes and analyze distribution
- **Gene Annotation**: WGAP gene annotations + custom ergosterol pathway annotation
- **Distance Categories**: within gene, 0-1kb, 1-5kb, 5-50kb, >50kb
- **Selection Analysis**: Observed/expected variant ratio in genes vs. non-genic regions
- **Statistical Test**: Binomial test for selection analysis
- **Scripts**: 
  - `scripts/gene_analysis/variation.py`
  - `scripts/general_gene_analysis/generate_gene_mapping_full.py`

### Functional Impact Analysis

- **Method**: Analyze functional consequences of variants
- **Impact Categories**: HIGH, MODERATE, LOW, MODIFIER (from snpEff)
- **Distance Analysis**: Variant impact distribution by distance from target genes
- **Annotation Integration**: Cross-reference with UniProt functional annotations
- **Scripts**:
  - `scripts/functional_impact/analyze_high_impact_variants.py`
  - `scripts/functional_impact/analyze_high_impact_variants_by_distance.py`
  - `scripts/variants/variant_proximity_impact_summary.py`

### Variant Proximity Impact Analysis

- **Method**: Analyze relationships between variant impact and distance to ERG genes
- **Distance Categories**: within gene, 0-1kb, 1-5kb, 5-10kb, 10-50kb, >50kb
- **Visualization**: Heatmaps, boxplots, and interactive HTML reports
- **Statistical Analysis**: Chi-square test for independence between impact and distance
- **Script**: `scripts/variants/variant_proximity_impact_summary.py`
- **Parameters**:
  - Maximum distance: 100kb from ergosterol genes
  - Impact categories: HIGH, MODERATE, LOW, MODIFIER
  - Treatment grouping: WT-37, WTA, CAS, STC
- **Outputs**: 
  - Variant count heatmap by ERG gene and treatment
  - Distance distribution heatmap
  - Interactive HTML report
  - Tabular data for further analysis

### Network Analysis

- **Method**: Construct and analyze gene interaction networks
- **Nodes**: Genes with variants or in ergosterol pathway
- **Edges**: Based on genetic distance and functional relationships
- **Metrics**: Degree centrality, betweenness centrality, clustering coefficient
- **Software**: NetworkX Python package
- **Script**: `scripts/functional_impact/build_extended_erg_network.py`

### OSH Gene Analysis

- **Method**: Map and analyze OSH gene family and relationship to ERG genes
- **Gene Identification**: Based on gene annotations and homology
- **Distance Calculation**: Compute genomic distances between OSH and ERG genes
- **Variant Analysis**: Examine variants near OSH genes across treatments
- **Software**: BioPython, Pandas, Matplotlib
- **Scripts**:
  - `scripts/osh_analysis/analyze_osh_genes.py`
  - `scripts/osh_analysis/osh_variants.py`
  - `scripts/osh_analysis/osh_erg_distance.py`
  - `scripts/osh_analysis/count_variant_locations.py`

### Enhanced Annotation Approach

- **Method**: Comprehensive annotation with genome-wide gene proximity
- **Proximity Categories**: Distance-based categorization to nearest genes
- **Gene Context**: Functional categories based on gene proximity
- **Integration**: Combines variant data with gene annotations, treatments, and impacts
- **Script**: `scripts/variants/extract_filtered_scaffold_variants.py`
- **Key Features**:
  - Maps variants to nearest genes and calculates distances
  - Enhances standard snpEff annotations with gene proximity data
  - Categorizes variants by distance from genes of interest
  - Creates expanded annotation for treatment-specific variants
  - Generates statistical summaries of proximity-impact relationships

## Statistical Methods

### Variant Statistics

- **Variant Calling**: Minimum QUAL ≥ 20, DP ≥ 10
- **Replication Consistency**: Variant presence in 1, 2, or 3 replicates
- **Treatment Comparison**: Chi-square test for independence
- **Multiple Testing Correction**: Benjamini-Hochberg FDR control
- **Significance Threshold**: Adjusted p < 0.05

### Enrichment Analysis

- **Method**: Binomial test comparing observed vs. expected counts
- **Expected Model**: 
  - Uniform distribution of variants across genome
  - Adjusted for mappable regions and GC content
- **Fold Enrichment**: Observed/expected ratio
- **Significance**: FDR-corrected p < 0.05
- **Visualization**: -log10(p-value) vs. log2(fold change)

### Correlation Analysis

- **Methods**:
  - Pearson correlation for continuous variables
  - Spearman rank correlation for ordinal data
  - Jaccard similarity for binary presence/absence data
- **Significance**: Permutation test with 10,000 iterations
- **Visualization**: Heatmaps with hierarchical clustering

### Dimensionality Reduction

- **PCA**:
  - Scaled and centered data
  - First 5 principal components analyzed
  - Variance explained reported for each component
- **MDS**:
  - Based on Jaccard distance matrix
  - Two dimensions for visualization
  - Stress value reported
- **Clustering**:
  - Hierarchical clustering (UPGMA method)
  - Dendrogram with bootstrap support values

## Visualization Methods

### Plot Generation

- **Software**: Python with matplotlib, seaborn, plotly
- **Color Schemes**: Consistent color mapping across analyses
- **Statistical Annotations**: p-values and significance stars (* p<0.05, ** p<0.01, *** p<0.001)
- **Output Formats**: PNG (publication quality), HTML (interactive)

### Interactive Dashboards

- **Framework**: Python with plotly and Bootstrap
- **Features**:
  - Interactive zooming and panning
  - Tooltips with detailed information
  - Filtering and selection controls
  - Tabbed organization of related visualizations
- **Generation Scripts**:
  - `scripts/utils/generate_ergosterol_variant_report.py`
  - `scripts/utils/generate_functional_impact_report.py`
  - `scripts/sterols/generate_html_report.py`
  - `scripts/variants/variant_proximity_impact_summary.py`
  - `scripts/variants/generate_filtered_variants_visualizations.py`

### Regulatory Analysis

- **Framework**: Custom Python scripts for regulatory region analysis
- **Features**:
  - Maps variants relative to transcription start sites (TSS)
  - Identifies variants in promoter regions
  - Analyzes potential transcription factor binding site impacts
  - Creates distribution visualizations for regulatory variants
- **Generation Scripts**:
  - `scripts/regulatory_analysis/regulatory_region_mapping.py`
  - `scripts/regulatory_analysis/regulatory_motif_analysis.py`
  - `scripts/regulatory_analysis/accurate_regulatory_analysis.py`

## File Formats

### Input Formats

- **Raw Sequencing**: FASTQ format
- **Reference Genome**: FASTA format
- **Genome Annotations**: GFF3 and GenBank formats
- **Sterol Data**: TSV with concentration measurements

### Intermediate Formats

- **Alignments**: BAM format with BAI indices
- **Variants**: VCF format with CSI indices
- **Mutation Data**: TSV with mutation categories and counts
- **Gene Mappings**: TSV with gene annotations and distances

### Output Formats

- **Reports**: Text, Markdown, and HTML formats
- **Visualizations**: PNG for static, HTML for interactive
- **Data Tables**: TSV and CSV formats
- **Network Data**: JSON format for node and edge data

## Hardware Requirements

The analysis pipeline was designed to run on a standard workstation with:

- **Processor**: Multi-core CPU (8+ cores recommended)
- **Memory**: 16GB RAM minimum, 32GB recommended
- **Storage**: 100GB free space for intermediate files
- **Runtime**: ~4-6 hours for complete pipeline execution

## Software Dependencies

### Python Packages

- **numpy**: Numerical computing
- **pandas**: Data manipulation
- **matplotlib**: Plotting and visualization
- **seaborn**: Statistical data visualization
- **scipy**: Scientific computing and statistics
- **biopython**: Biological sequence analysis
- **scikit-learn**: Machine learning for clustering and PCA
- **statsmodels**: Statistical models
- **networkx**: Network analysis
- **plotly**: Interactive visualizations
- **logomaker**: Sequence logo generation

### Command Line Tools

- **FastQC**: Quality control for sequencing data
- **BWA**: Burrows-Wheeler Aligner for read mapping
- **samtools**: SAM/BAM file manipulation
- **bcftools**: VCF file manipulation and variant calling
- **bedtools**: Genome arithmetic tools
- **snpEff**: Variant annotation and effect prediction

This technical reference provides a comprehensive overview of the methods, parameters, and tools used in the Yeast MSA project.