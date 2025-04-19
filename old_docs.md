# Yeast Mutation Accumulation and Adaptation Analysis Pipeline

## Project Overview
This repository implements a comprehensive end‑to‑end pipeline for analyzing mutation accumulation and adaptation experiments in yeast under various treatments. The goals of the project are:
- To preprocess and map sequencing reads to a high‑quality yeast reference genome
- To call, filter, merge, and compare genetic variants (VCFs) across treatment and control samples
- To annotate variants with functional and genomic context (e.g., gene, scaffold)
- To perform comparative analyses and generate summary statistics and visualizations:
  - Mutation spectrum analysis
  - Genomic context characterization
  - Mutational signature extraction
  - Population structure and variant sharing
  - Regional enrichment and hotspot identification
  - Scaffold distribution and density analysis
  - Statistical pattern exploration (correlation, PCA, regression)
  - Treatment vs. control comparisons and fold‑change analysis

## Experimental Design
- Biological samples: replicates (e.g., CAS‑55‑1, STC‑55‑2, WT‑37‑55‑3, WTA‑55‑1) and corresponding controls (e.g., CAS‑CTRL). 55 generations under different treatments (CAS, STC, WT‑37, WTA).
- Sequencing data: paired‑end reads (FASTQ) cleaned and processed.
- Treatments represent environmental or chemical conditions applied to yeast populations.

## Directory Structure
```
├── data/                   # Raw and cleaned sequencing reads
│   └── processed_data/     # Fastq.gz files after quality trimming and filtering
├── reference/              # Yeast reference genome in FASTA and GenBank formats
├── scripts/                # Analysis, annotation, and pipeline scripts
│   ├── pipeline/           # Main bash pipeline stages (01_*.sh … 12_run_pipeline.sh)
│   ├── phase1/             # Legacy or exploratory preprocessing scripts
│   ├── vcf/                # Variant file (VCF) manipulation and comparison utilities
│   ├── annotation/         # Variant annotation and gene/scaffold analysis scripts
│   ├── analysis/           # Python scripts for downstream statistical and graphical analyses
│   └── old/                # Archived previous versions of scripts
├── annotation/             # Annotated VCFs, reports, and enhanced annotation outputs
├── analysis/               # Intermediate analysis outputs (plots, summaries by category)
├── mutation_spectrum_analysis/  # Raw mutation count inputs for spectrum analysis
├── results/                # Final pipeline outputs: BAMs, merged VCFs, summaries, and plots
│   ├── bam/
│   ├── merged/
│   ├── gene_scaffold_analysis/
│   ├── analysis/           # Organized results matching the analysis/ directory
│   └── preprocessing/
└── old_version/            # Archived legacy pipeline, analyses, and conclusions
```

## Pipeline Workflow
The `scripts/pipeline/12_run_pipeline.sh` script orchestrates the following steps:
1. **Setup** (`01_setup.sh`): prepare directory structure, download or link reference files, build indices (BWA, samtools, snpEff).
2. **Filter & Normalize** (`02_filter_normalize.sh`): apply quality filters to VCFs, normalize variant representation.
3. **Fix Contigs** (`03_fix_contigs.sh`): correct contig naming and ordering for consistency.
4. **Merge VCFs** (`04_merge_vcfs.sh`): combine replicates into group and high‑confidence VCFs per treatment.
5. **Test Samples** (`05_test_samples.sh`): validate sample inclusion and check for outliers.
6. **Group Comparisons** (`06_compare_groups.sh`): perform treatment vs. control comparisons.
7. **Cross‑Treatment** (`07_cross_treatment.sh`): identify shared and unique variants across treatments.
8. **Identify Unique** (`08_identify_unique.sh`): extract variants specific to each treatment.
9. **Direct VCF Comparison** (`09_direct_comparison.sh`): pairwise VCF comparisons with controls.
10. **Consistency Analysis** (`10_consistency_analysis.sh`): measure within‑group replicate consistency (e.g., Jaccard index).
11. **Summary Report** (`11_summary_report.sh`): compile counts, rates, and summary tables into `results/merged/summary`.

To execute the full pipeline:
```bash
bash scripts/pipeline/12_run_pipeline.sh
```

## Variant Annotation
The `scripts/annotation/` directory contains tools to annotate VCFs with gene and scaffold context:
- **snpEff‑based annotation**: prepare and run snpEff databases, generate annotated VCFs and impact reports.
- **Custom annotators** (`variant_gene_annotator.py`, `jriu_based_annotator.py`, `enhanced_jriu_annotator.py`): map variants to genes of interest, validate with GenBank annotations.
- **Gene & scaffold analysis** (`gene_scaffold_analysis.py`): compute variant densities and enrichment across genes or scaffolds.
- **Extraction utilities**: scripts to isolate variants in user‑defined genes of interest (`extract_genes_of_interest.py`).

Annotation outputs are stored under `annotation/` and integrated into final reports (HTML, TSV).

## Downstream Analysis
Python scripts in `scripts/analysis/` process merged and annotated data to generate figures and statistical summaries. Key analyses include:
- **Genomic Context** (`genomic_context_analysis.py`): flanking base composition, homopolymer and GC content distributions.
- **Mutation Spectrum** (`mutation_spectrum_analysis.py`): 6‑type mutation frequency plots and Ti/Tv ratios.
- **Mutational Signatures** (`mutational_signature_analysis.py`): de novo signature extraction and context logos.
- **Population Structure** (`population_spectrum_analysis.py`): PCA, MDS, dendrograms, variant sharing heatmaps.
- **Regional Enrichment** (`regional_enrichment_analysis.py`): identify genomic regions enriched for variants per treatment.
- **Scaffold Distribution** (`scaffold_distribution_analysis.py`): density heatmaps and hotspot identification across scaffolds.
- **Statistical Patterns** (`statistical_pattern_analysis.py`): correlation matrices, regression analyses, PCA biplots.
- **Treatment vs Control** (`TC_Visualization.py`): fold‑change barplots and comparative statistics.
- **Utility** (`variation.py`): core functions for reading variant tables and summarizing counts.

Generated figures and tables are found in the top‑level `analysis/` directory, organized by result category.

## Results
The `results/` directory mirrors the analysis structure with final BAM files (`results/bam/`), merged VCFs (`results/merged/`), gene‑scaffold analysis outputs, and organized subfolders under `results/analysis/` for each analysis type. Key deliverables include:
- High‑confidence variant sets per treatment
- Treatment‑unique variant lists
- Mutation spectrum and signature plots
- PCA/MDS plots for population structure
- Regional and scaffold enrichment tables and heatmaps
- Statistical summary reports (text and HTML)

## Dependencies
- Unix‐like environment with Bash, GNU coreutils
- BWA, samtools, bcftools, bedtools
- snpEff for variant annotation
- Python 3.x with packages: pandas, numpy, matplotlib, seaborn, scipy
- R (optional) for RStudio project (`yeast_analysis.Rproj`)

## Reproducing the Analysis
1. Install software dependencies and clone this repository.
2. Place raw FASTQ files in `data/processed_data/` following naming conventions.
3. Update paths to reference genome and snpEff database in `scripts/pipeline/01_setup.sh`.
4. Run the pipeline: `bash scripts/pipeline/12_run_pipeline.sh`.
5. View final results in `results/` and visualization outputs in `analysis/`.
