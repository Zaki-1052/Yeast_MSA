#!/bin/bash

# Script to reorganize files into the new directory structure
# This will move existing analysis scripts into their appropriate subdirectories

# Mutation spectrum analysis
mv scripts/mutation_spectrum_analysis.py scripts/analysis/mutation_spectrum/
cp scripts/analysis/__init__.py scripts/analysis/mutation_spectrum/

# Genomic context analysis
mv scripts/genomic_context_analysis.py scripts/analysis/genomic_context/
cp scripts/analysis/__init__.py scripts/analysis/genomic_context/

# Mutational signature analysis
mv scripts/mutational_signature_analysis.py scripts/analysis/mutational_signatures/
cp scripts/analysis/__init__.py scripts/analysis/mutational_signatures/

# Population structure analysis
mv scripts/population_spectrum_analysis.py scripts/analysis/population_structure/
cp scripts/analysis/__init__.py scripts/analysis/population_structure/

# Regional enrichment analysis
mv scripts/regional_enrichment_analysis.py scripts/analysis/regional_enrichment/
cp scripts/analysis/__init__.py scripts/analysis/regional_enrichment/

# Scaffold distribution analysis
mv scripts/scaffold_distribution_analysis.py scripts/analysis/scaffold_distribution/
cp scripts/analysis/__init__.py scripts/analysis/scaffold_distribution/

# Statistical pattern analysis
mv scripts/statistical_pattern_analysis.py scripts/analysis/statistical_patterns/
cp scripts/analysis/__init__.py scripts/analysis/statistical_patterns/

# Move pipeline scripts
cp scripts/pipeline/*.sh scripts/pipeline/
cp scripts/phase1/preprocessing.sh scripts/pipeline/01_preprocessing.sh
cp scripts/phase1/alignment.sh scripts/pipeline/02_alignment.sh
cp scripts/vcf/variants.sh scripts/pipeline/03_variant_calling.sh
cp scripts/pipeline/02_filter_normalize.sh scripts/pipeline/04_filtering.sh
cp scripts/annotation/jriu_based_annotator.py scripts/pipeline/05_annotation.sh
cp scripts/vcf/comparisons.sh scripts/pipeline/06_comparisons.sh

# Create master pipeline script
cat > scripts/pipeline/07_run_full_pipeline.sh << 'EOF'
#!/bin/bash

# Master script to run the entire analysis pipeline
echo "Running Yeast MSA analysis pipeline..."

# Run preprocessing
echo "Step 1: Preprocessing raw sequencing data"
bash scripts/pipeline/01_preprocessing.sh

# Run alignment
echo "Step 2: Aligning reads to reference genome"
bash scripts/pipeline/02_alignment.sh

# Run variant calling
echo "Step 3: Calling genetic variants"
bash scripts/pipeline/03_variant_calling.sh

# Run filtering
echo "Step 4: Filtering and normalizing variants"
bash scripts/pipeline/04_filtering.sh

# Run annotation
echo "Step 5: Annotating variants"
bash scripts/pipeline/05_annotation.sh

# Run comparisons
echo "Step 6: Comparing variants across treatments"
bash scripts/pipeline/06_comparisons.sh

echo "Pipeline complete! Results are in the 'results' directory."
EOF

chmod +x scripts/pipeline/07_run_full_pipeline.sh

# Create utility scripts directory with common functions
mkdir -p scripts/utils
cat > scripts/utils/config.py << 'EOF'
#!/usr/bin/env python3

"""
Global configuration parameters for the Yeast MSA analysis.
"""

# Treatment information
TREATMENTS = ['WT-37', 'WTA', 'STC', 'CAS']

# Define treatment information for better biological context
TREATMENT_INFO = {
    'WT-37': {'description': 'Temperature-adapted wild type', 'adaptation': 'Temperature'},
    'WTA': {'description': 'Low oxygen-adapted wild type', 'adaptation': 'Low Oxygen'},
    'STC': {'description': 'STC gene with low oxygen adaptation', 'adaptation': 'Low Oxygen', 'gene': 'STC'},
    'CAS': {'description': 'CAS gene with temperature adaptation', 'adaptation': 'Temperature', 'gene': 'CAS'}
}

# Treatment colors for consistent visualization
TREATMENT_COLORS = {
    'WT-37': '#1b9e77',  # Temperature-adapted
    'WTA': '#d95f02',    # Low oxygen-adapted
    'STC': '#7570b3',    # STC gene + low oxygen
    'CAS': '#e7298a'     # CAS gene + temperature
}

# Reference genome path
REFERENCE_GENOME = "reference/yeast_w303.fasta"

# Output directories
OUTPUT_DIR = "results"
EOF

# Create README for scripts directory
cat > scripts/README.md << 'EOF'
# Yeast MSA Analysis Scripts

This directory contains all scripts used in the Yeast MSA analysis pipeline.

## Directory Structure

- `pipeline/`: Sequential scripts for running the complete analysis pipeline
- `analysis/`: Modules for different types of analyses
- `visualization/`: Scripts for generating visualizations
- `utils/`: Utility functions used across multiple scripts
- `reporting/`: Scripts for generating analysis reports

## Usage

For running the complete pipeline:

```bash
bash scripts/pipeline/07_run_full_pipeline.sh
```

For running specific analysis modules:

```bash
python scripts/analysis/mutation_spectrum/spectrum_analysis.py
```

See the project_organization.md file for more details on how the scripts and results are organized.
EOF

# Create visualization directory
mkdir -p scripts/visualization
cat > scripts/visualization/core_plots.R << 'EOF'
#!/usr/bin/env Rscript

# Core plotting functions for Yeast MSA analysis
# This script provides reusable plotting functions used across different analyses

library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Function to create a standardized mutation spectrum plot
create_mutation_spectrum <- function(data, output_file, title = "Mutation Spectrum") {
  p <- ggplot(data, aes(x = Mutation, y = Count, fill = Treatment)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(
      title = title,
      x = "Mutation Type",
      y = "Count"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
  
  ggsave(output_file, p, width = 10, height = 6, dpi = 300)
  return(p)
}

# Function to create a standardized heatmap
create_heatmap <- function(data, output_file, title = "Heatmap") {
  # Implementation depends on specific data structure
  # This is a placeholder for the actual implementation
  message("Heatmap function called with data for: ", title)
  message("Output will be saved to: ", output_file)
}

# Function to create a standardized PCA plot
create_pca_plot <- function(data, output_file, title = "PCA Analysis") {
  # Implementation depends on specific data structure
  # This is a placeholder for the actual implementation
  message("PCA plot function called with data for: ", title)
  message("Output will be saved to: ", output_file)
}
EOF

echo "Reorganization script created. Run with: bash scripts/reorganize.sh"
echo "This script will move existing files into the new directory structure."