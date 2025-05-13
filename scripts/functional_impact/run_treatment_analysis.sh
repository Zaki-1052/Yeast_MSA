#!/bin/bash

# run_treatment_analysis.sh
# Script to analyze treatment-specific variant patterns

# Set default values
VARIANTS_FILE="results/functional_impact/key_genomic_regions/region_variants.tsv"
OUTPUT_DIR="results/treatment_analysis"
CONTROL_GROUP="WT"
INCLUDE_CONTROLS="--include_controls"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run the analysis
echo "Running treatment-specific variant analysis..."
python scripts/functional_impact/analyze_treatment_specific_patterns.py \
  --variants_file "$VARIANTS_FILE" \
  --output_dir "$OUTPUT_DIR" \
  --control_group "$CONTROL_GROUP"

echo "Analysis complete! Results are in $OUTPUT_DIR"