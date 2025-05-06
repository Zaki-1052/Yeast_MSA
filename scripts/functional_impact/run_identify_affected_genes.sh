#!/bin/bash

# run_identify_affected_genes.sh - FIXED VERSION 2
# Script to identify and characterize genes affected by HIGH/MODERATE impact variants

# Set default values
VARIANTS_FILE="results/functional_impact/variants_by_distance/high_impact_variants_by_distance.tsv"
GENBANK_DIR="reference/w303_annotations"
MAPPING_FILE="reference/chromosome_mapping.tsv"
OUTPUT_DIR="results/functional_impact/affected_genes_fixed"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run the analysis with debugging enabled
echo "Identifying genes affected by HIGH/MODERATE impact variants..."
python scripts/functional_impact/identify_affected_genes.py \
  --variants_file "$VARIANTS_FILE" \
  --genbank_dir "$GENBANK_DIR" \
  --mapping_file "$MAPPING_FILE" \
  --output_dir "$OUTPUT_DIR" \
  --debug

echo "Analysis complete! Results are in $OUTPUT_DIR"