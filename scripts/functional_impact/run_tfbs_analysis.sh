#!/bin/bash

# run_tfbs_analysis.sh
# Script to analyze how variants might affect transcription factor binding sites

# Set default values
VARIANT_FILE="results/regulatory_analysis/promoters/variants_with_tss_distance.tsv"
GENOME_FILE="reference/w303_chromosomal.fasta"
MAPPING_FILE="reference/chromosome_mapping.tsv"
OUTPUT_DIR="results/regulatory_analysis/tfbs"
CONTEXT_SIZE=30
PREDICTION_METHOD="motif_match"  # Options: jaspar, motif_match, both

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run the TFBS analysis
echo "Running transcription factor binding site analysis..."
python scripts/functional_impact/analyze_tfbs.py \
  --variant_file "$VARIANT_FILE" \
  --genome_file "$GENOME_FILE" \
  --mapping_file "$MAPPING_FILE" \
  --output_dir "$OUTPUT_DIR" \
  --context_size "$CONTEXT_SIZE" \
  --prediction_method "$PREDICTION_METHOD"

echo "Analysis complete! Results are in $OUTPUT_DIR"