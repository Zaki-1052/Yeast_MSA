#!/bin/bash

# run_promoter_analysis.sh
# Script to run promoter element analysis on ergosterol pathway gene variants

# Set default values
VARIANT_FILE="results/gene_variants_expanded/all_gene_variants.tsv"
GENE_MAPPING="reference/genes_of_interest_mapping.tsv"
OUTPUT_DIR="results/regulatory_analysis/promoters"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run the promoter element analysis
echo "Running promoter element analysis..."
python scripts/functional_impact/analyze_promoter_elements.py \
  --variant_file "$VARIANT_FILE" \
  --gene_mapping "$GENE_MAPPING" \
  --output_dir "$OUTPUT_DIR"

echo "Analysis complete! Results are in $OUTPUT_DIR"