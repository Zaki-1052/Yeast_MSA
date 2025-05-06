#!/bin/bash

# run_extended_network_analysis.sh
# Script to build and analyze the extended ergosterol pathway network

# Set default values
ERG_GENES="reference/genes_of_interest_mapping.tsv"
AFFECTED_GENES="results/functional_impact/affected_genes_fixed/affected_genes.tsv"
VARIANTS="results/functional_impact/affected_genes_fixed/variants_with_genes.tsv" 
OUTPUT_DIR="results/network_analysis"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run the network analysis
echo "Building extended ergosterol pathway network..."
python scripts/functional_impact/build_extended_erg_network.py \
  --erg_genes "$ERG_GENES" \
  --affected_genes "$AFFECTED_GENES" \
  --variants "$VARIANTS" \
  --output_dir "$OUTPUT_DIR"

echo "Analysis complete! Results are in $OUTPUT_DIR"