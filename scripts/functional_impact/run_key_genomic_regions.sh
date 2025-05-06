#!/bin/bash

# run_key_genomic_regions.sh
# Script to analyze key genomic regions showing interesting variant patterns

# Set default values
VARIANTS_FILE="results/functional_impact/variants_by_distance/high_impact_variants_by_distance.tsv"
GENOME_FILE="reference/w303_chromosomal.fasta"
GENBANK_DIR="reference/w303_annotations"
GENE_MAPPING="reference/genes_of_interest_mapping.tsv"
MAPPING_FILE="reference/chromosome_mapping.tsv"
OUTPUT_DIR="results/functional_impact/key_genomic_regions"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run the analysis
echo "Running key genomic regions analysis..."
python scripts/functional_impact/analyze_key_genomic_regions.py \
  --variants_file "$VARIANTS_FILE" \
  --genome_file "$GENOME_FILE" \
  --genbank_dir "$GENBANK_DIR" \
  --gene_mapping "$GENE_MAPPING" \
  --mapping_file "$MAPPING_FILE" \
  --output_dir "$OUTPUT_DIR"

echo "Analysis complete! Results are in $OUTPUT_DIR"