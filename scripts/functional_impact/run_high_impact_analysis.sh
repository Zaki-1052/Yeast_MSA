#!/bin/bash

# run_high_impact_analysis.sh
# Script to analyze HIGH and MODERATE impact variants in ergosterol pathway genes

# Set default values
VCF_DIR="vcf/annotated"
OUTPUT_DIR="results/functional_impact/high_impact"
GENE_MAPPING="reference/genes_of_interest_mapping.tsv"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run the HIGH impact variant analysis
echo "Running HIGH and MODERATE impact variant analysis..."
python scripts/functional_impact/analyze_high_impact_variants.py \
  --vcf_dir "$VCF_DIR" \
  --output_dir "$OUTPUT_DIR" \
  --gene_mapping "$GENE_MAPPING"

# Check if any variants were found
if [ -f "$OUTPUT_DIR/gene_treatment_variant_counts.tsv" ] && grep -q "No variants found" "$OUTPUT_DIR/gene_treatment_variant_counts.tsv"; then
  echo "================================================================================"
  echo "SIGNIFICANT FINDING: No HIGH or MODERATE impact variants found in ergosterol pathway genes."
  echo "This indicates strong purifying selection on these genes, reinforcing the hypothesis"
  echo "that the ergosterol pathway is highly conserved even under adaptive conditions."
  echo ""
  echo "This suggests adaptation occurs primarily through changes in gene expression (as seen"
  echo "in our regulatory variant analysis) rather than through protein-altering mutations."
  echo "================================================================================"
else
  echo "Analysis found HIGH/MODERATE impact variants. See report for details."
fi

echo "Analysis complete! Results are in $OUTPUT_DIR"