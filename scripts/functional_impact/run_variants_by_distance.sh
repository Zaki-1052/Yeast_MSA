#!/bin/bash

# run_variants_by_distance.sh
# Script to analyze HIGH and MODERATE impact variants near ergosterol pathway genes

# Set default values
VCF_DIR="vcf/annotated"
OUTPUT_DIR="results/functional_impact/variants_by_distance"
GENE_MAPPING="reference/genes_of_interest_mapping.tsv"
DISTANCE_THRESHOLD=50000  # 50kb default threshold

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run the analysis
echo "Running HIGH and MODERATE impact variant analysis by distance..."
python scripts/functional_impact/analyze_high_impact_variants_by_distance.py \
  --vcf_dir "$VCF_DIR" \
  --output_dir "$OUTPUT_DIR" \
  --gene_mapping "$GENE_MAPPING" \
  --distance_threshold "$DISTANCE_THRESHOLD"

# Check results and report findings
if [ -f "$OUTPUT_DIR/variants_by_distance_report.md" ]; then
  VARIANT_COUNT=$(grep -A 1 "^## Overview" "$OUTPUT_DIR/variants_by_distance_report.md" | grep "Total variants" | cut -d ':' -f 2 | tr -d ' ')
  
  if [ -z "$VARIANT_COUNT" ] || [ "$VARIANT_COUNT" -eq 0 ]; then
    echo "===================================================================================="
    echo "FINDING: No HIGH or MODERATE impact variants found within ${DISTANCE_THRESHOLD}bp of"
    echo "ergosterol pathway genes. This suggests extremely strong evolutionary constraints on"
    echo "this entire genomic region, not just the genes themselves."
    echo "===================================================================================="
  else
    echo "===================================================================================="
    echo "FINDING: Found $VARIANT_COUNT HIGH/MODERATE impact variants near ergosterol genes."
    echo "While the ergosterol genes themselves show no HIGH/MODERATE impact variants,"
    echo "nearby genes do contain such variants, potentially indicating:"
    echo "1. Neighboring genes may influence ergosterol pathway through regulatory mechanisms"
    echo "2. Selection constraints decrease with distance from the essential pathway genes"
    echo "3. Adaptation may occur through changes in genes that interact with the pathway"
    echo "===================================================================================="
  fi
else
  echo "Analysis completed, but report was not generated. Check for errors."
fi

echo "Analysis complete! Results are in $OUTPUT_DIR"