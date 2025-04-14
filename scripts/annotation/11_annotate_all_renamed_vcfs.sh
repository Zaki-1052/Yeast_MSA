#!/bin/bash

# File: scripts/annotation/11_annotate_all_renamed_vcfs.sh
# Purpose: Annotate all renamed VCF files

echo "=== Annotating All Renamed VCF Files ==="
echo "Date: $(date)"
echo ""

# Define directories
RENAMED_VCF_DIR="annotation/vcf_renamed"
RESULTS_DIR="annotation/results_renamed"
STATS_DIR="annotation/stats_renamed"
SNPEFF_JAR="/Users/zakiralibhai/snpEff/snpEff.jar"

# Create output directories
mkdir -p "$RESULTS_DIR"
mkdir -p "$STATS_DIR"

# Find all renamed VCF files
VCF_FILES=$(find "$RENAMED_VCF_DIR" -name "*.renamed.vcf.gz")

# Process each file
for VCF_FILE in $VCF_FILES; do
    SAMPLE=$(basename "$VCF_FILE" .renamed.vcf.gz)
    echo "Processing $SAMPLE..."
    
    # Run SnpEff
    java -Xmx4g -jar "$SNPEFF_JAR" \
        -v \
        -stats "$STATS_DIR/${SAMPLE}.stats.html" \
        w303 \
        "$VCF_FILE" \
        > "$RESULTS_DIR/${SAMPLE}.snpeff.vcf"
    
    # Compress and index
    bgzip -f "$RESULTS_DIR/${SAMPLE}.snpeff.vcf"
    tabix -p vcf "$RESULTS_DIR/${SAMPLE}.snpeff.vcf.gz"
    
    echo "Done processing $SAMPLE"
    echo ""
done

echo "All VCF files have been annotated"
echo ""
echo "=== Annotation Complete ==="
