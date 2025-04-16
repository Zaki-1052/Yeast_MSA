#!/bin/bash

# File: scripts/annotation/03_run_snpeff_annotation.sh
# Purpose: Annotate VCF files using SnpEff with W303 genome database

echo "=== Running SnpEff Annotation ==="
echo "Date: $(date)"
echo ""

# Define directories
VCF_SOURCE="annotation/vcf_ready"
ANNO_DEST="annotation/results"
STATS_DIR="annotation/stats"

# Create output directories
mkdir -p "$ANNO_DEST"
mkdir -p "$STATS_DIR"

# Path to SnpEff - adjust if needed
SNPEFF_PATH="/Users/zakiralibhai/snpEff"
SNPEFF_JAR="$SNPEFF_PATH/snpEff.jar"

# Ensure SnpEff jar exists
if [ ! -f "$SNPEFF_JAR" ]; then
    echo "ERROR: SnpEff jar not found at $SNPEFF_JAR"
    echo "Please provide the correct path to your SnpEff installation"
    exit 1
fi

# Ensure w303 database is available
if ! java -jar "$SNPEFF_JAR" databases | grep -q "w303"; then
    echo "ERROR: w303 database not found in SnpEff"
    echo "Please ensure the w303 database has been properly built"
    exit 1
fi

# Find all sorted VCF files
echo "Scanning for sorted VCF files..."
VCF_FILES=$(find "$VCF_SOURCE" -name "*.sorted.vcf.gz")
VCF_COUNT=$(echo "$VCF_FILES" | wc -l)

if [ "$VCF_COUNT" -eq 0 ]; then
    echo "ERROR: No sorted VCF files found in $VCF_SOURCE"
    exit 1
fi

echo "Found $VCF_COUNT VCF files to annotate"
echo ""

# Create a log file
LOG_FILE="annotation/snpeff_annotation.log"
echo "SnpEff Annotation Log" > "$LOG_FILE"
echo "Date: $(date)" >> "$LOG_FILE"
echo "==============================================" >> "$LOG_FILE"
echo "SnpEff version: $(java -jar $SNPEFF_JAR -version 2>&1 | head -n 1)" >> "$LOG_FILE"
echo "Genome database: w303" >> "$LOG_FILE"
echo "" >> "$LOG_FILE"

# Process each VCF file
for VCF_FILE in $VCF_FILES; do
    FILENAME=$(basename "$VCF_FILE")
    SAMPLE=${FILENAME%.sorted.vcf.gz}
    OUTPUT_FILE="$ANNO_DEST/${SAMPLE}.snpeff.vcf"
    STATS_FILE="$STATS_DIR/${SAMPLE}.snpeff.stats.html"
    
    echo "Processing $SAMPLE..."
    echo "[$SAMPLE] Annotation started at $(date)" >> "$LOG_FILE"
    
    echo "- Running SnpEff annotation..."
    echo "  Input: $VCF_FILE"
    echo "  Output: $OUTPUT_FILE"
    echo "  Stats: $STATS_FILE"
    
    # Run SnpEff with detailed statistics
    java -Xmx4g -jar "$SNPEFF_JAR" \
        -v \
        -stats "$STATS_FILE" \
        w303 \
        "$VCF_FILE" \
        > "$OUTPUT_FILE" \
        2>> "$LOG_FILE"
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Annotation successful"
        echo "[$SAMPLE] Annotation successful" >> "$LOG_FILE"
        
        # Compress and index the annotated VCF
        echo "- Compressing annotated VCF file..."
        bgzip -f "$OUTPUT_FILE"
        
        if [ $? -eq 0 ]; then
            echo "  ✓ Compression successful"
            echo "[$SAMPLE] Compression successful" >> "$LOG_FILE"
            
            echo "- Indexing annotated VCF file..."
            tabix -p vcf "$OUTPUT_FILE.gz"
            
            if [ $? -eq 0 ]; then
                echo "  ✓ Indexing successful"
                echo "[$SAMPLE] Indexing successful" >> "$LOG_FILE"
            else
                echo "  ✗ Indexing failed!"
                echo "[$SAMPLE] ERROR: Indexing failed!" >> "$LOG_FILE"
            fi
        else
            echo "  ✗ Compression failed!"
            echo "[$SAMPLE] ERROR: Compression failed!" >> "$LOG_FILE"
        fi
    else
        echo "  ✗ Annotation failed!"
        echo "[$SAMPLE] ERROR: Annotation failed!" >> "$LOG_FILE"
    fi
    
    echo "Done processing $SAMPLE"
    echo "[$SAMPLE] Processing completed at $(date)" >> "$LOG_FILE"
    echo "" >> "$LOG_FILE"
    echo ""
done

# Count successfully annotated files
ANNO_COUNT=$(ls "$ANNO_DEST"/*.snpeff.vcf.gz 2>/dev/null | wc -l)
echo "Summary: $ANNO_COUNT of $VCF_COUNT files successfully annotated"
echo "Log file saved to: $LOG_FILE"
echo ""

# Create a consolidated summary report
SUMMARY_FILE="annotation/annotation_summary.txt"
echo "SnpEff Annotation Summary" > "$SUMMARY_FILE"
echo "Date: $(date)" >> "$SUMMARY_FILE"
echo "==============================================" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"
echo "Total files processed: $VCF_COUNT" >> "$SUMMARY_FILE"
echo "Successfully annotated: $ANNO_COUNT" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"
echo "Files with annotation:" >> "$SUMMARY_FILE"

for ANNO_FILE in "$ANNO_DEST"/*.snpeff.vcf.gz; do
    if [ -f "$ANNO_FILE" ]; then
        FILENAME=$(basename "$ANNO_FILE")
        SAMPLE=${FILENAME%.snpeff.vcf.gz}
        
        # Count variants in file
        VARIANT_COUNT=$(bcftools view -H "$ANNO_FILE" | wc -l)
        
        # Count annotations
        ANNO_COUNT=$(bcftools view -H "$ANNO_FILE" | grep -c "ANN=")
        
        # Add to summary
        echo "- $SAMPLE: $VARIANT_COUNT variants, $ANNO_COUNT with annotations" >> "$SUMMARY_FILE"
    fi
done

echo "" >> "$SUMMARY_FILE"
echo "Statistics reports are available in: $STATS_DIR" >> "$SUMMARY_FILE"

echo "Annotation summary saved to: $SUMMARY_FILE"
echo ""
echo "=== SnpEff Annotation complete ==="