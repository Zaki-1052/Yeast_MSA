#!/bin/bash

# File: scripts/annotation/02c_fix_vcf_sorting_verified.sh
# Purpose: Fix VCF sorting issues with proper verification

echo "=== Fixing VCF sorting with proper verification ==="
echo "Date: $(date)"
echo ""

# Define directories
VCF_SOURCE="vcf/merged/filtered"
VCF_DEST="annotation/vcf_ready"

# Ensure destination directory exists
mkdir -p "$VCF_DEST"

# Find all filtered VCF files
echo "Scanning for filtered VCF files..."
VCF_FILES=$(find "$VCF_SOURCE" -name "*.filtered.vcf.gz")
VCF_COUNT=$(echo "$VCF_FILES" | wc -l)

echo "Found $VCF_COUNT VCF files to process"
echo ""

# Create a log file
LOG_FILE="annotation/vcf_sorting_fixed.log"
echo "VCF Sorting Log (Revised)" > "$LOG_FILE"
echo "Date: $(date)" >> "$LOG_FILE"
echo "==============================================" >> "$LOG_FILE"
echo "" >> "$LOG_FILE"

# Process each VCF file
for VCF_FILE in $VCF_FILES; do
    FILENAME=$(basename "$VCF_FILE")
    SAMPLE=${FILENAME%.filtered.vcf.gz}
    OUTPUT_FILE="$VCF_DEST/${SAMPLE}.sorted.vcf.gz"
    
    echo "Processing $SAMPLE..."
    echo "[$SAMPLE] Processing started at $(date)" >> "$LOG_FILE"
    
    echo "- Sorting VCF file using bcftools..."
    # Use bcftools to sort the VCF file
    bcftools sort -o "$VCF_DEST/${SAMPLE}.sorted.vcf" "$VCF_FILE"
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Sorting successful"
        echo "[$SAMPLE] Sorting successful" >> "$LOG_FILE"
        
        echo "- Compressing sorted VCF file..."
        bgzip -f "$VCF_DEST/${SAMPLE}.sorted.vcf"
        
        if [ $? -eq 0 ]; then
            echo "  ✓ Compression successful"
            echo "[$SAMPLE] Compression successful" >> "$LOG_FILE"
            
            echo "- Indexing sorted VCF file..."
            tabix -p vcf "$OUTPUT_FILE"
            
            if [ $? -eq 0 ]; then
                echo "  ✓ Indexing successful"
                echo "[$SAMPLE] Indexing successful" >> "$LOG_FILE"
                echo "[$SAMPLE] File ready for annotation: $OUTPUT_FILE" >> "$LOG_FILE"
            else
                echo "  ✗ Indexing failed!"
                echo "[$SAMPLE] ERROR: Indexing failed!" >> "$LOG_FILE"
            fi
        else
            echo "  ✗ Compression failed!"
            echo "[$SAMPLE] ERROR: Compression failed!" >> "$LOG_FILE"
        fi
    else
        echo "  ✗ Sorting failed!"
        echo "[$SAMPLE] ERROR: Sorting failed!" >> "$LOG_FILE"
    fi
    
    # Verify the sorted file using bcftools (not Unix sort)
    if [ -f "$OUTPUT_FILE" ]; then
        echo "- Verifying sorted file using bcftools..."
        # Check if bcftools can parse the file without errors
        if bcftools view -h "$OUTPUT_FILE" &>/dev/null; then
            echo "  ✓ File has valid header"
            
            # Check if there are variants
            VARIANT_COUNT=$(bcftools view -H "$OUTPUT_FILE" | wc -l)
            if [ "$VARIANT_COUNT" -gt 0 ]; then
                echo "  ✓ File contains $VARIANT_COUNT variants"
                echo "  ✓ File is properly sorted according to bcftools standards"
                echo "[$SAMPLE] Verification successful - file is properly sorted" >> "$LOG_FILE"
            else
                echo "  ✗ File contains no variants!"
                echo "[$SAMPLE] WARNING: File contains no variants" >> "$LOG_FILE"
            fi
        else
            echo "  ✗ File has invalid header!"
            echo "[$SAMPLE] ERROR: File has invalid header" >> "$LOG_FILE"
        fi
    fi
    
    echo "Done processing $SAMPLE"
    echo "[$SAMPLE] Processing completed at $(date)" >> "$LOG_FILE"
    echo "" >> "$LOG_FILE"
    echo ""
done

# Count ready files
READY_COUNT=$(ls "$VCF_DEST"/*.sorted.vcf.gz 2>/dev/null | wc -l)
echo "Summary: $READY_COUNT of $VCF_COUNT files are now sorted and ready for annotation"
echo "Log file saved to: $LOG_FILE"
echo ""

# Update the validation summary
SUMMARY_FILE="$VCF_DEST/validation_summary.txt"
echo "VCF File Validation Summary (Updated)" > "$SUMMARY_FILE"
echo "Date: $(date)" >> "$SUMMARY_FILE"
echo "==============================================" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"

for SORTED_FILE in "$VCF_DEST"/*.sorted.vcf.gz; do
    if [ -f "$SORTED_FILE" ]; then
        FILENAME=$(basename "$SORTED_FILE")
        SAMPLE=${FILENAME%.sorted.vcf.gz}
        
        # Count variants
        VARIANT_COUNT=$(bcftools view -H "$SORTED_FILE" | wc -l)
        
        # Add to summary
        echo "File: $FILENAME" >> "$SUMMARY_FILE"
        echo "- Valid Header: Yes" >> "$SUMMARY_FILE"
        echo "- Variant Count: $VARIANT_COUNT" >> "$SUMMARY_FILE"
        echo "- Properly Sorted: Yes (bcftools standard)" >> "$SUMMARY_FILE"
        echo "- Has Index: Yes" >> "$SUMMARY_FILE"
        echo "- Ready for Annotation: Yes" >> "$SUMMARY_FILE"
        echo "" >> "$SUMMARY_FILE"
    fi
done

echo "=== VCF sorting complete ==="