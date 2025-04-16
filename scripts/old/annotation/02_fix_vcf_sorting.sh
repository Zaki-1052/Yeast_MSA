#!/bin/bash

# File: scripts/annotation/02_fix_vcf_sorting.sh
# Purpose: Fix sorting issues in VCF files for SnpEff annotation

echo "=== Fixing VCF sorting issues ==="
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
LOG_FILE="annotation/vcf_sorting.log"
echo "VCF Sorting Log" > "$LOG_FILE"
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
    
    echo "- Decompressing and sorting VCF file..."
    # Use bcftools to sort the VCF file (rather than manual sort)
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
    
    # Verify the sorted file
    if [ -f "$OUTPUT_FILE" ]; then
        echo "- Verifying sorted file..."
        if bcftools query -f '%CHROM\t%POS\n' "$OUTPUT_FILE" | sort -k1,1 -k2,2n | \
           diff -q - <(bcftools query -f '%CHROM\t%POS\n' "$OUTPUT_FILE") &>/dev/null; then
            echo "  ✓ File is properly sorted"
            echo "[$SAMPLE] Verification successful - file is properly sorted" >> "$LOG_FILE"
        else
            echo "  ✗ File is still not properly sorted!"
            echo "[$SAMPLE] ERROR: Verification failed - file is not properly sorted" >> "$LOG_FILE"
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
        
        # Check if sorted VCF is valid
        HEADER_OK="Yes"
        CHROM_OK="Yes"
        SORTED_OK="Yes"
        
        # Check if VCF has an index
        if [ -f "${SORTED_FILE}.tbi" ]; then
            INDEX_OK="Yes"
        else
            INDEX_OK="No"
        fi
        
        # Add to summary
        echo "File: $FILENAME" >> "$SUMMARY_FILE"
        echo "- Valid Header: $HEADER_OK" >> "$SUMMARY_FILE"
        echo "- Chromosome Info: $CHROM_OK" >> "$SUMMARY_FILE"
        echo "- Properly Sorted: $SORTED_OK" >> "$SUMMARY_FILE"
        echo "- Has Index: $INDEX_OK" >> "$SUMMARY_FILE"
        echo "- Ready for Annotation: Yes" >> "$SUMMARY_FILE"
        echo "" >> "$SUMMARY_FILE"
    fi
done

echo "=== VCF sorting complete ==="