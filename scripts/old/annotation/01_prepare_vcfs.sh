#!/bin/bash

# File: scripts/annotation/01_prepare_vcfs.sh
# Purpose: Prepare and validate VCF files for SnpEff annotation

echo "=== Preparing VCF files for annotation ==="
echo "Date: $(date)"
echo ""

# Create output directory for prepared files
mkdir -p annotation/vcf_ready

# Define source directory for filtered VCFs
VCF_SOURCE="vcf/merged/filtered"
VCF_DEST="annotation/vcf_ready"

# Check if source directory exists
if [ ! -d "$VCF_SOURCE" ]; then
    echo "ERROR: Source directory $VCF_SOURCE not found!"
    exit 1
fi

# Find all filtered VCF files
echo "Scanning for filtered VCF files..."
VCF_FILES=$(find "$VCF_SOURCE" -name "*.filtered.vcf.gz")
VCF_COUNT=$(echo "$VCF_FILES" | wc -l)

if [ "$VCF_COUNT" -eq 0 ]; then
    echo "ERROR: No filtered VCF files found in $VCF_SOURCE"
    exit 1
fi

echo "Found $VCF_COUNT VCF files to process"
echo ""

# Create a validation summary file
SUMMARY_FILE="$VCF_DEST/validation_summary.txt"
echo "VCF File Validation Summary" > "$SUMMARY_FILE"
echo "Date: $(date)" >> "$SUMMARY_FILE"
echo "==============================================" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"

# Process each VCF file
for VCF_FILE in $VCF_FILES; do
    FILENAME=$(basename "$VCF_FILE")
    SAMPLE=${FILENAME%.filtered.vcf.gz}
    
    echo "Processing $SAMPLE..."
    echo "- File: $VCF_FILE"
    
    # Check if VCF is valid
    echo "- Validating VCF format..."
    if bcftools view -h "$VCF_FILE" &>/dev/null; then
        echo "  ✓ VCF header is valid"
        HEADER_OK="Yes"
    else
        echo "  ✗ VCF header is invalid"
        HEADER_OK="No"
    fi
    
    # Check if VCF has chromosome information
    echo "- Checking chromosome information..."
    CHROM_COUNT=$(bcftools view -h "$VCF_FILE" | grep -c "##contig=")
    if [ "$CHROM_COUNT" -gt 0 ]; then
        echo "  ✓ Found $CHROM_COUNT chromosome/contig definitions"
        CHROM_OK="Yes"
    else
        echo "  ✗ No chromosome/contig definitions found"
        CHROM_OK="No"
    fi
    
    # Check if VCF is sorted
    echo "- Checking if VCF is sorted..."
    if bcftools query -f '%CHROM\t%POS\n' "$VCF_FILE" | sort -k1,1 -k2,2n | \
       diff -q - <(bcftools query -f '%CHROM\t%POS\n' "$VCF_FILE") &>/dev/null; then
        echo "  ✓ VCF is properly sorted"
        SORTED_OK="Yes"
    else
        echo "  ✗ VCF is not properly sorted"
        SORTED_OK="No"
    fi
    
    # Check if VCF has an index
    echo "- Checking for index file..."
    if [ -f "${VCF_FILE}.tbi" ]; then
        echo "  ✓ Index file exists"
        INDEX_OK="Yes"
    else
        echo "  ✗ No index file found"
        echo "  - Creating index..."
        tabix -p vcf "$VCF_FILE"
        if [ $? -eq 0 ]; then
            echo "  ✓ Index created successfully"
            INDEX_OK="Yes (Created)"
        else
            echo "  ✗ Failed to create index"
            INDEX_OK="No (Failed)"
        fi
    fi
    
    # Copy file to destination if everything is OK
    if [ "$HEADER_OK" = "Yes" ] && [ "$CHROM_OK" = "Yes" ] && [ "$SORTED_OK" = "Yes" ]; then
        echo "- Copying validated file to annotation directory..."
        cp "$VCF_FILE" "$VCF_DEST/"
        cp "${VCF_FILE}.tbi" "$VCF_DEST/" 2>/dev/null || true
        COPY_OK="Yes"
    else
        echo "- File needs repair before annotation"
        COPY_OK="No"
    fi
    
    # Add to summary
    echo "File: $FILENAME" >> "$SUMMARY_FILE"
    echo "- Valid Header: $HEADER_OK" >> "$SUMMARY_FILE"
    echo "- Chromosome Info: $CHROM_OK" >> "$SUMMARY_FILE"
    echo "- Properly Sorted: $SORTED_OK" >> "$SUMMARY_FILE"
    echo "- Has Index: $INDEX_OK" >> "$SUMMARY_FILE"
    echo "- Ready for Annotation: $COPY_OK" >> "$SUMMARY_FILE"
    echo "" >> "$SUMMARY_FILE"
    
    echo "Done processing $SAMPLE"
    echo ""
done

# Count ready files
READY_COUNT=$(grep -c "Ready for Annotation: Yes" "$SUMMARY_FILE")
echo "Summary: $READY_COUNT of $VCF_COUNT files are ready for annotation"
echo "Validation summary saved to: $SUMMARY_FILE"
echo ""
echo "=== VCF preparation complete ==="