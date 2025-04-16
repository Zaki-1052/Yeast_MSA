#!/bin/bash

# File: scripts/annotation/09_test_annotation.sh
# Purpose: Test SnpEff annotation with the rebuilt database

echo "=== Testing SnpEff Annotation ==="
echo "Date: $(date)"
echo ""

# Define directories
SNPEFF_DIR="/Users/zakiralibhai/snpEff"
VCF_SOURCE="annotation/vcf_ready"
TEST_DIR="annotation/test_annotation"

# Create test directory
mkdir -p "$TEST_DIR"

# Select first VCF file for testing
TEST_VCF=$(find "$VCF_SOURCE" -name "*.sorted.vcf.gz" | head -n 1)
if [ -z "$TEST_VCF" ]; then
    echo "ERROR: No VCF files found for testing"
    exit 1
fi

TEST_SAMPLE=$(basename "$TEST_VCF" .sorted.vcf.gz)
echo "Selected test file: $TEST_SAMPLE"
echo ""

# Run SnpEff annotation
echo "Running SnpEff annotation..."
TEST_OUTPUT="${TEST_DIR}/${TEST_SAMPLE}.test.vcf"
TEST_STATS="${TEST_DIR}/${TEST_SAMPLE}.test.stats.html"

java -Xmx4g -jar "${SNPEFF_DIR}/snpEff.jar" \
    -v \
    -stats "$TEST_STATS" \
    w303 \
    "$TEST_VCF" \
    > "$TEST_OUTPUT"

# Check for errors
ERROR_COUNT=$(grep -c "ERROR_CHROMOSOME_NOT_FOUND" "$TEST_OUTPUT")
if [ "$ERROR_COUNT" -gt 0 ]; then
    echo "❌ Test failed: Found $ERROR_COUNT chromosome errors"
    echo "First few errors:"
    grep "ERROR_CHROMOSOME_NOT_FOUND" "$TEST_OUTPUT" | head -5
    
    # Check chromosome names in the file
    echo ""
    echo "Checking chromosome names in VCF..."
    CHROMS=$(bcftools view -h "$TEST_VCF" | grep '##contig=' | head -5)
    echo "$CHROMS"
    
    # Check chromosome names in database
    echo ""
    echo "Checking available chromosomes in SnpEff database..."
    cd "$SNPEFF_DIR"
    CHROMS_DB=$(java -jar snpEff.jar dump w303 | grep 'Chromosomes' -A 10)
    echo "$CHROMS_DB"
    cd - > /dev/null
    
    echo ""
    echo "Possible issues:"
    echo "1. The chromosome names in VCF still don't match the database"
    echo "2. The database build might not have processed all chromosomes"
    echo ""
    echo "Next steps:"
    echo "1. Try modifying the VCF chromosome names instead"
    echo "2. Check the SnpEff build logs for more details"
else
    echo "✅ Test successful: No chromosome errors found"
    
    # Check for annotations
    ANN_COUNT=$(grep -c "ANN=" "$TEST_OUTPUT")
    if [ "$ANN_COUNT" -gt 0 ]; then
        echo "Found $ANN_COUNT variants with annotations"
        echo "First few annotated variants:"
        grep "ANN=" "$TEST_OUTPUT" | head -3
        
        # Compress and index
        bgzip -f "$TEST_OUTPUT"
        tabix -p vcf "${TEST_OUTPUT}.gz"
        
        echo ""
        echo "Next steps:"
        echo "1. Proceed with annotating all VCF files using the same approach"
        echo "2. Re-run the target gene extraction script"
    else
        echo "⚠️ No annotations found. This is unusual."
        echo "The database may not be properly identifying genes."
    fi
fi

echo ""
echo "Test results saved to: $TEST_DIR"
echo ""
echo "=== SnpEff Annotation Test Complete ==="
