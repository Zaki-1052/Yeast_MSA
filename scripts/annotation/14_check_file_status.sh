#!/bin/bash

# File: scripts/annotation/14_check_file_status.sh
# Purpose: Check status of files after annotation attempts

echo "=== Checking File Status ==="
echo "Date: $(date)"
echo ""

# Create diagnosis directory
mkdir -p annotation/diagnosis

# Check directory structure and file existence
echo "1. Checking directory structure..."
echo "Renamed VCF files:"
ls -l annotation/vcf_renamed/
echo ""

echo "Annotation results directory:"
ls -l annotation/results_renamed/ 2>/dev/null || echo "  Directory empty or does not exist"
echo ""

echo "Statistics directory:"
ls -l annotation/stats_renamed/ 2>/dev/null || echo "  Directory empty or does not exist"
echo ""

# Check for annotation logs
echo "2. Checking for annotation logs..."
if [ -f "annotation/snpeff_annotation.log" ]; then
    echo "Found annotation log:"
    tail -n 20 annotation/snpeff_annotation.log
else
    echo "No annotation log found"
fi
echo ""

# Check SnpEff database
echo "3. Checking SnpEff database..."
SNPEFF_DIR="/Users/zakiralibhai/snpEff"
DATA_DIR="${SNPEFF_DIR}/data/w303"

echo "Data directory contents:"
ls -l "$DATA_DIR" 2>/dev/null || echo "  Directory empty or does not exist"
echo ""

echo "Looking for snpEffectPredictor.bin:"
if [ -f "${DATA_DIR}/snpEffectPredictor.bin" ]; then
    echo "  ✓ File exists ($(stat -f%z "${DATA_DIR}/snpEffectPredictor.bin") bytes)"
else
    echo "  ✗ File does not exist"
fi
echo ""

# Run basic database dump
echo "4. Testing basic SnpEff database dump..."
cd "$SNPEFF_DIR"
echo "Running: java -jar snpEff.jar databases | grep -i w303"
java -jar snpEff.jar databases | grep -i w303
echo ""

echo "Running: java -jar snpEff.jar dump w303 | head -20"
java -jar snpEff.jar dump w303 | head -20
cd - > /dev/null

# Run a very basic test
echo ""
echo "5. Running minimal test annotation..."
# Create a tiny test VCF
TEST_DIR="annotation/diagnosis"
TEST_VCF="${TEST_DIR}/tiny_test.vcf"

cat > "$TEST_VCF" << 'VCFEND'
##fileformat=VCFv4.2
##contig=<ID=w303_scaffold_1>
##reference=w303_scaffold_1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
w303_scaffold_1	100	.	A	T	100	PASS	.
VCFEND

echo "Created tiny test VCF with a single variant"

# Run SnpEff on this minimal file
java -jar "${SNPEFF_DIR}/snpEff.jar" -v w303 "$TEST_VCF" > "${TEST_DIR}/tiny_test_annotated.vcf"

echo "Test annotation result:"
cat "${TEST_DIR}/tiny_test_annotated.vcf"

echo ""
echo "=== File Status Check Complete ==="
