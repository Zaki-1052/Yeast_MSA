#!/bin/bash

# File: scripts/annotation/02b_debug_vcf_sorting.sh
# Purpose: Debug why VCF sorting verification is failing

echo "=== Debugging VCF Sorting Issues ==="
echo "Date: $(date)"
echo ""

# Define test file (use first file for debugging)
VCF_SOURCE="annotation/vcf_ready"
DEBUG_DIR="annotation/debug"
mkdir -p "$DEBUG_DIR"

# Get first sorted VCF file for debugging
TEST_FILE=$(find "$VCF_SOURCE" -name "*.sorted.vcf.gz" | head -n 1)
if [ -z "$TEST_FILE" ]; then
    echo "ERROR: No sorted VCF files found!"
    exit 1
fi

FILENAME=$(basename "$TEST_FILE")
SAMPLE=${FILENAME%.sorted.vcf.gz}
echo "Using $FILENAME for debugging..."
echo ""

# Extract chromosome information from header
echo "Extracting contig information from header..."
bcftools view -h "$TEST_FILE" | grep "##contig=" > "$DEBUG_DIR/contig_info.txt"
CONTIG_COUNT=$(wc -l < "$DEBUG_DIR/contig_info.txt")
echo "Found $CONTIG_COUNT contigs in header"
echo ""

# Extract first few chromosome names from header
echo "First 5 contigs from header:"
head -n 5 "$DEBUG_DIR/contig_info.txt"
echo "..."
echo ""

# Extract chromosome and position data
echo "Extracting chromosome and position data..."
bcftools query -f '%CHROM\t%POS\n' "$TEST_FILE" > "$DEBUG_DIR/chrom_pos_original.txt"
VARIANT_COUNT=$(wc -l < "$DEBUG_DIR/chrom_pos_original.txt")
echo "File contains $VARIANT_COUNT variants"
echo ""

# Sort data using Unix sort
echo "Sorting data with Unix sort..."
sort -k1,1 -k2,2n "$DEBUG_DIR/chrom_pos_original.txt" > "$DEBUG_DIR/chrom_pos_unixsort.txt"

# Check if already sorted
if diff -q "$DEBUG_DIR/chrom_pos_original.txt" "$DEBUG_DIR/chrom_pos_unixsort.txt" >/dev/null; then
    echo "✓ The file is already sorted according to Unix sort"
else
    echo "✗ The file is NOT sorted according to Unix sort"
fi
echo ""

# Extract first few lines of original and sorted data
echo "First 10 lines of original data:"
head -n 10 "$DEBUG_DIR/chrom_pos_original.txt"
echo ""
echo "First 10 lines of Unix-sorted data:"
head -n 10 "$DEBUG_DIR/chrom_pos_unixsort.txt"
echo ""

# Get unique chromosomes in order of appearance
echo "Extracting unique chromosome names in order of appearance..."
cut -f1 "$DEBUG_DIR/chrom_pos_original.txt" | uniq > "$DEBUG_DIR/chrom_order_original.txt"
CHROM_COUNT=$(wc -l < "$DEBUG_DIR/chrom_order_original.txt")
echo "File contains $CHROM_COUNT unique chromosomes"
echo ""

echo "First 10 chromosomes in original order:"
head -n 10 "$DEBUG_DIR/chrom_order_original.txt"
echo ""

# Get same list sorted alphabetically
echo "Sorting chromosome names alphabetically..."
sort "$DEBUG_DIR/chrom_order_original.txt" > "$DEBUG_DIR/chrom_order_sorted.txt"

# Check if already in alphabetical order
if diff -q "$DEBUG_DIR/chrom_order_original.txt" "$DEBUG_DIR/chrom_order_sorted.txt" >/dev/null; then
    echo "✓ Chromosomes are already in alphabetical order"
else
    echo "✗ Chromosomes are NOT in alphabetical order"
fi
echo ""

echo "First 10 chromosomes in alphabetical order:"
head -n 10 "$DEBUG_DIR/chrom_order_sorted.txt"
echo ""

# Test sorting with bedtools
echo "Checking if bedtools is available for alternative sorting..."
if command -v bedtools &> /dev/null; then
    echo "✓ bedtools is available, using for alternative sorting"
    bedtools sort -i "$DEBUG_DIR/chrom_pos_original.txt" > "$DEBUG_DIR/chrom_pos_bedtools.txt"
    
    if diff -q "$DEBUG_DIR/chrom_pos_original.txt" "$DEBUG_DIR/chrom_pos_bedtools.txt" >/dev/null; then
        echo "✓ The file is already sorted according to bedtools sort"
    else
        echo "✗ The file is NOT sorted according to bedtools sort"
    fi
else
    echo "✗ bedtools not found, skipping alternative sorting"
fi
echo ""

echo "Debug information saved to $DEBUG_DIR directory"
echo ""
echo "Conclusion:"
echo "============"
echo "The most likely issue is that bcftools is sorting based on the contig order in the header,"
echo "while our verification test is sorting alphabetically."
echo ""
echo "To fix this, we should either:"
echo "1. Extract the contig order from the header and use it to verify sorting"
echo "2. Skip the verification step and trust bcftools' sorting"
echo "3. Use bcftools itself to verify the sorting"
echo ""
echo "=== Debugging complete ==="