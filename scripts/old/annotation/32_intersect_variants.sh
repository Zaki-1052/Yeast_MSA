#!/bin/bash

# File: scripts/annotation/32_intersect_variants.sh
# Purpose: Intersect renamed VCFs with target gene BED file

# Exit immediately if a command exits with a non-zero status.
set -e

echo "=== Intersecting Variants with Target Gene Regions ==="
echo "Date: $(date)"
echo ""

# Define directories
VCF_DIR="annotation/vcf_fixed_direct" # Use the VCFs with renamed chromosomes (w303_scaffold_XXX)
BED_FILE="annotation/gene_coordinates/target_genes_robust.bed" # Use the newly created robust BED file
OUTPUT_DIR="annotation/gene_variants_intersected"

# Create output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/by_gene"

# Check if BED file exists
if [ ! -f "$BED_FILE" ]; then
    echo "ERROR: BED file not found: $BED_FILE"
    exit 1
fi

# Find VCF files
# Use find directly in bash for better handling of paths and spaces
VCF_FILES=$(find "$VCF_DIR" -name '*.fixed.vcf.gz')
if [ -z "$VCF_FILES" ]; then
    echo "ERROR: No fixed VCF files found in $VCF_DIR"
    exit 1
fi

VCF_COUNT=$(echo "$VCF_FILES" | wc -l | xargs) # Get count properly
echo "Found $VCF_COUNT VCF files to process."

# Process each VCF file
COMBINED_OUTPUT="$OUTPUT_DIR/all_target_variants.vcf"

# Initialize combined output file with header from the first VCF
FIRST_VCF=$(echo "$VCF_FILES" | head -n 1)
if [ -z "$FIRST_VCF" ]; then echo "ERROR: Could not determine first VCF file."; exit 1; fi

echo "Initializing combined output file with header from $FIRST_VCF..."
# Ensure header is written correctly, overwriting if file exists
bcftools view -h "$FIRST_VCF" > "$COMBINED_OUTPUT"
if [ $? -ne 0 ]; then echo "ERROR: Failed to get header from $FIRST_VCF"; exit 1; fi

echo "Processing VCF files..."
SUCCESS_COUNT=0
FAILURE_COUNT=0
for VCF_FILE in $VCF_FILES; do
    SAMPLE=$(basename "$VCF_FILE" .fixed.vcf.gz)
    echo -n "Processing $SAMPLE... "

    # Intersect using bcftools view -R and append data
    # Use temporary file to check if variants were found before appending
    TEMP_VARIANTS=$(mktemp)
    bcftools view -H -R "$BED_FILE" "$VCF_FILE" > "$TEMP_VARIANTS"
    BCF_EXIT_CODE=$?

    if [ $BCF_EXIT_CODE -ne 0 ] && [ $BCF_EXIT_CODE -ne 1 ]; then
        echo "ERROR: bcftools failed for $SAMPLE with exit code $BCF_EXIT_CODE"
        FAILURE_COUNT=$((FAILURE_COUNT + 1))
    else
        # Check if any variants were actually found
        if [ -s "$TEMP_VARIANTS" ]; then
            cat "$TEMP_VARIANTS" >> "$COMBINED_OUTPUT"
            echo "Done (variants found)."
            SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
        else
            echo "Done (no variants in regions)."
            SUCCESS_COUNT=$((SUCCESS_COUNT + 1)) # Count as success even if no variants found
        fi
    fi
    rm -f "$TEMP_VARIANTS" # Clean up temp file
done

echo "Finished processing VCFs. Success: $SUCCESS_COUNT, Failures: $FAILURE_COUNT."

# Check if combined file has variants before proceeding
VARIANT_LINES=$(grep -vc '^#' "$COMBINED_OUTPUT")
if [ "$VARIANT_LINES" -eq 0 ]; then
    echo "No variants found in target regions across all samples. Skipping compression and splitting."
    # Optionally remove the empty combined file: rm -f "$COMBINED_OUTPUT"
else
    echo "Combined intersected variants saved to: $COMBINED_OUTPUT"
    echo "Compressing and indexing combined file..."
    bgzip -f "$COMBINED_OUTPUT"
    tabix -p vcf "${COMBINED_OUTPUT}.gz"
    COMBINED_GZ="${COMBINED_OUTPUT}.gz" # Reference the gzipped file

    # Split combined file by gene based on BED coordinates
    echo "Splitting combined file by gene based on coordinates..."
    # Extract target gene names from BED file, ensure unique names
    TARGET_GENES=$(tail -n +2 "$BED_FILE" | cut -f4 | sort | uniq)
    TOTAL_SPLIT_VARIANTS=0
    for GENE in $TARGET_GENES; do
        echo "  Extracting variants for $GENE..."
        GENE_OUTPUT_VCF="$OUTPUT_DIR/by_gene/${GENE}_variants.vcf"
        GENE_OUTPUT_GZ="${GENE_OUTPUT_VCF}.gz"

        # Get region for this gene from BED file
        # Use awk for robustness and handle potential multiple entries for a gene (take first)
        REGION=$(awk -v gene="$GENE" '$4 == gene {{printf "%s:%s-%s", $1, $2+1, $3; exit}}' "$BED_FILE")
        if [ -z "$REGION" ]; then
            echo "    ERROR: Could not find region for $GENE in BED file."
            continue
        fi
        echo "    Region: $REGION"

        # Extract variants for this specific region, include header
        # Handle potential errors if no variants are found in the region
        { bcftools view -h "$COMBINED_GZ"; bcftools view -H -r "$REGION" "$COMBINED_GZ"; } > "$GENE_OUTPUT_VCF" 2>/dev/null
        COUNT=$(grep -cv '^#' "$GENE_OUTPUT_VCF") # Use grep -c -v
        echo "    Found $COUNT variants in region for $GENE"
        # Only keep file if variants were found
        if [ "$COUNT" -gt 0 ]; then
            TOTAL_SPLIT_VARIANTS=$((TOTAL_SPLIT_VARIANTS + COUNT))
            bgzip -f "$GENE_OUTPUT_VCF"
            tabix -p vcf "$GENE_OUTPUT_GZ"
        else
            rm -f "$GENE_OUTPUT_VCF" # Remove empty file
        fi
    done

    echo "Total variants assigned to specific gene regions: $TOTAL_SPLIT_VARIANTS"
    echo "Combined results: ${COMBINED_GZ}"
    echo "Results split by gene coordinates: $OUTPUT_DIR/by_gene/"
    echo "NOTE: Review the combined file for any variants missed by the split."

fi # End check for variants in combined file

echo -e "\n=== Variant Intersection Complete ===" # Use -e for newline