#!/bin/bash

# File: scripts/annotation/30_verify_identifier_presence.sh
# Purpose: Directly verify the presence of gene identifiers in VCF files using shell commands.

# Exit immediately if a command exits with a non-zero status.
set -e

echo "=== Verifying Identifier Presence in VCF Files ==="
echo "Date: $(date)"
echo ""

# Define directories where annotated files SHOULD be
# Prioritize the ones potentially created by the chromosome renaming step
VCF_DIRS=(
    "annotation/vcf_fixed_direct"  # Output of script 22 (after renaming)
    "annotation/results_fixed"     # Output of script 23 (SnpEff on fixed VCFs)
    "annotation/results_renamed"   # Output of script 11 (SnpEff on maybe-renamed VCFs)
    "annotation/results"           # Output of script 03 (original SnpEff run)
    "annotation/vcf_ready"         # Output of script 02c (sorted VCFs before renaming)
)

# Define a few specific identifiers to search for
SEARCH_TERMS=(
    "YHR007C"           # SGD ID
    "W303_0EI00110"     # W303 ID (from mapping)
    "W3030EI00110"      # Gene ID (from mapping)
    "ERG11"             # Common Name (from mapping)
    "w303_scaffold_139" # Scaffold associated with YHR007C
)

# Create a log file
LOG_FILE="annotation/identifier_verification.log"
echo "Identifier Verification Log" > "$LOG_FILE"
echo "Date: $(date)" >> "$LOG_FILE"
echo "==============================================" >> "$LOG_FILE"
echo "" >> "$LOG_FILE"

FOUND_ANY=0

# Loop through each potential directory
for VCF_DIR in "${VCF_DIRS[@]}"; do
    echo ">>> Checking directory: $VCF_DIR" | tee -a "$LOG_FILE"
    if [ ! -d "$VCF_DIR" ]; then
        echo "    Directory does not exist. Skipping." | tee -a "$LOG_FILE"
        echo "" | tee -a "$LOG_FILE"
        continue
    fi

    # Find VCF files in this directory
    VCF_FILES=$(find "$VCF_DIR" -name '*.vcf.gz' -print -quit) # Find only the first VCF file to test

    if [ -z "$VCF_FILES" ]; then
        echo "    No VCF files found in this directory." | tee -a "$LOG_FILE"
        echo "" | tee -a "$LOG_FILE"
        continue
    fi

    # Use the first found VCF file for testing this directory
    TEST_VCF="$VCF_FILES"
    echo "    Using sample VCF: $TEST_VCF" | tee -a "$LOG_FILE"

    # Test access method again - prioritize gunzip, then bcftools, then zcat
    ACCESS_CMD=""
    echo "    Testing access methods..." | tee -a "$LOG_FILE"
    if gunzip -c "$TEST_VCF" > /dev/null 2>&1; then
        ACCESS_CMD="gunzip -c"
        echo "      Using gunzip -c" | tee -a "$LOG_FILE"
    elif bcftools view "$TEST_VCF" > /dev/null 2>&1; then
        ACCESS_CMD="bcftools view"
         echo "      Using bcftools view" | tee -a "$LOG_FILE"
    elif zcat "$TEST_VCF" > /dev/null 2>&1; then
         ACCESS_CMD="zcat"
         echo "      Using zcat (might be problematic)" | tee -a "$LOG_FILE"
    else
         echo "      ERROR: Cannot access VCF content with gunzip, bcftools, or zcat." | tee -a "$LOG_FILE"
         continue # Skip to next directory if we can't access the file
    fi


    # Search for each term in the test VCF
    echo "    Searching for identifiers..." | tee -a "$LOG_FILE"
    for TERM in "${SEARCH_TERMS[@]}"; do
        echo "      Searching for '$TERM'..." | tee -a "$LOG_FILE"
        # Use grep -c for count, add -w for word boundaries
        COUNT=$($ACCESS_CMD "$TEST_VCF" | grep -c -w "$TERM" || true) # Allow grep to return 1 (no match) without exiting script

        if [ "$COUNT" -gt 0 ]; then
            echo "        ✓ Found $COUNT occurrences of '$TERM'" | tee -a "$LOG_FILE"
            echo "          Example lines:" | tee -a "$LOG_FILE"
            # Show first 3 matching lines
            EXAMPLES=$($ACCESS_CMD "$TEST_VCF" | grep -w "$TERM" | head -n 3)
            echo "$EXAMPLES" | sed 's/^/            /' | tee -a "$LOG_FILE"
            FOUND_ANY=1
        else
            echo "        ✗ No occurrences found for '$TERM'" | tee -a "$LOG_FILE"
        fi
    done
    echo "" | tee -a "$LOG_FILE"
done

echo "=== Verification Summary ===" | tee -a "$LOG_FILE"
if [ "$FOUND_ANY" -eq 1 ]; then
    echo "✓ At least one identifier was found in one of the checked VCF files." | tee -a "$LOG_FILE"
    echo "  Please review the log file ($LOG_FILE) for details." | tee -a "$LOG_FILE"
    echo "  The issue likely lies in how the Python search script (29) executes the search." | tee -a "$LOG_FILE"
else
    echo "✗ None of the specific identifiers were found in the first VCF file of any checked directory." | tee -a "$LOG_FILE"
    echo "  This suggests the gene identifiers are truly absent or in an unexpected format in the VCF files." | tee -a "$LOG_FILE"
    echo "  Possible causes: Annotation failed to add gene info, or file corruption." | tee -a "$LOG_FILE"
    echo "  Next steps: Manually inspect VCF files, re-run annotation carefully checking logs." | tee -a "$LOG_FILE"
fi

echo "" | tee -a "$LOG_FILE"
echo "=== Verification Complete ==="
