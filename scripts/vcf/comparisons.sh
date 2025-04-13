#!/bin/bash

# Create directories for results
mkdir -p results/vcf/direct_comparison

# Function to compare VCFs from treatment samples against a control
compare_treatment_to_control() {
    GROUP=$1
    CONTROL_VCF=$2
    TREATMENT_VCF_PREFIX=$3
    
    echo "Comparing $GROUP treatments vs control using direct VCF comparison"
    
    # Create output directory
    OUT_DIR="results/vcf/direct_comparison/${GROUP}"
    mkdir -p "$OUT_DIR"
    
    # Find treatment VCF files
    TREATMENT_VCFS=$(ls results/vcf/filtered/${TREATMENT_VCF_PREFIX}*.filtered.vcf.gz)
    
    # Compare each treatment sample with control
    for TREAT_VCF in $TREATMENT_VCFS; do
        SAMPLE=$(basename "$TREAT_VCF" .filtered.vcf.gz)
        echo "  Processing $SAMPLE"
        
        # Create a directory for this comparison
        ISEC_DIR="${OUT_DIR}/${SAMPLE}_vs_control"
        mkdir -p "$ISEC_DIR"
        
        # Find variants present in treatment but not in control
        bcftools isec -p "$ISEC_DIR" \
                     -n =1 -c none \
                     "$TREAT_VCF" "$CONTROL_VCF"
        
        # Count treatment-specific variants
        if [ -f "${ISEC_DIR}/0000.vcf" ]; then
            COUNT=$(grep -v "^#" "${ISEC_DIR}/0000.vcf" | wc -l)
            echo "  $SAMPLE specific variants: $COUNT"
            
            # Copy the result file to a clearer name
            cp "${ISEC_DIR}/0000.vcf" "${OUT_DIR}/${SAMPLE}_specific.vcf"
        else
            echo "  No treatment-specific variants found for $SAMPLE"
        fi
    done
    
    # Now find high-confidence variants (present in at least 2 of 3 replicates)
    if ls "${OUT_DIR}/"*_specific.vcf > /dev/null 2>&1; then
        echo "  Finding high-confidence variants for $GROUP"
        
        # Create a list of treatment-specific variant files
        SPECIFIC_VCFS="${OUT_DIR}/"*_specific.vcf
        
        # Use a more direct approach - extract positions and count occurrences
        cat $SPECIFIC_VCFS | grep -v "^#" | awk '{print $1":"$2}' | sort | uniq -c | sort -nr > "${OUT_DIR}/variant_counts.txt"
        
        # Extract variants that appear in at least 2 samples
        HIGH_CONF=$(cat "${OUT_DIR}/variant_counts.txt" | awk '$1 >= 2 {print}' | wc -l)
        echo "$GROUP high-confidence treatment-specific variants: $HIGH_CONF"
    else
        echo "  No treatment-specific variants to analyze for $GROUP"
    fi
}

# Proceed with direct comparisons
# WT group
compare_treatment_to_control "WT" "results/vcf/filtered/WT-CTRL.filtered.vcf.gz" "WT-37-55"

# STC group
compare_treatment_to_control "STC" "results/vcf/filtered/STC-CTRL.filtered.vcf.gz" "STC-55"

# CAS group
compare_treatment_to_control "CAS" "results/vcf/filtered/CAS-CTRL.filtered.vcf.gz" "CAS-55"

# WTA group (compare to WT control)
compare_treatment_to_control "WTA" "results/vcf/filtered/WT-CTRL.filtered.vcf.gz" "WTA-55"

echo "Direct VCF comparison complete."