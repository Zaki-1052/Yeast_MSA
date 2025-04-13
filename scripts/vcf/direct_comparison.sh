#!/bin/bash

# Create directories for results
mkdir -p results/vcf/final_comparison

# Function to compare VCFs from treatment samples against a control
compare_direct() {
    GROUP=$1
    CONTROL_VCF=$2
    TREATMENT_VCF_PREFIX=$3
    
    echo "Comparing $GROUP treatments vs control using direct VCF comparison"
    
    # Create output directory
    OUT_DIR="results/vcf/final_comparison/${GROUP}"
    mkdir -p "$OUT_DIR"
    
    # Find treatment VCF files
    TREATMENT_VCFS=$(ls results/vcf/filtered/${TREATMENT_VCF_PREFIX}*.filtered.vcf.gz)
    
    # Compare each treatment sample with control
    for TREAT_VCF in $TREATMENT_VCFS; do
        SAMPLE=$(basename "$TREAT_VCF" .filtered.vcf.gz)
        echo "  Processing $SAMPLE"
        
        # Find variants present in treatment but not in control
        bcftools isec -p "${OUT_DIR}/${SAMPLE}_vs_control" \
                     -n =1 -c any \
                     "$TREAT_VCF" "$CONTROL_VCF"
        
        # Count treatment-specific variants
        if [ -f "${OUT_DIR}/${SAMPLE}_vs_control/0000.vcf" ]; then
            COUNT=$(grep -v "^#" "${OUT_DIR}/${SAMPLE}_vs_control/0000.vcf" | wc -l)
            echo "  $SAMPLE specific variants: $COUNT"
            cp "${OUT_DIR}/${SAMPLE}_vs_control/0000.vcf" "${OUT_DIR}/${SAMPLE}_specific.vcf"
        else
            echo "  No treatment-specific variants found for $SAMPLE"
        fi
    done
    
    # Find high-confidence variants (in multiple replicates)
    find_high_confidence_variants "$GROUP" "$OUT_DIR"
}

# Extract high-confidence variants
find_high_confidence_variants() {
    GROUP=$1
    OUT_DIR=$2
    
    echo "  Finding high-confidence variants for $GROUP"
    
    # Process all variant files to find overlaps
    cat ${OUT_DIR}/*_specific.vcf | grep -v "^#" | \
      awk '{print $1"\t"$2"\t"$4"\t"$5}' | sort | uniq -c | sort -nr > "${OUT_DIR}/variant_counts.txt"
    
    # Count high-confidence variants (≥2 replicates)
    HIGH_CONF=$(awk '$1>=2' "${OUT_DIR}/variant_counts.txt" | wc -l)
    echo "$GROUP high-confidence treatment-specific variants: $HIGH_CONF"
    
    # Create detailed high-confidence list
    echo -e "Count\tCHROM\tPOS\tREF\tALT" > "${OUT_DIR}/highconf_variants.txt"
    awk '$1>=2 {print $0}' "${OUT_DIR}/variant_counts.txt" >> "${OUT_DIR}/highconf_variants.txt"
}

# Run comparisons
compare_direct "WT" "results/vcf/filtered/WT-CTRL.filtered.vcf.gz" "WT-37-55"
compare_direct "STC" "results/vcf/filtered/STC-CTRL.filtered.vcf.gz" "STC-55"
compare_direct "CAS" "results/vcf/filtered/CAS-CTRL.filtered.vcf.gz" "CAS-55"
compare_direct "WTA" "results/vcf/filtered/WT-CTRL.filtered.vcf.gz" "WTA-55"

echo "Direct comparison complete!"