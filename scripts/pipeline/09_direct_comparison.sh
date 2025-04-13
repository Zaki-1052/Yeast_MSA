#!/bin/bash

echo "=== Performing direct VCF comparisons ==="

# Load the force flag if needed
if [ -f "results/merged/fixed/force_flag.txt" ]; then
    source results/merged/fixed/force_flag.txt
else
    FORCE_FLAG=""
fi

# Function for direct comparison of VCFs
compare_treatment_to_control_direct() {
    GROUP=$1
    CONTROL_VCF=$2
    TREATMENT_VCF_PREFIX=$3
    
    echo "  Direct comparison for $GROUP treatments vs control"
    
    # Create output directory
    OUT_DIR="results/merged/direct_comparison/${GROUP}"
    mkdir -p "$OUT_DIR"
    
    # Find treatment VCF files
    TREATMENT_VCFS=$(ls results/merged/fixed/${TREATMENT_VCF_PREFIX}*.fixed.vcf.gz 2>/dev/null)
    
    # Check if we found any files
    if [ -z "$TREATMENT_VCFS" ]; then
        echo "    No treatment VCF files found for pattern: results/merged/fixed/${TREATMENT_VCF_PREFIX}*.fixed.vcf.gz"
        mkdir -p "results/merged/summary/${GROUP}"
        echo "0" > "results/merged/summary/${GROUP}/direct_highconf_count.txt"
        return
    fi
    
    # Correct the control VCF path if needed
    if [ ! -f "$CONTROL_VCF" ]; then
        # Try to find the control in the fixed directory
        CONTROL_BASE=$(basename "$CONTROL_VCF" .filtered.vcf.gz)
        if [ -f "results/merged/fixed/${CONTROL_BASE}.fixed.vcf.gz" ]; then
            CONTROL_VCF="results/merged/fixed/${CONTROL_BASE}.fixed.vcf.gz"
        else
            echo "    Control VCF not found: $CONTROL_VCF"
            echo "    Skipping direct comparison for $GROUP"
            mkdir -p "results/merged/summary/${GROUP}"
            echo "0" > "results/merged/summary/${GROUP}/direct_highconf_count.txt"
            return
        fi
    fi
    
    echo "    Using control: $CONTROL_VCF"
    echo "    Using treatment files: $TREATMENT_VCFS"
    
    # Compare each treatment sample with control
    for TREAT_VCF in $TREATMENT_VCFS; do
        SAMPLE=$(basename "$TREAT_VCF" | sed 's/\.fixed\.vcf\.gz$//')
        echo "    Processing $SAMPLE"
        
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
            echo "    $SAMPLE specific variants: $COUNT"
            
            # Copy the result file to a clearer name
            cp "${ISEC_DIR}/0000.vcf" "${OUT_DIR}/${SAMPLE}_specific.vcf"
        else
            echo "    No treatment-specific variants found for $SAMPLE"
        fi
    done
    
    # Find high-confidence variants (present in at least 2 of 3 replicates)
    if ls "${OUT_DIR}/"*_specific.vcf > /dev/null 2>&1; then
        echo "    Finding high-confidence variants for $GROUP"
        
        # Extract positions and count occurrences
        cat "${OUT_DIR}/"*_specific.vcf | grep -v "^#" | awk '{print $1":"$2":"$4":"$5}' | sort | uniq -c | sort -nr > "${OUT_DIR}/variant_counts.txt"
        
        # Extract variants that appear in at least 2 samples
        HIGH_CONF=$(awk '$1 >= 2' "${OUT_DIR}/variant_counts.txt" | wc -l)
        echo "    $GROUP high-confidence direct-comparison variants: $HIGH_CONF"
        
        # Save to summary
        mkdir -p "results/merged/summary/${GROUP}"
        echo "$HIGH_CONF" > "results/merged/summary/${GROUP}/direct_highconf_count.txt"
    else
        echo "    No treatment-specific variants to analyze for $GROUP"
        mkdir -p "results/merged/summary/${GROUP}"
        echo "0" > "results/merged/summary/${GROUP}/direct_highconf_count.txt"
    fi
}

# Run direct comparisons with biologically correct groupings
compare_treatment_to_control_direct "WT-37" "results/merged/fixed/WT-CTRL.fixed.vcf.gz" "WT-37-55"
compare_treatment_to_control_direct "WTA" "results/merged/fixed/WT-CTRL.fixed.vcf.gz" "WTA-55"
compare_treatment_to_control_direct "STC" "results/merged/fixed/WT-CTRL.fixed.vcf.gz" "STC-55"
compare_treatment_to_control_direct "CAS" "results/merged/fixed/WT-CTRL.fixed.vcf.gz" "CAS-55"

# Also keep original control comparisons
compare_treatment_to_control_direct "STC-vs-STCCTRL" "results/merged/fixed/STC-CTRL.fixed.vcf.gz" "STC-55"
compare_treatment_to_control_direct "CAS-vs-CASCTRL" "results/merged/fixed/CAS-CTRL.fixed.vcf.gz" "CAS-55"