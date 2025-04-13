#!/bin/bash

echo "=== Performing treatment vs control comparisons ==="

# Load the force flag if needed
if [ -f "results/merged/fixed/force_flag.txt" ]; then
    source results/merged/fixed/force_flag.txt
else
    FORCE_FLAG=""
fi

# Function to compare treatment samples with control
compare_treatment_control() {
    local GROUP="$1"
    local CONTROL="$2"
    local TREATMENT_SAMPLES="$3"
    
    echo "  Comparing $GROUP treatments vs control ($CONTROL)..."
    
    # Create output directory
    mkdir -p "results/merged/analysis/${GROUP}"
    
    # Verify variables are defined before using
    if [ -z "$GROUP" ] || [ -z "$CONTROL" ] || [ -z "$TREATMENT_SAMPLES" ]; then
        echo "ERROR: Missing parameters for compare_treatment_control function"
        echo "  GROUP: $GROUP"
        echo "  CONTROL: $CONTROL"
        echo "  TREATMENT_SAMPLES: $TREATMENT_SAMPLES"
        return 1
    fi
    
    # Extract just this treatment group with its control
    echo "    Extracting $GROUP group data..."
    bcftools view $FORCE_FLAG -s "$CONTROL,$TREATMENT_SAMPLES" results/merged/fixed/all_samples.vcf.gz \
      -Oz -o "results/merged/analysis/${GROUP}/group.vcf.gz"
    tabix -p vcf "results/merged/analysis/${GROUP}/group.vcf.gz"
    
    # Get control index
    CONTROL_IDX=$(bcftools query -l "results/merged/analysis/${GROUP}/group.vcf.gz" | grep -n "$CONTROL" | cut -d: -f1)
    CONTROL_IDX=$((CONTROL_IDX - 1))  # Convert to 0-based index
    
    # Get treatment samples indices
    TREATMENTS=$(bcftools query -l "results/merged/analysis/${GROUP}/group.vcf.gz" | grep -v "$CONTROL" | tr '\n' ',')
    TREATMENTS=${TREATMENTS%,}  # Remove trailing comma
    
    # Create a dynamic filter expression for any number of treatment samples
    FILTER_EXPR="GT[$CONTROL_IDX]=\".\""
    for SAMPLE in $(echo $TREATMENTS | tr ',' ' '); do
        SAMPLE_IDX=$(bcftools query -l "results/merged/analysis/${GROUP}/group.vcf.gz" | grep -n "$SAMPLE" | cut -d: -f1)
        SAMPLE_IDX=$((SAMPLE_IDX - 1))
        FILTER_EXPR="$FILTER_EXPR || GT[$SAMPLE_IDX]=\"1\""
    done
    FILTER_EXPR="$FILTER_EXPR)"
    FILTER_EXPR="GT[$CONTROL_IDX]=\".\" && ($FILTER_EXPR"
    
    # Find variants present in at least one treatment sample but not in control
    echo "    Identifying treatment-specific variants..."
    bcftools view -i "$FILTER_EXPR" "results/merged/analysis/${GROUP}/group.vcf.gz" \
      -Oz -o "results/merged/analysis/${GROUP}/specific.vcf.gz"
    tabix -p vcf "results/merged/analysis/${GROUP}/specific.vcf.gz"
    
    # Count treatment-specific variants
    COUNT=$(bcftools view -H "results/merged/analysis/${GROUP}/specific.vcf.gz" | wc -l)
    echo "    Found $COUNT treatment-specific variants"
    
    # Find high-confidence variants (in multiple replicates)
    echo "    Identifying high-confidence variants..."
    
    # Count samples to set threshold
    NUM_SAMPLES=$(echo "$TREATMENT_SAMPLES" | tr -cd ',' | wc -c)
    NUM_SAMPLES=$((NUM_SAMPLES + 1))
    THRESHOLD=$((NUM_SAMPLES - 1))  # n-1 replicates for high confidence
    
    bcftools view $FORCE_FLAG -s "$TREATMENT_SAMPLES" -c $THRESHOLD "results/merged/analysis/${GROUP}/specific.vcf.gz" \
      -Oz -o "results/merged/analysis/${GROUP}/highconf.vcf.gz"
    tabix -p vcf "results/merged/analysis/${GROUP}/highconf.vcf.gz"
    
    HC_COUNT=$(bcftools view -H "results/merged/analysis/${GROUP}/highconf.vcf.gz" | wc -l)
    echo "    Found $HC_COUNT high-confidence variants"
    
    # Generate summary statistics for this group
    mkdir -p "results/merged/summary/${GROUP}"
    bcftools view -H "results/merged/analysis/${GROUP}/specific.vcf.gz" | wc -l > \
      "results/merged/summary/${GROUP}/specific_count.txt"
    bcftools view -H "results/merged/analysis/${GROUP}/highconf.vcf.gz" | wc -l > \
      "results/merged/summary/${GROUP}/highconf_count.txt"
}

# New function to compare two treatment groups directly
compare_treatment_groups() {
    local COMPARISON_NAME="$1"
    local GROUP1_SAMPLES="$2"
    local GROUP2_SAMPLES="$3"
    
    echo "  Comparing treatment groups: $COMPARISON_NAME"
    echo "    Group 1: $GROUP1_SAMPLES"
    echo "    Group 2: $GROUP2_SAMPLES"
    
    # Create output directory
    OUT_DIR="results/merged/group_comparisons/${COMPARISON_NAME}"
    mkdir -p "$OUT_DIR"
    
    # Extract variants from each group
    bcftools view $FORCE_FLAG -s $GROUP1_SAMPLES results/merged/fixed/all_samples.vcf.gz -c 2 -Oz -o "${OUT_DIR}/group1.vcf.gz"
    bcftools view $FORCE_FLAG -s $GROUP2_SAMPLES results/merged/fixed/all_samples.vcf.gz -c 2 -Oz -o "${OUT_DIR}/group2.vcf.gz"
    
    # Index files
    bcftools index "${OUT_DIR}/group1.vcf.gz"
    bcftools index "${OUT_DIR}/group2.vcf.gz"
    
    # Find variants unique to each group (using isec for more reliable comparison)
    mkdir -p "${OUT_DIR}/isec"
    bcftools isec -p "${OUT_DIR}/isec" "${OUT_DIR}/group1.vcf.gz" "${OUT_DIR}/group2.vcf.gz"
    
    # Count and report differences
    G1_SPECIFIC=$(grep -v "^#" "${OUT_DIR}/isec/0000.vcf" | wc -l)
    G2_SPECIFIC=$(grep -v "^#" "${OUT_DIR}/isec/0001.vcf" | wc -l)
    SHARED=$(grep -v "^#" "${OUT_DIR}/isec/0002.vcf" | wc -l)
    
    echo "    Group 1 specific: $G1_SPECIFIC variants"
    echo "    Group 2 specific: $G2_SPECIFIC variants"
    echo "    Shared: $SHARED variants"
    
    # Save results to summary
    mkdir -p "results/merged/summary/${COMPARISON_NAME}"
    echo "Group 1 specific: $G1_SPECIFIC variants" > "results/merged/summary/${COMPARISON_NAME}/comparison.txt"
    echo "Group 2 specific: $G2_SPECIFIC variants" >> "results/merged/summary/${COMPARISON_NAME}/comparison.txt"
    echo "Shared: $SHARED variants" >> "results/merged/summary/${COMPARISON_NAME}/comparison.txt"
}

# Primary comparisons against WT-CTRL (biologically correct)
compare_treatment_control "WT-37" "WT-CTRL" "WT-37-55-1,WT-37-55-2,WT-37-55-3"
compare_treatment_control "WTA" "WT-CTRL" "WTA-55-1,WTA-55-2,WTA-55-3"
compare_treatment_control "STC" "WT-CTRL" "STC-55-1,STC-55-2,STC-55-3"
compare_treatment_control "CAS" "WT-CTRL" "CAS-55-1,CAS-55-2,CAS-55-3"

# Original control comparisons
compare_treatment_control "STC-vs-STCCTRL" "STC-CTRL" "STC-55-1,STC-55-2,STC-55-3"
compare_treatment_control "CAS-vs-CASCTRL" "CAS-CTRL" "CAS-55-1,CAS-55-2,CAS-55-3"

# Biologically relevant direct comparisons
# Gene effects under adaptation conditions
compare_treatment_groups "STC-vs-WTA" "STC-55-1,STC-55-2,STC-55-3" "WTA-55-1,WTA-55-2,WTA-55-3"
compare_treatment_groups "CAS-vs-WT37" "CAS-55-1,CAS-55-2,CAS-55-3" "WT-37-55-1,WT-37-55-2,WT-37-55-3"

# Different adaptation mechanisms
compare_treatment_groups "WT37-vs-WTA" "WT-37-55-1,WT-37-55-2,WT-37-55-3" "WTA-55-1,WTA-55-2,WTA-55-3"