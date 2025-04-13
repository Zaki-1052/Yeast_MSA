#!/bin/bash

# Create directory for comparison results
mkdir -p results/vcf/fixed_comparison

# First, check sample names in fixed VCF file
echo "Using sample names from fixed VCF..."
bcftools query -l results/vcf/merged/all_samples_fixed.vcf.gz > results/vcf/fixed_sample_names.txt

# Define sample groups
WT_CTRL="WT-CTRL"
WT_TREAT="WT-37-55-1,WT-37-55-2,WT-37-55-3"
STC_CTRL="STC-CTRL"
STC_TREAT="STC-55-1,STC-55-2,STC-55-3"
CAS_CTRL="CAS-CTRL"
CAS_TREAT="CAS-55-1,CAS-55-2,CAS-55-3"
WTA_TREAT="WTA-55-1,WTA-55-2,WTA-55-3"

echo "Using these sample groups:"
echo "WT Control: $WT_CTRL"
echo "WT Treatment: $WT_TREAT"
echo "STC Control: $STC_CTRL"
echo "STC Treatment: $STC_TREAT"
echo "CAS Control: $CAS_CTRL"
echo "CAS Treatment: $CAS_TREAT"
echo "WTA Treatment: $WTA_TREAT"

# Function to compare a treatment group with its control
compare_group() {
    GROUP=$1
    TREATMENT_SAMPLES=$2
    CONTROL_SAMPLE=$3
    
    echo "Comparing $GROUP treatments vs control"
    
    # Extract just this group with --force-samples flag
    bcftools view --force-samples -s ${CONTROL_SAMPLE},${TREATMENT_SAMPLES} \
      results/vcf/merged/all_samples_fixed.vcf.gz -Oz -o results/vcf/fixed_comparison/${GROUP}_group.vcf.gz
    bcftools index results/vcf/fixed_comparison/${GROUP}_group.vcf.gz
    
    # Find variants present in treatment but not in control
    bcftools view --force-samples -s ${TREATMENT_SAMPLES} -c 1 results/vcf/fixed_comparison/${GROUP}_group.vcf.gz | \
    bcftools view --force-samples -s ${CONTROL_SAMPLE} -c 0 -Oz -o results/vcf/fixed_comparison/${GROUP}_specific.vcf.gz
    
    # Count treatment-specific variants
    COUNT=$(bcftools view -H results/vcf/fixed_comparison/${GROUP}_specific.vcf.gz | wc -l)
    echo "$GROUP treatment-specific variants: $COUNT"
    
    # Calculate threshold for high-confidence variants
    NUM_SAMPLES=$(echo "$TREATMENT_SAMPLES" | tr -cd ',' | wc -c)
    NUM_SAMPLES=$((NUM_SAMPLES + 1))
    THRESHOLD=$((NUM_SAMPLES - 1))  # Allow for one missing sample
    
    # Find high-confidence variants
    bcftools view --force-samples -s ${TREATMENT_SAMPLES} -c $THRESHOLD \
      results/vcf/fixed_comparison/${GROUP}_specific.vcf.gz -Oz -o results/vcf/fixed_comparison/${GROUP}_specific_highconf.vcf.gz
    
    HIGH_COUNT=$(bcftools view -H results/vcf/fixed_comparison/${GROUP}_specific_highconf.vcf.gz | wc -l)
    echo "$GROUP high-confidence treatment-specific variants: $HIGH_COUNT"
    
    # Index result files
    bcftools index results/vcf/fixed_comparison/${GROUP}_specific.vcf.gz
    bcftools index results/vcf/fixed_comparison/${GROUP}_specific_highconf.vcf.gz
}

# Compare each treatment group
compare_group "WT" "$WT_TREAT" "$WT_CTRL"
compare_group "STC" "$STC_TREAT" "$STC_CTRL" 
compare_group "CAS" "$CAS_TREAT" "$CAS_CTRL"
compare_group "WTA" "$WTA_TREAT" "$WT_CTRL"

echo "Fixed VCF comparison complete!"