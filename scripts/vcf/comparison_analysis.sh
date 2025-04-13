#!/bin/bash

# Create directory for comparison results
mkdir -p results/vcf/comparison

# First, let's check what sample names we actually have
echo "Checking actual sample names in VCF file..."
bcftools query -l results/vcf/merged/all_samples.vcf.gz > results/vcf/sample_names.txt
cat results/vcf/sample_names.txt

# Define a function to identify correct sample names based on patterns
find_samples() {
    PATTERN=$1
    grep "$PATTERN" results/vcf/sample_names.txt | tr '\n' ',' | sed 's/,$//'
}

# Get sample names for each group
WT_CTRL=$(find_samples "WT-CTRL")
WT_TREAT=$(find_samples "WT-37-55")
STC_CTRL=$(find_samples "STC-CTRL")
STC_TREAT=$(find_samples "STC-55")
CAS_CTRL=$(find_samples "CAS-CTRL")
CAS_TREAT=$(find_samples "CAS-55")
WTA_TREAT=$(find_samples "WTA-55")

echo "Using these sample names:"
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
    
    # Skip if any sample list is empty
    if [ -z "$TREATMENT_SAMPLES" ] || [ -z "$CONTROL_SAMPLE" ]; then
        echo "Skipping $GROUP comparison - missing sample names"
        return
    fi
    
    echo "Comparing $GROUP treatments vs control"
    
    # Extract just this group (control + treatments)
    bcftools view -s ${CONTROL_SAMPLE},${TREATMENT_SAMPLES} results/vcf/merged/all_samples.vcf.gz -Oz -o results/vcf/comparison/${GROUP}_group.vcf.gz
    bcftools index results/vcf/comparison/${GROUP}_group.vcf.gz
    
    # Find variants present in at least one treatment sample but not in control
    bcftools view -s ${TREATMENT_SAMPLES} -c 1 results/vcf/comparison/${GROUP}_group.vcf.gz | \
    bcftools view -s ${CONTROL_SAMPLE} -c 0 -Oz -o results/vcf/comparison/${GROUP}_specific.vcf.gz
    
    # Count treatment-specific variants
    COUNT=$(bcftools view -H results/vcf/comparison/${GROUP}_specific.vcf.gz | wc -l)
    echo "$GROUP treatment-specific variants: $COUNT"
    
    # Extract high-confidence treatment-specific variants (present in most treatment replicates)
    # Count how many treatment samples we have by counting commas
    NUM_SAMPLES=$(echo "$TREATMENT_SAMPLES" | tr -cd ',' | wc -c)
    NUM_SAMPLES=$((NUM_SAMPLES + 1))
    THRESHOLD=$((NUM_SAMPLES - 1))  # Allow for one missing sample
    
    bcftools view -s ${TREATMENT_SAMPLES} -c $THRESHOLD results/vcf/comparison/${GROUP}_specific.vcf.gz -Oz -o results/vcf/comparison/${GROUP}_specific_highconf.vcf.gz
    HIGH_COUNT=$(bcftools view -H results/vcf/comparison/${GROUP}_specific_highconf.vcf.gz | wc -l)
    echo "$GROUP high-confidence treatment-specific variants: $HIGH_COUNT"
}

# Compare each treatment group
[ ! -z "$WT_CTRL" ] && [ ! -z "$WT_TREAT" ] && compare_group "WT" "$WT_TREAT" "$WT_CTRL"
[ ! -z "$STC_CTRL" ] && [ ! -z "$STC_TREAT" ] && compare_group "STC" "$STC_TREAT" "$STC_CTRL" 
[ ! -z "$CAS_CTRL" ] && [ ! -z "$CAS_TREAT" ] && compare_group "CAS" "$CAS_TREAT" "$CAS_CTRL"

# For WTA, compare against WT control (since no WTA control exists)
[ ! -z "$WT_CTRL" ] && [ ! -z "$WTA_TREAT" ] && compare_group "WTA" "$WTA_TREAT" "$WT_CTRL"

echo "Comparison analysis complete."