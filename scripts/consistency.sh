#!/bin/bash

# Analyze consistency of variants within replicates
mkdir -p results/vcf/consistency

# Function to analyze group consistency
analyze_group() {
    GROUP=$1
    SAMPLES=$2
    
    echo "Analyzing consistency for $GROUP group"
    
    # Extract just these samples
    bcftools view -s $SAMPLES results/vcf/merged/all_samples.vcf.gz -Oz -o results/vcf/consistency/${GROUP}_samples.vcf.gz
    
    # Count variants present in all replicates
    bcftools view -s $SAMPLES -c 3 results/vcf/consistency/${GROUP}_samples.vcf.gz -Oz -o results/vcf/consistency/${GROUP}_all_replicates.vcf.gz
    ALL=$(bcftools view -H results/vcf/consistency/${GROUP}_all_replicates.vcf.gz | wc -l)
    
    # Count variants present in at least 2 replicates
    bcftools view -s $SAMPLES -c 2 results/vcf/consistency/${GROUP}_samples.vcf.gz -Oz -o results/vcf/consistency/${GROUP}_at_least_2.vcf.gz
    AT_LEAST_2=$(bcftools view -H results/vcf/consistency/${GROUP}_at_least_2.vcf.gz | wc -l)
    
    # Count variants present in only 1 replicate
    bcftools view -s $SAMPLES -c 1 -C 2 results/vcf/consistency/${GROUP}_samples.vcf.gz -Oz -o results/vcf/consistency/${GROUP}_only_1.vcf.gz
    ONLY_1=$(bcftools view -H results/vcf/consistency/${GROUP}_only_1.vcf.gz | wc -l)
    
    echo "$GROUP consistency:"
    echo "  Variants in all 3 replicates: $ALL"
    echo "  Variants in exactly 2 replicates: $((AT_LEAST_2 - ALL))"
    echo "  Variants in only 1 replicate: $ONLY_1"
    echo "  Total variants: $((ALL + (AT_LEAST_2 - ALL) + ONLY_1))"
}

# Analyze each treatment group
analyze_group "WT" "WT-37-55-1,WT-37-55-2,WT-37-55-3"
analyze_group "STC" "STC-55-1,STC-55-2,STC-55-3"
analyze_group "CAS" "CAS-55-1,CAS-55-2,CAS-55-3"
analyze_group "WTA" "WTA-55-1,WTA-55-2,WTA-55-3"