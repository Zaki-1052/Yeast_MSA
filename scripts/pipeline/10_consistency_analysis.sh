#!/bin/bash

echo "=== Analyzing consistency within treatment replicates ==="

# Load the force flag if needed
if [ -f "results/merged/fixed/force_flag.txt" ]; then
    source results/merged/fixed/force_flag.txt
else
    FORCE_FLAG=""
fi

# Function to analyze replicate consistency
analyze_group_consistency() {
    GROUP=$1
    SAMPLES=$2
    
    echo "  Analyzing consistency for $GROUP group"
    
    # Create output directory
    OUT_DIR="results/merged/consistency/${GROUP}"
    mkdir -p "$OUT_DIR"
    
    # Extract samples from merged VCF
    bcftools view $FORCE_FLAG -s $SAMPLES results/merged/fixed/all_samples.vcf.gz -Oz -o "${OUT_DIR}/samples.vcf.gz"
    bcftools index "${OUT_DIR}/samples.vcf.gz"
    
    # Count variants present in all replicates
    NUM_SAMPLES=$(echo "$SAMPLES" | tr -cd ',' | wc -c)
    NUM_SAMPLES=$((NUM_SAMPLES + 1))
    
    bcftools view -s $SAMPLES -c $NUM_SAMPLES "${OUT_DIR}/samples.vcf.gz" -Oz -o "${OUT_DIR}/all_replicates.vcf.gz"
    ALL=$(bcftools view -H "${OUT_DIR}/all_replicates.vcf.gz" | wc -l)
    
    # Count variants present in at least 2 replicates
    AT_LEAST_2=$(bcftools view -s $SAMPLES -c 2 "${OUT_DIR}/samples.vcf.gz" | bcftools view -H | wc -l)
    
    # Count variants present in only 1 replicate
    ONLY_1=$(bcftools view -s $SAMPLES -c 1 -C 2 "${OUT_DIR}/samples.vcf.gz" | bcftools view -H | wc -l)
    
    echo "  $GROUP consistency results:"
    echo "    Variants in all $NUM_SAMPLES replicates: $ALL"
    echo "    Variants in at least 2 replicates: $AT_LEAST_2"
    echo "    Variants in only 1 replicate: $ONLY_1"
    echo "    Total variants: $((ALL + ONLY_1))"
    
    # Save results to summary
    mkdir -p "results/merged/summary/${GROUP}"
    echo "$GROUP consistency:" > "results/merged/summary/${GROUP}/consistency.txt"
    echo "  Variants in all $NUM_SAMPLES replicates: $ALL" >> "results/merged/summary/${GROUP}/consistency.txt"
    echo "  Variants in exactly 2 replicates: $((AT_LEAST_2 - ALL))" >> "results/merged/summary/${GROUP}/consistency.txt"
    echo "  Variants in only 1 replicate: $ONLY_1" >> "results/merged/summary/${GROUP}/consistency.txt"
    echo "  Total variants: $((AT_LEAST_2 + ONLY_1))" >> "results/merged/summary/${GROUP}/consistency.txt"
}

# Analyze consistency for each treatment group with correct biological groupings
analyze_group_consistency "WT-37" "WT-37-55-1,WT-37-55-2,WT-37-55-3"
analyze_group_consistency "WTA" "WTA-55-1,WTA-55-2,WTA-55-3"
analyze_group_consistency "STC" "STC-55-1,STC-55-2,STC-55-3"
analyze_group_consistency "CAS" "CAS-55-1,CAS-55-2,CAS-55-3"