#!/bin/bash

echo "=== Direct Single-Step Comparison ==="

mkdir -p results/fixed_comparison

# Function to perform direct comparison
compare_direct() {
    GROUP=$1
    CONTROL=$2
    TREATMENT_SAMPLES=$3
    
    echo "Processing $GROUP comparison..."
    
    # Create a file with sample names for treatment
    echo "$TREATMENT_SAMPLES" | tr ',' '\n' > results/fixed_comparison/${GROUP}_treatment.txt
    
    # Create a file with the control sample
    echo "$CONTROL" > results/fixed_comparison/${GROUP}_control.txt
    
    # Directly extract variants that are in treatments but not control in one step
    echo "  Extracting treatment-specific variants..."
    bcftools view -S results/fixed_comparison/${GROUP}_treatment.txt --force-samples results/merged/fixed/all_samples.vcf.gz | \
      bcftools view -S results/fixed_comparison/${GROUP}_control.txt --force-samples -c 0 -Oz -o results/fixed_comparison/${GROUP}_specific.vcf.gz
    
    # Count variants
    COUNT=$(bcftools view -H results/fixed_comparison/${GROUP}_specific.vcf.gz | wc -l)
    echo "  $GROUP treatment-specific variants: $COUNT"
    
    # Find high-confidence variants
    echo "  Identifying high-confidence variants..."
    bcftools view -S results/fixed_comparison/${GROUP}_treatment.txt --force-samples results/fixed_comparison/${GROUP}_specific.vcf.gz -c 2 -Oz -o results/fixed_comparison/${GROUP}_highconf.vcf.gz
    
    # Count high-confidence variants
    HC_COUNT=$(bcftools view -H results/fixed_comparison/${GROUP}_highconf.vcf.gz | wc -l)
    echo "  $GROUP high-confidence variants: $HC_COUNT"
}

# Run comparisons
compare_direct "WT" "WT-CTRL" "WT-37-55-1,WT-37-55-2,WT-37-55-3"
compare_direct "STC" "STC-CTRL" "STC-55-1,STC-55-2,STC-55-3"
compare_direct "CAS" "CAS-CTRL" "CAS-55-1,CAS-55-2,CAS-55-3"
compare_direct "WTA" "WT-CTRL" "WTA-55-1,WTA-55-2,WTA-55-3"

echo "Direct comparison complete!"