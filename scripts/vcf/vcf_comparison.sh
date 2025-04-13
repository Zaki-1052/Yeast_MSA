#!/bin/bash

echo "=== Method Comparison: Information Gain Analysis ==="
mkdir -p method_comparison

echo "1. Running original direct comparison method..."
mkdir -p method_comparison/original

# Original method function (simulating what we did before)
original_method() {
    GROUP=$1
    CONTROL=$2
    TREATMENT_PREFIX=$3
    
    # Find treatment-specific variants
    echo "  Processing $GROUP with original method"
    
    # Compare files using bcftools isec directly
    mkdir -p method_comparison/original/$GROUP
    
    # Run comparison for each sample
    count=0
    for TREAT_VCF in results/vcf/filtered/${TREATMENT_PREFIX}*.filtered.vcf.gz; do
        bcftools isec -p method_comparison/original/$GROUP/temp \
            -n =1 -c none \
            "$TREAT_VCF" "results/vcf/filtered/${CONTROL}.filtered.vcf.gz"
            
        # Count variants
        if [ -f "method_comparison/original/$GROUP/temp/0000.vcf" ]; then
            this_count=$(grep -v "^#" "method_comparison/original/$GROUP/temp/0000.vcf" | wc -l)
            count=$((count + this_count))
        fi
    done
    
    # Generate position-only counts
    mkdir -p method_comparison/original/$GROUP/positions
    for TREAT_VCF in results/vcf/filtered/${TREATMENT_PREFIX}*.filtered.vcf.gz; do
        SAMPLE=$(basename "$TREAT_VCF" .filtered.vcf.gz)
        bcftools isec -p method_comparison/original/$GROUP/positions/$SAMPLE \
            -n =1 -c none \
            "$TREAT_VCF" "results/vcf/filtered/${CONTROL}.filtered.vcf.gz"
            
        # Extract just positions
        if [ -f "method_comparison/original/$GROUP/positions/$SAMPLE/0000.vcf" ]; then
            grep -v "^#" "method_comparison/original/$GROUP/positions/$SAMPLE/0000.vcf" | \
            awk '{print $1"_"$2}' > method_comparison/original/$GROUP/positions/${SAMPLE}_pos.txt
        fi
    done
    
    # Original method for finding shared positions (imprecise)
    cat method_comparison/original/$GROUP/positions/*_pos.txt | sort | uniq -c | \
    sort -nr > method_comparison/original/$GROUP/shared_positions.txt
    
    HC_COUNT=$(awk '$1 >= 2' method_comparison/original/$GROUP/shared_positions.txt | wc -l)
    
    echo "  $GROUP: $count total variants across samples"
    echo "  $GROUP: $HC_COUNT high-confidence positions (old method)"
}

# Run original method
original_method "WT" "WT-CTRL" "WT-37-55"
original_method "STC" "STC-CTRL" "STC-55"
original_method "CAS" "CAS-CTRL" "CAS-55"
original_method "WTA" "WT-CTRL" "WTA-55"

echo -e "\n2. Running improved method with better variant tracking..."
mkdir -p method_comparison/improved

improved_method() {
    GROUP=$1
    CONTROL=$2
    TREATMENT_PREFIX=$3
    
    echo "  Processing $GROUP with improved method"
    mkdir -p method_comparison/improved/$GROUP
    
    # Process each sample, extracting clean variant identifiers
    for TREAT_VCF in results/vcf/filtered/${TREATMENT_PREFIX}*.filtered.vcf.gz; do
        SAMPLE=$(basename "$TREAT_VCF" .filtered.vcf.gz)
        
        # Run bcftools isec
        bcftools isec -p method_comparison/improved/$GROUP/$SAMPLE \
            -n =1 -c all \
            "$TREAT_VCF" "results/vcf/filtered/${CONTROL}.filtered.vcf.gz"
        
        # Extract full variant info (chr_pos_ref_alt)
        if [ -f "method_comparison/improved/$GROUP/$SAMPLE/0000.vcf" ]; then
            grep -v "^#" "method_comparison/improved/$GROUP/$SAMPLE/0000.vcf" | \
            awk '{print $1"_"$2"_"$4"_"$5}' > method_comparison/improved/$GROUP/${SAMPLE}_variants.txt
        fi
    done
    
    # Find shared variants (with allele specificity)
    cat method_comparison/improved/$GROUP/*_variants.txt | sort | uniq -c | \
    sort -nr > method_comparison/improved/$GROUP/shared_variants.txt
    
    # Categorize variants by how many samples they appear in
    TOTAL=$(cat method_comparison/improved/$GROUP/*_variants.txt | wc -l)
    ONE_SAMPLE=$(awk '$1 == 1' method_comparison/improved/$GROUP/shared_variants.txt | wc -l)
    TWO_SAMPLES=$(awk '$1 == 2' method_comparison/improved/$GROUP/shared_variants.txt | wc -l)
    THREE_SAMPLES=$(awk '$1 == 3' method_comparison/improved/$GROUP/shared_variants.txt | wc -l)
    HC_TOTAL=$((TWO_SAMPLES + THREE_SAMPLES))
    
    echo "  $GROUP: $TOTAL total variants across samples"
    echo "  $GROUP: $HC_TOTAL high-confidence variants (improved method)"
    echo "    - $ONE_SAMPLE variants in 1 sample only"
    echo "    - $TWO_SAMPLES variants in exactly 2 samples"
    echo "    - $THREE_SAMPLES variants in all 3 samples"
    
    # Output top shared variants for examination
    echo "  Top variants:" > method_comparison/improved/$GROUP/summary.txt
    head -5 method_comparison/improved/$GROUP/shared_variants.txt >> method_comparison/improved/$GROUP/summary.txt
}

# Run improved method
improved_method "WT" "WT-CTRL" "WT-37-55"
improved_method "STC" "STC-CTRL" "STC-55"
improved_method "CAS" "CAS-CTRL" "CAS-55"
improved_method "WTA" "WT-CTRL" "WTA-55"

# Compare results - information gain analysis
echo -e "\n3. Information Gain Analysis by Group:"
for GROUP in WT STC CAS WTA; do
    # Get counts from both methods
    OLD_HC=$(awk '$1 >= 2' method_comparison/original/$GROUP/shared_positions.txt | wc -l)
    NEW_HC=$(awk '$1 >= 2' method_comparison/improved/$GROUP/shared_variants.txt | wc -l)
    
    # Calculate information gain
    VARIANT_GAIN=$((NEW_HC - OLD_HC))
    if [ $OLD_HC -eq 0 ]; then
        PERCENT="∞%"  # Avoid division by zero
    else
        PERCENT=$(echo "scale=1; 100 * $VARIANT_GAIN / $OLD_HC" | bc)"%"
    fi
    
    # Check for shared variants that would be missed
    THREE_SAMPLES=$(awk '$1 == 3' method_comparison/improved/$GROUP/shared_variants.txt | wc -l)
    
    echo "$GROUP Group:"
    echo "  - Original method: $OLD_HC high-confidence positions"
    echo "  - Improved method: $NEW_HC high-confidence variants"
    echo "  - Information gain: $VARIANT_GAIN more variants ($PERCENT improvement)"
    echo "  - $THREE_SAMPLES variants present in ALL replicates (highest confidence)"
    
    # Show example of variants we would have missed
    if [ $THREE_SAMPLES -gt 0 ]; then
        echo "  - Example high-quality variants we would have missed:"
        awk '$1 == 3 {print "    * "$2}' method_comparison/improved/$GROUP/shared_variants.txt | head -3
    fi
    echo ""
done

echo "=== Complete Comparison Analysis ==="