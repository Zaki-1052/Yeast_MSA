#!/bin/bash

echo "=== Simplified VCF Fix and Direct Comparison ==="
mkdir -p final_analysis

# 1. Since merging is problematic, let's use the direct comparison approach that works
echo "Using direct comparison approach (which we've confirmed works)"

# Function to compare samples
compare_samples() {
    GROUP=$1
    CONTROL=$2
    TREATMENT_PREFIX=$3
    
    echo "Analyzing $GROUP group..."
    mkdir -p final_analysis/$GROUP
    
    # Find treatment-specific variants for each sample
    for TREAT_VCF in results/vcf/filtered/${TREATMENT_PREFIX}*.filtered.vcf.gz; do
        SAMPLE=$(basename "$TREAT_VCF" .filtered.vcf.gz)
        echo "  Processing $SAMPLE"
        
        # Create temp directory for bcftools isec
        mkdir -p final_analysis/$GROUP/${SAMPLE}_vs_control
        
        # Use appropriate comparison mode - exact allele matches
        bcftools isec -p final_analysis/$GROUP/${SAMPLE}_vs_control \
                     -n =1 -c all \
                     "$TREAT_VCF" "results/vcf/filtered/${CONTROL}.filtered.vcf.gz"
        
        # Move the result to a cleaner location
        if [ -f "final_analysis/$GROUP/${SAMPLE}_vs_control/0000.vcf" ]; then
            cp "final_analysis/$GROUP/${SAMPLE}_vs_control/0000.vcf" \
               "final_analysis/$GROUP/${SAMPLE}_specific.vcf"
               
            COUNT=$(grep -v "^#" "final_analysis/$GROUP/${SAMPLE}_specific.vcf" | wc -l)
            echo "    $SAMPLE has $COUNT specific variants"
        else
            echo "    No treatment-specific variants found for $SAMPLE"
        fi
    done
    
    # Find high-confidence variants (present in at least 2 samples)
    if ls final_analysis/$GROUP/*_specific.vcf > /dev/null 2>&1; then
        echo "  Finding high-confidence variants for $GROUP"
        
        # Extract variant positions and count occurrences
        grep -v "^#" final_analysis/$GROUP/*_specific.vcf | \
        awk '{print $1":"$2":"$4":"$5}' | sort | uniq -c | sort -nr > \
        final_analysis/$GROUP/variant_counts.txt
        
        # Count variants present in at least 2 samples
        HIGH_CONF=$(awk '$1 >= 2' final_analysis/$GROUP/variant_counts.txt | wc -l)
        echo "  $GROUP high-confidence treatment-specific variants: $HIGH_CONF"
        
        # Create a detailed high-confidence variant list
        echo -e "Count\tCHROM\tPOS\tREF\tALT" > final_analysis/$GROUP/highconf_details.txt
        awk '$1 >= 2 {split($2,a,":"); print $1"\t"a[1]"\t"a[2]"\t"a[3]"\t"a[4]}' \
            final_analysis/$GROUP/variant_counts.txt >> \
            final_analysis/$GROUP/highconf_details.txt
    else
        echo "  No variants to analyze for $GROUP"
    fi
}

# Run analysis for each group
compare_samples "WT" "WT-CTRL" "WT-37-55"
compare_samples "STC" "STC-CTRL" "STC-55"
compare_samples "CAS" "CAS-CTRL" "CAS-55"
compare_samples "WTA" "WT-CTRL" "WTA-55"

# Create a cross-treatment comparison
echo "Creating cross-treatment comparison matrix..."

# Collect all high-confidence variant positions
cat final_analysis/*/highconf_details.txt | grep -v "^Count" | \
awk '{print $2":"$3":"$4":"$5}' | sort | uniq > \
final_analysis/all_highconf_variants.txt

# Create comparison matrix header
echo -e "CHROM\tPOS\tREF\tALT\tWT\tSTC\tCAS\tWTA" > \
final_analysis/treatment_comparison.txt

# Fill in the matrix
cat final_analysis/all_highconf_variants.txt | while read VARIANT; do
    IFS=':' read -r CHROM POS REF ALT <<< "$VARIANT"
    
    # Start building the line
    LINE="$CHROM\t$POS\t$REF\t$ALT"
    
    # Check presence in each treatment group
    for GROUP in WT STC CAS WTA; do
        if grep -q "$VARIANT" final_analysis/$GROUP/variant_counts.txt; then
            COUNT=$(grep "$VARIANT" final_analysis/$GROUP/variant_counts.txt | awk '{print $1}')
            LINE="$LINE\t$COUNT"
        else
            LINE="$LINE\t0"
        fi
    done
    
    echo -e "$LINE" >> final_analysis/treatment_comparison.txt
done

# Generate summary
echo "Summary of treatment-specific variants:"
for GROUP in WT STC CAS WTA; do
    if [ -f "final_analysis/$GROUP/highconf_details.txt" ]; then
        COUNT=$(grep -v "^Count" final_analysis/$GROUP/highconf_details.txt | wc -l)
        echo "$GROUP high-confidence variants: $COUNT"
    else
        echo "$GROUP high-confidence variants: 0"
    fi
done

echo "Analysis complete! Results in final_analysis/"

# Add this code block to the end of the script to debug the high-confidence issue:

echo "Debugging high-confidence variants issue:"
SAMPLE_DIR="final_analysis/WT"

# Examine the variant counts file
echo "Checking variant counts format:"
head -5 "$SAMPLE_DIR/variant_counts.txt"

# Check if awk is processing the file correctly 
echo "Testing awk filter:"
awk '$1 >= 2 {print $0}' "$SAMPLE_DIR/variant_counts.txt" | head -5

# Try a simpler approach to count high-confidence variants
echo "Using grep to find patterns:"
HIGH_CONF_COUNT=$(grep -E "^[ ]*[2-9]|^[ ]*[0-9]{2,}" "$SAMPLE_DIR/variant_counts.txt" | wc -l)
echo "High-confidence variants (grep method): $HIGH_CONF_COUNT"

# Add to the same script
echo -e "\n=== Debugging High-Confidence Variants ==="

# 1. Check the raw output of grep for variants
echo "Checking raw variant files for identical variants:"
mkdir -p debug_variants

# Extract clean variant strings from each file
for SAMPLE in WT-37-55-1 WT-37-55-2 WT-37-55-3; do
  grep -v "^#" final_analysis/WT/${SAMPLE}_specific.vcf | \
  awk '{print $1":"$2":"$4":"$5}' > debug_variants/${SAMPLE}_variants.txt
done

# Look for any variants that might be in multiple files
cat debug_variants/*_variants.txt | sort | uniq -c | sort -rn | head -10 > debug_variants/all_counts.txt
echo "Top variant counts:"
cat debug_variants/all_counts.txt

# 2. Try a more direct comparison of actual variant records
echo -e "\nTrying more direct variant comparison:"

# Extract chromosome and position only 
for SAMPLE in WT-37-55-1 WT-37-55-2 WT-37-55-3; do
  grep -v "^#" final_analysis/WT/${SAMPLE}_specific.vcf | \
  awk '{print $1"\t"$2}' > debug_variants/${SAMPLE}_positions.txt
done

# Count shared positions
cat debug_variants/*_positions.txt | sort | uniq -c | sort -rn | head -10 > debug_variants/position_counts.txt
echo "Top position counts (ignoring exact alleles):"
cat debug_variants/position_counts.txt

echo "=== Debugging Complete ==="