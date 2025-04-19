#!/bin/bash

echo "=== Direct Comparison Workflow ==="

# Create organized directory structure
mkdir -p results/direct/{filtered,comparison,highconf,summary}

# 1. Filter, compress, and index each normalized VCF
echo "Step 1: Creating filtered VCF files..."

for NORM_VCF in results/vcf/individual/*.norm.vcf; do
    SAMPLE=$(basename "$NORM_VCF" .norm.vcf)
    echo "  Processing $SAMPLE..."
    
    # Filter by quality and read depth (FORMAT/DP)
    bcftools filter -i 'QUAL>=20 && FORMAT/DP>=10' \
      -o "results/direct/filtered/${SAMPLE}.filtered.vcf" \
      "$NORM_VCF"
    
    # Compress and index
    bgzip -f "results/direct/filtered/${SAMPLE}.filtered.vcf"
    tabix -p vcf "results/direct/filtered/${SAMPLE}.filtered.vcf.gz"
    
    # Generate basic statistics
    bcftools stats "results/direct/filtered/${SAMPLE}.filtered.vcf.gz" > \
      "results/direct/filtered/${SAMPLE}.stats.txt"
done

# Create a sample summary
echo -e "Sample\tTotal_Variants\tSNPs\tIndels" > results/direct/summary/variant_counts.txt
for STATS in results/direct/filtered/*.stats.txt; do
    SAMPLE=$(basename "$STATS" .stats.txt)
    TOTAL=$(bcftools view -H "results/direct/filtered/${SAMPLE}.filtered.vcf.gz" | wc -l)
    SNPS=$(bcftools view -H -v snps "results/direct/filtered/${SAMPLE}.filtered.vcf.gz" | wc -l)
    INDELS=$(bcftools view -H -v indels "results/direct/filtered/${SAMPLE}.filtered.vcf.gz" | wc -l)
    
    echo -e "$SAMPLE\t$TOTAL\t$SNPS\t$INDELS" >> results/direct/summary/variant_counts.txt
done

# 2. Direct treatment vs control comparison
echo -e "\nStep 2: Performing direct treatment vs control comparisons..."

# Function for direct comparison
compare_treatment_control() {
    GROUP=$1
    CONTROL=$2
    TREATMENT_PREFIX=$3
    
    echo "  Comparing $GROUP treatments vs control..."
    mkdir -p "results/direct/comparison/$GROUP"
    
    # Find all treatment samples
    TREATMENT_VCFS=$(ls results/direct/filtered/${TREATMENT_PREFIX}*.filtered.vcf.gz)
    
    # Process each treatment sample
    for TREAT_VCF in $TREATMENT_VCFS; do
        SAMPLE=$(basename "$TREAT_VCF" .filtered.vcf.gz)
        echo "    Processing $SAMPLE..."
        
        # Create directory for this comparison
        mkdir -p "results/direct/comparison/$GROUP/${SAMPLE}_vs_control"
        
        # Find variants present in treatment but not in control
        # Using -c all to match exact alleles, not just positions
        bcftools isec -p "results/direct/comparison/$GROUP/${SAMPLE}_vs_control" \
                    -n =1 -c all \
                    "$TREAT_VCF" "results/direct/filtered/${CONTROL}.filtered.vcf.gz"
        
        # Check if any treatment-specific variants were found
        if [ -f "results/direct/comparison/$GROUP/${SAMPLE}_vs_control/0000.vcf" ]; then
            # Copy to a clean location with clear naming
            cp "results/direct/comparison/$GROUP/${SAMPLE}_vs_control/0000.vcf" \
               "results/direct/comparison/$GROUP/${SAMPLE}_specific.vcf"
               
            COUNT=$(grep -v "^#" "results/direct/comparison/$GROUP/${SAMPLE}_specific.vcf" | wc -l)
            echo "    $SAMPLE has $COUNT specific variants"
        else
            echo "    No treatment-specific variants found for $SAMPLE"
        fi
    done
}

# Run comparisons for each treatment group
compare_treatment_control "WT" "WT-CTRL" "WT-37-55"
compare_treatment_control "STC" "STC-CTRL" "STC-55"
compare_treatment_control "CAS" "CAS-CTRL" "CAS-55"
compare_treatment_control "WTA" "WT-CTRL" "WTA-55"  # WTA uses WT control

# 3. Find high-confidence variants (present in multiple replicates)
echo -e "\nStep 3: Identifying high-confidence variants across replicates..."

find_high_confidence() {
    GROUP=$1
    
    echo "  Finding high-confidence variants for $GROUP..."
    mkdir -p "results/direct/highconf/$GROUP"
    
    # Check if specific variant files exist
    if ls "results/direct/comparison/$GROUP/"*_specific.vcf >/dev/null 2>&1; then
        # Extract variants with allele specificity
        for SPEC_VCF in "results/direct/comparison/$GROUP/"*_specific.vcf; do
            SAMPLE=$(basename "$SPEC_VCF" _specific.vcf)
            
            grep -v "^#" "$SPEC_VCF" | \
            awk '{print $1"_"$2"_"$4"_"$5}' > \
            "results/direct/highconf/$GROUP/${SAMPLE}_variants.txt"
        done
        
        # Find shared variants
        cat "results/direct/highconf/$GROUP/"*_variants.txt | \
        sort | uniq -c | sort -nr > \
        "results/direct/highconf/$GROUP/variant_counts.txt"
        
        # Extract high-confidence variants (in 2+ samples)
        awk '$1 >= 2 {print}' "results/direct/highconf/$GROUP/variant_counts.txt" > \
        "results/direct/highconf/$GROUP/highconf_variants.txt"
        
        # Create detailed table
        echo -e "Count\tChrom\tPos\tRef\tAlt" > "results/direct/highconf/$GROUP/highconf_details.tsv"
        awk '$1 >= 2 {
            split($2, a, "_");
            print $1"\t"a[1]"\t"a[2]"\t"a[3]"\t"a[4]
        }' "results/direct/highconf/$GROUP/variant_counts.txt" >> \
        "results/direct/highconf/$GROUP/highconf_details.tsv"
        
        # Count categories
        ALL_COUNT=$(cat "results/direct/highconf/$GROUP/variant_counts.txt" | wc -l)
        HC_COUNT=$(cat "results/direct/highconf/$GROUP/highconf_variants.txt" | wc -l)
        THREE_COUNT=$(awk '$1 == 3' "results/direct/highconf/$GROUP/variant_counts.txt" | wc -l)
        TWO_COUNT=$(awk '$1 == 2' "results/direct/highconf/$GROUP/variant_counts.txt" | wc -l)
        ONE_COUNT=$(awk '$1 == 1' "results/direct/highconf/$GROUP/variant_counts.txt" | wc -l)
        
        echo "  $GROUP variant statistics:"
        echo "    Total variants: $ALL_COUNT"
        echo "    High-confidence variants (2+ samples): $HC_COUNT"
        echo "    In all 3 samples: $THREE_COUNT"
        echo "    In exactly 2 samples: $TWO_COUNT"
        echo "    In only 1 sample: $ONE_COUNT"
    else
        echo "  No specific variants found for $GROUP"
    fi
}

# Process each treatment group
find_high_confidence "WT"
find_high_confidence "STC"
find_high_confidence "CAS"
find_high_confidence "WTA"

# 4. Create cross-treatment comparison
echo -e "\nStep 4: Creating cross-treatment comparison..."

# Collect all high-confidence variants
mkdir -p "results/direct/summary"
echo -e "Chrom\tPos\tRef\tAlt\tWT\tSTC\tCAS\tWTA" > "results/direct/summary/treatment_comparison.tsv"

# Create master list of all high-confidence variants
for GROUP in WT STC CAS WTA; do
    if [ -f "results/direct/highconf/$GROUP/highconf_variants.txt" ]; then
        awk '{print $2}' "results/direct/highconf/$GROUP/highconf_variants.txt"
    fi
done | sort | uniq > "results/direct/summary/all_highconf_variants.txt"

# Fill in the comparison matrix
while read VARIANT; do
    IFS='_' read -r CHROM POS REF ALT <<< "$VARIANT"
    LINE="$CHROM\t$POS\t$REF\t$ALT"
    
    # Check presence in each treatment group
    for GROUP in WT STC CAS WTA; do
        if [ -f "results/direct/highconf/$GROUP/variant_counts.txt" ] && \
           grep -q "$VARIANT" "results/direct/highconf/$GROUP/variant_counts.txt"; then
            COUNT=$(grep "$VARIANT" "results/direct/highconf/$GROUP/variant_counts.txt" | awk '{print $1}')
            LINE="$LINE\t$COUNT"
        else
            LINE="$LINE\t0"
        fi
    done
    
    echo -e "$LINE" >> "results/direct/summary/treatment_comparison.tsv"
done < "results/direct/summary/all_highconf_variants.txt"

# 5. Identify unique variants per treatment
echo -e "\nStep 5: Identifying treatment-unique variants..."

echo -e "Treatment\tUnique_Variants" > "results/direct/summary/unique_variants.txt"

# Count unique variants per treatment
awk 'NR>1 && ($5>=2 && $6==0 && $7==0 && $8==0) {count++} END {print "WT\t"count}' \
  "results/direct/summary/treatment_comparison.tsv" >> "results/direct/summary/unique_variants.txt"
  
awk 'NR>1 && ($5==0 && $6>=2 && $7==0 && $8==0) {count++} END {print "STC\t"count}' \
  "results/direct/summary/treatment_comparison.tsv" >> "results/direct/summary/unique_variants.txt"
  
awk 'NR>1 && ($5==0 && $6==0 && $7>=2 && $8==0) {count++} END {print "CAS\t"count}' \
  "results/direct/summary/treatment_comparison.tsv" >> "results/direct/summary/unique_variants.txt"
  
awk 'NR>1 && ($5==0 && $6==0 && $7==0 && $8>=2) {count++} END {print "WTA\t"count}' \
  "results/direct/summary/treatment_comparison.tsv" >> "results/direct/summary/unique_variants.txt"

echo "Direct comparison workflow complete! Results in results/direct/ directory."