#!/bin/bash

echo "=== Merged VCF Approach with Contig Fix ==="

# Create organized directory structure
mkdir -p results/merged/{filtered,fixed,analysis,summary}

# 1. Filter, compress, and index each normalized VCF (if not already done)
echo "Step 1: Creating filtered VCF files..."

for NORM_VCF in results/vcf/individual/*.norm.vcf; do
    SAMPLE=$(basename "$NORM_VCF" .norm.vcf)
    echo "  Processing $SAMPLE..."
    
    # Skip if already filtered
    if [ -f "results/merged/filtered/${SAMPLE}.filtered.vcf.gz" ]; then
        echo "    Already filtered, skipping..."
        continue
    fi
    
    # Filter by quality and read depth
    bcftools filter -i 'QUAL>=20 && FORMAT/DP>=10' \
      -o "results/merged/filtered/${SAMPLE}.filtered.vcf" \
      "$NORM_VCF"
    
    # Compress and index
    bgzip -f "results/merged/filtered/${SAMPLE}.filtered.vcf"
    tabix -p vcf "results/merged/filtered/${SAMPLE}.filtered.vcf.gz"
done

# 2. Extract contig information from reference genome
echo -e "\nStep 2: Extracting contig information from reference..."

# Generate FAI index if it doesn't exist
if [ ! -f "reference/yeast_w303.fasta.fai" ]; then
    echo "  Creating reference genome index..."
    samtools faidx reference/yeast_w303.fasta
fi

# Create a proper contig header from the reference index
echo "  Creating contig header lines..."
awk 'BEGIN {print "##fileformat=VCFv4.2"} 
     {print "##contig=<ID="$1",length="$2">"}' reference/yeast_w303.fasta.fai > \
     results/merged/fixed/contig_header.txt

# 3. Add proper contig definitions to each filtered VCF
echo -e "\nStep 3: Adding proper contig definitions to VCFs..."

for FILT_VCF in results/merged/filtered/*.filtered.vcf.gz; do
    SAMPLE=$(basename "$FILT_VCF" .filtered.vcf.gz)
    echo "  Fixing header for $SAMPLE..."
    
    # Extract the existing header without contig lines
    bcftools view -h "$FILT_VCF" | grep -v "^##contig" | grep -v "^##fileformat" > \
      results/merged/fixed/existing_header_${SAMPLE}.txt
    
    # Combine the headers
    cat results/merged/fixed/contig_header.txt \
        results/merged/fixed/existing_header_${SAMPLE}.txt > \
        results/merged/fixed/complete_header_${SAMPLE}.txt
    
    # Replace the header
    bcftools reheader -h results/merged/fixed/complete_header_${SAMPLE}.txt \
      -o results/merged/fixed/${SAMPLE}.fixed.vcf.gz \
      "$FILT_VCF"
    
    # Index the fixed VCF
    tabix -p vcf results/merged/fixed/${SAMPLE}.fixed.vcf.gz
done

# 4. Create a list of fixed VCFs to merge
echo -e "\nStep 4: Merging VCFs with proper contig definitions..."
ls results/merged/fixed/*.fixed.vcf.gz > results/merged/fixed/vcf_list.txt

# Merge the VCFs with proper headers
bcftools merge -l results/merged/fixed/vcf_list.txt \
  -Oz -o results/merged/fixed/all_samples.vcf.gz

# Index the merged file
tabix -p vcf results/merged/fixed/all_samples.vcf.gz

# 5. Verify the fix worked by testing sample selection
echo -e "\nStep 5: Verifying sample selection works with fixed VCF..."

# Check actual sample names
bcftools query -l results/merged/fixed/all_samples.vcf.gz > \
  results/merged/fixed/sample_names.txt
echo "  Available samples:"
cat results/merged/fixed/sample_names.txt

# Test sample selection with a control sample
echo "  Testing sample selection..."
bcftools view -s WT-CTRL \
  -o results/merged/fixed/test_selection.vcf \
  results/merged/fixed/all_samples.vcf.gz
  
if [ -s results/merged/fixed/test_selection.vcf ]; then
    echo "  SUCCESS: Sample selection works with fixed VCF!"
    
    # Test multiple sample selection
    bcftools view -s WT-CTRL,WT-37-55-1 \
      -o results/merged/fixed/test_multi_selection.vcf \
      results/merged/fixed/all_samples.vcf.gz
      
    if [ -s results/merged/fixed/test_multi_selection.vcf ]; then
        echo "  SUCCESS: Multi-sample selection works!"
        FORCE_FLAG=""  # No force flag needed
    else
        echo "  WARNING: Multi-sample selection failed, using --force-samples"
        FORCE_FLAG="--force-samples"
    fi
else
    echo "  WARNING: Sample selection not working, using --force-samples"
    FORCE_FLAG="--force-samples"
fi

# 6. Perform treatment vs control comparisons using the fixed merged VCF
echo -e "\nStep 6: Performing treatment vs control comparisons..."

# Function to compare treatment samples with control
# FIXED: Consistently define the function parameters
compare_treatment_control() {
    local GROUP="$1"
    local CONTROL="$2"
    local TREATMENT_SAMPLES="$3"
    
    echo "  Comparing $GROUP treatments vs control..."
    echo "  DEBUG: Checking if ${GROUP}_group.vcf.gz has variants..."
    
    # Check if the group VCF has any variants
    TOTAL_VARIANTS=$(bcftools view -H "results/merged/analysis/${GROUP}_group.vcf.gz" | wc -l)
    echo "  DEBUG: ${GROUP}_group.vcf.gz has $TOTAL_VARIANTS total variants"
    
    # Check actual sample names in the file
    echo "  DEBUG: Sample names in ${GROUP}_group.vcf.gz:"
    bcftools query -l "results/merged/analysis/${GROUP}_group.vcf.gz"
    
    # Verify variables are defined before using
    if [ -z "$GROUP" ] || [ -z "$CONTROL" ] || [ -z "$TREATMENT_SAMPLES" ]; then
        echo "ERROR: Missing parameters for compare_treatment_control function"
        echo "  GROUP: $GROUP"
        echo "  CONTROL: $CONTROL"
        echo "  TREATMENT_SAMPLES: $TREATMENT_SAMPLES"
        return 1
    fi

    echo "STRING: '$CONTROL,$TREATMENT_SAMPLES'"
    
    # Extract just this treatment group with its control
    echo "    Extracting $GROUP group data..."
    bcftools view $FORCE_FLAG -s "$CONTROL,$TREATMENT_SAMPLES" results/merged/fixed/all_samples.vcf.gz -Oz -o "results/merged/analysis/${GROUP}_group.vcf.gz"
    tabix -p vcf "results/merged/analysis/${GROUP}_group.vcf.gz"
    
    # Find variants present in at least one treatment sample but not in control
echo "    Identifying treatment-specific variants..."
# Step 1: Find variants in treatment samples
#bcftools view $FORCE_FLAG -s "$TREATMENT_SAMPLES" -c 1 "results/merged/analysis/${GROUP}_group.vcf.gz" -o "results/merged/analysis/${GROUP}_treatments_temp.vcf"

# Step 2: From the original file, extract those exact variant positions
#bcftools view $FORCE_FLAG -s "$CONTROL" -c 0 -R "results/merged/analysis/${GROUP}_treatments_temp.vcf" "results/merged/analysis/${GROUP}_group.vcf.gz" -Oz -o "results/merged/analysis/${GROUP}_specific.vcf.gz"

# Get sample indices (positions in the VCF file)
# Get control index
CONTROL_IDX=$(bcftools query -l "results/merged/analysis/${GROUP}_group.vcf.gz" | grep -n "$CONTROL" | cut -d: -f1)
CONTROL_IDX=$((CONTROL_IDX - 1))  # Convert to 0-based index

# Get treatment samples indices
TREATMENTS=$(bcftools query -l "results/merged/analysis/${GROUP}_group.vcf.gz" | grep -v "$CONTROL" | tr '\n' ',')
TREATMENTS=${TREATMENTS%,}  # Remove trailing comma

# Create a dynamic filter expression for any number of treatment samples
FILTER_EXPR="GT[$CONTROL_IDX]=\".\""
for SAMPLE in $(echo $TREATMENTS | tr ',' ' '); do
    SAMPLE_IDX=$(bcftools query -l "results/merged/analysis/${GROUP}_group.vcf.gz" | grep -n "$SAMPLE" | cut -d: -f1)
    SAMPLE_IDX=$((SAMPLE_IDX - 1))
    FILTER_EXPR="$FILTER_EXPR || GT[$SAMPLE_IDX]=\"1\""
done
FILTER_EXPR="$FILTER_EXPR)"
FILTER_EXPR="GT[$CONTROL_IDX]=\".\" && ($FILTER_EXPR"

# Use the dynamic filter expression
bcftools view -i "$FILTER_EXPR" "results/merged/analysis/${GROUP}_group.vcf.gz" -Oz -o "results/merged/analysis/${GROUP}_specific.vcf.gz"

# Use direct filtering with sample indices
#bcftools view -i 'GT[0]="." && (GT[1]="1" || GT[2]="1" || GT[3]="1")' "results/merged/analysis/${GROUP}_group.vcf.gz" -Oz -o "results/merged/analysis/${GROUP}_specific.vcf.gz"    
    # Count treatment-specific variants
    COUNT=$(bcftools view -H "results/merged/analysis/${GROUP}_specific.vcf.gz" | wc -l)
    echo "    Found $COUNT treatment-specific variants"
    
    # Find high-confidence variants (in multiple replicates)
    echo "    Identifying high-confidence variants..."
    
    # Count samples to set threshold
    NUM_SAMPLES=$(echo "$TREATMENT_SAMPLES" | tr -cd ',' | wc -c)
    NUM_SAMPLES=$((NUM_SAMPLES + 1))
    THRESHOLD=$((NUM_SAMPLES - 1))  # n-1 replicates for high confidence
    
    bcftools view $FORCE_FLAG -s "$TREATMENT_SAMPLES" -c $THRESHOLD "results/merged/analysis/${GROUP}_specific.vcf.gz" -Oz -o "results/merged/analysis/${GROUP}_highconf.vcf.gz"
    tabix -p vcf "results/merged/analysis/${GROUP}_highconf.vcf.gz"
    
    HC_COUNT=$(bcftools view -H "results/merged/analysis/${GROUP}_highconf.vcf.gz" | wc -l)
    echo "    Found $HC_COUNT high-confidence variants"
    
    # Generate summary statistics for this group
    bcftools view -H "results/merged/analysis/${GROUP}_specific.vcf.gz" | wc -l > \
      "results/merged/summary/${GROUP}_specific_count.txt"
    bcftools view -H "results/merged/analysis/${GROUP}_highconf.vcf.gz" | wc -l > \
      "results/merged/summary/${GROUP}_highconf_count.txt"
}

# Run comparisons for each treatment group with explicit parameter passing
# Order is GROUP, CONTROL, TREATMENT_SAMPLES
compare_treatment_control "WT" "WT-CTRL" "WT-37-55-1,WT-37-55-2,WT-37-55-3"
compare_treatment_control "STC" "STC-CTRL" "STC-55-1,STC-55-2,STC-55-3"
compare_treatment_control "CAS" "CAS-CTRL" "CAS-55-1,CAS-55-2,CAS-55-3"
compare_treatment_control "WTA" "WT-CTRL" "WTA-55-1,WTA-55-2,WTA-55-3"

# 7. Create cross-treatment comparison
echo -e "\nStep 7: Creating cross-treatment comparison..."

# Extract variant details for all high-confidence variants
for GROUP in WT STC CAS WTA; do
    if [ -f "results/merged/analysis/${GROUP}_highconf.vcf.gz" ]; then
        echo "  Extracting variant details for $GROUP..."
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
          "results/merged/analysis/${GROUP}_highconf.vcf.gz" > \
          "results/merged/summary/${GROUP}_highconf_variants.txt"
    else
        echo "  WARNING: No high-confidence variants file for $GROUP"
        touch "results/merged/summary/${GROUP}_highconf_variants.txt"
    fi
done

# Create combined variant list
echo -e "Chrom\tPos\tRef\tAlt\tWT\tSTC\tCAS\tWTA" > \
  "results/merged/summary/treatment_comparison.tsv"

# Create a master list of all high-confidence variant positions
cat results/merged/summary/*_highconf_variants.txt | \
  sort -k1,1 -k2,2n | uniq > \
  "results/merged/summary/all_highconf_positions.txt"

# Process each variant position
while read CHROM POS REF ALT; do
    if [ -z "$CHROM" ] || [ -z "$POS" ] || [ -z "$REF" ] || [ -z "$ALT" ]; then
        continue  # Skip empty lines
    fi
    
    LINE="$CHROM\t$POS\t$REF\t$ALT"
    
    # Check each treatment group
    for GROUP in WT STC CAS WTA; do
        if [ -f "results/merged/summary/${GROUP}_highconf_variants.txt" ] && \
           grep -q "^$CHROM[[:space:]]$POS[[:space:]]$REF[[:space:]]$ALT$" \
           "results/merged/summary/${GROUP}_highconf_variants.txt"; then
            LINE="$LINE\t1"
        else
            LINE="$LINE\t0"
        fi
    done
    
    echo -e "$LINE" >> "results/merged/summary/treatment_comparison.tsv"
done < "results/merged/summary/all_highconf_positions.txt"

# 8. Identify unique variants per treatment group
echo -e "\nStep 8: Identifying treatment-unique variants..."

echo -e "Treatment\tUnique_Variants" > "results/merged/summary/unique_variants.txt"

# Count unique variants for each treatment
awk 'NR>1 && ($5==1 && $6==0 && $7==0 && $8==0) {count++} \
     END {print "WT\t"count}' \
     "results/merged/summary/treatment_comparison.tsv" >> \
     "results/merged/summary/unique_variants.txt"
     
awk 'NR>1 && ($5==0 && $6==1 && $7==0 && $8==0) {count++} \
     END {print "STC\t"count}' \
     "results/merged/summary/treatment_comparison.tsv" >> \
     "results/merged/summary/unique_variants.txt"
     
awk 'NR>1 && ($5==0 && $6==0 && $7==1 && $8==0) {count++} \
     END {print "CAS\t"count}' \
     "results/merged/summary/treatment_comparison.tsv" >> \
     "results/merged/summary/unique_variants.txt"
     
awk 'NR>1 && ($5==0 && $6==0 && $7==0 && $8==1) {count++} \
     END {print "WTA\t"count}' \
     "results/merged/summary/treatment_comparison.tsv" >> \
     "results/merged/summary/unique_variants.txt"

# Function for direct comparison of VCFs (more reliable approach)
compare_treatment_to_control_direct() {
    GROUP=$1
    CONTROL_VCF=$2
    TREATMENT_VCF_PREFIX=$3
    
    echo "  Direct comparison for $GROUP treatments vs control"
    
    # Create output directory
    OUT_DIR="results/merged/direct_comparison/${GROUP}"
    mkdir -p "$OUT_DIR"
    
    # Find treatment VCF files - CORRECTED PATH
    TREATMENT_VCFS=$(ls results/merged/filtered/${TREATMENT_VCF_PREFIX}*.filtered.vcf.gz 2>/dev/null)
    
    # Check if we found any files
    if [ -z "$TREATMENT_VCFS" ]; then
        echo "    No treatment VCF files found for pattern: results/merged/filtered/${TREATMENT_VCF_PREFIX}*.filtered.vcf.gz"
        echo "    Trying alternative path: results/merged/fixed/${TREATMENT_VCF_PREFIX}*.fixed.vcf.gz"
        
        # Try alternative path with fixed files
        TREATMENT_VCFS=$(ls results/merged/fixed/${TREATMENT_VCF_PREFIX}*.fixed.vcf.gz 2>/dev/null)
        
        if [ -z "$TREATMENT_VCFS" ]; then
            echo "    Still no files found. Skipping direct comparison for $GROUP"
            echo "0" > "results/merged/summary/${GROUP}_direct_highconf_count.txt"
            return
        fi
    fi
    
    # Also correct the control VCF path
    if [ ! -f "$CONTROL_VCF" ]; then
        # Try to find the control in the filtered directory
        ALT_CONTROL=$(basename "$CONTROL_VCF")
        if [ -f "results/merged/filtered/$ALT_CONTROL" ]; then
            CONTROL_VCF="results/merged/filtered/$ALT_CONTROL"
        elif [ -f "results/merged/fixed/${ALT_CONTROL%.filtered.vcf.gz}.fixed.vcf.gz" ]; then
            CONTROL_VCF="results/merged/fixed/${ALT_CONTROL%.filtered.vcf.gz}.fixed.vcf.gz"
        else
            echo "    Control VCF not found: $CONTROL_VCF"
            echo "    Skipping direct comparison for $GROUP"
            echo "0" > "results/merged/summary/${GROUP}_direct_highconf_count.txt"
            return
        fi
    fi
    
    echo "    Using control: $CONTROL_VCF"
    echo "    Using treatment files: $TREATMENT_VCFS"
    
    # Compare each treatment sample with control
    for TREAT_VCF in $TREATMENT_VCFS; do
        SAMPLE=$(basename "$TREAT_VCF" | sed 's/\.filtered\.vcf\.gz$\|\.fixed\.vcf\.gz$//')
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
        echo "$HIGH_CONF" > "results/merged/summary/${GROUP}_direct_highconf_count.txt"
    else
        echo "    No treatment-specific variants to analyze for $GROUP"
        echo "0" > "results/merged/summary/${GROUP}_direct_highconf_count.txt"
    fi
}

# Run direct comparisons for all groups - WITH CORRECTED PATHS
compare_treatment_to_control_direct "WT" "results/merged/filtered/WT-CTRL.filtered.vcf.gz" "WT-37-55"
compare_treatment_to_control_direct "STC" "results/merged/filtered/STC-CTRL.filtered.vcf.gz" "STC-55"
compare_treatment_to_control_direct "CAS" "results/merged/filtered/CAS-CTRL.filtered.vcf.gz" "CAS-55"
compare_treatment_to_control_direct "WTA" "results/merged/filtered/WT-CTRL.filtered.vcf.gz" "WTA-55"
# 8b. Analyze consistency of variants within replicates
echo -e "\nStep 8b: Analyzing consistency within treatment replicates..."

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
    echo "$GROUP consistency:" > "results/merged/summary/${GROUP}_consistency.txt"
    echo "  Variants in all $NUM_SAMPLES replicates: $ALL" >> "results/merged/summary/${GROUP}_consistency.txt"
    echo "  Variants in exactly 2 replicates: $((AT_LEAST_2 - ALL))" >> "results/merged/summary/${GROUP}_consistency.txt"
    echo "  Variants in only 1 replicate: $ONLY_1" >> "results/merged/summary/${GROUP}_consistency.txt"
    echo "  Total variants: $((AT_LEAST_2 + ONLY_1))" >> "results/merged/summary/${GROUP}_consistency.txt"
}

# Analyze consistency for each treatment group
analyze_group_consistency "WT" "WT-37-55-1,WT-37-55-2,WT-37-55-3"
analyze_group_consistency "STC" "STC-55-1,STC-55-2,STC-55-3"
analyze_group_consistency "CAS" "CAS-55-1,CAS-55-2,CAS-55-3"
analyze_group_consistency "WTA" "WTA-55-1,WTA-55-2,WTA-55-3"

# 8c. Add method comparison to the summary
echo -e "\nStep 8c: Comparing analysis methods..."

echo "Method Comparison:" > "results/merged/summary/method_comparison.txt"
for GROUP in WT STC CAS WTA; do
    if [ -f "results/merged/summary/${GROUP}_highconf_count.txt" ] && [ -f "results/merged/summary/${GROUP}_direct_highconf_count.txt" ]; then
        MERGED_COUNT=$(cat "results/merged/summary/${GROUP}_highconf_count.txt")
        DIRECT_COUNT=$(cat "results/merged/summary/${GROUP}_direct_highconf_count.txt")
        echo "  $GROUP: $MERGED_COUNT high-conf variants (merged), $DIRECT_COUNT high-conf variants (direct)" >> "results/merged/summary/method_comparison.txt"
    else
        echo "  $GROUP: Incomplete data" >> "results/merged/summary/method_comparison.txt"
    fi
done

# 9. Generate a summary report
echo -e "\nStep 9: Generating final summary report..."

echo "=== Variant Analysis Summary Report ===" > "results/merged/summary/analysis_report.txt"
echo "Date: $(date)" >> "results/merged/summary/analysis_report.txt"
echo "" >> "results/merged/summary/analysis_report.txt"

echo "Treatment-Specific Variants:" >> "results/merged/summary/analysis_report.txt"
for GROUP in WT STC CAS WTA; do
    if [ -f "results/merged/summary/${GROUP}_specific_count.txt" ] && [ -f "results/merged/summary/${GROUP}_highconf_count.txt" ]; then
        TOTAL=$(cat "results/merged/summary/${GROUP}_specific_count.txt")
        HC=$(cat "results/merged/summary/${GROUP}_highconf_count.txt")
        echo "  $GROUP: $TOTAL total, $HC high-confidence" >> "results/merged/summary/analysis_report.txt"
    else
        echo "  $GROUP: No data available" >> "results/merged/summary/analysis_report.txt"
    fi
done

echo "" >> "results/merged/summary/analysis_report.txt"
echo "Treatment-Unique Variants:" >> "results/merged/summary/analysis_report.txt"
cat "results/merged/summary/unique_variants.txt" | sed 's/^/  /' >> "results/merged/summary/analysis_report.txt"

# Add consistency information to the report
echo "" >> "results/merged/summary/analysis_report.txt"
echo "Variant Consistency Within Replicates:" >> "results/merged/summary/analysis_report.txt"
for GROUP in WT STC CAS WTA; do
    if [ -f "results/merged/summary/${GROUP}_consistency.txt" ]; then
        cat "results/merged/summary/${GROUP}_consistency.txt" | sed 's/^/  /' >> "results/merged/summary/analysis_report.txt"
        echo "" >> "results/merged/summary/analysis_report.txt"
    fi
done

# Add method comparison to the report
echo "Variant Analysis Method Comparison:" >> "results/merged/summary/analysis_report.txt"
cat "results/merged/summary/method_comparison.txt" | sed 's/^/  /' >> "results/merged/summary/analysis_report.txt"

echo -e "\nMerged VCF analysis complete with contig fix! Results in results/merged/ directory."
echo "Summary report available at: results/merged/summary/analysis_report.txt"