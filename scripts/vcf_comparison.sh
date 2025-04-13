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

echo -e "\nMerged VCF analysis complete with contig fix! Results in results/merged/ directory."
echo "Summary report available at: results/merged/summary/analysis_report.txt"