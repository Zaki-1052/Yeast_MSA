#!/bin/bash

# Create directories for functional analysis
mkdir -p results/functional
mkdir -p results/visualization

# 1. Extract high-confidence variants for detailed analysis
for GROUP in WT STC CAS WTA; do
    # Create analysis directory for this group
    mkdir -p "results/functional/${GROUP}"
    
    # Extract variants by location (for easier analysis)
    echo "Extracting variant details for $GROUP"
    
    # Get list of variants from the variant_counts file
    if [ -f "results/vcf/direct_comparison/${GROUP}/variant_counts.txt" ]; then
        # Extract variants that appear in at least 2 samples (high confidence)
        cat "results/vcf/direct_comparison/${GROUP}/variant_counts.txt" | awk '$1 >= 2 {print $2}' > "results/functional/${GROUP}/highconf_variants.txt"
        
        # Count variants per scaffold
        cut -d: -f1 "results/functional/${GROUP}/highconf_variants.txt" | sort | uniq -c | sort -nr > "results/functional/${GROUP}/scaffold_distribution.txt"
        
        # Extract detailed variant info from VCF files
        echo -e "Scaffold\tPosition\tRef\tAlt\tQuality\tDepth" > "results/functional/${GROUP}/variant_details.tsv"
        
        # Process each specific variant file to extract details
        for VCF in results/vcf/direct_comparison/${GROUP}/*_specific.vcf; do
            SAMPLE=$(basename "$VCF" _specific.vcf)
            grep -v "^#" "$VCF" | awk -v sample="$SAMPLE" '{
                split($8, info, ";");
                dp = "";
                for (i in info) {
                    if (info[i] ~ /^DP=/) {
                        split(info[i], dp_val, "=");
                        dp = dp_val[2];
                        break;
                    }
                }
                print $1"\t"$2"\t"$4"\t"$5"\t"$6"\t"dp;
            }' >> "results/functional/${GROUP}/variant_details.tsv"
        done
    else
        echo "No variant count data found for $GROUP"
    fi
done

# 2. Generate comparative statistics
echo -e "Group\tVariant_Count\tSNPs\tIndels\tAvg_Quality\tMedian_Depth" > results/functional/group_comparison.tsv

for GROUP in WT STC CAS WTA; do
    if [ -f "results/functional/${GROUP}/variant_details.tsv" ]; then
        # Count total variants (excluding header)
        TOTAL=$(wc -l < "results/functional/${GROUP}/variant_details.tsv")
        TOTAL=$((TOTAL - 1))
        
        # Count SNPs vs indels (simple heuristic: ref and alt are same length for SNPs)
        SNPS=$(awk 'NR>1 && length($3) == length($4) {count++} END {print count}' "results/functional/${GROUP}/variant_details.tsv")
        INDELS=$((TOTAL - SNPS))
        
        # Calculate average quality
        AVG_QUAL=$(awk 'NR>1 {sum+=$5; count++} END {print sum/count}' "results/functional/${GROUP}/variant_details.tsv")
        
        # Get median depth (simplified - not a true median calculation)
        MEDIAN_DP=$(awk 'NR>1 && $6!="" {print $6}' "results/functional/${GROUP}/variant_details.tsv" | sort -n | awk '{a[NR]=$1} END {print a[int(NR/2)]}')
        
        echo -e "${GROUP}\t${TOTAL}\t${SNPS}\t${INDELS}\t${AVG_QUAL}\t${MEDIAN_DP}" >> results/functional/group_comparison.tsv
    else
        echo -e "${GROUP}\tNA\tNA\tNA\tNA\tNA" >> results/functional/group_comparison.tsv
    fi
done

# 3. Find shared and unique variants between treatment groups
mkdir -p results/functional/overlap

# Create sets of variant locations for comparison
for GROUP in WT STC CAS WTA; do
    if [ -f "results/functional/${GROUP}/highconf_variants.txt" ]; then
        cp "results/functional/${GROUP}/highconf_variants.txt" "results/functional/overlap/${GROUP}_variants.txt"
    fi
done

# Find overlaps between pairs of treatments
for GROUP1 in WT STC CAS WTA; do
    for GROUP2 in WT STC CAS WTA; do
        if [ "$GROUP1" != "$GROUP2" ] && [ -f "results/functional/overlap/${GROUP1}_variants.txt" ] && [ -f "results/functional/overlap/${GROUP2}_variants.txt" ]; then
            # Find common variants
            comm -12 <(sort "results/functional/overlap/${GROUP1}_variants.txt") <(sort "results/functional/overlap/${GROUP2}_variants.txt") > "results/functional/overlap/${GROUP1}_${GROUP2}_common.txt"
            
            # Count common variants
            COUNT=$(wc -l < "results/functional/overlap/${GROUP1}_${GROUP2}_common.txt")
            echo "${GROUP1}-${GROUP2} shared variants: $COUNT"
        fi
    done
done

echo "Functional analysis complete. Results available in results/functional/ directory."