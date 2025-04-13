#!/bin/bash

echo "=== Complete Fix for Both Issues ==="

# Part 1: Fix the high-confidence variant counting
echo "Fixing high-confidence variant counting..."
mkdir -p fixed_results

# Process each treatment group
for GROUP in WT STC CAS WTA; do
  echo "Processing $GROUP high-confidence variants..."
  
  # Create directories
  mkdir -p fixed_results/$GROUP
  
  # Extract positions from each sample's specific variants
  grep -v "^#" final_analysis/$GROUP/*_specific.vcf | \
  awk '{print $1":"$2":"$4":"$5}' | sort | uniq -c | sort -nr > \
  fixed_results/$GROUP/variant_counts.txt
  
  # Count variants present in at least 2 samples (count ≥ 2)
  HC_COUNT=$(awk '$1 >= 2' fixed_results/$GROUP/variant_counts.txt | wc -l)
  echo "  $GROUP high-confidence variants: $HC_COUNT"
  
  # Output detailed list of high-confidence variants
  awk '$1 >= 2 {
    split($2,parts,":");
    print parts[1] "\t" parts[2] "\t" parts[3] "\t" parts[4] "\t" $1
  }' fixed_results/$GROUP/variant_counts.txt > fixed_results/$GROUP/high_conf_variants.tsv
done

# Part 2: Fix the VCF selection and comparison
echo -e "\nFixing merged VCF sample selection with --force-samples..."

# Check sample names again
bcftools query -l fixed_vcfs/properly_merged.vcf.gz > fixed_results/sample_names.txt

# Try WT test selection with --force-samples
echo "Testing WT treatment vs control comparison with --force-samples:"
bcftools view --force-samples -s WT-CTRL,WT-37-55-1,WT-37-55-2,WT-37-55-3 \
  fixed_vcfs/properly_merged.vcf.gz -Oz -o fixed_results/wt_group.vcf.gz
  
tabix -p vcf fixed_results/wt_group.vcf.gz

# Find variants in treatment but not control using --force-samples
bcftools view --force-samples -s WT-37-55-1,WT-37-55-2,WT-37-55-3 -c 1 \
  fixed_results/wt_group.vcf.gz | \
bcftools view --force-samples -s WT-CTRL -c 0 -Oz -o fixed_results/wt_specific.vcf.gz

# Index the specific variants file
tabix -p vcf fixed_results/wt_specific.vcf.gz

# Find high-confidence variants using --force-samples
bcftools view --force-samples -s WT-37-55-1,WT-37-55-2,WT-37-55-3 -c 2 \
  fixed_results/wt_specific.vcf.gz -Oz -o fixed_results/wt_highconf.vcf.gz

# Count results
TOTAL=$(bcftools view -H fixed_results/wt_specific.vcf.gz | wc -l)
HIGH=$(bcftools view -H fixed_results/wt_highconf.vcf.gz | wc -l)

echo "WT treatment-specific variants: $TOTAL"
echo "WT high-confidence variants: $HIGH"

echo "=== Complete Fix Done ==="