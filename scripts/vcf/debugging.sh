#!/bin/bash

# Create a debug directory
mkdir -p debug

# 1. Extract the exact VCF header to examine how samples are defined
bcftools view -h results/vcf/merged/all_samples.vcf.gz > debug/header.txt

# 2. Check for hidden characters in sample names
bcftools query -l results/vcf/merged/all_samples.vcf.gz | xxd > debug/sample_hex.txt

# 3. Try with --force-samples flag
echo "Testing with --force-samples flag on WT group..."
bcftools view --force-samples -s "WT-CTRL,WT-37-55-1,WT-37-55-2,WT-37-55-3" \
  results/vcf/merged/all_samples.vcf.gz -Oz -o debug/wt_group_force.vcf.gz

# 4. Check if any output was created
if [ -s debug/wt_group_force.vcf.gz ]; then
  echo "  Success! Output file created with --force-samples"
  bcftools query -l debug/wt_group_force.vcf.gz
else
  echo "  Failed! No output created even with --force-samples"
fi

# 5. Try creating a file with sample renaming
echo "Testing with explicit sample renaming..."
bcftools reheader -s <(bcftools query -l results/vcf/merged/all_samples.vcf.gz) \
  -o debug/renamed.vcf.gz results/vcf/merged/all_samples.vcf.gz