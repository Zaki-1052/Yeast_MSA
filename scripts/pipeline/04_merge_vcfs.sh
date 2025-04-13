#!/bin/bash

echo "=== Merging VCFs with proper contig definitions ==="

# Create a list of fixed VCFs to merge
ls results/merged/fixed/*.fixed.vcf.gz > results/merged/fixed/vcf_list.txt

# Merge the VCFs with proper headers
bcftools merge -l results/merged/fixed/vcf_list.txt \
  -Oz -o results/merged/fixed/all_samples.vcf.gz

# Index the merged file
tabix -p vcf results/merged/fixed/all_samples.vcf.gz