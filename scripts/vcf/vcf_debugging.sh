#!/bin/bash

echo "=== Methodical Debugging of VCF Issues ==="

# 1. Examine the exact sample names in our merged file
echo "Examining merged file sample names:"
bcftools query -l fixed_vcfs/properly_merged.vcf.gz > sample_list.txt
cat sample_list.txt

# 2. Check header structure of first few lines
echo -e "\nExamining merged file header structure:"
bcftools view -h fixed_vcfs/properly_merged.vcf.gz | tail -5

# 3. Check if there are quote issues with sample names
echo -e "\nChecking for potential quote issues in sample names:"
bcftools view -h fixed_vcfs/properly_merged.vcf.gz | grep "^#CHROM" | xxd

# 4. Try with --force-samples to see if that helps
echo -e "\nTesting with --force-samples:"
bcftools view --force-samples -s WT-CTRL fixed_vcfs/properly_merged.vcf.gz | head -5 > forced_test.txt
if [ $? -eq 0 ]; then
  echo "SUCCESS: --force-samples allows selection"
else
  echo "FAILED: Even force-samples doesn't work"
fi

# 5. Check for invisible characters or encoding issues
echo -e "\nChecking for encoding issues:"
echo "WT-CTRL" > expected.txt
bcftools query -l fixed_vcfs/properly_merged.vcf.gz | grep "WT-CTRL" > actual.txt
diff -u expected.txt actual.txt
hexdump -C expected.txt
hexdump -C actual.txt