#!/bin/bash

echo "Fresh merge attempt..."

# 1. Create a simple sample name list (no special characters)
mkdir -p fresh_merge
for SAMPLE in $(bcftools query -l results/vcf/merged/all_samples_fixed.vcf.gz); do
    echo "$SAMPLE" >> fresh_merge/sample_names.txt
done

# 2. Create a brand new merged VCF with explicit sample naming
bcftools merge --force-samples results/vcf/filtered/*.filtered.vcf.gz -Oz -o fresh_merge/freshly_merged.vcf.gz

# 3. Test with a simple sample selection
echo "Testing sample selection on fresh merge:"
bcftools view -s WT-CTRL fresh_merge/freshly_merged.vcf.gz -Ov -o /dev/null
if [ $? -eq 0 ]; then
    echo "SUCCESS: Sample selection works on fresh merge!"
else
    echo "FAILED: Even a fresh merge has sample selection issues"
fi

echo "Fresh merge test complete"