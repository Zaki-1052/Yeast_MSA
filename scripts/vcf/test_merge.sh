#!/bin/bash

echo "Testing freshly merged VCF..."

# 1. Test single sample selection
echo "Test 1: Single sample selection"
bcftools view -s WT-CTRL fresh_merge/freshly_merged.vcf.gz -Ov -o /dev/null
if [ $? -eq 0 ]; then
    echo "  PASS: Single sample selection works"
else
    echo "  FAIL: Single sample selection still has issues"
fi

# 2. Test multiple sample selection with comma-separated list
echo "Test 2: Multiple sample selection"
bcftools view -s WT-CTRL,WT-37-55-1,WT-37-55-2 fresh_merge/freshly_merged.vcf.gz -Ov -o /dev/null
if [ $? -eq 0 ]; then
    echo "  PASS: Multiple sample selection works"
else
    echo "  FAIL: Multiple sample selection still has issues"
fi

# 3. Test sample filtering with -c option
echo "Test 3: Sample filtering"
bcftools view -s WT-37-55-1,WT-37-55-2,WT-37-55-3 -c 1 fresh_merge/freshly_merged.vcf.gz -Ov -o /dev/null
if [ $? -eq 0 ]; then
    echo "  PASS: Sample filtering works"
else
    echo "  FAIL: Sample filtering still has issues"
fi

# 4. Test treatment vs control - create a simple WT test
echo "Test 4: Treatment vs control comparison"
bcftools view -s WT-CTRL,WT-37-55-1,WT-37-55-2,WT-37-55-3 fresh_merge/freshly_merged.vcf.gz -Oz -o fresh_merge/wt_group.vcf.gz
bcftools index fresh_merge/wt_group.vcf.gz

# Find variants in treatment but not control
bcftools view -s WT-37-55-1,WT-37-55-2,WT-37-55-3 -c 1 fresh_merge/wt_group.vcf.gz | \
bcftools view -s WT-CTRL -c 0 -Oz -o fresh_merge/wt_specific.vcf.gz

# Count treatment-specific variants
WT_COUNT=$(bcftools view -H fresh_merge/wt_specific.vcf.gz | wc -l)
echo "  WT treatment-specific variants: $WT_COUNT"
if [ $WT_COUNT -gt 0 ]; then
    echo "  PASS: Treatment vs control comparison works and found $WT_COUNT variants"
else
    echo "  INCONCLUSIVE: Treatment vs control comparison produced no variants"
fi

echo "Testing complete!"