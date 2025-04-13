#!/bin/bash

echo "Testing alternative sample selection methods..."

# Create sample lists as files
echo "WT-CTRL" > wt_ctrl.txt
echo "WT-37-55-1" > wt_treat1.txt
echo "WT-37-55-2" >> wt_treat1.txt
echo "WT-37-55-3" >> wt_treat1.txt

# 1. Test with a sample file instead of comma-separated list
echo "Testing with sample file:"
bcftools view -S wt_ctrl.txt results/vcf/merged/all_samples_fixed.vcf.gz -Ov -o /dev/null
if [ $? -eq 0 ]; then
    echo "  Success: Sample file selection works"
else
    echo "  Failed: Sample file selection doesn't work"
fi

# 2. Try with --force-samples to see if that helps
echo "Testing with --force-samples:"
bcftools view --force-samples -s WT-CTRL,WT-37-55-1 results/vcf/merged/all_samples_fixed.vcf.gz -Ov -o /dev/null
if [ $? -eq 0 ]; then
    echo "  Success: Force-samples option works"
else
    echo "  Failed: Even force-samples doesn't work"
fi

# 3. Examine the header directly
echo "Examining fixed VCF header:"
bcftools view -h results/vcf/merged/all_samples_fixed.vcf.gz | grep "^#CHROM" | sed 's/\t/|/g'
echo ""

# 4. Try with a completely different sample selection approach
echo "Testing alternative selection method:"
bcftools view results/vcf/merged/all_samples_fixed.vcf.gz | \
  awk -v OFS="\t" '{if($1 ~ /^#CHROM/) {for(i=1;i<=NF;i++) if($i=="WT-CTRL") ctrlcol=i; print} else if(!($1 ~ /^#/)) {print}}' | \
  head -20 > test_awkselect.vcf

WT_AWKCOUNT=$(grep -v "^#" test_awkselect.vcf | wc -l)
echo "  Direct extraction sample count: $WT_AWKCOUNT"

# Clean up
rm -f wt_ctrl.txt wt_treat1.txt test_awkselect.vcf

echo "Advanced testing complete."