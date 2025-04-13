#!/bin/bash

echo "Testing individual filtered VCFs..."

# Pick a sample to test
TEST_VCF="results/vcf/filtered/WT-CTRL.filtered.vcf.gz"

# 1. Check if the header has sample information
echo "Examining header of individual filtered VCF:"
bcftools view -h $TEST_VCF | grep "^#CHROM" | sed 's/\t/|/g'

# 2. Check if sample selection works on an individual filtered VCF
echo "Testing sample selection on individual filtered VCF:"
bcftools view -s WT-CTRL $TEST_VCF -Ov -o /dev/null
if [ $? -eq 0 ]; then
    echo "  PASS: Sample selection works on individual filtered VCF"
else
    echo "  FAIL: Sample selection fails even on individual filtered VCF"
fi

# 3. Check FORMAT fields to see if DP is present
echo "Checking FORMAT fields in the filtered VCF:"
bcftools view $TEST_VCF | grep -v "^#" | head -1 | cut -f9- | sed 's/\t/|/g'

echo "Individual VCF test complete"