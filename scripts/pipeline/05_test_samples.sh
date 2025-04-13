#!/bin/bash

echo "=== Verifying sample selection works with fixed VCF ==="

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
        echo "FORCE_FLAG=\"\"" > results/merged/fixed/force_flag.txt
    else
        echo "  WARNING: Multi-sample selection failed, using --force-samples"
        echo "FORCE_FLAG=\"--force-samples\"" > results/merged/fixed/force_flag.txt
    fi
else
    echo "  WARNING: Sample selection not working, using --force-samples"
    echo "FORCE_FLAG=\"--force-samples\"" > results/merged/fixed/force_flag.txt
fi