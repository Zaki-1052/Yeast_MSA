#!/bin/bash

echo "=== Using BCFtools Native Header Creation ==="

mkdir -p results/native_fix

# 1. Create a properly formatted sequence dictionary
echo "Step 1: Creating proper sequence dictionary..."
samtools dict reference/yeast_w303.fasta -o results/native_fix/yeast_w303.dict

# 2. Use bcftools to add headers correctly
echo "Step 2: Adding reference dictionary to VCFs..."

for VCF in results/merged/filtered/*.filtered.vcf.gz; do
    SAMPLE=$(basename "$VCF" .filtered.vcf.gz)
    echo "  Processing $SAMPLE..."
    
    # Use bcftools annotate to add contig lines properly
    bcftools annotate --fai reference/yeast_w303.fasta.fai \
      --header-lines results/native_fix/yeast_w303.dict \
      "$VCF" -Oz -o results/native_fix/${SAMPLE}.fixed.vcf.gz
    
    # Index
    tabix -p vcf results/native_fix/${SAMPLE}.fixed.vcf.gz
done

# 3. Create a list of properly fixed VCFs
echo "Step 3: Merging properly fixed VCFs..."
ls results/native_fix/*.fixed.vcf.gz > results/native_fix/vcf_list.txt

# 4. Merge using bcftools native operations
bcftools merge -l results/native_fix/vcf_list.txt \
  -Oz -o results/native_fix/all_samples.vcf.gz

# 5. Index the merged file
tabix -p vcf results/native_fix/all_samples.vcf.gz

# 6. Test sample selection
echo "Step 6: Testing sample selection..."
bcftools view -s WT-CTRL results/native_fix/all_samples.vcf.gz \
  -Ov -o results/native_fix/test_selection.vcf

# Check if file was created and has content
if [ -s results/native_fix/test_selection.vcf ]; then
    echo "  SUCCESS: Sample selection works!"
    
    # Test multi-sample selection
    echo "  Testing multi-sample selection..."
    bcftools view -s WT-CTRL,WT-37-55-1 results/native_fix/all_samples.vcf.gz \
      -Ov -o results/native_fix/test_multi.vcf
    
    if [ -s results/native_fix/test_multi.vcf ]; then
        echo "  SUCCESS: Multi-sample selection works!"
        echo -e "\nReady to proceed with full comparison using native tools."
    else
        echo "  FAILED: Multi-sample selection failed."
    fi
else
    echo "  FAILED: Sample selection still not working."
fi

echo "=== Header Fixing Complete ==="