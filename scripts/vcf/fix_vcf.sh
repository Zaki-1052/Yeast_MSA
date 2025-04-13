#!/bin/bash

echo "=== Contig Order Fixing and Proper Merging ==="
mkdir -p vcf_fixed

# 1. Create a master sequence dictionary from the reference
echo "Creating master sequence dictionary..."
bcftools view -h results/vcf/filtered/WT-CTRL.filtered.vcf.gz | grep "^##contig" > vcf_fixed/master_contigs.txt

# 2. Reorder each VCF to match this master ordering
echo "Reordering VCFs to consistent contig ordering..."

# Process each VCF file
for VCF in results/vcf/filtered/*.filtered.vcf.gz; do
  SAMPLE=$(basename "$VCF" .filtered.vcf.gz)
  echo "  Processing $SAMPLE..."
  
  # Create a new header with consistent contig ordering
  bcftools view -h "$VCF" | grep -v "^##contig" > vcf_fixed/${SAMPLE}_header.txt
  cat vcf_fixed/master_contigs.txt >> vcf_fixed/${SAMPLE}_header.txt
  bcftools view -h "$VCF" | grep "^#CHROM" >> vcf_fixed/${SAMPLE}_header.txt
  
  # Apply the fixed header and sort
  bcftools reheader -h vcf_fixed/${SAMPLE}_header.txt "$VCF" | \
  bcftools sort -Oz -o vcf_fixed/${SAMPLE}.vcf.gz
  
  # Index the fixed VCF
  bcftools index vcf_fixed/${SAMPLE}.vcf.gz
done

# 3. Create a list of reordered VCFs
ls vcf_fixed/*.vcf.gz > vcf_fixed/vcf_list.txt

# 4. Merge the reordered VCFs
echo "Merging VCFs with consistent contig ordering..."
bcftools merge -l vcf_fixed/vcf_list.txt -Oz -o vcf_fixed/all_samples_fixed.vcf.gz

# 5. Index the merged file
bcftools index vcf_fixed/all_samples_fixed.vcf.gz

# 6. Test a simple sample selection on merged file
echo "Testing sample selection on properly merged file:"
bcftools view -s WT-CTRL,WT-37-55-1 vcf_fixed/all_samples_fixed.vcf.gz -Ov -o /dev/null
if [ $? -eq 0 ]; then
  echo "  SUCCESS: Sample selection works on properly merged file!"
else
  echo "  FAILED: Sample selection still fails on merged file"
fi

# 7. Final verification - try a comparison
echo "Testing treatment vs control comparison:"
bcftools view -s WT-CTRL,WT-37-55-1,WT-37-55-2,WT-37-55-3 vcf_fixed/all_samples_fixed.vcf.gz -Oz -o vcf_fixed/wt_group.vcf.gz
bcftools index vcf_fixed/wt_group.vcf.gz

# Find variants in treatment but not control
bcftools view -s WT-37-55-1,WT-37-55-2,WT-37-55-3 -c 1 vcf_fixed/wt_group.vcf.gz | \
bcftools view -s WT-CTRL -c 0 -Oz -o vcf_fixed/wt_specific.vcf.gz

# Count treatment-specific variants
WT_COUNT=$(bcftools view -H vcf_fixed/wt_specific.vcf.gz | wc -l)
echo "  WT treatment-specific variants: $WT_COUNT"

echo "=== VCF Fix Complete ==="