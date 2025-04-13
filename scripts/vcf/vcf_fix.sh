#!/bin/bash

echo "=== Fixing VCF Header Format - Final Attempt ==="
mkdir -p fixed_vcfs

# 1. Extract contig lines in the correct format
echo "Extracting properly formatted contig lines..."
bcftools view -h results/vcf/filtered/WT-CTRL.filtered.vcf.gz | grep "^##contig" > fixed_vcfs/contig_lines.txt

# 2. Fix each VCF file with a properly formatted header
for VCF in results/vcf/filtered/*.filtered.vcf.gz; do
  SAMPLE=$(basename "$VCF" .filtered.vcf.gz)
  echo "Fixing $SAMPLE..."
  
  # Extract header without contig lines
  bcftools view -h "$VCF" | grep -v "^##contig" > fixed_vcfs/${SAMPLE}_header_part1.txt
  
  # Extract #CHROM line
  bcftools view -h "$VCF" | grep "^#CHROM" > fixed_vcfs/${SAMPLE}_header_part3.txt
  
  # Combine to create a fixed header
  cat fixed_vcfs/${SAMPLE}_header_part1.txt fixed_vcfs/contig_lines.txt fixed_vcfs/${SAMPLE}_header_part3.txt > fixed_vcfs/${SAMPLE}_fixed_header.txt
  
  # Get data portion of the VCF
  bcftools view -H "$VCF" > fixed_vcfs/${SAMPLE}_data.vcf
  
  # Combine header with data
  cat fixed_vcfs/${SAMPLE}_fixed_header.txt fixed_vcfs/${SAMPLE}_data.vcf > fixed_vcfs/${SAMPLE}_fixed.vcf
  
  # Compress and index
  bgzip -f fixed_vcfs/${SAMPLE}_fixed.vcf
  tabix -p vcf fixed_vcfs/${SAMPLE}_fixed.vcf.gz
done

# 3. Create list of fixed VCFs
ls fixed_vcfs/*_fixed.vcf.gz > fixed_vcfs/vcf_list.txt

# 4. Merge fixed VCFs
echo "Merging fixed VCFs..."
bcftools merge -l fixed_vcfs/vcf_list.txt -Oz -o fixed_vcfs/properly_merged.vcf.gz

# 5. Index the merged file
tabix -p vcf fixed_vcfs/properly_merged.vcf.gz

# 6. Test sample selection on merged file
echo "Testing sample selection on properly merged file:"
bcftools view -s WT-CTRL,WT-37-55-1 fixed_vcfs/properly_merged.vcf.gz -Ov -o /dev/null
if [ $? -eq 0 ]; then
  echo "SUCCESS: Sample selection works on properly merged file!"
  
  # 7. Test high-confidence variant identification
  echo "Testing treatment vs control comparison:"
  bcftools view -s WT-CTRL,WT-37-55-1,WT-37-55-2,WT-37-55-3 fixed_vcfs/properly_merged.vcf.gz -Oz -o fixed_vcfs/wt_group.vcf.gz
  tabix -p vcf fixed_vcfs/wt_group.vcf.gz
  
  # Find variants present in at least one treatment but not in control
  bcftools view -s WT-37-55-1,WT-37-55-2,WT-37-55-3 -c 1 fixed_vcfs/wt_group.vcf.gz | \
    bcftools view -s WT-CTRL -c 0 -Oz -o fixed_vcfs/wt_specific.vcf.gz
  
  # Find high-confidence variants (present in at least 2 treatment samples)
  bcftools view -s WT-37-55-1,WT-37-55-2,WT-37-55-3 -c 2 fixed_vcfs/wt_specific.vcf.gz -Oz -o fixed_vcfs/wt_highconf.vcf.gz
  
  # Count results
  TOTAL=$(bcftools view -H fixed_vcfs/wt_specific.vcf.gz | wc -l)
  HIGH=$(bcftools view -H fixed_vcfs/wt_highconf.vcf.gz | wc -l)
  
  echo "WT treatment-specific variants: $TOTAL"
  echo "WT high-confidence variants: $HIGH"
else
  echo "FAILED: Sample selection still fails on merged file"
fi

echo "=== Fix Complete ==="