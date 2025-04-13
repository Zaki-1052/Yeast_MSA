#!/bin/bash

echo "=== VCF Format Deep Dive ==="
mkdir -p vcf_debug

# 1. Verify basic file access
echo "Testing basic file access with bcftools:"
bcftools view -h results/vcf/filtered/WT-CTRL.filtered.vcf.gz > vcf_debug/test_header.txt
if [ $? -eq 0 ]; then
  echo "  SUCCESS: Can read header with bcftools"
  head -5 vcf_debug/test_header.txt
else
  echo "  FAILED: Cannot read header with bcftools"
fi

# 2. Check if we can decompress properly
echo -e "\nTesting decompression:"
gzip -dc results/vcf/filtered/WT-CTRL.filtered.vcf.gz > vcf_debug/test_decompressed.vcf 2>vcf_debug/decompress_error.txt
if [ $? -eq 0 ]; then
  echo "  SUCCESS: File decompressed correctly"
  head -5 vcf_debug/test_decompressed.vcf
else
  echo "  FAILED: Decompression error"
  cat vcf_debug/decompress_error.txt
fi

# 3. Check sample names in each file
echo -e "\nChecking sample names in each VCF:"
for VCF in results/vcf/filtered/*.filtered.vcf.gz; do
  SAMPLE=$(basename "$VCF" .filtered.vcf.gz)
  HEADER_SAMPLE=$(bcftools query -l "$VCF")
  echo "$SAMPLE → $HEADER_SAMPLE" >> vcf_debug/sample_names.txt
done
cat vcf_debug/sample_names.txt

# 4. Inspect contig lines in headers
echo -e "\nInspecting contig lines in headers:"
for VCF in $(ls results/vcf/filtered/*.filtered.vcf.gz | head -3); do
  SAMPLE=$(basename "$VCF" .filtered.vcf.gz)
  echo "=== $SAMPLE contigs ===" >> vcf_debug/contig_lines.txt
  bcftools view -h "$VCF" | grep "##contig" | head -3 >> vcf_debug/contig_lines.txt
done
head -10 vcf_debug/contig_lines.txt

# 5. Try a simpler merge approach with just 2 files
echo -e "\nTrying a simple merge with 2 files:"
bcftools merge results/vcf/filtered/WT-CTRL.filtered.vcf.gz results/vcf/filtered/WT-37-55-1.filtered.vcf.gz -Oz -o vcf_debug/two_sample_merge.vcf.gz
if [ $? -eq 0 ]; then
  echo "  SUCCESS: Two-file merge worked"
  bcftools query -l vcf_debug/two_sample_merge.vcf.gz
else
  echo "  FAILED: Even simple merge failed"
fi

# 6. Try using bgzip directly to ensure proper compression
echo -e "\nTesting proper BGZF compression:"
bcftools view results/vcf/filtered/WT-CTRL.filtered.vcf.gz | bgzip -c > vcf_debug/recompressed.vcf.gz
if [ $? -eq 0 ]; then
  echo "  SUCCESS: Recompression successful"
  tabix -p vcf vcf_debug/recompressed.vcf.gz
  if [ $? -eq 0 ]; then
    echo "  SUCCESS: Indexing successful"
  else
    echo "  FAILED: Couldn't index recompressed file"
  fi
else
  echo "  FAILED: Recompression failed"
fi

# 7. Check if there are contigs with different ordering
echo -e "\nChecking contig ordering across files:"
for VCF in $(ls results/vcf/filtered/*.filtered.vcf.gz | head -3); do
  SAMPLE=$(basename "$VCF" .filtered.vcf.gz)
  echo "=== $SAMPLE first 5 records ===" >> vcf_debug/record_order.txt
  bcftools view "$VCF" | grep -v "^#" | head -5 >> vcf_debug/record_order.txt
done
head -20 vcf_debug/record_order.txt

echo "=== VCF Format Investigation Complete ==="