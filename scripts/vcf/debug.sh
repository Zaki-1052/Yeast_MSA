#!/bin/bash

echo "Detailed VCF header analysis..."

# Create debug directory
mkdir -p debug

# 1. Extract raw header line and show byte-by-byte hex dump
echo "Raw header bytes (hex dump):"
bcftools view -h results/vcf/merged/all_samples_fixed.vcf.gz | grep "^#CHROM" | xxd > debug/header_hex.txt
cat debug/header_hex.txt

# 2. Count fields in the header line
echo -e "\nCounting fields in header:"
FIELD_COUNT=$(bcftools view -h results/vcf/merged/all_samples_fixed.vcf.gz | grep "^#CHROM" | tr '\t' '\n' | wc -l)
echo "Field count: $FIELD_COUNT"

# 3. Extract sample names manually
echo -e "\nExtracting column names manually:"
bcftools view -h results/vcf/merged/all_samples_fixed.vcf.gz | grep "^#CHROM" | tr '\t' '\n' | tail -n +10 > debug/column_names.txt
cat debug/column_names.txt

# 4. Try alternative manual solution
echo -e "\nTrying direct extraction with cut:"
LINE=$(bcftools view -h results/vcf/merged/all_samples_fixed.vcf.gz | grep "^#CHROM")
echo "Sample names for manual solution:"
for i in {10..24}; do
    COL=$(echo "$LINE" | cut -f $i)
    echo "$i: $COL"
done > debug/manual_columns.txt
cat debug/manual_columns.txt

# 5. Test if sample file selection actually works with a different approach
echo -e "\nTesting sample file selection with sample count check:"
echo "WT-CTRL" > debug/wt_ctrl.txt
bcftools view -S debug/wt_ctrl.txt results/vcf/merged/all_samples_fixed.vcf.gz -Ov -o debug/wt_ctrl_subset.vcf
SAMPLE_COUNT=$(bcftools query -l debug/wt_ctrl_subset.vcf | wc -l)
echo "Number of samples in output VCF: $SAMPLE_COUNT"

echo "Detailed analysis complete."