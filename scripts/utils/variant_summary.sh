#!/bin/bash

# Create fixed variant summary
echo -e "Sample\tTotal_Variants\tSNPs\tIndels" > results/vcf/variant_summary_fixed.txt

for VCF in results/vcf/filtered/*.filtered.vcf.gz; do
    SAMPLE=$(basename "$VCF" .filtered.vcf.gz)
    
    # Count variants directly from VCF files
    TOTAL=$(bcftools view -H "$VCF" | wc -l)
    SNPS=$(bcftools view -H -v snps "$VCF" | wc -l)
    INDELS=$(bcftools view -H -v indels "$VCF" | wc -l)
    
    echo -e "$SAMPLE\t$TOTAL\t$SNPS\t$INDELS" >> results/vcf/variant_summary_fixed.txt
done

cat results/vcf/variant_summary_fixed.txt