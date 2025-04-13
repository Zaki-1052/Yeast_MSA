#!/bin/bash

echo "=== Creating filtered VCF files ==="

# Filter, compress, and index each normalized VCF
for NORM_VCF in results/vcf/individual/*.norm.vcf; do
    SAMPLE=$(basename "$NORM_VCF" .norm.vcf)
    echo "  Processing $SAMPLE..."
    
    # Skip if already filtered
    if [ -f "results/merged/filtered/${SAMPLE}.filtered.vcf.gz" ]; then
        echo "    Already filtered, skipping..."
        continue
    fi
    
    # Filter by quality and read depth
    bcftools filter -i 'QUAL>=20 && FORMAT/DP>=10' \
      -o "results/merged/filtered/${SAMPLE}.filtered.vcf" \
      "$NORM_VCF"
    
    # Compress and index
    bgzip -f "results/merged/filtered/${SAMPLE}.filtered.vcf"
    tabix -p vcf "results/merged/filtered/${SAMPLE}.filtered.vcf.gz"
done