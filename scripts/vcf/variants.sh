#!/bin/bash

# Create directories for variant calls
mkdir -p results/vcf/individual
mkdir -p results/vcf/filtered
mkdir -p results/vcf/merged

# Path to reference genome
REF="reference/yeast_w303.fasta"

# 1. Call variants for each sample individually
for BAM in results/bam/*.final.bam; do
    SAMPLE=$(basename "$BAM" .final.bam)
    echo "Calling variants for sample: $SAMPLE ($(date))"
    
    # Generate mpileup and call variants (using haploid model for yeast)
    bcftools mpileup -f "$REF" -a "DP,AD" "$BAM" | \
    bcftools call -mv --ploidy 1 -Ov -o "results/vcf/individual/${SAMPLE}.vcf"
    
    # Normalize indel representations
    bcftools norm -f "$REF" "results/vcf/individual/${SAMPLE}.vcf" -Ov -o "results/vcf/individual/${SAMPLE}.norm.vcf"
    
    # Filter variants by quality and depth - using FORMAT/DP to be specific
    bcftools filter -i 'QUAL>=20 && FORMAT/DP>=10' \
    -o "results/vcf/filtered/${SAMPLE}.filtered.vcf" \
    "results/vcf/individual/${SAMPLE}.norm.vcf"
    
    # Compress and index VCF files
    bgzip -f "results/vcf/filtered/${SAMPLE}.filtered.vcf"
    tabix -p vcf "results/vcf/filtered/${SAMPLE}.filtered.vcf.gz"
    
    # Generate per-sample variant statistics
    bcftools stats "results/vcf/filtered/${SAMPLE}.filtered.vcf.gz" > "results/vcf/filtered/${SAMPLE}.stats.txt"
    
    echo "Completed variant calling for $SAMPLE ($(date))"
done

# 2. Create a merged VCF file with all samples
echo "Merging VCF files from all samples ($(date))"

# Create a list of VCF files to merge
ls results/vcf/filtered/*.filtered.vcf.gz > results/vcf/vcf_list.txt

# Merge VCFs
bcftools merge -l results/vcf/vcf_list.txt -Oz -o "results/vcf/merged/all_samples.vcf.gz"

# Index the merged file
tabix -p vcf "results/vcf/merged/all_samples.vcf.gz"

# 3. Generate comprehensive variant statistics
bcftools stats "results/vcf/merged/all_samples.vcf.gz" > "results/vcf/merged/variant_stats.txt"

# 4. Basic summary counts
echo -e "Sample\tTotal_Variants\tSNPs\tIndels" > "results/vcf/variant_summary.txt"
for STATS in results/vcf/filtered/*.stats.txt; do
    SAMPLE=$(basename "$STATS" .stats.txt)
    VARIANTS=$(grep "number of records:" "$STATS" | awk '{print $4}')
    SNPS=$(grep "number of SNPs:" "$STATS" | awk '{print $4}')
    INDELS=$(grep "number of indels:" "$STATS" | awk '{print $4}')
    
    echo -e "${SAMPLE%.filtered}\t$VARIANTS\t$SNPS\t$INDELS" >> "results/vcf/variant_summary.txt"
done

echo "Variant calling complete. Results available in results/vcf/ directory."
echo "Merged VCF file: results/vcf/merged/all_samples.vcf.gz"
echo "Summary statistics: results/vcf/variant_summary.txt"