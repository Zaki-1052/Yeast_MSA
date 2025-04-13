#!/bin/bash

echo "=== Fixing contig definitions in VCF headers ==="

# Extract contig information from reference genome
echo "  Extracting contig information from reference..."

# Generate FAI index if it doesn't exist
if [ ! -f "reference/yeast_w303.fasta.fai" ]; then
    echo "  Creating reference genome index..."
    samtools faidx reference/yeast_w303.fasta
fi

# Create a proper contig header from the reference index
echo "  Creating contig header lines..."
awk 'BEGIN {print "##fileformat=VCFv4.2"} 
     {print "##contig=<ID="$1",length="$2">"}' reference/yeast_w303.fasta.fai > \
     results/merged/fixed/contig_header.txt

# Add proper contig definitions to each filtered VCF
echo "  Adding proper contig definitions to VCFs..."

for FILT_VCF in results/merged/filtered/*.filtered.vcf.gz; do
    SAMPLE=$(basename "$FILT_VCF" .filtered.vcf.gz)
    echo "  Fixing header for $SAMPLE..."
    
    # Extract the existing header without contig lines
    bcftools view -h "$FILT_VCF" | grep -v "^##contig" | grep -v "^##fileformat" > \
      results/merged/fixed/existing_header_${SAMPLE}.txt
    
    # Combine the headers
    cat results/merged/fixed/contig_header.txt \
        results/merged/fixed/existing_header_${SAMPLE}.txt > \
        results/merged/fixed/complete_header_${SAMPLE}.txt
    
    # Replace the header
    bcftools reheader -h results/merged/fixed/complete_header_${SAMPLE}.txt \
      -o results/merged/fixed/${SAMPLE}.fixed.vcf.gz \
      "$FILT_VCF"
    
    # Index the fixed VCF
    tabix -p vcf results/merged/fixed/${SAMPLE}.fixed.vcf.gz
done