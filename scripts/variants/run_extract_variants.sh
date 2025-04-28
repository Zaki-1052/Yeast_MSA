#!/bin/bash

# run_extract_variants.sh - Extract variants affecting ergosterol pathway genes

# Check if the required directories exist
if [ ! -d "vcf/annotated" ]; then
    echo "Error: Directory vcf/annotated does not exist"
    exit 1
fi

if [ ! -f "reference/genes_of_interest_mapping.tsv" ]; then
    echo "Error: File reference/genes_of_interest_mapping.tsv does not exist"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p results/gene_variants

# Run the Python script
echo "Extracting variants affecting ergosterol pathway genes..."
python scripts/extract_gene_variants.py \
  --vcf_dir vcf/annotated \
  --gene_mapping reference/genes_of_interest_mapping.tsv \
  --output_dir results/gene_variants

# Check if the extraction was successful
if [ $? -eq 0 ]; then
    echo "Extraction complete. Results are in results/gene_variants/"
else
    echo "Error: Extraction failed"
    exit 1
fi

# Display summary info
if [ -f "results/gene_variants/gene_summary.tsv" ]; then
    echo -e "\nGene Summary:"
    column -t -s $'\t' results/gene_variants/gene_summary.tsv
fi

echo -e "\nDone!"