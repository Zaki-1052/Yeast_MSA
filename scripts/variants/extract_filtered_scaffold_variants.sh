#!/bin/bash

# extract_filtered_scaffold_variants.sh - Extract treatment-specific variants from scaffolds containing ERG genes

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
mkdir -p results/filtered_scaffold_variants

# Run the Python script
echo "Extracting treatment-specific variants from scaffolds containing ergosterol pathway genes..."
python3 scripts/variants/extract_filtered_scaffold_variants.py \
  --vcf_dir vcf/annotated \
  --gene_mapping reference/genes_of_interest_mapping.tsv \
  --output_dir results/filtered_scaffold_variants \
  --distance_bins 0,500,1000,5000,10000,50000 \
  --max_distance 50000

# Check if the extraction was successful
if [ $? -eq 0 ]; then
    echo "Extraction complete. Results are in results/filtered_scaffold_variants/"
    
    # Display summary
    if [ -f "results/filtered_scaffold_variants/treatment_specific_scaffold_variant_summary.txt" ]; then
        echo -e "\nSummary:"
        head -n 20 results/filtered_scaffold_variants/treatment_specific_scaffold_variant_summary.txt
        echo -e "\n(See full report for more details)"
    fi
else
    echo "Error: Extraction failed"
    exit 1
fi

echo -e "\nDone!"