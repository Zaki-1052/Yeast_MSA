#!/bin/bash

# verify_gene_coordinates.sh - Verify gene annotation coordinates

# Check if the required directories exist
if [ ! -d "reference/w303_annotations" ]; then
    echo "Error: Directory reference/w303_annotations does not exist"
    exit 1
fi

if [ ! -f "reference/genes_of_interest_mapping.tsv" ]; then
    echo "Error: File reference/genes_of_interest_mapping.tsv does not exist"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p results/gene_verification

# Run the Python script
echo "Verifying gene coordinates..."
python scripts/verify_gene_coordinates.py \
  --gene_mapping reference/genes_of_interest_mapping.tsv \
  --genbank_dir reference/w303_annotations \
  --output_dir results/gene_verification

# Check if the verification was successful
if [ $? -eq 0 ]; then
    echo "Verification complete. Results are in results/gene_verification/"
    
    # Display summary
    if [ -f "results/gene_verification/gene_verification_summary.txt" ]; then
        echo -e "\nSummary:"
        head -n 10 results/gene_verification/gene_verification_summary.txt
        echo -e "\n(See full report for more details)"
    fi
else
    echo "Error: Verification failed"
    exit 1
fi

echo -e "\nDone!"