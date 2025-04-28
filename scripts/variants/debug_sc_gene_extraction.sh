#!/bin/bash

# debug_sc_gene_extraction.sh - Debug SC gene ID extraction from GenBank files

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
mkdir -p results/debugging

# Run the Python script
echo "Debugging SC gene ID extraction from GenBank files..."
python scripts/debug_sc_gene_extraction.py \
  --genbank_dir reference/w303_annotations \
  --gene_mapping reference/genes_of_interest_mapping.tsv \
  --output results/debugging/sc_gene_extraction_debug.txt

# Check if the debugging was successful
if [ $? -eq 0 ]; then
    echo "Debugging complete. Results are in results/debugging/sc_gene_extraction_debug.txt"
    
    # Display a summary
    echo -e "\nSummary from debugging:"
    echo "========================="
    
    # Extract and display found/missing gene counts
    echo -e "\nGene matching status:"
    grep -A 1 "Total SC gene IDs found:" results/debugging/sc_gene_extraction_debug.txt
    
    # Extract and display inference patterns
    echo -e "\nCommon inference patterns:"
    grep -A 3 "Common inference patterns:" results/debugging/sc_gene_extraction_debug.txt | tail -n 3
    
    # Extract and display recommendations
    echo -e "\nRecommended fix:"
    grep -A 5 "Recommended fixes for verify_gene_coordinates.py:" results/debugging/sc_gene_extraction_debug.txt | tail -n 5
else
    echo "Error: Debugging failed"
    exit 1
fi

echo -e "\nDone!"