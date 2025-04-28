#!/bin/bash

# debug_genbank_scaffolds.sh - Run the GenBank scaffold debugging script

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
echo "Debugging GenBank scaffold naming issues..."
python scripts/variants/debug_genbank_scaffolds.py \
  --genbank_dir reference/w303_annotations \
  --gene_mapping reference/genes_of_interest_mapping.tsv \
  --output results/debugging/genbank_scaffold_debug.txt

# Check if the debugging was successful
if [ $? -eq 0 ]; then
    echo "Debugging complete. Results are in results/debugging/genbank_scaffold_debug.txt"
    
    # Display a summary
    echo -e "\nSummary from debugging:"
    echo "========================="
    echo -e "\nProbable issue:"
    grep -A 10 "Diagnosis and Fix Recommendation" results/debugging/genbank_scaffold_debug.txt
    
    echo -e "\nSuggested fix:"
    echo "1. Modify the verify_gene_coordinates.py script to extract scaffold IDs correctly"
    echo "2. Update the parse_genbank_files function to use the correct source for scaffold names"
    echo "3. Consider extracting scaffold info from filenames if not available in record fields"
else
    echo "Error: Debugging failed"
    exit 1
fi

echo -e "\nDone!"