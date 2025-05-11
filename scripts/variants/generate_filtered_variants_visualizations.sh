#!/bin/bash

# generate_filtered_variants_visualizations.sh - Create visualizations for treatment-specific variants

# Check if the required directories and files exist
if [ ! -f "results/filtered_scaffold_variants/treatment_specific_scaffold_variants.tsv" ]; then
    echo "Error: File results/filtered_scaffold_variants/treatment_specific_scaffold_variants.tsv does not exist"
    echo "Please run extract_filtered_scaffold_variants.sh first to generate the filtered variants data"
    exit 1
fi

if [ ! -f "reference/genes_of_interest_mapping.tsv" ]; then
    echo "Error: File reference/genes_of_interest_mapping.tsv does not exist"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p results/filtered_scaffold_variants/visualizations

# Run the Python script
echo "Generating visualizations for treatment-specific variants near ERG genes..."
source venv/bin/activate
python3 scripts/variants/generate_filtered_variants_visualizations.py \
  --input_file results/filtered_scaffold_variants/treatment_specific_scaffold_variants.tsv \
  --gene_mapping reference/genes_of_interest_mapping.tsv \
  --output_dir results/filtered_scaffold_variants/visualizations

# Check if the visualization generation was successful
if [ $? -eq 0 ]; then
    echo "Visualization complete. Results are in results/filtered_scaffold_variants/visualizations/"
    echo "HTML report available at: results/filtered_scaffold_variants/visualizations/filtered_variants_report.html"
else
    echo "Error: Visualization generation failed"
    exit 1
fi

echo -e "\nDone!"