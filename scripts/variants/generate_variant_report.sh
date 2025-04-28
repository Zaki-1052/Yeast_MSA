#!/bin/bash

# generate_variant_report.sh - Create a visual report of the scaffold variant analysis

# Check if the required directories exist
if [ ! -d "results/scaffold_variants" ]; then
    echo "Error: Directory results/scaffold_variants does not exist"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p results/reports

# Run the Python script
echo "Generating variant analysis report..."
python scripts/generate_variant_report.py \
  --data_dir results/scaffold_variants \
  --output results/reports/ergosterol_variant_analysis.html

# Check if the report generation was successful
if [ $? -eq 0 ]; then
    echo "Report generated successfully: results/reports/ergosterol_variant_analysis.html"
    
    # Try to open the report in the default browser
    if [ "$(uname)" == "Darwin" ]; then
        # macOS
        open results/reports/ergosterol_variant_analysis.html
    elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
        # Linux
        if command -v xdg-open > /dev/null; then
            xdg-open results/reports/ergosterol_variant_analysis.html
        fi
    fi
else
    echo "Error: Report generation failed"
    exit 1
fi

echo -e "\nDone!"