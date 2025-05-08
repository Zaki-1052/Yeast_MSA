#!/bin/bash
# Run script for sterol analysis pipeline

# Set up error handling
set -e
set -o pipefail

echo "===== Yeast MSA Sterol Analysis Pipeline ====="
echo "Starting analysis at $(date)"

# Ensure we have necessary directories
mkdir -p sterol_data
mkdir -p results/sterol_analysis/{basic_stats,comparative,pathway,correlation,visualizations}

# Check if sterol data file exists
if [ ! -f "sterol_data/sterol_data_with_sd.csv" ]; then
    echo "Copying sterol data file to sterol_data directory..."
    cp sterol_data_with_sd.csv sterol_data/
fi

# Execute sterol preprocessing
echo "===== Step 1: Sterol Data Preprocessing ====="
echo "Running preprocessing script..."
python scripts/sterols/preprocess_sterols.py
echo "Preprocessing complete!"

# Execute comparative analysis
echo "===== Step 2: Comparative Sterol Analysis ====="
echo "Running comparative analysis script..."
python scripts/sterols/sterol_analysis.py
echo "Comparative analysis complete!"

# Execute pathway analysis
echo "===== Step 3: Ergosterol Pathway Analysis ====="
echo "Running pathway analysis script..."
python scripts/sterols/sterol_pathway.py
echo "Pathway analysis complete!"

# Execute integration analysis
echo "===== Step 4: Sterol-Genomic Integration ====="
echo "Running integration analysis script..."
python scripts/sterols/sterol_integration.py
echo "Integration analysis complete!"

echo "===== All sterol analyses completed successfully! ====="
echo "Analysis finished at $(date)"
echo ""
echo "Results are available in:"
echo "- Basic statistics: results/sterol_analysis/basic_stats/"
echo "- Comparative analysis: results/sterol_analysis/comparative/"
echo "- Pathway analysis: results/sterol_analysis/pathway/"
echo "- Integration analysis: results/sterol_analysis/correlation/"
echo "- Visualizations: results/sterol_analysis/visualizations/"
echo ""
echo "Final integrated report: results/sterol_analysis/correlation/integrated_findings_report.md"