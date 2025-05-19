#!/bin/bash
# /Users/zakiralibhai/Documents/GitHub/Yeast_MSA/scripts/specific_regulatory_analysis/run_specific_analysis.sh

# This script runs the complete specific regulatory analysis pipeline
# which focuses on gene variants (close to genes) rather than scaffold variants

# Ensure we're in the right directory
cd /Users/zakiralibhai/Documents/GitHub/Yeast_MSA

# Activate virtual environment
source venv/bin/activate

# Step 1: Fix variant annotations
echo "===== Step 1: Fixing variant annotations ====="
python3 scripts/specific_regulatory_analysis/fix_variant_annotations.py

# Step 2: Run variant regulatory mapping
echo "===== Step 2: Running variant regulatory mapping ====="
python3 scripts/specific_regulatory_analysis/variant_regulatory_mapping.py

# Step 3: Run treatment regulatory analysis
echo "===== Step 3: Running treatment regulatory analysis ====="
python3 scripts/specific_regulatory_analysis/treatment_regulatory_analysis.py

echo "===== Specific Regulatory Analysis Complete ====="
echo "Results are available in results/specific_regulatory_analysis/"