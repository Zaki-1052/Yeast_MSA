#!/bin/bash

echo "Testing variation.py with chromosome ID mapping fix"
echo "=================================================="

# Clear any existing log
> analysis.log

# Run the variation.py script with detailed output
python3 scripts/gene_analysis/variation.py

# Show gene mapping section of the log
echo ""
echo "Showing gene mapping output from the log:"
echo "----------------------------------------"
grep -A 20 "Loaded .* genes from reference/gene_mapping.tsv" analysis.log || cat analysis.log

# Show mapping statistics
echo ""
echo "Showing mapping statistics:"
echo "-------------------------"
grep "Mapped .* out of .* variants to genes" analysis.log

echo ""
echo "Test complete"