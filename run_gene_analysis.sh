#!/bin/bash

echo "Running Gene-Specific Analysis Scripts"
echo "====================================="

# Create output directory if it doesn't exist
mkdir -p analysis/genes_of_interest

# Function to run a script and time it
run_script() {
    script_name=$1
    echo ""
    echo "Running $script_name..."
    start_time=$(date +%s)
    python3 scripts/gene_analysis/$script_name
    end_time=$(date +%s)
    runtime=$((end_time - start_time))
    echo "$script_name completed in $runtime seconds"
}

# Run all gene analysis scripts
echo "Starting analysis pipeline..."

# Step 1: Mutation Spectrum Analysis
run_script "mutation_spectrum_analysis.py"

# Step 2: Genomic Context Analysis
run_script "genomic_context_analysis.py"

# Step 3: Mutational Signature Analysis
run_script "mutational_signature_analysis.py"

# Step 4: Population Spectrum Analysis
run_script "population_spectrum_analysis.py"

# Step 5: Regional Enrichment Analysis
run_script "regional_enrichment_analysis.py"

# Step 6: Scaffold Distribution Analysis
run_script "scaffold_distribution_analysis.py"

# Step 7: Statistical Pattern Analysis
run_script "statistical_pattern_analysis.py"

# Step 8: Variation Analysis
run_script "variation.py"

# Step 9: TC Visualization
run_script "TC_Visualization.py"

# Final Step: Check outputs
echo ""
echo "Analysis complete! Checking outputs..."
find analysis/genes_of_interest -type f | wc -l
echo "Output files generated in analysis/genes_of_interest/"