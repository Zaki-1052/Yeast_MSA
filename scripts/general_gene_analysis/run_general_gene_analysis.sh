#!/bin/bash

echo "Running General Gene Analysis Scripts with Updated Gene Mapping"
echo "=============================================================="

# Create output directories if they don't exist
mkdir -p analysis/general_gene_analysis
mkdir -p analysis/general_gene_analysis/mutation_spectrum_results
mkdir -p analysis/general_gene_analysis/gene_mutation_spectrum_results
mkdir -p analysis/general_gene_analysis/population_structure_results
mkdir -p analysis/general_gene_analysis/gene_population_structure_results
mkdir -p analysis/general_gene_analysis/genomic_context_results
mkdir -p analysis/general_gene_analysis/gene_genomic_context_results
mkdir -p analysis/general_gene_analysis/statistical_pattern_results
mkdir -p analysis/general_gene_analysis/gene_statistical_pattern_results
mkdir -p analysis/general_gene_analysis/regional_enrichment_results
mkdir -p analysis/general_gene_analysis/gene_regional_enrichment_results
mkdir -p analysis/general_gene_analysis/treatment_control_analysis
mkdir -p analysis/general_gene_analysis/gene_treatment_control_analysis

# Function to run a script and time it
run_script() {
    script_name=$1
    echo ""
    echo "Running $script_name..."
    start_time=$(date +%s)
    # Use environment variable for OUTPUT_DIR to ensure all scripts use the same output location
    OUTPUT_DIR="analysis/general_gene_analysis" python3 scripts/general_gene_analysis/$script_name
    end_time=$(date +%s)
    runtime=$((end_time - start_time))
    echo "$script_name completed in $runtime seconds"
}

# Run all gene analysis scripts
echo "Starting analysis pipeline..."

# Step 1: First run the gene mapping full generator to ensure latest gene mapping
echo "Generating comprehensive gene mapping..."
run_script "generate_gene_mapping_full.py"

# Step 2: Mutation Spectrum Analysis
run_script "mutation_spectrum_analysis.py"

# Step 3: Variation Analysis
run_script "variation.py"

# Step 4: TC Visualization
run_script "TC_Visualization.py"

# Final Step: Check outputs
echo ""
echo "Analysis complete! Checking outputs..."
echo "Main output directory:"
find analysis/general_gene_analysis -type f | wc -l
echo "Output files generated in analysis/general_gene_analysis/"

echo ""
echo "All related output directories:"
echo "Mutation spectrum results:"
find analysis/general_gene_analysis/mutation_spectrum_results -type f | wc -l
echo "Gene mutation spectrum results:"
find analysis/general_gene_analysis/gene_mutation_spectrum_results -type f | wc -l
echo "Population structure results:"
find analysis/general_gene_analysis/population_structure_results -type f | wc -l
echo "Gene population structure results:"
find analysis/general_gene_analysis/gene_population_structure_results -type f | wc -l
echo "Genomic context results:"
find analysis/general_gene_analysis/genomic_context_results -type f | wc -l
echo "Gene genomic context results:" 
find analysis/general_gene_analysis/gene_genomic_context_results -type f | wc -l
echo "Statistical pattern results:"
find analysis/general_gene_analysis/statistical_pattern_results -type f | wc -l
echo "Gene statistical pattern results:"
find analysis/general_gene_analysis/gene_statistical_pattern_results -type f | wc -l
echo "Regional enrichment results:"
find analysis/general_gene_analysis/regional_enrichment_results -type f | wc -l
echo "Gene regional enrichment results:"
find analysis/general_gene_analysis/gene_regional_enrichment_results -type f | wc -l
echo "Treatment control analysis results:"
find analysis/general_gene_analysis/treatment_control_analysis -type f | wc -l
echo "Gene treatment control analysis results:"
find analysis/general_gene_analysis/gene_treatment_control_analysis -type f | wc -l