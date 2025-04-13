#!/bin/bash

echo "===== Starting Yeast Variant Analysis Pipeline ====="
echo "Running pipeline with biological groupings"
echo "Date: $(date)"
echo ""

# Create directories for scripts if they don't exist
mkdir -p scripts/pipeline

# Make all scripts executable
chmod +x scripts/pipeline/*.sh

# Run each step of the pipeline
echo "Step 1: Setting up directory structure"
bash scripts/pipeline/01_setup.sh

echo "Step 2: Filtering and normalizing VCF files"
bash scripts/pipeline/02_filter_normalize.sh

echo "Step 3: Fixing contig definitions"
bash scripts/pipeline/03_fix_contigs.sh

echo "Step 4: Merging VCF files"
bash scripts/pipeline/04_merge_vcfs.sh

echo "Step 5: Testing sample selection"
bash scripts/pipeline/05_test_samples.sh

echo "Step 6: Performing treatment vs control comparisons"
bash scripts/pipeline/06_compare_groups.sh

echo "Step 7: Creating cross-treatment comparison"
bash scripts/pipeline/07_cross_treatment.sh

echo "Step 8: Identifying treatment-unique variants"
bash scripts/pipeline/08_identify_unique.sh

echo "Step 9: Running direct VCF comparisons"
bash scripts/pipeline/09_direct_comparison.sh

echo "Step 10: Analyzing within-group consistency"
bash scripts/pipeline/10_consistency_analysis.sh

echo "Step 11: Generating final summary report"
bash scripts/pipeline/11_summary_report.sh

echo ""
echo "===== Yeast Variant Analysis Pipeline Complete ====="
echo "Results available in results/merged/ directory"
echo "Summary report: results/merged/summary/analysis_report.txt"