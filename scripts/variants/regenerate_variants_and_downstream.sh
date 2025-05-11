#!/bin/bash
# regenerate_variants_and_downstream.sh
#
# This script regenerates gene variant files using the fixed script and then
# reruns all downstream analyses that depended on the flawed data.
#
# It should be run from the project root directory.

set -e  # Exit on any error

# Activate the virtual environment if it exists
if [ -d "venv" ]; then
  echo "Activating virtual environment..."
  source venv/bin/activate
fi

echo "======================================================================"
echo "Regenerating Gene Variants and Downstream Analyses"
echo "======================================================================"

# Define important directories
PROJECT_DIR=$(pwd)
SCRIPTS_DIR="$PROJECT_DIR/scripts"
RESULTS_DIR="$PROJECT_DIR/results"
VCF_DIR="$PROJECT_DIR/vcf/annotated"
GENE_MAPPING="$PROJECT_DIR/reference/genes_of_interest_mapping.tsv"

# Step 1: Backup the existing gene_variants directories
echo "Backing up existing gene_variants directories..."
timestamp=$(date +%Y%m%d_%H%M%S)
mkdir -p "$RESULTS_DIR/backups"
backup_dir="$RESULTS_DIR/backups/$timestamp"
mkdir -p "$backup_dir"

# Only backup if the directories exist
if [ -d "$RESULTS_DIR/gene_variants" ]; then
  echo "Backing up $RESULTS_DIR/gene_variants"
  cp -r "$RESULTS_DIR/gene_variants" "$backup_dir/gene_variants"
fi

if [ -d "$RESULTS_DIR/gene_variants_expanded" ]; then
  echo "Backing up $RESULTS_DIR/gene_variants_expanded"
  cp -r "$RESULTS_DIR/gene_variants_expanded" "$backup_dir/gene_variants_expanded"
fi

echo "Backups created in $backup_dir"

# Step 2: Regenerate the gene variants files
echo "Regenerating gene variants files..."

# a) Standard gene variants (1000bp upstream)
echo "Generating standard gene variants (1000bp upstream)..."
rm -rf "$RESULTS_DIR/gene_variants"  # Remove existing directory
mkdir -p "$RESULTS_DIR/gene_variants"
python3 "$SCRIPTS_DIR/variants/extract_gene_variants_fixed.py" \
  --vcf_dir "$VCF_DIR" \
  --gene_mapping "$GENE_MAPPING" \
  --output_dir "$RESULTS_DIR/gene_variants" \
  --exclude_controls  # Optional: exclude variants present in control samples from treatment samples

# b) Expanded gene variants (5000bp upstream)
echo "Generating expanded gene variants (5000bp upstream)..."
rm -rf "$RESULTS_DIR/gene_variants_expanded"  # Remove existing directory
mkdir -p "$RESULTS_DIR/gene_variants_expanded"
python3 "$SCRIPTS_DIR/variants/extract_gene_variants_fixed.py" \
  --vcf_dir "$VCF_DIR" \
  --gene_mapping "$GENE_MAPPING" \
  --output_dir "$RESULTS_DIR/gene_variants_expanded" \
  --upstream_distance 5000 \
  --exclude_controls  # Optional: exclude variants present in control samples from treatment samples

# Step 3: Regenerate all downstream analyses
echo "Regenerating downstream analyses..."

# a) Regenerate regulatory analysis reports
echo "Regenerating regulatory analysis..."
if [ -d "$SCRIPTS_DIR/regulatory_analysis" ]; then
  # Clear existing results
  rm -rf "$RESULTS_DIR/regulatory_analysis"
  mkdir -p "$RESULTS_DIR/regulatory_analysis/accurate_mapping"
  
  # Run the accurate regulatory analysis
  python3 "$SCRIPTS_DIR/regulatory_analysis/accurate_regulatory_analysis.py" \
    --variants "$RESULTS_DIR/gene_variants_expanded/all_gene_variants.tsv" \
    --output-dir "$RESULTS_DIR/regulatory_analysis/accurate_mapping"

  # Run the regulatory motif analysis
  mkdir -p "$RESULTS_DIR/regulatory_analysis/motif_analysis"
  python3 "$SCRIPTS_DIR/regulatory_analysis/regulatory_motif_analysis.py" \
    --variants "$RESULTS_DIR/gene_variants_expanded/all_gene_variants.tsv" \
    --output-dir "$RESULTS_DIR/regulatory_analysis/motif_analysis"
    
  echo "Regulatory analysis regenerated"
fi

# b) Regenerate satellite gene analysis
echo "Regenerating satellite gene analysis..."
if [ -f "$SCRIPTS_DIR/satellite_genes/run_satellite_analysis.sh" ]; then
  bash "$SCRIPTS_DIR/satellite_genes/run_satellite_analysis.sh"
  echo "Satellite gene analysis regenerated"
else
  echo "Satellite gene analysis script not found, skipping"
fi

# c) Regenerate OSH gene analysis
echo "Regenerating OSH gene analysis..."
if [ -f "$SCRIPTS_DIR/osh_analysis/run_osh_analysis.sh" ]; then
  bash "$SCRIPTS_DIR/osh_analysis/run_osh_analysis.sh"
  echo "OSH gene analysis regenerated"
else
  echo "OSH gene analysis script not found, skipping"
fi

# d) Regenerate promoter analysis
echo "Regenerating promoter analysis..."
if [ -f "$SCRIPTS_DIR/functional_impact/run_promoter_analysis.sh" ]; then
  bash "$SCRIPTS_DIR/functional_impact/run_promoter_analysis.sh"
  echo "Promoter analysis regenerated"
else
  echo "Promoter analysis script not found, skipping"
fi

# e) Regenerate the ergosterol variant report
echo "Regenerating ergosterol variant report..."
if [ -f "$SCRIPTS_DIR/utils/generate_ergosterol_variant_report.py" ]; then
  # This script doesn't accept command line arguments, so we need to edit it directly
  # to point to the correct paths before running
  TMP_SCRIPT="/tmp/generate_ergosterol_variant_report_tmp.py"
  cp "$SCRIPTS_DIR/utils/generate_ergosterol_variant_report.py" "$TMP_SCRIPT"
  sed -i '' "s|output_path = \".*\"|output_path = \"$RESULTS_DIR/reports/ergosterol_variant_analysis.html\"|" "$TMP_SCRIPT"
  python3 "$TMP_SCRIPT"
  rm "$TMP_SCRIPT"
  echo "Ergosterol variant report regenerated"
else
  echo "Ergosterol variant report script not found, skipping"
fi

# f) Check if functional impact report script exists and run it if it does
echo "Checking for functional impact report script..."
if [ -f "$SCRIPTS_DIR/utils/generate_functional_impact_report.py" ]; then
  # Use Task to analyze how to run this script
  echo "Generating functional impact report..."
  python3 "$SCRIPTS_DIR/utils/generate_functional_impact_report.py"
  echo "Functional impact report regenerated"
else
  echo "Functional impact report script not found, skipping"
fi

echo "======================================================================"
echo "Regeneration Complete!"
echo "======================================================================"
echo "All gene variant files and downstream analyses have been regenerated."
echo "Original data was backed up to $backup_dir"
echo ""
echo "IMPORTANT: Please review the new results to ensure they are correct."
echo "Particularly look for variant count differences between treatments,"
echo "which should now show biological patterns rather than identical variants."
echo "======================================================================"