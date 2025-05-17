#!/bin/bash
# /Users/zakiralibhai/Documents/GitHub/Yeast_MSA/scripts/osh_analysis/run_osh_analysis.sh

# This script runs the complete OSH gene family analysis
# It executes all three analysis scripts in sequence:
# 1. analyze_osh_genes.py - Maps OSH family genes in the reference genome
# 2. osh_variants.py - Analyzes variants in and around OSH genes
# 3. osh_erg_distance.py - Calculates genomic distances between OSH and ERG genes

# Set base directories
SCRIPT_DIR="$(dirname "$(realpath "$0")")"
PROJECT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"
OUTPUT_DIR="$PROJECT_DIR/results/osh_analysis"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

echo "========================================================"
echo "OSH Gene Family Analysis"
echo "========================================================"
echo "Project directory: $PROJECT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "========================================================"

# Step 1: Map OSH genes in the reference genome
echo
echo "Step 1: Mapping OSH family genes in the reference genome..."
echo "========================================================"
python3 "$SCRIPT_DIR/analyze_osh_genes.py" \
  --ref-genome "$PROJECT_DIR/reference/w303_chromosomal.fasta" \
  --gene-mapping "$PROJECT_DIR/reference/gene_mapping_full.tsv" \
  --output-dir "$OUTPUT_DIR"

if [ $? -ne 0 ]; then
  echo "Error in OSH gene mapping analysis. Exiting."
  exit 1
fi

# Step 2: Analyze variants in and around OSH genes
echo
echo "Step 2: Analyzing variants in and around OSH genes..."
echo "========================================================"
python3 "$SCRIPT_DIR/osh_variants.py" \
  --gene-mapping "$PROJECT_DIR/reference/gene_mapping_full.tsv" \
  --variant-dir "$PROJECT_DIR/results/gene_variants" \
  --output-dir "$OUTPUT_DIR" \
  --distance-threshold 25000

if [ $? -ne 0 ]; then
  echo "Error in OSH variant analysis. Exiting."
  exit 1
fi

# Step 3: Calculate distances between OSH and ERG genes
echo
echo "Step 3: Calculating genomic distances between OSH and ERG genes..."
echo "========================================================"
python3 "$SCRIPT_DIR/osh_erg_distance.py" \
  --gene-mapping "$PROJECT_DIR/reference/gene_mapping_full.tsv" \
  --genome-file "$PROJECT_DIR/reference/w303_chromosomal.fasta" \
  --variant-dir "$PROJECT_DIR/results/gene_variants" \
  --output-dir "$OUTPUT_DIR"

if [ $? -ne 0 ]; then
  echo "Error in OSH-ERG distance analysis. Exiting."
  exit 1
fi

echo
echo "========================================================"
echo "OSH Gene Family Analysis Complete!"
echo "Results are available in: $OUTPUT_DIR"
echo "========================================================"

# Create a brief summary report
SUMMARY_FILE="$OUTPUT_DIR/osh_analysis_summary.txt"
echo "OSH Gene Family Analysis Summary" > "$SUMMARY_FILE"
echo "===============================" >> "$SUMMARY_FILE"
echo "Date: $(date)" >> "$SUMMARY_FILE"
echo >> "$SUMMARY_FILE"

echo "OSH Gene Mapping:" >> "$SUMMARY_FILE"
if [ -f "$OUTPUT_DIR/osh_gene_summary.tsv" ]; then
  OSH_COUNT=$(wc -l < "$OUTPUT_DIR/osh_gene_summary.tsv")
  OSH_COUNT=$((OSH_COUNT - 1))  # Subtract header line
  echo "  Found $OSH_COUNT OSH family genes" >> "$SUMMARY_FILE"
else
  echo "  OSH gene mapping file not found" >> "$SUMMARY_FILE"
fi

echo >> "$SUMMARY_FILE"
echo "OSH Variant Analysis:" >> "$SUMMARY_FILE"
if [ -f "$OUTPUT_DIR/osh_variants.tsv" ]; then
  VARIANT_COUNT=$(wc -l < "$OUTPUT_DIR/osh_variants.tsv")
  VARIANT_COUNT=$((VARIANT_COUNT - 1))  # Subtract header line
  echo "  Found $VARIANT_COUNT variants near OSH genes" >> "$SUMMARY_FILE"
else
  echo "  OSH variants file not found" >> "$SUMMARY_FILE"
fi

echo >> "$SUMMARY_FILE"
echo "OSH-ERG Distance Analysis:" >> "$SUMMARY_FILE"
if [ -f "$OUTPUT_DIR/osh_erg_distances.tsv" ]; then
  DISTANCE_COUNT=$(wc -l < "$OUTPUT_DIR/osh_erg_distances.tsv")
  DISTANCE_COUNT=$((DISTANCE_COUNT - 1))  # Subtract header line
  echo "  Analyzed $DISTANCE_COUNT OSH-ERG gene relationships" >> "$SUMMARY_FILE"
else
  echo "  OSH-ERG distance file not found" >> "$SUMMARY_FILE"
fi

echo >> "$SUMMARY_FILE"
echo "For detailed results, please review the individual analysis files in $OUTPUT_DIR" >> "$SUMMARY_FILE"

echo "Summary report created: $SUMMARY_FILE"
exit 0