#!/bin/bash
# /Users/zakiralibhai/Documents/GitHub/Yeast_MSA/scripts/satellite_genes/run_satellite_analysis.sh

# This script runs the complete satellite gene analysis pipeline
# It executes all three analysis scripts in sequence:
# 1. satellite_gene_identification.py - Identifies genes in the satellite zone of ERG genes
# 2. satellite_annotation.py - Gathers functional annotations for satellite genes
# 3. satellite_variant_profiling.py - Analyzes variant patterns in satellite genes

# Set base directories
SCRIPT_DIR="$(dirname "$(realpath "$0")")"
PROJECT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"
OUTPUT_DIR="$PROJECT_DIR/results/satellite_genes"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

echo "========================================================"
echo "Satellite Gene Analysis Pipeline"
echo "========================================================"
echo "Project directory: $PROJECT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "========================================================"

# Step 1: Identify genes in the satellite zone of ergosterol pathway genes
echo
echo "Step 1: Identifying satellite genes around ergosterol pathway genes..."
echo "========================================================"
python3 "$SCRIPT_DIR/satellite_gene_identification.py" \
  --gene-mapping "$PROJECT_DIR/reference/gene_mapping_full.tsv" \
  --output-dir "$OUTPUT_DIR" \
  --min-distance 50000 \
  --max-distance 100000

if [ $? -ne 0 ]; then
  echo "Error in satellite gene identification. Exiting."
  exit 1
fi

# Step 2: Annotate satellite genes with functional information
echo
echo "Step 2: Annotating satellite genes with functional information..."
echo "========================================================"
python3 "$SCRIPT_DIR/satellite_annotation.py" \
  --satellite-genes "$OUTPUT_DIR/satellite_genes.tsv" \
  --genbank-dir "$PROJECT_DIR/reference/w303_annotations" \
  --gene-mapping "$PROJECT_DIR/reference/gene_mapping_full.tsv" \
  --output-dir "$OUTPUT_DIR"

if [ $? -ne 0 ]; then
  echo "Error in satellite gene annotation. Exiting."
  exit 1
fi

# Step 3: Analyze variant patterns in satellite genes
echo
echo "Step 3: Analyzing variant patterns in satellite genes..."
echo "========================================================"
python3 "$SCRIPT_DIR/satellite_variant_profiling.py" \
  --satellite-genes "$OUTPUT_DIR/satellite_genes_annotated.tsv" \
  --variant-dir "$PROJECT_DIR/results/gene_variants_expanded" \
  --output-dir "$OUTPUT_DIR"

if [ $? -ne 0 ]; then
  echo "Error in satellite variant profiling. Exiting."
  exit 1
fi

echo
echo "========================================================"
echo "Satellite Gene Analysis Complete!"
echo "Results are available in: $OUTPUT_DIR"
echo "========================================================"

# Create a brief summary report
SUMMARY_FILE="$OUTPUT_DIR/satellite_analysis_summary.txt"
echo "Satellite Gene Analysis Summary" > "$SUMMARY_FILE"
echo "=============================" >> "$SUMMARY_FILE"
echo "Date: $(date)" >> "$SUMMARY_FILE"
echo >> "$SUMMARY_FILE"

echo "1. Satellite Gene Identification:" >> "$SUMMARY_FILE"
if [ -f "$OUTPUT_DIR/satellite_genes.tsv" ]; then
  SATELLITE_COUNT=$(wc -l < "$OUTPUT_DIR/satellite_genes.tsv")
  SATELLITE_COUNT=$((SATELLITE_COUNT - 1))  # Subtract header line
  echo "  - Found $SATELLITE_COUNT genes in the satellite zone (50-100kb from ERG genes)" >> "$SUMMARY_FILE"
else
  echo "  - Satellite gene file not found" >> "$SUMMARY_FILE"
fi

echo >> "$SUMMARY_FILE"
echo "2. Satellite Gene Annotation:" >> "$SUMMARY_FILE"
if [ -f "$OUTPUT_DIR/satellite_annotation_summary.txt" ]; then
  # Extract key metrics from annotation summary
  echo "  - Annotation metrics:" >> "$SUMMARY_FILE"
  grep -A 10 "Functional Categories:" "$OUTPUT_DIR/satellite_annotation_summary.txt" | grep -v "^For" >> "$SUMMARY_FILE"
else
  echo "  - Satellite gene annotation summary not found" >> "$SUMMARY_FILE"
fi

echo >> "$SUMMARY_FILE"
echo "3. Satellite Variant Profiling:" >> "$SUMMARY_FILE"
if [ -f "$OUTPUT_DIR/satellite_variants.tsv" ]; then
  VARIANT_COUNT=$(wc -l < "$OUTPUT_DIR/satellite_variants.tsv")
  VARIANT_COUNT=$((VARIANT_COUNT - 1))  # Subtract header line
  echo "  - Found $VARIANT_COUNT variants in satellite genes" >> "$SUMMARY_FILE"
  
  # Extract key metrics from variant profiling report
  if [ -f "$OUTPUT_DIR/satellite_variant_profiling_report.txt" ]; then
    echo "  - Variant distribution:" >> "$SUMMARY_FILE"
    grep -A 5 "Variants by Impact:" "$OUTPUT_DIR/satellite_variant_profiling_report.txt" | grep -v "^[0-9]" >> "$SUMMARY_FILE"
    
    echo "  - Adaptation-specific patterns:" >> "$SUMMARY_FILE"
    grep -A 2 "Temperature adaptation" "$OUTPUT_DIR/satellite_variant_profiling_report.txt" | grep -v "^Impact" >> "$SUMMARY_FILE"
  fi
else
  echo "  - Satellite variant file not found" >> "$SUMMARY_FILE"
fi

echo >> "$SUMMARY_FILE"
echo "For detailed results, please review the individual analysis files in $OUTPUT_DIR" >> "$SUMMARY_FILE"

echo "Summary report created: $SUMMARY_FILE"

# Make the script executable
chmod +x "$SCRIPT_DIR/satellite_gene_identification.py"
chmod +x "$SCRIPT_DIR/satellite_annotation.py"
chmod +x "$SCRIPT_DIR/satellite_variant_profiling.py"

exit 0