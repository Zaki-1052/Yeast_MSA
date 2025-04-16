#!/bin/bash

echo "=== W303 Yeast Variant Analysis Data Exploration ==="
echo "Running data gathering script on $(date)"
echo

# 1. Exploring directory structure
echo "=== DIRECTORY STRUCTURE ==="
echo "Checking high-confidence variant files location:"
find results/merged/analysis -name "*highconf.vcf.gz" | sort
echo

echo "Checking summary files location:"
find results/merged/summary -type f | sort
echo

# 2. Examine VCF structure and count variants
echo "=== VCF FILE ANALYSIS ==="
for GROUP in WT STC CAS WTA; do
  echo "--- $GROUP Treatment Group ---"
  
  # Check if the high-confidence VCF file exists
  HC_FILE="results/merged/analysis/${GROUP}_highconf.vcf.gz"
  if [ -f "$HC_FILE" ]; then
    echo "File exists: $HC_FILE"
    
    # Count total variants
    echo "Total high-confidence variants:"
    bcftools view -H "$HC_FILE" | wc -l
    
    # List samples in the VCF
    echo "Samples in VCF:"
    bcftools query -l "$HC_FILE"
    
    # Get file size
    echo "File size:"
    ls -lh "$HC_FILE"
    
    # Show header structure
    echo "VCF header structure (first 5 lines):"
    bcftools view -h "$HC_FILE" | head -5
    
    # Sample of variants (first 5)
    echo "Sample variants (first 5):"
    bcftools view -H "$HC_FILE" | head -5
    
    # Count mutation types (simple version)
    echo "Initial mutation type counts:"
    bcftools query -f '%REF>%ALT\n' "$HC_FILE" | sort | uniq -c | sort -nr
  else
    echo "File not found: $HC_FILE"
  fi
  echo
done

# 3. Check summary files content
echo "=== SUMMARY FILES ANALYSIS ==="
echo "Treatment comparison file (first 10 lines):"
if [ -f "results/merged/summary/treatment_comparison.tsv" ]; then
  head -10 "results/merged/summary/treatment_comparison.tsv"
else
  echo "Treatment comparison file not found"
fi
echo

echo "Unique variants summary:"
if [ -f "results/merged/summary/unique_variants.txt" ]; then
  cat "results/merged/summary/unique_variants.txt"
else
  echo "Unique variants file not found"
fi
echo

echo "Analysis report summary:"
if [ -f "results/merged/summary/analysis_report.txt" ]; then
  cat "results/merged/summary/analysis_report.txt"
else
  echo "Analysis report not found"
fi
echo

# 4. Create list of necessary files for further analysis
echo "=== PREPARING DATA FOR MUTATION SPECTRUM ANALYSIS ==="
echo "Creating list of necessary VCF files:"

ANALYSIS_DIR="mutation_spectrum_analysis"
mkdir -p "$ANALYSIS_DIR"

for GROUP in WT STC CAS WTA; do
  HC_FILE="results/merged/analysis/${GROUP}_highconf.vcf.gz"
  if [ -f "$HC_FILE" ]; then
    echo "$HC_FILE" >> "$ANALYSIS_DIR/highconf_vcf_files.txt"
    
    # Extract mutation types for quick check
    echo "Extracting mutation data for $GROUP..."
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$HC_FILE" > "$ANALYSIS_DIR/${GROUP}_mutations.txt"
    
    # Count how many of each type
    echo "Mutation types in $GROUP:"
    cut -f 3,4 "$ANALYSIS_DIR/${GROUP}_mutations.txt" | sort | uniq -c | sort -nr > "$ANALYSIS_DIR/${GROUP}_mutation_counts.txt"
    cat "$ANALYSIS_DIR/${GROUP}_mutation_counts.txt"
  else
    echo "Warning: $HC_FILE not found, cannot extract mutation data for $GROUP"
  fi
done

echo
echo "Data exploration complete. Results saved to $ANALYSIS_DIR directory."