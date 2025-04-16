#!/bin/bash

# File: scripts/annotation/28_search_gene_variants_directly.sh
# Purpose: Search directly for gene variants using grep

echo "=== Direct Gene Variant Search ==="
echo "Date: $(date)"
echo ""

# Define directories
OUTPUT_DIR="annotation/gene_variants_grep"
VCF_DIR="annotation/vcf_ready"

# Create output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/by_gene"

# Define genes to search for
GENES=(
  "YHR190W" "W303_0BY00190" "W3030BY00190" "ERG9"
  "YGR175C" "W303_0AJ00440" "W3030AJ00440" "ERG1"
  "YHR072W" "W303_0CB00150" "W3030CB00150" "ERG7"
  "YHR007C" "W303_0EI00110" "W3030EI00110" "ERG11"
  "YNL280C" "W303_0O00140" "W3030O00140" "ERG24"
  "YGR060W" "W303_0W00270" "W3030W00270" "ERG25"
  "YML008C" "W303_0S00270" "W3030S00270" "ERG6"
  "YMR202W" "W303_0AD00260" "W3030AD00260" "ERG2"
  "YLR056W" "W303_0E01010" "W3030E01010" "ERG3"
  "YMR015C" "W303_0S00490" "W3030S00490" "ERG5"
  "YGL012W" "W303_0Y00390" "W3030Y00390" "ERG4"
)

# Find VCF files
VCF_FILES=$(find "$VCF_DIR" -name "*.sorted.vcf.gz")

# For each SGD gene
for ((i=0; i<${#GENES[@]}; i+=4)); do
  SGD_GENE="${GENES[$i]}"
  W303_ID="${GENES[$i+1]}"
  GENE_ID="${GENES[$i+2]}"
  ERG_NAME="${GENES[$i+3]}"
  
  echo "Searching for $SGD_GENE ($ERG_NAME)..."
  
  # Output file
  OUTPUT_FILE="$OUTPUT_DIR/by_gene/${SGD_GENE}_variants.txt"
  echo "# Variants for $SGD_GENE ($ERG_NAME)" > "$OUTPUT_FILE"
  echo "# W303 ID: $W303_ID" >> "$OUTPUT_FILE"
  echo "# Gene ID: $GENE_ID" >> "$OUTPUT_FILE"
  echo "# Date: $(date)" >> "$OUTPUT_FILE"
  echo "#" >> "$OUTPUT_FILE"
  echo "# Sample Chromosome Position Ref Alt Info" >> "$OUTPUT_FILE"
  
  # Search for all identifiers in all VCF files
  for VCF_FILE in $VCF_FILES; do
    SAMPLE=$(basename "$VCF_FILE" .sorted.vcf.gz)
    echo "  Processing $SAMPLE..."
    
    # Search for all possible identifiers
    {
      zcat "$VCF_FILE" | grep -v "^#" | grep -e "$SGD_GENE" -e "$W303_ID" -e "$GENE_ID" -e "$ERG_NAME" || true
    } | while read -r LINE; do
      # Extract fields
      CHROM=$(echo "$LINE" | awk '{print $1}')
      POS=$(echo "$LINE" | awk '{print $2}')
      REF=$(echo "$LINE" | awk '{print $4}')
      ALT=$(echo "$LINE" | awk '{print $5}')
      INFO=$(echo "$LINE" | awk '{print $8}')
      
      # Add to output
      echo "$SAMPLE $CHROM $POS $REF $ALT $INFO" >> "$OUTPUT_FILE"
    done
  done
  
  # Count variants
  COUNT=$(wc -l < "$OUTPUT_FILE")
  COUNT=$((COUNT - 6))  # Subtract header lines
  if [ $COUNT -gt 0 ]; then
    echo "  Found $COUNT variants for $SGD_GENE"
  else
    echo "  No variants found for $SGD_GENE"
  fi
  echo ""
done

echo "=== Direct Gene Search Complete ==="
echo "Results saved to $OUTPUT_DIR"
