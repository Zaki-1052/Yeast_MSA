#!/bin/bash
# simple_signatures.sh

mkdir -p results/analysis/signatures

# Define top scaffolds for each treatment (hardcoded)
echo "Analyzing WT-specific hotspot: JRIU01000157.1"
grep "JRIU01000157.1" results/functional/WT/variant_details.tsv > results/analysis/signatures/WT_specific_variants.tsv

echo "Analyzing STC-specific hotspot: JRIU01000117.1"
grep "JRIU01000117.1" results/functional/STC/variant_details.tsv > results/analysis/signatures/STC_specific_variants.tsv

echo "Analyzing CAS-specific hotspot: JRIU01000289.1"
grep "JRIU01000289.1" results/functional/CAS/variant_details.tsv > results/analysis/signatures/CAS_specific_variants.tsv

echo "Analyzing WTA-specific hotspot: JRIU01000341.1"
grep "JRIU01000341.1" results/functional/WTA/variant_details.tsv > results/analysis/signatures/WTA_specific_variants.tsv

# Create mutation type summary
echo "=== Treatment-Specific Mutation Patterns ===" > results/analysis/signatures/mutation_patterns.txt

for GROUP in WT STC CAS WTA; do
  VARIANTS_FILE="results/analysis/signatures/${GROUP}_specific_variants.tsv"
  
  if [ -s "$VARIANTS_FILE" ]; then
    echo "" >> results/analysis/signatures/mutation_patterns.txt
    echo "$GROUP specific hotspot mutations:" >> results/analysis/signatures/mutation_patterns.txt
    
    # Count total variants
    TOTAL=$(wc -l < $VARIANTS_FILE)
    echo "Total variants: $TOTAL" >> results/analysis/signatures/mutation_patterns.txt
    
    # Count SNPs and indels
    SNPS=$(awk 'length($3)==1 && length($4)==1' $VARIANTS_FILE | wc -l)
    INDELS=$(awk 'length($3)!=1 || length($4)!=1' $VARIANTS_FILE | wc -l)
    echo "SNPs: $SNPS, Indels: $INDELS" >> results/analysis/signatures/mutation_patterns.txt
    
    # For SNPs, analyze transition vs transversion
    if [ $SNPS -gt 0 ]; then
      echo "Mutation types:" >> results/analysis/signatures/mutation_patterns.txt
      awk 'length($3)==1 && length($4)==1 {print $3">"$4}' $VARIANTS_FILE | sort | uniq -c | sort -nr >> results/analysis/signatures/mutation_patterns.txt
    fi
  fi
done