#!/bin/bash

echo "=== High-Confidence Variant Analysis: Working Solution ==="
mkdir -p results/vcf/final_results

# Part 1: Extract high-confidence variants from direct comparison
echo "Extracting high-confidence variants from direct comparison:"

for GROUP in WT STC CAS WTA; do
  echo "Processing $GROUP..."
  mkdir -p results/vcf/final_results/$GROUP
  
  # Extract variants from each sample's file
  for SAMPLE_FILE in results/vcf/final_analysis/$GROUP/*_specific.vcf; do
    SAMPLE=$(basename "$SAMPLE_FILE" _specific.vcf)
    grep -v "^#" "$SAMPLE_FILE" | awk '{print $1"_"$2"_"$4"_"$5}' > "results/vcf/final_results/$GROUP/${SAMPLE}_variants.txt"
  done
  
  # Find shared variants
  cat results/vcf/final_results/$GROUP/*_variants.txt | sort | uniq -c | sort -nr > results/vcf/final_results/$GROUP/shared_variants.txt
  
  # Count high-confidence variants (in 2+ samples)
  HC_COUNT=$(awk '$1 >= 2' results/vcf/final_results/$GROUP/shared_variants.txt | wc -l)
  echo "  $GROUP high-confidence variants: $HC_COUNT"
  
  # Create detailed list
  echo -e "Count\tChrom\tPos\tRef\tAlt" > results/vcf/final_results/$GROUP/high_conf_variants.tsv
  awk '$1 >= 2 {
    split($2,parts,"_");
    print $1"\t"parts[1]"\t"parts[2]"\t"parts[3]"\t"parts[4]
  }' results/vcf/final_results/$GROUP/shared_variants.txt >> results/vcf/final_results/$GROUP/high_conf_variants.tsv
  
  # Print top 5 variants
  echo "  Top variants:"
  head -5 results/vcf/final_results/$GROUP/shared_variants.txt
done

# Part 2: Create cross-treatment comparison
echo -e "\nCreating cross-treatment comparison:"
# Collect all high-confidence variants
echo -e "Chrom\tPos\tRef\tAlt\tWT\tSTC\tCAS\tWTA" > results/vcf/final_results/treatment_comparison.tsv

# Create a master list of variants
for GROUP in WT STC CAS WTA; do
  awk '$1 >= 2 {print $2}' results/vcf/final_results/$GROUP/shared_variants.txt
done | sort | uniq > results/vcf/final_results/all_variants.txt

# Fill in the comparison matrix
while read VARIANT; do
  IFS='_' read -r CHROM POS REF ALT <<< "$VARIANT"
  LINE="$CHROM\t$POS\t$REF\t$ALT"
  
  # Check presence in each group
  for GROUP in WT STC CAS WTA; do
    if grep -q "$VARIANT" results/vcf/final_results/$GROUP/shared_variants.txt; then
      COUNT=$(grep "$VARIANT" results/vcf/final_results/$GROUP/shared_variants.txt | awk '{print $1}')
      LINE="$LINE\t$COUNT"
    else
      LINE="$LINE\t0"
    fi
  done
  
  echo -e "$LINE" >> results/vcf/final_results/treatment_comparison.tsv
done < results/vcf/final_results/all_variants.txt

# Part 3: Summary statistics
echo -e "\nSummary of high-confidence variants by treatment:"
for GROUP in WT STC CAS WTA; do
  COUNT=$(awk '$1 >= 2' results/vcf/final_results/$GROUP/shared_variants.txt | wc -l)
  echo "$GROUP: $COUNT high-confidence variants"
done

# Count unique variants per treatment
echo -e "\nUnique high-confidence variants by treatment:"
awk 'NR>1 && ($5>=2 && $6==0 && $7==0 && $8==0) {count++} END {print "WT: "count" unique variants"}' results/vcf/final_results/treatment_comparison.tsv
awk 'NR>1 && ($5==0 && $6>=2 && $7==0 && $8==0) {count++} END {print "STC: "count" unique variants"}' results/vcf/final_results/treatment_comparison.tsv
awk 'NR>1 && ($5==0 && $6==0 && $7>=2 && $8==0) {count++} END {print "CAS: "count" unique variants"}' results/vcf/final_results/treatment_comparison.tsv
awk 'NR>1 && ($5==0 && $6==0 && $7==0 && $8>=2) {count++} END {print "WTA: "count" unique variants"}' results/vcf/final_results/treatment_comparison.tsv

echo "Analysis complete! Results in results/vcf/final_results/"