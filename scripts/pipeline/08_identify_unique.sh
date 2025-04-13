#!/bin/bash

echo "=== Identifying treatment-unique variants ==="

echo -e "Treatment\tUnique_Variants" > "results/merged/summary/unique_variants.txt"

# Count unique variants for each treatment (corrected for WT-37)
awk 'NR>1 && ($5==1 && $6==0 && $7==0 && $8==0) {count++} \
     END {print "WT-37\t"count}' \
     "results/merged/summary/treatment_comparison.tsv" >> \
     "results/merged/summary/unique_variants.txt"
     
awk 'NR>1 && ($5==0 && $6==1 && $7==0 && $8==0) {count++} \
     END {print "STC\t"count}' \
     "results/merged/summary/treatment_comparison.tsv" >> \
     "results/merged/summary/unique_variants.txt"
     
awk 'NR>1 && ($5==0 && $6==0 && $7==1 && $8==0) {count++} \
     END {print "CAS\t"count}' \
     "results/merged/summary/treatment_comparison.tsv" >> \
     "results/merged/summary/unique_variants.txt"
     
awk 'NR>1 && ($5==0 && $6==0 && $7==0 && $8==1) {count++} \
     END {print "WTA\t"count}' \
     "results/merged/summary/treatment_comparison.tsv" >> \
     "results/merged/summary/unique_variants.txt"