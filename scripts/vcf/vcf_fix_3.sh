#!/bin/bash

echo "=== Explicit File Path Debugging ==="

# 1. Check if the specific variant files exist
echo "Checking if variant files exist:"
for GROUP in WT STC CAS WTA; do
  echo -n "$GROUP: "
  ls final_analysis/$GROUP/*_specific.vcf 2>/dev/null || echo "No files found!"
done

# 2. Examine content of a specific file
echo -e "\nExamining content of a specific file:"
FILE=$(ls final_analysis/WT/*_specific.vcf 2>/dev/null | head -1)
if [ -n "$FILE" ]; then
  echo "First 3 lines of $FILE:"
  head -3 "$FILE"
  echo "Variant count: $(grep -v '^#' "$FILE" | wc -l)"
else
  echo "No variant files found to examine"
fi

# 3. Try extremely explicit variant counting
echo -e "\nExplicit variant counting for WT group:"

# Extract chromosomes and positions explicitly
for SAMPLE in WT-37-55-1 WT-37-55-2 WT-37-55-3; do
  VARFILE="final_analysis/WT/${SAMPLE}_specific.vcf"
  if [ -f "$VARFILE" ]; then
    grep -v "^#" "$VARFILE" | awk '{print $1"_"$2"_"$4"_"$5}' > "debug_${SAMPLE}_variants.txt"
    echo "$SAMPLE: $(wc -l < debug_${SAMPLE}_variants.txt) variants"
  else
    echo "$SAMPLE: File not found!"
  fi
done

# Count shared variants explicitly
echo -e "\nFinding shared variants:"
cat debug_WT*.txt | sort | uniq -c | sort -nr > debug_all_variants.txt
SHARED=$(awk '$1 >= 2' debug_all_variants.txt | wc -l)
echo "Variants in ≥2 samples: $SHARED"
echo "Top 5 shared variants:"
head -5 debug_all_variants.txt

# 4. Try a new approach for the merged VCF
echo -e "\nTesting alternative merged VCF approach:"
# Let's try using a specific merge approach
bcftools merge --force-samples fixed_vcfs/WT-CTRL_fixed.vcf.gz fixed_vcfs/WT-37-55-1_fixed.vcf.gz -Oz -o debug_test_merge.vcf.gz
tabix -p vcf debug_test_merge.vcf.gz

# Check if this works better
echo "Testing sample selection on this merge:"
bcftools query -l debug_test_merge.vcf.gz

echo "=== Debugging Complete ==="