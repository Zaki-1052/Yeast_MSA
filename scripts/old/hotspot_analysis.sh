#!/bin/bash
# hotspot_analysis.sh - Characterize variant hotspots across treatments

# Create output directory for hotspot analysis
mkdir -p results/analysis/hotspots

# 1. Get scaffold lengths from reference genome
echo "Extracting scaffold lengths from reference..."
samtools faidx reference/yeast_w303.fasta
cut -f1,2 reference/yeast_w303.fasta.fai > results/analysis/hotspots/scaffold_lengths.txt

# 2. Calculate variant density for each treatment
echo "Calculating variant density across scaffolds..."

# Process each treatment group
for GROUP in WT STC CAS WTA; do
  if [ -f "results/functional/${GROUP}/scaffold_distribution.txt" ]; then
    echo "Processing ${GROUP} variants..."
    
    # Create output file for variant density
    echo -e "Scaffold\tVariants\tLength_kb\tDensity_per_kb" > results/analysis/hotspots/${GROUP}_density.tsv
    
    # Join scaffold counts with lengths and calculate density
    while read COUNT SCAFFOLD; do
      # Get scaffold length from the scaffold lengths file
      LENGTH=$(grep -w "^${SCAFFOLD}" results/analysis/hotspots/scaffold_lengths.txt | cut -f2)
      
      if [ ! -z "$LENGTH" ]; then
        # Calculate density (variants per kb)
        KB_LENGTH=$(echo "scale=3; ${LENGTH}/1000" | bc)
        DENSITY=$(echo "scale=3; ${COUNT}/${KB_LENGTH}" | bc)
        
        # Write to output file
        echo -e "${SCAFFOLD}\t${COUNT}\t${KB_LENGTH}\t${DENSITY}" >> results/analysis/hotspots/${GROUP}_density.tsv
      fi
    done < results/functional/${GROUP}/scaffold_distribution.txt
    
    # Sort by density (descending) and get top 10
    echo "Top 10 variant-dense scaffolds for ${GROUP}:" > results/analysis/hotspots/${GROUP}_top_hotspots.txt
    sort -k4,4nr results/analysis/hotspots/${GROUP}_density.tsv | head -11 >> results/analysis/hotspots/${GROUP}_top_hotspots.txt
  fi
done

# 3. Find common hotspots across treatments
echo "Identifying common hotspots across treatments..."
# Extract scaffold lists from each treatment
for GROUP in WT STC CAS WTA; do
  if [ -f "results/analysis/hotspots/${GROUP}_density.tsv" ]; then
    # Get top 20 scaffolds by density
    tail -n +2 results/analysis/hotspots/${GROUP}_density.tsv | sort -k4,4nr | head -20 | cut -f1 > results/analysis/hotspots/${GROUP}_top20.txt
  fi
done

# Find scaffolds that appear as hotspots in multiple treatments
cat results/analysis/hotspots/*_top20.txt | sort | uniq -c | sort -nr > results/analysis/hotspots/shared_hotspots.txt

# 4. Analyze variant positions within top hotspot (JRIU01000031.1)
echo "Analyzing variant positions in top hotspot scaffold..."
mkdir -p results/analysis/hotspots/position_analysis

for GROUP in WT STC CAS WTA; do
  if [ -f "results/functional/${GROUP}/variant_details.tsv" ]; then
    # Extract positions in JRIU01000031.1
    grep "JRIU01000031.1" results/functional/${GROUP}/variant_details.tsv | cut -f1,2,3,4,5,6 > results/analysis/hotspots/position_analysis/${GROUP}_JRIU01000031.1_variants.tsv
    
    # Count SNPs vs indels in this scaffold
    TOTAL=$(grep "JRIU01000031.1" results/functional/${GROUP}/variant_details.tsv | wc -l)
    SNPS=$(grep "JRIU01000031.1" results/functional/${GROUP}/variant_details.tsv | awk 'length($3)==1 && length($4)==1' | wc -l)
    INDELS=$(grep "JRIU01000031.1" results/functional/${GROUP}/variant_details.tsv | awk 'length($3)!=1 || length($4)!=1' | wc -l)
    
    echo "${GROUP} variants in JRIU01000031.1: Total=${TOTAL}, SNPs=${SNPS}, Indels=${INDELS}" >> results/analysis/hotspots/position_analysis/mutation_types.txt
  fi
done

# 5. Generate summary of hotspot analysis
echo "Creating hotspot analysis summary..."
{
  echo "===== Variant Hotspot Analysis ====="
  echo "Date: $(date)"
  echo ""
  echo "Top shared hotspot scaffolds across treatments:"
  cat results/analysis/hotspots/shared_hotspots.txt | head -10
  echo ""
  echo "Variant counts in top hotspot scaffold (JRIU01000031.1):"
  cat results/analysis/hotspots/position_analysis/mutation_types.txt
  echo ""
  echo "Variant density statistics:"
  for GROUP in WT STC CAS WTA; do
    if [ -f "results/analysis/hotspots/${GROUP}_density.tsv" ]; then
      MAX=$(tail -n +2 results/analysis/hotspots/${GROUP}_density.tsv | sort -k4,4nr | head -1)
      echo "${GROUP} max density: ${MAX}"
    fi
  done
} > results/analysis/hotspots/hotspot_summary.txt

echo "Hotspot analysis complete. Results in results/analysis/hotspots/"