#!/bin/bash
# position_analysis.sh

mkdir -p results/analysis/position_clustering

# Extract position data from JRIU01000031.1
for GROUP in WT STC CAS WTA; do
  echo "Processing $GROUP..."
  
  # Extract positions
  grep "JRIU01000031.1" results/functional/${GROUP}/variant_details.tsv | cut -f2 > results/analysis/position_clustering/${GROUP}_positions.txt
  
  # Sort positions
  sort -n results/analysis/position_clustering/${GROUP}_positions.txt > results/analysis/position_clustering/${GROUP}_positions_sorted.txt
  
  # Calculate distances between consecutive positions
  awk 'NR>1 {print $1 - prev} {prev=$1}' results/analysis/position_clustering/${GROUP}_positions_sorted.txt > results/analysis/position_clustering/${GROUP}_distances.txt
  
  # Basic statistics about distances
  echo "Statistics for $GROUP distances between variants:" > results/analysis/position_clustering/${GROUP}_stats.txt
  echo "Total variants: $(wc -l < results/analysis/position_clustering/${GROUP}_positions.txt)" >> results/analysis/position_clustering/${GROUP}_stats.txt
  
  if [ -s results/analysis/position_clustering/${GROUP}_distances.txt ]; then
    echo "Minimum distance: $(sort -n results/analysis/position_clustering/${GROUP}_distances.txt | head -1)" >> results/analysis/position_clustering/${GROUP}_stats.txt
    echo "Maximum distance: $(sort -n results/analysis/position_clustering/${GROUP}_distances.txt | tail -1)" >> results/analysis/position_clustering/${GROUP}_stats.txt
    
    # Count small distances (potential clusters, <10bp apart)
    SMALL=$(awk '$1 < 10 {count++} END {print count}' results/analysis/position_clustering/${GROUP}_distances.txt)
    echo "Variants <10bp apart: $SMALL" >> results/analysis/position_clustering/${GROUP}_stats.txt
  else
    echo "Not enough variants for distance calculation" >> results/analysis/position_clustering/${GROUP}_stats.txt
  fi
done

# Create summary
echo "=== Position Clustering Analysis for JRIU01000031.1 ===" > results/analysis/position_clustering/summary.txt
cat results/analysis/position_clustering/*_stats.txt >> results/analysis/position_clustering/summary.txt