# Create a script to analyze position clustering
   #!/bin/bash
   # hotspot_position_analysis.sh
   
   mkdir -p results/analysis/position_clustering
   
   # Extract position data from JRIU01000031.1
   for GROUP in WT STC CAS WTA; do
     POSITIONS_FILE="results/analysis/position_clustering/${GROUP}_JRIU01000031.1_positions.txt"
     grep "JRIU01000031.1" results/functional/${GROUP}/variant_details.tsv | cut -f2 > $POSITIONS_FILE
     
     # Sort positions numerically
     sort -n $POSITIONS_FILE > "${POSITIONS_FILE}.sorted"
     mv "${POSITIONS_FILE}.sorted" $POSITIONS_FILE
     
     # Calculate distances between consecutive positions
     awk 'NR>1 {print $1 - prev} {prev=$1}' $POSITIONS_FILE > "results/analysis/position_clustering/${GROUP}_distances.txt"
   done