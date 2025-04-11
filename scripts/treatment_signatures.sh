# Script to compare variant contexts
   #!/bin/bash
   # treatment_signatures.sh
   
   mkdir -p results/analysis/signatures
   
   # Analyze top treatment-specific scaffolds
   declare -A TOP_SCAFFOLDS=(
     ["WT"]="JRIU01000157.1"
     ["STC"]="JRIU01000117.1"
     ["CAS"]="JRIU01000289.1"
     ["WTA"]="JRIU01000341.1"
   )
   
   for GROUP in "${!TOP_SCAFFOLDS[@]}"; do
     SCAFFOLD=${TOP_SCAFFOLDS[$GROUP]}
     echo "Analyzing $GROUP-specific hotspot: $SCAFFOLD"
     
     # Extract variant details
     grep "$SCAFFOLD" results/functional/${GROUP}/variant_details.tsv > "results/analysis/signatures/${GROUP}_specific_variants.tsv"
     
     # Count transition vs transversion (for SNPs)
     grep "$SCAFFOLD" results/functional/${GROUP}/variant_details.tsv | \
       awk 'length($3)==1 && length($4)==1' | \
       awk '{print $3">"$4}' | sort | uniq -c > "results/analysis/signatures/${GROUP}_${SCAFFOLD}_mutations.txt"
   done