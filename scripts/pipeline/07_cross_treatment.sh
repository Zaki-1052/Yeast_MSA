#!/bin/bash

echo "=== Creating cross-treatment comparison ==="

# Create combined variant list
echo -e "Chrom\tPos\tRef\tAlt\tWT-37\tSTC\tCAS\tWTA" > \
  "results/merged/summary/treatment_comparison.tsv"

# Extract variant details for all high-confidence variants from primary comparisons
for GROUP in WT-37 STC CAS WTA; do
    if [ -f "results/merged/analysis/${GROUP}/highconf.vcf.gz" ]; then
        echo "  Extracting variant details for $GROUP..."
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
          "results/merged/analysis/${GROUP}/highconf.vcf.gz" > \
          "results/merged/summary/${GROUP}/highconf_variants.txt"
    else
        echo "  WARNING: No high-confidence variants file for $GROUP"
        mkdir -p "results/merged/summary/${GROUP}"
        touch "results/merged/summary/${GROUP}/highconf_variants.txt"
    fi
done

# Create a master list of all high-confidence variant positions
cat results/merged/summary/*/highconf_variants.txt | \
  sort -k1,1 -k2,2n | uniq > \
  "results/merged/summary/all_highconf_positions.txt"

# Process each variant position
while read CHROM POS REF ALT; do
    if [ -z "$CHROM" ] || [ -z "$POS" ] || [ -z "$REF" ] || [ -z "$ALT" ]; then
        continue  # Skip empty lines
    fi
    
    LINE="$CHROM\t$POS\t$REF\t$ALT"
    
    # Check each treatment group
    for GROUP in WT-37 STC CAS WTA; do
        if [ -f "results/merged/summary/${GROUP}/highconf_variants.txt" ] && \
           grep -q "^$CHROM[[:space:]]$POS[[:space:]]$REF[[:space:]]$ALT$" \
           "results/merged/summary/${GROUP}/highconf_variants.txt"; then
            LINE="$LINE\t1"
        else
            LINE="$LINE\t0"
        fi
    done
    
    echo -e "$LINE" >> "results/merged/summary/treatment_comparison.tsv"
done < "results/merged/summary/all_highconf_positions.txt"