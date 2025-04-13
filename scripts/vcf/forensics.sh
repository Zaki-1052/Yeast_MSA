#!/bin/bash

# Create directory for deeper analysis
mkdir -p results/vcf/forensic

# Examine the exact VCF header format
echo "Examining VCF header format..."
bcftools view -h results/vcf/merged/all_samples.vcf.gz | grep "^#CHROM" > results/vcf/forensic/header_line.txt
cat results/vcf/forensic/header_line.txt

# Try a direct approach using grep and awk to filter variants
# This bypasses the sample selection mechanism entirely

echo "Using direct filtering approach for WT group..."

# 1. Extract just the variant lines and crucial columns
bcftools view -H results/vcf/merged/all_samples.vcf.gz | \
  awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24}' > results/vcf/forensic/all_variants.txt

# Create a WT-specific header
echo -e "CHROM\tPOS\tREF\tALT\tFORMAT\tCAS-55-1\tCAS-55-2\tCAS-55-3\tCAS-CTRL\tSTC-55-1\tSTC-55-2\tSTC-55-3\tSTC-CTRL\tWT-37-55-1\tWT-37-55-2\tWT-37-55-3\tWT-CTRL\tWTA-55-1\tWTA-55-2\tWTA-55-3" > results/vcf/forensic/header.txt

# 2. Write a Python script to handle the filtering directly
cat << 'EOF' > results/vcf/forensic/filter_variants.py
#!/usr/bin/env python3
import sys

# Define sample indices (0-based, after the first 5 columns)
sample_indices = {
    "CAS-55-1": 5, "CAS-55-2": 6, "CAS-55-3": 7, "CAS-CTRL": 8,
    "STC-55-1": 9, "STC-55-2": 10, "STC-55-3": 11, "STC-CTRL": 12,
    "WT-37-55-1": 13, "WT-37-55-2": 14, "WT-37-55-3": 15, "WT-CTRL": 16,
    "WTA-55-1": 17, "WTA-55-2": 18, "WTA-55-3": 19
}

# Get group definition from command line
group = sys.argv[1]
if group == "WT":
    treat_samples = ["WT-37-55-1", "WT-37-55-2", "WT-37-55-3"]
    ctrl_samples = ["WT-CTRL"]
elif group == "STC":
    treat_samples = ["STC-55-1", "STC-55-2", "STC-55-3"]
    ctrl_samples = ["STC-CTRL"]
elif group == "CAS":
    treat_samples = ["CAS-55-1", "CAS-55-2", "CAS-55-3"]
    ctrl_samples = ["CAS-CTRL"]
elif group == "WTA":
    treat_samples = ["WTA-55-1", "WTA-55-2", "WTA-55-3"]
    ctrl_samples = ["WT-CTRL"]

# Get high-confidence threshold from command line
threshold = int(sys.argv[2]) if len(sys.argv) > 2 else 1

# Process each line
treatment_specific = []
for line in sys.stdin:
    if line.startswith("CHROM"):
        continue  # Skip header
        
    fields = line.strip().split("\t")
    
    # Check if any treatment sample has variant
    treat_has_variant = False
    treat_count = 0
    for sample in treat_samples:
        idx = sample_indices[sample]
        if idx < len(fields) and fields[idx].startswith("1"):
            treat_has_variant = True
            treat_count += 1
            
    # Check if any control sample has variant
    ctrl_has_variant = False
    for sample in ctrl_samples:
        idx = sample_indices[sample]
        if idx < len(fields) and fields[idx].startswith("1"):
            ctrl_has_variant = True
            
    # Output treatment-specific variants meeting threshold
    if treat_has_variant and not ctrl_has_variant and treat_count >= threshold:
        print('\t'.join(fields[:5]))  # Just output variant info columns

# Print summary
print(f"Found {len(treatment_specific)} {group} variants with {threshold}+ samples", file=sys.stderr)
EOF

chmod +x results/vcf/forensic/filter_variants.py

# 3. Run the filtering script for each group
echo "Filtering variants directly for WT group..."
cat results/vcf/forensic/all_variants.txt | results/vcf/forensic/filter_variants.py WT 1 > results/vcf/forensic/WT_specific.txt
cat results/vcf/forensic/all_variants.txt | results/vcf/forensic/filter_variants.py WT 2 > results/vcf/forensic/WT_highconf.txt

echo "Filtering variants directly for STC group..."
cat results/vcf/forensic/all_variants.txt | results/vcf/forensic/filter_variants.py STC 1 > results/vcf/forensic/STC_specific.txt
cat results/vcf/forensic/all_variants.txt | results/vcf/forensic/filter_variants.py STC 2 > results/vcf/forensic/STC_highconf.txt

echo "Filtering variants directly for CAS group..."
cat results/vcf/forensic/all_variants.txt | results/vcf/forensic/filter_variants.py CAS 1 > results/vcf/forensic/CAS_specific.txt
cat results/vcf/forensic/all_variants.txt | results/vcf/forensic/filter_variants.py CAS 2 > results/vcf/forensic/CAS_highconf.txt

echo "Filtering variants directly for WTA group..."
cat results/vcf/forensic/all_variants.txt | results/vcf/forensic/filter_variants.py WTA 1 > results/vcf/forensic/WTA_specific.txt
cat results/vcf/forensic/all_variants.txt | results/vcf/forensic/filter_variants.py WTA 2 > results/vcf/forensic/WTA_highconf.txt

# 4. Count results
for GROUP in WT STC CAS WTA; do
    SPECIFIC=$(wc -l < results/vcf/forensic/${GROUP}_specific.txt)
    HIGHCONF=$(wc -l < results/vcf/forensic/${GROUP}_highconf.txt)
    echo "$GROUP treatment-specific variants: $SPECIFIC"
    echo "$GROUP high-confidence treatment-specific variants: $HIGHCONF"
done

echo "Direct filtering analysis complete."