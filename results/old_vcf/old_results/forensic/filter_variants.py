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
