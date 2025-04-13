#!/usr/bin/env python3
import sys
import gzip
import os

# Get all filtered VCF files
filtered_dir = "results/vcf/filtered"
vcf_files = [os.path.join(filtered_dir, f) for f in os.listdir(filtered_dir) if f.endswith('.filtered.vcf.gz')]

# Extract sample names
samples = [os.path.basename(f).replace('.filtered.vcf.gz', '') for f in vcf_files]

# Process first 100 lines from each file to build a merged VCF
variants = {}
for i, vcf_file in enumerate(vcf_files):
    with gzip.open(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom, pos, id, ref, alt = fields[0:5]
            key = f"{chrom}_{pos}_{ref}_{alt}"
            if key not in variants:
                variants[key] = fields[0:8] + ["GT"] + ["." for _ in range(len(samples))]
            # Add genotype for this sample
            variants[key][9 + i] = fields[9]
            # Only process first 100 variants
            if len(variants) >= 100:
                break

# Write header
with open("vcf_debug/new_header.txt", 'r') as f:
    header = f.read()
with open("vcf_debug/manual_merged.vcf", 'w') as out:
    out.write(header)
    # Write variants
    for key in sorted(variants.keys()):
        out.write('\t'.join(variants[key]) + '\n')
