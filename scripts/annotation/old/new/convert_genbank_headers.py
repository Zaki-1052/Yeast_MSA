#!/usr/bin/env python3
"""
convert_genbank_headers.py - Convert YGAP GenBank headers to match VCF chromosome naming
"""

import os
import re
import sys

# Configuration
INPUT_DIR = "annotation/reference/w303_scaffolds"
OUTPUT_DIR = "annotation/reference/w303_scaffolds_converted"
MAPPING_FILE = "annotation/reference/header_mapping.txt"

def main():
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Read mapping file
    mapping = {}
    reverse_mapping = {}
    with open(MAPPING_FILE, 'r') as f:
        # Skip header line
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                original = parts[0].strip('>').split()[0]  # Extract CM007XXX.1
                ygap = parts[1].strip('>')                # Extract I, II, etc.
                mapping[ygap] = original
                reverse_mapping[original] = ygap
    
    print(f"Loaded {len(mapping)} chromosome mappings")
    
    # Process each GenBank file
    gb_files = [f for f in os.listdir(INPUT_DIR) if f.endswith('.gb') or f.endswith('.gbk')]
    print(f"Found {len(gb_files)} GenBank files to process")
    
    for gb_file in gb_files:
        input_path = os.path.join(INPUT_DIR, gb_file)
        output_path = os.path.join(OUTPUT_DIR, gb_file)
        
        print(f"Processing {gb_file}...")
        
        with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
            # Read the content
            content = infile.read()
            
            # Find chromosome identifier (could be Roman numeral or other format)
            # This is a simplistic approach - might need refinement based on actual file format
            chrom_match = re.search(r'LOCUS\s+(\S+)', content)
            if chrom_match:
                ygap_chrom = chrom_match.group(1)
                
                # Try to match against our mapping
                for key in mapping:
                    if key in ygap_chrom or ygap_chrom.endswith(key):
                        # Replace with CM007XXX.1 format
                        new_content = content.replace(ygap_chrom, mapping[key])
                        outfile.write(new_content)
                        print(f"  Converted {ygap_chrom} to {mapping[key]}")
                        break
                else:
                    # No replacement found, write original
                    outfile.write(content)
                    print(f"  WARNING: No mapping found for {ygap_chrom}")
            else:
                # No LOCUS line found, write original
                outfile.write(content)
                print(f"  WARNING: No chromosome identifier found in {gb_file}")
    
    print("Conversion complete!")

if __name__ == "__main__":
    main()