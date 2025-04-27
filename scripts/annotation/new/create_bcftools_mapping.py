#!/usr/bin/env python3
"""
create_bcftools_mapping.py

This script creates a chromosome mapping file for bcftools rename operation.
"""

import os
import argparse
import pandas as pd

def create_mapping_file(input_file, output_file, reverse=False):
    """Create a mapping file for bcftools from our chromosome mapping."""
    print(f"Reading mapping from {input_file}")
    df = pd.read_csv(input_file, sep='\t')
    
    # Check column names
    if 'chromosome_id' not in df.columns or 'w303_scaffold' not in df.columns:
        print(f"ERROR: Expected columns 'chromosome_id' and 'w303_scaffold' in {input_file}")
        return False
    
    # Create the mapping (CM007XXX.1 -> w303_scaffold_X)
    if reverse:
        mapping = pd.DataFrame({
            'old': df['w303_scaffold'],
            'new': df['chromosome_id']
        })
    else:
        mapping = pd.DataFrame({
            'old': df['chromosome_id'],
            'new': df['w303_scaffold']
        })
    
    # Save to file without header
    mapping.to_csv(output_file, sep='\t', index=False, header=False)
    
    print(f"Created bcftools mapping file: {output_file}")
    print(f"First few mappings:")
    print(mapping.head().to_string(index=False))
    
    return True

def main():
    parser = argparse.ArgumentParser(description="Create a mapping file for bcftools rename operation")
    parser.add_argument("--input", required=True, help="Input chromosome mapping file (TSV)")
    parser.add_argument("--output", required=True, help="Output bcftools mapping file")
    parser.add_argument("--reverse", action="store_true", help="Create reverse mapping (w303_scaffold_X -> CM007XXX.1)")
    args = parser.parse_args()
    
    create_mapping_file(args.input, args.output, args.reverse)

if __name__ == "__main__":
    main()