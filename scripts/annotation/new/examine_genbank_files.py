#!/usr/bin/env python3
"""
locate_genbank_files.py - Find and examine GenBank files
"""

import os
import re
import sys

def main():
    # Possible locations to check
    locations = [
        "annotation/reference/w303_scaffolds",
        "annotation/reference",
        "annotation",
        "."
    ]
    
    genbank_files = []
    
    # Search for GenBank files
    for loc in locations:
        if os.path.exists(loc):
            for root, dirs, files in os.walk(loc):
                gb_files = [os.path.join(root, f) for f in files 
                           if f.endswith('.gb') or f.endswith('.gbk') or f.endswith('.genbank')]
                genbank_files.extend(gb_files)
    
    print(f"Found {len(genbank_files)} GenBank files")
    
    # Examine first few files to determine format
    if genbank_files:
        print("\nExamining first 3 GenBank files:")
        
        for i, gb_file in enumerate(genbank_files[:3]):
            print(f"\nFile {i+1}: {gb_file}")
            
            with open(gb_file, 'r') as f:
                # Read first 10 lines to check format
                lines = [f.readline().strip() for _ in range(10)]
                
                # Look for chromosome identifiers
                locus_line = next((line for line in lines if line.startswith('LOCUS')), None)
                if locus_line:
                    print(f"  LOCUS line: {locus_line}")
                
                # Check for CM007XXX.1 format
                cm_format = any('CM007' in line for line in lines)
                roman_format = any(re.search(r'\b(I|II|III|IV|V|VI|VII|VIII|IX|X|XI|XII|XIII|XIV|XV|XVI)\b', line) for line in lines)
                
                print(f"  Contains CM007XXX.1 format: {cm_format}")
                print(f"  Contains Roman numeral format: {roman_format}")
    else:
        print("No GenBank files found. Please check the OneDrive folder from Israel.")

if __name__ == "__main__":
    main()