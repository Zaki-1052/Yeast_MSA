#!/usr/bin/env python3
"""
create_chromosome_mapping.py

This script creates a bidirectional mapping between chromosome identifiers 
(CM007XXX.1 or LYZE01XXXXXX.1) and w303_scaffold_X identifiers by parsing GenBank annotation files.

Usage:
    python create_chromosome_mapping.py --genbank_dir <genbank_directory> --output_dir <output_directory>
"""

import os
import argparse
from Bio import SeqIO
import pandas as pd
import re

def parse_genbank_files(genbank_dir):
    """
    Parse GenBank files to extract chromosome mapping information.
    
    Args:
        genbank_dir (str): Directory containing GenBank files
        
    Returns:
        list: List of dictionaries with mapping information
    """
    mapping_data = []
    
    # List all GenBank files
    genbank_files = [f for f in os.listdir(genbank_dir) if f.endswith('.genbank')]
    
    if not genbank_files:
        print(f"No GenBank files found in {genbank_dir}")
        return mapping_data
    
    print(f"Found {len(genbank_files)} GenBank files to process")
    
    # Process each GenBank file
    for genbank_file in sorted(genbank_files):
        file_path = os.path.join(genbank_dir, genbank_file)
        print(f"Processing {genbank_file}...")
        
        try:
            # Parse GenBank file
            record = SeqIO.read(file_path, "genbank")
            
            # Extract scaffold name from LOCUS line
            scaffold_name = record.name  # This should be w303_scaffold_X
            
            # Extract chromosome identifier from /note field
            chromosome_id = None
            scaffold_number = None
            
            for feature in record.features:
                if feature.type == "source":
                    notes = feature.qualifiers.get("note", [])
                    
                    for note in notes:
                        # Look for CM007XXX.1 or LYZE01XXXXXX.1 pattern
                        if (note.startswith("CM007") and "." in note) or (note.startswith("LYZE") and "." in note):
                            chromosome_id = note.strip()
                        # Look for scaffold number
                        elif note.startswith("scaffold"):
                            scaffold_number = note.strip().split()[1]
            
            # Add to mapping data
            if chromosome_id and scaffold_name:
                entry = {
                    "w303_scaffold": scaffold_name,
                    "scaffold_number": scaffold_number,
                    "chromosome_id": chromosome_id,
                    "length": len(record),
                    "description": record.description
                }
                mapping_data.append(entry)
                print(f"  Mapped {chromosome_id} -> {scaffold_name} (Length: {len(record)} bp)")
            else:
                print(f"  WARNING: Could not extract mapping information from {genbank_file}")
                
        except Exception as e:
            print(f"  ERROR processing {genbank_file}: {str(e)}")
    
    return mapping_data

def create_mapping_files(mapping_data, output_dir):
    """
    Create mapping files from the extracted data.
    
    Args:
        mapping_data (list): List of dictionaries with mapping information
        output_dir (str): Directory to save output files
    """
    if not mapping_data:
        print("No mapping data to save")
        return
    
    # Create DataFrame
    df = pd.DataFrame(mapping_data)
    
    # Sort by scaffold number
    df['scaffold_num'] = df['scaffold_number'].astype(int)
    df = df.sort_values('scaffold_num')
    df = df.drop('scaffold_num', axis=1)
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Create forward mapping (chromosome_id -> w303_scaffold_X)
    forward_mapping = df[['chromosome_id', 'w303_scaffold']].copy()
    forward_mapping_file = os.path.join(output_dir, "chromosome_mapping.tsv")
    forward_mapping.to_csv(forward_mapping_file, sep='\t', index=False)
    print(f"Created forward mapping file: {forward_mapping_file}")
    
    # Create reverse mapping (w303_scaffold_X -> chromosome_id)
    reverse_mapping = df[['w303_scaffold', 'chromosome_id']].copy()
    reverse_mapping_file = os.path.join(output_dir, "chromosome_mapping_reverse.tsv")
    reverse_mapping.to_csv(reverse_mapping_file, sep='\t', index=False)
    print(f"Created reverse mapping file: {reverse_mapping_file}")
    
    # Create comprehensive summary
    summary_file = os.path.join(output_dir, "chromosome_summary.tsv")
    df.to_csv(summary_file, sep='\t', index=False)
    print(f"Created chromosome summary file: {summary_file}")
    
    # Create SnpEff synonym file format
    snpeff_synonyms = df[['chromosome_id', 'w303_scaffold']].copy()
    snpeff_file = os.path.join(output_dir, "snpeff_chromosome_synonyms.txt")
    snpeff_synonyms.to_csv(snpeff_file, sep='\t', index=False, header=False)
    print(f"Created SnpEff chromosome synonyms file: {snpeff_file}")
    
    # Print summary
    print(f"\nSummary: Mapped {len(df)} chromosomes")
    print("\nFirst few mappings:")
    print(df[['chromosome_id', 'w303_scaffold', 'length']].head().to_string(index=False))

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Create chromosome mapping files from GenBank annotations")
    parser.add_argument("--genbank_dir", required=True, help="Directory containing GenBank files")
    parser.add_argument("--output_dir", required=True, help="Directory to save output files")
    args = parser.parse_args()
    
    # Process GenBank files
    mapping_data = parse_genbank_files(args.genbank_dir)
    
    # Create mapping files
    create_mapping_files(mapping_data, args.output_dir)

if __name__ == "__main__":
    main()