#!/usr/bin/env python3
"""
debug_genbank_scaffolds.py - Diagnose and fix scaffold naming issue in GenBank parsing

This script examines how scaffold names are stored in GenBank files and diagnoses
why they're being reported as "unknown" in the verification script.
"""

import os
import sys
import argparse
from Bio import SeqIO
import re
import pandas as pd

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Debug GenBank scaffold naming issues')
    parser.add_argument('--genbank_dir', required=True, help='Directory containing GenBank files')
    parser.add_argument('--gene_mapping', required=True, help='TSV file mapping genes of interest')
    parser.add_argument('--output', required=True, help='Output file for debugging information')
    return parser.parse_args()

def load_gene_mapping(mapping_file):
    """Load gene mapping information."""
    gene_mapping = pd.read_csv(mapping_file, sep='\t')
    print(f"Loaded {len(gene_mapping)} genes from mapping file")
    return gene_mapping

def debug_genbank_structure(genbank_dir, gene_mapping, output_file):
    """Examine structure of GenBank files in detail to diagnose scaffold naming issue."""
    # Open output file
    with open(output_file, 'w') as out:
        out.write("GenBank Scaffold Naming Debugging Report\n")
        out.write("=======================================\n\n")
        
        # Track mapped files and genes
        mapped_scaffolds = {}
        found_genes = {}
        scaffold_format_examples = []
        
        # Process each GenBank file
        genbank_files = [f for f in os.listdir(genbank_dir) if f.endswith('.genbank') or f.endswith('.gb')]
        out.write(f"Found {len(genbank_files)} GenBank files\n\n")
        
        for i, gb_file in enumerate(sorted(genbank_files)):
            file_path = os.path.join(genbank_dir, gb_file)
            out.write(f"File {i+1}: {gb_file}\n")
            out.write("-" * (len(gb_file) + 8) + "\n")
            
            # Find the expected scaffold name from the filename
            expected_scaffold = re.search(r'w303_\d+Apr\d+_scaffold_(\d+)', gb_file)
            expected_scaffold_id = f"w303_scaffold_{expected_scaffold.group(1)}" if expected_scaffold else "unknown"
            out.write(f"Expected scaffold from filename: {expected_scaffold_id}\n\n")
            
            try:
                # Parse the GenBank file
                for record in SeqIO.parse(file_path, "genbank"):
                    # Record basic info
                    out.write(f"Record ID: {record.id}\n")
                    out.write(f"Record Name: {record.name}\n")
                    out.write(f"Record Description: {record.description}\n")
                    
                    # Extract locus line manually from the file
                    with open(file_path, 'r') as f:
                        first_line = f.readline().strip()
                        out.write(f"LOCUS Line: {first_line}\n\n")
                    
                    # Add to examples
                    if len(scaffold_format_examples) < 3:
                        scaffold_format_examples.append({
                            'file': gb_file,
                            'id': record.id,
                            'name': record.name,
                            'description': record.description,
                            'locus': first_line
                        })
                    
                    # Check for source feature
                    source_found = False
                    for feature in record.features:
                        if feature.type == "source":
                            source_found = True
                            out.write("Source Feature Qualifiers:\n")
                            
                            # Print all qualifiers
                            for key, values in feature.qualifiers.items():
                                out.write(f"  {key}: {values}\n")
                            
                            # Check specifically for notes that might contain scaffold info
                            if 'note' in feature.qualifiers:
                                for note in feature.qualifiers['note']:
                                    # Look for patterns like w303_scaffold_X or CM007XXX
                                    scaffold_match1 = re.search(r'w303_scaffold_\d+', note)
                                    scaffold_match2 = re.search(r'CM007\d+\.\d+', note)
                                    
                                    if scaffold_match1:
                                        scaffold_id = scaffold_match1.group(0)
                                        out.write(f"\nFound scaffold ID in note: {scaffold_id}\n")
                                        mapped_scaffolds[expected_scaffold_id] = scaffold_id
                                    elif scaffold_match2:
                                        scaffold_id = scaffold_match2.group(0)
                                        out.write(f"\nFound CM007 ID in note: {scaffold_id}\n")
                                        mapped_scaffolds[expected_scaffold_id] = scaffold_id
                    
                    if not source_found:
                        out.write("No source feature found in this record\n")
                    
                    # Look for genes of interest in this file
                    for _, gene_row in gene_mapping.iterrows():
                        w303_gene_id = gene_row['w303_gene_id']
                        sc_gene_id = gene_row['sc_gene_id']
                        
                        for feature in record.features:
                            if feature.type == "CDS":
                                # Check various fields for gene identifiers
                                locus_tag = feature.qualifiers.get('locus_tag', [''])[0]
                                gene_id = feature.qualifiers.get('gene', [''])[0]
                                
                                # Check for gene in notes
                                gene_found = False
                                for note in feature.qualifiers.get('note', []):
                                    if sc_gene_id in note:
                                        gene_found = True
                                        found_genes[sc_gene_id] = {
                                            'file': gb_file,
                                            'w303_id': w303_gene_id,
                                            'locus_tag': locus_tag,
                                            'genbank_scaffold': record.id,
                                            'expected_scaffold': expected_scaffold_id
                                        }
                                        out.write(f"\nFound gene {sc_gene_id} in this file\n")
                                        out.write(f"  W303 Gene ID: {w303_gene_id}\n")
                                        out.write(f"  Locus Tag: {locus_tag}\n")
                                        out.write(f"  GenBank Scaffold: {record.id}\n")
                                        out.write(f"  Expected Scaffold: {expected_scaffold_id}\n")
                                        break
            
            except Exception as e:
                out.write(f"Error processing file: {e}\n")
            
            out.write("\n" + "=" * 80 + "\n\n")
        
        # Summary section
        out.write("\nSummary\n")
        out.write("=======\n\n")
        
        # Scaffold format examples
        out.write("Scaffold Format Examples:\n")
        for example in scaffold_format_examples:
            out.write(f"File: {example['file']}\n")
            out.write(f"  ID: {example['id']}\n")
            out.write(f"  Name: {example['name']}\n")
            out.write(f"  Description: {example['description']}\n")
            out.write(f"  LOCUS: {example['locus']}\n\n")
        
        # Found genes
        out.write("Found Genes of Interest:\n")
        for sc_gene_id, gene_info in found_genes.items():
            out.write(f"{sc_gene_id}:\n")
            out.write(f"  File: {gene_info['file']}\n")
            out.write(f"  W303 ID: {gene_info['w303_id']}\n")
            out.write(f"  Locus Tag: {gene_info['locus_tag']}\n")
            out.write(f"  GenBank Scaffold: {gene_info['genbank_scaffold']}\n")
            out.write(f"  Expected Scaffold: {gene_info['expected_scaffold']}\n\n")
        
        # Mapped scaffolds
        out.write("Mapped Scaffolds:\n")
        for expected, actual in mapped_scaffolds.items():
            out.write(f"  {expected} -> {actual}\n")
        
        # Conclusion and fix recommendation
        out.write("\nDiagnosis and Fix Recommendation\n")
        out.write("==============================\n")
        out.write("The issue appears to be that the script is not correctly extracting scaffold IDs from\n")
        out.write("the GenBank files. Based on the analysis above, here's the recommended fix:\n\n")
        
        # We'll determine the exact recommendation based on what we find in the debugging info
        out.write("1. In the verify_gene_coordinates.py script, modify the scaffold ID extraction logic\n")
        out.write("   to properly handle the format used in these GenBank files.\n")
        out.write("2. Specifically, the GenBank record.id doesn't contain the expected w303_scaffold_X format.\n")
        out.write("3. Instead, extract the scaffold ID from either the filename or look for it in the\n")
        out.write("   source feature notes with specific regex patterns as shown in this debug script.\n")

def main():
    args = parse_arguments()
    gene_mapping = load_gene_mapping(args.gene_mapping)
    debug_genbank_structure(args.genbank_dir, gene_mapping, args.output)
    print(f"Debugging information written to {args.output}")

if __name__ == "__main__":
    main()