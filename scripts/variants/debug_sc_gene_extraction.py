#!/usr/bin/env python3
"""
debug_sc_gene_extraction.py - Debug issue with SC gene ID extraction from GenBank files

This script performs detailed analysis of how SC gene IDs (like YHR190W) appear in GenBank files
and why they might not be correctly matched in the verification script.
"""

import os
import sys
import argparse
import re
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Debug SC gene ID extraction from GenBank files')
    parser.add_argument('--genbank_dir', required=True, help='Directory containing GenBank files')
    parser.add_argument('--gene_mapping', required=True, help='TSV file mapping genes of interest')
    parser.add_argument('--output', required=True, help='Output file for debugging information')
    return parser.parse_args()

def load_gene_mapping(mapping_file):
    """Load gene mapping information."""
    gene_mapping = pd.read_csv(mapping_file, sep='\t')
    print(f"Loaded {len(gene_mapping)} genes from mapping file")
    
    # Extract SC gene IDs for easier lookup
    sc_gene_ids = set(gene_mapping['sc_gene_id'].values)
    print(f"Target SC gene IDs: {', '.join(sorted(sc_gene_ids))}")
    
    return gene_mapping, sc_gene_ids

def analyze_genbank_fields(genbank_dir, sc_gene_ids, output_file):
    """Analyze how SC gene IDs appear in GenBank files."""
    # Data collection structures
    gene_matches = defaultdict(list)  # Tracks all potential matches for each SC gene ID
    inference_formats = set()         # Tracks unique formats of inference field
    note_formats = set()              # Tracks unique formats of note field
    found_genes = {}                  # Successfully identified genes
    
    # Open output file
    with open(output_file, 'w') as out:
        out.write("SC Gene ID Extraction Debugging Report\n")
        out.write("====================================\n\n")
        
        # Process each GenBank file
        genbank_files = [f for f in os.listdir(genbank_dir) if f.endswith('.genbank') or f.endswith('.gb')]
        out.write(f"Found {len(genbank_files)} GenBank files\n\n")
        
        for gb_file in sorted(genbank_files):
            file_path = os.path.join(genbank_dir, gb_file)
            out.write(f"Processing {gb_file}...\n")
            
            # Extract scaffold name from filename
            scaffold_match = re.search(r'w303_\d+Apr\d+_scaffold_(\d+)', gb_file)
            scaffold_id = f"w303_scaffold_{scaffold_match.group(1)}" if scaffold_match else "unknown"
            
            # Track potential matches in this file
            file_matches = []
            
            try:
                # Parse the GenBank file
                for record in SeqIO.parse(file_path, "genbank"):
                    # Process each feature
                    for feature in record.features:
                        if feature.type == "CDS":
                            # Extract all relevant fields
                            w303_gene_id = feature.qualifiers.get('gene', [''])[0]
                            locus_tag = feature.qualifiers.get('locus_tag', [''])[0]
                            
                            # Analyze inference field
                            inferences = feature.qualifiers.get('inference', [])
                            for inf in inferences:
                                # Add to unique formats
                                inference_formats.add(inf)
                                
                                # Look for SC gene ID pattern in inference
                                sc_match = re.search(r'Y[A-Z]{2}\d{3}[WC]', inf)
                                if sc_match:
                                    potential_sc_id = sc_match.group(0)
                                    if potential_sc_id in sc_gene_ids:
                                        match_info = {
                                            'file': gb_file,
                                            'scaffold': scaffold_id,
                                            'w303_id': w303_gene_id,
                                            'locus_tag': locus_tag,
                                            'source': 'inference',
                                            'pattern': inf,
                                            'location': feature.location
                                        }
                                        gene_matches[potential_sc_id].append(match_info)
                                        file_matches.append((potential_sc_id, match_info))
                                        found_genes[potential_sc_id] = match_info
                            
                            # Analyze note field
                            notes = feature.qualifiers.get('note', [])
                            for note in notes:
                                # Add to unique formats
                                if len(note) < 100:  # Avoid overly long notes
                                    note_formats.add(note)
                                    
                                # Look for "similar to Saccharomyces cerevisiae X (YYY)" pattern
                                sc_match = re.search(r'similar to Saccharomyces cerevisiae\s+(\S+)\s+\(([A-Z]{3}\d{3}[WC])\)', note)
                                if sc_match:
                                    gene_name = sc_match.group(1)
                                    potential_sc_id = sc_match.group(2)
                                    if potential_sc_id in sc_gene_ids:
                                        match_info = {
                                            'file': gb_file,
                                            'scaffold': scaffold_id,
                                            'w303_id': w303_gene_id,
                                            'locus_tag': locus_tag,
                                            'source': 'note',
                                            'pattern': note,
                                            'gene_name': gene_name,
                                            'location': feature.location
                                        }
                                        gene_matches[potential_sc_id].append(match_info)
                                        file_matches.append((potential_sc_id, match_info))
                                        found_genes[potential_sc_id] = match_info
                                
                                # Alternative: direct lookup for SC gene ID pattern
                                sc_match = re.search(r'Y[A-Z]{2}\d{3}[WC]', note)
                                if sc_match and not potential_sc_id in sc_gene_ids:
                                    potential_sc_id = sc_match.group(0)
                                    if potential_sc_id in sc_gene_ids:
                                        match_info = {
                                            'file': gb_file,
                                            'scaffold': scaffold_id,
                                            'w303_id': w303_gene_id,
                                            'locus_tag': locus_tag,
                                            'source': 'note_direct',
                                            'pattern': note,
                                            'location': feature.location
                                        }
                                        gene_matches[potential_sc_id].append(match_info)
                                        file_matches.append((potential_sc_id, match_info))
                                        found_genes[potential_sc_id] = match_info
            
                # Report matches in this file
                if file_matches:
                    out.write(f"  Found {len(file_matches)} potential matches in this file:\n")
                    for sc_id, match in file_matches:
                        src = match['source']
                        out.write(f"    {sc_id} - {match['w303_id']} ({src}): {match['pattern'][:60]}...\n")
                else:
                    out.write("  No matches found in this file\n")
                
            except Exception as e:
                out.write(f"  Error processing file: {e}\n")
            
            out.write("\n")
        
        # Summary section
        out.write("\nSummary\n")
        out.write("=======\n\n")
        
        # SC Gene ID matches
        out.write("SC Gene ID Matches:\n")
        for sc_id in sorted(sc_gene_ids):
            matches = gene_matches[sc_id]
            if matches:
                out.write(f"{sc_id}: {len(matches)} matches\n")
                for i, match in enumerate(matches, 1):
                    out.write(f"  Match {i}:\n")
                    out.write(f"    File: {match['file']}\n")
                    out.write(f"    Scaffold: {match['scaffold']}\n")
                    out.write(f"    W303 ID: {match['w303_id']}\n")
                    out.write(f"    Locus Tag: {match['locus_tag']}\n")
                    out.write(f"    Source: {match['source']}\n")
                    if 'gene_name' in match:
                        out.write(f"    Gene Name: {match['gene_name']}\n")
                    out.write(f"    Pattern: {match['pattern'][:60]}...\n")
                    out.write(f"    Location: {match['location']}\n\n")
            else:
                out.write(f"{sc_id}: No matches found\n\n")
        
        # Unique inference formats
        out.write("\nUnique Inference Field Formats:\n")
        for fmt in sorted(inference_formats):
            out.write(f"  {fmt}\n")
        
        # Unique note formats relevant to genes
        out.write("\nSample Note Field Formats:\n")
        for fmt in sorted(list(note_formats)[:10]):
            out.write(f"  {fmt}\n")
        
        # Analysis and recommendations
        out.write("\nAnalysis and Recommendations\n")
        out.write("===========================\n")
        
        # Count successful matches
        found_count = len(found_genes)
        out.write(f"Total SC gene IDs found: {found_count}/{len(sc_gene_ids)}\n\n")
        
        if found_count < len(sc_gene_ids):
            out.write("Missing SC gene IDs:\n")
            for sc_id in sorted(sc_gene_ids):
                if sc_id not in found_genes:
                    out.write(f"  {sc_id}\n")
            out.write("\n")
        
        # Analyze common patterns
        inference_patterns = defaultdict(int)
        note_patterns = defaultdict(int)
        
        for matches in gene_matches.values():
            for match in matches:
                if match['source'] == 'inference':
                    pattern = re.sub(r'Y[A-Z]{2}\d{3}[WC]', 'YXXXNNNX', match['pattern'])
                    inference_patterns[pattern] += 1
                elif match['source'] == 'note':
                    pattern = re.sub(r'Y[A-Z]{2}\d{3}[WC]', 'YXXXNNNX', match['pattern'])
                    pattern = re.sub(r'\S+\s+\(', 'GENE_NAME (', pattern)
                    note_patterns[pattern] += 1
        
        out.write("Common inference patterns:\n")
        for pattern, count in sorted(inference_patterns.items(), key=lambda x: x[1], reverse=True):
            out.write(f"  {pattern} ({count} occurrences)\n")
        
        out.write("\nCommon note patterns:\n")
        for pattern, count in sorted(note_patterns.items(), key=lambda x: x[1], reverse=True):
            out.write(f"  {pattern} ({count} occurrences)\n")
        
        # Recommendations
        out.write("\nRecommended fixes for verify_gene_coordinates.py:\n")
        out.write("1. Update the SC gene ID extraction logic to check both inference and note fields\n")
        out.write("2. Use regular expressions to extract the SC gene IDs instead of simple word matching\n")
        out.write("3. Specific regex patterns to use:\n")
        out.write("   - For inference field: r'Y[A-Z]{2}\\d{3}[WC]'\n")
        out.write("   - For note field: r'similar to Saccharomyces cerevisiae\\s+(\\S+)\\s+\\(([A-Z]{3}\\d{3}[WC])\\)'\n")
        out.write("4. Modify the parse_genbank_files function to use these patterns\n")

    return found_genes, gene_matches, inference_formats, note_formats

def main():
    args = parse_arguments()
    gene_mapping, sc_gene_ids = load_gene_mapping(args.gene_mapping)
    found_genes, gene_matches, inference_formats, note_formats = analyze_genbank_fields(
        args.genbank_dir, sc_gene_ids, args.output)
    
    print(f"Found {len(found_genes)}/{len(sc_gene_ids)} SC gene IDs")
    print(f"Detailed report written to {args.output}")

if __name__ == "__main__":
    main()