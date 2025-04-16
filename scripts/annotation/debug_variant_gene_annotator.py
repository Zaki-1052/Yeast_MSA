#!/usr/bin/env python3
"""
debug_variant_gene_annotator.py - Debug version to identify annotation issues
"""

import os
import sys
import csv
import argparse
from collections import defaultdict, Counter
import gzip
from Bio import SeqIO

def debug_variant_gene_matching(vcf_file, scaffold_mapping_file, genes_of_interest_file, output_file):
    """
    Debug function to trace the variant-gene matching process
    """
    # Load scaffold mapping
    jriu_to_scaffold = {}
    scaffold_to_jriu = {}
    
    print("Loading scaffold mapping...")
    with open(scaffold_mapping_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            scaffold_id = row['scaffold_id']
            jriu_id = row['jriu_id']
            
            if jriu_id and jriu_id != 'unknown':
                jriu_to_scaffold[jriu_id] = scaffold_id
                scaffold_to_jriu[scaffold_id] = jriu_id
    
    print(f"Loaded {len(jriu_to_scaffold)} JRIU to scaffold mappings")
    
    # Load genes of interest
    genes_of_interest = {}
    
    print("Loading genes of interest...")
    with open(genes_of_interest_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            sc_gene_id = row['sc_gene_id']
            w303_gene_id = row['w303_gene_id']
            scaffold_id = row['scaffold_id']
            jriu_id = row['jriu_id']
            start = int(row['start'])
            end = int(row['end'])
            
            genes_of_interest[w303_gene_id] = {
                'sc_gene_id': sc_gene_id,
                'scaffold_id': scaffold_id,
                'jriu_id': jriu_id,
                'start': start,
                'end': end
            }
    
    print(f"Loaded {len(genes_of_interest)} genes of interest")
    
    # Show a few genes of interest for debugging
    print("\nSample genes of interest:")
    for i, (gene_id, info) in enumerate(list(genes_of_interest.items())[:3]):
        print(f"  {gene_id}: {info['scaffold_id']} ({info['jriu_id']}) at {info['start']}-{info['end']}")
    
    # Process VCF file
    is_gzipped = vcf_file.endswith('.gz')
    opener = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'
    
    # Extract sample name
    sample_name = os.path.basename(vcf_file).split('.')[0]
    
    print(f"\nProcessing VCF file: {vcf_file}")
    
    # Debug counters
    total_variants = 0
    mapped_variants = 0
    goi_jriu_matches = 0
    goi_scaffold_matches = 0
    potential_overlaps = 0
    confirmed_overlaps = 0
    
    # Store matching variants
    matching_variants = []
    
    # Create inverse lookup from JRIU to gene
    jriu_to_genes = defaultdict(list)
    for gene_id, info in genes_of_interest.items():
        jriu_id = info['jriu_id']
        if jriu_id:
            jriu_to_genes[jriu_id].append(gene_id)
    
    print(f"Created lookup for {len(jriu_to_genes)} JRIUs containing genes of interest")
    
    # Process VCF
    with opener(vcf_file, mode) as f:
        # Skip header lines
        for line in f:
            if line.startswith('#'):
                continue
            
            total_variants += 1
            fields = line.strip().split('\t')
            
            if len(fields) < 8:
                continue
            
            # Extract variant information
            chrom = fields[0]  # JRIU ID
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            
            # Check if this JRIU is in our mapping
            scaffold_id = jriu_to_scaffold.get(chrom)
            if scaffold_id:
                mapped_variants += 1
            
            # Check if this JRIU has any genes of interest
            if chrom in jriu_to_genes:
                goi_jriu_matches += 1
                print(f"\nFound variant on JRIU with genes of interest: {chrom} at position {pos}")
                
                # Check each gene of interest on this JRIU
                for gene_id in jriu_to_genes[chrom]:
                    gene_info = genes_of_interest[gene_id]
                    gene_start = gene_info['start']
                    gene_end = gene_info['end']
                    
                    print(f"  Checking against {gene_id}: {gene_start}-{gene_end}")
                    
                    # Check if position falls within gene boundaries
                    if gene_start <= pos <= gene_end:
                        confirmed_overlaps += 1
                        print(f"  MATCH! Position {pos} is within gene {gene_id} ({gene_start}-{gene_end})")
                        
                        matching_variants.append({
                            'chrom': chrom,
                            'pos': pos,
                            'ref': ref,
                            'alt': alt,
                            'gene_id': gene_id,
                            'rel_position': pos - gene_start
                        })
            
            # Check by scaffold ID too
            elif scaffold_id and scaffold_id in [info['scaffold_id'] for info in genes_of_interest.values()]:
                goi_scaffold_matches += 1
                print(f"\nFound variant on scaffold with genes of interest: {scaffold_id} (from {chrom}) at position {pos}")
                
                # Check each gene on this scaffold
                for gene_id, gene_info in genes_of_interest.items():
                    if gene_info['scaffold_id'] == scaffold_id:
                        gene_start = gene_info['start']
                        gene_end = gene_info['end']
                        
                        print(f"  Checking against {gene_id}: {gene_start}-{gene_end}")
                        
                        # Check if position falls within gene boundaries
                        if gene_start <= pos <= gene_end:
                            confirmed_overlaps += 1
                            print(f"  MATCH! Position {pos} is within gene {gene_id} ({gene_start}-{gene_end})")
                            
                            matching_variants.append({
                                'chrom': chrom,
                                'pos': pos,
                                'ref': ref,
                                'alt': alt,
                                'gene_id': gene_id,
                                'rel_position': pos - gene_start
                            })
            
            # Early exit after processing a reasonable number of variants for debugging
            if total_variants >= 10000:
                print("\nProcessed 10,000 variants for debugging, stopping early...")
                break
    
    # Print summary
    print("\nDebug Summary:")
    print(f"Total variants processed: {total_variants}")
    print(f"Variants with mapped scaffold: {mapped_variants} ({mapped_variants/total_variants*100:.2f}%)")
    print(f"Variants on JRIU with genes of interest: {goi_jriu_matches}")
    print(f"Variants on scaffold with genes of interest: {goi_scaffold_matches}")
    print(f"Confirmed overlapping variants: {confirmed_overlaps}")
    
    # Write matching variants to output file
    if matching_variants:
        with open(output_file, 'w') as f:
            f.write("CHROM\tPOS\tREF\tALT\tGENE_ID\tREL_POSITION\n")
            for var in matching_variants:
                f.write(f"{var['chrom']}\t{var['pos']}\t{var['ref']}\t{var['alt']}\t{var['gene_id']}\t{var['rel_position']}\n")
        
        print(f"\nWrote {len(matching_variants)} matching variants to {output_file}")
    else:
        print("\nNo matching variants found.")
    
    # Suggest fixes
    print("\nPossible issues and fixes:")
    
    if mapped_variants == 0:
        print("- CRITICAL: No variants could be mapped to scaffolds. Check the scaffold mapping file format.")
    elif mapped_variants < total_variants * 0.5:
        print("- WARNING: Less than 50% of variants could be mapped to scaffolds.")
    
    if goi_jriu_matches == 0 and goi_scaffold_matches == 0:
        print("- CRITICAL: No variants found on scaffolds or JRIUs containing genes of interest.")
        print("  Suggestion: Verify that the genes of interest are on scaffolds present in the VCF.")
    
    if goi_jriu_matches > 0 and confirmed_overlaps == 0:
        print("- WARNING: Found variants on JRIUs with genes of interest, but no position overlaps.")
        print("  Suggestion: Check for coordinate system mismatches (0-based vs 1-based).")
    
    # Print all JRIUs from VCF that contain genes of interest
    jriu_with_goi = set(jriu_to_genes.keys())
    print(f"\nThere are {len(jriu_with_goi)} JRIUs containing genes of interest:")
    for jriu in sorted(list(jriu_with_goi)[:10]):
        print(f"  {jriu} -> {', '.join(jriu_to_genes[jriu])}")
    if len(jriu_with_goi) > 10:
        print(f"  ... and {len(jriu_with_goi) - 10} more")

def main():
    parser = argparse.ArgumentParser(description='Debug variant gene matching')
    parser.add_argument('--vcf', required=True, help='VCF file to analyze')
    parser.add_argument('--scaffold-mapping', required=True, help='Scaffold to JRIU mapping file')
    parser.add_argument('--genes-of-interest', required=True, help='Genes of interest file')
    parser.add_argument('--output', required=True, help='Output file for matching variants')
    
    args = parser.parse_args()
    debug_variant_gene_matching(args.vcf, args.scaffold_mapping, args.genes_of_interest, args.output)

if __name__ == '__main__':
    main()