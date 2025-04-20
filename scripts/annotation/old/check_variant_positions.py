#!/usr/bin/env python3
"""
check_variant_positions.py - Directly compare variant positions with gene positions
"""

import os
import sys
import csv
import argparse
import gzip
from Bio import SeqIO
from collections import defaultdict

def check_variant_positions(vcf_file, genbank_dir, genes_of_interest_file, output_file):
    """
    Directly check if variants in VCF file overlap with genes of interest by position
    """
    # Load genes of interest
    genes_by_jriu = defaultdict(list)
    all_genes = []
    
    with open(genes_of_interest_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene_id = row.get('w303_gene_id')
            jriu_id = row.get('jriu_id')
            
            if gene_id and jriu_id:
                # Store gene info
                gene_info = {
                    'gene_id': gene_id,
                    'sc_gene_id': row.get('sc_gene_id'),
                    'jriu_id': jriu_id,
                    'start': int(row.get('start', 0)),
                    'end': int(row.get('end', 0)),
                    'strand': row.get('strand', '+')
                }
                genes_by_jriu[jriu_id].append(gene_info)
                all_genes.append(gene_info)
    
    print(f"Loaded {len(all_genes)} genes of interest on {len(genes_by_jriu)} JRIU IDs")
    
    # Load GenBank files for JRIUs with genes of interest
    genbank_by_jriu = {}
    for jriu_id in genes_by_jriu:
        for filename in os.listdir(genbank_dir):
            if not filename.endswith(('.genbank', '.gb', '.gbk')):
                continue
            
            filepath = os.path.join(genbank_dir, filename)
            try:
                record = SeqIO.read(filepath, "genbank")
                for feature in record.features:
                    if feature.type == "source":
                        for note in feature.qualifiers.get('note', []):
                            if note == jriu_id:
                                genbank_by_jriu[jriu_id] = filepath
                                break
            except Exception as e:
                pass
    
    print(f"Found GenBank files for {len(genbank_by_jriu)}/{len(genes_by_jriu)} JRIU IDs")
    
    # Extract variants from VCF file
    is_gzipped = vcf_file.endswith('.gz')
    opener = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'
    
    variants_by_jriu = defaultdict(list)
    total_variants = 0
    
    with opener(vcf_file, mode) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            
            chrom = fields[0]  # JRIU ID
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            
            variants_by_jriu[chrom].append({
                'pos': pos,
                'ref': ref,
                'alt': alt
            })
            
            total_variants += 1
    
    print(f"Loaded {total_variants} variants on {len(variants_by_jriu)} JRIU IDs")
    
    # Find overlapping JRIUs
    overlapping_jrius = set(genes_by_jriu.keys()) & set(variants_by_jriu.keys())
    print(f"Found {len(overlapping_jrius)} JRIU IDs with both genes and variants")
    
    # Check each overlapping JRIU for position-based overlaps
    overlapping_variants = []
    jriu_results = {}
    
    for jriu_id in overlapping_jrius:
        genes = genes_by_jriu[jriu_id]
        variants = variants_by_jriu[jriu_id]
        
        # Sort variants by position
        variants.sort(key=lambda v: v['pos'])
        
        # Check for overlaps
        jriu_overlaps = []
        
        for gene in genes:
            gene_start = gene['start']
            gene_end = gene['end']
            
            # Count variants in different regions
            within_gene = 0
            upstream_gene = 0
            downstream_gene = 0
            
            for variant in variants:
                pos = variant['pos']
                
                # Direct position check
                if gene_start <= pos <= gene_end:
                    within_gene += 1
                    # Remember this overlapping variant
                    overlapping_variants.append({
                        'jriu_id': jriu_id,
                        'pos': pos,
                        'ref': variant['ref'],
                        'alt': variant['alt'],
                        'gene_id': gene['gene_id'],
                        'sc_gene_id': gene['sc_gene_id'],
                        'gene_start': gene_start,
                        'gene_end': gene_end,
                        'relative_pos': pos - gene_start
                    })
                elif pos < gene_start:
                    upstream_gene += 1
                else:  # pos > gene_end
                    downstream_gene += 1
            
            jriu_overlaps.append({
                'gene_id': gene['gene_id'],
                'sc_gene_id': gene['sc_gene_id'],
                'gene_start': gene_start,
                'gene_end': gene_end,
                'within': within_gene,
                'upstream': upstream_gene,
                'downstream': downstream_gene
            })
        
        jriu_results[jriu_id] = {
            'genes': len(genes),
            'variants': len(variants),
            'variant_range': f"{variants[0]['pos']}-{variants[-1]['pos']}" if variants else "N/A",
            'overlaps': jriu_overlaps
        }
    
    # Also check GenBank files directly for genes
    print("\nChecking GenBank files directly for gene positions...\n")
    
    for jriu_id in overlapping_jrius:
        # Get GenBank file for this JRIU
        genbank_file = genbank_by_jriu.get(jriu_id)
        if not genbank_file:
            continue
        
        try:
            # Parse GenBank file
            record = SeqIO.read(genbank_file, "genbank")
            
            # Find all gene features
            gene_features = {}
            for feature in record.features:
                if feature.type == "gene":
                    gene_id = None
                    if 'gene' in feature.qualifiers:
                        gene_id = feature.qualifiers['gene'][0]
                    elif 'locus_tag' in feature.qualifiers:
                        gene_id = feature.qualifiers['locus_tag'][0]
                    
                    if gene_id:
                        gene_features[gene_id] = {
                            'start': int(feature.location.start),
                            'end': int(feature.location.end),
                            'strand': '+' if feature.location.strand == 1 else '-'
                        }
            
            # Compare with our gene info
            for gene in genes_by_jriu[jriu_id]:
                gene_id = gene['gene_id']
                gene_start = gene['start']
                gene_end = gene['end']
                
                if gene_id in gene_features:
                    gb_start = gene_features[gene_id]['start']
                    gb_end = gene_features[gene_id]['end']
                    
                    # Check if positions match
                    start_diff = gene_start - gb_start
                    end_diff = gene_end - gb_end
                    
                    print(f"Gene {gene_id} on {jriu_id}:")
                    print(f"  Our coords: {gene_start}-{gene_end}")
                    print(f"  GenBank coords: {gb_start}-{gb_end}")
                    print(f"  Difference: start={start_diff}, end={end_diff}")
                    
                    # If there's a consistent offset, check if variants would overlap with adjusted coordinates
                    if start_diff == end_diff and start_diff != 0:
                        adjusted_overlaps = 0
                        for variant in variants_by_jriu[jriu_id]:
                            pos = variant['pos']
                            if gb_start <= pos <= gb_end:
                                adjusted_overlaps += 1
                        
                        print(f"  Variants overlapping with GenBank coords: {adjusted_overlaps}")
                    
                    print("")
        except Exception as e:
            print(f"Error processing {genbank_file}: {str(e)}")
    
    # Print summary of results
    print("\n=== Results Summary ===\n")
    total_overlaps = 0
    
    for jriu_id, result in jriu_results.items():
        print(f"JRIU {jriu_id}:")
        print(f"  {result['genes']} genes, {result['variants']} variants")
        print(f"  Variant position range: {result['variant_range']}")
        
        for overlap in result['overlaps']:
            gene_id = overlap['gene_id']
            sc_gene_id = overlap['sc_gene_id']
            within = overlap['within']
            total_overlaps += within
            
            print(f"  Gene {gene_id} ({sc_gene_id}):")
            print(f"    Position: {overlap['gene_start']}-{overlap['gene_end']}")
            print(f"    Variants: {within} within, {overlap['upstream']} upstream, {overlap['downstream']} downstream")
    
    print(f"\nTotal overlapping variants found: {total_overlaps}")
    
    # Output to file
    with open(output_file, 'w') as f:
        f.write("=== Variant-Gene Position Check Results ===\n\n")
        f.write(f"VCF File: {vcf_file}\n")
        f.write(f"Total variants: {total_variants}\n")
        f.write(f"Total genes of interest: {len(all_genes)}\n")
        f.write(f"JRIUs with both genes and variants: {len(overlapping_jrius)}\n")
        f.write(f"Total overlapping variants found: {total_overlaps}\n\n")
        
        f.write("Gene-specific results:\n")
        for jriu_id, result in sorted(jriu_results.items()):
            f.write(f"JRIU {jriu_id}:\n")
            f.write(f"  Variants: {result['variants']} (Range: {result['variant_range']})\n")
            
            for overlap in result['overlaps']:
                f.write(f"  Gene {overlap['gene_id']} ({overlap['sc_gene_id']}):\n")
                f.write(f"    Position: {overlap['gene_start']}-{overlap['gene_end']}\n")
                f.write(f"    Overlaps: {overlap['within']} within, {overlap['upstream']} upstream, {overlap['downstream']} downstream\n")
            
            f.write("\n")
        
        if overlapping_variants:
            f.write("\nDetailed overlapping variants:\n")
            for i, ov in enumerate(overlapping_variants):
                f.write(f"{i+1}. {ov['jriu_id']}:{ov['pos']} {ov['ref']}>{ov['alt']} in gene {ov['gene_id']} ({ov['sc_gene_id']})\n")
                f.write(f"   Gene range: {ov['gene_start']}-{ov['gene_end']}, Relative position: {ov['relative_pos']}\n")
        else:
            f.write("\nNo overlapping variants found.\n")
            
        # If no overlaps were found, provide possible explanations
        if total_overlaps == 0:
            f.write("\nPossible explanations for no overlaps:\n")
            f.write("1. Coordinate system mismatch (0-based vs 1-based)\n")
            f.write("2. Different reference versions for VCF and gene annotations\n")
            f.write("3. Off-by-one error in position comparison\n")
            f.write("4. Issues in parsing gene positions from GenBank files\n")
    
    print(f"\nDetailed results written to {output_file}")
    return total_overlaps

def main():
    parser = argparse.ArgumentParser(description='Check variant positions against gene positions')
    parser.add_argument('--vcf', required=True, help='VCF file to analyze')
    parser.add_argument('--genbank-dir', required=True, help='Directory containing GenBank files')
    parser.add_argument('--genes-of-interest', required=True, help='Genes of interest file')
    parser.add_argument('--output', required=True, help='Output file for results')
    
    args = parser.parse_args()
    check_variant_positions(args.vcf, args.genbank_dir, args.genes_of_interest, args.output)

if __name__ == '__main__':
    main()