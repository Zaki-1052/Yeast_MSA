#!/usr/bin/env python3
"""
verify_gene_coordinates.py - Verify gene annotations by comparing mapping file with GenBank files

This script compares the gene coordinates in our mapping file with the original GenBank annotations
to ensure we're targeting the correct regions in our analysis.

Usage:
    python verify_gene_coordinates.py --gene_mapping <gene_mapping_file> --genbank_dir <genbank_directory> --output_dir <output_directory>
"""
import re
import os
import sys
import argparse
import csv
from collections import defaultdict
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Verify gene coordinates by comparing mapping with GenBank files')
    parser.add_argument('--gene_mapping', required=True, help='TSV file mapping genes of interest')
    parser.add_argument('--genbank_dir', required=True, help='Directory containing GenBank annotation files')
    parser.add_argument('--output_dir', required=True, help='Directory for output files')
    return parser.parse_args()

def load_gene_mapping(mapping_file):
    """
    Load the gene mapping information for the target genes.
    
    Returns:
        dict: Dictionary mapping SC gene ID to gene information
        dict: Dictionary mapping W303 gene ID to gene information
        dict: Dictionary mapping scaffold to list of genes on that scaffold
    """
    sc_genes = {}
    w303_genes = {}
    scaffold_genes = defaultdict(list)
    
    with open(mapping_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            # Store gene information
            sc_gene_id = row['sc_gene_id']
            w303_gene_id = row['w303_gene_id']
            
            gene_info = {
                'sc_gene_id': sc_gene_id,
                'erg_name': row['erg_name'],
                'w303_gene_id': w303_gene_id,
                'locus_tag': row['locus_tag'],
                'w303_scaffold': row['w303_scaffold'],
                'start': int(row['start']),
                'end': int(row['end']),
                'strand': row['strand'],
                'product': row['product']
            }
            
            sc_genes[sc_gene_id] = gene_info
            w303_genes[w303_gene_id] = gene_info
            scaffold_genes[row['w303_scaffold']].append(w303_gene_id)
    
    return sc_genes, w303_genes, scaffold_genes

def parse_genbank_files(genbank_dir, sc_genes):
    """
    Parse GenBank files to extract gene information with improved SC gene ID extraction.
    
    Returns:
        dict: Dictionary mapping SC gene ID to GenBank gene information
        dict: Dictionary mapping W303 gene ID to GenBank gene information
        dict: Dictionary mapping scaffold to info about the scaffold
    """
    genbank_sc_genes = {}
    genbank_w303_genes = {}
    scaffolds_info = {}
    
    # Process each GenBank file
    genbank_files = [f for f in os.listdir(genbank_dir) if f.endswith('.genbank') or f.endswith('.gb')]
    print(f"Found {len(genbank_files)} GenBank files")
    
    sc_gene_ids = set(sc_genes.keys())
    target_w303_ids = set()
    for gene_info in sc_genes.values():
        target_w303_ids.add(gene_info['w303_gene_id'])
    
    for gb_file in sorted(genbank_files):
        file_path = os.path.join(genbank_dir, gb_file)
        print(f"Processing {gb_file}...")
        
        # Parse the GenBank file
        for record in SeqIO.parse(file_path, "genbank"):
            # Use record.name for scaffold ID (fixed in previous step)
            scaffold_id = record.name
            
            # Store scaffold information
            scaffolds_info[scaffold_id] = {
                'length': len(record.seq),
                'description': record.description
            }
            
            # Process each feature
            for feature in record.features:
                # We're looking for CDS features with information about our target genes
                if feature.type == "CDS":
                    # Extract gene information
                    gene_id = None
                    w303_id = None
                    locus_tag = None
                    sc_id = None
                    product = "Unknown"
                    
                    # Extract locus tag
                    if 'locus_tag' in feature.qualifiers:
                        locus_tag = feature.qualifiers['locus_tag'][0]
                        w303_id = locus_tag.replace('W303_0', 'W3030')
                    
                    # Extract gene ID
                    if 'gene' in feature.qualifiers:
                        gene_id = feature.qualifiers['gene'][0]
                    
                    # NEW: Check inference field for SC gene ID using regex
                    inferences = feature.qualifiers.get('inference', [])
                    for inference in inferences:
                        sc_match = re.search(r'Y[A-Z]{2}\d{3}[WC](-A)?', inference)
                        if sc_match:
                            potential_sc_id = sc_match.group(0)
                            # Remove -A suffix if present (e.g., YHR072W-A)
                            base_sc_id = potential_sc_id.split('-')[0]
                            if base_sc_id in sc_gene_ids:
                                sc_id = base_sc_id
                                break
                    
                    # Check note field if SC gene ID not found in inference
                    if not sc_id:
                        notes = feature.qualifiers.get('note', [])
                        for note in notes:
                            # Look for pattern: similar to Saccharomyces cerevisiae GENE (YXX123W)
                            sc_match = re.search(r'similar to Saccharomyces cerevisiae\s+\S+\s+\(([A-Z]{3}\d{3}[WC])\)', note)
                            if sc_match:
                                potential_sc_id = sc_match.group(1)
                                if potential_sc_id in sc_gene_ids:
                                    sc_id = potential_sc_id
                                    break
                            
                            # Fallback: direct search for SC gene ID pattern
                            if not sc_id:
                                sc_match = re.search(r'Y[A-Z]{2}\d{3}[WC]', note)
                                if sc_match:
                                    potential_sc_id = sc_match.group(0)
                                    if potential_sc_id in sc_gene_ids:
                                        sc_id = potential_sc_id
                                        break
                    
                    # Extract product
                    if 'product' in feature.qualifiers:
                        product = feature.qualifiers['product'][0]
                    
                    # Skip if we can't identify the gene
                    if not (sc_id or w303_id in target_w303_ids):
                        continue
                    
                    # Get coordinates
                    if feature.location:
                        start = int(feature.location.start)
                        end = int(feature.location.end)
                        strand = '+' if feature.location.strand == 1 else '-'
                    else:
                        continue
                    
                    # Create a record for this gene
                    gene_record = {
                        'sc_gene_id': sc_id,
                        'w303_gene_id': w303_id,
                        'locus_tag': locus_tag,
                        'scaffold': scaffold_id,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'product': product,
                        'source': 'genbank'
                    }
                    
                    # Store by both IDs if available
                    if sc_id:
                        genbank_sc_genes[sc_id] = gene_record
                    
                    if w303_id and w303_id in target_w303_ids:
                        genbank_w303_genes[w303_id] = gene_record
    
    print(f"Found {len(genbank_sc_genes)} genes by SC ID and {len(genbank_w303_genes)} genes by W303 ID in GenBank files")
    return genbank_sc_genes, genbank_w303_genes, scaffolds_info

def compare_coordinates(sc_genes, w303_genes, genbank_sc_genes, genbank_w303_genes):
    """
    Compare gene coordinates between mapping file and GenBank records.
    
    Returns:
        list: List of discrepancies
    """
    discrepancies = []
    
    # Check each gene in our mapping
    for sc_id, gene_info in sc_genes.items():
        w303_id = gene_info['w303_gene_id']
        
        mapping_record = {
            'sc_gene_id': sc_id,
            'w303_gene_id': w303_id,
            'scaffold': gene_info['w303_scaffold'],
            'start': gene_info['start'],
            'end': gene_info['end'],
            'strand': gene_info['strand'],
            'product': gene_info['product'],
            'source': 'mapping'
        }
        
        # Check if gene exists in GenBank records
        genbank_record_sc = genbank_sc_genes.get(sc_id)
        genbank_record_w303 = genbank_w303_genes.get(w303_id)
        
        # If found by SC ID, compare coordinates
        if genbank_record_sc:
            if (genbank_record_sc['start'] != mapping_record['start'] or
                genbank_record_sc['end'] != mapping_record['end'] or
                genbank_record_sc['strand'] != mapping_record['strand'] or
                genbank_record_sc['scaffold'] != mapping_record['scaffold']):
                
                discrepancies.append({
                    'gene_id': sc_id,
                    'w303_id': w303_id,
                    'mapping': mapping_record,
                    'genbank': genbank_record_sc
                })
        
        # If found by W303 ID but not SC ID, compare coordinates
        elif genbank_record_w303:
            if (genbank_record_w303['start'] != mapping_record['start'] or
                genbank_record_w303['end'] != mapping_record['end'] or
                genbank_record_w303['strand'] != mapping_record['strand'] or
                genbank_record_w303['scaffold'] != mapping_record['scaffold']):
                
                discrepancies.append({
                    'gene_id': sc_id,
                    'w303_id': w303_id,
                    'mapping': mapping_record,
                    'genbank': genbank_record_w303
                })
        
        # If not found in GenBank records
        else:
            discrepancies.append({
                'gene_id': sc_id,
                'w303_id': w303_id,
                'mapping': mapping_record,
                'genbank': None,
                'issue': 'Gene not found in GenBank records'
            })
    
    return discrepancies

def generate_gene_visualization(sc_genes, w303_genes, genbank_sc_genes, genbank_w303_genes, output_dir, scaffolds_info):
    """Generate visualization of gene locations for verification."""
    # Create output directory for visualizations
    viz_dir = os.path.join(output_dir, 'visualizations')
    os.makedirs(viz_dir, exist_ok=True)
    
    # Group genes by scaffold
    scaffold_genes_mapping = defaultdict(list)
    scaffold_genes_genbank = defaultdict(list)
    
    for gene_info in w303_genes.values():
        scaffold = gene_info['w303_scaffold']
        scaffold_genes_mapping[scaffold].append(gene_info)
    
    for gene_info in genbank_w303_genes.values():
        scaffold = gene_info['scaffold']
        scaffold_genes_genbank[scaffold].append(gene_info)
    
    # For each scaffold with our genes, create a visualization
    for scaffold in set(scaffold_genes_mapping.keys()) | set(scaffold_genes_genbank.keys()):
        # Get scaffold length
        scaffold_length = scaffolds_info.get(scaffold, {}).get('length', 0)
        if not scaffold_length:
            continue
        
        # Create figure
        plt.figure(figsize=(15, 8))
        plt.title(f"Gene Locations on {scaffold}")
        plt.xlabel("Position (bp)")
        plt.yticks([])
        
        # Set axis limits
        plt.xlim(0, scaffold_length)
        
        # Draw chromosome
        plt.axhline(y=1, xmin=0, xmax=1, color='black', linewidth=2)
        
        # Track for y-position offset to avoid overlaps
        y_positions = {'mapping': 1.2, 'genbank': 0.8}
        
        # Plot genes from mapping
        mapping_genes = scaffold_genes_mapping.get(scaffold, [])
        for gene_info in mapping_genes:
            gene_id = gene_info['sc_gene_id']
            start = gene_info['start']
            end = gene_info['end']
            strand = gene_info['strand']
            
            # Create arrow for gene
            arrow_length = max(end - start, scaffold_length * 0.01)  # Ensure visibility
            
            if strand == '+':
                plt.arrow(start, y_positions['mapping'], arrow_length, 0, 
                          head_width=0.05, head_length=scaffold_length * 0.01, 
                          fc='blue', ec='blue', length_includes_head=True)
            else:
                plt.arrow(end, y_positions['mapping'], -arrow_length, 0, 
                          head_width=0.05, head_length=scaffold_length * 0.01, 
                          fc='blue', ec='blue', length_includes_head=True)
            
            # Add gene label
            plt.text((start + end) / 2, y_positions['mapping'] + 0.1, 
                     f"{gene_id} ({gene_info['w303_gene_id']})", 
                     ha='center', va='bottom', fontsize=8)
            
            # Update y-position for next gene
            y_positions['mapping'] += 0.25
        
        # Plot genes from GenBank
        genbank_genes = scaffold_genes_genbank.get(scaffold, [])
        for gene_info in genbank_genes:
            gene_id = gene_info['sc_gene_id'] or gene_info['w303_gene_id']
            start = gene_info['start']
            end = gene_info['end']
            strand = gene_info['strand']
            
            # Create arrow for gene
            arrow_length = max(end - start, scaffold_length * 0.01)  # Ensure visibility
            
            if strand == '+':
                plt.arrow(start, y_positions['genbank'], arrow_length, 0, 
                          head_width=0.05, head_length=scaffold_length * 0.01, 
                          fc='red', ec='red', length_includes_head=True)
            else:
                plt.arrow(end, y_positions['genbank'], -arrow_length, 0, 
                          head_width=0.05, head_length=scaffold_length * 0.01, 
                          fc='red', ec='red', length_includes_head=True)
            
            # Add gene label
            plt.text((start + end) / 2, y_positions['genbank'] - 0.1, 
                     f"{gene_id}", 
                     ha='center', va='top', fontsize=8)
            
            # Update y-position for next gene
            y_positions['genbank'] -= 0.25
        
        # Add legend
        blue_patch = mpatches.Patch(color='blue', label='Mapping File')
        red_patch = mpatches.Patch(color='red', label='GenBank')
        plt.legend(handles=[blue_patch, red_patch], loc='upper right')
        
        # Save figure
        plt.savefig(os.path.join(viz_dir, f"{scaffold}_genes.png"), dpi=300, bbox_inches='tight')
        plt.close()

def generate_detailed_report(sc_genes, w303_genes, genbank_sc_genes, genbank_w303_genes, discrepancies, output_dir):
    """Generate detailed report on gene coordinates."""
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Write summary report
    with open(os.path.join(output_dir, 'gene_verification_summary.txt'), 'w') as f:
        f.write("Gene Coordinate Verification Summary\n")
        f.write("===================================\n\n")
        
        f.write(f"Total genes in mapping: {len(sc_genes)}\n")
        f.write(f"Genes found in GenBank by SC ID: {len(genbank_sc_genes)}\n")
        f.write(f"Genes found in GenBank by W303 ID: {len(genbank_w303_genes)}\n")
        f.write(f"Discrepancies found: {len(discrepancies)}\n\n")
        
        if discrepancies:
            f.write("Discrepancies:\n")
            for i, disc in enumerate(discrepancies, 1):
                gene_id = disc['gene_id']
                w303_id = disc['w303_id']
                mapping = disc['mapping']
                genbank = disc.get('genbank')
                
                f.write(f"{i}. Gene: {gene_id} ({w303_id})\n")
                
                if genbank:
                    f.write(f"   Issue: Coordinate mismatch\n")
                    f.write(f"   Mapping: {mapping['scaffold']}:{mapping['start']}-{mapping['end']} ({mapping['strand']})\n")
                    f.write(f"   GenBank: {genbank['scaffold']}:{genbank['start']}-{genbank['end']} ({genbank['strand']})\n")
                else:
                    f.write(f"   Issue: {disc.get('issue', 'Unknown issue')}\n")
                
                f.write("\n")
    
    # Write detailed gene information
    with open(os.path.join(output_dir, 'gene_details.tsv'), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow([
            'SC_Gene_ID', 'ERG_Name', 'W303_Gene_ID', 'Scaffold', 
            'Start_Mapping', 'End_Mapping', 'Strand_Mapping',
            'Start_GenBank', 'End_GenBank', 'Strand_GenBank',
            'Found_in_GenBank', 'Coordinate_Match'
        ])
        
        for sc_id, gene_info in sorted(sc_genes.items()):
            w303_id = gene_info['w303_gene_id']
            
            # Get GenBank info if available
            genbank_info = genbank_sc_genes.get(sc_id) or genbank_w303_genes.get(w303_id)
            
            writer.writerow([
                sc_id,
                gene_info.get('erg_name', ''),
                w303_id,
                gene_info['w303_scaffold'],
                gene_info['start'],
                gene_info['end'],
                gene_info['strand'],
                genbank_info['start'] if genbank_info else '',
                genbank_info['end'] if genbank_info else '',
                genbank_info['strand'] if genbank_info else '',
                'Yes' if genbank_info else 'No',
                'Yes' if genbank_info and (
                    genbank_info['start'] == gene_info['start'] and
                    genbank_info['end'] == gene_info['end'] and
                    genbank_info['strand'] == gene_info['strand']
                ) else 'No' if genbank_info else 'N/A'
            ])

def main():
    """Main function to verify gene coordinates."""
    args = parse_arguments()
    
    # Load gene mapping
    print(f"Loading gene mapping from {args.gene_mapping}")
    sc_genes, w303_genes, scaffold_genes = load_gene_mapping(args.gene_mapping)
    print(f"Loaded information for {len(sc_genes)} target genes")
    
    # Parse GenBank files
    print(f"Parsing GenBank files from {args.genbank_dir}")
    genbank_sc_genes, genbank_w303_genes, scaffolds_info = parse_genbank_files(args.genbank_dir, sc_genes)
    
    # Compare coordinates
    print("Comparing gene coordinates...")
    discrepancies = compare_coordinates(sc_genes, w303_genes, genbank_sc_genes, genbank_w303_genes)
    
    # Generate reports
    print("Generating reports...")
    generate_detailed_report(sc_genes, w303_genes, genbank_sc_genes, genbank_w303_genes, discrepancies, args.output_dir)
    
    # Generate visualizations
    print("Generating visualizations...")
    generate_gene_visualization(sc_genes, w303_genes, genbank_sc_genes, genbank_w303_genes, args.output_dir, scaffolds_info)
    
    print(f"Verification complete. Results in {args.output_dir}")
    print(f"Found {len(discrepancies)} discrepancies")
    
    # Print immediate summary
    perfect_matches = sum(1 for sc_id in sc_genes 
                         if (genbank_sc_genes.get(sc_id) or genbank_w303_genes.get(sc_genes[sc_id]['w303_gene_id'])) 
                         and sc_id not in [d['gene_id'] for d in discrepancies])
    
    print(f"Perfect matches: {perfect_matches}/{len(sc_genes)}")
    print(f"Detailed report written to {os.path.join(args.output_dir, 'gene_verification_summary.txt')}")

if __name__ == "__main__":
    main()