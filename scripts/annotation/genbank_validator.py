#!/usr/bin/env python3
"""
genbank_validator.py - Validate GenBank annotation files for consistency
"""

import os
import sys
from collections import Counter, defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def validate_genbank_files(directory):
    """
    Validate all GenBank files in the specified directory for consistency.
    
    Args:
        directory (str): Path to directory containing GenBank files
    
    Returns:
        dict: Summary statistics and validation results
    """
    # Initialize counters and containers
    stats = {
        'total_files': 0,
        'valid_files': 0,
        'invalid_files': 0,
        'total_genes': 0,
        'total_cds': 0,
        'total_trna': 0,
        'scaffolds': [],
        'original_ids': [],
        'issues': [],
        'gene_counts': [],
        'scaffold_lengths': [],
        'gene_names': set(),
        'scaffold_to_jriu': {},
    }
    
    # Process each GenBank file
    for filename in sorted(os.listdir(directory)):
        if not filename.endswith(('.genbank', '.gb', '.gbk')):
            continue
            
        filepath = os.path.join(directory, filename)
        stats['total_files'] += 1
        
        try:
            # Parse the GenBank file
            record = SeqIO.read(filepath, "genbank")
            
            # Extract basic information
            scaffold_id = record.id
            stats['scaffolds'].append(scaffold_id)
            stats['scaffold_lengths'].append(len(record.seq))
            
            # Look for original ID in the source feature notes
            original_id = None
            for feature in record.features:
                if feature.type == "source":
                    for note in feature.qualifiers.get('note', []):
                        if note.startswith('JRIU'):
                            original_id = note
                            stats['original_ids'].append(original_id)
                            stats['scaffold_to_jriu'][scaffold_id] = original_id
            
            # Count features
            gene_count = 0
            cds_count = 0
            trna_count = 0
            
            for feature in record.features:
                if feature.type == "gene":
                    gene_count += 1
                    # Extract gene name if available
                    if 'gene' in feature.qualifiers:
                        stats['gene_names'].add(feature.qualifiers['gene'][0])
                elif feature.type == "CDS":
                    cds_count += 1
                elif feature.type == "tRNA":
                    trna_count += 1
            
            stats['total_genes'] += gene_count
            stats['total_cds'] += cds_count
            stats['total_trna'] += trna_count
            stats['gene_counts'].append(gene_count)
            
            # Validation successful
            stats['valid_files'] += 1
            
        except Exception as e:
            # Record validation error
            stats['invalid_files'] += 1
            stats['issues'].append(f"Error in {filename}: {str(e)}")
    
    # Calculate summary statistics
    if stats['scaffold_lengths']:
        stats['min_length'] = min(stats['scaffold_lengths'])
        stats['max_length'] = max(stats['scaffold_lengths'])
        stats['total_length'] = sum(stats['scaffold_lengths'])
    
    if stats['gene_counts']:
        stats['min_genes'] = min(stats['gene_counts'])
        stats['max_genes'] = max(stats['gene_counts'])
        stats['mean_genes'] = sum(stats['gene_counts']) / len(stats['gene_counts'])
    
    return stats

def print_validation_report(stats):
    """Print a formatted report of validation results"""
    print("\n=== GenBank Validation Report ===\n")
    
    print(f"Total files processed: {stats['total_files']}")
    print(f"Valid files: {stats['valid_files']}")
    print(f"Invalid files: {stats['invalid_files']}")
    
    if stats['invalid_files'] > 0:
        print("\nIssues found:")
        for issue in stats['issues']:
            print(f"  - {issue}")
    
    print("\nAnnotation Statistics:")
    print(f"  Total scaffolds: {len(stats['scaffolds'])}")
    print(f"  Total genome length: {stats.get('total_length', 'N/A')} bp")
    print(f"  Scaffold length range: {stats.get('min_length', 'N/A')} - {stats.get('max_length', 'N/A')} bp")
    print(f"  Total genes: {stats['total_genes']}")
    print(f"  Total CDS features: {stats['total_cds']}")
    print(f"  Total tRNA features: {stats['total_trna']}")
    print(f"  Genes per scaffold range: {stats.get('min_genes', 'N/A')} - {stats.get('max_genes', 'N/A')}")
    print(f"  Unique gene names: {len(stats['gene_names'])}")
    
    print("\nScaffold mapping:")
    print(f"  Scaffolds with JRIU mapping: {len(stats['scaffold_to_jriu'])}")
    print(f"  Scaffolds missing mapping: {len(stats['scaffolds']) - len(stats['scaffold_to_jriu'])}")
    
    # List a few example mappings
    if stats['scaffold_to_jriu']:
        print("\nExample scaffold mappings:")
        count = 0
        for scaffold, jriu in list(stats['scaffold_to_jriu'].items())[:5]:
            print(f"  {scaffold} → {jriu}")
            count += 1
        if count < len(stats['scaffold_to_jriu']):
            print(f"  ... and {len(stats['scaffold_to_jriu']) - count} more")

def create_scaffold_mapping_file(stats, output_file):
    """Create a TSV file with scaffold to original ID mapping"""
    with open(output_file, 'w') as f:
        f.write("scaffold_id\toriginal_id\tscaffold_length\tgene_count\n")
        
        # Get gene counts by scaffold
        scaffold_to_gene_count = defaultdict(int)
        for i, scaffold in enumerate(stats['scaffolds']):
            if i < len(stats['gene_counts']):
                scaffold_to_gene_count[scaffold] = stats['gene_counts'][i]
        
        # Get scaffold lengths
        scaffold_to_length = defaultdict(int)
        for i, scaffold in enumerate(stats['scaffolds']):
            if i < len(stats['scaffold_lengths']):
                scaffold_to_length[scaffold] = stats['scaffold_lengths'][i]
        
        # Write mapping data
        for scaffold in sorted(stats['scaffolds']):
            original_id = stats['scaffold_to_jriu'].get(scaffold, "unknown")
            length = scaffold_to_length.get(scaffold, 0)
            gene_count = scaffold_to_gene_count.get(scaffold, 0)
            f.write(f"{scaffold}\t{original_id}\t{length}\t{gene_count}\n")

def main():
    if len(sys.argv) < 2:
        print("Usage: python genbank_validator.py <genbank_directory> [output_mapping_file]")
        sys.exit(1)
    
    genbank_dir = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else "scaffold_mapping.tsv"
    
    if not os.path.isdir(genbank_dir):
        print(f"Error: {genbank_dir} is not a directory")
        sys.exit(1)
    
    print(f"Validating GenBank files in {genbank_dir}...")
    stats = validate_genbank_files(genbank_dir)
    print_validation_report(stats)
    
    print(f"\nCreating scaffold mapping file: {output_file}")
    create_scaffold_mapping_file(stats, output_file)
    print("Done!")

if __name__ == "__main__":
    main()