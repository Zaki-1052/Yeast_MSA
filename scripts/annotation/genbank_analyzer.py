#!/usr/bin/env python3
"""
genbank_analyzer.py - Extract scaffold mappings and gene information from GenBank files
"""

import os
import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def analyze_genbank_files(directory, genes_of_interest=None):
    """
    Analyze all GenBank files in the specified directory and extract mappings
    and gene information.
    
    Args:
        directory (str): Path to directory containing GenBank files
        genes_of_interest (list): List of gene names to specially track
        
    Returns:
        dict: Dictionary containing extracted data and statistics
    """
    # Initialize data containers
    data = {
        'total_files': 0,
        'valid_files': 0,
        'invalid_files': 0,
        'total_genes': 0,
        'total_cds': 0,
        'total_trna': 0,
        'scaffold_info': {},
        'gene_index': {},
        'genes_of_interest_data': {},
        'issues': [],
    }
    
    if genes_of_interest is None:
        genes_of_interest = []
    genes_of_interest = [g.upper() for g in genes_of_interest]
    
    # Process each GenBank file
    for filename in sorted(os.listdir(directory)):
        if not filename.endswith(('.genbank', '.gb', '.gbk')):
            continue
            
        filepath = os.path.join(directory, filename)
        data['total_files'] += 1
        
        try:
            # Parse the GenBank file
            record = SeqIO.read(filepath, "genbank")
            
            # Extract scaffold ID - use ID from record or from filename
            scaffold_id = record.id
            if scaffold_id == "unknown":
                # Try to get scaffold number from filename
                if "_scaffold_" in filename:
                    scaffold_num = filename.split("_scaffold_")[1].split(".")[0]
                    scaffold_id = f"scaffold_{scaffold_num}"
            
            # Initialize scaffold info
            data['scaffold_info'][scaffold_id] = {
                'length': len(record.seq),
                'jriu_id': None,
                'gene_count': 0,
                'cds_count': 0,
                'trna_count': 0,
                'genes': [],
                'filepath': filepath
            }
            
            # Extract JRIU ID from source feature notes
            for feature in record.features:
                if feature.type == "source":
                    for note in feature.qualifiers.get('note', []):
                        if note.startswith('JRIU'):
                            data['scaffold_info'][scaffold_id]['jriu_id'] = note
                        elif note.startswith('scaffold '):
                            scaffold_num = note.split('scaffold ')[1]
                            alt_id = f"scaffold_{scaffold_num}"
                            if scaffold_id == "unknown":
                                scaffold_id = alt_id
                                data['scaffold_info'][scaffold_id] = data['scaffold_info'].pop("unknown")
            
            # Process gene features
            gene_features = defaultdict(dict)
            for feature in record.features:
                if feature.type == "gene":
                    gene_id = None
                    if 'gene' in feature.qualifiers:
                        gene_id = feature.qualifiers['gene'][0]
                    elif 'locus_tag' in feature.qualifiers:
                        gene_id = feature.qualifiers['locus_tag'][0]
                    
                    if gene_id:
                        data['total_genes'] += 1
                        data['scaffold_info'][scaffold_id]['gene_count'] += 1
                        data['scaffold_info'][scaffold_id]['genes'].append(gene_id)
                        
                        # Store gene location info
                        location = feature.location
                        strand = location.strand
                        gene_features[gene_id]['location'] = (
                            int(location.start),
                            int(location.end), 
                            '+' if strand == 1 else '-'
                        )
                        gene_features[gene_id]['scaffold_id'] = scaffold_id
                        
                        # Add to gene index
                        data['gene_index'][gene_id] = {
                            'scaffold_id': scaffold_id,
                            'location': gene_features[gene_id]['location'],
                            'product': None,
                            'note': None
                        }
                        
                        # Check if this is a gene of interest
                        sc_gene_id = None
                        if 'inference' in feature.qualifiers:
                            for inf in feature.qualifiers['inference']:
                                if 'similar to' in inf and ':' in inf:
                                    parts = inf.split(':')
                                    if len(parts) > 1:
                                        sc_genes = parts[-1].split(',')
                                        for gene in sc_genes:
                                            # Remove any parentheses and get S. cerevisiae gene ID
                                            sc_gene = gene.strip().split('(')[-1].split(')')[0]
                                            if sc_gene.upper() in genes_of_interest:
                                                sc_gene_id = sc_gene.upper()
                        
                        # If note contains gene ID, check if it matches genes of interest
                        if 'note' in feature.qualifiers:
                            for note in feature.qualifiers['note']:
                                gene_features[gene_id]['note'] = note
                                data['gene_index'][gene_id]['note'] = note
                                
                                # Check for S. cerevisiae genes in note
                                if 'similar to Saccharomyces cerevisiae' in note:
                                    for goi in genes_of_interest:
                                        if goi in note:
                                            sc_gene_id = goi
                        
                        # If this is a gene of interest, store additional details
                        if sc_gene_id:
                            if sc_gene_id not in data['genes_of_interest_data']:
                                data['genes_of_interest_data'][sc_gene_id] = []
                            
                            data['genes_of_interest_data'][sc_gene_id].append({
                                'gene_id': gene_id,
                                'scaffold_id': scaffold_id,
                                'location': gene_features[gene_id]['location'],
                                'note': gene_features[gene_id].get('note')
                            })
                
                elif feature.type == "CDS":
                    data['total_cds'] += 1
                    data['scaffold_info'][scaffold_id]['cds_count'] += 1
                    
                    cds_gene_id = None
                    if 'gene' in feature.qualifiers:
                        cds_gene_id = feature.qualifiers['gene'][0]
                    elif 'locus_tag' in feature.qualifiers:
                        cds_gene_id = feature.qualifiers['locus_tag'][0]
                    
                    if cds_gene_id and cds_gene_id in gene_features:
                        if 'product' in feature.qualifiers:
                            product = feature.qualifiers['product'][0]
                            gene_features[cds_gene_id]['product'] = product
                            if cds_gene_id in data['gene_index']:
                                data['gene_index'][cds_gene_id]['product'] = product
                    
                elif feature.type == "tRNA":
                    data['total_trna'] += 1
                    data['scaffold_info'][scaffold_id]['trna_count'] += 1
                    
            data['valid_files'] += 1
            
        except Exception as e:
            # Record validation error
            data['invalid_files'] += 1
            data['issues'].append(f"Error in {filename}: {str(e)}")
    
    # Calculate summary statistics
    scaffold_lengths = [info['length'] for info in data['scaffold_info'].values()]
    if scaffold_lengths:
        data['min_length'] = min(scaffold_lengths)
        data['max_length'] = max(scaffold_lengths)
        data['total_length'] = sum(scaffold_lengths)
    
    gene_counts = [info['gene_count'] for info in data['scaffold_info'].values()]
    if gene_counts:
        data['min_genes'] = min(gene_counts)
        data['max_genes'] = max(gene_counts)
        data['avg_genes'] = sum(gene_counts) / len(gene_counts)
    
    # Count scaffolds with JRIU IDs
    jriu_count = sum(1 for info in data['scaffold_info'].values() if info.get('jriu_id'))
    data['scaffolds_with_jriu'] = jriu_count
    
    return data

def print_analysis_report(data):
    """Print a formatted report of analysis results"""
    print("\n=== GenBank Analysis Report ===\n")
    
    print(f"Total files processed: {data['total_files']}")
    print(f"Valid files: {data['valid_files']}")
    print(f"Invalid files: {data['invalid_files']}")
    
    if data['invalid_files'] > 0:
        print("\nIssues found:")
        for issue in data['issues']:
            print(f"  - {issue}")
    
    print("\nAnnotation Statistics:")
    print(f"  Total scaffolds: {len(data['scaffold_info'])}")
    print(f"  Total genome length: {data.get('total_length', 'N/A')} bp")
    print(f"  Scaffold length range: {data.get('min_length', 'N/A')} - {data.get('max_length', 'N/A')} bp")
    print(f"  Total genes: {data['total_genes']}")
    print(f"  Total CDS features: {data['total_cds']}")
    print(f"  Total tRNA features: {data['total_trna']}")
    print(f"  Genes per scaffold range: {data.get('min_genes', 'N/A')} - {data.get('max_genes', 'N/A')}")
    print(f"  Average genes per scaffold: {data.get('avg_genes', 'N/A'):.2f}")
    
    print("\nScaffold mapping:")
    print(f"  Scaffolds with JRIU mapping: {data['scaffolds_with_jriu']}")
    print(f"  Scaffolds missing mapping: {len(data['scaffold_info']) - data['scaffolds_with_jriu']}")
    
    # Print some example mappings
    if data['scaffold_info']:
        print("\nExample scaffold mappings:")
        count = 0
        for scaffold_id, info in sorted(data['scaffold_info'].items())[:5]:
            if info.get('jriu_id'):
                print(f"  {scaffold_id} → {info['jriu_id']}")
                count += 1
        if count < data['scaffolds_with_jriu']:
            print(f"  ... and {data['scaffolds_with_jriu'] - count} more")
    
    # Print gene index stats
    print(f"\nGene index entries: {len(data['gene_index'])}")
    
    # Print genes of interest data
    if data['genes_of_interest_data']:
        print("\nGenes of interest found:")
        for sc_gene_id, matches in sorted(data['genes_of_interest_data'].items()):
            print(f"  {sc_gene_id}: {len(matches)} matching gene(s)")
            for i, match in enumerate(matches):
                gene_id = match['gene_id']
                scaffold = match['scaffold_id']
                loc = match['location']
                print(f"    {i+1}. {gene_id} on {scaffold} at {loc[0]}-{loc[1]} ({loc[2]})")

def create_scaffold_mapping_file(data, output_file):
    """Create a TSV file with scaffold to JRIU mapping"""
    with open(output_file, 'w') as f:
        f.write("scaffold_id\tjriu_id\tscaffold_length\tgene_count\tcds_count\ttrna_count\n")
        
        # Write mapping data
        for scaffold_id, info in sorted(data['scaffold_info'].items()):
            jriu_id = info.get('jriu_id', 'unknown')
            length = info.get('length', 0)
            gene_count = info.get('gene_count', 0)
            cds_count = info.get('cds_count', 0)
            trna_count = info.get('trna_count', 0)
            f.write(f"{scaffold_id}\t{jriu_id}\t{length}\t{gene_count}\t{cds_count}\t{trna_count}\n")

def create_gene_index_file(data, output_file):
    """Create a TSV file with gene index information"""
    with open(output_file, 'w') as f:
        f.write("gene_id\tscaffold_id\tstart\tend\tstrand\tproduct\tnote\n")
        
        # Write gene data
        for gene_id, info in sorted(data['gene_index'].items()):
            scaffold_id = info.get('scaffold_id', 'unknown')
            location = info.get('location', (0, 0, '+'))
            start = location[0]
            end = location[1]
            strand = location[2]
            product = info.get('product', '')
            note = info.get('note', '')
            
            # Escape tabs and newlines in long text fields
            if product:
                product = product.replace('\t', ' ').replace('\n', ' ')
            if note:
                note = note.replace('\t', ' ').replace('\n', ' ')
                
            f.write(f"{gene_id}\t{scaffold_id}\t{start}\t{end}\t{strand}\t{product}\t{note}\n")

def create_genes_of_interest_file(data, output_file):
    """Create a TSV file with detailed information about genes of interest"""
    with open(output_file, 'w') as f:
        f.write("sc_gene_id\tw303_gene_id\tscaffold_id\tjriu_id\tstart\tend\tstrand\tproduct\tnote\n")
        
        # Write gene of interest data
        for sc_gene_id, matches in sorted(data['genes_of_interest_data'].items()):
            for match in matches:
                gene_id = match['gene_id']
                scaffold_id = match['scaffold_id']
                jriu_id = data['scaffold_info'][scaffold_id].get('jriu_id', 'unknown')
                location = match['location']
                start = location[0]
                end = location[1]
                strand = location[2]
                
                # Get product and note
                product = ''
                note = match.get('note', '')
                if gene_id in data['gene_index']:
                    product = data['gene_index'][gene_id].get('product', '')
                
                # Escape tabs and newlines in long text fields
                if product:
                    product = product.replace('\t', ' ').replace('\n', ' ')
                if note:
                    note = note.replace('\t', ' ').replace('\n', ' ')
                
                f.write(f"{sc_gene_id}\t{gene_id}\t{scaffold_id}\t{jriu_id}\t{start}\t{end}\t{strand}\t{product}\t{note}\n")

def main():
    if len(sys.argv) < 2:
        print("Usage: python genbank_analyzer.py <genbank_directory> [output_directory] [genes_of_interest_file]")
        sys.exit(1)
    
    genbank_dir = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "."
    genes_of_interest_file = sys.argv[3] if len(sys.argv) > 3 else None
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Get list of genes of interest
    genes_of_interest = []
    if genes_of_interest_file and os.path.isfile(genes_of_interest_file):
        with open(genes_of_interest_file, 'r') as f:
            genes_of_interest = [line.strip() for line in f if line.strip()]
    else:
        # Default list of genes of interest from the project
        genes_of_interest = [
            "YHR190W", "YGR175C", "YHR072W", "YHR007C", "YNL280C", 
            "YGR060W", "YML008C", "YMR202W", "YLR056W", "YMR015C", "YGL012W"
        ]
    
    if not os.path.isdir(genbank_dir):
        print(f"Error: {genbank_dir} is not a directory")
        sys.exit(1)
    
    print(f"Analyzing GenBank files in {genbank_dir}...")
    print(f"Looking for the following genes of interest: {', '.join(genes_of_interest)}")
    
    data = analyze_genbank_files(genbank_dir, genes_of_interest)
    print_analysis_report(data)
    
    scaffold_mapping_file = os.path.join(output_dir, "scaffold_mapping.tsv")
    print(f"\nCreating scaffold mapping file: {scaffold_mapping_file}")
    create_scaffold_mapping_file(data, scaffold_mapping_file)
    
    gene_index_file = os.path.join(output_dir, "gene_index.tsv")
    print(f"Creating gene index file: {gene_index_file}")
    create_gene_index_file(data, gene_index_file)
    
    genes_of_interest_file = os.path.join(output_dir, "genes_of_interest.tsv")
    print(f"Creating genes of interest file: {genes_of_interest_file}")
    create_genes_of_interest_file(data, genes_of_interest_file)
    
    print("Done!")

if __name__ == "__main__":
    main()