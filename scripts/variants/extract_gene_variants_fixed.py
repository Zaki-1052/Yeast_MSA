#!/usr/bin/env python3
"""
extract_gene_variants_fixed.py - Extract variants affecting ergosterol pathway genes from annotated VCF files

This script parses annotated VCF files, extracts variants that affect the 11 genes of interest
in the ergosterol biosynthesis pathway, and generates comprehensive tables for analysis.

The script has been fixed to correctly handle variant inclusion, ensuring that control samples
are properly used for comparison and that variants are not erroneously duplicated across treatments.

Usage:
    python extract_gene_variants_fixed.py --vcf_dir <vcf_directory> --gene_mapping <gene_mapping_file> --output_dir <output_directory>
"""

import os
import sys
import argparse
import csv
import gzip
import re
from collections import defaultdict

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Extract variants affecting ergosterol pathway genes from annotated VCF files')
    parser.add_argument('--vcf_dir', required=True, help='Directory containing annotated VCF files')
    parser.add_argument('--gene_mapping', required=True, help='TSV file mapping genes of interest')
    parser.add_argument('--output_dir', required=True, help='Directory for output files')
    parser.add_argument('--upstream_distance', type=int, default=1000, 
                        help='Distance upstream to include for regulatory variants (default: 1000bp)')
    parser.add_argument('--exclude_controls', action='store_true', 
                        help='Exclude variants present in control samples from treatment samples')
    return parser.parse_args()

def load_gene_mapping(mapping_file):
    """
    Load the gene mapping information for the 11 target genes.
    
    Returns:
        dict: Dictionary mapping w303_gene_id to gene information
        dict: Dictionary mapping scaffold to list of genes on that scaffold
        dict: Dictionary mapping sample name to treatment group
    """
    genes = {}
    scaffold_genes = defaultdict(list)
    
    with open(mapping_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            # Store gene information
            gene_id = row['w303_gene_id']
            genes[gene_id] = {
                'sc_gene_id': row['sc_gene_id'],
                'erg_name': row['erg_name'],
                'locus_tag': row['locus_tag'],
                'w303_scaffold': row['w303_scaffold'],
                'start': int(row['start']),
                'end': int(row['end']),
                'strand': row['strand'],
                'product': row['product']
            }
            
            # Store genes by scaffold for faster lookup
            scaffold_genes[row['w303_scaffold']].append(gene_id)
    
    # Define sample to treatment mapping
    sample_treatment = {
        'WT-CTRL': 'WT-CTRL',
        'WT-37-55-1': 'WT-37',
        'WT-37-55-2': 'WT-37',
        'WT-37-55-3': 'WT-37',
        'WTA-55-1': 'WTA',
        'WTA-55-2': 'WTA',
        'WTA-55-3': 'WTA',
        'STC-55-1': 'STC',
        'STC-55-2': 'STC',
        'STC-55-3': 'STC',
        'STC-CTRL': 'STC-CTRL',
        'CAS-55-1': 'CAS',
        'CAS-55-2': 'CAS',
        'CAS-55-3': 'CAS',
        'CAS-CTRL': 'CAS-CTRL'
    }
    
    return genes, scaffold_genes, sample_treatment

def create_output_dirs(output_dir):
    """Create the necessary output directories."""
    os.makedirs(os.path.join(output_dir, 'gene_specific'), exist_ok=True)
    os.makedirs(os.path.join(output_dir, 'treatment_specific'), exist_ok=True)
    
    return {
        'master': os.path.join(output_dir, 'all_gene_variants.tsv'),
        'gene_specific': os.path.join(output_dir, 'gene_specific'),
        'treatment_specific': os.path.join(output_dir, 'treatment_specific')
    }

def initialize_output_files(output_paths, genes, sample_treatment):
    """Initialize all output files with headers."""
    # Define the header for all files
    header = ['Sample', 'Treatment', 'Scaffold', 'Position', 'Ref', 'Alt', 'Quality', 
              'Gene_ID', 'ERG_Name', 'SC_Gene_ID', 'Effect', 'Impact', 'Codon_Change', 
              'Amino_Acid_Change', 'Distance', 'Variant_Type']
    
    # Initialize master file
    with open(output_paths['master'], 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(header)
    
    # Initialize gene-specific files
    for gene_id in genes:
        gene_file = os.path.join(output_paths['gene_specific'], f"{gene_id}.tsv")
        with open(gene_file, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(header)
    
    # Initialize treatment-specific files
    treatments = set(sample_treatment.values())
    for treatment in treatments:
        treatment_file = os.path.join(output_paths['treatment_specific'], f"{treatment}.tsv")
        with open(treatment_file, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(header)

def is_within_gene_region(position, gene_info, upstream_distance):
    """Check if a position is within a gene region (including upstream regulatory region)."""
    # For + strand genes
    if gene_info['strand'] == '+':
        return (gene_info['start'] - upstream_distance) <= position <= gene_info['end']
    # For - strand genes
    else:
        return gene_info['start'] <= position <= (gene_info['end'] + upstream_distance)

def get_variant_type(ref, alt):
    """Determine the variant type based on ref and alt alleles."""
    if len(ref) == 1 and len(alt) == 1:
        return 'SNV'
    elif len(ref) > len(alt):
        return 'DELETION'
    elif len(ref) < len(alt):
        return 'INSERTION'
    else:
        return 'COMPLEX'

def parse_ann_field(ann_str):
    """Parse the ANN field from VCF into a usable format."""
    # ANN format: A|missense_variant|MODERATE|W3030G02340|Gene_469848_471270|transcript|W303_0G02340|protein_coding||c.123A>T|p.Gln41His|123|123|41|Q/H|cAa/cTa||
    annotations = []
    
    # Split multiple annotations
    for ann in ann_str.split(','):
        fields = ann.split('|')
        if len(fields) < 10:  # Basic sanity check
            continue
        
        annotation = {
            'alt': fields[0],
            'effect': fields[1],
            'impact': fields[2],
            'gene_id': fields[3],
            'codon_change': fields[9],
            'aa_change': fields[10],
            'distance': fields[14] if len(fields) > 14 and fields[14].isdigit() else ''
        }
        annotations.append(annotation)
    
    return annotations

def extract_all_variants(vcf_dir, genes, scaffold_genes, sample_treatment, upstream_distance):
    """
    Extract all variants from VCF files first, then filter and organize them.
    This allows for proper treatment-control comparisons.
    """
    all_variants = []
    control_variants = set()  # Set of (scaffold, pos, ref, alt, gene_id) tuples for control samples
    
    # Process each VCF file
    print(f"Processing VCF files from {vcf_dir}")
    vcf_files = [os.path.join(vcf_dir, f) for f in os.listdir(vcf_dir) 
                 if f.endswith('.vcf') or f.endswith('.vcf.gz')]
    
    # First, collect all variants
    for vcf_file in sorted(vcf_files):
        # Extract sample name from filename
        sample_name = os.path.basename(vcf_file).split('.')[0]
        if sample_name not in sample_treatment:
            print(f"Warning: Unknown sample {sample_name}, skipping")
            continue
        
        treatment = sample_treatment[sample_name]
        print(f"Processing {sample_name} (Treatment: {treatment})")
        
        # Flag if this is a control sample
        is_control = 'CTRL' in sample_name
        
        # Open VCF file (handle both gzipped and regular)
        opener = gzip.open if vcf_file.endswith('.gz') else open
        with opener(vcf_file, 'rt') as f:
            for line in f:
                # Skip header lines
                if line.startswith('#'):
                    continue
                
                # Parse VCF line
                fields = line.strip().split('\t')
                if len(fields) < 8:
                    continue
                    
                chrom, pos, _, ref, alt, qual, filter_val, info = fields[:8]
                pos = int(pos)
                
                # Skip if not on a scaffold with our genes of interest
                if chrom not in scaffold_genes:
                    continue
                
                # Find ANN field
                ann_match = re.search(r'ANN=([^;]+)', info)
                if not ann_match:
                    continue
                    
                ann_str = ann_match.group(1)
                annotations = parse_ann_field(ann_str)
                
                # Check each annotation for our genes of interest
                for annotation in annotations:
                    gene_id = annotation['gene_id']
                    
                    # Skip if not one of our target genes
                    if gene_id not in genes:
                        continue
                    
                    # Check if position is within gene region (including upstream)
                    gene_info = genes[gene_id]
                    if not is_within_gene_region(pos, gene_info, upstream_distance):
                        continue
                    
                    # Create variant record
                    variant_record = {
                        'sample': sample_name,
                        'treatment': treatment,
                        'scaffold': chrom,
                        'position': pos,
                        'ref': ref,
                        'alt': alt,
                        'quality': qual,
                        'gene_id': gene_id,
                        'erg_name': gene_info['erg_name'],
                        'sc_gene_id': gene_info['sc_gene_id'],
                        'effect': annotation['effect'],
                        'impact': annotation['impact'],
                        'codon_change': annotation['codon_change'],
                        'aa_change': annotation['aa_change'],
                        'distance': annotation['distance'],
                        'variant_type': get_variant_type(ref, alt)
                    }
                    
                    # Add to all variants
                    all_variants.append(variant_record)
                    
                    # Add to control variants if this is a control sample
                    if is_control:
                        control_variants.add((chrom, pos, ref, alt, gene_id))
    
    return all_variants, control_variants

def write_variants_to_files(all_variants, control_variants, output_paths, exclude_controls):
    """Write processed variants to output files, applying any filtering."""
    # Initialize counters
    gene_counts = defaultdict(int)
    treatment_counts = defaultdict(int)
    sample_counts = defaultdict(int)
    effect_counts = defaultdict(int)
    impact_counts = defaultdict(int)
    
    # Write to files
    for variant in all_variants:
        # Skip if this variant is present in control samples and we're excluding control variants
        variant_key = (
            variant['scaffold'], 
            variant['position'], 
            variant['ref'], 
            variant['alt'], 
            variant['gene_id']
        )
        if exclude_controls and variant_key in control_variants and 'CTRL' not in variant['treatment']:
            continue
        
        # Prepare the variant record for writing
        variant_record = [
            variant['sample'],
            variant['treatment'],
            variant['scaffold'],
            variant['position'],
            variant['ref'],
            variant['alt'],
            variant['quality'],
            variant['gene_id'],
            variant['erg_name'],
            variant['sc_gene_id'],
            variant['effect'],
            variant['impact'],
            variant['codon_change'],
            variant['aa_change'],
            variant['distance'],
            variant['variant_type']
        ]
        
        # Write to master file
        with open(output_paths['master'], 'a') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(variant_record)
        
        # Update counters
        gene_counts[variant['gene_id']] += 1
        treatment_counts[variant['treatment']] += 1
        sample_counts[variant['sample']] += 1
        effect_counts[variant['effect']] += 1
        impact_counts[variant['impact']] += 1
        
        # Write to gene-specific file
        gene_file = os.path.join(output_paths['gene_specific'], f"{variant['gene_id']}.tsv")
        with open(gene_file, 'a') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(variant_record)
        
        # Write to treatment-specific file
        treatment_file = os.path.join(output_paths['treatment_specific'], f"{variant['treatment']}.tsv")
        with open(treatment_file, 'a') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(variant_record)
    
    return gene_counts, treatment_counts, sample_counts, effect_counts, impact_counts

def generate_summary_statistics(output_dir, genes, gene_counts, treatment_counts, sample_counts, effect_counts, impact_counts):
    """Generate summary statistics for the extracted variants."""
    print("\nGenerating summary statistics...")
    
    # Calculate total variants
    total_variants = sum(gene_counts.values())
    
    # Write gene summary
    with open(os.path.join(output_dir, 'gene_summary.tsv'), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Gene_ID', 'ERG_Name', 'SC_Gene_ID', 'Variant_Count', 'Percentage'])
        for gene_id in genes:
            count = gene_counts.get(gene_id, 0)
            percentage = (count / total_variants * 100) if total_variants > 0 else 0
            writer.writerow([
                gene_id, 
                genes[gene_id]['erg_name'], 
                genes[gene_id]['sc_gene_id'], 
                count, 
                f"{percentage:.2f}%"
            ])
    
    # Write treatment summary
    with open(os.path.join(output_dir, 'treatment_summary.tsv'), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Treatment', 'Variant_Count', 'Percentage'])
        for treatment, count in sorted(treatment_counts.items()):
            percentage = (count / total_variants * 100) if total_variants > 0 else 0
            writer.writerow([treatment, count, f"{percentage:.2f}%"])
    
    # Write sample summary
    with open(os.path.join(output_dir, 'sample_summary.tsv'), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Sample', 'Variant_Count'])
        for sample, count in sorted(sample_counts.items()):
            writer.writerow([sample, count])
    
    # Write effect summary
    with open(os.path.join(output_dir, 'effect_summary.tsv'), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Effect', 'Variant_Count', 'Percentage'])
        for effect, count in sorted(effect_counts.items(), key=lambda x: x[1], reverse=True):
            percentage = (count / total_variants * 100) if total_variants > 0 else 0
            writer.writerow([effect, count, f"{percentage:.2f}%"])
    
    # Write impact summary
    with open(os.path.join(output_dir, 'impact_summary.tsv'), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Impact', 'Variant_Count', 'Percentage'])
        for impact, count in sorted(impact_counts.items(), key=lambda x: x[1], reverse=True):
            percentage = (count / total_variants * 100) if total_variants > 0 else 0
            writer.writerow([impact, count, f"{percentage:.2f}%"])
    
    print(f"Total variants in target genes: {total_variants}")
    print(f"Summary files written to {output_dir}")

def main():
    """Main function to extract variants affecting ergosterol pathway genes."""
    args = parse_arguments()
    
    # Load gene mapping
    print(f"Loading gene mapping from {args.gene_mapping}")
    genes, scaffold_genes, sample_treatment = load_gene_mapping(args.gene_mapping)
    print(f"Loaded information for {len(genes)} target genes")
    
    # Create output directories
    print(f"Creating output directories in {args.output_dir}")
    output_paths = create_output_dirs(args.output_dir)
    
    # Initialize output files
    print("Initializing output files")
    initialize_output_files(output_paths, genes, sample_treatment)
    
    # Extract all variants first
    all_variants, control_variants = extract_all_variants(
        args.vcf_dir, genes, scaffold_genes, sample_treatment, args.upstream_distance
    )
    
    # Print some statistics about control variants
    print(f"Found {len(control_variants)} unique variants in control samples")
    print(f"Found {len(all_variants)} total variant records across all samples")
    
    # Write variants to files, filtering as needed
    gene_counts, treatment_counts, sample_counts, effect_counts, impact_counts = write_variants_to_files(
        all_variants, control_variants, output_paths, args.exclude_controls
    )
    
    # Generate summary statistics
    generate_summary_statistics(
        args.output_dir, genes, gene_counts, treatment_counts, sample_counts, effect_counts, impact_counts
    )
    
    print("Done!")

if __name__ == "__main__":
    main()