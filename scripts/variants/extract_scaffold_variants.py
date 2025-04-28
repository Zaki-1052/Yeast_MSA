#!/usr/bin/env python3
"""
extract_scaffold_variants.py - Extract variants from scaffolds containing ergosterol pathway genes

This script extracts all variants from the scaffolds containing our genes of interest,
calculates their distance to the nearest target gene, and generates statistics and 
visualizations on variant distribution patterns.

Usage:
    python extract_scaffold_variants.py --vcf_dir <vcf_directory> --gene_mapping <gene_mapping_file> --output_dir <output_directory>
"""

import os
import sys
import argparse
import csv
import gzip
import re
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Extract variants from scaffolds containing ergosterol pathway genes')
    parser.add_argument('--vcf_dir', required=True, help='Directory containing annotated VCF files')
    parser.add_argument('--gene_mapping', required=True, help='TSV file mapping genes of interest')
    parser.add_argument('--output_dir', required=True, help='Directory for output files')
    parser.add_argument('--distance_bins', type=str, default='0,500,1000,5000,10000,50000', 
                        help='Comma-separated list of distance bins for categorization (default: 0,500,1000,5000,10000,50000)')
    parser.add_argument('--max_distance', type=int, default=50000,
                        help='Maximum distance to consider for gene proximity (default: 50000bp)')
    return parser.parse_args()

def load_gene_mapping(mapping_file):
    """
    Load the gene mapping information for the target genes.
    
    Returns:
        dict: Dictionary mapping gene ID to gene information
        set: Set of scaffolds containing target genes
        dict: Dictionary mapping scaffold to list of genes on that scaffold
        dict: Dictionary mapping sample name to treatment group
    """
    genes = {}
    scaffolds_of_interest = set()
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
                'scaffold': row['w303_scaffold'],
                'start': int(row['start']),
                'end': int(row['end']),
                'strand': row['strand'],
                'product': row['product']
            }
            
            # Store scaffold information
            scaffolds_of_interest.add(row['w303_scaffold'])
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
    
    return genes, scaffolds_of_interest, scaffold_genes, sample_treatment

def create_output_dirs(output_dir):
    """Create the necessary output directories."""
    os.makedirs(os.path.join(output_dir, 'scaffold_variants'), exist_ok=True)
    os.makedirs(os.path.join(output_dir, 'gene_proximity'), exist_ok=True)
    os.makedirs(os.path.join(output_dir, 'visualizations'), exist_ok=True)
    
    return {
        'base': output_dir,
        'scaffold_variants': os.path.join(output_dir, 'scaffold_variants'),
        'gene_proximity': os.path.join(output_dir, 'gene_proximity'),
        'visualizations': os.path.join(output_dir, 'visualizations')
    }

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

def parse_vcf_file(vcf_file, scaffolds_of_interest, sample_treatment):
    """
    Parse a VCF file and extract variants on scaffolds of interest.
    
    Returns:
        list: List of variant dictionaries for the sample
        str: Sample name
        str: Treatment group
    """
    # Extract sample name from filename
    sample_name = os.path.basename(vcf_file).split('.')[0]
    if sample_name not in sample_treatment:
        print(f"Warning: Unknown sample {sample_name}, skipping")
        return [], None, None
    
    treatment = sample_treatment[sample_name]
    print(f"Processing {sample_name} (Treatment: {treatment})")
    
    variants = []
    
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
                
            chrom, pos, id_, ref, alt, qual, filter_val, info = fields[:8]
            pos = int(pos)
            
            # Skip if not on a scaffold of interest
            if chrom not in scaffolds_of_interest:
                continue
            
            # Extract variant information
            variant = {
                'scaffold': chrom,
                'position': pos,
                'ref': ref,
                'alt': alt,
                'quality': qual,
                'filter': filter_val,
                'variant_type': get_variant_type(ref, alt),
                'sample': sample_name,
                'treatment': treatment
            }
            
            # Parse ANN field if present
            ann_match = re.search(r'ANN=([^;]+)', info)
            if ann_match:
                ann_str = ann_match.group(1)
                annotations = []
                
                # Split multiple annotations
                for ann in ann_str.split(','):
                    fields = ann.split('|')
                    if len(fields) < 3:  # Basic sanity check
                        continue
                    
                    annotation = {
                        'alt': fields[0],
                        'effect': fields[1],
                        'impact': fields[2],
                        'gene_id': fields[3] if len(fields) > 3 else ''
                    }
                    annotations.append(annotation)
                
                variant['annotations'] = annotations
                
                # Set primary effect and impact (first annotation)
                if annotations:
                    variant['effect'] = annotations[0]['effect']
                    variant['impact'] = annotations[0]['impact']
                    variant['gene_id'] = annotations[0]['gene_id']
                else:
                    variant['effect'] = ''
                    variant['impact'] = ''
                    variant['gene_id'] = ''
            else:
                variant['annotations'] = []
                variant['effect'] = ''
                variant['impact'] = ''
                variant['gene_id'] = ''
            
            variants.append(variant)
    
    return variants, sample_name, treatment

def calculate_distance_to_genes(variants, genes, scaffold_genes):
    """
    Calculate the distance of each variant to the nearest target gene.
    
    Args:
        variants: List of variant dictionaries
        genes: Dictionary mapping gene ID to gene information
        scaffold_genes: Dictionary mapping scaffold to list of genes on that scaffold
        
    Returns:
        list: Updated variant list with distance information
    """
    for variant in variants:
        scaffold = variant['scaffold']
        position = variant['position']
        
        min_distance = float('inf')
        nearest_gene = None
        position_relative = None
        
        # Check each gene on this scaffold
        for gene_id in scaffold_genes.get(scaffold, []):
            gene_info = genes[gene_id]
            gene_start = gene_info['start']
            gene_end = gene_info['end']
            
            # Check if variant is within gene boundaries
            if gene_start <= position <= gene_end:
                distance = 0
                position_relative = 'within'
            # Check if variant is upstream of gene
            elif position < gene_start:
                distance = gene_start - position
                position_relative = 'upstream'
            # Check if variant is downstream of gene
            else:  # position > gene_end
                distance = position - gene_end
                position_relative = 'downstream'
            
            # Update nearest gene if this distance is smaller
            if distance < min_distance:
                min_distance = distance
                nearest_gene = gene_id
                variant['position_relative'] = position_relative
        
        variant['distance_to_nearest_gene'] = min_distance
        variant['nearest_gene'] = nearest_gene
        
        # Add nearest gene info
        if nearest_gene and nearest_gene in genes:
            variant['nearest_gene_sc_id'] = genes[nearest_gene]['sc_gene_id']
            variant['nearest_gene_erg_name'] = genes[nearest_gene]['erg_name']
        else:
            variant['nearest_gene_sc_id'] = ''
            variant['nearest_gene_erg_name'] = ''
    
    return variants

def categorize_variants_by_distance(variants, distance_bins):
    """
    Categorize variants by distance to nearest gene.
    
    Args:
        variants: List of variant dictionaries with distance information
        distance_bins: List of distance thresholds for binning
        
    Returns:
        dict: Dictionary mapping distance category to list of variants
    """
    categories = defaultdict(list)
    
    for variant in variants:
        distance = variant['distance_to_nearest_gene']
        
        # Find the appropriate bin
        for i in range(len(distance_bins) - 1):
            if distance_bins[i] <= distance < distance_bins[i+1]:
                category = f"{distance_bins[i]}-{distance_bins[i+1]}"
                categories[category].append(variant)
                variant['distance_category'] = category
                break
        else:
            # If distance is greater than the last bin
            if distance >= distance_bins[-1]:
                category = f">{distance_bins[-1]}"
                categories[category].append(variant)
                variant['distance_category'] = category
    
    return categories

def write_variant_files(variants, categories, scaffold_genes, output_paths):
    """Write variant data to output files."""
    # Write all scaffold variants to a single file
    with open(os.path.join(output_paths['base'], 'all_scaffold_variants.tsv'), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow([
            'Sample', 'Treatment', 'Scaffold', 'Position', 'Ref', 'Alt', 'Quality', 'Variant_Type',
            'Effect', 'Impact', 'Gene_ID', 'Nearest_Gene', 'Nearest_Gene_SC_ID', 'Nearest_Gene_ERG',
            'Distance', 'Position_Relative', 'Distance_Category'
        ])
        
        for variant in variants:
            writer.writerow([
                variant['sample'],
                variant['treatment'],
                variant['scaffold'],
                variant['position'],
                variant['ref'],
                variant['alt'],
                variant['quality'],
                variant['variant_type'],
                variant['effect'],
                variant['impact'],
                variant['gene_id'],
                variant.get('nearest_gene', ''),
                variant.get('nearest_gene_sc_id', ''),
                variant.get('nearest_gene_erg_name', ''),
                variant.get('distance_to_nearest_gene', ''),
                variant.get('position_relative', ''),
                variant.get('distance_category', '')
            ])
    
    # Write variants for each scaffold
    for scaffold, genes_list in scaffold_genes.items():
        scaffold_variants = [v for v in variants if v['scaffold'] == scaffold]
        
        if not scaffold_variants:
            continue
        
        with open(os.path.join(output_paths['scaffold_variants'], f"{scaffold}_variants.tsv"), 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow([
                'Sample', 'Treatment', 'Position', 'Ref', 'Alt', 'Quality', 'Variant_Type',
                'Effect', 'Impact', 'Gene_ID', 'Nearest_Gene', 'Nearest_Gene_SC_ID', 'Nearest_Gene_ERG',
                'Distance', 'Position_Relative'
            ])
            
            for variant in scaffold_variants:
                writer.writerow([
                    variant['sample'],
                    variant['treatment'],
                    variant['position'],
                    variant['ref'],
                    variant['alt'],
                    variant['quality'],
                    variant['variant_type'],
                    variant['effect'],
                    variant['impact'],
                    variant['gene_id'],
                    variant.get('nearest_gene', ''),
                    variant.get('nearest_gene_sc_id', ''),
                    variant.get('nearest_gene_erg_name', ''),
                    variant.get('distance_to_nearest_gene', ''),
                    variant.get('position_relative', '')
                ])
    
    # Write variants for each distance category
    for category, category_variants in categories.items():
        with open(os.path.join(output_paths['gene_proximity'], f"distance_{category.replace('>', 'over_')}.tsv"), 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow([
                'Sample', 'Treatment', 'Scaffold', 'Position', 'Ref', 'Alt', 'Quality', 'Variant_Type',
                'Effect', 'Impact', 'Gene_ID', 'Nearest_Gene', 'Nearest_Gene_SC_ID', 'Nearest_Gene_ERG',
                'Distance', 'Position_Relative'
            ])
            
            for variant in category_variants:
                writer.writerow([
                    variant['sample'],
                    variant['treatment'],
                    variant['scaffold'],
                    variant['position'],
                    variant['ref'],
                    variant['alt'],
                    variant['quality'],
                    variant['variant_type'],
                    variant['effect'],
                    variant['impact'],
                    variant['gene_id'],
                    variant.get('nearest_gene', ''),
                    variant.get('nearest_gene_sc_id', ''),
                    variant.get('nearest_gene_erg_name', ''),
                    variant.get('distance_to_nearest_gene', ''),
                    variant.get('position_relative', '')
                ])

def generate_summary_statistics(variants, categories, genes, scaffold_genes, output_paths):
    """Generate summary statistics on variant distribution."""
    # Count variants by scaffold
    scaffold_counts = defaultdict(int)
    for variant in variants:
        scaffold_counts[variant['scaffold']] += 1
    
    # Count variants by distance category
    category_counts = {category: len(variants_list) for category, variants_list in categories.items()}
    
    # Count variants by gene proximity
    gene_proximity_counts = defaultdict(int)
    for variant in variants:
        if 'nearest_gene' in variant and variant['nearest_gene']:
            gene_proximity_counts[variant['nearest_gene']] += 1
    
    # Count variants by treatment
    treatment_counts = defaultdict(int)
    for variant in variants:
        treatment_counts[variant['treatment']] += 1
    
    # Count variants by effect
    effect_counts = defaultdict(int)
    for variant in variants:
        if variant['effect']:
            effect_counts[variant['effect']] += 1
    
    # Count variants by impact
    impact_counts = defaultdict(int)
    for variant in variants:
        if variant['impact']:
            impact_counts[variant['impact']] += 1
    
    # Count variants by type
    type_counts = defaultdict(int)
    for variant in variants:
        if 'variant_type' in variant:
            type_counts[variant['variant_type']] += 1
    
    # Write summary statistics
    with open(os.path.join(output_paths['base'], 'scaffold_variant_summary.txt'), 'w') as f:
        f.write("Scaffold Variant Analysis Summary\n")
        f.write("===============================\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write(f"Total variants on target scaffolds: {len(variants)}\n")
        f.write(f"Target scaffolds: {len(scaffold_counts)}\n")
        f.write(f"Target genes: {len(genes)}\n\n")
        
        f.write("Variants by scaffold:\n")
        for scaffold, count in sorted(scaffold_counts.items(), key=lambda x: x[1], reverse=True):
            genes_on_scaffold = scaffold_genes.get(scaffold, [])
            gene_names = [f"{genes[gene_id]['sc_gene_id']} ({genes[gene_id]['erg_name']})" for gene_id in genes_on_scaffold]
            f.write(f"  {scaffold}: {count} variants, Genes: {', '.join(gene_names)}\n")
        f.write("\n")
        
        f.write("Variants by distance category:\n")
        for category, count in sorted(category_counts.items(), key=lambda x: x[0]):
            percentage = (count / len(variants) * 100) if variants else 0
            f.write(f"  {category}: {count} ({percentage:.2f}%)\n")
        f.write("\n")
        
        f.write("Variants by nearest gene:\n")
        for gene_id, count in sorted(gene_proximity_counts.items(), key=lambda x: x[1], reverse=True):
            if gene_id in genes:
                gene_info = genes[gene_id]
                percentage = (count / len(variants) * 100) if variants else 0
                f.write(f"  {gene_info['sc_gene_id']} ({gene_info['erg_name']}): {count} ({percentage:.2f}%)\n")
        f.write("\n")
        
        f.write("Variants by treatment:\n")
        for treatment, count in sorted(treatment_counts.items(), key=lambda x: x[1], reverse=True):
            percentage = (count / len(variants) * 100) if variants else 0
            f.write(f"  {treatment}: {count} ({percentage:.2f}%)\n")
        f.write("\n")
        
        f.write("Variants by type:\n")
        for variant_type, count in sorted(type_counts.items(), key=lambda x: x[1], reverse=True):
            percentage = (count / len(variants) * 100) if variants else 0
            f.write(f"  {variant_type}: {count} ({percentage:.2f}%)\n")
        f.write("\n")
        
        f.write("Top 10 variant effects:\n")
        for effect, count in sorted(effect_counts.items(), key=lambda x: x[1], reverse=True)[:10]:
            percentage = (count / len(variants) * 100) if variants else 0
            f.write(f"  {effect}: {count} ({percentage:.2f}%)\n")
        f.write("\n")
        
        f.write("Variants by impact:\n")
        for impact, count in sorted(impact_counts.items(), key=lambda x: x[1], reverse=True):
            percentage = (count / len(variants) * 100) if variants else 0
            f.write(f"  {impact}: {count} ({percentage:.2f}%)\n")
    
    # Write detailed gene proximity statistics
    with open(os.path.join(output_paths['base'], 'gene_proximity_summary.tsv'), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Gene_ID', 'ERG_Name', 'SC_Gene_ID', 'Variants_Within', 'Variants_Upstream_1kb', 
                         'Variants_Downstream_1kb', 'Variants_Upstream_5kb', 'Variants_Downstream_5kb', 
                         'Total_Within_5kb'])
        
        gene_position_counts = defaultdict(lambda: defaultdict(int))
        for variant in variants:
            if 'nearest_gene' in variant and variant['nearest_gene']:
                gene_id = variant['nearest_gene']
                position_relative = variant.get('position_relative', 'unknown')
                distance = variant.get('distance_to_nearest_gene', float('inf'))
                
                if position_relative == 'within':
                    gene_position_counts[gene_id]['within'] += 1
                elif position_relative == 'upstream' and distance <= 1000:
                    gene_position_counts[gene_id]['upstream_1kb'] += 1
                elif position_relative == 'downstream' and distance <= 1000:
                    gene_position_counts[gene_id]['downstream_1kb'] += 1
                elif position_relative == 'upstream' and distance <= 5000:
                    gene_position_counts[gene_id]['upstream_5kb'] += 1
                elif position_relative == 'downstream' and distance <= 5000:
                    gene_position_counts[gene_id]['downstream_5kb'] += 1
        
        for gene_id in sorted(genes.keys()):
            if gene_id in gene_position_counts:
                counts = gene_position_counts[gene_id]
                within = counts.get('within', 0)
                upstream_1kb = counts.get('upstream_1kb', 0)
                downstream_1kb = counts.get('downstream_1kb', 0)
                upstream_5kb = counts.get('upstream_5kb', 0)
                downstream_5kb = counts.get('downstream_5kb', 0)
                total_5kb = within + upstream_1kb + downstream_1kb + upstream_5kb + downstream_5kb
                
                writer.writerow([
                    gene_id,
                    genes[gene_id]['erg_name'],
                    genes[gene_id]['sc_gene_id'],
                    within,
                    upstream_1kb,
                    downstream_1kb,
                    upstream_5kb,
                    downstream_5kb,
                    total_5kb
                ])
            else:
                writer.writerow([
                    gene_id,
                    genes[gene_id]['erg_name'],
                    genes[gene_id]['sc_gene_id'],
                    0, 0, 0, 0, 0, 0
                ])
    
    # Write treatment comparison statistics
    with open(os.path.join(output_paths['base'], 'treatment_comparison.tsv'), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Treatment', 'Total_Variants', 'Variants_Within_Genes', 'Variants_Within_1kb', 
                         'Variants_Within_5kb', 'High_Impact', 'Moderate_Impact', 'Low_Impact', 'Modifier_Impact'])
        
        treatment_detail_counts = defaultdict(lambda: defaultdict(int))
        for variant in variants:
            treatment = variant['treatment']
            position_relative = variant.get('position_relative', 'unknown')
            distance = variant.get('distance_to_nearest_gene', float('inf'))
            impact = variant.get('impact', '')
            
            treatment_detail_counts[treatment]['total'] += 1
            
            if position_relative == 'within':
                treatment_detail_counts[treatment]['within'] += 1
            
            if position_relative in ['within', 'upstream', 'downstream'] and distance <= 1000:
                treatment_detail_counts[treatment]['within_1kb'] += 1
            
            if position_relative in ['within', 'upstream', 'downstream'] and distance <= 5000:
                treatment_detail_counts[treatment]['within_5kb'] += 1
            
            if impact:
                treatment_detail_counts[treatment][impact.lower()] += 1
        
        for treatment in sorted(treatment_detail_counts.keys()):
            counts = treatment_detail_counts[treatment]
            writer.writerow([
                treatment,
                counts['total'],
                counts['within'],
                counts['within_1kb'],
                counts['within_5kb'],
                counts['high'],
                counts['moderate'],
                counts['low'],
                counts['modifier']
            ])

def generate_visualizations(variants, genes, scaffold_genes, categories, output_paths, max_distance):
    """Generate visualizations of variant distribution."""
    # 1. Distance distribution histogram
    plt.figure(figsize=(12, 6))
    distances = [min(variant.get('distance_to_nearest_gene', 0), max_distance) 
                 for variant in variants if 'distance_to_nearest_gene' in variant]
    plt.hist(distances, bins=50, alpha=0.75)
    plt.xlabel('Distance to Nearest Gene (bp)')
    plt.ylabel('Number of Variants')
    plt.title('Distribution of Variants by Distance to Nearest Gene')
    plt.grid(True, alpha=0.3)
    plt.savefig(os.path.join(output_paths['visualizations'], 'distance_distribution.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Variant count by distance category
    plt.figure(figsize=(14, 6))
    categories_sorted = sorted(categories.keys(), key=lambda x: (
        int(x.split('-')[0].replace('>', '0')) if x.split('-')[0].replace('>', '0').isdigit() else 0, 
        int(x.split('-')[1]) if '-' in x and x.split('-')[1].isdigit() else float('inf')
    ))
    counts = [len(categories[cat]) for cat in categories_sorted]
    
    # Create better labels for categories
    def prettify_category(cat):
        if cat.startswith('>'):
            return f">{cat[1:]} bp"
        else:
            parts = cat.split('-')
            return f"{parts[0]}-{parts[1]} bp"
    
    labels = [prettify_category(cat) for cat in categories_sorted]
    
    plt.bar(labels, counts, alpha=0.7)
    plt.xlabel('Distance Category')
    plt.ylabel('Number of Variants')
    plt.title('Variant Count by Distance to Nearest Gene')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.grid(True, alpha=0.3, axis='y')
    plt.savefig(os.path.join(output_paths['visualizations'], 'distance_category_counts.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Variant count by scaffold and gene
    scaffolds = sorted(scaffold_genes.keys())
    gene_counts = defaultdict(lambda: defaultdict(int))
    
    for variant in variants:
        if 'nearest_gene' in variant and variant['nearest_gene'] and variant.get('distance_to_nearest_gene', float('inf')) <= 5000:
            scaffold = variant['scaffold']
            gene_id = variant['nearest_gene']
            gene_counts[scaffold][gene_id] += 1
    
    # Create stacked bar chart for each scaffold
    plt.figure(figsize=(14, 8))
    
    x = np.arange(len(scaffolds))
    bottom = np.zeros(len(scaffolds))
    
    # Create a color map for genes
    gene_ids = sorted(genes.keys())
    colors = plt.cm.tab20(np.linspace(0, 1, len(gene_ids)))
    gene_colors = {gene_id: colors[i] for i, gene_id in enumerate(gene_ids)}
    
    # Plot each gene's variants as a stack
    for gene_id in gene_ids:
        gene_data = [gene_counts[scaffold][gene_id] for scaffold in scaffolds]
        if sum(gene_data) > 0:  # Only plot genes with variants
            plt.bar(x, gene_data, bottom=bottom, width=0.8, 
                    label=f"{genes[gene_id]['sc_gene_id']} ({genes[gene_id]['erg_name']})",
                    color=gene_colors[gene_id], alpha=0.7)
            bottom += np.array(gene_data)
    
    plt.xlabel('Scaffold')
    plt.ylabel('Number of Variants (within 5kb of gene)')
    plt.title('Variants by Scaffold and Nearest Gene')
    plt.xticks(x, scaffolds, rotation=45, ha='right')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(os.path.join(output_paths['visualizations'], 'scaffold_gene_counts.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 4. Gene proximity heatmap
    gene_treatment_counts = defaultdict(lambda: defaultdict(int))
    for variant in variants:
        if 'nearest_gene' in variant and variant['nearest_gene'] and variant.get('distance_to_nearest_gene', float('inf')) < 5000:
            treatment = variant['treatment']
            gene_id = variant['nearest_gene']
            gene_treatment_counts[gene_id][treatment] += 1
    
    # Create data for heatmap
    gene_ids = sorted(gene_treatment_counts.keys())
    treatments = sorted(set(treatment for counts in gene_treatment_counts.values() for treatment in counts.keys()))
    
    data = np.zeros((len(gene_ids), len(treatments)))
    for i, gene_id in enumerate(gene_ids):
        for j, treatment in enumerate(treatments):
            data[i, j] = gene_treatment_counts[gene_id][treatment]
    
    plt.figure(figsize=(12, 10))
    plt.imshow(data, aspect='auto', cmap='YlOrRd')
    
    # Add gene labels with ERG names
    gene_labels = [f"{genes[gene_id]['sc_gene_id']} ({genes[gene_id]['erg_name']})" for gene_id in gene_ids]
    
    plt.yticks(range(len(gene_ids)), gene_labels)
    plt.xticks(range(len(treatments)), treatments, rotation=45, ha='right')
    plt.colorbar(label='Number of Variants')
    plt.title('Variants Within 5kb of Each Gene by Treatment')
    plt.tight_layout()
    plt.savefig(os.path.join(output_paths['visualizations'], 'gene_treatment_heatmap.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 5. Variant effect distribution
    effect_counts = defaultdict(int)
    for variant in variants:
        if variant['effect']:
            effect_counts[variant['effect']] += 1
    
    # Only plot top 10 effects
    top_effects = sorted(effect_counts.items(), key=lambda x: x[1], reverse=True)[:10]
    effect_labels = [effect for effect, count in top_effects]
    effect_values = [count for effect, count in top_effects]
    
    plt.figure(figsize=(12, 6))
    plt.bar(effect_labels, effect_values, alpha=0.7)
    plt.xlabel('Variant Effect')
    plt.ylabel('Count')
    plt.title('Top 10 Variant Effects')
    plt.xticks(rotation=45, ha='right')
    plt.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(os.path.join(output_paths['visualizations'], 'variant_effect_distribution.png'), dpi=300, bbox_inches='tight')
    plt.close()

def main():
    """Main function to extract and analyze scaffold variants."""
    args = parse_arguments()
    
    # Parse distance bins
    distance_bins = [int(x) for x in args.distance_bins.split(',')]
    
    # Load gene mapping
    print(f"Loading gene mapping from {args.gene_mapping}")
    genes, scaffolds_of_interest, scaffold_genes, sample_treatment = load_gene_mapping(args.gene_mapping)
    print(f"Loaded information for {len(genes)} target genes on {len(scaffolds_of_interest)} scaffolds")
    
    # Create output directories
    print(f"Creating output directories in {args.output_dir}")
    output_paths = create_output_dirs(args.output_dir)
    
    # Process each VCF file
    print(f"Processing VCF files from {args.vcf_dir}")
    vcf_files = [os.path.join(args.vcf_dir, f) for f in os.listdir(args.vcf_dir) 
                 if f.endswith('.vcf') or f.endswith('.vcf.gz')]
    
    all_variants = []
    for vcf_file in sorted(vcf_files):
        sample_variants, sample_name, treatment = parse_vcf_file(
            vcf_file, scaffolds_of_interest, sample_treatment)
        
        if sample_variants:
            # Calculate distance to nearest gene for each variant
            sample_variants = calculate_distance_to_genes(
                sample_variants, genes, scaffold_genes)
            
            all_variants.extend(sample_variants)
    
    print(f"Extracted {len(all_variants)} variants from target scaffolds")
    
    # Categorize variants by distance
    print("Categorizing variants by distance to genes...")
    categories = categorize_variants_by_distance(all_variants, distance_bins)
    
    # Write variant files
    print("Writing variant files...")
    write_variant_files(all_variants, categories, scaffold_genes, output_paths)
    
    # Generate summary statistics
    print("Generating summary statistics...")
    generate_summary_statistics(all_variants, categories, genes, scaffold_genes, output_paths)
    
    # Generate visualizations
    print("Generating visualizations...")
    generate_visualizations(all_variants, genes, scaffold_genes, categories, output_paths, args.max_distance)
    
    print(f"Analysis complete. Results are in {args.output_dir}")

if __name__ == "__main__":
    main()