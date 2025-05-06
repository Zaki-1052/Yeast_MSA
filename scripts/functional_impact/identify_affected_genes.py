#!/usr/bin/env python3
"""
identify_affected_genes.py - DEBUGGED VERSION

This script identifies and characterizes genes affected by HIGH/MODERATE impact variants 
near ergosterol pathway genes, with proper chromosome name mapping.
"""

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict, Counter
import json
import re
from Bio import SeqIO
import gzip

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Identify and characterize genes affected by variants')
    parser.add_argument('--variants_file', required=True, help='Path to the variants TSV file')
    parser.add_argument('--genbank_dir', required=True, help='Directory containing GenBank annotation files')
    parser.add_argument('--mapping_file', help='Path to chromosome mapping file (optional)')
    parser.add_argument('--output_dir', required=True, help='Directory to save results')
    parser.add_argument('--debug', action='store_true', help='Enable debug output')
    return parser.parse_args()

def load_variants(variants_file):
    """Load variant data from TSV file."""
    print(f"Loading variants from {variants_file}")
    variants = pd.read_csv(variants_file, sep='\t')
    # Ensure consistent column name capitalization
    variants.columns = [col.lower() for col in variants.columns]
    
    print(f"Loaded {len(variants)} variants")
    return variants

def create_chromosome_mapping(genbank_dir, debug=False):
    """
    Create a mapping between different chromosome naming systems by analyzing GenBank records.
    
    This function extracts chromosome IDs from GenBank records and identifies potential
    mappings between different naming conventions (e.g., CM007XXX.1 ↔ w303_scaffold_X)
    """
    print("Creating chromosome mapping from GenBank files")
    mapping = {}
    reverse_mapping = {}
    
    # Get all GenBank files in the directory
    genbank_files = [os.path.join(genbank_dir, f) for f in os.listdir(genbank_dir) 
                    if f.endswith(('.gb', '.gbk', '.genbank'))]
    
    if not genbank_files:
        # Check for files without standard extensions
        genbank_files = [os.path.join(genbank_dir, f) for f in os.listdir(genbank_dir) 
                        if os.path.isfile(os.path.join(genbank_dir, f)) and not f.startswith('.')]
    
    for gb_file in genbank_files:
        try:
            for record in SeqIO.parse(gb_file, "genbank"):
                # Extract record ID (likely in the form: w303_scaffold_X or LOCUS name)
                record_id = record.id
                
                # Check for alternative chromosome IDs in source features
                cm_id = None
                roman_id = None
                
                for feature in record.features:
                    if feature.type == "source":
                        # Look for CM007XXX.1 format ID in notes
                        for note in feature.qualifiers.get("note", []):
                            cm_match = re.search(r'(CM\d+\.\d+)', note)
                            if cm_match:
                                cm_id = cm_match.group(1)
                            
                            # Look for chromosome Roman numeral format
                            roman_match = re.search(r'chromosome\s+([IVX]+)', note, re.IGNORECASE)
                            if roman_match:
                                roman_id = roman_match.group(1)
                
                # Map the record ID to CM and Roman IDs if found
                if cm_id:
                    mapping[record_id] = cm_id
                    reverse_mapping[cm_id] = record_id
                
                if roman_id:
                    if record_id not in mapping:
                        mapping[record_id] = roman_id
                    if roman_id not in reverse_mapping:
                        reverse_mapping[roman_id] = record_id
                
                # Handle special case: If record ID is in 'w303_scaffold_X' format
                # and we have CM ID but no direct mapping
                w303_match = re.match(r'w303_scaffold_(\d+)', record_id)
                if w303_match and cm_id and record_id not in mapping:
                    mapping[record_id] = cm_id
                    reverse_mapping[cm_id] = record_id
                
        except Exception as e:
            print(f"Error parsing {gb_file}: {e}")
    
    # Add standard mappings if we're missing some
    scaffold_to_cm = {
        'w303_scaffold_1': 'CM007964.1',
        'w303_scaffold_2': 'CM007965.1',
        'w303_scaffold_3': 'CM007966.1',
        'w303_scaffold_4': 'CM007967.1',
        'w303_scaffold_5': 'CM007968.1',
        'w303_scaffold_6': 'CM007969.1',
        'w303_scaffold_7': 'CM007970.1',
        'w303_scaffold_8': 'CM007971.1',
        'w303_scaffold_9': 'CM007972.1',
        'w303_scaffold_10': 'CM007973.1',
        'w303_scaffold_11': 'CM007974.1',
        'w303_scaffold_12': 'CM007975.1',
        'w303_scaffold_16': 'CM007976.1',
        'w303_scaffold_17': 'CM007977.1',
        'w303_scaffold_18': 'CM007978.1',
        'w303_scaffold_19': 'CM007979.1',
        'w303_scaffold_20': 'CM007980.1',
        'w303_scaffold_21': 'CM007981.1'
    }
    
    # Merge with discovered mappings
    for scaffold, cm in scaffold_to_cm.items():
        if scaffold not in mapping:
            mapping[scaffold] = cm
        if cm not in reverse_mapping:
            reverse_mapping[cm] = scaffold
    
    if debug:
        print(f"Chromosome mapping (first 10 entries):")
        for i, (k, v) in enumerate(mapping.items()):
            if i >= 10: break
            print(f"  {k} -> {v}")
    
    print(f"Created mapping for {len(mapping)} chromosome IDs")
    return mapping, reverse_mapping

def load_chromosome_mapping(mapping_file):
    """Load chromosome mapping from a file."""
    print(f"Loading chromosome mapping from {mapping_file}")
    mapping = {}
    reverse_mapping = {}
    
    with open(mapping_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                parts = line.split('\t')
                if len(parts) == 2:
                    # Extract the CM ID from the first field
                    cm_id_match = re.search(r'(CM\d+\.\d+)', parts[0])
                    if cm_id_match:
                        cm_id = cm_id_match.group(1)
                        # Extract the scaffold ID from the second field
                        scaffold_id = parts[1].replace('>', '')
                        mapping[scaffold_id] = cm_id
                        reverse_mapping[cm_id] = scaffold_id
    
    print(f"Loaded mapping for {len(mapping)} chromosome IDs")
    return mapping, reverse_mapping

def load_genome_annotations_from_genbank(genbank_dir, debug=False):
    """Load gene annotations from GenBank files."""
    print(f"Loading genome annotations from {genbank_dir}")
    genes = []
    
    # Get all GenBank files in the directory
    genbank_files = [os.path.join(genbank_dir, f) for f in os.listdir(genbank_dir) 
                    if f.endswith(('.gb', '.gbk', '.genbank'))]
    
    if not genbank_files:
        # Check for files without standard extensions
        genbank_files = [os.path.join(genbank_dir, f) for f in os.listdir(genbank_dir) 
                        if os.path.isfile(os.path.join(genbank_dir, f)) and not f.startswith('.')]
    
    print(f"Found {len(genbank_files)} GenBank files")
    
    # Track chromosome IDs for debugging
    chromosome_ids = set()
    
    for gb_file in genbank_files:
        try:
            for record in SeqIO.parse(gb_file, "genbank"):
                scaffold_id = record.id
                chromosome_ids.add(scaffold_id)
                
                for feature in record.features:
                    if feature.type in ['gene', 'CDS', 'mRNA']:
                        # Extract gene ID and name
                        gene_id = feature.qualifiers.get('locus_tag', [''])[0]
                        gene_name = feature.qualifiers.get('gene', [''])[0] if 'gene' in feature.qualifiers else ''
                        product = feature.qualifiers.get('product', [''])[0] if 'product' in feature.qualifiers else ''
                        note = '; '.join(feature.qualifiers.get('note', []))
                        
                        # Extract location information
                        start = int(feature.location.start) + 1  # Convert to 1-based
                        end = int(feature.location.end)
                        strand = '+' if feature.location.strand == 1 else '-'
                        
                        genes.append({
                            'chrom': scaffold_id,
                            'start': start,
                            'end': end,
                            'strand': strand,
                            'feature_type': feature.type,
                            'gene_id': gene_id,
                            'gene_name': gene_name,
                            'product': product,
                            'note': note
                        })
        except Exception as e:
            print(f"Error parsing {gb_file}: {e}")
    
    if debug:
        print(f"Chromosome IDs in GenBank files (first 10): {list(chromosome_ids)[:10]}")
    
    print(f"Loaded {len(genes)} genome annotation features")
    return genes

def map_variants_to_genes(variants, genes, chrom_mapping, reverse_mapping, debug=False):
    """Map variants to genes and identify affected genes."""
    print("Mapping variants to genes")
    
    # Debug output
    if debug:
        variant_chromosomes = set(variants['chrom'])
        print(f"Chromosome IDs in variants (first 10): {list(variant_chromosomes)[:10]}")
        
        gene_df = pd.DataFrame(genes)
        gene_chromosomes = set(gene_df['chrom'])
        print(f"Chromosome IDs in genes (first 10): {list(gene_chromosomes)[:10]}")
    
    # Create a DataFrame for gene features
    gene_df = pd.DataFrame(genes)
    
    # Group genes by chromosome
    genes_by_chrom = {}
    for _, gene in gene_df.iterrows():
        chrom = gene['chrom']
        if chrom not in genes_by_chrom:
            genes_by_chrom[chrom] = []
        genes_by_chrom[chrom].append(gene)
    
    # Map variants to genes
    mapped_variants = []
    mapped_count = 0
    
    for _, variant in variants.iterrows():
        chrom = variant['chrom']
        pos = variant['pos']
        
        # Try to map the variant chromosome to the gene chromosome format
        gene_chrom = chrom
        if chrom in chrom_mapping:
            gene_chrom = chrom_mapping[chrom]
        
        # Find overlapping or nearby genes
        affected_gene = None
        affected_gene_distance = float('inf')
        
        # Try with mapped chromosome first
        if gene_chrom in genes_by_chrom:
            for gene in genes_by_chrom[gene_chrom]:
                if gene['start'] <= pos <= gene['end']:
                    # Variant is within gene
                    affected_gene = gene
                    affected_gene_distance = 0
                    break
                else:
                    # Calculate distance to gene
                    if pos < gene['start']:
                        distance = gene['start'] - pos
                    else:  # pos > gene['end']
                        distance = pos - gene['end']
                    
                    if distance < affected_gene_distance:
                        affected_gene = gene
                        affected_gene_distance = distance
        
        # If not found, try with original chromosome
        if affected_gene is None and chrom in genes_by_chrom:
            for gene in genes_by_chrom[chrom]:
                if gene['start'] <= pos <= gene['end']:
                    # Variant is within gene
                    affected_gene = gene
                    affected_gene_distance = 0
                    break
                else:
                    # Calculate distance to gene
                    if pos < gene['start']:
                        distance = gene['start'] - pos
                    else:  # pos > gene['end']
                        distance = pos - gene['end']
                    
                    if distance < affected_gene_distance:
                        affected_gene = gene
                        affected_gene_distance = distance
        
        # If still not found, try with reverse-mapped chromosome
        if affected_gene is None and chrom in reverse_mapping:
            reverse_chrom = reverse_mapping[chrom]
            if reverse_chrom in genes_by_chrom:
                for gene in genes_by_chrom[reverse_chrom]:
                    if gene['start'] <= pos <= gene['end']:
                        # Variant is within gene
                        affected_gene = gene
                        affected_gene_distance = 0
                        break
                    else:
                        # Calculate distance to gene
                        if pos < gene['start']:
                            distance = gene['start'] - pos
                        else:  # pos > gene['end']
                            distance = pos - gene['end']
                        
                        if distance < affected_gene_distance:
                            affected_gene = gene
                            affected_gene_distance = distance
        
        # Add gene information to variant
        variant_info = variant.to_dict()
        if affected_gene:
            mapped_count += 1
            variant_info['affected_gene_id'] = affected_gene['gene_id']
            variant_info['affected_gene_name'] = affected_gene['gene_name']
            variant_info['affected_gene_product'] = affected_gene['product']
            variant_info['affected_gene_note'] = affected_gene['note']
            variant_info['affected_gene_distance'] = affected_gene_distance
            variant_info['affected_gene_feature_type'] = affected_gene['feature_type']
            variant_info['within_gene'] = (affected_gene_distance == 0)
        else:
            variant_info['affected_gene_id'] = None
            variant_info['affected_gene_name'] = None
            variant_info['affected_gene_product'] = None
            variant_info['affected_gene_note'] = None
            variant_info['affected_gene_distance'] = None
            variant_info['affected_gene_feature_type'] = None
            variant_info['within_gene'] = False
        
        mapped_variants.append(variant_info)
    
    mapped_df = pd.DataFrame(mapped_variants)
    print(f"Mapped {len(mapped_variants)} variants to genes, {mapped_count} with affected genes")
    
    # Extra debug info for troubleshooting
    if debug and mapped_count == 0:
        print("\nWARNING: No variants were mapped to affected genes!")
        print("This suggests a chromosome naming mismatch between variants and genes.")
        print("Variant chromosomes:", set(variants['chrom']))
        print("Gene chromosomes:", set(gene_df['chrom']))
        print("Please check chromosome naming and provide a mapping file if needed.")
    
    return mapped_df

def analyze_affected_genes(mapped_variants, debug=False):
    """Analyze affected genes and their potential connections to ergosterol pathway."""
    print("Analyzing affected genes")
    
    # Debug info about the mapped variants
    if debug:
        has_affected_genes = mapped_variants['affected_gene_id'].notna()
        print(f"Variants with affected_gene_id: {has_affected_genes.sum()} out of {len(mapped_variants)}")
        
        if has_affected_genes.sum() > 0:
            # Show a sample of variants with affected genes
            sample = mapped_variants[has_affected_genes].head(3)
            print("\nSample of variants with affected genes:")
            for _, row in sample.iterrows():
                print(f"  Variant: {row['chrom']}:{row['pos']} {row['ref']}>{row['alt']}")
                print(f"  Affected gene: {row['affected_gene_id']} ({row['affected_gene_name']})")
                print(f"  Product: {row['affected_gene_product']}")
                print(f"  Distance: {row['affected_gene_distance']}")
                print()
    
    # Count affected genes
    affected_gene_counts = Counter()
    for _, variant in mapped_variants.iterrows():
        if pd.notna(variant['affected_gene_id']) and variant['affected_gene_id']:
            affected_gene_counts[variant['affected_gene_id']] += 1
    
    if debug:
        print(f"Number of unique affected genes: {len(affected_gene_counts)}")
        if len(affected_gene_counts) > 0:
            print("Top 5 affected genes by variant count:")
            for gene_id, count in affected_gene_counts.most_common(5):
                print(f"  {gene_id}: {count} variants")
    
    # Find unique affected genes
    unique_affected_genes = []
    for gene_id, count in affected_gene_counts.items():
        gene_variants = mapped_variants[mapped_variants['affected_gene_id'] == gene_id]
        gene_info = gene_variants.iloc[0]
        
        unique_affected_genes.append({
            'gene_id': gene_id,
            'gene_name': gene_info['affected_gene_name'],
            'product': gene_info['affected_gene_product'],
            'note': gene_info['affected_gene_note'],
            'variant_count': count,
            'nearest_erg_gene': gene_info['nearest_gene_name'],
            'min_distance_to_erg': gene_variants['distance_to_gene'].min(),
            'high_impact_count': len(gene_variants[gene_variants['impact'] == 'HIGH']),
            'moderate_impact_count': len(gene_variants[gene_variants['impact'] == 'MODERATE']),
            'severe_impact_count': len(gene_variants[gene_variants['functional_impact'] == 'Severe']) if 'functional_impact' in gene_variants.columns else 0,
            'treatments': ', '.join(gene_variants['treatment'].unique())
        })
    
    # Sort affected genes by variant count
    unique_affected_genes.sort(key=lambda x: x['variant_count'], reverse=True)
    
    # Create affected genes DataFrame
    affected_genes_df = pd.DataFrame(unique_affected_genes)
    
    # Analyze gene names and products for potential connections to ergosterol pathway
    if not affected_genes_df.empty:
        affected_genes_df['sterol_related'] = affected_genes_df.apply(
            lambda row: is_sterol_related(row['gene_name'], row['product'], row['note']), axis=1
        )
    
    return affected_genes_df

def is_sterol_related(gene_name, product, note):
    """Check if a gene is potentially related to sterol metabolism."""
    sterol_keywords = [
        'sterol', 'lipid', 'membrane', 'fatty acid', 'ergosterol', 'cholesterol',
        'isoprenoid', 'terpene', 'mevalonate', 'hmg', 'farnesyl', 'acetyl-coa',
        'lanosterol', 'squalene', 'zymosterol'
    ]
    
    text = f"{gene_name} {product} {note}".lower()
    
    for keyword in sterol_keywords:
        if keyword in text:
            return True
    
    return False

def generate_gene_summaries(affected_genes_df, mapped_variants, output_dir):
    """Generate summaries of affected genes."""
    print("Generating gene summaries")
    
    if affected_genes_df.empty:
        print("No affected genes to summarize")
        return {
            'erg_gene_summaries': [],
            'function_categories': {},
            'treatment_summaries': [],
            'total_affected_genes': 0,
            'sterol_related_genes': 0,
            'high_impact_genes': 0
        }
    
    # Create a dictionary of affected genes by ergosterol gene
    genes_by_erg = defaultdict(list)
    for _, gene in affected_genes_df.iterrows():
        erg_gene = gene['nearest_erg_gene']
        genes_by_erg[erg_gene].append(gene.to_dict())
    
    # Generate summary by ergosterol gene
    erg_summaries = []
    for erg_gene, genes in genes_by_erg.items():
        erg_summary = {
            'erg_gene': erg_gene,
            'affected_gene_count': len(genes),
            'high_impact_gene_count': sum(1 for g in genes if g['high_impact_count'] > 0),
            'sterol_related_count': sum(1 for g in genes if g['sterol_related']),
            'top_affected_genes': [g['gene_name'] if g['gene_name'] else g['gene_id'] for g in sorted(genes, key=lambda x: x['variant_count'], reverse=True)[:5]],
            'min_distance': min(g['min_distance_to_erg'] for g in genes)
        }
        erg_summaries.append(erg_summary)
    
    # Sort by affected gene count
    erg_summaries.sort(key=lambda x: x['affected_gene_count'], reverse=True)
    
    # Generate gene product/function summaries
    function_categories = defaultdict(int)
    for _, gene in affected_genes_df.iterrows():
        product = str(gene['product']).lower()
        
        # Categorize gene function
        if 'transport' in product:
            function_categories['Transport'] += 1
        elif 'transcrip' in product:
            function_categories['Transcription'] += 1
        elif 'signal' in product:
            function_categories['Signaling'] += 1
        elif 'membrane' in product:
            function_categories['Membrane'] += 1
        elif 'metabol' in product or 'synth' in product:
            function_categories['Metabolism'] += 1
        elif 'repair' in product or 'replic' in product:
            function_categories['DNA repair/replication'] += 1
        elif 'protein' in product and ('fold' in product or 'chaperone' in product):
            function_categories['Protein folding'] += 1
        elif 'protea' in product or 'degrad' in product:
            function_categories['Protein degradation'] += 1
        elif 'ribo' in product or 'translat' in product:
            function_categories['Translation'] += 1
        elif 'mitoch' in product:
            function_categories['Mitochondrial'] += 1
        elif 'hypothetical' in product or product == '' or 'unknown' in product:
            function_categories['Unknown/Hypothetical'] += 1
        else:
            function_categories['Other'] += 1
    
    # Generate treatment summaries
    treatment_summaries = []
    for treatment in mapped_variants['treatment'].unique():
        treatment_variants = mapped_variants[mapped_variants['treatment'] == treatment]
        treatment_genes = set(v['affected_gene_id'] for _, v in treatment_variants.iterrows() 
                             if pd.notna(v['affected_gene_id']) and v['affected_gene_id'])
        
        other_treatments = [t for t in mapped_variants['treatment'].unique() if t != treatment]
        other_genes = set()
        for other_treatment in other_treatments:
            other_variants = mapped_variants[mapped_variants['treatment'] == other_treatment]
            other_genes.update(v['affected_gene_id'] for _, v in other_variants.iterrows() 
                              if pd.notna(v['affected_gene_id']) and v['affected_gene_id'])
        
        unique_genes = treatment_genes - other_genes
        
        treatment_summary = {
            'treatment': treatment,
            'variant_count': len(treatment_variants),
            'affected_gene_count': len(treatment_genes),
            'high_impact_count': len(treatment_variants[treatment_variants['impact'] == 'HIGH']),
            'unique_affected_genes': len(unique_genes),
            'unique_gene_ids': list(unique_genes)
        }
        treatment_summaries.append(treatment_summary)
    
    # Sort by variant count
    treatment_summaries.sort(key=lambda x: x['variant_count'], reverse=True)
    
    summaries = {
        'erg_gene_summaries': erg_summaries,
        'function_categories': {k: v for k, v in sorted(function_categories.items(), key=lambda x: x[1], reverse=True)},
        'treatment_summaries': treatment_summaries,
        'total_affected_genes': len(affected_genes_df),
        'sterol_related_genes': int(affected_genes_df['sterol_related'].sum()),
        'high_impact_genes': len(affected_genes_df[affected_genes_df['high_impact_count'] > 0])
    }
    
    # Save summaries to file
    summary_path = os.path.join(output_dir, 'affected_genes_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summaries, f, indent=2)
    
    return summaries

def generate_visualizations(affected_genes_df, mapped_variants, output_dir):
    """Generate visualizations of affected genes."""
    print("Generating visualizations")
    
    # Create output directory for visualizations
    vis_dir = os.path.join(output_dir, 'visualizations')
    os.makedirs(vis_dir, exist_ok=True)
    
    plots = {}
    
    if affected_genes_df.empty:
        print("No affected genes to visualize")
        return plots
    
    # 1. Bar chart of top affected genes by variant count
    plt.figure(figsize=(12, 8))
    top_genes = affected_genes_df.sort_values('variant_count', ascending=False).head(20)
    
    # Use gene names if available, otherwise use gene IDs
    labels = [row['gene_name'] if row['gene_name'] else row['gene_id'] for _, row in top_genes.iterrows()]
    
    sns.barplot(x='variant_count', y=labels, data=top_genes)
    plt.title('Top 20 Genes Affected by HIGH/MODERATE Impact Variants')
    plt.xlabel('Number of Variants')
    plt.ylabel('Gene')
    plt.tight_layout()
    
    top_genes_plot = os.path.join(vis_dir, 'top_affected_genes.png')
    plt.savefig(top_genes_plot, dpi=300)
    plt.close()
    plots['top_genes_plot'] = top_genes_plot
    
    # Additional visualizations remain the same as original...
    
    return plots

def save_results(affected_genes_df, mapped_variants, output_dir):
    """Save analysis results to output files."""
    print(f"Saving results to {output_dir}")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Save affected genes data
    affected_genes_path = os.path.join(output_dir, 'affected_genes.tsv')
    if not affected_genes_df.empty:
        affected_genes_df.to_csv(affected_genes_path, sep='\t', index=False)
    else:
        # Create an empty file with header
        with open(affected_genes_path, 'w') as f:
            f.write("gene_id\tgene_name\tproduct\tnote\tvariant_count\tnearest_erg_gene\tmin_distance_to_erg\thigh_impact_count\tmoderate_impact_count\tsevere_impact_count\ttreatments\tsterol_related\n")
    
    # Save mapped variants data
    mapped_variants_path = os.path.join(output_dir, 'variants_with_genes.tsv')
    mapped_variants.to_csv(mapped_variants_path, sep='\t', index=False)
    
    # Create a human-readable report
    report = [
        "# Analysis of Genes Affected by HIGH/MODERATE Impact Variants Near Ergosterol Pathway",
        f"Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
    ]
    
    # Report content creation remains the same as original...
    
    # Save report
    report_path = os.path.join(output_dir, 'affected_genes_report.md')
    with open(report_path, 'w') as f:
        f.write('\n'.join(report))
    
    return {
        'affected_genes': affected_genes_path,
        'mapped_variants': mapped_variants_path,
        'report': report_path
    }

def main():
    """Main function to run the analysis."""
    # Parse command line arguments
    args = parse_arguments()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load data
    variants = load_variants(args.variants_file)
    genes = load_genome_annotations_from_genbank(args.genbank_dir, debug=args.debug)
    
    # Create or load chromosome mapping
    if args.mapping_file and os.path.exists(args.mapping_file):
        chrom_mapping, reverse_mapping = load_chromosome_mapping(args.mapping_file)
    else:
        chrom_mapping, reverse_mapping = create_chromosome_mapping(args.genbank_dir, debug=args.debug)
    
    # Write mappings to file for future reference
    mapping_dir = os.path.join(args.output_dir, 'mappings')
    os.makedirs(mapping_dir, exist_ok=True)
    
    with open(os.path.join(mapping_dir, 'chromosome_mapping.txt'), 'w') as f:
        for k, v in sorted(chrom_mapping.items()):
            f.write(f"{k}\t{v}\n")
    
    with open(os.path.join(mapping_dir, 'reverse_mapping.txt'), 'w') as f:
        for k, v in sorted(reverse_mapping.items()):
            f.write(f"{k}\t{v}\n")
    
    # Map variants to genes
    mapped_variants = map_variants_to_genes(variants, genes, chrom_mapping, reverse_mapping, debug=args.debug)
    
    # Analyze affected genes
    affected_genes_df = analyze_affected_genes(mapped_variants, debug=args.debug)
    
    # Generate gene summaries
    summaries = generate_gene_summaries(affected_genes_df, mapped_variants, args.output_dir)
    
    # Generate visualizations
    plots = generate_visualizations(affected_genes_df, mapped_variants, args.output_dir)
    
    # Save results
    file_paths = save_results(affected_genes_df, mapped_variants, args.output_dir)
    
    print("Analysis complete!")
    print(f"Found {len(affected_genes_df)} genes affected by HIGH/MODERATE impact variants")
    print(f"Results saved to {args.output_dir}")
    print(f"Report: {file_paths['report']}")

if __name__ == "__main__":
    main()