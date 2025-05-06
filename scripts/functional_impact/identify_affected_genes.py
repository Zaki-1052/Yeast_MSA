#!/usr/bin/env python3
"""
identify_affected_genes.py

This script identifies and characterizes genes affected by HIGH/MODERATE impact variants 
near ergosterol pathway genes, with fixes for GenBank annotation parsing.
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

def load_tabular_mapping(mapping_file):
    """
    Load chromosome mapping from a tabular file in format:
    chromosome_id  w303_scaffold
    CM007964.1     w303_scaffold_1
    ...
    """
    print(f"Loading chromosome mapping from {mapping_file}")
    try:
        mapping = pd.read_csv(mapping_file, sep='\t')
        
        # Create the mappings
        cm_to_scaffold = dict(zip(mapping['chromosome_id'], mapping['w303_scaffold']))
        scaffold_to_cm = dict(zip(mapping['w303_scaffold'], mapping['chromosome_id']))
        
        print(f"Loaded mapping for {len(cm_to_scaffold)} chromosome IDs")
        print(f"Sample mappings (first 3):")
        for i, (k, v) in enumerate(cm_to_scaffold.items()):
            if i >= 3: break
            print(f"  {k} -> {v}")
            
        return scaffold_to_cm, cm_to_scaffold
        
    except Exception as e:
        print(f"Error loading mapping file: {e}")
        print("Continuing with empty mapping...")
        return {}, {}

def extract_cm_id_from_genbank(record):
    """
    Extract the CM identifier from a GenBank record.
    Look in the /note field of source features for CM007XXX.1 format IDs.
    """
    # Look for source features
    for feature in record.features:
        if feature.type == "source":
            # Check the note qualifiers
            for note in feature.qualifiers.get("note", []):
                # Look for CM007XXX.1 pattern
                cm_match = re.search(r'(CM\d+\.\d+)', note)
                if cm_match:
                    return cm_match.group(1)
    
    # If not found, return None
    return None

def extract_gene_id(feature, debug=False):
    """
    Extract gene ID from GenBank feature, trying multiple qualifier fields.
    
    Args:
        feature: BioPython SeqFeature object
        debug: Whether to print debug info
        
    Returns:
        str: Gene ID
    """
    # Try locus_tag first (more specific)
    if 'locus_tag' in feature.qualifiers and feature.qualifiers['locus_tag']:
        return feature.qualifiers['locus_tag'][0]
    
    # Try gene next (more commonly present)
    if 'gene' in feature.qualifiers and feature.qualifiers['gene']:
        return feature.qualifiers['gene'][0]
    
    # If no gene ID found
    if debug:
        print(f"Warning: No gene ID found for feature at {feature.location}")
        print(f"Available qualifiers: {list(feature.qualifiers.keys())}")
    
    return ""

def load_genome_annotations_from_genbank(genbank_dir, scaffold_to_cm=None, debug=False):
    """Load gene annotations from GenBank files with proper chromosome ID extraction."""
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
    
    # Create a mapping to track locus -> scaffold relationships
    locus_to_scaffold = {}
    
    # First pass: extract chromosome/scaffold IDs
    for gb_file in genbank_files:
        try:
            for record in SeqIO.parse(gb_file, "genbank"):
                # Get locus name from LOCUS line
                locus_name = record.name
                
                # Extract CM ID from source feature notes
                cm_id = extract_cm_id_from_genbank(record)
                
                if cm_id and locus_name:
                    # This maps the locus name to the CM ID
                    locus_to_scaffold[locus_name] = cm_id
                    
                    if debug:
                        print(f"Mapped locus {locus_name} to {cm_id}")
        except Exception as e:
            if debug:
                print(f"Error in first pass parsing {gb_file}: {e}")
    
    if debug:
        print(f"Extracted {len(locus_to_scaffold)} locus mappings")
        if locus_to_scaffold:
            print("Sample mappings (first 3):")
            for i, (k, v) in enumerate(locus_to_scaffold.items()):
                if i >= 3: break
                print(f"  {k} -> {v}")
    
    # Second pass: extract gene annotations using correct chromosome IDs
    for gb_file in genbank_files:
        try:
            for record in SeqIO.parse(gb_file, "genbank"):
                # Get locus name from LOCUS line
                locus_name = record.name
                
                # Use our mapping to get proper chromosome ID
                chromosome_id = locus_name  # default to locus name
                
                if scaffold_to_cm and locus_name in scaffold_to_cm:
                    # If we have a mapping from the mapping file, use it
                    chromosome_id = scaffold_to_cm[locus_name]
                elif locus_name in locus_to_scaffold:
                    # If we extracted a CM ID in the first pass, use it
                    chromosome_id = locus_to_scaffold[locus_name]
                
                for feature in record.features:
                    if feature.type in ['gene', 'CDS', 'mRNA']:
                        # Extract gene ID using our improved function
                        gene_id = extract_gene_id(feature, debug)
                        
                        # Extract other annotation info
                        gene_name = feature.qualifiers.get('gene', [''])[0] if 'gene' in feature.qualifiers else ''
                        product = feature.qualifiers.get('product', [''])[0] if 'product' in feature.qualifiers else ''
                        note = '; '.join(feature.qualifiers.get('note', []))
                        
                        # Skip features without a gene ID
                        if not gene_id:
                            continue
                        
                        # Extract location information
                        start = int(feature.location.start) + 1  # Convert to 1-based
                        end = int(feature.location.end)
                        strand = '+' if feature.location.strand == 1 else '-'
                        
                        genes.append({
                            'chrom': locus_name,  # Use the original locus name (w303_scaffold_X)
                            'chrom_id': chromosome_id,  # Store CM ID for reference
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
            if debug:
                print(f"Error in second pass parsing {gb_file}: {e}")
    
    # Track unique chromosome IDs for debugging
    chromosome_ids = set(g['chrom'] for g in genes)
    
    if debug:
        print(f"Chromosome IDs in GenBank files (first 10): {list(chromosome_ids)[:10]}")
        
        # Print some sample gene entries
        if genes:
            print("\nSample gene entries:")
            for i, gene in enumerate(genes[:3]):
                print(f"Gene {i+1}:")
                print(f"  Chromosome: {gene['chrom']}")
                print(f"  Gene ID: {gene['gene_id']}")
                print(f"  Gene Name: {gene['gene_name']}")
                print(f"  Location: {gene['start']}-{gene['end']} ({gene['strand']})")
                print(f"  Product: {gene['product']}")
                print()
    
    print(f"Loaded {len(genes)} genome annotation features")
    return genes

def map_variants_to_genes(variants, genes, debug=False):
    """Map variants to genes and identify affected genes."""
    print("Mapping variants to genes")
    
    # Debug output
    if debug:
        variant_chromosomes = set(variants['chrom'])
        print(f"Variant chromosomes (first 10): {list(variant_chromosomes)[:10]}")
        
        gene_chromosomes = set(g['chrom'] for g in genes)
        print(f"Gene chromosomes (first 10): {list(gene_chromosomes)[:10]}")
    
    # Group genes by chromosome
    genes_by_chrom = defaultdict(list)
    for gene in genes:
        chrom = gene['chrom']
        genes_by_chrom[chrom].append(gene)
    
    # Map variants to genes
    mapped_variants = []
    mapped_count = 0
    
    for _, variant in variants.iterrows():
        chrom = variant['chrom']
        pos = variant['pos']
        
        # Find overlapping or nearby genes
        affected_gene = None
        affected_gene_distance = float('inf')
        
        if chrom in genes_by_chrom:
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
        
        if variants.shape[0] > 0:
            print(f"Sample variant 1: {variants.iloc[0]['chrom']}:{variants.iloc[0]['pos']}")
            
        if len(genes) > 0:
            print(f"Sample gene 1: {genes[0]['chrom']}:{genes[0]['start']}-{genes[0]['end']}")
        
        print(f"Variant chromosomes: {set(variants['chrom'])}")
        print(f"Gene chromosomes: {set(g['chrom'] for g in genes)}")
    
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
        if pd.notna(variant['affected_gene_id']) and variant['affected_gene_id'] and variant['affected_gene_id'].strip():
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
                             if pd.notna(v['affected_gene_id']) and v['affected_gene_id'] and v['affected_gene_id'].strip())
        
        other_treatments = [t for t in mapped_variants['treatment'].unique() if t != treatment]
        other_genes = set()
        for other_treatment in other_treatments:
            other_variants = mapped_variants[mapped_variants['treatment'] == other_treatment]
            other_genes.update(v['affected_gene_id'] for _, v in other_variants.iterrows() 
                              if pd.notna(v['affected_gene_id']) and v['affected_gene_id'] and v['affected_gene_id'].strip())
        
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
    
    # 2. Pie chart of gene functions
    function_counts = Counter()
    for _, gene in affected_genes_df.iterrows():
        product = str(gene['product']).lower()
        
        # Categorize gene function
        if 'transport' in product:
            function_counts['Transport'] += 1
        elif 'transcrip' in product:
            function_counts['Transcription'] += 1
        elif 'signal' in product:
            function_counts['Signaling'] += 1
        elif 'membrane' in product:
            function_counts['Membrane'] += 1
        elif 'metabol' in product or 'synth' in product:
            function_counts['Metabolism'] += 1
        elif 'repair' in product or 'replic' in product:
            function_counts['DNA repair/replication'] += 1
        elif 'protein' in product and ('fold' in product or 'chaperone' in product):
            function_counts['Protein folding'] += 1
        elif 'protea' in product or 'degrad' in product:
            function_counts['Protein degradation'] += 1
        elif 'ribo' in product or 'translat' in product:
            function_counts['Translation'] += 1
        elif 'mitoch' in product:
            function_counts['Mitochondrial'] += 1
        elif 'hypothetical' in product or product == '' or 'unknown' in product:
            function_counts['Unknown/Hypothetical'] += 1
        else:
            function_counts['Other'] += 1
    
    if function_counts:
        plt.figure(figsize=(10, 10))
        plt.pie(
            function_counts.values(), 
            labels=function_counts.keys(), 
            autopct='%1.1f%%',
            colors=sns.color_palette('pastel', len(function_counts))
        )
        plt.title('Functional Categories of Affected Genes')
        plt.tight_layout()
        
        function_plot = os.path.join(vis_dir, 'gene_functions.png')
        plt.savefig(function_plot, dpi=300)
        plt.close()
        plots['function_plot'] = function_plot
    
    # 3. Heatmap of gene-by-treatment
    if 'affected_gene_name' in mapped_variants.columns:
        # Create a copy to avoid SettingWithCopyWarning
        mapped_variants_copy = mapped_variants.copy()
        
        # Fill NaN affected_gene_name with affected_gene_id
        mapped_variants_copy.loc[mapped_variants_copy['affected_gene_name'].isna(), 'affected_gene_name'] = \
            mapped_variants_copy.loc[mapped_variants_copy['affected_gene_name'].isna(), 'affected_gene_id']
        
        # Count variants by gene and treatment
        gene_treatment_counts = mapped_variants_copy.groupby(['affected_gene_name', 'treatment']).size().unstack(fill_value=0)
        
        if not gene_treatment_counts.empty:
            # Select top genes
            top_genes_for_heatmap = affected_genes_df.sort_values('variant_count', ascending=False).head(15)['gene_name']
            if 'gene_name' in top_genes_for_heatmap:
                filtered_heatmap = gene_treatment_counts.loc[gene_treatment_counts.index.isin(top_genes_for_heatmap)]
                
                if not filtered_heatmap.empty:
                    plt.figure(figsize=(12, 10))
                    sns.heatmap(filtered_heatmap, annot=True, cmap='YlGnBu', fmt='d')
                    plt.title('Number of Variants by Gene and Treatment')
                    plt.ylabel('Gene')
                    plt.xlabel('Treatment')
                    plt.tight_layout()
                    
                    heatmap_plot = os.path.join(vis_dir, 'gene_treatment_heatmap.png')
                    plt.savefig(heatmap_plot, dpi=300)
                    plt.close()
                    plots['heatmap_plot'] = heatmap_plot
    
    # 4. Bar chart of gene count by nearest ergosterol gene
    plt.figure(figsize=(12, 8))
    erg_gene_counts = affected_genes_df['nearest_erg_gene'].value_counts()
    
    sns.barplot(x=erg_gene_counts.index, y=erg_gene_counts.values)
    plt.title('Affected Genes by Nearest Ergosterol Gene')
    plt.xlabel('Ergosterol Gene')
    plt.ylabel('Number of Affected Genes')
    plt.xticks(rotation=45)
    plt.tight_layout()
    
    erg_gene_plot = os.path.join(vis_dir, 'affected_genes_by_erg.png')
    plt.savefig(erg_gene_plot, dpi=300)
    plt.close()
    plots['erg_gene_plot'] = erg_gene_plot
    
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
    
    if affected_genes_df.empty:
        report.extend([
            "## Overview",
            "No genes were found to be affected by the HIGH/MODERATE impact variants.",
            "This could be due to the variants being in intergenic regions or due to incomplete gene annotations.",
            "",
            "## Biological Implications",
            "",
            "The absence of identifiable genes affected by these variants suggests that the HIGH/MODERATE impact",
            "variants near ergosterol pathway genes might be affecting intergenic regions or unannotated genes.",
            "This is consistent with our previous finding that the ergosterol pathway itself is under strong",
            "purifying selection, extending to a considerable distance around these genes."
        ])
    else:
        report.extend([
            "## Overview",
            f"Total affected genes: {len(affected_genes_df)}",
            f"Sterol-related genes: {int(affected_genes_df['sterol_related'].sum())}",
            f"Genes with HIGH impact variants: {len(affected_genes_df[affected_genes_df['high_impact_count'] > 0])}",
            "",
            "## Top Affected Genes",
        ])
        
        # Add top genes to report
        top_genes = affected_genes_df.sort_values('variant_count', ascending=False).head(15)
        for _, gene in top_genes.iterrows():
            gene_name = gene['gene_name'] if gene['gene_name'] else gene['gene_id']
            report.append(f"### {gene_name}")
            report.append(f"- Product: {gene['product']}")
            report.append(f"- Variant count: {gene['variant_count']}")
            report.append(f"- Nearest ergosterol gene: {gene['nearest_erg_gene']}")
            report.append(f"- Distance to ergosterol gene: {gene['min_distance_to_erg']} bp")
            report.append(f"- HIGH impact variants: {gene['high_impact_count']}")
            report.append(f"- MODERATE impact variants: {gene['moderate_impact_count']}")
            report.append(f"- Found in treatments: {gene['treatments']}")
            report.append(f"- Potential sterol relation: {'Yes' if gene['sterol_related'] else 'No'}")
            report.append("")
        
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
    
    # Load variant data
    variants = load_variants(args.variants_file)
    
    # Load mapping file if provided
    scaffold_to_cm = {}
    cm_to_scaffold = {}
    if args.mapping_file and os.path.exists(args.mapping_file):
        scaffold_to_cm, cm_to_scaffold = load_tabular_mapping(args.mapping_file)
    
    # Load genome annotations from GenBank files with proper chromosome ID extraction
    genes = load_genome_annotations_from_genbank(args.genbank_dir, scaffold_to_cm, debug=args.debug)
    
    # Map variants to genes
    mapped_variants = map_variants_to_genes(variants, genes, debug=args.debug)
    
    # Analyze affected genes
    affected_genes_df = analyze_affected_genes(mapped_variants, debug=args.debug)
    
    # Generate gene summaries
    summaries = generate_gene_summaries(affected_genes_df, mapped_variants, args.output_dir)
    
    # Generate visualizations
    plots = generate_visualizations(affected_genes_df, mapped_variants, args.output_dir)
    
    # Save results
    file_paths = save_results(affected_genes_df, mapped_variants, args.output_dir)
    
    # Output results
    print("\nAnalysis complete!")
    print(f"Found {len(affected_genes_df)} genes affected by HIGH/MODERATE impact variants")
    print(f"Results saved to {args.output_dir}")
    print(f"Report: {file_paths['report']}")

if __name__ == "__main__":
    main()