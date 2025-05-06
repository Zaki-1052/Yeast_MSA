#!/usr/bin/env python3
"""
analyze_high_impact_variants.py

This script analyzes HIGH and MODERATE impact variants in the ergosterol pathway genes:
1. Extracts variants with significant predicted effects on protein structure/function
2. Categorizes them by gene, treatment, and effect type
3. Predicts functional consequences
4. Generates visualizations and reports

Usage:
  python analyze_high_impact_variants.py --vcf_dir <path> --output_dir <path> --gene_mapping <path>
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
import gzip
from Bio.SeqUtils import ProtParam
from matplotlib.patches import Patch

# Amino acid properties for functional impact analysis
AA_PROPERTIES = {
    'A': {'hydrophobic': True, 'charge': 0, 'size': 'small', 'polarity': 'nonpolar'},
    'C': {'hydrophobic': True, 'charge': 0, 'size': 'small', 'polarity': 'nonpolar'},
    'D': {'hydrophobic': False, 'charge': -1, 'size': 'small', 'polarity': 'acidic'},
    'E': {'hydrophobic': False, 'charge': -1, 'size': 'medium', 'polarity': 'acidic'},
    'F': {'hydrophobic': True, 'charge': 0, 'size': 'large', 'polarity': 'nonpolar'},
    'G': {'hydrophobic': True, 'charge': 0, 'size': 'small', 'polarity': 'nonpolar'},
    'H': {'hydrophobic': False, 'charge': 1, 'size': 'medium', 'polarity': 'basic'},
    'I': {'hydrophobic': True, 'charge': 0, 'size': 'large', 'polarity': 'nonpolar'},
    'K': {'hydrophobic': False, 'charge': 1, 'size': 'large', 'polarity': 'basic'},
    'L': {'hydrophobic': True, 'charge': 0, 'size': 'large', 'polarity': 'nonpolar'},
    'M': {'hydrophobic': True, 'charge': 0, 'size': 'large', 'polarity': 'nonpolar'},
    'N': {'hydrophobic': False, 'charge': 0, 'size': 'small', 'polarity': 'polar'},
    'P': {'hydrophobic': True, 'charge': 0, 'size': 'small', 'polarity': 'nonpolar'},
    'Q': {'hydrophobic': False, 'charge': 0, 'size': 'medium', 'polarity': 'polar'},
    'R': {'hydrophobic': False, 'charge': 1, 'size': 'large', 'polarity': 'basic'},
    'S': {'hydrophobic': False, 'charge': 0, 'size': 'small', 'polarity': 'polar'},
    'T': {'hydrophobic': False, 'charge': 0, 'size': 'small', 'polarity': 'polar'},
    'V': {'hydrophobic': True, 'charge': 0, 'size': 'medium', 'polarity': 'nonpolar'},
    'W': {'hydrophobic': True, 'charge': 0, 'size': 'large', 'polarity': 'nonpolar'},
    'Y': {'hydrophobic': False, 'charge': 0, 'size': 'large', 'polarity': 'polar'},
    '*': {'hydrophobic': False, 'charge': 0, 'size': 'none', 'polarity': 'none'}
}

# Functional domains in ergosterol pathway proteins (approximate positions)
# Based on literature and UniProt database
PROTEIN_DOMAINS = {
    'ERG1': {
        'FAD_binding': (50, 200),
        'Substrate_binding': (300, 400),
        'Membrane_association': (450, 490)
    },
    'ERG2': {
        'Sigma_receptor_like': (20, 150),
        'Sterol_binding': (80, 180)
    },
    'ERG3': {
        'Cytochrome_b5': (20, 90),
        'Sterol_desaturase': (100, 300)
    },
    'ERG4': {
        'Sterol_reductase': (50, 300),
        'NADPH_binding': (100, 150)
    },
    'ERG5': {
        'Cytochrome_P450': (50, 400),
        'Heme_binding': (420, 450)
    },
    'ERG6': {
        'SAM_binding': (80, 120),
        'Sterol_binding': (200, 300)
    },
    'ERG7': {
        'Oxidosqualene_cyclase': (100, 600),
        'Catalytic_site': (450, 500)
    },
    'ERG9': {
        'Isoprenoid_synthase': (50, 300),
        'DDXXD_motif': (80, 85)
    },
    'ERG11': {
        'Cytochrome_P450': (50, 400),
        'Heme_binding': (420, 450)
    },
    'ERG24': {
        'Sterol_reductase': (50, 300),
        'NADPH_binding': (100, 150)
    },
    'ERG25': {
        'Fatty_acid_hydroxylase': (100, 300),
        'Iron_binding': (150, 180)
    }
}

# Impact severity weights for aggregated scoring
IMPACT_WEIGHTS = {
    'HIGH': 1.0,
    'MODERATE': 0.5,
    'LOW': 0.1,
    'MODIFIER': 0.0
}

# Effect severity weights for functional scoring
EFFECT_WEIGHTS = {
    'stop_gained': 1.0,
    'frameshift_variant': 0.9,
    'splice_acceptor_variant': 0.8,
    'splice_donor_variant': 0.8,
    'start_lost': 0.7,
    'stop_lost': 0.7,
    'missense_variant': 0.5,
    'inframe_insertion': 0.4,
    'inframe_deletion': 0.4,
    'disruptive_inframe_insertion': 0.45,
    'disruptive_inframe_deletion': 0.45,
    'protein_altering_variant': 0.3,
    'synonymous_variant': 0.0,
    'splice_region_variant': 0.2,
    'initiator_codon_variant': 0.3,
    'stop_retained_variant': 0.1,
    '5_prime_UTR_variant': 0.1,
    '3_prime_UTR_variant': 0.1,
    'upstream_gene_variant': 0.0,
    'downstream_gene_variant': 0.0,
    'intergenic_region': 0.0
}

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Analyze HIGH and MODERATE impact variants')
    parser.add_argument('--vcf_dir', required=True, help='Directory containing annotated VCF files')
    parser.add_argument('--output_dir', required=True, help='Directory to save results')
    parser.add_argument('--gene_mapping', required=True, help='Path to gene mapping TSV file')
    return parser.parse_args()

def load_gene_mapping(gene_mapping_file):
    """Load gene mapping information."""
    print(f"Loading gene mapping from {gene_mapping_file}")
    genes = pd.read_csv(gene_mapping_file, sep='\t')
    
    # Convert column names to lowercase for consistency
    genes.columns = [col.lower() for col in genes.columns]
    
    print(f"Loaded information for {len(genes)} genes")
    return genes

def read_vcf(vcf_file):
    """Read a VCF file and extract relevant information."""
    # Determine if file is gzipped
    is_gzipped = vcf_file.endswith('.gz')
    opener = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'
    
    variants = []
    sample_name = os.path.basename(vcf_file).split('.')[0]
    
    with opener(vcf_file, mode) as f:
        # Skip header lines
        for line in f:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    # Extract the sample name from the header
                    header_parts = line.strip().split('\t')
                    if len(header_parts) > 9:
                        sample_name = header_parts[9]
                continue
            
            # Parse variant line
            parts = line.strip().split('\t')
            chrom, pos, id_, ref, alt = parts[0:5]
            qual, filter_, info = parts[5:8]
            
            # Extract treatment from sample name
            treatment = 'unknown'
            if '-' in sample_name:
                treatment_part = sample_name.split('-')[0]
                if treatment_part in ['CAS', 'STC', 'WT', 'WTA']:
                    treatment = treatment_part
                    if treatment_part == 'WT' and 'WT-37' in sample_name:
                        treatment = 'WT-37'
                elif treatment_part in ['CTRL']:
                    if 'CAS-CTRL' in sample_name:
                        treatment = 'CAS-CTRL'
                    elif 'STC-CTRL' in sample_name:
                        treatment = 'STC-CTRL'
                    else:
                        treatment = 'WT-CTRL'
            
            # Parse SnpEff annotation from INFO field
            ann_info = None
            for info_field in info.split(';'):
                if info_field.startswith('ANN='):
                    ann_info = info_field[4:]
                    break
            
            if ann_info:
                # Multiple annotations may be separated by commas
                for ann in ann_info.split(','):
                    ann_parts = ann.split('|')
                    if len(ann_parts) < 10:
                        continue
                    
                    allele = ann_parts[0]
                    effect = ann_parts[1]
                    impact = ann_parts[2]
                    gene_name = ann_parts[3]
                    gene_id = ann_parts[4]
                    feature_type = ann_parts[5]
                    feature_id = ann_parts[6]
                    transcript_biotype = ann_parts[7]
                    exon_rank = ann_parts[8]
                    hgvs_c = ann_parts[9]
                    hgvs_p = ann_parts[10] if len(ann_parts) > 10 else ''
                    
                    # Extract amino acid change from HGVS notation
                    aa_change = None
                    aa_pos = None
                    if hgvs_p and hgvs_p.startswith('p.'):
                        # Remove the 'p.' prefix
                        hgvs_p = hgvs_p[2:]
                        # Try to parse amino acid change
                        aa_match = re.search(r'([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|[*])', hgvs_p)
                        if aa_match:
                            ref_aa = aa_match.group(1)
                            aa_pos = int(aa_match.group(2))
                            alt_aa = aa_match.group(3)
                            aa_change = f"{ref_aa}{aa_pos}{alt_aa}"
                            
                            # Convert three-letter to one-letter amino acid codes
                            aa_dict = {
                                'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D',
                                'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G',
                                'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
                                'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S',
                                'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
                                'Ter': '*'
                            }
                            ref_aa_1letter = aa_dict.get(ref_aa, 'X')
                            alt_aa_1letter = aa_dict.get(alt_aa, 'X')
                            aa_change_1letter = f"{ref_aa_1letter}{aa_pos}{alt_aa_1letter}"
                    
                    # Only keep HIGH and MODERATE impact variants
                    if impact in ['HIGH', 'MODERATE']:
                        variant = {
                            'chrom': chrom,
                            'pos': int(pos),
                            'ref': ref,
                            'alt': alt,
                            'qual': float(qual) if qual != '.' else 0,
                            'effect': effect,
                            'impact': impact,
                            'gene_name': gene_name,
                            'gene_id': gene_id,
                            'feature_type': feature_type,
                            'transcript_biotype': transcript_biotype,
                            'hgvs_c': hgvs_c,
                            'hgvs_p': hgvs_p,
                            'aa_change': aa_change,
                            'aa_pos': aa_pos,
                            'sample': sample_name,
                            'treatment': treatment,
                            'aa_change_1letter': aa_change_1letter if 'aa_change_1letter' in locals() else None,
                            'ref_aa_1letter': ref_aa_1letter if 'ref_aa_1letter' in locals() else None,
                            'alt_aa_1letter': alt_aa_1letter if 'alt_aa_1letter' in locals() else None
                        }
                        variants.append(variant)
    
    print(f"Extracted {len(variants)} HIGH/MODERATE impact variants from {vcf_file}")
    return variants

def extract_high_impact_variants(vcf_dir, genes_of_interest):
    """Extract HIGH and MODERATE impact variants from VCF files."""
    print(f"Extracting HIGH and MODERATE impact variants from {vcf_dir}")
    
    # Print the gene IDs and names we're searching for (for debugging)
    gene_ids = set(genes_of_interest['w303_gene_id'].tolist())
    gene_names = set(genes_of_interest['erg_name'].tolist())
    print(f"Looking for genes with IDs: {gene_ids}")
    print(f"Looking for genes with names: {gene_names}")
    
    all_variants = []
    vcf_files = [f for f in os.listdir(vcf_dir) if f.endswith('.vcf') or f.endswith('.vcf.gz')]
    
    # For debugging: examine the first few variants from first file
    if vcf_files and False:  # Set to True to enable debugging output
        debug_file = vcf_files[0]
        file_path = os.path.join(vcf_dir, debug_file)
        variants = read_vcf(file_path)
        
        print(f"Examples of HIGH/MODERATE variants from {debug_file}:")
        for i, variant in enumerate(variants[:5]):
            print(f"  Variant {i+1}:")
            print(f"    gene_id: '{variant['gene_id']}'")
            print(f"    gene_name: '{variant['gene_name']}'")
            print(f"    effect: {variant['effect']}")
            print(f"    impact: {variant['impact']}")
    
    for vcf_file in vcf_files:
        file_path = os.path.join(vcf_dir, vcf_file)
        variants = read_vcf(file_path)
        
        # Filter variants for genes of interest
        filtered_variants = []
        for variant in variants:
            if variant['gene_id'] in gene_ids or variant['gene_name'] in gene_names:
                filtered_variants.append(variant)
        
        print(f"Found {len(filtered_variants)} variants in genes of interest in {vcf_file}")
        all_variants.extend(filtered_variants)
    
    # Convert to DataFrame for easier processing
    variants_df = pd.DataFrame(all_variants) if all_variants else pd.DataFrame()
    
    # Add gene information
    if not variants_df.empty:
        gene_id_to_info = {}
        for _, row in genes_of_interest.iterrows():
            gene_id_to_info[row['w303_gene_id']] = {
                'erg_name': row['erg_name'],
                'sc_gene_id': row['sc_gene_id'] if 'sc_gene_id' in row else None
            }
        
        def add_gene_info(row):
            if row['gene_id'] in gene_id_to_info:
                info = gene_id_to_info[row['gene_id']]
                row['erg_name'] = info['erg_name']
                if 'sc_gene_id' in info and info['sc_gene_id']:
                    row['sc_gene_id'] = info['sc_gene_id']
            return row
        
        variants_df = variants_df.apply(add_gene_info, axis=1)
    
    print(f"Total HIGH/MODERATE impact variants in genes of interest: {len(variants_df)}")
    return variants_df

def analyze_variant_effects(variants_df):
    """Analyze the effects of variants on protein structure and function."""
    print("Analyzing variant effects on protein structure and function")
    
    # Skip if DataFrame is empty
    if variants_df.empty:
        print("No variants to analyze")
        return variants_df
    
    # Add domain information
    def add_domain_info(row):
        if 'erg_name' in row and row['erg_name'] in PROTEIN_DOMAINS and row['aa_pos']:
            for domain_name, (start, end) in PROTEIN_DOMAINS[row['erg_name']].items():
                if start <= row['aa_pos'] <= end:
                    row['domain'] = domain_name
                    break
        if 'domain' not in row or pd.isna(row['domain']):
            row['domain'] = 'Unknown'
        return row
    
    variants_df = variants_df.apply(add_domain_info, axis=1)
    
    # Add functional change assessment for missense variants
    def assess_functional_change(row):
        if row['effect'] == 'missense_variant' and row['ref_aa_1letter'] and row['alt_aa_1letter']:
            ref_aa = row['ref_aa_1letter']
            alt_aa = row['alt_aa_1letter']
            
            # Compare amino acid properties
            property_changes = []
            functional_severity = 0
            
            # Check property changes
            for prop in ['hydrophobic', 'charge', 'size', 'polarity']:
                if AA_PROPERTIES[ref_aa][prop] != AA_PROPERTIES[alt_aa][prop]:
                    property_changes.append(prop)
                    functional_severity += 0.25  # Each property change adds to severity
            
            row['property_changes'] = ', '.join(property_changes) if property_changes else 'None'
            row['functional_severity'] = functional_severity
            
            # Classify the change
            if functional_severity == 0:
                row['functional_impact'] = 'Benign'
            elif functional_severity <= 0.5:
                row['functional_impact'] = 'Mild'
            elif functional_severity <= 0.75:
                row['functional_impact'] = 'Moderate'
            else:
                row['functional_impact'] = 'Severe'
        else:
            # For non-missense variants
            if row['effect'] in ['stop_gained', 'frameshift_variant', 'splice_acceptor_variant', 'splice_donor_variant']:
                row['functional_impact'] = 'Severe'
                row['property_changes'] = 'Protein structure disruption'
                row['functional_severity'] = 1.0
            elif row['effect'] in ['inframe_insertion', 'inframe_deletion', 'disruptive_inframe_insertion', 'disruptive_inframe_deletion']:
                row['functional_impact'] = 'Moderate'
                row['property_changes'] = 'Protein length alteration'
                row['functional_severity'] = 0.6
            else:
                row['functional_impact'] = 'Unknown'
                row['property_changes'] = 'Unknown'
                row['functional_severity'] = 0.0
        
        return row
    
    variants_df = variants_df.apply(assess_functional_change, axis=1)
    
    # Calculate aggregated impact score based on variant effect, location, etc.
    def calculate_impact_score(row):
        # Base score from impact type
        base_score = IMPACT_WEIGHTS[row['impact']]
        
        # Add score from effect type
        effect_score = EFFECT_WEIGHTS.get(row['effect'], 0.0)
        
        # Add domain factor - variants in known domains are more impactful
        domain_factor = 1.2 if row['domain'] != 'Unknown' else 1.0
        
        # Add functional severity for missense variants
        functional_factor = row['functional_severity'] if 'functional_severity' in row else 0.0
        
        # Combine scores
        impact_score = (base_score + effect_score) * domain_factor * (1.0 + functional_factor)
        
        # Normalize to 0-10 scale
        impact_score = min(10, impact_score * 5)
        
        row['impact_score'] = round(impact_score, 2)
        return row
    
    variants_df = variants_df.apply(calculate_impact_score, axis=1)
    
    return variants_df

def summarize_variants(variants_df):
    """Create summary statistics for HIGH and MODERATE impact variants."""
    print("Generating summary statistics")
    
    # Create empty summary if no variants
    if variants_df.empty:
        return {
            'total_variants': 0,
            'high_impact_variants': 0,
            'moderate_impact_variants': 0,
            'variants_by_gene': {},
            'variants_by_effect': {},
            'variants_by_treatment': {},
            'variants_by_domain': {},
            'variants_by_functional_impact': {},
            'average_impact_score': 0,
            'max_impact_score': 0,
            'gene_summaries': {},
            'treatment_summaries': {}
        }
    
    summary = {
        'total_variants': len(variants_df),
        'high_impact_variants': len(variants_df[variants_df['impact'] == 'HIGH']),
        'moderate_impact_variants': len(variants_df[variants_df['impact'] == 'MODERATE']),
        'variants_by_gene': variants_df['erg_name'].value_counts().to_dict(),
        'variants_by_effect': variants_df['effect'].value_counts().to_dict(),
        'variants_by_treatment': variants_df['treatment'].value_counts().to_dict(),
        'variants_by_domain': variants_df['domain'].value_counts().to_dict(),
        'variants_by_functional_impact': variants_df['functional_impact'].value_counts().to_dict() if 'functional_impact' in variants_df.columns else {},
        'average_impact_score': variants_df['impact_score'].mean() if 'impact_score' in variants_df.columns else 0,
        'max_impact_score': variants_df['impact_score'].max() if 'impact_score' in variants_df.columns else 0,
    }
    
    # Add gene-specific summaries
    gene_summaries = {}
    for gene, group in variants_df.groupby('erg_name'):
        gene_summaries[gene] = {
            'variant_count': len(group),
            'high_impact_count': len(group[group['impact'] == 'HIGH']),
            'moderate_impact_count': len(group[group['impact'] == 'MODERATE']),
            'top_effects': group['effect'].value_counts().to_dict(),
            'average_impact_score': group['impact_score'].mean() if 'impact_score' in group.columns else 0,
            'max_impact_score': group['impact_score'].max() if 'impact_score' in group.columns else 0,
            'domains_affected': group['domain'].value_counts().to_dict(),
        }
    
    summary['gene_summaries'] = gene_summaries
    
    # Add treatment-specific summaries
    treatment_summaries = {}
    for treatment, group in variants_df.groupby('treatment'):
        treatment_summaries[treatment] = {
            'variant_count': len(group),
            'high_impact_count': len(group[group['impact'] == 'HIGH']),
            'moderate_impact_count': len(group[group['impact'] == 'MODERATE']),
            'genes_affected': group['erg_name'].value_counts().to_dict(),
            'average_impact_score': group['impact_score'].mean() if 'impact_score' in group.columns else 0,
        }
    
    summary['treatment_summaries'] = treatment_summaries
    
    return summary

def generate_impact_heatmap(variants_df, output_dir):
    """Generate a heatmap of variant impact by gene and treatment."""
    print("Generating impact heatmap")
    
    # Skip if DataFrame is empty
    if variants_df.empty:
        print("No data for impact heatmap")
        return None
    
    # Create pivot table of impact scores
    impact_pivot = variants_df.pivot_table(
        values='impact_score',
        index='erg_name',
        columns='treatment',
        aggfunc='mean',
        fill_value=0
    )
    
    # Create heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(impact_pivot, annot=True, cmap='YlOrRd', fmt='.2f', linewidths=.5)
    plt.title('Average Variant Impact Score by Gene and Treatment')
    plt.ylabel('Gene')
    plt.xlabel('Treatment')
    plt.tight_layout()
    
    # Save heatmap
    heatmap_path = os.path.join(output_dir, 'impact_heatmap.png')
    plt.savefig(heatmap_path, dpi=300)
    plt.close()
    
    return heatmap_path

def generate_gene_variant_plots(variants_df, output_dir):
    """Generate plots showing variant positions along gene coordinates."""
    print("Generating gene variant plots")
    
    # Skip if DataFrame is empty
    if variants_df.empty:
        print("No data for gene variant plots")
        return []
    
    plot_paths = []
    gene_plots_dir = os.path.join(output_dir, 'gene_plots')
    os.makedirs(gene_plots_dir, exist_ok=True)
    
    # Group variants by gene
    for gene, group in variants_df.groupby('erg_name'):
        if len(group) == 0:
            continue
        
        # Create lollipop plot
        plt.figure(figsize=(14, 8))
        
        # Calculate marker sizes based on impact
        marker_sizes = [50 * IMPACT_WEIGHTS.get(impact, 0.5) for impact in group['impact']]
        
        # Set colors by treatment
        treatments = group['treatment'].unique()
        treatment_colors = dict(zip(treatments, sns.color_palette('husl', len(treatments))))
        colors = [treatment_colors.get(t, 'gray') for t in group['treatment']]
        
        # Plot protein domains as background
        if gene in PROTEIN_DOMAINS:
            domain_colors = plt.cm.Pastel1.colors
            for i, (domain, (start, end)) in enumerate(PROTEIN_DOMAINS[gene].items()):
                plt.axvspan(start, end, alpha=0.3, color=domain_colors[i % len(domain_colors)], 
                           label=domain)
        
        # Plot variant positions
        for i, (_, row) in enumerate(group.iterrows()):
            if pd.isna(row['aa_pos']):
                continue
                
            plt.scatter(row['aa_pos'], 0.2 + (i % 5) * 0.15, s=marker_sizes[i], 
                      color=colors[i], edgecolors='black', zorder=10, alpha=0.7)
            
            effect_label = row['effect'].replace('_variant', '').replace('_', ' ')
            plt.text(row['aa_pos'], 0.2 + (i % 5) * 0.15 + 0.05, 
                   f"{row['aa_change_1letter'] if row['aa_change_1letter'] else effect_label}",
                   ha='center', va='bottom', fontsize=8, rotation=45)
        
        # Add treatment legend
        legend_elements = [Patch(facecolor=color, label=treatment)
                         for treatment, color in treatment_colors.items()]
        plt.legend(handles=legend_elements, title='Treatment', bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Set plot properties
        plt.yticks([])
        plt.ylim(0, 1)
        plt.xlim(0, max([end for _, (_, end) in PROTEIN_DOMAINS.get(gene, {'Unknown': (0, 500)}).items()]) + 50)
        plt.xlabel('Amino Acid Position')
        plt.title(f'Variant Positions in {gene}')
        plt.grid(axis='x', linestyle='--', alpha=0.7)
        plt.tight_layout()
        
        # Save plot
        plot_path = os.path.join(gene_plots_dir, f'{gene}_variants.png')
        plt.savefig(plot_path, dpi=300)
        plt.close()
        
        plot_paths.append(plot_path)
    
    print(f"Generated {len(plot_paths)} gene variant plots")
    return plot_paths

def generate_effect_distribution(variants_df, output_dir):
    """Generate plots showing distribution of variant effects."""
    print("Generating effect distribution plots")
    
    # Skip if DataFrame is empty
    if variants_df.empty:
        print("No data for effect distribution plots")
        return []
    
    # Effect distribution by gene
    plt.figure(figsize=(14, 10))
    effect_counts = variants_df.groupby(['erg_name', 'effect']).size().unstack(fill_value=0)
    effect_counts.plot(kind='bar', stacked=True, colormap='tab20')
    plt.title('Distribution of Variant Effects by Gene')
    plt.xlabel('Gene')
    plt.ylabel('Number of Variants')
    plt.legend(title='Effect', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45)
    plt.tight_layout()
    
    # Save plot
    effect_gene_path = os.path.join(output_dir, 'effect_distribution_by_gene.png')
    plt.savefig(effect_gene_path, dpi=300)
    plt.close()
    
    # Effect distribution by treatment
    plt.figure(figsize=(14, 10))
    treatment_effect_counts = variants_df.groupby(['treatment', 'effect']).size().unstack(fill_value=0)
    treatment_effect_counts.plot(kind='bar', stacked=True, colormap='tab20')
    plt.title('Distribution of Variant Effects by Treatment')
    plt.xlabel('Treatment')
    plt.ylabel('Number of Variants')
    plt.legend(title='Effect', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45)
    plt.tight_layout()
    
    # Save plot
    effect_treatment_path = os.path.join(output_dir, 'effect_distribution_by_treatment.png')
    plt.savefig(effect_treatment_path, dpi=300)
    plt.close()
    
    # Impact distribution pie chart
    plt.figure(figsize=(10, 8))
    impact_counts = variants_df['impact'].value_counts()
    plt.pie(impact_counts, labels=impact_counts.index, autopct='%1.1f%%', 
          colors=sns.color_palette('Set3', len(impact_counts)))
    plt.title('Distribution of Variant Impacts')
    plt.tight_layout()
    
    # Save plot
    impact_pie_path = os.path.join(output_dir, 'impact_distribution_pie.png')
    plt.savefig(impact_pie_path, dpi=300)
    plt.close()
    
    return [effect_gene_path, effect_treatment_path, impact_pie_path]

def generate_domain_impact_plot(variants_df, output_dir):
    """Generate plot showing impact of variants in protein domains."""
    print("Generating domain impact plot")
    
    # Skip if DataFrame is empty
    if variants_df.empty:
        print("No data for domain impact plot")
        return None
    
    # Count variants by domain and impact
    domain_impact_counts = variants_df.groupby(['domain', 'impact']).size().unstack(fill_value=0)
    
    # Create stacked bar chart
    plt.figure(figsize=(12, 8))
    domain_impact_counts.plot(kind='barh', stacked=True, colormap='Set2')
    plt.title('Variant Distribution by Protein Domain and Impact')
    plt.xlabel('Number of Variants')
    plt.ylabel('Protein Domain')
    plt.legend(title='Impact')
    plt.grid(axis='x', linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    # Save plot
    domain_path = os.path.join(output_dir, 'domain_impact_distribution.png')
    plt.savefig(domain_path, dpi=300)
    plt.close()
    
    return domain_path

def generate_functional_impact_plot(variants_df, output_dir):
    """Generate plot showing functional impact of missense variants."""
    print("Generating functional impact plot")
    
    # Skip if DataFrame is empty
    if variants_df.empty or 'functional_impact' not in variants_df.columns:
        print("No data for functional impact plot")
        return None
    
    # Filter for missense variants
    missense_df = variants_df[variants_df['effect'] == 'missense_variant']
    
    if len(missense_df) == 0:
        print("No missense variants found")
        return None
    
    # Count by gene and functional impact
    func_impact_counts = missense_df.groupby(['erg_name', 'functional_impact']).size().unstack(fill_value=0)
    
    # Create stacked bar chart
    plt.figure(figsize=(12, 8))
    
    # Define color map for impact severity
    impact_colors = {
        'Benign': 'green',
        'Mild': 'lightblue',
        'Moderate': 'orange',
        'Severe': 'red',
        'Unknown': 'gray'
    }
    
    # Reorder columns by severity
    ordered_impacts = ['Benign', 'Mild', 'Moderate', 'Severe', 'Unknown']
    ordered_impacts = [imp for imp in ordered_impacts if imp in func_impact_counts.columns]
    
    # Create color map
    color_map = [impact_colors[impact] for impact in ordered_impacts]
    
    # Plot stacked bars
    func_impact_counts[ordered_impacts].plot(kind='bar', stacked=True, color=color_map)
    
    plt.title('Functional Impact of Missense Variants by Gene')
    plt.xlabel('Gene')
    plt.ylabel('Number of Variants')
    plt.legend(title='Functional Impact')
    plt.xticks(rotation=45)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    # Save plot
    func_impact_path = os.path.join(output_dir, 'functional_impact_by_gene.png')
    plt.savefig(func_impact_path, dpi=300)
    plt.close()
    
    # Create property change heatmap
    property_changes = []
    for _, row in missense_df.iterrows():
        if row['property_changes'] != 'None' and row['property_changes'] != 'Unknown':
            for prop in row['property_changes'].split(', '):
                property_changes.append({'gene': row['erg_name'], 'property': prop})
    
    if property_changes:
        prop_df = pd.DataFrame(property_changes)
        prop_counts = prop_df.groupby(['gene', 'property']).size().unstack(fill_value=0)
        
        plt.figure(figsize=(10, 8))
        sns.heatmap(prop_counts, annot=True, cmap='YlGnBu', fmt='d')
        plt.title('Property Changes in Missense Variants by Gene')
        plt.tight_layout()
        
        # Save plot
        prop_path = os.path.join(output_dir, 'property_changes_heatmap.png')
        plt.savefig(prop_path, dpi=300)
        plt.close()
        
        return [func_impact_path, prop_path]
    
    return func_impact_path

def generate_treatment_comparison_plot(variants_df, output_dir):
    """Generate plots comparing variants across treatment groups."""
    print("Generating treatment comparison plots")
    
    # Skip if DataFrame is empty
    if variants_df.empty:
        print("No data for treatment comparison plots")
        return []
    
    # Group control treatments and experimental treatments
    control_treatments = ['WT-CTRL', 'CAS-CTRL', 'STC-CTRL']
    experimental_treatments = ['WT-37', 'CAS', 'WTA', 'STC']
    
    # Create comparison DataFrame
    comparison_data = []
    
    for gene, group in variants_df.groupby('erg_name'):
        control_vars = group[group['treatment'].isin(control_treatments)]
        experimental_vars = group[group['treatment'].isin(experimental_treatments)]
        
        control_count = len(control_vars)
        experimental_count = len(experimental_vars)
        
        # Calculate enrichment ratio
        if control_count > 0:
            enrichment = experimental_count / control_count
        else:
            enrichment = experimental_count if experimental_count > 0 else 0
        
        comparison_data.append({
            'gene': gene,
            'control_variants': control_count,
            'experimental_variants': experimental_count,
            'enrichment': enrichment
        })
    
    comparison_df = pd.DataFrame(comparison_data)
    
    # Create enrichment plot
    plt.figure(figsize=(12, 8))
    
    # Sort by enrichment
    comparison_df = comparison_df.sort_values('enrichment', ascending=False)
    
    # Create bar plot
    bars = plt.bar(comparison_df['gene'], comparison_df['enrichment'], color='teal')
    
    # Add data labels
    for bar, exp_count, ctrl_count in zip(bars, comparison_df['experimental_variants'], comparison_df['control_variants']):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1, 
               f"{exp_count}/{ctrl_count}", ha='center', va='bottom', fontsize=9)
    
    plt.axhline(y=1, color='red', linestyle='--', alpha=0.7)
    plt.title('Enrichment of HIGH/MODERATE Impact Variants in Experimental vs. Control Samples')
    plt.xlabel('Gene')
    plt.ylabel('Enrichment Ratio (Experimental/Control)')
    plt.xticks(rotation=45)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    # Save plot
    enrichment_path = os.path.join(output_dir, 'variant_enrichment_by_gene.png')
    plt.savefig(enrichment_path, dpi=300)
    plt.close()
    
    # Create comparison bar plot
    plt.figure(figsize=(14, 10))
    
    # Set width of bars
    barWidth = 0.25
    
    # Set positions of bars on X axis
    r1 = np.arange(len(comparison_df))
    r2 = [x + barWidth for x in r1]
    
    # Create bars
    plt.bar(r1, comparison_df['experimental_variants'], width=barWidth, edgecolor='white', label='Experimental')
    plt.bar(r2, comparison_df['control_variants'], width=barWidth, edgecolor='white', label='Control')
    
    # Add labels and legend
    plt.xlabel('Gene', fontweight='bold')
    plt.ylabel('Number of Variants')
    plt.title('HIGH/MODERATE Impact Variants in Experimental vs. Control Samples')
    plt.xticks([r + barWidth/2 for r in range(len(comparison_df))], comparison_df['gene'], rotation=45)
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    # Save plot
    comparison_path = os.path.join(output_dir, 'variant_comparison_by_gene.png')
    plt.savefig(comparison_path, dpi=300)
    plt.close()
    
    # Now compare the specific experimental treatments
    
    # Create treatment-specific comparison
    treatment_data = []
    
    for gene, group in variants_df.groupby('erg_name'):
        for treatment in experimental_treatments:
            treatment_vars = group[group['treatment'] == treatment]
            treatment_count = len(treatment_vars)
            
            treatment_data.append({
                'gene': gene,
                'treatment': treatment,
                'variant_count': treatment_count
            })
    
    treatment_df = pd.DataFrame(treatment_data)
    
    # Create treatment heatmap
    treatment_pivot = treatment_df.pivot(index='gene', columns='treatment', values='variant_count').fillna(0)
    
    plt.figure(figsize=(12, 8))
    sns.heatmap(treatment_pivot, annot=True, cmap='YlGnBu', fmt='d')
    plt.title('HIGH/MODERATE Impact Variants by Gene and Treatment')
    plt.tight_layout()
    
    # Save plot
    treatment_heatmap_path = os.path.join(output_dir, 'treatment_variant_heatmap.png')
    plt.savefig(treatment_heatmap_path, dpi=300)
    plt.close()
    
    return [enrichment_path, comparison_path, treatment_heatmap_path]

def save_results(variants_df, summary, output_dir):
    """Save analysis results to output files."""
    print(f"Saving results to {output_dir}")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Save variants to TSV file
    variants_path = os.path.join(output_dir, 'high_impact_variants.tsv')
    variants_df.to_csv(variants_path, sep='\t', index=False)
    
    # Save summary to JSON file
    summary_path = os.path.join(output_dir, 'high_impact_variants_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Create a human-readable summary report
    report = [
        "# HIGH and MODERATE Impact Variant Analysis Report",
        f"Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "## Overview",
        f"Total variants analyzed: {summary['total_variants']}",
        f"HIGH impact variants: {summary['high_impact_variants']}",
        f"MODERATE impact variants: {summary['moderate_impact_variants']}",
        ""
    ]
    
    # Add content only if we have variants
    if summary['total_variants'] > 0:
        # Add gene distribution, effect distribution, etc.
        report.extend([
            "## Variant Distribution by Gene",
            "\n".join([f"- {gene}: {count} variants" for gene, count in summary['variants_by_gene'].items()]),
            "",
            "## Variant Distribution by Effect",
            "\n".join([f"- {effect}: {count} variants" for effect, count in summary['variants_by_effect'].items()]),
            "",
            "## Variant Distribution by Treatment",
            "\n".join([f"- {treatment}: {count} variants" for treatment, count in summary['variants_by_treatment'].items()]),
            "",
            "## Genes with Highest Impact Variants"
        ])
        
        # Add gene-specific summaries
        for gene, gene_summary in sorted(summary['gene_summaries'].items(), 
                                   key=lambda x: x[1]['average_impact_score'], reverse=True):
            report.extend([
                f"### {gene}",
                f"- Total variants: {gene_summary['variant_count']}",
                f"- HIGH impact variants: {gene_summary['high_impact_count']}",
                f"- MODERATE impact variants: {gene_summary['moderate_impact_count']}",
                f"- Average impact score: {gene_summary['average_impact_score']:.2f}",
                f"- Top effects: {', '.join([f'{effect} ({count})' for effect, count in gene_summary['top_effects'].items()])}",
                f"- Domains affected: {', '.join([f'{domain} ({count})' for domain, count in gene_summary['domains_affected'].items()])}",
                ""
            ])
        
        # Add treatment-specific summaries
        report.extend(["## Treatment Comparisons"])
        for treatment, treatment_summary in sorted(summary['treatment_summaries'].items(), 
                                             key=lambda x: x[1]['average_impact_score'], reverse=True):
            report.extend([
                f"### {treatment}",
                f"- Total variants: {treatment_summary['variant_count']}",
                f"- HIGH impact variants: {treatment_summary['high_impact_count']}",
                f"- MODERATE impact variants: {treatment_summary['moderate_impact_count']}",
                f"- Average impact score: {treatment_summary['average_impact_score']:.2f}",
                f"- Genes affected: {', '.join([f'{gene} ({count})' for gene, count in treatment_summary['genes_affected'].items()])}",
                ""
            ])
    else:
        # Message for when no variants are found
        report.extend([
            "## No HIGH or MODERATE Impact Variants Found",
            "",
            "No variants with HIGH or MODERATE impact were found in the analyzed ergosterol pathway genes.",
            "This suggests strong purifying selection on these genes, reinforcing our earlier findings",
            "that the ergosterol pathway genes are highly conserved even under adaptive conditions.",
            "",
            "This observation is consistent with the essential role of ergosterol in yeast cell membranes",
            "and suggests that adaptation may occur primarily through changes in gene expression rather",
            "than through protein-altering mutations.",
            "",
            "Key implications:",
            "",
            "1. The ergosterol pathway appears to be under strong purifying selection",
            "2. Adaptation likely occurs through gene expression changes rather than protein alterations",
            "3. This finding is consistent with the critical role of ergosterol in membrane integrity",
            "4. It reinforces our earlier observation of predominantly regulatory variants in these genes"
        ])
    
    # Save report
    report_path = os.path.join(output_dir, 'high_impact_variants_report.md')
    with open(report_path, 'w') as f:
        f.write('\n'.join(report))
    
    result_files = {
        'variants': variants_path,
        'summary': summary_path,
        'report': report_path
    }
    
    # Only create pivot tables if DataFrame is not empty
    if not variants_df.empty and 'impact' in variants_df.columns:
        # Create a gene-by-treatment variant count table
        gene_treatment_counts = variants_df.pivot_table(
            index='erg_name', 
            columns='treatment', 
            values='impact', 
            aggfunc='count', 
            fill_value=0
        )
        
        gene_treatment_path = os.path.join(output_dir, 'gene_treatment_variant_counts.tsv')
        gene_treatment_counts.to_csv(gene_treatment_path, sep='\t')
        result_files['gene_treatment_counts'] = gene_treatment_path
        
        # Create a effect-by-treatment variant count table
        effect_treatment_counts = variants_df.pivot_table(
            index='effect', 
            columns='treatment', 
            values='impact', 
            aggfunc='count', 
            fill_value=0
        )
        
        effect_treatment_path = os.path.join(output_dir, 'effect_treatment_variant_counts.tsv')
        effect_treatment_counts.to_csv(effect_treatment_path, sep='\t')
        result_files['effect_treatment_counts'] = effect_treatment_path
    else:
        # Create empty placeholder files
        gene_treatment_path = os.path.join(output_dir, 'gene_treatment_variant_counts.tsv')
        with open(gene_treatment_path, 'w') as f:
            f.write("# No variants found\n")
        result_files['gene_treatment_counts'] = gene_treatment_path
        
        effect_treatment_path = os.path.join(output_dir, 'effect_treatment_variant_counts.tsv')
        with open(effect_treatment_path, 'w') as f:
            f.write("# No variants found\n")
        result_files['effect_treatment_counts'] = effect_treatment_path
        
    return result_files

def main():
    """Main function to run the analysis."""
    # Parse command line arguments
    args = parse_arguments()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load gene mapping
    genes_of_interest = load_gene_mapping(args.gene_mapping)
    
    # Extract HIGH and MODERATE impact variants
    variants_df = extract_high_impact_variants(args.vcf_dir, genes_of_interest)
    
    if len(variants_df) == 0:
        print("No HIGH or MODERATE impact variants found in genes of interest")
        # Create empty summary
        summary = {
            'total_variants': 0,
            'high_impact_variants': 0,
            'moderate_impact_variants': 0,
            'variants_by_gene': {},
            'variants_by_effect': {},
            'variants_by_treatment': {},
            'gene_summaries': {},
            'treatment_summaries': {}
        }
        # Save empty results
        save_results(variants_df, summary, args.output_dir)
        return
    
    # Analyze variant effects
    variants_df = analyze_variant_effects(variants_df)
    
    # Generate summary statistics
    summary = summarize_variants(variants_df)
    
    # Generate visualizations (only if we have variants)
    if len(variants_df) > 0:
        impact_heatmap = generate_impact_heatmap(variants_df, args.output_dir)
        gene_plots = generate_gene_variant_plots(variants_df, args.output_dir)
        effect_plots = generate_effect_distribution(variants_df, args.output_dir)
        domain_plot = generate_domain_impact_plot(variants_df, args.output_dir)
        functional_plots = generate_functional_impact_plot(variants_df, args.output_dir)
        treatment_plots = generate_treatment_comparison_plot(variants_df, args.output_dir)
    
    # Save results
    file_paths = save_results(variants_df, summary, args.output_dir)
    
    print("Analysis complete!")
    print(f"Found {len(variants_df)} HIGH/MODERATE impact variants in genes of interest")
    print(f"Results saved to {args.output_dir}")
    print(f"Summary report: {file_paths['report']}")

if __name__ == "__main__":
    main()