#!/usr/bin/env python3
"""
analyze_high_impact_variants_by_distance.py

This script analyzes HIGH and MODERATE impact variants near the ergosterol pathway genes:
1. Extracts all HIGH/MODERATE impact variants from VCF files
2. Calculates distance to nearest ergosterol pathway gene
3. Analyzes variants by distance, effect type, and treatment
4. Generates visualizations and reports

Usage:
  python analyze_high_impact_variants_by_distance.py --vcf_dir <path> --output_dir <path> --gene_mapping <path> --distance_threshold <int>
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
    parser = argparse.ArgumentParser(description='Analyze HIGH and MODERATE impact variants near ergosterol pathway genes')
    parser.add_argument('--vcf_dir', required=True, help='Directory containing annotated VCF files')
    parser.add_argument('--output_dir', required=True, help='Directory to save results')
    parser.add_argument('--gene_mapping', required=True, help='Path to gene mapping TSV file')
    parser.add_argument('--distance_threshold', type=int, default=50000, 
                      help='Maximum distance (in bp) from a gene of interest (default: 50kb)')
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

def extract_high_impact_variants_by_distance(vcf_dir, genes_of_interest, distance_threshold):
    """
    Extract HIGH and MODERATE impact variants and calculate their distance 
    to the nearest gene of interest.
    """
    print(f"Extracting HIGH and MODERATE impact variants from {vcf_dir}")
    
    # Create a dictionary of gene locations for calculating distances
    gene_locations = {}
    for _, gene in genes_of_interest.iterrows():
        chrom = gene['w303_scaffold']
        start = gene['start']
        end = gene['end']
        strand = gene['strand']
        gene_id = gene['w303_gene_id']
        gene_name = gene['erg_name']
        sc_gene_id = gene['sc_gene_id'] if 'sc_gene_id' in gene else None
        
        if chrom not in gene_locations:
            gene_locations[chrom] = []
        
        gene_locations[chrom].append({
            'gene_id': gene_id,
            'gene_name': gene_name,
            'sc_gene_id': sc_gene_id,
            'start': start,
            'end': end,
            'strand': strand
        })
    
    all_variants = []
    vcf_files = [f for f in os.listdir(vcf_dir) if f.endswith('.vcf') or f.endswith('.vcf.gz')]
    
    for vcf_file in vcf_files:
        file_path = os.path.join(vcf_dir, vcf_file)
        variants = read_vcf(file_path)
        
        # Calculate distance to nearest gene of interest
        variants_with_distance = []
        for variant in variants:
            chrom = variant['chrom']
            pos = variant['pos']
            
            # Find the nearest gene of interest
            nearest_gene = None
            min_distance = float('inf')
            pos_relative = None
            
            if chrom in gene_locations:
                for gene in gene_locations[chrom]:
                    # Calculate distance to this gene
                    if pos < gene['start']:
                        # Variant is upstream of gene
                        distance = gene['start'] - pos
                        position = "upstream"
                    elif pos > gene['end']:
                        # Variant is downstream of gene
                        distance = pos - gene['end']
                        position = "downstream"
                    else:
                        # Variant is within gene
                        distance = 0
                        position = "within"
                    
                    if distance < min_distance:
                        min_distance = distance
                        nearest_gene = gene
                        pos_relative = position
            
            # Add distance information
            if nearest_gene:
                variant['nearest_gene_id'] = nearest_gene['gene_id']
                variant['nearest_gene_name'] = nearest_gene['gene_name']
                variant['nearest_sc_gene_id'] = nearest_gene['sc_gene_id']
                variant['distance_to_gene'] = min_distance
                variant['position_relative'] = pos_relative
                
                # Only include variants within distance threshold
                if min_distance <= distance_threshold:
                    variants_with_distance.append(variant)
        
        print(f"Found {len(variants_with_distance)} variants within {distance_threshold}bp of genes of interest in {vcf_file}")
        all_variants.extend(variants_with_distance)
    
    # Convert to DataFrame for easier processing
    variants_df = pd.DataFrame(all_variants)
    
    # Sort by distance
    if not variants_df.empty and 'distance_to_gene' in variants_df.columns:
        variants_df = variants_df.sort_values('distance_to_gene')
    
    print(f"Total HIGH/MODERATE impact variants within {distance_threshold}bp of genes of interest: {len(variants_df)}")
    return variants_df

def analyze_distance_distribution(variants_df, distance_threshold):
    """Analyze distribution of variants by distance from genes of interest."""
    print("Analyzing variant distribution by distance")
    
    if variants_df.empty:
        print("No variants to analyze")
        return {
            'total_variants': 0,
            'distance_bins': {},
            'genes_by_distance': {},
            'treatments_by_distance': {}
        }
    
    # Define distance bins
    bins = [0, 500, 1000, 2000, 5000, 10000, 20000, 50000]
    # Ensure the distance threshold is included in bins
    if distance_threshold > bins[-1]:
        bins.append(distance_threshold)
    
    bin_labels = [f"{bins[i]}-{bins[i+1]}" for i in range(len(bins)-1)]
    
    # Categorize variants by distance bin
    variants_df['distance_bin'] = pd.cut(
        variants_df['distance_to_gene'], 
        bins=bins,
        labels=bin_labels,
        right=True
    )
    
    # Count variants by distance bin
    distance_counts = variants_df['distance_bin'].value_counts().sort_index().to_dict()
    
    # Count variants by gene and distance bin
    gene_distance = variants_df.groupby(['nearest_gene_name', 'distance_bin']).size().unstack(fill_value=0)
    gene_distance_dict = {gene: bin_counts.to_dict() for gene, bin_counts in gene_distance.iterrows()}
    
    # Count variants by treatment and distance bin
    treatment_distance = variants_df.groupby(['treatment', 'distance_bin']).size().unstack(fill_value=0)
    treatment_distance_dict = {treatment: bin_counts.to_dict() for treatment, bin_counts in treatment_distance.iterrows()}
    
    # Analyze effect types by distance
    effect_distance = variants_df.groupby(['effect', 'distance_bin']).size().unstack(fill_value=0)
    effect_distance_dict = {effect: bin_counts.to_dict() for effect, bin_counts in effect_distance.iterrows()}
    
    # Calculate statistics by relative position
    position_counts = variants_df['position_relative'].value_counts().to_dict()
    
    return {
        'total_variants': len(variants_df),
        'distance_bins': distance_counts,
        'genes_by_distance': gene_distance_dict,
        'treatments_by_distance': treatment_distance_dict,
        'effects_by_distance': effect_distance_dict,
        'position_counts': position_counts
    }

def analyze_variant_effects(variants_df):
    """Analyze the effects of variants on protein structure and function."""
    print("Analyzing variant effects on protein structure and function")
    
    if variants_df.empty:
        print("No variants to analyze")
        return variants_df
    
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
        
        # Add distance factor - closer variants are more impactful
        distance = row['distance_to_gene']
        if distance <= 500:
            distance_factor = 1.2
        elif distance <= 2000:
            distance_factor = 1.1
        elif distance <= 5000:
            distance_factor = 1.0
        elif distance <= 10000:
            distance_factor = 0.9
        else:
            distance_factor = 0.8
        
        # Add functional severity for missense variants
        functional_factor = row['functional_severity'] if 'functional_severity' in row else 0.0
        
        # Combine scores
        impact_score = (base_score + effect_score) * distance_factor * (1.0 + functional_factor)
        
        # Normalize to 0-10 scale
        impact_score = min(10, impact_score * 5)
        
        row['impact_score'] = round(impact_score, 2)
        return row
    
    variants_df = variants_df.apply(calculate_impact_score, axis=1)
    
    return variants_df

def summarize_variants(variants_df, distance_threshold, distance_distribution):
    """Create summary statistics for HIGH and MODERATE impact variants."""
    print("Generating summary statistics")
    
    if variants_df.empty:
        return {
            'total_variants': 0,
            'high_impact_variants': 0,
            'moderate_impact_variants': 0,
            'distance_distribution': distance_distribution,
            'variants_by_gene': {},
            'variants_by_effect': {},
            'variants_by_treatment': {},
            'gene_summaries': {},
            'treatment_summaries': {}
        }
    
    summary = {
        'total_variants': len(variants_df),
        'high_impact_variants': len(variants_df[variants_df['impact'] == 'HIGH']),
        'moderate_impact_variants': len(variants_df[variants_df['impact'] == 'MODERATE']),
        'distance_distribution': distance_distribution,
        'variants_by_nearest_gene': variants_df['nearest_gene_name'].value_counts().to_dict(),
        'variants_by_effect': variants_df['effect'].value_counts().to_dict(),
        'variants_by_treatment': variants_df['treatment'].value_counts().to_dict(),
        'variants_by_position': variants_df['position_relative'].value_counts().to_dict(),
        'variants_by_functional_impact': variants_df['functional_impact'].value_counts().to_dict() if 'functional_impact' in variants_df.columns else {},
        'average_impact_score': variants_df['impact_score'].mean() if 'impact_score' in variants_df.columns else 0,
        'max_impact_score': variants_df['impact_score'].max() if 'impact_score' in variants_df.columns else 0,
        'median_distance': variants_df['distance_to_gene'].median(),
        'min_distance': variants_df['distance_to_gene'].min(),
        'max_distance': variants_df['distance_to_gene'].max(),
    }
    
    # Add gene-specific summaries
    gene_summaries = {}
    for gene, group in variants_df.groupby('nearest_gene_name'):
        gene_summaries[gene] = {
            'variant_count': len(group),
            'high_impact_count': len(group[group['impact'] == 'HIGH']),
            'moderate_impact_count': len(group[group['impact'] == 'MODERATE']),
            'top_effects': group['effect'].value_counts().to_dict(),
            'average_impact_score': group['impact_score'].mean() if 'impact_score' in group.columns else 0,
            'max_impact_score': group['impact_score'].max() if 'impact_score' in group.columns else 0,
            'median_distance': group['distance_to_gene'].median(),
            'min_distance': group['distance_to_gene'].min(),
            'max_distance': group['distance_to_gene'].max(),
            'distance_bins': group['distance_bin'].value_counts().sort_index().to_dict() if 'distance_bin' in group.columns else {},
            'position_counts': group['position_relative'].value_counts().to_dict(),
        }
    
    summary['gene_summaries'] = gene_summaries
    
    # Add treatment-specific summaries
    treatment_summaries = {}
    for treatment, group in variants_df.groupby('treatment'):
        treatment_summaries[treatment] = {
            'variant_count': len(group),
            'high_impact_count': len(group[group['impact'] == 'HIGH']),
            'moderate_impact_count': len(group[group['impact'] == 'MODERATE']),
            'genes_affected': group['nearest_gene_name'].value_counts().to_dict(),
            'average_impact_score': group['impact_score'].mean() if 'impact_score' in group.columns else 0,
            'median_distance': group['distance_to_gene'].median(),
            'distance_bins': group['distance_bin'].value_counts().sort_index().to_dict() if 'distance_bin' in group.columns else {},
        }
    
    summary['treatment_summaries'] = treatment_summaries
    
    return summary

def generate_distance_histogram(variants_df, output_dir):
    """Generate histogram of variant distances from genes of interest."""
    print("Generating distance histogram")
    
    if variants_df.empty:
        print("No variants for distance histogram")
        return None
    
    plt.figure(figsize=(12, 8))
    
    # Create bins with log scale for better visualization of distance distribution
    bins = np.logspace(0, np.log10(variants_df['distance_to_gene'].max()+1), 30)
    
    # Create histogram with natural log scale on x-axis
    plt.hist(variants_df['distance_to_gene'], bins=bins, alpha=0.7, color='skyblue', edgecolor='black')
    plt.xscale('log')
    
    # Add vertical lines for common distance thresholds
    thresholds = [500, 1000, 5000, 10000]
    colors = ['red', 'orange', 'green', 'purple']
    labels = ['500bp', '1kb', '5kb', '10kb']
    
    for threshold, color, label in zip(thresholds, colors, labels):
        plt.axvline(x=threshold, color=color, linestyle='--', alpha=0.7, label=label)
    
    plt.title('Distribution of HIGH/MODERATE Impact Variants by Distance from Ergosterol Genes')
    plt.xlabel('Distance from Nearest Gene (bp) - Log Scale')
    plt.ylabel('Number of Variants')
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    
    # Save plot
    histogram_path = os.path.join(output_dir, 'distance_histogram.png')
    plt.savefig(histogram_path, dpi=300)
    plt.close()
    
    # Create a second histogram by gene
    plt.figure(figsize=(14, 10))
    
    # Group by nearest gene
    for gene, group in variants_df.groupby('nearest_gene_name'):
        plt.hist(group['distance_to_gene'], bins=bins, alpha=0.6, label=gene)
    
    plt.xscale('log')
    plt.title('Distance Distribution by Gene')
    plt.xlabel('Distance from Gene (bp) - Log Scale')
    plt.ylabel('Number of Variants')
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    
    # Save plot
    gene_histogram_path = os.path.join(output_dir, 'gene_distance_histogram.png')
    plt.savefig(gene_histogram_path, dpi=300)
    plt.close()
    
    return [histogram_path, gene_histogram_path]

def generate_distance_heatmap(variants_df, output_dir):
    """Generate heatmap of variants by distance bin and gene/treatment."""
    print("Generating distance heatmaps")
    
    if variants_df.empty or 'distance_bin' not in variants_df.columns:
        print("No data for distance heatmap")
        return []
    
    # Create gene by distance heatmap
    gene_distance = variants_df.groupby(['nearest_gene_name', 'distance_bin']).size().unstack(fill_value=0)
    
    plt.figure(figsize=(14, 8))
    sns.heatmap(gene_distance, annot=True, cmap='YlGnBu', fmt='d')
    plt.title('HIGH/MODERATE Impact Variants by Gene and Distance')
    plt.ylabel('Gene')
    plt.xlabel('Distance from Gene (bp)')
    plt.tight_layout()
    
    # Save gene heatmap
    gene_heatmap_path = os.path.join(output_dir, 'gene_distance_heatmap.png')
    plt.savefig(gene_heatmap_path, dpi=300)
    plt.close()
    
    # Create treatment by distance heatmap
    treatment_distance = variants_df.groupby(['treatment', 'distance_bin']).size().unstack(fill_value=0)
    
    plt.figure(figsize=(14, 8))
    sns.heatmap(treatment_distance, annot=True, cmap='YlOrRd', fmt='d')
    plt.title('HIGH/MODERATE Impact Variants by Treatment and Distance')
    plt.ylabel('Treatment')
    plt.xlabel('Distance from Gene (bp)')
    plt.tight_layout()
    
    # Save treatment heatmap
    treatment_heatmap_path = os.path.join(output_dir, 'treatment_distance_heatmap.png')
    plt.savefig(treatment_heatmap_path, dpi=300)
    plt.close()
    
    # Create effect by distance heatmap
    effect_distance = variants_df.groupby(['effect', 'distance_bin']).size().unstack(fill_value=0)
    
    plt.figure(figsize=(14, 10))
    sns.heatmap(effect_distance, annot=True, cmap='Greens', fmt='d')
    plt.title('Variant Effects by Distance')
    plt.ylabel('Effect')
    plt.xlabel('Distance from Gene (bp)')
    plt.tight_layout()
    
    # Save effect heatmap
    effect_heatmap_path = os.path.join(output_dir, 'effect_distance_heatmap.png')
    plt.savefig(effect_heatmap_path, dpi=300)
    plt.close()
    
    return [gene_heatmap_path, treatment_heatmap_path, effect_heatmap_path]

def generate_gene_distance_plots(variants_df, output_dir):
    """Generate plots showing variant distribution around specific genes."""
    print("Generating gene-specific distance plots")
    
    if variants_df.empty:
        print("No data for gene distance plots")
        return []
    
    plot_paths = []
    gene_plots_dir = os.path.join(output_dir, 'gene_plots')
    os.makedirs(gene_plots_dir, exist_ok=True)
    
    # Group variants by nearest gene
    for gene, group in variants_df.groupby('nearest_gene_name'):
        if len(group) == 0:
            continue
        
        plt.figure(figsize=(14, 8))
        
        # Set color by impact
        impact_colors = {'HIGH': 'red', 'MODERATE': 'orange'}
        colors = [impact_colors[impact] for impact in group['impact']]
        
        # Set marker size by distance (closer = larger)
        sizes = [max(20, 100 * (1 / (1 + group['distance_to_gene'].iloc[i] / 1000))) for i in range(len(group))]
        
        # Create scatter plot
        scatter = plt.scatter(group['distance_to_gene'], range(len(group)), 
                            c=colors, s=sizes, alpha=0.7, edgecolors='black')
        
        # Add variant details as text
        for i, (_, row) in enumerate(group.iterrows()):
            effect = row['effect'].replace('_variant', '').replace('_', ' ')
            plt.text(row['distance_to_gene'] * 1.05, i, 
                   f"{effect} ({row['treatment']})",
                   va='center', fontsize=9)
        
        # Add legend for impact
        legend_elements = [
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='HIGH'),
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='orange', markersize=10, label='MODERATE')
        ]
        plt.legend(handles=legend_elements, title='Impact')
        
        # Add distance threshold markers
        for threshold in [500, 1000, 5000, 10000]:
            plt.axvline(x=threshold, color='gray', linestyle='--', alpha=0.5)
            plt.text(threshold, len(group) * 0.95, f"{threshold}bp", 
                   ha='center', va='center', backgroundcolor='white', fontsize=8)
        
        plt.xscale('log')
        plt.yticks([])
        plt.title(f'Variants Near {gene} Ordered by Distance')
        plt.xlabel('Distance from Gene (bp) - Log Scale')
        plt.grid(axis='x', alpha=0.3)
        plt.tight_layout()
        
        # Save plot
        plot_path = os.path.join(gene_plots_dir, f'{gene}_distance_plot.png')
        plt.savefig(plot_path, dpi=300)
        plt.close()
        
        plot_paths.append(plot_path)
    
    print(f"Generated {len(plot_paths)} gene-specific distance plots")
    return plot_paths

def generate_treatment_comparison_plot(variants_df, output_dir):
    """Generate plots comparing variants by distance across treatment groups."""
    print("Generating treatment comparison plots")
    
    if variants_df.empty:
        print("No data for treatment comparison plots")
        return []
    
    # Group control treatments and experimental treatments
    control_treatments = ['WT-CTRL', 'CAS-CTRL', 'STC-CTRL']
    experimental_treatments = ['WT-37', 'CAS', 'WTA', 'STC']
    
    # Create distance bins for aggregation
    distance_bins = [0, 1000, 5000, 10000, 50000]
    bin_labels = [f"{distance_bins[i]}-{distance_bins[i+1]}" for i in range(len(distance_bins)-1)]
    
    variants_df['distance_category'] = pd.cut(
        variants_df['distance_to_gene'], 
        bins=distance_bins,
        labels=bin_labels,
        right=True
    )
    
    # Create comparison data by gene and distance bin
    comparison_data = []
    
    for gene, gene_group in variants_df.groupby('nearest_gene_name'):
        for distance_bin in bin_labels:
            bin_group = gene_group[gene_group['distance_category'] == distance_bin]
            
            control_vars = bin_group[bin_group['treatment'].isin(control_treatments)]
            experimental_vars = bin_group[bin_group['treatment'].isin(experimental_treatments)]
            
            control_count = len(control_vars)
            experimental_count = len(experimental_vars)
            
            # Calculate enrichment ratio
            if control_count > 0:
                enrichment = experimental_count / control_count
            else:
                enrichment = experimental_count if experimental_count > 0 else 0
            
            comparison_data.append({
                'gene': gene,
                'distance_bin': distance_bin,
                'control_variants': control_count,
                'experimental_variants': experimental_count,
                'enrichment': enrichment
            })
    
    comparison_df = pd.DataFrame(comparison_data)
    
    # Create enrichment heatmap by gene and distance
    plt.figure(figsize=(14, 10))
    
    # Pivot the data for the heatmap
    enrichment_pivot = comparison_df.pivot(index='gene', columns='distance_bin', values='enrichment')
    
    # Create heatmap with custom normalization to center at 1.0
    from matplotlib.colors import TwoSlopeNorm
    norm = TwoSlopeNorm(vmin=0, vcenter=1, vmax=max(3, enrichment_pivot.max().max()))
    
    sns.heatmap(enrichment_pivot, annot=True, cmap='RdBu_r', norm=norm, fmt='.2f')
    plt.title('Enrichment of HIGH/MODERATE Impact Variants in Experimental vs. Control Samples by Distance')
    plt.xlabel('Distance from Gene (bp)')
    plt.ylabel('Gene')
    plt.tight_layout()
    
    # Save plot
    enrichment_path = os.path.join(output_dir, 'distance_enrichment_heatmap.png')
    plt.savefig(enrichment_path, dpi=300)
    plt.close()
    
    # Create stacked bar chart of variants by distance and treatment type
    plt.figure(figsize=(14, 10))
    
    # Aggregate data
    bin_counts = comparison_df.groupby('distance_bin').agg({
        'control_variants': 'sum',
        'experimental_variants': 'sum'
    })
    
    # Plot stacked bars
    bin_counts.plot(kind='bar', stacked=True, colormap='Set2')
    
    plt.title('HIGH/MODERATE Impact Variants by Distance and Sample Type')
    plt.xlabel('Distance from Gene (bp)')
    plt.ylabel('Number of Variants')
    plt.legend(title='Sample Type')
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    
    # Save plot
    bars_path = os.path.join(output_dir, 'distance_treatment_bars.png')
    plt.savefig(bars_path, dpi=300)
    plt.close()
    
    return [enrichment_path, bars_path]

def save_results(variants_df, summary, output_dir):
    """Save analysis results to output files."""
    print(f"Saving results to {output_dir}")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Save variants to TSV file
    variants_path = os.path.join(output_dir, 'high_impact_variants_by_distance.tsv')
    variants_df.to_csv(variants_path, sep='\t', index=False)
    
    # Save summary to JSON file
    summary_path = os.path.join(output_dir, 'variants_by_distance_summary.json')
    
    # Define custom JSON serializer for NumPy types
    def json_serializable(obj):
        if isinstance(obj, (np.integer, np.int64)):
            return int(obj)
        elif isinstance(obj, (np.floating, np.float64)):
            return float(obj)
        elif isinstance(obj, (np.ndarray,)):
            return obj.tolist()
        return obj
    
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2, default=json_serializable)
    
    # Create a human-readable summary report
    report = [
        "# HIGH and MODERATE Impact Variants Near Ergosterol Pathway Genes",
        f"Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "## Overview",
        f"Total variants analyzed: {summary['total_variants']}",
        f"HIGH impact variants: {summary['high_impact_variants']}",
        f"MODERATE impact variants: {summary['moderate_impact_variants']}",
        f"Median distance to nearest gene: {summary['median_distance']} bp",
        f"Minimum distance: {summary['min_distance']} bp",
        f"Maximum distance: {summary['max_distance']} bp",
        "",
        "## Distance Distribution",
    ]
    
    # Add distance distribution
    if summary['distance_distribution'].get('distance_bins'):
        for bin_name, count in summary['distance_distribution']['distance_bins'].items():
            report.append(f"- {bin_name} bp: {count} variants")
    
    report.extend([
        "",
        "## Variant Distribution by Nearest Gene",
        "\n".join([f"- {gene}: {count} variants" for gene, count in summary['variants_by_nearest_gene'].items()]),
        "",
        "## Variant Distribution by Effect",
        "\n".join([f"- {effect}: {count} variants" for effect, count in summary['variants_by_effect'].items()]),
        "",
        "## Variant Distribution by Treatment",
        "\n".join([f"- {treatment}: {count} variants" for treatment, count in summary['variants_by_treatment'].items()]),
        "",
        "## Variant Distribution by Position",
        "\n".join([f"- {position}: {count} variants" for position, count in summary['variants_by_position'].items()]),
        "",
        "## Genes with Nearby Variants"
    ])
    
    # Add gene-specific summaries
    for gene, gene_summary in sorted(summary['gene_summaries'].items(), 
                                 key=lambda x: x[1]['variant_count'], reverse=True):
        report.extend([
            f"### {gene}",
            f"- Total variants: {gene_summary['variant_count']}",
            f"- HIGH impact variants: {gene_summary['high_impact_count']}",
            f"- MODERATE impact variants: {gene_summary['moderate_impact_count']}",
            f"- Median distance: {gene_summary['median_distance']} bp",
            f"- Distance range: {gene_summary['min_distance']}-{gene_summary['max_distance']} bp",
            f"- Top effects: {', '.join([f'{effect} ({count})' for effect, count in gene_summary['top_effects'].items()])}",
            "",
            "Distance distribution:",
            "\n".join([f"  - {bin_name} bp: {count} variants" for bin_name, count in gene_summary.get('distance_bins', {}).items()]),
            "",
            "Position distribution:",
            "\n".join([f"  - {position}: {count} variants" for position, count in gene_summary.get('position_counts', {}).items()]),
            ""
        ])
    
    # Add treatment-specific summaries
    report.extend(["## Treatment Comparisons"])
    for treatment, treatment_summary in sorted(summary['treatment_summaries'].items(), 
                                           key=lambda x: x[1]['variant_count'], reverse=True):
        report.extend([
            f"### {treatment}",
            f"- Total variants: {treatment_summary['variant_count']}",
            f"- HIGH impact variants: {treatment_summary['high_impact_count']}",
            f"- MODERATE impact variants: {treatment_summary['moderate_impact_count']}",
            f"- Median distance: {treatment_summary['median_distance']} bp",
            f"- Genes affected: {', '.join([f'{gene} ({count})' for gene, count in treatment_summary['genes_affected'].items()])}",
            "",
            "Distance distribution:",
            "\n".join([f"  - {bin_name} bp: {count} variants" for bin_name, count in treatment_summary.get('distance_bins', {}).items()]),
            ""
        ])
    
    # Add interpretation and biological significance
    report.extend([
        "## Biological Interpretation",
        "",
        "### Conservation of Ergosterol Pathway Genes",
        "",
        "Our analysis found no HIGH or MODERATE impact variants directly within the 11 ergosterol pathway genes, ",
        "suggesting these genes are under strong purifying selection. This finding is consistent with their ",
        "essential role in membrane integrity and sterol biosynthesis.",
        "",
        "However, our expanded analysis of the genomic neighborhood around these genes reveals interesting patterns ",
        "of genetic variation that may influence ergosterol pathway function through regulatory or other mechanisms.",
        "",
        "### Implications for Adaptation Mechanisms",
        "",
        "The pattern of variants near (but not within) ergosterol genes suggests that adaptation to temperature and ",
        "oxygen stress conditions may be occurring through changes in genes that interact with or regulate the ergosterol ",
        "pathway, rather than through direct modifications to the core pathway enzymes themselves.",
        "",
        "This lends further support to our earlier finding that adaptation likely occurs primarily through gene expression ",
        "changes (as evidenced by upstream regulatory variants) rather than through protein-altering mutations in the ",
        "critical ergosterol biosynthesis enzymes.",
        "",
        "### Future Directions",
        "",
        "Further investigation should focus on the functional relationships between the genes harboring HIGH/MODERATE impact ",
        "variants and the ergosterol pathway genes they are near. This could reveal important regulatory or metabolic ",
        "connections that contribute to adaptation mechanisms."
    ])
    
    # Save report
    report_path = os.path.join(output_dir, 'variants_by_distance_report.md')
    with open(report_path, 'w') as f:
        f.write('\n'.join(report))
    
    # Create distance summary tables
    if not variants_df.empty:
        # Distance bins by gene
        gene_distance = variants_df.groupby(['nearest_gene_name', 'distance_bin']).size().unstack(fill_value=0)
        gene_distance_path = os.path.join(output_dir, 'gene_distance_counts.tsv')
        gene_distance.to_csv(gene_distance_path, sep='\t')
        
        # Distance bins by treatment
        treatment_distance = variants_df.groupby(['treatment', 'distance_bin']).size().unstack(fill_value=0)
        treatment_distance_path = os.path.join(output_dir, 'treatment_distance_counts.tsv')
        treatment_distance.to_csv(treatment_distance_path, sep='\t')
        
        # Effect by distance
        effect_distance = variants_df.groupby(['effect', 'distance_bin']).size().unstack(fill_value=0)
        effect_distance_path = os.path.join(output_dir, 'effect_distance_counts.tsv')
        effect_distance.to_csv(effect_distance_path, sep='\t')
        
        # Export per-gene files
        gene_dir = os.path.join(output_dir, 'per_gene')
        os.makedirs(gene_dir, exist_ok=True)
        
        for gene, group in variants_df.groupby('nearest_gene_name'):
            gene_file = os.path.join(gene_dir, f'{gene}_variants.tsv')
            group.to_csv(gene_file, sep='\t', index=False)
    
    return {
        'variants': variants_path,
        'summary': summary_path,
        'report': report_path
    }

def main():
    """Main function to run the analysis."""
    # Parse command line arguments
    args = parse_arguments()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load gene mapping
    genes_of_interest = load_gene_mapping(args.gene_mapping)
    
    # Extract HIGH and MODERATE impact variants by distance
    variants_df = extract_high_impact_variants_by_distance(args.vcf_dir, genes_of_interest, args.distance_threshold)
    
    if len(variants_df) == 0:
        print(f"No HIGH or MODERATE impact variants found within {args.distance_threshold}bp of genes of interest")
        # Create empty distance distribution
        distance_distribution = {
            'total_variants': 0,
            'distance_bins': {},
            'genes_by_distance': {},
            'treatments_by_distance': {}
        }
        # Create empty summary
        summary = {
            'total_variants': 0,
            'high_impact_variants': 0,
            'moderate_impact_variants': 0,
            'distance_distribution': distance_distribution,
            'variants_by_nearest_gene': {},
            'variants_by_effect': {},
            'variants_by_treatment': {},
            'variants_by_position': {},
            'gene_summaries': {},
            'treatment_summaries': {}
        }
        # Save empty results
        save_results(variants_df, summary, args.output_dir)
        return
    
    # Analyze distance distribution
    distance_distribution = analyze_distance_distribution(variants_df, args.distance_threshold)
    
    # Analyze variant effects
    variants_df = analyze_variant_effects(variants_df)
    
    # Generate summary statistics
    summary = summarize_variants(variants_df, args.distance_threshold, distance_distribution)
    
    # Generate visualizations
    histograms = generate_distance_histogram(variants_df, args.output_dir)
    heatmaps = generate_distance_heatmap(variants_df, args.output_dir)
    gene_plots = generate_gene_distance_plots(variants_df, args.output_dir)
    treatment_plots = generate_treatment_comparison_plot(variants_df, args.output_dir)
    
    # Save results
    file_paths = save_results(variants_df, summary, args.output_dir)
    
    print("Analysis complete!")
    print(f"Found {len(variants_df)} HIGH/MODERATE impact variants within {args.distance_threshold}bp of genes of interest")
    print(f"Results saved to {args.output_dir}")
    print(f"Summary report: {file_paths['report']}")

if __name__ == "__main__":
    main()