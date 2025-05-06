#!/usr/bin/env python3
"""
analyze_promoter_elements.py

This script analyzes the upstream variants in the ergosterol pathway genes to:
1. Map their precise positions relative to transcription start sites
2. Analyze distribution patterns of regulatory variants
3. Compare variant locations across genes of interest

Usage:
  python analyze_promoter_elements.py --variant_file <path> --gene_mapping <path> --output_dir <path>
"""

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict, Counter
import json

# Known regulatory motifs in ergosterol pathway genes
# Source: Literature on ergosterol regulation in yeast
ERGOSTEROL_REGULATORY_MOTIFS = {
    'SRE': 'TCGTATA',           # Sterol Regulatory Element (Upc2p/Ecm22p binding)
    'AR1': 'CGGNNNTANCGGN',     # Anaerobic Response Element (Hap1p)
    'ROX1': 'YSYATTGTTCTC',     # Rox1p binding site (oxygen regulation)
    'MOT3': 'AAGGKA',           # Mot3p binding site (hypoxic repression)
    'PDR': 'TCCGYGGR',          # Pleiotropic Drug Response Element
}

# Define sterol regulatory element (SRE) positions
# These are approximate promoter regions (distance from TSS)
# Based on literature for ERG genes
SRE_POSITIONS = {
    'YGL012W': (-150, -100),   # ERG4
    'YGR060W': (-200, -150),   # ERG25
    'YGR175C': (-170, -120),   # ERG1
    'YHR007C': (-180, -130),   # ERG11
    'YHR072W': (-160, -110),   # ERG7
    'YHR190W': (-190, -140),   # ERG9
    'YLR056W': (-175, -125),   # ERG3
    'YML008C': (-165, -115),   # ERG6
    'YMR015C': (-155, -105),   # ERG5
    'YMR202W': (-155, -105),   # ERG2
    'YNL280C': (-160, -110),   # ERG24
}

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Analyze promoter elements for ergosterol pathway genes')
    parser.add_argument('--variant_file', required=True, help='Path to the variants TSV file')
    parser.add_argument('--gene_mapping', required=True, help='Path to the gene mapping TSV file')
    parser.add_argument('--output_dir', required=True, help='Directory to save results')
    parser.add_argument('--upstream_distance', type=int, default=1000, 
                        help='Distance upstream of TSS to analyze (default: 1000bp)')
    return parser.parse_args()

def load_variants(variant_file):
    """Load variant data from TSV file."""
    print(f"Loading variants from {variant_file}")
    variants = pd.read_csv(variant_file, sep='\t')
    # Ensure consistent column name capitalization
    variants.columns = [col.lower() for col in variants.columns]
    
    # Keep only upstream variants
    upstream_variants = variants[variants['effect'] == 'upstream_gene_variant']
    print(f"Loaded {len(variants)} variants, {len(upstream_variants)} are upstream variants")
    return upstream_variants

def load_gene_mapping(gene_mapping_file):
    """Load gene mapping information."""
    print(f"Loading gene mapping from {gene_mapping_file}")
    genes = pd.read_csv(gene_mapping_file, sep='\t')
    # Ensure consistent column name capitalization
    genes.columns = [col.lower() for col in genes.columns]
    
    print(f"Loaded information for {len(genes)} genes")
    print(f"Column names: {', '.join(genes.columns)}")
    return genes

def calculate_tss_distance(variant_row, gene_row):
    """
    Calculate the distance from a variant to the transcription start site.
    
    Args:
        variant_row: A pandas Series containing variant information
        gene_row: A pandas Series containing gene information
    
    Returns:
        int: Distance to TSS (negative values are upstream)
    """
    variant_pos = variant_row['position']
    
    # TSS depends on strand
    if gene_row['strand'] == '+':
        tss = gene_row['start']
        distance = variant_pos - tss  # Positive is downstream, negative is upstream
    else:  # '-' strand
        tss = gene_row['end']
        distance = tss - variant_pos  # Positive is downstream, negative is upstream
    
    return distance

def map_variants_to_tss(variants, genes):
    """Map variants to transcription start sites."""
    print("Mapping variant positions relative to TSS")
    
    # Create a dictionary mapping gene IDs to gene info (using lowercase column names)
    sc_gene_dict = {row['sc_gene_id']: row for _, row in genes.iterrows()}
    w303_gene_dict = {row['w303_gene_id']: row for _, row in genes.iterrows()}
    
    # Add TSS distance to each variant
    tss_distances = []
    for _, variant in variants.iterrows():
        gene_id = variant['sc_gene_id'] if 'sc_gene_id' in variant else None
        w303_id = variant['gene_id'] if 'gene_id' in variant else None
        
        # Try SC gene ID first, then W303 ID
        if gene_id and gene_id in sc_gene_dict:
            distance = calculate_tss_distance(variant, sc_gene_dict[gene_id])
            tss_distances.append(distance)
        elif w303_id and w303_id in w303_gene_dict:
            distance = calculate_tss_distance(variant, w303_gene_dict[w303_id])
            tss_distances.append(distance)
        else:
            tss_distances.append(None)
    
    variants_with_tss = variants.copy()
    variants_with_tss['tss_distance'] = tss_distances
    
    # Check if variant overlaps with known regulatory regions
    sre_overlap = []
    for _, variant in variants_with_tss.iterrows():
        gene_id = variant['sc_gene_id'] if 'sc_gene_id' in variant else None
        tss_dist = variant['tss_distance']
        
        if pd.isna(tss_dist) or not gene_id or gene_id not in SRE_POSITIONS:
            sre_overlap.append(False)
            continue
        
        sre_start, sre_end = SRE_POSITIONS[gene_id]
        overlap = sre_start <= tss_dist <= sre_end
        sre_overlap.append(overlap)
    
    variants_with_tss['sre_overlap'] = sre_overlap
    
    print(f"Successfully mapped {sum(pd.notnull(variants_with_tss['tss_distance']))} variants to TSS")
    print(f"Found {sum(variants_with_tss['sre_overlap'])} variants overlapping SRE regions")
    
    return variants_with_tss

def analyze_variant_distributions(variants):
    """Analyze the distribution patterns of variants."""
    print("Analyzing variant distribution patterns")
    
    # Filter variants with valid TSS distances
    valid_variants = variants[pd.notnull(variants['tss_distance'])]
    
    # Analyze distribution by gene
    gene_distributions = valid_variants.groupby('erg_name')['tss_distance'].agg(['count', 'min', 'max', 'mean', 'median'])
    
    # Analyze distribution by distance bins
    distance_bins = [-1000, -750, -500, -250, -100, -50, 0, 50, 100, 250, 500, 750, 1000]
    bin_labels = [f"{distance_bins[i]} to {distance_bins[i+1]}" for i in range(len(distance_bins)-1)]
    
    valid_variants['distance_bin'] = pd.cut(
        valid_variants['tss_distance'], 
        bins=distance_bins,
        labels=bin_labels
    )
    
    distance_distribution = valid_variants.groupby(['erg_name', 'distance_bin']).size().unstack(fill_value=0)
    
    # Analyze SRE overlap
    sre_overlap = valid_variants.groupby('erg_name')['sre_overlap'].sum()
    
    # Analyze distribution by treatment
    treatment_distribution = valid_variants.groupby(['treatment', 'erg_name']).size().unstack(fill_value=0)
    
    return {
        'gene_distributions': gene_distributions,
        'distance_distribution': distance_distribution,
        'sre_overlap': sre_overlap,
        'treatment_distribution': treatment_distribution
    }

def categorize_variant_positions(variants):
    """Categorize variants by their positions relative to TSS."""
    # Define position categories
    categories = {
        'Proximal Promoter': (-250, 0),
        'Core Promoter': (-50, 0),
        'TATA Box Region': (-150, -50),
        'Upstream Activating Sequence': (-500, -250),
        'Far Upstream': (-1000, -500),
        'Downstream': (0, 1000)
    }
    
    # Categorize each variant
    for category, (start, end) in categories.items():
        variants[f'in_{category.replace(" ", "_").lower()}'] = (
            variants['tss_distance'].between(start, end)
        )
    
    # Create a summary of categories
    category_counts = {
        category: sum(variants[f'in_{category.replace(" ", "_").lower()}'])
        for category in categories.keys()
    }
    
    return variants, category_counts

def generate_tss_distribution_plot(variants, output_dir):
    """Generate distribution plot of variant distances from TSS."""
    print("Generating TSS distance distribution plot")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Filter variants with valid TSS distances
    valid_variants = variants[pd.notnull(variants['tss_distance'])]
    
    # Create a filtered version with distances between -1000 and 200
    promoter_variants = valid_variants[
        (valid_variants['tss_distance'] >= -1000) & 
        (valid_variants['tss_distance'] <= 200)
    ]
    
    # Create plot
    plt.figure(figsize=(12, 6))
    
    # Plot histogram of distances
    sns.histplot(data=promoter_variants, x='tss_distance', hue='erg_name', 
                 bins=50, kde=True, element='step')
    
    # Add vertical lines for typical promoter elements
    plt.axvspan(-150, -50, alpha=0.2, color='lightgreen', label='TATA Box Region')
    plt.axvspan(-50, 0, alpha=0.2, color='lightblue', label='Core Promoter')
    
    plt.title('Distribution of Variant Distances from Transcription Start Sites')
    plt.xlabel('Distance from TSS (bp) - Negative values are upstream')
    plt.ylabel('Count')
    plt.legend(title='Gene', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    
    # Save plot
    plot_path = os.path.join(output_dir, 'tss_distance_distribution.png')
    plt.savefig(plot_path, dpi=300)
    plt.close()
    
    print(f"Saved TSS distribution plot to {plot_path}")
    
    # Create a second plot showing the entire distribution
    plt.figure(figsize=(12, 6))
    
    # Create bins with higher resolution near TSS
    custom_bins = list(range(-1000, -500, 100)) + list(range(-500, -100, 50)) + list(range(-100, 101, 10)) + list(range(150, 1001, 100))
    
    # Plot histogram of distances
    sns.histplot(data=valid_variants, x='tss_distance', bins=custom_bins, kde=False)
    
    plt.title('Full Distribution of Variant Distances from TSS')
    plt.xlabel('Distance from TSS (bp) - Negative values are upstream')
    plt.ylabel('Count')
    plt.axvline(x=0, color='red', linestyle='--', label='TSS')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    
    # Save plot
    full_plot_path = os.path.join(output_dir, 'full_tss_distance_distribution.png')
    plt.savefig(full_plot_path, dpi=300)
    plt.close()
    
    return plot_path, full_plot_path

def generate_gene_specific_plots(variants, output_dir):
    """Generate gene-specific variant position plots."""
    print("Generating gene-specific variant plots")
    
    # Create output directory if it doesn't exist
    gene_plots_dir = os.path.join(output_dir, 'gene_plots')
    os.makedirs(gene_plots_dir, exist_ok=True)
    
    # Filter variants with valid TSS distances
    valid_variants = variants[pd.notnull(variants['tss_distance'])]
    
    # Get list of genes
    genes = valid_variants['erg_name'].unique()
    
    plot_paths = []
    for gene in genes:
        # Filter variants for this gene
        gene_variants = valid_variants[valid_variants['erg_name'] == gene]
        
        # Skip if no variants
        if len(gene_variants) == 0:
            continue
            
        # Get SC gene ID for this gene
        sc_gene_id = gene_variants['sc_gene_id'].iloc[0] if 'sc_gene_id' in gene_variants.columns else None
        
        # Create plot
        plt.figure(figsize=(10, 6))
        
        # Plot histogram of distances
        sns.histplot(data=gene_variants, x='tss_distance', bins=30, kde=True)
        
        # Add SRE region if available
        if sc_gene_id and sc_gene_id in SRE_POSITIONS:
            sre_start, sre_end = SRE_POSITIONS[sc_gene_id]
            plt.axvspan(sre_start, sre_end, alpha=0.3, color='orange', 
                       label=f'SRE ({sre_start} to {sre_end}bp)')
        
        # Add standard promoter regions
        plt.axvspan(-150, -50, alpha=0.2, color='lightgreen', label='TATA Box Region')
        plt.axvspan(-50, 0, alpha=0.2, color='lightblue', label='Core Promoter')
        
        plt.title(f'Variant Positions Relative to TSS - {gene} ({sc_gene_id or "Unknown SC ID"})')
        plt.xlabel('Distance from TSS (bp) - Negative values are upstream')
        plt.ylabel('Count')
        plt.legend()
        plt.xlim(-1000, 200)  # Focus on promoter region
        plt.grid(alpha=0.3)
        plt.tight_layout()
        
        # Save plot
        plot_path = os.path.join(gene_plots_dir, f'{gene}_tss_distances.png')
        plt.savefig(plot_path, dpi=300)
        plt.close()
        
        plot_paths.append(plot_path)
    
    print(f"Generated {len(plot_paths)} gene-specific plots")
    return plot_paths

def generate_heatmap(variants, output_dir):
    """Generate heatmap of variant positions by gene and treatment."""
    print("Generating variant position heatmap")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Filter variants with valid TSS distances
    valid_variants = variants[pd.notnull(variants['tss_distance'])]
    
    # Define distance bins
    bins = [-1000, -500, -250, -100, -50, 0, 50, 100, 250, 500]
    bin_labels = [f"{bins[i]} to {bins[i+1]}" for i in range(len(bins)-1)]
    
    # Create bin column
    valid_variants['position_bin'] = pd.cut(
        valid_variants['tss_distance'], 
        bins=bins,
        labels=bin_labels
    )
    
    # Count variants by gene, bin and treatment
    gene_bin_counts = valid_variants.groupby(['erg_name', 'position_bin']).size().unstack(fill_value=0)
    treatment_bin_counts = valid_variants.groupby(['treatment', 'position_bin']).size().unstack(fill_value=0)
    
    # Create heatmaps
    plt.figure(figsize=(12, 8))
    sns.heatmap(gene_bin_counts, annot=True, cmap='YlGnBu', fmt='d')
    plt.title('Variant Distribution by Gene and Position')
    plt.ylabel('Gene')
    plt.xlabel('Position Relative to TSS')
    plt.tight_layout()
    
    # Save gene heatmap
    gene_heatmap_path = os.path.join(output_dir, 'gene_position_heatmap.png')
    plt.savefig(gene_heatmap_path, dpi=300)
    plt.close()
    
    # Create treatment heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(treatment_bin_counts, annot=True, cmap='YlOrRd', fmt='d')
    plt.title('Variant Distribution by Treatment and Position')
    plt.ylabel('Treatment')
    plt.xlabel('Position Relative to TSS')
    plt.tight_layout()
    
    # Save treatment heatmap
    treatment_heatmap_path = os.path.join(output_dir, 'treatment_position_heatmap.png')
    plt.savefig(treatment_heatmap_path, dpi=300)
    plt.close()
    
    return gene_heatmap_path, treatment_heatmap_path

def save_results(variants, distributions, category_counts, output_dir):
    """Save analysis results to output directory."""
    print(f"Saving results to {output_dir}")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Save annotated variants
    variants_path = os.path.join(output_dir, 'variants_with_tss_distance.tsv')
    variants.to_csv(variants_path, sep='\t', index=False)
    
    # Save distribution statistics
    for key, df in distributions.items():
        if isinstance(df, pd.DataFrame):
            df.to_csv(os.path.join(output_dir, f'{key}.tsv'), sep='\t')
        else:
            df.to_csv(os.path.join(output_dir, f'{key}.tsv'), sep='\t', header=True)
    
    # Save category counts
    with open(os.path.join(output_dir, 'position_category_counts.json'), 'w') as f:
        json.dump(category_counts, f, indent=2)
    
    # Create a summary report
    valid_variants = variants[pd.notnull(variants['tss_distance'])]
    
    report = {
        'total_variants': len(variants),
        'variants_with_tss_distance': len(valid_variants),
        'genes_analyzed': valid_variants['erg_name'].nunique(),
        'variants_in_sre': sum(variants['sre_overlap']) if 'sre_overlap' in variants.columns else 0,
        'position_categories': category_counts,
        'gene_counts': valid_variants['erg_name'].value_counts().to_dict(),
        'treatment_counts': valid_variants['treatment'].value_counts().to_dict(),
    }
    
    report_path = os.path.join(output_dir, 'promoter_analysis_summary.json')
    with open(report_path, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"Results saved to {output_dir}")
    return report

def main():
    """Main function to run the analysis."""
    # Parse command line arguments
    args = parse_arguments()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load data
    variants = load_variants(args.variant_file)
    genes = load_gene_mapping(args.gene_mapping)
    
    # Map variants to TSS
    variants_with_tss = map_variants_to_tss(variants, genes)
    
    # Categorize variant positions
    variants_categorized, category_counts = categorize_variant_positions(variants_with_tss)
    
    # Analyze variant distributions
    distributions = analyze_variant_distributions(variants_categorized)
    
    # Generate visualizations
    dist_plot, full_plot = generate_tss_distribution_plot(variants_categorized, args.output_dir)
    gene_plots = generate_gene_specific_plots(variants_categorized, args.output_dir)
    gene_heatmap, treatment_heatmap = generate_heatmap(variants_categorized, args.output_dir)
    
    # Save results
    report = save_results(variants_categorized, distributions, category_counts, args.output_dir)
    
    print("Analysis complete!")
    print(f"Found {report['variants_with_tss_distance']} variants with TSS distance")
    print(f"Found {report['variants_in_sre']} variants in SRE regions")
    print(f"Generated {len(gene_plots) + 4} visualization plots")

if __name__ == "__main__":
    main()