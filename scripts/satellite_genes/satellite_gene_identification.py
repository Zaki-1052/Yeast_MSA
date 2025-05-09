#!/usr/bin/env python3
# /Users/zakiralibhai/Documents/GitHub/Yeast_MSA/scripts/satellite_genes/satellite_gene_identification.py

"""
Script to systematically identify genes in the satellite zone of each ERG gene.
This analysis builds on the four-zone conservation architecture discovered in the project,
focusing on genes in the satellite zone (50-100kb distance from ERG genes).
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import argparse
import sys
import re

# Add the parent directory to the path to access shared modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.tools import ensure_dir, save_tsv, load_tsv, setup_plotting_style, save_plot

# Ergosterol pathway genes
ERG_GENES = {
    'ERG1': 'YGR175C',
    'ERG2': 'YMR202W',
    'ERG3': 'YLR056W',
    'ERG4': 'YGL012W',
    'ERG5': 'YMR015C',
    'ERG6': 'YML008C',
    'ERG7': 'YHR072W',
    'ERG9': 'YHR190W',
    'ERG11': 'YHR007C',
    'ERG24': 'YNL280C',
    'ERG25': 'YGR060W',
}

# Conservation zone boundaries (in bp)
ZONE_BOUNDARIES = {
    'core': 0,        # Within gene
    'buffer': 5000,   # 0-5kb
    'intermediate': 50000,  # 5-50kb
    'satellite': 100000,    # 50-100kb
    'distant': float('inf')  # >100kb
}

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Identify genes in the satellite zone of ERG genes")
    parser.add_argument("--gene-mapping", default="/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/reference/gene_mapping_full.tsv",
                      help="Path to gene mapping file")
    parser.add_argument("--output-dir", default="/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/satellite_genes",
                      help="Directory to store output files")
    parser.add_argument("--min-distance", type=int, default=50000,
                      help="Minimum distance (in bp) for satellite zone")
    parser.add_argument("--max-distance", type=int, default=100000,
                      help="Maximum distance (in bp) for satellite zone")
    return parser.parse_args()

def load_gene_mapping(gene_mapping_file):
    """Load gene mapping information from file"""
    if not os.path.exists(gene_mapping_file):
        print(f"ERROR: Gene mapping file not found: {gene_mapping_file}")
        sys.exit(1)
    
    try:
        gene_df = pd.read_csv(gene_mapping_file, sep="\t")
        print(f"Loaded gene mapping with {len(gene_df)} genes")
        return gene_df
    except Exception as e:
        print(f"ERROR: Failed to load gene mapping file: {e}")
        sys.exit(1)

def extract_erg_genes(gene_df):
    """Extract ergosterol pathway gene information from gene mapping"""
    # Find ergosterol genes in sc_gene_id column
    erg_df = gene_df[gene_df['sc_gene_id'].isin(ERG_GENES.values())]
    
    # If we don't find all genes, try alternate approaches
    if len(erg_df) < len(ERG_GENES):
        # Try looking in std_gene_name column
        std_name_matches = gene_df[gene_df['std_gene_name'].isin(ERG_GENES.values())]
        if not std_name_matches.empty:
            erg_df = pd.concat([erg_df, std_name_matches])
        
        # Try looking for ERG names in other columns
        for erg_name, systematic_id in ERG_GENES.items():
            # If we didn't find this gene by systematic ID, try by ERG name
            if systematic_id not in erg_df['sc_gene_id'].values:
                # Try looking in erg_name column if it exists
                if 'erg_name' in gene_df.columns:
                    name_matches = gene_df[gene_df['erg_name'] == erg_name]
                    if not name_matches.empty:
                        erg_df = pd.concat([erg_df, name_matches])
                
                # Try looking in product column
                if 'product' in gene_df.columns:
                    product_matches = gene_df[gene_df['product'].str.contains(erg_name, case=False, na=False)]
                    if not product_matches.empty:
                        erg_df = pd.concat([erg_df, product_matches])
    
    # Reset index and drop duplicates
    erg_df = erg_df.drop_duplicates().reset_index(drop=True)
    
    return erg_df

def identify_satellite_genes(gene_df, erg_df, min_distance, max_distance):
    """Identify genes in the satellite zone around ergosterol pathway genes"""
    satellite_genes = []
    
    # For each ERG gene
    for erg_idx, erg_gene in erg_df.iterrows():
        erg_scaffold = erg_gene.get('w303_scaffold', '')
        erg_start = erg_gene.get('start', 0)
        erg_end = erg_gene.get('end', 0)
        
        # Get ERG gene name (prefer systematic ID, then standard name)
        erg_systematic_id = erg_gene.get('sc_gene_id', '')
        erg_name = None
        for name, sid in ERG_GENES.items():
            if sid == erg_systematic_id:
                erg_name = name
                break
        
        if not erg_name:
            erg_name = erg_gene.get('std_gene_name', erg_systematic_id)
        
        if not erg_scaffold or erg_start == 0 or erg_end == 0:
            print(f"WARNING: Incomplete data for ERG gene {erg_name}, skipping")
            continue
        
        # Filter genes on the same scaffold
        scaffold_genes = gene_df[gene_df['w303_scaffold'] == erg_scaffold].copy()
        
        # Calculate distance to the ERG gene for each gene
        for gene_idx, gene in scaffold_genes.iterrows():
            gene_start = gene.get('start', 0)
            gene_end = gene.get('end', 0)
            gene_id = gene.get('sc_gene_id', '')
            
            # Skip the ERG gene itself
            if gene_id == erg_systematic_id:
                continue
            
            # Calculate minimum distance between genes
            if gene_end < erg_start:
                distance = erg_start - gene_end  # Gene is upstream of ERG
                relative_position = "upstream"
            elif erg_end < gene_start:
                distance = gene_start - erg_end  # Gene is downstream of ERG
                relative_position = "downstream"
            else:
                distance = 0  # Genes overlap
                relative_position = "overlapping"
            
            # Check if this gene is in the satellite zone
            if min_distance <= distance <= max_distance:
                # Get gene information
                gene_name = gene.get('std_gene_name', gene_id)
                gene_product = gene.get('product', 'Unknown')
                
                satellite_genes.append({
                    'satellite_gene_id': gene_id,
                    'satellite_gene_name': gene_name,
                    'erg_gene_id': erg_systematic_id,
                    'erg_gene_name': erg_name,
                    'scaffold': erg_scaffold,
                    'distance': distance,
                    'distance_kb': round(distance / 1000, 2),
                    'relative_position': relative_position,
                    'satellite_start': gene_start,
                    'satellite_end': gene_end,
                    'erg_start': erg_start,
                    'erg_end': erg_end,
                    'satellite_product': gene_product,
                    'zone': 'satellite'
                })
    
    # Convert to DataFrame
    satellite_df = pd.DataFrame(satellite_genes)
    return satellite_df

def visualize_satellite_genes(satellite_df, output_dir):
    """Create visualizations of satellite gene distribution"""
    if satellite_df.empty:
        print("ERROR: No satellite genes found, skipping visualizations")
        return
    
    setup_plotting_style()
    
    # Create output directory for plots
    plot_dir = os.path.join(output_dir, "visualizations")
    ensure_dir(plot_dir)
    
    # 1. Distribution of satellite genes by ERG gene
    plt.figure(figsize=(12, 6))
    counts = satellite_df.groupby('erg_gene_name').size().sort_values(ascending=False)
    ax = sns.barplot(x=counts.index, y=counts.values)
    plt.title('Number of Satellite Genes by ERG Gene')
    plt.ylabel('Number of Satellite Genes')
    plt.xlabel('ERG Gene')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    save_plot(plt, os.path.join(plot_dir, "satellite_genes_by_erg.png"))
    plt.close()
    
    # 2. Distance distribution of satellite genes
    plt.figure(figsize=(10, 6))
    sns.histplot(data=satellite_df, x='distance_kb', bins=20, kde=True)
    plt.axvline(x=50, color='red', linestyle='--', alpha=0.7, label='Min Satellite Zone (50kb)')
    plt.axvline(x=100, color='red', linestyle='--', alpha=0.7, label='Max Satellite Zone (100kb)')
    plt.title('Distribution of Distances to ERG Genes')
    plt.xlabel('Distance (kb)')
    plt.ylabel('Number of Genes')
    plt.legend()
    plt.tight_layout()
    save_plot(plt, os.path.join(plot_dir, "satellite_gene_distances.png"))
    plt.close()
    
    # 3. Upstream vs downstream distribution
    plt.figure(figsize=(8, 6))
    pos_counts = satellite_df['relative_position'].value_counts()
    ax = sns.barplot(x=pos_counts.index, y=pos_counts.values)
    plt.title('Distribution of Satellite Genes by Position')
    plt.ylabel('Number of Genes')
    plt.xlabel('Position Relative to ERG Gene')
    plt.tight_layout()
    save_plot(plt, os.path.join(plot_dir, "satellite_gene_positions.png"))
    plt.close()
    
    # 4. Heatmap of satellite gene counts by ERG gene and position
    plt.figure(figsize=(12, 8))
    heatmap_data = pd.crosstab(satellite_df['erg_gene_name'], satellite_df['relative_position'])
    sns.heatmap(heatmap_data, annot=True, fmt='d', cmap='YlGnBu')
    plt.title('Distribution of Satellite Genes by ERG Gene and Position')
    plt.tight_layout()
    save_plot(plt, os.path.join(plot_dir, "satellite_gene_heatmap.png"))
    plt.close()

def main():
    # Parse command line arguments
    args = parse_args()
    
    # Ensure output directory exists
    ensure_dir(args.output_dir)
    
    print("======================================================")
    print("Satellite Gene Identification")
    print("======================================================")
    print(f"Minimum distance: {args.min_distance} bp")
    print(f"Maximum distance: {args.max_distance} bp")
    print("======================================================")
    
    # Load gene mapping data
    gene_df = load_gene_mapping(args.gene_mapping)
    
    # Extract ERG gene information
    print("Extracting ERG gene information...")
    erg_df = extract_erg_genes(gene_df)
    print(f"Found {len(erg_df)} ERG genes out of {len(ERG_GENES)} expected")
    
    # Identify satellite genes around ERG genes
    print("Identifying satellite genes around ERG genes...")
    satellite_df = identify_satellite_genes(gene_df, erg_df, args.min_distance, args.max_distance)
    print(f"Found {len(satellite_df)} satellite genes across all ERG genes")
    
    # Save satellite gene data to file
    satellite_file = os.path.join(args.output_dir, "satellite_genes.tsv")
    save_tsv(satellite_df, satellite_file)
    print(f"Satellite gene data saved to {satellite_file}")
    
    # Create summary statistics
    if not satellite_df.empty:
        print("\nSummary Statistics:")
        print(f"Total satellite genes: {len(satellite_df)}")
        per_erg = satellite_df.groupby('erg_gene_name').size()
        print(f"Average satellite genes per ERG gene: {per_erg.mean():.2f}")
        print(f"Median satellite genes per ERG gene: {per_erg.median():.2f}")
        
        upstream = len(satellite_df[satellite_df['relative_position'] == 'upstream'])
        downstream = len(satellite_df[satellite_df['relative_position'] == 'downstream'])
        print(f"Upstream satellite genes: {upstream} ({upstream/len(satellite_df)*100:.1f}%)")
        print(f"Downstream satellite genes: {downstream} ({downstream/len(satellite_df)*100:.1f}%)")
        
        # Create visualizations
        print("\nCreating visualizations...")
        visualize_satellite_genes(satellite_df, args.output_dir)
    
    # Write summary report
    summary_file = os.path.join(args.output_dir, "satellite_gene_identification_summary.txt")
    with open(summary_file, 'w') as f:
        f.write("Satellite Gene Identification Summary\n")
        f.write("====================================\n\n")
        
        f.write(f"Analysis Parameters:\n")
        f.write(f"- Minimum distance: {args.min_distance} bp\n")
        f.write(f"- Maximum distance: {args.max_distance} bp\n\n")
        
        f.write(f"Results Summary:\n")
        f.write(f"- Total ERG genes analyzed: {len(erg_df)}\n")
        f.write(f"- Total satellite genes identified: {len(satellite_df)}\n")
        
        if not satellite_df.empty:
            f.write(f"- Average satellite genes per ERG gene: {per_erg.mean():.2f}\n")
            f.write(f"- Median satellite genes per ERG gene: {per_erg.median():.2f}\n")
            f.write(f"- Upstream satellite genes: {upstream} ({upstream/len(satellite_df)*100:.1f}%)\n")
            f.write(f"- Downstream satellite genes: {downstream} ({downstream/len(satellite_df)*100:.1f}%)\n\n")
            
            f.write("Satellite Genes per ERG Gene:\n")
            for erg_name, count in per_erg.items():
                f.write(f"- {erg_name}: {count} genes\n")
        
        f.write("\nFor detailed results, please see the satellite_genes.tsv file.\n")
    
    print(f"Summary report saved to {summary_file}")
    print("Satellite gene identification complete!")

if __name__ == "__main__":
    main()