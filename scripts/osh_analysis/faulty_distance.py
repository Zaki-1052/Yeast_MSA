#!/usr/bin/env python3
# /Users/zakiralibhai/Documents/GitHub/Yeast_MSA/scripts/osh_analysis/osh_erg_distance.py

"""
Script to calculate genomic distances between OSH genes and ergosterol pathway genes.
This analysis will help understand potential regulatory relationships and implications
for sterol transport and adaptation.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import sys
from collections import defaultdict

# Add the parent directory to the path to access shared modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.tools import ensure_dir

# OSH gene names and standard names
OSH_GENES = {
    'OSH2': 'YDL019C',
    'OSH3': 'YHR073W',
    'OSH5': 'YOR237W',
    'OSH6': 'YKR003W',
    'OSH7': 'YHR001W'
}

# Ergosterol pathway genes (for comparison)
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

# Conservation zones definitions (based on the 4-zone architecture)
ZONES = {
    'Core': (0, 0),           # Within gene
    'Buffer': (1, 5000),      # 1-5000bp
    'Intermediate': (5001, 50000),  # 5-50kb
    'Satellite': (50001, float('inf'))  # >50kb
}

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Calculate genomic distances between OSH and ERG genes")
    parser.add_argument("--gene-mapping", default="/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/reference/gene_mapping_full.tsv",
                      help="Path to gene mapping file")
    parser.add_argument("--output-dir", default="/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/osh_analysis",
                      help="Directory to store output files")
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

def extract_genes_of_interest(gene_df, gene_list, gene_type_name):
    """Extract specific genes from the gene mapping dataframe"""
    # Print available columns in gene_df for debugging
    print(f"Available columns in gene_df: {gene_df.columns.tolist()}")
    
    # Try to find genes of interest using multiple approaches
    genes_of_interest = pd.DataFrame()
    
    # Approach 1: Try matching on sc_gene_id
    if 'sc_gene_id' in gene_df.columns:
        sc_matches = gene_df[gene_df['sc_gene_id'].isin(gene_list.values())]
        if not sc_matches.empty:
            print(f"Found {len(sc_matches)} {gene_type_name} genes by matching sc_gene_id")
            genes_of_interest = pd.concat([genes_of_interest, sc_matches])
    
    # Approach 2: Try matching on std_gene_name
    if 'std_gene_name' in gene_df.columns:
        std_matches = gene_df[gene_df['std_gene_name'].isin(gene_list.values())]
        if not std_matches.empty:
            print(f"Found {len(std_matches)} {gene_type_name} genes by matching std_gene_name")
            genes_of_interest = pd.concat([genes_of_interest, std_matches])
    
    # Approach 3: Try matching gene names against any column that might contain them
    for col in gene_df.columns:
        if col.endswith('_name') or col.endswith('_id') or col == 'description' or col == 'product':
            try:
                # Try matching against gene list keys (gene names)
                name_matches = gene_df[gene_df[col].isin(gene_list.keys())]
                if not name_matches.empty:
                    print(f"Found {len(name_matches)} {gene_type_name} genes by matching gene names in {col}")
                    genes_of_interest = pd.concat([genes_of_interest, name_matches])
            except Exception as e:
                print(f"Error checking {col}: {e}")
    
    # Approach 4: Try fuzzy matching
    for value in list(gene_list.values()) + list(gene_list.keys()):
        for col in gene_df.columns:
            if col.endswith('_name') or col.endswith('_id') or col == 'description' or col == 'product':
                try:
                    # Handle different column types
                    if gene_df[col].dtype == 'object':  # String column
                        fuzzy_matches = gene_df[gene_df[col].astype(str).str.contains(value, case=False, na=False)]
                        if not fuzzy_matches.empty:
                            print(f"Found {len(fuzzy_matches)} {gene_type_name} genes by fuzzy matching '{value}' in {col}")
                            genes_of_interest = pd.concat([genes_of_interest, fuzzy_matches])
                except Exception as e:
                    pass  # Silently continue if error
    
    # Remove duplicates
    if not genes_of_interest.empty:
        genes_of_interest = genes_of_interest.drop_duplicates().reset_index(drop=True)
    
    if genes_of_interest.empty:
        print(f"Warning: No {gene_type_name} genes found in the gene mapping.")
        return pd.DataFrame()
    
    print(f"Found {len(genes_of_interest)} {gene_type_name} genes for analysis")
    return genes_of_interest

def standardize_gene_info(gene_df):
    """Standardize gene information for distance calculation"""
    # Create a standardized dataframe with consistent column names
    std_genes = []
    
    for _, gene in gene_df.iterrows():
        # Get gene name (try multiple fields)
        gene_name = None
        for field in ['std_gene_name', 'erg_name', 'osh_name']:
            if field in gene and not pd.isna(gene[field]) and gene[field]:
                gene_name = gene[field]
                break
        if not gene_name:
            # As a last resort, use w303_gene_id or sc_gene_id
            gene_name = gene.get('w303_gene_id', gene.get('sc_gene_id', 'Unknown'))
        
        # Get other gene information
        scaffold = gene.get('w303_scaffold', '')
        start = gene.get('start', 0)
        end = gene.get('end', 0)
        
        if not scaffold or start == 0 or end == 0:
            print(f"Warning: Incomplete gene information for {gene_name}, skipping")
            continue
        
        # Determine gene type
        gene_type = 'Other'
        if gene_name in OSH_GENES.values() or any(gene.get(col) in OSH_GENES.keys() for col in gene.index if not pd.isna(gene.get(col))):
            gene_type = 'OSH'
        elif gene_name in ERG_GENES.values() or any(gene.get(col) in ERG_GENES.keys() for col in gene.index if not pd.isna(gene.get(col))):
            gene_type = 'ERG'
        
        # Add standardized gene info
        std_genes.append({
            'gene_name': gene_name,
            'scaffold': scaffold,
            'start': start,
            'end': end,
            'gene_type': gene_type,
            'w303_gene_id': gene.get('w303_gene_id', ''),
            'sc_gene_id': gene.get('sc_gene_id', ''),
            'length': end - start + 1
        })
    
    return pd.DataFrame(std_genes)

def calculate_gene_distances(osh_genes, erg_genes):
    """Calculate distances between OSH and ERG genes"""
    # Create a list to store distance information
    distances = []
    
    # Iterate over each OSH gene
    for _, osh in osh_genes.iterrows():
        osh_name = osh['gene_name']
        osh_scaffold = osh['scaffold']
        osh_start = osh['start']
        osh_end = osh['end']
        
        # Find ERG genes on the same scaffold
        same_scaffold_erg = erg_genes[erg_genes['scaffold'] == osh_scaffold].copy()
        
        # Calculate distance to each ERG gene on the same scaffold
        for _, erg in same_scaffold_erg.iterrows():
            erg_name = erg['gene_name']
            erg_start = erg['start']
            erg_end = erg['end']
            
            # Calculate distance between genes
            if osh_end < erg_start:  # OSH is upstream of ERG
                distance = erg_start - osh_end
                relative_position = "upstream"
            elif erg_end < osh_start:  # OSH is downstream of ERG
                distance = osh_start - erg_end
                relative_position = "downstream"
            else:  # Genes overlap
                distance = 0
                relative_position = "overlapping"
            
            # Determine which conservation zone this falls into
            zone = "Unknown"
            for zone_name, (min_dist, max_dist) in ZONES.items():
                if min_dist <= distance <= max_dist:
                    zone = zone_name
                    break
            
            # Add to distances list
            distances.append({
                'osh_gene': osh_name,
                'erg_gene': erg_name,
                'scaffold': osh_scaffold,
                'distance': distance,
                'relative_position': relative_position,
                'zone': zone,
                'osh_start': osh_start,
                'osh_end': osh_end,
                'erg_start': erg_start,
                'erg_end': erg_end
            })
    
    # Convert to DataFrame
    distances_df = pd.DataFrame(distances)
    
    # Calculate additional statistics if we have data
    if not distances_df.empty:
        # Sort by distance
        distances_df = distances_df.sort_values('distance')
        
        # Count number of pairs in each zone
        zone_counts = distances_df['zone'].value_counts()
        print("\nOSH-ERG gene pairs by conservation zone:")
        for zone, count in zone_counts.items():
            print(f"  {zone}: {count} pairs")
        
        # Calculate minimum distance for each OSH gene to any ERG gene
        min_distances = distances_df.groupby('osh_gene')['distance'].min().reset_index()
        min_distances = min_distances.rename(columns={'distance': 'min_distance'})
        print("\nMinimum distances from OSH genes to nearest ERG gene:")
        for _, row in min_distances.iterrows():
            print(f"  {row['osh_gene']}: {row['min_distance']} bp")
        
        # Calculate minimum distance for each ERG gene to any OSH gene
        min_distances_erg = distances_df.groupby('erg_gene')['distance'].min().reset_index()
        min_distances_erg = min_distances_erg.rename(columns={'distance': 'min_distance'})
        print("\nMinimum distances from ERG genes to nearest OSH gene:")
        for _, row in min_distances_erg.iterrows():
            print(f"  {row['erg_gene']}: {row['min_distance']} bp")
    
    return distances_df

def visualize_distances(distances_df, output_dir):
    """Create visualizations for the OSH-ERG gene distances"""
    if distances_df.empty:
        print("No distance data to visualize")
        return
    
    # Create a histogram of distances
    plt.figure(figsize=(12, 6))
    sns.histplot(data=distances_df, x='distance', bins=20, hue='zone', multiple='stack')
    plt.title('Distribution of Distances Between OSH and ERG Genes')
    plt.xlabel('Distance (bp)')
    plt.ylabel('Count')
    plt.tight_layout()
    
    # Save the figure
    output_file = os.path.join(output_dir, "osh_erg_distance_histogram.png")
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Distance histogram saved to {output_file}")
    
    # Create a heatmap of distances between OSH and ERG genes
    pivot_df = distances_df.pivot_table(
        index='osh_gene', 
        columns='erg_gene', 
        values='distance',
        aggfunc='min'  # Use minimum distance if multiple pairs
    )
    
    # Create a mask for missing values to display them as white
    mask = pivot_df.isnull()
    
    plt.figure(figsize=(14, 10))
    ax = sns.heatmap(
        pivot_df,
        cmap="YlGnBu",
        mask=mask,
        annot=True,
        fmt=".0f",
        cbar_kws={'label': 'Distance (bp)'}
    )
    
    plt.title('Distances Between OSH and ERG Genes')
    plt.tight_layout()
    
    # Save the figure
    output_file = os.path.join(output_dir, "osh_erg_distance_heatmap.png")
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Distance heatmap saved to {output_file}")
    
    # Create a barplot of minimum distances for each OSH gene
    min_distances = distances_df.groupby('osh_gene')['distance'].min().reset_index()
    
    plt.figure(figsize=(10, 6))
    sns.barplot(data=min_distances, x='osh_gene', y='distance')
    plt.title('Minimum Distance from OSH Genes to ERG Genes')
    plt.xlabel('OSH Gene')
    plt.ylabel('Minimum Distance (bp)')
    plt.xticks(rotation=45)
    plt.tight_layout()
    
    # Save the figure
    output_file = os.path.join(output_dir, "osh_min_distance_barplot.png")
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Minimum distance barplot saved to {output_file}")
    
    # Create a zone distribution plot
    plt.figure(figsize=(10, 6))
    zone_counts = distances_df['zone'].value_counts().reset_index()
    zone_counts.columns = ['Zone', 'Count']
    
    # Ensure zones are in the right order
    zone_order = ['Core', 'Buffer', 'Intermediate', 'Satellite']
    zone_counts['Zone'] = pd.Categorical(zone_counts['Zone'], categories=zone_order, ordered=True)
    zone_counts = zone_counts.sort_values('Zone')
    
    sns.barplot(data=zone_counts, x='Zone', y='Count')
    plt.title('Distribution of OSH-ERG Gene Pairs by Conservation Zone')
    plt.xlabel('Conservation Zone')
    plt.ylabel('Number of Gene Pairs')
    plt.tight_layout()
    
    # Save the figure
    output_file = os.path.join(output_dir, "osh_erg_zone_distribution.png")
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Zone distribution plot saved to {output_file}")

def generate_summary_report(distances_df, osh_genes, erg_genes, output_dir):
    """Generate a summary report of OSH-ERG gene distances"""
    # Create a summary report
    summary = [
        "OSH-ERG Gene Distance Analysis Summary",
        "=====================================",
        f"OSH genes analyzed: {len(osh_genes)}",
        f"ERG genes analyzed: {len(erg_genes)}",
        f"Total OSH-ERG gene pairs analyzed: {len(distances_df)}",
        "",
        "Zone Definitions:",
        "  Core: 0 bp (overlapping genes)",
        "  Buffer: 1-5,000 bp",
        "  Intermediate: 5,001-50,000 bp",
        "  Satellite: >50,000 bp",
        "",
        "OSH-ERG Gene Pairs by Zone:",
    ]
    
    if not distances_df.empty:
        # Count gene pairs by zone
        zone_counts = distances_df['zone'].value_counts().sort_index()
        for zone, count in zone_counts.items():
            summary.append(f"  {zone}: {count} pairs ({100 * count / len(distances_df):.1f}%)")
        
        # Calculate statistics
        summary.extend([
            "",
            "Distance Statistics:",
            f"  Minimum distance: {distances_df['distance'].min()} bp",
            f"  Maximum distance: {distances_df['distance'].max()} bp",
            f"  Median distance: {distances_df['distance'].median()} bp",
            f"  Mean distance: {distances_df['distance'].mean():.1f} bp",
            "",
            "Top 5 Closest OSH-ERG Gene Pairs:",
        ])
        
        # Add top 5 closest pairs
        closest_pairs = distances_df.sort_values('distance').head(5)
        for i, (_, pair) in enumerate(closest_pairs.iterrows(), 1):
            summary.append(f"  {i}. {pair['osh_gene']} - {pair['erg_gene']}: {pair['distance']} bp ({pair['relative_position']})")
        
        # Add OSH genes with no nearby ERG genes
        if len(osh_genes) > distances_df['osh_gene'].nunique():
            summary.append("")
            summary.append("OSH genes with no ERG genes on same scaffold:")
            for _, osh in osh_genes.iterrows():
                if osh['gene_name'] not in distances_df['osh_gene'].values:
                    summary.append(f"  {osh['gene_name']} (scaffold {osh['scaffold']})")
    else:
        summary.append("  No OSH-ERG gene pairs found on the same scaffolds")
    
    # Write summary to file
    summary_file = os.path.join(output_dir, "osh_erg_distance_summary.txt")
    with open(summary_file, 'w') as f:
        f.write('\n'.join(summary))
    
    print(f"Summary report saved to {summary_file}")

def main():
    # Parse command line arguments
    args = parse_args()
    
    # Ensure output directory exists
    ensure_dir(args.output_dir)
    
    # Load gene mapping data
    gene_df = load_gene_mapping(args.gene_mapping)
    
    # Extract OSH genes
    print("Extracting OSH gene information...")
    osh_genes_df = extract_genes_of_interest(gene_df, OSH_GENES, "OSH")
    
    # Extract ERG genes
    print("Extracting ERG gene information...")
    erg_genes_df = extract_genes_of_interest(gene_df, ERG_GENES, "ERG")
    
    # Standardize gene information
    print("Standardizing gene information...")
    osh_genes = standardize_gene_info(osh_genes_df)
    erg_genes = standardize_gene_info(erg_genes_df)
    
    # Calculate distances between OSH and ERG genes
    print("Calculating distances between OSH and ERG genes...")
    distances_df = calculate_gene_distances(osh_genes, erg_genes)
    
    # Save distance data to file
    distances_file = os.path.join(args.output_dir, "osh_erg_distances.tsv")
    if not distances_df.empty:
        distances_df.to_csv(distances_file, sep='\t', index=False)
        print(f"OSH-ERG distances saved to {distances_file}")
    else:
        print("No distances to save")
    
    # Visualize distances
    print("Generating visualizations...")
    visualize_distances(distances_df, args.output_dir)
    
    # Generate summary report
    print("Generating summary report...")
    generate_summary_report(distances_df, osh_genes, erg_genes, args.output_dir)
    
    print("OSH-ERG distance analysis complete!")

if __name__ == "__main__":
    main()