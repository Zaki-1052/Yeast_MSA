#!/usr/bin/env python3
# /Users/zakiralibhai/Documents/GitHub/Yeast_MSA/scripts/osh_analysis/osh_erg_distance.py

"""
Script to calculate genomic distances between OSH genes and ergosterol pathway genes,
exploring their spatial relationships and potential regulatory connections.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
import argparse
import sys
from collections import defaultdict
import networkx as nx
from scipy import stats

# Add the parent directory to the path to access shared modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.tools import ensure_dir

# OSH gene names and standard names
# Note: Some OSH genes might not be present in the W303 genome annotation
# We identify them by systematic ID in the std_gene_name column
OSH_GENES = {
    'OSH1': 'YAR042W',
    'OSH2': 'YDL019C',
    'OSH3': 'YHR073W',
    'OSH4': 'YPL145C',
    'OSH5': 'YOR237W',
    'OSH6': 'YKR003W',
    'OSH7': 'YHR001W'
}

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

# Define distance thresholds for proximity zones
DISTANCE_ZONES = {
    'Core': 0,
    'Buffer': 5000,
    'Intermediate': 25000,
    'Satellite': 100000
}

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Calculate distances between OSH and ERG genes")
    parser.add_argument("--gene-mapping", default="/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/reference/gene_mapping_full.tsv",
                      help="Path to gene mapping file")
    parser.add_argument("--output-dir", default="/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/osh_analysis",
                      help="Directory to store output files")
    parser.add_argument("--genome-file", default="/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/reference/w303_chromosomal.fasta",
                      help="Path to reference genome FASTA file")
    parser.add_argument("--variant-dir", default="/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/gene_variants_expanded",
                      help="Directory containing variant data")
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

def extract_osh_genes(gene_df):
    """Extract OSH gene family information from gene mapping"""
    print("Looking for OSH genes with the following IDs:", list(OSH_GENES.values()))

    # First try exact matches for standard gene names in std_gene_name column
    osh_df = gene_df[gene_df['std_gene_name'].isin(OSH_GENES.values())]
    if not osh_df.empty:
        print(f"Found {len(osh_df)} OSH genes by std_gene_name match")

    # Also try sc_gene_id column
    sc_matches = gene_df[gene_df['sc_gene_id'].isin(OSH_GENES.values())]
    if not sc_matches.empty:
        print(f"Found {len(sc_matches)} OSH genes by sc_gene_id match")
        osh_df = pd.concat([osh_df, sc_matches])

    # If we don't find all genes, try matching on gene name
    if len(osh_df) < len(OSH_GENES):
        for osh_name, systematic_id in OSH_GENES.items():
            # Try all possible ways to find the gene
            if systematic_id not in osh_df['std_gene_name'].values and systematic_id not in osh_df['sc_gene_id'].values:
                # Try direct match on various columns
                for col in ['std_gene_name', 'sc_gene_id', 'w303_gene_id', 'product']:
                    if col in gene_df.columns:
                        try:
                            # Look for either the name or the ID
                            for search_term in [osh_name, systematic_id]:
                                partial_matches = gene_df[gene_df[col].astype(str).str.contains(search_term, case=False, na=False)]
                                if not partial_matches.empty:
                                    print(f"Found {len(partial_matches)} matches for {osh_name} ({systematic_id}) in {col}")
                                    osh_df = pd.concat([osh_df, partial_matches])
                        except Exception as e:
                            print(f"Error searching for {systematic_id} in {col}: {e}")

    # If we still don't have all genes, try fuzzy search with OSH names
    if len(osh_df) < len(OSH_GENES):
        # Try to find any gene with "OSH" in the name or product
        for search_term in ["OSH", "oxysterol"]:
            for col in ['std_gene_name', 'sc_gene_id', 'w303_gene_id', 'product']:
                if col in gene_df.columns:
                    try:
                        osh_matches = gene_df[gene_df[col].astype(str).str.contains(search_term, case=False, na=False)]
                        if not osh_matches.empty:
                            print(f"Found {len(osh_matches)} matches with '{search_term}' in {col}")
                            osh_df = pd.concat([osh_df, osh_matches])
                    except Exception as e:
                        print(f"Error searching for {search_term} in {col}: {e}")

    # Reset index and drop duplicates
    osh_df = osh_df.drop_duplicates().reset_index(drop=True)

    return osh_df

def extract_erg_genes(gene_df):
    """Extract ergosterol pathway gene information from gene mapping"""
    print("Looking for ERG genes with the following IDs:", list(ERG_GENES.values()))

    # First look for genes where erg_name column matches ERG gene names
    erg_df = gene_df[gene_df['erg_name'].isin(ERG_GENES.keys())]
    if not erg_df.empty:
        print(f"Found {len(erg_df)} ERG genes by matching erg_name column")

    # Also try looking for gene IDs in the std_gene_name column
    std_matches = gene_df[gene_df['std_gene_name'].isin(ERG_GENES.values())]
    if not std_matches.empty:
        print(f"Found {len(std_matches)} ERG genes by matching std_gene_name column")
        erg_df = pd.concat([erg_df, std_matches])

    # Also try looking in sc_gene_id column
    sc_matches = gene_df[gene_df['sc_gene_id'].isin(ERG_GENES.values())]
    if not sc_matches.empty:
        print(f"Found {len(sc_matches)} ERG genes by matching sc_gene_id column")
        erg_df = pd.concat([erg_df, sc_matches])

    # Reset index and drop duplicates
    erg_df = erg_df.drop_duplicates().reset_index(drop=True)

    return erg_df

def load_variant_data(variant_dir):
    """Load variant data from variant files"""
    # Check if the directory exists
    if not os.path.exists(variant_dir):
        print(f"ERROR: Variant directory not found: {variant_dir}")
        sys.exit(1)

    # Load all gene variants file
    all_variants_file = os.path.join(variant_dir, "all_gene_variants.tsv")
    if not os.path.exists(all_variants_file):
        print(f"ERROR: All gene variants file not found: {all_variants_file}")
        sys.exit(1)

    try:
        variants_df = pd.read_csv(all_variants_file, sep="\t")
        print(f"Loaded variant data with {len(variants_df)} variants")
        return variants_df
    except Exception as e:
        print(f"ERROR: Failed to load variant data: {e}")
        sys.exit(1)

def calculate_gene_distances(osh_df, erg_df):
    """Calculate distances between OSH and ERG genes"""
    # Create a list to store results
    distance_data = []
    
    # Loop through OSH genes
    for _, osh_gene in osh_df.iterrows():
        osh_scaffold = osh_gene.get('w303_scaffold', '')
        osh_start = osh_gene.get('start', 0)
        osh_end = osh_gene.get('end', 0)
        
        # Get OSH gene name
        if 'sc_gene_id' in osh_gene and pd.notna(osh_gene['sc_gene_id']):
            # Try to find the gene name from our mapping
            osh_id = osh_gene['sc_gene_id']
            osh_name = None
            for gene_name, gene_id in OSH_GENES.items():
                if gene_id == osh_id:
                    osh_name = gene_name
                    break
            if not osh_name:
                osh_name = osh_id
        else:
            osh_name = osh_gene.get('w303_gene_id', 'Unknown')
            osh_id = ''
        
        # Loop through ERG genes
        for _, erg_gene in erg_df.iterrows():
            erg_scaffold = erg_gene.get('w303_scaffold', '')
            erg_start = erg_gene.get('start', 0)
            erg_end = erg_gene.get('end', 0)
            
            # Get ERG gene name
            erg_name = erg_gene.get('erg_name', '')
            if not erg_name or pd.isna(erg_name):
                if 'sc_gene_id' in erg_gene and pd.notna(erg_gene['sc_gene_id']):
                    erg_id = erg_gene['sc_gene_id']
                    for gene_name, gene_id in ERG_GENES.items():
                        if gene_id == erg_id:
                            erg_name = gene_name
                            break
                    if not erg_name:
                        erg_name = erg_id
                else:
                    erg_name = erg_gene.get('w303_gene_id', 'Unknown')
                    erg_id = ''
            else:
                erg_id = erg_gene.get('sc_gene_id', '')
            
            # Calculate distance only for genes on the same scaffold
            if osh_scaffold == erg_scaffold:
                # Calculate the distance between genes
                if osh_end < erg_start:
                    distance = erg_start - osh_end
                    direction = 'upstream'  # OSH is upstream of ERG
                elif erg_end < osh_start:
                    distance = osh_start - erg_end
                    direction = 'downstream'  # OSH is downstream of ERG
                else:
                    # Genes overlap
                    distance = 0
                    if osh_start <= erg_start and osh_end >= erg_end:
                        direction = 'contains'  # OSH contains ERG
                    elif erg_start <= osh_start and erg_end >= osh_end:
                        direction = 'within'  # OSH is within ERG
                    else:
                        direction = 'overlaps'  # Partial overlap

                # Determine distance zone
                if distance == 0:
                    zone = 'Core'
                elif distance <= DISTANCE_ZONES['Buffer']:
                    zone = 'Buffer'
                elif distance <= DISTANCE_ZONES['Intermediate']:
                    zone = 'Intermediate'
                elif distance <= DISTANCE_ZONES['Satellite']:
                    zone = 'Satellite'
                else:
                    zone = 'Distant'

                # Add the data
                distance_data.append({
                    'osh_gene': osh_name,
                    'osh_id': osh_id,
                    'erg_gene': erg_name,
                    'erg_id': erg_id,
                    'scaffold': osh_scaffold,
                    'distance': distance,
                    'distance_kb': round(distance / 1000, 2),
                    'direction': direction,
                    'zone': zone
                })
            else:
                # Genes on different scaffolds
                distance_data.append({
                    'osh_gene': osh_name,
                    'osh_id': osh_id,
                    'erg_gene': erg_name,
                    'erg_id': erg_id,
                    'scaffold': f"{osh_scaffold}/{erg_scaffold}",
                    'distance': None,
                    'distance_kb': None,
                    'direction': 'different_scaffold',
                    'zone': 'Different Scaffold'
                })

    # Convert to DataFrame
    distance_df = pd.DataFrame(distance_data)
    return distance_df

def create_proximity_heatmap(distance_df, output_dir):
    """Create a heatmap showing the proximity of OSH genes to ERG genes"""
    # Filter to only include genes on the same scaffold
    same_scaffold = distance_df[distance_df['distance'].notna()].copy()

    if same_scaffold.empty:
        print("No OSH and ERG genes found on the same scaffolds.")
        return

    # Create a pivot table of distances (in kb)
    pivot_df = same_scaffold.pivot_table(
        index='osh_gene',
        columns='erg_gene',
        values='distance_kb',
        aggfunc='first'
    ).fillna(1000)  # Use a large value for empty cells

    # Create a custom colormap for distance zones
    cmap = sns.diverging_palette(220, 20, as_cmap=True)

    # Create the heatmap
    plt.figure(figsize=(14, 10))

    # Define a custom normalization to highlight different distance zones
    norm = plt.Normalize(0, 100)  # Adjust as needed

    ax = sns.heatmap(
        pivot_df,
        cmap=cmap,
        norm=norm,
        annot=True,
        fmt=".1f",
        cbar_kws={'label': 'Distance (kb)'}
    )

    plt.title("Genomic Proximity Between OSH and Ergosterol Pathway Genes")
    plt.tight_layout()

    # Save the figure
    output_file = os.path.join(output_dir, "osh_erg_proximity_heatmap.png")
    plt.savefig(output_file, dpi=300)
    plt.close()

    print(f"Proximity heatmap saved to {output_file}")

def create_zone_distribution_plot(distance_df, output_dir):
    """Create a bar plot showing the distribution of OSH genes across distance zones"""
    # Filter to only include genes on the same scaffold
    same_scaffold = distance_df[distance_df['distance'].notna()].copy()

    if same_scaffold.empty:
        print("No OSH and ERG genes found on the same scaffolds.")
        return

    # Count the number of relationships in each zone
    zone_counts = same_scaffold['zone'].value_counts().reset_index()
    zone_counts.columns = ['Zone', 'Count']

    # Ensure all zones are represented (even if count is 0)
    all_zones = pd.DataFrame({
        'Zone': ['Core', 'Buffer', 'Intermediate', 'Satellite', 'Distant'],
        'Order': [1, 2, 3, 4, 5]
    })

    zone_counts = pd.merge(all_zones, zone_counts, on='Zone', how='left').fillna(0)
    zone_counts = zone_counts.sort_values('Order')

    # Create a bar plot
    plt.figure(figsize=(12, 6))

    # Define custom colors for zones
    zone_colors = {
        'Core': '#1f77b4',  # Blue
        'Buffer': '#ff7f0e',  # Orange
        'Intermediate': '#2ca02c',  # Green
        'Satellite': '#d62728',  # Red
        'Distant': '#9467bd'  # Purple
    }

    # Create the bar plot with custom colors
    bars = plt.bar(
        zone_counts['Zone'],
        zone_counts['Count'],
        color=[zone_colors.get(zone, '#777777') for zone in zone_counts['Zone']]
    )

    # Add count labels above bars
    for bar in bars:
        height = bar.get_height()
        plt.text(
            bar.get_x() + bar.get_width()/2.,
            height + 0.1,
            f'{int(height)}',
            ha='center',
            va='bottom'
        )

    plt.title("Distribution of OSH-ERG Gene Relationships by Distance Zone")
    plt.xlabel("Distance Zone")
    plt.ylabel("Number of Gene Relationships")
    plt.tight_layout()

    # Save the figure
    output_file = os.path.join(output_dir, "osh_erg_zone_distribution.png")
    plt.savefig(output_file, dpi=300)
    plt.close()

    print(f"Zone distribution plot saved to {output_file}")

def create_gene_network_graph(distance_df, output_dir):
    """Create a network graph visualizing the relationships between OSH and ERG genes"""
    # Filter to only include genes on the same scaffold
    same_scaffold = distance_df[distance_df['distance'].notna()].copy()

    if same_scaffold.empty:
        print("No OSH and ERG genes found on the same scaffolds.")
        return

    # Create a networkx graph
    G = nx.Graph()

    # Add nodes
    osh_genes = same_scaffold['osh_gene'].unique()
    erg_genes = same_scaffold['erg_gene'].unique()

    for gene in osh_genes:
        G.add_node(gene, type='OSH')

    for gene in erg_genes:
        G.add_node(gene, type='ERG')

    # Add edges based on proximity
    close_relationships = same_scaffold[same_scaffold['zone'].isin(['Core', 'Buffer', 'Intermediate'])]

    for _, row in close_relationships.iterrows():
        G.add_edge(
            row['osh_gene'],
            row['erg_gene'],
            distance=row['distance_kb'],
            zone=row['zone']
        )

    # Create the plot
    plt.figure(figsize=(14, 10))

    # Define positions using spring layout
    pos = nx.spring_layout(G, seed=42)

    # Define node colors based on type
    node_colors = [
        '#1f77b4' if G.nodes[node]['type'] == 'OSH' else '#ff7f0e'
        for node in G.nodes
    ]

    # Define edge colors based on zone
    edge_colors = []
    for u, v, data in G.edges(data=True):
        if data['zone'] == 'Core':
            edge_colors.append('#e41a1c')  # Red
        elif data['zone'] == 'Buffer':
            edge_colors.append('#377eb8')  # Blue
        else:  # Intermediate
            edge_colors.append('#4daf4a')  # Green

    # Draw the network
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=700, alpha=0.8)
    nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold')
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=2, alpha=0.7)

    # Add a legend
    legend_elements = [
        Patch(facecolor='#1f77b4', label='OSH Gene'),
        Patch(facecolor='#ff7f0e', label='ERG Gene'),
        Patch(facecolor='#e41a1c', label='Core Zone'),
        Patch(facecolor='#377eb8', label='Buffer Zone'),
        Patch(facecolor='#4daf4a', label='Intermediate Zone')
    ]

    plt.legend(handles=legend_elements, loc='upper right')
    plt.title("Network of OSH-ERG Gene Relationships")
    plt.axis('off')
    plt.tight_layout()

    # Save the figure
    output_file = os.path.join(output_dir, "osh_erg_network.png")
    plt.savefig(output_file, dpi=300)
    plt.close()

    print(f"Network graph saved to {output_file}")

def analyze_variants_by_distance(variants_df, distance_df, output_dir):
    """Analyze how variants are distributed relative to OSH-ERG gene distances"""
    # Create a dictionary to map scaffold and position to OSH-ERG relationships
    relationship_map = defaultdict(list)

    # Process only relationships on the same scaffold
    same_scaffold = distance_df[distance_df['distance'].notna()].copy()

    # This is a simplified approach - in reality, you'd need more complex
    # logic to accurately map variants to gene relationships, considering
    # gene boundaries and orientations
    for _, row in same_scaffold.iterrows():
        scaffold = row['scaffold']
        key = (scaffold, row['osh_gene'], row['erg_gene'])
        relationship_map[key].append(row)

    # Process variants
    variant_zones = []

    for _, variant in variants_df.iterrows():
        scaffold = variant.get('scaffold', '')
        position = variant.get('position', 0)

        # Check all relationships for this scaffold
        for key, relationships in relationship_map.items():
            if key[0] == scaffold:
                for rel in relationships:
                    # This is a simplified approach - in reality, you'd need to 
                    # consider exact genomic coordinates and gene directions
                    variant_zones.append({
                        'variant_id': variant.get('variant_id', ''),
                        'scaffold': scaffold,
                        'position': position,
                        'osh_gene': rel['osh_gene'],
                        'erg_gene': rel['erg_gene'],
                        'distance_kb': rel['distance_kb'],
                        'zone': rel['zone'],
                        'treatment': variant.get('treatment', ''),
                        'impact': variant.get('impact', '')
                    })

    # Convert to DataFrame
    if variant_zones:
        variant_zones_df = pd.DataFrame(variant_zones)
        
        # Save to file
        output_file = os.path.join(output_dir, "variants_by_osh_erg_distance.tsv")
        variant_zones_df.to_csv(output_file, sep='\t', index=False)
        print(f"Variant zone analysis saved to {output_file}")
        
        # Create visualizations
        create_variant_zone_plots(variant_zones_df, output_dir)
        
        return variant_zones_df
    else:
        print("No variants found in OSH-ERG relationship zones.")
        return pd.DataFrame()

def create_variant_zone_plots(variant_zones_df, output_dir):
    """Create plots showing variant distribution across OSH-ERG distance zones"""
    # Count variants by zone
    zone_counts = variant_zones_df.groupby('zone').size().reset_index(name='count')
    
    # Ensure all zones are represented
    all_zones = pd.DataFrame({
        'zone': ['Core', 'Buffer', 'Intermediate', 'Satellite', 'Distant'],
        'order': [1, 2, 3, 4, 5]
    })
    
    zone_counts = pd.merge(all_zones, zone_counts, on='zone', how='left').fillna(0)
    zone_counts = zone_counts.sort_values('order')
    
    # Create a bar plot
    plt.figure(figsize=(12, 6))
    
    # Define custom colors for zones
    zone_colors = {
        'Core': '#1f77b4',  # Blue
        'Buffer': '#ff7f0e',  # Orange
        'Intermediate': '#2ca02c',  # Green
        'Satellite': '#d62728',  # Red
        'Distant': '#9467bd'  # Purple
    }
    
    # Create the bar plot with custom colors
    bars = plt.bar(
        zone_counts['zone'], 
        zone_counts['count'],
        color=[zone_colors.get(zone, '#777777') for zone in zone_counts['zone']]
    )
    
    # Add count labels above bars
    for bar in bars:
        height = bar.get_height()
        plt.text(
            bar.get_x() + bar.get_width()/2., 
            height + 0.1,
            f'{int(height)}',
            ha='center', 
            va='bottom'
        )
    
    plt.title("Variant Distribution by OSH-ERG Distance Zone")
    plt.xlabel("Distance Zone")
    plt.ylabel("Number of Variants")
    plt.tight_layout()
    
    # Save the figure
    output_file = os.path.join(output_dir, "variants_by_osh_erg_zone.png")
    plt.savefig(output_file, dpi=300)
    plt.close()
    
    print(f"Variant zone distribution plot saved to {output_file}")
    
    # Create a heatmap of variants by zone and treatment
    zone_treatment_counts = variant_zones_df.groupby(['zone', 'treatment']).size().reset_index(name='count')
    
    # Create a pivot table
    pivot_df = zone_treatment_counts.pivot_table(
        index='zone',
        columns='treatment',
        values='count',
        aggfunc='sum'
    ).fillna(0)
    
    # Ensure proper zone order
    if not all(zone in pivot_df.index for zone in ['Core', 'Buffer', 'Intermediate', 'Satellite', 'Distant']):
        # Create a complete index with all zones
        new_index = pd.Index(['Core', 'Buffer', 'Intermediate', 'Satellite', 'Distant'])
        pivot_df = pivot_df.reindex(new_index, fill_value=0)
    
    # Create the heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(
        pivot_df, 
        annot=True, 
        fmt="g", 
        cmap="YlGnBu",
        cbar_kws={'label': 'Variant Count'}
    )
    
    plt.title("Variant Distribution by OSH-ERG Distance Zone and Treatment")
    plt.tight_layout()
    
    # Save the figure
    output_file = os.path.join(output_dir, "zone_treatment_heatmap.png")
    plt.savefig(output_file, dpi=300)
    plt.close()
    
    print(f"Zone-treatment heatmap saved to {output_file}")
    
    # Create a heatmap of variants by zone and impact
    zone_impact_counts = variant_zones_df.groupby(['zone', 'impact']).size().reset_index(name='count')
    
    # Create a pivot table
    pivot_df = zone_impact_counts.pivot_table(
        index='zone',
        columns='impact',
        values='count',
        aggfunc='sum'
    ).fillna(0)
    
    # Ensure proper zone order
    if not all(zone in pivot_df.index for zone in ['Core', 'Buffer', 'Intermediate', 'Satellite', 'Distant']):
        # Create a complete index with all zones
        new_index = pd.Index(['Core', 'Buffer', 'Intermediate', 'Satellite', 'Distant'])
        pivot_df = pivot_df.reindex(new_index, fill_value=0)
    
    # Create the heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(
        pivot_df, 
        annot=True, 
        fmt="g", 
        cmap="YlGnBu",
        cbar_kws={'label': 'Variant Count'}
    )
    
    plt.title("Variant Distribution by OSH-ERG Distance Zone and Impact")
    plt.tight_layout()
    
    # Save the figure
    output_file = os.path.join(output_dir, "zone_impact_heatmap.png")
    plt.savefig(output_file, dpi=300)
    plt.close()
    
    print(f"Zone-impact heatmap saved to {output_file}")

def create_distance_summary(distance_df, output_dir):
    """Create a summary of distances between OSH and ERG genes"""
    # Filter to only include genes on the same scaffold
    same_scaffold = distance_df[distance_df['distance'].notna()].copy()
    
    if same_scaffold.empty:
        print("No OSH and ERG genes found on the same scaffolds.")
        return
    
    # Calculate statistics for each OSH gene
    osh_stats = []
    
    for osh_gene in same_scaffold['osh_gene'].unique():
        osh_data = same_scaffold[same_scaffold['osh_gene'] == osh_gene]
        
        # Calculate closest ERG gene
        if not osh_data.empty:
            min_distance_row = osh_data.loc[osh_data['distance'].idxmin()]
            
            osh_stats.append({
                'OSH Gene': osh_gene,
                'Number of ERG Relationships': len(osh_data),
                'Closest ERG Gene': min_distance_row['erg_gene'],
                'Minimum Distance (kb)': min_distance_row['distance_kb'],
                'Average Distance (kb)': osh_data['distance_kb'].mean(),
                'Core Zone Relationships': sum(osh_data['zone'] == 'Core'),
                'Buffer Zone Relationships': sum(osh_data['zone'] == 'Buffer'),
                'Intermediate Zone Relationships': sum(osh_data['zone'] == 'Intermediate'),
                'Satellite Zone Relationships': sum(osh_data['zone'] == 'Satellite'),
                'Distant Zone Relationships': sum(osh_data['zone'] == 'Distant')
            })
    
    # Convert to DataFrame
    if osh_stats:
        osh_stats_df = pd.DataFrame(osh_stats)
        
        # Save to file
        output_file = os.path.join(output_dir, "osh_erg_distance_summary.tsv")
        osh_stats_df.to_csv(output_file, sep='\t', index=False)
        
        print(f"Distance summary saved to {output_file}")
        
        # Also save a text report
        report_lines = [
            "OSH-ERG Distance Analysis Summary",
            "=================================",
            "",
            f"Total OSH genes analyzed: {len(osh_stats_df)}",
            f"Total ERG genes analyzed: {len(same_scaffold['erg_gene'].unique())}",
            f"Total gene relationships: {len(same_scaffold)}",
            "",
            "Distance Zones:",
            f"  Core (0bp): {sum(same_scaffold['zone'] == 'Core')} relationships",
            f"  Buffer (1-{DISTANCE_ZONES['Buffer']}bp): {sum(same_scaffold['zone'] == 'Buffer')} relationships",
            f"  Intermediate ({DISTANCE_ZONES['Buffer']+1}-{DISTANCE_ZONES['Intermediate']}bp): {sum(same_scaffold['zone'] == 'Intermediate')} relationships",
            f"  Satellite ({DISTANCE_ZONES['Intermediate']+1}-{DISTANCE_ZONES['Satellite']}bp): {sum(same_scaffold['zone'] == 'Satellite')} relationships",
            f"  Distant (>{DISTANCE_ZONES['Satellite']}bp): {sum(same_scaffold['zone'] == 'Distant')} relationships",
            "",
            "OSH Gene Details:",
        ]
        
        for _, row in osh_stats_df.iterrows():
            report_lines.append(f"  {row['OSH Gene']}:")
            report_lines.append(f"    Closest to {row['Closest ERG Gene']} ({row['Minimum Distance (kb)']:.1f} kb)")
            report_lines.append(f"    Average distance to ERG genes: {row['Average Distance (kb)']:.1f} kb")
            report_lines.append(f"    Relationships: {row['Number of ERG Relationships']} ERG genes")
            report_lines.append(f"    Zone breakdown: Core={row['Core Zone Relationships']}, Buffer={row['Buffer Zone Relationships']}, Intermediate={row['Intermediate Zone Relationships']}, Satellite={row['Satellite Zone Relationships']}, Distant={row['Distant Zone Relationships']}")
            report_lines.append("")
        
        # Write the report
        report_file = os.path.join(output_dir, "osh_erg_distance_report.txt")
        with open(report_file, 'w') as f:
            f.write('\n'.join(report_lines))
        
        print(f"Distance report saved to {report_file}")

def main():
    # Parse command line arguments
    args = parse_args()
    
    # Ensure output directory exists
    ensure_dir(args.output_dir)
    
    # Load gene mapping data
    gene_df = load_gene_mapping(args.gene_mapping)
    
    # Extract OSH gene information
    print("Extracting OSH gene family information...")
    osh_df = extract_osh_genes(gene_df)
    print(f"Found {len(osh_df)} OSH family genes out of {len(OSH_GENES)} expected")
    
    # Extract ergosterol pathway gene information
    print("Extracting ergosterol pathway gene information...")
    erg_df = extract_erg_genes(gene_df)
    print(f"Found {len(erg_df)} ergosterol pathway genes out of {len(ERG_GENES)} expected")
    
    # Calculate distances between OSH and ERG genes
    print("Calculating distances between OSH and ERG genes...")
    distance_df = calculate_gene_distances(osh_df, erg_df)
    
    # Save distance data to file
    distance_file = os.path.join(args.output_dir, "osh_erg_distances.tsv")
    distance_df.to_csv(distance_file, sep='\t', index=False)
    print(f"Distance data saved to {distance_file}")
    
    # Create a proximity heatmap
    print("Creating proximity heatmap...")
    create_proximity_heatmap(distance_df, args.output_dir)
    
    # Create a zone distribution plot
    print("Creating zone distribution plot...")
    create_zone_distribution_plot(distance_df, args.output_dir)
    
    # Create a gene network graph
    print("Creating gene network graph...")
    create_gene_network_graph(distance_df, args.output_dir)
    
    # Create a distance summary
    print("Creating distance summary...")
    create_distance_summary(distance_df, args.output_dir)
    
    # Load variant data
    print("Loading variant data...")
    variants_df = load_variant_data(args.variant_dir)
    
    # Analyze variants by distance
    print("Analyzing variants by distance...")
    analyze_variants_by_distance(variants_df, distance_df, args.output_dir)
    
    print("OSH-ERG distance analysis complete!")

if __name__ == "__main__":
    main()