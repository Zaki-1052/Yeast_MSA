#!/usr/bin/env python3
# /Users/zakiralibhai/Documents/GitHub/Yeast_MSA/scripts/osh_analysis/analyze_osh_genes.py

"""
Script to identify and map OSH (OxySterol binding Homology) gene family members
in the reference genome, analyze their sequence characteristics, and identify
their relationship to ergosterol pathway genes.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
import re
import argparse
import sys

# Add the parent directory to the path to access shared modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.tools import ensure_dir

# OSH gene names and standard names
# Values must match the sc_gene_id column in the gene mapping file
OSH_GENES = {
    'OSH2': 'YDL019C',
    'OSH3': 'YHR073W',
    'OSH5': 'YOR237W',
    'OSH6': 'YKR003W',
    'OSH7': 'YHR001W'
}

# Ergosterol pathway genes (from the biochemistry.md documentation)
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

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Map OSH gene family in the reference genome")
    parser.add_argument("--ref-genome", default="/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/reference/w303_chromosomal.fasta",
                      help="Path to reference genome FASTA file")
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

def extract_osh_genes(gene_df):
    """Extract OSH gene family information from gene mapping"""
    # First try exact matches for standard gene names in std_gene_name column
    osh_df = gene_df[gene_df['std_gene_name'].isin(OSH_GENES.values())]

    # If we don't find all genes, try matching on standard gene name
    if len(osh_df) < len(OSH_GENES):
        for osh_name, systematic_id in OSH_GENES.items():
            # If we didn't find this gene by standard gene name, try by name
            if systematic_id not in osh_df['std_gene_name'].values:
                # Try looking by standard gene name
                name_matches = gene_df[gene_df['std_gene_name'] == systematic_id]
                if not name_matches.empty:
                    osh_df = pd.concat([osh_df, name_matches])

                # Also try looking in sc_gene_id column
                id_matches = gene_df[gene_df['sc_gene_id'] == systematic_id]
                if not id_matches.empty:
                    osh_df = pd.concat([osh_df, id_matches])

    # If we still don't have all genes, try fuzzy matching on product
    if len(osh_df) < len(OSH_GENES):
        for osh_name, std_gene_name in OSH_GENES.items():
            # Check if we've already found this gene
            if not any((osh_df['std_gene_name'] == std_gene_name) |
                      (osh_df['std_gene_name'] == osh_name) |
                      (osh_df['sc_gene_id'] == std_gene_name)):
                # Look for oxysterol-binding in product description
                if 'product' in gene_df.columns:
                    desc_matches = gene_df[gene_df['product'].str.contains('oxysterol', case=False, na=False)]
                    if not desc_matches.empty:
                        osh_df = pd.concat([osh_df, desc_matches])

                # Also try partial matching on standard gene name
                for col in ['std_gene_name', 'sc_gene_id', 'w303_gene_id', 'product']:
                    if col in gene_df.columns:
                        # Look for the OSH name or the standard gene name
                        for search_term in [osh_name, std_gene_name]:
                            partial_matches = gene_df[gene_df[col].astype(str).str.contains(search_term, na=False, case=False)]
                            if not partial_matches.empty:
                                print(f"Found {len(partial_matches)} partial matches for {osh_name} in {col}")
                                osh_df = pd.concat([osh_df, partial_matches])

    # Reset index and drop duplicates
    osh_df = osh_df.drop_duplicates().reset_index(drop=True)

    return osh_df

def extract_erg_genes(gene_df):
    """Extract ergosterol pathway gene information from gene mapping"""
    # Find ergosterol genes using sc_gene_id column
    erg_df = gene_df[gene_df['sc_gene_id'].isin(ERG_GENES.values())]

    # If we don't find all genes, try matching on ERG name
    if len(erg_df) < len(ERG_GENES):
        for erg_name, systematic_id in ERG_GENES.items():
            # If we didn't find this gene by systematic ID, try by ERG name
            if systematic_id not in erg_df['sc_gene_id'].values:
                name_matches = gene_df[gene_df['erg_name'] == erg_name]
                if not name_matches.empty:
                    erg_df = pd.concat([erg_df, name_matches])

    # Also find genes where erg_name column matches the ERG gene names
    erg_name_matches = gene_df[gene_df['erg_name'].isin(ERG_GENES.keys())]
    if not erg_name_matches.empty:
        erg_df = pd.concat([erg_df, erg_name_matches])

    # Reset index and drop duplicates
    erg_df = erg_df.drop_duplicates().reset_index(drop=True)

    return erg_df

def calculate_genomic_distances(osh_df, erg_df):
    """Calculate genomic distances between OSH and ergosterol pathway genes"""
    # Create an empty DataFrame to store distances
    distance_data = []

    # Loop through each OSH gene
    for _, osh_row in osh_df.iterrows():
        osh_scaffold = osh_row.get('w303_scaffold', '')
        osh_start = osh_row.get('start', 0)
        osh_end = osh_row.get('end', 0)
        osh_name = osh_row.get('std_gene_name', osh_row.get('sc_gene_id', 'Unknown'))
        if not osh_name or pd.isna(osh_name):
            # If std_gene_name is empty, try to get the name from our mapping
            for gene_name, gene_id in OSH_GENES.items():
                if gene_id == osh_row.get('sc_gene_id', ''):
                    osh_name = gene_name
                    break
            # If still empty, use w303 gene ID
            if not osh_name or pd.isna(osh_name):
                osh_name = osh_row.get('w303_gene_id', 'Unknown')

        # Loop through each ERG gene
        for _, erg_row in erg_df.iterrows():
            erg_scaffold = erg_row.get('w303_scaffold', '')
            erg_start = erg_row.get('start', 0)
            erg_end = erg_row.get('end', 0)
            erg_name = erg_row.get('erg_name', '')
            if not erg_name or pd.isna(erg_name):
                erg_name = erg_row.get('std_gene_name', erg_row.get('sc_gene_id', 'Unknown'))
                if not erg_name or pd.isna(erg_name):
                    # If std_gene_name is empty, try to get the name from our mapping
                    for gene_name, gene_id in ERG_GENES.items():
                        if gene_id == erg_row.get('sc_gene_id', ''):
                            erg_name = gene_name
                            break
                    # If still empty, use w303 gene ID
                    if not erg_name or pd.isna(erg_name):
                        erg_name = erg_row.get('w303_gene_id', 'Unknown')

            # Skip if not on the same scaffold/chromosome
            if osh_scaffold != erg_scaffold:
                distance = None  # Different chromosomes
                distance_kb = None
            else:
                # Calculate distance between genes (minimum distance between any parts)
                if osh_end < erg_start:
                    distance = erg_start - osh_end  # OSH is upstream of ERG
                elif erg_end < osh_start:
                    distance = osh_start - erg_end  # OSH is downstream of ERG
                else:
                    distance = 0  # Genes overlap

                distance_kb = round(distance / 1000, 2) if distance is not None else None

            # Store the data
            distance_data.append({
                'osh_gene': osh_name,
                'osh_systematic_id': osh_row.get('sc_gene_id', ''),
                'erg_gene': erg_name,
                'erg_systematic_id': erg_row.get('sc_gene_id', ''),
                'scaffold': osh_scaffold,
                'distance_bp': distance,
                'distance_kb': distance_kb,
                'same_chromosome': osh_scaffold == erg_scaffold
            })

    # Convert to DataFrame
    distance_df = pd.DataFrame(distance_data)
    return distance_df

def create_distance_heatmap(distance_df, output_dir):
    """Create a heatmap visualization of the distances between OSH and ERG genes"""
    # Filter only genes on the same chromosome for the heatmap
    if 'same_chromosome' in distance_df.columns:
        same_chrom_df = distance_df[distance_df['same_chromosome']].copy()
    elif 'same_scaffold' in distance_df.columns:
        same_chrom_df = distance_df[distance_df['same_scaffold']].copy()
    else:
        print("Warning: No column indicating chromosome/scaffold relationship found.")
        same_chrom_df = pd.DataFrame()

    if same_chrom_df.empty:
        print("No OSH and ERG genes found on the same chromosomes.")
        return

    # Create a pivot table for the heatmap
    pivot_df = same_chrom_df.pivot_table(
        index='osh_gene',
        columns='erg_gene',
        values='distance_kb',
        aggfunc='first'
    )
    
    # Create the heatmap
    plt.figure(figsize=(12, 8))
    mask = pivot_df.isna()
    ax = sns.heatmap(
        pivot_df, 
        annot=True, 
        fmt=".1f", 
        cmap="YlGnBu", 
        mask=mask,
        cbar_kws={'label': 'Distance (kb)'}
    )
    
    plt.title("Genomic Distances Between OSH and Ergosterol Pathway Genes (kb)")
    plt.tight_layout()
    
    # Save the figure
    output_file = os.path.join(output_dir, "osh_erg_distance_heatmap.png")
    plt.savefig(output_file, dpi=300)
    plt.close()
    
    print(f"Heatmap saved to {output_file}")

def analyze_osh_gene_conservation(osh_df, output_dir):
    """Analyze basic conservation characteristics of OSH genes"""
    # Create a summary DataFrame
    summary_data = []

    for _, gene in osh_df.iterrows():
        # Get gene name (prefer standard name, then sc_gene_id, then w303 ID)
        gene_name = gene.get('std_gene_name', '')
        if not gene_name or pd.isna(gene_name):
            # Try to get the name from our mapping
            sc_gene_id = gene.get('sc_gene_id', '')
            for osh_name, osh_id in OSH_GENES.items():
                if osh_id == sc_gene_id:
                    gene_name = osh_name
                    break
            # If still empty, use w303 gene ID
            if not gene_name or pd.isna(gene_name):
                gene_name = gene.get('w303_gene_id', 'Unknown')

        # Get basic gene properties
        summary_data.append({
            'Gene': gene_name,
            'Systematic_ID': gene.get('sc_gene_id', ''),
            'W303_Gene_ID': gene.get('w303_gene_id', ''),
            'Chromosome': gene.get('w303_scaffold', ''),
            'Chromosome_ID': gene.get('chromosome_id', ''),
            'Length': gene.get('end', 0) - gene.get('start', 0),
            'Product': gene.get('product', 'No description available'),
            'Position': f"{gene.get('start', 0)}-{gene.get('end', 0)}",
            'Strand': gene.get('strand', '')
        })

    # Convert to DataFrame
    summary_df = pd.DataFrame(summary_data)

    # Save the summary data
    output_file = os.path.join(output_dir, "osh_gene_summary.tsv")
    summary_df.to_csv(output_file, sep='\t', index=False)

    print(f"OSH gene summary saved to {output_file}")
    return summary_df

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
    
    # Analyze OSH gene conservation
    print("Analyzing OSH gene conservation...")
    osh_summary_df = analyze_osh_gene_conservation(osh_df, args.output_dir)
    
    # Calculate genomic distances between OSH and ERG genes
    print("Calculating genomic distances between OSH and ergosterol pathway genes...")
    distance_df = calculate_genomic_distances(osh_df, erg_df)
    
    # Save distance DataFrame to file
    distance_file = os.path.join(args.output_dir, "osh_erg_distances.tsv")
    distance_df.to_csv(distance_file, sep='\t', index=False)
    print(f"Distance data saved to {distance_file}")
    
    # Create distance heatmap visualization
    print("Creating distance heatmap visualization...")
    create_distance_heatmap(distance_df, args.output_dir)
    
    # Print a summary of the closest relationships
    if 'same_chromosome' in distance_df.columns:
        same_chrom_df = distance_df[distance_df['same_chromosome']].copy()
    elif 'same_scaffold' in distance_df.columns:
        same_chrom_df = distance_df[distance_df['same_scaffold']].copy()
    else:
        print("Warning: No column indicating chromosome/scaffold relationship found.")
        same_chrom_df = pd.DataFrame()

    if not same_chrom_df.empty:
        same_chrom_df = same_chrom_df.sort_values('distance_kb')
        print("\nClosest OSH-ERG gene pairs:")
        for i, (_, row) in enumerate(same_chrom_df.iterrows()):
            if i < 5:  # Show top 5 closest pairs
                scaffold = row.get('scaffold', row.get('osh_scaffold', 'unknown'))
                print(f"  {row['osh_gene']} and {row['erg_gene']}: {row['distance_kb']} kb apart on {scaffold}")
    
    print("OSH gene analysis complete!")

if __name__ == "__main__":
    main()