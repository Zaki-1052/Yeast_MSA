#!/usr/bin/env python3
# /Users/zakiralibhai/Documents/GitHub/Yeast_MSA/scripts/osh_analysis/osh_variants.py

"""
Script to analyze variants in and around OSH genes in all treatment conditions.
This analysis will help understand if OSH genes follow the same conservation pattern
as ergosterol pathway genes, and identify treatment-specific patterns.
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
from scipy import stats

# Add the parent directory to the path to access shared modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.tools import ensure_dir

# OSH gene names and standard names
# Values must match the sc_gene_id column in gene_mapping_full.tsv
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

# Treatment groups
TREATMENTS = {
    'WT-37': 'Temperature-adapted wild type',
    'WTA': 'Low oxygen-adapted wild type',
    'STC': 'STC gene-modified with low oxygen adaptation',
    'CAS': 'CAS gene-modified with temperature adaptation'
}

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Analyze variants in OSH genes")
    parser.add_argument("--gene-mapping", default="/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/reference/gene_mapping_full.tsv",
                      help="Path to gene mapping file")
    parser.add_argument("--variant-dir", default="/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/gene_variants_expanded",
                      help="Directory containing variant data")
    parser.add_argument("--output-dir", default="/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/osh_analysis",
                      help="Directory to store output files")
    parser.add_argument("--distance-threshold", type=int, default=25000,
                      help="Distance threshold for considering variants near genes (bp)")
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

    # If we don't find all genes, try matching on gene ID
    if len(osh_df) < len(OSH_GENES):
        for osh_name, std_gene_name in OSH_GENES.items():
            # Try all possible ways to find the gene
            if std_gene_name not in osh_df['std_gene_name'].values:
                # Try direct match on std_gene_name
                direct_matches = gene_df[gene_df['std_gene_name'] == std_gene_name]
                if not direct_matches.empty:
                    osh_df = pd.concat([osh_df, direct_matches])
                    continue

                # Try partial match on various columns
                for col in ['std_gene_name', 'sc_gene_id', 'w303_gene_id', 'product']:
                    if col in gene_df.columns:
                        try:
                            partial_matches = gene_df[gene_df[col].astype(str).str.contains(std_gene_name, case=False, na=False)]
                            if not partial_matches.empty:
                                print(f"Found {len(partial_matches)} matches for {osh_name} ({std_gene_name}) in {col}")
                                osh_df = pd.concat([osh_df, partial_matches])
                        except Exception as e:
                            print(f"Error searching for {std_gene_name} in {col}: {e}")

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
    # Find ergosterol genes using sc_gene_id column
    erg_df = gene_df[gene_df['sc_gene_id'].isin(ERG_GENES.values())]

    # If we don't find all genes, try matching on ERG name
    if len(erg_df) < len(ERG_GENES):
        for erg_name, systematic_id in ERG_GENES.items():
            # If we didn't find this gene by systematic ID, try by name
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

def load_variant_data(variant_dir):
    """Load variant data from variant files"""
    # Check if the directory exists
    if not os.path.exists(variant_dir):
        print(f"ERROR: Variant directory not found: {variant_dir}")
        sys.exit(1)

    # Define project root directory
    project_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    # Load all scaffold variants file instead of gene-filtered variants
    all_variants_file = os.path.join(project_dir, "results/scaffold_variants/all_scaffold_variants.tsv")
    if not os.path.exists(all_variants_file):
        print(f"ERROR: All scaffold variants file not found: {all_variants_file}")
        print(f"Falling back to gene variants file")

        # Fallback to gene variants file
        all_variants_file = os.path.join(variant_dir, "all_gene_variants.tsv")
        if not os.path.exists(all_variants_file):
            print(f"ERROR: All gene variants file not found: {all_variants_file}")
            sys.exit(1)

    try:
        variants_df = pd.read_csv(all_variants_file, sep="\t")
        print(f"Loaded variant data with {len(variants_df)} variants")

        # Debug: Print column names
        print(f"Variant data columns: {', '.join(variants_df.columns)}")

        # Check for scaffold columns (case-insensitive)
        scaffold_found = False

        # Check for uppercase 'Scaffold' column
        if 'Scaffold' in variants_df.columns:
            print(f"Using existing 'Scaffold' column for scaffold information")
            variants_df['scaffold'] = variants_df['Scaffold']
            scaffold_found = True
        # Check for lowercase 'scaffold' column
        elif 'scaffold' in variants_df.columns:
            print(f"Using existing 'scaffold' column for scaffold information")
            scaffold_found = True
        # Check for w303_scaffold column
        elif 'w303_scaffold' in variants_df.columns:
            print("Using 'w303_scaffold' column for scaffold information")
            variants_df['scaffold'] = variants_df['w303_scaffold']
            scaffold_found = True
        # Try chromosome_id column
        elif 'chromosome_id' in variants_df.columns:
            print("Using 'chromosome_id' column as scaffold information")
            variants_df['scaffold'] = variants_df['chromosome_id']
            scaffold_found = True

        # If no scaffold column found, show a warning
        if not scaffold_found:
            print("WARNING: No scaffold information found in variant data. Using placeholder.")
            variants_df['scaffold'] = "unknown"

        # Standardize column names (adapt scaffold variants file columns to match expected format)
        # Handle Position column
        if 'Position' not in variants_df.columns and 'position' in variants_df.columns:
            variants_df['Position'] = variants_df['position']
        elif 'Position' not in variants_df.columns:
            print("WARNING: No position column found. Using placeholder.")
            variants_df['Position'] = 0

        # Handle Ref column
        if 'Ref' not in variants_df.columns and 'ref' in variants_df.columns:
            variants_df['Ref'] = variants_df['ref']
        elif 'Ref' not in variants_df.columns:
            print("WARNING: No reference column found. Using placeholder.")
            variants_df['Ref'] = ""

        # Handle Alt column
        if 'Alt' not in variants_df.columns and 'alt' in variants_df.columns:
            variants_df['Alt'] = variants_df['alt']
        elif 'Alt' not in variants_df.columns:
            print("WARNING: No alternate column found. Using placeholder.")
            variants_df['Alt'] = ""

        # Handle Treatment column
        if 'Treatment' not in variants_df.columns and 'treatment' in variants_df.columns:
            variants_df['Treatment'] = variants_df['treatment']

        # Handle Impact column
        if 'Impact' not in variants_df.columns and 'Impact' in variants_df.columns:
            variants_df['Impact'] = variants_df['impact']

        # Handle Effect column
        if 'Effect' not in variants_df.columns and 'effect' in variants_df.columns:
            variants_df['Effect'] = variants_df['effect']

        # Print sample values to verify
        if 'scaffold' in variants_df.columns:
            unique_scaffolds = variants_df['scaffold'].unique()
            print(f"Found these scaffold values: {unique_scaffolds[:5]}...")

        return variants_df
    except Exception as e:
        print(f"ERROR: Failed to load variant data: {e}")
        sys.exit(1)

def find_gene_variants(variants_df, gene_df, gene_list, distance_threshold=25000):
    """Find variants in or near specific genes"""
    # Create an empty list to store results
    gene_variants = []

    # Print available columns in gene_df for debugging
    print(f"Available columns in gene_df: {gene_df.columns.tolist()}")

    # Print some sample values from key columns
    for col in ['std_gene_name', 'sc_gene_id', 'gene_id', 'w303_gene_id', 'gene_name']:
        if col in gene_df.columns:
            print(f"Sample values in {col} column: {gene_df[col].iloc[:5].tolist()}")

    # Print values we're searching for
    print(f"Searching for gene IDs: {list(gene_list.values())}")
    print(f"Searching for gene names: {list(gene_list.keys())}")

    # Try to find genes of interest using multiple approaches
    genes_of_interest = pd.DataFrame()

    # Approach 1: Try matching on sc_gene_id
    if 'sc_gene_id' in gene_df.columns:
        sc_matches = gene_df[gene_df['sc_gene_id'].isin(gene_list.values())]
        if not sc_matches.empty:
            print(f"Found {len(sc_matches)} genes by matching sc_gene_id")
            genes_of_interest = pd.concat([genes_of_interest, sc_matches])

    # Approach 2: Try matching on std_gene_name
    if 'std_gene_name' in gene_df.columns:
        std_matches = gene_df[gene_df['std_gene_name'].isin(gene_list.values())]
        if not std_matches.empty:
            print(f"Found {len(std_matches)} genes by matching std_gene_name")
            genes_of_interest = pd.concat([genes_of_interest, std_matches])

    # Approach 3: Try matching gene names against any column that might contain them
    for col in gene_df.columns:
        if col.endswith('_name') or col.endswith('_id') or col == 'description' or col == 'product':
            try:
                # Try matching against gene list keys (gene names)
                name_matches = gene_df[gene_df[col].isin(gene_list.keys())]
                if not name_matches.empty:
                    print(f"Found {len(name_matches)} genes by matching gene names in {col}")
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
                            print(f"Found {len(fuzzy_matches)} genes by fuzzy matching '{value}' in {col}")
                            genes_of_interest = pd.concat([genes_of_interest, fuzzy_matches])
                except Exception as e:
                    pass  # Silently continue if error

    # Remove duplicates
    if not genes_of_interest.empty:
        genes_of_interest = genes_of_interest.drop_duplicates().reset_index(drop=True)

    if genes_of_interest.empty:
        print(f"Warning: No genes of interest found in the gene mapping.")
        return pd.DataFrame()

    print(f"Found {len(genes_of_interest)} genes of interest for analysis")

    # Loop through each gene
    for _, gene in genes_of_interest.iterrows():
        # Get gene name (try multiple fields)
        gene_name = None
        for field in ['std_gene_name', 'erg_name']:
            if field in gene and not pd.isna(gene[field]) and gene[field]:
                gene_name = gene[field]
                break
        if not gene_name:
            # As a last resort, use w303_gene_id or sc_gene_id
            gene_name = gene.get('w303_gene_id', gene.get('sc_gene_id', 'Unknown'))

        sc_gene_id = gene.get('sc_gene_id', '')
        scaffold = gene.get('w303_scaffold', '')
        gene_start = gene.get('start', 0)
        gene_end = gene.get('end', 0)

        # Debug info
        print(f"Checking gene {gene_name} on scaffold {scaffold} from {gene_start}-{gene_end}")

        if not scaffold:
            print(f"Skipping gene {gene_name} with no scaffold information")
            continue

        # Debug what scaffold values are present in the variant data
        # Use a simplified approach for scaffold matching - check both uppercase and lowercase scaffold columns
        mask = False

        # List scaffold columns to check
        scaffold_cols = ['Scaffold', 'scaffold', 'w303_scaffold', 'chromosome_id']

        for scaffold_col in scaffold_cols:
            if scaffold_col in variants_df.columns:
                # Try direct match with the scaffold name
                current_mask = variants_df[scaffold_col] == scaffold
                if any(current_mask):
                    print(f"Found {sum(current_mask)} variants with {scaffold_col}={scaffold}")
                    mask = mask | current_mask

        # If no matches found, print debug info
        if isinstance(mask, bool) and not mask or isinstance(mask, pd.Series) and not any(mask):
            print(f"No variants found on scaffold {scaffold}")
            # Print sample of available scaffolds for debugging
            for scaffold_col in scaffold_cols:
                if scaffold_col in variants_df.columns:
                    print(f"Sample values in {scaffold_col} column:")
                    print(variants_df[scaffold_col].value_counts().head())

        # Get variants matching the scaffold
        if isinstance(mask, bool):
            scaffold_variants = pd.DataFrame() if not mask else variants_df.copy()
        else:
            scaffold_variants = variants_df[mask].copy() if any(mask) else pd.DataFrame()
        print(f"Found {len(scaffold_variants)} variants on scaffold {scaffold}")

        for _, variant in scaffold_variants.iterrows():
            # Get position from Position (uppercase) or position (lowercase)
            variant_pos = variant.get('Position', variant.get('position', 0))

            # Calculate distance to gene
            if variant_pos < gene_start:
                distance = gene_start - variant_pos
                position = "upstream"
            elif variant_pos > gene_end:
                distance = variant_pos - gene_end
                position = "downstream"
            else:
                distance = 0
                position = "within"

            # Check if variant is within threshold distance
            if distance <= distance_threshold:
                gene_variants.append({
                    'gene_name': gene_name,
                    'sc_gene_id': sc_gene_id,
                    'scaffold': scaffold,
                    'gene_start': gene_start,
                    'gene_end': gene_end,
                    'variant_id': variant.get('variant_id', ''),
                    'position': variant_pos,
                    'distance': distance,
                    'location': position,
                    'ref': variant.get('Ref', variant.get('ref', '')),
                    'alt': variant.get('Alt', variant.get('alt', '')),
                    'treatment': variant.get('Treatment', variant.get('treatment', '')),
                    'impact': variant.get('Impact', variant.get('impact', '')),
                    'effect': variant.get('Effect', variant.get('effect', ''))
                })

    # Convert to DataFrame
    if gene_variants:
        result_df = pd.DataFrame(gene_variants)
        return result_df
    else:
        return pd.DataFrame()

def analyze_variant_distribution(osh_variants_df, erg_variants_df, output_dir):
    """Analyze and compare variant distribution between OSH and ERG genes"""
    # Combine data with a gene type column
    osh_variants_df['gene_type'] = 'OSH'
    erg_variants_df['gene_type'] = 'ERG'
    
    combined_df = pd.concat([osh_variants_df, erg_variants_df])
    
    # Count variants by gene type and treatment
    treatment_counts = combined_df.groupby(['gene_type', 'treatment']).size().reset_index(name='variant_count')
    
    # Calculate statistics for OSH vs ERG genes
    osh_count = len(osh_variants_df)
    erg_count = len(erg_variants_df)
    osh_genes_count = len(OSH_GENES)
    erg_genes_count = len(ERG_GENES)
    
    # Average variants per gene
    osh_avg = osh_count / osh_genes_count if osh_genes_count > 0 else 0
    erg_avg = erg_count / erg_genes_count if erg_genes_count > 0 else 0
    
    # Create a bar plot of variant counts by gene type and treatment
    plt.figure(figsize=(12, 6))
    ax = sns.barplot(x='treatment', y='variant_count', hue='gene_type', data=treatment_counts)
    
    # Add counts above bars
    for p in ax.patches:
        ax.annotate(f'{int(p.get_height())}', 
                   (p.get_x() + p.get_width() / 2., p.get_height()), 
                   ha = 'center', va = 'bottom', 
                   xytext = (0, 5), textcoords = 'offset points')
    
    plt.title('Variant Counts by Gene Type and Treatment')
    plt.xlabel('Treatment')
    plt.ylabel('Variant Count')
    plt.tight_layout()
    
    # Save the figure
    output_file = os.path.join(output_dir, "osh_erg_variant_comparison.png")
    plt.savefig(output_file, dpi=300)
    plt.close()
    
    print(f"Variant comparison plot saved to {output_file}")
    
    # Create a summary report
    summary = [
        "OSH vs ERG Gene Variant Analysis",
        "==============================",
        f"Total OSH gene variants: {osh_count}",
        f"Total ERG gene variants: {erg_count}",
        f"OSH genes analyzed: {osh_genes_count}",
        f"ERG genes analyzed: {erg_genes_count}",
        f"Average variants per OSH gene: {osh_avg:.2f}",
        f"Average variants per ERG gene: {erg_avg:.2f}",
        "",
        "Variants by treatment:",
    ]
    
    for gene_type in ['OSH', 'ERG']:
        summary.append(f"\n{gene_type} gene variants:")
        for treatment in TREATMENTS:
            count = treatment_counts[
                (treatment_counts['gene_type'] == gene_type) & 
                (treatment_counts['treatment'] == treatment)
            ]['variant_count'].values
            
            if len(count) > 0:
                summary.append(f"  {treatment}: {count[0]}")
            else:
                summary.append(f"  {treatment}: 0")
    
    # Perform statistical test to compare OSH vs ERG conservation
    if osh_count > 0 and erg_count > 0:
        # Create contingency table
        contingency = np.array([[osh_count, osh_genes_count], 
                               [erg_count, erg_genes_count]])
        
        # Perform Fisher's exact test
        odds_ratio, p_value = stats.fisher_exact(contingency)
        
        summary.append("\nStatistical comparison:")
        summary.append(f"Fisher's exact test p-value: {p_value:.4f}")
        summary.append(f"Odds ratio: {odds_ratio:.4f}")
        
        if p_value < 0.05:
            if odds_ratio > 1:
                summary.append("Conclusion: OSH genes have significantly MORE variants than ERG genes")
            else:
                summary.append("Conclusion: OSH genes have significantly FEWER variants than ERG genes")
        else:
            summary.append("Conclusion: No significant difference in variant counts between OSH and ERG genes")
    
    # Write summary to file
    summary_file = os.path.join(output_dir, "osh_variant_analysis_summary.txt")
    with open(summary_file, 'w') as f:
        f.write('\n'.join(summary))
    
    print(f"Summary report saved to {summary_file}")
    
    return combined_df

def analyze_variant_distance_distribution(combined_df, output_dir):
    """Analyze the distribution of variant distances from genes"""
    # Check if the dataframe is empty
    if combined_df.empty:
        print("No data available for distance distribution analysis.")
        return

    # Create bins for distance
    bins = [0, 1, 1000, 5000, 10000, 25000]
    labels = ['Within gene', '1-1000bp', '1001-5000bp', '5001-10000bp', '10001-25000bp']

    # Create distance category
    try:
        combined_df['distance_category'] = pd.cut(
            combined_df['distance'],
            bins=bins,
            labels=labels,
            include_lowest=True
        )
    except Exception as e:
        print(f"Error creating distance categories: {e}")
        combined_df['distance_category'] = 'Unknown'
    
    # Count variants by gene type, distance category, and treatment
    distance_counts = combined_df.groupby(['gene_type', 'distance_category', 'treatment']).size().reset_index(name='count')
    
    # Plot variant distances
    plt.figure(figsize=(14, 8))
    g = sns.catplot(
        data=distance_counts,
        x='treatment',
        y='count',
        hue='distance_category',
        col='gene_type',
        kind='bar',
        height=6,
        aspect=0.8
    )
    
    g.set_axis_labels("Treatment", "Variant Count")
    g.set_titles("{col_name} Genes")
    g.fig.suptitle('Variant Distance Distribution by Gene Type and Treatment', y=1.05)
    plt.tight_layout()
    
    # Save the figure
    output_file = os.path.join(output_dir, "osh_erg_distance_distribution.png")
    plt.savefig(output_file, dpi=300)
    plt.close()
    
    print(f"Distance distribution plot saved to {output_file}")
    
    # Calculate distance statistics
    within_gene = combined_df[combined_df['distance'] == 0]
    within_1kb = combined_df[combined_df['distance'] <= 1000]
    within_5kb = combined_df[combined_df['distance'] <= 5000]
    
    # Summarize by gene type
    summary_data = []
    
    for gene_type in ['OSH', 'ERG']:
        type_data = combined_df[combined_df['gene_type'] == gene_type]
        
        if not type_data.empty:
            total = len(type_data)
            
            summary_data.append({
                'Gene Type': gene_type,
                'Total Variants': total,
                'Within Gene': len(type_data[type_data['distance'] == 0]),
                'Within Gene %': 100 * len(type_data[type_data['distance'] == 0]) / total if total > 0 else 0,
                'Within 1kb': len(type_data[type_data['distance'] <= 1000]),
                'Within 1kb %': 100 * len(type_data[type_data['distance'] <= 1000]) / total if total > 0 else 0,
                'Within 5kb': len(type_data[type_data['distance'] <= 5000]),
                'Within 5kb %': 100 * len(type_data[type_data['distance'] <= 5000]) / total if total > 0 else 0,
                'Upstream Variants': len(type_data[type_data['location'] == 'upstream']),
                'Upstream %': 100 * len(type_data[type_data['location'] == 'upstream']) / total if total > 0 else 0,
                'Downstream Variants': len(type_data[type_data['location'] == 'downstream']),
                'Downstream %': 100 * len(type_data[type_data['location'] == 'downstream']) / total if total > 0 else 0
            })
    
    # Create a DataFrame and save to file
    summary_df = pd.DataFrame(summary_data)
    summary_file = os.path.join(output_dir, "osh_erg_distance_summary.tsv")
    summary_df.to_csv(summary_file, sep='\t', index=False)
    
    print(f"Distance summary saved to {summary_file}")

def analyze_variant_impact(combined_df, output_dir):
    """Analyze the impact of variants on OSH and ERG genes"""
    # Count variants by gene type, impact, and treatment
    impact_counts = combined_df.groupby(['gene_type', 'impact', 'treatment']).size().reset_index(name='count')
    
    # Plot impact distribution
    plt.figure(figsize=(14, 8))
    g = sns.catplot(
        data=impact_counts,
        x='treatment',
        y='count',
        hue='impact',
        col='gene_type',
        kind='bar',
        height=6,
        aspect=0.8
    )
    
    g.set_axis_labels("Treatment", "Variant Count")
    g.set_titles("{col_name} Genes")
    g.fig.suptitle('Variant Impact Distribution by Gene Type and Treatment', y=1.05)
    plt.tight_layout()
    
    # Save the figure
    output_file = os.path.join(output_dir, "osh_erg_impact_distribution.png")
    plt.savefig(output_file, dpi=300)
    plt.close()
    
    print(f"Impact distribution plot saved to {output_file}")
    
    # Calculate impact statistics
    high_impact = combined_df[combined_df['impact'] == 'HIGH']
    moderate_impact = combined_df[combined_df['impact'] == 'MODERATE']
    low_impact = combined_df[combined_df['impact'] == 'LOW']
    modifier_impact = combined_df[combined_df['impact'] == 'MODIFIER']
    
    # Summarize by gene type
    summary_data = []
    
    for gene_type in ['OSH', 'ERG']:
        type_data = combined_df[combined_df['gene_type'] == gene_type]
        
        if not type_data.empty:
            total = len(type_data)
            
            summary_data.append({
                'Gene Type': gene_type,
                'Total Variants': total,
                'HIGH Impact': len(type_data[type_data['impact'] == 'HIGH']),
                'HIGH Impact %': 100 * len(type_data[type_data['impact'] == 'HIGH']) / total if total > 0 else 0,
                'MODERATE Impact': len(type_data[type_data['impact'] == 'MODERATE']),
                'MODERATE Impact %': 100 * len(type_data[type_data['impact'] == 'MODERATE']) / total if total > 0 else 0,
                'LOW Impact': len(type_data[type_data['impact'] == 'LOW']),
                'LOW Impact %': 100 * len(type_data[type_data['impact'] == 'LOW']) / total if total > 0 else 0,
                'MODIFIER Impact': len(type_data[type_data['impact'] == 'MODIFIER']),
                'MODIFIER Impact %': 100 * len(type_data[type_data['impact'] == 'MODIFIER']) / total if total > 0 else 0
            })
    
    # Create a DataFrame and save to file
    summary_df = pd.DataFrame(summary_data)
    summary_file = os.path.join(output_dir, "osh_erg_impact_summary.tsv")
    summary_df.to_csv(summary_file, sep='\t', index=False)
    
    print(f"Impact summary saved to {summary_file}")

def analyze_treatment_effects(combined_df, output_dir):
    """Analyze treatment-specific effects on OSH genes"""
    # Check if dataframe is empty
    if combined_df.empty:
        print("No data available for treatment effect analysis")
        # Create empty placeholder files
        summary_file = os.path.join(output_dir, "osh_treatment_summary.tsv")
        with open(summary_file, 'w') as f:
            f.write("No data available for analysis\n")
        print(f"Empty treatment summary saved to {summary_file}")
        return

    # Get counts by gene, gene type, and treatment
    gene_treatment_counts = combined_df.groupby(['gene_name', 'gene_type', 'treatment']).size().reset_index(name='count')

    # Check if we have any OSH data
    osh_data = gene_treatment_counts[gene_treatment_counts['gene_type'] == 'OSH']
    if osh_data.empty:
        print("No OSH gene data available for heatmap")
        # Create a simple placeholder figure
        plt.figure(figsize=(10, 6))
        plt.text(0.5, 0.5, "No OSH gene variants found",
                 horizontalalignment='center', verticalalignment='center',
                 fontsize=14)
        plt.axis('off')
    else:
        # Pivot the data for heatmap
        osh_pivot = osh_data.pivot(
            index='gene_name',
            columns='treatment',
            values='count'
        ).fillna(0)

        # Create a heatmap
        plt.figure(figsize=(12, 8))
        ax = sns.heatmap(
            osh_pivot,
            annot=True,
            fmt="g",
            cmap="YlGnBu",
            cbar_kws={'label': 'Variant Count'}
        )

        plt.title("OSH Gene Variants by Treatment")

    plt.tight_layout()

    # Save the figure
    output_file = os.path.join(output_dir, "osh_treatment_heatmap.png")
    plt.savefig(output_file, dpi=300)
    plt.close()

    print(f"Treatment heatmap saved to {output_file}")

    # Calculate treatment-specific statistics - handle empty data gracefully
    temperature_adaptations = combined_df[combined_df['treatment'].isin(['WT-37', 'CAS'])] if not combined_df.empty else pd.DataFrame()
    low_oxygen_adaptations = combined_df[combined_df['treatment'].isin(['WTA', 'STC'])] if not combined_df.empty else pd.DataFrame()

    gene_modified = combined_df[combined_df['treatment'].isin(['CAS', 'STC'])] if not combined_df.empty else pd.DataFrame()
    non_modified = combined_df[combined_df['treatment'].isin(['WT-37', 'WTA'])] if not combined_df.empty else pd.DataFrame()

    # Summarize by adaptation type
    summary_data = []

    # Add OSH data by adaptation type
    osh_data = combined_df[combined_df['gene_type'] == 'OSH'] if not combined_df.empty else pd.DataFrame()
    erg_data = combined_df[combined_df['gene_type'] == 'ERG'] if not combined_df.empty else pd.DataFrame()

    if not osh_data.empty:
        total_osh = len(osh_data)

        # Define groups to check
        groups = [
            ('Adaptation Type', 'Temperature (WT-37, CAS)', ['WT-37', 'CAS']),
            ('Adaptation Type', 'Low Oxygen (WTA, STC)', ['WTA', 'STC']),
            ('Gene Modification', 'Modified (CAS, STC)', ['CAS', 'STC']),
            ('Gene Modification', 'Non-modified (WT-37, WTA)', ['WT-37', 'WTA'])
        ]

        for category, group_name, treatments in groups:
            osh_in_group = len(osh_data[osh_data['treatment'].isin(treatments)])
            erg_in_group = len(erg_data[erg_data['treatment'].isin(treatments)]) if not erg_data.empty else 0

            summary_data.append({
                'Category': category,
                'Group': group_name,
                'OSH Variants': osh_in_group,
                'OSH %': 100 * osh_in_group / total_osh if total_osh > 0 else 0,
                'ERG Variants': erg_in_group,
                'Ratio OSH/ERG': osh_in_group / erg_in_group if erg_in_group > 0 else 0
            })
    else:
        # Add placeholder data
        summary_data.append({
            'Category': 'Analysis',
            'Group': 'All Data',
            'OSH Variants': 0,
            'OSH %': 0,
            'ERG Variants': len(erg_data),
            'Ratio OSH/ERG': 0
        })

    # Create a DataFrame and save to file
    summary_df = pd.DataFrame(summary_data)
    summary_file = os.path.join(output_dir, "osh_treatment_summary.tsv")
    summary_df.to_csv(summary_file, sep='\t', index=False)

    print(f"Treatment summary saved to {summary_file}")

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
    
    # Extract ergosterol pathway gene information for comparison
    print("Extracting ergosterol pathway gene information...")
    erg_df = extract_erg_genes(gene_df)
    print(f"Found {len(erg_df)} ergosterol pathway genes out of {len(ERG_GENES)} expected")
    
    # Load variant data
    print("Loading variant data...")
    variants_df = load_variant_data(args.variant_dir)
    
    # Find variants in/near OSH genes
    print(f"Finding variants within {args.distance_threshold}bp of OSH genes...")
    osh_variants_df = find_gene_variants(
        variants_df, 
        osh_df, 
        OSH_GENES, 
        distance_threshold=args.distance_threshold
    )
    print(f"Found {len(osh_variants_df)} variants near OSH genes")
    
    # Find variants in/near ERG genes for comparison
    print(f"Finding variants within {args.distance_threshold}bp of ergosterol pathway genes...")
    erg_variants_df = find_gene_variants(
        variants_df, 
        erg_df, 
        ERG_GENES, 
        distance_threshold=args.distance_threshold
    )
    print(f"Found {len(erg_variants_df)} variants near ergosterol pathway genes")
    
    # Save variant data to files
    osh_variants_file = os.path.join(args.output_dir, "osh_variants.tsv")
    osh_variants_df.to_csv(osh_variants_file, sep='\t', index=False)
    print(f"OSH variants saved to {osh_variants_file}")
    
    # Analyze variant distribution
    print("Analyzing variant distribution between OSH and ERG genes...")
    combined_df = analyze_variant_distribution(osh_variants_df, erg_variants_df, args.output_dir)
    
    # Analyze variant distance distribution
    print("Analyzing variant distance distribution...")
    analyze_variant_distance_distribution(combined_df, args.output_dir)
    
    # Analyze variant impact
    print("Analyzing variant impact...")
    analyze_variant_impact(combined_df, args.output_dir)
    
    # Analyze treatment effects
    print("Analyzing treatment-specific effects...")
    analyze_treatment_effects(combined_df, args.output_dir)
    
    print("OSH variant analysis complete!")

if __name__ == "__main__":
    main()