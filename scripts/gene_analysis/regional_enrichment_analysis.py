#!/usr/bin/env python3

'''
Regional Enrichment Analysis with Gene Mapping

This script analyzes the regional enrichment of mutations in yeast strains,
with a focus on gene-specific mapping and analysis. It identifies regions with
statistically significant enrichment of variants, maps them to genes, and
produces visualizations and reports of the enrichment patterns with gene context.

The script adds gene-specific functionality to the original regional_enrichment_analysis.py,
allowing for the identification of enriched regions within genes, particularly
focusing on genes involved in the ergosterol biosynthesis pathway which may be
under purifying selection.
'''

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict, Counter
from scipy.stats import poisson, ttest_ind
from scipy.cluster import hierarchy
import subprocess
import re
import warnings
import logging
warnings.filterwarnings('ignore')

# Set up logging for file only, not console
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("regional_enrichment_debug.log")
    ]
)

# Create scaffold matching dictionary to handle ID format differences
def create_scaffold_mapping(mutation_chroms, scaffold_ids):
    """
    Create a mapping between mutation chromosome IDs and reference scaffold IDs.
    
    Args:
        mutation_chroms: List of chromosome IDs from mutation data
        scaffold_ids: List of scaffold IDs from reference
    
    Returns:
        Dictionary mapping mutation chromosome IDs to reference scaffold IDs
    """
    mapping = {}
    
    # Try exact matching first
    for chrom in mutation_chroms:
        if chrom in scaffold_ids:
            mapping[chrom] = chrom
    
    # If no exact matches, try more flexible matching
    if not mapping:
        logging.warning("No exact matches between mutation chromosomes and reference scaffolds")
        logging.info("Attempting flexible matching...")
        
        # Try to extract patterns from chromosome IDs
        chrom_pattern = None
        if len(mutation_chroms) > 0:
            sample_chrom = mutation_chroms[0]
            # Check if it matches a known pattern (e.g., JRIU01000xxx.1)
            match = re.match(r'([A-Za-z]+\d+)(\d+)(\.\d+)?', sample_chrom)
            if match:
                chrom_pattern = match.group(1)
                logging.info(f"Detected chromosome pattern: {chrom_pattern}")
        
        # Try to extract patterns from scaffold IDs
        scaffold_pattern = None
        if len(scaffold_ids) > 0:
            sample_scaffold = scaffold_ids[0]
            # Check if it matches a known pattern
            match = re.match(r'([A-Za-z]+\d+)(\d+)(\.\d+)?', sample_scaffold)
            if match:
                scaffold_pattern = match.group(1)
                logging.info(f"Detected scaffold pattern: {scaffold_pattern}")
        
        # If both patterns are detected and they are different, try mapping
        if chrom_pattern and scaffold_pattern and chrom_pattern != scaffold_pattern:
            logging.info(f"Attempting to map between patterns: {chrom_pattern} -> {scaffold_pattern}")
            
            # For each chromosome, try to find a matching scaffold
            for chrom in mutation_chroms:
                # Extract the numeric part
                chrom_match = re.match(r'[A-Za-z]+\d+(\d+)(\.\d+)?', chrom)
                if chrom_match:
                    chrom_num = chrom_match.group(1)
                    # Look for a scaffold with the same numeric part
                    for scaffold in scaffold_ids:
                        scaffold_match = re.match(r'[A-Za-z]+\d+(\d+)(\.\d+)?', scaffold)
                        if scaffold_match and scaffold_match.group(1) == chrom_num:
                            mapping[chrom] = scaffold
                            break
        
        # If still no matches, try ignoring any suffix like ".1"
        if not mapping:
            logging.info("Attempting to match by ignoring suffixes...")
            for chrom in mutation_chroms:
                base_chrom = chrom.split('.')[0]
                for scaffold in scaffold_ids:
                    base_scaffold = scaffold.split('.')[0]
                    if base_chrom == base_scaffold:
                        mapping[chrom] = scaffold
                        break
    
    logging.info(f"Created mapping for {len(mapping)} out of {len(mutation_chroms)} chromosomes")
    return mapping

# Set plotting style
plt.style.use('ggplot')
sns.set(font_scale=1.2)
sns.set_style("whitegrid")

# Define output directories
OUTPUT_DIR = "analysis/regional_enrichment_results"
GENE_OUTPUT_DIR = "analysis/genes_of_interest/regional_enrichment_results"
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(GENE_OUTPUT_DIR, exist_ok=True)

# Initialize gene data structures
GENE_DATA = {}  # Dictionary mapping gene IDs to their details
SCAFFOLD_GENES = defaultdict(list)  # Dictionary mapping scaffolds to lists of genes
GENES_OF_INTEREST = set()  # Set of gene IDs involved in the ergosterol pathway

# Updated biologically correct treatment groups
TREATMENTS = ['WT-37', 'WTA', 'STC', 'CAS']

# Adaptation types
ADAPTATION_TYPES = ['Temperature', 'Low Oxygen']

# Define adaptation colors for consistent visualization
ADAPTATION_COLORS = {
    'Temperature': '#1b9e77',  # Temperature adaptation
    'Low Oxygen': '#d95f02',  # Low oxygen adaptation
}

# Define treatment information for better biological context
TREATMENT_INFO = {
    'WT-37': {'description': 'Temperature-adapted wild type', 'adaptation': 'Temperature'},
    'WTA': {'description': 'Low oxygen-adapted wild type', 'adaptation': 'Low Oxygen'},
    'STC': {'description': 'STC gene with low oxygen adaptation', 'adaptation': 'Low Oxygen', 'gene': 'STC'},
    'CAS': {'description': 'CAS gene with temperature adaptation', 'adaptation': 'Temperature', 'gene': 'CAS'}
}

# Treatment colors for consistent visualization
TREATMENT_COLORS = {
    'WT-37': '#1b9e77',  # Temperature-adapted
    'WTA': '#d95f02',    # Low oxygen-adapted
    'STC': '#7570b3',    # STC gene + low oxygen
    'CAS': '#e7298a'     # CAS gene + temperature
}

# Function to load gene mapping data
def load_gene_mapping():
    """
    Load gene mapping data from reference files.
    
    This function loads gene data from the reference directory, including:
    1. Gene coordinates and information from gene_mapping.tsv
    2. Genes of interest (ergosterol pathway) from genes_of_interest_mapping.tsv
    
    The function populates three global data structures:
    - GENE_DATA: Dictionary mapping gene IDs to their details
    - SCAFFOLD_GENES: Dictionary mapping scaffolds to lists of genes
    - GENES_OF_INTEREST: Set of ergosterol pathway gene IDs
    
    Returns:
        bool: True if data was loaded successfully, False otherwise
    """
    global GENE_DATA, SCAFFOLD_GENES, GENES_OF_INTEREST
    
    # Clear existing data
    GENE_DATA.clear()
    SCAFFOLD_GENES.clear()
    GENES_OF_INTEREST.clear()
    
    # Define possible file paths for gene mapping data
    gene_mapping_paths = [
        "reference/gene_mapping.tsv",
        "reference/w303_annotations/gene_mapping.tsv"
    ]
    
    # Load gene mapping data
    gene_mapping_file = None
    for path in gene_mapping_paths:
        if os.path.exists(path):
            gene_mapping_file = path
            break
    
    if gene_mapping_file:
        try:
            # Load gene mapping data
            gene_df = pd.read_csv(gene_mapping_file, sep='\t')
            print(f"Loaded {len(gene_df)} genes from {gene_mapping_file}")
            logging.info(f"Loaded {len(gene_df)} genes from {gene_mapping_file}")
            
            # Process each gene
            for _, row in gene_df.iterrows():
                gene_id = row['w303_gene_id']
                scaffold = row['w303_scaffold']
                
                # Store gene data
                GENE_DATA[gene_id] = {
                    'gene_id': gene_id,
                    'locus_tag': row['locus_tag'] if 'locus_tag' in row else None,
                    'sc_gene_id': row['sc_gene_id'] if 'sc_gene_id' in row else None,
                    'erg_name': row['erg_name'] if 'erg_name' in row else None,
                    'scaffold': scaffold,
                    'start': int(row['start']),
                    'end': int(row['end']),
                    'strand': row['strand'] if 'strand' in row else None,
                    'product': row['product'] if 'product' in row else None
                }
                
                # Map scaffold to genes
                SCAFFOLD_GENES[scaffold].append(gene_id)
            
            # Load genes of interest (ergosterol pathway genes)
            goi_paths = [
                "reference/genes_of_interest_mapping.tsv",
                "reference/w303_annotations/genes_of_interest.tsv"
            ]
            
            goi_file = None
            for path in goi_paths:
                if os.path.exists(path):
                    goi_file = path
                    break
            
            if goi_file:
                try:
                    goi_df = pd.read_csv(goi_file, sep='\t')
                    print(f"Loaded {len(goi_df)} genes of interest from {goi_file}")
                    logging.info(f"Loaded {len(goi_df)} genes of interest from {goi_file}")
                    
                    # Add to our set of genes of interest
                    for _, row in goi_df.iterrows():
                        if 'w303_gene_id' in row:
                            GENES_OF_INTEREST.add(row['w303_gene_id'])
                except Exception as e:
                    print(f"Error loading genes of interest: {e}")
                    logging.error(f"Error loading genes of interest: {e}")
            else:
                print("No genes of interest file found. Using empty set.")
                logging.warning("No genes of interest file found.")
            
            return True
        except Exception as e:
            print(f"Error loading gene mapping data: {e}")
            logging.error(f"Error loading gene mapping data: {e}")
            return False
    else:
        print("No gene mapping file found. Gene mapping will not be available.")
        logging.warning("No gene mapping file found. Gene mapping will not be available.")
        return False

# Function to map variants to genes
def map_variants_to_genes(variant_df):
    """
    Map variants to genes based on their genomic coordinates.
    
    This function takes a DataFrame of variants with scaffold and position information
    and maps each variant to the corresponding gene if it falls within a gene's coordinates.
    Variants are also categorized as being in genes of interest or not.
    
    Args:
        variant_df (pandas.DataFrame): DataFrame containing variant information with at minimum
                                      'CHROM' and 'POS' columns
    
    Returns:
        pandas.DataFrame: The original DataFrame with additional gene-related columns:
                         - in_gene: Boolean indicating if variant is in a gene
                         - gene_id: Gene identifier (if in_gene)
                         - gene_name: Gene name (if available)
                         - gene_type: 'ergosterol' or 'other' (if in_gene)
                         - gene_product: Gene product description (if available)
    """
    if len(GENE_DATA) == 0 or len(SCAFFOLD_GENES) == 0:
        print("Gene mapping data not loaded. Cannot map variants to genes.")
        logging.warning("Gene mapping data not loaded. Cannot map variants to genes.")
        return variant_df
    
    # Create a copy of the input dataframe
    result_df = variant_df.copy()
    
    # Initialize gene-related columns
    result_df['in_gene'] = False
    result_df['gene_id'] = None
    result_df['gene_name'] = None
    result_df['gene_type'] = None
    result_df['gene_product'] = None
    
    # Map each variant to genes
    for idx, row in result_df.iterrows():
        scaffold = row['CHROM']
        position = row['POS']
        
        # Skip if scaffold has no mapped genes
        if scaffold not in SCAFFOLD_GENES:
            continue
        
        # Check each gene in this scaffold
        for gene_id in SCAFFOLD_GENES[scaffold]:
            gene_data = GENE_DATA[gene_id]
            
            # Check if position falls within gene coordinates
            if gene_data['start'] <= position <= gene_data['end']:
                result_df.at[idx, 'in_gene'] = True
                result_df.at[idx, 'gene_id'] = gene_id
                
                # Add gene name if available
                if gene_data['erg_name']:
                    result_df.at[idx, 'gene_name'] = gene_data['erg_name']
                elif gene_data['sc_gene_id']:
                    result_df.at[idx, 'gene_name'] = gene_data['sc_gene_id']
                
                # Set gene type based on presence in genes of interest
                if gene_id in GENES_OF_INTEREST:
                    result_df.at[idx, 'gene_type'] = 'ergosterol'
                else:
                    result_df.at[idx, 'gene_type'] = 'other'
                
                # Add gene product description if available
                if gene_data['product']:
                    result_df.at[idx, 'gene_product'] = gene_data['product']
                
                # Break since we found a matching gene
                break
    
    # Log the mapping results
    in_gene_count = sum(result_df['in_gene'])
    ergosterol_count = sum(result_df['gene_type'] == 'ergosterol')
    
    print(f"Mapped {in_gene_count} out of {len(result_df)} variants to genes")
    print(f"Found {ergosterol_count} variants in ergosterol pathway genes")
    
    logging.info(f"Mapped {in_gene_count} out of {len(result_df)} variants to genes")
    logging.info(f"Found {ergosterol_count} variants in ergosterol pathway genes")
    
    return result_df

# Function to analyze gene-specific enrichment or depletion patterns
def analyze_gene_specific_enrichment(data, scaffold_info):
    """
    Analyze enrichment or depletion of variants within genes, with special consideration for genes of interest.
    
    This function analyzes the distribution of variants within genes to identify
    genes with mutation rates that differ from expected. Special attention is given
    to ergosterol pathway genes, which may be under purifying selection (showing negative
    enrichment/depletion of variants).
    
    Args:
        data (pandas.DataFrame): DataFrame containing variant information with
                               gene mapping information
        scaffold_info (dict): Dictionary mapping scaffolds to lengths
    
    Returns:
        dict: Dictionary with gene-specific enrichment results:
              - 'all_genes': DataFrame of all genes with their enrichment values
              - 'enriched_genes': DataFrame of significantly enriched genes
              - 'depleted_genes': DataFrame of significantly depleted genes
              - 'erg_genes': DataFrame of ergosterol pathway genes
              - 'non_erg_genes': DataFrame of non-ergosterol genes
              - 'top_enriched_genes': DataFrame of top enriched genes
              - 'top_depleted_genes': DataFrame of top depleted genes
    """
    if data is None or not scaffold_info or 'in_gene' not in data.columns:
        logging.error("Cannot perform gene-specific enrichment analysis: missing data or gene mapping")
        return {}
    
    # Calculate genome-wide mutation rate
    total_mutations = len(data)
    total_length = sum(scaffold_info.values())
    genome_wide_rate = total_mutations / total_length
    
    print(f"Analyzing gene-specific enrichment using {len(data)} variants")
    print(f"Genome-wide mutation rate: {genome_wide_rate:.8f} mutations per base")
    logging.info(f"Analyzing gene-specific enrichment using {len(data)} variants")
    logging.info(f"Genome-wide mutation rate: {genome_wide_rate:.8f} mutations per base")
    
    # Get gene-specific variant counts
    gene_variants = data[data['in_gene']].groupby('gene_id').size().reset_index(name='variant_count')
    
    # Calculate enrichment or depletion for each gene
    analyzed_genes = []
    
    # Also analyze genes with no variants (for potential purifying selection)
    # First, create a set of all genes that have variants
    genes_with_variants = set(gene_variants['gene_id'])
    
    # Analyze both genes with variants and genes without variants
    for gene_id in GENE_DATA.keys():
        gene_info = GENE_DATA[gene_id]
        gene_length = gene_info['end'] - gene_info['start'] + 1
        expected_variants = genome_wide_rate * gene_length
        
        # Get variant count (0 if gene has no variants)
        if gene_id in genes_with_variants:
            variant_count = gene_variants[gene_variants['gene_id'] == gene_id]['variant_count'].iloc[0]
        else:
            variant_count = 0
        
        # Calculate fold enrichment/depletion
        # Values > 1 indicate enrichment, values < 1 indicate depletion/purifying selection
        fold_enrichment = variant_count / expected_variants if expected_variants > 0 else 0
        
        # Calculate log2 fold change for better visualization of both enrichment and depletion
        if fold_enrichment > 0:
            log2_fold_change = np.log2(fold_enrichment)
        else:
            log2_fold_change = float('-inf')
        
        # Calculate p-value using Poisson distribution
        # For variant_count > expected: 1 - CDF of (variant_count - 1) gives p-value for enrichment
        # For variant_count < expected: CDF gives p-value for depletion
        if variant_count > expected_variants:
            # Testing for enrichment (more variants than expected)
            p_value = 1 - poisson.cdf(variant_count - 1, expected_variants)
            direction = 'enriched'
        else:
            # Testing for depletion (fewer variants than expected - purifying selection)
            p_value = poisson.cdf(variant_count, expected_variants)
            direction = 'depleted'
        
        # Extract treatment info for variants in this gene
        gene_data = data[data['gene_id'] == gene_id]
        if len(gene_data) > 0:
            treatments = gene_data['Treatment'].value_counts().to_dict()
            adaptations = gene_data['Adaptation'].value_counts().to_dict()
            
            # Format as strings for easier reporting
            treatment_str = '; '.join([f"{t}: {count} ({count / len(gene_data) * 100:.1f}%)" for t, count in treatments.items()])
            adaptation_str = '; '.join([f"{a}: {count} ({count / len(gene_data) * 100:.1f}%)" for a, count in adaptations.items()])
        else:
            treatment_str = ""
            adaptation_str = ""
        
        # Create gene record
        gene_record = {
            'gene_id': gene_id,
            'gene_name': gene_info.get('erg_name') or gene_info.get('sc_gene_id') or gene_id,
            'scaffold': gene_info['scaffold'],
            'start': gene_info['start'],
            'end': gene_info['end'],
            'length': gene_length,
            'variant_count': variant_count,
            'expected_variants': expected_variants,
            'fold_enrichment': fold_enrichment,
            'log2_fold_change': log2_fold_change,
            'p_value': p_value,
            'direction': direction,
            'is_erg_gene': gene_id in GENES_OF_INTEREST,
            'treatments': treatment_str,
            'adaptations': adaptation_str
        }
        
        analyzed_genes.append(gene_record)
    
    # Convert to DataFrame
    if analyzed_genes:
        enriched_df = pd.DataFrame(analyzed_genes)
        
        # Multiple testing correction using FDR
        if len(enriched_df) > 1:
            from statsmodels.stats.multitest import multipletests
            try:
                _, corrected_pvals, _, _ = multipletests(
                    enriched_df['p_value'], 
                    method='fdr_bh',
                    alpha=0.5  # Using a lenient threshold 
                )
                enriched_df['q_value'] = corrected_pvals
            except Exception as e:
                logging.warning(f"Error in multiple testing correction: {e}. Using uncorrected p-values.")
                enriched_df['q_value'] = enriched_df['p_value']
        else:
            enriched_df['q_value'] = enriched_df['p_value']
        
        # Get significantly enriched genes
        sig_enriched = enriched_df[(enriched_df['direction'] == 'enriched') & (enriched_df['q_value'] < 0.5)].copy()
        
        # Get significantly depleted genes (under purifying selection)
        sig_depleted = enriched_df[(enriched_df['direction'] == 'depleted') & (enriched_df['q_value'] < 0.5)].copy()
        
        # Get ergosterol genes
        erg_genes = enriched_df[enriched_df['is_erg_gene']].copy()
        
        # Get non-ergosterol genes
        non_erg_genes = enriched_df[~enriched_df['is_erg_gene']].copy()
        
        # Get top 20 most enriched and depleted genes
        top_enriched = enriched_df[enriched_df['direction'] == 'enriched'].sort_values('fold_enrichment', ascending=False).head(20).copy()
        top_depleted = enriched_df[enriched_df['direction'] == 'depleted'].sort_values('fold_enrichment').head(20).copy()
        
        print(f"Found {len(enriched_df)} total genes analyzed")
        print(f"Found {len(sig_enriched)} significantly enriched genes")
        print(f"Found {len(sig_depleted)} significantly depleted genes (under purifying selection)")
        print(f"Found {len(erg_genes)} ergosterol pathway genes with data")
        
        logging.info(f"Found {len(enriched_df)} total genes analyzed")
        logging.info(f"Found {len(sig_enriched)} significantly enriched genes")
        logging.info(f"Found {len(sig_depleted)} significantly depleted genes (under purifying selection)")
        logging.info(f"Found {len(erg_genes)} ergosterol pathway genes with data")
        
        return {
            'all_genes': enriched_df,
            'enriched_genes': sig_enriched,
            'depleted_genes': sig_depleted,
            'erg_genes': erg_genes,
            'non_erg_genes': non_erg_genes,
            'top_enriched_genes': top_enriched,
            'top_depleted_genes': top_depleted
        }
    else:
        print("No genes with variants found")
        logging.warning("No genes with variants found")
        return {}

# Function to plot gene-specific enrichment patterns
def plot_gene_enrichment(gene_results, output_dir):
    """
    Generate visualizations of gene-specific enrichment patterns.
    
    Args:
        gene_results (dict): Dictionary containing gene-specific enrichment results
        output_dir (str): Directory to save the visualizations
    """
    if not gene_results or not gene_results.get('all_genes') or len(gene_results['all_genes']) == 0:
        logging.warning("No gene enrichment data to visualize")
        return
    
    all_genes_df = gene_results['all_genes']
    
    # Plot 1: Log2 fold change by gene type (shows both enrichment and depletion)
    plt.figure(figsize=(12, 8))
    
    # Create a new column for gene type
    all_genes_df['gene_category'] = all_genes_df['is_erg_gene'].apply(
        lambda x: 'Ergosterol Pathway' if x else 'Other Genes'
    )
    
    # Create boxplot with swarmplot overlay for log2 fold change by gene type
    ax = sns.boxplot(
        x='gene_category', 
        y='log2_fold_change',
        data=all_genes_df,
        palette={'Ergosterol Pathway': '#1b9e77', 'Other Genes': '#d95f02'}
    )
    
    # Add swarmplot for individual genes
    sns.swarmplot(
        x='gene_category', 
        y='log2_fold_change',
        data=all_genes_df,
        color='black',
        alpha=0.5,
        size=4
    )
    
    # Add a horizontal line at y=0 (no change)
    plt.axhline(y=0, color='red', linestyle='--', alpha=0.7)
    
    plt.title('Log2 Fold Change by Gene Type (Positive = Enriched, Negative = Depleted)')
    plt.xlabel('Gene Type')
    plt.ylabel('Log2 Fold Change')
    
    # Add text annotation for purifying selection
    plt.text(
        0.02, 0.02, 
        "Negative values indicate potential purifying selection", 
        transform=plt.gca().transAxes,
        fontsize=10, 
        bbox=dict(facecolor='white', alpha=0.7)
    )
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'gene_type_log2fc.png'), dpi=300)
    plt.close()
    
    # Plot 2: Top genes by fold enrichment
    plt.figure(figsize=(12, 8))
    
    # Check if we have enriched genes
    if 'top_enriched_genes' in gene_results and len(gene_results['top_enriched_genes']) > 0:
        top_genes = gene_results['top_enriched_genes']
        
        # Create bar plot with gene type coloring
        colors = ['#1b9e77' if is_erg else '#d95f02' for is_erg in top_genes['is_erg_gene']]
        
        ax = sns.barplot(
            x='gene_name', 
            y='fold_enrichment', 
            data=top_genes,
            palette=colors
        )
        
        # Add variant count as text on bars
        for i, (_, row) in enumerate(top_genes.iterrows()):
            ax.text(
                i, 
                row['fold_enrichment'] + 0.2, 
                f"{int(row['variant_count'])}", 
                ha='center'
            )
        
        plt.xticks(rotation=90)
        plt.title('Top Genes by Fold Enrichment')
        plt.xlabel('Gene')
        plt.ylabel('Fold Enrichment')
        
        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#1b9e77', label='Ergosterol Pathway'),
            Patch(facecolor='#d95f02', label='Other Genes')
        ]
        plt.legend(handles=legend_elements)
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'top_enriched_genes.png'), dpi=300)
        plt.close()
    
    # Plot 3: Top genes showing depletion (potential purifying selection)
    plt.figure(figsize=(12, 8))
    
    # Check if we have depleted genes
    if 'top_depleted_genes' in gene_results and len(gene_results['top_depleted_genes']) > 0:
        top_depleted = gene_results['top_depleted_genes']
        
        # Create bar plot with gene type coloring for depletion (using negative log2 fold change)
        colors = ['#1b9e77' if is_erg else '#d95f02' for is_erg in top_depleted['is_erg_gene']]
        
        ax = sns.barplot(
            x='gene_name', 
            y='log2_fold_change', 
            data=top_depleted,
            palette=colors
        )
        
        # Add variant count as text on bars
        for i, (_, row) in enumerate(top_depleted.iterrows()):
            ax.text(
                i, 
                row['log2_fold_change'] - 0.2, 
                f"{int(row['variant_count'])}", 
                ha='center',
                va='top',
                color='white' if row['log2_fold_change'] < -2 else 'black'
            )
        
        plt.xticks(rotation=90)
        plt.title('Top Genes with Potential Purifying Selection (Depleted Variants)')
        plt.xlabel('Gene')
        plt.ylabel('Log2 Fold Change')
        
        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#1b9e77', label='Ergosterol Pathway'),
            Patch(facecolor='#d95f02', label='Other Genes')
        ]
        plt.legend(handles=legend_elements)
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'genes_under_purifying_selection.png'), dpi=300)
        plt.close()
    
    # Plot 4: Ergosterol pathway genes comparison
    if 'erg_genes' in gene_results and len(gene_results['erg_genes']) > 0:
        erg_genes_df = gene_results['erg_genes']
        
        plt.figure(figsize=(14, 8))
        
        # Sort by log2 fold change to show both enrichment and depletion
        erg_genes_df = erg_genes_df.sort_values('log2_fold_change', ascending=False)
        
        # Create horizontal bar plot of log2 fold change
        ax = sns.barplot(
            x='log2_fold_change',
            y='gene_name',
            data=erg_genes_df,
            palette=erg_genes_df['log2_fold_change'].apply(
                lambda x: '#1b9e77' if x > 0 else '#d95f02'  # Green for enriched, orange for depleted
            )
        )
        
        # Add variant count as text on bars
        for i, (_, row) in enumerate(erg_genes_df.iterrows()):
            text_color = 'black'
            ha_position = 'left' if row['log2_fold_change'] < 0 else 'right'
            x_offset = -0.2 if row['log2_fold_change'] < 0 else 0.2
            
            ax.text(
                row['log2_fold_change'] + x_offset, 
                i,
                f"{int(row['variant_count'])} variants", 
                va='center',
                ha=ha_position,
                color=text_color
            )
        
        # Add a vertical line at x=0 (no change)
        plt.axvline(x=0, color='black', linestyle='--', alpha=0.7)
        
        plt.title('Ergosterol Pathway Genes: Enrichment vs. Depletion')
        plt.xlabel('Log2 Fold Change')
        plt.ylabel('Gene')
        
        # Add explanatory text
        plt.figtext(
            0.5, 0.01, 
            "Positive values indicate enrichment, negative values indicate depletion (potential purifying selection)", 
            ha='center',
            fontsize=10,
            bbox=dict(facecolor='white', alpha=0.7)
        )
        
        plt.tight_layout(rect=[0, 0.03, 1, 0.97])  # Adjust for the text at the bottom
        plt.savefig(os.path.join(output_dir, 'ergosterol_genes_enrichment_depletion.png'), dpi=300)
        plt.close()
        
        # Plot 5: Statistical comparison between ERG and non-ERG genes
        plt.figure(figsize=(10, 6))
        
        # Prepare data for comparison
        erg_log2fc = erg_genes_df['log2_fold_change'].dropna()
        
        if 'non_erg_genes' in gene_results and len(gene_results['non_erg_genes']) > 0:
            non_erg_genes_df = gene_results['non_erg_genes']
            non_erg_log2fc = non_erg_genes_df['log2_fold_change'].dropna()
            
            # Create violin plot for comparison
            comparison_data = pd.DataFrame({
                'Gene Type': ['Ergosterol Pathway'] * len(erg_log2fc) + ['Other Genes'] * len(non_erg_log2fc),
                'Log2 Fold Change': list(erg_log2fc) + list(non_erg_log2fc)
            })
            
            sns.violinplot(
                x='Gene Type',
                y='Log2 Fold Change',
                data=comparison_data,
                palette={'Ergosterol Pathway': '#1b9e77', 'Other Genes': '#d95f02'},
                inner='quartile'
            )
            
            # Add points for individual genes
            sns.stripplot(
                x='Gene Type',
                y='Log2 Fold Change',
                data=comparison_data,
                color='black',
                alpha=0.3,
                jitter=True
            )
            
            # Add a horizontal line at y=0 (no change)
            plt.axhline(y=0, color='red', linestyle='--', alpha=0.7)
            
            # Perform statistical test
            t_stat, p_val = ttest_ind(erg_log2fc, non_erg_log2fc, equal_var=False)
            
            # Add statistical annotation
            plt.title(f'Comparison: Ergosterol vs Other Genes (p-value: {p_val:.4f})')
            
            # Add more detailed statistics
            erg_mean = erg_log2fc.mean()
            non_erg_mean = non_erg_log2fc.mean()
            
            plt.figtext(
                0.5, 0.01,
                f"Mean Log2FC: ERG={erg_mean:.2f}, Non-ERG={non_erg_mean:.2f}, Difference={erg_mean-non_erg_mean:.2f}",
                ha='center',
                fontsize=10,
                bbox=dict(facecolor='white', alpha=0.7)
            )
            
            plt.xlabel('Gene Type')
            plt.ylabel('Log2 Fold Change')
            plt.tight_layout(rect=[0, 0.03, 1, 0.97])  # Adjust for the text at the bottom
            plt.savefig(os.path.join(output_dir, 'erg_vs_nonerg_statistical_comparison.png'), dpi=300)
            plt.close()

# Function to create gene-specific reports
def create_gene_specific_reports(gene_results, output_dir):
    """
    Create detailed reports of gene-specific enrichment analysis.
    
    Args:
        gene_results (dict): Dictionary containing gene-specific enrichment results
        output_dir (str): Directory to save the reports
    """
    if not gene_results or not gene_results.get('all_genes') or len(gene_results['all_genes']) == 0:
        logging.warning("No gene enrichment data for reporting")
        return
    
    # Save all genes enrichment data
    all_genes_df = gene_results['all_genes']
    all_genes_df.to_csv(os.path.join(output_dir, 'all_genes_analysis.csv'), index=False)
    
    # Save significantly enriched genes
    if 'enriched_genes' in gene_results and len(gene_results['enriched_genes']) > 0:
        gene_results['enriched_genes'].to_csv(os.path.join(output_dir, 'significantly_enriched_genes.csv'), index=False)
    
    # Save significantly depleted genes (under purifying selection)
    if 'depleted_genes' in gene_results and len(gene_results['depleted_genes']) > 0:
        gene_results['depleted_genes'].to_csv(os.path.join(output_dir, 'genes_under_purifying_selection.csv'), index=False)
    
    # Save ergosterol pathway genes
    if 'erg_genes' in gene_results and len(gene_results['erg_genes']) > 0:
        gene_results['erg_genes'].to_csv(os.path.join(output_dir, 'ergosterol_pathway_genes.csv'), index=False)
    
    # Save non-ergosterol genes
    if 'non_erg_genes' in gene_results and len(gene_results['non_erg_genes']) > 0:
        gene_results['non_erg_genes'].to_csv(os.path.join(output_dir, 'non_ergosterol_genes.csv'), index=False)
    
    # Create a comprehensive text report
    with open(os.path.join(output_dir, 'gene_enrichment_summary.txt'), 'w') as f:
        f.write("Gene-Specific Enrichment Analysis Summary\n")
        f.write("=========================================\n\n")
        
        # Overall statistics
        f.write("Overall Statistics:\n")
        f.write("-----------------\n")
        
        total_genes = len(all_genes_df)
        enriched_genes = len(gene_results.get('enriched_genes', pd.DataFrame()))
        depleted_genes = len(gene_results.get('depleted_genes', pd.DataFrame()))
        erg_genes = len(gene_results.get('erg_genes', pd.DataFrame()))
        
        f.write(f"Total genes analyzed: {total_genes}\n")
        f.write(f"Significantly enriched genes: {enriched_genes}\n")
        f.write(f"Significantly depleted genes (potential purifying selection): {depleted_genes}\n")
        f.write(f"Ergosterol pathway genes analyzed: {erg_genes}\n\n")
        
        # Ergosterol pathway genes statistics
        if 'erg_genes' in gene_results and len(gene_results['erg_genes']) > 0:
            erg_genes_df = gene_results['erg_genes']
            
            # Count enriched and depleted ERG genes
            erg_enriched = len(erg_genes_df[erg_genes_df['log2_fold_change'] > 0])
            erg_depleted = len(erg_genes_df[erg_genes_df['log2_fold_change'] < 0])
            
            f.write("Ergosterol Pathway Gene Statistics:\n")
            f.write("--------------------------------\n")
            f.write(f"Ergosterol genes with enrichment: {erg_enriched}\n")
            f.write(f"Ergosterol genes with depletion: {erg_depleted}\n")
            f.write(f"Average log2 fold change: {erg_genes_df['log2_fold_change'].mean():.2f}\n\n")
            
            # Statistical comparison with non-ERG genes if available
            if 'non_erg_genes' in gene_results and len(gene_results['non_erg_genes']) > 0:
                non_erg_genes_df = gene_results['non_erg_genes']
                
                # Perform t-test
                erg_log2fc = erg_genes_df['log2_fold_change'].dropna()
                non_erg_log2fc = non_erg_genes_df['log2_fold_change'].dropna()
                
                if len(erg_log2fc) > 0 and len(non_erg_log2fc) > 0:
                    t_stat, p_val = ttest_ind(erg_log2fc, non_erg_log2fc, equal_var=False)
                    
                    f.write("Statistical Comparison (ERG vs. Non-ERG genes):\n")
                    f.write("-----------------------------------------\n")
                    f.write(f"ERG genes average log2FC: {erg_log2fc.mean():.2f}\n")
                    f.write(f"Non-ERG genes average log2FC: {non_erg_log2fc.mean():.2f}\n")
                    f.write(f"Difference: {erg_log2fc.mean() - non_erg_log2fc.mean():.2f}\n")
                    f.write(f"T-test p-value: {p_val:.4f}\n\n")
                    
                    # Interpretation
                    if p_val < 0.05:
                        if erg_log2fc.mean() < non_erg_log2fc.mean():
                            f.write("Interpretation: Ergosterol genes show significantly lower mutation rates\n")
                            f.write("compared to other genes, suggesting purifying selection.\n\n")
                        else:
                            f.write("Interpretation: Ergosterol genes show significantly higher mutation rates\n")
                            f.write("compared to other genes, suggesting potential selection pressure.\n\n")
                    else:
                        f.write("Interpretation: No significant difference in mutation rates between\n")
                        f.write("ergosterol genes and other genes.\n\n")
        
        # Top enriched genes
        if 'top_enriched_genes' in gene_results and len(gene_results['top_enriched_genes']) > 0:
            f.write("Top 10 Enriched Genes:\n")
            f.write("------------------\n")
            f.write(f"{'Gene Name':<15}{'Variants':<10}{'Expected':<10}{'Fold Enrichment':<20}{'Log2FC':<10}{'Q-value':<15}\n")
            f.write(f"{'-' * 70}\n")
            
            for _, row in gene_results['top_enriched_genes'].head(10).iterrows():
                f.write(f"{row['gene_name']:<15}{int(row['variant_count']):<10}{row['expected_variants']:.2f:<10}{row['fold_enrichment']:<20.2f}{row['log2_fold_change']:<10.2f}{row['q_value']:<15.2e}\n")
            
            f.write("\n")
        
        # Top depleted genes (purifying selection)
        if 'top_depleted_genes' in gene_results and len(gene_results['top_depleted_genes']) > 0:
            f.write("Top 10 Genes Under Purifying Selection (Depleted):\n")
            f.write("------------------------------------------\n")
            f.write(f"{'Gene Name':<15}{'Variants':<10}{'Expected':<10}{'Fold Enrichment':<20}{'Log2FC':<10}{'Q-value':<15}\n")
            f.write(f"{'-' * 70}\n")
            
            for _, row in gene_results['top_depleted_genes'].head(10).iterrows():
                f.write(f"{row['gene_name']:<15}{int(row['variant_count']):<10}{row['expected_variants']:.2f:<10}{row['fold_enrichment']:<20.2f}{row['log2_fold_change']:<10.2f}{row['q_value']:<15.2e}\n")
            
            f.write("\n")
        
        # Ergosterol genes detailed analysis
        if 'erg_genes' in gene_results and len(gene_results['erg_genes']) > 0:
            f.write("Ergosterol Pathway Genes Analysis:\n")
            f.write("-----------------------------\n")
            f.write(f"{'Gene Name':<15}{'Variants':<10}{'Expected':<10}{'Fold Enrichment':<20}{'Log2FC':<10}{'Status':<15}\n")
            f.write(f"{'-' * 70}\n")
            
            for _, row in gene_results['erg_genes'].iterrows():
                status = "Enriched" if row['log2_fold_change'] > 0 else "Depleted"
                f.write(f"{row['gene_name']:<15}{int(row['variant_count']):<10}{row['expected_variants']:.2f:<10}{row['fold_enrichment']:<20.2f}{row['log2_fold_change']:<10.2f}{status:<15}\n")
            
            f.write("\n")
        
        # Main conclusions
        f.write("Main Conclusions:\n")
        f.write("---------------\n")
        
        # Analyze ERG genes for overall trend
        if 'erg_genes' in gene_results and len(gene_results['erg_genes']) > 0:
            erg_genes_df = gene_results['erg_genes']
            erg_mean_log2fc = erg_genes_df['log2_fold_change'].mean()
            
            if erg_mean_log2fc < 0:
                f.write("1. Ergosterol pathway genes show evidence of purifying selection with\n")
                f.write("   overall lower mutation rates than expected by chance.\n")
            elif erg_mean_log2fc > 0:
                f.write("1. Ergosterol pathway genes show higher mutation rates than expected,\n")
                f.write("   suggesting potential selection pressure on this pathway.\n")
            else:
                f.write("1. Ergosterol pathway genes show mutation rates similar to the genomic average.\n")
            
            # Check if there's a mix of enrichment and depletion
            erg_enriched = len(erg_genes_df[erg_genes_df['log2_fold_change'] > 0])
            erg_depleted = len(erg_genes_df[erg_genes_df['log2_fold_change'] < 0])
            
            if erg_enriched > 0 and erg_depleted > 0:
                f.write("2. There is heterogeneity within the ergosterol pathway, with some genes\n")
                f.write("   showing enrichment and others showing depletion of variants.\n")
            
            # Compare with other genes
            if 'non_erg_genes' in gene_results and len(gene_results['non_erg_genes']) > 0:
                non_erg_genes_df = gene_results['non_erg_genes']
                non_erg_mean_log2fc = non_erg_genes_df['log2_fold_change'].mean()
                
                f.write(f"3. Compared to other genes (mean log2FC: {non_erg_mean_log2fc:.2f}), ergosterol genes\n")
                f.write(f"   (mean log2FC: {erg_mean_log2fc:.2f}) show ")
                
                if erg_mean_log2fc < non_erg_mean_log2fc:
                    f.write("stronger purifying selection.\n")
                else:
                    f.write("less purifying selection.\n")
        else:
            f.write("1. Limited data available for ergosterol pathway genes analysis.\n")
        
        f.write("4. For detailed gene-level analysis, refer to the CSV files and visualizations\n")
        f.write("   in the output directory.\n")

# Function to create gene-specific summary report
def create_gene_specific_summary(gene_results, data, output_dir):
    """
    Create a comprehensive gene-specific enrichment analysis summary.
    
    Args:
        gene_results (dict): Dictionary containing gene-specific enrichment results
        data (pandas.DataFrame): The original variant data with gene mapping
        output_dir (str): Directory to save the summary report
    """
    # Create report file
    with open(os.path.join(output_dir, 'gene_analysis_summary.txt'), 'w') as f:
        f.write("Gene-Level Mutation Analysis Summary\n")
        f.write("=================================\n\n")
        
        # Overall statistics
        total_variants = len(data)
        in_gene_variants = sum(data['in_gene'] if 'in_gene' in data.columns else 0)
        erg_gene_variants = sum(data['gene_type'] == 'ergosterol' if 'gene_type' in data.columns else 0)
        
        f.write("Overall Statistics:\n")
        f.write("-----------------\n")
        f.write(f"Total variants analyzed: {total_variants}\n")
        
        if 'in_gene' in data.columns:
            pct_in_genes = in_gene_variants / total_variants * 100 if total_variants > 0 else 0
            f.write(f"Variants in genes: {in_gene_variants} ({pct_in_genes:.1f}%)\n")
            
            if 'gene_type' in data.columns:
                pct_in_erg = erg_gene_variants / total_variants * 100 if total_variants > 0 else 0
                f.write(f"Variants in ergosterol pathway genes: {erg_gene_variants} ({pct_in_erg:.1f}%)\n")
                
                if in_gene_variants > 0:
                    pct_of_gene_variants = erg_gene_variants / in_gene_variants * 100
                    f.write(f"Percent of genic variants in ergosterol pathway: {pct_of_gene_variants:.1f}%\n")
        
        f.write("\n")
        
        # Treatment-specific statistics
        f.write("Treatment-Specific Statistics:\n")
        f.write("--------------------------\n")
        
        for treatment in TREATMENTS:
            treatment_data = data[data['Treatment'] == treatment]
            treatment_count = len(treatment_data)
            if treatment_count == 0:
                continue
                
            treatment_in_gene = sum(treatment_data['in_gene'] if 'in_gene' in treatment_data.columns else 0)
            treatment_in_erg = sum(treatment_data['gene_type'] == 'ergosterol' if 'gene_type' in treatment_data.columns else 0)
            
            pct_in_genes = treatment_in_gene / treatment_count * 100 if treatment_count > 0 else 0
            pct_in_erg = treatment_in_erg / treatment_count * 100 if treatment_count > 0 else 0
            
            f.write(f"\n{treatment} Treatment:\n")
            f.write(f"  Total variants: {treatment_count}\n")
            
            if 'in_gene' in treatment_data.columns:
                f.write(f"  Variants in genes: {treatment_in_gene} ({pct_in_genes:.1f}%)\n")
                
                if 'gene_type' in treatment_data.columns:
                    f.write(f"  Variants in ergosterol genes: {treatment_in_erg} ({pct_in_erg:.1f}%)\n")
                    
                    if treatment_in_gene > 0:
                        pct_of_gene_variants = treatment_in_erg / treatment_in_gene * 100
                        f.write(f"  Percent of genic variants in ergosterol pathway: {pct_of_gene_variants:.1f}%\n")
        
        f.write("\n")
        
        # Gene enrichment statistics
        if gene_results and 'all_genes' in gene_results and len(gene_results['all_genes']) > 0:
            all_genes_df = gene_results['all_genes']
            erg_genes_df = gene_results.get('erg_genes', pd.DataFrame())
            
            f.write("Gene Enrichment Statistics:\n")
            f.write("------------------------\n")
            f.write(f"Genes analyzed: {len(all_genes_df)}\n")
            
            # Count genes with enrichment vs depletion
            enriched_count = len(all_genes_df[all_genes_df['log2_fold_change'] > 0])
            depleted_count = len(all_genes_df[all_genes_df['log2_fold_change'] < 0])
            pct_enriched = enriched_count / len(all_genes_df) * 100 if len(all_genes_df) > 0 else 0
            pct_depleted = depleted_count / len(all_genes_df) * 100 if len(all_genes_df) > 0 else 0
            
            f.write(f"Genes with enrichment (log2FC > 0): {enriched_count} ({pct_enriched:.1f}%)\n")
            f.write(f"Genes with depletion (log2FC < 0): {depleted_count} ({pct_depleted:.1f}%)\n")
            
            if 'enriched_genes' in gene_results:
                f.write(f"Significantly enriched genes: {len(gene_results['enriched_genes'])}\n")
            
            if 'depleted_genes' in gene_results:
                f.write(f"Significantly depleted genes: {len(gene_results['depleted_genes'])}\n")
            
            if len(erg_genes_df) > 0:
                f.write(f"Ergosterol pathway genes analyzed: {len(erg_genes_df)}\n")
                
                # Count ERG genes with enrichment vs depletion
                erg_enriched = len(erg_genes_df[erg_genes_df['log2_fold_change'] > 0])
                erg_depleted = len(erg_genes_df[erg_genes_df['log2_fold_change'] < 0])
                
                f.write(f"Ergosterol genes with enrichment: {erg_enriched}\n")
                f.write(f"Ergosterol genes with depletion: {erg_depleted}\n")
                
                # Average fold changes
                avg_erg_log2fc = erg_genes_df['log2_fold_change'].mean()
                f.write(f"Average log2 fold change in ergosterol genes: {avg_erg_log2fc:.2f}\n")
                
                # Compare with non-ergosterol genes
                non_erg_df = all_genes_df[~all_genes_df['is_erg_gene']]
                if len(non_erg_df) > 0:
                    avg_non_erg_log2fc = non_erg_df['log2_fold_change'].mean()
                    f.write(f"Average log2 fold change in non-ergosterol genes: {avg_non_erg_log2fc:.2f}\n")
                    
                    # Statistical comparison
                    if len(erg_genes_df) > 5 and len(non_erg_df) > 5:
                        t_stat, p_val = ttest_ind(
                            erg_genes_df['log2_fold_change'].dropna(), 
                            non_erg_df['log2_fold_change'].dropna(),
                            equal_var=False
                        )
                        f.write(f"Statistical comparison (t-test): p-value = {p_val:.4f}\n")
                        
                        # Interpretation
                        if p_val < 0.05:
                            if avg_erg_log2fc < avg_non_erg_log2fc:
                                f.write("  Interpretation: Ergosterol genes show significantly stronger purifying selection\n")
                                f.write("  compared to other genes.\n")
                            else:
                                f.write("  Interpretation: Ergosterol genes show significantly less purifying selection\n")
                                f.write("  compared to other genes.\n")
            
            f.write("\n")
            
            # Top enriched and depleted genes
            if 'top_enriched_genes' in gene_results and len(gene_results['top_enriched_genes']) > 0:
                f.write("Top 5 Most Enriched Genes:\n")
                for _, row in gene_results['top_enriched_genes'].head(5).iterrows():
                    f.write(f"  {row['gene_name']} - {int(row['variant_count'])} variants, {row['fold_enrichment']:.2f}x enrichment")
                    if row['is_erg_gene']:
                        f.write(" (ergosterol pathway)")
                    f.write("\n")
                
                f.write("\n")
            
            if 'top_depleted_genes' in gene_results and len(gene_results['top_depleted_genes']) > 0:
                f.write("Top 5 Genes Under Strongest Purifying Selection:\n")
                for _, row in gene_results['top_depleted_genes'].head(5).iterrows():
                    f.write(f"  {row['gene_name']} - {int(row['variant_count'])} variants, {row['log2_fold_change']:.2f} log2FC")
                    if row['is_erg_gene']:
                        f.write(" (ergosterol pathway)")
                    f.write("\n")
                
                f.write("\n")
        
        # Main conclusions
        f.write("Main Conclusions:\n")
        f.write("---------------\n")
        f.write("1. This analysis reveals the distribution of variants within genes, with special\n")
        f.write("   attention to purifying selection in ergosterol pathway genes.\n")
        
        if 'in_gene' in data.columns and total_variants > 0:
            if pct_in_genes > 50:
                f.write("2. The majority of variants are located within gene regions.\n")
            else:
                f.write("2. Most variants are located in non-genic regions.\n")
        
        if gene_results and 'erg_genes' in gene_results and len(gene_results['erg_genes']) > 0:
            erg_genes_df = gene_results['erg_genes']
            avg_erg_log2fc = erg_genes_df['log2_fold_change'].mean()
            
            # Conclusion about ergosterol genes
            if avg_erg_log2fc < 0:
                f.write("3. Ergosterol pathway genes show evidence of purifying selection, with lower\n")
                f.write("   mutation rates than expected by chance.\n")
            else:
                f.write("3. Ergosterol pathway genes do not show clear evidence of purifying selection\n")
                f.write("   across the entire pathway.\n")
            
            # Check for mixed pattern
            erg_enriched = len(erg_genes_df[erg_genes_df['log2_fold_change'] > 0])
            erg_depleted = len(erg_genes_df[erg_genes_df['log2_fold_change'] < 0])
            
            if erg_enriched > 0 and erg_depleted > 0:
                f.write("4. There is heterogeneity within the ergosterol pathway, with some genes showing\n")
                f.write("   enrichment and others showing depletion of variants.\n")
        
        # Look for treatment patterns in genes
        if 'in_gene' in data.columns and len(TREATMENTS) > 1:
            treatment_gene_counts = {}
            treatment_erg_counts = {}
            
            for treatment in TREATMENTS:
                treatment_data = data[data['Treatment'] == treatment]
                if len(treatment_data) > 0:
                    treatment_gene_counts[treatment] = sum(treatment_data['in_gene']) / len(treatment_data)
                    
                    if 'gene_type' in treatment_data.columns:
                        erg_count = sum(treatment_data['gene_type'] == 'ergosterol')
                        treatment_erg_counts[treatment] = erg_count / len(treatment_data)
            
            # See if there's variation between treatments
            if treatment_gene_counts and max(treatment_gene_counts.values()) - min(treatment_gene_counts.values()) > 0.1:
                f.write("5. Different treatments show varying proportions of mutations in genes,\n")
                f.write("   suggesting treatment-specific effects on mutation localization.\n")
            
            # Check if there are differences in ERG gene targeting
            if treatment_erg_counts and max(treatment_erg_counts.values()) - min(treatment_erg_counts.values()) > 0.05:
                highest_treatment = max(treatment_erg_counts, key=treatment_erg_counts.get)
                lowest_treatment = min(treatment_erg_counts, key=treatment_erg_counts.get)
                
                f.write(f"6. Treatment {highest_treatment} shows the highest proportion of mutations in ergosterol genes,\n")
                f.write(f"   while treatment {lowest_treatment} shows the lowest, suggesting differential\n")
                f.write("   effects on this pathway.\n")
        
        f.write("7. For detailed gene-level analysis, refer to the visualization files and CSV reports\n")
        f.write("   in the output directory.\n")

# Function to load mutation data
def load_data():
    """
    Load mutation data from pre-processed files or directly from VCF files.
    
    Returns:
        pandas.DataFrame: A DataFrame containing mutation data
    """
    all_data = []
    
    # Try to load from mutation_spectrum_analysis files first
    for treatment in TREATMENTS:
        file_patterns = [
            f"mutation_spectrum_analysis/{treatment}_mutations.txt",
            f"analysis/mutation_spectrum_analysis/{treatment}_mutations.txt"
        ]
        
        file_found = False
        for file_path in file_patterns:
            if os.path.exists(file_path):
                try:
                    print(f"Found file for {treatment} at {file_path}")
                    
                    # Try to determine the file structure
                    with open(file_path, 'r') as f:
                        first_line = f.readline().strip()
                        columns = first_line.split('\t')
                        if len(columns) == 5:
                            print(f"File structure: 5 columns in first line")
                            # Typical format: CHROM, POS, REF, ALT, Treatment
                            data = pd.read_csv(file_path, sep='\t', header=None,
                                               names=['CHROM', 'POS', 'REF', 'ALT', 'Treatment'])
                        else:
                            print(f"File structure: {len(columns)} columns in first line")
                            # Try to adapt to whatever format is there
                            data = pd.read_csv(file_path, sep='\t')
                    
                    print(f"Loaded {len(data)} mutations for {treatment}")
                    
                    # Ensure we use the correct treatment name regardless of file name
                    data['Treatment'] = treatment
                    all_data.append(data)
                    file_found = True
                    break
                except Exception as e:
                    print(f"Error reading {file_path}: {e}")
                    continue
        
        if not file_found:
            print(f"No mutation data found for {treatment}")
    
    if all_data:
        combined_data = pd.concat(all_data, ignore_index=True)
        print(f"Loaded {len(combined_data)} mutations across {len(all_data)} treatments")
        return combined_data
    else:
        print("No mutation data found.")
        return None

def main():
    # Load data
    print("Loading mutation data...")
    logging.info("Loading mutation data...")
    mutation_data = load_data()
    
    print("Loading scaffold information...")
    logging.info("Loading scaffold information...")
    scaffold_info = load_scaffold_info()
    
    print("Loading gene mapping data...")
    logging.info("Loading gene mapping data...")
    gene_mapping_loaded = load_gene_mapping()
    
    if mutation_data is None or not scaffold_info:
        print("Cannot proceed with analysis due to missing data")
        logging.error("Cannot proceed with analysis due to missing data")
        return
    
    # Log key statistics
    logging.info(f"Loaded {len(mutation_data)} mutations across {len(mutation_data['Treatment'].unique())} treatments")
    logging.info(f"Loaded {len(scaffold_info)} scaffolds from reference")
    if gene_mapping_loaded:
        logging.info(f"Loaded {len(GENE_DATA)} genes and {len(GENES_OF_INTEREST)} genes of interest")
    
    # Map variants to genes if gene mapping was loaded
    if gene_mapping_loaded:
        print("Mapping variants to genes...")
        logging.info("Mapping variants to genes...")
        mutation_data = map_variants_to_genes(mutation_data)
    
    # Perform treatment-specific enrichment analysis with adaptive parameters
    print("\nAnalyzing regional enrichment by treatment...")
    logging.info("\nAnalyzing regional enrichment by treatment...")
    treatment_results = {}
    
    for treatment in TREATMENTS:
        print(f"\nAnalyzing enrichment for {treatment} treatment...")
        logging.info(f"\nAnalyzing enrichment for {treatment} treatment...")
        
        treatment_data = mutation_data[mutation_data['Treatment'] == treatment]
        if len(treatment_data) == 0:
            print(f"No data available for {treatment}")
            logging.warning(f"No data available for {treatment}")
            continue
        
        # Use adaptive analysis
        enriched = adaptive_enrichment_analysis(treatment_data, scaffold_info)
        
        if enriched is not None and len(enriched) > 0:
            treatment_results[treatment] = enriched
    
    # Perform adaptation-specific enrichment analysis
    print("\nAnalyzing regional enrichment by adaptation type...")
    logging.info("\nAnalyzing regional enrichment by adaptation type...")
    adaptation_results = {}
    
    # Group by adaptation type
    for adaptation in mutation_data['Adaptation'].unique():
        print(f"\nAnalyzing enrichment for {adaptation} adaptation...")
        logging.info(f"\nAnalyzing enrichment for {adaptation} adaptation...")
        
        adaptation_data = mutation_data[mutation_data['Adaptation'] == adaptation]
        if len(adaptation_data) == 0:
            print(f"No data available for {adaptation} adaptation")
            logging.warning(f"No data available for {adaptation} adaptation")
            continue
        
        # Use adaptive analysis
        enriched = adaptive_enrichment_analysis(adaptation_data, scaffold_info)
        
        if enriched is not None and len(enriched) > 0:
            adaptation_results[adaptation] = enriched
    
    # Perform gene-modification enrichment analysis
    print("\nAnalyzing regional enrichment by gene modification status...")
    logging.info("\nAnalyzing regional enrichment by gene modification status...")
    gene_modification_results = {}
    
    # Split by gene modification status
    for gene_status in mutation_data['Has_Gene'].unique():
        print(f"\nAnalyzing enrichment for {'gene-modified' if gene_status == 'Yes' else 'non-modified'} strains...")
        logging.info(f"\nAnalyzing enrichment for {'gene-modified' if gene_status == 'Yes' else 'non-modified'} strains...")
        
        gene_data = mutation_data[mutation_data['Has_Gene'] == gene_status]
        if len(gene_data) == 0:
            continue
        
        # Use adaptive analysis
        enriched = adaptive_enrichment_analysis(gene_data, scaffold_info)
        
        if enriched is not None and len(enriched) > 0:
            gene_modification_results[gene_status] = enriched
    
    # Find consistently enriched regions across treatments
    print("\nFinding regions consistently enriched across treatments...")
    logging.info("\nFinding regions consistently enriched across treatments...")
    shared_regions = find_consistently_enriched_regions(treatment_results)
    
    # Analyze adaptation-specific patterns
    print("\nAnalyzing adaptation-specific enrichment patterns...")
    logging.info("\nAnalyzing adaptation-specific enrichment patterns...")
    adaptation_specific = analyze_adaptation_patterns(treatment_results)
    
    # Gene-specific enrichment analysis (new in gene_analysis version)
    gene_analysis_results = {}
    if gene_mapping_loaded and 'in_gene' in mutation_data.columns:
        print("\nPerforming gene-specific enrichment analysis...")
        logging.info("\nPerforming gene-specific enrichment analysis...")
        gene_analysis_results = analyze_gene_specific_enrichment(mutation_data, scaffold_info)
    
    # Log summary of the results
    logging.info("\nAnalysis Results Summary:")
    logging.info(f"Treatment-specific enriched regions: {sum(len(results) for results in treatment_results.values()) if treatment_results else 0}")
    logging.info(f"Adaptation-specific enriched regions: {sum(len(results) for results in adaptation_results.values()) if adaptation_results else 0}")
    logging.info(f"Gene-status-specific enriched regions: {sum(len(results) for results in gene_modification_results.values()) if gene_modification_results else 0}")
    logging.info(f"Shared enriched regions: {len(shared_regions) if shared_regions is not None else 0}")
    
    if gene_mapping_loaded and gene_analysis_results:
        logging.info(f"Gene-specific enrichment results: {len(gene_analysis_results.get('all_genes', [])) if 'all_genes' in gene_analysis_results else 0} genes analyzed")
    
    # Generate visualizations if we have results
    if treatment_results:
        print("\nGenerating visualizations...")
        logging.info("\nGenerating visualizations...")
        plot_enrichment_patterns(treatment_results, OUTPUT_DIR)
        
        # Create detailed reports
        print("\nCreating detailed enrichment reports...")
        logging.info("\nCreating detailed enrichment reports...")
        create_enrichment_reports(
            treatment_results, 
            adaptation_results, 
            gene_modification_results, 
            shared_regions, 
            adaptation_specific, 
            OUTPUT_DIR
        )
        
        # Create summary report
        print("\nCreating summary report...")
        logging.info("\nCreating summary report...")
        create_summary_report(
            treatment_results,
            adaptation_results,
            gene_modification_results,
            shared_regions,
            adaptation_specific,
            OUTPUT_DIR
        )
        
        print(f"\nRegional enrichment analysis complete! Results saved to {OUTPUT_DIR}/")
    else:
        logging.error("No enriched regions found after adaptive analysis.")
        print("\nNo enriched regions found for regional enrichment analysis. Please check regional_enrichment_debug.log for details.")
        
        # Create placeholder summary file to explain the issue
        with open(os.path.join(OUTPUT_DIR, "regional_enrichment_summary.txt"), 'w') as f:
            f.write("Regional Enrichment Analysis Summary\n")
            f.write("==================================\n\n")
            f.write("No statistically significant enriched regions were found in this analysis.\n\n")
            f.write("Potential reasons:\n")
            f.write("1. Mutations are too sparsely distributed across the genome\n")
            f.write("2. There is no significant clustering of mutations in specific regions\n")
            f.write("3. The mutation rate is too low to detect statistically significant enrichment\n\n")
            f.write("Recommendations:\n")
            f.write("1. Consider using a smaller window size for the analysis\n")
            f.write("2. Consider relaxing the statistical threshold (p-value)\n")
            f.write("3. Check if the mutation data and reference genome scaffolds match correctly\n")
            f.write("4. Consider analyzing the data at a different genomic resolution\n")
    
    # Perform gene-specific analysis if gene mapping was loaded
    if gene_mapping_loaded and 'in_gene' in mutation_data.columns:
        if gene_analysis_results and 'all_genes' in gene_analysis_results and len(gene_analysis_results['all_genes']) > 0:
            print("\nGenerating gene-specific visualizations...")
            logging.info("\nGenerating gene-specific visualizations...")
            plot_gene_enrichment(gene_analysis_results, GENE_OUTPUT_DIR)
            
            print("\nCreating gene-specific reports...")
            logging.info("\nCreating gene-specific reports...")
            create_gene_specific_reports(gene_analysis_results, GENE_OUTPUT_DIR)
            
            print("\nCreating gene-specific summary...")
            logging.info("\nCreating gene-specific summary...")
            create_gene_specific_summary(gene_analysis_results, mutation_data, GENE_OUTPUT_DIR)
            
            print(f"\nGene-specific analysis complete! Results saved to {GENE_OUTPUT_DIR}/")
        else:
            logging.warning("No gene-specific enrichment results found.")
            print("\nNo gene-specific enrichment results. Please check regional_enrichment_debug.log for details.")
            
            # Create placeholder gene summary file
            with open(os.path.join(GENE_OUTPUT_DIR, "gene_analysis_summary.txt"), 'w') as f:
                f.write("Gene-Level Mutation Analysis Summary\n")
                f.write("=================================\n\n")
                f.write("No significant gene-specific enrichment patterns were found.\n\n")
                f.write("Potential reasons:\n")
                f.write("1. Too few variants mapped to genes\n")
                f.write("2. Variants are widely distributed across genes without clear enrichment\n")
                f.write("3. The gene-specific mutation rate is too low for statistical significance\n\n")
                f.write("Recommendations:\n")
                f.write("1. Consider focusing on specific gene families or pathways\n")
                f.write("2. Consider examining gene functional categories instead of individual genes\n")
                f.write("3. Check if the gene mapping information is correct\n")
    
    print("\nAll analyses complete!")
    logging.info("All analyses complete!")

if __name__ == "__main__":
    main()