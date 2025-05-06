#!/usr/bin/env python3

'''
Variation Analysis with Gene Mapping

This script analyzes the statistical significance of treatment vs control differences
in variant patterns, with specific focus on gene-level analysis. It enhances the original
variation.py script by adding gene-specific functionality, allowing for the analysis of
variant patterns within genes, particularly focusing on genes involved in the ergosterol
biosynthesis pathway which may be under purifying selection.
'''

import os
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact, chi2_contingency, ttest_ind, mannwhitneyu, poisson
import statsmodels.stats.multitest as mt
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("analysis.log", mode='a'),
        logging.StreamHandler()
    ]
)

# Set matplotlib style
plt.style.use('ggplot')
sns.set(font_scale=1.2)
sns.set_style("whitegrid")

# Define output directories
OUTPUT_DIR = "analysis/treatment_control_analysis"
GENE_OUTPUT_DIR = "analysis/genes_of_interest/treatment_control_analysis"
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(GENE_OUTPUT_DIR, exist_ok=True)

# Define biologically correct treatment groups
TREATMENTS = ['WT-37', 'WTA', 'STC', 'CAS']

# Define treatment information for better biological context
TREATMENT_INFO = {
    'WT-37': {'description': 'Temperature-adapted wild type', 'adaptation': 'Temperature'},
    'WTA': {'description': 'Low oxygen-adapted wild type', 'adaptation': 'Low Oxygen'},
    'STC': {'description': 'STC gene with low oxygen adaptation', 'adaptation': 'Low Oxygen', 'gene': 'STC'},
    'CAS': {'description': 'CAS gene with temperature adaptation', 'adaptation': 'Temperature', 'gene': 'CAS'}
}

# Define treatment colors for consistent visualizations
TREATMENT_COLORS = {
    'WT-37': '#1b9e77',  # Temperature-adapted
    'WTA': '#d95f02',    # Low oxygen-adapted
    'STC': '#7570b3',    # STC gene + low oxygen
    'CAS': '#e7298a',    # CAS gene + temperature
    'STC-vs-STCCTRL': '#66a61e',  # STC with original control
    'CAS-vs-CASCTRL': '#e6ab02'   # CAS with original control
}

# Gene status colors
GENE_COLORS = {
    'ERG': '#2ca02c',    # Ergosterol pathway genes
    'Non-ERG': '#1f77b4', # Non-ergosterol genes
    'No Gene': '#7f7f7f'  # No gene
}

# Initialize gene data structures
GENE_DATA = {}  # Dictionary mapping gene IDs to their details
SCAFFOLD_GENES = defaultdict(list)  # Dictionary mapping scaffolds to lists of genes
GENES_OF_INTEREST = set()  # Set of gene IDs involved in the ergosterol pathway

def find_vcf_file(base_name, possible_locations):
    """Find a VCF file by trying multiple possible locations and formats."""
    for location in possible_locations:
        if os.path.exists(location.format(base_name)):
            return location.format(base_name)
    return None

# Function to load gene mapping data
def load_gene_mapping():
    """
    Load gene mapping data from reference files.
    
    This function loads gene data from the reference directory, including:
    1. Gene coordinates and information from gene_mapping.tsv
    2. Genes of interest (ergosterol pathway) from genes_of_interest_mapping.tsv
    
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
                
                # Check if chromosome_id is directly available in the gene_mapping.tsv file
                chromosome_id = None
                if 'chromosome_id' in row and not pd.isna(row['chromosome_id']):
                    chromosome_id = row['chromosome_id']
                
                # Store gene data
                gene_data = {
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
                
                # Add chromosome_id if it was found
                if chromosome_id:
                    gene_data['chromosome_id'] = chromosome_id
                
                GENE_DATA[gene_id] = gene_data
                
                # Map scaffold to genes
                if scaffold not in SCAFFOLD_GENES:
                    SCAFFOLD_GENES[scaffold] = []
                SCAFFOLD_GENES[scaffold].append(gene_id)
                
                # Also map chromosome_id to genes if available
                if chromosome_id:
                    if chromosome_id not in SCAFFOLD_GENES:
                        SCAFFOLD_GENES[chromosome_id] = []
                    SCAFFOLD_GENES[chromosome_id].append(gene_id)
            
            # Check if we need to load chromosome mapping separately
            # (only needed if the gene mapping doesn't have chromosome_id)
            if not any('chromosome_id' in gene_data for gene_data in GENE_DATA.values()):
                print("No chromosome_id found in gene mapping, attempting to load chromosome mapping files...")
                
                # Look for chromosome mapping file to map between scaffold IDs and chromosome IDs
                chromosome_mapping = {}
                mapping_paths = [
                    "reference/chromosome_mapping.tsv",
                    "reference/chromosome_mapping_reverse.tsv"
                ]
                
                for mapping_path in mapping_paths:
                    if os.path.exists(mapping_path):
                        try:
                            mapping_df = pd.read_csv(mapping_path, sep='\t')
                            print(f"Loaded chromosome mapping from {mapping_path}")
                            
                            # Check columns and create mapping
                            if mapping_path.endswith('chromosome_mapping.tsv') and 'w303_scaffold' in mapping_df.columns:
                                for _, row in mapping_df.iterrows():
                                    chromosome_id = row['chromosome_id']
                                    scaffold = row['w303_scaffold']
                                    chromosome_mapping[scaffold] = chromosome_id
                            elif mapping_path.endswith('chromosome_mapping_reverse.tsv') and 'w303_scaffold' in mapping_df.columns:
                                for _, row in mapping_df.iterrows():
                                    chromosome_id = row['chromosome_id']
                                    scaffold = row['w303_scaffold']
                                    chromosome_mapping[scaffold] = chromosome_id
                            elif 'scaffold' in mapping_df.columns and 'chromosome_id' in mapping_df.columns:
                                for _, row in mapping_df.iterrows():
                                    scaffold = row['scaffold']
                                    chromosome_id = row['chromosome_id']
                                    chromosome_mapping[scaffold] = chromosome_id
                            else:
                                print(f"Warning: Unexpected columns in {mapping_path}")
                                print(f"Available columns: {mapping_df.columns.tolist()}")
                        except Exception as e:
                            print(f"Error loading chromosome mapping: {e}")
                
                # Apply chromosome mapping to genes
                if chromosome_mapping:
                    print(f"Adding chromosome mapping to {len(chromosome_mapping)} genes...")
                    
                    for gene_id, gene_data in GENE_DATA.items():
                        scaffold = gene_data['scaffold']
                        if scaffold in chromosome_mapping:
                            chromosome_id = chromosome_mapping[scaffold]
                            GENE_DATA[gene_id]['chromosome_id'] = chromosome_id
                            
                            # Add to SCAFFOLD_GENES
                            if chromosome_id not in SCAFFOLD_GENES:
                                SCAFFOLD_GENES[chromosome_id] = []
                            if gene_id not in SCAFFOLD_GENES[chromosome_id]:
                                SCAFFOLD_GENES[chromosome_id].append(gene_id)
            
            # Debug: Print statistics about loaded data
            genes_with_chroms = sum(1 for gene_data in GENE_DATA.values() if 'chromosome_id' in gene_data)
            print(f"Loaded {len(GENE_DATA)} genes, {genes_with_chroms} have chromosome IDs")
            print(f"Mapped genes to {len(SCAFFOLD_GENES)} scaffolds/chromosomes")
            
            # Debug: Print sample of chromosome IDs and scaffolds
            chrom_ids = set(gene_data.get('chromosome_id') for gene_data in GENE_DATA.values() 
                         if 'chromosome_id' in gene_data)
            scaffold_ids = set(gene_data['scaffold'] for gene_data in GENE_DATA.values())
            
            print(f"Sample chromosome IDs: {list(chrom_ids)[:5]}")
            print(f"Sample scaffold IDs: {list(scaffold_ids)[:5]}")
            print(f"Sample SCAFFOLD_GENES keys: {list(SCAFFOLD_GENES.keys())[:10]}")
            
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

# Function to extract variants from a VCF file and map them to genes
def extract_variants_from_vcf(vcf_file):
    """
    Extract variants from a VCF file and map them to genes.
    
    Args:
        vcf_file (str): Path to the VCF file
        
    Returns:
        pandas.DataFrame: DataFrame containing variant information mapped to genes,
                         or None if extraction failed
    """
    if not vcf_file or not os.path.exists(vcf_file):
        print(f"Error: VCF file not found: {vcf_file}")
        return None
    
    try:
        # Extract variant information from VCF
        cmd = f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' {vcf_file}"
        output = subprocess.check_output(cmd, shell=True, universal_newlines=True)
        
        # Parse the output into a DataFrame
        lines = output.strip().split('\n')
        if not lines or not lines[0]:
            return pd.DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT'])
            
        data = [line.split('\t') for line in lines if line]
        variant_df = pd.DataFrame(data, columns=['CHROM', 'POS', 'REF', 'ALT'])
        
        # Convert position to numeric
        variant_df['POS'] = pd.to_numeric(variant_df['POS'])
        
        # Add variant ID column
        variant_df['Variant_ID'] = variant_df.apply(
            lambda row: f"{row['CHROM']}_{row['POS']}_{row['REF']}_{row['ALT']}", axis=1
        )
        
        # Debug: Print unique chromosome IDs in variant data
        unique_chroms = variant_df['CHROM'].unique()
        print(f"Found {len(unique_chroms)} unique chromosomes in variant data")
        print(f"First 5 chromosome IDs in variants: {list(unique_chroms)[:5]}")
        
        # Map variants to genes if gene data is available
        if GENE_DATA and SCAFFOLD_GENES:
            # Initialize gene-related columns
            variant_df['in_gene'] = False
            variant_df['gene_id'] = None
            variant_df['gene_name'] = None
            variant_df['gene_type'] = None
            variant_df['gene_product'] = None
            
            # Debug: Print scaffold information
            print(f"Total scaffold entries in SCAFFOLD_GENES: {len(SCAFFOLD_GENES)}")
            scaffold_keys = list(SCAFFOLD_GENES.keys())
            print(f"Sample scaffold keys: {scaffold_keys[:5]}")
            
            # Count matches and misses for debugging
            matches = 0
            misses = 0
            
            # Map each variant to genes
            for idx, row in variant_df.iterrows():
                chrom = row['CHROM']
                position = row['POS']
                
                # Skip if chromosome has no mapped genes
                if chrom not in SCAFFOLD_GENES:
                    misses += 1
                    
                    # Only print first few misses to avoid flooding the log
                    if misses <= 5:
                        print(f"Debug: No genes mapped for chromosome {chrom}")
                    continue
                else:
                    matches += 1
                
                # Check each gene in this chromosome
                for gene_id in SCAFFOLD_GENES[chrom]:
                    gene_data = GENE_DATA[gene_id]
                    
                    # Check if position falls within gene coordinates
                    if gene_data['start'] <= position <= gene_data['end']:
                        variant_df.at[idx, 'in_gene'] = True
                        variant_df.at[idx, 'gene_id'] = gene_id
                        
                        # Add gene name if available
                        if gene_data['erg_name']:
                            variant_df.at[idx, 'gene_name'] = gene_data['erg_name']
                        elif gene_data['sc_gene_id']:
                            variant_df.at[idx, 'gene_name'] = gene_data['sc_gene_id']
                        
                        # Set gene type based on presence in genes of interest
                        if gene_id in GENES_OF_INTEREST:
                            variant_df.at[idx, 'gene_type'] = 'ergosterol'
                        else:
                            variant_df.at[idx, 'gene_type'] = 'other'
                        
                        # Add gene product description if available
                        if gene_data['product']:
                            variant_df.at[idx, 'gene_product'] = gene_data['product']
                        
                        # Break since we found a matching gene
                        break
            
            # Log the mapping results
            in_gene_count = sum(variant_df['in_gene'])
            ergosterol_count = sum(variant_df['gene_type'] == 'ergosterol')
            
            print(f"Mapped {in_gene_count} out of {len(variant_df)} variants to genes")
            print(f"Found {ergosterol_count} variants in ergosterol pathway genes")
            print(f"Chromosome matches: {matches}, misses: {misses}")
            
            logging.info(f"Mapped {in_gene_count} out of {len(variant_df)} variants to genes")
            logging.info(f"Found {ergosterol_count} variants in ergosterol pathway genes")
        
        return variant_df
    
    except Exception as e:
        print(f"Error extracting variants from {vcf_file}: {e}")
        logging.error(f"Error extracting variants from {vcf_file}: {e}")
        return None

def load_vcf_counts(vcf_file):
    """Extract variant counts from a VCF file."""
    if not vcf_file or not os.path.exists(vcf_file):
        print(f"  Error: VCF file not found: {vcf_file}")
        return 0
    
    try:
        # Count total variants
        cmd = f"bcftools view -H {vcf_file} | wc -l"
        variant_count = int(subprocess.check_output(cmd, shell=True))
        return variant_count
    except Exception as e:
        print(f"  Error processing {vcf_file}: {e}")
        return 0

def analyze_gene_specific_patterns(treatment_variants, control_variants, treatment_name):
    """
    Analyze gene-specific patterns between treatment and control.
    
    This function analyzes the distribution of variants within genes, especially
    focusing on genes involved in the ergosterol pathway. It detects patterns of
    enrichment or depletion (purifying selection) within genes.
    
    Args:
        treatment_variants (pandas.DataFrame): DataFrame with treatment variants
        control_variants (pandas.DataFrame): DataFrame with control variants
        treatment_name (str): Name of the treatment
        
    Returns:
        dict: Dictionary containing gene-specific analysis results
    """
    results = {}
    
    # Skip if either DataFrame is empty
    if treatment_variants is None or control_variants is None:
        return results
    if len(treatment_variants) == 0 or len(control_variants) == 0:
        return results
    
    # Calculate gene-level statistics
    # Count variants in each gene for treatment and control
    treatment_gene_counts = defaultdict(int)
    control_gene_counts = defaultdict(int)
    gene_lengths = {}
    
    # Calculate total counts for each gene status
    treatment_erg_count = 0
    treatment_non_erg_count = 0
    treatment_non_gene_count = 0
    control_erg_count = 0
    control_non_erg_count = 0
    control_non_gene_count = 0
    
    # Process treatment variants
    for _, row in treatment_variants.iterrows():
        if row['in_gene']:
            gene_id = row['gene_id']
            treatment_gene_counts[gene_id] += 1
            
            # Store gene length if not already stored
            if gene_id not in gene_lengths and gene_id in GENE_DATA:
                gene_lengths[gene_id] = GENE_DATA[gene_id]['end'] - GENE_DATA[gene_id]['start'] + 1
            
            # Count by gene type
            if row['gene_type'] == 'ergosterol':
                treatment_erg_count += 1
            else:
                treatment_non_erg_count += 1
        else:
            treatment_non_gene_count += 1
    
    # Process control variants
    for _, row in control_variants.iterrows():
        if row['in_gene']:
            gene_id = row['gene_id']
            control_gene_counts[gene_id] += 1
            
            # Store gene length if not already stored
            if gene_id not in gene_lengths and gene_id in GENE_DATA:
                gene_lengths[gene_id] = GENE_DATA[gene_id]['end'] - GENE_DATA[gene_id]['start'] + 1
            
            # Count by gene type
            if row['gene_type'] == 'ergosterol':
                control_erg_count += 1
            else:
                control_non_erg_count += 1
        else:
            control_non_gene_count += 1
    
    # Store the gene status counts
    results['gene_status_counts'] = {
        'treatment': {
            'ergosterol': treatment_erg_count,
            'non_ergosterol': treatment_non_erg_count,
            'non_gene': treatment_non_gene_count
        },
        'control': {
            'ergosterol': control_erg_count,
            'non_ergosterol': control_non_erg_count,
            'non_gene': control_non_gene_count
        }
    }
    
    # Calculate the fold change for each gene type
    erg_fold_change = (treatment_erg_count / max(1, control_erg_count))
    non_erg_fold_change = (treatment_non_erg_count / max(1, control_non_erg_count))
    non_gene_fold_change = (treatment_non_gene_count / max(1, control_non_gene_count))
    
    results['gene_status_fold_change'] = {
        'ergosterol': erg_fold_change,
        'non_ergosterol': non_erg_fold_change,
        'non_gene': non_gene_fold_change
    }
    
    # Calculate gene-specific fold changes and statistical significance
    gene_results = []
    
    # Get all genes with variants in either treatment or control
    all_genes = set(treatment_gene_counts.keys()) | set(control_gene_counts.keys())
    
    # Calculate expected variants per position
    treatment_total = len(treatment_variants)
    control_total = len(control_variants)
    genome_size = 12000000  # Approximate S. cerevisiae genome size
    
    treatment_rate = treatment_total / genome_size
    control_rate = control_total / genome_size
    
    # Calculate enrichment/depletion for each gene
    for gene_id in all_genes:
        # Get gene information
        gene_info = GENE_DATA.get(gene_id, {})
        gene_name = gene_info.get('erg_name') or gene_info.get('sc_gene_id') or gene_id
        gene_type = 'ergosterol' if gene_id in GENES_OF_INTEREST else 'other'
        gene_length = gene_lengths.get(gene_id, 0)
        
        # Skip genes with unknown length
        if gene_length == 0:
            continue
        
        # Get observed variant counts
        treatment_count = treatment_gene_counts.get(gene_id, 0)
        control_count = control_gene_counts.get(gene_id, 0)
        
        # Calculate expected counts based on genome-wide rates
        expected_treatment = treatment_rate * gene_length
        expected_control = control_rate * gene_length
        
        # Calculate enrichment (> 1) or depletion (< 1) fold change
        treatment_fold_change = treatment_count / max(1, expected_treatment)
        control_fold_change = control_count / max(1, expected_control)
        
        # Calculate differential fold change (treatment vs control)
        differential_fold_change = treatment_fold_change / max(0.001, control_fold_change)
        
        # Calculate statistical significance using Poisson distribution
        # For treatment count compared to expected
        if treatment_count > expected_treatment:
            # Testing for enrichment
            treatment_p_value = 1 - poisson.cdf(treatment_count - 1, expected_treatment)
            treatment_direction = 'enriched'
        else:
            # Testing for depletion (purifying selection)
            treatment_p_value = poisson.cdf(treatment_count, expected_treatment)
            treatment_direction = 'depleted'
        
        # For control count compared to expected
        if control_count > expected_control:
            control_p_value = 1 - poisson.cdf(control_count - 1, expected_control)
            control_direction = 'enriched'
        else:
            control_p_value = poisson.cdf(control_count, expected_control)
            control_direction = 'depleted'
        
        # Store gene-specific results
        gene_results.append({
            'gene_id': gene_id,
            'gene_name': gene_name,
            'gene_type': gene_type,
            'gene_length': gene_length,
            'treatment_count': treatment_count,
            'control_count': control_count,
            'expected_treatment': expected_treatment,
            'expected_control': expected_control,
            'treatment_fold_change': treatment_fold_change,
            'control_fold_change': control_fold_change,
            'differential_fold_change': differential_fold_change,
            'treatment_p_value': treatment_p_value,
            'treatment_direction': treatment_direction,
            'control_p_value': control_p_value,
            'control_direction': control_direction
        })
    
    # Convert to DataFrame
    gene_df = pd.DataFrame(gene_results)
    
    # Apply multiple testing correction if there are results
    if len(gene_df) > 0:
        # Correct p-values for treatment
        gene_df['treatment_q_value'] = mt.multipletests(gene_df['treatment_p_value'], method='fdr_bh')[1]
        
        # Correct p-values for control
        gene_df['control_q_value'] = mt.multipletests(gene_df['control_p_value'], method='fdr_bh')[1]
        
        # Add log2 fold change for better visualization
        gene_df['log2_treatment_fold_change'] = np.log2(gene_df['treatment_fold_change'])
        gene_df['log2_control_fold_change'] = np.log2(gene_df['control_fold_change'])
        gene_df['log2_differential_fold_change'] = np.log2(gene_df['differential_fold_change'])
        
        # Sort by significance and fold change
        gene_df = gene_df.sort_values(['gene_type', 'treatment_p_value'])
    
    # Store gene DataFrame in results
    results['gene_analysis'] = gene_df
    
    # Compare ergosterol vs. non-ergosterol genes
    if len(gene_df) > 0:
        # Get ergosterol and non-ergosterol genes
        erg_genes = gene_df[gene_df['gene_type'] == 'ergosterol']
        non_erg_genes = gene_df[gene_df['gene_type'] == 'other']
        
        # Compare fold changes if we have enough data
        if len(erg_genes) > 0 and len(non_erg_genes) > 0:
            # Parametric test (t-test)
            t_stat, t_pvalue = ttest_ind(
                erg_genes['treatment_fold_change'], 
                non_erg_genes['treatment_fold_change'],
                equal_var=False
            )
            
            # Non-parametric test (Mann-Whitney U)
            u_stat, u_pvalue = mannwhitneyu(
                erg_genes['treatment_fold_change'], 
                non_erg_genes['treatment_fold_change'],
                alternative='two-sided'
            )
            
            results['erg_vs_non_erg'] = {
                'mean_erg_fold_change': erg_genes['treatment_fold_change'].mean(),
                'mean_non_erg_fold_change': non_erg_genes['treatment_fold_change'].mean(),
                'median_erg_fold_change': erg_genes['treatment_fold_change'].median(),
                'median_non_erg_fold_change': non_erg_genes['treatment_fold_change'].median(),
                't_statistic': float(t_stat),
                't_pvalue': float(t_pvalue),
                'u_statistic': float(u_stat),
                'u_pvalue': float(u_pvalue)
            }
    
    return results

def create_gene_specific_visualizations(all_gene_results, results_df):
    """
    Create gene-specific visualizations from the analysis results.
    
    Args:
        all_gene_results (dict): Dictionary mapping treatment names to gene analysis results
        results_df (pandas.DataFrame): DataFrame with treatment vs control statistics
    """
    # Create output directory
    os.makedirs(GENE_OUTPUT_DIR, exist_ok=True)
    
    # Prepare combined data
    gene_status_data = []
    erg_gene_data = []
    
    # Process each treatment
    for treatment_name, gene_results in all_gene_results.items():
        # Skip if no gene status counts
        if 'gene_status_counts' not in gene_results:
            continue
        
        # Extract gene status counts
        status_counts = gene_results['gene_status_counts']['treatment']
        
        # Add to combined data
        gene_status_data.append({
            'Treatment': treatment_name,
            'Gene Status': 'ERG',
            'Variants': status_counts['ergosterol']
        })
        gene_status_data.append({
            'Treatment': treatment_name,
            'Gene Status': 'Non-ERG',
            'Variants': status_counts['non_ergosterol']
        })
        gene_status_data.append({
            'Treatment': treatment_name,
            'Gene Status': 'No Gene',
            'Variants': status_counts['non_gene']
        })
        
        # Add ergosterol gene details if gene analysis is available
        if 'gene_analysis' in gene_results and len(gene_results['gene_analysis']) > 0:
            gene_df = gene_results['gene_analysis']
            
            # Get ergosterol genes
            erg_genes = gene_df[gene_df['gene_type'] == 'ergosterol']
            
            # Add each gene's data
            for _, row in erg_genes.iterrows():
                erg_gene_data.append({
                    'Treatment': treatment_name,
                    'Gene': row['gene_name'],
                    'Variants': row['treatment_count'],
                    'Fold Change': row['treatment_fold_change']
                })
    
    # Convert to DataFrames
    status_df = pd.DataFrame(gene_status_data)
    erg_df = pd.DataFrame(erg_gene_data)
    
    # 1. Gene Status Distribution
    if len(status_df) > 0:
        plt.figure(figsize=(14, 8))
        
        # Create a grouped bar chart
        ax = sns.barplot(
            data=status_df,
            x='Treatment',
            y='Variants',
            hue='Gene Status',
            palette=GENE_COLORS
        )
        
        # Add value labels
        for container in ax.containers:
            ax.bar_label(container, fmt='%d')
        
        # Customize plot
        plt.title('Distribution of Variants by Gene Status Across Treatments', fontsize=14)
        plt.xlabel('Treatment')
        plt.ylabel('Number of Variants')
        plt.legend(title='Gene Status')
        
        # Save figure
        plt.tight_layout()
        plt.savefig(f'{GENE_OUTPUT_DIR}/gene_status_distribution.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # 2. Stacked Gene Status Proportion
    if len(status_df) > 0:
        plt.figure(figsize=(14, 8))
        
        # Pivot data for stacked plot
        pivot_df = status_df.pivot(index='Treatment', columns='Gene Status', values='Variants').fillna(0)
        
        # Calculate proportions
        pivot_df = pivot_df.div(pivot_df.sum(axis=1), axis=0) * 100
        
        # Plot stacked bars
        pivot_df.plot(
            kind='bar', 
            stacked=True,
            color=[GENE_COLORS[status] for status in pivot_df.columns],
            figsize=(14, 8)
        )
        
        # Customize plot
        plt.title('Proportion of Variants by Gene Status Across Treatments', fontsize=14)
        plt.xlabel('Treatment')
        plt.ylabel('Percentage of Variants')
        plt.legend(title='Gene Status')
        
        # Format percentages
        for container in plt.gca().containers:
            plt.gca().bar_label(container, label_type='center', fmt='%.1f%%')
        
        # Save figure
        plt.tight_layout()
        plt.savefig(f'{GENE_OUTPUT_DIR}/gene_status_proportion.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # 3. Ergosterol Gene Distribution (if we have data)
    if len(erg_df) > 0:
        plt.figure(figsize=(14, 10))
        
        # Create a grouped bar chart
        ax = sns.barplot(
            data=erg_df,
            x='Gene',
            y='Variants',
            hue='Treatment',
            palette=TREATMENT_COLORS
        )
        
        # Customize plot
        plt.title('Ergosterol Pathway Gene Variant Distribution by Treatment', fontsize=14)
        plt.xlabel('Gene')
        plt.ylabel('Number of Variants')
        plt.xticks(rotation=45, ha='right')
        
        # Add value labels
        for container in ax.containers:
            ax.bar_label(container, fmt='%d')
        
        # Save figure
        plt.tight_layout()
        plt.savefig(f'{GENE_OUTPUT_DIR}/erg_gene_distribution.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # 4. Fold Change Visualization
    fold_change_data = []
    
    # Process each treatment
    for treatment_name, gene_results in all_gene_results.items():
        # Skip if no gene analysis
        if 'gene_analysis' not in gene_results or len(gene_results['gene_analysis']) == 0:
            continue
        
        # Get treatment row from results_df
        treatment_row = results_df[results_df['Treatment'] == treatment_name]
        if len(treatment_row) == 0:
            continue
            
        treatment_fold_change = treatment_row['Fold_Change'].values[0]
        
        # Extract fold changes by gene type
        gene_df = gene_results['gene_analysis']
        
        # Get mean fold change for ergosterol and non-ergosterol genes
        erg_fold = gene_df[gene_df['gene_type'] == 'ergosterol']['treatment_fold_change'].mean()
        non_erg_fold = gene_df[gene_df['gene_type'] == 'other']['treatment_fold_change'].mean()
        
        # Add to data
        fold_change_data.append({
            'Treatment': treatment_name,
            'Category': 'Overall',
            'Fold Change': treatment_fold_change
        })
        
        if not np.isnan(erg_fold):
            fold_change_data.append({
                'Treatment': treatment_name,
                'Category': 'ERG Genes',
                'Fold Change': erg_fold
            })
        
        if not np.isnan(non_erg_fold):
            fold_change_data.append({
                'Treatment': treatment_name,
                'Category': 'Non-ERG Genes',
                'Fold Change': non_erg_fold
            })
    
    # Create fold change comparison chart
    if len(fold_change_data) > 0:
        fold_df = pd.DataFrame(fold_change_data)
        
        plt.figure(figsize=(14, 8))
        
        # Create a grouped bar chart
        ax = sns.barplot(
            data=fold_df,
            x='Treatment',
            y='Fold Change',
            hue='Category',
            palette={
                'Overall': '#ff7f0e',
                'ERG Genes': '#2ca02c',
                'Non-ERG Genes': '#1f77b4'
            }
        )
        
        # Add horizontal line at fold change = 1
        plt.axhline(y=1, color='black', linestyle='--', alpha=0.5)
        
        # Customize plot
        plt.title('Fold Change Comparison by Gene Category', fontsize=14)
        plt.xlabel('Treatment')
        plt.ylabel('Fold Change (Treatment/Control)')
        plt.legend(title='Category')
        
        # Add value labels
        for container in ax.containers:
            ax.bar_label(container, fmt='%.2f')
        
        # Save figure
        plt.tight_layout()
        plt.savefig(f'{GENE_OUTPUT_DIR}/fold_change_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # 5. Purifying Selection Visualization
    # For treatments with gene analysis, create volcano plots
    for treatment_name, gene_results in all_gene_results.items():
        # Skip if no gene analysis
        if 'gene_analysis' not in gene_results or len(gene_results['gene_analysis']) == 0:
            continue
        
        gene_df = gene_results['gene_analysis']
        
        # Create volcano plot
        plt.figure(figsize=(12, 10))
        
        # Plot all genes
        plt.scatter(
            gene_df[gene_df['gene_type'] == 'other']['log2_treatment_fold_change'],
            -np.log10(gene_df[gene_df['gene_type'] == 'other']['treatment_p_value']),
            alpha=0.7,
            s=50,
            color=GENE_COLORS['Non-ERG'],
            label='Non-ERG Genes'
        )
        
        # Highlight ergosterol genes
        plt.scatter(
            gene_df[gene_df['gene_type'] == 'ergosterol']['log2_treatment_fold_change'],
            -np.log10(gene_df[gene_df['gene_type'] == 'ergosterol']['treatment_p_value']),
            alpha=0.9,
            s=100,
            color=GENE_COLORS['ERG'],
            label='ERG Genes'
        )
        
        # Add gene labels for significant ergosterol genes
        sig_erg_genes = gene_df[
            (gene_df['gene_type'] == 'ergosterol') & 
            (gene_df['treatment_q_value'] < 0.1)
        ]
        
        for _, row in sig_erg_genes.iterrows():
            plt.text(
                row['log2_treatment_fold_change'],
                -np.log10(row['treatment_p_value']),
                row['gene_name'],
                fontsize=10,
                ha='center',
                va='bottom'
            )
        
        # Add reference lines
        plt.axvline(x=0, color='gray', linestyle='--', alpha=0.5)  # Vertical line at log2FC = 0
        plt.axhline(y=-np.log10(0.05), color='gray', linestyle='--', alpha=0.5)  # Horizontal line at p = 0.05
        
        # Add regions of interest
        plt.text(2, 0.5, 'Enriched\n(more variants)', ha='center', fontsize=12, color='darkred')
        plt.text(-2, 0.5, 'Depleted\n(purifying selection)', ha='center', fontsize=12, color='darkblue')
        
        # Customize plot
        plt.title(f'Gene Variant Enrichment/Depletion in {treatment_name}', fontsize=14)
        plt.xlabel('Log2 Fold Change (observed/expected variants)', fontsize=12)
        plt.ylabel('-Log10 P-value', fontsize=12)
        plt.legend(title='Gene Type')
        plt.grid(True, alpha=0.3)
        
        # Set axis limits
        x_min = min(-3, gene_df['log2_treatment_fold_change'].min() - 0.5)
        x_max = max(3, gene_df['log2_treatment_fold_change'].max() + 0.5)
        plt.xlim(x_min, x_max)
        
        # Save figure
        plt.tight_layout()
        plt.savefig(f'{GENE_OUTPUT_DIR}/{treatment_name}_purifying_selection.png', dpi=300, bbox_inches='tight')
        plt.close()

def create_gene_specific_report(all_gene_results, results_df):
    """
    Create a gene-specific report with analysis results.
    
    Args:
        all_gene_results (dict): Dictionary mapping treatment names to gene analysis results
        results_df (pandas.DataFrame): DataFrame with treatment vs control statistics
    """
    # Create output directory
    os.makedirs(GENE_OUTPUT_DIR, exist_ok=True)
    
    # Prepare report
    report_lines = [
        "# Gene-Specific Analysis of Treatment vs Control Data",
        "\n## Overview\n",
        "This report provides an analysis of gene-specific variation patterns across treatments,",
        "with a focus on detecting purifying selection in ergosterol pathway genes.\n"
    ]
    
    # Add summary table
    report_lines.extend([
        "## Treatment Summary\n",
        "| Treatment | Total Variants | ERG Genes | Non-ERG Genes | Non-Genic | Overall Fold Change |",
        "| --- | --- | --- | --- | --- | --- |"
    ])
    
    # Add row for each treatment
    for treatment_name, gene_results in all_gene_results.items():
        # Get treatment row from results_df
        treatment_row = results_df[results_df['Treatment'] == treatment_name]
        if len(treatment_row) == 0 or 'gene_status_counts' not in gene_results:
            continue
            
        # Get fold change
        fold_change = treatment_row['Fold_Change'].values[0]
        
        # Get gene status counts
        status_counts = gene_results['gene_status_counts']['treatment']
        
        # Add to table
        report_lines.append(
            f"| {treatment_name} | {sum(status_counts.values())} | "
            f"{status_counts['ergosterol']} | {status_counts['non_ergosterol']} | "
            f"{status_counts['non_gene']} | {fold_change:.2f} |"
        )
    
    # Add detailed purifying selection analysis
    report_lines.extend([
        "\n## Purifying Selection Analysis\n",
        "Purifying selection can be detected by comparing the observed number of variants in a gene",
        "to the expected number based on gene length and genome-wide mutation rate.",
        "A fold change < 1 indicates fewer variants than expected (purifying selection),",
        "while a fold change > 1 indicates more variants than expected (positive selection).\n"
    ])
    
    # Add treatment-specific analysis
    for treatment_name, gene_results in all_gene_results.items():
        # Skip if no gene analysis or erg vs non-erg comparison
        if 'gene_analysis' not in gene_results or 'erg_vs_non_erg' not in gene_results:
            continue
        
        # Get statistics
        stats = gene_results['erg_vs_non_erg']
        
        # Add treatment section
        report_lines.extend([
            f"### {treatment_name} Treatment\n",
            f"- Mean fold change in ERG genes: {stats['mean_erg_fold_change']:.2f}",
            f"- Mean fold change in non-ERG genes: {stats['mean_non_erg_fold_change']:.2f}",
            f"- Median fold change in ERG genes: {stats['median_erg_fold_change']:.2f}",
            f"- Median fold change in non-ERG genes: {stats['median_non_erg_fold_change']:.2f}",
            f"- t-test p-value: {stats['t_pvalue']:.4f}",
            f"- Mann-Whitney U p-value: {stats['u_pvalue']:.4f}",
            ""
        ])
        
        # Interpret results
        erg_ratio = stats['mean_erg_fold_change'] / stats['mean_non_erg_fold_change']
        if erg_ratio < 0.75:
            report_lines.append("**Interpretation:** Strong evidence of purifying selection in ergosterol pathway genes (significantly fewer mutations than in other genes).")
        elif erg_ratio < 0.9:
            report_lines.append("**Interpretation:** Moderate evidence of purifying selection in ergosterol pathway genes (fewer mutations than in other genes).")
        elif erg_ratio < 1.1:
            report_lines.append("**Interpretation:** No evidence of differential selection between ergosterol and other genes.")
        else:
            report_lines.append("**Interpretation:** Evidence of increased mutation rate in ergosterol pathway genes compared to other genes.")
        
        report_lines.append("")
        
        # Add significant gene table
        gene_df = gene_results['gene_analysis']
        
        # Filter for significant genes (q-value < 0.1)
        sig_genes = gene_df[gene_df['treatment_q_value'] < 0.1].sort_values('treatment_p_value')
        
        if len(sig_genes) > 0:
            report_lines.extend([
                f"#### Significantly Enriched/Depleted Genes in {treatment_name}\n",
                "| Gene | Type | Observed | Expected | Fold Change | Direction | Q-value |",
                "| --- | --- | --- | --- | --- | --- | --- |"
            ])
            
            # Add each significant gene
            for _, row in sig_genes.iterrows():
                gene_type = 'ERG' if row['gene_type'] == 'ergosterol' else 'Non-ERG'
                
                report_lines.append(
                    f"| {row['gene_name']} | {gene_type} | {row['treatment_count']} | "
                    f"{row['expected_treatment']:.1f} | {row['treatment_fold_change']:.2f} | "
                    f"{row['treatment_direction']} | {row['treatment_q_value']:.3e} |"
                )
            
            report_lines.append("")
    
    # Write report to file
    with open(f'{GENE_OUTPUT_DIR}/gene_specific_report.md', 'w') as f:
        f.write('\n'.join(report_lines))
    
    print(f"Gene-specific report written to {GENE_OUTPUT_DIR}/gene_specific_report.md")

def analyze_treatment_control_differences():
    """Analyze statistical significance of treatment vs control differences with gene-specific analysis."""
    
    # Load gene mapping data if available
    gene_mapping_loaded = load_gene_mapping()
    
    # Define treatments and their controls
    # Updated to reflect the correct biological grouping
    treatments = {
        # Primary comparisons with WT-CTRL as baseline
        'WT-37': {'treatment': 'WT-37-55', 'control': 'WT-CTRL', 'description': 'Temperature-adapted wild type'},
        'WTA': {'treatment': 'WTA-55', 'control': 'WT-CTRL', 'description': 'Low oxygen-adapted wild type'},
        'STC': {'treatment': 'STC-55', 'control': 'WT-CTRL', 'description': 'STC gene with low oxygen adaptation'},
        'CAS': {'treatment': 'CAS-55', 'control': 'WT-CTRL', 'description': 'CAS gene with temperature adaptation'},
        
        # Original control comparisons (preserved for completeness)
        'STC-vs-STCCTRL': {'treatment': 'STC-55', 'control': 'STC-CTRL', 'description': 'STC vs STC control'},
        'CAS-vs-CASCTRL': {'treatment': 'CAS-55', 'control': 'CAS-CTRL', 'description': 'CAS vs CAS control'}
    }
    
    # Define possible treatment VCF locations to check
    treatment_locations = [
        "results/merged/analysis/{}/highconf.vcf.gz",
        "results/merged/analysis/{}_highconf.vcf.gz",
        "results/merged/analysis/{}/specific.vcf.gz",
        "results/merged/analysis/{}_specific.vcf.gz"
    ]
    
    # Define possible control VCF locations to check
    control_locations = [
        # Try individual directory with different extensions
        "results/vcf/individual/{}.vcf.gz",
        "results/vcf/individual/{}.vcf",
        "results/vcf/individual/{}.norm.vcf",
        # Try merged filtered directory
        "results/vcf/merged/filtered/{}.filtered.vcf.gz",
        # Try merged fixed directory
        "results/vcf/merged/fixed/{}.fixed.vcf.gz"
    ]
    
    results = []
    all_gene_results = {}  # Store gene-specific results for each treatment
    
    # Process each treatment
    for treatment_name, info in treatments.items():
        print(f"Treatment: {treatment_name}")
        
        # Find treatment VCF
        if treatment_name in ['WT-37', 'WTA', 'STC', 'CAS']:
            # For main treatment groups
            treatment_vcf = find_vcf_file(treatment_name, treatment_locations)
        else:
            # For the original control comparisons
            base_treatment = treatment_name.split('-vs-')[0]
            treatment_vcf = find_vcf_file(base_treatment, treatment_locations)
        
        # Find control VCF
        control_vcf = find_vcf_file(info['control'], control_locations)
        
        print(f"  Treatment VCF: {treatment_vcf}")
        print(f"  Control VCF: {control_vcf}")
        
        # Get variant counts
        treatment_count = load_vcf_counts(treatment_vcf)
        control_count = load_vcf_counts(control_vcf)
        
        print(f"  Treatment count: {treatment_count}")
        print(f"  Control count: {control_count}")
        
        # Create contingency table
        # Using genome size as background (total possible positions)
        genome_size = 12000000  # Approximate S. cerevisiae genome size
        
        contingency_table = np.array([
            [treatment_count, genome_size - treatment_count],
            [control_count, genome_size - control_count]
        ])
        
        # Perform Fisher's exact test (stable and doesn't require chi-square assumptions)
        odds_ratio, p_value = fisher_exact(contingency_table)
        
        # Store results
        results.append({
            'Treatment': treatment_name,
            'Description': info['description'],
            'Treatment_Variants': treatment_count,
            'Control': info['control'],
            'Control_Variants': control_count,
            'Odds_Ratio': odds_ratio,
            'P_Value': p_value
        })
        
        # Perform gene-specific analysis if gene mapping was loaded
        if gene_mapping_loaded:
            print(f"  Performing gene-specific analysis for {treatment_name}...")
            
            # Extract variants and map to genes
            treatment_variants = extract_variants_from_vcf(treatment_vcf)
            control_variants = extract_variants_from_vcf(control_vcf)
            
            # Add treatment information
            if treatment_variants is not None:
                treatment_variants['Treatment'] = treatment_name
            
            if control_variants is not None:
                control_variants['Treatment'] = info['control']
            
            # Analyze gene-specific patterns
            gene_results = analyze_gene_specific_patterns(
                treatment_variants, 
                control_variants, 
                treatment_name
            )
            
            # Store gene-specific results
            all_gene_results[treatment_name] = gene_results
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results)
    
    # Perform multiple testing correction
    results_df['Q_Value'] = mt.multipletests(results_df['P_Value'], method='fdr_bh')[1]
    
    # Add fold change
    results_df['Fold_Change'] = results_df['Treatment_Variants'] / results_df['Control_Variants'].replace(0, 1)
    
    # Sort by significance
    results_df = results_df.sort_values('P_Value')
    
    # Create output directory if it doesn't exist
    os.makedirs('analysis/treatment_control_analysis', exist_ok=True)
    
    # Save results
    results_df.to_csv('analysis/treatment_control_analysis/treatment_vs_control_statistics.csv', index=False)
    
    # Create a detailed report
    with open('analysis/treatment_control_analysis/statistical_analysis_report.txt', 'w') as f:
        f.write("Treatment vs Control Statistical Analysis\n")
        f.write("=====================================\n\n")
        
        for _, row in results_df.iterrows():
            f.write(f"{row['Treatment']} Treatment Analysis:\n")
            f.write("-" * (len(row['Treatment']) + 20) + "\n")
            f.write(f"Description: {row['Description']}\n")
            f.write(f"Control: {row['Control']}\n")
            f.write(f"Treatment variants: {row['Treatment_Variants']}\n")
            f.write(f"Control variants: {row['Control_Variants']}\n")
            f.write(f"Fold change: {row['Fold_Change']:.2f}\n")
            f.write(f"Odds ratio: {row['Odds_Ratio']:.2f}\n")
            f.write(f"P-value: {row['P_Value']:.2e}\n")
            f.write(f"Q-value (FDR-corrected): {row['Q_Value']:.2e}\n")
            f.write(f"Statistical significance: {'***' if row['Q_Value'] < 0.001 else '**' if row['Q_Value'] < 0.01 else '*' if row['Q_Value'] < 0.05 else 'ns'}\n\n")
    
    # Create gene-specific visualizations and report if gene mapping was loaded
    if gene_mapping_loaded and all_gene_results:
        create_gene_specific_visualizations(all_gene_results, results_df)
        create_gene_specific_report(all_gene_results, results_df)
        print("Gene-specific analysis complete! Results saved to genes_of_interest/treatment_control_analysis/")
    
    print("Analysis complete! Results saved to treatment_control_analysis/")
    return results_df

if __name__ == "__main__":
    # Load gene mapping data for gene-specific analysis
    print("Yeast MSA Variation Analysis with Gene Mapping")
    print("=============================================")
    
    # Run analysis
    results = analyze_treatment_control_differences()
    
    # Print summary
    print("\nSummary of Results:")
    print(results[['Treatment', 'Control', 'Treatment_Variants', 'Control_Variants', 'Fold_Change', 'Q_Value']])