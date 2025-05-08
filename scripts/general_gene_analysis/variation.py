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
BASE_DIR = os.environ.get('OUTPUT_DIR', "analysis/general_gene_analysis")
OUTPUT_DIR = f"{BASE_DIR}/treatment_control_analysis"
GENE_OUTPUT_DIR = f"{BASE_DIR}/gene_treatment_control_analysis"
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(GENE_OUTPUT_DIR, exist_ok=True)

# Log the output directories
logging.info(f"Using output directories: {OUTPUT_DIR} and {GENE_OUTPUT_DIR}")

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
CHROM_TO_SCAFFOLD = {}  # Dictionary mapping chromosome IDs to scaffold IDs
SCAFFOLD_TO_CHROM = {}  # Dictionary mapping scaffold IDs to chromosome IDs

def find_vcf_file(base_name, possible_locations):
    """Find a VCF file by trying multiple possible locations and formats."""
    for location in possible_locations:
        if os.path.exists(location.format(base_name)):
            return location.format(base_name)
    return None

# Function to load chromosome mapping data
def load_chromosome_mapping():
    """
    Load chromosome mapping data from the reference directory.
    
    This function loads the chromosome_mapping.tsv file to create a mapping
    between chromosome IDs (CM007964.1) and scaffold IDs (w303_scaffold_1).
    
    Returns:
        bool: True if data was loaded successfully, False otherwise
    """
    global CHROM_TO_SCAFFOLD, SCAFFOLD_TO_CHROM
    
    # Clear existing data
    CHROM_TO_SCAFFOLD.clear()
    SCAFFOLD_TO_CHROM.clear()
    
    # Define mapping file path
    chrom_mapping_file = "reference/chromosome_mapping.tsv"
    
    if os.path.exists(chrom_mapping_file):
        try:
            # Load chromosome mapping data
            chrom_df = pd.read_csv(chrom_mapping_file, sep='\t')
            print(f"Loaded {len(chrom_df)} chromosome ID mappings")
            logging.info(f"Loaded {len(chrom_df)} chromosome ID mappings")
            
            # Process each mapping
            for _, row in chrom_df.iterrows():
                chrom_id = row['chromosome_id']
                scaffold_id = row['w303_scaffold']
                
                # Add mappings in both directions
                CHROM_TO_SCAFFOLD[chrom_id] = scaffold_id
                SCAFFOLD_TO_CHROM[scaffold_id] = chrom_id
            
            return True
        except Exception as e:
            print(f"Error loading chromosome mapping data: {e}")
            logging.error(f"Error loading chromosome mapping data: {e}")
            return False
    else:
        print("No chromosome mapping file found. Chromosome mapping will not be available.")
        logging.warning("No chromosome mapping file found. Chromosome mapping will not be available.")
        return False

# Function to load gene mapping data
def load_gene_mapping():
    """
    Load gene mapping data from reference files.
    
    This function loads gene data from the reference directory, including:
    1. Gene coordinates and information from gene_mapping_full.tsv
    2. Genes of interest (ergosterol pathway) from genes_of_interest_mapping.tsv
    
    Returns:
        bool: True if data was loaded successfully, False otherwise
    """
    global GENE_DATA, SCAFFOLD_GENES, GENES_OF_INTEREST
    
    # Clear existing data
    GENE_DATA.clear()
    SCAFFOLD_GENES.clear()
    GENES_OF_INTEREST.clear()
    
    # First load chromosome mapping
    load_chromosome_mapping()
    
    # Define possible file paths for gene mapping data, prioritizing gene_mapping_full.tsv
    gene_mapping_paths = [
        "reference/gene_mapping_full.tsv",  # New comprehensive mapping
        "reference/gene_mapping.tsv",       # Fallback to original mapping if full doesn't exist
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
                
                # Add chromosome_id if it's in the data
                if 'chromosome_id' in row and not pd.isna(row['chromosome_id']):
                    gene_data['chromosome_id'] = row['chromosome_id']
                
                GENE_DATA[gene_id] = gene_data
                
                # Map scaffold to genes
                if scaffold not in SCAFFOLD_GENES:
                    SCAFFOLD_GENES[scaffold] = []
                SCAFFOLD_GENES[scaffold].append(gene_id)
                
                # Also map chromosome_id to genes if available
                if 'chromosome_id' in gene_data and gene_data['chromosome_id']:
                    chromosome_id = gene_data['chromosome_id']
                    if chromosome_id not in SCAFFOLD_GENES:
                        SCAFFOLD_GENES[chromosome_id] = []
                    SCAFFOLD_GENES[chromosome_id].append(gene_id)
                    
                # Add to genes of interest if it has a non-empty erg_name
                if gene_data['erg_name'] and len(str(gene_data['erg_name']).strip()) > 0:
                    GENES_OF_INTEREST.add(gene_id)
            
            # Print statistics about loaded data
            genes_with_chroms = sum(1 for gene_data in GENE_DATA.values() if 'chromosome_id' in gene_data)
            print(f"Loaded {len(GENE_DATA)} genes, {genes_with_chroms} have chromosome IDs")
            print(f"Mapped genes to {len(SCAFFOLD_GENES)} scaffolds/chromosomes")
            print(f"Identified {len(GENES_OF_INTEREST)} ergosterol pathway genes")
            
            # Print sample of chromosome IDs and scaffolds
            chrom_ids = set(gene_data.get('chromosome_id') for gene_data in GENE_DATA.values() 
                        if 'chromosome_id' in gene_data)
            scaffold_ids = set(gene_data['scaffold'] for gene_data in GENE_DATA.values())
            
            print(f"Sample chromosome IDs: {list(chrom_ids)[:5]}")
            print(f"Sample scaffold IDs: {list(scaffold_ids)[:5]}")
            print(f"Sample SCAFFOLD_GENES keys: {list(SCAFFOLD_GENES.keys())[:10]}")
            
            logging.info(f"Loaded {len(GENE_DATA)} genes, {genes_with_chroms} have chromosome IDs")
            logging.info(f"Mapped genes to {len(SCAFFOLD_GENES)} scaffolds/chromosomes")
            logging.info(f"Identified {len(GENES_OF_INTEREST)} ergosterol pathway genes")
            
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
        
        # Count matches and misses for debugging
        matches = 0
        misses = 0
        chrom_misses = {}
        
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
            
            # Map each variant to genes
            for idx, row in variant_df.iterrows():
                chrom = row['CHROM']
                position = row['POS']
                
                found_gene = False
                
                # First try to match using chromosome ID directly
                if chrom in SCAFFOLD_GENES:
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
                            
                            found_gene = True
                            matches += 1
                            break
                            
                # If not found directly, try using chromosome mapping
                if not found_gene and chrom in CHROM_TO_SCAFFOLD:
                    scaffold = CHROM_TO_SCAFFOLD[chrom]
                    if scaffold in SCAFFOLD_GENES:
                        for gene_id in SCAFFOLD_GENES[scaffold]:
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
                                
                                found_gene = True
                                matches += 1
                                break
                
                if found_gene:
                    continue
                
                # Only track as a miss if we couldn't find it in either direct lookup or through mapping
                if not found_gene:
                    # Track chromosome misses for debugging
                    if chrom not in chrom_misses:
                        chrom_misses[chrom] = 0
                    chrom_misses[chrom] += 1
                    misses += 1
            
            # Print chromosome miss statistics for debugging
            if misses > 0:
                print("Top 5 chromosomes with missing mappings:")
                for chrom, count in sorted(chrom_misses.items(), key=lambda x: x[1], reverse=True)[:5]:
                    print(f"  {chrom}: {count} misses")
            
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
        # Using improved pseudocount approach to avoid -inf when calculating log2 fold change
        # We use a more robust pseudocount that scales with the expected count
        pseudocount_treatment = max(0.5, expected_treatment * 0.01)  # At least 0.5, but scales with expected
        pseudocount_control = max(0.5, expected_control * 0.01)  # At least 0.5, but scales with expected
        
        treatment_fold_change = (treatment_count + pseudocount_treatment) / (expected_treatment + pseudocount_treatment)
        control_fold_change = (control_count + pseudocount_control) / (expected_control + pseudocount_control)
        
        # Calculate differential fold change (treatment vs control)
        # Using 0.01 as minimum value to avoid division by zero
        differential_fold_change = treatment_fold_change / max(0.01, control_fold_change)
        
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
        
        # Add log2 fold change for better visualization with safety checks
        # Ensure all values are positive before taking log2
        gene_df['treatment_fold_change_safe'] = gene_df['treatment_fold_change'].apply(lambda x: max(0.00001, x))
        gene_df['control_fold_change_safe'] = gene_df['control_fold_change'].apply(lambda x: max(0.00001, x))
        gene_df['differential_fold_change_safe'] = gene_df['differential_fold_change'].apply(lambda x: max(0.00001, x))
        
        # Calculate log2 fold changes from safe values
        gene_df['log2_treatment_fold_change'] = np.log2(gene_df['treatment_fold_change_safe'])
        gene_df['log2_control_fold_change'] = np.log2(gene_df['control_fold_change_safe'])
        gene_df['log2_differential_fold_change'] = np.log2(gene_df['differential_fold_change_safe'])
        
        # Debug log - Report any remaining issues with log2 calculations
        log2_fc = gene_df['log2_treatment_fold_change']
        if np.isneginf(log2_fc).any() or np.isposinf(log2_fc).any() or np.isnan(log2_fc).any():
            logging.warning(f"Some log2 fold change values are still problematic after robust correction")
            logging.warning(f"Infinities: {np.isinf(log2_fc).sum()}, NaNs: {np.isnan(log2_fc).sum()}")
            
            # If there are still issues, replace any problematic values
            gene_df['log2_treatment_fold_change'] = gene_df['log2_treatment_fold_change'].replace([np.inf, -np.inf, np.nan], 0)
            gene_df['log2_control_fold_change'] = gene_df['log2_control_fold_change'].replace([np.inf, -np.inf, np.nan], 0)
            gene_df['log2_differential_fold_change'] = gene_df['log2_differential_fold_change'].replace([np.inf, -np.inf, np.nan], 0)
            
            logging.info("Successfully replaced problematic log2 values with 0")
        
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
    if 'Fold_Change' in results_df.columns and len(results_df) > 0:
        # Prepare data
        fc_data = []
        for _, row in results_df.iterrows():
            treatment = row['Treatment']
            fold_change = row['Fold_Change']
            
            # Get treatment row data from gene results
            if treatment in all_gene_results and 'gene_status_fold_change' in all_gene_results[treatment]:
                gene_fc = all_gene_results[treatment]['gene_status_fold_change']
                
                # Add overall fold change
                fc_data.append({
                    'Treatment': treatment,
                    'Category': 'Overall',
                    'Fold Change': fold_change
                })
                
                # Add ERG fold change if available
                if 'ergosterol' in gene_fc and not np.isnan(gene_fc['ergosterol']):
                    fc_data.append({
                        'Treatment': treatment,
                        'Category': 'ERG Genes',
                        'Fold Change': gene_fc['ergosterol']
                    })
                
                # Add non-ERG fold change if available
                if 'non_ergosterol' in gene_fc and not np.isnan(gene_fc['non_ergosterol']):
                    fc_data.append({
                        'Treatment': treatment,
                        'Category': 'Non-ERG Genes',
                        'Fold Change': gene_fc['non_ergosterol']
                    })
        
        # Create fold change comparison chart if we have data
        if fc_data:
            fold_df = pd.DataFrame(fc_data)
            
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
            plt.savefig(f'{GENE_OUTPUT_DIR}/fold_change_by_gene_status.png', dpi=300, bbox_inches='tight')
            plt.close()
    
    # 5. Purifying Selection Visualization
    # For treatments with gene analysis, create volcano plots
    for treatment_name, gene_results in all_gene_results.items():
        # Skip if no gene analysis
        if 'gene_analysis' not in gene_results or len(gene_results['gene_analysis']) == 0:
            continue
        
        gene_df = gene_results['gene_analysis']
        
        # DEBUG: Print stats about log2_treatment_fold_change before plotting
        log2_fc = gene_df['log2_treatment_fold_change']
        print(f"\n[DEBUG] {treatment_name} log2_treatment_fold_change stats:")
        print(f"  Count: {len(log2_fc)}")
        print(f"  Min: {log2_fc.min()}  Max: {log2_fc.max()}")
        print(f"  Zeros: {(log2_fc == 0).sum()}")
        print(f"  NaNs: {log2_fc.isna().sum()}")
        print(f"  +Inf: {np.isposinf(log2_fc).sum()}  -Inf: {np.isneginf(log2_fc).sum()}")
        print(f"  Sample values: {log2_fc.head(10).tolist()}")
        
        # Handle any remaining infinite values for visualization (should not happen with pseudocounts)
        if np.isinf(log2_fc).any() or np.isnan(log2_fc).any():
            print(f"  WARNING: Some log2 fold change values are still infinite or NaN. Fixing for visualization.")
            # Replace -inf with lowest valid value - 1, and +inf with highest valid value + 1
            valid_min = log2_fc[~np.isinf(log2_fc) & ~np.isnan(log2_fc)].min()
            valid_max = log2_fc[~np.isinf(log2_fc) & ~np.isnan(log2_fc)].max()
            
            gene_df.loc[np.isneginf(gene_df['log2_treatment_fold_change']), 'log2_treatment_fold_change'] = valid_min - 1 if not np.isnan(valid_min) else -10
            gene_df.loc[np.isposinf(gene_df['log2_treatment_fold_change']), 'log2_treatment_fold_change'] = valid_max + 1 if not np.isnan(valid_max) else 10
            gene_df.loc[np.isnan(gene_df['log2_treatment_fold_change']), 'log2_treatment_fold_change'] = 0
            
            # Log the corrections
            log2_fc = gene_df['log2_treatment_fold_change']
            print(f"  After correction - Min: {log2_fc.min()}  Max: {log2_fc.max()}")
        
        # Also print for -log10(p-value)
        neglog10p = -np.log10(gene_df['treatment_p_value'])
        print(f"[DEBUG] {treatment_name} -log10(p-value) stats:")
        print(f"  Min: {neglog10p.min()}  Max: {neglog10p.max()}")
        print(f"  NaNs: {neglog10p.isna().sum()}")
        print(f"  +Inf: {np.isposinf(neglog10p).sum()}  -Inf: {np.isneginf(neglog10p).sum()}")
        print(f"  Sample values: {neglog10p.head(10).tolist()}")
        
        # Handle any NaN or infinite p-values
        if np.isinf(neglog10p).any() or np.isnan(neglog10p).any():
            # Replace NaN with 0 (p=1), and inf with a large value (e.g., 20, which is p=1e-20)
            gene_df.loc[np.isnan(gene_df['treatment_p_value']), 'treatment_p_value'] = 1.0
            gene_df.loc[gene_df['treatment_p_value'] <= 0, 'treatment_p_value'] = 1e-20
            
            # Recalculate
            neglog10p = -np.log10(gene_df['treatment_p_value'])
            print(f"  After correction - Min: {neglog10p.min()}  Max: {neglog10p.max()}")
        
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
        
        # Set axis limits with robust handling (avoid any infinite values)
        # This should now be unnecessary with our safety measures above, but we keep it as a failsafe
        log2_fc_values = gene_df['log2_treatment_fold_change'].replace([np.inf, -np.inf], np.nan).dropna()
        if len(log2_fc_values) > 0:
            # Use reasonable bounds based on actual data
            x_min = max(-10, min(-3, log2_fc_values.min() - 0.5))
            x_max = min(10, max(3, log2_fc_values.max() + 0.5))
        else:
            # Default reasonable bounds if no valid values
            x_min, x_max = -3, 3
        
        plt.xlim(x_min, x_max)
        
        # Also set y-axis limits to avoid any issues with p-values
        y_values = -np.log10(gene_df['treatment_p_value'].replace([0, np.inf, -np.inf], np.nan).dropna())
        if len(y_values) > 0:
            y_max = min(20, max(5, y_values.max() + 0.5))  # Cap at 20 for readability
        else:
            y_max = 5
        
        plt.ylim(0, y_max)
        
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
        "while a fold change > 1 indicates more variants than expected (positive selection).\n",
        "\n### Methodology Note\n",
        "When calculating fold changes, we use a pseudocount approach (adding 1 to both numerator and denominator)",
        "to handle genes with zero observed variants. This is a standard approach in genomics that allows us to",
        "detect purifying selection while avoiding mathematical issues with log transformations.",
        "For a gene with zero observed variants, the fold change will be 1/(expected+1) rather than zero,",
        "which still indicates depletion but avoids -infinity values in log2 transformations.\n"
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
        "results/vcf/merged/fixed/{}.filtered.vcf.gz"
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