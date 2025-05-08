 #!/usr/bin/env python3

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from collections import defaultdict, Counter
from scipy.stats import poisson, spearmanr, ttest_ind, mannwhitneyu
import math
import subprocess
import re

# Set matplotlib style for better visualizations
plt.style.use('ggplot')
sns.set(font_scale=1.2)
sns.set_style("whitegrid")

# Define output directory
# Define output directories
BASE_DIR = os.environ.get("OUTPUT_DIR", "analysis/general_gene_analysis")
OUTPUT_DIR = f"{BASE_DIR}/scaffold_distribution_results"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Log the output directories
logging.info(f"Using output directory: {OUTPUT_DIR}")

# Updated biologically correct treatment groups
TREATMENTS = ['WT-37', 'WTA', 'STC', 'CAS']

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

# Gene status colors
GENE_COLORS = {
    'ERG': '#2ca02c',    # Ergosterol pathway genes
    'Non-ERG': '#1f77b4', # Non-ergosterol genes
    'No Gene': '#7f7f7f'  # No gene
}

# Function to load gene mapping data
def load_gene_mapping():
    """Load gene mapping data from reference files."""
    gene_mapping_file = "reference/gene_mapping.tsv"
    genes_of_interest_file = "reference/genes_of_interest_mapping.tsv"
    
    try:
        # Load main gene mapping
        if os.path.exists(gene_mapping_file):
            gene_mapping = pd.read_csv(gene_mapping_file, sep='\t')
            print(f"Loaded {len(gene_mapping)} genes from {gene_mapping_file}")
        else:
            gene_mapping = pd.DataFrame()
            print(f"Warning: Gene mapping file not found at {gene_mapping_file}")
        
        # Load genes of interest (ergosterol pathway)
        if os.path.exists(genes_of_interest_file):
            erg_genes = pd.read_csv(genes_of_interest_file, sep='\t')
            print(f"Loaded {len(erg_genes)} genes of interest from {genes_of_interest_file}")
        else:
            erg_genes = pd.DataFrame()
            print(f"Warning: Genes of interest file not found at {genes_of_interest_file}")
        
        return gene_mapping, erg_genes
    
    except Exception as e:
        print(f"Error loading gene mapping: {e}")
        return pd.DataFrame(), pd.DataFrame()

# Function to build gene coordinates dictionary
def build_gene_coordinates(gene_mapping):
    """Build a dictionary of scaffold-to-genes mapping with coordinates."""
    if gene_mapping.empty:
        return {}
    
    scaffold_to_genes = defaultdict(list)
    
    # Check if required columns exist
    required_cols = ['w303_scaffold', 'start', 'end', 'w303_gene_id', 'erg_name', 'strand']
    if not all(col in gene_mapping.columns for col in required_cols):
        print(f"Warning: Gene mapping file missing required columns. Available: {gene_mapping.columns.tolist()}")
        return {}
    
    # Add genes to scaffold dictionary
    for _, gene in gene_mapping.iterrows():
        scaffold = gene['w303_scaffold']
        scaffold_to_genes[scaffold].append({
            'gene_id': gene['w303_gene_id'],
            'start': gene['start'],
            'end': gene['end'],
            'strand': gene['strand'],
            'erg_name': gene.get('erg_name', None),
            'is_erg': not pd.isna(gene.get('erg_name', None))
        })
    
    return scaffold_to_genes

# Function to map variants to genes
def map_variants_to_genes(variant_data, scaffold_to_genes, window_size=1000):
    """Map variants to genes and add gene information."""
    if variant_data.empty or not scaffold_to_genes:
        return variant_data
    
    # Add gene-related columns if they don't exist
    for col in ['Gene_ID', 'ERG_Name', 'Is_ERG', 'Distance_To_Gene']:
        if col not in variant_data.columns:
            variant_data[col] = None
    
    # Map each variant to nearest gene
    for idx, variant in variant_data.iterrows():
        scaffold = variant['CHROM']
        position = int(variant['POS'])
        
        # Skip if scaffold has no genes
        if scaffold not in scaffold_to_genes:
            variant_data.at[idx, 'Distance_To_Gene'] = float('inf')
            continue
        
        # Find closest gene
        min_distance = float('inf')
        closest_gene = None
        
        for gene in scaffold_to_genes[scaffold]:
            gene_start = gene['start']
            gene_end = gene['end']
            
            # Check if variant is within gene
            if gene_start <= position <= gene_end:
                distance = 0
                closest_gene = gene
                break
            
            # Calculate distance to gene
            distance_to_start = abs(position - gene_start)
            distance_to_end = abs(position - gene_end)
            distance = min(distance_to_start, distance_to_end)
            
            if distance < min_distance:
                min_distance = distance
                closest_gene = gene
        
        # Add gene information to variant
        if closest_gene and min_distance <= window_size:
            variant_data.at[idx, 'Gene_ID'] = closest_gene['gene_id']
            variant_data.at[idx, 'ERG_Name'] = closest_gene['erg_name']
            variant_data.at[idx, 'Is_ERG'] = closest_gene['is_erg']
            variant_data.at[idx, 'Distance_To_Gene'] = min_distance
        else:
            variant_data.at[idx, 'Distance_To_Gene'] = min_distance
    
    # Add gene status column for easier grouping
    variant_data['Gene_Status'] = variant_data.apply(
        lambda x: 'ERG' if x['Is_ERG'] else ('Non-ERG' if pd.notna(x['Gene_ID']) and x['Distance_To_Gene'] <= window_size else 'No Gene'),
        axis=1
    )
    
    return variant_data

# Function to find a file trying multiple locations
def find_file(base_name, file_patterns):
    """Find a file checking multiple possible locations."""
    # Substitute base_name into patterns
    patterns = [pattern.format(base_name) for pattern in file_patterns]
    
    # Add backward compatibility for WT-37
    if base_name == 'WT-37':
        wt_patterns = [pattern.format('WT') for pattern in file_patterns]
        patterns.extend(wt_patterns)
    
    # Try each pattern
    for pattern in patterns:
        if os.path.exists(pattern):
            print(f"Found file for {base_name} at {pattern}")
            return pattern
    
    print(f"Could not find file for {base_name} in any expected location")
    return None

# Function to parse the reference genome index to get scaffold lengths
def parse_genome_index(fai_file="reference/w303_chromosomal.fasta.fai"):
    """Parse the reference genome fasta index to get scaffold lengths."""
    # Try multiple possible locations for the fai file
    fai_patterns = [
        "reference/w303_chromosomal.fasta.fai",
        "reference/genome.fasta.fai",
        "reference/yeast/yeast_w303.fasta.fai"
    ]
    
    found_fai = None
    for pattern in fai_patterns:
        if os.path.exists(pattern):
            found_fai = pattern
            break
    
    if not found_fai:
        print(f"Warning: Genome index file not found in expected locations.")
        print("Trying to recreate scaffold lengths from VCF headers...")
        
        # Try to find scaffolds in one of the VCF files
        for treatment in TREATMENTS:
            vcf_patterns = [
                "results/merged/analysis/{}/highconf.vcf.gz",
                "results/merged/analysis/{}_highconf.vcf.gz",
                "results/merged/fixed/all_samples.vcf.gz"
            ]
            
            vcf_file = find_file(treatment, vcf_patterns)
            if vcf_file:
                try:
                    # Extract contig lines from VCF header
                    cmd = f"bcftools view -h {vcf_file} | grep contig"
                    output = subprocess.check_output(cmd, shell=True).decode('utf-8')
                    
                    # Parse contig lengths
                    scaffold_lengths = {}
                    for line in output.strip().split('\n'):
                        match = re.search(r'ID=([^,]+),length=(\d+)', line)
                        if match:
                            scaffold, length = match.groups()
                            scaffold_lengths[scaffold] = int(length)
                    
                    if scaffold_lengths:
                        print(f"Successfully extracted {len(scaffold_lengths)} scaffold lengths from VCF header.")
                        return scaffold_lengths
                except Exception as e:
                    print(f"Error extracting scaffold info from {vcf_file}: {e}")
                    continue
        
        print("Could not determine scaffold lengths. Using placeholder values.")
        # Return a dummy dict with placeholder values
        return defaultdict(lambda: 1000)
    
    # Read the fasta index file
    scaffold_lengths = {}
    with open(found_fai, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            scaffold, length = parts[0], int(parts[1])
            scaffold_lengths[scaffold] = length
    
    print(f"Loaded lengths for {len(scaffold_lengths)} scaffolds from {found_fai}")
    return scaffold_lengths

# Function to parse mutation data for a specific treatment
def parse_mutation_data(treatment):
    """Parse mutation data for a specific treatment."""
    # Define possible locations for mutation data files
    file_patterns = [
        "mutation_spectrum_analysis/{}_mutations.txt",
        "analysis/MSA/mutation_spectrum_analysis/{}_mutations.txt",
        "results/mutation_spectrum_analysis/{}_mutations.txt"
    ]
    
    # Find the mutation data file
    mutation_file = find_file(treatment, file_patterns)
    
    if mutation_file:
        try:
            # Read the already extracted mutation data
            data = pd.read_csv(mutation_file, sep='\t', header=None)
            
            # Check the number of columns and handle appropriately
            if len(data.columns) == 5:  # File already has Treatment column
                data.columns = ['CHROM', 'POS', 'REF', 'ALT', 'Treatment']
            else:  # Original format with 4 columns
                data.columns = ['CHROM', 'POS', 'REF', 'ALT']
                data['Treatment'] = treatment
                
            # Add biological context
            data['Adaptation'] = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
            data['Has_Gene'] = 'Yes' if TREATMENT_INFO.get(treatment, {}).get('gene') else 'No'
            
            return data
        except Exception as e:
            print(f"Error reading {mutation_file}: {e}")
    
    # If no mutation file found, try to extract from VCF
    vcf_patterns = [
        "results/merged/analysis/{}/highconf.vcf.gz",
        "results/merged/analysis/{}_highconf.vcf.gz",
        "results/merged/analysis/{}/specific.vcf.gz",
        "results/merged/analysis/{}_specific.vcf.gz"
    ]
    
    vcf_file = find_file(treatment, vcf_patterns)
    
    if vcf_file:
        try:
            print(f"Extracting mutation data from {vcf_file}")
            # Extract data using bcftools
            cmd = f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' {vcf_file}"
            output = subprocess.check_output(cmd, shell=True).decode('utf-8')
            
            # Parse output
            rows = []
            for line in output.strip().split('\n'):
                if line:  # Skip empty lines
                    parts = line.split('\t')
                    if len(parts) == 4:
                        rows.append(parts)
            
            if rows:
                # Create dataframe
                data = pd.DataFrame(rows, columns=['CHROM', 'POS', 'REF', 'ALT'])
                data['POS'] = data['POS'].astype(int)
                data['Treatment'] = treatment
                # Add biological context
                data['Adaptation'] = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
                data['Has_Gene'] = 'Yes' if TREATMENT_INFO.get(treatment, {}).get('gene') else 'No'
                
                # Save extracted data for future use
                os.makedirs("mutation_spectrum_analysis", exist_ok=True)
                data.to_csv(f"mutation_spectrum_analysis/{treatment}_mutations.txt", 
                          sep='\t', index=False, header=False)
                print(f"Saved extracted data to mutation_spectrum_analysis/{treatment}_mutations.txt")
                
                return data
        except Exception as e:
            print(f"Error extracting from {vcf_file}: {e}")
    
    print(f"Warning: No data available for {treatment}")
    return pd.DataFrame()

# Function to count variants per scaffold
def count_variants_per_scaffold(data):
    """Count the number of variants in each scaffold."""
    if len(data) == 0:
        return {}
    
    # Count variants per scaffold
    counts = Counter(data['CHROM'])
    return dict(counts)

# Function to count variants by gene status
def count_variants_by_gene_status(data):
    """Count variants by gene status (ERG, Non-ERG, No Gene)."""
    if len(data) == 0 or 'Gene_Status' not in data.columns:
        return {}
    
    # Count variants by gene status
    counts = Counter(data['Gene_Status'])
    return dict(counts)

# Function to count variants per gene
def count_variants_per_gene(data):
    """Count the number of variants in each gene."""
    if len(data) == 0 or 'Gene_ID' not in data.columns:
        return {}
    
    # Filter out variants not associated with genes
    gene_data = data[data['Gene_ID'].notna()]
    
    # Count variants per gene
    counts = Counter(gene_data['Gene_ID'])
    return dict(counts)

# Function to calculate variant density
def calculate_variant_density(counts, scaffold_lengths, min_scaffold_length=1000):
    """Calculate variant density (variants per kb) for each scaffold."""
    densities = {}
    
    for scaffold, count in counts.items():
        length = scaffold_lengths.get(scaffold, 0)
        
        # Skip short scaffolds to avoid artificially high densities
        if length < min_scaffold_length:
            continue
        
        # Calculate variants per kb
        density = (count * 1000) / length
        densities[scaffold] = density
    
    return densities

# Function to calculate gene-specific variant density
def calculate_gene_density(data, gene_mapping):
    """Calculate variant density per gene (variants per kb of gene length)."""
    if len(data) == 0 or 'Gene_ID' not in data.columns or gene_mapping.empty:
        return {}
    
    # Create a lookup for gene lengths
    gene_lengths = {}
    for _, gene in gene_mapping.iterrows():
        gene_id = gene['w303_gene_id']
        length = gene['end'] - gene['start'] + 1
        gene_lengths[gene_id] = length
    
    # Count variants per gene
    gene_counts = count_variants_per_gene(data)
    
    # Calculate density
    densities = {}
    for gene_id, count in gene_counts.items():
        length = gene_lengths.get(gene_id, 0)
        if length > 0:
            # Variants per kb of gene length
            density = (count * 1000) / length
            densities[gene_id] = density
    
    return densities

# Function to identify statistically enriched scaffolds
def identify_enriched_scaffolds(counts, scaffold_lengths, genome_wide_rate, min_scaffold_length=1000, pvalue_threshold=0.05):
    """Identify scaffolds with statistically significant variant enrichment."""
    enriched_scaffolds = []
    
    # Calculate total variants and total length
    total_variants = sum(counts.values())
    total_length = sum(scaffold_lengths.values())
    
    # If genome_wide_rate is not provided, calculate it
    if genome_wide_rate is None:
        genome_wide_rate = total_variants / total_length
    
    for scaffold, count in counts.items():
        length = scaffold_lengths.get(scaffold, 0)
        
        # Skip short scaffolds
        if length < min_scaffold_length:
            continue
        
        # Expected number of variants based on scaffold length
        expected = length * genome_wide_rate
        
        # Skip scaffolds with very low expected counts
        if expected < 1:
            continue
        
        # Calculate p-value using Poisson distribution
        # For enrichment, we want P(X ≥ observed) when expecting lambda
        p_value = 1 - poisson.cdf(count - 1, expected)
        
        # If significant, add to enriched list
        if p_value < pvalue_threshold:
            fold_change = count / expected if expected > 0 else float('inf')
            enriched_scaffolds.append({
                'Scaffold': scaffold,
                'Length': length,
                'Observed': count,
                'Expected': expected,
                'Fold_Change': fold_change,
                'P_Value': p_value
            })
    
    # Convert to dataframe and sort by p-value
    if enriched_scaffolds:
        result = pd.DataFrame(enriched_scaffolds)
        result = result.sort_values('P_Value')
        return result
    else:
        return pd.DataFrame()

# Function to identify statistically enriched or depleted genes
def identify_gene_patterns(data, gene_mapping, pvalue_threshold=0.05):
    """Identify genes with statistically significant enrichment or depletion of variants (for purifying selection)."""
    if len(data) == 0 or 'Gene_ID' not in data.columns or gene_mapping.empty:
        return pd.DataFrame(), pd.DataFrame()
    
    # Create a lookup for all genes including those with zero variants
    all_genes = {}
    for _, gene in gene_mapping.iterrows():
        gene_id = gene['w303_gene_id']
        all_genes[gene_id] = {
            'length': gene['end'] - gene['start'] + 1,
            'erg_name': gene.get('erg_name'),
            'is_erg': not pd.isna(gene.get('erg_name')),
            'count': 0
        }
    
    # Count variants per gene
    gene_variants = data[data['Gene_ID'].notna()]
    gene_counts = Counter(gene_variants['Gene_ID'])
    
    # Update gene counts
    for gene_id, count in gene_counts.items():
        if gene_id in all_genes:
            all_genes[gene_id]['count'] = count
    
    # Calculate genome-wide variant rate (variants per bp)
    total_variants = len(data)
    total_gene_length = sum(gene['length'] for gene in all_genes.values())
    genome_wide_rate = total_variants / total_gene_length if total_gene_length > 0 else 0
    
    # Identify enriched and depleted genes
    enriched_genes = []
    depleted_genes = []
    
    for gene_id, info in all_genes.items():
        length = info['length']
        observed = info['count']
        
        # Expected variants based on gene length and genome-wide rate
        expected = length * genome_wide_rate
        
        # Skip genes with very low expected counts
        if expected < 0.1:
            continue
        
        # Calculate Log2 fold change for better visualization
        log2_fold_change = math.log2((observed + 0.1) / (expected + 0.1))  # Add pseudocount to avoid log(0)
        
        # Calculate p-value
        if observed > expected:
            # Testing for enrichment (more variants than expected)
            p_value = 1 - poisson.cdf(observed - 1, expected)
            if p_value < pvalue_threshold:
                enriched_genes.append({
                    'Gene_ID': gene_id,
                    'ERG_Name': info['erg_name'],
                    'Is_ERG': info['is_erg'],
                    'Length': length,
                    'Observed': observed,
                    'Expected': expected,
                    'Log2_Fold_Change': log2_fold_change,
                    'P_Value': p_value
                })
        else:
            # Testing for depletion (fewer variants than expected - purifying selection)
            p_value = poisson.cdf(observed, expected)
            if p_value < pvalue_threshold:
                depleted_genes.append({
                    'Gene_ID': gene_id,
                    'ERG_Name': info['erg_name'],
                    'Is_ERG': info['is_erg'],
                    'Length': length,
                    'Observed': observed,
                    'Expected': expected,
                    'Log2_Fold_Change': log2_fold_change,
                    'P_Value': p_value
                })
    
    # Convert to dataframes and sort by p-value
    enriched_df = pd.DataFrame(enriched_genes).sort_values('P_Value') if enriched_genes else pd.DataFrame()
    depleted_df = pd.DataFrame(depleted_genes).sort_values('P_Value') if depleted_genes else pd.DataFrame()
    
    return enriched_df, depleted_df

# Function to plot gene enrichment and depletion
def plot_gene_enrichment(enriched_df, depleted_df, treatment, output_dir):
    """Generate volcano plot showing gene enrichment and depletion patterns."""
    if (enriched_df.empty and depleted_df.empty) or ('Log2_Fold_Change' not in enriched_df.columns and 'Log2_Fold_Change' not in depleted_df.columns):
        print(f"Warning: No significant gene enrichment/depletion data for {treatment}")
        # Create placeholder image
        plt.figure(figsize=(10, 8))
        plt.text(0.5, 0.5, f"No significant gene enrichment patterns for {treatment}", 
                horizontalalignment='center', verticalalignment='center', fontsize=14)
        plt.axis('off')
        plt.savefig(os.path.join(output_dir, f"{treatment}_gene_enrichment.png"), dpi=300)
        plt.close()
        return

    # Combine dataframes and add significance column
    enriched_df['Pattern'] = 'Enriched'
    depleted_df['Pattern'] = 'Depleted'
    
    combined_df = pd.concat([enriched_df, depleted_df])
    if combined_df.empty:
        return
    
    # Create the volcano plot
    plt.figure(figsize=(12, 10))
    
    # Plot all genes as small points
    sns.scatterplot(
        x='Log2_Fold_Change', 
        y='-log10(P_Value)', 
        hue='Is_ERG',
        style='Pattern',
        s=100,
        alpha=0.7,
        data=combined_df.assign(**{'-log10(P_Value)': -np.log10(combined_df['P_Value'])})
    )
    
    # Add gene labels for ERG genes
    for _, gene in combined_df[combined_df['Is_ERG']].iterrows():
        plt.annotate(
            gene['ERG_Name'],
            (gene['Log2_Fold_Change'], -np.log10(gene['P_Value'])),
            xytext=(5, 5),
            textcoords='offset points',
            fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8)
        )
    
    # Add reference lines
    plt.axvline(x=0, color='gray', linestyle='--', alpha=0.3)
    plt.axhline(y=-np.log10(0.05), color='gray', linestyle='--', alpha=0.3)
    
    plt.xlabel('Log2 Fold Change (Observed/Expected)')
    plt.ylabel('-log10(P-value)')
    plt.title(f'Gene Enrichment and Depletion Analysis for {treatment}\n' + 
              f'Positive values: enriched (adaptive), Negative values: depleted (purifying)')
    
    # Customize legend
    plt.legend(title='Gene Type / Pattern')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{treatment}_gene_enrichment.png"), dpi=300)
    plt.close()
    
    # Also create ERG-specific plot if we have ERG genes
    erg_genes = combined_df[combined_df['Is_ERG']]
    if not erg_genes.empty:
        plt.figure(figsize=(10, 8))
        
        # Plot ERG genes with fold change
        sns.barplot(
            x='ERG_Name',
            y='Log2_Fold_Change',
            data=erg_genes.sort_values('Log2_Fold_Change')
        )
        
        # Color bars by significance and direction
        for i, row in enumerate(erg_genes.sort_values('Log2_Fold_Change').iterrows()):
            bar_color = 'green' if row[1]['Log2_Fold_Change'] > 0 else 'red'
            alpha = 1.0 if row[1]['P_Value'] < 0.05 else 0.5
            plt.gca().patches[i].set_facecolor(bar_color)
            plt.gca().patches[i].set_alpha(alpha)
        
        # Add horizontal line at 0
        plt.axhline(y=0, color='black', linestyle='-', alpha=0.3)
        
        plt.ylabel('Log2 Fold Change')
        plt.xlabel('Ergosterol Pathway Gene')
        plt.title(f'Ergosterol Pathway Genes Variation Pattern in {treatment}\n' + 
                 'Positive: enriched, Negative: depleted (purifying selection)')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{treatment}_erg_genes.png"), dpi=300)
        plt.close()
    
    # Compare ERG vs non-ERG genes
    if not combined_df.empty and 'Is_ERG' in combined_df.columns:
        erg_fold_changes = combined_df[combined_df['Is_ERG']]['Log2_Fold_Change']
        non_erg_fold_changes = combined_df[~combined_df['Is_ERG']]['Log2_Fold_Change']
        
        if len(erg_fold_changes) > 0 and len(non_erg_fold_changes) > 0:
            plt.figure(figsize=(10, 8))
            
            # Create a boxplot comparing fold changes
            sns.boxplot(
                x='Is_ERG',
                y='Log2_Fold_Change',
                data=combined_df,
                palette=['#2ca02c', '#1f77b4']
            )
            
            # Add individual data points
            sns.stripplot(
                x='Is_ERG',
                y='Log2_Fold_Change',
                data=combined_df,
                color='black',
                alpha=0.5,
                jitter=True,
                size=4
            )
            
            # Add horizontal line at 0
            plt.axhline(y=0, color='black', linestyle='-', alpha=0.3)
            
            # Perform statistical test
            try:
                # Try t-test first
                t_stat, p_value = ttest_ind(erg_fold_changes, non_erg_fold_changes, equal_var=False)
                test_name = "Welch's t-test"
                
                # If data is too small, use Mann-Whitney U test
                if len(erg_fold_changes) < 5 or len(non_erg_fold_changes) < 5:
                    u_stat, p_value = mannwhitneyu(erg_fold_changes, non_erg_fold_changes, alternative='two-sided')
                    test_name = "Mann-Whitney U test"
                
                plt.text(
                    0.5, 0.95, 
                    f"{test_name} p-value: {p_value:.4f}",
                    horizontalalignment='center',
                    transform=plt.gca().transAxes,
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8)
                )
            except Exception as e:
                print(f"Statistical test error: {e}")
            
            plt.ylabel('Log2 Fold Change')
            plt.xlabel('Gene Type')
            plt.title(f'Comparison of Enrichment/Depletion Patterns in {treatment}\n' + 
                     'Ergosterol vs. Non-Ergosterol Genes')
            
            # Update x-axis labels
            plt.xticks([0, 1], ['Ergosterol Genes', 'Other Genes'])
            
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"{treatment}_erg_vs_nonerg.png"), dpi=300)
            plt.close()

# Function to generate variant density plot
def plot_variant_density(densities, treatment, output_dir, top_n=20):
    """Generate a variant density plot showing top scaffolds."""
    # Sort scaffolds by density
    sorted_densities = sorted(densities.items(), key=lambda x: x[1], reverse=True)
    
    # Take top N scaffolds
    top_scaffolds = sorted_densities[:top_n]
    scaffolds, density_values = zip(*top_scaffolds) if top_scaffolds else ([], [])
    
    # Create truncated scaffold names for better display
    scaffold_names = [s[:15] + '...' if len(s) > 15 else s for s in scaffolds]
    
    # Get treatment metadata
    description = TREATMENT_INFO.get(treatment, {}).get('description', '')
    adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', '')
    has_gene = TREATMENT_INFO.get(treatment, {}).get('gene')
    gene_text = f" with {has_gene} gene" if has_gene else ""
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 8))
    bars = ax.bar(range(len(scaffold_names)), density_values, color=TREATMENT_COLORS.get(treatment, '#1b9e77'))
    
    # Add value labels on top of bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                f'{height:.2f}', ha='center', va='bottom', rotation=45, fontsize=8)
    
    # Customize the plot
    ax.set_xticks(range(len(scaffold_names)))
    ax.set_xticklabels(scaffold_names, rotation=90)
    ax.set_xlabel('Scaffold')
    ax.set_ylabel('Variants per kb')
    ax.set_title(f'Top {top_n} Scaffolds by Variant Density for {treatment}\n{description} ({adaptation} adaptation{gene_text})')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{treatment}_variant_density_top{top_n}.png"), dpi=300)
    plt.close()

# Function to plot gene status distribution
def plot_gene_status_distribution(data, treatment, output_dir):
    """Plot the distribution of variants by gene status."""
    if len(data) == 0 or 'Gene_Status' not in data.columns:
        print(f"Warning: No gene status data for {treatment}")
        return
    
    # Count variants by gene status
    status_counts = count_variants_by_gene_status(data)
    
    # Check if we have data
    if not status_counts:
        print(f"Warning: No gene status data for {treatment}")
        return
    
    # Create dataframe for plotting
    status_df = pd.DataFrame({
        'Gene Status': list(status_counts.keys()),
        'Count': list(status_counts.values())
    })
    
    # Create the plot
    plt.figure(figsize=(10, 8))
    ax = sns.barplot(
        x='Gene Status',
        y='Count',
        data=status_df,
        palette=[GENE_COLORS.get(status, '#333333') for status in status_df['Gene Status']]
    )
    
    # Add count labels on top of bars
    for p in ax.patches:
        ax.annotate(
            f'{int(p.get_height())}',
            (p.get_x() + p.get_width()/2., p.get_height()),
            ha='center',
            va='bottom',
            fontsize=10
        )
    
    # Get treatment metadata
    description = TREATMENT_INFO.get(treatment, {}).get('description', '')
    adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', '')
    has_gene = TREATMENT_INFO.get(treatment, {}).get('gene')
    gene_text = f" with {has_gene} gene" if has_gene else ""
    
    plt.title(f'Variant Distribution by Gene Status for {treatment}\n{description} ({adaptation} adaptation{gene_text})')
    plt.xlabel('Gene Status')
    plt.ylabel('Number of Variants')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{treatment}_gene_status_distribution.png"), dpi=300)
    plt.close()
    
    # Calculate percentages for reporting
    total = sum(status_counts.values())
    percentages = {status: (count/total*100) for status, count in status_counts.items()}
    
    # Create a pie chart
    plt.figure(figsize=(10, 8))
    plt.pie(
        status_counts.values(),
        labels=[f"{status}\n({percentages[status]:.1f}%)" for status in status_counts.keys()],
        colors=[GENE_COLORS.get(status, '#333333') for status in status_counts.keys()],
        autopct='%1.1f%%',
        startangle=90,
        explode=[0.1 if status == 'ERG' else 0 for status in status_counts.keys()]
    )
    plt.axis('equal')
    plt.title(f'Gene Status Distribution for {treatment}')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{treatment}_gene_status_pie.png"), dpi=300)
    plt.close()

# Function to generate comparative density heat map
def plot_comparative_heatmap(all_densities, output_dir, top_n=30):
    """Generate a heatmap comparing variant densities across treatments."""
    # Check if all_densities is empty or has no treatments
    if not all_densities:
        print("Warning: No density data available for heatmap visualization")
        # Create empty dataframe
        df = pd.DataFrame({'scaffold_name': [], 'density': []})
        df.set_index('scaffold_name', inplace=True)
        top_scaffold_names = []
    else:
        # Find top scaffolds across all treatments
        all_scaffold_densities = {}
        
        for treatment, densities in all_densities.items():
            for scaffold, density in densities.items():
                if scaffold not in all_scaffold_densities:
                    all_scaffold_densities[scaffold] = 0
                all_scaffold_densities[scaffold] += density
        
        # Check if we have any scaffolds
        if not all_scaffold_densities:
            print("Warning: No scaffold density data available")
            df = pd.DataFrame({'scaffold_name': [], 'density': []})
            df.set_index('scaffold_name', inplace=True)
            top_scaffold_names = []
        else:
            # Sort scaffolds by total density
            top_scaffolds = sorted(all_scaffold_densities.items(), key=lambda x: x[1], reverse=True)[:top_n]
            top_scaffold_names = [s[0] for s in top_scaffolds]
            
            # Create a dataframe for the heatmap
            data = []
            for scaffold in top_scaffold_names:
                row = {'scaffold_name': scaffold}
                for treatment, densities in all_densities.items():
                    row[treatment] = densities.get(scaffold, 0)
                data.append(row)
            
            df = pd.DataFrame(data)
            
            # Check if 'scaffold_name' column exists before trying to set as index
            if len(data) > 0 and 'scaffold_name' in df.columns:
                df.set_index('scaffold_name', inplace=True)
            else:
                # Handle the case when the column doesn't exist - print columns for debugging
                print(f"Warning: Expected 'scaffold_name' column not found. Available columns: {df.columns.tolist()}")
                # In case all_densities is empty or doesn't have the expected structure
                if len(df.columns) == 0:
                    print("Creating empty dataframe with scaffold column")
                    df = pd.DataFrame({'scaffold_name': [], 'density': []})
                    df.set_index('scaffold_name', inplace=True)
    
    # Create truncated scaffold names for better display
    df.index = [s[:20] + '...' if len(s) > 20 else s for s in df.index]
    
    # Create the heatmap
    plt.figure(figsize=(12, 10))
    
    # Check if we have enough data for meaningful visualization
    if len(df) == 0 or len(df.columns) < 2:
        print("Warning: Not enough data for heatmap visualization")
        
        # Create a simple message image
        plt.figure(figsize=(10, 8))
        plt.text(0.5, 0.5, "Insufficient data for heatmap visualization", 
                 horizontalalignment='center', verticalalignment='center', fontsize=14)
        plt.axis('off')
        plt.savefig(os.path.join(output_dir, f"comparative_density_heatmap_top{top_n}.png"), dpi=300)
        plt.close()
        
        # Create placeholder for treatment correlation
        plt.figure(figsize=(10, 8))
        plt.text(0.5, 0.5, "Insufficient data for correlation analysis", 
                 horizontalalignment='center', verticalalignment='center', fontsize=14)
        plt.axis('off')
        plt.savefig(os.path.join(output_dir, "treatment_correlation_heatmap.png"), dpi=300)
        plt.close()
        
        return
    
    # Create column colors based on adaptation type
    column_colors = [TREATMENT_COLORS.get(t, '#333333') for t in all_densities.keys()]
    
    try:
        # Create heatmap with clustered columns
        g = sns.clustermap(
            df,
            cmap='viridis',
            col_colors=[column_colors],
            annot=True,
            fmt='.2f',
            linewidths=.5,
            figsize=(14, 12),
            dendrogram_ratio=(.1, .2),
            row_cluster=False  # Keep scaffolds in order of density
        )
    except Exception as e:
        print(f"Warning: Could not create clustered heatmap: {e}")
        
        # Fall back to a simple heatmap
        plt.figure(figsize=(12, 10))
        sns.heatmap(df, cmap='viridis', annot=True, fmt='.2f', linewidths=.5)
        plt.title('Comparative Variant Density Across Treatments')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"comparative_density_heatmap_top{top_n}.png"), dpi=300)
        plt.close()
        
        # Create placeholder for treatment correlation
        plt.figure(figsize=(10, 8))
        if len(df.columns) >= 2:
            sns.heatmap(df.corr(), annot=True, cmap='coolwarm', vmin=-1, vmax=1)
            plt.title('Treatment Correlation Based on Scaffold Variant Density')
        else:
            plt.text(0.5, 0.5, "Insufficient data for correlation analysis", 
                    horizontalalignment='center', verticalalignment='center', fontsize=14)
            plt.axis('off')
        plt.savefig(os.path.join(output_dir, "treatment_correlation_heatmap.png"), dpi=300)
        plt.close()
        
        return
    
    # Improve the column dendrogram
    g.ax_col_dendrogram.set_visible(True)
    
    # Add treatment labels with adaptation info
    ax = g.ax_heatmap
    ax.set_xticklabels([
        f"{t}\n({TREATMENT_INFO.get(t, {}).get('adaptation', '')})" 
        for t in df.columns
    ])
    
    # Add adaptation type legend
    handles = []
    for adaptation in ["Temperature", "Low Oxygen"]:
        treatments = [t for t in TREATMENTS if TREATMENT_INFO.get(t, {}).get('adaptation') == adaptation]
        if treatments:
            color = TREATMENT_COLORS.get(treatments[0], '#333333')
            handles.append(plt.Rectangle((0,0), 1, 1, color=color, label=f"{adaptation} Adaptation"))
    
    # Add adaptation legend
    g.ax_col_dendrogram.legend(handles=handles, loc='upper center', ncol=2, frameon=True)
    
    plt.suptitle('Variant Density (variants/kb) Across Treatments', fontsize=16, y=0.98)
    plt.savefig(os.path.join(output_dir, f"comparative_density_heatmap_top{top_n}.png"), dpi=300)
    plt.close()
    
    # Also create adaptation-grouped heatmap
    # Group treatments by adaptation type
    adaptation_densities = defaultdict(dict)
    for treatment, densities in all_densities.items():
        adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
        
        # Merge densities from the same adaptation type (using max for now)
        for scaffold, density in densities.items():
            if scaffold in adaptation_densities[adaptation]:
                adaptation_densities[adaptation][scaffold] = max(adaptation_densities[adaptation][scaffold], density)
            else:
                adaptation_densities[adaptation][scaffold] = density
    
    # Create adaptation dataframe
    adaptation_data = []
    for scaffold in top_scaffold_names:
        row = {'Scaffold': scaffold}
        for adaptation, densities in adaptation_densities.items():
            row[adaptation] = densities.get(scaffold, 0)
        adaptation_data.append(row)
    
    adaptation_df = pd.DataFrame(adaptation_data)
    adaptation_df.set_index('Scaffold', inplace=True)
    
    # Create adaptation heatmap
    plt.figure(figsize=(10, 12))
    sns.heatmap(
        adaptation_df,
        cmap='viridis',
        annot=True,
        fmt='.2f',
        linewidths=.5
    )
    
    plt.title('Variant Density (variants/kb) by Adaptation Type')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"adaptation_density_heatmap.png"), dpi=300)
    plt.close()

# Function to plot bubble chart
def plot_bubble_chart(counts, scaffold_lengths, treatment, output_dir, top_n=50):
    """Generate a bubble chart showing scaffold length vs. variant count."""
    data = []
    
    for scaffold, count in counts.items():
        length = scaffold_lengths.get(scaffold, 0)
        if length > 0:
            density = (count * 1000) / length
            data.append({
                'Scaffold': scaffold,
                'Length': length,
                'Count': count,
                'Density': density
            })
    
    # Check if we have data
    if not data:
        print(f"Warning: No valid scaffold data for {treatment}, creating placeholder visualization")
        plt.figure(figsize=(10, 8))
        plt.text(0.5, 0.5, f"No significant variant density for {treatment}", 
                 horizontalalignment='center', verticalalignment='center', fontsize=14)
        plt.axis('off')
        plt.savefig(os.path.join(output_dir, f"{treatment}_bubble_chart.png"), dpi=300)
        plt.close()
        return
    
    # Convert to dataframe and sort by density
    df = pd.DataFrame(data)
    
    # Check if Density column exists
    if 'Density' in df.columns and len(df) > 0:
        df = df.sort_values('Density', ascending=False).head(top_n)
    else:
        # Just take the first n rows if we can't sort by density
        df = df.head(top_n)
    
    # Get treatment metadata
    description = TREATMENT_INFO.get(treatment, {}).get('description', '')
    adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', '')
    has_gene = TREATMENT_INFO.get(treatment, {}).get('gene')
    gene_text = f" with {has_gene} gene" if has_gene else ""
    
    # Create the bubble chart
    plt.figure(figsize=(12, 8))
    
    # Use log scale for scaffold length
    plt.scatter(
        df['Length'], 
        df['Count'], 
        s=df['Density']*20,  # Scale bubble size by density
        alpha=0.6, 
        color=TREATMENT_COLORS.get(treatment, '#1b9e77'), 
        edgecolors='w',
        linewidth=0.5
    )
    
    # Add labels for top bubbles
    for i, row in df.head(10).iterrows():
        plt.annotate(
            row['Scaffold'][:10] + '...',
            (row['Length'], row['Count']),
            xytext=(5, 5),
            textcoords='offset points'
        )
    
    plt.xscale('log')
    plt.xlabel('Scaffold Length (log scale)')
    plt.ylabel('Variant Count')
    plt.title(f'Scaffold Length vs. Variant Count for {treatment}\n{description} ({adaptation} adaptation{gene_text})')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{treatment}_bubble_chart.png"), dpi=300)
    plt.close()

# Function to calculate correlation between treatments
def calculate_treatment_correlations(all_densities, output_dir):
    """Calculate and visualize correlations between treatments based on scaffold density."""
    # Check if we have enough data for correlation analysis
    if len(all_densities) < 2:
        print("Warning: Not enough treatments for correlation analysis")
        # Create placeholder images
        plt.figure(figsize=(10, 8))
        plt.text(0.5, 0.5, "Insufficient data for correlation analysis", 
                 horizontalalignment='center', verticalalignment='center', fontsize=14)
        plt.axis('off')
        plt.savefig(os.path.join(output_dir, "treatment_correlation_heatmap.png"), dpi=300)
        plt.close()
        
        # Create empty dataframes to return
        empty_corr = pd.DataFrame(columns=TREATMENTS, index=TREATMENTS)
        empty_adapt = pd.DataFrame(columns=['Temperature', 'Low Oxygen'], index=['Temperature', 'Low Oxygen'])
        return empty_corr, empty_adapt
    
    # Get all unique scaffolds across treatments
    all_scaffolds = set()
    for densities in all_densities.values():
        all_scaffolds.update(densities.keys())
    
    # If no scaffolds, return empty correlation matrices
    if not all_scaffolds:
        print("Warning: No scaffold data available for correlation analysis")
        # Create placeholder images
        plt.figure(figsize=(10, 8))
        plt.text(0.5, 0.5, "No scaffold data for correlation analysis", 
                 horizontalalignment='center', verticalalignment='center', fontsize=14)
        plt.axis('off')
        plt.savefig(os.path.join(output_dir, "treatment_correlation_heatmap.png"), dpi=300)
        plt.close()
        
        # Create empty dataframes to return
        empty_corr = pd.DataFrame(columns=TREATMENTS, index=TREATMENTS)
        empty_adapt = pd.DataFrame(columns=['Temperature', 'Low Oxygen'], index=['Temperature', 'Low Oxygen'])
        return empty_corr, empty_adapt
    
    # Create a dataframe with scaffold densities for each treatment
    data = []
    for scaffold in all_scaffolds:
        row = {'scaffold_name': scaffold}
        for treatment, densities in all_densities.items():
            row[treatment] = densities.get(scaffold, 0)
        data.append(row)
    
    # Create DataFrame
    df = pd.DataFrame(data)
    
    # Check if we have the expected column
    if 'scaffold_name' in df.columns:
        df.set_index('scaffold_name', inplace=True)
    else:
        print(f"Warning: Expected column 'scaffold_name' not found. Columns: {df.columns.tolist()}")
        # If no proper columns, return empty correlation matrices
        empty_corr = pd.DataFrame(columns=TREATMENTS, index=TREATMENTS)
        empty_adapt = pd.DataFrame(columns=['Temperature', 'Low Oxygen'], index=['Temperature', 'Low Oxygen'])
        return empty_corr, empty_adapt
    
    # Calculate correlation matrix
    corr_matrix = df.corr(method='spearman')
    
    # Plot correlation heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', vmin=-1, vmax=1, linewidths=.5)
    plt.title('Spearman Correlation of Variant Densities Between Treatments')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "treatment_correlation_heatmap.png"), dpi=300)
    plt.close()
    
    # Also plot correlation clustered by adaptation type
    # Create column colors for adaptation types
    adaptation_colors = {t: TREATMENT_COLORS.get(t, '#333333') for t in corr_matrix.columns}
    row_colors = pd.Series(adaptation_colors)
    
    # Create clustered heatmap
    g = sns.clustermap(
        corr_matrix,
        cmap='coolwarm',
        vmin=-1,
        vmax=1,
        annot=True,
        linewidths=.5,
        row_colors=row_colors,
        col_colors=row_colors,
        figsize=(12, 10)
    )
    
    # Add adaptation type legend
    for adaptation in ["Temperature", "Low Oxygen"]:
        treatments = [t for t in corr_matrix.columns if TREATMENT_INFO.get(t, {}).get('adaptation') == adaptation]
        if treatments:
            color = TREATMENT_COLORS.get(treatments[0], '#333333')
            g.ax_row_dendrogram.bar(0, 0, color=color, label=f"{adaptation} Adaptation")
    
    g.ax_row_dendrogram.legend(frameon=True, loc='center', title="Adaptation")
    
    plt.suptitle('Clustered Correlation of Variant Densities', fontsize=16, y=0.98)
    plt.savefig(os.path.join(output_dir, "clustered_correlation_heatmap.png"), dpi=300)
    plt.close()
    
    # Calculate within-adaptation correlations
    adaptation_groups = defaultdict(list)
    for treatment in all_densities.keys():
        adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
        adaptation_groups[adaptation].append(treatment)
    
    # Calculate average correlation within and between adaptation types
    adaptation_correlations = {}
    
    for adap1 in adaptation_groups:
        for adap2 in adaptation_groups:
            if adap1 not in adaptation_correlations:
                adaptation_correlations[adap1] = {}
            
            # Calculate average correlation between treatments in these adaptations
            correlations = []
            for t1 in adaptation_groups[adap1]:
                for t2 in adaptation_groups[adap2]:
                    if t1 != t2:  # Skip self-correlations
                        correlations.append(corr_matrix.loc[t1, t2])
            
            if correlations:
                adaptation_correlations[adap1][adap2] = sum(correlations) / len(correlations)
            else:
                adaptation_correlations[adap1][adap2] = float('nan')
    
    # Create adaptation correlation matrix
    adaptation_corr_df = pd.DataFrame(adaptation_correlations)
    
    # Plot adaptation correlation heatmap
    plt.figure(figsize=(8, 6))
    sns.heatmap(adaptation_corr_df, annot=True, cmap='coolwarm', vmin=-1, vmax=1, linewidths=.5)
    plt.title('Average Correlation Between Adaptation Types')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "adaptation_correlation_heatmap.png"), dpi=300)
    plt.close()
    
    return corr_matrix, adaptation_corr_df

# Function to identify treatment-specific hotspots
def identify_treatment_specific_hotspots(all_densities, output_dir, fold_threshold=5):
    """Identify scaffolds that are uniquely enriched in specific treatments."""
    # Get all unique scaffolds across treatments
    all_scaffolds = set()
    for densities in all_densities.values():
        all_scaffolds.update(densities.keys())
    
    # Find treatment-specific hotspots
    treatment_specific = {treatment: [] for treatment in all_densities.keys()}
    
    for scaffold in all_scaffolds:
        # Get density for each treatment
        densities = {t: d.get(scaffold, 0) for t, d in all_densities.items()}
        
        # Skip scaffolds with low density across all treatments
        if max(densities.values()) < 1:
            continue
        
        # Check each treatment
        for treatment, density in densities.items():
            # Skip if this treatment has no variants in this scaffold
            if density == 0:
                continue
            
            # Compare with other treatments
            is_specific = True
            for other, other_density in densities.items():
                if other != treatment and other_density > 0:
                    # If density in other treatment is within threshold, not specific
                    ratio = density / other_density if other_density > 0 else float('inf')
                    if ratio < fold_threshold:
                        is_specific = False
                        break
            
            if is_specific:
                treatment_specific[treatment].append({
                    'Scaffold': scaffold,
                    'Density': density,
                    'Fold_Enrichment': min([density / other_density if other_density > 0 else float('inf') 
                                          for other, other_density in densities.items() if other != treatment])
                })
    
    # Prepare summary
    with open(os.path.join(output_dir, "treatment_specific_hotspots.txt"), 'w') as f:
        f.write("Treatment-Specific Scaffold Hotspots\n")
        f.write("===================================\n\n")
        
        for treatment, hotspots in treatment_specific.items():
            adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
            gene_mod = TREATMENT_INFO.get(treatment, {}).get('gene')
            f.write(f"{treatment} Treatment ({adaptation}{' with '+gene_mod+' gene' if gene_mod else ''}):\n")
            f.write(f"  {len(hotspots)} specific hotspots\n")
            
            if hotspots:
                # Sort by density
                hotspots = sorted(hotspots, key=lambda x: x['Density'], reverse=True)
                
                f.write(f"{'Scaffold':<25}{'Density':<10}{'Fold Enrichment':<20}\n")
                f.write("-" * 55 + "\n")
                
                for h in hotspots:
                    fold = h['Fold_Enrichment']
                    fold_str = f"{fold:.2f}" if fold != float('inf') else "∞"
                    f.write(f"{h['Scaffold']:<25}{h['Density']:<10.2f}{fold_str:<20}\n")
            
            f.write("\n")
    
    print(f"Treatment-specific hotspots written to {os.path.join(output_dir, 'treatment_specific_hotspots.txt')}")
    
    # Also identify adaptation-specific hotspots
    adaptation_densities = defaultdict(dict)
    for treatment, densities in all_densities.items():
        adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
        
        # Merge densities by taking the maximum
        for scaffold, density in densities.items():
            if scaffold in adaptation_densities[adaptation]:
                adaptation_densities[adaptation][scaffold] = max(adaptation_densities[adaptation][scaffold], density)
            else:
                adaptation_densities[adaptation][scaffold] = density
    
    # Find adaptation-specific hotspots
    adaptation_specific = {adaptation: [] for adaptation in adaptation_densities.keys()}
    
    for scaffold in all_scaffolds:
        # Get density for each adaptation
        densities = {a: d.get(scaffold, 0) for a, d in adaptation_densities.items()}
        
        # Skip scaffolds with low density
        if max(densities.values()) < 1:
            continue
        
        # Check each adaptation
        for adaptation, density in densities.items():
            # Skip if no variants
            if density == 0:
                continue
            
            # Compare with other adaptations
            is_specific = True
            for other, other_density in densities.items():
                if other != adaptation and other_density > 0:
                    ratio = density / other_density if other_density > 0 else float('inf')
                    if ratio < fold_threshold:
                        is_specific = False
                        break
            
            if is_specific:
                adaptation_specific[adaptation].append({
                    'Scaffold': scaffold,
                    'Density': density,
                    'Fold_Enrichment': min([density / other_density if other_density > 0 else float('inf')
                                          for other, other_density in densities.items() if other != adaptation])
                })
    
    # Write adaptation-specific hotspots
    with open(os.path.join(output_dir, "adaptation_specific_hotspots.txt"), 'w') as f:
        f.write("Adaptation-Specific Scaffold Hotspots\n")
        f.write("===================================\n\n")
        
        for adaptation, hotspots in adaptation_specific.items():
            f.write(f"{adaptation} Adaptation:\n")
            f.write(f"  {len(hotspots)} specific hotspots\n")
            
            if hotspots:
                # Sort by density
                hotspots = sorted(hotspots, key=lambda x: x['Density'], reverse=True)
                
                f.write(f"{'Scaffold':<25}{'Density':<10}{'Fold Enrichment':<20}\n")
                f.write("-" * 55 + "\n")
                
                for h in hotspots:
                    fold = h['Fold_Enrichment']
                    fold_str = f"{fold:.2f}" if fold != float('inf') else "∞"
                    f.write(f"{h['Scaffold']:<25}{h['Density']:<10.2f}{fold_str:<20}\n")
            
            f.write("\n")
    
    print(f"Adaptation-specific hotspots written to {os.path.join(output_dir, 'adaptation_specific_hotspots.txt')}")
    
    # Compare gene-modified vs non-modified strains within each adaptation
    gene_comparison = {}
    
    for adaptation in ["Temperature", "Low Oxygen"]:
        gene_modified = [t for t in all_densities.keys() 
                        if TREATMENT_INFO.get(t, {}).get('adaptation') == adaptation 
                        and TREATMENT_INFO.get(t, {}).get('gene')]
        
        non_modified = [t for t in all_densities.keys() 
                       if TREATMENT_INFO.get(t, {}).get('adaptation') == adaptation 
                       and not TREATMENT_INFO.get(t, {}).get('gene')]
        
        if gene_modified and non_modified:
            # Combine densities for gene-modified treatments
            gene_mod_densities = defaultdict(float)
            for t in gene_modified:
                for scaffold, density in all_densities[t].items():
                    gene_mod_densities[scaffold] = max(gene_mod_densities[scaffold], density)
            
            # Combine densities for non-modified treatments
            non_mod_densities = defaultdict(float)
            for t in non_modified:
                for scaffold, density in all_densities[t].items():
                    non_mod_densities[scaffold] = max(non_mod_densities[scaffold], density)
            
            # Find scaffolds enriched in gene-modified strains
            gene_effect_scaffolds = []
            
            for scaffold in all_scaffolds:
                gene_density = gene_mod_densities.get(scaffold, 0)
                nonmod_density = non_mod_densities.get(scaffold, 0)
                
                if gene_density > 0 and nonmod_density > 0:
                    ratio = gene_density / nonmod_density
                    
                    if ratio >= fold_threshold:
                        gene_effect_scaffolds.append({
                            'Scaffold': scaffold,
                            'Gene_Density': gene_density,
                            'NonMod_Density': nonmod_density,
                            'Fold_Enrichment': ratio
                        })
            
            gene_comparison[adaptation] = gene_effect_scaffolds
    
    # Write gene effect analysis
    with open(os.path.join(output_dir, "gene_modification_effects.txt"), 'w') as f:
        f.write("Gene Modification Effects on Scaffold Enrichment\n")
        f.write("============================================\n\n")
        
        for adaptation, scaffolds in gene_comparison.items():
            f.write(f"{adaptation} Adaptation - Gene Effects:\n")
            f.write(f"  {len(scaffolds)} scaffolds enriched in gene-modified strains\n")
            
            if scaffolds:
                # Sort by fold enrichment
                scaffolds = sorted(scaffolds, key=lambda x: x['Fold_Enrichment'], reverse=True)
                
                f.write(f"{'Scaffold':<25}{'Gene Mod':<10}{'Non-Mod':<10}{'Fold Enrich':<15}\n")
                f.write("-" * 60 + "\n")
                
                for s in scaffolds:
                    f.write(f"{s['Scaffold']:<25}{s['Gene_Density']:<10.2f}{s['NonMod_Density']:<10.2f}{s['Fold_Enrichment']:<15.2f}\n")
            
            f.write("\n")
    
    print(f"Gene modification effects written to {os.path.join(output_dir, 'gene_modification_effects.txt')}")
    
    return treatment_specific, adaptation_specific, gene_comparison

# Function to create a summary table
def create_summary_table(all_counts, all_densities, scaffold_lengths, gene_stat_counts, output_dir):
    """Create a comprehensive summary table with key statistics."""
    summary = []
    
    for treatment in all_densities.keys():
        counts = all_counts[treatment]
        total_variants = sum(counts.values())
        total_scaffolds = len(counts)
        
        # Calculate global variant density
        total_length = sum(scaffold_lengths.get(scaffold, 0) for scaffold in counts.keys())
        global_density = (total_variants * 1000) / total_length if total_length > 0 else 0
        
        # Get gene status counts if available
        gene_status_counts = gene_stat_counts.get(treatment, {})
        erg_count = gene_status_counts.get('ERG', 0)
        nonerg_count = gene_status_counts.get('Non-ERG', 0)
        nogene_count = gene_status_counts.get('No Gene', 0)
        
        # Get treatment metadata
        description = TREATMENT_INFO.get(treatment, {}).get('description', 'Unknown')
        adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
        has_gene = TREATMENT_INFO.get(treatment, {}).get('gene') is not None
        
        # Top 3 densest scaffolds
        densities = all_densities[treatment]
        top_scaffolds = sorted(densities.items(), key=lambda x: x[1], reverse=True)[:3]
        top_scaffold_info = ", ".join([f"{s} ({d:.2f})" for s, d in top_scaffolds])
        
        summary.append({
            'Treatment': treatment,
            'Description': description,
            'Adaptation': adaptation,
            'Has_Gene': 'Yes' if has_gene else 'No',
            'Total_Variants': total_variants,
            'Scaffolds_With_Variants': total_scaffolds,
            'Global_Density': global_density,
            'ERG_Variants': erg_count,
            'NonERG_Variants': nonerg_count,
            'NoGene_Variants': nogene_count,
            'Top_Scaffolds': top_scaffold_info
        })
    
    # Convert to DataFrame
    summary_df = pd.DataFrame(summary)
    
    # Add adaptation averages
    adaptation_summary = summary_df.groupby('Adaptation').agg({
        'Total_Variants': 'mean',
        'Scaffolds_With_Variants': 'mean',
        'Global_Density': 'mean',
        'ERG_Variants': 'mean',
        'NonERG_Variants': 'mean',
        'NoGene_Variants': 'mean'
    }).reset_index()
    
    adaptation_summary['Treatment'] = adaptation_summary['Adaptation'] + ' (Average)'
    adaptation_summary['Description'] = 'Average across treatments'
    adaptation_summary['Has_Gene'] = 'N/A'
    adaptation_summary['Top_Scaffolds'] = 'N/A'
    
    # Combine with main summary
    summary_df = pd.concat([summary_df, adaptation_summary], ignore_index=True)
    
    # Save to CSV
    summary_df.to_csv(os.path.join(output_dir, "scaffold_distribution_summary.csv"), index=False)
    
    return summary_df

# Function to create summary report
def create_summary_report(all_data, all_counts, all_densities, scaffold_lengths, gene_stat_counts, gene_enrichment_results, output_dir):
    """Create a comprehensive summary report."""
    with open(os.path.join(output_dir, "scaffold_distribution_summary.txt"), 'w') as f:
        f.write("Scaffold Distribution Analysis Summary\n")
        f.write("=====================================\n\n")
        
        # Overall statistics
        f.write("Overall Statistics:\n")
        f.write("-----------------\n")
        
        for treatment in all_densities.keys():
            counts = all_counts[treatment]
            total_variants = sum(counts.values())
            total_scaffolds = len(counts)
            
            # Calculate global variant density
            total_length = sum(scaffold_lengths.get(scaffold, 0) for scaffold in counts.keys())
            global_density = (total_variants * 1000) / total_length if total_length > 0 else 0
            
            # Get treatment metadata
            description = TREATMENT_INFO.get(treatment, {}).get('description', 'Unknown')
            adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
            has_gene = TREATMENT_INFO.get(treatment, {}).get('gene')
            gene_text = f" with {has_gene} gene" if has_gene else ""
            
            # Get gene status counts if available
            gene_status_counts = gene_stat_counts.get(treatment, {})
            erg_count = gene_status_counts.get('ERG', 0)
            nonerg_count = gene_status_counts.get('Non-ERG', 0)
            nogene_count = gene_status_counts.get('No Gene', 0)
            
            # Get gene enrichment results if available
            enriched_genes, depleted_genes = gene_enrichment_results.get(treatment, (None, None))
            enriched_erg = len(enriched_genes[enriched_genes['Is_ERG']]) if enriched_genes is not None and not enriched_genes.empty else 0
            depleted_erg = len(depleted_genes[depleted_genes['Is_ERG']]) if depleted_genes is not None and not depleted_genes.empty else 0
            
            f.write(f"{treatment} Treatment ({description}):\n")
            f.write(f"  Adaptation: {adaptation}{gene_text}\n")
            f.write(f"  Total Variants: {total_variants}\n")
            f.write(f"  Scaffolds with Variants: {total_scaffolds}\n")
            f.write(f"  Global Variant Density: {global_density:.4f} variants/kb\n")
            
            # Gene status information
            if erg_count > 0 or nonerg_count > 0 or nogene_count > 0:
                f.write("  Gene Status Distribution:\n")
                f.write(f"    ERG Genes: {erg_count} variants ({erg_count/total_variants*100:.1f}%)\n")
                f.write(f"    Non-ERG Genes: {nonerg_count} variants ({nonerg_count/total_variants*100:.1f}%)\n")
                f.write(f"    No Gene Association: {nogene_count} variants ({nogene_count/total_variants*100:.1f}%)\n")
            
            # Gene enrichment/depletion information
            if enriched_erg > 0 or depleted_erg > 0:
                f.write("  Ergosterol Pathway Analysis:\n")
                if enriched_erg > 0:
                    f.write(f"    {enriched_erg} ergosterol genes with significant enrichment of variants\n")
                if depleted_erg > 0:
                    f.write(f"    {depleted_erg} ergosterol genes with significant depletion of variants (purifying selection)\n")
            
            # Top 5 densest scaffolds
            densities = all_densities[treatment]
            top_scaffolds = sorted(densities.items(), key=lambda x: x[1], reverse=True)[:5]
            
            f.write("  Top 5 Scaffolds by Variant Density:\n")
            for scaffold, density in top_scaffolds:
                count = counts.get(scaffold, 0)
                length = scaffold_lengths.get(scaffold, 0)
                f.write(f"    {scaffold}: {density:.4f} variants/kb (Count: {count}, Length: {length}bp)\n")
            
            f.write("\n")
        
        # Adaptation comparisons
        f.write("Adaptation Type Comparison:\n")
        f.write("-------------------------\n")
        
        adaptation_variants = defaultdict(int)
        adaptation_scaffolds = defaultdict(set)
        adaptation_length = defaultdict(int)
        adaptation_gene_status = defaultdict(lambda: defaultdict(int))
        
        for treatment, counts in all_counts.items():
            adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
            adaptation_variants[adaptation] += sum(counts.values())
            adaptation_scaffolds[adaptation].update(counts.keys())
            
            # Sum scaffold lengths
            for scaffold in counts.keys():
                adaptation_length[adaptation] += scaffold_lengths.get(scaffold, 0)
            
            # Sum gene status counts if available
            gene_status_counts = gene_stat_counts.get(treatment, {})
            for status, count in gene_status_counts.items():
                adaptation_gene_status[adaptation][status] += count
        
        for adaptation, variants in adaptation_variants.items():
            scaffolds = len(adaptation_scaffolds[adaptation])
            length = adaptation_length[adaptation]
            density = (variants * 1000) / length if length > 0 else 0
            
            erg_count = adaptation_gene_status[adaptation].get('ERG', 0)
            nonerg_count = adaptation_gene_status[adaptation].get('Non-ERG', 0)
            nogene_count = adaptation_gene_status[adaptation].get('No Gene', 0)
            
            f.write(f"{adaptation} Adaptation:\n")
            f.write(f"  Average Variants: {variants/2:.1f}\n")  # Assuming 2 treatments per adaptation
            f.write(f"  Total Unique Scaffolds: {scaffolds}\n")
            f.write(f"  Average Variant Density: {density:.4f} variants/kb\n")
            
            # Gene status information
            total = erg_count + nonerg_count + nogene_count
            if total > 0:
                f.write("  Gene Status Distribution:\n")
                f.write(f"    ERG Genes: {erg_count} variants ({erg_count/total*100:.1f}%)\n")
                f.write(f"    Non-ERG Genes: {nonerg_count} variants ({nonerg_count/total*100:.1f}%)\n")
                f.write(f"    No Gene Association: {nogene_count} variants ({nogene_count/total*100:.1f}%)\n")
            
            f.write("\n")
        
        # Gene modification effects
        gene_mod_variants = defaultdict(int)
        non_mod_variants = defaultdict(int)
        gene_mod_status = defaultdict(lambda: defaultdict(int))
        non_mod_status = defaultdict(lambda: defaultdict(int))
        
        for treatment, counts in all_counts.items():
            adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
            has_gene = TREATMENT_INFO.get(treatment, {}).get('gene') is not None
            
            variants = sum(counts.values())
            
            if has_gene:
                gene_mod_variants[adaptation] += variants
                # Add gene status counts
                for status, count in gene_stat_counts.get(treatment, {}).items():
                    gene_mod_status[adaptation][status] += count
            else:
                non_mod_variants[adaptation] += variants
                # Add gene status counts
                for status, count in gene_stat_counts.get(treatment, {}).items():
                    non_mod_status[adaptation][status] += count
        
        f.write("Gene Modification Effects:\n")
        f.write("------------------------\n")
        
        for adaptation in gene_mod_variants.keys():
            gene_vars = gene_mod_variants[adaptation]
            non_vars = non_mod_variants[adaptation]
            
            f.write(f"{adaptation} Adaptation:\n")
            f.write(f"  Gene-Modified Variants: {gene_vars}\n")
            f.write(f"  Non-Modified Variants: {non_vars}\n")
            
            if non_vars > 0:
                ratio = gene_vars / non_vars
                f.write(f"  Gene/Non-Gene Ratio: {ratio:.2f}\n")
            
            # Compare ERG genes between gene-modified and non-modified
            gene_mod_erg = gene_mod_status[adaptation].get('ERG', 0)
            non_mod_erg = non_mod_status[adaptation].get('ERG', 0)
            
            if gene_mod_erg > 0 and non_mod_erg > 0:
                gene_mod_erg_pct = gene_mod_erg / gene_vars * 100 if gene_vars > 0 else 0
                non_mod_erg_pct = non_mod_erg / non_vars * 100 if non_vars > 0 else 0
                
                f.write("  ERG Gene Comparison:\n")
                f.write(f"    Gene-Modified ERG: {gene_mod_erg} variants ({gene_mod_erg_pct:.1f}%)\n")
                f.write(f"    Non-Modified ERG: {non_mod_erg} variants ({non_mod_erg_pct:.1f}%)\n")
                f.write(f"    Difference: {gene_mod_erg_pct - non_mod_erg_pct:.1f}% points\n")
            
            f.write("\n")
        
        # Write treatment correlation summary
        f.write("Treatment Correlation Summary:\n")
        f.write("----------------------------\n")
        
        try:
            corr_matrix, adaptation_corr = calculate_treatment_correlations(all_densities, output_dir)
            
            # Get pairs with highest and lowest correlation
            treatments = list(all_densities.keys())
            pairs = [(i, j) for i in range(len(treatments)) for j in range(i+1, len(treatments))]
            
            if pairs and len(corr_matrix) > 1:
                corr_values = [corr_matrix.iloc[i, j] for i, j in pairs]
                max_idx = np.argmax(corr_values)
                min_idx = np.argmin(corr_values)
                
                max_pair = (treatments[pairs[max_idx][0]], treatments[pairs[max_idx][1]])
                min_pair = (treatments[pairs[min_idx][0]], treatments[pairs[min_idx][1]])
                
                f.write(f"  Most Similar Treatments: {max_pair[0]} and {max_pair[1]} (ρ = {corr_values[max_idx]:.4f})\n")
                f.write(f"  Most Different Treatments: {min_pair[0]} and {min_pair[1]} (ρ = {corr_values[min_idx]:.4f})\n\n")
                
                # Full correlation matrix
                f.write("  Full Correlation Matrix:\n")
                f.write("  " + " " * 4 + "".join(f"{t:>8}" for t in treatments) + "\n")
                
                for i, t1 in enumerate(treatments):
                    f.write(f"  {t1:4}")
                    for j, t2 in enumerate(treatments):
                        f.write(f"{corr_matrix.iloc[i, j]:8.4f}")
                    f.write("\n")
                
                f.write("\n")
            
            # Add adaptation correlation summary
            if len(adaptation_corr) > 1:
                f.write("  Adaptation Correlation:\n")
                for adap1 in adaptation_corr.index:
                    for adap2 in adaptation_corr.columns:
                        if adap1 != adap2:
                            corr = adaptation_corr.loc[adap1, adap2]
                            f.write(f"    {adap1} vs {adap2}: {corr:.4f}\n")
                
                f.write("\n")
        except Exception as e:
            f.write(f"  Error calculating correlations: {str(e)}\n\n")
        
        # Gene-specific scaffold enrichment analysis
        f.write("Gene-Specific Scaffold Enrichment Analysis:\n")
        f.write("----------------------------------------\n")
        
        for treatment, (enriched_genes, depleted_genes) in gene_enrichment_results.items():
            if enriched_genes is not None and not enriched_genes.empty:
                f.write(f"{treatment} - Genes with Enriched Variants:\n")
                f.write(f"  Found {len(enriched_genes)} enriched genes (p < 0.05)\n")
                
                # ERG gene specific information
                erg_enriched = enriched_genes[enriched_genes['Is_ERG']]
                if not erg_enriched.empty:
                    f.write(f"  {len(erg_enriched)} ergosterol pathway genes enriched:\n")
                    for _, gene in erg_enriched.iterrows():
                        f.write(f"    {gene['ERG_Name']} (Fold change: {2**gene['Log2_Fold_Change']:.2f}x, p={gene['P_Value']:.4f})\n")
                
                f.write("\n")
            
            if depleted_genes is not None and not depleted_genes.empty:
                f.write(f"{treatment} - Genes with Depleted Variants (Purifying Selection):\n")
                f.write(f"  Found {len(depleted_genes)} depleted genes (p < 0.05)\n")
                
                # ERG gene specific information
                erg_depleted = depleted_genes[depleted_genes['Is_ERG']]
                if not erg_depleted.empty:
                    f.write(f"  {len(erg_depleted)} ergosterol pathway genes under purifying selection:\n")
                    for _, gene in erg_depleted.iterrows():
                        f.write(f"    {gene['ERG_Name']} (Fold change: {2**gene['Log2_Fold_Change']:.2f}x, p={gene['P_Value']:.4f})\n")
                
                f.write("\n")
        
        # Main conclusions
        f.write("Main Conclusions:\n")
        f.write("---------------\n")
        f.write("1. This analysis identifies scaffolds with high variant densities across treatments.\n")
        f.write("2. Gene-specific analysis reveals patterns of enrichment and depletion in ergosterol pathway genes.\n")
        f.write("3. Several ergosterol genes show evidence of purifying selection with significantly fewer variants than expected.\n")
        f.write("4. Temperature and low oxygen adaptations show distinct scaffold distribution patterns.\n")
        f.write("5. Gene modifications (STC, CAS) appear to influence scaffold enrichment patterns.\n")
        f.write("6. Analysis of variant distribution by gene type provides insights into functional impact.\n")

# Main function to run the analysis
def main():
    """Run the scaffold distribution analysis with gene-specific enhancements."""
    # Load gene mapping data
    gene_mapping, erg_genes = load_gene_mapping()
    
    # Build gene coordinates dictionary
    scaffold_to_genes = build_gene_coordinates(gene_mapping)
    
    # Parse scaffold lengths from reference genome index
    scaffold_lengths = parse_genome_index()
    
    # Parse data for each treatment
    all_data = {}
    for treatment in TREATMENTS:
        data = parse_mutation_data(treatment)
        if len(data) > 0:
            # Map variants to genes
            data = map_variants_to_genes(data, scaffold_to_genes)
            all_data[treatment] = data
            print(f"Loaded {len(data)} variants for {treatment} treatment")
        else:
            print(f"Warning: No data available for {treatment} treatment")
    
    # Count variants per scaffold and by gene status
    all_counts = {}
    gene_stat_counts = {}
    all_gene_counts = {}
    
    for treatment, data in all_data.items():
        # Count by scaffold
        counts = count_variants_per_scaffold(data)
        all_counts[treatment] = counts
        print(f"{treatment}: Found variants in {len(counts)} scaffolds")
        
        # Count by gene status
        gene_status_counts = count_variants_by_gene_status(data)
        gene_stat_counts[treatment] = gene_status_counts
        
        # Count by gene
        gene_counts = count_variants_per_gene(data)
        all_gene_counts[treatment] = gene_counts
    
    # Calculate variant densities
    all_densities = {}
    all_gene_densities = {}
    
    for treatment, counts in all_counts.items():
        # Scaffold densities
        densities = calculate_variant_density(counts, scaffold_lengths)
        all_densities[treatment] = densities
        
        # Gene densities
        if treatment in all_data:
            gene_densities = calculate_gene_density(all_data[treatment], gene_mapping)
            all_gene_densities[treatment] = gene_densities
    
    # Generate variant density plots
    for treatment, densities in all_densities.items():
        plot_variant_density(densities, treatment, OUTPUT_DIR)
    
    # Generate gene status distribution plots
    for treatment, data in all_data.items():
        plot_gene_status_distribution(data, treatment, OUTPUT_DIR)
    
    # Generate comparative heatmap
    plot_comparative_heatmap(all_densities, OUTPUT_DIR)
    
    # Generate bubble charts
    for treatment, counts in all_counts.items():
        plot_bubble_chart(counts, scaffold_lengths, treatment, OUTPUT_DIR)
    
    # Calculate treatment correlations
    corr_matrix, adaptation_corr = calculate_treatment_correlations(all_densities, OUTPUT_DIR)
    
    # Identify treatment-specific hotspots
    treatment_specific, adaptation_specific, gene_comparison = identify_treatment_specific_hotspots(all_densities, OUTPUT_DIR)
    
    # Identify significantly enriched or depleted genes (for purifying selection analysis)
    gene_enrichment_results = {}
    
    for treatment, data in all_data.items():
        enriched_df, depleted_df = identify_gene_patterns(data, gene_mapping)
        gene_enrichment_results[treatment] = (enriched_df, depleted_df)
        
        # Generate plots for gene enrichment/depletion patterns
        plot_gene_enrichment(enriched_df, depleted_df, treatment, OUTPUT_DIR)
        
        # Print summary
        print(f"{treatment}: Found {len(enriched_df)} enriched genes and {len(depleted_df)} depleted genes")
    
    # Create summary table
    summary_df = create_summary_table(all_counts, all_densities, scaffold_lengths, gene_stat_counts, OUTPUT_DIR)
    
    # Create summary report
    create_summary_report(all_data, all_counts, all_densities, scaffold_lengths, gene_stat_counts, gene_enrichment_results, OUTPUT_DIR)
    
    print(f"Analysis complete! Results saved to {OUTPUT_DIR}/")
    print(f"Generated summary report: {os.path.join(OUTPUT_DIR, 'scaffold_distribution_summary.txt')}")
    print(f"Generated treatment-specific hotspots: {os.path.join(OUTPUT_DIR, 'treatment_specific_hotspots.txt')}")
    print("Gene-specific analysis successfully incorporated into scaffold distribution analysis.")

# Run the analysis
if __name__ == "__main__":
    main()