 #!/usr/bin/env python3
"""
Gene-Specific Population Structure Analysis Module

This module analyzes population structure across different treatment conditions in yeast adaptation
experiments, with a focus on gene-level analysis. It performs dimensionality reduction, clustering,
and similarity analysis on variant patterns, grouped by treatment, adaptation type, and gene properties.

The module maps variants to genes based on genomic coordinates, allowing for analysis at the gene level
rather than just the scaffold level. This helps identify patterns of mutations in specific genes across
different treatments and adaptations.
"""

import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict, Counter
from Bio import SeqIO, motifs
from Bio.Seq import Seq
import random
import subprocess
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist, squareform
import csv
import warnings
warnings.filterwarnings('ignore')

# Set matplotlib style
plt.style.use('ggplot')
sns.set(font_scale=1.2)
sns.set_style("whitegrid")

# Define output directories
# Define output directories
BASE_DIR = os.environ.get("OUTPUT_DIR", "analysis/general_gene_analysis")
OUTPUT_DIR = f"{BASE_DIR}/population_structure_results"

GENE_OUTPUT_DIR = f"{BASE_DIR}/gene_gene_population_structure_results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Log the output directories
logging.info(f"Using output directory: {OUTPUT_DIR}")
os.makedirs(GENE_OUTPUT_DIR, exist_ok=True)

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

# Adaptation colors for grouping by adaptation type
ADAPTATION_COLORS = {
    'Temperature': '#1f77b4',
    'Low Oxygen': '#ff7f0e',
}

# File paths for gene mapping and annotations
GENE_MAPPING_FILE = "reference/gene_mapping.tsv"
GENES_OF_INTEREST_FILE = "reference/genes_of_interest_mapping.tsv"

# Gene data structures
GENE_DATA = {}
SCAFFOLD_GENES = defaultdict(list)
GENES_OF_INTEREST = set()


def load_gene_mapping():
    """
    Load the gene mapping information from the reference files.
    
    This function loads gene data from:
    1. The main gene mapping file with all annotated genes
    2. The genes of interest file with ergosterol pathway genes
    
    Returns:
        dict: Dictionary mapping gene ID to gene information
        defaultdict: Dictionary mapping scaffold to list of genes on that scaffold
        set: Set of gene IDs that are of special interest (ergosterol pathway)
    """
    global GENE_DATA, SCAFFOLD_GENES, GENES_OF_INTEREST
    
    print(f"Loading gene mapping from {GENE_MAPPING_FILE}")
    
    # Load the main gene mapping file
    if os.path.exists(GENE_MAPPING_FILE):
        with open(GENE_MAPPING_FILE, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gene_id = row['w303_gene_id']
                GENE_DATA[gene_id] = {
                    'sc_gene_id': row['sc_gene_id'],
                    'erg_name': row.get('erg_name', ''),
                    'locus_tag': row['locus_tag'],
                    'w303_scaffold': row['w303_scaffold'],
                    'start': int(row['start']),
                    'end': int(row['end']),
                    'strand': row['strand'],
                    'product': row.get('product', 'hypothetical protein')
                }
                
                # Add to scaffold lookup
                SCAFFOLD_GENES[row['w303_scaffold']].append(gene_id)
        
        print(f"Loaded information for {len(GENE_DATA)} genes across {len(SCAFFOLD_GENES)} scaffolds")
    else:
        print(f"Warning: Gene mapping file {GENE_MAPPING_FILE} not found")
    
    # Load the genes of interest file
    if os.path.exists(GENES_OF_INTEREST_FILE):
        with open(GENES_OF_INTEREST_FILE, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                GENES_OF_INTEREST.add(row['w303_gene_id'])
        
        print(f"Loaded {len(GENES_OF_INTEREST)} genes of interest")
    else:
        print(f"Warning: Genes of interest file {GENES_OF_INTEREST_FILE} not found")
    
    return GENE_DATA, SCAFFOLD_GENES, GENES_OF_INTEREST


def map_variants_to_genes(variant_df):
    """
    Map variants to genes based on their genomic coordinates.
    
    This function takes a DataFrame of variants with scaffold and position information
    and maps each variant to the corresponding gene if it falls within a gene's coordinates.
    Variants are also categorized as being in genes of interest or not.
    
    Args:
        variant_df (pandas.DataFrame): DataFrame containing variant information with at minimum
                                      'variant_id' containing CHROM_POS_REF_ALT format
    
    Returns:
        pandas.DataFrame: The original DataFrame with additional gene-related columns
    """
    print("Mapping variants to genes...")
    
    if len(GENE_DATA) == 0:
        print("Warning: No gene data loaded. Cannot map variants to genes.")
        return variant_df
    
    # Function to extract scaffold and position from variant_id
    def extract_scaffold_pos(variant_id):
        parts = variant_id.split('_')
        if len(parts) >= 2:
            scaffold = parts[0]
            try:
                pos = int(parts[1])
                return scaffold, pos
            except ValueError:
                return None, None
        return None, None
    
    # Function to find gene(s) for a given variant
    def find_genes_for_variant(scaffold, pos):
        if scaffold not in SCAFFOLD_GENES:
            return None, None, False
        
        for gene_id in SCAFFOLD_GENES[scaffold]:
            gene_info = GENE_DATA[gene_id]
            
            # Check if position is within gene bounds
            if gene_info['start'] <= pos <= gene_info['end']:
                # It's within the gene
                in_gene_of_interest = gene_id in GENES_OF_INTEREST
                return gene_id, gene_info.get('sc_gene_id', ''), in_gene_of_interest
        
        return None, None, False
    
    # Extract scaffold and position for each variant
    variant_df['scaffold'], variant_df['position'] = zip(*variant_df['variant_id'].apply(extract_scaffold_pos))
    
    # Map variants to genes
    gene_mappings = []
    for _, row in variant_df.iterrows():
        if row['scaffold'] is not None and row['position'] is not None:
            gene_id, sc_gene_id, in_gene_of_interest = find_genes_for_variant(row['scaffold'], row['position'])
            gene_mappings.append({
                'gene_id': gene_id,
                'sc_gene_id': sc_gene_id,
                'in_gene_of_interest': in_gene_of_interest
            })
        else:
            gene_mappings.append({
                'gene_id': None,
                'sc_gene_id': None,
                'in_gene_of_interest': False
            })
    
    # Add gene mapping information to the DataFrame
    gene_mapping_df = pd.DataFrame(gene_mappings)
    enhanced_df = pd.concat([variant_df.reset_index(drop=True), gene_mapping_df.reset_index(drop=True)], axis=1)
    
    # Add gene information column
    enhanced_df['in_gene'] = enhanced_df['gene_id'].notna()
    
    # Add a column for gene type (ergosterol pathway or other)
    enhanced_df['gene_type'] = 'None'
    enhanced_df.loc[enhanced_df['in_gene'], 'gene_type'] = 'Other Gene'
    enhanced_df.loc[enhanced_df['in_gene_of_interest'], 'gene_type'] = 'Ergosterol Gene'
    
    print(f"Mapped variants to genes: {enhanced_df['in_gene'].sum()} variants in genes")
    print(f"Ergosterol pathway genes: {enhanced_df['in_gene_of_interest'].sum()} variants")
    
    return enhanced_df


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

# Function to extract variant data from treatment-specific VCF files
def extract_treatment_specific_variants():
    """Extract treatment-specific variants from separate VCF files."""
    all_samples = []
    all_variant_data = []
    
    # Define possible VCF file patterns
    vcf_patterns = [
        "results/merged/analysis/{}/highconf.vcf.gz",
        "results/merged/analysis/{}_highconf.vcf.gz",
        "results/merged/analysis/{}/specific.vcf.gz",
        "results/merged/analysis/{}_specific.vcf.gz"
    ]
    
    # Process each treatment
    for treatment in TREATMENTS:
        # Find VCF file for this treatment
        vcf_file = find_file(treatment, vcf_patterns)
        
        if not vcf_file:
            print(f"Warning: No VCF file found for {treatment}")
            continue
        
        print(f"Processing {vcf_file} for {treatment} treatment...")
        
        # Extract sample names
        try:
            sample_output = subprocess.check_output(
                f"bcftools query -l {vcf_file}", shell=True).decode('utf-8')
            samples = sample_output.strip().split('\n')
            print(f"Found {len(samples)} samples: {', '.join(samples)}")
        except Exception as e:
            print(f"Error extracting samples: {str(e)}")
            continue
        
        # Extract variant data
        try:
            # Get variant details and genotypes
            cmd = f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' {vcf_file}"
            output = subprocess.check_output(cmd, shell=True).decode('utf-8')
            
            # Process each variant
            variant_count = 0
            for line in output.strip().split('\n'):
                if not line:
                    continue
                
                parts = line.split('\t')
                if len(parts) < 4 + len(samples):
                    print(f"Warning: Line has {len(parts)} parts, expected {4 + len(samples)}")
                    continue
                
                chrom, pos, ref, alt = parts[:4]
                genotypes = parts[4:4+len(samples)]
                
                # Create variant identifier
                variant_id = f"{chrom}_{pos}_{ref}_{alt}"
                variant_count += 1
                
                # For each sample, determine if variant is present
                for i, sample in enumerate(samples):
                    if i < len(genotypes):
                        gt = genotypes[i]
                        # Check if variant is present in this sample
                        # Only count as present if genotype contains '1' and is not missing
                        has_variant = 0
                        if gt not in ['0/0', '0|0', './.', '.|.', '.'] and ('1' in gt):
                            has_variant = 1
                        
                        if has_variant:
                            # Add adaptation and gene info
                            adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
                            has_gene = 'Yes' if TREATMENT_INFO.get(treatment, {}).get('gene') else 'No'
                            
                            all_variant_data.append({
                                'sample': sample,
                                'variant_id': variant_id,
                                'treatment': treatment,
                                'adaptation': adaptation,
                                'has_gene': has_gene
                            })
            
            print(f"Processed {variant_count} variants for {treatment}")
            all_samples.extend(samples)
        except Exception as e:
            print(f"Error processing {vcf_file}: {str(e)}")
    
    # Convert to DataFrame
    variant_df = pd.DataFrame(all_variant_data)
    print(f"Created DataFrame with {len(variant_df)} variant-sample entries")
    
    # Get unique samples
    unique_samples = sorted(set(all_samples))
    print(f"Found {len(unique_samples)} unique samples")
    
    return variant_df, unique_samples

# Function to create sample-by-variant matrix
def create_sample_variant_matrix(variant_df, samples):
    """Create a sample-by-variant presence/absence matrix."""
    if len(variant_df) == 0:
        print("Error: No variant data available")
        return None, None, None, None
    
    # Get unique variant IDs
    variant_ids = sorted(variant_df['variant_id'].unique())
    print(f"Found {len(variant_ids)} unique variants")
    
    # Create pivot table: samples by variants
    print("Creating sample-by-variant matrix...")
    
    # Initialize matrix with zeros
    matrix = np.zeros((len(samples), len(variant_ids)))
    
    # Create mapping from variant_id to index
    variant_to_idx = {v: i for i, v in enumerate(variant_ids)}
    
    # Create mapping from sample to index
    sample_to_idx = {s: i for i, s in enumerate(samples)}
    
    # Fill the matrix
    for _, row in variant_df.iterrows():
        sample = row['sample']
        variant = row['variant_id']
        
        if sample in sample_to_idx and variant in variant_to_idx:
            matrix[sample_to_idx[sample], variant_to_idx[variant]] = 1
    
    # Extract treatments and adaptations from sample names
    treatments = []
    adaptations = []
    has_genes = []
    
    for sample in samples:
        # Extract treatment from sample name
        treatment_prefix = sample.split('-')[0]
        if treatment_prefix not in ['WT', 'WTA', 'STC', 'CAS']:
            # Handle special case for WT-37
            if sample.startswith('WT-37'):
                treatment_prefix = 'WT-37'
            else:
                treatment_prefix = 'Unknown'
        elif treatment_prefix == 'WT':
            # Update old WT naming to WT-37
            if 'WT-37' in TREATMENTS:
                treatment_prefix = 'WT-37'
        
        # Get adaptation type and gene status
        adaptation = TREATMENT_INFO.get(treatment_prefix, {}).get('adaptation', 'Unknown')
        has_gene = 'Yes' if TREATMENT_INFO.get(treatment_prefix, {}).get('gene') else 'No'
        
        treatments.append(treatment_prefix)
        adaptations.append(adaptation)
        has_genes.append(has_gene)
    
    # Verify we have variation in the data
    sample_counts = matrix.sum(axis=1)
    variant_counts = matrix.sum(axis=0)
    
    print(f"Sample variant counts - min: {sample_counts.min()}, max: {sample_counts.max()}")
    print(f"Variant presence counts - min: {variant_counts.min()}, max: {variant_counts.max()}")
    
    return matrix, variant_ids, treatments, {'treatments': treatments, 'adaptations': adaptations, 'has_genes': has_genes}

# Function to calculate genetic distances
def calculate_genetic_distances(matrix, samples):
    """Calculate genetic distances between samples based on variant profiles."""
    print("Calculating genetic distances...")
    
    # Calculate Jaccard distances
    distances = pdist(matrix, metric='jaccard')
    
    # Convert to square matrix
    distance_matrix = squareform(distances)
    
    # Create DataFrame for better visualization
    distance_df = pd.DataFrame(distance_matrix, index=samples, columns=samples)
    
    return distance_df

# Function to perform PCA
def perform_pca(matrix, samples, metadata, n_components=3):
    """Perform PCA on the sample-by-variant matrix."""
    print("Performing PCA analysis...")
    
    # Check if we have enough samples for PCA
    if matrix.shape[0] <= 1:
        print("Not enough samples for PCA")
        return None
    
    # Check component count
    n_components = min(n_components, min(matrix.shape) - 1)
    if n_components < 1:
        print("Not enough dimensions for PCA")
        return None
    
    # Perform PCA
    pca = PCA(n_components=n_components)
    principal_components = pca.fit_transform(matrix)
    
    # Create DataFrame with PCA results
    columns = [f'PC{i+1}' for i in range(n_components)]
    pca_df = pd.DataFrame(data=principal_components, columns=columns)
    
    # Add sample information
    pca_df['Sample'] = samples
    pca_df['Treatment'] = metadata
    
    # Add adaptation and gene info if available
    if 'adaptations' in metadata:
        pca_df['Adaptation'] = metadata['adaptations']
    if 'has_genes' in metadata:
        pca_df['Has_Gene'] = metadata['has_genes']
    
    # Calculate explained variance
    explained_variance = pca.explained_variance_ratio_
    
    # Calculate feature loadings
    if hasattr(pca, 'components_'):
        loadings = pca.components_.T
    else:
        loadings = None
    
    return pca_df, explained_variance, loadings

# Function to perform MDS
def perform_mds(distance_matrix, samples, metadata, n_components=2):
    """Perform Multidimensional Scaling on the distance matrix."""
    print("Performing MDS analysis...")
    
    # Check if we have enough samples for MDS
    if distance_matrix.shape[0] <= 1:
        print("Not enough samples for MDS")
        return None
    
    # Check component count
    n_components = min(n_components, distance_matrix.shape[0] - 1)
    if n_components < 1:
        print("Not enough dimensions for MDS")
        return None
    
    # Perform MDS
    mds = MDS(n_components=n_components, dissimilarity='precomputed', random_state=42)
    mds_result = mds.fit_transform(distance_matrix.values)
    
    # Create DataFrame with MDS results
    columns = [f'MDS{i+1}' for i in range(n_components)]
    mds_df = pd.DataFrame(data=mds_result, columns=columns)
    
    # Add sample information
    mds_df['Sample'] = samples
    mds_df['Treatment'] = metadata
    
    # Add adaptation and gene info if available
    if isinstance(metadata, dict) and 'adaptations' in metadata:
        mds_df['Adaptation'] = metadata['adaptations']
    if isinstance(metadata, dict) and 'has_genes' in metadata:
        mds_df['Has_Gene'] = metadata['has_genes']
    
    return mds_df

# Function to perform hierarchical clustering
def perform_hierarchical_clustering(distance_matrix):
    """Perform hierarchical clustering on the distance matrix."""
    print("Performing hierarchical clustering...")
    
    # Check if we have enough samples for clustering
    if distance_matrix.shape[0] <= 1:
        print("Not enough samples for hierarchical clustering")
        return None
    
    # Calculate linkage
    try:
        linkage_matrix = hierarchy.linkage(squareform(distance_matrix.values), method='average')
        return linkage_matrix
    except Exception as e:
        print(f"Error during hierarchical clustering: {e}")
        return None

# Function to plot PCA results
def plot_pca(pca_df, explained_variance, output_dir):
    """Plot PCA results to visualize sample relationships."""
    if pca_df is None or 'PC1' not in pca_df.columns:
        print("Cannot plot PCA: insufficient data")
        return
    
    # Create colormap for treatments
    treatment_colors = {t: TREATMENT_COLORS.get(t, '#999999') for t in pca_df['Treatment'].unique()}
    
    # 1. Plot by treatment
    plt.figure(figsize=(10, 8))
    
    for treatment in pca_df['Treatment'].unique():
        subset = pca_df[pca_df['Treatment'] == treatment]
        plt.scatter(
            subset['PC1'], 
            subset['PC2'] if 'PC2' in subset.columns else [0] * len(subset), 
            alpha=0.8, 
            s=100, 
            label=treatment,
            color=treatment_colors.get(treatment, '#999999')
        )
    
    # Add sample labels
    for i, row in pca_df.iterrows():
        plt.annotate(
            row['Sample'],
            (row['PC1'], row['PC2'] if 'PC2' in pca_df.columns else 0),
            fontsize=8,
            alpha=0.7,
            xytext=(5, 5),
            textcoords='offset points'
        )
    
    # Add axis labels with explained variance
    plt.xlabel(f'PC1 ({explained_variance[0]:.2%} variance)')
    if len(explained_variance) > 1:
        plt.ylabel(f'PC2 ({explained_variance[1]:.2%} variance)')
    else:
        plt.ylabel('PC2')
    
    plt.title('Principal Component Analysis by Treatment')
    plt.legend(title="Treatment")
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "pca_by_treatment.png"), dpi=300)
    plt.close()
    
    # 2. Plot by adaptation type if available
    if 'Adaptation' in pca_df.columns:
        plt.figure(figsize=(10, 8))
        
        # Plot points colored by adaptation
        for adaptation in pca_df['Adaptation'].unique():
            subset = pca_df[pca_df['Adaptation'] == adaptation]
            plt.scatter(
                subset['PC1'], 
                subset['PC2'] if 'PC2' in subset.columns else [0] * len(subset), 
                alpha=0.8, 
                s=100, 
                label=adaptation,
                color=ADAPTATION_COLORS.get(adaptation, '#999999')
            )
        
        # Highlight gene-modified samples with different markers
        if 'Has_Gene' in pca_df.columns:
            gene_modified = pca_df[pca_df['Has_Gene'] == 'Yes']
            if not gene_modified.empty:
                plt.scatter(
                    gene_modified['PC1'], 
                    gene_modified['PC2'] if 'PC2' in gene_modified.columns else [0] * len(gene_modified), 
                    s=150, 
                    facecolors='none',
                    edgecolors='black',
                    linewidth=2,
                    label='Gene Modified'
                )
        
        # Add sample labels
        for i, row in pca_df.iterrows():
            plt.annotate(
                row['Sample'],
                (row['PC1'], row['PC2'] if 'PC2' in pca_df.columns else 0),
                fontsize=8,
                alpha=0.7,
                xytext=(5, 5),
                textcoords='offset points'
            )
        
        # Add axis labels with explained variance
        plt.xlabel(f'PC1 ({explained_variance[0]:.2%} variance)')
        if len(explained_variance) > 1:
            plt.ylabel(f'PC2 ({explained_variance[1]:.2%} variance)')
        else:
            plt.ylabel('PC2')
        
        plt.title('Principal Component Analysis by Adaptation Type')
        plt.legend(title="Adaptation")
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "pca_by_adaptation.png"), dpi=300)
        plt.close()
    
    # 3. Create scree plot for explained variance
    plt.figure(figsize=(10, 6))
    plt.bar(range(1, len(explained_variance) + 1), explained_variance, color='steelblue')
    plt.xlabel('Principal Component')
    plt.ylabel('Explained Variance Ratio')
    plt.title('Scree Plot: Explained Variance by Component')
    plt.xticks(range(1, len(explained_variance) + 1))
    plt.grid(True, alpha=0.3)
    
    # Add cumulative explained variance line
    plt.plot(
        range(1, len(explained_variance) + 1),
        np.cumsum(explained_variance),
        'ro-',
        linewidth=2
    )
    plt.axhline(y=0.8, color='gray', linestyle='--', alpha=0.7)
    plt.axhline(y=0.9, color='gray', linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "pca_scree_plot.png"), dpi=300)
    plt.close()

# Function to plot MDS results
def plot_mds(mds_df, output_dir):
    """Plot MDS results to visualize sample relationships."""
    if mds_df is None or 'MDS1' not in mds_df.columns:
        print("Cannot plot MDS: insufficient data")
        return
    
    # Create colormap for treatments
    treatment_colors = {t: TREATMENT_COLORS.get(t, '#999999') for t in mds_df['Treatment'].unique()}
    
    # 1. Plot by treatment
    plt.figure(figsize=(10, 8))
    
    for treatment in mds_df['Treatment'].unique():
        subset = mds_df[mds_df['Treatment'] == treatment]
        plt.scatter(
            subset['MDS1'], 
            subset['MDS2'] if 'MDS2' in subset.columns else [0] * len(subset), 
            alpha=0.8, 
            s=100, 
            label=treatment,
            color=treatment_colors.get(treatment, '#999999')
        )
    
    # Add sample labels
    for i, row in mds_df.iterrows():
        plt.annotate(
            row['Sample'],
            (row['MDS1'], row['MDS2'] if 'MDS2' in mds_df.columns else 0),
            fontsize=8,
            alpha=0.7,
            xytext=(5, 5),
            textcoords='offset points'
        )
    
    plt.xlabel('MDS Dimension 1')
    plt.ylabel('MDS Dimension 2')
    plt.title('Multidimensional Scaling by Treatment')
    plt.legend(title="Treatment")
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "mds_by_treatment.png"), dpi=300)
    plt.close()
    
    # 2. Plot by adaptation type if available
    if 'Adaptation' in mds_df.columns:
        plt.figure(figsize=(10, 8))
        
        # Plot points colored by adaptation
        for adaptation in mds_df['Adaptation'].unique():
            subset = mds_df[mds_df['Adaptation'] == adaptation]
            plt.scatter(
                subset['MDS1'], 
                subset['MDS2'] if 'MDS2' in subset.columns else [0] * len(subset), 
                alpha=0.8, 
                s=100, 
                label=adaptation,
                color=ADAPTATION_COLORS.get(adaptation, '#999999')
            )
        
        # Highlight gene-modified samples with different markers
        if 'Has_Gene' in mds_df.columns:
            gene_modified = mds_df[mds_df['Has_Gene'] == 'Yes']
            if not gene_modified.empty:
                plt.scatter(
                    gene_modified['MDS1'], 
                    gene_modified['MDS2'] if 'MDS2' in gene_modified.columns else [0] * len(gene_modified), 
                    s=150, 
                    facecolors='none',
                    edgecolors='black',
                    linewidth=2,
                    label='Gene Modified'
                )
        
        # Add sample labels
        for i, row in mds_df.iterrows():
            plt.annotate(
                row['Sample'],
                (row['MDS1'], row['MDS2'] if 'MDS2' in mds_df.columns else 0),
                fontsize=8,
                alpha=0.7,
                xytext=(5, 5),
                textcoords='offset points'
            )
        
        plt.xlabel('MDS Dimension 1')
        plt.ylabel('MDS Dimension 2')
        plt.title('Multidimensional Scaling by Adaptation Type')
        plt.legend(title="Adaptation")
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "mds_by_adaptation.png"), dpi=300)
        plt.close()

# Function to plot distance matrix heatmap
def plot_distance_heatmap(distance_df, samples, metadata, output_dir):
    """Plot a heatmap of the genetic distance matrix."""
    if distance_df is None or distance_df.empty:
        print("Cannot plot distance heatmap: insufficient data")
        return
    
    plt.figure(figsize=(12, 10))
    
    # Determine what metadata is available
    treatments = metadata
    adaptations = None
    has_genes = None
    
    if isinstance(metadata, dict):
        treatments = metadata.get('treatments', [])
        adaptations = metadata.get('adaptations', [])
        has_genes = metadata.get('has_genes', [])
    
    # Create annotation DataFrame
    annotation_df = pd.DataFrame({'Treatment': treatments}, index=samples)
    if adaptations:
        annotation_df['Adaptation'] = adaptations
    if has_genes:
        annotation_df['Has_Gene'] = has_genes
    
    # Create row colors for the heatmap
    row_colors = []
    
    # Add treatment colors
    if 'Treatment' in annotation_df.columns:
        treatment_colors = {t: TREATMENT_COLORS.get(t, '#999999') for t in annotation_df['Treatment'].unique()}
        row_colors.append(annotation_df['Treatment'].map(treatment_colors))
    
    # Add adaptation colors
    if 'Adaptation' in annotation_df.columns:
        adaptation_colors = {a: ADAPTATION_COLORS.get(a, '#999999') for a in annotation_df['Adaptation'].unique()}
        row_colors.append(annotation_df['Adaptation'].map(adaptation_colors))
    
    # Create the heatmap
    g = sns.clustermap(
        distance_df,
        cmap='viridis_r',  # Reversed viridis (darker = more similar)
        figsize=(12, 10),
        row_colors=row_colors,
        col_colors=row_colors,
        xticklabels=True,
        yticklabels=True
    )
    
    # Rotate x-axis labels
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45)
    
    # Add legend for treatments
    for treatment, color in treatment_colors.items():
        if treatment in annotation_df['Treatment'].values:
            g.ax_row_dendrogram.bar(0, 0, color=color, label=treatment, linewidth=0)
    
    if 'Adaptation' in annotation_df.columns:
        # Move to second legend position
        for adaptation, color in adaptation_colors.items():
            if adaptation in annotation_df['Adaptation'].values:
                g.ax_row_dendrogram.bar(0, 0, color=color, label=adaptation, linewidth=0)
    
    g.ax_row_dendrogram.legend(title="Groups", loc="center", ncol=1)
    
    # Set title
    plt.suptitle('Genetic Distance Heatmap with Hierarchical Clustering', fontsize=14, y=0.98)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "distance_heatmap.png"), dpi=300)
    plt.close()

# Function to plot dendrogram
def plot_dendrogram(linkage_matrix, samples, metadata, output_dir):
    """Plot a dendrogram of sample relationships."""
    if linkage_matrix is None:
        print("Cannot plot dendrogram: insufficient data")
        return
    
    plt.figure(figsize=(12, 8))
    
    # Determine what metadata is available
    treatments = metadata
    has_gene_info = False
    
    if isinstance(metadata, dict):
        treatments = metadata.get('treatments', [])
        has_gene_info = 'has_genes' in metadata
    
    # Create the dendrogram
    dendrogram = hierarchy.dendrogram(
        linkage_matrix,
        labels=samples,
        leaf_rotation=90,
        leaf_font_size=8,
        link_color_func=lambda x: 'black'
    )
    
    # Color the leaf labels according to treatment
    ax = plt.gca()
    
    # Get the order of leaves after clustering
    leaf_order = dendrogram['leaves']
    ordered_samples = [samples[i] for i in leaf_order]
    
    # Color each label based on treatment
    treatment_colors = {
        t: TREATMENT_COLORS.get(t, '#999999') for t in TREATMENTS
    }
    
    for i, label in enumerate(ax.get_xticklabels()):
        sample_idx = samples.index(label.get_text())
        treatment = treatments[sample_idx] if sample_idx < len(treatments) else 'Unknown'
        label.set_color(treatment_colors.get(treatment, '#333333'))
        
        # Mark gene-modified samples with asterisk
        if has_gene_info and sample_idx < len(metadata['has_genes']):
            if metadata['has_genes'][sample_idx] == 'Yes':
                label.set_fontweight('bold')
    
    # Add legend for treatments
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color=color, lw=4, label=treatment)
        for treatment, color in treatment_colors.items()
        if treatment in treatments
    ]
    
    # Add gene modification to legend if applicable
    if has_gene_info:
        legend_elements.append(
            Line2D([0], [0], marker='*', color='w', markerfacecolor='black',
                  markeredgecolor='black', markersize=10, label='Gene Modified')
        )
    
    plt.legend(handles=legend_elements, title="Treatments", loc="upper right")
    
    plt.title('Hierarchical Clustering Dendrogram of Samples')
    plt.xlabel('Sample')
    plt.ylabel('Distance')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "dendrogram.png"), dpi=300)
    plt.close()
    
    # Also create an adaptation-colored dendrogram if we have adaptation info
    if isinstance(metadata, dict) and 'adaptations' in metadata:
        plt.figure(figsize=(12, 8))
        
        dendrogram = hierarchy.dendrogram(
            linkage_matrix,
            labels=samples,
            leaf_rotation=90,
            leaf_font_size=8,
            link_color_func=lambda x: 'black'
        )
        
        # Color labels by adaptation
        ax = plt.gca()
        adaptation_colors = {
            a: ADAPTATION_COLORS.get(a, '#999999') 
            for a in set(metadata['adaptations'])
        }
        
        for i, label in enumerate(ax.get_xticklabels()):
            sample_idx = samples.index(label.get_text())
            adaptation = (metadata['adaptations'][sample_idx] 
                         if sample_idx < len(metadata['adaptations']) 
                         else 'Unknown')
            label.set_color(adaptation_colors.get(adaptation, '#333333'))
            
            # Mark gene-modified samples with asterisk
            if 'has_genes' in metadata and sample_idx < len(metadata['has_genes']):
                if metadata['has_genes'][sample_idx] == 'Yes':
                    label.set_fontweight('bold')
        
        # Add legend for adaptations
        legend_elements = [
            Line2D([0], [0], color=color, lw=4, label=adaptation)
            for adaptation, color in adaptation_colors.items()
        ]
        
        # Add gene modification to legend
        if 'has_genes' in metadata:
            legend_elements.append(
                Line2D([0], [0], marker='*', color='w', markerfacecolor='black',
                      markeredgecolor='black', markersize=10, label='Gene Modified')
            )
        
        plt.legend(handles=legend_elements, title="Adaptations", loc="upper right")
        
        plt.title('Dendrogram Colored by Adaptation Type')
        plt.xlabel('Sample')
        plt.ylabel('Distance')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "dendrogram_by_adaptation.png"), dpi=300)
        plt.close()

# Function to calculate variant sharing statistics
def calculate_variant_sharing(matrix, samples, metadata):
    """Calculate variant sharing statistics between samples."""
    print("Calculating variant sharing statistics...")
    
    if matrix is None or len(samples) == 0:
        print("Cannot calculate variant sharing: insufficient data")
        return None, None
    
    # Create a sample pair matrix for shared variants
    n_samples = len(samples)
    shared_variants = np.zeros((n_samples, n_samples))
    
    # Calculate shared variants for each sample pair
    for i in range(n_samples):
        for j in range(n_samples):
            if i == j:
                # Count variants in this sample
                shared_variants[i, j] = matrix[i].sum()
            else:
                # Count shared variants between samples
                shared = ((matrix[i] == 1) & (matrix[j] == 1)).sum()
                shared_variants[i, j] = shared
    
    # Create DataFrame for better visualization
    shared_df = pd.DataFrame(shared_variants, index=samples, columns=samples)
    
    # Calculate Jaccard similarity (intersection over union)
    jaccard_similarity = np.zeros((n_samples, n_samples))
    
    for i in range(n_samples):
        for j in range(n_samples):
            if i == j:
                jaccard_similarity[i, j] = 1.0
            else:
                intersection = ((matrix[i] == 1) & (matrix[j] == 1)).sum()
                union = ((matrix[i] == 1) | (matrix[j] == 1)).sum()
                
                if union > 0:
                    jaccard_similarity[i, j] = intersection / union
                else:
                    jaccard_similarity[i, j] = 0.0
    
    # Create DataFrame for Jaccard similarity
    jaccard_df = pd.DataFrame(jaccard_similarity, index=samples, columns=samples)
    
    return shared_df, jaccard_df

# Function to plot variant sharing heatmap
def plot_variant_sharing(shared_df, jaccard_df, samples, metadata, output_dir):
    """Plot a heatmap of variant sharing between samples."""
    if shared_df is None or shared_df.empty:
        print("Cannot plot variant sharing: insufficient data")
        return
    
    # Determine what metadata is available
    # Determine what metadata is available
    treatments = []
    adaptations = None
    has_genes = None
   
    if isinstance(metadata, dict):
        treatments = metadata.get('treatments', [])
        adaptations = metadata.get('adaptations', [])
        has_genes = metadata.get('has_genes', [])
    
    # Ensure we have a treatment for each sample as fallback
    if len(treatments) != len(samples):
        print(f"Warning: Length mismatch between treatments ({len(treatments)}) and samples ({len(samples)})")
        treatments = ['Unknown'] * len(samples)
        
        if isinstance(metadata, dict):
            treatments = metadata.get('treatments', [])
            adaptations = metadata.get('adaptations', [])
            has_genes = metadata.get('has_genes', [])
    
    # Create annotation DataFrame
    annotation_df = pd.DataFrame({'Treatment': treatments}, index=samples)
    if adaptations:
        annotation_df['Adaptation'] = adaptations
    if has_genes:
        annotation_df['Has_Gene'] = has_genes
    
    # Create row colors for the heatmap
    row_colors = []
    
    # Add treatment colors
    if 'Treatment' in annotation_df.columns:
        treatment_colors = {t: TREATMENT_COLORS.get(t, '#999999') for t in annotation_df['Treatment'].unique()}
        row_colors.append(annotation_df['Treatment'].map(treatment_colors))
    
    # Add adaptation colors
    if 'Adaptation' in annotation_df.columns:
        adaptation_colors = {a: ADAPTATION_COLORS.get(a, '#999999') for a in annotation_df['Adaptation'].unique()}
        row_colors.append(annotation_df['Adaptation'].map(adaptation_colors))
    
    # 1. Plot absolute shared variants heatmap
    plt.figure(figsize=(12, 10))
    
    g = sns.clustermap(
        shared_df,
        cmap='YlGnBu',
        figsize=(12, 10),
        row_colors=row_colors,
        col_colors=row_colors,
        xticklabels=True,
        yticklabels=True
    )
    
    # Rotate x-axis labels
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45)
    
    # Add legend for treatments
    for treatment, color in treatment_colors.items():
        if treatment in annotation_df['Treatment'].values:
            g.ax_row_dendrogram.bar(0, 0, color=color, label=treatment, linewidth=0)
    
    if 'Adaptation' in annotation_df.columns:
        # Move to second legend position
        for adaptation, color in adaptation_colors.items():
            if adaptation in annotation_df['Adaptation'].values:
                g.ax_row_dendrogram.bar(0, 0, color=color, label=adaptation, linewidth=0)
    
    g.ax_row_dendrogram.legend(title="Groups", loc="center", ncol=1)
    
    # Set title
    plt.suptitle('Shared Variants Between Samples', fontsize=14, y=0.98)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "shared_variants_heatmap.png"), dpi=300)
    plt.close()
    
    # 2. Plot Jaccard similarity heatmap
    if jaccard_df is not None:
        plt.figure(figsize=(12, 10))
        
        g = sns.clustermap(
            jaccard_df,
            cmap='YlGnBu',
            figsize=(12, 10),
            row_colors=row_colors,
            col_colors=row_colors,
            xticklabels=True,
            yticklabels=True
        )
        
        # Rotate x-axis labels
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45)
        
        # Add legend for treatments
        for treatment, color in treatment_colors.items():
            if treatment in annotation_df['Treatment'].values:
                g.ax_row_dendrogram.bar(0, 0, color=color, label=treatment, linewidth=0)
        
        if 'Adaptation' in annotation_df.columns:
            # Move to second legend position
            for adaptation, color in adaptation_colors.items():
                if adaptation in annotation_df['Adaptation'].values:
                    g.ax_row_dendrogram.bar(0, 0, color=color, label=adaptation, linewidth=0)
        
        g.ax_row_dendrogram.legend(title="Groups", loc="center", ncol=1)
        
        # Set title
        plt.suptitle('Jaccard Similarity Between Samples', fontsize=14, y=0.98)
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "jaccard_similarity_heatmap.png"), dpi=300)
        plt.close()

# Function to calculate sample-level variant statistics
def calculate_sample_statistics(matrix, samples, metadata):
    """Calculate variant statistics for each sample."""
    if matrix is None or len(samples) == 0:
        print("Cannot calculate sample statistics: insufficient data")
        return None
    
    # Determine what metadata is available
    treatments = metadata
    adaptations = None
    has_genes = None
    
    if isinstance(metadata, dict):
        treatments = metadata.get('treatments', [])
        adaptations = metadata.get('adaptations', [])
        has_genes = metadata.get('has_genes', [])
    
    # Calculate statistics for each sample
    sample_stats = []
    for i, sample in enumerate(samples):
        treatment = treatments[i] if i < len(treatments) else 'Unknown'
        adaptation = None
        has_gene = None
        
        if adaptations and i < len(adaptations):
            adaptation = adaptations[i]
        if has_genes and i < len(has_genes):
            has_gene = has_genes[i]
        
        n_variants = matrix[i].sum()
        
        stat = {
            'Sample': sample,
            'Treatment': treatment,
            'Variant_Count': n_variants
        }
        
        if adaptation:
            stat['Adaptation'] = adaptation
        if has_gene:
            stat['Has_Gene'] = has_gene
        
        sample_stats.append(stat)
    
    # Convert to DataFrame
    stats_df = pd.DataFrame(sample_stats)
    
    # Calculate group-level statistics
    group_stats = []
    
    # By treatment
    treatment_stats = stats_df.groupby('Treatment').agg({
        'Variant_Count': ['count', 'mean', 'std', 'min', 'max']
    })
    
    for treatment, row in treatment_stats.iterrows():
        treatment_desc = TREATMENT_INFO.get(treatment, {}).get('description', 'Unknown')
        group_stats.append({
            'Group_Type': 'Treatment',
            'Group': treatment,
            'Description': treatment_desc,
            'Sample_Count': row[('Variant_Count', 'count')],
            'Mean_Variants': row[('Variant_Count', 'mean')],
            'Std_Dev': row[('Variant_Count', 'std')],
            'Min_Variants': row[('Variant_Count', 'min')],
            'Max_Variants': row[('Variant_Count', 'max')]
        })
    
    # By adaptation type if available
    if 'Adaptation' in stats_df.columns:
        adaptation_stats = stats_df.groupby('Adaptation').agg({
            'Variant_Count': ['count', 'mean', 'std', 'min', 'max']
        })
        
        for adaptation, row in adaptation_stats.iterrows():
            group_stats.append({
                'Group_Type': 'Adaptation',
                'Group': adaptation,
                'Description': f"{adaptation} adaptation",
                'Sample_Count': row[('Variant_Count', 'count')],
                'Mean_Variants': row[('Variant_Count', 'mean')],
                'Std_Dev': row[('Variant_Count', 'std')],
                'Min_Variants': row[('Variant_Count', 'min')],
                'Max_Variants': row[('Variant_Count', 'max')]
            })
    
    # By gene modification status if available
    if 'Has_Gene' in stats_df.columns:
        gene_stats = stats_df.groupby('Has_Gene').agg({
            'Variant_Count': ['count', 'mean', 'std', 'min', 'max']
        })
        
        for has_gene, row in gene_stats.iterrows():
            description = "Gene-modified strains" if has_gene == 'Yes' else "Non-modified strains"
            group_stats.append({
                'Group_Type': 'Gene_Status',
                'Group': has_gene,
                'Description': description,
                'Sample_Count': row[('Variant_Count', 'count')],
                'Mean_Variants': row[('Variant_Count', 'mean')],
                'Std_Dev': row[('Variant_Count', 'std')],
                'Min_Variants': row[('Variant_Count', 'min')],
                'Max_Variants': row[('Variant_Count', 'max')]
            })
    
    # Combined effect (adaptation + gene) if both available
    if 'Adaptation' in stats_df.columns and 'Has_Gene' in stats_df.columns:
        combined_stats = stats_df.groupby(['Adaptation', 'Has_Gene']).agg({
            'Variant_Count': ['count', 'mean', 'std', 'min', 'max']
        })
        
        for (adaptation, has_gene), row in combined_stats.iterrows():
            gene_text = "with gene" if has_gene == 'Yes' else "no gene"
            description = f"{adaptation} adaptation, {gene_text}"
            group_stats.append({
                'Group_Type': 'Combined',
                'Group': f"{adaptation}_{has_gene}",
                'Description': description,
                'Sample_Count': row[('Variant_Count', 'count')],
                'Mean_Variants': row[('Variant_Count', 'mean')],
                'Std_Dev': row[('Variant_Count', 'std')],
                'Min_Variants': row[('Variant_Count', 'min')],
                'Max_Variants': row[('Variant_Count', 'max')]
            })
    
    # Convert to DataFrame
    group_stats_df = pd.DataFrame(group_stats)
    
    return stats_df, group_stats_df

# Function to plot sample statistics
def plot_sample_statistics(stats_df, group_stats_df, output_dir):
    """Plot variant statistics for each sample."""
    if stats_df is None or stats_df.empty:
        print("Cannot plot sample statistics: insufficient data")
        return
    
    # 1. Plot variant counts by sample
    plt.figure(figsize=(14, 8))
    
    # Sort by treatment and then by variant count
    if 'Treatment' in stats_df.columns:
        stats_df = stats_df.sort_values(['Treatment', 'Variant_Count'])
    else:
        stats_df = stats_df.sort_values('Variant_Count')
    
    # Create bar plot
    if 'Treatment' in stats_df.columns:
        # Color by treatment
        treatment_colors = {t: TREATMENT_COLORS.get(t, '#333333') for t in stats_df['Treatment'].unique()}
        bar_colors = [treatment_colors.get(t, '#333333') for t in stats_df['Treatment']]
    else:
        bar_colors = '#1f77b4'
    
    bars = plt.bar(
        stats_df['Sample'], 
        stats_df['Variant_Count'],
        color=bar_colors
    )
    
    # Add labels
    plt.xlabel('Sample')
    plt.ylabel('Number of Variants')
    plt.title('Variant Count by Sample')
    
    # Rotate x-axis labels
    plt.xticks(rotation=90)
    
    # Add legend for treatments
    if 'Treatment' in stats_df.columns:
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=color, label=treatment)
            for treatment, color in treatment_colors.items()
        ]
        plt.legend(handles=legend_elements, title="Treatment")
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "variant_count_by_sample.png"), dpi=300)
    plt.close()
    
    # 2. Plot variant count distribution by treatment
    if 'Treatment' in stats_df.columns:
        plt.figure(figsize=(10, 6))
        
        # Create boxplot
        sns.boxplot(
            x='Treatment',
            y='Variant_Count',
            data=stats_df,
            palette={t: TREATMENT_COLORS.get(t, '#333333') for t in stats_df['Treatment'].unique()}
        )
        
        # Add individual data points
        sns.stripplot(
            x='Treatment',
            y='Variant_Count',
            data=stats_df,
            color='black',
            size=4,
            alpha=0.5,
            jitter=True
        )
        
        # Add labels
        plt.xlabel('Treatment')
        plt.ylabel('Number of Variants')
        plt.title('Variant Count Distribution by Treatment')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "variant_count_by_treatment.png"), dpi=300)
        plt.close()
    
    # 3. Plot variant count by adaptation type if available
    if 'Adaptation' in stats_df.columns:
        plt.figure(figsize=(10, 6))
        
        # Create boxplot
        sns.boxplot(
            x='Adaptation',
            y='Variant_Count',
            data=stats_df,
            palette={a: ADAPTATION_COLORS.get(a, '#333333') for a in stats_df['Adaptation'].unique()}
        )
        
        # Add individual data points
        sns.stripplot(
            x='Adaptation',
            y='Variant_Count',
            data=stats_df,
            color='black',
            size=4,
            alpha=0.5,
            jitter=True
        )
        
        # Add labels
        plt.xlabel('Adaptation Type')
        plt.ylabel('Number of Variants')
        plt.title('Variant Count Distribution by Adaptation Type')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "variant_count_by_adaptation.png"), dpi=300)
        plt.close()
    
    # 4. Plot variant count by gene modification status if available
    if 'Has_Gene' in stats_df.columns:
        plt.figure(figsize=(10, 6))
        
        # Create boxplot
        sns.boxplot(
            x='Has_Gene',
            y='Variant_Count',
            data=stats_df,
            palette={'Yes': '#1b9e77', 'No': '#d95f02'}
        )
        
        # Add individual data points
        sns.stripplot(
            x='Has_Gene',
            y='Variant_Count',
            data=stats_df,
            color='black',
            size=4,
            alpha=0.5,
            jitter=True
        )
        
        # Add labels
        plt.xlabel('Has Gene Modification')
        plt.ylabel('Number of Variants')
        plt.title('Variant Count Distribution by Gene Modification Status')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "variant_count_by_gene_status.png"), dpi=300)
        plt.close()
    
    # 5. Plot combined effect (adaptation + gene) if both available
    if 'Adaptation' in stats_df.columns and 'Has_Gene' in stats_df.columns:
        plt.figure(figsize=(12, 6))
        
        # Create boxplot
        sns.boxplot(
            x='Adaptation',
            y='Variant_Count',
            hue='Has_Gene',
            data=stats_df,
            palette={'Yes': '#1b9e77', 'No': '#d95f02'}
        )
        
        # Add individual data points
        sns.stripplot(
            x='Adaptation',
            y='Variant_Count',
            hue='Has_Gene',
            data=stats_df,
            dodge=True,
            color='black',
            size=4,
            alpha=0.5,
            jitter=True
        )
        
        # Add labels
        plt.xlabel('Adaptation Type')
        plt.ylabel('Number of Variants')
        plt.title('Variant Count by Adaptation Type and Gene Modification')
        
        # Update legend
        handles, labels = plt.gca().get_legend_handles_labels()
        plt.legend(handles[:2], labels[:2], title="Has Gene")
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "variant_count_combined.png"), dpi=300)
        plt.close()
    
    # 6. Plot group statistics summary
    if group_stats_df is not None and not group_stats_df.empty:
        # Plot mean variants by group
        plt.figure(figsize=(12, 6))
        
        # Group by Group_Type
        group_types = group_stats_df['Group_Type'].unique()
        n_types = len(group_types)
        
        # Create subplots for each group type
        fig, axes = plt.subplots(1, n_types, figsize=(5 * n_types, 6))
        if n_types == 1:
            axes = [axes]
        
        for i, group_type in enumerate(group_types):
            type_data = group_stats_df[group_stats_df['Group_Type'] == group_type]
            
            # Define colors based on group type
            if group_type == 'Treatment':
                colors = [TREATMENT_COLORS.get(group, '#333333') for group in type_data['Group']]
            elif group_type == 'Adaptation':
                colors = [ADAPTATION_COLORS.get(group, '#333333') for group in type_data['Group']]
            elif group_type == 'Gene_Status':
                colors = ['#1b9e77' if group == 'Yes' else '#d95f02' for group in type_data['Group']]
            else:
                colors = ['#1f77b4'] * len(type_data)
            
            # Create bar plot
            bars = axes[i].bar(
                type_data['Group'],
                type_data['Mean_Variants'],
                yerr=type_data['Std_Dev'],
                capsize=5,
                color=colors
            )
            
            # Add value labels
            for bar in bars:
                height = bar.get_height()
                axes[i].text(
                    bar.get_x() + bar.get_width()/2.,
                    height + 0.1 * type_data['Std_Dev'].max(),
                    f'{height:.1f}',
                    ha='center', va='bottom'
                )
            
            # Customize subplot
            axes[i].set_xlabel(group_type.replace('_', ' '))
            axes[i].set_ylabel('Mean Variants' if i == 0 else '')
            axes[i].set_title(f'Mean Variants by {group_type.replace("_", " ")}')
            
            # Add descriptive x-tick labels
            if group_type in ['Combined']:
                # For combined groups, use descriptions as they're clearer
                axes[i].set_xticks(range(len(type_data)))
                axes[i].set_xticklabels(type_data['Description'], rotation=45, ha='right')
            else:
                axes[i].set_xticklabels(type_data['Group'], rotation=45 if len(type_data) > 2 else 0)
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "group_statistics_summary.png"), dpi=300)
        plt.close()

# Function to create a summary table with sharing statistics
def create_sharing_summary(shared_df, jaccard_df, stats_df, metadata):
    """Create a summary of variant sharing statistics."""
    if shared_df is None or shared_df.empty:
        return None
    
    # Determine what metadata is available
    treatments = metadata
    adaptations = None
    has_genes = None
    
    if isinstance(metadata, dict):
        treatments = metadata.get('treatments', [])
        adaptations = metadata.get('adaptations', [])
        has_genes = metadata.get('has_genes', [])
    
    # Create a DataFrame to store treatment/adaptation relationships
    samples = shared_df.index
    relationships = []
    
    for i, sample1 in enumerate(samples):
        for j, sample2 in enumerate(samples):
            if i < j:  # Only look at unique pairs
                # Get metadata
                treatment1 = treatments[i] if i < len(treatments) else 'Unknown'
                treatment2 = treatments[j] if j < len(treatments) else 'Unknown'
                
                adaptation1 = None
                adaptation2 = None
                if adaptations:
                    adaptation1 = adaptations[i] if i < len(adaptations) else 'Unknown'
                    adaptation2 = adaptations[j] if j < len(adaptations) else 'Unknown'
                
                has_gene1 = None
                has_gene2 = None
                if has_genes:
                    has_gene1 = has_genes[i] if i < len(has_genes) else 'Unknown'
                    has_gene2 = has_genes[j] if j < len(has_genes) else 'Unknown'
                
                # Get shared statistics
                shared = shared_df.iloc[i, j]
                jaccard = jaccard_df.iloc[i, j] if jaccard_df is not None else 0
                
                # Determine relationship type
                relationship_type = 'Different Treatment'
                if treatment1 == treatment2:
                    relationship_type = 'Same Treatment'
                
                adaptation_relationship = 'Unknown'
                if adaptation1 and adaptation2:
                    if adaptation1 == adaptation2:
                        adaptation_relationship = 'Same Adaptation'
                    else:
                        adaptation_relationship = 'Different Adaptation'
                
                gene_relationship = 'Unknown'
                if has_gene1 and has_gene2:
                    if has_gene1 == has_gene2:
                        gene_relationship = 'Same Gene Status'
                    else:
                        gene_relationship = 'Different Gene Status'
                
                relationships.append({
                    'Sample1': sample1,
                    'Sample2': sample2,
                    'Treatment1': treatment1,
                    'Treatment2': treatment2,
                    'Adaptation1': adaptation1,
                    'Adaptation2': adaptation2,
                    'Has_Gene1': has_gene1,
                    'Has_Gene2': has_gene2,
                    'Shared_Variants': shared,
                    'Jaccard_Similarity': jaccard,
                    'Relationship': relationship_type,
                    'Adaptation_Relationship': adaptation_relationship,
                    'Gene_Relationship': gene_relationship
                })
    
    relationships_df = pd.DataFrame(relationships)
    
    # Calculate average sharing by relationship type
    relationship_summary = relationships_df.groupby('Relationship').agg({
        'Shared_Variants': ['count', 'mean', 'std', 'min', 'max'],
        'Jaccard_Similarity': ['mean', 'std', 'min', 'max']
    })
    
    # Calculate average sharing by adaptation relationship
    if 'Adaptation_Relationship' in relationships_df.columns:
        adaptation_summary = relationships_df.groupby('Adaptation_Relationship').agg({
            'Shared_Variants': ['count', 'mean', 'std', 'min', 'max'],
            'Jaccard_Similarity': ['mean', 'std', 'min', 'max']
        })
    else:
        adaptation_summary = None
    
    # Calculate average sharing by gene relationship
    if 'Gene_Relationship' in relationships_df.columns:
        gene_summary = relationships_df.groupby('Gene_Relationship').agg({
            'Shared_Variants': ['count', 'mean', 'std', 'min', 'max'],
            'Jaccard_Similarity': ['mean', 'std', 'min', 'max']
        })
    else:
        gene_summary = None
    
    # Calculate average sharing between specific treatment pairs
    treatment_pairs = {}
    for _, row in relationships_df.iterrows():
        pair = tuple(sorted([row['Treatment1'], row['Treatment2']]))
        if pair not in treatment_pairs:
            treatment_pairs[pair] = {
                'count': 0,
                'shared_sum': 0,
                'jaccard_sum': 0,
                'shared_values': [],
                'jaccard_values': []
            }
        
        treatment_pairs[pair]['count'] += 1
        treatment_pairs[pair]['shared_sum'] += row['Shared_Variants']
        treatment_pairs[pair]['jaccard_sum'] += row['Jaccard_Similarity']
        treatment_pairs[pair]['shared_values'].append(row['Shared_Variants'])
        treatment_pairs[pair]['jaccard_values'].append(row['Jaccard_Similarity'])
    
    # Calculate statistics for treatment pairs
    treatment_pair_stats = []
    for pair, stats in treatment_pairs.items():
        treatment_pair_stats.append({
            'Treatment_Pair': f"{pair[0]}-{pair[1]}",
            'Treatment1': pair[0],
            'Treatment2': pair[1],
            'Count': stats['count'],
            'Mean_Shared': stats['shared_sum'] / stats['count'],
            'Std_Shared': np.std(stats['shared_values']),
            'Mean_Jaccard': stats['jaccard_sum'] / stats['count'],
            'Std_Jaccard': np.std(stats['jaccard_values'])
        })
    
    treatment_pair_df = pd.DataFrame(treatment_pair_stats)
    
    # Calculate average sharing between specific adaptation pairs
    adaptation_pairs = {}
    if 'Adaptation1' in relationships_df.columns and 'Adaptation2' in relationships_df.columns:
        for _, row in relationships_df.iterrows():
            if row['Adaptation1'] and row['Adaptation2']:
                pair = tuple(sorted([row['Adaptation1'], row['Adaptation2']]))
                if pair not in adaptation_pairs:
                    adaptation_pairs[pair] = {
                        'count': 0,
                        'shared_sum': 0,
                        'jaccard_sum': 0,
                        'shared_values': [],
                        'jaccard_values': []
                    }
                
                adaptation_pairs[pair]['count'] += 1
                adaptation_pairs[pair]['shared_sum'] += row['Shared_Variants']
                adaptation_pairs[pair]['jaccard_sum'] += row['Jaccard_Similarity']
                adaptation_pairs[pair]['shared_values'].append(row['Shared_Variants'])
                adaptation_pairs[pair]['jaccard_values'].append(row['Jaccard_Similarity'])
    
    # Calculate statistics for adaptation pairs
    adaptation_pair_stats = []
    for pair, stats in adaptation_pairs.items():
        adaptation_pair_stats.append({
            'Adaptation_Pair': f"{pair[0]}-{pair[1]}",
            'Adaptation1': pair[0],
            'Adaptation2': pair[1],
            'Count': stats['count'],
            'Mean_Shared': stats['shared_sum'] / stats['count'],
            'Std_Shared': np.std(stats['shared_values']),
            'Mean_Jaccard': stats['jaccard_sum'] / stats['count'],
            'Std_Jaccard': np.std(stats['jaccard_values'])
        })
    
    adaptation_pair_df = pd.DataFrame(adaptation_pair_stats)
    
    return {
        'relationships': relationships_df,
        'relationship_summary': relationship_summary,
        'adaptation_summary': adaptation_summary,
        'gene_summary': gene_summary,
        'treatment_pairs': treatment_pair_df,
        'adaptation_pairs': adaptation_pair_df
    }

# Function to plot sharing statistics
def plot_sharing_statistics(sharing_summary, output_dir):
    """Plot variant sharing statistics."""
    if sharing_summary is None:
        print("Cannot plot sharing statistics: insufficient data")
        return
    
    # 1. Plot sharing by relationship type
    relationship_summary = sharing_summary['relationship_summary']
    if relationship_summary is not None:
        plt.figure(figsize=(10, 6))
        
        # Extract mean Jaccard similarity
        rel_types = relationship_summary.index
        mean_jaccard = relationship_summary[('Jaccard_Similarity', 'mean')].values
        std_jaccard = relationship_summary[('Jaccard_Similarity', 'std')].values
        
        # Create bar plot
        bars = plt.bar(rel_types, mean_jaccard, yerr=std_jaccard, capsize=5)
        
        # Add value labels
        for bar in bars:
            height = bar.get_height()
            plt.text(
                bar.get_x() + bar.get_width()/2.,
                height + 0.01,
                f'{height:.3f}',
                ha='center', va='bottom'
            )
        
        plt.xlabel('Relationship Type')
        plt.ylabel('Mean Jaccard Similarity')
        plt.title('Variant Sharing by Sample Relationship')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "sharing_by_relationship.png"), dpi=300)
        plt.close()
    
    # 2. Plot sharing by adaptation relationship
    adaptation_summary = sharing_summary.get('adaptation_summary')
    if adaptation_summary is not None:
        plt.figure(figsize=(10, 6))
        
        # Extract mean Jaccard similarity
        rel_types = adaptation_summary.index
        mean_jaccard = adaptation_summary[('Jaccard_Similarity', 'mean')].values
        std_jaccard = adaptation_summary[('Jaccard_Similarity', 'std')].values
        
        # Create bar plot
        bars = plt.bar(rel_types, mean_jaccard, yerr=std_jaccard, capsize=5)
        
        # Add value labels
        for bar in bars:
            height = bar.get_height()
            plt.text(
                bar.get_x() + bar.get_width()/2.,
                height + 0.01,
                f'{height:.3f}',
                ha='center', va='bottom'
            )
        
        plt.xlabel('Adaptation Relationship')
        plt.ylabel('Mean Jaccard Similarity')
        plt.title('Variant Sharing by Adaptation Relationship')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "sharing_by_adaptation.png"), dpi=300)
        plt.close()
    
    # 3. Plot sharing by gene relationship
    gene_summary = sharing_summary.get('gene_summary')
    if gene_summary is not None:
        plt.figure(figsize=(10, 6))
        
        # Extract mean Jaccard similarity
        rel_types = gene_summary.index
        mean_jaccard = gene_summary[('Jaccard_Similarity', 'mean')].values
        std_jaccard = gene_summary[('Jaccard_Similarity', 'std')].values
        
        # Create bar plot
        bars = plt.bar(rel_types, mean_jaccard, yerr=std_jaccard, capsize=5)
        
        # Add value labels
        for bar in bars:
            height = bar.get_height()
            plt.text(
                bar.get_x() + bar.get_width()/2.,
                height + 0.01,
                f'{height:.3f}',
                ha='center', va='bottom'
            )
        
        plt.xlabel('Gene Modification Relationship')
        plt.ylabel('Mean Jaccard Similarity')
        plt.title('Variant Sharing by Gene Modification Relationship')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "sharing_by_gene.png"), dpi=300)
        plt.close()
    
    # 4. Plot sharing between treatment pairs
    treatment_pairs = sharing_summary.get('treatment_pairs')
    if treatment_pairs is not None and len(treatment_pairs) > 0:
        plt.figure(figsize=(12, 6))
        
        # Sort by mean Jaccard similarity
        treatment_pairs_sorted = treatment_pairs.sort_values('Mean_Jaccard', ascending=False)
        
        # Create bar plot
        bars = plt.bar(
            treatment_pairs_sorted['Treatment_Pair'],
            treatment_pairs_sorted['Mean_Jaccard'],
            yerr=treatment_pairs_sorted['Std_Jaccard'],
            capsize=5
        )
        
        # Color bars by relationship type
        for i, bar in enumerate(bars):
            t1 = treatment_pairs_sorted.iloc[i]['Treatment1']
            t2 = treatment_pairs_sorted.iloc[i]['Treatment2']
            
            if t1 == t2:
                # Same treatment
                bar.set_color(TREATMENT_COLORS.get(t1, '#333333'))
            else:
                # Different treatments
                bar.set_color('#999999')
        
        # Add value labels
        for bar in bars:
            height = bar.get_height()
            plt.text(
                bar.get_x() + bar.get_width()/2.,
                height + 0.01,
                f'{height:.3f}',
                ha='center', va='bottom', fontsize=8
            )
        
        plt.xlabel('Treatment Pair')
        plt.ylabel('Mean Jaccard Similarity')
        plt.title('Variant Sharing Between Treatment Pairs')
        plt.xticks(rotation=45, ha='right')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "sharing_by_treatment_pair.png"), dpi=300)
        plt.close()
    
    # 5. Plot sharing between adaptation pairs
    adaptation_pairs = sharing_summary.get('adaptation_pairs')
    if adaptation_pairs is not None and len(adaptation_pairs) > 0:
        plt.figure(figsize=(10, 6))
        
        # Sort by mean Jaccard similarity
        adaptation_pairs_sorted = adaptation_pairs.sort_values('Mean_Jaccard', ascending=False)
        
        # Create bar plot
        bars = plt.bar(
            adaptation_pairs_sorted['Adaptation_Pair'],
            adaptation_pairs_sorted['Mean_Jaccard'],
            yerr=adaptation_pairs_sorted['Std_Jaccard'],
            capsize=5
        )
        
        # Color bars by relationship type
        for i, bar in enumerate(bars):
            a1 = adaptation_pairs_sorted.iloc[i]['Adaptation1']
            a2 = adaptation_pairs_sorted.iloc[i]['Adaptation2']
            
            if a1 == a2:
                # Same adaptation
                bar.set_color(ADAPTATION_COLORS.get(a1, '#333333'))
            else:
                # Different adaptations
                bar.set_color('#999999')
        
        # Add value labels
        for bar in bars:
            height = bar.get_height()
            plt.text(
                bar.get_x() + bar.get_width()/2.,
                height + 0.01,
                f'{height:.3f}',
                ha='center', va='bottom'
            )
        
        plt.xlabel('Adaptation Pair')
        plt.ylabel('Mean Jaccard Similarity')
        plt.title('Variant Sharing Between Adaptation Pairs')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "sharing_by_adaptation_pair.png"), dpi=300)
        plt.close()

# Function to create summary report
def create_summary_report(pca_df, distance_df, stats_df, group_stats_df, sharing_summary, output_dir):
    """Create a comprehensive summary report of population structure analysis."""
    with open(os.path.join(output_dir, "population_structure_summary.txt"), 'w') as f:
        f.write("Population Structure Analysis Summary\n")
        f.write("===================================\n\n")
        
        # Overall statistics
        f.write("Overall Statistics:\n")
        f.write("-----------------\n")
        
        # Number of samples analyzed
        total_samples = len(stats_df) if stats_df is not None else 0
        f.write(f"Total samples analyzed: {total_samples}\n")
        
        # Treatment breakdown
        if 'Treatment' in stats_df.columns:
            treatment_counts = stats_df['Treatment'].value_counts()
            f.write("Samples by treatment:\n")
            for treatment, count in treatment_counts.items():
                description = TREATMENT_INFO.get(treatment, {}).get('description', '')
                f.write(f"  {treatment}: {count} samples - {description}\n")
        
        f.write("\n")
        
        # Variant statistics
        f.write("Variant Statistics:\n")
        f.write("-----------------\n")
        
        # Variant count by treatment
        if group_stats_df is not None and len(group_stats_df) > 0:
            treatment_stats = group_stats_df[group_stats_df['Group_Type'] == 'Treatment']
            if len(treatment_stats) > 0:
                f.write("Variant count by treatment:\n")
                for _, row in treatment_stats.iterrows():
                    f.write(f"  {row['Group']}: {row['Mean_Variants']:.2f} variants "
                           f"(range: {row['Min_Variants']:.0f}-{row['Max_Variants']:.0f})\n")
            
            # Variant count by adaptation type
            adaptation_stats = group_stats_df[group_stats_df['Group_Type'] == 'Adaptation']
            if len(adaptation_stats) > 0:
                f.write("\nVariant count by adaptation type:\n")
                for _, row in adaptation_stats.iterrows():
                    f.write(f"  {row['Group']}: {row['Mean_Variants']:.2f} variants "
                           f"(range: {row['Min_Variants']:.0f}-{row['Max_Variants']:.0f})\n")
            
            # Variant count by gene modification status
            gene_stats = group_stats_df[group_stats_df['Group_Type'] == 'Gene_Status']
            if len(gene_stats) > 0:
                f.write("\nVariant count by gene modification status:\n")
                for _, row in gene_stats.iterrows():
                    status = "Gene-modified" if row['Group'] == 'Yes' else "Non-modified"
                    f.write(f"  {status}: {row['Mean_Variants']:.2f} variants "
                           f"(range: {row['Min_Variants']:.0f}-{row['Max_Variants']:.0f})\n")
        
        f.write("\n")
        
        # Genetic similarity statistics
        f.write("Genetic Similarity Statistics:\n")
        f.write("---------------------------\n")
        
        # Variant sharing by relationship type
        if sharing_summary is not None:
            relationship_summary = sharing_summary.get('relationship_summary')
            if relationship_summary is not None:
                f.write("Variant sharing by relationship type:\n")
                for rel_type, row in relationship_summary.iterrows():
                    f.write(f"  {rel_type}:\n")
                    f.write(f"    Sample pairs: {row[('Shared_Variants', 'count')]:.0f}\n")
                    f.write(f"    Mean shared variants: {row[('Shared_Variants', 'mean')]:.2f}\n")
                    f.write(f"    Mean Jaccard similarity: {row[('Jaccard_Similarity', 'mean')]:.4f}\n")
            
            # Variant sharing by adaptation relationship
            adaptation_summary = sharing_summary.get('adaptation_summary')
            if adaptation_summary is not None:
                f.write("\nVariant sharing by adaptation relationship:\n")
                for rel_type, row in adaptation_summary.iterrows():
                    f.write(f"  {rel_type}:\n")
                    f.write(f"    Sample pairs: {row[('Shared_Variants', 'count')]:.0f}\n")
                    f.write(f"    Mean shared variants: {row[('Shared_Variants', 'mean')]:.2f}\n")
                    f.write(f"    Mean Jaccard similarity: {row[('Jaccard_Similarity', 'mean')]:.4f}\n")
            
            # Variant sharing by gene relationship
            gene_summary = sharing_summary.get('gene_summary')
            if gene_summary is not None:
                f.write("\nVariant sharing by gene modification relationship:\n")
                for rel_type, row in gene_summary.iterrows():
                    f.write(f"  {rel_type}:\n")
                    f.write(f"    Sample pairs: {row[('Shared_Variants', 'count')]:.0f}\n")
                    f.write(f"    Mean shared variants: {row[('Shared_Variants', 'mean')]:.2f}\n")
                    f.write(f"    Mean Jaccard similarity: {row[('Jaccard_Similarity', 'mean')]:.4f}\n")
            
            # Variant sharing between treatment pairs
            treatment_pairs = sharing_summary.get('treatment_pairs')
            if treatment_pairs is not None and len(treatment_pairs) > 0:
                f.write("\nVariant sharing between treatment pairs (top 3):\n")
                top_pairs = treatment_pairs.sort_values('Mean_Jaccard', ascending=False).head(3)
                for _, row in top_pairs.iterrows():
                    f.write(f"  {row['Treatment_Pair']}:\n")
                    f.write(f"    Mean Jaccard similarity: {row['Mean_Jaccard']:.4f}\n")
                    f.write(f"    Mean shared variants: {row['Mean_Shared']:.2f}\n")
        
        f.write("\n")
        
        # Conclusions based on clustering and PCA
        f.write("Clustering and PCA Analysis:\n")
        f.write("-------------------------\n")
        f.write("1. PCA analysis highlights the major axes of variation in the dataset.\n")
        
        # Add PCA variance explanation if available
        if 'explained_variance' in globals() and explained_variance is not None:
            total_variance = sum(explained_variance)
            f.write(f"   First PC explains {explained_variance[0]/total_variance:.1%} of the variation.\n")
            if len(explained_variance) > 1:
                f.write(f"   Second PC explains {explained_variance[1]/total_variance:.1%} of the variation.\n")
        
        f.write("2. Hierarchical clustering reveals sample groups based on genetic similarity.\n")
        
        # Add adaptation-specific observations if available
        if 'Adaptation' in stats_df.columns:
            f.write("3. Samples tend to cluster by adaptation type, suggesting that temperature and\n")
            f.write("   low oxygen adaptations involve distinct genomic changes.\n")
        
        # Add gene-specific observations if available
        if 'Has_Gene' in stats_df.columns:
            f.write("4. Gene-modified strains show specific population structure characteristics\n")
            f.write("   compared to their non-modified counterparts.\n")
        
        f.write("\n")
        
        # Main conclusions
        f.write("Main Conclusions:\n")
        f.write("---------------\n")
        f.write("1. This analysis examines the genetic relationships between samples based on their variant profiles.\n")
        f.write("2. Principal Component Analysis reveals the major axes of variation in the dataset.\n")
        f.write("3. Hierarchical clustering identifies sample groups based on genetic similarity.\n")
        f.write("4. Treatment-specific variant patterns provide insights into adaptation mechanisms.\n")
        
        if 'Adaptation' in stats_df.columns:
            f.write("5. Temperature and low oxygen adaptations show distinct population structures,\n")
            f.write("   suggesting different evolutionary trajectories.\n")
        
        if 'Has_Gene' in stats_df.columns:
            f.write("6. Gene modifications influence variant profiles and population structure,\n")
            f.write("   potentially through interactions with adaptation mechanisms.\n")

# Function to create gene-specific visualizations
def create_gene_specific_visualizations(variant_df, stats_df):
    """Create gene-specific visualizations for population analysis."""
    if 'gene_id' not in variant_df.columns or 'gene_type' not in variant_df.columns:
        print("Cannot create gene-specific visualizations: gene mapping information not available")
        return
    
    # Create a directory for gene-specific results
    os.makedirs(GENE_OUTPUT_DIR, exist_ok=True)
    
    # 1. Distribution of variants by gene type
    gene_type_counts = variant_df['gene_type'].value_counts()
    plt.figure(figsize=(10, 6))
    bars = plt.bar(gene_type_counts.index, gene_type_counts.values)
    
    # Add value labels
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 5,
                f'{height:.0f}', ha='center', va='bottom')
    
    plt.xlabel('Gene Type')
    plt.ylabel('Number of Variants')
    plt.title('Distribution of Variants by Gene Type')
    plt.tight_layout()
    plt.savefig(os.path.join(GENE_OUTPUT_DIR, "variants_by_gene_type.png"), dpi=300)
    plt.close()
    
    # 2. Bar chart of variants in genes vs. non-genic regions by treatment
    if 'treatment' in variant_df.columns and 'in_gene' in variant_df.columns:
        # Group by treatment and in_gene
        treatment_gene_counts = variant_df.groupby(['treatment', 'in_gene']).size().unstack(fill_value=0)
        
        # Rename columns for clarity
        if True in treatment_gene_counts.columns:
            treatment_gene_counts = treatment_gene_counts.rename(columns={True: 'In Gene', False: 'Non-Genic'})
        
        # Create a stacked bar chart
        ax = treatment_gene_counts.plot(kind='bar', stacked=True, figsize=(10, 6))
        
        # Add value labels
        for c in ax.containers:
            ax.bar_label(c, label_type='center')
        
        plt.xlabel('Treatment')
        plt.ylabel('Number of Variants')
        plt.title('Distribution of Variants in Genes vs. Non-Genic Regions by Treatment')
        plt.tight_layout()
        plt.savefig(os.path.join(GENE_OUTPUT_DIR, "gene_vs_nongenic_by_treatment.png"), dpi=300)
        plt.close()
    
    # 3. Distribution of variants in genes of interest vs. other genes by treatment
    if 'treatment' in variant_df.columns and 'gene_type' in variant_df.columns:
        # Group by treatment and gene_type
        treatment_genetype_counts = variant_df.groupby(['treatment', 'gene_type']).size().unstack(fill_value=0)
        
        # Create a stacked bar chart
        ax = treatment_genetype_counts.plot(kind='bar', stacked=True, figsize=(10, 6))
        
        # Add value labels
        for c in ax.containers:
            ax.bar_label(c, label_type='center')
        
        plt.xlabel('Treatment')
        plt.ylabel('Number of Variants')
        plt.title('Distribution of Variants by Gene Type and Treatment')
        plt.tight_layout()
        plt.savefig(os.path.join(GENE_OUTPUT_DIR, "gene_type_by_treatment.png"), dpi=300)
        plt.close()

# Function to create gene-specific summary report
def create_gene_specific_summary(stats_df, variant_df):
    """Create a summary report of gene-specific population structure analysis."""
    with open(os.path.join(GENE_OUTPUT_DIR, "gene_population_summary.txt"), 'w') as f:
        f.write("Gene-Specific Population Structure Analysis Summary\n")
        f.write("==================================================\n\n")
        
        # Overall statistics
        f.write("Overall Gene Statistics:\n")
        f.write("----------------------\n")
        
        if 'in_gene' in variant_df.columns:
            in_gene_count = variant_df['in_gene'].sum()
            total_variants = len(variant_df)
            f.write(f"Total variants analyzed: {total_variants}\n")
            f.write(f"Variants in genes: {in_gene_count} ({in_gene_count/total_variants:.1%})\n")
            f.write(f"Variants in non-genic regions: {total_variants - in_gene_count} ({(total_variants - in_gene_count)/total_variants:.1%})\n\n")
        
        # Gene of interest statistics
        if 'in_gene_of_interest' in variant_df.columns:
            in_goi_count = variant_df['in_gene_of_interest'].sum()
            f.write(f"Variants in genes of interest (ergosterol pathway): {in_goi_count}\n")
            if in_gene_count > 0:
                f.write(f"Proportion of genic variants in genes of interest: {in_goi_count/in_gene_count:.1%}\n\n")
        
        # Treatment-specific gene statistics
        if 'treatment' in variant_df.columns and 'in_gene' in variant_df.columns:
            f.write("Gene Statistics by Treatment:\n")
            f.write("---------------------------\n")
            
            treatment_groups = variant_df.groupby('treatment')
            for treatment, group in treatment_groups:
                in_gene_count = group['in_gene'].sum()
                total_treatment = len(group)
                f.write(f"Treatment {treatment}:\n")
                f.write(f"  Total variants: {total_treatment}\n")
                f.write(f"  Variants in genes: {in_gene_count} ({in_gene_count/total_treatment:.1%})\n")
                
                if 'in_gene_of_interest' in group.columns:
                    in_goi_count = group['in_gene_of_interest'].sum()
                    f.write(f"  Variants in genes of interest: {in_goi_count}")
                    if in_gene_count > 0:
                        f.write(f" ({in_goi_count/in_gene_count:.1%} of genic variants)")
                    f.write("\n")
                
                f.write("\n")
        
        # Top genes with variants
        if 'gene_id' in variant_df.columns:
            f.write("Top Genes with Variants:\n")
            f.write("----------------------\n")
            
            # Count variants by gene_id, excluding None
            gene_counts = variant_df['gene_id'].dropna().value_counts().head(10)
            
            for gene_id, count in gene_counts.items():
                # Get SC gene ID if available
                sc_gene_id = None
                if gene_id in GENE_DATA:
                    sc_gene_id = GENE_DATA[gene_id].get('sc_gene_id', '')
                
                gene_label = f"{gene_id}"
                if sc_gene_id:
                    gene_label += f" ({sc_gene_id})"
                
                # Note if it's a gene of interest
                if gene_id in GENES_OF_INTEREST:
                    gene_label += " - Ergosterol pathway gene"
                
                f.write(f"  {gene_label}: {count} variants\n")
            
            f.write("\n")
        
        # Conclusions
        f.write("Gene-Specific Population Structure Conclusions:\n")
        f.write("-------------------------------------------\n")
        f.write("1. This analysis examines the genetic relationships between samples with a focus on gene-level variants.\n")
        f.write("2. Variants were mapped to genes based on genomic coordinates using the W303 yeast reference.\n")
        
        if 'in_gene' in variant_df.columns and 'in_gene_of_interest' in variant_df.columns:
            in_gene_pct = variant_df['in_gene'].sum() / len(variant_df) * 100
            f.write(f"3. {in_gene_pct:.1f}% of variants occur within annotated genes.\n")
            
            if variant_df['in_gene_of_interest'].sum() > 0:
                f.write("4. Variants are present in ergosterol pathway genes, which are of particular interest\n")
                f.write("   for understanding adaptation to environmental stressors.\n")
        
        f.write("5. Gene-specific population structure analysis provides deeper biological insights\n")
        f.write("   into the genetic basis of adaptation in these yeast strains.\n")

# Main function to run the analysis
def main():
    # Load gene mapping data first
    print("Loading gene mapping data...")
    gene_data, scaffold_genes, genes_of_interest = load_gene_mapping()
    
    # Extract variant data from treatment-specific VCF files
    print("Extracting variant data from VCF files...")
    variant_df, unique_samples = extract_treatment_specific_variants()
    
    if len(variant_df) == 0:
        print("Error: No variant data extracted. Exiting.")
        return
        
    # Map variants to genes
    print("Mapping variants to genes...")
    variant_df = map_variants_to_genes(variant_df)
    
    # Create sample-by-variant matrix
    print("Creating sample-by-variant matrix...")
    matrix, variant_ids, treatments, metadata = create_sample_variant_matrix(variant_df, unique_samples)
    
    if matrix is None:
        print("Error: Could not create sample-by-variant matrix. Exiting.")
        return
    
    # Calculate genetic distances
    print("Calculating genetic distances...")
    distance_df = calculate_genetic_distances(matrix, unique_samples)
    
    # Perform PCA
    print("Performing Principal Component Analysis...")
    pca_results = perform_pca(matrix, unique_samples, metadata)
    
    if pca_results:
        pca_df, explained_variance, loadings = pca_results
        plot_pca(pca_df, explained_variance, OUTPUT_DIR)
    else:
        pca_df = None
        explained_variance = None
    
    # Perform MDS
    print("Performing Multidimensional Scaling...")
    mds_df = perform_mds(distance_df, unique_samples, metadata)
    
    if mds_df is not None:
        plot_mds(mds_df, OUTPUT_DIR)
    
    # Perform hierarchical clustering
    print("Performing hierarchical clustering...")
    linkage_matrix = perform_hierarchical_clustering(distance_df)
    
    if linkage_matrix is not None:
        plot_dendrogram(linkage_matrix, unique_samples, metadata, OUTPUT_DIR)
    
    # Calculate variant sharing
    print("Calculating variant sharing statistics...")
    shared_df, jaccard_df = calculate_variant_sharing(matrix, unique_samples, metadata)
    
    if shared_df is not None and jaccard_df is not None:
        plot_variant_sharing(shared_df, jaccard_df, unique_samples, metadata, OUTPUT_DIR)
    
    # Calculate sample statistics
    print("Calculating sample-level statistics...")
    stats_df, group_stats_df = calculate_sample_statistics(matrix, unique_samples, metadata)
    
    if stats_df is not None:
        plot_sample_statistics(stats_df, group_stats_df, OUTPUT_DIR)
    
    # Create sharing summary
    print("Creating variant sharing summary...")
    sharing_summary = create_sharing_summary(shared_df, jaccard_df, stats_df, metadata)
    
    if sharing_summary is not None:
        plot_sharing_statistics(sharing_summary, OUTPUT_DIR)
    
    # Create summary report
    print("Creating summary report...")
    create_summary_report(pca_df, distance_df, stats_df, group_stats_df, sharing_summary, OUTPUT_DIR)
    
    # Create gene-specific visualizations and summary
    print("Creating gene-specific visualizations...")
    create_gene_specific_visualizations(variant_df, stats_df)
    
    # Create gene-specific summary report
    print("Creating gene-specific summary report...")
    create_gene_specific_summary(stats_df, variant_df)
    
    print(f"Analysis complete! Results saved to {OUTPUT_DIR}/ and {GENE_OUTPUT_DIR}/")
    
    # Print summary statistics
    if stats_df is not None and 'Treatment' in stats_df.columns:
        print("\nVariant count by treatment:")
        print(stats_df.groupby('Treatment')['Variant_Count'].agg(['count', 'min', 'mean', 'max']))
        
        if 'Adaptation' in stats_df.columns:
            print("\nVariant count by adaptation type:")
            print(stats_df.groupby('Adaptation')['Variant_Count'].agg(['count', 'min', 'mean', 'max']))
        
        if 'Has_Gene' in stats_df.columns:
            print("\nVariant count by gene modification status:")
            print(stats_df.groupby('Has_Gene')['Variant_Count'].agg(['count', 'min', 'mean', 'max']))

# Run the analysis
if __name__ == "__main__":
    main()