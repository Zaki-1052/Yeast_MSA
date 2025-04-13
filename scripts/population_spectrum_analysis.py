#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist, squareform
import subprocess
from collections import defaultdict, Counter
import re
import warnings
warnings.filterwarnings('ignore')

# Set matplotlib style
plt.style.use('ggplot')
sns.set(font_scale=1.2)
sns.set_style("whitegrid")

# Define output directory
OUTPUT_DIR = "population_structure_results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Function to extract variant data from treatment-specific VCF files
def extract_treatment_specific_variants():
    """Extract treatment-specific variants from separate VCF files."""
    treatments = ['WT', 'STC', 'CAS', 'WTA']
    all_samples = []
    all_variant_data = []
    
    # Process each treatment
    for treatment in treatments:
        # Use treatment-specific VCF files
        vcf_file = f"results/merged/analysis/{treatment}_specific.vcf.gz"
        if not os.path.exists(vcf_file):
            print(f"Warning: File {vcf_file} not found")
            # Try alternative file name pattern
            vcf_file = f"results/merged/analysis/{treatment}_highconf.vcf.gz"
            if not os.path.exists(vcf_file):
                print(f"Warning: Alternative file {vcf_file} not found")
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
                            all_variant_data.append({
                                'sample': sample,
                                'variant_id': variant_id,
                                'treatment': treatment
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
        return None, None, None
    
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
    
    # Extract treatments from sample names
    treatments = [s.split('-')[0] if '-' in s else 'Unknown' for s in samples]
    
    # Verify we have variation in the data
    sample_counts = matrix.sum(axis=1)
    variant_counts = matrix.sum(axis=0)
    
    print(f"Sample variant counts - min: {sample_counts.min()}, max: {sample_counts.max()}")
    print(f"Variant presence counts - min: {variant_counts.min()}, max: {variant_counts.max()}")
    
    return matrix, variant_ids, treatments

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
def perform_pca(matrix, samples, n_components=3):
    """Perform PCA on the sample-by-variant matrix."""
    # Transpose matrix to get samples x variants
    sample_matrix = matrix.T
    
    # Perform PCA
    pca = PCA(n_components=n_components)
    principal_components = pca.fit_transform(sample_matrix)
    
    # Create DataFrame with PCA results
    columns = [f'PC{i+1}' for i in range(n_components)]
    pca_df = pd.DataFrame(data=principal_components, columns=columns, index=samples)
    
    # Add sample information
    pca_df['Sample'] = pca_df.index
    
    # Extract treatment information from sample names
    pca_df['Treatment'] = pca_df['Sample'].apply(
        lambda x: x.split('-')[0] if '-' in x else 'Unknown')
    
    # Calculate explained variance
    explained_variance = pca.explained_variance_ratio_
    
    return pca_df, explained_variance

# Function to perform MDS
def perform_mds(distance_matrix, n_components=2):
    """Perform Multidimensional Scaling on the distance matrix."""
    # Perform MDS
    mds = MDS(n_components=n_components, dissimilarity='precomputed', random_state=42)
    mds_result = mds.fit_transform(distance_matrix.values)
    
    # Create DataFrame with MDS results
    columns = [f'MDS{i+1}' for i in range(n_components)]
    mds_df = pd.DataFrame(data=mds_result, columns=columns, index=distance_matrix.index)
    
    # Add sample information
    mds_df['Sample'] = mds_df.index
    
    # Extract treatment information from sample names
    mds_df['Treatment'] = mds_df['Sample'].apply(
        lambda x: x.split('-')[0] if '-' in x else 'Unknown')
    
    return mds_df

# Function to perform hierarchical clustering
def perform_hierarchical_clustering(distance_matrix):
    """Perform hierarchical clustering on the distance matrix."""
    # Calculate linkage
    linkage_matrix = linkage(squareform(distance_matrix.values), method='average')
    
    return linkage_matrix

# Function to plot PCA results
def plot_pca(pca_df, explained_variance, output_dir):
    """Plot PCA results to visualize sample relationships."""
    # Define colors for treatments
    treatment_colors = {
        'WT': '#1b9e77',
        'STC': '#d95f02',
        'CAS': '#7570b3',
        'WTA': '#e7298a',
        'Unknown': '#999999'
    }
    
    # Plot PC1 vs PC2
    plt.figure(figsize=(10, 8))
    
    for treatment in pca_df['Treatment'].unique():
        subset = pca_df[pca_df['Treatment'] == treatment]
        plt.scatter(
            subset['PC1'], 
            subset['PC2'], 
            alpha=0.8, 
            s=100, 
            label=treatment,
            color=treatment_colors.get(treatment, '#999999')
        )
    
    # Add sample labels
    for i, row in pca_df.iterrows():
        plt.annotate(
            row['Sample'],
            (row['PC1'], row['PC2']),
            fontsize=8,
            alpha=0.7,
            xytext=(5, 5),
            textcoords='offset points'
        )
    
    # Add axis labels with explained variance
    plt.xlabel(f'PC1 ({explained_variance[0]:.2%} variance)')
    plt.ylabel(f'PC2 ({explained_variance[1]:.2%} variance)')
    
    plt.title('Principal Component Analysis of Sample Variant Profiles')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "pca_plot.png"), dpi=300)
    plt.close()
    
    # Plot PC1 vs PC3 if we have at least 3 components
    if 'PC3' in pca_df.columns:
        plt.figure(figsize=(10, 8))
        
        for treatment in pca_df['Treatment'].unique():
            subset = pca_df[pca_df['Treatment'] == treatment]
            plt.scatter(
                subset['PC1'], 
                subset['PC3'], 
                alpha=0.8, 
                s=100, 
                label=treatment,
                color=treatment_colors.get(treatment, '#999999')
            )
        
        # Add sample labels
        for i, row in pca_df.iterrows():
            plt.annotate(
                row['Sample'],
                (row['PC1'], row['PC3']),
                fontsize=8,
                alpha=0.7,
                xytext=(5, 5),
                textcoords='offset points'
            )
        
        # Add axis labels with explained variance
        plt.xlabel(f'PC1 ({explained_variance[0]:.2%} variance)')
        plt.ylabel(f'PC3 ({explained_variance[2]:.2%} variance)')
        
        plt.title('Principal Component Analysis (PC1 vs PC3)')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "pca_plot_pc1_pc3.png"), dpi=300)
        plt.close()
    
    # Plot the scree plot (explained variance by component)
    plt.figure(figsize=(10, 6))
    
    plt.bar(
        range(1, len(explained_variance) + 1), 
        explained_variance, 
        alpha=0.6, 
        color='steelblue'
    )
    
    plt.xlabel('Principal Component')
    plt.ylabel('Explained Variance Ratio')
    plt.title('Scree Plot: Explained Variance by Principal Component')
    plt.xticks(range(1, len(explained_variance) + 1))
    plt.grid(True, alpha=0.3)
    
    # Add cumulative explained variance line
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    ax2.plot(
        range(1, len(explained_variance) + 1),
        np.cumsum(explained_variance),
        color='red',
        marker='o'
    )
    ax2.set_ylabel('Cumulative Explained Variance')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "pca_scree_plot.png"), dpi=300)
    plt.close()

# Function to plot MDS results
def plot_mds(mds_df, output_dir):
    """Plot MDS results to visualize sample relationships."""
    # Define colors for treatments
    treatment_colors = {
        'WT': '#1b9e77',
        'STC': '#d95f02',
        'CAS': '#7570b3',
        'WTA': '#e7298a',
        'Unknown': '#999999'
    }
    
    # Plot MDS1 vs MDS2
    plt.figure(figsize=(10, 8))
    
    for treatment in mds_df['Treatment'].unique():
        subset = mds_df[mds_df['Treatment'] == treatment]
        plt.scatter(
            subset['MDS1'], 
            subset['MDS2'], 
            alpha=0.8, 
            s=100, 
            label=treatment,
            color=treatment_colors.get(treatment, '#999999')
        )
    
    # Add sample labels
    for i, row in mds_df.iterrows():
        plt.annotate(
            row['Sample'],
            (row['MDS1'], row['MDS2']),
            fontsize=8,
            alpha=0.7,
            xytext=(5, 5),
            textcoords='offset points'
        )
    
    plt.xlabel('MDS Dimension 1')
    plt.ylabel('MDS Dimension 2')
    plt.title('Multidimensional Scaling of Sample Genetic Distances')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "mds_plot.png"), dpi=300)
    plt.close()

# Function to plot distance matrix heatmap
def plot_distance_heatmap(distance_df, output_dir):
    """Plot a heatmap of the genetic distance matrix."""
    plt.figure(figsize=(12, 10))
    
    # Extract treatment info from sample names
    samples = distance_df.index
    treatments = [s.split('-')[0] if '-' in s else 'Unknown' for s in samples]
    
    # Create annotation DataFrame
    annotation_df = pd.DataFrame({'Treatment': treatments}, index=samples)
    
    # Define colors for treatments
    treatment_colors = {
        'WT': '#1b9e77',
        'STC': '#d95f02',
        'CAS': '#7570b3',
        'WTA': '#e7298a',
        'Unknown': '#999999'
    }
    
    treatment_palette = [treatment_colors.get(t, '#999999') for t in annotation_df['Treatment']]
    
    # Create the heatmap
    g = sns.clustermap(
        distance_df,
        cmap='viridis_r',  # Reversed viridis (darker = more similar)
        figsize=(12, 10),
        row_colors=treatment_palette,
        col_colors=treatment_palette,
        xticklabels=True,
        yticklabels=True
    )
    
    # Rotate x-axis labels
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45)
    
    # Add legend for treatments
    for treatment, color in treatment_colors.items():
        if treatment in treatments:
            g.ax_row_dendrogram.bar(0, 0, color=color, label=treatment, linewidth=0)
    g.ax_row_dendrogram.legend(title="Treatment", loc="center", ncol=1)
    
    # Set title
    plt.suptitle('Genetic Distance Heatmap with Hierarchical Clustering', fontsize=14, y=0.98)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "distance_heatmap.png"), dpi=300)
    plt.close()

# Function to plot dendrogram
def plot_dendrogram(linkage_matrix, samples, output_dir):
    """Plot a dendrogram of sample relationships."""
    plt.figure(figsize=(12, 8))
    
    # Extract treatment info from sample names
    treatments = [s.split('-')[0] if '-' in s else 'Unknown' for s in samples]
    
    # Define colors for treatments
    treatment_colors = {
        'WT': '#1b9e77',
        'STC': '#d95f02',
        'CAS': '#7570b3',
        'WTA': '#e7298a',
        'Unknown': '#999999'
    }
    
    # Map treatments to colors
    leaf_colors = [treatment_colors.get(t, '#999999') for t in treatments]
    
    # Create the dendrogram
    dendrogram(
        linkage_matrix,
        labels=samples,
        leaf_rotation=90,
        leaf_font_size=8,
        link_color_func=lambda x: 'black'
    )
    
    # Color the leaf labels according to treatment
    ax = plt.gca()
    for i, (color, label) in enumerate(zip(leaf_colors, ax.get_xticklabels())):
        label.set_color(color)
    
    # Add legend for treatments
    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0], [0], color=color, lw=4, label=treatment)
                      for treatment, color in treatment_colors.items()
                      if treatment in treatments]
    plt.legend(handles=legend_elements, title="Treatment", loc="upper right")
    
    plt.title('Hierarchical Clustering Dendrogram of Samples')
    plt.xlabel('Sample')
    plt.ylabel('Distance')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "dendrogram.png"), dpi=300)
    plt.close()

# Function to calculate variant sharing statistics
def calculate_variant_sharing(df, samples):
    """Calculate variant sharing statistics between samples."""
    if df is None or samples is None:
        return None
    
    # Get variant columns (sample names)
    sample_cols = [col for col in df.columns if col in samples]
    
    # Create a sample pair matrix for shared variants
    n_samples = len(sample_cols)
    shared_variants = np.zeros((n_samples, n_samples))
    
    # Calculate shared variants for each sample pair
    for i, sample1 in enumerate(sample_cols):
        for j, sample2 in enumerate(sample_cols):
            if i == j:
                # Count variants in this sample
                shared_variants[i, j] = df[sample1].sum()
            else:
                # Count shared variants between samples
                shared = ((df[sample1] == 1) & (df[sample2] == 1)).sum()
                shared_variants[i, j] = shared
    
    # Create DataFrame for better visualization
    shared_df = pd.DataFrame(shared_variants, index=sample_cols, columns=sample_cols)
    
    return shared_df

# Function to plot variant sharing heatmap
def plot_variant_sharing(shared_df, output_dir):
    """Plot a heatmap of variant sharing between samples."""
    plt.figure(figsize=(12, 10))
    
    # Extract treatment info from sample names
    samples = shared_df.index
    treatments = [s.split('-')[0] if '-' in s else 'Unknown' for s in samples]
    
    # Create annotation DataFrame
    annotation_df = pd.DataFrame({'Treatment': treatments}, index=samples)
    
    # Define colors for treatments
    treatment_colors = {
        'WT': '#1b9e77',
        'STC': '#d95f02',
        'CAS': '#7570b3',
        'WTA': '#e7298a',
        'Unknown': '#999999'
    }
    
    treatment_palette = [treatment_colors.get(t, '#999999') for t in annotation_df['Treatment']]
    
    # Create the heatmap
    g = sns.clustermap(
        shared_df,
        cmap='YlGnBu',
        figsize=(12, 10),
        row_colors=treatment_palette,
        col_colors=treatment_palette,
        xticklabels=True,
        yticklabels=True
    )
    
    # Rotate x-axis labels
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45)
    
    # Add legend for treatments
    for treatment, color in treatment_colors.items():
        if treatment in treatments:
            g.ax_row_dendrogram.bar(0, 0, color=color, label=treatment, linewidth=0)
    g.ax_row_dendrogram.legend(title="Treatment", loc="center", ncol=1)
    
    # Set title
    plt.suptitle('Variant Sharing Between Samples', fontsize=14, y=0.98)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "variant_sharing_heatmap.png"), dpi=300)
    plt.close()

# Function to calculate sample-level variant statistics
def calculate_sample_statistics(df, samples):
    """Calculate variant statistics for each sample."""
    if df is None or samples is None:
        return None
    
    # Get variant columns (sample names)
    sample_cols = [col for col in df.columns if col in samples]
    
    # Calculate statistics for each sample
    sample_stats = []
    for sample in sample_cols:
        n_variants = df[sample].sum()
        sample_stats.append({
            'Sample': sample,
            'Treatment': sample.split('-')[0] if '-' in sample else 'Unknown',
            'Variant_Count': n_variants
        })
    
    # Convert to DataFrame
    stats_df = pd.DataFrame(sample_stats)
    
    return stats_df

# Function to plot sample statistics
def plot_sample_statistics(stats_df, output_dir):
    """Plot variant statistics for each sample."""
    # Define colors for treatments
    treatment_colors = {
        'WT': '#1b9e77',
        'STC': '#d95f02',
        'CAS': '#7570b3',
        'WTA': '#e7298a',
        'Unknown': '#999999'
    }
    
    # Plot variant counts by sample
    plt.figure(figsize=(14, 8))
    
    # Sort by treatment and then by variant count
    stats_df = stats_df.sort_values(['Treatment', 'Variant_Count'])
    
    # Create bar plot
    bars = plt.bar(
        stats_df['Sample'], 
        stats_df['Variant_Count'],
        color=[treatment_colors.get(t, '#999999') for t in stats_df['Treatment']]
    )
    
    # Add labels
    plt.xlabel('Sample')
    plt.ylabel('Number of Variants')
    plt.title('Variant Count by Sample')
    
    # Rotate x-axis labels
    plt.xticks(rotation=90)
    
    # Add legend for treatments
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=color, label=treatment)
                      for treatment, color in treatment_colors.items()
                      if treatment in stats_df['Treatment'].values]
    plt.legend(handles=legend_elements, title="Treatment")
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "variant_count_by_sample.png"), dpi=300)
    plt.close()
    
    # Plot variant count distribution by treatment
    plt.figure(figsize=(10, 6))
    
    # Create boxplot
    sns.boxplot(
        x='Treatment',
        y='Variant_Count',
        data=stats_df,
        palette=treatment_colors
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
    plt.savefig(os.path.join(output_dir, "variant_count_boxplot.png"), dpi=300)
    plt.close()

# Function to create summary report
def create_summary_report(pca_df, distance_df, stats_df, output_dir):
    """Create a comprehensive summary report of population structure analysis."""
    with open(os.path.join(output_dir, "population_structure_summary.txt"), 'w') as f:
        f.write("Population Structure Analysis Summary\n")
        f.write("===================================\n\n")
        
        # Overall statistics
        f.write("Overall Statistics:\n")
        f.write("-----------------\n")
        
        # Number of samples analyzed
        total_samples = len(pca_df)
        f.write(f"Total samples analyzed: {total_samples}\n")
        
        # Treatment breakdown
        treatments = pca_df['Treatment'].unique()
        f.write("Samples by treatment:\n")
        for treatment in treatments:
            count = len(pca_df[pca_df['Treatment'] == treatment])
            f.write(f"  {treatment}: {count} samples\n")
        
        f.write("\n")
        
        # Variant statistics
        f.write("Variant Statistics:\n")
        f.write("-----------------\n")
        
        # Variant count by treatment
        f.write("Average variant count by treatment:\n")
        for treatment in treatments:
            treatment_stats = stats_df[stats_df['Treatment'] == treatment]
            avg_variants = treatment_stats['Variant_Count'].mean()
            min_variants = treatment_stats['Variant_Count'].min()
            max_variants = treatment_stats['Variant_Count'].max()
            f.write(f"  {treatment}: {avg_variants:.2f} variants (range: {min_variants}-{max_variants})\n")
        
        f.write("\n")
        
        # Sample distances
        f.write("Genetic Distance Statistics:\n")
        f.write("-------------------------\n")
        
        # Calculate average within-treatment and between-treatment distances
        treatment_groups = {t: [] for t in treatments}
        for sample, treatment in zip(pca_df['Sample'], pca_df['Treatment']):
            treatment_groups[treatment].append(sample)
        
        # Within-treatment distances
        f.write("Average within-treatment genetic distances:\n")
        for treatment in treatments:
            samples = treatment_groups[treatment]
            if len(samples) < 2:
                f.write(f"  {treatment}: N/A (only one sample)\n")
                continue
                
            # Calculate average distance between all sample pairs in this treatment
            within_distances = []
            for i, sample1 in enumerate(samples):
                for sample2 in samples[i+1:]:
                    within_distances.append(distance_df.loc[sample1, sample2])
            
            avg_distance = np.mean(within_distances)
            f.write(f"  {treatment}: {avg_distance:.4f}\n")
        
        f.write("\n")
        
        # Between-treatment distances
        f.write("Average between-treatment genetic distances:\n")
        for i, t1 in enumerate(treatments):
            for t2 in treatments[i+1:]:
                samples1 = treatment_groups[t1]
                samples2 = treatment_groups[t2]
                
                # Calculate average distance between all sample pairs across treatments
                between_distances = []
                for sample1 in samples1:
                    for sample2 in samples2:
                        between_distances.append(distance_df.loc[sample1, sample2])
                
                if between_distances:
                    avg_distance = np.mean(between_distances)
                    f.write(f"  {t1} vs {t2}: {avg_distance:.4f}\n")
        
        f.write("\n")
        
        # Main conclusions
        f.write("Main Conclusions:\n")
        f.write("---------------\n")
        f.write("1. This analysis examines the genetic relationships between samples based on their variant profiles.\n")
        f.write("2. Principal Component Analysis (PCA) reveals the major axes of variation in the dataset.\n")
        f.write("3. Hierarchical clustering identifies sample groups based on genetic similarity.\n")
        f.write("4. The distance heatmap visualizes the genetic relationships between all samples.\n")
        f.write("5. Treatment-specific patterns may provide insights into the effects of different conditions.\n")

# Main function to run the analysis
def main():
    # Extract variant data from treatment-specific VCF files
    variant_df, unique_samples = extract_treatment_specific_variants()
    
    if len(variant_df) == 0:
        print("Error: No variant data extracted. Exiting.")
        return
    
    # Create sample-by-variant matrix
    matrix, variant_ids, treatment_labels = create_sample_variant_matrix(variant_df, unique_samples)
    
    if matrix is None:
        print("Error: Could not create sample-by-variant matrix. Exiting.")
        return
    
    # Create treatment mapping for each sample
    sample_treatments = {sample: treatment for sample, treatment in zip(unique_samples, treatment_labels)}
    
    # Add sample count statistics
    sample_stats = []
    for i, sample in enumerate(unique_samples):
        sample_stats.append({
            'Sample': sample,
            'Treatment': treatment_labels[i],
            'Variant_Count': int(matrix[i].sum())
        })
    stats_df = pd.DataFrame(sample_stats)
    
    # Calculate genetic distances
    distance_df = calculate_genetic_distances(matrix, unique_samples)
    
    # Perform PCA
    print("Performing Principal Component Analysis...")
    n_components = min(3, min(matrix.shape[0], matrix.shape[1]) - 1)
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(matrix)
    
    # Create PCA DataFrame
    pca_df = pd.DataFrame(
        data=pca_result,
        columns=[f'PC{i+1}' for i in range(n_components)],
        index=unique_samples
    )
    pca_df['Sample'] = pca_df.index
    pca_df['Treatment'] = treatment_labels
    
    # Get explained variance
    explained_variance = pca.explained_variance_ratio_
    
    # Perform MDS
    print("Performing Multidimensional Scaling...")
    n_components_mds = min(2, distance_df.shape[0] - 1)
    if n_components_mds > 0:
        mds = MDS(n_components=n_components_mds, dissimilarity='precomputed', random_state=42)
        mds_result = mds.fit_transform(distance_df.values)
        
        # Create MDS DataFrame
        mds_df = pd.DataFrame(
            data=mds_result,
            columns=[f'MDS{i+1}' for i in range(n_components_mds)],
            index=unique_samples
        )
        mds_df['Sample'] = mds_df.index
        mds_df['Treatment'] = treatment_labels
    else:
        print("Warning: Not enough samples for MDS analysis")
        mds_df = None
    
    # Perform hierarchical clustering
    print("Performing hierarchical clustering...")
    if len(unique_samples) > 1:
        linkage_matrix = linkage(squareform(distance_df.values), method='average')
    else:
        print("Warning: Not enough samples for hierarchical clustering")
        linkage_matrix = None
    
    # Calculate variant sharing statistics
    print("Calculating variant sharing statistics...")
    shared_matrix = np.zeros((len(unique_samples), len(unique_samples)))
    
    # Group variant data by sample
    sample_variants = {sample: set() for sample in unique_samples}
    for _, row in variant_df.iterrows():
        sample_variants[row['sample']].add(row['variant_id'])
    
    # Calculate shared variants
    for i, sample1 in enumerate(unique_samples):
        vars1 = sample_variants[sample1]
        for j, sample2 in enumerate(unique_samples):
            vars2 = sample_variants[sample2]
            if i == j:
                shared_matrix[i, j] = len(vars1)  # Number of variants in this sample
            else:
                shared_matrix[i, j] = len(vars1.intersection(vars2))  # Shared variants
    
    shared_df = pd.DataFrame(shared_matrix, index=unique_samples, columns=unique_samples)
    
    # Generate visualizations
    print("Generating visualizations...")
    
    # Plot PCA results
    if 'PC2' in pca_df.columns:
        plot_pca(pca_df, explained_variance, OUTPUT_DIR)
    else:
        print("Warning: Not enough components for PCA visualization")
    
    # Plot MDS results
    if mds_df is not None and 'MDS2' in mds_df.columns:
        plot_mds(mds_df, OUTPUT_DIR)
    else:
        print("Warning: Not enough components for MDS visualization")
    
    # Plot distance heatmap
    plot_distance_heatmap(distance_df, OUTPUT_DIR)
    
    # Plot dendrogram
    if linkage_matrix is not None and len(unique_samples) > 1:
        plot_dendrogram(linkage_matrix, unique_samples, OUTPUT_DIR)
    else:
        print("Warning: Not enough samples for dendrogram visualization")
    
    # Plot variant sharing
    plot_variant_sharing(shared_df, OUTPUT_DIR)
    
    # Plot sample statistics
    plot_sample_statistics(stats_df, OUTPUT_DIR)
    
    # Create summary report
    create_summary_report(pca_df, distance_df, stats_df, OUTPUT_DIR)
    
    print(f"Analysis complete! Results saved to {OUTPUT_DIR}/")
    
    # Print summary statistics
    print("\nVariant count by treatment:")
    print(stats_df.groupby('Treatment')['Variant_Count'].agg(['count', 'min', 'mean', 'max']))
    
    print("\nSample genetic distances (min, max, mean):")
    non_zero_distances = distance_df.values[~np.eye(distance_df.shape[0], dtype=bool)]
    if len(non_zero_distances) > 0:
        min_dist = non_zero_distances.min()
        max_dist = non_zero_distances.max()
        mean_dist = non_zero_distances.mean()
        print(f"Min: {min_dist:.4f}, Max: {max_dist:.4f}, Mean: {mean_dist:.4f}")
    else:
        print("No distance data available")
        
# Run the analysis
if __name__ == "__main__":
    main()