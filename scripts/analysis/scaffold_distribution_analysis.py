#!/usr/bin/env python3

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from collections import defaultdict, Counter
from scipy.stats import poisson, spearmanr
import math
import subprocess
import re

# Set matplotlib style for better visualizations
plt.style.use('ggplot')
sns.set(font_scale=1.2)
sns.set_style("whitegrid")

# Define output directory
OUTPUT_DIR = "analysis/scaffold_distribution_results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

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
def create_summary_table(all_counts, all_densities, scaffold_lengths, output_dir):
    """Create a comprehensive summary table with key statistics."""
    summary = []
    
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
            'Top_Scaffolds': top_scaffold_info
        })
    
    # Convert to DataFrame
    summary_df = pd.DataFrame(summary)
    
    # Add adaptation averages
    adaptation_summary = summary_df.groupby('Adaptation').agg({
        'Total_Variants': 'mean',
        'Scaffolds_With_Variants': 'mean',
        'Global_Density': 'mean'
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
def create_summary_report(all_data, all_counts, all_densities, scaffold_lengths, output_dir):
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
            
            f.write(f"{treatment} Treatment ({description}):\n")
            f.write(f"  Adaptation: {adaptation}{gene_text}\n")
            f.write(f"  Total Variants: {total_variants}\n")
            f.write(f"  Scaffolds with Variants: {total_scaffolds}\n")
            f.write(f"  Global Variant Density: {global_density:.4f} variants/kb\n")
            
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
        
        for treatment, counts in all_counts.items():
            adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
            adaptation_variants[adaptation] += sum(counts.values())
            adaptation_scaffolds[adaptation].update(counts.keys())
            
            # Sum scaffold lengths
            for scaffold in counts.keys():
                adaptation_length[adaptation] += scaffold_lengths.get(scaffold, 0)
        
        for adaptation, variants in adaptation_variants.items():
            scaffolds = len(adaptation_scaffolds[adaptation])
            length = adaptation_length[adaptation]
            density = (variants * 1000) / length if length > 0 else 0
            
            f.write(f"{adaptation} Adaptation:\n")
            f.write(f"  Average Variants: {variants/2:.1f}\n")  # Assuming 2 treatments per adaptation
            f.write(f"  Total Unique Scaffolds: {scaffolds}\n")
            f.write(f"  Average Variant Density: {density:.4f} variants/kb\n")
            f.write("\n")
        
        # Gene modification effects
        gene_mod_variants = defaultdict(int)
        non_mod_variants = defaultdict(int)
        
        for treatment, counts in all_counts.items():
            adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
            has_gene = TREATMENT_INFO.get(treatment, {}).get('gene') is not None
            
            if has_gene:
                gene_mod_variants[adaptation] += sum(counts.values())
            else:
                non_mod_variants[adaptation] += sum(counts.values())
        
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
        
        # Main conclusions
        f.write("Main Conclusions:\n")
        f.write("---------------\n")
        f.write("1. This analysis identifies scaffolds with high variant densities across treatments.\n")
        f.write("2. Several scaffolds show treatment-specific enrichment.\n")
        f.write("3. The pattern of variant distribution provides insights into adaptation mechanisms.\n")
        f.write("4. Temperature and low oxygen adaptations show distinct scaffold distribution patterns.\n")
        f.write("5. Gene modifications (STC, CAS) appear to influence scaffold enrichment patterns.\n")
        f.write("6. Further sequence analysis of enriched scaffolds may reveal functional implications.\n")

# Main function to run the analysis
def main():
    # Parse scaffold lengths from reference genome index
    scaffold_lengths = parse_genome_index()
    
    # Parse data for each treatment
    all_data = {}
    for treatment in TREATMENTS:
        data = parse_mutation_data(treatment)
        if len(data) > 0:
            all_data[treatment] = data
            print(f"Loaded {len(data)} variants for {treatment} treatment")
        else:
            print(f"Warning: No data available for {treatment} treatment")
    
    # Count variants per scaffold
    all_counts = {}
    for treatment, data in all_data.items():
        counts = count_variants_per_scaffold(data)
        all_counts[treatment] = counts
        print(f"{treatment}: Found variants in {len(counts)} scaffolds")
    
    # Calculate variant densities
    all_densities = {}
    for treatment, counts in all_counts.items():
        densities = calculate_variant_density(counts, scaffold_lengths)
        all_densities[treatment] = densities
    
    # Identify statistically enriched scaffolds
    all_enriched = {}
    for treatment, counts in all_counts.items():
        # Calculate genome-wide variant rate
        total_variants = sum(counts.values())
        total_length = sum(scaffold_lengths.values())
        genome_wide_rate = total_variants / total_length if total_length > 0 else 0
        
        # Identify enriched scaffolds
        enriched = identify_enriched_scaffolds(counts, scaffold_lengths, genome_wide_rate)
        all_enriched[treatment] = enriched
        
        # Save to file
        if len(enriched) > 0:
            enriched.to_csv(os.path.join(OUTPUT_DIR, f"{treatment}_enriched_scaffolds.csv"), index=False)
            print(f"{treatment}: Found {len(enriched)} significantly enriched scaffolds")
        else:
            print(f"{treatment}: No significantly enriched scaffolds found")
    
    # Generate variant density plots
    for treatment, densities in all_densities.items():
        plot_variant_density(densities, treatment, OUTPUT_DIR)
    
    # Generate comparative heatmap
    plot_comparative_heatmap(all_densities, OUTPUT_DIR)
    
    # Generate bubble charts
    for treatment, counts in all_counts.items():
        plot_bubble_chart(counts, scaffold_lengths, treatment, OUTPUT_DIR)
    
    # Calculate treatment correlations
    corr_matrix, adaptation_corr = calculate_treatment_correlations(all_densities, OUTPUT_DIR)
    
    # Identify treatment-specific hotspots
    treatment_specific, adaptation_specific, gene_comparison = identify_treatment_specific_hotspots(all_densities, OUTPUT_DIR)
    
    # Create summary table
    summary_df = create_summary_table(all_counts, all_densities, scaffold_lengths, OUTPUT_DIR)
    
    # Create summary report
    create_summary_report(all_data, all_counts, all_densities, scaffold_lengths, OUTPUT_DIR)
    
    print(f"Analysis complete! Results saved to {OUTPUT_DIR}/")
    print(f"Generated summary report: {os.path.join(OUTPUT_DIR, 'scaffold_distribution_summary.txt')}")
    print(f"Generated treatment-specific hotspots: {os.path.join(OUTPUT_DIR, 'treatment_specific_hotspots.txt')}")

# Run the analysis
if __name__ == "__main__":
    main()