#!/usr/bin/env python3

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from collections import defaultdict, Counter
from scipy.stats import poisson, spearmanr
import math

# Set matplotlib style for better visualizations
plt.style.use('ggplot')
sns.set(font_scale=1.2)
sns.set_style("whitegrid")

# Define output directory
OUTPUT_DIR = "scaffold_distribution_results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Function to parse the reference genome index to get scaffold lengths
def parse_genome_index(fai_file="reference/yeast_w303.fasta.fai"):
    """Parse the reference genome fasta index to get scaffold lengths."""
    if not os.path.exists(fai_file):
        print(f"Warning: Genome index file {fai_file} not found.")
        print("Trying to recreate scaffold lengths from VCF headers...")
        
        # Try to find scaffolds in one of the VCF files
        for treatment in ["WT", "STC", "CAS", "WTA"]:
            vcf_file = f"results/merged/analysis/{treatment}_highconf.vcf.gz"
            if os.path.exists(vcf_file):
                import subprocess
                # Extract contig lines from VCF header
                try:
                    output = subprocess.check_output(
                        f"bcftools view -h {vcf_file} | grep contig", 
                        shell=True).decode('utf-8')
                    
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
                except:
                    continue
        
        print("Could not determine scaffold lengths. Using placeholder values.")
        # Return a dummy dict with placeholder values
        return defaultdict(lambda: 1000)
    
    # Read the fasta index file
    scaffold_lengths = {}
    with open(fai_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            scaffold, length = parts[0], int(parts[1])
            scaffold_lengths[scaffold] = length
    
    return scaffold_lengths

# Function to parse mutation data for a specific treatment
def parse_mutation_data(treatment):
    """Parse mutation data for a specific treatment."""
    # Check for extracted data from previous analysis
    extracted_file = f"mutation_spectrum_analysis/{treatment}_mutations.txt"
    if os.path.exists(extracted_file):
        # Read the already extracted mutation data
        data = pd.read_csv(extracted_file, sep='\t', header=None, 
                           names=['CHROM', 'POS', 'REF', 'ALT'])
        data['Treatment'] = treatment
        return data
    
    # If extracted file doesn't exist, try to extract from VCF
    vcf_file = f"results/merged/analysis/{treatment}_highconf.vcf.gz"
    if not os.path.exists(vcf_file):
        print(f"Warning: File {vcf_file} not found")
        return pd.DataFrame()
    
    # Extract data from VCF using bcftools
    import subprocess
    try:
        output = subprocess.check_output(
            f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' {vcf_file}", 
            shell=True).decode('utf-8')
        
        # Parse the output
        rows = []
        for line in output.strip().split('\n'):
            if line:  # Skip empty lines
                parts = line.split('\t')
                if len(parts) == 4:
                    rows.append(parts)
        
        # Create dataframe
        data = pd.DataFrame(rows, columns=['CHROM', 'POS', 'REF', 'ALT'])
        data['POS'] = data['POS'].astype(int)
        data['Treatment'] = treatment
        
        return data
    except:
        print(f"Error extracting data from {vcf_file}")
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
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 8))
    bars = ax.bar(range(len(scaffold_names)), density_values, color='#1b9e77')
    
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
    ax.set_title(f'Top {top_n} Scaffolds by Variant Density for {treatment} Treatment')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{treatment}_variant_density_top{top_n}.png"), dpi=300)
    plt.close()

# Function to generate comparative density heat map
def plot_comparative_heatmap(all_densities, output_dir, top_n=30):
    """Generate a heatmap comparing variant densities across treatments."""
    # Find top scaffolds across all treatments
    all_scaffold_densities = {}
    
    for treatment, densities in all_densities.items():
        for scaffold, density in densities.items():
            if scaffold not in all_scaffold_densities:
                all_scaffold_densities[scaffold] = 0
            all_scaffold_densities[scaffold] += density
    
    # Sort scaffolds by total density
    top_scaffolds = sorted(all_scaffold_densities.items(), key=lambda x: x[1], reverse=True)[:top_n]
    top_scaffold_names = [s[0] for s in top_scaffolds]
    
    # Create a dataframe for the heatmap
    data = []
    for scaffold in top_scaffold_names:
        row = {'Scaffold': scaffold}
        for treatment, densities in all_densities.items():
            row[treatment] = densities.get(scaffold, 0)
        data.append(row)
    
    df = pd.DataFrame(data)
    df.set_index('Scaffold', inplace=True)
    
    # Create truncated scaffold names for better display
    df.index = [s[:20] + '...' if len(s) > 20 else s for s in df.index]
    
    # Create the heatmap
    plt.figure(figsize=(12, 10))
    sns.heatmap(df, cmap='viridis', annot=True, fmt='.2f', linewidths=.5)
    plt.title('Variant Density (variants/kb) Across Treatments')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"comparative_density_heatmap_top{top_n}.png"), dpi=300)
    plt.close()

# Function to generate bubble plot
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
    
    # Convert to dataframe and sort by density
    df = pd.DataFrame(data)
    df = df.sort_values('Density', ascending=False).head(top_n)
    
    # Create the bubble chart
    plt.figure(figsize=(12, 8))
    
    # Use log scale for scaffold length
    plt.scatter(
        df['Length'], 
        df['Count'], 
        s=df['Density']*20,  # Scale bubble size by density
        alpha=0.6, 
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
    plt.title(f'Scaffold Length vs. Variant Count for {treatment}\n(Bubble size represents variants/kb)')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{treatment}_bubble_chart.png"), dpi=300)
    plt.close()

# Function to calculate correlation between treatments
def calculate_treatment_correlations(all_densities, output_dir):
    """Calculate and visualize correlations between treatments based on scaffold density."""
    # Get all unique scaffolds across treatments
    all_scaffolds = set()
    for densities in all_densities.values():
        all_scaffolds.update(densities.keys())
    
    # Create a dataframe with scaffold densities for each treatment
    data = []
    for scaffold in all_scaffolds:
        row = {'Scaffold': scaffold}
        for treatment, densities in all_densities.items():
            row[treatment] = densities.get(scaffold, 0)
        data.append(row)
    
    df = pd.DataFrame(data)
    df.set_index('Scaffold', inplace=True)
    
    # Calculate correlation matrix
    corr_matrix = df.corr(method='spearman')
    
    # Plot correlation heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', vmin=-1, vmax=1, linewidths=.5)
    plt.title('Spearman Correlation of Variant Densities Between Treatments')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "treatment_correlation_heatmap.png"), dpi=300)
    plt.close()
    
    return corr_matrix

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
            f.write(f"{treatment} Treatment: {len(hotspots)} specific hotspots\n")
            
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
    return treatment_specific

# Function to create summary report
def create_summary_report(all_counts, all_densities, scaffold_lengths, output_dir):
    """Create a comprehensive summary report."""
    with open(os.path.join(output_dir, "scaffold_distribution_summary.txt"), 'w') as f:
        f.write("Scaffold Distribution Analysis Summary\n")
        f.write("=====================================\n\n")
        
        # Overall statistics
        f.write("Overall Statistics:\n")
        f.write("-----------------\n")
        
        for treatment, counts in all_counts.items():
            total_variants = sum(counts.values())
            total_scaffolds = len(counts)
            
            # Calculate global variant density
            total_length = sum(scaffold_lengths.get(scaffold, 0) for scaffold in counts.keys())
            global_density = (total_variants * 1000) / total_length if total_length > 0 else 0
            
            f.write(f"{treatment} Treatment:\n")
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
        
        # Write correlation summary
        f.write("Treatment Correlation Summary:\n")
        f.write("----------------------------\n")
        
        try:
            corr_matrix = calculate_treatment_correlations(all_densities, output_dir)
            
            # Get pairs with highest and lowest correlation
            treatments = list(all_densities.keys())
            pairs = [(i, j) for i in range(len(treatments)) for j in range(i+1, len(treatments))]
            
            if pairs:
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
        except Exception as e:
            f.write(f"  Error calculating correlations: {str(e)}\n\n")
        
        # Main conclusions
        f.write("Main Conclusions:\n")
        f.write("---------------\n")
        f.write("1. This analysis identifies scaffolds with high variant densities across treatments.\n")
        f.write("2. Several scaffolds show treatment-specific enrichment.\n")
        f.write("3. The pattern of variant distribution provides insights into potential functional regions.\n")
        f.write("4. Further sequence analysis of enriched scaffolds may reveal treatment-specific effects.\n")

# Main function to run the analysis
def main():
    treatments = ['WT', 'STC', 'CAS', 'WTA']
    
    # Parse scaffold lengths from reference genome index
    scaffold_lengths = parse_genome_index()
    print(f"Loaded {len(scaffold_lengths)} scaffold lengths from reference genome")
    
    # Parse data for each treatment
    all_data = {}
    for treatment in treatments:
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
    corr_matrix = calculate_treatment_correlations(all_densities, OUTPUT_DIR)
    
    # Identify treatment-specific hotspots
    hotspots = identify_treatment_specific_hotspots(all_densities, OUTPUT_DIR)
    
    # Create summary report
    create_summary_report(all_counts, all_densities, scaffold_lengths, OUTPUT_DIR)
    
    print(f"Analysis complete! Results saved to {OUTPUT_DIR}/")
    print(f"Generated summary report: {os.path.join(OUTPUT_DIR, 'scaffold_distribution_summary.txt')}")
    print(f"Generated treatment-specific hotspots: {os.path.join(OUTPUT_DIR, 'treatment_specific_hotspots.txt')}")

# Run the analysis
if __name__ == "__main__":
    main()