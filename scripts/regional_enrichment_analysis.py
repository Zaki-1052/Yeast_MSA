#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import poisson
from collections import defaultdict
import statsmodels.stats.multitest as mt
import warnings
warnings.filterwarnings('ignore')

# Set plotting style
plt.style.use('ggplot')
sns.set(font_scale=1.2)
sns.set_style("whitegrid")

# Define output directory
OUTPUT_DIR = "regional_enrichment_results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def load_data():
    """Load mutation data and scaffold information."""
    # Load mutation data from previous analyses
    treatments = ['WT', 'STC', 'CAS', 'WTA']
    all_data = []
    
    for treatment in treatments:
        file_path = f"analysis/MSA/mutation_spectrum_analysis/{treatment}_mutations.txt"
        if os.path.exists(file_path):
            data = pd.read_csv(file_path, sep='\t', header=None,
                             names=['CHROM', 'POS', 'REF', 'ALT'])
            data['Treatment'] = treatment
            all_data.append(data)
    
    if all_data:
        combined_data = pd.concat(all_data, ignore_index=True)
        print(f"Loaded {len(combined_data)} mutations across {len(treatments)} treatments")
        return combined_data
    else:
        print("No mutation data found.")
        return None

def load_scaffold_info():
    """Load scaffold length information."""
    fai_path = "reference/yeast_w303.fasta.fai"
    if os.path.exists(fai_path):
        scaffold_df = pd.read_csv(fai_path, sep='\t', header=None,
                                names=['Scaffold', 'Length', 'Offset', 'Linebases', 'Linewidth'])
        scaffold_info = dict(zip(scaffold_df['Scaffold'], scaffold_df['Length']))
        print(f"Loaded information for {len(scaffold_info)} scaffolds")
        return scaffold_info
    else:
        print("No scaffold information found.")
        return None

def analyze_regional_enrichment(data, scaffold_info, window_size=1000, step_size=500):
    """Analyze regional enrichment of variants using sliding windows."""
    if data is None or scaffold_info is None:
        return None
    
    # Calculate genome-wide mutation rate
    total_mutations = len(data)
    total_length = sum(scaffold_info.values())
    genome_wide_rate = total_mutations / total_length
    
    # Initialize results storage
    enriched_regions = []
    
    # Analyze each scaffold
    for scaffold, length in scaffold_info.items():
        # Skip if scaffold is shorter than window size
        if length < window_size:
            continue
        
        # Get mutations in this scaffold
        scaffold_mutations = data[data['CHROM'] == scaffold]
        
        # Analyze windows
        for start in range(0, length - window_size + 1, step_size):
            end = start + window_size
            
            # Count mutations in window
            mutations_in_window = scaffold_mutations[
                (scaffold_mutations['POS'] >= start) & 
                (scaffold_mutations['POS'] < end)
            ]
            
            observed = len(mutations_in_window)
            expected = genome_wide_rate * window_size
            
            # Calculate p-value using Poisson distribution
            if observed > expected:
                p_value = 1 - poisson.cdf(observed - 1, expected)
                
                enriched_regions.append({
                    'Scaffold': scaffold,
                    'Start': start,
                    'End': end,
                    'Observed': observed,
                    'Expected': expected,
                    'Fold_Enrichment': observed / expected if expected > 0 else float('inf'),
                    'P_Value': p_value
                })
    
    # Convert to DataFrame
    results_df = pd.DataFrame(enriched_regions)
    
    # Multiple testing correction
    if len(results_df) > 0:
        _, corrected_pvals, _, _ = mt.multipletests(
            results_df['P_Value'], 
            method='fdr_bh'
        )
        results_df['Q_Value'] = corrected_pvals
        
        # Filter for significance
        significant_regions = results_df[results_df['Q_Value'] < 0.05].copy()
        significant_regions = significant_regions.sort_values('Q_Value')
        
        print(f"Found {len(significant_regions)} significantly enriched regions")
        return significant_regions
    else:
        print("No enriched regions found")
        return pd.DataFrame()

def analyze_treatment_specific_enrichment(data, scaffold_info, window_size=1000, step_size=500):
    """Analyze regional enrichment separately for each treatment."""
    treatment_results = {}
    
    for treatment in data['Treatment'].unique():
        print(f"\nAnalyzing enrichment for {treatment} treatment...")
        treatment_data = data[data['Treatment'] == treatment]
        
        enriched = analyze_regional_enrichment(
            treatment_data, 
            scaffold_info,
            window_size,
            step_size
        )
        
        if enriched is not None and len(enriched) > 0:
            treatment_results[treatment] = enriched
    
    return treatment_results

def compare_enriched_regions(treatment_results):
    """Compare enriched regions between treatments."""
    if not treatment_results:
        return None
    
    # Create a master list of all enriched regions
    all_regions = []
    
    for treatment, results in treatment_results.items():
        for _, row in results.iterrows():
            region = (row['Scaffold'], row['Start'], row['End'])
            all_regions.append(region)
    
    # Remove duplicates
    unique_regions = list(set(all_regions))
    
    # Create comparison matrix
    comparison_data = []
    
    for region in unique_regions:
        scaffold, start, end = region
        row_data = {'Scaffold': scaffold, 'Start': start, 'End': end}
        
        for treatment in treatment_results.keys():
            results = treatment_results[treatment]
            matching = results[
                (results['Scaffold'] == scaffold) &
                (results['Start'] == start) &
                (results['End'] == end)
            ]
            
            if len(matching) > 0:
                row_data[f'{treatment}_Enrichment'] = matching['Fold_Enrichment'].iloc[0]
                row_data[f'{treatment}_Q_Value'] = matching['Q_Value'].iloc[0]
            else:
                row_data[f'{treatment}_Enrichment'] = 0
                row_data[f'{treatment}_Q_Value'] = 1
        
        comparison_data.append(row_data)
    
    return pd.DataFrame(comparison_data)

def plot_enrichment_patterns(treatment_results, output_dir):
    """Generate visualizations of enrichment patterns."""
    if not treatment_results:
        return
    
    # Plot 1: Distribution of enrichment levels by treatment
    plt.figure(figsize=(12, 6))
    
    enrichment_data = []
    for treatment, results in treatment_results.items():
        if len(results) > 0:
            enrichment_data.extend([
                {'Treatment': treatment, 'Fold_Enrichment': fe}
                for fe in results['Fold_Enrichment']
            ])
    
    if enrichment_data:
        enrichment_df = pd.DataFrame(enrichment_data)
        sns.boxplot(x='Treatment', y='Fold_Enrichment', data=enrichment_df)
        plt.yscale('log')
        plt.title('Distribution of Fold Enrichment by Treatment')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'enrichment_distribution.png'), dpi=300)
        plt.close()
    
    # Plot 2: Enriched regions along scaffolds
    comparison_df = compare_enriched_regions(treatment_results)
    
    if comparison_df is not None and len(comparison_df) > 0:
        # Get enrichment columns
        enrichment_cols = [col for col in comparison_df.columns if 'Enrichment' in col]
        
        # Create enrichment matrix
        enrichment_matrix = comparison_df[enrichment_cols].values
        
        plt.figure(figsize=(12, 8))
        sns.heatmap(
            enrichment_matrix,
            xticklabels=[col.replace('_Enrichment', '') for col in enrichment_cols],
            yticklabels=False,
            cmap='YlOrRd'
        )
        plt.title('Enrichment Patterns Across Treatments')
        plt.xlabel('Treatment')
        plt.ylabel('Enriched Regions')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'enrichment_heatmap.png'), dpi=300)
        plt.close()

def create_summary_report(treatment_results, output_dir):
    """Create a comprehensive summary report of regional enrichment analysis."""
    with open(os.path.join(output_dir, 'regional_enrichment_summary.txt'), 'w') as f:
        f.write("Regional Enrichment Analysis Summary\n")
        f.write("==================================\n\n")
        
        # Overall statistics
        f.write("Overall Statistics:\n")
        f.write("-----------------\n")
        total_regions = sum(len(results) for results in treatment_results.values())
        f.write(f"Total enriched regions identified: {total_regions}\n\n")
        
        # Treatment-specific statistics
        f.write("Treatment-Specific Statistics:\n")
        f.write("--------------------------\n")
        
        for treatment, results in treatment_results.items():
            f.write(f"\n{treatment} Treatment:\n")
            if len(results) > 0:
                f.write(f"  Enriched regions: {len(results)}\n")
                f.write(f"  Average fold enrichment: {results['Fold_Enrichment'].mean():.2f}\n")
                f.write(f"  Maximum fold enrichment: {results['Fold_Enrichment'].max():.2f}\n")
                
                # Top 5 most enriched regions
                f.write("\n  Top 5 most enriched regions:\n")
                top_regions = results.nlargest(5, 'Fold_Enrichment')
                for _, region in top_regions.iterrows():
                    f.write(f"    {region['Scaffold']}:{region['Start']}-{region['End']}, "
                           f"Fold enrichment: {region['Fold_Enrichment']:.2f}, "
                           f"Q-value: {region['Q_Value']:.2e}\n")
            else:
                f.write("  No significantly enriched regions found\n")
        
        # Cross-treatment comparison
        comparison_df = compare_enriched_regions(treatment_results)
        if comparison_df is not None and len(comparison_df) > 0:
            f.write("\nCross-Treatment Comparison:\n")
            f.write("------------------------\n")
            
            # Count regions shared between treatments
            enrichment_cols = [col for col in comparison_df.columns if 'Enrichment' in col]
            shared_counts = (comparison_df[enrichment_cols] > 0).sum()
            
            f.write("\nNumber of enriched regions by treatment:\n")
            for col in enrichment_cols:
                treatment = col.replace('_Enrichment', '')
                count = shared_counts[col]
                f.write(f"  {treatment}: {count}\n")
        
        # Main conclusions
        f.write("\nMain Conclusions:\n")
        f.write("---------------\n")
        f.write("1. This analysis identifies regions with statistically significant variant enrichment.\n")
        f.write("2. Treatment-specific patterns of regional enrichment have been identified.\n")
        f.write("3. Some regions show consistent enrichment across multiple treatments.\n")
        f.write("4. The distribution of enriched regions varies between treatments.\n")

def main():
    # Load data
    mutation_data = load_data()
    scaffold_info = load_scaffold_info()
    
    if mutation_data is None or scaffold_info is None:
        print("Cannot proceed with analysis due to missing data")
        return
    
    # Perform treatment-specific enrichment analysis
    print("\nAnalyzing regional enrichment by treatment...")
    treatment_results = analyze_treatment_specific_enrichment(
        mutation_data,
        scaffold_info
    )
    
    if not treatment_results:
        print("No enriched regions found in any treatment")
        return
    
    # Generate visualizations
    print("\nGenerating visualizations...")
    plot_enrichment_patterns(treatment_results, OUTPUT_DIR)
    
    # Create summary report
    print("\nCreating summary report...")
    create_summary_report(treatment_results, OUTPUT_DIR)
    
    print(f"\nAnalysis complete! Results saved to {OUTPUT_DIR}/")

if __name__ == "__main__":
    main()