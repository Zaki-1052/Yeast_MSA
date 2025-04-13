#!/usr/bin/env python3

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from collections import defaultdict, Counter
import re
from scipy.stats import chi2_contingency

# Set matplotlib style for better visualizations
plt.style.use('ggplot')
sns.set(font_scale=1.2)
sns.set_style("whitegrid")

# Define output directory
OUTPUT_DIR = "mutation_spectrum_results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Define complementary base pairs and mutation categories
COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
TRANSITIONS = [('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')]
TRANSVERSIONS = [('A', 'C'), ('C', 'A'), ('A', 'T'), ('T', 'A'), 
                ('G', 'T'), ('T', 'G'), ('G', 'C'), ('C', 'G')]

# All possible single nucleotide substitutions
ALL_SUBSTITUTIONS = TRANSITIONS + TRANSVERSIONS

# Standardized substitution representation (use pyrimidine as reference)
STD_SUBSTITUTIONS = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']

# Function to parse mutation data files
def parse_mutation_data(treatment):
    """Parse mutation data for a specific treatment."""
    filename = f"mutation_spectrum_analysis/{treatment}_mutations.txt"
    if not os.path.exists(filename):
        print(f"Warning: File {filename} not found")
        return pd.DataFrame()
    
    # Read the mutation data file
    data = pd.read_csv(filename, sep='\t', header=None, 
                       names=['CHROM', 'POS', 'REF', 'ALT'])
    
    # Add treatment column
    data['Treatment'] = treatment
    
    return data

# Function to filter data for single nucleotide variants
def filter_snvs(data):
    """Filter data to include only single nucleotide variants."""
    # Keep only rows where REF and ALT are single nucleotides
    snv_data = data[(data['REF'].str.len() == 1) & (data['ALT'].str.len() == 1)]
    
    # Keep only ACGT bases (filter out N or other ambiguous bases)
    valid_bases = snv_data['REF'].isin(['A', 'C', 'G', 'T']) & snv_data['ALT'].isin(['A', 'C', 'G', 'T'])
    snv_data = snv_data[valid_bases]
    
    return snv_data

# Function to classify mutations
def classify_mutations(data):
    """Classify each mutation as transition or transversion."""
    if len(data) == 0:
        return data
    
    # Create mutation type column
    data['Mutation'] = data['REF'] + '>' + data['ALT']
    
    # Classify as transition or transversion
    data['Class'] = 'Unknown'
    for ref, alt in TRANSITIONS:
        mask = (data['REF'] == ref) & (data['ALT'] == alt)
        data.loc[mask, 'Class'] = 'Transition'
    
    for ref, alt in TRANSVERSIONS:
        mask = (data['REF'] == ref) & (data['ALT'] == alt)
        data.loc[mask, 'Class'] = 'Transversion'
    
    # Standardize mutation representation (pyrimidine-based)
    data['Std_Mutation'] = data.apply(standardize_mutation, axis=1)
    
    return data

# Function to standardize mutation representation
def standardize_mutation(row):
    """Convert mutation to standardized format with pyrimidine as reference."""
    ref, alt = row['REF'], row['ALT']
    
    # If reference is a purine (A or G), convert to pyrimidine-based
    if ref in ['A', 'G']:
        ref = COMPLEMENT[ref]
        alt = COMPLEMENT[alt]
    
    return f"{ref}>{alt}"

# Function to calculate transition/transversion ratio
def calculate_ti_tv_ratio(data):
    """Calculate the transition/transversion ratio."""
    if len(data) == 0:
        return 0
    
    transitions = len(data[data['Class'] == 'Transition'])
    transversions = len(data[data['Class'] == 'Transversion'])
    
    if transversions == 0:
        return float('inf')  # Avoid division by zero
    
    return transitions / transversions

# Function to count mutation types
def count_mutation_types(data):
    """Count occurrences of each mutation type."""
    if len(data) == 0:
        return {}
    
    # Count standard mutations
    counts = Counter(data['Std_Mutation'])
    
    # Ensure all possible substitutions are represented
    for sub in STD_SUBSTITUTIONS:
        if sub not in counts:
            counts[sub] = 0
    
    return counts

# Function to generate mutation spectrum plot
def plot_mutation_spectrum(mutation_counts, treatment, output_dir):
    """Generate mutation spectrum plot for a treatment."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Prepare data for plotting
    categories = STD_SUBSTITUTIONS
    values = [mutation_counts.get(cat, 0) for cat in categories]
    
    # Define colors for different mutation types
    colors = ['#2166ac', '#4393c3', '#92c5de', '#d6604d', '#f4a582', '#fddbc7']
    
    # Plot the bars
    bars = ax.bar(range(len(categories)), values, color=colors)
    
    # Add value labels on top of bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                f'{height}', ha='center', va='bottom')
    
    # Customize the plot
    ax.set_xticks(range(len(categories)))
    ax.set_xticklabels(categories, rotation=45)
    ax.set_xlabel('Mutation Type')
    ax.set_ylabel('Count')
    ax.set_title(f'Mutation Spectrum for {treatment} Treatment')
    
    # Add transition/transversion annotations
    ax.text(0.02, 0.95, 'Transversions', transform=ax.transAxes, 
            fontsize=12, va='top', color='#2166ac')
    ax.text(0.5, 0.95, 'Transitions', transform=ax.transAxes, 
            fontsize=12, va='top', color='#d6604d')
    
    # Add vertical lines to separate mutation types
    ax.axvline(x=2.5, color='black', linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{treatment}_mutation_spectrum.png"), dpi=300)
    plt.close()

# Function to plot comparative mutation spectrum
def plot_comparative_spectrum(all_counts, output_dir):
    """Generate comparative mutation spectrum plot for all treatments."""
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Number of categories and treatments
    categories = STD_SUBSTITUTIONS
    treatments = list(all_counts.keys())
    n_cats = len(categories)
    n_treatments = len(treatments)
    
    # Width of bars
    width = 0.8 / n_treatments
    
    # Define colors for treatments
    treatment_colors = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a']
    
    # Plot grouped bars
    for i, treatment in enumerate(treatments):
        values = [all_counts[treatment].get(cat, 0) for cat in categories]
        positions = [j + (i - n_treatments/2 + 0.5) * width for j in range(n_cats)]
        bars = ax.bar(positions, values, width, label=treatment, color=treatment_colors[i])
    
    # Customize the plot
    ax.set_xticks(range(n_cats))
    ax.set_xticklabels(categories, rotation=45)
    ax.set_xlabel('Mutation Type')
    ax.set_ylabel('Count')
    ax.set_title('Comparative Mutation Spectrum Across Treatments')
    ax.legend()
    
    # Add vertical line to separate transitions and transversions
    ax.axvline(x=2.5, color='black', linestyle='--', alpha=0.3)
    ax.text(0.02, 0.95, 'Transversions', transform=ax.transAxes, 
            fontsize=12, va='top')
    ax.text(0.5, 0.95, 'Transitions', transform=ax.transAxes, 
            fontsize=12, va='top')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "comparative_mutation_spectrum.png"), dpi=300)
    plt.close()

# Function to plot transition/transversion ratios
def plot_ti_tv_ratios(ratios, output_dir):
    """Plot transition/transversion ratios for all treatments."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot the bars
    treatments = list(ratios.keys())
    values = list(ratios.values())
    bars = ax.bar(treatments, values, color='#5ab4ac')
    
    # Add value labels on top of bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.05,
                f'{height:.2f}', ha='center', va='bottom')
    
    # Customize the plot
    ax.set_xlabel('Treatment')
    ax.set_ylabel('Ti/Tv Ratio')
    ax.set_title('Transition/Transversion Ratio by Treatment')
    ax.set_ylim(0, max(values) * 1.2)  # Add some space for labels
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "ti_tv_ratios.png"), dpi=300)
    plt.close()

# Function to perform statistical test on mutation patterns
def test_mutation_differences(all_counts):
    """Perform chi-square test to check if mutation patterns differ significantly."""
    # Prepare data for chi-square test
    treatments = list(all_counts.keys())
    categories = STD_SUBSTITUTIONS
    
    # Create contingency table
    table = []
    for treatment in treatments:
        row = [all_counts[treatment].get(cat, 0) for cat in categories]
        table.append(row)
    
    # Perform chi-square test
    chi2, p, dof, expected = chi2_contingency(table)
    
    return {
        'chi2': chi2,
        'p_value': p,
        'degrees_of_freedom': dof,
        'expected_counts': expected
    }

# Function to create a summary table
def create_summary_table(all_data, all_counts, ti_tv_ratios):
    """Create a summary table with key statistics."""
    summary = []
    
    for treatment in all_data.keys():
        data = all_data[treatment]
        total_snvs = len(data)
        transitions = len(data[data['Class'] == 'Transition'])
        transversions = len(data[data['Class'] == 'Transversion'])
        
        # Most common mutation
        std_counts = {k: v for k, v in all_counts[treatment].items()}
        most_common = max(std_counts.items(), key=lambda x: x[1]) if std_counts else ('N/A', 0)
        
        summary.append({
            'Treatment': treatment,
            'Total SNVs': total_snvs,
            'Transitions': transitions,
            'Transversions': transversions,
            'Ti/Tv Ratio': ti_tv_ratios[treatment],
            'Most Common': most_common[0],
            'Most Common Count': most_common[1]
        })
    
    return pd.DataFrame(summary)

# Main function to run the analysis
def main():
    treatments = ['WT', 'STC', 'CAS', 'WTA']
    
    # Parse data for each treatment
    all_raw_data = {}
    for treatment in treatments:
        all_raw_data[treatment] = parse_mutation_data(treatment)
    
    # Filter for SNVs and classify mutations
    all_data = {}
    for treatment, data in all_raw_data.items():
        snv_data = filter_snvs(data)
        all_data[treatment] = classify_mutations(snv_data)
    
    # Calculate transition/transversion ratios
    ti_tv_ratios = {}
    for treatment, data in all_data.items():
        ti_tv_ratios[treatment] = calculate_ti_tv_ratio(data)
    
    # Count mutation types
    all_counts = {}
    for treatment, data in all_data.items():
        all_counts[treatment] = count_mutation_types(data)
    
    # Generate mutation spectrum plots
    for treatment, counts in all_counts.items():
        plot_mutation_spectrum(counts, treatment, OUTPUT_DIR)
    
    # Generate comparative plot
    plot_comparative_spectrum(all_counts, OUTPUT_DIR)
    
    # Plot transition/transversion ratios
    plot_ti_tv_ratios(ti_tv_ratios, OUTPUT_DIR)
    
    # Perform statistical test
    test_results = test_mutation_differences(all_counts)
    
    # Create summary table
    summary_table = create_summary_table(all_data, all_counts, ti_tv_ratios)
    
    # Save summary table
    summary_table.to_csv(os.path.join(OUTPUT_DIR, "mutation_spectrum_summary.csv"), index=False)
    
    # Save test results
    with open(os.path.join(OUTPUT_DIR, "statistical_test_results.txt"), 'w') as f:
        f.write("Chi-square test for differences in mutation patterns:\n")
        f.write(f"Chi-square value: {test_results['chi2']:.4f}\n")
        f.write(f"p-value: {test_results['p_value']:.4f}\n")
        f.write(f"Degrees of freedom: {test_results['degrees_of_freedom']}\n")
        f.write("\nInterpretation:\n")
        if test_results['p_value'] < 0.05:
            f.write("The mutation patterns differ significantly between treatments (p < 0.05).\n")
        else:
            f.write("No significant difference detected in mutation patterns between treatments (p >= 0.05).\n")
    
    print(f"Analysis complete! Results saved to {OUTPUT_DIR}/")
    print(f"Summary table saved as {OUTPUT_DIR}/mutation_spectrum_summary.csv")
    print(f"Statistical test results saved as {OUTPUT_DIR}/statistical_test_results.txt")

# Run the analysis
if __name__ == "__main__":
    main()