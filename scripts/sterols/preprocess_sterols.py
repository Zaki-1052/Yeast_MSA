#!/usr/bin/env python3
"""
Sterol data preprocessing script for Yeast MSA project.
This script performs initial processing and normalization of sterol profile data.

Input:
- sterol_data_with_sd.csv: Raw sterol data with standard deviations

Output:
- sterol_data_processed.csv: Processed and normalized sterol data
- basic statistics files in results/sterol_analysis/basic_stats/
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set paths
INPUT_FILE = 'sterol_data/sterol_data_with_sd.csv'
OUTPUT_FILE = 'sterol_data/sterol_data_processed.csv'
RESULTS_DIR = 'results/sterol_analysis'
STATS_DIR = f'{RESULTS_DIR}/basic_stats'
VIS_DIR = f'{RESULTS_DIR}/visualizations'

def ensure_directories():
    """Ensure all required directories exist."""
    for directory in [STATS_DIR, VIS_DIR]:
        Path(directory).mkdir(parents=True, exist_ok=True)

def load_sterol_data(file_path=INPUT_FILE):
    """Load sterol data from CSV file."""
    print(f"Loading sterol data from {file_path}")
    df = pd.read_csv(file_path)
    print(f"Loaded {len(df)} sterol measurements across {df['sample'].nunique()} samples")
    return df

def parse_sample_metadata(df):
    """Parse sample names into treatment, generation, and condition."""
    # Extract metadata from sample names
    # Format: [Treatment]_[Generation]_[Condition]
    
    # Initialize columns
    df['treatment'] = ''
    df['generation'] = ''
    df['condition'] = 'Unknown'
    df['adaptation_type'] = ''
    
    for idx, row in df.iterrows():
        sample = row['sample']
        parts = sample.split('_')
        
        # Treatment is always the first part
        df.at[idx, 'treatment'] = parts[0]
        
        # Second part is generation
        if len(parts) > 1:
            df.at[idx, 'generation'] = parts[1]
        
        # Condition is the third part if it exists
        if len(parts) > 2:
            df.at[idx, 'condition'] = parts[2]
        
        # Determine adaptation type based on condition or treatment
        if 'MA' in sample:
            df.at[idx, 'adaptation_type'] = 'Low Oxygen'
        elif '37C' in sample or df.at[idx, 'treatment'] in ['WT-37', 'CAS']:
            df.at[idx, 'adaptation_type'] = 'Temperature'
        elif df.at[idx, 'treatment'] in ['STC', 'WTA']:
            df.at[idx, 'adaptation_type'] = 'Low Oxygen'
    
    # Mark modified genes
    df['gene_modified'] = df['treatment'].apply(lambda x: x in ['CAS', 'STC'])
    
    return df

def normalize_sterol_data(df):
    """Normalize sterol concentrations within each sample."""
    # Create a copy of the dataframe
    normalized_df = df.copy()
    
    # For each sample, calculate total sterol content and normalize
    sample_totals = df.groupby('sample')['concentration'].sum()
    
    # Calculate relative abundance for each sterol
    for idx, row in normalized_df.iterrows():
        sample = row['sample']
        total = sample_totals[sample]
        normalized_df.at[idx, 'relative_abundance'] = row['concentration'] / total * 100
    
    return normalized_df

def calculate_basic_statistics(df):
    """Calculate basic statistics for sterol profiles."""
    # Sterol statistics by treatment
    sterol_by_treatment = df.pivot_table(
        index='treatment', 
        columns='sterol', 
        values='concentration', 
        aggfunc='mean'
    ).fillna(0)
    
    # Sterol statistics by adaptation type
    sterol_by_adaptation = df.pivot_table(
        index='adaptation_type', 
        columns='sterol', 
        values='concentration', 
        aggfunc='mean'
    ).fillna(0)
    
    # Sterol statistics by generation
    sterol_by_generation = df.pivot_table(
        index=['treatment', 'generation'], 
        columns='sterol', 
        values='concentration', 
        aggfunc='mean'
    ).fillna(0)
    
    # Coefficient of variation for each sterol
    cv_by_sterol = df.groupby('sterol').apply(
        lambda x: np.std(x['concentration']) / np.mean(x['concentration']) * 100
    ).sort_values(ascending=False)
    
    # Save statistics to files
    sterol_by_treatment.to_csv(f'{STATS_DIR}/sterol_by_treatment.csv')
    sterol_by_adaptation.to_csv(f'{STATS_DIR}/sterol_by_adaptation.csv')
    sterol_by_generation.to_csv(f'{STATS_DIR}/sterol_by_generation.csv')
    cv_by_sterol.to_csv(f'{STATS_DIR}/cv_by_sterol.csv', header=['cv_percent'])
    
    # Create summary statistics file
    with open(f'{STATS_DIR}/summary_statistics.txt', 'w') as f:
        f.write("# Sterol Profile Basic Statistics\n\n")
        
        f.write("## Overview\n")
        f.write(f"- Total samples: {df['sample'].nunique()}\n")
        f.write(f"- Total sterols measured: {df['sterol'].nunique()}\n")
        f.write(f"- Treatments: {', '.join(df['treatment'].unique())}\n")
        f.write(f"- Adaptation types: {', '.join(df['adaptation_type'].unique())}\n\n")
        
        f.write("## Ergosterol Levels\n")
        ergosterol_df = df[df['sterol'] == 'Ergosterol'].copy()
        for treatment in ergosterol_df['treatment'].unique():
            treatment_data = ergosterol_df[ergosterol_df['treatment'] == treatment]
            f.write(f"- {treatment}: {treatment_data['concentration'].mean():.2f} ± {treatment_data['std_dev'].mean():.2f}\n")
        f.write("\n")
        
        f.write("## Unique Sterols by Treatment\n")
        for treatment in df['treatment'].unique():
            treatment_sterols = df[df['treatment'] == treatment]['sterol'].unique()
            f.write(f"- {treatment}: {', '.join(treatment_sterols)}\n")
        f.write("\n")
        
        f.write("## Coefficient of Variation\n")
        f.write("Sterols sorted by variability (highest first):\n")
        for sterol, cv in cv_by_sterol.items():
            f.write(f"- {sterol}: {cv:.2f}%\n")
    
    return {
        'sterol_by_treatment': sterol_by_treatment,
        'sterol_by_adaptation': sterol_by_adaptation,
        'sterol_by_generation': sterol_by_generation,
        'cv_by_sterol': cv_by_sterol
    }

def create_basic_visualizations(df, stats):
    """Create basic visualizations for sterol data."""
    # Set up plotting style
    sns.set(style="whitegrid")
    plt.rcParams['figure.figsize'] = (12, 8)
    
    # 1. Heatmap of sterol profiles by treatment
    plt.figure(figsize=(14, 10))
    sterol_matrix = df.pivot_table(
        index='sample', 
        columns='sterol', 
        values='concentration', 
        aggfunc='mean'
    ).fillna(0)
    
    # Get sample order by treatment and generation
    sample_metadata = df.drop_duplicates('sample')[['sample', 'treatment', 'adaptation_type', 'generation']]
    sample_order = sample_metadata.sort_values(['adaptation_type', 'treatment', 'generation'])['sample'].tolist()
    
    # Reorder matrix
    sterol_matrix = sterol_matrix.reindex(sample_order)
    
    # Create heatmap
    sns.heatmap(sterol_matrix, cmap="YlGnBu", linewidths=0.5, annot=False)
    plt.title("Sterol Profiles Across Samples")
    plt.ylabel("Sample")
    plt.xlabel("Sterol")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f'{VIS_DIR}/sterol_heatmap.png', dpi=300)
    plt.close()
    
    # 2. Bar chart for ergosterol levels by treatment
    plt.figure(figsize=(12, 6))
    ergosterol_df = df[df['sterol'] == 'Ergosterol'].copy()
    
    # Sort by adaptation type and treatment
    order = ergosterol_df.sort_values(['adaptation_type', 'treatment', 'generation'])['sample'].unique()
    
    chart = sns.barplot(
        data=ergosterol_df,
        x='sample',
        y='concentration',
        order=order,
        palette='viridis'
    )
    
    # Add error bars
    for i, row in enumerate(ergosterol_df.iterrows()):
        _, data = row
        if pd.notna(data['std_dev']):
            plt.errorbar(
                i, data['concentration'], 
                yerr=data['std_dev'], 
                fmt='none', capsize=5, 
                color='black', alpha=0.7
            )
    
    chart.set_xticklabels(chart.get_xticklabels(), rotation=45, ha='right')
    plt.title("Ergosterol Levels Across Samples")
    plt.ylabel("Concentration")
    plt.tight_layout()
    plt.savefig(f'{VIS_DIR}/ergosterol_levels.png', dpi=300)
    plt.close()
    
    # 3. Sterol composition by treatment
    plt.figure(figsize=(15, 10))
    
    # Get treatment colors
    treatments = df['treatment'].unique()
    n_treatments = len(treatments)
    colors = plt.cm.viridis(np.linspace(0, 1, n_treatments))
    
    # Prepare sterol counts by treatment
    treatment_sterols = {}
    for i, treatment in enumerate(treatments):
        treatment_df = df[df['treatment'] == treatment]
        treatment_sterols[treatment] = treatment_df['sterol'].value_counts()
    
    # Set up bar positions
    x = np.arange(len(df['sterol'].unique()))
    width = 0.8 / n_treatments
    offsets = np.linspace(-(0.8/2) + (width/2), (0.8/2) - (width/2), n_treatments)
    
    # Plot bars
    for i, treatment in enumerate(treatments):
        sterol_counts = treatment_sterols[treatment]
        plt.bar(
            x + offsets[i], 
            [sterol_counts.get(sterol, 0) for sterol in df['sterol'].unique()], 
            width, 
            label=treatment,
            color=colors[i],
            alpha=0.7
        )
    
    plt.xticks(x, df['sterol'].unique(), rotation=45, ha='right')
    plt.legend(title='Treatment')
    plt.title('Sterol Diversity by Treatment')
    plt.ylabel('Count')
    plt.xlabel('Sterol')
    plt.tight_layout()
    plt.savefig(f'{VIS_DIR}/sterol_diversity.png', dpi=300)
    plt.close()
    
    # 4. Adaptation type comparison
    plt.figure(figsize=(12, 6))
    adaptation_ergosterol = df[
        (df['sterol'] == 'Ergosterol') & 
        (df['adaptation_type'].isin(['Temperature', 'Low Oxygen']))
    ].copy()
    
    sns.boxplot(
        data=adaptation_ergosterol,
        x='adaptation_type',
        y='concentration',
        palette='viridis'
    )
    
    # Add individual points
    sns.stripplot(
        data=adaptation_ergosterol,
        x='adaptation_type',
        y='concentration',
        color='black',
        alpha=0.7,
        jitter=True
    )
    
    plt.title("Ergosterol Levels by Adaptation Type")
    plt.tight_layout()
    plt.savefig(f'{VIS_DIR}/adaptation_comparison.png', dpi=300)
    plt.close()
    
    # 5. Generation comparison
    plt.figure(figsize=(14, 8))
    
    # Prepare data
    gen_data = df[df['sterol'] == 'Ergosterol'].copy()
    
    # Plot generation comparison for each treatment
    sns.catplot(
        data=gen_data,
        x='treatment',
        y='concentration',
        hue='generation',
        kind='bar',
        palette='viridis',
        height=6,
        aspect=1.5
    )
    
    plt.title("Ergosterol Levels by Generation and Treatment")
    plt.tight_layout()
    plt.savefig(f'{VIS_DIR}/generation_comparison.png', dpi=300)
    plt.close()
    
    print(f"Basic visualizations saved to {VIS_DIR}")

def main():
    """Main processing function."""
    ensure_directories()
    
    # Load data
    df = load_sterol_data()
    
    # Parse sample metadata
    df = parse_sample_metadata(df)
    
    # Normalize sterol data
    normalized_df = normalize_sterol_data(df)
    
    # Calculate basic statistics
    stats = calculate_basic_statistics(normalized_df)
    
    # Create basic visualizations
    create_basic_visualizations(normalized_df, stats)
    
    # Save processed data
    normalized_df.to_csv(OUTPUT_FILE, index=False)
    print(f"Processed sterol data saved to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()