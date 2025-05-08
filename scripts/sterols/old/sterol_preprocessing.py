import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load the sterol data
def load_sterol_data(filepath):
    """Load sterol data from CSV file and perform initial preprocessing."""
    df = pd.read_csv(filepath)
    
    # Parse sample names to extract treatment conditions
    df['treatment'] = df['sample'].apply(lambda x: x.split('_')[0])
    df['temperature'] = df['sample'].apply(lambda x: '55' if '55' in x else '5')
    df['condition'] = df['sample'].apply(lambda x: 'MA' if 'MA' in x else '37C')
    
    # Normalize concentrations if needed
    # df['concentration_normalized'] = ...
    
    return df

# Exploratory analysis
def explore_sterol_profiles(df):
    """Perform exploratory analysis of sterol profiles."""
    # Summary statistics
    summary = df.groupby(['treatment', 'sterol']).agg({
        'concentration': ['mean', 'std', 'min', 'max'],
        'std_dev': 'mean'
    }).reset_index()
    
    # Count unique sterols per treatment
    sterol_counts = df.groupby('treatment')['sterol'].nunique()
    
    # Distribution of sterol types
    sterol_dist = df.groupby('sterol')['concentration'].sum().sort_values(ascending=False)
    
    return {
        'summary': summary,
        'sterol_counts': sterol_counts,
        'sterol_dist': sterol_dist
    }

# Basic visualizations
def visualize_sterol_distributions(df):
    """Create basic visualizations of sterol distributions."""
    # Set up plotting
    plt.figure(figsize=(12, 8))
    
    # Ergosterol levels by treatment
    ergosterol_df = df[df['sterol'] == 'Ergosterol']
    sns.barplot(x='treatment', y='concentration', data=ergosterol_df)
    plt.title('Ergosterol Levels by Treatment')
    plt.savefig('results/sterol_analysis/ergosterol_by_treatment.png')
    
    # Sterol composition by treatment
    plt.figure(figsize=(14, 10))
    sns.barplot(x='treatment', y='concentration', hue='sterol', data=df)
    plt.title('Sterol Composition by Treatment')
    plt.xticks(rotation=45)
    plt.savefig('results/sterol_analysis/sterol_composition_by_treatment.png')