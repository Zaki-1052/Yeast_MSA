import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

def load_genomic_data():
    """Load the relevant genomic data for correlation with sterol profiles."""
    # This is a placeholder - actual implementation would load the variant data
    # from the previous genomic analyses
    
    ergosterol_genes = {
        'ERG1': {'chromosome': 'CM007967.1', 'position': 175000, 'gene_id': 'YGR175C'},
        'ERG2': {'chromosome': 'CM007973.1', 'position': 202000, 'gene_id': 'YMR202W'},
        'ERG3': {'chromosome': 'CM007971.1', 'position': 56000, 'gene_id': 'YLR056W'},
        'ERG4': {'chromosome': 'CM007969.1', 'position': 109000, 'gene_id': 'YGL012W'},
        'ERG5': {'chromosome': 'CM007968.1', 'position': 40000, 'gene_id': 'YMR015C'},
        'ERG6': {'chromosome': 'CM007972.1', 'position': 108000, 'gene_id': 'YML008C'},
        'ERG7': {'chromosome': 'CM007968.1', 'position': 258000, 'gene_id': 'YHR072W'},
        'ERG9': {'chromosome': 'CM007969.1', 'position': 18000, 'gene_id': 'YGR157W'},
        'ERG11': {'chromosome': 'CM007968.1', 'position': 7000, 'gene_id': 'YHR007C'},
        'ERG24': {'chromosome': 'CM007963.1', 'position': 302000, 'gene_id': 'YNL280C'},
        'ERG25': {'chromosome': 'CM007967.1', 'position': 60000, 'gene_id': 'YGR060W'}
    }
    
    satellite_genes = {
        'W3030H00610': {'chromosome': 'CM007968.1', 'position': 15149, 'distance_to': 'ERG11', 'distance': 8149},
        'W3030G02910': {'chromosome': 'CM007967.1', 'position': 44051, 'distance_to': 'ERG25', 'distance': 15949},
        'W3030G02200': {'chromosome': 'CM007969.1', 'position': 82870, 'distance_to': 'ERG4', 'distance': 26130},
        'W3030G03230': {'chromosome': 'CM007967.1', 'position': 100586, 'distance_to': 'ERG25', 'distance': 40586},
        'W3030L01080': {'chromosome': 'CM007971.1', 'position': 8394, 'distance_to': 'ERG3', 'distance': 47606},
        'W3030H01660': {'chromosome': 'CM007968.1', 'position': 305676, 'distance_to': 'ERG7', 'distance': 47676}
    }
    
    # Create a combined dataset for genetic context
    genetic_data = {
        'ergosterol_genes': ergosterol_genes,
        'satellite_genes': satellite_genes
    }
    
    return genetic_data

def map_sample_to_treatment(sample):
    """Map sample names to treatment conditions for correlation analysis."""
    if 'CAS_5' in sample:
        return 'CAS'
    elif 'CAS_55' in sample:
        return 'CAS'
    elif 'STC_5' in sample:
        return 'STC'
    elif 'STC_55' in sample:
        return 'STC'
    elif 'WT_5_37C' in sample:
        return 'WT-CTRL'
    elif 'WT_5_MA' in sample:
        return 'WTA'
    elif 'WT_55_37C' in sample:
        return 'WT-37'
    elif 'WT_55_MA' in sample:
        return 'WTA'
    else:
        return 'Unknown'

def correlate_sterols_with_variants(sterol_df, variant_counts):
    """Correlate sterol levels with variant counts across samples."""
    # Map samples to treatments for correlation
    sterol_df['genomic_treatment'] = sterol_df['sample'].apply(map_sample_to_treatment)
    
    # Aggregate sterol data by treatment
    sterol_by_treatment = sterol_df.groupby(['genomic_treatment', 'sterol'])['concentration'].mean().reset_index()
    
    # Example variant count data (placeholder)
    variant_data = {
        'WT-CTRL': 4,
        'WT-37': 12, 
        'WTA': 12,
        'STC': 16,
        'CAS': 16
    }
    
    # Add variant counts to sterol data
    sterol_by_treatment['variant_count'] = sterol_by_treatment['genomic_treatment'].map(variant_data)
    
    # Calculate correlation for each sterol
    correlation_results = []
    
    for sterol in sterol_by_treatment['sterol'].unique():
        sterol_data = sterol_by_treatment[sterol_by_treatment['sterol'] == sterol]
        
        if len(sterol_data) >= 3:  # Need at least 3 points for correlation
            corr, p_value = stats.pearsonr(
                sterol_data['variant_count'], 
                sterol_data['concentration']
            )
            
            correlation_results.append({
                'sterol': sterol,
                'pearson_correlation': corr,
                'p_value': p_value,
                'significant': p_value < 0.05
            })
    
    return pd.DataFrame(correlation_results)

def visualize_genomic_sterol_relationship(sterol_df, genetic_data, correlation_df):
    """Create visualizations showing relationships between genomic patterns and sterol profiles."""
    # Plot correlations
    plt.figure(figsize=(10, 6))
    sns.barplot(x='sterol', y='pearson_correlation', data=correlation_df)
    plt.axhline(0, color='black', linestyle='-', alpha=0.3)
    plt.xticks(rotation=45, ha='right')
    plt.title('Correlation Between Variant Counts and Sterol Levels')
    plt.tight_layout()
    plt.savefig('results/sterol_analysis/correlation/variant_sterol_correlation.png')
    
    # Ergosterol content vs variant counts plot
    ergosterol_df = sterol_df[sterol_df['sterol'] == 'Ergosterol'].copy()
    ergosterol_df['genomic_treatment'] = ergosterol_df['sample'].apply(map_sample_to_treatment)
    
    # Example variant count data (placeholder)
    variant_data = {
        'WT-CTRL': 4,
        'WT-37': 12, 
        'WTA': 12,
        'STC': 16,
        'CAS': 16
    }
    
    ergosterol_df['variant_count'] = ergosterol_df['genomic_treatment'].map(variant_data)
    
    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        x='variant_count', 
        y='concentration',
        hue='genomic_treatment',
        size='concentration',
        data=ergosterol_df,
        alpha=0.8
    )
    
    plt.title('Ergosterol Content vs. Variant Count')
    plt.xlabel('Variant Count')
    plt.ylabel('Ergosterol Concentration')
    plt.savefig('results/sterol_analysis/correlation/ergosterol_vs_variants.png')
    
    # Create heatmap showing sterol profiles vs. adaptation types
    sterol_by_adaptation = sterol_df.copy()
    
    # Map to adaptation type
    adaptation_map = {
        'CAS_5_37C': 'Temperature',
        'CAS_55_37C': 'Temperature',
        'STC_5': 'Low Oxygen',
        'STC_55': 'Low Oxygen',
        'WT_5_37C': 'None',
        'WT_5_MA': 'Low Oxygen',
        'WT_55_37C': 'Temperature',
        'WT_55_MA': 'Low Oxygen'
    }
    
    sterol_by_adaptation['adaptation'] = sterol_by_adaptation['sample'].map(adaptation_map)
    
    # Create pivot table
    pivot_data = sterol_by_adaptation.pivot_table(
        index='sterol',
        columns='adaptation',
        values='concentration',
        aggfunc='mean'
    )
    
    # Visualize as heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(pivot_data, annot=True, cmap='viridis')
    plt.title('Sterol Profile by Adaptation Type')
    plt.savefig('results/sterol_analysis/correlation/sterol_by_adaptation_heatmap.png')