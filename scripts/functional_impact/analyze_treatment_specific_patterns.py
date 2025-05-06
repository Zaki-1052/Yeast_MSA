#!/usr/bin/env python3
"""
analyze_treatment_specific_patterns.py

This script performs comparative analysis of variant patterns across different treatment groups.
It identifies treatment-specific enrichment of variants and generates visualizations.

Usage:
  python analyze_treatment_specific_patterns.py --variants_file <variants_file> --output_dir <output_dir>
"""

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
import warnings
warnings.filterwarnings('ignore')  # Suppress matplotlib warnings

# Define treatment groups
TREATMENT_GROUPS = {
    'CAS': ['CAS-55-1', 'CAS-55-2', 'CAS-55-3'],
    'STC': ['STC-55-1', 'STC-55-2', 'STC-55-3'],
    'WT-37': ['WT-37-55-1', 'WT-37-55-2', 'WT-37-55-3'],
    'WTA': ['WTA-55-1', 'WTA-55-2', 'WTA-55-3']
}

CONTROL_GROUPS = {
    'WT': ['WT-CTRL'],
    'CAS-CTRL': ['CAS-CTRL'],
    'STC-CTRL': ['STC-CTRL']
}

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Analyze treatment-specific variant patterns')
    parser.add_argument('--variants_file', required=True, help='Path to the variants TSV file')
    parser.add_argument('--output_dir', required=True, help='Directory to save results')
    parser.add_argument('--control_group', default='WT', choices=['WT', 'CAS-CTRL', 'STC-CTRL', 'all'], 
                        help='Control group to use for comparison (default: WT)')
    parser.add_argument('--include_controls', action='store_true', 
                        help='Include control samples in the analysis')
    return parser.parse_args()

def load_variants(variants_file):
    """Load variants from TSV file."""
    print(f"Loading variants from {variants_file}")
    variants = pd.read_csv(variants_file, sep='\t')
    print(f"Loaded {len(variants)} variants")
    return variants

def assign_treatment_groups(variants):
    """Assign treatment groups to variants based on sample names."""
    print("Assigning treatment groups")
    
    # Create a mapping from sample to treatment
    sample_to_treatment = {}
    for treatment, samples in TREATMENT_GROUPS.items():
        for sample in samples:
            sample_to_treatment[sample] = treatment
    
    for control_group, samples in CONTROL_GROUPS.items():
        for sample in samples:
            sample_to_treatment[sample] = control_group
    
    # Assign treatment group to each variant
    if 'sample' in variants.columns:
        variants['treatment_group'] = variants['sample'].map(sample_to_treatment)
    else:
        print("Warning: 'sample' column not found, assuming 'treatment' column exists")
    
    # Check if any samples are missing treatment assignments
    if 'treatment_group' in variants.columns:
        missing_treatments = variants[variants['treatment_group'].isna()]
        if len(missing_treatments) > 0:
            print(f"Warning: {len(missing_treatments)} variants have unknown treatment groups")
            print(f"Unique samples without treatment: {missing_treatments['sample'].unique()}")
    
    return variants

def calculate_variant_frequencies(variants):
    """Calculate variant frequencies by gene, region, and treatment."""
    print("Calculating variant frequencies")
    
    # Identify the treatment column
    if 'treatment_group' in variants.columns:
        treatment_column = 'treatment_group'
    elif 'treatment' in variants.columns:
        treatment_column = 'treatment'
    else:
        print("Error: No treatment column found")
        return pd.DataFrame()
    
    # Identify possible grouping columns
    grouping_columns = []
    for col in ['gene_id', 'nearest_gene_name', 'effect', 'impact', 'region', 'erg_gene']:
        if col in variants.columns:
            grouping_columns.append(col)
    
    if 'chrom' in variants.columns and 'position_relative' in variants.columns:
        variants['region'] = variants['chrom'] + '_' + variants['position_relative']
        if 'region' not in grouping_columns:
            grouping_columns.append('region')
    
    if not grouping_columns:
        print("Warning: No suitable grouping columns found")
        if 'chrom' in variants.columns:
            grouping_columns = ['chrom']
        else:
            return pd.DataFrame()
    
    # Calculate variant counts by group and treatment
    result_frames = []
    
    for group_col in grouping_columns:
        # Count variants by group and treatment
        counts = variants.groupby([group_col, treatment_column]).size().reset_index(name='count')
        
        # Add group type column
        counts['group_type'] = group_col
        
        # Append to result frames
        result_frames.append(counts)
    
    # Combine results
    if result_frames:
        result = pd.concat(result_frames, ignore_index=True)
    else:
        result = pd.DataFrame()
    
    return result

def perform_statistical_tests(variants, control_group='WT'):
    """Perform statistical tests to identify significant treatment-specific patterns."""
    print("Performing statistical tests")
    
    # Get treatment groups for the analysis
    treatments = list(TREATMENT_GROUPS.keys())
    
    # Identify the treatment column
    if 'treatment_group' in variants.columns:
        treatment_column = 'treatment_group'
    elif 'treatment' in variants.columns:
        treatment_column = 'treatment'
    else:
        print("Error: No treatment column found")
        return pd.DataFrame()
    
    # Identify possible grouping columns
    grouping_columns = []
    for col in ['gene_id', 'nearest_gene_name', 'effect', 'impact', 'region', 'erg_gene']:
        if col in variants.columns:
            grouping_columns.append(col)
    
    if not grouping_columns:
        print("Warning: No suitable grouping columns found")
        if 'chrom' in variants.columns:
            grouping_columns = ['chrom']
        else:
            return pd.DataFrame()
    
    # Prepare results
    stat_results = []
    
    for group_col in grouping_columns:
        print(f"Analyzing {group_col}")
        
        # Get unique group values
        group_values = variants[group_col].dropna().unique()
        
        for group_value in group_values:
            # For each treatment
            for treatment in treatments:
                # Get variants for this treatment
                treatment_variants = variants[variants[treatment_column] == treatment]
                if len(treatment_variants) == 0:
                    continue
                
                # Get variants in this group for this treatment
                group_treatment_variants = treatment_variants[treatment_variants[group_col] == group_value]
                
                # Get control variants
                if control_group == 'all':
                    control_variants = variants[variants[treatment_column].isin(CONTROL_GROUPS.keys())]
                else:
                    control_variants = variants[variants[treatment_column] == control_group]
                
                if len(control_variants) == 0:
                    continue
                
                # Get variants in this group for control
                group_control_variants = control_variants[control_variants[group_col] == group_value]
                
                # Create contingency table
                # [group_treatment, not_group_treatment]
                # [group_control, not_group_control]
                contingency = np.array([
                    [len(group_treatment_variants), len(treatment_variants) - len(group_treatment_variants)],
                    [len(group_control_variants), len(control_variants) - len(group_control_variants)]
                ])
                
                # Skip if any expected count is too small
                if np.min(contingency) < 1:
                    continue
                
                # Fisher's exact test
                try:
                    odds_ratio, p_value = fisher_exact(contingency)
                except ValueError:
                    odds_ratio, p_value = 0, 1.0
                
                # Calculate enrichment
                treatment_frac = len(group_treatment_variants) / len(treatment_variants) if len(treatment_variants) > 0 else 0
                control_frac = len(group_control_variants) / len(control_variants) if len(control_variants) > 0 else 0
                enrichment = treatment_frac / control_frac if control_frac > 0 else 0
                
                # Store results
                stat_results.append({
                    'group_type': group_col,
                    'group_value': group_value,
                    'treatment': treatment,
                    'treatment_count': len(group_treatment_variants),
                    'control_count': len(group_control_variants),
                    'total_treatment_variants': len(treatment_variants),
                    'total_control_variants': len(control_variants),
                    'treatment_fraction': treatment_frac,
                    'control_fraction': control_frac,
                    'odds_ratio': odds_ratio,
                    'p_value': p_value,
                    'enrichment': enrichment
                })
    
    # Convert to DataFrame
    results_df = pd.DataFrame(stat_results)
    
    # Apply multiple testing correction if we have results
    if len(results_df) > 0:
        # Apply FDR correction
        _, corrected_pvals, _, _ = multipletests(results_df['p_value'].fillna(1), method='fdr_bh')
        results_df['q_value'] = corrected_pvals
        
        # Mark significant results
        results_df['significant'] = results_df['q_value'] < 0.05
    
    return results_df

def analyze_variant_patterns(variants):
    """Analyze patterns of variants across treatments."""
    print("Analyzing variant patterns")
    
    # Identify the treatment column
    if 'treatment_group' in variants.columns:
        treatment_column = 'treatment_group'
    elif 'treatment' in variants.columns:
        treatment_column = 'treatment'
    else:
        print("Error: No treatment column found")
        return {}
    
    # Get unique treatments
    treatments = variants[treatment_column].unique()
    print(f"Found {len(treatments)} unique treatments: {', '.join(treatments)}")
    
    # Initialize patterns dictionary
    patterns = {}
    
    # Count variants by treatment
    patterns['variants_by_treatment'] = variants[treatment_column].value_counts().to_dict()
    
    # If we have gene information, analyze by gene
    if 'gene_id' in variants.columns:
        # Count variants by gene and treatment
        gene_treatment_counts = variants.groupby(['gene_id', treatment_column]).size().unstack(fill_value=0)
        
        # Convert to dictionary format
        patterns['variants_by_gene_treatment'] = {}
        for gene in gene_treatment_counts.index:
            patterns['variants_by_gene_treatment'][gene] = gene_treatment_counts.loc[gene].to_dict()
    
    # If we have effect information, analyze by effect
    if 'effect' in variants.columns:
        # Count variants by effect and treatment
        effect_treatment_counts = variants.groupby(['effect', treatment_column]).size().unstack(fill_value=0)
        
        # Convert to dictionary format
        patterns['variants_by_effect_treatment'] = {}
        for effect in effect_treatment_counts.index:
            patterns['variants_by_effect_treatment'][effect] = effect_treatment_counts.loc[effect].to_dict()
    
    # If we have impact information, analyze by impact
    if 'impact' in variants.columns:
        # Count variants by impact and treatment
        impact_treatment_counts = variants.groupby(['impact', treatment_column]).size().unstack(fill_value=0)
        
        # Convert to dictionary format
        patterns['variants_by_impact_treatment'] = {}
        for impact in impact_treatment_counts.index:
            patterns['variants_by_impact_treatment'][impact] = impact_treatment_counts.loc[impact].to_dict()
    
    # If we have ERG gene information, analyze by ERG gene
    if 'nearest_gene_name' in variants.columns:
        # Count variants by ERG gene and treatment
        erg_treatment_counts = variants.groupby(['nearest_gene_name', treatment_column]).size().unstack(fill_value=0)
        
        # Convert to dictionary format
        patterns['variants_by_erg_treatment'] = {}
        for erg in erg_treatment_counts.index:
            patterns['variants_by_erg_treatment'][erg] = erg_treatment_counts.loc[erg].to_dict()
    
    # If we have erg_gene column, analyze specifically by that
    if 'erg_gene' in variants.columns:
        # Count variants by ERG gene and treatment
        erg_treatment_counts = variants.groupby(['erg_gene', treatment_column]).size().unstack(fill_value=0)
        
        # Convert to dictionary format
        patterns['variants_by_erg_gene_treatment'] = {}
        for erg in erg_treatment_counts.index:
            patterns['variants_by_erg_gene_treatment'][erg] = erg_treatment_counts.loc[erg].to_dict()
    
    return patterns

def generate_visualizations(variants, variant_frequencies, stat_results, patterns, output_dir):
    """Generate visualizations of treatment-specific patterns."""
    print("Generating visualizations")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Identify the treatment column
    if 'treatment_group' in variants.columns:
        treatment_column = 'treatment_group'
    elif 'treatment' in variants.columns:
        treatment_column = 'treatment'
    else:
        print("Error: No treatment column found")
        return
    
    # 1. Create a bar plot of variant counts by treatment
    plt.figure(figsize=(10, 6))
    treatment_counts = variants[treatment_column].value_counts()
    colors = ['blue' if treatment in TREATMENT_GROUPS else 'gray' for treatment in treatment_counts.index]
    treatment_counts.sort_values().plot(kind='barh', color=colors)
    plt.title('Variant Counts by Treatment')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'variant_counts_by_treatment.png'), dpi=300)
    plt.close()
    
    # 2. Create a heatmap of variant counts by gene/group and treatment
    if len(variant_frequencies) > 0:
        # Process each group type
        for group_type in variant_frequencies['group_type'].unique():
            # Filter for this group type
            group_data = variant_frequencies[variant_frequencies['group_type'] == group_type]
            
            # Get the actual column name
            group_col = group_type
            
            # Get top groups by variant count
            top_groups = group_data.groupby(group_col)['count'].sum().nlargest(20).index.tolist()
            
            if top_groups:
                # Filter for top groups
                plot_data = group_data[group_data[group_col].isin(top_groups)]
                
                # Create pivot table
                pivot_data = plot_data.pivot(index=group_col, columns=treatment_column, values='count').fillna(0)
                
                # Sort by total count
                pivot_data['total'] = pivot_data.sum(axis=1)
                pivot_data = pivot_data.sort_values('total', ascending=False)
                pivot_data = pivot_data.drop(columns=['total'])
                
                # Create heatmap
                plt.figure(figsize=(12, len(top_groups) * 0.4 + 2))
                sns.heatmap(pivot_data, cmap='YlGnBu', annot=True, fmt='.0f', linewidths=0.5)
                plt.title(f'Variant Counts by {group_type} and Treatment')
                plt.tight_layout()
                plt.savefig(os.path.join(output_dir, f'{group_type}_treatment_heatmap.png'), dpi=300)
                plt.close()
    
    # 3. Create a heatmap of enrichment values
    if len(stat_results) > 0:
        # Process each group type
        for group_type in stat_results['group_type'].unique():
            # Filter for this group type
            group_results = stat_results[stat_results['group_type'] == group_type]
            
            # Get top groups by significance
            top_sig_groups = group_results.sort_values('q_value').drop_duplicates('group_value')['group_value'].head(20).tolist()
            
            if top_sig_groups:
                # Filter for top groups
                plot_data = group_results[group_results['group_value'].isin(top_sig_groups)]
                
                # Create pivot table
                pivot_data = plot_data.pivot(index='group_value', columns='treatment', values='enrichment').fillna(1)
                
                # Create heatmap
                plt.figure(figsize=(10, len(top_sig_groups) * 0.4 + 2))
                sns.heatmap(pivot_data, cmap='RdBu_r', center=1, annot=True, fmt='.2f', linewidths=0.5)
                plt.title(f'Variant Enrichment by {group_type} and Treatment')
                plt.tight_layout()
                plt.savefig(os.path.join(output_dir, f'{group_type}_enrichment_heatmap.png'), dpi=300)
                plt.close()
    
    # 4. Create a bar plot of significant results by treatment
    if len(stat_results) > 0 and 'significant' in stat_results.columns:
        sig_results = stat_results[stat_results['significant']]
        
        if len(sig_results) > 0:
            treatment_sig_counts = sig_results['treatment'].value_counts()
            
            plt.figure(figsize=(10, 6))
            treatment_sig_counts.plot(kind='bar', color='lightblue')
            plt.title('Number of Significant Associations by Treatment')
            plt.ylabel('Count')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'significant_associations_by_treatment.png'), dpi=300)
            plt.close()
    
    # 5. Create visualizations for specific attributes
    for attr in ['effect', 'impact', 'nearest_gene_name', 'erg_gene']:
        if attr in variants.columns:
            # Count variants by attribute and treatment
            attr_counts = variants.groupby([attr, treatment_column]).size().reset_index(name='count')
            
            # Create pivot table
            pivot_data = attr_counts.pivot(index=attr, columns=treatment_column, values='count').fillna(0)
            
            # Sort by total count
            pivot_data['total'] = pivot_data.sum(axis=1)
            pivot_data = pivot_data.sort_values('total', ascending=False)
            pivot_data = pivot_data.drop(columns=['total'])
            
            # Create bar plot
            plt.figure(figsize=(12, 8))
            pivot_data.plot(kind='bar', stacked=True)
            plt.title(f'Variants by {attr} and Treatment')
            plt.ylabel('Count')
            plt.legend(title='Treatment')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'variants_by_{attr}.png'), dpi=300)
            plt.close()
            
            # Create heatmap
            plt.figure(figsize=(10, len(pivot_data) * 0.4 + 2))
            sns.heatmap(pivot_data, cmap='YlGnBu', annot=True, fmt='.0f', linewidths=0.5)
            plt.title(f'Variant Counts by {attr} and Treatment')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'{attr}_heatmap.png'), dpi=300)
            plt.close()
    
    # 6. Create distribution comparisons by treatment
    if 'aa_change' in variants.columns:
        # Count amino acid changes by treatment
        aa_counts = variants.groupby(['aa_change', treatment_column]).size().reset_index(name='count')
        
        # Get top 20 amino acid changes
        top_aa = aa_counts.groupby('aa_change')['count'].sum().nlargest(20).index.tolist()
        
        if top_aa:
            # Filter for top amino acid changes
            plot_data = aa_counts[aa_counts['aa_change'].isin(top_aa)]
            
            # Create pivot table
            pivot_data = plot_data.pivot(index='aa_change', columns=treatment_column, values='count').fillna(0)
            
            # Create heatmap
            plt.figure(figsize=(10, len(top_aa) * 0.4 + 2))
            sns.heatmap(pivot_data, cmap='YlGnBu', annot=True, fmt='.0f', linewidths=0.5)
            plt.title('Amino Acid Changes by Treatment')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'aa_change_heatmap.png'), dpi=300)
            plt.close()
    
    # 7. Create a PCA plot if possible
    try:
        # Check if we have sample column
        if 'sample' in variants.columns:
            # Create feature columns based on available attributes
            feature_cols = []
            for attr in ['gene_id', 'effect', 'impact', 'nearest_gene_name', 'erg_gene']:
                if attr in variants.columns:
                    feature_cols.append(attr)
            
            if feature_cols:
                for feature in feature_cols:
                    # Create a pivot table of samples by feature values
                    try:
                        pivot = pd.crosstab(variants['sample'], variants[feature])
                        
                        # Apply PCA if possible
                        if pivot.shape[0] >= 3 and pivot.shape[1] >= 3:
                            # Perform PCA
                            pca = PCA(n_components=2)
                            pca_result = pca.fit_transform(pivot)
                            
                            # Create a DataFrame with PCA results
                            pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])
                            pca_df['sample'] = pivot.index
                            
                            # Add treatment information
                            sample_to_treatment = {}
                            for treatment, samples in TREATMENT_GROUPS.items():
                                for sample in samples:
                                    sample_to_treatment[sample] = treatment
                            
                            for control_group, samples in CONTROL_GROUPS.items():
                                for sample in samples:
                                    sample_to_treatment[sample] = control_group
                            
                            pca_df['treatment'] = pca_df['sample'].map(sample_to_treatment)
                            
                            # Create PCA plot
                            plt.figure(figsize=(10, 8))
                            treatments = pca_df['treatment'].dropna().unique()
                            colors = plt.cm.tab10(np.linspace(0, 1, len(treatments)))
                            
                            for i, treatment in enumerate(treatments):
                                subset = pca_df[pca_df['treatment'] == treatment]
                                plt.scatter(subset['PC1'], subset['PC2'], c=[colors[i]], label=treatment, s=100)
                            
                            plt.title(f'PCA of Samples Based on {feature} Patterns')
                            plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%})')
                            plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%})')
                            plt.legend()
                            plt.grid(alpha=0.3)
                            plt.tight_layout()
                            plt.savefig(os.path.join(output_dir, f'pca_by_{feature}.png'), dpi=300)
                            plt.close()
                    except Exception as e:
                        print(f"Could not create PCA plot for {feature}: {e}")
    except Exception as e:
        print(f"Could not create PCA plot: {e}")
    
    print(f"Visualizations saved to {output_dir}")

def save_results(variant_frequencies, stat_results, patterns, output_dir):
    """Save analysis results to output files."""
    print("Saving results")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Save variant frequencies
    if len(variant_frequencies) > 0:
        variant_frequencies.to_csv(os.path.join(output_dir, 'variant_frequencies.tsv'), sep='\t', index=False)
    
    # Save statistical results
    if len(stat_results) > 0:
        stat_results.to_csv(os.path.join(output_dir, 'statistical_results.tsv'), sep='\t', index=False)
        
        # Save significant results separately
        if 'significant' in stat_results.columns:
            significant_results = stat_results[stat_results['significant']]
            if len(significant_results) > 0:
                significant_results.to_csv(os.path.join(output_dir, 'significant_results.tsv'), sep='\t', index=False)
    
    # Save patterns as JSON
    if patterns:
        try:
            import json
            with open(os.path.join(output_dir, 'variant_patterns.json'), 'w') as f:
                json.dump(patterns, f, indent=2, default=str)
        except Exception as e:
            print(f"Could not save patterns as JSON: {e}")
            # Try saving as TSV instead
            try:
                # Convert nested dictionary to dataframe
                pattern_rows = []
                for pattern_type, pattern_data in patterns.items():
                    if isinstance(pattern_data, dict):
                        for key, value in pattern_data.items():
                            if isinstance(value, dict):
                                for subkey, count in value.items():
                                    pattern_rows.append({
                                        'pattern_type': pattern_type,
                                        'key': key,
                                        'treatment': subkey,
                                        'count': count
                                    })
                            else:
                                pattern_rows.append({
                                    'pattern_type': pattern_type,
                                    'key': key,
                                    'count': value
                                })
                
                if pattern_rows:
                    pd.DataFrame(pattern_rows).to_csv(os.path.join(output_dir, 'variant_patterns.tsv'), sep='\t', index=False)
            except:
                print("Could not save patterns in any format")
    
    # Create a summary report
    with open(os.path.join(output_dir, 'analysis_summary.txt'), 'w') as f:
        f.write("Treatment-Specific Variant Analysis Summary\n")
        f.write("========================================\n\n")
        
        f.write("1. Variant Counts by Treatment\n")
        f.write("----------------------------\n")
        for treatment, count in patterns.get('variants_by_treatment', {}).items():
            f.write(f"{treatment}: {count} variants\n")
        f.write("\n")
        
        if len(stat_results) > 0 and 'significant' in stat_results.columns:
            f.write("2. Statistical Analysis\n")
            f.write("---------------------\n")
            f.write(f"Total comparisons: {len(stat_results)}\n")
            f.write(f"Significant associations (q < 0.05): {len(stat_results[stat_results['significant']])}\n\n")
            
            if len(stat_results[stat_results['significant']]) > 0:
                f.write("Top 10 Significant Associations:\n")
                top_results = stat_results[stat_results['significant']].sort_values('q_value').head(10)
                for _, row in top_results.iterrows():
                    f.write(f"- {row['group_value']} in {row['treatment']}: ")
                    f.write(f"Enrichment: {row['enrichment']:.2f}, p-value: {row['p_value']:.2e}, q-value: {row['q_value']:.2e}\n")
                f.write("\n")
        
        f.write("3. Treatment-Specific Patterns\n")
        f.write("---------------------------\n")
        
        if 'variants_by_effect_treatment' in patterns:
            f.write("Variant Effects by Treatment:\n")
            for effect, treatments in patterns['variants_by_effect_treatment'].items():
                f.write(f"- {effect}: ")
                treatment_strs = [f"{treatment}: {count}" for treatment, count in treatments.items()]
                f.write(", ".join(treatment_strs) + "\n")
            f.write("\n")
        
        if 'variants_by_impact_treatment' in patterns:
            f.write("Variant Impacts by Treatment:\n")
            for impact, treatments in patterns['variants_by_impact_treatment'].items():
                f.write(f"- {impact}: ")
                treatment_strs = [f"{treatment}: {count}" for treatment, count in treatments.items()]
                f.write(", ".join(treatment_strs) + "\n")
            f.write("\n")
        
        if 'variants_by_erg_treatment' in patterns:
            f.write("Variants by ERG Gene and Treatment:\n")
            for erg, treatments in patterns['variants_by_erg_treatment'].items():
                f.write(f"- {erg}: ")
                treatment_strs = [f"{treatment}: {count}" for treatment, count in treatments.items()]
                f.write(", ".join(treatment_strs) + "\n")
            f.write("\n")
            
        if 'variants_by_erg_gene_treatment' in patterns:
            f.write("Variants by ERG Gene Column and Treatment:\n")
            for erg, treatments in patterns['variants_by_erg_gene_treatment'].items():
                f.write(f"- {erg}: ")
                treatment_strs = [f"{treatment}: {count}" for treatment, count in treatments.items()]
                f.write(", ".join(treatment_strs) + "\n")
            f.write("\n")
    
    print(f"Results saved to {output_dir}")

def main():
    """Main function to run the analysis."""
    # Parse command line arguments
    args = parse_arguments()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load variants
    variants = load_variants(args.variants_file)
    
    # Check if we need to assign treatment groups
    if 'treatment' not in variants.columns and 'treatment_group' not in variants.columns:
        variants = assign_treatment_groups(variants)
    
    # Calculate variant frequencies
    variant_frequencies = calculate_variant_frequencies(variants)
    
    # Analyze variant patterns
    patterns = analyze_variant_patterns(variants)
    
    # Perform statistical tests
    stat_results = perform_statistical_tests(variants, args.control_group)
    
    # Generate visualizations
    generate_visualizations(variants, variant_frequencies, stat_results, patterns, args.output_dir)
    
    # Save results
    save_results(variant_frequencies, stat_results, patterns, args.output_dir)
    
    print("Analysis complete!")

if __name__ == "__main__":
    main()