#!/usr/bin/env python3
"""
Core sterol analysis script for Yeast MSA project.
This script performs comparative analysis of sterol profiles across treatments.

Input:
- sterol_data_processed.csv: Processed sterol data

Output:
- Various statistical analysis files in results/sterol_analysis/comparative/
- Visualization files in results/sterol_analysis/visualizations/
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import scipy.stats

# Set paths
INPUT_FILE = 'sterol_data/sterol_data_processed.csv'
RESULTS_DIR = 'results/sterol_analysis'
COMP_DIR = f'{RESULTS_DIR}/comparative'
VIS_DIR = f'{RESULTS_DIR}/visualizations'

def ensure_directories():
    """Ensure all required directories exist."""
    for directory in [COMP_DIR, VIS_DIR]:
        Path(directory).mkdir(parents=True, exist_ok=True)

def load_processed_data(file_path=INPUT_FILE):
    """Load processed sterol data."""
    print(f"Loading processed sterol data from {file_path}")
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Processed sterol data file not found: {file_path}")
    
    df = pd.read_csv(file_path)
    print(f"Loaded {len(df)} processed sterol measurements")
    return df

def perform_statistical_tests(df):
    """Perform statistical tests on sterol data."""
    results = {}
    
    # 1. Adaptation type comparison (t-test for ergosterol)
    temp_erg = df[(df['sterol'] == 'Ergosterol') & (df['adaptation_type'] == 'Temperature')]['concentration']
    low_ox_erg = df[(df['sterol'] == 'Ergosterol') & (df['adaptation_type'] == 'Low Oxygen')]['concentration']
    
    if len(temp_erg) > 0 and len(low_ox_erg) > 0:
        t_stat, p_val = scipy.stats.ttest_ind(temp_erg, low_ox_erg, equal_var=False)
        results['ergosterol_by_adaptation'] = {
            'test': 'Independent t-test',
            'statistic': t_stat,
            'p_value': p_val,
            'temperature_mean': temp_erg.mean(),
            'low_oxygen_mean': low_ox_erg.mean(),
            'significant': p_val < 0.05
        }
    
    # 2. Gene modification comparison (t-test for ergosterol)
    mod_erg = df[(df['sterol'] == 'Ergosterol') & (df['gene_modified'])]['concentration']
    nonmod_erg = df[(df['sterol'] == 'Ergosterol') & (~df['gene_modified'])]['concentration']
    
    if len(mod_erg) > 0 and len(nonmod_erg) > 0:
        t_stat, p_val = scipy.stats.ttest_ind(mod_erg, nonmod_erg, equal_var=False)
        results['ergosterol_by_gene_modification'] = {
            'test': 'Independent t-test',
            'statistic': t_stat,
            'p_value': p_val,
            'modified_mean': mod_erg.mean(),
            'nonmodified_mean': nonmod_erg.mean(),
            'significant': p_val < 0.05
        }
    
    # 3. Treatment comparison (ANOVA for ergosterol)
    treatments = df['treatment'].unique()
    if len(treatments) > 2:
        anova_groups = []
        for treatment in treatments:
            group = df[(df['sterol'] == 'Ergosterol') & (df['treatment'] == treatment)]['concentration']
            if len(group) > 0:
                anova_groups.append(group)
        
        if len(anova_groups) > 1:
            f_stat, p_val = scipy.stats.f_oneway(*anova_groups)
            results['ergosterol_by_treatment'] = {
                'test': 'One-way ANOVA',
                'statistic': f_stat,
                'p_value': p_val,
                'significant': p_val < 0.05
            }
            
            # If ANOVA is significant, perform post-hoc tests
            if p_val < 0.05:
                # Flatten data for tukey test
                treatment_data = df[df['sterol'] == 'Ergosterol'].copy()
                from statsmodels.stats.multicomp import pairwise_tukeyhsd
                tukey = pairwise_tukeyhsd(
                    treatment_data['concentration'], 
                    treatment_data['treatment'],
                    alpha=0.05
                )
                # Convert tukey results to DataFrame
                tukey_df = pd.DataFrame(
                    data=tukey._results_table.data[1:],
                    columns=tukey._results_table.data[0]
                )
                tukey_df.to_csv(f'{COMP_DIR}/tukey_treatment_comparison.csv', index=False)
                results['tukey_results'] = 'Saved to tukey_treatment_comparison.csv'
    
    # 4. Generation comparison (paired t-test where possible)
    gen_tests = {}
    for treatment in treatments:
        gen5 = df[(df['sterol'] == 'Ergosterol') & (df['treatment'] == treatment) & (df['generation'] == '5')]['concentration']
        gen55 = df[(df['sterol'] == 'Ergosterol') & (df['treatment'] == treatment) & (df['generation'] == '55')]['concentration']
        
        if len(gen5) > 0 and len(gen55) > 0:
            # Use independent t-test since we don't have paired data
            t_stat, p_val = scipy.stats.ttest_ind(gen5, gen55, equal_var=False)
            gen_tests[treatment] = {
                'test': 'Independent t-test',
                'statistic': t_stat,
                'p_value': p_val,
                'gen5_mean': gen5.mean(),
                'gen55_mean': gen55.mean(),
                'change': gen55.mean() - gen5.mean(),
                'percent_change': (gen55.mean() - gen5.mean()) / gen5.mean() * 100 if gen5.mean() != 0 else np.nan,
                'significant': p_val < 0.05
            }
    
    results['generation_tests'] = gen_tests
    
    # 5. Sterol diversity analysis
    diversity_stats = {}
    for treatment in df['treatment'].unique():
        treatment_sterols = df[df['treatment'] == treatment]['sterol'].unique()
        diversity_stats[treatment] = {
            'sterol_count': len(treatment_sterols),
            'sterols': list(treatment_sterols)
        }
    
    # Compare sterol diversity by gene modification
    mod_diversity = np.mean([stat_info['sterol_count'] for treatment, stat_info in diversity_stats.items() 
                          if any(df[(df['treatment'] == treatment) & df['gene_modified']].index)])
    nonmod_diversity = np.mean([stat_info['sterol_count'] for treatment, stat_info in diversity_stats.items() 
                             if any(df[(df['treatment'] == treatment) & ~df['gene_modified']].index)])
    
    if not np.isnan(mod_diversity) and not np.isnan(nonmod_diversity):
        results['sterol_diversity_by_modification'] = {
            'modified_mean_diversity': mod_diversity,
            'nonmodified_mean_diversity': nonmod_diversity,
            'difference': mod_diversity - nonmod_diversity,
            'fold_difference': mod_diversity / nonmod_diversity if nonmod_diversity > 0 else np.nan
        }
    
    # Compare sterol diversity by adaptation type
    temp_diversity = np.mean([stat_info['sterol_count'] for treatment, stat_info in diversity_stats.items() 
                           if any(df[(df['treatment'] == treatment) & (df['adaptation_type'] == 'Temperature')].index)])
    lowox_diversity = np.mean([stat_info['sterol_count'] for treatment, stat_info in diversity_stats.items() 
                            if any(df[(df['treatment'] == treatment) & (df['adaptation_type'] == 'Low Oxygen')].index)])
    
    if not np.isnan(temp_diversity) and not np.isnan(lowox_diversity):
        results['sterol_diversity_by_adaptation'] = {
            'temperature_mean_diversity': temp_diversity,
            'low_oxygen_mean_diversity': lowox_diversity,
            'difference': temp_diversity - lowox_diversity,
            'fold_difference': temp_diversity / lowox_diversity if lowox_diversity > 0 else np.nan
        }
    
    # 6. Analysis of unique sterols
    # Get all sterols present in each treatment
    treatment_sterols = {}
    for treatment in df['treatment'].unique():
        treatment_sterols[treatment] = set(df[df['treatment'] == treatment]['sterol'].unique())
    
    # Identify unique sterols for each treatment
    unique_sterols = {}
    for treatment, sterols in treatment_sterols.items():
        other_treatments = set(df['treatment'].unique()) - {treatment}
        other_sterols = set()
        for other in other_treatments:
            other_sterols.update(treatment_sterols[other])
        
        unique_sterols[treatment] = sterols - other_sterols
    
    results['unique_sterols'] = unique_sterols
    results['treatment_diversity'] = diversity_stats
    
    # Save statistical test results
    with open(f'{COMP_DIR}/statistical_test_results.txt', 'w') as f:
        f.write("# Statistical Tests for Sterol Analysis\n\n")
        
        f.write("## Ergosterol by Adaptation Type\n")
        if 'ergosterol_by_adaptation' in results:
            r = results['ergosterol_by_adaptation']
            f.write(f"Test: {r['test']}\n")
            f.write(f"Statistic: {r['statistic']:.4f}\n")
            f.write(f"P-value: {r['p_value']:.4f}\n")
            f.write(f"Temperature Mean: {r['temperature_mean']:.2f}\n")
            f.write(f"Low Oxygen Mean: {r['low_oxygen_mean']:.2f}\n")
            f.write(f"Significant: {'Yes' if r['significant'] else 'No'}\n\n")
        else:
            f.write("Not enough data for this test.\n\n")
        
        f.write("## Ergosterol by Gene Modification\n")
        if 'ergosterol_by_gene_modification' in results:
            r = results['ergosterol_by_gene_modification']
            f.write(f"Test: {r['test']}\n")
            f.write(f"Statistic: {r['statistic']:.4f}\n")
            f.write(f"P-value: {r['p_value']:.4f}\n")
            f.write(f"Modified Mean: {r['modified_mean']:.2f}\n")
            f.write(f"Non-modified Mean: {r['nonmodified_mean']:.2f}\n")
            f.write(f"Significant: {'Yes' if r['significant'] else 'No'}\n\n")
        else:
            f.write("Not enough data for this test.\n\n")
        
        f.write("## Ergosterol by Treatment\n")
        if 'ergosterol_by_treatment' in results:
            r = results['ergosterol_by_treatment']
            f.write(f"Test: {r['test']}\n")
            f.write(f"Statistic: {r['statistic']:.4f}\n")
            f.write(f"P-value: {r['p_value']:.4f}\n")
            f.write(f"Significant: {'Yes' if r['significant'] else 'No'}\n")
            if 'tukey_results' in results:
                f.write(f"Post-hoc test: {results['tukey_results']}\n\n")
        else:
            f.write("Not enough data for this test.\n\n")
        
        f.write("## Generation Comparison Tests\n")
        if results['generation_tests']:
            for treatment, r in results['generation_tests'].items():
                f.write(f"### {treatment}\n")
                f.write(f"Test: {r['test']}\n")
                f.write(f"Statistic: {r['statistic']:.4f}\n")
                f.write(f"P-value: {r['p_value']:.4f}\n")
                f.write(f"Gen 5 Mean: {r['gen5_mean']:.2f}\n")
                f.write(f"Gen 55 Mean: {r['gen55_mean']:.2f}\n")
                f.write(f"Change: {r['change']:.2f}\n")
                if not np.isnan(r['percent_change']):
                    f.write(f"Percent Change: {r['percent_change']:.2f}%\n")
                f.write(f"Significant: {'Yes' if r['significant'] else 'No'}\n\n")
        else:
            f.write("Not enough data for these tests.\n\n")
            
        # Add sterol diversity analysis results
        f.write("## Sterol Diversity Analysis\n\n")
        
        f.write("### Diversity by Treatment\n")
        for treatment, stats in results['treatment_diversity'].items():
            f.write(f"**{treatment}**: {stats['sterol_count']} sterols\n")
            f.write(f"- Sterols: {', '.join(stats['sterols'])}\n\n")
        
        if 'sterol_diversity_by_modification' in results:
            r = results['sterol_diversity_by_modification']
            f.write("### Diversity by Gene Modification\n")
            f.write(f"- Modified strains: {r['modified_mean_diversity']:.2f} sterols (mean)\n")
            f.write(f"- Non-modified strains: {r['nonmodified_mean_diversity']:.2f} sterols (mean)\n")
            f.write(f"- Difference: {r['difference']:.2f} more sterols in modified strains\n")
            if not np.isnan(r['fold_difference']):
                f.write(f"- Fold difference: {r['fold_difference']:.2f}x\n\n")
        
        if 'sterol_diversity_by_adaptation' in results:
            r = results['sterol_diversity_by_adaptation']
            f.write("### Diversity by Adaptation Type\n")
            f.write(f"- Temperature adaptation: {r['temperature_mean_diversity']:.2f} sterols (mean)\n")
            f.write(f"- Low oxygen adaptation: {r['low_oxygen_mean_diversity']:.2f} sterols (mean)\n")
            f.write(f"- Difference: {r['difference']:.2f}\n")
            if not np.isnan(r['fold_difference']):
                f.write(f"- Fold difference: {r['fold_difference']:.2f}x\n\n")
        
        # Add unique sterol analysis
        f.write("## Unique Sterol Analysis\n\n")
        for treatment, unique in results['unique_sterols'].items():
            if unique:
                f.write(f"### {treatment}\n")
                f.write(f"Unique sterols: {', '.join(unique)}\n\n")
            else:
                f.write(f"### {treatment}\n")
                f.write(f"No unique sterols detected.\n\n")
    
    # Also save as CSV for easier programmatic access
    results_df = []
    
    if 'ergosterol_by_adaptation' in results:
        r = results['ergosterol_by_adaptation']
        results_df.append({
            'comparison': 'Ergosterol by Adaptation',
            'test': r['test'],
            'statistic': r['statistic'],
            'p_value': r['p_value'],
            'group1': 'Temperature',
            'group1_mean': r['temperature_mean'],
            'group2': 'Low Oxygen',
            'group2_mean': r['low_oxygen_mean'],
            'significant': r['significant']
        })
    
    if 'ergosterol_by_gene_modification' in results:
        r = results['ergosterol_by_gene_modification']
        results_df.append({
            'comparison': 'Ergosterol by Gene Modification',
            'test': r['test'],
            'statistic': r['statistic'],
            'p_value': r['p_value'],
            'group1': 'Modified',
            'group1_mean': r['modified_mean'],
            'group2': 'Non-modified',
            'group2_mean': r['nonmodified_mean'],
            'significant': r['significant']
        })
    
    if results['generation_tests']:
        for treatment, r in results['generation_tests'].items():
            results_df.append({
                'comparison': f'Generation Effect in {treatment}',
                'test': r['test'],
                'statistic': r['statistic'],
                'p_value': r['p_value'],
                'group1': 'Gen 5',
                'group1_mean': r['gen5_mean'],
                'group2': 'Gen 55',
                'group2_mean': r['gen55_mean'],
                'change': r['change'],
                'percent_change': r['percent_change'],
                'significant': r['significant']
            })
    
    pd.DataFrame(results_df).to_csv(f'{COMP_DIR}/statistical_results_summary.csv', index=False)
    
    # Save sterol diversity analysis to CSV
    diversity_df = []
    
    for treatment, stat_info in results['treatment_diversity'].items():
        diversity_df.append({
            'treatment': treatment,
            'sterol_count': stat_info['sterol_count'],
            'sterols': ', '.join(stat_info['sterols'])
        })
    
    pd.DataFrame(diversity_df).to_csv(f'{COMP_DIR}/sterol_diversity_by_treatment.csv', index=False)
    
    # Save unique sterols analysis to CSV
    unique_df = []
    
    for treatment, unique_set in results['unique_sterols'].items():
        unique_df.append({
            'treatment': treatment,
            'unique_sterol_count': len(unique_set),
            'unique_sterols': ', '.join(unique_set) if unique_set else 'None'
        })
    
    pd.DataFrame(unique_df).to_csv(f'{COMP_DIR}/unique_sterols_by_treatment.csv', index=False)
    
    print(f"Statistical test results saved to {COMP_DIR}")
    return results

def calculate_fold_changes(df):
    """Calculate fold changes between conditions."""
    # Calculate fold changes for ergosterol
    erg_df = df[df['sterol'] == 'Ergosterol'].copy()
    
    # Create pivot table for easier comparison
    erg_pivot = erg_df.pivot_table(
        index=['treatment', 'generation'],
        values='concentration',
        aggfunc='mean'
    ).reset_index()
    
    # Calculate fold changes
    fold_changes = []
    
    # Temperature vs Low Oxygen adaptation fold changes
    temp_samples = erg_df[erg_df['adaptation_type'] == 'Temperature']
    lowox_samples = erg_df[erg_df['adaptation_type'] == 'Low Oxygen']
    
    if not temp_samples.empty and not lowox_samples.empty:
        temp_mean = temp_samples['concentration'].mean()
        lowox_mean = lowox_samples['concentration'].mean()
        
        if lowox_mean != 0:
            fold_change = temp_mean / lowox_mean
            fold_changes.append({
                'comparison': 'Temperature vs Low Oxygen',
                'group1': 'Temperature',
                'group1_mean': temp_mean,
                'group2': 'Low Oxygen',
                'group2_mean': lowox_mean,
                'fold_change': fold_change,
                'log2_fold_change': np.log2(fold_change) if fold_change > 0 else np.nan
            })
    
    # Gene-modified vs non-modified fold changes
    mod_samples = erg_df[erg_df['gene_modified']]
    nonmod_samples = erg_df[~erg_df['gene_modified']]
    
    if not mod_samples.empty and not nonmod_samples.empty:
        mod_mean = mod_samples['concentration'].mean()
        nonmod_mean = nonmod_samples['concentration'].mean()
        
        if nonmod_mean != 0:
            fold_change = mod_mean / nonmod_mean
            fold_changes.append({
                'comparison': 'Modified vs Non-modified',
                'group1': 'Modified',
                'group1_mean': mod_mean,
                'group2': 'Non-modified',
                'group2_mean': nonmod_mean,
                'fold_change': fold_change,
                'log2_fold_change': np.log2(fold_change) if fold_change > 0 else np.nan
            })
    
    # Generation fold changes (within treatments)
    for treatment in erg_df['treatment'].unique():
        gen5 = erg_df[(erg_df['treatment'] == treatment) & (erg_df['generation'] == '5')]['concentration']
        gen55 = erg_df[(erg_df['treatment'] == treatment) & (erg_df['generation'] == '55')]['concentration']
        
        if not gen5.empty and not gen55.empty:
            gen5_mean = gen5.mean()
            gen55_mean = gen55.mean()
            
            if gen5_mean != 0:
                fold_change = gen55_mean / gen5_mean
                fold_changes.append({
                    'comparison': f'Generation Effect in {treatment}',
                    'group1': 'Gen 55',
                    'group1_mean': gen55_mean,
                    'group2': 'Gen 5',
                    'group2_mean': gen5_mean,
                    'fold_change': fold_change,
                    'log2_fold_change': np.log2(fold_change) if fold_change > 0 else np.nan
                })
    
    # Treatment comparisons (pairwise)
    treatments = erg_df['treatment'].unique()
    for i, t1 in enumerate(treatments):
        for t2 in treatments[i+1:]:
            t1_samples = erg_df[erg_df['treatment'] == t1]['concentration']
            t2_samples = erg_df[erg_df['treatment'] == t2]['concentration']
            
            if not t1_samples.empty and not t2_samples.empty:
                t1_mean = t1_samples.mean()
                t2_mean = t2_samples.mean()
                
                if t2_mean != 0:
                    fold_change = t1_mean / t2_mean
                    fold_changes.append({
                        'comparison': f'{t1} vs {t2}',
                        'group1': t1,
                        'group1_mean': t1_mean,
                        'group2': t2,
                        'group2_mean': t2_mean,
                        'fold_change': fold_change,
                        'log2_fold_change': np.log2(fold_change) if fold_change > 0 else np.nan
                    })
    
    # Save fold changes to CSV
    fold_changes_df = pd.DataFrame(fold_changes)
    fold_changes_df.to_csv(f'{COMP_DIR}/ergosterol_fold_changes.csv', index=False)
    
    # Create a summary text file
    with open(f'{COMP_DIR}/fold_change_analysis.txt', 'w') as f:
        f.write("# Ergosterol Fold Change Analysis\n\n")
        
        f.write("## Summary of Key Fold Changes\n")
        for fc in fold_changes:
            f.write(f"### {fc['comparison']}\n")
            f.write(f"- {fc['group1']}: {fc['group1_mean']:.2f}\n")
            f.write(f"- {fc['group2']}: {fc['group2_mean']:.2f}\n")
            f.write(f"- Fold change ({fc['group1']}/{fc['group2']}): {fc['fold_change']:.2f}\n")
            if 'log2_fold_change' in fc and not pd.isna(fc['log2_fold_change']):
                f.write(f"- Log2 fold change: {fc['log2_fold_change']:.2f}\n")
            f.write("\n")
    
    print(f"Fold change analysis saved to {COMP_DIR}")
    return fold_changes_df

def create_comparative_visualizations(df, stats_results, fold_changes):
    """Create visualizations for comparative analysis."""
    # Set up plotting style
    sns.set(style="whitegrid")
    plt.rcParams['figure.figsize'] = (12, 8)
    
    # 1. Statistical comparison plot
    if 'ergosterol_by_adaptation' in stats_results:
        plt.figure(figsize=(14, 8))
        
        # Adaptation type comparison
        sns.barplot(
            data=df[df['sterol'] == 'Ergosterol'],
            x='adaptation_type',
            y='concentration',
            palette='viridis',
            capsize=0.1,
            errwidth=1.5
        )
        
        # Add p-value annotation
        r = stats_results['ergosterol_by_adaptation']
        p_val = r['p_value']
        
        if p_val < 0.001:
            p_text = 'p < 0.001'
        elif p_val < 0.01:
            p_text = 'p < 0.01'
        elif p_val < 0.05:
            p_text = 'p < 0.05'
        else:
            p_text = f'p = {p_val:.3f}'
        
        plt.title(f"Ergosterol Levels by Adaptation Type ({p_text})")
        plt.ylabel("Ergosterol Concentration")
        plt.tight_layout()
        plt.savefig(f'{VIS_DIR}/statistical_adaptation_comparison.png', dpi=300)
        plt.close()
    
    # 2. Fold change plot
    plt.figure(figsize=(12, 6))
    
    # Sort fold changes for better visualization
    sorted_fc = fold_changes.sort_values('fold_change', ascending=False)
    
    # Color bars based on fold change magnitude
    colors = []
    for fc in sorted_fc['fold_change']:
        if fc >= 2:
            colors.append('forestgreen')
        elif fc >= 1:
            colors.append('lightgreen')
        elif fc >= 0.5:
            colors.append('lightsalmon')
        else:
            colors.append('firebrick')
    
    plt.bar(
        range(len(sorted_fc)),
        sorted_fc['fold_change'],
        color=colors
    )
    
    # Add line at fold change = 1
    plt.axhline(y=1, color='black', linestyle='--', alpha=0.7)
    
    plt.xticks(
        range(len(sorted_fc)),
        sorted_fc['comparison'],
        rotation=45,
        ha='right'
    )
    plt.title("Ergosterol Fold Changes Between Conditions")
    plt.ylabel("Fold Change")
    plt.tight_layout()
    plt.savefig(f'{VIS_DIR}/ergosterol_fold_changes.png', dpi=300)
    plt.close()
    
    # 3. Generation effect plot
    plt.figure(figsize=(14, 8))
    
    # Prepare generation data
    gen_data = []
    for treatment in df['treatment'].unique():
        for gen in ['5', '55']:
            gen_samples = df[(df['treatment'] == treatment) & 
                           (df['generation'] == gen) & 
                           (df['sterol'] == 'Ergosterol')]
            
            if not gen_samples.empty:
                for _, row in gen_samples.iterrows():
                    gen_data.append({
                        'treatment': treatment,
                        'generation': gen,
                        'concentration': row['concentration']
                    })
    
    if gen_data:
        gen_df = pd.DataFrame(gen_data)
        
        # Check if we have enough data for a visualization
        if len(gen_df) > 1:
            # Create connected scatter plot
            plt.figure(figsize=(12, 6))
            
            # Group by treatment and generation
            gen_means = gen_df.groupby(['treatment', 'generation'])['concentration'].mean().reset_index()
            
            # Create a pivot table for easier plotting
            gen_pivot = gen_means.pivot(index='treatment', columns='generation', values='concentration')
            
            # Plot connected points
            for treatment in gen_pivot.index:
                if '5' in gen_pivot.columns and '55' in gen_pivot.columns:
                    values = [gen_pivot.loc[treatment, '5'], gen_pivot.loc[treatment, '55']]
                    plt.plot([0, 1], values, 'o-', linewidth=2, markersize=10, label=treatment)
            
            plt.xticks([0, 1], ['Generation 5', 'Generation 55'])
            plt.legend(title='Treatment')
            plt.title("Ergosterol Changes Across Generations")
            plt.ylabel("Ergosterol Concentration")
            plt.tight_layout()
            plt.savefig(f'{VIS_DIR}/generation_effect.png', dpi=300)
            plt.close()
    
    # 4. Heatmap of all fold changes
    plt.figure(figsize=(10, 8))
    
    # Create matrix of fold changes
    treatments = sorted(df['treatment'].unique())
    n_treatments = len(treatments)
    fc_matrix = np.ones((n_treatments, n_treatments))
    
    for i, t1 in enumerate(treatments):
        for j, t2 in enumerate(treatments):
            if i != j:
                # Find this comparison in fold_changes
                fc_row = fold_changes[
                    ((fold_changes['group1'] == t1) & (fold_changes['group2'] == t2)) |
                    ((fold_changes['group1'] == t2) & (fold_changes['group2'] == t1))
                ]
                
                if not fc_row.empty:
                    fc = fc_row.iloc[0]['fold_change']
                    # If t2/t1 instead of t1/t2, invert the fold change
                    if fc_row.iloc[0]['group1'] == t2 and fc_row.iloc[0]['group2'] == t1:
                        fc = 1 / fc
                    fc_matrix[i, j] = fc
    
    # Use log2 for better visualization
    log2_fc_matrix = np.log2(fc_matrix)
    
    # Create custom diverging colormap centered at 0
    cmap = sns.diverging_palette(240, 10, as_cmap=True)
    
    sns.heatmap(
        log2_fc_matrix,
        annot=True,
        fmt=".2f",
        cmap=cmap,
        center=0,
        xticklabels=treatments,
        yticklabels=treatments
    )
    
    plt.title("Log2 Fold Change in Ergosterol Between Treatments")
    plt.tight_layout()
    plt.savefig(f'{VIS_DIR}/treatment_fold_change_heatmap.png', dpi=300)
    plt.close()
    
    # 5. Radar chart for sterol profiles
    # Get unique sterols
    sterols = df['sterol'].unique()
    if len(sterols) >= 3:  # Need at least 3 points for a radar chart
        plt.figure(figsize=(12, 12))
        
        # Number of variables
        N = len(sterols)
        
        # Create angles for radar chart
        angles = np.linspace(0, 2*np.pi, N, endpoint=False).tolist()
        angles += angles[:1]  # Close the loop
        
        # Create figure
        ax = plt.subplot(111, polar=True)
        
        # Draw one axis per variable and add labels
        plt.xticks(angles[:-1], sterols, size=10)
        
        # Draw y-axis labels (assuming scale from 0 to 100%)
        plt.yticks([0.25, 0.5, 0.75, 1], ['25%', '50%', '75%', '100%'], color='grey', size=8)
        plt.ylim(0, 1)
        
        # Group by treatment and calculate mean relative abundance
        treatment_profiles = {}
        for treatment in df['treatment'].unique():
            treatment_df = df[df['treatment'] == treatment]
            profile = []
            for sterol in sterols:
                sterol_rows = treatment_df[treatment_df['sterol'] == sterol]
                if not sterol_rows.empty:
                    mean_rel = sterol_rows['relative_abundance'].mean() / 100.0  # Convert to 0-1 scale
                    profile.append(mean_rel)
                else:
                    profile.append(0)
            
            # Close the loop
            profile.append(profile[0])
            treatment_profiles[treatment] = profile
        
        # Colors for treatments
        colors = plt.cm.viridis(np.linspace(0, 1, len(treatment_profiles)))
        
        # Plot each treatment
        for i, (treatment, profile) in enumerate(treatment_profiles.items()):
            ax.plot(angles, profile, 'o-', linewidth=2, label=treatment, color=colors[i])
            ax.fill(angles, profile, alpha=0.1, color=colors[i])
        
        # Add legend
        plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))
        
        plt.title("Sterol Composition Profiles by Treatment")
        plt.tight_layout()
        plt.savefig(f'{VIS_DIR}/sterol_radar_chart.png', dpi=300)
        plt.close()
    
    # 6. NEW: Sterol diversity bar chart
    plt.figure(figsize=(12, 6))
    
    # Create dataframe of treatment sterol counts
    diversity_data = []
    for treatment, stat_info in stats_results['treatment_diversity'].items():
        diversity_data.append({
            'treatment': treatment,
            'sterol_count': stat_info['sterol_count'],
            # Add relevant metadata for coloring
            'gene_modified': any(df[(df['treatment'] == treatment) & df['gene_modified']].index),
            'adaptation_type': df[df['treatment'] == treatment]['adaptation_type'].iloc[0] if any(df[df['treatment'] == treatment].index) else 'Unknown'
        })
    
    diversity_df = pd.DataFrame(diversity_data)
    
    # Sort by sterol count for better visualization
    diversity_df = diversity_df.sort_values('sterol_count', ascending=False)
    
    # Create grouped bar chart
    ax = sns.barplot(
        data=diversity_df,
        x='treatment',
        y='sterol_count',
        hue='adaptation_type',
        palette='viridis'
    )
    
    # Add gene modification indicator
    for i, row in enumerate(diversity_df.itertuples()):
        if row.gene_modified:
            ax.patches[i].set_edgecolor('red')
            ax.patches[i].set_linewidth(2)
    
    # Add direct labels on the bars
    for i, patch in enumerate(ax.patches):
        height = patch.get_height()
        ax.text(
            patch.get_x() + patch.get_width()/2,
            height + 0.1,
            f"{height:.0f}",
            ha='center'
        )
    
    plt.title("Sterol Diversity by Treatment")
    plt.ylabel("Number of Unique Sterols")
    plt.legend(title="Adaptation Type")
    
    # Add note about red outline for gene modified strains
    plt.figtext(0.5, 0.01, "Red outline indicates gene-modified strain", 
               ha='center', fontsize=10, color='red')
    
    plt.tight_layout()
    plt.savefig(f'{VIS_DIR}/sterol_diversity_barchart.png', dpi=300)
    plt.close()
    
    # 7. NEW: Unique sterol visualization
    plt.figure(figsize=(12, 8))
    
    # Create dataframe for unique sterols
    unique_data = []
    for treatment, unique_set in stats_results['unique_sterols'].items():
        if unique_set:  # Only include treatments with unique sterols
            for sterol in unique_set:
                # Find the concentration of this sterol
                conc = df[(df['treatment'] == treatment) & (df['sterol'] == sterol)]['concentration'].mean()
                
                unique_data.append({
                    'treatment': treatment,
                    'unique_sterol': sterol,
                    'concentration': conc,
                    'gene_modified': any(df[(df['treatment'] == treatment) & df['gene_modified']].index),
                    'adaptation_type': df[df['treatment'] == treatment]['adaptation_type'].iloc[0] if any(df[df['treatment'] == treatment].index) else 'Unknown'
                })
    
    if unique_data:
        unique_df = pd.DataFrame(unique_data)
        
        # Create grouped bar chart
        ax = sns.barplot(
            data=unique_df,
            x='treatment',
            y='concentration',
            hue='unique_sterol',
            palette='Set2'
        )
        
        plt.title("Unique Sterols by Treatment")
        plt.ylabel("Concentration")
        plt.legend(title="Unique Sterol")
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(f'{VIS_DIR}/unique_sterols_barchart.png', dpi=300)
        plt.close()
    
    # 8. NEW: Relative abundance heatmap
    plt.figure(figsize=(14, 10))
    
    # Create pivot table of relative abundances
    rel_abundance = df.pivot_table(
        index='sterol',
        columns='treatment',
        values='relative_abundance',
        aggfunc='mean'
    ).fillna(0)
    
    # Create heatmap
    sns.heatmap(
        rel_abundance,
        cmap='YlGnBu',
        annot=True,
        fmt='.1f',
        linewidths=0.5
    )
    
    plt.title("Relative Abundance (%) of Sterols by Treatment")
    plt.tight_layout()
    plt.savefig(f'{VIS_DIR}/relative_abundance_heatmap.png', dpi=300)
    plt.close()
    
    # 9. NEW: Adaptation type sterol comparison
    # Group by adaptation type and sterol
    adaptation_sterols = df.groupby(['adaptation_type', 'sterol'])['concentration'].mean().reset_index()
    
    # Pivot for easier plotting
    adaptation_pivot = adaptation_sterols.pivot(
        index='sterol',
        columns='adaptation_type',
        values='concentration'
    ).fillna(0)
    
    if not adaptation_pivot.empty:
        plt.figure(figsize=(12, 8))
        
        # Get sterols present in both groups
        common_sterols = []
        temp_only = []
        lowox_only = []
        
        for sterol in adaptation_pivot.index:
            temp_val = adaptation_pivot.loc[sterol, 'Temperature'] if 'Temperature' in adaptation_pivot.columns else 0
            lowox_val = adaptation_pivot.loc[sterol, 'Low Oxygen'] if 'Low Oxygen' in adaptation_pivot.columns else 0
            
            if temp_val > 0 and lowox_val > 0:
                common_sterols.append(sterol)
            elif temp_val > 0:
                temp_only.append(sterol)
            elif lowox_val > 0:
                lowox_only.append(sterol)
        
        # Plot common sterols
        if common_sterols:
            common_df = adaptation_pivot.loc[common_sterols]
            ax = common_df.plot(kind='bar', figsize=(12, 6))
            plt.title("Sterols Present in Both Adaptation Types")
            plt.ylabel("Mean Concentration")
            plt.xticks(rotation=45, ha='right')
            plt.legend(title="Adaptation Type")
            plt.tight_layout()
            plt.savefig(f'{VIS_DIR}/common_sterols_by_adaptation.png', dpi=300)
            plt.close()
        
        # Create bar chart for adaptation-specific sterols
        plt.figure(figsize=(12, 6))
        
        # Combine adaptation-specific sterols
        specific_data = []
        
        for sterol in temp_only:
            specific_data.append({
                'adaptation_type': 'Temperature',
                'sterol': sterol,
                'concentration': adaptation_pivot.loc[sterol, 'Temperature']
            })
            
        for sterol in lowox_only:
            specific_data.append({
                'adaptation_type': 'Low Oxygen',
                'sterol': sterol,
                'concentration': adaptation_pivot.loc[sterol, 'Low Oxygen']
            })
        
        if specific_data:
            specific_df = pd.DataFrame(specific_data)
            
            sns.barplot(
                data=specific_df,
                x='sterol',
                y='concentration',
                hue='adaptation_type',
                palette={'Temperature': 'firebrick', 'Low Oxygen': 'steelblue'}
            )
            
            plt.title("Adaptation-Specific Sterols")
            plt.ylabel("Mean Concentration")
            plt.xticks(rotation=45, ha='right')
            plt.legend(title="Adaptation Type")
            plt.tight_layout()
            plt.savefig(f'{VIS_DIR}/adaptation_specific_sterols.png', dpi=300)
            plt.close()
    
    print(f"Comparative visualizations saved to {VIS_DIR}")

def main():
    """Main analysis function."""
    ensure_directories()
    
    # Load processed data
    df = load_processed_data()
    
    # Perform statistical tests
    stats_results = perform_statistical_tests(df)
    
    # Calculate fold changes
    fold_changes = calculate_fold_changes(df)
    
    # Create comparative visualizations
    create_comparative_visualizations(df, stats_results, fold_changes)
    
    print("Comparative sterol analysis completed.")

if __name__ == "__main__":
    main()