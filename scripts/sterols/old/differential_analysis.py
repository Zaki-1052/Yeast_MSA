import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

def perform_statistical_tests(df):
    """Perform statistical tests to identify significant differences in sterol profiles."""
    results = []
    
    # Pairwise comparisons for Ergosterol across treatments
    ergosterol_df = df[df['sterol'] == 'Ergosterol']
    treatments = ergosterol_df['treatment'].unique()
    
    for i in range(len(treatments)):
        for j in range(i+1, len(treatments)):
            t1 = treatments[i]
            t2 = treatments[j]
            
            t1_data = ergosterol_df[ergosterol_df['treatment'] == t1]['concentration']
            t2_data = ergosterol_df[ergosterol_df['treatment'] == t2]['concentration']
            
            t_stat, p_val = stats.ttest_ind(t1_data, t2_data, equal_var=False)
            
            results.append({
                'sterol': 'Ergosterol',
                'treatment1': t1,
                'treatment2': t2,
                't_statistic': t_stat,
                'p_value': p_val,
                'significant': p_val < 0.05
            })
    
    # Compare treatment types (Temperature vs Low Oxygen)
    # Group samples by adaptation type
    temp_adaptation = df[df['treatment'].isin(['WT-37', 'CAS'])]
    oxygen_adaptation = df[df['treatment'].isin(['WTA', 'STC'])]
    
    for sterol in df['sterol'].unique():
        temp_data = temp_adaptation[temp_adaptation['sterol'] == sterol]['concentration']
        oxygen_data = oxygen_adaptation[oxygen_adaptation['sterol'] == sterol]['concentration']
        
        if len(temp_data) > 0 and len(oxygen_data) > 0:
            t_stat, p_val = stats.ttest_ind(temp_data, oxygen_data, equal_var=False)
            
            results.append({
                'sterol': sterol,
                'treatment1': 'Temperature Adaptation',
                'treatment2': 'Low Oxygen Adaptation',
                't_statistic': t_stat,
                'p_value': p_val,
                'significant': p_val < 0.05
            })
    
    return pd.DataFrame(results)

def calculate_fold_changes(df):
    """Calculate fold changes in sterol concentrations between treatments."""
    # Use WT_5_37C as baseline
    baseline = df[df['sample'] == 'WT_5_37C'].copy()
    baseline_vals = baseline.set_index('sterol')['concentration'].to_dict()
    
    results = []
    
    for sample in df['sample'].unique():
        if sample == 'WT_5_37C':
            continue
            
        sample_data = df[df['sample'] == sample]
        
        for _, row in sample_data.iterrows():
            sterol = row['sterol']
            
            if sterol in baseline_vals and baseline_vals[sterol] > 0:
                fold_change = row['concentration'] / baseline_vals[sterol]
                log2_fold_change = np.log2(fold_change)
                
                results.append({
                    'sample': sample,
                    'treatment': row['treatment'],
                    'sterol': sterol,
                    'concentration': row['concentration'],
                    'baseline_concentration': baseline_vals.get(sterol, np.nan),
                    'fold_change': fold_change,
                    'log2_fold_change': log2_fold_change
                })
            
    return pd.DataFrame(results)

def visualize_differential_results(stat_results, fold_changes):
    """Visualize the results of differential analysis."""
    # Volcano plot
    plt.figure(figsize=(10, 8))
    plt.scatter(
        fold_changes['log2_fold_change'],
        -np.log10(stat_results['p_value']),
        alpha=0.7
    )
    plt.axhline(-np.log10(0.05), color='red', linestyle='--')
    plt.axvline(-1, color='blue', linestyle='--')
    plt.axvline(1, color='blue', linestyle='--')
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10 P-value')
    plt.title('Volcano Plot of Sterol Changes')
    plt.savefig('results/sterol_analysis/differential/volcano_plot.png')
    
    # Heatmap of fold changes
    pivot_df = fold_changes.pivot_table(
        index='sterol',
        columns='sample',
        values='fold_change'
    )
    
    plt.figure(figsize=(12, 10))
    sns.heatmap(pivot_df, annot=True, cmap='RdBu_r', center=1)
    plt.title('Fold Changes in Sterol Concentrations Relative to WT_5_37C')
    plt.savefig('results/sterol_analysis/differential/fold_change_heatmap.png')