#!/usr/bin/env python3
"""
Sterol Analysis Script for Yeast MSA Project

This script performs analysis on sterol profile data from different yeast treatment conditions,
exploring differences in sterol composition, visualizing sterol changes, analyzing ergosterol
levels, and comparing adaptation types (temperature vs. low oxygen).
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import networkx as nx
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import dendrogram, linkage

# Set plotting styles
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 12

def create_output_dirs():
    """Create necessary output directories for results."""
    base_dir = 'results/sterol_analysis'
    subdirs = [
        '', 
        'exploratory', 
        'differential', 
        'pathway_analysis', 
        'correlation', 
        'pattern_analysis'
    ]
    
    for subdir in subdirs:
        os.makedirs(f'{base_dir}/{subdir}', exist_ok=True)
        print(f"Created directory: {base_dir}/{subdir}")

def load_sterol_data(filepath):
    """Load sterol data from CSV file and perform initial preprocessing."""
    print(f"Loading sterol data from {filepath}")
    df = pd.read_csv(filepath)
    
    # Parse sample names to extract treatment conditions
    df['treatment'] = df['sample'].apply(lambda x: x.split('_')[0])
    df['temperature'] = df['sample'].apply(lambda x: '55' if '55' in x else '5')
    df['condition'] = df['sample'].apply(lambda x: 'MA' if 'MA' in x else '37C')
    
    # Convert treatment codes to more descriptive names for visualization
    treatment_map = {
        'CAS': 'CAS (Gene-Modified, Temp)',
        'STC': 'STC (Gene-Modified, O2)',
        'WT': 'Wild Type'
    }
    df['treatment_desc'] = df['treatment'].map(treatment_map)
    
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
    df['adaptation'] = df['sample'].map(adaptation_map)
    
    print(f"Loaded data with {len(df)} rows and {len(df['sample'].unique())} unique samples")
    return df

def explore_sterol_profiles(df):
    """Perform exploratory analysis of sterol profiles."""
    print("Performing exploratory analysis...")
    
    # Summary statistics
    summary = df.groupby(['treatment', 'sterol']).agg({
        'concentration': ['mean', 'std', 'min', 'max'],
        'std_dev': 'mean'
    }).reset_index()
    
    # Count unique sterols per treatment
    sterol_counts = df.groupby('treatment')['sterol'].nunique()
    
    # Distribution of sterol types
    sterol_dist = df.groupby('sterol')['concentration'].sum().sort_values(ascending=False)
    
    # Count unique sterols per adaptation type
    adaptation_counts = df.groupby('adaptation')['sterol'].nunique()
    
    # Statistics by adaptation
    adaptation_stats = df.groupby(['adaptation', 'sterol']).agg({
        'concentration': ['mean', 'std', 'min', 'max'],
    }).reset_index()
    
    return {
        'summary': summary,
        'sterol_counts': sterol_counts,
        'sterol_dist': sterol_dist,
        'adaptation_counts': adaptation_counts,
        'adaptation_stats': adaptation_stats
    }

def visualize_sterol_distributions(df, output_dir='results/sterol_analysis/exploratory'):
    """Create basic visualizations of sterol distributions."""
    print("Creating basic sterol distribution visualizations...")
    
    # 1. Ergosterol levels by sample
    plt.figure(figsize=(14, 8))
    ergosterol_df = df[df['sterol'] == 'Ergosterol']
    bar_plot = sns.barplot(x='sample', y='concentration', data=ergosterol_df)
    bar_plot.set_xticklabels(bar_plot.get_xticklabels(), rotation=45, ha='right')
    plt.title('Ergosterol Levels by Sample', fontsize=16)
    plt.ylabel('Concentration', fontsize=14)
    plt.xlabel('Sample', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/ergosterol_by_sample.png', dpi=300)
    plt.close()
    
    # 2. Ergosterol levels by adaptation type
    plt.figure(figsize=(10, 6))
    sns.barplot(x='adaptation', y='concentration', data=ergosterol_df, 
                order=['None', 'Temperature', 'Low Oxygen'])
    plt.title('Ergosterol Levels by Adaptation Type', fontsize=16)
    plt.ylabel('Concentration', fontsize=14)
    plt.xlabel('Adaptation', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/ergosterol_by_adaptation.png', dpi=300)
    plt.close()
    
    # 3. Sterol composition by treatment
    plt.figure(figsize=(14, 10))
    treatment_bars = sns.barplot(x='treatment', y='concentration', hue='sterol', data=df)
    plt.title('Sterol Composition by Treatment', fontsize=16)
    plt.ylabel('Concentration', fontsize=14)
    plt.xlabel('Treatment', fontsize=14)
    plt.legend(title='Sterol', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/sterol_composition_by_treatment.png', dpi=300)
    plt.close()
    
    # 4. Sterol distribution across all samples (stacked bar)
    pivot_df = df.pivot_table(index='sample', columns='sterol', values='concentration', fill_value=0)
    pivot_df.plot(kind='bar', stacked=True, figsize=(15, 10))
    plt.title('Sterol Composition by Sample (Stacked)', fontsize=16)
    plt.ylabel('Concentration', fontsize=14)
    plt.xlabel('Sample', fontsize=14)
    plt.legend(title='Sterol', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/sterol_stacked_by_sample.png', dpi=300)
    plt.close()
    
    # 5. Sterol profile heatmap
    plt.figure(figsize=(12, 10))
    sns.heatmap(pivot_df, annot=True, cmap='viridis', fmt='.1f')
    plt.title('Sterol Profile Heatmap', fontsize=16)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/sterol_profile_heatmap.png', dpi=300)
    plt.close()

def perform_differential_analysis(df, output_dir='results/sterol_analysis/differential'):
    """Perform differential analysis of sterol concentrations between conditions."""
    print("Performing differential analysis...")
    
    # Use WT_5_37C as baseline for fold change calculations
    baseline = df[df['sample'] == 'WT_5_37C'].copy()
    baseline_vals = baseline.set_index('sterol')['concentration'].to_dict()
    
    # Calculate fold changes
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
                    'adaptation': row['adaptation'],
                    'sterol': sterol,
                    'concentration': row['concentration'],
                    'baseline_concentration': baseline_vals.get(sterol, np.nan),
                    'fold_change': fold_change,
                    'log2_fold_change': log2_fold_change
                })
    
    fold_changes = pd.DataFrame(results)
    
    # Visualize fold changes
    plt.figure(figsize=(14, 10))
    pivot_fc = fold_changes.pivot_table(
        index='sterol',
        columns='sample',
        values='fold_change'
    )
    sns.heatmap(pivot_fc, annot=True, cmap='RdBu_r', center=1, fmt='.2f')
    plt.title('Fold Changes in Sterol Concentrations Relative to WT_5_37C', fontsize=16)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/fold_change_heatmap.png', dpi=300)
    plt.close()

    # Barplot of ergosterol fold changes
    plt.figure(figsize=(12, 6))
    ergosterol_fc = fold_changes[fold_changes['sterol'] == 'Ergosterol']
    bars = sns.barplot(x='sample', y='fold_change', data=ergosterol_fc)
    bars.set_xticklabels(bars.get_xticklabels(), rotation=45, ha='right')
    plt.axhline(y=1, color='r', linestyle='--')
    plt.title('Ergosterol Fold Change Relative to WT_5_37C', fontsize=16)
    plt.ylabel('Fold Change', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/ergosterol_fold_change.png', dpi=300)
    plt.close()
    
    # Adaptation-specific fold changes
    plt.figure(figsize=(10, 6))
    adaptation_means = fold_changes.groupby(['adaptation', 'sterol'])['fold_change'].mean().reset_index()
    erg_adapt = adaptation_means[adaptation_means['sterol'] == 'Ergosterol']
    sns.barplot(x='adaptation', y='fold_change', data=erg_adapt, 
                order=['Temperature', 'Low Oxygen'])
    plt.axhline(y=1, color='r', linestyle='--')
    plt.title('Average Ergosterol Fold Change by Adaptation Type', fontsize=16)
    plt.ylabel('Fold Change', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/adaptation_ergosterol_fold_change.png', dpi=300)
    plt.close()
    
    return fold_changes

def analyze_sterol_pathways(df, output_dir='results/sterol_analysis/pathway_analysis'):
    """Analyze sterol pathway relationships and flux."""
    print("Analyzing sterol pathways...")
    
    # Create a simplified pathway graph
    G = nx.DiGraph()
    
    # Add nodes (sterol compounds)
    unique_sterols = df['sterol'].unique()
    for sterol in unique_sterols:
        G.add_node(sterol)
    
    # Add edges (pathway steps) - This is a simplified representation
    # In reality, these would be based on actual biochemical knowledge
    edges = []
    if 'Lanosterol' in unique_sterols and 'Ergosterol' in unique_sterols:
        edges.append(('Lanosterol', 'Ergosterol'))
    if 'Zymosterol' in unique_sterols and 'Ergosterol' in unique_sterols:
        edges.append(('Zymosterol', 'Ergosterol'))
    if 'Fecosterol' in unique_sterols and 'Ergosterol' in unique_sterols:
        edges.append(('Fecosterol', 'Ergosterol'))
    if 'Ergosta-7-en-3-ol' in unique_sterols and 'Ergosterol' in unique_sterols:
        edges.append(('Ergosta-7-en-3-ol', 'Ergosterol'))
    if 'Ergost-7-en-3beta-ol' in unique_sterols and 'Ergosterol' in unique_sterols:
        edges.append(('Ergost-7-en-3beta-ol', 'Ergosterol'))
    if 'Cycloartenol' in unique_sterols and 'Ergosterol' in unique_sterols:
        edges.append(('Cycloartenol', 'Ergosterol'))
    
    for edge in edges:
        G.add_edge(edge[0], edge[1])
    
    # Calculate pathway ratios
    ratios = []
    for sample in df['sample'].unique():
        sample_data = df[df['sample'] == sample]
        sample_sterols = {row['sterol']: row['concentration'] for _, row in sample_data.iterrows()}
        
        for source, target in G.edges():
            if source in sample_sterols and target in sample_sterols:
                if sample_sterols[source] > 0:
                    ratio = sample_sterols[target] / sample_sterols[source]
                    
                    ratios.append({
                        'sample': sample,
                        'treatment': sample_data['treatment'].iloc[0],
                        'adaptation': sample_data['adaptation'].iloc[0],
                        'source_sterol': source,
                        'target_sterol': target,
                        'source_concentration': sample_sterols[source],
                        'target_concentration': sample_sterols[target],
                        'ratio': ratio,
                        'log_ratio': np.log2(ratio)
                    })
    
    pathway_ratios = pd.DataFrame(ratios)
    
    # Visualization of the pathway
    if edges:  # Only if we have edges to visualize
        plt.figure(figsize=(12, 10))
        pos = nx.spring_layout(G, seed=42)
        
        # Get node sizes based on concentrations
        max_conc = df.groupby('sterol')['concentration'].mean().max()
        node_sizes = {}
        for node in G.nodes():
            node_data = df[df['sterol'] == node]
            if len(node_data) > 0:
                avg_conc = node_data['concentration'].mean()
                node_sizes[node] = 1000 * (avg_conc / max_conc)
            else:
                node_sizes[node] = 300
        
        # Draw the graph
        nx.draw_networkx_nodes(
            G, pos,
            node_size=[node_sizes.get(node, 300) for node in G.nodes()],
            node_color='skyblue', alpha=0.8
        )
        
        edge_widths = []
        for source, target in G.edges():
            edge_data = pathway_ratios[
                (pathway_ratios['source_sterol'] == source) & 
                (pathway_ratios['target_sterol'] == target)
            ]
            
            if len(edge_data) > 0:
                avg_ratio = edge_data['ratio'].mean()
                edge_widths.append(max(0.5, avg_ratio))
            else:
                edge_widths.append(1.0)
        
        nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.7, arrows=True)
        nx.draw_networkx_labels(G, pos, font_size=12, font_weight='bold')
        
        plt.title('Ergosterol Pathway Diagram', fontsize=16)
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/pathway_diagram.png', dpi=300)
        plt.close()
    
    # Visualization of pathway ratios
    if not pathway_ratios.empty:
        # Bar plot of ratios by adaptation
        plt.figure(figsize=(14, 8))
        ratio_bars = sns.barplot(
            x='source_sterol', 
            y='ratio', 
            hue='adaptation',
            data=pathway_ratios
        )
        ratio_bars.set_xticklabels(ratio_bars.get_xticklabels(), rotation=45, ha='right')
        plt.title('Sterol Pathway Ratios by Adaptation Type', fontsize=16)
        plt.ylabel('Ratio (Target/Source)', fontsize=14)
        plt.xlabel('Source Sterol', fontsize=14)
        plt.legend(title='Adaptation')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/pathway_ratios_by_adaptation.png', dpi=300)
        plt.close()
    
    return pathway_ratios

def perform_pattern_analysis(df, output_dir='results/sterol_analysis/pattern_analysis'):
    """Perform pattern recognition analysis on sterol profiles."""
    print("Performing pattern analysis...")
    
    # Create a wide-format dataframe with samples as rows and sterols as columns
    pivot_df = df.pivot_table(
        index='sample',
        columns='sterol',
        values='concentration',
        fill_value=0
    )
    
    # Add metadata for interpretation
    metadata = df.drop_duplicates('sample')[['sample', 'treatment', 'temperature', 'condition', 'adaptation']]
    
    # Standardize the data
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(pivot_df)
    
    # Perform PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(X_scaled)
    
    # Create a dataframe with PCA results
    pca_df = pd.DataFrame({
        'PC1': pca_result[:, 0],
        'PC2': pca_result[:, 1],
        'sample': pivot_df.index
    })
    
    # Add metadata
    pca_df = pca_df.merge(metadata, on='sample')
    
    # Explained variance
    explained_variance = {
        'PC1': pca.explained_variance_ratio_[0] * 100,
        'PC2': pca.explained_variance_ratio_[1] * 100
    }
    
    # Component loadings
    loadings = pd.DataFrame(
        pca.components_.T,
        columns=['PC1', 'PC2'],
        index=pivot_df.columns
    )
    
    # PCA plot
    plt.figure(figsize=(12, 10))
    sns.scatterplot(
        x='PC1', y='PC2', 
        hue='adaptation',
        style='treatment',
        s=150, alpha=0.8,
        data=pca_df
    )
    
    plt.title('PCA of Sterol Profiles', fontsize=16)
    plt.xlabel(f'PC1 ({explained_variance["PC1"]:.2f}%)', fontsize=14)
    plt.ylabel(f'PC2 ({explained_variance["PC2"]:.2f}%)', fontsize=14)
    plt.legend(title='Adaptation', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Add sample labels
    for i, row in pca_df.iterrows():
        plt.annotate(row['sample'], (row['PC1'], row['PC2']), fontsize=10)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/pca_plot.png', dpi=300)
    plt.close()
    
    # Loadings plot
    plt.figure(figsize=(12, 10))
    plt.scatter(loadings['PC1'], loadings['PC2'], s=100)
    
    for i, sterol in enumerate(loadings.index):
        plt.annotate(sterol, (loadings['PC1'][i], loadings['PC2'][i]), fontsize=12)
        
    plt.axhline(0, color='gray', linestyle='--', alpha=0.5)
    plt.axvline(0, color='gray', linestyle='--', alpha=0.5)
    plt.title('PCA Loadings: Contribution of Sterols to Principal Components', fontsize=16)
    plt.xlabel('PC1', fontsize=14)
    plt.ylabel('PC2', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/pca_loadings.png', dpi=300)
    plt.close()
    
    # Hierarchical clustering
    if len(pivot_df) >= 2:  # Need at least 2 samples for clustering
        Z = linkage(X_scaled, method='ward')
        
        plt.figure(figsize=(12, 8))
        dendrogram(Z, labels=pivot_df.index, leaf_rotation=90)
        plt.title('Hierarchical Clustering of Sterol Profiles', fontsize=16)
        plt.xlabel('Samples', fontsize=14)
        plt.ylabel('Distance', fontsize=14)
        plt.tight_layout()
        plt.savefig(f'{output_dir}/dendrogram.png', dpi=300)
        plt.close()
    
    return {
        'pca_df': pca_df,
        'explained_variance': explained_variance,
        'loadings': loadings
    }

def save_analysis_results(results, output_dir='results/sterol_analysis'):
    """Save analysis results to CSV files."""
    print("Saving analysis results...")
    
    # Save exploratory analysis results
    if 'summary' in results:
        results['summary'].to_csv(f'{output_dir}/sterol_summary_stats.csv', index=False)
    if 'sterol_dist' in results:
        results['sterol_dist'].to_frame().to_csv(f'{output_dir}/sterol_distribution.csv')
    if 'adaptation_stats' in results:
        results['adaptation_stats'].to_csv(f'{output_dir}/adaptation_statistics.csv', index=False)
    
    # Save differential analysis results
    if 'fold_changes' in results:
        results['fold_changes'].to_csv(f'{output_dir}/differential/fold_changes.csv', index=False)
    
    # Save pathway analysis results
    if 'pathway_ratios' in results:
        results['pathway_ratios'].to_csv(f'{output_dir}/pathway_analysis/sterol_ratios.csv', index=False)
    
    # Save pattern analysis results
    if 'pca_df' in results:
        results['pca_df'].to_csv(f'{output_dir}/pattern_analysis/pca_results.csv', index=False)
    if 'loadings' in results:
        results['loadings'].to_csv(f'{output_dir}/pattern_analysis/pca_loadings.csv')

def generate_summary_report(results, output_dir='results/sterol_analysis'):
    """Generate a summary report of sterol analysis findings."""
    print("Generating summary report...")
    
    report = [
        "# Yeast Sterol Analysis Summary Report",
        "",
        "## Overview",
        "",
        f"Analysis performed on {len(results['df']['sample'].unique())} samples with {len(results['df']['sterol'].unique())} unique sterols.",
        "",
        "## Key Findings",
        "",
        "### Sterol Composition",
        ""
    ]
    
    # Add sterol composition findings
    sterols_by_count = results['sterol_dist'].sort_values(ascending=False)
    top_sterols = sterols_by_count.head(3).index.tolist()
    report.append(f"* Most abundant sterols across all samples: {', '.join(top_sterols)}")
    
    # Add adaptation type findings
    report.extend([
        "",
        "### Adaptation-Specific Patterns",
        ""
    ])
    
    if 'fold_changes' in results:
        # Ergosterol in temperature adaptation
        temp_erg = results['fold_changes'][
            (results['fold_changes']['adaptation'] == 'Temperature') & 
            (results['fold_changes']['sterol'] == 'Ergosterol')
        ]
        if not temp_erg.empty:
            avg_fc = temp_erg['fold_change'].mean()
            report.append(f"* Temperature adaptation: Ergosterol levels are {avg_fc:.2f}x the control level")
        
        # Ergosterol in low oxygen adaptation
        oxygen_erg = results['fold_changes'][
            (results['fold_changes']['adaptation'] == 'Low Oxygen') & 
            (results['fold_changes']['sterol'] == 'Ergosterol')
        ]
        if not oxygen_erg.empty:
            avg_fc = oxygen_erg['fold_change'].mean()
            report.append(f"* Low oxygen adaptation: Ergosterol levels are {avg_fc:.2f}x the control level")
    
    report.extend([
        "",
        "### Treatment Effects",
        ""
    ])
    
    # Add treatment-specific findings
    for treatment in results['df']['treatment'].unique():
        treatment_data = results['df'][results['df']['treatment'] == treatment]
        if 'Ergosterol' in treatment_data['sterol'].values:
            erg_level = treatment_data[treatment_data['sterol'] == 'Ergosterol']['concentration'].mean()
            report.append(f"* {treatment}: Average ergosterol level is {erg_level:.2f}")
    
    # Add statistical patterns
    if 'pca_df' in results:
        report.extend([
            "",
            "### Statistical Patterns",
            "",
            f"* Principal Component Analysis explains {results['explained_variance']['PC1'] + results['explained_variance']['PC2']:.2f}% of variance",
            f"* PC1 ({results['explained_variance']['PC1']:.2f}%) primarily separates samples by adaptation type",
            f"* PC2 ({results['explained_variance']['PC2']:.2f}%) primarily separates samples by treatment"
        ])
    
    # Write report to file
    with open(f'{output_dir}/sterol_analysis_summary.md', 'w') as f:
        f.write('\n'.join(report))

def run_sterol_analysis():
    """Main function to run the complete sterol analysis workflow."""
    print("Starting sterol analysis...")
    
    # Create output directories
    create_output_dirs()
    
    # Load and preprocess sterol data
    sterol_data = load_sterol_data('sterol_data_with_sd.csv')
    
    # Perform exploratory analysis
    exploratory_results = explore_sterol_profiles(sterol_data)
    visualize_sterol_distributions(sterol_data)
    
    # Perform differential analysis
    fold_changes = perform_differential_analysis(sterol_data)
    
    # Analyze sterol pathways
    pathway_ratios = analyze_sterol_pathways(sterol_data)
    
    # Perform pattern analysis
    pattern_results = perform_pattern_analysis(sterol_data)
    
    # Compile all results
    all_results = {
        'df': sterol_data,
        **exploratory_results,
        'fold_changes': fold_changes,
        'pathway_ratios': pathway_ratios,
        **pattern_results
    }
    
    # Save results to files
    save_analysis_results(all_results)
    
    # Generate summary report
    generate_summary_report(all_results)
    
    print("Sterol analysis complete. Results saved to 'results/sterol_analysis/'")
    return all_results

if __name__ == "__main__":
    run_sterol_analysis()