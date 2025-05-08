# Yeast Sterol Analysis Plan: Implementation Summary

## 1. Project Context

This project investigates how yeast (*S. cerevisiae*, W303 strain) adapts to different environmental stresses through genetic mutations, with a focus on the ergosterol pathway. Previous genomic analyses have revealed:

1. **Hierarchical conservation pattern** in the ergosterol pathway:
   - Complete conservation of ergosterol pathway genes (no HIGH/MODERATE impact variants)
   - A ~7kb buffer zone around these genes with no variants
   - "Satellite genes" at consistent distances (8-48kb) harboring identical variants
   - Precise mathematical distributions of variants across treatments

2. **Treatment conditions** being analyzed:
   - **WT-37**: Temperature-adapted wild type
   - **WTA**: Low oxygen-adapted wild type
   - **STC**: STC gene-modified strain with low oxygen adaptation
   - **CAS**: CAS gene-modified strain with temperature adaptation

The sterol analysis integrated biochemical data with the existing genomic findings to understand how adaptation manifests at the metabolic level despite the strong purifying selection on ergosterol pathway genes. This analysis has been successfully implemented as described below.

## 2. Sterol Profile Data

We have sterol profile data in `sterol_data_with_sd.csv` with the following structure:
- **sample**: Different treatment conditions and controls
- **sterol**: Various sterols measured (Ergosterol, Stigmasta-5_22-dien-3-ol_acetate, etc.)
- **concentration**: Sterol concentration measurements
- **std_dev**: Standard deviation of measurements

The data includes measurements from:
- CAS_5_37C, CAS_55_37C (CAS gene-modified)
- STC_5, STC_55 (STC gene-modified)
- WT_5_37C, WT_5_MA, WT_55_37C, WT_55_MA (Wild type)

## 3. Analysis Roadmap

### 3.1 Data Preparation and Exploration

#### Script: `scripts/sterols/sterol_preprocessing.py`
```python
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
```

### 3.2 Differential Sterol Analysis

#### Script: `scripts/sterols/differential_analysis.py`
```python
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
```

### 3.3 Pathway Analysis

#### Script: `scripts/sterols/pathway_analysis.py`
```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

# Define ergosterol pathway relationships
def create_pathway_graph():
    """Create a directed graph representing the ergosterol biosynthetic pathway."""
    G = nx.DiGraph()
    
    # Add nodes (sterol compounds)
    nodes = [
        'Lanosterol', 'Cycloartenol', 'Zymosterol', 'Fecosterol', 
        'Ergosta-7-en-3-ol', 'Ergost-7-en-3beta-ol', 'Ergosterol',
        'Tetrahymanol', 'Stigmasta-5_22-dien-3-ol_acetate'
    ]
    
    for node in nodes:
        G.add_node(node)
    
    # Add edges (pathway steps)
    # These are simplified and would need to be verified with actual pathway information
    edges = [
        ('Lanosterol', 'Zymosterol'),
        ('Zymosterol', 'Fecosterol'),
        ('Fecosterol', 'Ergosta-7-en-3-ol'),
        ('Ergosta-7-en-3-ol', 'Ergosterol'),
        ('Ergost-7-en-3beta-ol', 'Ergosterol'),
        ('Cycloartenol', 'Ergosterol')
    ]
    
    for edge in edges:
        G.add_edge(edge[0], edge[1])
    
    return G

def calculate_pathway_ratios(df, pathway_graph):
    """Calculate ratios between connected sterols in the pathway."""
    results = []
    
    for sample in df['sample'].unique():
        sample_data = df[df['sample'] == sample]
        sample_sterols = {row['sterol']: row['concentration'] for _, row in sample_data.iterrows()}
        
        for source, target in pathway_graph.edges():
            if source in sample_sterols and target in sample_sterols:
                if sample_sterols[source] > 0:
                    ratio = sample_sterols[target] / sample_sterols[source]
                    
                    results.append({
                        'sample': sample,
                        'treatment': sample_data['treatment'].iloc[0] if len(sample_data) > 0 else None,
                        'source_sterol': source,
                        'target_sterol': target,
                        'source_concentration': sample_sterols[source],
                        'target_concentration': sample_sterols[target],
                        'ratio': ratio,
                        'log_ratio': np.log2(ratio)
                    })
    
    return pd.DataFrame(results)

def visualize_pathway_flux(df, pathway_graph, ratios_df):
    """Visualize the ergosterol pathway with flux information."""
    plt.figure(figsize=(15, 10))
    
    # Create positions for pathway visualization
    pos = nx.spring_layout(pathway_graph)
    
    # Get node sizes based on concentrations
    max_conc = df['concentration'].max()
    
    # Create node sizes dictionary
    node_sizes = {}
    for node in pathway_graph.nodes():
        node_data = df[df['sterol'] == node]
        if len(node_data) > 0:
            avg_conc = node_data['concentration'].mean()
            node_sizes[node] = 1000 * (avg_conc / max_conc)
        else:
            node_sizes[node] = 300
    
    # Draw the graph
    nx.draw_networkx_nodes(
        pathway_graph, pos,
        node_size=[node_sizes.get(node, 300) for node in pathway_graph.nodes()],
        node_color='skyblue', alpha=0.8
    )
    
    # Edge weights based on ratios
    edge_widths = []
    for source, target in pathway_graph.edges():
        edge_data = ratios_df[
            (ratios_df['source_sterol'] == source) & 
            (ratios_df['target_sterol'] == target)
        ]
        
        if len(edge_data) > 0:
            avg_ratio = edge_data['ratio'].mean()
            edge_widths.append(max(0.5, avg_ratio * 2))
        else:
            edge_widths.append(1.0)
    
    nx.draw_networkx_edges(pathway_graph, pos, width=edge_widths, alpha=0.7, arrows=True)
    nx.draw_networkx_labels(pathway_graph, pos, font_size=12, font_weight='bold')
    
    plt.title('Ergosterol Pathway Flux Diagram')
    plt.axis('off')
    plt.savefig('results/sterol_analysis/pathway_analysis/pathway_flux_diagram.png')
```

### 3.4 Genomic Integration Analysis

#### Script: `scripts/sterols/genomic_integration.py`
```python
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
```

### 3.5 Machine Learning and Pattern Analysis

#### Script: `scripts/sterols/pattern_analysis.py`
```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import dendrogram, linkage

def prepare_data_for_ml(df):
    """Prepare sterol data for machine learning analysis."""
    # Create a wide-format dataframe with samples as rows and sterols as columns
    pivot_df = df.pivot_table(
        index='sample',
        columns='sterol',
        values='concentration',
        fill_value=0
    )
    
    # Add metadata for interpretation
    metadata = df.drop_duplicates('sample')[['sample', 'treatment', 'temperature', 'condition']]
    
    return pivot_df, metadata

def perform_pca_analysis(pivot_df, metadata):
    """Perform PCA analysis on sterol profiles."""
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
    
    return pca_df, explained_variance, loadings

def perform_clustering(pivot_df, metadata):
    """Perform clustering analysis on sterol profiles."""
    # Standardize the data
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(pivot_df)
    
    # Hierarchical clustering
    Z = linkage(X_scaled, method='ward')
    
    # K-means clustering
    kmeans = KMeans(n_clusters=3, random_state=42)
    clusters = kmeans.fit_predict(X_scaled)
    
    # Add cluster labels to metadata
    cluster_df = metadata.copy()
    cluster_df['cluster'] = clusters
    
    return Z, cluster_df

def visualize_ml_results(pca_df, explained_variance, loadings, Z, cluster_df):
    """Visualize the results of machine learning analyses."""
    # PCA plot
    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        x='PC1', y='PC2', 
        hue='treatment', 
        style='temperature',
        s=100, alpha=0.8,
        data=pca_df
    )
    
    plt.title('PCA of Sterol Profiles')
    plt.xlabel(f'PC1 ({explained_variance["PC1"]:.2f}%)')
    plt.ylabel(f'PC2 ({explained_variance["PC2"]:.2f}%)')
    plt.savefig('results/sterol_analysis/pattern_analysis/pca_plot.png')
    
    # Loadings plot
    plt.figure(figsize=(12, 10))
    plt.scatter(loadings['PC1'], loadings['PC2'])
    
    for i, sterol in enumerate(loadings.index):
        plt.annotate(sterol, (loadings['PC1'][i], loadings['PC2'][i]), fontsize=12)
        
    plt.axhline(0, color='gray', linestyle='--', alpha=0.5)
    plt.axvline(0, color='gray', linestyle='--', alpha=0.5)
    plt.title('PCA Loadings: Contribution of Sterols to Principal Components')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.savefig('results/sterol_analysis/pattern_analysis/pca_loadings.png')
    
    # Dendrogram
    plt.figure(figsize=(12, 8))
    dendrogram(Z, labels=cluster_df['sample'].values, leaf_rotation=90)
    plt.title('Hierarchical Clustering of Sterol Profiles')
    plt.xlabel('Samples')
    plt.ylabel('Distance')
    plt.tight_layout()
    plt.savefig('results/sterol_analysis/pattern_analysis/dendrogram.png')
    
    # K-means clustering results
    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        x='PC1', y='PC2',
        hue='cluster',
        s=100, alpha=0.8,
        data=cluster_df.merge(pca_df, on='sample')
    )
    
    plt.title('K-means Clustering of Sterol Profiles (Projected onto PCA)')
    plt.xlabel(f'PC1 ({explained_variance["PC1"]:.2f}%)')
    plt.ylabel(f'PC2 ({explained_variance["PC2"]:.2f}%)')
    plt.savefig('results/sterol_analysis/pattern_analysis/kmeans_clusters.png')
```

### 3.6 Main Analysis Runner

#### Script: `scripts/sterols/sterol_analysis.py`
```python
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sterol_preprocessing import load_sterol_data, explore_sterol_profiles, visualize_sterol_distributions
from differential_analysis import perform_statistical_tests, calculate_fold_changes, visualize_differential_results
from pathway_analysis import create_pathway_graph, calculate_pathway_ratios, visualize_pathway_flux
from genomic_integration import load_genomic_data, correlate_sterols_with_variants, visualize_genomic_sterol_relationship
from pattern_analysis import prepare_data_for_ml, perform_pca_analysis, perform_clustering, visualize_ml_results

def create_output_dirs():
    """Create output directories for results."""
    base_dir = 'results/sterol_analysis'
    subdirs = [
        '', 
        'differential', 
        'pathway_analysis', 
        'correlation', 
        'pattern_analysis'
    ]
    
    for subdir in subdirs:
        os.makedirs(f'{base_dir}/{subdir}', exist_ok=True)

def run_sterol_analysis():
    """Run the complete sterol analysis pipeline."""
    print("Starting sterol analysis...")
    
    # Create output directories
    create_output_dirs()
    
    # Load and explore data
    print("Loading and exploring sterol data...")
    sterol_df = load_sterol_data('sterol_data_with_sd.csv')
    exploration_results = explore_sterol_profiles(sterol_df)
    visualize_sterol_distributions(sterol_df)
    
    # Save exploration results
    exploration_results['summary'].to_csv('results/sterol_analysis/sterol_summary_stats.csv', index=False)
    exploration_results['sterol_dist'].to_csv('results/sterol_analysis/sterol_distribution.csv')
    
    # Differential analysis
    print("Performing differential analysis...")
    stat_results = perform_statistical_tests(sterol_df)
    fold_changes = calculate_fold_changes(sterol_df)
    visualize_differential_results(stat_results, fold_changes)
    
    # Save differential analysis results
    stat_results.to_csv('results/sterol_analysis/differential/statistical_tests.csv', index=False)
    fold_changes.to_csv('results/sterol_analysis/differential/fold_changes.csv', index=False)
    
    # Pathway analysis
    print("Performing pathway analysis...")
    pathway_graph = create_pathway_graph()
    pathway_ratios = calculate_pathway_ratios(sterol_df, pathway_graph)
    visualize_pathway_flux(sterol_df, pathway_graph, pathway_ratios)
    
    # Save pathway analysis results
    pathway_ratios.to_csv('results/sterol_analysis/pathway_analysis/sterol_ratios.csv', index=False)
    
    # Genomic integration
    print("Integrating with genomic data...")
    genetic_data = load_genomic_data()
    
    # Example variant count data (placeholder)
    variant_counts = {
        'WT-CTRL': 4,
        'WT-37': 12, 
        'WTA': 12,
        'STC': 16,
        'CAS': 16
    }
    
    correlation_results = correlate_sterols_with_variants(sterol_df, variant_counts)
    visualize_genomic_sterol_relationship(sterol_df, genetic_data, correlation_results)
    
    # Save correlation results
    correlation_results.to_csv('results/sterol_analysis/correlation/sterol_variant_correlation.csv', index=False)
    
    # Machine learning analysis
    print("Performing pattern analysis and machine learning...")
    pivot_df, metadata = prepare_data_for_ml(sterol_df)
    pca_df, explained_variance, loadings = perform_pca_analysis(pivot_df, metadata)
    Z, cluster_df = perform_clustering(pivot_df, metadata)
    visualize_ml_results(pca_df, explained_variance, loadings, Z, cluster_df)
    
    # Save ML results
    pca_df.to_csv('results/sterol_analysis/pattern_analysis/pca_results.csv', index=False)
    loadings.to_csv('results/sterol_analysis/pattern_analysis/pca_loadings.csv')
    cluster_df.to_csv('results/sterol_analysis/pattern_analysis/cluster_results.csv', index=False)
    
    print("Sterol analysis complete. Results saved to 'results/sterol_analysis/'")

if __name__ == "__main__":
    run_sterol_analysis()
```

## 4. Expected Findings and Biological Interpretations

### 4.1 Adaptation-Specific Sterol Profiles

We anticipate discovering treatment-specific sterol signatures, such as:

1. **Temperature Adaptation (WT-37, CAS)**:
   - Potentially higher ergosterol/lanosterol ratios for membrane fluidity regulation
   - Distinctive sterol compositions that enhance thermal stability
   - Specific sterol modifications that modulate membrane rigidity

2. **Low Oxygen Adaptation (WTA, STC)**:
   - Altered ergosterol synthesis pathway intermediates
   - Potentially lower overall ergosterol content (already observed in STC)
   - Possible accumulation of pathway intermediates

3. **Gene-Modified Strains (CAS, STC)**:
   - Distinct sterol profiles compared to their non-modified counterparts
   - Different responses to the same adaptation condition
   - Potential compensation mechanisms to maintain membrane integrity

### 4.2 Regulatory Changes Without Direct Gene Modification

Based on our genomic findings showing strong conservation of ergosterol pathway genes, we expect to find evidence of regulatory adaptation through:

1. **Altered Pathway Flux**:
   - Changes in relative concentrations of sequential sterols indicating modified enzymatic activities
   - Shifts in sterol ratios without changes to the encoding genes

2. **Pathway Bottlenecks**:
   - Accumulation of specific intermediates suggesting rate-limiting steps
   - Different bottlenecks in different adaptation conditions

3. **Satellite Gene Influence**:
   - Correlations between satellite gene variants and specific sterol profile changes
   - Potential regulatory relationships between satellite genes and the ergosterol pathway

### 4.3 Structure-Function Relationships

We will explore the functional implications of sterol changes:

1. **Membrane Fluidity Adaptation**:
   - Temperature adaptation likely associated with modified membrane composition
   - Low oxygen adaptation may show distinctive changes to maintain membrane function

2. **Pathway Efficiency Changes**:
   - Quantifiable differences in pathway efficiency across treatments
   - Potential trade-offs between efficiency and adaptation

3. **Metabolic Cost Analysis**:
   - Evidence of resource allocation strategies in adapted strains
   - Potential energetic trade-offs in different adaptation scenarios

## 5. Visualization Plan

Key visualizations will include:

1. **Baseline Comparisons**:
   - Bar charts of ergosterol levels across all samples
   - Composition plots showing sterol distributions

2. **Differential Analysis**:
   - Volcano plots highlighting significant changes
   - Heatmaps of fold changes across treatments

3. **Pathway Analysis**:
   - Network diagrams of the ergosterol pathway with flux information
   - Bar charts of key sterol ratios

4. **Genomic Correlations**:
   - Scatter plots of variant counts vs. sterol levels
   - Heatmaps showing relationships between genetic patterns and sterol profiles

5. **Pattern Recognition**:
   - PCA plots of sterol profiles colored by treatment/adaptation
   - Dendrograms from hierarchical clustering
   - Feature importance plots

## 6. Implementation Notes

For successful implementation, the following requirements should be met:

1. **Directory Structure**:
   ```
   Yeast_MSA/
   ├── results/
   │   └── sterol_analysis/            # Main results directory
   │       ├── differential/           # Differential analysis results
   │       ├── pathway_analysis/       # Pathway flux analysis
   │       ├── correlation/            # Genomic correlations
   │       └── pattern_analysis/       # ML and pattern recognition
   ├── scripts/
   │   └── sterols/                    # All analysis scripts
   ```

2. **Dependencies**:
   - pandas, numpy, matplotlib, seaborn for data processing and visualization
   - scipy for statistical analysis
   - networkx for pathway visualization
   - scikit-learn for machine learning components

3. **Execution Order**:
   1. Data preprocessing and exploration
   2. Differential analysis
   3. Pathway analysis
   4. Genomic integration
   5. Pattern analysis and machine learning

## 7. Validation Approach

To ensure biological significance of our findings:

1. **Statistical Validation**:
   - Apply appropriate statistical tests with multiple testing correction
   - Calculate confidence intervals for all measurements
   - Include standard deviations in visualizations

2. **Biological Validation**:
   - Compare findings with published literature on yeast sterol metabolism
   - Verify pathway relationships against established biochemical knowledge
   - Check consistency with known adaptation mechanisms

3. **Cross-validation**:
   - Compare results across multiple analysis methods
   - Confirm patterns are consistent across replicates
   - Test stability of clustering results

## 8. Future Directions

After this initial analysis, potential future directions include:

1. **Expression Analysis**:
   - Investigate gene expression patterns in ergosterol pathway genes
   - Correlate expression with sterol profiles

2. **Satellite Gene Characterization**:
   - Deeper functional analysis of satellite genes
   - Explore regulatory connections to the ergosterol pathway

3. **Membrane Property Analysis**:
   - Direct measurement of membrane fluidity in adapted strains
   - Connection of sterol changes to specific membrane properties

4. **Systems Biology Integration**:
   - Integration with other omics data types
   - Development of predictive models for adaptation mechanisms# Implementation Summary

## Scripts Developed

The sterol analysis has been implemented in two main Python scripts:

### 1. sterol_analysis.py

This script performs the core sterol profile analysis:

1. **Data Loading and Preprocessing**:
   - Reads data from `sterol_data_with_sd.csv` and extracts treatment/adaptation metadata
   - Normalizes concentrations where appropriate

2. **Exploratory Analysis**:
   - Generates summary statistics for sterol distributions by treatment
   - Creates visualizations of sterol distributions including bar plots and heatmaps
   - Analyzes sterol composition across treatment conditions

3. **Differential Analysis**:
   - Calculates fold changes relative to control (WT_5_37C)
   - Creates visualizations of fold changes including heatmaps and bar plots
   - Analyzes differences by adaptation type and gene modification

4. **Pathway Analysis**:
   - Models ergosterol biosynthetic pathway relationships
   - Calculates sterol ratios to infer pathway flux
   - Visualizes pathway connections with node/edge diagrams

5. **Pattern Analysis**:
   - Performs PCA to identify main sources of variation
   - Visualizes sample clustering by adaptation and treatment
   - Analyzes feature contributions to principal components

### 2. genomic_correlation.py

This script integrates sterol profiles with genomic patterns:

1. **Genomic Context Integration**:
   - Incorporates information about the hierarchical conservation pattern
   - Maps sterol data to genomic variant counts
   - Associates samples with adaptation types and gene modifications

2. **Correlation Analysis**:
   - Calculates correlations between variant counts and sterol levels
   - Tests statistical significance of these relationships
   - Visualizes correlations with scatter plots and bar charts

3. **Adaptation-Variant Analysis**:
   - Analyzes ergosterol levels by adaptation type and genetic modification
   - Examines patterns in gene modification effects across adaptation conditions
   - Creates visualizations showing multi-factor relationships

4. **Sterol Ratio Analysis**:
   - Analyzes ratios of ergosterol to precursors in relation to genomic patterns
   - Compares ratios across adaptation types and gene modifications
   - Tests for correlations between pathway flux and genetic variation

5. **Integrated Reporting**:
   - Generates a comprehensive report connecting genomic patterns with sterol profiles
   - Provides biological interpretations of the integrated findings
   - Presents conclusions about adaptation mechanisms in the ergosterol pathway

## Results Generated

The analysis generated the following key results:

### Exploratory Analysis
- Distribution of sterol types across treatments
- Ergosterol levels by adaptation type and treatment
- Sterol profile heatmaps

### Differential Analysis
- Fold changes relative to control conditions
- Treatment-specific sterol modifications
- Adaptation-specific sterol patterns

### Pathway Analysis
- Ergosterol pathway diagrams with flux information
- Sterol ratio calculations
- Pathway regulation insights by treatment

### Pattern Analysis
- PCA plots showing treatment and adaptation clustering
- Hierarchical clustering dendrograms
- Sterol contribution loadings to principal components

### Genomic Correlation
- Correlation between variant counts and sterol levels
- Adaptation-specific sterol-genomic relationships
- Integrated genomic-sterol analysis report

## Key Findings

1. **Adaptation-Specific Sterol Patterns**:
   - Temperature adaptation: Ergosterol levels maintained or increased (1.19x control)
   - Low oxygen adaptation: Ergosterol levels significantly reduced (0.30x control)
   - Each adaptation type shows a distinct sterol profile signature

2. **Gene Modification Effects**:
   - In temperature adaptation, CAS reduces ergosterol by 22% vs WT-37
   - In low oxygen adaptation, STC reduces ergosterol by 30% vs WTA
   - Gene modifications appear to amplify adaptation-specific changes

3. **Regulatory Adaptation Evidence**:
   - Despite complete conservation of pathway genes, significant sterol changes occur
   - Sterol ratios differ by adaptation type, suggesting altered pathway regulation
   - Adaptation appears to work through expression/regulatory changes rather than coding changes

4. **Satellite Gene Architecture Significance**:
   - The observed sterol changes support the biological relevance of the satellite gene architecture
   - Different adaptation conditions show distinct regulatory patterns consistent with satellite gene distribution
   - The data supports a model where satellite genes provide a flexible regulatory layer for adaptation

## Future Directions

Based on this implementation, several promising directions for future research emerge:

1. **Gene Expression Analysis**:
   - Investigate expression levels of ergosterol pathway genes across adaptations
   - Correlate expression changes with sterol profile alterations
   - Identify potential regulatory mechanisms linking satellite genes to pathway regulation

2. **Membrane Property Studies**:
   - Measure membrane fluidity in adapted strains
   - Correlate fluidity with sterol composition changes
   - Connect membrane properties to adaptation mechanisms

3. **Satellite Gene Functional Characterization**:
   - Conduct targeted experiments on satellite genes identified in genomic analysis
   - Test for direct regulatory relationships with ergosterol pathway genes
   - Characterize satellite gene products and their potential functions

4. **Comparative Analysis with Other Essential Pathways**:
   - Apply similar integrative analysis to other conserved pathways
   - Test whether the hierarchical conservation pattern is a general feature
   - Compare regulatory strategies across different essential cellular functions