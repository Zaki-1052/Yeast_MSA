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