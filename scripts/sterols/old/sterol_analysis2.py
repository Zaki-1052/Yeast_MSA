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
    base_dir = 'results/sterol_analysis2'
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
    exploration_results['summary'].to_csv('results/sterol_analysis2/sterol_summary_stats.csv', index=False)
    exploration_results['sterol_dist'].to_csv('results/sterol_analysis2/sterol_distribution.csv')
    
    # Differential analysis
    print("Performing differential analysis...")
    stat_results = perform_statistical_tests(sterol_df)
    fold_changes = calculate_fold_changes(sterol_df)
    visualize_differential_results(stat_results, fold_changes)
    
    # Save differential analysis results
    stat_results.to_csv('results/sterol_analysis2/differential/statistical_tests.csv', index=False)
    fold_changes.to_csv('results/sterol_analysis2/differential/fold_changes.csv', index=False)
    
    # Pathway analysis
    print("Performing pathway analysis...")
    pathway_graph = create_pathway_graph()
    pathway_ratios = calculate_pathway_ratios(sterol_df, pathway_graph)
    visualize_pathway_flux(sterol_df, pathway_graph, pathway_ratios)
    
    # Save pathway analysis results
    pathway_ratios.to_csv('results/sterol_analysis2/pathway_analysis/sterol_ratios.csv', index=False)
    
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
    correlation_results.to_csv('results/sterol_analysis2/correlation/sterol_variant_correlation.csv', index=False)
    
    # Machine learning analysis
    print("Performing pattern analysis and machine learning...")
    pivot_df, metadata = prepare_data_for_ml(sterol_df)
    pca_df, explained_variance, loadings = perform_pca_analysis(pivot_df, metadata)
    Z, cluster_df = perform_clustering(pivot_df, metadata)
    visualize_ml_results(pca_df, explained_variance, loadings, Z, cluster_df)
    
    # Save ML results
    pca_df.to_csv('results/sterol_analysis2/pattern_analysis/pca_results.csv', index=False)
    loadings.to_csv('results/sterol_analysis2/pattern_analysis/pca_loadings.csv')
    cluster_df.to_csv('results/sterol_analysis2/pattern_analysis/cluster_results.csv', index=False)
    
    print("Sterol analysis complete. Results saved to 'results/sterol_analysis2/'")

if __name__ == "__main__":
    run_sterol_analysis()