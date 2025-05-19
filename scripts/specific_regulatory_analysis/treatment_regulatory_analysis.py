#!/usr/bin/env python3
# /Users/zakiralibhai/Documents/GitHub/Yeast_MSA/scripts/specific_regulatory_analysis/treatment_regulatory_analysis.py

"""
Treatment-Specific Regulatory Pattern Analysis for the Yeast MSA project.

This script identifies regulatory patterns specific to different adaptation conditions,
compares regulatory changes between temperature and low oxygen adaptation, and
analyzes the effect of gene modifications (CAS, STC) on regulatory patterns.

Part of the New Regulatory Analysis Planning implementation.
"""

import os
import sys
import json
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict, Counter
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


class TreatmentRegulatoryAnalyzer:
    """Analyzes treatment-specific regulatory patterns with statistical validation"""
    
    def __init__(self, mapped_variants_file, gene_regulatory_map=None, output_dir=None):
        """Initialize with necessary files and configurations"""
        # Setup directories
        self.output_dir = self._ensure_dir(output_dir or '/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/specific_regulatory_analysis')
        self.data_dir = self._ensure_dir(os.path.join(self.output_dir, 'data'))
        self.plot_dir = self._ensure_dir(os.path.join(self.output_dir, 'plots'))
        
        # Load variant data
        print("Loading mapped variant data...")
        try:
            self.variants = pd.read_csv(mapped_variants_file, sep='\t')
            print(f"Loaded {len(self.variants)} mapped variants")
        except Exception as e:
            print(f"Error loading mapped variants file: {e}")
            sys.exit(1)
        
        # Check for required columns
        required_cols = ['Regulatory_Region', 'Treatment', 'Scaffold', 'Position', 'Ref', 'Alt']
        missing_cols = [col for col in required_cols if col not in self.variants.columns]
        if missing_cols:
            print(f"Error: Missing required columns: {missing_cols}")
            print("These columns are necessary for treatment-specific pattern analysis.")
            sys.exit(1)
            
        # Identify treatment samples (excluding any control samples if present)
        treatment_samples = [t for t in self.variants['Treatment'].unique() if not t.endswith('-CTRL')]
        control_samples = [t for t in self.variants['Treatment'].unique() if t.endswith('-CTRL')]
        
        print(f"Found {len(treatment_samples)} treatment groups: {', '.join(treatment_samples)}")
        if control_samples:
            print(f"Found {len(control_samples)} control groups: {', '.join(control_samples)}")
            
        # Save a backup of full dataset
        self.full_variants = self.variants.copy()
        
        # Create a unique identifier for each variant
        self.variants['variant_key'] = self.variants.apply(
            lambda row: f"{row['Scaffold']}:{row['Position']}:{row['Ref']}:{row['Alt']}", axis=1
        )
        
        # For the fixed variant annotations, we're already using treatment-specific variants
        # Mark all variants as treatment_specific for consistency with the rest of the code
        self.variants['treatment_specific'] = True
        
        # Show variant count by treatment
        treatment_counts = self.variants['Treatment'].value_counts()
        print("\nVariant counts by treatment:")
        for treatment, count in treatment_counts.items():
            print(f"  {treatment}: {count} variants")
            
        # Create a column for treatment group
        if not 'Treatment_Group' in self.variants.columns:
            self.variants['Treatment_Group'] = self.variants['Treatment'].apply(
                lambda x: self._get_treatment_group(x)
            )
        
        # Load gene regulatory map if provided
        self.gene_map = None
        if gene_regulatory_map and os.path.exists(gene_regulatory_map):
            try:
                with open(gene_regulatory_map, 'r') as f:
                    self.gene_map = json.load(f)
                print(f"Loaded gene regulatory map with {len(self.gene_map)} genes")
            except Exception as e:
                print(f"Error loading gene regulatory map: {e}")
                print("Continuing without gene regulatory map information")
        
        # Define treatment groups
        self.treatment_groups = {
            'Temperature': ['WT-37', 'CAS'],
            'Low Oxygen': ['WTA', 'STC'],
            'Gene Modified': ['CAS', 'STC'],
            'Non-Modified': ['WT-37', 'WTA']
        }
        
        # Define control groups
        self.control_groups = {
            'WT-37': 'WT-CTRL',
            'WTA': 'WT-CTRL',
            'CAS': 'CAS-CTRL',
            'STC': 'STC-CTRL'
        }
        
        # Analysis results storage
        self.treatment_patterns = {}
        self.enrichment_results = {}
        self.treatment_comparison = {}
        self.effect_sizes = {}
    
    def _ensure_dir(self, directory):
        """Create directory if it doesn't exist"""
        if not os.path.exists(directory):
            os.makedirs(directory)
            print(f"Created directory: {directory}")
        return directory
    
    def analyze_treatment_distribution(self):
        """Analyze the distribution of regulatory regions by treatment"""
        print("\nAnalyzing distribution of regulatory regions by treatment...")
        
        # Get treatment information (excluding controls)
        treatments = sorted([t for t in self.variants['Treatment'].unique() if not t.endswith('-CTRL')])
        print(f"Found {len(treatments)} treatment groups: {', '.join(treatments)}")
        
        # Filter for variants with regulatory region information
        reg_variants = self.variants[self.variants['Regulatory_Region'] != 'unknown'].copy()
        
        if len(reg_variants) == 0:
            print("No variants with regulatory region information.")
            return None
            
        # Add message about treatment-specific variants
        if 'treatment_specific' in reg_variants.columns:
            specific_count = reg_variants['treatment_specific'].sum()
            print(f"Analyzing {len(reg_variants)} treatment-specific variants with regulatory region information")
        
        # Create a cross-tabulation of treatment vs. regulatory region
        # Raw counts
        count_table = pd.crosstab(
            reg_variants['Treatment'],
            reg_variants['Regulatory_Region']
        )
        
        # Normalize by treatment (percentage within each treatment)
        percent_table = pd.crosstab(
            reg_variants['Treatment'],
            reg_variants['Regulatory_Region'],
            normalize='index'
        ) * 100
        
        # Save tables
        count_file = os.path.join(self.data_dir, 'treatment_region_counts.tsv')
        count_table.to_csv(count_file, sep='\t')
        
        percent_file = os.path.join(self.data_dir, 'treatment_region_percentages.tsv')
        percent_table.to_csv(percent_file, sep='\t')
        
        print(f"Saved treatment-region distribution to {count_file} and {percent_file}")
        
        # Visualize the distribution
        plt.figure(figsize=(14, 10))
        sns.heatmap(percent_table, annot=True, fmt='.1f', cmap='viridis',
                   linewidths=0.5, cbar_kws={'label': 'Percentage of Variants'})
        plt.title('Regulatory Region Distribution by Treatment', fontsize=14)
        plt.xlabel('Regulatory Region', fontsize=12)
        plt.ylabel('Treatment', fontsize=12)
        plt.tight_layout()
        plt.savefig(os.path.join(self.plot_dir, 'treatment_region_heatmap.png'), dpi=300)
        plt.close()
        
        # Analyze TFBS impacts by treatment if available
        if 'TFBS' in reg_variants.columns:
            tfbs_variants = reg_variants[reg_variants['TFBS'].notna()].copy()
            
            if len(tfbs_variants) > 0:
                print(f"\nAnalyzing TFBS impacts for {len(tfbs_variants)} variants across treatments")
                
                # Create a dictionary to store TF counts by treatment
                treatment_tf_counts = {}
                
                for treatment in treatments:
                    treatment_variants = tfbs_variants[tfbs_variants['Treatment'] == treatment]
                    
                    if len(treatment_variants) == 0:
                        continue
                    
                    # Extract all TFs for this treatment
                    treatment_tfs = []
                    for tfs in treatment_variants['TFBS'].dropna():
                        if isinstance(tfs, str):
                            treatment_tfs.extend(tfs.split(','))
                    
                    # Store counts
                    treatment_tf_counts[treatment] = Counter(treatment_tfs)
                
                # Create a combined dataframe for visualization
                all_tfs = sorted(set(tf for counts in treatment_tf_counts.values() for tf in counts.keys()))
                tf_matrix = pd.DataFrame(0, index=treatments, columns=all_tfs)
                
                for treatment in treatments:
                    if treatment in treatment_tf_counts:
                        for tf, count in treatment_tf_counts[treatment].items():
                            tf_matrix.at[treatment, tf] = count
                
                # Normalize by treatment total variants
                for treatment in treatments:
                    treatment_total = len(reg_variants[reg_variants['Treatment'] == treatment])
                    if treatment_total > 0:
                        tf_matrix.loc[treatment] = tf_matrix.loc[treatment] / treatment_total * 100
                
                # Save TFBS impact matrix
                tf_file = os.path.join(self.data_dir, 'treatment_tfbs_impacts.tsv')
                tf_matrix.to_csv(tf_file, sep='\t')
                
                # Visualize TFBS impacts by treatment
                plt.figure(figsize=(14, 10))
                sns.heatmap(tf_matrix, annot=True, fmt='.1f', cmap='YlOrRd',
                           linewidths=0.5, cbar_kws={'label': 'Percentage of Treatment Variants'})
                plt.title('Transcription Factor Impacts by Treatment', fontsize=14)
                plt.xlabel('Transcription Factor', fontsize=12)
                plt.ylabel('Treatment', fontsize=12)
                plt.tight_layout()
                plt.savefig(os.path.join(self.plot_dir, 'treatment_tfbs_heatmap.png'), dpi=300)
                plt.close()
        
        # Store treatment distribution results
        self.treatment_patterns = {
            'count_table': count_table.to_dict(),
            'percent_table': percent_table.to_dict()
        }
        
        return {
            'count_table': count_table,
            'percent_table': percent_table
        }
    
    def analyze_treatment_groups(self):
        """Analyze regulatory patterns by treatment groups (Temperature, Low Oxygen, etc.)"""
        print("\nAnalyzing regulatory patterns by treatment groups...")
        
        # Get all treatment groups
        group_comparisons = [
            ('Temperature', 'Low Oxygen'),
            ('Gene Modified', 'Non-Modified')
        ]
        
        # Work with all variants, not just those with regulatory regions
        # We'll filter based on Effect column instead
        variants = self.variants.copy()
        
        # Add treatment group column
        variants['Treatment_Group'] = variants['Treatment'].apply(
            lambda x: self._get_treatment_groups(x)
        )
        
        # Ensure all treatments are assigned to a group
        unknown_treatments = variants[variants['Treatment_Group'] == 'Unknown']['Treatment'].unique()
        if len(unknown_treatments) > 0:
            print(f"Warning: Unknown treatment groups for treatments: {', '.join(unknown_treatments)}")
        
        # Create effect type distributions for each treatment group
        group_effect_counts = {}
        group_effect_percents = {}
        
        # Define effect categories
        effect_types = {
            'coding_missense': 'missense_variant',
            'coding_synonymous': 'synonymous_variant',
            'coding_frameshift': 'frameshift',  # Will use str.contains for this
            'distal_regulatory': 'upstream_gene_variant'
        }
        
        for group_name in set(group for groups in self.treatment_groups.values() for group in groups):
            group_variants = variants[variants['Treatment_Group'].str.contains(group_name)]
            
            if len(group_variants) > 0:
                # Count by effect type
                effect_counts = {}
                for effect_name, effect_pattern in effect_types.items():
                    if effect_name == 'coding_frameshift':
                        # Use string contains for frameshift
                        count = len(group_variants[group_variants['Effect'].str.contains(effect_pattern, na=False)])
                    else:
                        # Exact match for other types
                        count = len(group_variants[group_variants['Effect'] == effect_pattern])
                    effect_counts[effect_name] = count
                
                # Calculate percentages
                effect_percent = {k: (v / len(group_variants) * 100).round(1) for k, v in effect_counts.items()}
                
                # Store results
                group_effect_counts[group_name] = effect_counts
                group_effect_percents[group_name] = effect_percent
                
                print(f"\nEffect distribution for {group_name} (n={len(group_variants)}):")
                for effect, count in effect_counts.items():
                    if count > 0:
                        print(f"  {effect}: {count} variants ({effect_percent[effect]}%)")
        
        # Analyze pairwise group comparisons
        group_comparison_results = {}
        
        for group1, group2 in group_comparisons:
            print(f"\nComparing {group1} vs {group2}")
            
            # Get variants for each group
            group1_variants = variants[variants['Treatment_Group'].str.contains(group1)]
            group2_variants = variants[variants['Treatment_Group'].str.contains(group2)]
            
            if len(group1_variants) == 0 or len(group2_variants) == 0:
                print(f"  Insufficient data for comparison ({len(group1_variants)} vs {len(group2_variants)} variants)")
                continue
            
            print(f"  Comparing {len(group1_variants)} vs {len(group2_variants)} variants")
            
            # Compare effect distributions
            effect_comparison = {}
            for effect_name, effect_pattern in effect_types.items():
                # Calculate counts and percentages for each group
                if effect_name == 'coding_frameshift':
                    # Use string contains for frameshift
                    g1_count = len(group1_variants[group1_variants['Effect'].str.contains(effect_pattern, na=False)])
                    g2_count = len(group2_variants[group2_variants['Effect'].str.contains(effect_pattern, na=False)])
                else:
                    # Exact match for other types
                    g1_count = len(group1_variants[group1_variants['Effect'] == effect_pattern])
                    g2_count = len(group2_variants[group2_variants['Effect'] == effect_pattern])
                
                # Calculate percentages
                g1_percent = (g1_count / len(group1_variants) * 100) if len(group1_variants) > 0 else 0
                g2_percent = (g2_count / len(group2_variants) * 100) if len(group2_variants) > 0 else 0
                
                # Store result
                effect_comparison[effect_name] = {
                    'Group1_Count': g1_count,
                    'Group1_Percent': g1_percent,
                    'Group2_Count': g2_count,
                    'Group2_Percent': g2_percent,
                    'Difference': abs(g1_percent - g2_percent)
                }
            
            # Convert to dataframe for easier handling
            comparison_df = pd.DataFrame.from_dict(effect_comparison, orient='index')
            
            # Save comparison
            comparison_file = os.path.join(self.data_dir, f'{group1}_vs_{group2}_effects.tsv')
            comparison_df.to_csv(comparison_file, sep='\t')
            
            # Plot comparison
            plt.figure(figsize=(12, 6))
            plt.bar(comparison_df.index, comparison_df['Group1_Percent'], alpha=0.7, label=group1)
            plt.bar(comparison_df.index, comparison_df['Group2_Percent'], alpha=0.7, label=group2)
            plt.title(f'Effect Type Comparison: {group1} vs {group2}', fontsize=14)
            plt.xlabel('Effect Type', fontsize=12)
            plt.ylabel('Percentage of Variants', fontsize=12)
            plt.xticks(rotation=45, ha='right')
            plt.legend(title='Treatment Group')
            plt.tight_layout()
            plt.savefig(os.path.join(self.plot_dir, f'{group1}_vs_{group2}_effects.png'), dpi=300)
            plt.close()
            
            # Statistical testing
            stat_results = []
            
            # Loop through each effect type with sufficient data
            for effect_name, effect_data in effect_comparison.items():
                g1_count = effect_data['Group1_Count']
                g2_count = effect_data['Group2_Count']
                
                # Skip if either count is zero
                if g1_count == 0 or g2_count == 0:
                    continue
                
                # Get totals
                g1_total = len(group1_variants)
                g2_total = len(group2_variants)
                
                # Create contingency table
                table = np.array([
                    [g1_count, g1_total - g1_count],
                    [g2_count, g2_total - g2_count]
                ])
                
                # Skip if expected counts are too low
                if np.any(table < 5):
                    continue
                
                # Fisher's exact test
                odds_ratio, p_value = stats.fisher_exact(table)
                
                # Calculate effect size (Phi coefficient)
                chi2, p, dof, expected = stats.chi2_contingency(table, correction=False)
                n = table.sum()
                phi = np.sqrt(chi2 / n)
                
                # Store result
                stat_results.append({
                    'Effect': effect_name,
                    'Group1_Count': g1_count,
                    'Group1_Percent': effect_data['Group1_Percent'],
                    'Group2_Count': g2_count,
                    'Group2_Percent': effect_data['Group2_Percent'],
                    'Odds_Ratio': odds_ratio,
                    'P_Value': p_value,
                    'Effect_Size': phi
                })
            
            # Check if we have any results
            if not stat_results:
                print(f"  No statistical tests could be performed due to insufficient data")
                continue
            
            # Convert to dataframe and apply multiple testing correction
            stats_df = pd.DataFrame(stat_results)
            if len(stats_df) > 0:
                _, corrected_pvals, _, _ = multipletests(stats_df['P_Value'], method='fdr_bh')
                stats_df['Adjusted_P_Value'] = corrected_pvals
                stats_df['Significant'] = stats_df['Adjusted_P_Value'] < 0.05
            
                # Save statistical results
                stats_file = os.path.join(self.data_dir, f'{group1}_vs_{group2}_effect_statistics.tsv')
                stats_df.to_csv(stats_file, sep='\t', index=False)
                
                # Store comparison results
                group_comparison_results[f"{group1}_vs_{group2}"] = {
                    'comparison_df': comparison_df,
                    'statistics': stats_df
                }
                
                # Print significant differences
                sig_results = stats_df[stats_df['Significant']]
                if len(sig_results) > 0:
                    print(f"  Significant differences in effect types ({len(sig_results)}):")
                    for _, row in sig_results.iterrows():
                        direction = f"Higher in {group1}" if row['Group1_Percent'] > row['Group2_Percent'] else f"Higher in {group2}"
                        print(f"    {row['Effect']}: {abs(row['Group1_Percent'] - row['Group2_Percent']):.1f}% difference ({direction}, p={row['Adjusted_P_Value']:.4f})")
                else:
                    print("  No significant differences in effect types after multiple testing correction")
        
        # Store group analysis results
        self.treatment_comparison = {
            'group_effect_counts': group_effect_counts,
            'group_effect_percents': group_effect_percents,
            'group_comparisons': {
                k: {
                    'comparison': v['comparison_df'].to_dict(),
                    'statistics': v['statistics'].to_dict('records') if len(v['statistics']) > 0 else []
                }
                for k, v in group_comparison_results.items()
            }
        }
        
        return {
            'group_effect_counts': group_effect_counts,
            'group_effect_percents': group_effect_percents,
            'group_comparisons': group_comparison_results
        }
    
    def _get_treatment_groups(self, treatment):
        """Get all treatment groups that a treatment belongs to"""
        # Handle treatments with sample numbers (e.g., WT-37-55-1)
        base_treatment = treatment
        if '-55-' in treatment:
            # Extract the base treatment type (e.g., WT-37 from WT-37-55-1)
            base_treatment = treatment.split('-55-')[0]
        
        groups = []
        for group, treatments in self.treatment_groups.items():
            if base_treatment in treatments:
                groups.append(group)
        
        if not groups:
            # Handle control samples
            if treatment.endswith('-CTRL'):
                return "Control"
            return "Unknown"
        
        return ",".join(groups)
        
    def _get_treatment_group(self, treatment):
        """Get single primary treatment group for a treatment"""
        # Handle treatments with sample numbers
        base_treatment = treatment
        if '-55-' in treatment:
            # Extract the base treatment type
            base_treatment = treatment.split('-55-')[0]
            
        # Temperature group
        if base_treatment in ['WT-37', 'CAS']:
            return 'Temperature'
        # Low Oxygen group
        elif base_treatment in ['WTA', 'STC']:
            return 'Low Oxygen'
        # Control group
        elif treatment.endswith('-CTRL'):
            return 'Control'
        else:
            return 'Unknown'
    
    def analyze_enriched_regulatory_elements(self):
        """Analyze enrichment of specific regulatory elements in different treatments"""
        print("\nAnalyzing enrichment of specific regulatory elements in treatments...")
        
        # Filter for variants with regulatory information
        reg_variants = self.variants[self.variants['Regulatory_Region'] != 'unknown'].copy()
        
        if len(reg_variants) == 0:
            print("No variants with regulatory region information.")
            return None
        
        # Get treatments excluding controls
        treatments = [t for t in reg_variants['Treatment'].unique() if not t.endswith('-CTRL')]
        
        if len(treatments) == 0:
            print("No non-control treatments found.")
            return None
        
        print(f"Analyzing enrichment for treatments: {', '.join(treatments)}")
        
        # Analyze region enrichment for each treatment
        enrichment_results = {}
        
        for treatment in treatments:
            print(f"\nAnalyzing enrichment for {treatment}...")
            
            # Get treatment's control if available
            control = self.control_groups.get(treatment)
            
            if not control or control not in reg_variants['Treatment'].unique():
                print(f"  No control group found for {treatment}, skipping treatment vs control enrichment")
                continue
            
            # Get treatment and control variants
            treatment_vars = reg_variants[reg_variants['Treatment'] == treatment]
            control_vars = reg_variants[reg_variants['Treatment'] == control]
            
            if len(treatment_vars) == 0 or len(control_vars) == 0:
                print(f"  Insufficient data for {treatment} vs {control} enrichment analysis")
                continue
            
            print(f"  Comparing {len(treatment_vars)} {treatment} variants vs {len(control_vars)} {control} variants")
            
            # Count regions in treatment and control
            treatment_regions = treatment_vars['Regulatory_Region'].value_counts()
            control_regions = control_vars['Regulatory_Region'].value_counts()
            
            # Calculate enrichment for each region
            enrichment_data = []
            
            for region in set(list(treatment_regions.index) + list(control_regions.index)):
                # Get counts (default to 0 if not present)
                treatment_count = treatment_regions.get(region, 0)
                control_count = control_regions.get(region, 0)
                
                # Skip if both are 0
                if treatment_count == 0 and control_count == 0:
                    continue
                
                # Calculate frequencies
                treatment_freq = treatment_count / len(treatment_vars)
                control_freq = control_count / len(control_vars)
                
                # Calculate fold change (handle division by zero)
                if control_freq > 0:
                    fold_change = treatment_freq / control_freq
                else:
                    fold_change = float('inf') if treatment_freq > 0 else 1.0
                
                # Statistical testing (Fisher's exact test)
                table = np.array([
                    [treatment_count, len(treatment_vars) - treatment_count],
                    [control_count, len(control_vars) - control_count]
                ])
                
                odds_ratio, p_value = stats.fisher_exact(table)
                
                # Calculate effect size (Phi coefficient)
                if np.min(table) >= 5:  # Check if counts are sufficient
                    chi2, p, dof, expected = stats.chi2_contingency(table, correction=False)
                    n = table.sum()
                    phi = np.sqrt(chi2 / n)
                else:
                    phi = None
                
                # Store result
                enrichment_data.append({
                    'Treatment': treatment,
                    'Control': control,
                    'Region': region,
                    'Treatment_Count': treatment_count,
                    'Control_Count': control_count,
                    'Treatment_Frequency': treatment_freq,
                    'Control_Frequency': control_freq,
                    'Fold_Change': fold_change,
                    'P_Value': p_value,
                    'Odds_Ratio': odds_ratio,
                    'Effect_Size': phi
                })
            
            # Apply multiple testing correction if we have results
            if enrichment_data:
                enrichment_df = pd.DataFrame(enrichment_data)
                _, corrected_pvals, _, _ = multipletests(enrichment_df['P_Value'], method='fdr_bh')
                enrichment_df['Adjusted_P_Value'] = corrected_pvals
                enrichment_df['Significant'] = enrichment_df['Adjusted_P_Value'] < 0.05
                
                # Save enrichment results
                enrichment_file = os.path.join(self.data_dir, f'{treatment}_region_enrichment.tsv')
                enrichment_df.to_csv(enrichment_file, sep='\t', index=False)
                
                # Store results
                enrichment_results[treatment] = enrichment_df
                
                # Report significant enrichments
                sig_enrichment = enrichment_df[enrichment_df['Significant']]
                if len(sig_enrichment) > 0:
                    print(f"  Significantly enriched regions in {treatment} vs {control} ({len(sig_enrichment)}):")
                    for _, row in sig_enrichment.iterrows():
                        direction = "enriched" if row['Fold_Change'] > 1 else "depleted"
                        print(f"    {row['Region']}: {row['Fold_Change']:.2f}-fold {direction} (p={row['Adjusted_P_Value']:.4f})")
                else:
                    print(f"  No significantly enriched regions in {treatment} vs {control}")
            else:
                print(f"  No enrichment data available for {treatment}")
        
        # Analyze TFBS enrichment if data is available
        if 'TFBS' in reg_variants.columns:
            tfbs_variants = reg_variants[reg_variants['TFBS'].notna()]
            
            if len(tfbs_variants) > 0:
                print("\nAnalyzing transcription factor binding site enrichment...")
                
                tf_enrichment_results = {}
                
                for treatment in treatments:
                    # Get treatment's control if available
                    control = self.control_groups.get(treatment)
                    
                    if not control or control not in tfbs_variants['Treatment'].unique():
                        continue
                    
                    # Get treatment and control TFBS variants
                    treatment_tfbs = tfbs_variants[tfbs_variants['Treatment'] == treatment]
                    control_tfbs = tfbs_variants[tfbs_variants['Treatment'] == control]
                    
                    if len(treatment_tfbs) == 0 or len(control_tfbs) == 0:
                        continue
                    
                    print(f"  Analyzing TFBS enrichment for {treatment} vs {control}")
                    
                    # Extract all TFs for treatment and control
                    treatment_tfs = []
                    for tfs in treatment_tfbs['TFBS'].dropna():
                        if isinstance(tfs, str):
                            treatment_tfs.extend(tfs.split(','))
                    
                    control_tfs = []
                    for tfs in control_tfbs['TFBS'].dropna():
                        if isinstance(tfs, str):
                            control_tfs.extend(tfs.split(','))
                    
                    # Count TF occurrences
                    treatment_tf_counts = Counter(treatment_tfs)
                    control_tf_counts = Counter(control_tfs)
                    
                    # Calculate enrichment for each TF
                    tf_enrichment_data = []
                    
                    for tf in set(list(treatment_tf_counts.keys()) + list(control_tf_counts.keys())):
                        # Get counts
                        treatment_count = treatment_tf_counts.get(tf, 0)
                        control_count = control_tf_counts.get(tf, 0)
                        
                        # Skip if both are 0
                        if treatment_count == 0 and control_count == 0:
                            continue
                        
                        # Calculate frequencies
                        treatment_freq = treatment_count / len(treatment_tfbs)
                        control_freq = control_count / len(control_tfbs)
                        
                        # Calculate fold change (handle division by zero)
                        if control_freq > 0:
                            fold_change = treatment_freq / control_freq
                        else:
                            fold_change = float('inf') if treatment_freq > 0 else 1.0
                        
                        # Statistical testing
                        table = np.array([
                            [treatment_count, len(treatment_tfbs) - treatment_count],
                            [control_count, len(control_tfbs) - control_count]
                        ])
                        
                        odds_ratio, p_value = stats.fisher_exact(table)
                        
                        # Store result
                        tf_enrichment_data.append({
                            'Treatment': treatment,
                            'Control': control,
                            'Transcription_Factor': tf,
                            'Treatment_Count': treatment_count,
                            'Control_Count': control_count,
                            'Treatment_Frequency': treatment_freq,
                            'Control_Frequency': control_freq,
                            'Fold_Change': fold_change,
                            'P_Value': p_value,
                            'Odds_Ratio': odds_ratio
                        })
                    
                    # Apply multiple testing correction if we have results
                    if tf_enrichment_data:
                        tf_enrichment_df = pd.DataFrame(tf_enrichment_data)
                        _, corrected_pvals, _, _ = multipletests(tf_enrichment_df['P_Value'], method='fdr_bh')
                        tf_enrichment_df['Adjusted_P_Value'] = corrected_pvals
                        tf_enrichment_df['Significant'] = tf_enrichment_df['Adjusted_P_Value'] < 0.05
                        
                        # Save results
                        tf_enrichment_file = os.path.join(self.data_dir, f'{treatment}_tfbs_enrichment.tsv')
                        tf_enrichment_df.to_csv(tf_enrichment_file, sep='\t', index=False)
                        
                        # Store results
                        tf_enrichment_results[treatment] = tf_enrichment_df
                        
                        # Report significant enrichments
                        sig_tf_enrichment = tf_enrichment_df[tf_enrichment_df['Significant']]
                        if len(sig_tf_enrichment) > 0:
                            print(f"  Significantly enriched TFs in {treatment} vs {control} ({len(sig_tf_enrichment)}):")
                            for _, row in sig_tf_enrichment.iterrows():
                                direction = "enriched" if row['Fold_Change'] > 1 else "depleted"
                                print(f"    {row['Transcription_Factor']}: {row['Fold_Change']:.2f}-fold {direction} (p={row['Adjusted_P_Value']:.4f})")
                        else:
                            print(f"  No significantly enriched TFs in {treatment} vs {control}")
                
                # Store TF enrichment results
                if tf_enrichment_results:
                    self.enrichment_results['tf_enrichment'] = {
                        treatment: df.to_dict('records') 
                        for treatment, df in tf_enrichment_results.items()
                    }
        
        # Store enrichment results
        self.enrichment_results['region_enrichment'] = {
            treatment: df.to_dict('records') 
            for treatment, df in enrichment_results.items()
        }
        
        return enrichment_results
    
    def calculate_effect_sizes(self):
        """Calculate effect sizes for treatment-specific regulatory patterns"""
        print("\nCalculating effect sizes for treatment-specific regulatory patterns...")
        
        # Filter for variants with regulatory information
        reg_variants = self.variants[self.variants['Regulatory_Region'] != 'unknown'].copy()
        
        if len(reg_variants) == 0:
            print("No variants with regulatory region information.")
            return None
        
        # Define effect measures to calculate
        effect_measures = {
            'Odds_Ratio': "The odds of a variant being in a specific regulatory region in treatment vs control",
            'Phi_Coefficient': "Correlation coefficient between treatment and regulatory region",
            'Percent_Difference': "Absolute percentage point difference between treatment and control"
        }
        
        # Get treatments and their controls
        treatment_controls = {}
        for treatment in reg_variants['Treatment'].unique():
            if not treatment.endswith('-CTRL'):
                control = self.control_groups.get(treatment)
                if control and control in reg_variants['Treatment'].unique():
                    treatment_controls[treatment] = control
        
        if not treatment_controls:
            print("No valid treatment-control pairs found.")
            return None
        
        print(f"Calculating effect sizes for {len(treatment_controls)} treatment-control pairs")
        
        # Calculate effect sizes for each treatment-control pair
        effect_size_results = {}
        
        for treatment, control in treatment_controls.items():
            print(f"\nCalculating effect sizes for {treatment} vs {control}")
            
            # Get treatment and control variants
            treatment_vars = reg_variants[reg_variants['Treatment'] == treatment]
            control_vars = reg_variants[reg_variants['Treatment'] == control]
            
            if len(treatment_vars) == 0 or len(control_vars) == 0:
                print(f"  Insufficient data for {treatment} vs {control}")
                continue
            
            # Get all regions
            all_regions = sorted(set(
                list(treatment_vars['Regulatory_Region'].unique()) + 
                list(control_vars['Regulatory_Region'].unique())
            ))
            
            # Calculate effect sizes for each region
            region_effects = []
            
            for region in all_regions:
                if region == 'unknown':
                    continue
                
                # Get counts
                t_region_count = len(treatment_vars[treatment_vars['Regulatory_Region'] == region])
                c_region_count = len(control_vars[control_vars['Regulatory_Region'] == region])
                
                t_other_count = len(treatment_vars) - t_region_count
                c_other_count = len(control_vars) - c_region_count
                
                # Skip if insufficient data
                if t_region_count == 0 and c_region_count == 0:
                    continue
                
                # Calculate percentages
                t_percent = t_region_count / len(treatment_vars) * 100 if len(treatment_vars) > 0 else 0
                c_percent = c_region_count / len(control_vars) * 100 if len(control_vars) > 0 else 0
                
                # Effect size calculations
                effect_data = {
                    'Treatment': treatment,
                    'Control': control,
                    'Region': region,
                    'Treatment_Count': t_region_count,
                    'Control_Count': c_region_count,
                    'Treatment_Percent': t_percent,
                    'Control_Percent': c_percent,
                    'Percent_Difference': abs(t_percent - c_percent)
                }
                
                # Odds ratio
                if t_other_count > 0 and c_other_count > 0 and t_region_count > 0 and c_region_count > 0:
                    odds_ratio = (t_region_count / t_other_count) / (c_region_count / c_other_count)
                    effect_data['Odds_Ratio'] = odds_ratio
                else:
                    effect_data['Odds_Ratio'] = None
                
                # Phi coefficient
                contingency_table = np.array([
                    [t_region_count, t_other_count],
                    [c_region_count, c_other_count]
                ])
                
                if np.min(contingency_table) >= 5:  # Check if counts are sufficient
                    chi2, _, _, _ = stats.chi2_contingency(contingency_table)
                    n = contingency_table.sum()
                    phi = np.sqrt(chi2 / n)
                    effect_data['Phi_Coefficient'] = phi
                else:
                    effect_data['Phi_Coefficient'] = None
                
                # Add to results
                region_effects.append(effect_data)
            
            # Convert to dataframe
            if region_effects:
                effects_df = pd.DataFrame(region_effects)
                
                # Save effect sizes
                effects_file = os.path.join(self.data_dir, f'{treatment}_vs_{control}_effect_sizes.tsv')
                effects_df.to_csv(effects_file, sep='\t', index=False)
                
                # Store results
                effect_size_results[f"{treatment}_vs_{control}"] = effects_df
                
                # Report largest effect sizes
                print(f"  Largest effect sizes by percent difference:")
                for _, row in effects_df.nlargest(3, 'Percent_Difference').iterrows():
                    print(f"    {row['Region']}: {row['Percent_Difference']:.1f}% difference")
                
                if 'Phi_Coefficient' in effects_df.columns:
                    print(f"  Largest effect sizes by phi coefficient:")
                    phi_df = effects_df.dropna(subset=['Phi_Coefficient'])
                    if not phi_df.empty:
                        for _, row in phi_df.nlargest(3, 'Phi_Coefficient').iterrows():
                            print(f"    {row['Region']}: phi = {row['Phi_Coefficient']:.3f}")
            else:
                print(f"  No effect size data available for {treatment} vs {control}")
        
        # Store effect size results
        self.effect_sizes = {
            pair: df.to_dict('records')
            for pair, df in effect_size_results.items()
        }
        
        return effect_size_results
    
    def generate_treatment_patterns_report(self):
        """Generate a comprehensive report of treatment-specific regulatory patterns"""
        print("\nGenerating comprehensive treatment-specific regulatory patterns report...")
        
        # Create report content
        report = ["# Treatment-Specific Regulatory Pattern Analysis Report\n"]
        
        # 1. Overview section
        report.append("## Overview\n")
        
        treatments = sorted(self.variants['Treatment'].unique())
        report.append(f"Total variants analyzed: {len(self.variants)}")
        report.append(f"Number of treatments: {len(treatments)}")
        report.append(f"Treatments: {', '.join(treatments)}\n")
        
        # Extract data directly from the Effect column
        effect_counts = {
            'coding_missense': len(self.variants[self.variants['Effect'] == 'missense_variant']),
            'coding_synonymous': len(self.variants[self.variants['Effect'] == 'synonymous_variant']),
            'coding_frameshift': len(self.variants[self.variants['Effect'].str.contains('frameshift', na=False)]),
            'distal_regulatory': len(self.variants[self.variants['Effect'] == 'upstream_gene_variant'])
        }
        
        # Remove empty categories
        effect_counts = {k: v for k, v in effect_counts.items() if v > 0}
        
        # Calculate percentages
        effect_percents = {k: (v / len(self.variants) * 100) for k, v in effect_counts.items()}
        
        # Report effect distribution
        report.append("## Variant Effect Distribution\n")
        report.append("| Effect Type | Count | Percentage |")
        report.append("|-------------|-------|------------|")
        for effect, count in effect_counts.items():
            report.append(f"| {effect} | {count} | {effect_percents[effect]:.1f}% |")
        report.append("")
        
        # Calculate distance statistics if available
        if 'Distance' in self.variants.columns:
            variants_with_distance = self.variants.dropna(subset=['Distance'])
            if len(variants_with_distance) > 0:
                mean_distance = variants_with_distance['Distance'].abs().mean()
                median_distance = variants_with_distance['Distance'].abs().median()
                min_distance = variants_with_distance['Distance'].abs().min()
                max_distance = variants_with_distance['Distance'].abs().max()
                
                report.append("## Distance Relationships\n")
                report.append(f"Mean distance from gene: {mean_distance:.1f} bp")
                report.append(f"Median distance from gene: {median_distance:.1f} bp")
                report.append(f"Minimum distance from gene: {min_distance:.1f} bp")
                report.append(f"Maximum distance from gene: {max_distance:.1f} bp")
                report.append("")
        
        # 2. Treatment Distribution section
        if self.treatment_patterns:
            report.append("## Effect Distribution by Treatment\n")
            
            # Non-control treatments
            non_control_treatments = [t for t in treatments if not t.endswith('-CTRL')]
            
            # Create table header
            report.append("| Treatment | Coding Missense | Coding Synonymous | Coding Frameshift | Distal Regulatory |")
            report.append("|-----------|----------------|-------------------|-------------------|-------------------|")
            
            # Add row for each treatment
            for treatment in non_control_treatments:
                treatment_variants = self.variants[self.variants['Treatment'] == treatment]
                if len(treatment_variants) == 0:
                    continue
                
                # Calculate counts by effect
                t_missense = len(treatment_variants[treatment_variants['Effect'] == 'missense_variant'])
                t_synonymous = len(treatment_variants[treatment_variants['Effect'] == 'synonymous_variant'])
                t_frameshift = len(treatment_variants[treatment_variants['Effect'].str.contains('frameshift', na=False)])
                t_regulatory = len(treatment_variants[treatment_variants['Effect'] == 'upstream_gene_variant'])
                
                # Calculate percentages
                total = len(treatment_variants)
                if total > 0:
                    p_missense = (t_missense / total * 100)
                    p_synonymous = (t_synonymous / total * 100)
                    p_frameshift = (t_frameshift / total * 100)
                    p_regulatory = (t_regulatory / total * 100)
                    
                    # Add to table
                    report.append(f"| {treatment} | {p_missense:.1f}% | {p_synonymous:.1f}% | {p_frameshift:.1f}% | {p_regulatory:.1f}% |")
            
            report.append("")
        
        # 3. Treatment Group Comparison section
        if self.treatment_comparison and 'group_comparisons' in self.treatment_comparison:
            report.append("## Treatment Group Comparisons\n")
            
            for comparison, data in self.treatment_comparison['group_comparisons'].items():
                # Skip if no statistical results
                if 'statistics' not in data or not data['statistics']:
                    continue
                
                # Extract group names
                group1, group2 = comparison.split('_vs_')
                report.append(f"### {group1} vs {group2}\n")
                
                # Add effect comparison data
                if 'comparison' in data:
                    report.append("#### Effect Type Distribution\n")
                    report.append("| Effect Type | " + group1 + " | " + group2 + " | Difference |")
                    report.append("|------------|----------|----------|------------|")
                    
                    # Get the comparison data as a dataframe
                    comparison_data = pd.DataFrame.from_dict(data['comparison'], orient='index')
                    
                    # Sort by difference (largest first)
                    sorted_effects = comparison_data.sort_values('Difference', ascending=False).index
                    
                    for effect in sorted_effects:
                        row = comparison_data.loc[effect]
                        g1_pct = row['Group1_Percent']
                        g2_pct = row['Group2_Percent']
                        diff = abs(g1_pct - g2_pct)
                        direction = f"Higher in {group1}" if g1_pct > g2_pct else f"Higher in {group2}"
                        
                        report.append(f"| {effect} | {g1_pct:.1f}% | {g2_pct:.1f}% | {diff:.1f}% ({direction}) |")
                    
                    report.append("")
                
                # Add statistical results
                significant_results = [result for result in data['statistics'] if result.get('Significant', False)]
                
                if significant_results:
                    report.append(f"#### Significant Differences\n")
                    report.append(f"Found {len(significant_results)} significantly different effect types after correction for multiple testing:\n")
                    
                    report.append("| Effect Type | " + group1 + " | " + group2 + " | Difference | P-value | Effect Size |")
                    report.append("|------------|----------|----------|------------|---------|------------|")
                    
                    for result in sorted(significant_results, key=lambda x: x['Adjusted_P_Value']):
                        g1_pct = result['Group1_Percent']
                        g2_pct = result['Group2_Percent']
                        diff = abs(g1_pct - g2_pct)
                        direction = f"Higher in {group1}" if g1_pct > g2_pct else f"Higher in {group2}"
                        
                        report.append(f"| {result['Effect']} | {g1_pct:.1f}% | {g2_pct:.1f}% | {diff:.1f}% ({direction}) | {result['Adjusted_P_Value']:.4f} | {result['Effect_Size']:.3f} |")
                    
                    report.append("")
                else:
                    report.append("No significant differences found after correction for multiple testing.\n")
        
        # 4. Enrichment Analysis section
        if self.enrichment_results and 'region_enrichment' in self.enrichment_results:
            report.append("## Regulatory Region Enrichment Analysis\n")
            
            for treatment, results in self.enrichment_results['region_enrichment'].items():
                # Find significant enrichments
                significant_results = [result for result in results if result.get('Significant', False)]
                
                if significant_results:
                    control = significant_results[0]['Control']  # Get control name from first result
                    report.append(f"### {treatment} vs {control}\n")
                    
                    report.append(f"Found {len(significant_results)} significantly enriched regulatory regions:\n")
                    
                    report.append("| Region | Fold Change | Treatment | Control | P-value |")
                    report.append("|--------|-------------|-----------|---------|---------|")
                    
                    for result in sorted(significant_results, key=lambda x: x['Adjusted_P_Value']):
                        direction = "enriched" if result['Fold_Change'] > 1 else "depleted"
                        report.append(f"| {result['Region']} | {result['Fold_Change']:.2f}-fold {direction} | {result['Treatment_Count']} ({result['Treatment_Frequency']*100:.1f}%) | {result['Control_Count']} ({result['Control_Frequency']*100:.1f}%) | {result['Adjusted_P_Value']:.4f} |")
                    
                    report.append("")
            
            # TFBS enrichment results
            if 'tf_enrichment' in self.enrichment_results:
                report.append("## Transcription Factor Binding Site Enrichment\n")
                
                for treatment, results in self.enrichment_results['tf_enrichment'].items():
                    # Find significant enrichments
                    significant_results = [result for result in results if result.get('Significant', False)]
                    
                    if significant_results:
                        control = significant_results[0]['Control']  # Get control name from first result
                        report.append(f"### {treatment} vs {control}\n")
                        
                        report.append(f"Found {len(significant_results)} significantly enriched transcription factors:\n")
                        
                        report.append("| Transcription Factor | Fold Change | Treatment | Control | P-value |")
                        report.append("|----------------------|-------------|-----------|---------|---------|")
                        
                        for result in sorted(significant_results, key=lambda x: x['Adjusted_P_Value']):
                            direction = "enriched" if result['Fold_Change'] > 1 else "depleted"
                            report.append(f"| {result['Transcription_Factor']} | {result['Fold_Change']:.2f}-fold {direction} | {result['Treatment_Count']} ({result['Treatment_Frequency']*100:.1f}%) | {result['Control_Count']} ({result['Control_Frequency']*100:.1f}%) | {result['Adjusted_P_Value']:.4f} |")
                        
                        report.append("")
        
        # 5. Effect Size Analysis section
        if self.effect_sizes:
            report.append("## Effect Size Analysis\n")
            
            for comparison, results in self.effect_sizes.items():
                # Extract treatment and control names
                treatment, control = comparison.split('_vs_')
                report.append(f"### {treatment} vs {control}\n")
                
                # Sort by different effect measures
                if results and 'Percent_Difference' in results[0]:
                    # Sort by percent difference
                    sorted_by_diff = sorted(results, key=lambda x: x.get('Percent_Difference', 0), reverse=True)
                    
                    report.append("#### Top Regions by Percentage Difference\n")
                    
                    report.append("| Region | Treatment | Control | Difference |")
                    report.append("|--------|-----------|---------|------------|")
                    
                    for result in sorted_by_diff[:5]:  # Show top 5
                        t_pct = result['Treatment_Percent']
                        c_pct = result['Control_Percent']
                        diff = result['Percent_Difference']
                        direction = "higher in treatment" if t_pct > c_pct else "higher in control"
                        
                        report.append(f"| {result['Region']} | {t_pct:.1f}% | {c_pct:.1f}% | {diff:.1f}% ({direction}) |")
                    
                    report.append("")
                
                if results and 'Phi_Coefficient' in results[0]:
                    # Sort by phi coefficient
                    valid_phi = [r for r in results if r.get('Phi_Coefficient') is not None]
                    if valid_phi:
                        sorted_by_phi = sorted(valid_phi, key=lambda x: x.get('Phi_Coefficient', 0), reverse=True)
                        
                        report.append("#### Top Regions by Phi Coefficient\n")
                        
                        report.append("| Region | Phi Coefficient | Treatment | Control |")
                        report.append("|--------|----------------|-----------|---------|")
                        
                        for result in sorted_by_phi[:5]:  # Show top 5
                            report.append(f"| {result['Region']} | {result['Phi_Coefficient']:.3f} | {result['Treatment_Count']} ({result['Treatment_Percent']:.1f}%) | {result['Control_Count']} ({result['Control_Percent']:.1f}%) |")
                        
                        report.append("")
        
        # 6. Biological Interpretation section
        report.append("## Biological Interpretation\n")
        
        # Get treatment-specific variants (excluding controls)
        treatments = [t for t in self.variants['Treatment'].unique() if not t.endswith('-CTRL')]
        
        # Calculate overall statistics
        variants_total = len(self.variants[self.variants['Treatment'].isin(treatments)])
        coding_missense = len(self.variants[(self.variants['Treatment'].isin(treatments)) & (self.variants['Effect'] == 'missense_variant')])
        coding_synonymous = len(self.variants[(self.variants['Treatment'].isin(treatments)) & (self.variants['Effect'] == 'synonymous_variant')])
        coding_frameshift = len(self.variants[(self.variants['Treatment'].isin(treatments)) & (self.variants['Effect'].str.contains('frameshift', na=False))])
        distal_regulatory = len(self.variants[(self.variants['Treatment'].isin(treatments)) & (self.variants['Effect'] == 'upstream_gene_variant')])
        
        coding_variants = coding_missense + coding_synonymous + coding_frameshift
        
        # Calculate temperature vs low oxygen statistics
        temp_treatments = ['WT-37', 'CAS']
        oxy_treatments = ['WTA', 'STC']
        
        # Filter to just base treatments, handling sample replicate naming
        temp_variants = self.variants[self.variants['Treatment'].apply(
            lambda x: any(t in x for t in temp_treatments) and not x.endswith('-CTRL')
        )]
        oxy_variants = self.variants[self.variants['Treatment'].apply(
            lambda x: any(t in x for t in oxy_treatments) and not x.endswith('-CTRL')
        )]
        
        # Missense proportions
        if len(temp_variants) > 0:
            temp_missense = len(temp_variants[temp_variants['Effect'] == 'missense_variant']) / len(temp_variants) * 100
        else:
            temp_missense = 0
            
        if len(oxy_variants) > 0:
            oxy_missense = len(oxy_variants[oxy_variants['Effect'] == 'missense_variant']) / len(oxy_variants) * 100
        else:
            oxy_missense = 0
        
        # Treatment-specific analyses
        treatment_analyses = {}
        for treatment in ['WT-37', 'CAS', 'STC', 'WTA']:
            # Get variants for this treatment, including replicates
            t_variants = self.variants[self.variants['Treatment'].apply(
                lambda x: treatment in x and not x.endswith('-CTRL')
            )]
            
            if len(t_variants) == 0:
                continue
                
            # Calculate effect distributions
            t_missense = len(t_variants[t_variants['Effect'] == 'missense_variant']) / len(t_variants) * 100
            t_synonymous = len(t_variants[t_variants['Effect'] == 'synonymous_variant']) / len(t_variants) * 100
            t_frameshift = len(t_variants[t_variants['Effect'].str.contains('frameshift', na=False)]) / len(t_variants) * 100
            t_regulatory = len(t_variants[t_variants['Effect'] == 'upstream_gene_variant']) / len(t_variants) * 100
            
            treatment_analyses[treatment] = {
                'missense': t_missense,
                'synonymous': t_synonymous,
                'frameshift': t_frameshift,
                'regulatory': t_regulatory
            }
        
        # Calculate variant distance range
        if 'Distance' in self.variants.columns:
            min_distance = self.variants['Distance'].abs().min()
            max_distance = self.variants['Distance'].abs().max()
            mean_distance = self.variants['Distance'].abs().mean()
        
        # Add main interpretations
        # Effect type distribution
        if coding_variants > 0 and distal_regulatory > 0:
            report.append(f"- **Effect Type Distribution**: The variants show a mix of coding changes ({coding_variants} variants, {coding_variants/variants_total*100:.1f}%) " +
                         f"and regulatory changes ({distal_regulatory} variants, {distal_regulatory/variants_total*100:.1f}%), " +
                         "highlighting how adaptation involves both protein-level and gene expression mechanisms.")
        
        # Temperature vs low oxygen adaptation strategies
        if temp_missense > oxy_missense:
            report.append(f"- **Temperature Adaptation Strategy**: Temperature adaptation shows a higher proportion " +
                        f"of coding missense variants ({temp_missense:.1f}% vs {oxy_missense:.1f}%), suggesting protein-level adaptation plays a stronger role " +
                        "in response to thermal stress compared to oxygen limitation.")
        elif oxy_missense < temp_missense:
            report.append(f"- **Low Oxygen Adaptation Strategy**: Low oxygen adaptation shows a lower proportion " +
                        f"of coding missense variants ({oxy_missense:.1f}% vs {temp_missense:.1f}%), suggesting adaptation may occur primarily through " +
                        "regulatory changes rather than protein modifications.")
        
        # Treatment-specific patterns
        report.append("- **Treatment-Specific Variant Patterns**:")
        for treatment, analysis in treatment_analyses.items():
            report.append(f"  - {treatment}: {analysis['missense']:.1f}% missense, {analysis['synonymous']:.1f}% synonymous, " +
                        f"{analysis['frameshift']:.1f}% frameshift, {analysis['regulatory']:.1f}% regulatory")
        
        # Distance relationships
        if 'Distance' in self.variants.columns:
            report.append(f"- **Variant Distribution Range**: Variants occur at distances ranging from {min_distance:.1f}bp to " +
                        f"{max_distance:.1f}bp from genes, with a mean distance of {mean_distance:.1f}bp, demonstrating a broad " +
                        "spatial distribution of genetic changes that may affect gene function.")
        
        # Hierarchical architecture
        report.append("- **Hierarchical Regulatory Architecture**: The distribution of variants across different effect types " +
                     "and distances from genes reveals a hierarchical structure that balances essential function " +
                     "preservation with adaptive flexibility.")
        
        # Adaptation mechanism
        if coding_variants < distal_regulatory:
            report.append("- **Regulatory-Focused Adaptation**: The dominance of regulatory variants over coding changes " +
                         "suggests that adaptation typically occurs through alterations in gene expression rather than " +
                         "through modifications to protein structure and function.")
        
        # Add TFBS interpretation if available
        if 'tf_enrichment' in self.enrichment_results:
            all_tfs = set()
            for results in self.enrichment_results['tf_enrichment'].values():
                significant_tfs = [r['Transcription_Factor'] for r in results if r.get('Significant', False) and r['Fold_Change'] > 1]
                all_tfs.update(significant_tfs)
            
            if all_tfs:
                report.append("- **Transcription Factor Binding Site Alterations**: The analysis identified significant alterations in binding sites for " +
                             f"{', '.join(sorted(list(all_tfs)))} transcription factors, suggesting specific regulatory network rewiring during adaptation.")
        
        # Write report to file
        report_path = os.path.join(self.output_dir, 'treatment_regulatory_patterns_report.md')
        with open(report_path, 'w') as f:
            f.write('\n'.join(report))
        
        print(f"Generated comprehensive report: {report_path}")
        
        # Save combined results as JSON
        results_file = os.path.join(self.data_dir, 'treatment_regulatory_patterns.json')
        
        # Combine all results
        combined_results = {
            'treatment_patterns': self.treatment_patterns,
            'treatment_comparison': self.treatment_comparison,
            'enrichment_results': self.enrichment_results,
            'effect_sizes': self.effect_sizes
        }
        
        # Convert to JSON-serializable format
        json_results = {}
        for category, cat_data in combined_results.items():
            if isinstance(cat_data, dict):
                json_results[category] = {}
                for k, v in cat_data.items():
                    if isinstance(v, dict):
                        json_results[category][k] = {str(kk): vv for kk, vv in v.items()}
                    elif isinstance(v, pd.DataFrame):
                        json_results[category][k] = v.to_dict('records')
                    else:
                        json_results[category][k] = v
            else:
                json_results[category] = cat_data
        
        with open(results_file, 'w') as f:
            json.dump(json_results, f, indent=2)
        
        print(f"Saved combined results to {results_file}")
        
        return report
    
    def run_analysis(self):
        """Run the complete treatment-specific regulatory pattern analysis pipeline"""
        print("\n=== Starting Treatment-Specific Regulatory Pattern Analysis ===\n")
        
        # Step 1: Analyze treatment distribution of regulatory regions
        self.analyze_treatment_distribution()
        
        # Step 2: Analyze treatment group patterns
        self.analyze_treatment_groups()
        
        # Step 3: Analyze enriched regulatory elements
        self.analyze_enriched_regulatory_elements()
        
        # Step 4: Calculate effect sizes
        self.calculate_effect_sizes()
        
        # Step 5: Generate comprehensive report
        self.generate_treatment_patterns_report()
        
        print("\n=== Treatment-Specific Regulatory Pattern Analysis Complete ===")
        
        # Return results for potential use by other scripts
        return {
            'treatment_patterns': self.treatment_patterns,
            'treatment_comparison': self.treatment_comparison,
            'enrichment_results': self.enrichment_results,
            'effect_sizes': self.effect_sizes
        }


def main():
    """Main function to parse arguments and run the analysis"""
    parser = argparse.ArgumentParser(description='Treatment-Specific Regulatory Pattern Analysis')
    
    # Required arguments
    parser.add_argument('--mapped-variants', 
                       default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/specific_regulatory_analysis/data/fixed_variant_regulatory_annotations.tsv',
                       help='Mapped variants TSV file')
    parser.add_argument('--output-dir', 
                       default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/specific_regulatory_analysis',
                       help='Output directory for results')
    
    # Optional arguments
    parser.add_argument('--gene-regulatory-map', 
                       default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/new_regulatory_analysis/data/gene_regulatory_map.json',
                       help='Gene regulatory map JSON file')
    
    args = parser.parse_args()
    
    # Initialize the analyzer
    analyzer = TreatmentRegulatoryAnalyzer(
        args.mapped_variants,
        args.gene_regulatory_map,
        args.output_dir
    )
    
    # Run the analysis
    analyzer.run_analysis()


if __name__ == "__main__":
    main()