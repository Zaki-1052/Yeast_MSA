#!/usr/bin/env python3
# /Users/zakiralibhai/Documents/GitHub/Yeast_MSA/scripts/regulatory_analysis/regulatory_region_mapping_v2.py

"""
Maps variants in regulatory regions to precise genomic features (promoters, UTRs, etc.)
Enhanced version with yeast-specific regulatory logic.
Part of the Comprehensive Regulatory Region Analysis in the Yeast MSA project.
"""

import os
import sys
import argparse
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

class RegulatoryRegionMapper:
    """Maps variants to precise regulatory regions and analyzes their distribution"""
    
    def __init__(self, output_dir, gene_mapping_file, variants_file, erg_genes_file):
        """Initialize with file paths and configurations"""
        # Setup directories
        self.output_dir = self._ensure_dir(output_dir)
        self.data_dir = self._ensure_dir(os.path.join(output_dir, 'data'))
        self.plot_dir = self._ensure_dir(os.path.join(output_dir, 'plots'))
        
        # Load required data
        print("Loading input files...")
        try:
            self.gene_mapping = pd.read_csv(gene_mapping_file, sep='\t')
            print(f"Loaded gene mapping file with {len(self.gene_mapping)} entries")
        except Exception as e:
            print(f"Error loading gene mapping file: {e}")
            sys.exit(1)
            
        try:
            self.variants = pd.read_csv(variants_file, sep='\t')
            print(f"Loaded variants file with {len(self.variants)} entries")
        except Exception as e:
            print(f"Error loading variants file: {e}")
            sys.exit(1)
            
        try:
            self.erg_genes = pd.read_csv(erg_genes_file, sep='\t')
            print(f"Loaded ERG genes file with {len(self.erg_genes)} entries")
        except Exception as e:
            print(f"Error loading ERG genes file: {e}")
            sys.exit(1)
        
        # Define yeast-specific regulatory regions (distances in bp)
        # Note: In yeast, gene-proximal regions are particularly important
        self.define_regulatory_regions()
        
        # Treatment groups
        self.treatments = ['WT-37', 'WTA', 'CAS', 'STC']
        self.treatment_groups = {
            'Temperature': ['WT-37', 'CAS'],
            'Low Oxygen': ['WTA', 'STC'],
            'Gene Modified': ['CAS', 'STC'],
            'Non-Modified': ['WT-37', 'WTA']
        }
        
        # Color scheme for consistent visualization
        self.colors = {
            'WT-37': '#1f77b4',   # Blue
            'WTA': '#ff7f0e',     # Orange
            'CAS': '#2ca02c',     # Green
            'STC': '#d62728',     # Red
            'Temperature': '#8c564b',  # Brown
            'Low Oxygen': '#9467bd',   # Purple
            'Gene Modified': '#e377c2', # Pink
            'Non-Modified': '#7f7f7f'   # Gray
        }
        
    def _ensure_dir(self, directory):
        """Create directory if it doesn't exist"""
        if not os.path.exists(directory):
            try:
                os.makedirs(directory)
                print(f"Created directory: {directory}")
            except Exception as e:
                print(f"Error creating directory {directory}: {e}")
                sys.exit(1)
        return directory
        
    def define_regulatory_regions(self):
        """Define yeast-specific regulatory region categories"""
        # These definitions are based on yeast literature
        self.regulatory_regions = {
            # Core promoter regions - most critical for gene expression
            'core_promoter': (0, 150),            # 0-150bp upstream
            'TATA_box_region': (40, 120),         # TATA box region (typically 40-120bp upstream)
            'proximal_UAS': (150, 500),           # Upstream Activating Sequence (proximal)
            'distal_UAS': (500, 1500),            # Distal regulatory elements
            
            # Extended regulatory regions
            'proximal_regulatory': (1500, 3000),  # Other proximal regulatory elements
            'distal_regulatory': (3000, 10000),   # Distal regulatory region
            
            # UTR regions
            'five_prime_UTR': (-60, 0),           # 5' UTR (usually ~20-60bp in yeast)
            'three_prime_UTR': (0, 120),          # 3' UTR (typically ~50-100bp in yeast)
            
            # Downstream elements
            'terminator': (0, 250),               # Terminator region
            'downstream_element': (250, 1000)     # Other downstream regulatory elements
        }
        
        # Define ERG gene proximity zones based on the four-zone architecture
        self.erg_proximity_zones = {
            'core_zone': (0, 0),                # The gene itself
            'buffer_zone': (0, 5000),           # Buffer zone (<5kb)
            'intermediate_zone': (5000, 50000), # Intermediate zone (5-50kb)
            'satellite_zone': (50000, 100000)   # Satellite zone (50-100kb)
        }
        
        # Define genetic context classifications
        self.genetic_contexts = {
            'divergent_promoter': 'Shared promoter region between divergent genes',
            'convergent_terminator': 'Shared terminator region between convergent genes',
            'tandem_intergenic': 'Intergenic region between tandem genes',
            'isolated_promoter': 'Promoter not shared with other genes',
            'isolated_terminator': 'Terminator not shared with other genes'
        }
    
    def prepare_data(self):
        """Prepare and preprocess data for analysis"""
        print("\nPreparing data for regulatory region mapping...")
        
        # Check variant effect annotation - see what types of variants we have
        if 'Effect' in self.variants.columns:
            effects = self.variants['Effect'].value_counts()
            print("\nVariant effects distribution:")
            for effect, count in effects.items():
                print(f"  {effect}: {count} variants")
                
        # Examine Distance values
        if 'Distance' in self.variants.columns:
            # Check range and distribution of distance values
            dist_min = self.variants['Distance'].min()
            dist_max = self.variants['Distance'].max()
            dist_median = self.variants['Distance'].median()
            dist_mean = self.variants['Distance'].mean()
            
            print(f"\nDistance statistics:")
            print(f"  Range: {dist_min} to {dist_max} bp")
            print(f"  Median: {dist_median} bp")
            print(f"  Mean: {dist_mean} bp")
            
            # Count positive and negative distances
            pos_dist = (self.variants['Distance'] >= 0).sum()
            neg_dist = (self.variants['Distance'] < 0).sum()
            print(f"  Positive distances: {pos_dist} variants")
            print(f"  Negative distances: {neg_dist} variants")
            
        # Show a few samples to understand the data structure
        print("\nSample variants for examination:")
        sample_vars = self.variants.head(5)
        pd.set_option('display.max_columns', None)
        print(sample_vars[['Gene_ID', 'Effect', 'Distance', 'Scaffold', 'Position']].to_string())
        
        # Filter for regulatory variants - these have potentially regulatory effects
        print("\nFiltering for regulatory variants...")
        regulatory_effects = [
            'upstream_gene_variant',
            'downstream_gene_variant',
            'intergenic_region',
            '5_prime_UTR_variant',
            '3_prime_UTR_variant'
        ]
        
        # Check if we have an Effect column - if not, classify based on Distance only
        if 'Effect' in self.variants.columns:
            # Filter on effect types
            reg_variants = self.variants[
                self.variants['Effect'].str.contains('|'.join(regulatory_effects), case=False, na=False)
            ].copy()
            print(f"Found {len(reg_variants)} variants with regulatory effects")
        else:
            # No Effect column - use all variants
            reg_variants = self.variants.copy()
            print(f"No Effect column found. Using all {len(reg_variants)} variants.")
            
        # Make sure distances are numeric
        if 'Distance' in reg_variants.columns:
            reg_variants['Distance'] = pd.to_numeric(reg_variants['Distance'], errors='coerce')
            
        # Track what data we're using
        self.reg_variants = reg_variants
        return reg_variants
    
    def classify_variants(self):
        """Classify variants by regulatory region and genomic context"""
        print("\nClassifying variants by regulatory region and genomic context...")
        
        # Create new columns for classifications
        self.reg_variants['Regulatory_Region'] = 'unknown'
        self.reg_variants['ERG_Zone'] = 'unknown'
        self.reg_variants['Genetic_Context'] = 'unknown'
        
        # Check if we have the required data for classification
        if 'Distance' not in self.reg_variants.columns:
            print("WARNING: Distance column not found. Cannot perform distance-based classification.")
            return self.reg_variants
            
        # Classify into regulatory regions based on distance
        for region, (min_dist, max_dist) in self.regulatory_regions.items():
            # Get absolute distance for region classification
            # In our variants file, Distance represents distance from gene (positive = upstream)
            self.reg_variants.loc[
                (self.reg_variants['Distance'] >= min_dist) & 
                (self.reg_variants['Distance'] < max_dist),
                'Regulatory_Region'
            ] = region
        
        # Classify by ERG zone if applicable
        # First, identify variants associated with ERG genes
        erg_gene_ids = set(self.erg_genes['w303_gene_id'].values)
        
        # Check for ERG gene association
        for zone, (min_dist, max_dist) in self.erg_proximity_zones.items():
            # Mark variants that are in specific zones around ERG genes
            self.reg_variants.loc[
                (self.reg_variants['Gene_ID'].isin(erg_gene_ids)) & 
                (self.reg_variants['Distance'] >= min_dist) & 
                (self.reg_variants['Distance'] < max_dist),
                'ERG_Zone'
            ] = zone
                
        # Calculate statistics
        reg_region_counts = self.reg_variants['Regulatory_Region'].value_counts()
        erg_zone_counts = self.reg_variants[
            self.reg_variants['ERG_Zone'] != 'unknown'
        ]['ERG_Zone'].value_counts()
        
        # Print statistics
        print("\nRegulatory Region Counts:")
        for region, count in reg_region_counts.items():
            print(f"  {region}: {count} variants ({count/len(self.reg_variants)*100:.1f}%)")
            
        if not erg_zone_counts.empty:
            print("\nERG Zone Counts:")
            for zone, count in erg_zone_counts.items():
                print(f"  {zone}: {count} variants")
        
        # Return the classified dataframe
        return self.reg_variants
    
    def analyze_erg_gene_regulation(self):
        """Analyze regulatory variants around ergosterol pathway genes"""
        print("\nAnalyzing regulatory variants around ergosterol pathway genes...")
        
        # Get ergosterol gene IDs
        erg_gene_ids = set(self.erg_genes['w303_gene_id'].values)
        
        # Filter for variants associated with ERG genes
        erg_variants = self.reg_variants[
            self.reg_variants['Gene_ID'].isin(erg_gene_ids)
        ].copy()
        
        print(f"Found {len(erg_variants)} variants associated with ergosterol pathway genes")
        
        if len(erg_variants) == 0:
            print("No ERG gene variants found. Skipping ERG gene analysis.")
            return
            
        # Analyze distribution by regulatory region and treatment
        if 'Treatment' in erg_variants.columns:
            reg_by_treatment = pd.crosstab(
                erg_variants['Treatment'],
                erg_variants['Regulatory_Region'],
                normalize='index'
            ) * 100  # Convert to percentages
            
            # Save results
            erg_variants.to_csv(os.path.join(self.data_dir, 'erg_regulatory_variants.tsv'), sep='\t', index=False)
            reg_by_treatment.to_csv(os.path.join(self.data_dir, 'erg_region_by_treatment.tsv'), sep='\t')
            
            # Plot distribution
            plt.figure(figsize=(12, 8))
            reg_by_treatment.plot(kind='bar', stacked=True, colormap='viridis')
            plt.title('Distribution of Variants in Regulatory Regions of ERG Genes', fontsize=14)
            plt.xlabel('Treatment', fontsize=12)
            plt.ylabel('Percentage of Variants', fontsize=12)
            plt.xticks(rotation=45)
            plt.legend(title='Regulatory Region', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.savefig(os.path.join(self.plot_dir, 'erg_region_by_treatment.png'), dpi=300)
            plt.close()
        
        # Analyze distance distribution in more detail
        if 'Distance' in erg_variants.columns:
            # Create distance bins - more granular for close distances
            distance_bins = [
                0, 50, 100, 150, 200, 300, 500, 1000, 
                1500, 2000, 3000, 5000, 10000, float('inf')
            ]
            bin_labels = [
                '0-50', '50-100', '100-150', '150-200', '200-300', '300-500', 
                '500-1000', '1000-1500', '1500-2000', '2000-3000', '3000-5000',
                '5000-10000', '>10000'
            ]
            
            erg_variants['Distance_Bin'] = pd.cut(
                erg_variants['Distance'],
                bins=distance_bins,
                labels=bin_labels,
                right=False
            )
            
            # Analyze distribution by distance bin and treatment
            if 'Treatment' in erg_variants.columns:
                dist_by_treatment = pd.crosstab(
                    erg_variants['Treatment'],
                    erg_variants['Distance_Bin'],
                    normalize='index'
                ) * 100  # Convert to percentages
                
                # Save results
                dist_by_treatment.to_csv(os.path.join(self.data_dir, 'erg_distance_by_treatment.tsv'), sep='\t')
                
                # Plot distance distribution
                plt.figure(figsize=(14, 8))
                dist_by_treatment.plot(kind='bar', stacked=True, colormap='plasma')
                plt.title('Distance Distribution of Variants from ERG Genes', fontsize=14)
                plt.xlabel('Treatment', fontsize=12)
                plt.ylabel('Percentage of Variants', fontsize=12)
                plt.xticks(rotation=45)
                plt.legend(title='Distance from Gene (bp)', bbox_to_anchor=(1.05, 1), loc='upper left')
                plt.tight_layout()
                plt.savefig(os.path.join(self.plot_dir, 'erg_distance_distribution.png'), dpi=300)
                plt.close()
                
                # Create a heatmap for distance vs treatment
                plt.figure(figsize=(14, 8))
                sns.heatmap(dist_by_treatment, annot=True, fmt='.1f', cmap='viridis',
                          linewidths=0.5, cbar_kws={'label': 'Percentage of Variants'})
                plt.title('Heatmap of Variant Distance from ERG Genes by Treatment', fontsize=14)
                plt.xlabel('Distance from Gene (bp)', fontsize=12)
                plt.ylabel('Treatment', fontsize=12)
                plt.tight_layout()
                plt.savefig(os.path.join(self.plot_dir, 'erg_distance_heatmap.png'), dpi=300)
                plt.close()
    
    def analyze_regulatory_distances(self):
        """Analyze the distribution of variants by distance from genes"""
        print("\nAnalyzing regulatory variant distances from genes...")
        
        if 'Distance' not in self.reg_variants.columns:
            print("WARNING: Distance column not found. Skipping distance analysis.")
            return
        
        # Create binned distances with more granularity near genes
        distance_bins = [
            0, 50, 100, 150, 200, 300, 500, 1000, 
            1500, 2000, 3000, 5000, 10000, float('inf')
        ]
        bin_labels = [
            '0-50', '50-100', '100-150', '150-200', '200-300', '300-500', 
            '500-1000', '1000-1500', '1500-2000', '2000-3000', '3000-5000',
            '5000-10000', '>10000'
        ]
        
        # Create a new column with binned distances
        self.reg_variants['Distance_Bin'] = pd.cut(
            self.reg_variants['Distance'],
            bins=distance_bins,
            labels=bin_labels,
            right=False
        )
        
        # Count variants by distance bin
        dist_counts = self.reg_variants['Distance_Bin'].value_counts().sort_index()
        
        # Calculate percentages
        dist_percent = (dist_counts / len(self.reg_variants) * 100).round(1)
        
        print("\nDistance distribution from genes:")
        for bin_label, percent in dist_percent.items():
            print(f"  {bin_label} bp: {percent}%")
        
        # Save results
        dist_counts_df = pd.DataFrame({
            'Distance_Bin': dist_counts.index,
            'Count': dist_counts.values,
            'Percent': dist_percent.values
        })
        dist_counts_df.to_csv(os.path.join(self.data_dir, 'distance_distribution.tsv'), sep='\t', index=False)
        
        # Plot distance distribution
        plt.figure(figsize=(12, 6))
        plt.bar(range(len(dist_counts)), dist_counts.values, tick_label=dist_counts.index)
        plt.title('Distribution of Variants by Distance from Genes', fontsize=14)
        plt.xlabel('Distance from Gene (bp)', fontsize=12)
        plt.ylabel('Number of Variants', fontsize=12)
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(self.plot_dir, 'distance_distribution.png'), dpi=300)
        plt.close()
        
        # Plot percentage distribution
        plt.figure(figsize=(12, 6))
        plt.bar(range(len(dist_percent)), dist_percent.values, tick_label=dist_percent.index)
        plt.title('Percentage Distribution of Variants by Distance from Genes', fontsize=14)
        plt.xlabel('Distance from Gene (bp)', fontsize=12)
        plt.ylabel('Percentage of Variants', fontsize=12)
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(self.plot_dir, 'distance_percent_distribution.png'), dpi=300)
        plt.close()
        
        # If we have treatment information, analyze distance by treatment
        if 'Treatment' in self.reg_variants.columns:
            dist_by_treatment = pd.crosstab(
                self.reg_variants['Treatment'],
                self.reg_variants['Distance_Bin'],
                normalize='index'
            ) * 100  # Convert to percentages
            
            # Save results
            dist_by_treatment.to_csv(os.path.join(self.data_dir, 'distance_by_treatment.tsv'), sep='\t')
            
            # Plot distance by treatment
            plt.figure(figsize=(14, 8))
            dist_by_treatment.plot(kind='bar', stacked=True, colormap='tab20c')
            plt.title('Distance Distribution by Treatment', fontsize=14)
            plt.xlabel('Treatment', fontsize=12)
            plt.ylabel('Percentage of Variants', fontsize=12)
            plt.xticks(rotation=45)
            plt.legend(title='Distance from Gene (bp)', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.savefig(os.path.join(self.plot_dir, 'distance_by_treatment.png'), dpi=300)
            plt.close()
    
    def analyze_treatment_differences(self):
        """Analyze differences in regulatory patterns between treatments"""
        print("\nAnalyzing treatment differences in regulatory patterns...")
        
        if 'Treatment' not in self.reg_variants.columns:
            print("WARNING: Treatment column not found. Skipping treatment analysis.")
            return
            
        if 'Regulatory_Region' not in self.reg_variants.columns:
            print("WARNING: Regulatory_Region column not found. Skipping regulatory region analysis.")
            return
        
        # Crosstab of Treatment vs Regulatory Region (normalized by treatment)
        region_by_treatment = pd.crosstab(
            self.reg_variants['Treatment'],
            self.reg_variants['Regulatory_Region'],
            normalize='index'
        ) * 100  # Convert to percentages
        
        # Crosstab of Treatment Group vs Regulatory Region
        if 'Treatment' in self.reg_variants.columns:
            # Map treatments to treatment groups
            treatment_to_group = {}
            for group, treatments in self.treatment_groups.items():
                for treatment in treatments:
                    treatment_to_group[treatment] = group
                    
            self.reg_variants['Treatment_Group'] = self.reg_variants['Treatment'].map(
                treatment_to_group
            )
            
            region_by_group = pd.crosstab(
                self.reg_variants['Treatment_Group'],
                self.reg_variants['Regulatory_Region'],
                normalize='index'
            ) * 100  # Convert to percentages
            
            # Save results
            region_by_treatment.to_csv(os.path.join(self.data_dir, 'region_by_treatment.tsv'), sep='\t')
            region_by_group.to_csv(os.path.join(self.data_dir, 'region_by_group.tsv'), sep='\t')
            
            # Plot region by treatment
            plt.figure(figsize=(14, 8))
            region_by_treatment.plot(kind='bar', stacked=True, colormap='viridis')
            plt.title('Regulatory Region Distribution by Treatment', fontsize=14)
            plt.xlabel('Treatment', fontsize=12)
            plt.ylabel('Percentage of Variants', fontsize=12)
            plt.xticks(rotation=45)
            plt.legend(title='Regulatory Region', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.savefig(os.path.join(self.plot_dir, 'region_by_treatment.png'), dpi=300)
            plt.close()
            
            # Plot region by treatment group
            plt.figure(figsize=(14, 8))
            region_by_group.plot(kind='bar', stacked=True, colormap='plasma')
            plt.title('Regulatory Region Distribution by Treatment Group', fontsize=14)
            plt.xlabel('Treatment Group', fontsize=12)
            plt.ylabel('Percentage of Variants', fontsize=12)
            plt.xticks(rotation=45)
            plt.legend(title='Regulatory Region', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.savefig(os.path.join(self.plot_dir, 'region_by_group.png'), dpi=300)
            plt.close()
            
            # Create heatmaps for clearer visualization
            plt.figure(figsize=(14, 8))
            sns.heatmap(region_by_treatment, annot=True, fmt='.1f', cmap='viridis',
                       linewidths=0.5, cbar_kws={'label': 'Percentage of Variants'})
            plt.title('Heatmap of Regulatory Region Distribution by Treatment', fontsize=14)
            plt.xlabel('Regulatory Region', fontsize=12)
            plt.ylabel('Treatment', fontsize=12)
            plt.tight_layout()
            plt.savefig(os.path.join(self.plot_dir, 'region_treatment_heatmap.png'), dpi=300)
            plt.close()
            
            plt.figure(figsize=(12, 8))
            sns.heatmap(region_by_group, annot=True, fmt='.1f', cmap='plasma',
                       linewidths=0.5, cbar_kws={'label': 'Percentage of Variants'})
            plt.title('Heatmap of Regulatory Region Distribution by Treatment Group', fontsize=14)
            plt.xlabel('Regulatory Region', fontsize=12)
            plt.ylabel('Treatment Group', fontsize=12)
            plt.tight_layout()
            plt.savefig(os.path.join(self.plot_dir, 'region_group_heatmap.png'), dpi=300)
            plt.close()
    
    def generate_summary_report(self):
        """Generate a comprehensive summary report of the regulatory analysis"""
        print("\nGenerating summary report...")
        
        # Gather summary statistics
        total_variants = len(self.variants)
        reg_variants = len(self.reg_variants)
        reg_percent = (reg_variants / total_variants * 100) if total_variants > 0 else 0
        
        # Regulatory region distribution
        if 'Regulatory_Region' in self.reg_variants.columns:
            region_counts = self.reg_variants['Regulatory_Region'].value_counts()
            region_percent = (region_counts / reg_variants * 100).round(1)
        else:
            region_counts = pd.Series()
            region_percent = pd.Series()
            
        # Treatment distribution
        if 'Treatment' in self.reg_variants.columns:
            treatment_counts = self.reg_variants['Treatment'].value_counts()
            treatment_percent = (treatment_counts / reg_variants * 100).round(1)
        else:
            treatment_counts = pd.Series()
            treatment_percent = pd.Series()
            
        # ERG gene variants
        erg_gene_ids = set(self.erg_genes['w303_gene_id'].values)
        erg_variants = len(self.reg_variants[self.reg_variants['Gene_ID'].isin(erg_gene_ids)])
        erg_percent = (erg_variants / reg_variants * 100) if reg_variants > 0 else 0
        
        # Create report content
        report = ["# Regulatory Region Mapping Analysis Report\n"]
        
        # Overview section
        report.append("## Overview\n")
        report.append(f"Total variants analyzed: {total_variants}")
        report.append(f"Regulatory variants: {reg_variants} ({reg_percent:.1f}%)")
        report.append(f"Ergosterol gene regulatory variants: {erg_variants} ({erg_percent:.1f}%)\n")
        
        # Regulatory region distribution
        if not region_counts.empty:
            report.append("## Regulatory Region Distribution\n")
            for region, count in region_counts.items():
                percent = region_percent[region]
                report.append(f"{region}: {count} variants ({percent}%)")
            report.append("")
            
        # Treatment distribution
        if not treatment_counts.empty:
            report.append("## Treatment Distribution\n")
            for treatment, count in treatment_counts.items():
                percent = treatment_percent[treatment]
                report.append(f"{treatment}: {count} variants ({percent}%)")
            report.append("")
                
        # Key findings and interpretation
        report.append("## Key Findings\n")
        
        # Note any standout regulatory regions
        if not region_percent.empty and len(region_percent) > 0:
            top_region = region_percent.idxmax()
            top_percent = region_percent.max()
            report.append(f"- Most common regulatory region: {top_region} ({top_percent}% of variants)")
        
        # Note core promoter stats
        if 'core_promoter' in region_counts:
            core_percent = region_percent.get('core_promoter', 0)
            report.append(f"- Core promoter variants: {core_percent}% of regulatory variants")
        else:
            report.append("- No core promoter variants identified")
            
        # Note ERG gene regulatory patterns
        report.append(f"- Ergosterol gene regulatory variants: {erg_percent:.1f}% of all regulatory variants")
        
        # Treatment-specific patterns
        if 'Treatment' in self.reg_variants.columns and 'Regulatory_Region' in self.reg_variants.columns:
            report.append("\n## Treatment-Specific Patterns\n")
            
            # Get region distribution for each treatment
            for treatment in self.treatments:
                if treatment in self.reg_variants['Treatment'].values:
                    treatment_data = self.reg_variants[self.reg_variants['Treatment'] == treatment]
                    if len(treatment_data) > 0:
                        report.append(f"### {treatment}")
                        
                        # Get top regulatory regions for this treatment
                        t_region_counts = treatment_data['Regulatory_Region'].value_counts()
                        t_region_percent = (t_region_counts / len(treatment_data) * 100).round(1)
                        
                        for region, percent in t_region_percent.nlargest(3).items():
                            report.append(f"- {region}: {percent}%")
                        
                        report.append("")
        
        # Write report to file
        report_path = os.path.join(self.output_dir, 'regulatory_mapping_report.txt')
        with open(report_path, 'w') as f:
            f.write('\n'.join(report))
            
        print(f"Summary report saved to {report_path}")
        
        # Save all data for further analysis
        all_data_path = os.path.join(self.data_dir, 'all_regulatory_variants.tsv')
        self.reg_variants.to_csv(all_data_path, sep='\t', index=False)
        print(f"All regulatory variant data saved to {all_data_path}")
        
        # Return summary stats
        return {
            'total_variants': total_variants,
            'regulatory_variants': reg_variants,
            'erg_variants': erg_variants,
            'region_distribution': region_counts.to_dict() if not region_counts.empty else {},
            'treatment_distribution': treatment_counts.to_dict() if not treatment_counts.empty else {}
        }
        
    def run_analysis(self):
        """Run the complete regulatory region mapping analysis pipeline"""
        # Step 1: Prepare and preprocess data
        self.prepare_data()
        
        # Step 2: Classify variants
        self.classify_variants()
        
        # Step 3: Analyze ERG gene regulation
        self.analyze_erg_gene_regulation()
        
        # Step 4: Analyze regulatory distances
        self.analyze_regulatory_distances()
        
        # Step 5: Analyze treatment differences
        self.analyze_treatment_differences()
        
        # Step 6: Generate summary report
        summary = self.generate_summary_report()
        
        print("\nRegulatory region mapping analysis complete!")
        return summary

def main():
    parser = argparse.ArgumentParser(description='Map and analyze regulatory regions around genes')
    parser.add_argument('--output-dir', default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/regulatory_analysis/regulatory_mapping',
                      help='Output directory for results')
    parser.add_argument('--gene-mapping', default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/reference/gene_mapping_full.tsv',
                      help='Gene mapping file with coordinates')
    parser.add_argument('--variants', default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/gene_variants_expanded/all_gene_variants.tsv',
                      help='Variants TSV file')
    parser.add_argument('--erg-genes', default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/reference/genes_of_interest_mapping.tsv',
                      help='ERG genes mapping file')
    
    args = parser.parse_args()
    
    print("\n=== Regulatory Region Mapping Analysis ===\n")
    
    # Initialize and run analysis
    mapper = RegulatoryRegionMapper(
        args.output_dir,
        args.gene_mapping,
        args.variants,
        args.erg_genes
    )
    
    # Run complete analysis
    mapper.run_analysis()

if __name__ == "__main__":
    main()