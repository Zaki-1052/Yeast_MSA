#!/usr/bin/env python3
# /Users/zakiralibhai/Documents/GitHub/Yeast_MSA/scripts/regulatory_analysis/regulatory_region_mapping_fixed.py

"""
Maps variants in regulatory regions to precise genomic features
Part of the Comprehensive Regulatory Region Analysis in the Yeast MSA project.

This script properly analyzes the Distance values in our dataset to classify
variants into biologically meaningful yeast-specific regulatory regions.
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
        
        # Load required data with error handling
        print("Loading input files...")
        self.gene_mapping = self._load_file(gene_mapping_file, "gene mapping")
        self.variants = self._load_file(variants_file, "variants")
        self.erg_genes = self._load_file(erg_genes_file, "ERG genes")
        
        # Define yeast-specific regulatory regions (distances in bp)
        self._define_regulatory_regions()
        
        # Treatment groups
        self.treatments = ['WT-37', 'WTA', 'CAS', 'STC', 
                          'WT-CTRL', 'CAS-CTRL', 'STC-CTRL']
        self.treatment_groups = {
            'Temperature': ['WT-37', 'CAS'],
            'Low Oxygen': ['WTA', 'STC'],
            'Gene Modified': ['CAS', 'STC'],
            'Non-Modified': ['WT-37', 'WTA'],
            'Control': ['WT-CTRL', 'CAS-CTRL', 'STC-CTRL']
        }
        
        # Visualization colors
        self.colors = {
            'WT-37': '#1f77b4',    # Blue
            'WTA': '#ff7f0e',      # Orange
            'CAS': '#2ca02c',      # Green
            'STC': '#d62728',      # Red
            'WT-CTRL': '#7f7f7f',  # Gray
            'CAS-CTRL': '#7f7f7f', # Gray 
            'STC-CTRL': '#7f7f7f'  # Gray
        }
    
    def _ensure_dir(self, directory):
        """Create directory if it doesn't exist"""
        if not os.path.exists(directory):
            os.makedirs(directory)
            print(f"Created directory: {directory}")
        return directory
    
    def _load_file(self, file_path, description):
        """Load a data file with error handling"""
        try:
            if file_path.endswith(('.tsv', '.txt')):
                data = pd.read_csv(file_path, sep='\t')
            elif file_path.endswith('.csv'):
                data = pd.read_csv(file_path)
            else:
                data = pd.read_csv(file_path, sep='\t')  # Default to TSV
            
            print(f"Loaded {description} file with {len(data)} entries")
            return data
        except Exception as e:
            print(f"Error loading {description} file ({file_path}): {e}")
            if description in ["gene mapping", "variants"]:
                print("Critical file missing. Cannot continue.")
                sys.exit(1)
            return pd.DataFrame()  # Return empty dataframe for non-critical files
    
    def _define_regulatory_regions(self):
        """Define yeast-specific regulatory regions with correct boundaries"""
        # IMPORTANT: In our variants, distance appears to be a POSITIVE value
        # representing base pairs upstream or downstream of gene
        
        # Upstream regions (promoter and regulatory elements)
        self.upstream_regions = {
            'core_promoter': (0, 150),          # Core promoter region
            'TATA_box_region': (40, 120),       # TATA box typically found here
            'proximal_UAS': (150, 500),         # Upstream Activating Sequence - proximal
            'distal_UAS': (500, 1500),          # Upstream Activating Sequence - distal
            'far_upstream': (1500, 10000)       # Distant regulatory elements
        }
        
        # Downstream regions (terminator and post-transcriptional elements)
        self.downstream_regions = {
            'five_prime_UTR': (-60, 0),         # 5' UTR
            'terminator': (0, 250),             # Termination region 
            'three_prime_UTR': (0, 120),        # 3' UTR
            'downstream_element': (250, 1000)   # Other downstream elements
        }
        
        # Ergosterol pathway gene proximity zones
        self.erg_zones = {
            'erg_core': (0, 0),               # Gene itself
            'erg_buffer': (0, 5000),          # Buffer zone (<5kb)
            'erg_intermediate': (5000, 50000), # Intermediate zone (5-50kb)
            'erg_satellite': (50000, 100000)  # Satellite zone (50-100kb)
        }
    
    def examine_variants(self):
        """Examine variant data to understand its structure"""
        print("\nExamining variant data structure...")
        
        # Check column names
        print(f"Columns in variants file: {list(self.variants.columns)}")
        
        # Check effect types if available
        if 'Effect' in self.variants.columns:
            effect_counts = self.variants['Effect'].value_counts()
            print("\nEffect types:")
            for effect, count in effect_counts.items():
                print(f"  {effect}: {count} variants")
        
        # Examine Distance values
        if 'Distance' in self.variants.columns:
            distances = self.variants['Distance'].dropna()
            print("\nDistance statistics:")
            print(f"  Range: {distances.min()} to {distances.max()} bp")
            print(f"  Mean: {distances.mean():.1f} bp")
            print(f"  Median: {distances.median()} bp")
            
            # Sample some distance values
            print("\nSample distance values:")
            print(self.variants[['Gene_ID', 'Effect', 'Distance']].head(10).to_string())
            
            # Check for negative distances
            neg_count = (distances < 0).sum()
            if neg_count > 0:
                print(f"\nFound {neg_count} negative distance values")
            else:
                print("\nAll distance values are positive")
        
        # Check effect/distance relationship for clearer understanding
        if 'Effect' in self.variants.columns and 'Distance' in self.variants.columns:
            print("\nEffect-Distance relationship:")
            print(self.variants.groupby('Effect')['Distance'].describe()[['count', 'min', 'max', 'mean']])
        
        # Determine if distances are upstream or downstream based on effect annotation
        if 'Effect' in self.variants.columns and 'Distance' in self.variants.columns:
            upstream_mask = self.variants['Effect'].str.contains('upstream', case=False, na=False)
            downstream_mask = self.variants['Effect'].str.contains('downstream', case=False, na=False)
            
            upstream_distances = self.variants.loc[upstream_mask, 'Distance']
            downstream_distances = self.variants.loc[downstream_mask, 'Distance']
            
            if not upstream_distances.empty:
                print(f"\nUpstream variant distances: {upstream_distances.min()} to {upstream_distances.max()} bp")
            
            if not downstream_distances.empty:
                print(f"Downstream variant distances: {downstream_distances.min()} to {downstream_distances.max()} bp")
    
    def prepare_variants(self):
        """Prepare variants for regulatory region analysis"""
        print("\nPreparing variants for regulatory region mapping...")
        
        # Filter for regulatory variants based on Effect if available
        if 'Effect' in self.variants.columns:
            # Common regulatory effects
            regulatory_effects = [
                'upstream_gene_variant',
                'downstream_gene_variant',
                'intergenic_region',
                '5_prime_UTR_variant',
                '3_prime_UTR_variant'
            ]
            
            # Filter for regulatory variants
            regulatory_mask = self.variants['Effect'].str.contains(
                '|'.join(regulatory_effects), case=False, na=False
            )
            reg_variants = self.variants[regulatory_mask].copy()
            print(f"Filtered {len(reg_variants)} regulatory variants from {len(self.variants)} total variants")
        else:
            # No Effect column - use all variants
            reg_variants = self.variants.copy()
            print(f"No Effect column found - using all {len(reg_variants)} variants")
        
        # Ensure Distance is numeric
        if 'Distance' in reg_variants.columns:
            reg_variants['Distance'] = pd.to_numeric(reg_variants['Distance'], errors='coerce')
        
        # Add classification columns
        reg_variants['Regulatory_Region'] = "unclassified"
        reg_variants['Region_Type'] = "unknown"
        reg_variants['ERG_Zone'] = "non_erg"
        
        # Store prepared variants
        self.reg_variants = reg_variants
        return reg_variants
    
    def classify_regulatory_regions(self):
        """Classify variants into specific yeast regulatory regions"""
        print("\nClassifying variants into regulatory regions...")
        
        if self.reg_variants is None or len(self.reg_variants) == 0:
            print("No variants to classify")
            return None
        
        # Ensure needed columns exist
        required_cols = ['Effect', 'Distance']
        missing_cols = [col for col in required_cols if col not in self.reg_variants.columns]
        if missing_cols:
            print(f"Warning: Missing required columns: {missing_cols}")
            return self.reg_variants
        
        # Process each variant
        for idx, row in self.reg_variants.iterrows():
            effect = row['Effect']
            distance = row['Distance']
            
            # Skip if missing data
            if pd.isna(distance) or pd.isna(effect):
                continue
            
            # Based on our examination of the data, all distances are positive
            # Classify based on the effect type and distance value
            
            # Upstream variants
            if 'upstream' in effect.lower():
                self.reg_variants.at[idx, 'Region_Type'] = 'upstream'
                
                # Classify into specific upstream regions
                for region, (min_dist, max_dist) in self.upstream_regions.items():
                    if min_dist <= distance < max_dist:
                        self.reg_variants.at[idx, 'Regulatory_Region'] = region
                        break
            
            # Downstream variants
            elif 'downstream' in effect.lower():
                self.reg_variants.at[idx, 'Region_Type'] = 'downstream'
                
                # Classify into specific downstream regions
                for region, (min_dist, max_dist) in self.downstream_regions.items():
                    if min_dist <= distance < max_dist:
                        self.reg_variants.at[idx, 'Regulatory_Region'] = region
                        break
            
            # UTR variants
            elif '5_prime_utr' in effect.lower():
                self.reg_variants.at[idx, 'Region_Type'] = 'utr'
                self.reg_variants.at[idx, 'Regulatory_Region'] = 'five_prime_UTR'
            
            elif '3_prime_utr' in effect.lower():
                self.reg_variants.at[idx, 'Region_Type'] = 'utr'
                self.reg_variants.at[idx, 'Regulatory_Region'] = 'three_prime_UTR'
            
            # Intergenic variants
            elif 'intergenic' in effect.lower():
                self.reg_variants.at[idx, 'Region_Type'] = 'intergenic'
                self.reg_variants.at[idx, 'Regulatory_Region'] = 'intergenic'
            
            # Classify ERG zone if near an ergosterol pathway gene
            if 'Gene_ID' in row and row['Gene_ID'] in self.erg_genes['w303_gene_id'].values:
                # It's an ERG gene - classify by distance
                for zone, (min_dist, max_dist) in self.erg_zones.items():
                    if min_dist <= distance < max_dist:
                        self.reg_variants.at[idx, 'ERG_Zone'] = zone
                        break
        
        # Calculate region statistics
        region_counts = self.reg_variants['Regulatory_Region'].value_counts()
        region_percent = (region_counts / len(self.reg_variants) * 100).round(1)
        
        # Display classification results
        print("\nRegulatory region classification:")
        for region, count in region_counts.items():
            percent = region_percent[region]
            print(f"  {region}: {count} variants ({percent}%)")
        
        # Calculate ERG zone statistics if any exist
        erg_variants = self.reg_variants[self.reg_variants['ERG_Zone'] != 'non_erg']
        if len(erg_variants) > 0:
            erg_zone_counts = erg_variants['ERG_Zone'].value_counts()
            print("\nERG zone classification:")
            for zone, count in erg_zone_counts.items():
                print(f"  {zone}: {count} variants")
        
        return self.reg_variants
    
    def analyze_erg_gene_regions(self):
        """Analyze regulatory regions specifically around ergosterol pathway genes"""
        print("\nAnalyzing regulatory regions around ergosterol pathway genes...")
        
        if self.reg_variants is None or len(self.reg_variants) == 0:
            print("No variants to analyze")
            return None
        
        # Get ergosterol gene IDs
        erg_gene_ids = set(self.erg_genes['w303_gene_id'].values)
        
        # Filter for variants associated with ERG genes
        erg_variants = self.reg_variants[self.reg_variants['Gene_ID'].isin(erg_gene_ids)].copy()
        
        if len(erg_variants) == 0:
            print("No variants found near ergosterol pathway genes")
            return None
        
        print(f"Found {len(erg_variants)} variants in regulatory regions of ergosterol genes")
        
        # Save ERG regulatory variants
        erg_variants.to_csv(os.path.join(self.data_dir, 'erg_regulatory_variants.tsv'), sep='\t', index=False)
        
        # Analyze regulatory region distribution
        region_counts = erg_variants['Regulatory_Region'].value_counts()
        region_percent = (region_counts / len(erg_variants) * 100).round(1)
        
        print("\nERG gene regulatory region distribution:")
        for region, count in region_counts.items():
            percent = region_percent[region]
            print(f"  {region}: {count} variants ({percent}%)")
        
        # Save region distribution
        region_distribution = pd.DataFrame({
            'Region': region_counts.index,
            'Count': region_counts.values,
            'Percent': region_percent.values
        })
        region_distribution.to_csv(os.path.join(self.data_dir, 'erg_region_distribution.tsv'), 
                                  sep='\t', index=False)
        
        # Plot region distribution
        plt.figure(figsize=(10, 6))
        plt.bar(region_counts.index, region_counts.values, color='darkblue')
        plt.title('Distribution of Regulatory Regions around ERG Genes', fontsize=14)
        plt.xlabel('Regulatory Region', fontsize=12)
        plt.ylabel('Number of Variants', fontsize=12)
        plt.xticks(rotation=45, ha='right')
        for i, v in enumerate(region_counts.values):
            plt.text(i, v + 0.5, str(v), ha='center')
        plt.tight_layout()
        plt.savefig(os.path.join(self.plot_dir, 'erg_region_distribution.png'), dpi=300)
        plt.close()
        
        # Analyze distribution by treatment if available
        if 'Treatment' in erg_variants.columns:
            treatment_region = pd.crosstab(
                erg_variants['Treatment'],
                erg_variants['Regulatory_Region'],
                normalize='index'
            ) * 100  # Convert to percentages
            
            # Save treatment-region distribution
            treatment_region.to_csv(os.path.join(self.data_dir, 'erg_treatment_region.tsv'), sep='\t')
            
            # Plot as heatmap
            plt.figure(figsize=(12, 8))
            sns.heatmap(treatment_region, annot=True, fmt='.1f', cmap='viridis',
                       linewidths=0.5, cbar_kws={'label': 'Percentage of Variants'})
            plt.title('Regulatory Region Distribution by Treatment (ERG Genes)', fontsize=14)
            plt.xlabel('Regulatory Region', fontsize=12)
            plt.ylabel('Treatment', fontsize=12)
            plt.tight_layout()
            plt.savefig(os.path.join(self.plot_dir, 'erg_treatment_region_heatmap.png'), dpi=300)
            plt.close()
            
            # Plot as stacked bar chart
            plt.figure(figsize=(12, 8))
            treatment_region.plot(kind='bar', stacked=True, colormap='viridis')
            plt.title('Distribution of Regulatory Regions by Treatment (ERG Genes)', fontsize=14)
            plt.xlabel('Treatment', fontsize=12)
            plt.ylabel('Percentage of Variants', fontsize=12)
            plt.xticks(rotation=45)
            plt.legend(title='Regulatory Region', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.savefig(os.path.join(self.plot_dir, 'erg_treatment_region_stacked.png'), dpi=300)
            plt.close()
        
        # Analyze distance distribution for ERG genes
        if 'Distance' in erg_variants.columns:
            # Create more precise bins for yeast promoter analysis
            distance_bins = [0, 50, 100, 150, 200, 300, 500, 1000, 1500, float('inf')]
            bin_labels = ['0-50', '50-100', '100-150', '150-200', '200-300', 
                         '300-500', '500-1000', '1000-1500', '>1500']
            
            erg_variants['Distance_Bin'] = pd.cut(
                erg_variants['Distance'], 
                bins=distance_bins,
                labels=bin_labels,
                right=False
            )
            
            # Count variants by distance bin
            distance_counts = erg_variants['Distance_Bin'].value_counts().sort_index()
            
            # Save distance distribution
            distance_distribution = pd.DataFrame({
                'Distance_Bin': distance_counts.index,
                'Count': distance_counts.values,
                'Percent': (distance_counts / len(erg_variants) * 100).round(1).values
            })
            distance_distribution.to_csv(os.path.join(self.data_dir, 'erg_distance_distribution.tsv'),
                                        sep='\t', index=False)
            
            # Plot distance distribution
            plt.figure(figsize=(10, 6))
            plt.bar(range(len(distance_counts)), distance_counts.values, tick_label=distance_counts.index)
            plt.title('Distribution of Variants by Distance from ERG Genes', fontsize=14)
            plt.xlabel('Distance from Gene (bp)', fontsize=12)
            plt.ylabel('Number of Variants', fontsize=12)
            plt.xticks(rotation=45)
            for i, v in enumerate(distance_counts.values):
                plt.text(i, v + 0.5, str(v), ha='center')
            plt.tight_layout()
            plt.savefig(os.path.join(self.plot_dir, 'erg_distance_distribution.png'), dpi=300)
            plt.close()
            
            # Analyze distance by treatment if available
            if 'Treatment' in erg_variants.columns:
                treatment_distance = pd.crosstab(
                    erg_variants['Treatment'],
                    erg_variants['Distance_Bin'],
                    normalize='index'
                ) * 100  # Convert to percentages
                
                # Save treatment-distance distribution
                treatment_distance.to_csv(os.path.join(self.data_dir, 'erg_treatment_distance.tsv'), sep='\t')
                
                # Plot as heatmap
                plt.figure(figsize=(12, 8))
                sns.heatmap(treatment_distance, annot=True, fmt='.1f', cmap='viridis',
                           linewidths=0.5, cbar_kws={'label': 'Percentage of Variants'})
                plt.title('Distance Distribution by Treatment (ERG Genes)', fontsize=14)
                plt.xlabel('Distance from Gene (bp)', fontsize=12)
                plt.ylabel('Treatment', fontsize=12)
                plt.tight_layout()
                plt.savefig(os.path.join(self.plot_dir, 'erg_treatment_distance_heatmap.png'), dpi=300)
                plt.close()
        
        return erg_variants
    
    def analyze_treatment_differences(self):
        """Analyze differences in regulatory patterns between treatments"""
        print("\nAnalyzing treatment differences in regulatory patterns...")
        
        if self.reg_variants is None or len(self.reg_variants) == 0:
            print("No variants to analyze")
            return None
        
        if 'Treatment' not in self.reg_variants.columns:
            print("Treatment information not available")
            return None
        
        # Count variants by treatment
        treatment_counts = self.reg_variants['Treatment'].value_counts()
        print("\nVariants by treatment:")
        for treatment, count in treatment_counts.items():
            print(f"  {treatment}: {count} variants")
        
        # Compare regulatory region distribution across treatments
        if 'Regulatory_Region' in self.reg_variants.columns:
            # Cross-tabulate regions by treatment
            treatment_region = pd.crosstab(
                self.reg_variants['Treatment'],
                self.reg_variants['Regulatory_Region'],
                normalize='index'
            ) * 100  # Convert to percentages
            
            # Save results
            treatment_region.to_csv(os.path.join(self.data_dir, 'treatment_region_distribution.tsv'), sep='\t')
            
            # Plot as heatmap
            plt.figure(figsize=(12, 8))
            sns.heatmap(treatment_region, annot=True, fmt='.1f', cmap='viridis',
                       linewidths=0.5, cbar_kws={'label': 'Percentage of Variants'})
            plt.title('Regulatory Region Distribution by Treatment', fontsize=14)
            plt.xlabel('Regulatory Region', fontsize=12)
            plt.ylabel('Treatment', fontsize=12)
            plt.tight_layout()
            plt.savefig(os.path.join(self.plot_dir, 'treatment_region_heatmap.png'), dpi=300)
            plt.close()
            
            # Plot as stacked bar chart
            plt.figure(figsize=(12, 8))
            treatment_region.plot(kind='bar', stacked=True, colormap='viridis')
            plt.title('Distribution of Regulatory Regions by Treatment', fontsize=14)
            plt.xlabel('Treatment', fontsize=12)
            plt.ylabel('Percentage of Variants', fontsize=12)
            plt.xticks(rotation=45)
            plt.legend(title='Regulatory Region', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.savefig(os.path.join(self.plot_dir, 'treatment_region_stacked.png'), dpi=300)
            plt.close()
        
        # Compare treatments grouped by adaptation type or gene modification
        treatment_to_group = {}
        for group, treatments in self.treatment_groups.items():
            for treatment in treatments:
                if treatment in treatment_to_group:
                    treatment_to_group[treatment] = f"{treatment_to_group[treatment]},{group}"
                else:
                    treatment_to_group[treatment] = group
        
        # Add group information to variants
        self.reg_variants['Treatment_Group'] = self.reg_variants['Treatment'].map(treatment_to_group)
        
        # Compare main adaptation types
        group_comparisons = [
            ('Temperature', 'Low Oxygen'),
            ('Gene Modified', 'Non-Modified')
        ]
        
        for group1, group2 in group_comparisons:
            # Define filtering function for group membership
            def in_group(x, group):
                if pd.isna(x):
                    return False
                return isinstance(x, str) and group in x.split(',')
            
            # Filter variants for these groups
            group1_variants = self.reg_variants[
                self.reg_variants['Treatment_Group'].apply(lambda x: in_group(x, group1))
            ]
            group2_variants = self.reg_variants[
                self.reg_variants['Treatment_Group'].apply(lambda x: in_group(x, group2))
            ]
            
            if len(group1_variants) == 0 or len(group2_variants) == 0:
                print(f"Insufficient data for {group1} vs {group2} comparison")
                continue
            
            print(f"\nComparing {group1} ({len(group1_variants)} variants) vs " + 
                  f"{group2} ({len(group2_variants)} variants)")
            
            # Compare regulatory region distributions
            if 'Regulatory_Region' in self.reg_variants.columns:
                group1_regions = group1_variants['Regulatory_Region'].value_counts(normalize=True) * 100
                group2_regions = group2_variants['Regulatory_Region'].value_counts(normalize=True) * 100
                
                # Combine into a dataframe
                region_comparison = pd.DataFrame({
                    group1: group1_regions,
                    group2: group2_regions
                }).fillna(0)
                
                # Save comparison
                region_comparison.to_csv(
                    os.path.join(self.data_dir, f'{group1}_vs_{group2}_regions.tsv'), sep='\t'
                )
                
                # Plot comparison
                plt.figure(figsize=(10, 6))
                region_comparison.plot(kind='bar')
                plt.title(f'Regulatory Region Comparison: {group1} vs {group2}', fontsize=14)
                plt.xlabel('Regulatory Region', fontsize=12)
                plt.ylabel('Percentage of Variants', fontsize=12)
                plt.xticks(rotation=45, ha='right')
                plt.legend(title='Treatment Group')
                plt.tight_layout()
                plt.savefig(os.path.join(self.plot_dir, f'{group1}_vs_{group2}_regions.png'), dpi=300)
                plt.close()
    
    def generate_summary_report(self):
        """Generate a comprehensive summary report of the analysis"""
        print("\nGenerating summary report...")
        
        if self.reg_variants is None or len(self.reg_variants) == 0:
            print("No variants to summarize")
            return None
        
        # Collect key statistics for report
        stats = {}
        
        # Basic counts
        stats['total_variants'] = len(self.variants)
        stats['regulatory_variants'] = len(self.reg_variants)
        stats['percent_regulatory'] = (stats['regulatory_variants'] / stats['total_variants'] * 100
                                     if stats['total_variants'] > 0 else 0)
        
        # Regulatory region distribution
        if 'Regulatory_Region' in self.reg_variants.columns:
            stats['region_counts'] = self.reg_variants['Regulatory_Region'].value_counts().to_dict()
            stats['region_percent'] = (self.reg_variants['Regulatory_Region'].value_counts(normalize=True) * 100
                                     ).round(1).to_dict()
        
        # Treatment distribution
        if 'Treatment' in self.reg_variants.columns:
            stats['treatment_counts'] = self.reg_variants['Treatment'].value_counts().to_dict()
            stats['treatment_percent'] = (self.reg_variants['Treatment'].value_counts(normalize=True) * 100
                                        ).round(1).to_dict()
        
        # ERG gene variants
        if 'Gene_ID' in self.reg_variants.columns:
            erg_gene_ids = set(self.erg_genes['w303_gene_id'].values)
            erg_variants = self.reg_variants[self.reg_variants['Gene_ID'].isin(erg_gene_ids)]
            stats['erg_variants'] = len(erg_variants)
            stats['percent_erg'] = (len(erg_variants) / len(self.reg_variants) * 100
                                  if len(self.reg_variants) > 0 else 0)
            
            # ERG region distribution
            if len(erg_variants) > 0 and 'Regulatory_Region' in erg_variants.columns:
                stats['erg_region_counts'] = erg_variants['Regulatory_Region'].value_counts().to_dict()
                stats['erg_region_percent'] = (erg_variants['Regulatory_Region'].value_counts(normalize=True) * 100
                                            ).round(1).to_dict()
        
        # Create report content
        report = ["# Regulatory Region Mapping Analysis Report\n"]
        
        # Overview section
        report.append("## Overview\n")
        report.append(f"Total variants analyzed: {stats['total_variants']}")
        report.append(f"Regulatory variants: {stats['regulatory_variants']} ({stats['percent_regulatory']:.1f}%)")
        if 'erg_variants' in stats:
            report.append(f"Ergosterol gene regulatory variants: {stats['erg_variants']} ({stats['percent_erg']:.1f}%)")
        report.append("")
        
        # Regulatory region distribution
        if 'region_counts' in stats:
            report.append("## Regulatory Region Distribution\n")
            for region, count in sorted(stats['region_counts'].items(), key=lambda x: x[1], reverse=True):
                percent = stats['region_percent'][region]
                report.append(f"{region}: {count} variants ({percent}%)")
            report.append("")
        
        # Treatment distribution
        if 'treatment_counts' in stats:
            report.append("## Treatment Distribution\n")
            for treatment, count in sorted(stats['treatment_counts'].items(), key=lambda x: x[1], reverse=True):
                percent = stats['treatment_percent'][treatment]
                report.append(f"{treatment}: {count} variants ({percent}%)")
            report.append("")
        
        # Key findings
        report.append("## Key Findings\n")
        
        # Add regulatory region observations
        if 'region_percent' in stats and stats['region_percent']:
            # Most common regulatory region
            top_region = max(stats['region_percent'].items(), key=lambda x: x[1])
            report.append(f"- Most common regulatory region: {top_region[0]} ({top_region[1]}% of variants)")
            
            # Core promoter statistics
            core_percent = stats['region_percent'].get('core_promoter', 0)
            report.append(f"- Core promoter variants: {core_percent}% of regulatory variants")
            
            # TATA box statistics
            tata_percent = stats['region_percent'].get('TATA_box_region', 0)
            if tata_percent > 0:
                report.append(f"- TATA box region variants: {tata_percent}% of regulatory variants")
            
            # UAS statistics
            proximal_uas = stats['region_percent'].get('proximal_UAS', 0)
            distal_uas = stats['region_percent'].get('distal_UAS', 0)
            uas_total = proximal_uas + distal_uas
            report.append(f"- UAS region variants: {uas_total}% of regulatory variants")
        
        # Add ERG gene observations
        if 'erg_variants' in stats and stats['erg_variants'] > 0:
            report.append(f"- Ergosterol gene regulatory variants: {stats['percent_erg']:.1f}% of all regulatory variants")
            
            if 'erg_region_percent' in stats and stats['erg_region_percent']:
                # Top ERG regulatory region
                top_erg_region = max(stats['erg_region_percent'].items(), key=lambda x: x[1])
                report.append(f"- Most common ERG gene regulatory region: {top_erg_region[0]} ({top_erg_region[1]}% of ERG variants)")
        
        # Treatment-specific patterns
        if 'Treatment' in self.reg_variants.columns and 'Regulatory_Region' in self.reg_variants.columns:
            report.append("\n## Treatment-Specific Patterns\n")
            
            # Get region distribution for each treatment
            for treatment in self.treatments:
                if treatment in self.reg_variants['Treatment'].values:
                    treatment_data = self.reg_variants[self.reg_variants['Treatment'] == treatment]
                    if len(treatment_data) > 0:
                        report.append(f"### {treatment}")
                        
                        # Top regulatory regions for this treatment
                        region_counts = treatment_data['Regulatory_Region'].value_counts()
                        if not region_counts.empty:
                            region_percent = (region_counts / len(treatment_data) * 100).round(1)
                            
                            for region, percent in region_percent.nlargest(3).items():
                                report.append(f"- {region}: {percent}%")
                        
                        report.append("")
        
        # Write report to file
        report_path = os.path.join(self.output_dir, 'regulatory_mapping_report.txt')
        with open(report_path, 'w') as f:
            f.write('\n'.join(report))
        
        print(f"Summary report saved to {report_path}")
        
        # Save all regulatory variants data
        variants_path = os.path.join(self.data_dir, 'all_regulatory_variants.tsv')
        self.reg_variants.to_csv(variants_path, sep='\t', index=False)
        print(f"All regulatory variants saved to {variants_path}")
        
        # Save statistics as JSON
        stats_path = os.path.join(self.data_dir, 'regulatory_mapping_statistics.json')
        with open(stats_path, 'w') as f:
            # Ensure all keys are strings for JSON serialization
            json_stats = {str(k): (v if not isinstance(v, dict) else {str(kk): vv for kk, vv in v.items()})
                         for k, v in stats.items()}
            json.dump(json_stats, f, indent=2)
        
        print(f"Statistics saved to {stats_path}")
        
        return stats
    
    def run_analysis(self):
        """Run the complete regulatory region mapping analysis pipeline"""
        print("\n=== Starting Regulatory Region Mapping Analysis ===\n")
        
        # Step 1: Examine variant structure
        self.examine_variants()
        
        # Step 2: Prepare variants for analysis
        self.prepare_variants()
        
        # Step 3: Classify into regulatory regions
        self.classify_regulatory_regions()
        
        # Step 4: Analyze ERG gene regulatory regions
        self.analyze_erg_gene_regions()
        
        # Step 5: Analyze treatment differences
        self.analyze_treatment_differences()
        
        # Step 6: Generate summary report
        self.generate_summary_report()
        
        print("\n=== Regulatory Region Mapping Analysis Complete ===")
        
        return self.reg_variants

def main():
    """Main function to run the analysis"""
    parser = argparse.ArgumentParser(description='Map and analyze regulatory regions around genes')
    
    # Required arguments
    parser.add_argument('--output-dir', 
                       default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/regulatory_analysis/regulatory_mapping_fixed',
                       help='Output directory for results')
    parser.add_argument('--gene-mapping', 
                       default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/reference/gene_mapping_full.tsv',
                       help='Gene mapping file with coordinates')
    parser.add_argument('--variants', 
                       default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/gene_variants_expanded/all_gene_variants.tsv',
                       help='Variants TSV file')
    parser.add_argument('--erg-genes', 
                       default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/reference/genes_of_interest_mapping.tsv',
                       help='ERG genes mapping file')
    
    args = parser.parse_args()
    
    # Initialize analyzer
    mapper = RegulatoryRegionMapper(
        args.output_dir,
        args.gene_mapping,
        args.variants,
        args.erg_genes
    )
    
    # Run analysis
    mapper.run_analysis()

if __name__ == "__main__":
    main()