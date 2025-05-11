#!/usr/bin/env python3
# /Users/zakiralibhai/Documents/GitHub/Yeast_MSA/scripts/regulatory_analysis/regulatory_region_mapping.py

"""
Maps variants in regulatory regions to precise genomic features (promoters, UTRs, etc.)
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

# Add project root to path for imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils.tools import ensure_dir, load_tsv

class RegulatoryRegionMapper:
    """Maps variants to precise regulatory regions and analyzes their distribution"""
    
    def __init__(self, output_dir, gene_mapping_file, variants_file, erg_genes_file):
        """Initialize with file paths and configurations"""
        self.output_dir = ensure_dir(output_dir)
        self.data_dir = ensure_dir(os.path.join(output_dir, 'data'))
        self.plot_dir = ensure_dir(os.path.join(output_dir, 'plots'))

        # Load required data with error checking
        try:
            self.gene_mapping = pd.read_csv(gene_mapping_file, sep='\t')
            print(f"Loaded gene mapping file with {len(self.gene_mapping)} entries")
        except Exception as e:
            print(f"Error loading gene mapping file: {e}")
            self.gene_mapping = None

        try:
            self.variants = pd.read_csv(variants_file, sep='\t')
            print(f"Loaded variants file with {len(self.variants)} entries")
        except Exception as e:
            print(f"Error loading variants file: {e}")
            self.variants = None

        try:
            self.erg_genes = pd.read_csv(erg_genes_file, sep='\t')
            print(f"Loaded ERG genes file with {len(self.erg_genes)} entries")
        except Exception as e:
            print(f"Error loading ERG genes file: {e}")
            self.erg_genes = None
        
        # Define yeast-specific regulatory regions (distances in bp)
        self.regulatory_regions = {
            'core_promoter': (-150, 0),         # Core promoter (0-150bp upstream)
            'TATA_box_region': (-120, -40),     # TATA box typically 40-120bp upstream in yeast
            'UAS_proximal': (-500, -150),       # Proximal UAS (Upstream Activating Sequence)
            'UAS_distal': (-1500, -500),        # Distal UAS regions
            'proximal_regulatory': (-3000, -1500), # Other proximal regulatory elements
            'distal_regulatory': (-10000, -3000),  # Distal regulatory region
            '5_prime_UTR': (0, 60),             # 5' UTR (usually short in yeast, ~20-60bp)
            '3_prime_UTR': (-120, 0),           # 3' UTR (typically 50-100bp in yeast)
            'terminator': (0, 250),             # Terminator region
            'downstream_regulatory': (250, 1000) # Downstream regulatory elements
        }

        # Define genetic context classifications
        self.genetic_contexts = {
            'divergent_promoter': 'Shared promoter region between divergent genes',
            'convergent_terminator': 'Shared terminator region between convergent genes',
            'tandem_intergenic': 'Intergenic region between tandem genes',
            'isolated_promoter': 'Promoter not shared with other genes',
            'isolated_terminator': 'Terminator not shared with other genes',
            'ORF_proximal': 'Near ORF but not in typical regulatory region',
            'erg_proximal': 'Near ergosterol pathway gene (<5kb)',
            'erg_distal': 'Distal to ergosterol pathway gene (5-50kb)',
            'erg_satellite': 'In satellite zone of ergosterol pathway gene (50-100kb)'
        }
        
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
    
    def prepare_data(self):
        """Prepare and preprocess data for analysis"""
        print("Preparing data for regulatory region mapping...")

        # Check if data was loaded successfully
        if self.variants is None or self.gene_mapping is None or self.erg_genes is None:
            print("ERROR: Required data files couldn't be loaded. Aborting analysis.")
            self.reg_variants = pd.DataFrame()
            return self.reg_variants

        # Filter variants for regulatory regions (non-coding)
        try:
            self.reg_variants = self.variants[
                (self.variants['Effect'].str.contains('upstream|downstream|intergenic|UTR')) |
                ((self.variants['Effect'].str.contains('intron')) &
                 (~self.variants['Effect'].str.contains('coding')))
            ].copy()

            # Add gene distance information where missing
            if 'Distance' not in self.reg_variants.columns or self.reg_variants['Distance'].isna().any():
                print("Computing missing distances...")
                self.compute_distances()

            # Add regulatory region classification
            self.classify_regulatory_regions()

            # Add treatment group information
            self.reg_variants['Treatment_Group'] = self.reg_variants['Treatment'].map(
                {t: g for g, treatments in self.treatment_groups.items() for t in treatments}
            )

            print(f"Prepared {len(self.reg_variants)} regulatory variants for analysis")
        except Exception as e:
            print(f"Error preparing regulatory variant data: {e}")
            self.reg_variants = pd.DataFrame()

        return self.reg_variants
    
    def compute_distances(self):
        """Compute distance from variant to gene start/end where missing"""
        # Extract gene coordinates
        gene_coords = self.gene_mapping[['w303_gene_id', 'start', 'end', 'strand']].drop_duplicates()
        gene_coords = gene_coords.rename(columns={'w303_gene_id': 'Gene_ID'})
        
        # Merge with variants
        variants_with_coords = pd.merge(
            self.reg_variants, gene_coords, on='Gene_ID', how='left'
        )
        
        # Compute distances
        def calculate_distance(row):
            if pd.notna(row['Distance']):
                return row['Distance']
                
            pos = row['Position']
            
            if row['strand'] == '+':
                if 'upstream' in row['Effect']:
                    return row['start'] - pos  # Upstream of start (negative value)
                elif 'downstream' in row['Effect']:
                    return pos - row['end']    # Downstream of end (positive value)
            else:  # Negative strand
                if 'upstream' in row['Effect']:
                    return pos - row['end']    # Upstream of end (negative value)
                elif 'downstream' in row['Effect']:
                    return row['start'] - pos  # Downstream of start (positive value)
            
            # For intergenic, calculate nearest gene
            return min(abs(pos - row['start']), abs(pos - row['end']))
            
        variants_with_coords['Distance'] = variants_with_coords.apply(calculate_distance, axis=1)
        
        # Update the original dataframe
        self.reg_variants = variants_with_coords.drop(columns=['start', 'end', 'strand'])
    
    def classify_regulatory_regions(self):
        """Classify variants into specific regulatory regions"""
        # Create additional columns for rich classification
        self.reg_variants['Regulatory_Region'] = 'unknown'
        self.reg_variants['Genetic_Context'] = 'unknown'
        self.reg_variants['Motif_Context'] = 'unknown'

        # First, determine basic region type based on distance and effect
        print("Classifying variants into regulatory regions...")

        def classify_primary_region(row):
            dist = row['Distance']
            effect = row['Effect']

            # Just check if it's an int/float
            if not isinstance(dist, (int, float)):
                return 'unknown'

            # Basic region classification based on distance
            if 'upstream' in effect:
                # Classify upstream variants by distance
                if -150 <= dist <= 0:
                    return 'core_promoter'
                elif -500 <= dist < -150:
                    return 'UAS_proximal'
                elif -1500 <= dist < -500:
                    return 'UAS_distal'
                elif -3000 <= dist < -1500:
                    return 'proximal_regulatory'
                elif -10000 <= dist < -3000:
                    return 'distal_regulatory'
                else:
                    return 'far_upstream'

            elif 'downstream' in effect:
                # Classify downstream variants
                if 0 <= dist <= 250:
                    return 'terminator'
                elif 250 < dist <= 1000:
                    return 'downstream_regulatory'
                else:
                    return 'far_downstream'

            # UTR variants
            elif '5_prime_UTR' in effect:
                return '5_prime_UTR'
            elif '3_prime_UTR' in effect:
                return '3_prime_UTR'

            # Intron variants
            elif 'intron' in effect:
                return 'intron'

            # Intergenic variants - classify based on distance to nearest genes
            elif 'intergenic' in effect:
                # We'll do more detailed classification in the next step
                return 'intergenic'

            return 'other_regulatory'

        # Apply basic classification
        self.reg_variants['Regulatory_Region'] = self.reg_variants.apply(classify_primary_region, axis=1)

        # Identify potential TATA box regions
        tata_mask = (
            (self.reg_variants['Regulatory_Region'] == 'core_promoter') &
            (self.reg_variants['Distance'] >= -120) &
            (self.reg_variants['Distance'] <= -40)
        )
        self.reg_variants.loc[tata_mask, 'Regulatory_Region'] = 'TATA_box_region'

        # Now add contextual classification based on gene relationships
        def classify_genetic_context(row):
            gene_id = row['Gene_ID']
            dist = row['Distance']
            effect = row['Effect']

            # Check if this is near an ERG gene
            if gene_id in self.erg_genes['w303_gene_id'].values:
                if abs(dist) < 5000:
                    return 'erg_proximal'
                elif abs(dist) < 50000:
                    return 'erg_distal'
                elif abs(dist) < 100000:
                    return 'erg_satellite'

            # For intergenic regions, determine if it's between convergent/divergent genes
            if 'intergenic' in effect:
                # Extract the gene names from the effect field if possible
                if '-' in effect:
                    try:
                        genes = effect.split('|')[2].split('-')
                        gene1, gene2 = genes[0], genes[1]

                        # Check if these genes exist in our mapping
                        gene1_data = self.gene_mapping[self.gene_mapping['w303_gene_id'] == gene1]
                        gene2_data = self.gene_mapping[self.gene_mapping['w303_gene_id'] == gene2]

                        if not gene1_data.empty and not gene2_data.empty:
                            # Check the strand to determine convergent/divergent
                            gene1_strand = gene1_data['strand'].iloc[0]
                            gene2_strand = gene2_data['strand'].iloc[0]

                            if gene1_strand == '+' and gene2_strand == '+':
                                return 'tandem_intergenic'
                            elif gene1_strand == '-' and gene2_strand == '-':
                                return 'tandem_intergenic'
                            elif gene1_strand == '+' and gene2_strand == '-':
                                return 'convergent_terminator'
                            elif gene1_strand == '-' and gene2_strand == '+':
                                return 'divergent_promoter'
                    except:
                        # If we can't parse it, just use default
                        pass

            # Determine if this is a shared or isolated regulatory region
            if 'upstream' in effect:
                # Check if this region is shared with another gene
                near_genes = self.gene_mapping[
                    (self.gene_mapping['scaffold'] == row['Scaffold']) &
                    (self.gene_mapping['w303_gene_id'] != row['Gene_ID'])
                ]

                # Look for genes that might share this promoter region
                for _, gene in near_genes.iterrows():
                    # Check distance to this gene
                    pos = row['Position']
                    if gene['strand'] == '+':
                        other_dist = gene['start'] - pos
                        if 0 < other_dist < 1000:  # Within 1kb of another gene start
                            return 'shared_promoter'
                    else:
                        other_dist = pos - gene['end']
                        if 0 < other_dist < 1000:  # Within 1kb of another gene end
                            return 'shared_promoter'

                return 'isolated_promoter'

            elif 'downstream' in effect:
                return 'isolated_terminator'  # Default for downstream

            return 'ORF_proximal'  # Default for others

        # Apply genetic context classification where we have data
        if not self.reg_variants.empty and 'Gene_ID' in self.reg_variants.columns:
            self.reg_variants['Genetic_Context'] = self.reg_variants.apply(classify_genetic_context, axis=1)

        print(f"Classified {len(self.reg_variants)} variants into regulatory regions")

        # Generate some statistics
        region_counts = self.reg_variants['Regulatory_Region'].value_counts()
        print("\nRegulatory Region Distribution:")
        for region, count in region_counts.items():
            print(f"  {region}: {count} variants ({count/len(self.reg_variants)*100:.1f}%)")

        context_counts = self.reg_variants['Genetic_Context'].value_counts()
        print("\nGenetic Context Distribution:")
        for context, count in context_counts.items():
            print(f"  {context}: {count} variants ({count/len(self.reg_variants)*100:.1f}%)")

        # Return the enriched dataframe
        return self.reg_variants
    
    def analyze_erg_regulatory_regions(self):
        """Specifically analyze regulatory regions around ergosterol pathway genes"""
        print("Analyzing regulatory regions around ergosterol pathway genes...")

        # Get list of ERG genes
        erg_gene_ids = self.erg_genes['w303_gene_id'].unique()

        # Filter variants near ERG genes (direct annotation or by genetic context)
        erg_variants = self.reg_variants[
            (self.reg_variants['Gene_ID'].isin(erg_gene_ids)) |
            (self.reg_variants['Genetic_Context'].str.contains('erg_'))
        ].copy()

        print(f"Found {len(erg_variants)} variants in regulatory regions of ERG genes")

        # Analyze both regulatory region distribution and genetic context distribution
        erg_region_counts = erg_variants.groupby(['Treatment', 'Regulatory_Region']).size().unstack(fill_value=0)
        erg_context_counts = erg_variants.groupby(['Treatment', 'Genetic_Context']).size().unstack(fill_value=0)

        # Normalize to see proportions
        erg_region_percent = erg_region_counts.div(erg_region_counts.sum(axis=1), axis=0) * 100
        erg_context_percent = erg_context_counts.div(erg_context_counts.sum(axis=1), axis=0) * 100

        # Save results
        erg_variants.to_csv(os.path.join(self.data_dir, 'erg_regulatory_variants.tsv'), sep='\t', index=False)
        erg_region_counts.to_csv(os.path.join(self.data_dir, 'erg_regulatory_region_counts.tsv'), sep='\t')
        erg_region_percent.to_csv(os.path.join(self.data_dir, 'erg_regulatory_region_percent.tsv'), sep='\t')
        erg_context_counts.to_csv(os.path.join(self.data_dir, 'erg_genetic_context_counts.tsv'), sep='\t')
        erg_context_percent.to_csv(os.path.join(self.data_dir, 'erg_genetic_context_percent.tsv'), sep='\t')

        # Plot regulatory region distribution
        plt.figure(figsize=(14, 10))
        erg_region_percent.plot(kind='bar', stacked=True, colormap='viridis')
        plt.title('Distribution of Variants in Regulatory Regions of Ergosterol Genes', fontsize=14)
        plt.xlabel('Treatment', fontsize=12)
        plt.ylabel('Percentage of Variants', fontsize=12)
        plt.xticks(rotation=45)
        plt.legend(title='Regulatory Region', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(os.path.join(self.plot_dir, 'erg_regulatory_distribution.png'), dpi=300)
        plt.close()

        # Plot genetic context distribution
        plt.figure(figsize=(14, 10))
        erg_context_percent.plot(kind='bar', stacked=True, colormap='magma')
        plt.title('Genetic Context of Variants Near Ergosterol Genes', fontsize=14)
        plt.xlabel('Treatment', fontsize=12)
        plt.ylabel('Percentage of Variants', fontsize=12)
        plt.xticks(rotation=45)
        plt.legend(title='Genetic Context', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(os.path.join(self.plot_dir, 'erg_genetic_context_distribution.png'), dpi=300)
        plt.close()

        # Distance-specific analysis for ergosterol genes
        if len(erg_variants) > 0:
            # Create distance bins specific to ERG genes
            bins = [-10000, -5000, -3000, -1500, -500, -150, -50, 0, 50, 150, 500, 1500, 3000, 5000, 10000]
            erg_variants['Distance_Bin'] = pd.cut(erg_variants['Distance'], bins=bins)

            # Group by treatment and distance bin
            distance_counts = erg_variants.groupby(['Treatment', 'Distance_Bin']).size().unstack(fill_value=0)
            distance_percent = distance_counts.div(distance_counts.sum(axis=1), axis=0) * 100

            # Save distance analysis
            distance_counts.to_csv(os.path.join(self.data_dir, 'erg_distance_bin_counts.tsv'), sep='\t')
            distance_percent.to_csv(os.path.join(self.data_dir, 'erg_distance_bin_percent.tsv'), sep='\t')

            # Plot distance distribution
            plt.figure(figsize=(16, 10))
            ax = distance_percent.plot(kind='bar', stacked=False, colormap='tab20')
            plt.title('Distribution of Variants by Distance from Ergosterol Genes', fontsize=14)
            plt.xlabel('Treatment', fontsize=12)
            plt.ylabel('Percentage of Variants', fontsize=12)
            plt.xticks(rotation=45)

            # Add a vertical line at distance 0 (gene boundary)
            # Find the position of the bin containing 0
            zero_bin_index = None
            for i, bin_edge in enumerate(bins):
                if bin_edge == 0:
                    zero_bin_index = i - 1
                    break

            if zero_bin_index is not None:
                # Add a vertical line at the right edge of this bin
                ax.axvline(x=zero_bin_index + 0.5, color='red', linestyle='--',
                           label='Gene Boundary', linewidth=2)

            plt.legend(title='Distance Bin (bp)', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.savefig(os.path.join(self.plot_dir, 'erg_distance_distribution.png'), dpi=300)
            plt.close()

            # Create a focused heatmap of distance vs treatment
            pivot_data = erg_variants.pivot_table(
                index='Treatment',
                columns='Distance_Bin',
                aggfunc='size',
                fill_value=0
            )

            # Normalize to percentages
            pivot_norm = pivot_data.div(pivot_data.sum(axis=1), axis=0) * 100

            plt.figure(figsize=(14, 8))
            sns.heatmap(pivot_norm, annot=True, fmt='.1f', cmap='viridis',
                       linewidths=0.5, cbar_kws={'label': 'Percentage of Variants'})
            plt.title('Heatmap of Variant Distribution by Distance from Ergosterol Genes', fontsize=14)
            plt.xlabel('Distance Bin (bp)', fontsize=12)
            plt.ylabel('Treatment', fontsize=12)
            plt.tight_layout()
            plt.savefig(os.path.join(self.plot_dir, 'erg_distance_heatmap.png'), dpi=300)
            plt.close()

        # Analyze enrichment of specific patterns by treatment
        # Look for enrichment of specific region/context combinations
        if len(erg_variants) > 5:  # Need enough data for meaningful stats
            # Create combined region+context patterns
            erg_variants['Pattern'] = erg_variants['Regulatory_Region'] + '_' + erg_variants['Genetic_Context']

            # Count patterns by treatment
            pattern_counts = erg_variants.groupby(['Treatment', 'Pattern']).size().unstack(fill_value=0)

            # Calculate enrichment compared to average
            pattern_avg = pattern_counts.mean(axis=0)
            pattern_enrichment = pattern_counts.div(pattern_avg, axis=1)

            # Save enrichment data
            pattern_enrichment.to_csv(os.path.join(self.data_dir, 'erg_pattern_enrichment.tsv'), sep='\t')

            # Plot top enriched patterns by treatment
            # Get top 10 patterns by variance across treatments
            pattern_var = pattern_enrichment.var(axis=0).sort_values(ascending=False)
            top_patterns = pattern_var.head(min(10, len(pattern_var))).index.tolist()

            if top_patterns:
                plt.figure(figsize=(14, 10))
                enrichment_subset = pattern_enrichment[top_patterns]

                # Create heatmap
                sns.heatmap(enrichment_subset, annot=True, fmt='.2f', cmap='RdBu_r',
                           center=1.0, linewidths=0.5,
                           cbar_kws={'label': 'Fold Enrichment vs Average'})
                plt.title('Treatment-Specific Enrichment of Regulatory Patterns', fontsize=14)
                plt.xlabel('Regulatory Pattern', fontsize=12)
                plt.ylabel('Treatment', fontsize=12)
                plt.xticks(rotation=45, ha='right')
                plt.tight_layout()
                plt.savefig(os.path.join(self.plot_dir, 'erg_pattern_enrichment.png'), dpi=300)
                plt.close()

        return erg_variants, erg_region_counts, erg_context_counts
    
    def analyze_regulatory_by_treatment(self):
        """Analyze regulatory region patterns by treatment"""
        print("Analyzing regulatory region patterns by treatment...")
        
        # Overall counts
        region_counts = self.reg_variants.groupby(['Treatment', 'Regulatory_Region']).size().unstack(fill_value=0)
        region_percent = region_counts.div(region_counts.sum(axis=1), axis=0) * 100
        
        # Adaptation type comparison (Temperature vs Low Oxygen)
        adaptation_group = self.reg_variants.groupby(['Treatment_Group', 'Regulatory_Region']).size().unstack(fill_value=0)
        adaptation_percent = adaptation_group.div(adaptation_group.sum(axis=1), axis=0) * 100
        
        # Save results
        region_counts.to_csv(os.path.join(self.data_dir, 'treatment_regulatory_counts.tsv'), sep='\t')
        region_percent.to_csv(os.path.join(self.data_dir, 'treatment_regulatory_percent.tsv'), sep='\t')
        adaptation_group.to_csv(os.path.join(self.data_dir, 'adaptation_regulatory_counts.tsv'), sep='\t')
        adaptation_percent.to_csv(os.path.join(self.data_dir, 'adaptation_regulatory_percent.tsv'), sep='\t')
        
        # Plot by treatment
        plt.figure(figsize=(12, 8))
        region_percent.plot(kind='bar', stacked=True, colormap='viridis')
        plt.title('Distribution of Variants in Regulatory Regions by Treatment', fontsize=14)
        plt.xlabel('Treatment', fontsize=12)
        plt.ylabel('Percentage of Variants', fontsize=12)
        plt.xticks(rotation=45)
        plt.legend(title='Regulatory Region', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(os.path.join(self.plot_dir, 'treatment_regulatory_distribution.png'), dpi=300)
        plt.close()
        
        # Plot by adaptation type
        plt.figure(figsize=(12, 8))
        adaptation_percent.plot(kind='bar', stacked=True, colormap='viridis')
        plt.title('Distribution of Variants in Regulatory Regions by Adaptation Type', fontsize=14)
        plt.xlabel('Adaptation Type', fontsize=12)
        plt.ylabel('Percentage of Variants', fontsize=12)
        plt.xticks(rotation=45)
        plt.legend(title='Regulatory Region', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(os.path.join(self.plot_dir, 'adaptation_regulatory_distribution.png'), dpi=300)
        plt.close()
        
        return region_counts, adaptation_group
    
    def analyze_distance_distribution(self):
        """Analyze the distance distribution of regulatory variants"""
        print("Analyzing distance distribution of regulatory variants...")
        
        # Filter to focus on proximal regulatory variants (<10kb)
        proximal_variants = self.reg_variants[
            (self.reg_variants['Distance'] > -10000) & 
            (self.reg_variants['Distance'] < 10000)
        ].copy()
        
        # Create binned distances
        bins = np.arange(-10000, 10001, 500)
        proximal_variants['Distance_Bin'] = pd.cut(
            proximal_variants['Distance'], 
            bins=bins,
            labels=[f"{bins[i]}-{bins[i+1]}" for i in range(len(bins)-1)]
        )
        
        # Group by treatment and distance bin
        dist_counts = proximal_variants.groupby(['Treatment', 'Distance_Bin']).size().unstack(fill_value=0)
        
        # Save results
        proximal_variants.to_csv(os.path.join(self.data_dir, 'distance_binned_variants.tsv'), sep='\t', index=False)
        dist_counts.to_csv(os.path.join(self.data_dir, 'distance_distribution.tsv'), sep='\t')
        
        # Plot distance distribution
        plt.figure(figsize=(14, 8))
        for treatment in self.treatments:
            if treatment in dist_counts.index:
                plt.plot(range(len(dist_counts.columns)), 
                        dist_counts.loc[treatment], 
                        label=treatment,
                        color=self.colors.get(treatment, 'gray'),
                        marker='o', markersize=4)
        
        plt.title('Distribution of Regulatory Variants by Distance from Gene', fontsize=14)
        plt.xlabel('Distance from Gene (bp)', fontsize=12)
        plt.ylabel('Number of Variants', fontsize=12)
        plt.xticks(range(0, len(dist_counts.columns), 4), 
                  [dist_counts.columns[i] for i in range(0, len(dist_counts.columns), 4)], 
                  rotation=45, ha='right')
        plt.legend(title='Treatment')
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        plt.savefig(os.path.join(self.plot_dir, 'distance_distribution.png'), dpi=300)
        plt.close()
        
        return proximal_variants, dist_counts
    
    def analyze_by_gene_type(self):
        """Analyze regulatory patterns by gene type (ERG vs non-ERG)"""
        print("Analyzing regulatory patterns by gene type...")
        
        # Get list of ERG genes
        erg_gene_ids = set(self.erg_genes['w303_gene_id'].unique())
        
        # Add gene type column
        self.reg_variants['Gene_Type'] = self.reg_variants['Gene_ID'].apply(
            lambda x: 'ERG' if x in erg_gene_ids else 'non-ERG'
        )
        
        # Group by gene type and regulatory region
        gene_type_counts = self.reg_variants.groupby(['Treatment', 'Gene_Type', 'Regulatory_Region']).size().unstack(fill_value=0)
        
        # Save results
        gene_type_counts.to_csv(os.path.join(self.data_dir, 'gene_type_regulatory_counts.tsv'), sep='\t')
        
        # Plot grouped by gene type
        for gene_type in ['ERG', 'non-ERG']:
            if gene_type in [idx[1] for idx in gene_type_counts.index]:
                # Filter for this gene type
                type_data = gene_type_counts.xs(gene_type, level='Gene_Type')
                
                # Calculate percentages
                type_percent = type_data.div(type_data.sum(axis=1), axis=0) * 100
                
                plt.figure(figsize=(12, 8))
                type_percent.plot(kind='bar', stacked=True, colormap='viridis')
                plt.title(f'Distribution of Variants in Regulatory Regions - {gene_type} Genes', fontsize=14)
                plt.xlabel('Treatment', fontsize=12)
                plt.ylabel('Percentage of Variants', fontsize=12)
                plt.xticks(rotation=45)
                plt.legend(title='Regulatory Region', bbox_to_anchor=(1.05, 1), loc='upper left')
                plt.tight_layout()
                plt.savefig(os.path.join(self.plot_dir, f'{gene_type.lower()}_regulatory_distribution.png'), dpi=300)
                plt.close()
        
        return gene_type_counts
    
    def compare_erg_vs_nonerg(self):
        """Compare regulatory patterns between ERG and non-ERG genes"""
        print("Comparing ERG vs non-ERG regulatory patterns...")
        
        # Get regulatory region proportions for both gene types
        erg_gene_ids = set(self.erg_genes['w303_gene_id'].unique())
        erg_variants = self.reg_variants[self.reg_variants['Gene_ID'].isin(erg_gene_ids)]
        nonerg_variants = self.reg_variants[~self.reg_variants['Gene_ID'].isin(erg_gene_ids)]
        
        # Calculate proportions
        erg_props = erg_variants.groupby('Regulatory_Region').size() / len(erg_variants) * 100
        nonerg_props = nonerg_variants.groupby('Regulatory_Region').size() / len(nonerg_variants) * 100
        
        # Combine for comparison
        comparison = pd.DataFrame({
            'ERG': erg_props,
            'non-ERG': nonerg_props
        }).fillna(0)
        
        # Save results
        comparison.to_csv(os.path.join(self.data_dir, 'erg_vs_nonerg_comparison.tsv'), sep='\t')
        
        # Plot comparison
        plt.figure(figsize=(10, 8))
        comparison.plot(kind='bar', width=0.8)
        plt.title('Comparison of Regulatory Region Distribution: ERG vs non-ERG Genes', fontsize=14)
        plt.xlabel('Regulatory Region', fontsize=12)
        plt.ylabel('Percentage of Variants', fontsize=12)
        plt.xticks(rotation=45, ha='right')
        plt.legend(title='Gene Type')
        plt.grid(True, linestyle='--', alpha=0.5, axis='y')
        plt.tight_layout()
        plt.savefig(os.path.join(self.plot_dir, 'erg_vs_nonerg_comparison.png'), dpi=300)
        plt.close()
        
        return comparison
    
    def generate_summary(self):
        """Generate a summary of the regulatory region mapping analysis"""
        print("Generating summary of regulatory region mapping...")
        
        summary = {
            "total_regulatory_variants": len(self.reg_variants),
            "erg_regulatory_variants": len(self.reg_variants[self.reg_variants['Gene_ID'].isin(self.erg_genes['w303_gene_id'])]),
            "regulatory_region_counts": self.reg_variants['Regulatory_Region'].value_counts().to_dict(),
            "treatment_counts": self.reg_variants['Treatment'].value_counts().to_dict(),
            "top_regulatory_regions": self.reg_variants['Regulatory_Region'].value_counts().nlargest(3).to_dict(),
            "core_promoter_percent": round(
                len(self.reg_variants[self.reg_variants['Regulatory_Region'] == 'core_promoter']) / 
                len(self.reg_variants) * 100, 2
            ),
            "region_by_treatment": {
                treatment: self.reg_variants[self.reg_variants['Treatment'] == treatment]['Regulatory_Region'].value_counts().to_dict()
                for treatment in self.treatments if treatment in self.reg_variants['Treatment'].unique()
            }
        }
        
        # Save summary
        with open(os.path.join(self.data_dir, 'regulatory_mapping_summary.json'), 'w') as f:
            json.dump(summary, f, indent=2)
        
        # Create summary text file
        with open(os.path.join(self.output_dir, 'regulatory_mapping_report.txt'), 'w') as f:
            f.write("# Regulatory Region Mapping Analysis Report\n\n")
            f.write(f"Total regulatory variants analyzed: {summary['total_regulatory_variants']}\n")
            f.write(f"Ergosterol gene regulatory variants: {summary['erg_regulatory_variants']}\n\n")
            
            f.write("## Regulatory Region Distribution\n\n")
            for region, count in sorted(summary['regulatory_region_counts'].items(), key=lambda x: x[1], reverse=True):
                percent = round(count / summary['total_regulatory_variants'] * 100, 2)
                f.write(f"{region}: {count} variants ({percent}%)\n")
            
            f.write("\n## Treatment Distribution\n\n")
            for treatment, count in sorted(summary['treatment_counts'].items()):
                percent = round(count / summary['total_regulatory_variants'] * 100, 2)
                f.write(f"{treatment}: {count} variants ({percent}%)\n")
            
            f.write("\n## Key Findings\n\n")
            f.write(f"- Core promoter variants represent {summary['core_promoter_percent']}% of all regulatory variants\n")
            f.write("- Top regulatory regions by frequency:\n")
            for region, count in summary['top_regulatory_regions'].items():
                percent = round(count / summary['total_regulatory_variants'] * 100, 2)
                f.write(f"  * {region}: {percent}%\n")
            
            f.write("\n## Treatment-Specific Patterns\n\n")
            for treatment, regions in summary['region_by_treatment'].items():
                f.write(f"### {treatment}\n")
                total = sum(regions.values())
                for region, count in sorted(regions.items(), key=lambda x: x[1], reverse=True)[:3]:
                    percent = round(count / total * 100, 2)
                    f.write(f"- {region}: {percent}%\n")
                f.write("\n")
        
        print(f"Summary saved to {os.path.join(self.output_dir, 'regulatory_mapping_report.txt')}")
        return summary
    
    def run_analysis(self):
        """Run the complete regulatory region mapping analysis pipeline"""
        # Prepare data
        self.prepare_data()

        # Check if we have data to analyze
        if self.reg_variants.empty:
            print("No regulatory variants to analyze. Exiting.")
            return None

        try:
            # Run analyses
            self.analyze_erg_regulatory_regions()
            self.analyze_regulatory_by_treatment()
            self.analyze_distance_distribution()
            self.analyze_by_gene_type()
            self.compare_erg_vs_nonerg()

            # Generate summary
            summary = self.generate_summary()

            print("Regulatory region mapping analysis complete!")
            return summary
        except Exception as e:
            print(f"Error during analysis: {e}")
            import traceback
            traceback.print_exc()
            return None


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
    
    print(f"Starting regulatory region mapping analysis...")
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