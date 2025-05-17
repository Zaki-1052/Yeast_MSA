#!/usr/bin/env python3
# /Users/zakiralibhai/Documents/GitHub/Yeast_MSA/scripts/regulatory_analysis/accurate_regulatory_analysis.py

"""
Accurate regulatory region analysis for yeast variants.
Properly analyzes the unique variants and maps them to regulatory regions.
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

class RegulatoryAnalyzer:
    """Accurately analyzes regulatory regions in yeast variants"""
    
    def __init__(self, output_dir, variants_file, gene_mapping_file, erg_genes_file):
        """Initialize with file paths and configurations"""
        # Setup directories
        self.output_dir = self._ensure_dir(output_dir)
        self.data_dir = self._ensure_dir(os.path.join(output_dir, 'data'))
        self.plots_dir = self._ensure_dir(os.path.join(output_dir, 'plots'))
        
        # Load data files
        print("Loading input files...")
        self.variants = self._load_file(variants_file, "variants")
        self.gene_mapping = self._load_file(gene_mapping_file, "gene mapping")
        self.erg_genes = self._load_file(erg_genes_file, "ERG genes")
        
        # Define yeast regulatory regions with accurate ranges
        self.define_regulatory_regions()
        
        # Store analysis results
        self.unique_variants = None
        self.classified_variants = None
        self.results = {}
    
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
            if description in ["variants", "gene mapping"]:
                print("Critical file missing. Exiting.")
                sys.exit(1)
            return pd.DataFrame()  # Empty dataframe for non-critical files
    
    def define_regulatory_regions(self):
        """Define yeast-specific regulatory regions with precise boundaries"""
        # Define regions based on distance from gene
        self.regulatory_regions = {
            'core_promoter': (0, 150),          # Core promoter (0-150bp upstream)
            'TATA_box_region': (40, 120),       # TATA box region
            'proximal_UAS': (150, 500),         # Proximal UAS (150-500bp upstream)
            'distal_UAS': (500, 1500),          # Distal UAS (500-1500bp upstream)
            'far_upstream': (1500, float('inf')) # Far upstream elements
        }
        
        # Define ERG gene zones based on distance
        self.erg_zones = {
            'core_zone': (0, 0),                  # The gene itself
            'buffer_zone': (0, 5000),             # Buffer zone (<5kb)
            'intermediate_zone': (5000, 50000),   # Intermediate zone (5-50kb)
            'satellite_zone': (50000, 100000)     # Satellite zone (50-100kb)
        }
    
    def extract_unique_variants(self):
        """Extract and analyze unique variants from the dataset"""
        print("\nExtracting unique variants...")
        
        if self.variants is None or len(self.variants) == 0:
            print("No variants to analyze")
            return None
        
        # Define columns that determine uniqueness
        uniqueness_columns = ['Scaffold', 'Position', 'Ref', 'Alt', 'Gene_ID']
        
        # Extract unique variants
        unique_variants = self.variants.drop_duplicates(subset=uniqueness_columns)
        
        print(f"Found {len(unique_variants)} unique variants out of {len(self.variants)} total variants")
        
        # Store unique variants
        self.unique_variants = unique_variants
        
        # Show details of unique variants
        print("\nUnique variants:")
        for _, row in unique_variants.iterrows():
            gene_id = row['Gene_ID']
            distance = row['Distance']
            effect = row['Effect']
            print(f"Gene: {gene_id}, Distance: {distance}bp, Effect: {effect}")
        
        # Analyze distance distribution
        distances = unique_variants['Distance'].dropna()
        print(f"\nDistance range: {distances.min()} to {distances.max()} bp")
        print(f"Unique distance values: {sorted(distances.unique())}")
        
        return unique_variants
    
    def classify_variants(self):
        """Classify variants into yeast-specific regulatory regions"""
        print("\nClassifying variants into yeast regulatory regions...")
        
        if self.unique_variants is None:
            print("No unique variants to classify. Run extract_unique_variants() first.")
            return None
        
        # Create a copy to work with
        classified = self.unique_variants.copy()
        
        # Add classification columns
        classified['Regulatory_Region'] = 'unknown'
        classified['ERG_Zone'] = 'unknown'
        
        # Classify each variant
        for idx, row in classified.iterrows():
            # Get distance and effect
            distance = row['Distance']
            effect = row['Effect']
            
            # Skip if missing data
            if pd.isna(distance) or pd.isna(effect):
                continue
                
            # For upstream variants, classify by distance
            if 'upstream' in effect.lower():
                for region, (min_dist, max_dist) in self.regulatory_regions.items():
                    if min_dist <= distance < max_dist:
                        classified.at[idx, 'Regulatory_Region'] = region
                        break
            
            # Check if this is an ERG gene variant
            gene_id = row['Gene_ID']
            if 'Gene_ID' in row and gene_id in self.erg_genes['w303_gene_id'].values:
                # It's an ERG gene - classify by distance
                for zone, (min_dist, max_dist) in self.erg_zones.items():
                    if min_dist <= distance < max_dist:
                        classified.at[idx, 'ERG_Zone'] = zone
                        break
        
        # Store classified variants
        self.classified_variants = classified
        
        # Show classification results
        print("\nRegulatory region classification:")
        region_counts = classified['Regulatory_Region'].value_counts()
        for region, count in region_counts.items():
            percent = count / len(classified) * 100
            print(f"  {region}: {count} variants ({percent:.1f}%)")
        
        # Show ERG zone classification if applicable
        erg_variants = classified[classified['ERG_Zone'] != 'unknown']
        if len(erg_variants) > 0:
            print("\nERG zone classification:")
            zone_counts = erg_variants['ERG_Zone'].value_counts()
            for zone, count in zone_counts.items():
                print(f"  {zone}: {count} variants")
        
        return classified
    
    def apply_classification_to_full_dataset(self):
        """Apply classification from unique variants to the full dataset"""
        print("\nApplying classification to full dataset...")
        
        if self.classified_variants is None:
            print("No classified variants. Run classify_variants() first.")
            return None
        
        # Create a mapping from variant key to classification
        variant_to_region = {}
        variant_to_erg_zone = {}
        
        for _, row in self.classified_variants.iterrows():
            # Create a key that uniquely identifies this variant
            key = (row['Scaffold'], row['Position'], row['Ref'], row['Alt'], row['Gene_ID'])
            variant_to_region[key] = row['Regulatory_Region']
            variant_to_erg_zone[key] = row['ERG_Zone']
        
        # Apply classification to full dataset
        full_dataset = self.variants.copy()
        
        # Add classification columns
        full_dataset['Regulatory_Region'] = 'unknown'
        full_dataset['ERG_Zone'] = 'unknown'
        
        # Classify each variant in the full dataset
        for idx, row in full_dataset.iterrows():
            key = (row['Scaffold'], row['Position'], row['Ref'], row['Alt'], row['Gene_ID'])
            if key in variant_to_region:
                full_dataset.at[idx, 'Regulatory_Region'] = variant_to_region[key]
            if key in variant_to_erg_zone:
                full_dataset.at[idx, 'ERG_Zone'] = variant_to_erg_zone[key]
        
        # Calculate region counts for full dataset
        region_counts = full_dataset['Regulatory_Region'].value_counts()
        region_percent = (region_counts / len(full_dataset) * 100).round(1)
        
        print("\nRegulatory region distribution in full dataset:")
        for region, count in region_counts.items():
            percent = region_percent[region]
            print(f"  {region}: {count} variants ({percent}%)")
        
        # Save full classified dataset
        full_dataset.to_csv(os.path.join(self.data_dir, 'classified_variants.tsv'), sep='\t', index=False)
        
        # Store region counts in results
        self.results['region_counts'] = region_counts.to_dict()
        self.results['region_percent'] = region_percent.to_dict()
        
        # Add full dataset to instance for further analysis
        self.full_classified = full_dataset
        
        return full_dataset
    
    def analyze_treatment_distribution(self):
        """Analyze the distribution of regulatory regions by treatment"""
        print("\nAnalyzing regulatory region distribution by treatment...")
        
        if not hasattr(self, 'full_classified') or self.full_classified is None:
            print("No classified full dataset. Run apply_classification_to_full_dataset() first.")
            return None
        
        if 'Treatment' not in self.full_classified.columns:
            print("Treatment information not available in the dataset.")
            return None
        
        # Create a crosstab of treatment vs. regulatory region
        treatment_region = pd.crosstab(
            self.full_classified['Treatment'],
            self.full_classified['Regulatory_Region'],
            normalize='index'
        ) * 100  # Convert to percentages
        
        # Store in results
        self.results['treatment_region'] = treatment_region.to_dict()
        
        # Save treatment-region distribution
        treatment_region.to_csv(os.path.join(self.data_dir, 'treatment_region_distribution.tsv'), sep='\t')
        
        # Show treatment distribution
        print("\nRegulatory region distribution by treatment:")
        print(treatment_region.round(1))
        
        # Plot as heatmap
        plt.figure(figsize=(12, 8))
        sns.heatmap(treatment_region, annot=True, fmt='.1f', cmap='viridis',
                   linewidths=0.5, cbar_kws={'label': 'Percentage of Variants'})
        plt.title('Regulatory Region Distribution by Treatment', fontsize=14)
        plt.xlabel('Regulatory Region', fontsize=12)
        plt.ylabel('Treatment', fontsize=12)
        plt.tight_layout()
        plt.savefig(os.path.join(self.plots_dir, 'treatment_region_heatmap.png'), dpi=300)
        plt.close()
        
        return treatment_region
    
    def analyze_erg_genes(self):
        """Analyze regulatory regions specifically around ERG genes"""
        print("\nAnalyzing regulatory regions around ERG genes...")
        
        if not hasattr(self, 'full_classified') or self.full_classified is None:
            print("No classified full dataset. Run apply_classification_to_full_dataset() first.")
            return None
        
        # Filter for variants near ERG genes
        erg_gene_ids = set(self.erg_genes['w303_gene_id'].values)
        erg_variants = self.full_classified[self.full_classified['Gene_ID'].isin(erg_gene_ids)]
        
        if len(erg_variants) == 0:
            print("No variants found near ERG genes.")
            return None
        
        print(f"Found {len(erg_variants)} variants near ERG genes")
        
        # Save ERG variants
        erg_variants.to_csv(os.path.join(self.data_dir, 'erg_gene_variants.tsv'), sep='\t', index=False)
        
        # Calculate region counts for ERG genes
        erg_region_counts = erg_variants['Regulatory_Region'].value_counts()
        erg_region_percent = (erg_region_counts / len(erg_variants) * 100).round(1)
        
        print("\nRegulatory region distribution around ERG genes:")
        for region, count in erg_region_counts.items():
            percent = erg_region_percent[region]
            print(f"  {region}: {count} variants ({percent}%)")
        
        # Store in results
        self.results['erg_region_counts'] = erg_region_counts.to_dict()
        self.results['erg_region_percent'] = erg_region_percent.to_dict()
        
        # Analyze by treatment if available
        if 'Treatment' in erg_variants.columns:
            erg_treatment_region = pd.crosstab(
                erg_variants['Treatment'],
                erg_variants['Regulatory_Region'],
                normalize='index'
            ) * 100  # Convert to percentages
            
            # Store in results
            self.results['erg_treatment_region'] = erg_treatment_region.to_dict()
            
            # Save ERG treatment-region distribution
            erg_treatment_region.to_csv(os.path.join(self.data_dir, 'erg_treatment_region.tsv'), sep='\t')
            
            # Plot as heatmap
            plt.figure(figsize=(12, 8))
            sns.heatmap(erg_treatment_region, annot=True, fmt='.1f', cmap='viridis',
                       linewidths=0.5, cbar_kws={'label': 'Percentage of Variants'})
            plt.title('Regulatory Region Distribution by Treatment (ERG Genes)', fontsize=14)
            plt.xlabel('Regulatory Region', fontsize=12)
            plt.ylabel('Treatment', fontsize=12)
            plt.tight_layout()
            plt.savefig(os.path.join(self.plots_dir, 'erg_treatment_region_heatmap.png'), dpi=300)
            plt.close()
        
        return erg_variants
    
    def generate_summary_report(self):
        """Generate a comprehensive summary report of findings"""
        print("\nGenerating summary report...")
        
        if not self.results:
            print("No analysis results to summarize.")
            return None
        
        # Create report content
        report = ["# Accurate Regulatory Region Analysis Report\n"]
        
        # Overview section
        report.append("## Overview\n")
        report.append(f"Total variants in dataset: {len(self.variants)}")
        report.append(f"Unique variants: {len(self.unique_variants)}")
        
        if 'erg_region_counts' in self.results:
            erg_count = sum(self.results['erg_region_counts'].values())
            report.append(f"Variants near ERG genes: {erg_count}")
        
        report.append("\n## Unique Regulatory Variants\n")
        for _, row in self.classified_variants.iterrows():
            gene_id = row['Gene_ID']
            distance = row['Distance']
            region = row['Regulatory_Region']
            report.append(f"- Gene: {gene_id}, Distance: {distance}bp -> {region}")
        
        # Regulatory region distribution
        if 'region_counts' in self.results:
            report.append("\n## Regulatory Region Distribution\n")
            for region, count in sorted(self.results['region_counts'].items(), key=lambda x: x[1], reverse=True):
                percent = self.results['region_percent'][region]
                report.append(f"{region}: {count} variants ({percent}%)")
        
        # Treatment distribution (treatment-wise)
        if 'treatment_region' in self.results:
            report.append("\n## Treatment-Specific Patterns\n")
            
            # Extract treatment-wise distributions
            treatments = self.full_classified['Treatment'].unique()
            
            for treatment in treatments:
                treatment_data = self.full_classified[self.full_classified['Treatment'] == treatment]
                if len(treatment_data) > 0:
                    report.append(f"### {treatment}")
                    
                    # Calculate region distribution for this treatment
                    region_counts = treatment_data['Regulatory_Region'].value_counts()
                    region_percent = (region_counts / len(treatment_data) * 100).round(1)
                    
                    for region, percent in region_percent.items():
                        report.append(f"- {region}: {percent}%")
                    
                    report.append("")  # Add blank line between treatments
        
        # ERG gene analysis
        if 'erg_region_counts' in self.results:
            report.append("\n## ERG Gene Regulatory Analysis\n")
            
            for region, count in sorted(self.results['erg_region_counts'].items(), key=lambda x: x[1], reverse=True):
                percent = self.results['erg_region_percent'][region]
                report.append(f"{region}: {count} variants ({percent}%)")
        
        # Biological interpretations
        report.append("\n## Biological Interpretations\n")
        
        # Add interpretations based on results
        if 'region_percent' in self.results:
            region_percent = self.results['region_percent']
            
            # Interpret core promoter vs. UAS
            core_percent = region_percent.get('core_promoter', 0)
            proximal_uas = region_percent.get('proximal_UAS', 0)
            distal_uas = region_percent.get('distal_UAS', 0)
            uas_total = proximal_uas + distal_uas
            
            if core_percent > uas_total:
                report.append("- **Core Promoter Preference**: A higher proportion of variants occur in core promoter regions, " +
                             "suggesting adaptation through changes in basal transcription machinery and RNA polymerase II recruitment.")
            else:
                report.append("- **UAS Region Preference**: A higher proportion of variants occur in Upstream Activating Sequence regions, " +
                             "suggesting adaptation through changes in transcription factor binding and regulatory protein interactions.")
            
            # Interpret TATA box region
            tata_percent = region_percent.get('TATA_box_region', 0)
            if tata_percent > 0:
                report.append("- **TATA Box Impacts**: Variants in the TATA box region (40-120bp upstream) may directly affect " +
                             "transcription factor IID (TFIID) binding and transcription initiation efficiency.")
            
            # Interpret distal UAS
            if distal_uas > 10:
                report.append("- **Distal Regulation**: The presence of variants in distal UAS regions suggests " +
                             "long-range regulatory mechanisms, potentially involving chromatin looping or insulator elements.")
        
        # Write report to file
        report_path = os.path.join(self.output_dir, 'accurate_regulatory_analysis_report.txt')
        with open(report_path, 'w') as f:
            f.write('\n'.join(report))
        
        print(f"Summary report saved to {report_path}")
        
        # Save results as JSON
        results_path = os.path.join(self.data_dir, 'regulatory_analysis_results.json')
        
        # Convert results to JSON-serializable format
        json_results = {}
        for key, value in self.results.items():
            if isinstance(value, dict):
                # Convert nested dictionaries
                json_results[key] = {str(k): v for k, v in value.items()}
            else:
                json_results[key] = value
        
        with open(results_path, 'w') as f:
            json.dump(json_results, f, indent=2)
        
        print(f"Analysis results saved to {results_path}")
        
        return report
    
    def run_analysis(self):
        """Run the complete regulatory region analysis pipeline"""
        print("\n=== Starting Accurate Regulatory Region Analysis ===\n")
        
        # Step 1: Extract unique variants
        self.extract_unique_variants()
        
        # Step 2: Classify unique variants
        self.classify_variants()
        
        # Step 3: Apply classification to full dataset
        self.apply_classification_to_full_dataset()
        
        # Step 4: Analyze treatment distribution
        self.analyze_treatment_distribution()
        
        # Step 5: Analyze ERG genes
        self.analyze_erg_genes()
        
        # Step 6: Generate summary report
        self.generate_summary_report()
        
        print("\n=== Accurate Regulatory Region Analysis Complete ===")
        
        return self.results

def main():
    """Main function to run the analysis"""
    parser = argparse.ArgumentParser(description='Accurate regulatory region mapping and analysis')
    
    # Required arguments
    parser.add_argument('--output-dir', 
                       default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/regulatory_analysis/accurate_mapping',
                       help='Output directory for results')
    parser.add_argument('--variants', 
                       default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/gene_variants/all_gene_variants.tsv',
                       help='Variants TSV file')
    parser.add_argument('--gene-mapping', 
                       default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/reference/gene_mapping_full.tsv',
                       help='Gene mapping file with coordinates')
    parser.add_argument('--erg-genes', 
                       default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/reference/genes_of_interest_mapping.tsv',
                       help='ERG genes mapping file')
    
    args = parser.parse_args()
    
    # Initialize analyzer
    analyzer = RegulatoryAnalyzer(
        args.output_dir,
        args.variants,
        args.gene_mapping,
        args.erg_genes
    )
    
    # Run complete analysis
    analyzer.run_analysis()

if __name__ == "__main__":
    main()