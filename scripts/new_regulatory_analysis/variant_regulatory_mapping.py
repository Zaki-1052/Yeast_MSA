#!/usr/bin/env python3
# /Users/zakiralibhai/Documents/GitHub/Yeast_MSA/scripts/new_regulatory_analysis/variant_regulatory_mapping.py

"""
Advanced Variant-to-Regulatory Feature Mapping System for the Yeast MSA project.

This script maps variants to specific regulatory features with higher precision,
categorizes variants by their potential regulatory impact, and
quantifies distance relationships to transcription start sites.

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
from Bio import SeqIO
import re


class VariantRegulatoryMapper:
    """Maps variants to regulatory features with high precision"""
    
    def __init__(self, variants_file, regulatory_map_file, reference_genome=None, 
                 regulatory_features_file=None, output_dir=None, motif_database=None):
        """Initialize with necessary files and configurations"""
        # Setup directories
        self.output_dir = self._ensure_dir(output_dir or '/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/new_regulatory_analysis')
        self.data_dir = self._ensure_dir(os.path.join(self.output_dir, 'data'))
        self.plot_dir = self._ensure_dir(os.path.join(self.output_dir, 'plots'))
        
        # Load input files
        print("Loading input files...")
        
        # Load variants
        self.variants = self._load_file(variants_file, "variants")
        
        # Load gene regulatory map
        try:
            with open(regulatory_map_file, 'r') as f:
                self.gene_regulatory_map = json.load(f)
            print(f"Loaded gene regulatory map with {len(self.gene_regulatory_map)} genes")
        except Exception as e:
            print(f"Error loading gene regulatory map: {e}")
            sys.exit(1)
        
        # Load regulatory features if provided
        self.regulatory_features = None
        if regulatory_features_file:
            try:
                with open(regulatory_features_file, 'r') as f:
                    self.regulatory_features = json.load(f)
                print(f"Loaded regulatory features with {len(self.regulatory_features)} categories")
            except Exception as e:
                print(f"Error loading regulatory features: {e}")
                print("Using default regulatory features from gene regulatory map")
        
        # Load reference genome for sequence context if provided
        self.genome_seq = None
        if reference_genome:
            try:
                self.genome_seq = SeqIO.to_dict(SeqIO.parse(reference_genome, "fasta"))
                print(f"Loaded reference genome with {len(self.genome_seq)} sequences")
            except Exception as e:
                print(f"Error loading reference genome: {e}")
                print("Sequence context analysis will be limited")
        
        # Load motif database if provided
        self.motifs = self._load_motifs(motif_database)
        
        # Results storage
        self.mapped_variants = None
        self.region_variants = {}
        self.tfbs_variants = {}
        self.mapping_statistics = {}
    
    def _ensure_dir(self, directory):
        """Create directory if it doesn't exist"""
        if not os.path.exists(directory):
            os.makedirs(directory)
            print(f"Created directory: {directory}")
        return directory
    
    def _load_file(self, file_path, description):
        """Load a data file with appropriate format handling"""
        try:
            if file_path.endswith(('.tsv', '.txt')):
                data = pd.read_csv(file_path, sep='\t')
            elif file_path.endswith('.csv'):
                data = pd.read_csv(file_path)
            else:
                # Default to TSV
                data = pd.read_csv(file_path, sep='\t')
            
            print(f"Loaded {description} file with {len(data)} entries")
            return data
        except Exception as e:
            print(f"Error loading {description} file ({file_path}): {e}")
            if description in ["variants"]:
                print("Critical file missing. Exiting.")
                sys.exit(1)
            return pd.DataFrame()
    
    def _load_motifs(self, motif_file):
        """Load transcription factor binding motifs"""
        motifs = {}
        
        if motif_file and os.path.exists(motif_file):
            try:
                with open(motif_file, 'r') as f:
                    motifs = json.load(f)
                print(f"Loaded motif database with {len(motifs)} TF binding motifs")
            except Exception as e:
                print(f"Error loading motif database: {e}")
                print("Using built-in set of common yeast motifs")
                motifs = self._get_default_motifs()
        else:
            print("Using built-in set of common yeast motifs")
            motifs = self._get_default_motifs()
        
        return motifs
    
    def _get_default_motifs(self):
        """Return a default set of common yeast transcription factor motifs"""
        # Common yeast motifs
        return {
            # Core promoter elements
            "TATA-box": "TATA[AT]A[AT]",       # TATA box consensus
            "Initiator": "YYANWYY",             # Initiator element
            
            # General transcription factors
            "GCN4": "TGA[CG]TCA",              # Amino acid biosynthesis
            "Rap1": "ACACCC[AG]TACATYW",       # Ribosomal protein genes regulator
            "Abf1": "TCACA[TC][CT][TA][AC]",   # ARS binding factor
            "Reb1": "TTACCCG",                 # RNA polymerase I enhancer binding protein
            "Cbf1": "TCACGTG",                 # Centromere binding factor
            "Gal4": "CGG[ACGT]{5}CCG",         # Galactose metabolism
            
            # Stress response elements
            "STRE": "AGGGG",                   # Stress response element
            "HSE": "NGAAN",                    # Heat shock element
            
            # Ergosterol pathway specific
            "Ecm22/Upc2": "TCGTATA",           # Sterol regulatory element (SRE)
            "Hap1": "CGGNNNTANCGG",            # Heme activator protein - oxygen sensing
            "Rox1": "YYNATTGTTY",              # Repressor of hypoxic genes
            
            # Other important TFs
            "Msn2/Msn4": "CCCCT",              # General stress response
            "Pdr1/Pdr3": "TCCGCGGA",           # Pleiotropic drug resistance
            "Yap1": "TTACTAA",                 # Oxidative stress response
            "Hsf1": "NGAANNTTCN",              # Heat shock response
            "Pho4": "CACGTK",                  # Phosphate response
            "Ste12": "TGAAACA",                # Mating response
            "Mcm1": "CCCWWWWWGG",              # Cell cycle and mating
            "Swi4/Swi6": "CACGAAA"             # Cell cycle regulation
        }
    
    def preprocess_variants(self):
        """Preprocess variants for regulatory mapping"""
        print("\nPreprocessing variants for regulatory mapping...")
        
        if len(self.variants) == 0:
            print("No variants to process.")
            return None
        
        # Check required columns
        required_cols = ['Scaffold', 'Position', 'Ref', 'Alt', 'Gene_ID']
        missing_cols = [col for col in required_cols if col not in self.variants.columns]
        
        if missing_cols:
            print(f"Error: Missing required columns: {missing_cols}")
            print("These columns are necessary for regulatory mapping.")
            return None
        
        # Display sample variants
        print("\nSample variants:")
        sample_cols = ['Gene_ID', 'Scaffold', 'Position', 'Ref', 'Alt']
        if 'Effect' in self.variants.columns:
            sample_cols.append('Effect')
        if 'Impact' in self.variants.columns:
            sample_cols.append('Impact')
        if 'Distance' in self.variants.columns:
            sample_cols.append('Distance')
        if 'Treatment' in self.variants.columns:
            sample_cols.append('Treatment')
        
        print(self.variants[sample_cols].head(5).to_string())
        
        # Basic statistics
        print("\nVariant counts:")
        print(f"Total variants: {len(self.variants)}")
        
        if 'Impact' in self.variants.columns:
            impact_counts = self.variants['Impact'].value_counts()
            print("\nImpact distribution:")
            for impact, count in impact_counts.items():
                print(f"  {impact}: {count} variants ({count/len(self.variants)*100:.1f}%)")
        
        if 'Effect' in self.variants.columns:
            effect_counts = self.variants['Effect'].value_counts()
            print("\nTop 5 effect types:")
            for effect, count in effect_counts.nlargest(5).items():
                print(f"  {effect}: {count} variants ({count/len(self.variants)*100:.1f}%)")
        
        if 'Treatment' in self.variants.columns:
            treatment_counts = self.variants['Treatment'].value_counts()
            print("\nTreatment distribution:")
            for treatment, count in treatment_counts.items():
                print(f"  {treatment}: {count} variants ({count/len(self.variants)*100:.1f}%)")
        
        # Copy variants
        processed_variants = self.variants.copy()
        
        # Ensure Distance is numeric
        if 'Distance' in processed_variants.columns:
            processed_variants['Distance'] = pd.to_numeric(processed_variants['Distance'], errors='coerce')
        
        # Create a unique identifier for each variant
        processed_variants['Variant_ID'] = processed_variants.apply(
            lambda row: f"{row['Scaffold']}:{row['Position']}:{row['Ref']}:{row['Alt']}", axis=1
        )
        
        # Store preprocessed variants
        return processed_variants
    
    def map_variants_to_regulatory_regions(self, processed_variants=None):
        """Map variants to regulatory regions and conservation zones"""
        print("\nMapping variants to regulatory regions and conservation zones...")
        
        if processed_variants is None:
            processed_variants = self.preprocess_variants()
            
        if processed_variants is None or len(processed_variants) == 0:
            print("No variants to map.")
            return None
        
        # Add regulatory columns
        processed_variants['Regulatory_Region'] = 'unknown'
        processed_variants['Regulatory_Type'] = 'unknown'
        processed_variants['Conservation_Zone'] = 'unknown'
        processed_variants['TFBS'] = None
        processed_variants['TFBS_Impact'] = None
        processed_variants['Motif_Disruption_Score'] = None
        
        # Track mapping statistics
        mapped_to_region = 0
        mapped_to_zone = 0
        
        # Map each variant
        print("Mapping variants to regulatory features...")
        for idx, variant in processed_variants.iterrows():
            # Get gene ID
            gene_id = variant['Gene_ID']
            
            # Skip if gene not in regulatory map
            if gene_id not in self.gene_regulatory_map:
                continue
            
            # Get gene info from regulatory map
            gene_info = self.gene_regulatory_map[gene_id]
            
            # Check regulatory regions
            for region, region_info in gene_info.get('regulatory_regions', {}).items():
                # Get region boundaries
                region_start = region_info.get('start', 0)
                region_end = region_info.get('end', 0)
                region_type = region_info.get('type', 'unknown')
                
                # Check if variant position is within region
                if region_start <= variant['Position'] <= region_end:
                    processed_variants.at[idx, 'Regulatory_Region'] = region
                    processed_variants.at[idx, 'Regulatory_Type'] = region_type
                    mapped_to_region += 1
                    break
            
            # Check conservation zones
            for zone, zone_info in gene_info.get('conservation_zones', {}).items():
                # Get zone boundaries
                zone_start = zone_info.get('start', 0)
                zone_end = zone_info.get('end', 0)
                
                # Check if variant position is within zone
                if zone_start <= variant['Position'] <= zone_end:
                    processed_variants.at[idx, 'Conservation_Zone'] = zone
                    mapped_to_zone += 1
                    break
        
        # Get transcription factor binding sites if reference genome is available
        if self.genome_seq is not None and self.motifs:
            print("Analyzing transcription factor binding sites...")
            processed_variants = self._analyze_tfbs(processed_variants)
        
        # Summarize mapping results
        total_variants = len(processed_variants)
        print(f"\nMapping results:")
        print(f"Total variants: {total_variants}")
        print(f"Mapped to regulatory regions: {mapped_to_region} ({mapped_to_region/total_variants*100:.1f}%)")
        print(f"Mapped to conservation zones: {mapped_to_zone} ({mapped_to_zone/total_variants*100:.1f}%)")
        
        # Count by region type
        region_counts = processed_variants['Regulatory_Region'].value_counts()
        region_percent = (region_counts / total_variants * 100).round(1)
        
        print("\nRegulatory region distribution:")
        for region, count in region_counts.nlargest(10).items():
            print(f"  {region}: {count} variants ({region_percent[region]}%)")
        
        # Count by conservation zone
        zone_counts = processed_variants['Conservation_Zone'].value_counts()
        zone_percent = (zone_counts / total_variants * 100).round(1)
        
        if 'unknown' in zone_counts:
            zone_counts = zone_counts.drop('unknown')
            zone_percent = zone_percent.drop('unknown')
        
        if not zone_counts.empty:
            print("\nConservation zone distribution:")
            for zone, count in zone_counts.items():
                print(f"  {zone}: {count} variants ({zone_percent[zone]}%)")
        
        # Save mapped variants
        self.mapped_variants = processed_variants
        
        # Save to file
        map_file = os.path.join(self.data_dir, 'variant_regulatory_annotations.tsv')
        processed_variants.to_csv(map_file, sep='\t', index=False)
        print(f"Saved mapped variants to {map_file}")
        
        # Store mapping statistics
        self.mapping_statistics = {
            'total_variants': total_variants,
            'mapped_to_region': mapped_to_region,
            'mapped_to_zone': mapped_to_zone,
            'region_counts': region_counts.to_dict(),
            'region_percent': region_percent.to_dict(),
            'zone_counts': zone_counts.to_dict(),
            'zone_percent': zone_percent.to_dict() if not zone_percent.empty else {}
        }
        
        # Save statistics
        stats_file = os.path.join(self.data_dir, 'variant_mapping_statistics.json')
        with open(stats_file, 'w') as f:
            # Convert statistics to JSON-serializable format
            json_stats = {}
            for category, cat_stats in self.mapping_statistics.items():
                if isinstance(cat_stats, dict):
                    json_stats[category] = {str(k): (float(v) if isinstance(v, np.number) else v)
                                          for k, v in cat_stats.items()}
                else:
                    json_stats[category] = cat_stats
            
            json.dump(json_stats, f, indent=2)
        
        print(f"Saved mapping statistics to {stats_file}")
        
        return processed_variants
    
    def _analyze_tfbs(self, variants):
        """Analyze transcription factor binding sites in variant contexts"""
        if self.genome_seq is None or not self.motifs:
            print("Reference genome or motif database not available. Skipping TFBS analysis.")
            return variants
        
        print("Analyzing transcription factor binding sites...")
        
        # Count variants with TFBS
        tfbs_count = 0
        
        for idx, variant in variants.iterrows():
            # Skip if not in a regulatory region
            if variant['Regulatory_Region'] == 'unknown':
                continue
            
            # Get scaffold and position
            scaffold = variant['Scaffold']
            position = variant['Position']
            
            # Skip if scaffold not in reference
            if scaffold not in self.genome_seq:
                continue
            
            try:
                # Extract sequence context (±20bp)
                seq_record = self.genome_seq[scaffold]
                start = max(0, position - 20)
                end = min(len(seq_record.seq), position + 20)
                
                # Skip if sequence bounds are invalid
                if start >= end or end - start < 10:
                    continue
                
                # Get reference sequence context
                ref_context = str(seq_record.seq[start:end])
                
                # Create alternative sequence with variant
                ref_allele = variant['Ref']
                alt_allele = variant['Alt']
                
                # Determine the relative position of the variant within the context
                relative_pos = position - start
                
                # Create the alternative sequence
                alt_context = ref_context[:relative_pos] + alt_allele + ref_context[relative_pos + len(ref_allele):]
                
                # Check for TF binding sites
                affected_tfs = []
                disruption_scores = {}
                
                for tf, motif in self.motifs.items():
                    ref_matches = list(re.finditer(motif, ref_context, re.IGNORECASE))
                    alt_matches = list(re.finditer(motif, alt_context, re.IGNORECASE))
                    
                    # Check if the variant creates or destroys a binding site
                    if len(ref_matches) != len(alt_matches):
                        affected_tfs.append(tf)
                        
                        # Calculate disruption score (simple version)
                        if len(ref_matches) > len(alt_matches):
                            # Site lost, negative score
                            disruption_scores[tf] = -1 * (len(ref_matches) - len(alt_matches))
                        else:
                            # Site gained, positive score
                            disruption_scores[tf] = len(alt_matches) - len(ref_matches)
                
                # If any TFs are affected, store the information
                if affected_tfs:
                    variants.at[idx, 'TFBS'] = ','.join(affected_tfs)
                    variants.at[idx, 'TFBS_Impact'] = 'affected'
                    
                    # Store the most significant disruption score
                    if disruption_scores:
                        variants.at[idx, 'Motif_Disruption_Score'] = max(disruption_scores.values(), key=abs)
                    
                    tfbs_count += 1
            
            except Exception as e:
                continue
        
        print(f"Found {tfbs_count} variants affecting transcription factor binding sites")
        
        return variants
    
    def analyze_region_distributions(self):
        """Analyze the distribution of variants across regulatory regions"""
        print("\nAnalyzing distribution of variants across regulatory regions...")
        
        if self.mapped_variants is None:
            print("No mapped variants. Run map_variants_to_regulatory_regions() first.")
            return
        
        # Distribution by regulatory region
        region_counts = self.mapped_variants['Regulatory_Region'].value_counts()
        region_percent = (region_counts / len(self.mapped_variants) * 100).round(1)
        
        # Save region counts to file
        region_file = os.path.join(self.data_dir, 'regulatory_region_distribution.tsv')
        pd.DataFrame({
            'Region': region_counts.index,
            'Count': region_counts.values,
            'Percentage': region_percent.values
        }).to_csv(region_file, sep='\t', index=False)
        
        # Plot region distribution
        plt.figure(figsize=(12, 6))
        sns.barplot(x=region_counts.index, y=region_counts.values, palette='viridis')
        plt.title('Distribution of Variants by Regulatory Region', fontsize=14)
        plt.xlabel('Regulatory Region', fontsize=12)
        plt.ylabel('Number of Variants', fontsize=12)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(os.path.join(self.plot_dir, 'regulatory_region_distribution.png'), dpi=300)
        plt.close()
        
        # Distribution by conservation zone
        zone_counts = self.mapped_variants['Conservation_Zone'].value_counts()
        zone_percent = (zone_counts / len(self.mapped_variants) * 100).round(1)
        
        # Remove 'unknown' for cleaner visualization
        if 'unknown' in zone_counts:
            unknown_count = zone_counts['unknown']
            unknown_percent = zone_percent['unknown']
            zone_counts = zone_counts.drop('unknown')
            zone_percent = zone_percent.drop('unknown')
            print(f"Note: {unknown_count} variants ({unknown_percent}%) not mapped to any conservation zone")
        
        if not zone_counts.empty:
            # Save zone counts to file
            zone_file = os.path.join(self.data_dir, 'conservation_zone_distribution.tsv')
            pd.DataFrame({
                'Zone': zone_counts.index,
                'Count': zone_counts.values,
                'Percentage': zone_percent.values
            }).to_csv(zone_file, sep='\t', index=False)
            
            # Plot zone distribution
            plt.figure(figsize=(12, 6))
            sns.barplot(x=zone_counts.index, y=zone_counts.values, palette='muted')
            plt.title('Distribution of Variants by Conservation Zone', fontsize=14)
            plt.xlabel('Conservation Zone', fontsize=12)
            plt.ylabel('Number of Variants', fontsize=12)
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(os.path.join(self.plot_dir, 'conservation_zone_distribution.png'), dpi=300)
            plt.close()
        
        # Treatment-specific distributions
        if 'Treatment' in self.mapped_variants.columns:
            # Create a crosstab of treatment vs. regulatory region
            treatment_region = pd.crosstab(
                self.mapped_variants['Treatment'],
                self.mapped_variants['Regulatory_Region'],
                normalize='index'
            ) * 100  # Convert to percentages
            
            # Save treatment-region distribution
            treatment_region.to_csv(os.path.join(self.data_dir, 'treatment_region_distribution.tsv'), sep='\t')
            
            # Plot as heatmap
            plt.figure(figsize=(14, 10))
            sns.heatmap(treatment_region, annot=True, fmt='.1f', cmap='viridis',
                       linewidths=0.5, cbar_kws={'label': 'Percentage of Variants'})
            plt.title('Regulatory Region Distribution by Treatment', fontsize=14)
            plt.xlabel('Regulatory Region', fontsize=12)
            plt.ylabel('Treatment', fontsize=12)
            plt.tight_layout()
            plt.savefig(os.path.join(self.plot_dir, 'treatment_region_heatmap.png'), dpi=300)
            plt.close()
        
        # TFBS analysis
        if 'TFBS' in self.mapped_variants.columns:
            # Filter for variants with TFBS impacts
            tfbs_variants = self.mapped_variants[self.mapped_variants['TFBS'].notna()].copy()
            
            if len(tfbs_variants) > 0:
                print(f"\nAnalyzing {len(tfbs_variants)} variants affecting transcription factor binding sites...")
                
                # Create a list of all affected TFs
                all_tfs = []
                for tfs in tfbs_variants['TFBS'].dropna():
                    all_tfs.extend(tfs.split(','))
                
                # Count TF occurrences
                tf_counts = Counter(all_tfs)
                
                # Save TF counts
                tf_file = os.path.join(self.data_dir, 'tfbs_impact_distribution.tsv')
                pd.DataFrame({
                    'Transcription_Factor': list(tf_counts.keys()),
                    'Count': list(tf_counts.values()),
                    'Percentage': [(count / len(tfbs_variants) * 100) for count in tf_counts.values()]
                }).to_csv(tf_file, sep='\t', index=False)
                
                # Plot TF distribution
                plt.figure(figsize=(12, 8))
                sns.barplot(x=list(tf_counts.keys()), y=list(tf_counts.values()), palette='deep')
                plt.title('Transcription Factors Affected by Variants', fontsize=14)
                plt.xlabel('Transcription Factor', fontsize=12)
                plt.ylabel('Number of Variants', fontsize=12)
                plt.xticks(rotation=45, ha='right')
                plt.tight_layout()
                plt.savefig(os.path.join(self.plot_dir, 'tfbs_impact_distribution.png'), dpi=300)
                plt.close()
                
                # If treatment information is available, analyze TF impacts by treatment
                if 'Treatment' in tfbs_variants.columns:
                    # Create a dictionary to store TF counts by treatment
                    treatment_tf_counts = {}
                    
                    for treatment in tfbs_variants['Treatment'].unique():
                        treatment_variants = tfbs_variants[tfbs_variants['Treatment'] == treatment]
                        
                        treatment_tfs = []
                        for tfs in treatment_variants['TFBS'].dropna():
                            treatment_tfs.extend(tfs.split(','))
                        
                        treatment_tf_counts[treatment] = Counter(treatment_tfs)
                    
                    # Create a combined dataframe for TF impacts by treatment
                    all_tfs = sorted(set(all_tfs))
                    treatments = sorted(tfbs_variants['Treatment'].unique())
                    
                    heatmap_data = []
                    for tf in all_tfs:
                        row = {'TF': tf}
                        for treatment in treatments:
                            row[treatment] = treatment_tf_counts[treatment].get(tf, 0)
                        heatmap_data.append(row)
                    
                    heatmap_df = pd.DataFrame(heatmap_data)
                    heatmap_df.set_index('TF', inplace=True)
                    
                    # Normalize by treatment totals
                    for treatment in treatments:
                        treatment_total = len(tfbs_variants[tfbs_variants['Treatment'] == treatment])
                        if treatment_total > 0:
                            heatmap_df[treatment] = heatmap_df[treatment] / treatment_total * 100
                    
                    # Save TF impacts by treatment
                    heatmap_df.to_csv(os.path.join(self.data_dir, 'tfbs_impact_by_treatment.tsv'), sep='\t')
                    
                    # Plot heatmap
                    plt.figure(figsize=(14, 10))
                    sns.heatmap(heatmap_df, annot=True, fmt='.1f', cmap='YlOrRd',
                               linewidths=0.5, cbar_kws={'label': 'Percentage of Treatment Variants'})
                    plt.title('Transcription Factor Impacts by Treatment', fontsize=14)
                    plt.xlabel('Treatment', fontsize=12)
                    plt.ylabel('Transcription Factor', fontsize=12)
                    plt.tight_layout()
                    plt.savefig(os.path.join(self.plot_dir, 'tfbs_impact_by_treatment_heatmap.png'), dpi=300)
                    plt.close()
                
                # Store TFBS variants for later use
                self.tfbs_variants = tfbs_variants
                tfbs_variants.to_csv(os.path.join(self.data_dir, 'tfbs_variants.tsv'), sep='\t', index=False)
        
        # Regulatory region variants by treatment
        if 'Treatment' in self.mapped_variants.columns:
            # Filter for variants mapped to regulatory regions
            reg_variants = self.mapped_variants[self.mapped_variants['Regulatory_Region'] != 'unknown'].copy()
            
            if len(reg_variants) > 0:
                # Store region variants for later use
                self.region_variants = reg_variants
                reg_variants.to_csv(os.path.join(self.data_dir, 'regulatory_region_variants.tsv'), sep='\t', index=False)
                
                # Create treatment-specific dataframes
                for treatment in reg_variants['Treatment'].unique():
                    treatment_variants = reg_variants[reg_variants['Treatment'] == treatment]
                    
                    # Save treatment-specific data
                    treatment_file = os.path.join(self.data_dir, f'{treatment}_regulatory_variants.tsv')
                    treatment_variants.to_csv(treatment_file, sep='\t', index=False)
                    
                    # Count by region
                    region_counts = treatment_variants['Regulatory_Region'].value_counts()
                    region_percent = (region_counts / len(treatment_variants) * 100).round(1)
                    
                    print(f"\nTop regulatory regions for {treatment}:")
                    for region, count in region_counts.nlargest(3).items():
                        print(f"  {region}: {count} variants ({region_percent[region]}%)")
    
    def analyze_distance_relationships(self):
        """Analyze distance relationships between variants and features"""
        print("\nAnalyzing distance relationships between variants and features...")
        
        if self.mapped_variants is None:
            print("No mapped variants. Run map_variants_to_regulatory_regions() first.")
            return
        
        # Check if Distance column is available
        if 'Distance' not in self.mapped_variants.columns:
            print("Distance information not available for analysis.")
            return
        
        # Filter for variants with distance information
        variants_with_distance = self.mapped_variants.dropna(subset=['Distance']).copy()
        
        if len(variants_with_distance) == 0:
            print("No variants with distance information available.")
            return
        
        print(f"Analyzing distance relationships for {len(variants_with_distance)} variants...")
        
        # Create distance bins for visualization
        bins = [0, 50, 100, 200, 500, 1000, 2000, 5000, 10000, float('inf')]
        labels = ['0-50', '51-100', '101-200', '201-500', '501-1000', 
                 '1001-2000', '2001-5000', '5001-10000', '>10000']
        
        variants_with_distance['Distance_Bin'] = pd.cut(
            variants_with_distance['Distance'].abs(), bins=bins, labels=labels, right=False
        )
        
        # Count variants by distance bin
        distance_counts = variants_with_distance['Distance_Bin'].value_counts().sort_index()
        distance_percent = (distance_counts / len(variants_with_distance) * 100).round(1)
        
        # Save distance distribution
        distance_file = os.path.join(self.data_dir, 'distance_distribution.tsv')
        pd.DataFrame({
            'Distance_Bin': distance_counts.index,
            'Count': distance_counts.values,
            'Percentage': distance_percent.values
        }).to_csv(distance_file, sep='\t', index=False)
        
        # Plot distance distribution
        plt.figure(figsize=(12, 6))
        sns.barplot(x=distance_counts.index, y=distance_counts.values, palette='Blues_d')
        plt.title('Distribution of Variants by Distance from Gene', fontsize=14)
        plt.xlabel('Distance from Gene (bp)', fontsize=12)
        plt.ylabel('Number of Variants', fontsize=12)
        plt.tight_layout()
        plt.savefig(os.path.join(self.plot_dir, 'distance_distribution.png'), dpi=300)
        plt.close()
        
        # Relationship between distance and regulatory region
        region_distance_data = []
        
        for region in variants_with_distance['Regulatory_Region'].unique():
            if region == 'unknown':
                continue
                
            region_variants = variants_with_distance[variants_with_distance['Regulatory_Region'] == region]
            
            if len(region_variants) == 0:
                continue
                
            region_distance_data.append({
                'Region': region,
                'Mean_Distance': region_variants['Distance'].abs().mean(),
                'Median_Distance': region_variants['Distance'].abs().median(),
                'Min_Distance': region_variants['Distance'].abs().min(),
                'Max_Distance': region_variants['Distance'].abs().max(),
                'Count': len(region_variants)
            })
        
        # Create DataFrame and save
        if region_distance_data:
            region_distance_df = pd.DataFrame(region_distance_data)
            region_distance_df.to_csv(os.path.join(self.data_dir, 'region_distance_statistics.tsv'), sep='\t', index=False)
            
            # Plot median distances by region
            plt.figure(figsize=(12, 6))
            sns.barplot(x='Region', y='Median_Distance', data=region_distance_df, palette='viridis')
            plt.title('Median Distance by Regulatory Region', fontsize=14)
            plt.xlabel('Regulatory Region', fontsize=12)
            plt.ylabel('Median Distance from Gene (bp)', fontsize=12)
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(os.path.join(self.plot_dir, 'region_median_distance.png'), dpi=300)
            plt.close()
        
        # Treatment-specific distance analysis
        if 'Treatment' in variants_with_distance.columns:
            # Create distance bins for each treatment
            treatment_distance_data = []
            
            for treatment in variants_with_distance['Treatment'].unique():
                treatment_variants = variants_with_distance[variants_with_distance['Treatment'] == treatment]
                
                for bin_label in labels:
                    bin_variants = treatment_variants[treatment_variants['Distance_Bin'] == bin_label]
                    
                    if len(bin_variants) > 0:
                        treatment_distance_data.append({
                            'Treatment': treatment,
                            'Distance_Bin': bin_label,
                            'Count': len(bin_variants),
                            'Percentage': len(bin_variants) / len(treatment_variants) * 100
                        })
            
            # Create DataFrame and save
            if treatment_distance_data:
                treatment_distance_df = pd.DataFrame(treatment_distance_data)
                treatment_distance_df.to_csv(os.path.join(self.data_dir, 'treatment_distance_distribution.tsv'), 
                                            sep='\t', index=False)
                
                # Create pivot table for heatmap
                pivot_df = treatment_distance_df.pivot(index='Treatment', columns='Distance_Bin', values='Percentage')
                pivot_df.fillna(0, inplace=True)
                
                # Plot heatmap
                plt.figure(figsize=(14, 8))
                sns.heatmap(pivot_df, annot=True, fmt='.1f', cmap='viridis',
                           linewidths=0.5, cbar_kws={'label': 'Percentage of Treatment Variants'})
                plt.title('Distance Distribution by Treatment', fontsize=14)
                plt.xlabel('Distance from Gene (bp)', fontsize=12)
                plt.ylabel('Treatment', fontsize=12)
                plt.tight_layout()
                plt.savefig(os.path.join(self.plot_dir, 'treatment_distance_heatmap.png'), dpi=300)
                plt.close()
    
    def generate_mapping_report(self):
        """Generate a comprehensive report of the variant regulatory mapping"""
        print("\nGenerating comprehensive variant regulatory mapping report...")
        
        if self.mapped_variants is None:
            print("No mapped variants. Run map_variants_to_regulatory_regions() first.")
            return
        
        # Create report content
        report = ["# Advanced Variant-to-Regulatory Feature Mapping Report\n"]
        
        # 1. Overview section
        report.append("## Overview\n")
        report.append(f"Total variants analyzed: {len(self.mapped_variants)}")
        
        mapped_regions = len(self.mapped_variants[self.mapped_variants['Regulatory_Region'] != 'unknown'])
        mapped_zones = len(self.mapped_variants[self.mapped_variants['Conservation_Zone'] != 'unknown'])
        
        report.append(f"Variants mapped to regulatory regions: {mapped_regions} ({mapped_regions/len(self.mapped_variants)*100:.1f}%)")
        report.append(f"Variants mapped to conservation zones: {mapped_zones} ({mapped_zones/len(self.mapped_variants)*100:.1f}%)")
        
        if 'TFBS' in self.mapped_variants.columns:
            tfbs_variants = len(self.mapped_variants[self.mapped_variants['TFBS'].notna()])
            report.append(f"Variants affecting transcription factor binding sites: {tfbs_variants} ({tfbs_variants/len(self.mapped_variants)*100:.1f}%)")
        
        report.append("")
        
        # 2. Regulatory Region Distribution
        report.append("## Regulatory Region Distribution\n")
        
        region_counts = self.mapped_variants['Regulatory_Region'].value_counts()
        region_percent = (region_counts / len(self.mapped_variants) * 100).round(1)
        
        report.append("| Region | Count | Percentage |")
        report.append("|--------|-------|------------|")
        
        for region, count in region_counts.items():
            if region == 'unknown':
                continue
            report.append(f"| {region} | {count} | {region_percent[region]}% |")
        
        report.append("")
        
        # 3. Conservation Zone Distribution
        zone_counts = self.mapped_variants['Conservation_Zone'].value_counts()
        zone_percent = (zone_counts / len(self.mapped_variants) * 100).round(1)
        
        if 'unknown' in zone_counts and zone_counts.drop('unknown').empty:
            report.append("## Conservation Zone Distribution\n")
            report.append("No variants mapped to specific conservation zones.")
        elif not zone_counts.drop('unknown', errors='ignore').empty:
            report.append("## Conservation Zone Distribution\n")
            
            report.append("| Zone | Count | Percentage |")
            report.append("|------|-------|------------|")
            
            for zone, count in zone_counts.items():
                if zone == 'unknown':
                    continue
                report.append(f"| {zone} | {count} | {zone_percent[zone]}% |")
            
            report.append("")
        
        # 4. Distance Relationships
        if 'Distance' in self.mapped_variants.columns:
            variants_with_distance = self.mapped_variants.dropna(subset=['Distance'])
            
            if len(variants_with_distance) > 0:
                report.append("## Distance Relationships\n")
                
                # Distance statistics
                mean_distance = variants_with_distance['Distance'].abs().mean()
                median_distance = variants_with_distance['Distance'].abs().median()
                min_distance = variants_with_distance['Distance'].abs().min()
                max_distance = variants_with_distance['Distance'].abs().max()
                
                report.append(f"Mean distance from gene: {mean_distance:.1f} bp")
                report.append(f"Median distance from gene: {median_distance:.1f} bp")
                report.append(f"Minimum distance from gene: {min_distance:.1f} bp")
                report.append(f"Maximum distance from gene: {max_distance:.1f} bp")
                report.append("")
                
                # Distance distribution by bin
                if 'Distance_Bin' in variants_with_distance.columns:
                    distance_counts = variants_with_distance['Distance_Bin'].value_counts().sort_index()
                    distance_percent = (distance_counts / len(variants_with_distance) * 100).round(1)
                    
                    report.append("### Distance Distribution\n")
                    report.append("| Distance Range (bp) | Count | Percentage |")
                    report.append("|---------------------|-------|------------|")
                    
                    for bin_label, count in distance_counts.items():
                        report.append(f"| {bin_label} | {count} | {distance_percent[bin_label]}% |")
                    
                    report.append("")
        
        # 5. Transcription Factor Binding Sites
        if 'TFBS' in self.mapped_variants.columns:
            tfbs_variants = self.mapped_variants[self.mapped_variants['TFBS'].notna()]
            
            if len(tfbs_variants) > 0:
                report.append("## Transcription Factor Binding Sites\n")
                
                # Extract all affected TFs
                all_tfs = []
                for tfs in tfbs_variants['TFBS'].dropna():
                    all_tfs.extend(tfs.split(','))
                
                # Count TF occurrences
                tf_counts = Counter(all_tfs)
                
                report.append("| Transcription Factor | Count | Percentage |")
                report.append("|----------------------|-------|------------|")
                
                for tf, count in tf_counts.most_common(10):
                    percent = count / len(tfbs_variants) * 100
                    report.append(f"| {tf} | {count} | {percent:.1f}% |")
                
                report.append("")
        
        # 6. Treatment-Specific Patterns
        if 'Treatment' in self.mapped_variants.columns:
            report.append("## Treatment-Specific Patterns\n")
            
            # Treatment counts
            treatment_counts = self.mapped_variants['Treatment'].value_counts()
            
            for treatment, count in treatment_counts.items():
                treatment_variants = self.mapped_variants[self.mapped_variants['Treatment'] == treatment]
                
                report.append(f"### {treatment}\n")
                report.append(f"Total variants: {count}\n")
                
                # Region distribution
                treatment_regions = treatment_variants['Regulatory_Region'].value_counts()
                treatment_region_percent = (treatment_regions / count * 100).round(1)
                
                report.append("#### Top Regulatory Regions\n")
                
                for region, region_count in treatment_regions.nlargest(5).items():
                    if region == 'unknown':
                        continue
                    report.append(f"- {region}: {region_count} variants ({treatment_region_percent[region]}%)")
                
                report.append("")
                
                # TFBS impacts
                if 'TFBS' in treatment_variants.columns:
                    tfbs_count = len(treatment_variants[treatment_variants['TFBS'].notna()])
                    if tfbs_count > 0:
                        report.append(f"#### Transcription Factor Binding Sites\n")
                        report.append(f"Variants affecting TFBS: {tfbs_count} ({tfbs_count/count*100:.1f}%)")
                        
                        # Extract TFs for this treatment
                        treatment_tfs = []
                        for tfs in treatment_variants['TFBS'].dropna():
                            treatment_tfs.extend(tfs.split(','))
                        
                        # Count and report most common TFs
                        if treatment_tfs:
                            tf_counts = Counter(treatment_tfs)
                            
                            report.append("\nMost affected transcription factors:")
                            for tf, tf_count in tf_counts.most_common(3):
                                percent = tf_count / tfbs_count * 100
                                report.append(f"- {tf}: {tf_count} variants ({percent:.1f}%)")
                        
                        report.append("")
        
        # 7. Biological Interpretations
        report.append("## Biological Interpretations\n")
        
        # Add interpretations based on the mapping results
        if 'Regulatory_Region' in self.mapped_variants.columns:
            core_promo_count = len(self.mapped_variants[self.mapped_variants['Regulatory_Region'] == 'core_promoter'])
            uas_count = len(self.mapped_variants[self.mapped_variants['Regulatory_Region'].isin(['UAS_proximal', 'UAS_distal'])])
            
            if core_promo_count > 0 or uas_count > 0:
                if core_promo_count > uas_count:
                    report.append("- **Core Promoter Focus**: A higher proportion of variants occur in core promoter regions, " +
                                 "suggesting adaptation through changes in basal transcription machinery recruitment and " +
                                 "transcription initiation.")
                else:
                    report.append("- **Upstream Activating Sequence Focus**: A higher proportion of variants occur in UAS regions, " +
                                 "suggesting adaptation through changes in binding of specific transcription factors and " +
                                 "condition-responsive gene regulation.")
        
        if 'TFBS' in self.mapped_variants.columns and len(self.mapped_variants[self.mapped_variants['TFBS'].notna()]) > 0:
            report.append("- **Transcription Factor Binding Site Alterations**: Variants affecting TF binding sites suggest " +
                         "adaptation through modified transcription factor recruitment and altered gene expression patterns.")
        
        if 'Treatment' in self.mapped_variants.columns:
            temp_treatments = ['WT-37', 'CAS']
            oxy_treatments = ['WTA', 'STC']
            
            temp_variants = self.mapped_variants[self.mapped_variants['Treatment'].isin(temp_treatments)]
            oxy_variants = self.mapped_variants[self.mapped_variants['Treatment'].isin(oxy_treatments)]
            
            if len(temp_variants) > 0 and len(oxy_variants) > 0:
                # Compare region distributions
                temp_regions = temp_variants['Regulatory_Region'].value_counts(normalize=True)
                oxy_regions = oxy_variants['Regulatory_Region'].value_counts(normalize=True)
                
                # Compare TFBS impacts
                if 'TFBS' in self.mapped_variants.columns:
                    temp_tfbs = len(temp_variants[temp_variants['TFBS'].notna()]) / len(temp_variants) if len(temp_variants) > 0 else 0
                    oxy_tfbs = len(oxy_variants[oxy_variants['TFBS'].notna()]) / len(oxy_variants) if len(oxy_variants) > 0 else 0
                    
                    if temp_tfbs > oxy_tfbs:
                        report.append("- **Temperature Adaptation TF Sensitivity**: Temperature adaptation shows a higher proportion " +
                                     "of variants affecting transcription factor binding sites, suggesting more extensive " +
                                     "rewiring of transcriptional networks in response to thermal stress.")
                    elif oxy_tfbs > temp_tfbs:
                        report.append("- **Oxygen Adaptation TF Sensitivity**: Low oxygen adaptation shows a higher proportion " +
                                     "of variants affecting transcription factor binding sites, suggesting more extensive " +
                                     "rewiring of transcriptional networks in response to oxygen limitation.")
        
        report.append("- **Hierarchical Regulatory Architecture**: The distribution of variants across different regulatory " +
                     "regions and conservation zones reveals a hierarchical structure that balances essential function " +
                     "preservation with adaptive flexibility.")
        
        # Write report to file
        report_path = os.path.join(self.output_dir, 'variant_regulatory_mapping_report.md')
        with open(report_path, 'w') as f:
            f.write('\n'.join(report))
        
        print(f"Generated comprehensive report: {report_path}")
        
        return report
    
    def run_analysis(self):
        """Run the complete variant regulatory mapping pipeline"""
        print("\n=== Starting Advanced Variant-to-Regulatory Feature Mapping ===\n")
        
        # Step 1: Preprocess variants
        processed_variants = self.preprocess_variants()
        
        # Step 2: Map variants to regulatory regions
        self.mapped_variants = self.map_variants_to_regulatory_regions(processed_variants)
        
        # Step 3: Analyze region distributions
        self.analyze_region_distributions()
        
        # Step 4: Analyze distance relationships
        self.analyze_distance_relationships()
        
        # Step 5: Generate mapping report
        self.generate_mapping_report()
        
        print("\n=== Advanced Variant-to-Regulatory Feature Mapping Complete ===")
        
        # Return mapped variants for potential use by other scripts
        return self.mapped_variants


def main():
    """Main function to parse arguments and run the analysis"""
    parser = argparse.ArgumentParser(description='Advanced Variant-to-Regulatory Feature Mapping System')
    
    # Required arguments
    parser.add_argument('--variants', 
                       default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/gene_variants_expanded/all_gene_variants.tsv',
                       help='Variants TSV file')
    parser.add_argument('--regulatory-map', 
                      default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/new_regulatory_analysis/data/gene_regulatory_map.json',
                      help='Gene regulatory map JSON file')
    parser.add_argument('--output-dir', 
                       default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/new_regulatory_analysis',
                       help='Output directory for results')
    
    # Optional arguments
    parser.add_argument('--reference-genome', 
                       default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/reference/w303_chromosomal.fasta',
                       help='Reference genome FASTA file')
    parser.add_argument('--regulatory-features', 
                      default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/scripts/new_regulatory_analysis/data/regulatory_features.json',
                      help='Regulatory features JSON file')
    parser.add_argument('--motif-database', 
                       default=None,
                       help='Motif database file (JSON)')
    
    args = parser.parse_args()
    
    # Initialize the mapper
    mapper = VariantRegulatoryMapper(
        args.variants,
        args.regulatory_map,
        args.reference_genome,
        args.regulatory_features,
        args.output_dir,
        args.motif_database
    )
    
    # Run the analysis
    mapper.run_analysis()


if __name__ == "__main__":
    main()