#!/usr/bin/env python3
# /Users/zakiralibhai/Documents/GitHub/Yeast_MSA/scripts/regulatory_analysis/regulatory_motif_analysis.py

"""
Identifies and analyzes transcription factor binding sites in regulatory regions
Part of the Comprehensive Regulatory Region Analysis in the Yeast MSA project.

This script:
1. Maps variants to precise regulatory regions based on yeast-specific promoter architecture
2. Identifies potential transcription factor binding sites that may be affected
3. Analyzes distribution of regulatory element variants by treatment
4. Generates visualizations of regulatory patterns across treatments

Yeast promoter architecture features:
- Core promoter (0-150bp upstream): Contains essential transcription elements
  - TATA box (40-120bp upstream): ~20% of yeast genes have canonical TATA elements
  - Initiator element (Inr): Located at the transcription start site
  - Transcription start site: Typically 40-120bp downstream of TATA box
- Upstream Activating Sequences (UAS): Enhancer-like elements often 100-700bp upstream
- Upstream Repressing Sequences (URS): Silencer elements that inhibit transcription
- 5' UTRs: Typically short in yeast (15-65bp)
- 3' UTRs and terminator regions: Post-transcriptional regulatory elements
"""

import os
import sys
import argparse
import json
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
import warnings

# Suppress FutureWarnings from pandas to clean output
warnings.simplefilter(action='ignore', category=FutureWarning)

class RegulatoryMotifAnalyzer:
    """Analyzes regulatory motifs and transcription factor binding sites in yeast"""
    
    def __init__(self, output_dir, gene_mapping_file, variants_file, erg_genes_file, 
                 reference_genome=None, motif_database=None):
        """Initialize with paths and configurations"""
        # Setup directories
        self.output_dir = self._ensure_dir(output_dir)
        self.data_dir = self._ensure_dir(os.path.join(output_dir, 'data'))
        self.plot_dir = self._ensure_dir(os.path.join(output_dir, 'plots'))
        self.motif_dir = self._ensure_dir(os.path.join(output_dir, 'motifs'))
        
        # Load required data
        print("Loading input files...")
        self.gene_mapping = self._load_file(gene_mapping_file, 'gene mapping')
        self.variants = self._load_file(variants_file, 'variants')
        self.erg_genes = self._load_file(erg_genes_file, 'ERG genes')
        
        # Optional reference genome for sequence context
        self.reference_genome = reference_genome
        self.genome_seq = None
        if reference_genome:
            try:
                self.genome_seq = self._load_reference_genome(reference_genome)
                print(f"Loaded reference genome: {reference_genome}")
            except Exception as e:
                print(f"Warning: Failed to load reference genome: {e}")
        
        # Motif database for TF binding sites
        self.motif_database = motif_database
        self.motifs = {}
        if motif_database:
            try:
                self.motifs = self._load_motif_database(motif_database)
                print(f"Loaded motif database with {len(self.motifs)} TF binding motifs")
            except Exception as e:
                print(f"Warning: Failed to load motif database: {e}")
        
        # Define yeast-specific regulatory regions with precise boundaries
        self._define_regulatory_regions()
        
        # Treatment groups
        self.treatments = ['WT-37', 'WTA', 'CAS', 'STC']
        self.treatment_groups = {
            'Temperature': ['WT-37', 'CAS'],
            'Low Oxygen': ['WTA', 'STC'],
            'Gene Modified': ['CAS', 'STC'],
            'Non-Modified': ['WT-37', 'WTA']
        }
        
        # For analysis results
        self.reg_variants = None
        self.tf_site_variants = None
    
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
            if description == 'gene mapping' or description == 'variants':
                print("Critical file missing. Exiting.")
                sys.exit(1)
            return pd.DataFrame()  # Return empty dataframe for non-critical files
    
    def _load_reference_genome(self, fasta_file):
        """Load reference genome from FASTA file"""
        return SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    
    def _load_motif_database(self, motif_file):
        """Load transcription factor binding motifs
        
        Expected format: JSON with TF names as keys and consensus sequences as values
        Example: {"GCN4": "RTGACTCAY", "HAP1": "CGGNNNTANCGG"}
        """
        if motif_file.endswith('.json'):
            with open(motif_file, 'r') as f:
                return json.load(f)
        else:
            # If not available, use a minimal set of known yeast motifs
            print("Using built-in set of common yeast motifs")
            return {
                # General transcription factors
                "TATA-box": "TATA[AT]A[AT]",       # TATA box consensus
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
                "Rox1": "YYNATTGTTY"               # Repressor of hypoxic genes
            }
    
    def _define_regulatory_regions(self):
        """Define detailed yeast-specific regulatory regions"""
        # Upstream regulatory regions (distance in bp upstream of ATG)
        self.upstream_regions = {
            'core_promoter': (0, 150),          # Core promoter contains essential elements
            'TATA_region': (40, 120),           # TATA box typically found here
            'tss_proximal': (0, 40),            # Region immediately upstream of TSS
            'UAS_proximal': (150, 500),         # Upstream activating sequences - proximal
            'UAS_distal': (500, 1500),          # Upstream activating sequences - distal
            'far_upstream': (1500, 10000)       # Distant regulatory elements
        }
        
        # Downstream and UTR regions (relative to coding sequence)
        self.downstream_regions = {
            'five_prime_UTR': (-60, 0),         # 5' UTR (negative values = before CDS start)
            'terminator': (0, 250),             # Transcription termination region
            'three_prime_UTR': (0, 120),        # 3' UTR region
            'downstream_element': (250, 1000)   # Other downstream control elements
        }
        
        # Special regulatory elements for yeast
        self.special_elements = {
            'TATA_box': {'consensus': 'TATA[AT]A[AT]', 'region': (40, 120)},
            'initiator': {'consensus': 'YYANWYY', 'region': (0, 40)},            # Initiator element
            'STRE': {'consensus': 'AGGGG', 'region': (100, 700)},                # Stress response
            'URS': {'consensus': '[ACG]CCCC[ACT]', 'region': (150, 800)},        # Upstream repressing
            'sterol_regulatory': {'consensus': 'TCGTATA', 'region': (100, 700)}  # Sterol regulation
        }
        
        # Known yeast promoter types
        self.promoter_classes = {
            'TATA_containing': 'Contains canonical TATA box (20% of yeast genes)',
            'TATA_less': 'No TATA box, often housekeeping genes (80% of yeast genes)',
            'constitutive': 'Consistently expressed genes like ribosomal proteins',
            'regulated': 'Genes regulated by environmental conditions',
            'stress_responsive': 'Induced under stress conditions'
        }
    
    def prepare_data(self):
        """Prepare and preprocess data for analysis"""
        print("\nPreparing data for regulatory motif analysis...")
        
        # Display sample variants to understand data
        print("\nSample variants:")
        if len(self.variants) > 0:
            sample = self.variants.head(5)
            display_cols = ['Gene_ID', 'Effect', 'Distance', 'Ref', 'Alt', 'Treatment']
            display_cols = [col for col in display_cols if col in sample.columns]
            print(sample[display_cols].to_string())
        
        # Check distance statistics if available
        if 'Distance' in self.variants.columns:
            distances = self.variants['Distance'].dropna()
            if len(distances) > 0:
                print(f"\nDistance statistics:")
                print(f"  Range: {distances.min()} to {distances.max()} bp")
                print(f"  Mean: {distances.mean():.1f} bp")
                print(f"  Median: {distances.median()} bp")
        
        # Filter for regulatory variants only
        if 'Effect' in self.variants.columns:
            # Common regulatory effect annotations
            regulatory_effects = [
                'upstream_gene_variant', 
                'downstream_gene_variant',
                'intergenic_region', 
                '5_prime_UTR_variant', 
                '3_prime_UTR_variant'
            ]
            
            # Filter for regulatory variants
            mask = self.variants['Effect'].str.contains(
                '|'.join(regulatory_effects), 
                case=False, 
                na=False
            )
            reg_variants = self.variants[mask].copy()
            print(f"\nFiltered {len(reg_variants)} regulatory variants from {len(self.variants)} total variants")
        else:
            # No Effect column - try to use Distance
            reg_variants = self.variants.copy()
        
        # Ensure distances are numeric
        if 'Distance' in reg_variants.columns:
            reg_variants['Distance'] = pd.to_numeric(reg_variants['Distance'], errors='coerce')
        
        # Store preprocessed data
        self.reg_variants = reg_variants
        return reg_variants
    
    def classify_regulatory_regions(self):
        """Classify variants into specific yeast regulatory regions"""
        print("\nClassifying variants into detailed yeast regulatory regions...")
        
        if self.reg_variants is None or len(self.reg_variants) == 0:
            print("No variants to classify. Run prepare_data() first.")
            return
        
        # Check essential columns
        required_cols = ['Effect', 'Distance']
        missing_cols = [col for col in required_cols if col not in self.reg_variants.columns]
        if missing_cols:
            print(f"Warning: Missing required columns: {missing_cols}")
            print("Classification will be limited.")
        
        # Add classification columns
        self.reg_variants['Regulatory_Region'] = 'unknown'
        self.reg_variants['Element_Type'] = 'unknown'
        
        # Extract effect and distance for classification
        if 'Effect' in self.reg_variants.columns and 'Distance' in self.reg_variants.columns:
            # Process each variant
            for idx, row in self.reg_variants.iterrows():
                effect = row['Effect']
                distance = row['Distance']
                
                # Skip if distance is not available
                if pd.isna(distance):
                    continue
                
                # Classify based on effect type and distance
                if 'upstream' in effect.lower():
                    # Upstream variants - classify into promoter regions
                    for region, (min_dist, max_dist) in self.upstream_regions.items():
                        if min_dist <= distance < max_dist:
                            self.reg_variants.at[idx, 'Regulatory_Region'] = region
                            break
                    
                    # Check for special elements based on region
                    region = self.reg_variants.at[idx, 'Regulatory_Region']
                    for element, props in self.special_elements.items():
                        element_region = props['region']
                        if element_region[0] <= distance < element_region[1]:
                            self.reg_variants.at[idx, 'Element_Type'] = element
                            break
                
                elif 'downstream' in effect.lower():
                    # Downstream variants - classify into terminator/downstream regions
                    for region, (min_dist, max_dist) in self.downstream_regions.items():
                        if min_dist <= distance < max_dist:
                            self.reg_variants.at[idx, 'Regulatory_Region'] = region
                            break
                
                elif '5_prime_utr' in effect.lower():
                    self.reg_variants.at[idx, 'Regulatory_Region'] = 'five_prime_UTR'
                
                elif '3_prime_utr' in effect.lower():
                    self.reg_variants.at[idx, 'Regulatory_Region'] = 'three_prime_UTR'
        
        # Summarize classification results
        region_counts = self.reg_variants['Regulatory_Region'].value_counts()
        print("\nRegulatory region distribution:")
        for region, count in region_counts.items():
            print(f"  {region}: {count} variants ({count/len(self.reg_variants)*100:.1f}%)")
        
        element_counts = self.reg_variants['Element_Type'].value_counts()
        if 'unknown' in element_counts:
            element_counts = element_counts.drop('unknown')
        
        if not element_counts.empty:
            print("\nRegulatory element type distribution:")
            for element, count in element_counts.items():
                print(f"  {element}: {count} variants")
        
        return self.reg_variants
    
    def identify_tf_binding_sites(self):
        """Identify potential transcription factor binding sites affected by variants"""
        print("\nIdentifying potential transcription factor binding sites...")
        
        if self.reg_variants is None or len(self.reg_variants) == 0:
            print("No variants to analyze. Run classify_regulatory_regions() first.")
            return
        
        # Check if reference genome is available
        if self.genome_seq is None:
            print("Reference genome not available. Cannot identify binding sites.")
            return
        
        # Add column for TF binding sites
        self.reg_variants['TF_Binding_Sites'] = None
        self.reg_variants['TF_Impact'] = None
        
        # Track affected binding sites
        tf_impacts = defaultdict(int)
        variant_contexts = {}
        
        # Process variants with sequence context
        variants_with_context = 0
        for idx, row in self.reg_variants.iterrows():
            # Skip if not in a relevant region
            if row['Regulatory_Region'] == 'unknown':
                continue
                
            # We need scaffold, position, and reference/alt alleles
            if not all(col in row for col in ['Scaffold', 'Position', 'Ref', 'Alt']):
                continue
            
            # Get sequence context (±10bp)
            scaffold = row['Scaffold']
            position = row['Position']
            
            try:
                # Extract sequence context (20bp around variant)
                if scaffold in self.genome_seq:
                    seq_record = self.genome_seq[scaffold]
                    # Adjust for 0-based indexing and sequence boundaries
                    start = max(0, position - 10)
                    end = min(len(seq_record.seq), position + 10)
                    context_seq = str(seq_record.seq[start:end])
                    
                    # Store sequence context for this variant
                    variant_contexts[idx] = context_seq
                    variants_with_context += 1
                    
                    # Check for TF binding motifs in this sequence
                    affected_tfs = []
                    for tf, motif in self.motifs.items():
                        # Look for pattern in the sequence context
                        if re.search(motif, context_seq, re.IGNORECASE):
                            affected_tfs.append(tf)
                            tf_impacts[tf] += 1
                    
                    # Store affected TFs for this variant
                    if affected_tfs:
                        self.reg_variants.at[idx, 'TF_Binding_Sites'] = ','.join(affected_tfs)
                        self.reg_variants.at[idx, 'TF_Impact'] = 'affected'
            
            except Exception as e:
                print(f"Error processing variant at {scaffold}:{position}: {e}")
        
        # Summarize TF binding site results
        if variants_with_context > 0:
            print(f"\nAnalyzed sequence context for {variants_with_context} variants")
            
            # Filter variants affecting TF binding sites
            tf_site_variants = self.reg_variants[self.reg_variants['TF_Impact'] == 'affected'].copy()
            self.tf_site_variants = tf_site_variants
            
            if not tf_site_variants.empty:
                print(f"Found {len(tf_site_variants)} variants potentially affecting TF binding sites")
                
                # Summarize affected TF types
                print("\nMost frequently affected transcription factors:")
                for tf, count in sorted(tf_impacts.items(), key=lambda x: x[1], reverse=True)[:5]:
                    print(f"  {tf}: {count} variants")
                
                # Save TF binding site variants
                tf_site_variants.to_csv(os.path.join(self.data_dir, 'tf_binding_site_variants.tsv'), 
                                       sep='\t', index=False)
            else:
                print("No variants affecting TF binding sites were identified")
        
        return tf_site_variants if 'tf_site_variants' in locals() else None
    
    def analyze_erg_regulatory_elements(self):
        """Analyze regulatory elements specifically around ERG genes"""
        print("\nAnalyzing regulatory elements around ergosterol pathway genes...")
        
        if self.reg_variants is None or len(self.reg_variants) == 0:
            print("No variants to analyze. Run classify_regulatory_regions() first.")
            return
        
        # Get ergosterol gene IDs
        erg_gene_ids = set(self.erg_genes['w303_gene_id'].values)
        
        # Filter for variants near ERG genes
        erg_variants = self.reg_variants[self.reg_variants['Gene_ID'].isin(erg_gene_ids)].copy()
        
        if len(erg_variants) == 0:
            print("No variants found near ergosterol pathway genes")
            return
        
        print(f"Found {len(erg_variants)} variants in regulatory regions of ergosterol genes")
        
        # Analyze distribution of regulatory regions by treatment
        if 'Treatment' in erg_variants.columns and 'Regulatory_Region' in erg_variants.columns:
            # Create region distributions for each treatment
            region_by_treatment = pd.crosstab(
                erg_variants['Treatment'], 
                erg_variants['Regulatory_Region'],
                normalize='index'
            ) * 100  # Convert to percentages
            
            # Save the data
            erg_variants.to_csv(os.path.join(self.data_dir, 'erg_regulatory_variants.tsv'), 
                               sep='\t', index=False)
            region_by_treatment.to_csv(os.path.join(self.data_dir, 'erg_region_by_treatment.tsv'), 
                                     sep='\t')
            
            # Plot distribution
            plt.figure(figsize=(12, 8))
            region_by_treatment.plot(kind='bar', stacked=True, colormap='viridis')
            plt.title('Distribution of Regulatory Regions around ERG Genes', fontsize=14)
            plt.xlabel('Treatment', fontsize=12)
            plt.ylabel('Percentage of Variants', fontsize=12)
            plt.xticks(rotation=45)
            plt.legend(title='Regulatory Region', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.savefig(os.path.join(self.plot_dir, 'erg_regulatory_distribution.png'), dpi=300)
            plt.close()
            
            # Create heatmap for better visualization
            plt.figure(figsize=(12, 8))
            sns.heatmap(region_by_treatment, annot=True, fmt='.1f', cmap='viridis',
                       linewidths=0.5, cbar_kws={'label': 'Percentage of Variants'})
            plt.title('Heatmap of Regulatory Elements around ERG Genes', fontsize=14)
            plt.xlabel('Regulatory Region', fontsize=12)
            plt.ylabel('Treatment', fontsize=12)
            plt.tight_layout()
            plt.savefig(os.path.join(self.plot_dir, 'erg_regulatory_heatmap.png'), dpi=300)
            plt.close()
        
        # Check for TF binding sites near ERG genes
        if 'TF_Binding_Sites' in erg_variants.columns:
            # Filter for variants affecting TF sites
            tf_variants = erg_variants[erg_variants['TF_Binding_Sites'].notna()]
            
            if len(tf_variants) > 0:
                print(f"Found {len(tf_variants)} variants affecting TF binding sites near ERG genes")
                
                # Extract all affected TFs and count occurrences
                all_tfs = []
                for tfs in tf_variants['TF_Binding_Sites'].dropna():
                    if isinstance(tfs, str):
                        all_tfs.extend(tfs.split(','))
                
                tf_counts = pd.Series(all_tfs).value_counts()
                
                # Save TF counts
                tf_counts.to_frame('Count').to_csv(os.path.join(self.data_dir, 'erg_tf_counts.tsv'), 
                                                sep='\t')
                
                # Plot TF distribution
                plt.figure(figsize=(10, 6))
                tf_counts.plot(kind='bar', color='darkblue')
                plt.title('Transcription Factors Affected Near ERG Genes', fontsize=14)
                plt.xlabel('Transcription Factor', fontsize=12)
                plt.ylabel('Number of Variants', fontsize=12)
                plt.xticks(rotation=45, ha='right')
                plt.tight_layout()
                plt.savefig(os.path.join(self.plot_dir, 'erg_tf_distribution.png'), dpi=300)
                plt.close()
            else:
                print("No TF binding sites affected near ERG genes")
        
        return erg_variants
    
    def analyze_treatment_differences(self):
        """Compare regulatory patterns between different treatments"""
        print("\nComparing regulatory patterns between treatments...")
        
        if self.reg_variants is None or len(self.reg_variants) == 0:
            print("No variants to analyze. Run classify_regulatory_regions() first.")
            return
        
        # Check if we have Treatment information
        if 'Treatment' not in self.reg_variants.columns:
            print("Treatment information not available for comparison")
            return
        
        # Count variants by treatment
        treatment_counts = self.reg_variants['Treatment'].value_counts()
        print("\nVariants by treatment:")
        for treatment, count in treatment_counts.items():
            print(f"  {treatment}: {count} variants")
        
        # Compare regulatory region distributions across treatments
        if 'Regulatory_Region' in self.reg_variants.columns:
            # Cross-tabulate regulatory regions by treatment
            region_by_treatment = pd.crosstab(
                self.reg_variants['Treatment'], 
                self.reg_variants['Regulatory_Region'],
                normalize='index'
            ) * 100  # Convert to percentages
            
            # Plot distribution
            plt.figure(figsize=(12, 8))
            region_by_treatment.plot(kind='bar', stacked=True, colormap='viridis')
            plt.title('Regulatory Region Distribution by Treatment', fontsize=14)
            plt.xlabel('Treatment', fontsize=12)
            plt.ylabel('Percentage of Variants', fontsize=12)
            plt.xticks(rotation=45)
            plt.legend(title='Regulatory Region', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.savefig(os.path.join(self.plot_dir, 'regulatory_region_by_treatment.png'), dpi=300)
            plt.close()
            
            # Create heatmap for clearer visualization
            plt.figure(figsize=(12, 8))
            sns.heatmap(region_by_treatment, annot=True, fmt='.1f', cmap='viridis',
                       linewidths=0.5, cbar_kws={'label': 'Percentage of Variants'})
            plt.title('Heatmap of Regulatory Region Distribution by Treatment', fontsize=14)
            plt.xlabel('Regulatory Region', fontsize=12)
            plt.ylabel('Treatment', fontsize=12)
            plt.tight_layout()
            plt.savefig(os.path.join(self.plot_dir, 'region_treatment_heatmap.png'), dpi=300)
            plt.close()
            
            # Save the data
            region_by_treatment.to_csv(os.path.join(self.data_dir, 'region_by_treatment.tsv'), sep='\t')
            
            # Compare treatment groups (Temperature vs Low Oxygen, Modified vs Non-Modified)
            # Add treatment group information
            treatment_to_group = {}
            for group, treatments in self.treatment_groups.items():
                for treatment in treatments:
                    if treatment in treatment_to_group:
                        treatment_to_group[treatment] = f"{treatment_to_group[treatment]},{group}"
                    else:
                        treatment_to_group[treatment] = group
            
            self.reg_variants['Treatment_Group'] = self.reg_variants['Treatment'].map(treatment_to_group)
            
            # Process treatment groups (handling multi-group assignments safely)
            for group in self.treatment_groups.keys():
                # Create a safer function to handle potential NaN values
                def is_in_group(x):
                    if pd.isna(x):
                        return 0
                    if isinstance(x, str) and group in x:
                        return 1
                    return 0

                self.reg_variants[f'Group_{group}'] = self.reg_variants['Treatment_Group'].apply(is_in_group)
            
            # Create separate analyses for main group comparisons
            group_comparisons = [
                ('Temperature', 'Low Oxygen'),
                ('Gene Modified', 'Non-Modified')
            ]
            
            for group1, group2 in group_comparisons:
                # Filter variants for these groups
                group1_variants = self.reg_variants[self.reg_variants[f'Group_{group1}'] == 1]
                group2_variants = self.reg_variants[self.reg_variants[f'Group_{group2}'] == 1]
                
                if len(group1_variants) == 0 or len(group2_variants) == 0:
                    print(f"Insufficient data for {group1} vs {group2} comparison")
                    continue
                
                print(f"\nComparing {group1} ({len(group1_variants)} variants) vs " + 
                      f"{group2} ({len(group2_variants)} variants)")
                
                # Compare region distributions
                group1_regions = group1_variants['Regulatory_Region'].value_counts(normalize=True) * 100
                group2_regions = group2_variants['Regulatory_Region'].value_counts(normalize=True) * 100
                
                # Combine into a single dataframe
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
                
                # Compare distances from gene
                if 'Distance' in self.reg_variants.columns:
                    # Create distance bins for visualization
                    bins = [0, 50, 100, 150, 200, 300, 500, 1000, 1500, 5000, float('inf')]
                    labels = ['0-50', '50-100', '100-150', '150-200', '200-300', '300-500', 
                             '500-1000', '1000-1500', '1500-5000', '>5000']
                    
                    group1_variants['Distance_Bin'] = pd.cut(
                        group1_variants['Distance'], bins=bins, labels=labels, right=False
                    )
                    group2_variants['Distance_Bin'] = pd.cut(
                        group2_variants['Distance'], bins=bins, labels=labels, right=False
                    )
                    
                    # Get distance distributions
                    group1_distances = group1_variants['Distance_Bin'].value_counts(normalize=True) * 100
                    group2_distances = group2_variants['Distance_Bin'].value_counts(normalize=True) * 100
                    
                    # Combine into a dataframe
                    distance_comparison = pd.DataFrame({
                        group1: group1_distances,
                        group2: group2_distances
                    }).fillna(0)
                    
                    # Save comparison
                    distance_comparison.to_csv(
                        os.path.join(self.data_dir, f'{group1}_vs_{group2}_distances.tsv'), sep='\t'
                    )
                    
                    # Plot comparison
                    plt.figure(figsize=(10, 6))
                    distance_comparison.plot(kind='bar')
                    plt.title(f'Distance Distribution: {group1} vs {group2}', fontsize=14)
                    plt.xlabel('Distance from Gene (bp)', fontsize=12)
                    plt.ylabel('Percentage of Variants', fontsize=12)
                    plt.xticks(rotation=45, ha='right')
                    plt.legend(title='Treatment Group')
                    plt.tight_layout()
                    plt.savefig(os.path.join(self.plot_dir, f'{group1}_vs_{group2}_distances.png'), dpi=300)
                    plt.close()
        
        # Compare TF binding site impacts between treatments
        if 'TF_Binding_Sites' in self.reg_variants.columns:
            # Filter for variants with TF binding sites
            tf_variants = self.reg_variants[self.reg_variants['TF_Binding_Sites'].notna()]
            
            if len(tf_variants) > 0:
                # Count TF occurrences by treatment
                tf_by_treatment = {}
                
                for treatment in self.treatments:
                    treatment_tf = tf_variants[tf_variants['Treatment'] == treatment]
                    
                    if len(treatment_tf) == 0:
                        continue
                    
                    # Extract all TFs for this treatment
                    all_tfs = []
                    for tfs in treatment_tf['TF_Binding_Sites'].dropna():
                        if isinstance(tfs, str):
                            all_tfs.extend(tfs.split(','))
                    
                    # Count occurrences
                    tf_by_treatment[treatment] = pd.Series(all_tfs).value_counts()
                
                # Plot TF distribution for each treatment with significant TF impacts
                for treatment, tf_counts in tf_by_treatment.items():
                    if len(tf_counts) > 0:
                        plt.figure(figsize=(10, 6))
                        tf_counts.plot(kind='bar', color=self.get_treatment_color(treatment))
                        plt.title(f'Transcription Factors Affected in {treatment}', fontsize=14)
                        plt.xlabel('Transcription Factor', fontsize=12)
                        plt.ylabel('Number of Variants', fontsize=12)
                        plt.xticks(rotation=45, ha='right')
                        plt.tight_layout()
                        plt.savefig(os.path.join(self.plot_dir, f'{treatment}_tf_distribution.png'), dpi=300)
                        plt.close()
                
                # Create a combined heatmap of TF impacts
                # First, create a combined dataframe with all TFs
                all_tfs = sorted(set(tf for counts in tf_by_treatment.values() for tf in counts.index))
                tf_matrix = pd.DataFrame(0, index=self.treatments, columns=all_tfs)
                
                for treatment, tf_counts in tf_by_treatment.items():
                    for tf, count in tf_counts.items():
                        tf_matrix.at[treatment, tf] = count
                
                # Normalize by total variants for each treatment
                for treatment in tf_matrix.index:
                    treatment_total = treatment_counts.get(treatment, 0)
                    if treatment_total > 0:
                        tf_matrix.loc[treatment] = tf_matrix.loc[treatment] / treatment_total * 100
                
                # Save the TF impact matrix
                tf_matrix.to_csv(os.path.join(self.data_dir, 'tf_impact_by_treatment.tsv'), sep='\t')
                
                # Plot the heatmap
                plt.figure(figsize=(14, 8))
                sns.heatmap(tf_matrix, cmap='YlOrRd', linewidths=0.5,
                           cbar_kws={'label': 'Percentage of Treatment Variants'})
                plt.title('Transcription Factor Impacts by Treatment', fontsize=14)
                plt.xlabel('Transcription Factor', fontsize=12)
                plt.ylabel('Treatment', fontsize=12)
                plt.xticks(rotation=45, ha='right')
                plt.tight_layout()
                plt.savefig(os.path.join(self.plot_dir, 'tf_impact_heatmap.png'), dpi=300)
                plt.close()
    
    def get_treatment_color(self, treatment):
        """Get consistent color for a treatment"""
        colors = {
            'WT-37': '#1f77b4',   # Blue
            'WTA': '#ff7f0e',     # Orange
            'CAS': '#2ca02c',     # Green
            'STC': '#d62728',     # Red
            'CAS-CTRL': '#7f7f7f', # Gray
            'STC-CTRL': '#7f7f7f', # Gray
            'WT-CTRL': '#7f7f7f'   # Gray
        }
        return colors.get(treatment, '#000000')  # Default to black
    
    def generate_summary_report(self):
        """Generate a comprehensive summary report of findings"""
        print("\nGenerating comprehensive summary report...")
        
        if self.reg_variants is None or len(self.reg_variants) == 0:
            print("No variants to summarize")
            return
        
        # Collect key statistics for the report
        stats = {}
        
        # Basic counts
        stats['total_variants'] = len(self.variants)
        stats['regulatory_variants'] = len(self.reg_variants)
        stats['percent_regulatory'] = (stats['regulatory_variants'] / stats['total_variants'] * 100) if stats['total_variants'] > 0 else 0
        
        # Regulatory region distribution
        if 'Regulatory_Region' in self.reg_variants.columns:
            stats['region_counts'] = self.reg_variants['Regulatory_Region'].value_counts().to_dict()
            stats['region_percent'] = (self.reg_variants['Regulatory_Region'].value_counts(normalize=True) * 100).to_dict()
        
        # Treatment distribution
        if 'Treatment' in self.reg_variants.columns:
            stats['treatment_counts'] = self.reg_variants['Treatment'].value_counts().to_dict()
            stats['treatment_percent'] = (self.reg_variants['Treatment'].value_counts(normalize=True) * 100).to_dict()
        
        # ERG gene variants
        erg_gene_ids = set(self.erg_genes['w303_gene_id'].values)
        erg_variants = self.reg_variants[self.reg_variants['Gene_ID'].isin(erg_gene_ids)]
        stats['erg_variants'] = len(erg_variants)
        stats['percent_erg'] = (len(erg_variants) / len(self.reg_variants) * 100) if len(self.reg_variants) > 0 else 0
        
        # TF binding site impacts
        stats['tf_variants'] = 0
        if 'TF_Binding_Sites' in self.reg_variants.columns:
            tf_variants = self.reg_variants[self.reg_variants['TF_Binding_Sites'].notna()]
            stats['tf_variants'] = len(tf_variants)
            stats['percent_tf'] = (len(tf_variants) / len(self.reg_variants) * 100) if len(self.reg_variants) > 0 else 0
            
            # Most common TFs
            if len(tf_variants) > 0:
                all_tfs = []
                for tfs in tf_variants['TF_Binding_Sites'].dropna():
                    if isinstance(tfs, str):
                        all_tfs.extend(tfs.split(','))
                
                if all_tfs:
                    stats['tf_counts'] = pd.Series(all_tfs).value_counts().to_dict()
        
        # Create report content
        report = ["# Regulatory Motif Analysis Report\n"]
        
        # Overview section
        report.append("## Overview\n")
        report.append(f"Total variants analyzed: {stats['total_variants']}")
        report.append(f"Regulatory variants: {stats['regulatory_variants']} ({stats['percent_regulatory']:.1f}%)")
        report.append(f"Ergosterol gene regulatory variants: {stats['erg_variants']} ({stats['percent_erg']:.1f}%)")
        report.append(f"Variants affecting TF binding sites: {stats['tf_variants']} ({stats.get('percent_tf', 0):.1f}%)\n")
        
        # Regulatory region distribution
        if 'region_counts' in stats:
            report.append("## Regulatory Region Distribution\n")
            for region, count in sorted(stats['region_counts'].items(), key=lambda x: x[1], reverse=True):
                percent = stats['region_percent'][region]
                report.append(f"{region}: {count} variants ({percent:.1f}%)")
            report.append("")
        
        # Treatment distribution
        if 'treatment_counts' in stats:
            report.append("## Treatment Distribution\n")
            for treatment, count in sorted(stats['treatment_counts'].items(), key=lambda x: x[1], reverse=True):
                percent = stats['treatment_percent'][treatment]
                report.append(f"{treatment}: {count} variants ({percent:.1f}%)")
            report.append("")
        
        # TF binding site impacts
        if 'tf_counts' in stats and stats['tf_counts']:
            report.append("## Transcription Factor Binding Sites\n")
            report.append(f"Most frequently affected transcription factors:")
            for tf, count in sorted(stats['tf_counts'].items(), key=lambda x: x[1], reverse=True)[:10]:
                report.append(f"- {tf}: {count} variants")
            report.append("")
        
        # Key findings section
        report.append("## Key Findings\n")
        
        # Add regulatory region observations
        if 'region_percent' in stats and stats['region_percent']:
            # Most common regulatory region
            top_region = max(stats['region_percent'].items(), key=lambda x: x[1])
            report.append(f"- Most common regulatory region: {top_region[0]} ({top_region[1]:.1f}% of variants)")
            
            # Note promoter element distribution
            promoters = ['core_promoter', 'TATA_region', 'tss_proximal', 'UAS_proximal', 'UAS_distal']
            promoter_percent = sum(stats['region_percent'].get(p, 0) for p in promoters)
            report.append(f"- Promoter elements account for {promoter_percent:.1f}% of regulatory variants")
            
            # Note UTR/terminator distribution
            downstream = ['five_prime_UTR', 'three_prime_UTR', 'terminator', 'downstream_element']
            downstream_percent = sum(stats['region_percent'].get(d, 0) for d in downstream)
            report.append(f"- UTR and terminator regions account for {downstream_percent:.1f}% of regulatory variants")
        
        # Add ERG gene and TF observations
        report.append(f"- Ergosterol gene regulatory variants: {stats['percent_erg']:.1f}% of all regulatory variants")
        if 'percent_tf' in stats and stats['percent_tf'] > 0:
            report.append(f"- Transcription factor binding site impacts: {stats['percent_tf']:.1f}% of regulatory variants")
        
        # Treatment-specific patterns
        if 'Treatment' in self.reg_variants.columns and 'Regulatory_Region' in self.reg_variants.columns:
            report.append("\n## Treatment-Specific Patterns\n")
            
            # Get region distribution for each treatment
            for treatment in sorted(stats.get('treatment_counts', {}).keys()):
                treatment_data = self.reg_variants[self.reg_variants['Treatment'] == treatment]
                if len(treatment_data) > 0:
                    report.append(f"### {treatment}")
                    
                    # Get top regulatory regions for this treatment
                    t_region_counts = treatment_data['Regulatory_Region'].value_counts()
                    if not t_region_counts.empty:
                        t_region_percent = (t_region_counts / len(treatment_data) * 100).round(1)
                        
                        for region, percent in t_region_percent.nlargest(3).items():
                            report.append(f"- {region}: {percent}%")
                    
                    # Check for TF impacts
                    if 'TF_Binding_Sites' in treatment_data.columns:
                        tf_variants = treatment_data[treatment_data['TF_Binding_Sites'].notna()]
                        if len(tf_variants) > 0:
                            all_tfs = []
                            for tfs in tf_variants['TF_Binding_Sites'].dropna():
                                if isinstance(tfs, str):
                                    all_tfs.extend(tfs.split(','))
                            
                            if all_tfs:
                                tf_counts = pd.Series(all_tfs).value_counts()
                                top_tf = tf_counts.index[0] if not tf_counts.empty else "None"
                                report.append(f"- Most affected TF: {top_tf}")
                    
                    report.append("")
        
        # Biological interpretation
        report.append("## Biological Interpretation\n")
        
        # Add interpretation based on results
        if 'region_percent' in stats:
            # Interpret upstream vs downstream distribution
            upstream_regions = ['core_promoter', 'TATA_region', 'tss_proximal', 'UAS_proximal', 
                              'UAS_distal', 'far_upstream']
            downstream_regions = ['five_prime_UTR', 'three_prime_UTR', 'terminator', 'downstream_element']
            
            upstream_percent = sum(stats['region_percent'].get(r, 0) for r in upstream_regions)
            downstream_percent = sum(stats['region_percent'].get(r, 0) for r in downstream_regions)
            
            if upstream_percent > downstream_percent:
                report.append("- **Promoter-Biased Regulation**: A higher proportion of variants occur in promoter/upstream regions, " +
                             "suggesting adaptation through changes in transcriptional regulation rather than post-transcriptional processes.")
            else:
                report.append("- **Post-Transcriptional Regulation**: A higher proportion of variants occur in UTR/terminator regions, " +
                             "suggesting adaptation through changes in RNA stability, translation efficiency, or transcript processing.")
            
            # Interpret TATA/UAS distribution
            tata_percent = stats['region_percent'].get('TATA_region', 0)
            uas_percent = stats['region_percent'].get('UAS_proximal', 0) + stats['region_percent'].get('UAS_distal', 0)
            
            if uas_percent > tata_percent:
                report.append("- **Distal Regulatory Preference**: Variants are more concentrated in UAS regions rather than core promoter elements, " +
                             "suggesting adaptation through modulation of condition-specific transcription factor binding rather than basal transcription machinery.")
            else:
                report.append("- **Core Promoter Modulation**: Variants are concentrated near core promoter elements, " +
                             "suggesting direct modulation of the basal transcription machinery.")
        
        # Interpret ERG gene patterns
        if stats['percent_erg'] > 50:
            report.append("- **Strong ERG Regulatory Focus**: The high proportion of variants associated with ergosterol pathway genes " +
                         "confirms the central role of sterol metabolism in the adaptive response to the experimental conditions.")
        
        # Interpret treatment patterns
        if 'treatment_counts' in stats and len(stats['treatment_counts']) >= 4:
            # Compare temperature vs low oxygen treatments
            temp_treatments = ['WT-37', 'CAS']
            oxy_treatments = ['WTA', 'STC']
            
            temp_variants = sum(stats['treatment_counts'].get(t, 0) for t in temp_treatments)
            oxy_variants = sum(stats['treatment_counts'].get(t, 0) for t in oxy_treatments)
            
            if temp_variants > oxy_variants:
                report.append("- **Temperature Adaptation Prominence**: Temperature adaptation conditions show a higher number of regulatory variants " +
                             "compared to low oxygen conditions, suggesting a more complex regulatory response to thermal stress.")
            elif oxy_variants > temp_variants:
                report.append("- **Oxygen Response Prominence**: Low oxygen adaptation conditions show a higher number of regulatory variants " +
                             "compared to temperature conditions, suggesting a more complex regulatory response to oxygen limitation.")
        
        # Write report to file
        report_path = os.path.join(self.output_dir, 'regulatory_motif_analysis_report.txt')
        with open(report_path, 'w') as f:
            f.write('\n'.join(report))
        
        print(f"Summary report saved to {report_path}")
        
        # Save key statistics as JSON
        stats_path = os.path.join(self.data_dir, 'regulatory_motif_statistics.json')
        with open(stats_path, 'w') as f:
            # Convert some stats to strings to ensure JSON serialization
            serializable_stats = {}
            for k, v in stats.items():
                if isinstance(v, dict):
                    serializable_stats[k] = {str(kk): vv for kk, vv in v.items()}
                else:
                    serializable_stats[k] = v
            
            json.dump(serializable_stats, f, indent=2)
        
        print(f"Statistics saved to {stats_path}")
        
        return stats
    
    def run_analysis(self):
        """Run the complete regulatory motif analysis pipeline"""
        print("\n=== Starting Regulatory Motif Analysis ===\n")
        
        # Step 1: Prepare data
        self.prepare_data()
        
        # Step 2: Classify regulatory regions
        self.classify_regulatory_regions()
        
        # Step 3: Identify TF binding sites (if reference genome available)
        if self.genome_seq:
            self.identify_tf_binding_sites()
        
        # Step 4: Analyze ERG gene regulatory elements
        self.analyze_erg_regulatory_elements()
        
        # Step 5: Analyze treatment differences
        self.analyze_treatment_differences()
        
        # Step 6: Generate summary report
        self.generate_summary_report()
        
        print("\n=== Regulatory Motif Analysis Complete ===")
        
        # Return the processed variants for further analysis
        return self.reg_variants

def main():
    """Main function to parse arguments and run the analysis"""
    parser = argparse.ArgumentParser(description='Analyze regulatory motifs and TF binding sites')
    
    # Required arguments
    parser.add_argument('--output-dir', default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/regulatory_analysis/motif_analysis',
                        help='Output directory for results')
    parser.add_argument('--gene-mapping', default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/reference/gene_mapping_full.tsv',
                        help='Gene mapping file with coordinates')
    parser.add_argument('--variants', default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/gene_variants_expanded/all_gene_variants.tsv',
                        help='Variants TSV file')
    parser.add_argument('--erg-genes', default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/reference/genes_of_interest_mapping.tsv',
                        help='ERG genes mapping file')
    
    # Optional arguments
    parser.add_argument('--reference-genome', 
                        default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/reference/w303_chromosomal.fasta',
                        help='Reference genome FASTA file')
    parser.add_argument('--motif-database', default=None,
                        help='JSON file with TF motif database')
    
    args = parser.parse_args()
    
    # Initialize analyzer
    analyzer = RegulatoryMotifAnalyzer(
        args.output_dir,
        args.gene_mapping,
        args.variants,
        args.erg_genes,
        args.reference_genome,
        args.motif_database
    )
    
    # Run analysis
    analyzer.run_analysis()

if __name__ == "__main__":
    main()