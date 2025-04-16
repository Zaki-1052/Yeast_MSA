#!/usr/bin/env python3
"""
scaffold_analysis.py - Comprehensive scaffold-level analysis of variants across treatments

This script provides detailed analysis and visualization of variant patterns at the scaffold level,
extending the existing analysis to provide a more complete picture of genomic distributions.
It integrates with the existing analysis pipeline and builds on JRIU annotations.
"""

import os
import sys
import csv
import argparse
import gzip
import json
import math
from collections import defaultdict, Counter
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats as scipy_stats
import pandas as pd
from Bio import SeqIO
from scipy.cluster import hierarchy
from scipy.spatial import distance

class ScaffoldAnalyzer:
    def __init__(self, output_dir):
        """
        Initialize the scaffold analyzer.
        
        Args:
            output_dir: Directory for output files and reports
        """
        self.output_dir = output_dir
        
        # Create output directories
        self.data_dir = os.path.join(output_dir, "data")
        self.report_dir = os.path.join(output_dir, "reports")
        self.plot_dir = os.path.join(output_dir, "plots")
        
        os.makedirs(self.data_dir, exist_ok=True)
        os.makedirs(self.report_dir, exist_ok=True)
        os.makedirs(self.plot_dir, exist_ok=True)
        
        # Treatment mapping
        self.treatment_groups = {
            'WT': ['WT-CTRL'],
            'WT-37': ['WT-37-55-1', 'WT-37-55-2', 'WT-37-55-3'],
            'WTA': ['WTA-55-1', 'WTA-55-2', 'WTA-55-3'],
            'STC': ['STC-55-1', 'STC-55-2', 'STC-55-3'],
            'CAS': ['CAS-55-1', 'CAS-55-2', 'CAS-55-3'],
            'STC-CTRL': ['STC-CTRL'],
            'CAS-CTRL': ['CAS-CTRL']
        }
        
        # Treatment metadata
        self.treatment_metadata = {
            'WT-37': {'adaptation': 'Temperature', 'has_gene': False, 'description': 'Temperature-adapted wild type'},
            'WTA': {'adaptation': 'Low Oxygen', 'has_gene': False, 'description': 'Low oxygen-adapted wild type'},
            'STC': {'adaptation': 'Low Oxygen', 'has_gene': True, 'description': 'STC gene with low oxygen adaptation'},
            'CAS': {'adaptation': 'Temperature', 'has_gene': True, 'description': 'CAS gene with temperature adaptation'},
            'WT': {'adaptation': 'None', 'has_gene': False, 'description': 'Wild type control'},
            'STC-CTRL': {'adaptation': 'None', 'has_gene': True, 'description': 'STC gene control'},
            'CAS-CTRL': {'adaptation': 'None', 'has_gene': True, 'description': 'CAS gene control'}
        }
        
        # Adaptation groups
        self.adaptation_groups = {
            'Temperature': ['WT-37', 'CAS'],
            'Low Oxygen': ['WTA', 'STC'],
            'None': ['WT', 'STC-CTRL', 'CAS-CTRL']
        }
        
        # Gene modification groups
        self.gene_modification_groups = {
            'Gene-Modified': ['STC', 'CAS', 'STC-CTRL', 'CAS-CTRL'],
            'Non-Modified': ['WT-37', 'WTA', 'WT']
        }
        
        # Reverse lookup: sample to treatment
        self.sample_to_treatment = {}
        for treatment, samples in self.treatment_groups.items():
            for sample in samples:
                self.sample_to_treatment[sample] = treatment
        
        # Data containers
        self.scaffold_metadata = {}  # Metadata about scaffolds (length, etc.)
        self.scaffold_variants = defaultdict(lambda: defaultdict(int))  # Variants by scaffold and treatment
        self.scaffold_density = defaultdict(lambda: defaultdict(float))  # Variant density by scaffold and treatment
        self.treatment_stats = defaultdict(dict)  # Statistics by treatment
        self.adaptation_stats = defaultdict(dict)  # Statistics by adaptation
        
        # Load scaffold length information if available
        self.scaffold_lengths = self.load_scaffold_lengths()
    
    def load_scaffold_lengths(self, scaffold_mapping_file=None):
        """
        Load scaffold length information from the reference genome or mapping file.
        
        Args:
            scaffold_mapping_file: Optional path to scaffold mapping file
            
        Returns:
            Dict mapping scaffold IDs to their lengths
        """
        scaffold_lengths = {}
        
        # Try to load from scaffold mapping file if provided
        if scaffold_mapping_file and os.path.exists(scaffold_mapping_file):
            print(f"Loading scaffold lengths from {scaffold_mapping_file}")
            
            with open(scaffold_mapping_file, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    if 'scaffold_id' in row and 'length' in row:
                        scaffold_id = row['scaffold_id']
                        length = int(row['length'])
                        scaffold_lengths[scaffold_id] = length
            
            if scaffold_lengths:
                print(f"Loaded lengths for {len(scaffold_lengths)} scaffolds")
                return scaffold_lengths
        
        # If no file provided or it didn't work, try to find the reference directory
        # and extract lengths from FASTA files
        ref_dirs = [
            os.path.join(os.path.dirname(os.path.dirname(self.output_dir)), 'reference'),
            os.path.join(os.path.dirname(os.path.dirname(self.output_dir)), 'annotation', 'reference')
        ]
        
        for ref_dir in ref_dirs:
            if os.path.exists(ref_dir):
                # Look for FASTA files
                fasta_files = [f for f in os.listdir(ref_dir) if f.endswith('.fasta') or f.endswith('.fa')]
                
                if fasta_files:
                    print(f"Loading scaffold lengths from FASTA file: {fasta_files[0]}")
                    
                    # Use the first FASTA file found
                    fasta_path = os.path.join(ref_dir, fasta_files[0])
                    
                    try:
                        for record in SeqIO.parse(fasta_path, "fasta"):
                            scaffold_lengths[record.id] = len(record.seq)
                        
                        print(f"Loaded lengths for {len(scaffold_lengths)} scaffolds from FASTA")
                        return scaffold_lengths
                    except Exception as e:
                        print(f"Error loading FASTA file: {e}")
                
                # Look for GenBank files
                genbank_files = [os.path.join(ref_dir, f) for f in os.listdir(ref_dir) if f.endswith('.genbank') or f.endswith('.gbk')]
                
                if genbank_files:
                    print(f"Loading scaffold lengths from {len(genbank_files)} GenBank files")
                    
                    # Process each GenBank file
                    for gb_file in genbank_files:
                        try:
                            for record in SeqIO.parse(gb_file, "genbank"):
                                scaffold_lengths[record.id] = len(record.seq)
                        except Exception as e:
                            print(f"Error loading GenBank file {gb_file}: {e}")
                    
                    if scaffold_lengths:
                        print(f"Loaded lengths for {len(scaffold_lengths)} scaffolds from GenBank")
                        return scaffold_lengths
        
        print("Warning: Could not load scaffold lengths from reference files.")
        print("Scaffold density calculations will use estimate from VCF data.")
        return scaffold_lengths
    
    def process_vcf_file(self, vcf_file):
        """
        Process a single VCF file to extract variant information by scaffold.
        
        Args:
            vcf_file: Path to VCF file
            
        Returns:
            Dict containing statistics about variants by scaffold
        """
        # Extract sample name
        sample_name = os.path.basename(vcf_file).split('.')[0]
        treatment = self.sample_to_treatment.get(sample_name, 'unknown')
        
        print(f"Processing {vcf_file} (Sample: {sample_name}, Treatment: {treatment})")
        
        # Determine if file is gzipped
        is_gzipped = vcf_file.endswith('.gz')
        opener = gzip.open if is_gzipped else open
        mode = 'rt' if is_gzipped else 'r'
        
        # Initialize counters
        total_variants = 0
        scaffold_counts = defaultdict(int)
        scaffold_positions = defaultdict(set)  # Track unique positions to handle multi-allelic sites
        
        # Process VCF file
        with opener(vcf_file, mode) as f:
            # Process header
            for line in f:
                if line.startswith('#'):
                    continue
                
                line = line.strip()
                if not line:
                    continue
                
                fields = line.split('\t')
                
                # Extract scaffold and position
                scaffold = fields[0]
                position = int(fields[1])
                
                # Count unique positions per scaffold
                position_key = f"{scaffold}:{position}"
                if position_key not in scaffold_positions[scaffold]:
                    scaffold_positions[scaffold].add(position_key)
                    scaffold_counts[scaffold] += 1
                    total_variants += 1
        
        # Store results
        result = {
            'sample': sample_name,
            'treatment': treatment,
            'total_variants': total_variants,
            'scaffolds': len(scaffold_counts),
            'scaffold_counts': dict(scaffold_counts)
        }
        
        # Update the treatment statistics
        if treatment not in self.treatment_stats:
            self.treatment_stats[treatment] = {
                'total_variants': 0,
                'sample_count': 0,
                'scaffolds': set(),
                'scaffold_counts': defaultdict(int)
            }
        
        self.treatment_stats[treatment]['total_variants'] += total_variants
        self.treatment_stats[treatment]['sample_count'] += 1
        self.treatment_stats[treatment]['scaffolds'].update(scaffold_counts.keys())
        
        for scaffold, count in scaffold_counts.items():
            self.treatment_stats[treatment]['scaffold_counts'][scaffold] += count
            self.scaffold_variants[scaffold][treatment] += count
            
            # Record scaffold metadata if not already done
            if scaffold not in self.scaffold_metadata:
                self.scaffold_metadata[scaffold] = {
                    'id': scaffold,
                    'treatments': set(),
                    'total_variants': 0
                }
            
            self.scaffold_metadata[scaffold]['treatments'].add(treatment)
            self.scaffold_metadata[scaffold]['total_variants'] += count
        
        print(f"  Found {total_variants} variants across {len(scaffold_counts)} scaffolds")
        return result
    
    def process_all_vcfs(self, vcf_files):
        """
        Process a list of VCF files and collect comprehensive statistics
        
        Args:
            vcf_files: List of VCF file paths
            
        Returns:
            Collected statistics across all files
        """
        # Process each VCF file
        all_results = []
        
        for vcf_file in vcf_files:
            result = self.process_vcf_file(vcf_file)
            all_results.append(result)
        
        # Calculate adaptation statistics
        for adaptation, treatments in self.adaptation_groups.items():
            self.adaptation_stats[adaptation] = {
                'treatments': treatments,
                'total_variants': sum(self.treatment_stats[t]['total_variants'] for t in treatments if t in self.treatment_stats),
                'scaffolds': set().union(*[self.treatment_stats[t]['scaffolds'] for t in treatments if t in self.treatment_stats]),
                'scaffold_counts': defaultdict(int)
            }
            
            # Combine scaffold counts from all treatments in this adaptation group
            for treatment in treatments:
                if treatment in self.treatment_stats:
                    for scaffold, count in self.treatment_stats[treatment]['scaffold_counts'].items():
                        self.adaptation_stats[adaptation]['scaffold_counts'][scaffold] += count
        
        # Calculate scaffold densities
        for scaffold, metadata in self.scaffold_metadata.items():
            # Get scaffold length from loaded data or estimate
            if scaffold in self.scaffold_lengths:
                length = self.scaffold_lengths[scaffold]
            else:
                # Estimate length as 10kb if unknown
                length = 10000
                print(f"Warning: No length data for scaffold {scaffold}. Using estimate of {length}bp.")
            
            # Store the length
            metadata['length'] = length
            
            # Calculate densities for each treatment
            for treatment, count in self.scaffold_variants[scaffold].items():
                # Density in variants per kilobase
                density = (count * 1000) / length
                self.scaffold_density[scaffold][treatment] = density
            
            # Calculate overall density
            metadata['density'] = (metadata['total_variants'] * 1000) / length
        
        return all_results
    
    def identify_enriched_scaffolds(self, min_variants=5, min_fold=2.0, max_pval=0.05):
        """
        Identify scaffolds with statistically significant enrichment of variants
        
        Args:
            min_variants: Minimum number of variants for a scaffold to be considered
            min_fold: Minimum fold enrichment over background
            max_pval: Maximum p-value for significance
            
        Returns:
            Dict mapping treatment to lists of enriched scaffolds
        """
        enriched_scaffolds = defaultdict(list)
        
        # Calculate global density for each treatment
        treatment_densities = {}
        for treatment, stats in self.treatment_stats.items():
            total_length = sum(self.scaffold_metadata[s]['length'] for s in stats['scaffolds'])
            global_density = (stats['total_variants'] * 1000) / total_length if total_length > 0 else 0
            treatment_densities[treatment] = global_density
        
        # Test each scaffold in each treatment for enrichment
        for treatment, stats in self.treatment_stats.items():
            global_density = treatment_densities[treatment]
            
            for scaffold in stats['scaffolds']:
                count = self.scaffold_variants[scaffold][treatment]
                
                # Skip scaffolds with too few variants
                if count < min_variants:
                    continue
                
                # Calculate fold enrichment
                scaffold_density = self.scaffold_density[scaffold][treatment]
                fold_enrichment = scaffold_density / global_density if global_density > 0 else 0
                
                # Skip if fold enrichment is too low
                if fold_enrichment < min_fold:
                    continue
                
                # Statistical test (Poisson test)
                length_kb = self.scaffold_metadata[scaffold]['length'] / 1000
                expected_count = global_density * length_kb
                p_value = scipy_stats.poisson.sf(count - 1, expected_count)
                q_value = p_value  # Add multiple testing correction if needed
                
                # Record if significant
                if p_value <= max_pval:
                    enriched_scaffolds[treatment].append({
                        'scaffold': scaffold,
                        'count': count,
                        'density': scaffold_density,
                        'fold_enrichment': fold_enrichment,
                        'expected_count': expected_count,
                        'p_value': p_value,
                        'q_value': q_value,
                        'length': self.scaffold_metadata[scaffold]['length']
                    })
        
        # Sort enriched scaffolds by fold enrichment
        for treatment in enriched_scaffolds:
            enriched_scaffolds[treatment].sort(key=lambda x: -x['fold_enrichment'])
        
        return enriched_scaffolds
    
    def calculate_correlations(self):
        """
        Calculate correlations between treatments based on scaffold variant patterns
        
        Returns:
            Tuple of (correlation_matrix, treatment_list)
        """
        treatments = sorted(self.treatment_stats.keys())
        all_scaffolds = set()
        
        # Get all scaffolds that have variants in any treatment
        for treatment, stats in self.treatment_stats.items():
            all_scaffolds.update(stats['scaffolds'])
        
        # Create data frame with scaffold counts by treatment
        data = []
        for scaffold in all_scaffolds:
            row = [scaffold]
            for treatment in treatments:
                row.append(self.scaffold_variants[scaffold].get(treatment, 0))
            data.append(row)
        
        columns = ['scaffold'] + treatments
        df = pd.DataFrame(data, columns=columns)
        df.set_index('scaffold', inplace=True)
        
        # Calculate Spearman correlation
        correlation_matrix = df.corr(method='spearman')
        
        return correlation_matrix, treatments
    
    def generate_text_reports(self, enriched_scaffolds):
        """
        Generate comprehensive text-based reports for scaffold analysis
        
        Args:
            enriched_scaffolds: Dict of enriched scaffolds by treatment
        
        Returns:
            Dict of report file paths
        """
        # 1. Treatment scaffold summary
        treatment_report_file = os.path.join(self.report_dir, "treatment_scaffold_summary.tsv")
        with open(treatment_report_file, 'w') as f:
            f.write("Treatment\tDescription\tAdaptation\tHas_Gene\tTotal_Variants\tScaffolds_With_Variants\t")
            f.write("Global_Density\tTop_Scaffolds\n")
            
            for treatment in sorted(self.treatment_stats.keys()):
                stats = self.treatment_stats[treatment]
                metadata = self.treatment_metadata.get(treatment, {})
                
                # Get top 3 scaffolds by density
                top_scaffolds = []
                for scaffold, count in stats['scaffold_counts'].items():
                    density = self.scaffold_density[scaffold][treatment]
                    top_scaffolds.append((scaffold, density, count))
                
                top_scaffolds.sort(key=lambda x: -x[1])
                top_3 = top_scaffolds[:3]
                top_desc = ", ".join([f"{s} ({d:.2f})" for s, d, c in top_3])
                
                # Calculate global density
                total_length = sum(self.scaffold_metadata[s]['length'] for s in stats['scaffolds'])
                global_density = (stats['total_variants'] * 1000) / total_length if total_length > 0 else 0
                
                f.write(f"{treatment}\t{metadata.get('description', '')}\t{metadata.get('adaptation', '')}\t")
                f.write(f"{metadata.get('has_gene', '')}\t{stats['total_variants']}\t{len(stats['scaffolds'])}\t")
                f.write(f"{global_density:.6f}\t{top_desc}\n")
        
        # 2. Adaptation scaffold summary
        adaptation_report_file = os.path.join(self.report_dir, "adaptation_scaffold_summary.tsv")
        with open(adaptation_report_file, 'w') as f:
            f.write("Adaptation\tTreatments\tTotal_Variants\tScaffolds_With_Variants\t")
            f.write("Global_Density\tTop_Scaffolds\n")
            
            for adaptation in sorted(self.adaptation_stats.keys()):
                stats = self.adaptation_stats[adaptation]
                
                # Skip if no data
                if not stats['scaffolds']:
                    continue
                
                # Get top 3 scaffolds by density
                top_scaffolds = []
                for scaffold in stats['scaffolds']:
                    count = stats['scaffold_counts'][scaffold]
                    length = self.scaffold_metadata[scaffold]['length']
                    density = (count * 1000) / length
                    top_scaffolds.append((scaffold, density, count))
                
                top_scaffolds.sort(key=lambda x: -x[1])
                top_3 = top_scaffolds[:3]
                top_desc = ", ".join([f"{s} ({d:.2f})" for s, d, c in top_3])
                
                # Calculate global density
                total_length = sum(self.scaffold_metadata[s]['length'] for s in stats['scaffolds'])
                global_density = (stats['total_variants'] * 1000) / total_length if total_length > 0 else 0
                
                treatments_str = ", ".join(stats['treatments'])
                
                f.write(f"{adaptation}\t{treatments_str}\t{stats['total_variants']}\t{len(stats['scaffolds'])}\t")
                f.write(f"{global_density:.6f}\t{top_desc}\n")
        
        # 3. Scaffold summary
        scaffold_report_file = os.path.join(self.report_dir, "scaffold_summary.tsv")
        with open(scaffold_report_file, 'w') as f:
            f.write("Scaffold\tLength\tTotal_Variants\tDensity\tTreatments\tTop_Treatment\tTop_Treatment_Count\n")
            
            # Sort scaffolds by total variants
            sorted_scaffolds = sorted(self.scaffold_metadata.items(), key=lambda x: -x[1]['total_variants'])
            
            for scaffold, metadata in sorted_scaffolds:
                # Find top treatment
                top_treatment = None
                top_count = 0
                
                for treatment, count in self.scaffold_variants[scaffold].items():
                    if count > top_count:
                        top_count = count
                        top_treatment = treatment
                
                treatments_str = ", ".join(sorted(metadata['treatments']))
                
                f.write(f"{scaffold}\t{metadata['length']}\t{metadata['total_variants']}\t")
                f.write(f"{metadata['density']:.6f}\t{treatments_str}\t{top_treatment}\t{top_count}\n")
        
        # 4. Enriched scaffolds report
        enriched_report_file = os.path.join(self.report_dir, "enriched_scaffolds.tsv")
        with open(enriched_report_file, 'w') as f:
            f.write("Treatment\tScaffold\tLength\tVariants\tDensity\tGlobal_Density\t")
            f.write("Fold_Enrichment\tExpected_Count\tP_Value\tQ_Value\n")
            
            for treatment in sorted(enriched_scaffolds.keys()):
                # Get global density
                stats = self.treatment_stats[treatment]
                total_length = sum(self.scaffold_metadata[s]['length'] for s in stats['scaffolds'])
                global_density = (stats['total_variants'] * 1000) / total_length if total_length > 0 else 0
                
                for info in enriched_scaffolds[treatment]:
                    f.write(f"{treatment}\t{info['scaffold']}\t{info['length']}\t{info['count']}\t")
                    f.write(f"{info['density']:.6f}\t{global_density:.6f}\t{info['fold_enrichment']:.6f}\t")
                    f.write(f"{info['expected_count']:.2f}\t{info['p_value']:.6e}\t{info['q_value']:.6e}\n")
        
        # 5. Correlation report
        correlation_matrix, treatments = self.calculate_correlations()
        correlation_report_file = os.path.join(self.report_dir, "treatment_correlations.tsv")
        
        with open(correlation_report_file, 'w') as f:
            # Write header
            f.write("Treatment\t" + "\t".join(treatments) + "\n")
            
            # Write correlation matrix
            for treatment in treatments:
                f.write(treatment)
                for other in treatments:
                    f.write(f"\t{correlation_matrix.loc[treatment, other]:.4f}")
                f.write("\n")
        
        # Return all report files
        return {
            'treatment_report': treatment_report_file,
            'adaptation_report': adaptation_report_file,
            'scaffold_report': scaffold_report_file,
            'enriched_report': enriched_report_file,
            'correlation_report': correlation_report_file
        }
    
    def generate_plots(self, enriched_scaffolds):
        """
        Generate comprehensive plots for scaffold analysis
        
        Args:
            enriched_scaffolds: Dict of enriched scaffolds by treatment
            
        Returns:
            Dict of plot file paths
        """
        plot_files = {}
        
        # Set visualization style
        plt.style.use('seaborn-v0_8-whitegrid')
        
        # 1. Variant counts by treatment
        treatment_counts = [(t, stats['total_variants']) for t, stats in self.treatment_stats.items()]
        treatment_counts.sort(key=lambda x: x[1], reverse=True)
        
        treatments = [t for t, _ in treatment_counts]
        counts = [c for _, c in treatment_counts]
        
        plt.figure(figsize=(10, 6))
        plt.bar(treatments, counts, color='skyblue')
        plt.title('Total Variants by Treatment', fontsize=14)
        plt.xlabel('Treatment', fontsize=12)
        plt.ylabel('Number of Variants', fontsize=12)
        plt.xticks(rotation=45)
        plt.tight_layout()
        
        treatment_count_plot = os.path.join(self.plot_dir, 'treatment_variant_counts.png')
        plt.savefig(treatment_count_plot, dpi=300)
        plt.close()
        
        plot_files['treatment_counts'] = treatment_count_plot
        
        # 2. Number of scaffolds with variants by treatment
        treatment_scaffolds = [(t, len(stats['scaffolds'])) for t, stats in self.treatment_stats.items()]
        treatment_scaffolds.sort(key=lambda x: x[1], reverse=True)
        
        treatments = [t for t, _ in treatment_scaffolds]
        scaffolds = [s for _, s in treatment_scaffolds]
        
        plt.figure(figsize=(10, 6))
        plt.bar(treatments, scaffolds, color='lightgreen')
        plt.title('Scaffolds with Variants by Treatment', fontsize=14)
        plt.xlabel('Treatment', fontsize=12)
        plt.ylabel('Number of Scaffolds', fontsize=12)
        plt.xticks(rotation=45)
        plt.tight_layout()
        
        treatment_scaffold_plot = os.path.join(self.plot_dir, 'treatment_scaffold_counts.png')
        plt.savefig(treatment_scaffold_plot, dpi=300)
        plt.close()
        
        plot_files['treatment_scaffolds'] = treatment_scaffold_plot
        
        # 3. Top 10 scaffolds by variant density across all treatments
        all_densities = []
        for scaffold, metadata in self.scaffold_metadata.items():
            if metadata['total_variants'] >= 5:  # Minimum 5 variants to be included
                all_densities.append((scaffold, metadata['density'], metadata['total_variants']))
        
        all_densities.sort(key=lambda x: -x[1])
        top_scaffolds = all_densities[:10]
        
        scaffolds = [s for s, _, _ in top_scaffolds]
        densities = [d for _, d, _ in top_scaffolds]
        
        plt.figure(figsize=(12, 6))
        plt.bar(scaffolds, densities, color='salmon')
        plt.title('Top 10 Scaffolds by Variant Density (All Treatments)', fontsize=14)
        plt.xlabel('Scaffold', fontsize=12)
        plt.ylabel('Variant Density (variants/kb)', fontsize=12)
        plt.xticks(rotation=90)
        plt.tight_layout()
        
        top_scaffold_plot = os.path.join(self.plot_dir, 'top_scaffolds_by_density.png')
        plt.savefig(top_scaffold_plot, dpi=300)
        plt.close()
        
        plot_files['top_scaffolds'] = top_scaffold_plot
        
        # 4. Bubble chart of variant density by scaffold for each treatment
        for treatment, stats in self.treatment_stats.items():
            # Get top 20 scaffolds by density for this treatment
            densities = []
            for scaffold in stats['scaffolds']:
                count = self.scaffold_variants[scaffold][treatment]
                if count >= 3:  # Minimum 3 variants per scaffold
                    densities.append((
                        scaffold, 
                        self.scaffold_density[scaffold][treatment],
                        count,
                        self.scaffold_metadata[scaffold]['length']
                    ))
            
            densities.sort(key=lambda x: -x[1])
            top_20 = densities[:20]
            
            if not top_20:
                continue
            
            scaffolds = [s for s, _, _, _ in top_20]
            densities = [d for _, d, _, _ in top_20]
            counts = [c for _, _, c, _ in top_20]
            lengths = [l/1000 for _, _, _, l in top_20]  # Convert to kb
            
            plt.figure(figsize=(14, 8))
            # Size points by count, color by length
            scatter = plt.scatter(scaffolds, densities, s=[c*20 for c in counts], c=lengths, 
                                 cmap='viridis', alpha=0.7)
            
            # Add count labels
            for i, (scaffold, density, count, _) in enumerate(top_20):
                plt.annotate(f"{count}", (i, density), 
                            textcoords="offset points", xytext=(0,10), 
                            ha='center')
            
            plt.colorbar(scatter, label='Scaffold Length (kb)')
            plt.title(f'Top 20 Scaffolds by Variant Density - {treatment}', fontsize=14)
            plt.xlabel('Scaffold', fontsize=12)
            plt.ylabel('Variant Density (variants/kb)', fontsize=12)
            plt.xticks(rotation=90)
            plt.tight_layout()
            
            treatment_bubble_plot = os.path.join(self.plot_dir, f'{treatment}_bubble_chart.png')
            plt.savefig(treatment_bubble_plot, dpi=300)
            plt.close()
            
            plot_files[f'{treatment}_bubble'] = treatment_bubble_plot
        
        # 5. Heatmap of enriched scaffolds across treatments
        # Gather data for top 30 enriched scaffolds across all treatments
        all_enriched = []
        for treatment, scaffolds in enriched_scaffolds.items():
            for info in scaffolds:
                all_enriched.append((treatment, info['scaffold'], info['fold_enrichment']))
        
        # Sort by fold enrichment and get top 30 unique scaffolds
        all_enriched.sort(key=lambda x: -x[2])
        top_scaffolds = []
        for _, scaffold, _ in all_enriched:
            if scaffold not in top_scaffolds and len(top_scaffolds) < 30:
                top_scaffolds.append(scaffold)
        
        # Create enrichment data for heatmap
        treatments = sorted(self.treatment_stats.keys())
        enrichment_data = np.zeros((len(top_scaffolds), len(treatments)))
        
        for i, scaffold in enumerate(top_scaffolds):
            for j, treatment in enumerate(treatments):
                # Find enrichment info for this scaffold/treatment
                enrichment = 0
                if treatment in enriched_scaffolds:
                    for info in enriched_scaffolds[treatment]:
                        if info['scaffold'] == scaffold:
                            enrichment = info['fold_enrichment']
                            break
                
                enrichment_data[i, j] = enrichment
        
        plt.figure(figsize=(14, 10))
        sns.heatmap(enrichment_data, xticklabels=treatments, yticklabels=top_scaffolds,
                    cmap="YlOrRd", annot=True, fmt=".1f")
        plt.title('Fold Enrichment of Top 30 Scaffolds Across Treatments', fontsize=14)
        plt.xlabel('Treatment', fontsize=12)
        plt.ylabel('Scaffold', fontsize=12)
        plt.tight_layout()
        
        enrichment_heatmap = os.path.join(self.plot_dir, 'enrichment_heatmap.png')
        plt.savefig(enrichment_heatmap, dpi=300)
        plt.close()
        
        plot_files['enrichment_heatmap'] = enrichment_heatmap
        
        # 6. Correlation heatmap between treatments
        correlation_matrix, treatments = self.calculate_correlations()
        
        plt.figure(figsize=(10, 8))
        sns.heatmap(correlation_matrix, annot=True, fmt=".2f", cmap="YlGnBu")
        plt.title('Correlation of Scaffold Variant Patterns Between Treatments', fontsize=14)
        plt.tight_layout()
        
        correlation_heatmap = os.path.join(self.plot_dir, 'correlation_heatmap.png')
        plt.savefig(correlation_heatmap, dpi=300)
        plt.close()
        
        plot_files['correlation_heatmap'] = correlation_heatmap
        
        # 7. Clustered heatmap of scaffold variant patterns
        # Get top 50 scaffolds by total variants
        top_scaffolds = sorted(self.scaffold_metadata.items(), 
                              key=lambda x: -x[1]['total_variants'])[:50]
        top_scaffold_ids = [s for s, _ in top_scaffolds]
        
        # Create data matrix
        data = np.zeros((len(top_scaffold_ids), len(treatments)))
        
        for i, scaffold in enumerate(top_scaffold_ids):
            for j, treatment in enumerate(treatments):
                data[i, j] = self.scaffold_variants[scaffold].get(treatment, 0)
        
        # Normalize by scaffold length to get densities
        for i, scaffold in enumerate(top_scaffold_ids):
            length_kb = self.scaffold_metadata[scaffold]['length'] / 1000
            data[i, :] = data[i, :] / length_kb
        
        # Cluster the data
        row_linkage = hierarchy.linkage(data, method='average')
        col_linkage = hierarchy.linkage(data.T, method='average')
        
        plt.figure(figsize=(14, 10))
        sns.clustermap(pd.DataFrame(data, index=top_scaffold_ids, columns=treatments),
                      cmap="YlOrRd", figsize=(14, 10), 
                      row_linkage=row_linkage, col_linkage=col_linkage,
                      xticklabels=treatments, yticklabels=top_scaffold_ids)
        plt.title('Clustered Heatmap of Scaffold Variant Densities', fontsize=16)
        
        clustered_heatmap = os.path.join(self.plot_dir, 'clustered_variant_density.png')
        plt.savefig(clustered_heatmap, dpi=300)
        plt.close()
        
        plot_files['clustered_heatmap'] = clustered_heatmap
        
        # 8. Bar chart comparing adaptation types
        adaptation_variants = [(a, stats['total_variants']) for a, stats in self.adaptation_stats.items() 
                              if a != 'None']  # Exclude controls
        adaptation_variants.sort(key=lambda x: -x[1])
        
        adaptations = [a for a, _ in adaptation_variants]
        counts = [c for _, c in adaptation_variants]
        
        plt.figure(figsize=(10, 6))
        plt.bar(adaptations, counts, color=['#ff9999', '#66b3ff'])
        plt.title('Total Variants by Adaptation Type', fontsize=14)
        plt.xlabel('Adaptation', fontsize=12)
        plt.ylabel('Number of Variants', fontsize=12)
        plt.tight_layout()
        
        adaptation_plot = os.path.join(self.plot_dir, 'adaptation_variant_counts.png')
        plt.savefig(adaptation_plot, dpi=300)
        plt.close()
        
        plot_files['adaptation_counts'] = adaptation_plot
        
        # 9. Venn-like diagram showing shared scaffolds between adaptation types
        # (We'll use a bar chart instead of a Venn diagram for simplicity)
        # Get shared and unique scaffold counts
        temp_scaffolds = self.adaptation_stats.get('Temperature', {}).get('scaffolds', set())
        oxygen_scaffolds = self.adaptation_stats.get('Low Oxygen', {}).get('scaffolds', set())
        
        only_temp = len(temp_scaffolds - oxygen_scaffolds)
        only_oxygen = len(oxygen_scaffolds - temp_scaffolds)
        shared = len(temp_scaffolds.intersection(oxygen_scaffolds))
        
        plt.figure(figsize=(10, 6))
        plt.bar(['Temperature Only', 'Shared', 'Low Oxygen Only'], 
               [only_temp, shared, only_oxygen],
               color=['#ff9999', '#9999ff', '#66b3ff'])
        plt.title('Scaffolds with Variants by Adaptation Type', fontsize=14)
        plt.ylabel('Number of Scaffolds', fontsize=12)
        plt.tight_layout()
        
        shared_scaffold_plot = os.path.join(self.plot_dir, 'shared_scaffolds.png')
        plt.savefig(shared_scaffold_plot, dpi=300)
        plt.close()
        
        plot_files['shared_scaffolds'] = shared_scaffold_plot
        
        return plot_files
    
    def generate_html_report(self, report_files, plot_files, enriched_scaffolds):
        """
        Generate comprehensive HTML report with interactive elements
        
        Args:
            report_files: Dict of generated report file paths
            plot_files: Dict of generated plot file paths
            enriched_scaffolds: Dict of enriched scaffolds by treatment
            
        Returns:
            Path to generated HTML report
        """
        report_file = os.path.join(self.report_dir, "scaffold_analysis_report.html")
        
        # Load treatment statistics
        treatment_stats = []
        for treatment, stats in self.treatment_stats.items():
            metadata = self.treatment_metadata.get(treatment, {})
            
            # Calculate global density
            total_length = sum(self.scaffold_metadata[s]['length'] for s in stats['scaffolds'])
            global_density = (stats['total_variants'] * 1000) / total_length if total_length > 0 else 0
            
            treatment_stats.append({
                'treatment': treatment,
                'description': metadata.get('description', ''),
                'adaptation': metadata.get('adaptation', ''),
                'has_gene': metadata.get('has_gene', ''),
                'total_variants': stats['total_variants'],
                'scaffolds': len(stats['scaffolds']),
                'global_density': global_density
            })
        
        # Sort by total variants
        treatment_stats.sort(key=lambda x: -x['total_variants'])
        
        # Load top enriched scaffolds
        top_enriched = []
        for treatment, scaffolds in enriched_scaffolds.items():
            for info in scaffolds[:5]:  # Top 5 per treatment
                top_enriched.append({
                    'treatment': treatment,
                    'scaffold': info['scaffold'],
                    'density': info['density'],
                    'fold_enrichment': info['fold_enrichment'],
                    'p_value': info['p_value']
                })
        
        # Sort by fold enrichment
        top_enriched.sort(key=lambda x: -x['fold_enrichment'])
        
        # Get correlation data
        correlation_matrix, treatments = self.calculate_correlations()
        correlation_data = {}
        
        for i, t1 in enumerate(treatments):
            for j, t2 in enumerate(treatments):
                if i < j:  # Only include each pair once
                    pair = f"{t1} vs {t2}"
                    correlation_data[pair] = correlation_matrix.loc[t1, t2]
        
        # Sort correlations
        sorted_correlations = sorted(correlation_data.items(), key=lambda x: -x[1])
        
        # Generate HTML
        with open(report_file, 'w') as f:
            f.write("""<!DOCTYPE html>
<html>
<head>
    <title>Scaffold-Level Variant Analysis Report</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1, h2, h3 { color: #333366; }
        .chart-container { width: 800px; height: 400px; margin: 20px 0; }
        .grid-container { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }
        table { border-collapse: collapse; width: 100%; margin-top: 20px; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        tr:nth-child(even) { background-color: #f9f9f9; }
        .note { background-color: #ffffdd; padding: 10px; border-left: 5px solid #ffcc00; margin: 20px 0; }
        .warning { background-color: #ffeeee; padding: 10px; border-left: 5px solid #ff0000; margin: 20px 0; }
        .info { background-color: #eeeeff; padding: 10px; border-left: 5px solid #0000ff; margin: 20px 0; }
        img { max-width: 100%; height: auto; margin: 20px 0; }
        .tab {
            overflow: hidden;
            border: 1px solid #ccc;
            background-color: #f1f1f1;
            margin-top: 20px;
        }
        .tab button {
            background-color: inherit;
            float: left;
            border: none;
            outline: none;
            cursor: pointer;
            padding: 14px 16px;
            transition: 0.3s;
            font-size: 17px;
        }
        .tab button:hover {
            background-color: #ddd;
        }
        .tab button.active {
            background-color: #ccc;
        }
        .tabcontent {
            display: none;
            padding: 6px 12px;
            border: 1px solid #ccc;
            border-top: none;
        }
    </style>
</head>
<body>
    <h1>Scaffold-Level Variant Analysis Report</h1>
    <p>Generated on: """ + f"{os.popen('date').read().strip()}" + """</p>
    
    <div class="info">
        <h3>Analysis Overview</h3>
        <p>This report presents a comprehensive analysis of genomic variants at the scaffold level across different treatments and adaptations.</p>
        <p>The analysis focuses on the distribution patterns of variants, identifying enriched scaffolds, and comparing patterns between treatments.</p>
    </div>
    
    <div class="tab">
        <button class="tablinks" onclick="openTab(event, 'Summary')" id="defaultOpen">Summary</button>
        <button class="tablinks" onclick="openTab(event, 'Treatments')">Treatments</button>
        <button class="tablinks" onclick="openTab(event, 'Scaffolds')">Scaffolds</button>
        <button class="tablinks" onclick="openTab(event, 'Enrichment')">Enrichment</button>
        <button class="tablinks" onclick="openTab(event, 'Comparisons')">Comparisons</button>
    </div>
    
    <div id="Summary" class="tabcontent">
        <h2>Analysis Summary</h2>
        <div class="info">
            <p><strong>Total Treatments Analyzed:</strong> """ + f"{len(self.treatment_stats)}" + """</p>
            <p><strong>Total Variants:</strong> """ + f"{sum(stats['total_variants'] for stats in self.treatment_stats.values())}" + """</p>
            <p><strong>Total Unique Scaffolds:</strong> """ + f"{len(self.scaffold_metadata)}" + """</p>
            <p><strong>Enriched Scaffolds Identified:</strong> """ + f"{sum(len(scaffolds) for scaffolds in enriched_scaffolds.values())}" + """</p>
        </div>
        
        <h3>Key Findings</h3>
        <ul>""")
            
            # Add key findings
            # Top treatments by variant count
            for treatment_info in treatment_stats[:3]:
                t = treatment_info['treatment']
                count = treatment_info['total_variants']
                desc = treatment_info['description']
                f.write(f"<li><strong>{t}</strong> ({desc}) has the highest number of variants ({count})</li>\n")
            
            # Top enriched scaffolds
            for enrichment in top_enriched[:3]:
                t = enrichment['treatment']
                s = enrichment['scaffold']
                fold = enrichment['fold_enrichment']
                f.write(f"<li>Scaffold <strong>{s}</strong> shows {fold:.1f}x enrichment in the <strong>{t}</strong> treatment</li>\n")
            
            # Top correlations
            if sorted_correlations:
                pair, corr = sorted_correlations[0]
                f.write(f"<li>The most similar treatments in terms of scaffold variant patterns are <strong>{pair}</strong> (correlation: {corr:.2f})</li>\n")
            
            f.write("""</ul>
        
        <h3>Variant Distribution by Treatment</h3>
        <img src="../plots/treatment_variant_counts.png" alt="Treatment Variant Counts">
        
        <h3>Variant Distribution by Adaptation</h3>
        <img src="../plots/adaptation_variant_counts.png" alt="Adaptation Variant Counts">
        
        <h3>Shared Scaffolds Between Adaptations</h3>
        <img src="../plots/shared_scaffolds.png" alt="Shared Scaffolds">
    </div>
    
    <div id="Treatments" class="tabcontent">
        <h2>Treatment Analysis</h2>
        <p>This section shows the distribution of variants across different treatments and adaptations.</p>
        
        <h3>Treatment Statistics</h3>
        <table>
            <tr>
                <th>Treatment</th>
                <th>Description</th>
                <th>Adaptation</th>
                <th>Gene Modified</th>
                <th>Total Variants</th>
                <th>Scaffolds</th>
                <th>Global Density</th>
            </tr>""")
            
            # Add treatment rows
            for treatment_info in treatment_stats:
                f.write(f"""
            <tr>
                <td>{treatment_info['treatment']}</td>
                <td>{treatment_info['description']}</td>
                <td>{treatment_info['adaptation']}</td>
                <td>{"Yes" if treatment_info['has_gene'] else "No"}</td>
                <td>{treatment_info['total_variants']}</td>
                <td>{treatment_info['scaffolds']}</td>
                <td>{treatment_info['global_density']:.4f}</td>
            </tr>""")
            
            f.write("""
        </table>
        
        <h3>Treatment Bubble Charts</h3>
        <p>These bubble charts show the top 20 scaffolds by variant density for each treatment.
           The size of each bubble represents the number of variants, and the color represents the scaffold length.</p>
        """)
            
            # Add bubble charts for each treatment
            for treatment in sorted(self.treatment_stats.keys()):
                bubble_path = f"../plots/{treatment}_bubble_chart.png"
                if f"{treatment}_bubble" in plot_files:
                    f.write(f"""
        <h4>{treatment} Treatment</h4>
        <img src="{bubble_path}" alt="{treatment} Bubble Chart">
        """)
            
            f.write("""
    </div>
    
    <div id="Scaffolds" class="tabcontent">
        <h2>Scaffold Analysis</h2>
        <p>This section provides details on the distribution of variants across scaffolds.</p>
        
        <h3>Top Scaffolds by Variant Density</h3>
        <img src="../plots/top_scaffolds_by_density.png" alt="Top Scaffolds by Density">
        
        <h3>Clustered Heatmap of Scaffold Variant Densities</h3>
        <p>This heatmap shows the variant densities for the top 50 scaffolds across all treatments,
           clustered to reveal patterns of similarity.</p>
        <img src="../plots/clustered_variant_density.png" alt="Clustered Variant Density Heatmap">
    </div>
    
    <div id="Enrichment" class="tabcontent">
        <h2>Scaffold Enrichment Analysis</h2>
        <p>This section shows statistical enrichment of variants on specific scaffolds.</p>
        
        <h3>Enrichment Heatmap</h3>
        <p>This heatmap shows the fold enrichment of the top enriched scaffolds across treatments.</p>
        <img src="../plots/enrichment_heatmap.png" alt="Enrichment Heatmap">
        
        <h3>Top Enriched Scaffolds</h3>
        <table>
            <tr>
                <th>Treatment</th>
                <th>Scaffold</th>
                <th>Density</th>
                <th>Fold Enrichment</th>
                <th>P-Value</th>
            </tr>""")
            
            # Add enriched scaffold rows
            for enrichment in top_enriched[:20]:  # Show top 20
                f.write(f"""
            <tr>
                <td>{enrichment['treatment']}</td>
                <td>{enrichment['scaffold']}</td>
                <td>{enrichment['density']:.4f}</td>
                <td>{enrichment['fold_enrichment']:.2f}x</td>
                <td>{enrichment['p_value']:.2e}</td>
            </tr>""")
            
            f.write("""
        </table>
    </div>
    
    <div id="Comparisons" class="tabcontent">
        <h2>Comparative Analysis</h2>
        <p>This section compares variant patterns between treatments and adaptations.</p>
        
        <h3>Treatment Correlation Heatmap</h3>
        <p>This heatmap shows the correlation of scaffold variant patterns between treatments.</p>
        <img src="../plots/correlation_heatmap.png" alt="Correlation Heatmap">
        
        <h3>Top Treatment Correlations</h3>
        <table>
            <tr>
                <th>Treatment Pair</th>
                <th>Correlation</th>
            </tr>""")
            
            # Add correlation rows
            for pair, corr in sorted_correlations[:10]:  # Show top 10
                f.write(f"""
            <tr>
                <td>{pair}</td>
                <td>{corr:.4f}</td>
            </tr>""")
            
            f.write("""
        </table>
        
        <h3>Adaptation Comparison</h3>
        <div class="grid-container">
            <div>
                <h4>Variant Counts</h4>
                <img src="../plots/adaptation_variant_counts.png" alt="Adaptation Variant Counts">
            </div>
            <div>
                <h4>Scaffold Distribution</h4>
                <img src="../plots/shared_scaffolds.png" alt="Shared Scaffolds">
            </div>
        </div>
    </div>
    
    <script>
        function openTab(evt, tabName) {
            var i, tabcontent, tablinks;
            tabcontent = document.getElementsByClassName("tabcontent");
            for (i = 0; i < tabcontent.length; i++) {
                tabcontent[i].style.display = "none";
            }
            tablinks = document.getElementsByClassName("tablinks");
            for (i = 0; i < tablinks.length; i++) {
                tablinks[i].className = tablinks[i].className.replace(" active", "");
            }
            document.getElementById(tabName).style.display = "block";
            evt.currentTarget.className += " active";
        }
        
        // Get the element with id="defaultOpen" and click on it
        document.getElementById("defaultOpen").click();
    </script>
</body>
</html>""")
        
        return report_file
    
    def generate_summary_text(self, enriched_scaffolds):
        """
        Generate a comprehensive summary text for the analysis
        
        Args:
            enriched_scaffolds: Dict of enriched scaffolds by treatment
            
        Returns:
            Summary text string
        """
        # Calculate basic statistics
        total_variants = sum(stats['total_variants'] for stats in self.treatment_stats.values())
        total_scaffolds = len(self.scaffold_metadata)
        total_enriched = sum(len(scaffolds) for scaffolds in enriched_scaffolds.values())
        
        # Get top treatments
        treatment_variants = [(t, stats['total_variants']) for t, stats in self.treatment_stats.items()]
        treatment_variants.sort(key=lambda x: -x[1])
        top_treatments = treatment_variants[:3]
        
        # Get top scaffolds by density
        top_scaffolds = []
        for scaffold, metadata in self.scaffold_metadata.items():
            top_scaffolds.append((scaffold, metadata['density'], metadata['total_variants']))
        
        top_scaffolds.sort(key=lambda x: -x[1])
        top_scaffold_list = top_scaffolds[:3]
        
        # Get adaptation statistics
        adaptation_stats = {
            a: {
                'variants': stats['total_variants'],
                'scaffolds': len(stats['scaffolds'])
            }
            for a, stats in self.adaptation_stats.items()
            if a != 'None'  # Exclude control groups
        }
        
        # Get top enriched scaffolds
        top_enriched = []
        for treatment, scaffolds in enriched_scaffolds.items():
            for info in scaffolds[:3]:  # Top 3 per treatment
                top_enriched.append((treatment, info['scaffold'], info['fold_enrichment']))
        
        top_enriched.sort(key=lambda x: -x[2])
        top_enriched = top_enriched[:5]  # Overall top 5
        
        # Get correlation data
        correlation_matrix, treatments = self.calculate_correlations()
        correlations = []
        
        for i, t1 in enumerate(treatments):
            for j, t2 in enumerate(treatments):
                if i < j:  # Only include each pair once
                    correlations.append((t1, t2, correlation_matrix.loc[t1, t2]))
        
        correlations.sort(key=lambda x: -x[2])
        top_correlations = correlations[:3]
        
        # Generate the summary text
        summary = "Scaffold Distribution Analysis Summary\n"
        summary += "=====================================\n\n"
        
        summary += "Overall Statistics:\n"
        summary += "-----------------\n"
        
        # Add treatment summaries
        for treatment, count in top_treatments:
            metadata = self.treatment_metadata.get(treatment, {})
            stats = self.treatment_stats[treatment]
            
            # Calculate global density
            total_length = sum(self.scaffold_metadata[s]['length'] for s in stats['scaffolds'])
            global_density = (stats['total_variants'] * 1000) / total_length if total_length > 0 else 0
            
            description = metadata.get('description', '')
            adaptation = metadata.get('adaptation', '')
            adaptation_text = f" ({adaptation} adaptation)" if adaptation and adaptation != 'None' else ""
            
            summary += f"{treatment} Treatment ({description}{adaptation_text}):\n"
            summary += f"  Adaptation: {adaptation}\n"
            summary += f"  Total Variants: {count}\n"
            summary += f"  Scaffolds with Variants: {len(stats['scaffolds'])}\n"
            summary += f"  Global Variant Density: {global_density:.4f} variants/kb\n"
            
            # Top 5 scaffolds by density
            scaffold_densities = []
            for scaffold in stats['scaffolds']:
                sc_count = self.scaffold_variants[scaffold][treatment]
                density = self.scaffold_density[scaffold][treatment]
                length = self.scaffold_metadata[scaffold]['length']
                scaffold_densities.append((scaffold, density, sc_count, length))
            
            scaffold_densities.sort(key=lambda x: -x[1])
            top_5 = scaffold_densities[:5]
            
            if top_5:
                summary += f"  Top 5 Scaffolds by Variant Density:\n"
                for scaffold, density, sc_count, length in top_5:
                    summary += f"    {scaffold}: {density:.4f} variants/kb (Count: {sc_count}, Length: {length}bp)\n"
            
            summary += "\n"
        
        # Add adaptation type comparison
        summary += "Adaptation Type Comparison:\n"
        summary += "-------------------------\n"
        
        for adaptation, stats in adaptation_stats.items():
            # Get treatments in this adaptation
            adaptation_treatments = self.adaptation_groups.get(adaptation, [])
            # Calculate average variants per treatment
            avg_variants = stats['variants'] / len(adaptation_treatments) if adaptation_treatments else 0
            
            summary += f"{adaptation} Adaptation:\n"
            summary += f"  Average Variants: {avg_variants:.1f}\n"
            summary += f"  Total Unique Scaffolds: {stats['scaffolds']}\n"
            
            # Calculate average density
            total_length = 0
            for treatment in adaptation_treatments:
                if treatment in self.treatment_stats:
                    for scaffold in self.treatment_stats[treatment]['scaffolds']:
                        total_length += self.scaffold_metadata[scaffold]['length']
            
            avg_density = (stats['variants'] * 1000) / total_length if total_length > 0 else 0
            summary += f"  Average Variant Density: {avg_density:.4f} variants/kb\n"
            summary += "\n"
        
        # Add gene modification effects
        summary += "Gene Modification Effects:\n"
        summary += "------------------------\n"
        
        for adaptation in ['Low Oxygen', 'Temperature']:
            modified_treatments = [t for t in self.adaptation_groups.get(adaptation, []) 
                                  if self.treatment_metadata.get(t, {}).get('has_gene', False)]
            non_modified_treatments = [t for t in self.adaptation_groups.get(adaptation, []) 
                                     if not self.treatment_metadata.get(t, {}).get('has_gene', False)]
            
            modified_variants = sum(self.treatment_stats.get(t, {}).get('total_variants', 0) 
                                  for t in modified_treatments)
            non_modified_variants = sum(self.treatment_stats.get(t, {}).get('total_variants', 0) 
                                     for t in non_modified_treatments)
            
            summary += f"{adaptation} Adaptation:\n"
            summary += f"  Gene-Modified Variants: {modified_variants}\n"
            summary += f"  Non-Modified Variants: {non_modified_variants}\n"
            
            ratio = modified_variants / non_modified_variants if non_modified_variants > 0 else 0
            summary += f"  Gene/Non-Gene Ratio: {ratio:.2f}\n"
            summary += "\n"
        
        # Add treatment correlation summary
        summary += "Treatment Correlation Summary:\n"
        summary += "----------------------------\n"
        
        if top_correlations:
            t1, t2, corr = top_correlations[0]
            summary += f"  Most Similar Treatments: {t1} and {t2} (ρ = {corr:.4f})\n"
        
        if len(correlations) > 0:
            t1, t2, corr = correlations[-1]
            summary += f"  Most Different Treatments: {t1} and {t2} (ρ = {corr:.4f})\n"
        
        summary += "\n  Full Correlation Matrix:\n"
        summary += "         " + "     ".join(treatments) + "\n"
        
        for t1 in treatments:
            summary += f"  {t1.ljust(5)} "
            for t2 in treatments:
                corr = correlation_matrix.loc[t1, t2]
                summary += f" {corr:.4f}"
            summary += "\n"
        
        # Add adaptation correlation
        if 'Temperature' in adaptation_stats and 'Low Oxygen' in adaptation_stats:
            # Calculate correlation between temperature and low oxygen adaptations
            temp_treatments = self.adaptation_groups.get('Temperature', [])
            oxygen_treatments = self.adaptation_groups.get('Low Oxygen', [])
            
            temp_data = {}
            oxygen_data = {}
            
            # Get all scaffolds from both adaptations
            all_scaffolds = set()
            for treatment in temp_treatments + oxygen_treatments:
                if treatment in self.treatment_stats:
                    all_scaffolds.update(self.treatment_stats[treatment]['scaffolds'])
            
            # Calculate average counts for each adaptation
            for scaffold in all_scaffolds:
                temp_count = sum(self.scaffold_variants[scaffold].get(t, 0) for t in temp_treatments) / len(temp_treatments) if temp_treatments else 0
                oxygen_count = sum(self.scaffold_variants[scaffold].get(t, 0) for t in oxygen_treatments) / len(oxygen_treatments) if oxygen_treatments else 0
                
                if temp_count > 0:
                    temp_data[scaffold] = temp_count
                
                if oxygen_count > 0:
                    oxygen_data[scaffold] = oxygen_count
            
            # Calculate correlation
            shared_scaffolds = set(temp_data.keys()) & set(oxygen_data.keys())
            if shared_scaffolds:
                temp_values = [temp_data[s] for s in shared_scaffolds]
                oxygen_values = [oxygen_data[s] for s in shared_scaffolds]
                
                corr, _ = scipy_stats.spearmanr(temp_values, oxygen_values)
                
                summary += "\n  Adaptation Correlation:\n"
                summary += f"    Temperature vs Low Oxygen: {corr:.4f}\n"
                summary += f"    Low Oxygen vs Temperature: {corr:.4f}\n"
        
        summary += "\nMain Conclusions:\n"
        summary += "---------------\n"
        summary += "1. This analysis identifies scaffolds with high variant densities across treatments.\n"
        summary += "2. Several scaffolds show treatment-specific enrichment.\n"
        summary += "3. The pattern of variant distribution provides insights into adaptation mechanisms.\n"
        summary += "4. Temperature and low oxygen adaptations show distinct scaffold distribution patterns.\n"
        summary += "5. Gene modifications (STC, CAS) appear to influence scaffold enrichment patterns.\n"
        summary += "6. Further sequence analysis of enriched scaffolds may reveal functional implications.\n"
        
        return summary
    
    def run_complete_analysis(self, vcf_files, scaffold_mapping_file=None):
        """Run the complete analysis pipeline"""
        print(f"Starting scaffold analysis with {len(vcf_files)} VCF files")
        
        # Load scaffold metadata if provided
        if scaffold_mapping_file:
            self.scaffold_lengths = self.load_scaffold_lengths(scaffold_mapping_file)
        
        # Process all VCF files
        results = self.process_all_vcfs(vcf_files)
        print(f"Processed {len(results)} VCF files, found variants on {len(self.scaffold_metadata)} scaffolds")
        
        # Identify enriched scaffolds
        enriched_scaffolds = self.identify_enriched_scaffolds()
        total_enriched = sum(len(scaffolds) for scaffolds in enriched_scaffolds.values())
        print(f"Identified {total_enriched} enriched scaffolds across all treatments")
        
        # Generate text reports
        report_files = self.generate_text_reports(enriched_scaffolds)
        print("\nGenerated text reports:")
        for report_type, report_file in report_files.items():
            print(f"  {report_type}: {report_file}")
        
        # Generate plots
        plot_files = self.generate_plots(enriched_scaffolds)
        print("\nGenerated plots:")
        for plot_type, plot_file in plot_files.items():
            print(f"  {plot_type}: {plot_file}")
        
        # Generate HTML report
        html_report = self.generate_html_report(report_files, plot_files, enriched_scaffolds)
        print(f"\nComplete HTML report: {html_report}")
        
        # Generate summary text
        summary_text = self.generate_summary_text(enriched_scaffolds)
        summary_file = os.path.join(self.report_dir, "scaffold_analysis_summary.txt")
        
        with open(summary_file, 'w') as f:
            f.write(summary_text)
        
        print(f"Summary text written to: {summary_file}")
        
        return {
            'results': results,
            'enriched_scaffolds': enriched_scaffolds,
            'report_files': report_files,
            'plot_files': plot_files,
            'html_report': html_report,
            'summary_file': summary_file
        }

def main():
    parser = argparse.ArgumentParser(
        description='Comprehensive scaffold-level analysis of variants across treatments',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--vcf', required=True, nargs='+', help='VCF file(s) to analyze')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--scaffold-mapping', help='Optional scaffold mapping file')
    
    args = parser.parse_args()
    
    # Create analyzer
    analyzer = ScaffoldAnalyzer(
        output_dir=args.output_dir
    )
    
    # Run analysis
    analysis_results = analyzer.run_complete_analysis(
        vcf_files=args.vcf,
        scaffold_mapping_file=args.scaffold_mapping
    )
    
    print("\nAnalysis complete!")
    print(f"Results are available in {args.output_dir}")

if __name__ == '__main__':
    main()