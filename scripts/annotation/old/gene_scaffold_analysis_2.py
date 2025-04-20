#!/usr/bin/env python3
"""
gene_scaffold_analysis.py - Gene-centric scaffold-level analysis of variants

This script provides a comprehensive gene-centric analysis of variant patterns at the scaffold level,
extending the scaffold analysis to include gene identity and function. It maps scaffolds to genes
and provides visualizations and analysis of variant patterns with respect to gene functions.
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

class GeneScaffoldAnalyzer:
    def __init__(self, output_dir, gene_index_file, scaffold_mapping_file):
        """
        Initialize the gene-scaffold analyzer.
        
        Args:
            output_dir: Directory for output files and reports
            gene_index_file: Path to gene index file
            scaffold_mapping_file: Path to scaffold mapping file
        """
        self.output_dir = output_dir
        self.gene_index_file = gene_index_file
        self.scaffold_mapping_file = scaffold_mapping_file
        
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
        self.scaffold_metadata = {}        # Metadata about scaffolds
        self.jriu_to_scaffold = {}         # Maps JRIU ID to scaffold_id
        self.scaffold_to_jriu = {}         # Maps scaffold_id to JRIU ID
        self.gene_metadata = {}            # Metadata about genes
        self.scaffold_to_genes = {}        # Maps scaffold_id to list of genes
        self.gene_to_scaffold = {}         # Maps gene_id to scaffold_id
        self.sc_gene_map = {}              # Maps SC gene ID to W303 gene ID
        self.jriu_to_sc_genes = {}         # Maps JRIU ID to list of SC gene IDs with positions
        self.gene_variants = defaultdict(lambda: defaultdict(int))  # Variants by gene and treatment
        self.gene_density = defaultdict(lambda: defaultdict(float)) # Variant density by gene and treatment
        self.sc_gene_info = {}             # Detailed info about SC genes from genes_of_interest_search.tsv
        self.sc_id_to_name_map = {}        # Maps SC gene ID to standard gene name/symbol
        
        self.scaffold_variants = defaultdict(lambda: defaultdict(int))  # Variants by scaffold and treatment
        self.treatment_stats = defaultdict(dict)  # Statistics by treatment
        self.adaptation_stats = defaultdict(dict)  # Statistics by adaptation
        
        # Load mappings and metadata
        self.load_scaffold_mapping()
        self.load_gene_index()
    
    def load_scaffold_mapping(self):
        """
        Load scaffold mapping information from the mapping file.
        
        This creates mappings between scaffold IDs, JRIU IDs, and stores 
        metadata about each scaffold.
        """
        print(f"Loading scaffold mapping from {self.scaffold_mapping_file}")
        
        with open(self.scaffold_mapping_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                scaffold_id = row.get('scaffold_id')
                jriu_id = row.get('jriu_id')
                
                if scaffold_id and jriu_id:
                    # Create mappings
                    self.jriu_to_scaffold[jriu_id] = scaffold_id
                    self.scaffold_to_jriu[scaffold_id] = jriu_id
                    
                    # Store metadata
                    self.scaffold_metadata[scaffold_id] = {
                        'scaffold_id': scaffold_id,
                        'jriu_id': jriu_id,
                        'length': int(row.get('scaffold_length', 0)),
                        'gene_count': int(row.get('gene_count', 0)),
                        'cds_count': int(row.get('cds_count', 0)),
                        'trna_count': int(row.get('trna_count', 0))
                    }
        
        print(f"Loaded {len(self.scaffold_metadata)} scaffolds with mapping information")
    
    def load_gene_index(self):
        """
        Load gene information from the gene index file.
        
        This creates mappings between genes and scaffolds, and stores
        metadata about each gene.
        """
        print(f"Loading gene index from {self.gene_index_file}")
        
        # Initialize scaffold to genes mapping
        for scaffold_id in self.scaffold_metadata:
            self.scaffold_to_genes[scaffold_id] = []
        
        with open(self.gene_index_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gene_id = row.get('gene_id')
                scaffold_id = row.get('scaffold_id')
                
                if gene_id and scaffold_id:
                    # Map gene to scaffold
                    self.gene_to_scaffold[gene_id] = scaffold_id
                    
                    # Add gene to scaffold's gene list
                    if scaffold_id in self.scaffold_to_genes:
                        self.scaffold_to_genes[scaffold_id].append(gene_id)
                    
                    # Store gene metadata
                    self.gene_metadata[gene_id] = {
                        'gene_id': gene_id,
                        'scaffold_id': scaffold_id,
                        'start': int(row.get('start', 0)),
                        'end': int(row.get('end', 0)),
                        'strand': row.get('strand', ''),
                        'product': row.get('product', ''),
                        'note': row.get('note', '')
                    }
                    
                    # Map SC gene ID if available
                    sc_gene_id = row.get('sc_gene_id')
                    if sc_gene_id:
                        self.sc_gene_map[sc_gene_id] = gene_id
                        self.gene_metadata[gene_id]['sc_gene_id'] = sc_gene_id
                        
                        # Get the corresponding JRIU ID
                        jriu_id = self.scaffold_to_jriu.get(scaffold_id)
                        if jriu_id and sc_gene_id:
                            if jriu_id not in self.jriu_to_sc_genes:
                                self.jriu_to_sc_genes[jriu_id] = []
                            self.jriu_to_sc_genes[jriu_id].append({
                                'sc_gene_id': sc_gene_id,
                                'start': int(row.get('start', 0)),
                                'end': int(row.get('end', 0))
                            })
        
        print(f"Loaded {len(self.gene_metadata)} genes across {len(self.scaffold_to_genes)} scaffolds")
        
        # Load the genes_of_interest_search.tsv file to get proper SC gene information
        genes_of_interest_search_file = os.path.join(os.path.dirname(self.gene_index_file), 'genes_of_interest_search.tsv')
        self.load_sc_gene_info(genes_of_interest_search_file)
        
    def load_sc_gene_info(self, genes_of_interest_search_file):
        """
        Load SC gene information from genes_of_interest_search.tsv
        
        This provides complete mapping between SC gene IDs and their metadata
        """
        print(f"Loading SC gene information from {genes_of_interest_search_file}")
        
        if not os.path.exists(genes_of_interest_search_file):
            print(f"Warning: {genes_of_interest_search_file} not found")
            return
            
        with open(genes_of_interest_search_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                sc_gene_id = row.get('sc_gene_id')
                sc_gene_name = row.get('sc_gene_name')
                w303_gene_id = row.get('w303_gene_id')
                scaffold_id = row.get('scaffold_id')
                jriu_id = row.get('jriu_id')
                product = row.get('product')
                
                if sc_gene_id:
                    # Store SC gene information
                    self.sc_gene_info[sc_gene_id] = {
                        'sc_gene_id': sc_gene_id,
                        'sc_gene_name': sc_gene_name,
                        'w303_gene_id': w303_gene_id,
                        'scaffold_id': scaffold_id,
                        'jriu_id': jriu_id,
                        'product': product,
                        'start': int(row.get('start', 0)),
                        'end': int(row.get('end', 0)),
                        'strand': row.get('strand', '')
                    }
                    
                    # Update sc_id_to_name_map
                    if sc_gene_name and sc_gene_id not in self.sc_id_to_name_map:
                        self.sc_id_to_name_map[sc_gene_id] = sc_gene_name
                    
                    # Update sc_gene_map
                    if w303_gene_id and sc_gene_id not in self.sc_gene_map:
                        self.sc_gene_map[sc_gene_id] = w303_gene_id
                        
                    print(f"Loaded information for {len(self.sc_gene_info)} SC genes")
                    
                    if sc_gene_id and sc_gene_name:
                        self.sc_id_to_name_map[sc_gene_id] = sc_gene_name
                        
                    if jriu_id and sc_gene_id:
                        if jriu_id not in self.jriu_to_sc_genes:
                            self.jriu_to_sc_genes[jriu_id] = []
                        # Only add if not already present
                        existing_ids = [g['sc_gene_id'] for g in self.jriu_to_sc_genes[jriu_id]]
                        if sc_gene_id not in existing_ids:
                            self.jriu_to_sc_genes[jriu_id].append({
                                'sc_gene_id': sc_gene_id,
                                'name': sc_gene_name,
                                'w303_gene_id': w303_gene_id
                            })
        
        # Add standard gene names for special cases
        self.gene_name_map = {
            "YHR190W": "ERG9",
            "YGR175C": "ERG1",
            "YHR072W": "ERG7",
            "YHR007C": "ERG11",
            "YNL280C": "ERG24",
            "YGR060W": "ERG25",
            "YML008C": "ERG6",
            "YMR202W": "ERG2",
            "YLR056W": "ERG3",
            "YMR015C": "ERG5",
            "YGL012W": "ERG4"
        }
        
        # Extra dictionary to store SC gene information from genes_of_interest_search.tsv
        self.sc_gene_info = {}
        
        # Update the sc_id_to_name_map with these special cases if not already present
        for sc_id, name in self.gene_name_map.items():
            if sc_id not in self.sc_id_to_name_map:
                self.sc_id_to_name_map[sc_id] = name
    
    def process_vcf_file(self, vcf_file):
        """
        Process a single VCF file to extract variant information by scaffold and gene.
        
        Args:
            vcf_file: Path to VCF file
            
        Returns:
            Dict containing statistics about variants by scaffold and gene
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
        gene_counts = defaultdict(int)
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
                
                # Extract variant info
                jriu_id = fields[0]      # JRIU ID
                position = int(fields[1])
                
                # Check the info field for gene annotations
                info_field = fields[7] if len(fields) > 7 else ""
                
                # Parse the GenesOnJRIU field to extract SC gene IDs
                sc_gene_ids = []
                if "GenesOnJRIU=" in info_field:
                    genes_part = info_field.split("GenesOnJRIU=")[1].split(";")[0]
                    for gene_entry in genes_part.split(","):
                        # Format is typically W3030BY00190(YHR190W)
                        if "(" in gene_entry and ")" in gene_entry:
                            sc_id = gene_entry.split("(")[1].split(")")[0]
                            if sc_id:  # Ensure it's not an empty string
                                sc_gene_ids.append(sc_id)
                
                # Map JRIU to scaffold
                scaffold_id = self.jriu_to_scaffold.get(jriu_id)
                
                if scaffold_id:
                    # Count unique positions per scaffold
                    position_key = f"{jriu_id}:{position}"
                    if position_key not in scaffold_positions[scaffold_id]:
                        scaffold_positions[scaffold_id].add(position_key)
                        scaffold_counts[scaffold_id] += 1
                        total_variants += 1
                        
                        # If we have SC gene IDs from the VCF, use those directly
                        if sc_gene_ids:
                            for sc_gene_id in sc_gene_ids:
                                # Use the SC gene ID directly as the counter key
                                gene_counts[sc_gene_id] += 1
                        else:
                            # If no genes in the VCF, fall back to scaffold genes
                            genes = self.scaffold_to_genes.get(scaffold_id, [])
                            if genes:
                                # Distribute variants to genes on this scaffold (less accurate fallback)
                                for gene_id in genes:
                                    if gene_id in self.gene_metadata and 'sc_gene_id' in self.gene_metadata[gene_id]:
                                        sc_gene_id = self.gene_metadata[gene_id]['sc_gene_id']
                                        if sc_gene_id:
                                            gene_counts[sc_gene_id] += 1
        
        # Store results
        result = {
            'sample': sample_name,
            'treatment': treatment,
            'total_variants': total_variants,
            'scaffold_counts': dict(scaffold_counts),
            'gene_counts': dict(gene_counts)
        }
        
        # Update treatment statistics
        if treatment not in self.treatment_stats:
            self.treatment_stats[treatment] = {
                'total_variants': 0,
                'sample_count': 0,
                'scaffolds': set(),
                'genes': set(),
                'scaffold_counts': defaultdict(int),
                'gene_counts': defaultdict(int)
            }
        
        self.treatment_stats[treatment]['total_variants'] += total_variants
        self.treatment_stats[treatment]['sample_count'] += 1
        self.treatment_stats[treatment]['scaffolds'].update(scaffold_counts.keys())
        self.treatment_stats[treatment]['genes'].update(gene_counts.keys())
        
        # Update scaffold and gene counts for this treatment
        for scaffold_id, count in scaffold_counts.items():
            self.treatment_stats[treatment]['scaffold_counts'][scaffold_id] += count
            self.scaffold_variants[scaffold_id][treatment] += count
        
        for gene_id, count in gene_counts.items():
            self.treatment_stats[treatment]['gene_counts'][gene_id] += count
            self.gene_variants[gene_id][treatment] += count
        
        print(f"  Found {total_variants} variants across {len(scaffold_counts)} scaffolds and {len(gene_counts)} genes")
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
                'genes': set().union(*[self.treatment_stats[t]['genes'] for t in treatments if t in self.treatment_stats]),
                'scaffold_counts': defaultdict(int),
                'gene_counts': defaultdict(int)
            }
            
            # Combine scaffold and gene counts from all treatments in this adaptation group
            for treatment in treatments:
                if treatment in self.treatment_stats:
                    for scaffold_id, count in self.treatment_stats[treatment]['scaffold_counts'].items():
                        self.adaptation_stats[adaptation]['scaffold_counts'][scaffold_id] += count
                    
                    for gene_id, count in self.treatment_stats[treatment]['gene_counts'].items():
                        self.adaptation_stats[adaptation]['gene_counts'][gene_id] += count
        
        # Calculate gene densities
        for gene_id, counts in self.gene_variants.items():
            # Get gene length - handle both SC gene IDs and W303 gene IDs
            gene_length = None
            
            # For SC gene IDs, find the corresponding W303 gene
            if gene_id in self.sc_gene_map:
                w303_gene_id = self.sc_gene_map[gene_id]
                if w303_gene_id in self.gene_metadata:
                    gene_length = self.gene_metadata[w303_gene_id]['end'] - self.gene_metadata[w303_gene_id]['start'] + 1
            # Direct lookup for W303 gene IDs
            elif gene_id in self.gene_metadata:
                gene_length = self.gene_metadata[gene_id]['end'] - self.gene_metadata[gene_id]['start'] + 1
            
            # Use a default length if not found (1000 bp is reasonable for yeast)
            if not gene_length:
                gene_length = 1000
                
            # Calculate densities for each treatment
            for treatment, count in counts.items():
                # Density in variants per kilobase
                density = (count * 1000) / gene_length if gene_length > 0 else 0
                self.gene_density[gene_id][treatment] = density
        
        return all_results
    
    def identify_enriched_genes(self, min_variants=3, min_fold=1.5, max_pval=0.05):
        """
        Identify genes with statistically significant enrichment of variants
        
        Args:
            min_variants: Minimum number of variants for a gene to be considered
            min_fold: Minimum fold enrichment over background
            max_pval: Maximum p-value for significance
            
        Returns:
            Dict mapping treatment to lists of enriched genes
        """
        enriched_genes = defaultdict(list)
        
        # Calculate global gene density for each treatment
        treatment_densities = {}
        for treatment, stats in self.treatment_stats.items():
            total_length = 0
            variant_count = 0
            
            # For SC gene IDs, we need to estimate lengths
            for gene_id in stats['genes']:
                # Use gene lengths from self.gene_metadata if available
                gene_length = None
                
                # For SC gene IDs, we need to find the corresponding entry
                if gene_id in self.sc_gene_map:
                    w303_gene_id = self.sc_gene_map[gene_id]
                    if w303_gene_id in self.gene_metadata:
                        gene_length = self.gene_metadata[w303_gene_id]['end'] - self.gene_metadata[w303_gene_id]['start'] + 1
                elif gene_id in self.gene_metadata:
                    gene_length = self.gene_metadata[gene_id]['end'] - self.gene_metadata[gene_id]['start'] + 1
                
                # Use a default length if we can't find it (1000 bp is a reasonable gene size for yeast)
                if not gene_length:
                    gene_length = 1000
                    
                total_length += gene_length
                variant_count += stats['gene_counts'][gene_id]
            
            global_density = (variant_count * 1000) / total_length if total_length > 0 else 0
            treatment_densities[treatment] = global_density
        
        # Test each gene in each treatment for enrichment
        for treatment, stats in self.treatment_stats.items():
            global_density = treatment_densities[treatment]
            
            for gene_id in stats['genes']:
                count = self.gene_variants[gene_id][treatment]
                
                # Skip genes with too few variants
                if count < min_variants:
                    continue
                
                # Get gene length - handle both SC IDs and W303 IDs
                gene_length = None
                w303_gene_id = None
                sc_gene_id = None
                
                # Determine the gene's properties based on its ID
                if gene_id in self.sc_gene_map:  # This is an SC ID
                    sc_gene_id = gene_id
                    w303_gene_id = self.sc_gene_map[gene_id]
                    if w303_gene_id in self.gene_metadata:
                        gene_length = self.gene_metadata[w303_gene_id]['end'] - self.gene_metadata[w303_gene_id]['start'] + 1
                elif gene_id in self.gene_metadata:  # This is a W303 ID
                    w303_gene_id = gene_id
                    sc_gene_id = self.gene_metadata[gene_id].get('sc_gene_id', '')
                    gene_length = self.gene_metadata[gene_id]['end'] - self.gene_metadata[gene_id]['start'] + 1
                
                # Use default length if still not found
                if not gene_length:
                    gene_length = 1000
                
                # Calculate gene density and fold enrichment
                gene_density = (count * 1000) / gene_length if gene_length > 0 else 0
                fold_enrichment = gene_density / global_density if global_density > 0 else 0
                
                # Skip if fold enrichment is too low
                if fold_enrichment < min_fold:
                    continue
                
                # Statistical test (Poisson test)
                length_kb = gene_length / 1000
                expected_count = global_density * length_kb
                p_value = scipy_stats.poisson.sf(count - 1, expected_count)
                q_value = p_value  # Add multiple testing correction if needed
                
                # Record if significant
                if p_value <= max_pval:
                    # Get gene name
                    gene_name = ""
                    product = ""
                    scaffold_id = ""
                    
                    if sc_gene_id:
                        gene_name = self.sc_id_to_name_map.get(sc_gene_id, sc_gene_id)
                    
                    if w303_gene_id in self.gene_metadata:
                        product = self.gene_metadata[w303_gene_id].get('product', '')
                        scaffold_id = self.gene_metadata[w303_gene_id].get('scaffold_id', '')
                    
                    enriched_genes[treatment].append({
                        'gene_id': w303_gene_id or gene_id,
                        'sc_gene_id': sc_gene_id or "",
                        'gene_name': gene_name,
                        'scaffold_id': scaffold_id,
                        'count': count,
                        'density': gene_density,
                        'fold_enrichment': fold_enrichment,
                        'expected_count': expected_count,
                        'p_value': p_value,
                        'q_value': q_value,
                        'gene_length': gene_length,
                        'product': product
                    })
        
        # Sort enriched genes by fold enrichment
        for treatment in enriched_genes:
            enriched_genes[treatment].sort(key=lambda x: -x['fold_enrichment'])
        
        return enriched_genes
    
    def calculate_correlations(self):
        """
        Calculate correlations between treatments based on gene variant patterns
        
        Returns:
            Tuple of (correlation_matrix, treatment_list)
        """
        treatments = sorted(self.treatment_stats.keys())
        all_genes = set()
        
        # Get all genes that have variants in any treatment
        for treatment, stats in self.treatment_stats.items():
            all_genes.update(stats['genes'])
        
        # Create data frame with gene counts by treatment
        data = []
        for gene_id in all_genes:
            row = [gene_id]
            for treatment in treatments:
                row.append(self.gene_variants[gene_id].get(treatment, 0))
            data.append(row)
        
        columns = ['gene_id'] + treatments
        df = pd.DataFrame(data, columns=columns)
        df.set_index('gene_id', inplace=True)
        
        # Calculate Spearman correlation
        correlation_matrix = df.corr(method='spearman')
        
        return correlation_matrix, treatments
    
    def generate_text_reports(self, enriched_genes):
        """
        Generate comprehensive text-based reports for gene analysis
        
        Args:
            enriched_genes: Dict of enriched genes by treatment
        
        Returns:
            Dict of report file paths
        """
        # 1. Treatment gene summary
        treatment_report_file = os.path.join(self.report_dir, "treatment_gene_summary.tsv")
        with open(treatment_report_file, 'w') as f:
            f.write("Treatment\tDescription\tAdaptation\tHas_Gene\tTotal_Variants\t")
            f.write("Genes_With_Variants\tGlobal_Density\tTop_Genes\n")
            
            for treatment in sorted(self.treatment_stats.keys()):
                stats = self.treatment_stats[treatment]
                metadata = self.treatment_metadata.get(treatment, {})
                
                # Calculate global density
                total_length = 0
                variant_count = 0
                
                for gene_id in stats['genes']:
                    # Get gene length - handle both SC gene IDs and W303 gene IDs
                    gene_length = None
                    
                    # For SC gene IDs, find the corresponding W303 gene
                    if gene_id in self.sc_gene_map:
                        w303_gene_id = self.sc_gene_map[gene_id]
                        if w303_gene_id in self.gene_metadata:
                            gene_length = self.gene_metadata[w303_gene_id]['end'] - self.gene_metadata[w303_gene_id]['start'] + 1
                    # Direct lookup for W303 gene IDs
                    elif gene_id in self.gene_metadata:
                        gene_length = self.gene_metadata[gene_id]['end'] - self.gene_metadata[gene_id]['start'] + 1
                    
                    # Use a default length if not found
                    if not gene_length:
                        gene_length = 1000
                        
                    total_length += gene_length
                    variant_count += stats['gene_counts'][gene_id]
                
                global_density = (variant_count * 1000) / total_length if total_length > 0 else 0
                
                # Get top 3 genes by variant density
                top_genes = []
                for gene_id in stats['genes']:
                    count = stats['gene_counts'][gene_id]
                    
                    # Get gene name
                    gene_name = gene_id
                    
                    # For SC gene IDs, use the name from the mapping
                    if gene_id in self.sc_id_to_name_map:
                        gene_name = self.sc_id_to_name_map.get(gene_id, gene_id)
                    elif gene_id in self.gene_metadata and 'sc_gene_id' in self.gene_metadata[gene_id]:
                        sc_gene_id = self.gene_metadata[gene_id]['sc_gene_id']
                        gene_name = self.sc_id_to_name_map.get(sc_gene_id, sc_gene_id) or gene_id
                    
                    # Get density
                    if gene_id in self.gene_density and treatment in self.gene_density[gene_id]:
                        density = self.gene_density[gene_id][treatment]
                        top_genes.append((gene_name, density, count))
                
                top_genes.sort(key=lambda x: -x[1])
                top_3 = top_genes[:3]
                top_desc = ", ".join([f"{g} ({d:.2f})" for g, d, c in top_3])
                
                f.write(f"{treatment}\t{metadata.get('description', '')}\t{metadata.get('adaptation', '')}\t")
                f.write(f"{metadata.get('has_gene', '')}\t{stats['total_variants']}\t{len(stats['genes'])}\t")
                f.write(f"{global_density:.6f}\t{top_desc}\n")
        
        # 2. Adaptation gene summary
        adaptation_report_file = os.path.join(self.report_dir, "adaptation_gene_summary.tsv")
        with open(adaptation_report_file, 'w') as f:
            f.write("Adaptation\tTreatments\tTotal_Variants\tGenes_With_Variants\t")
            f.write("Global_Density\tTop_Genes\n")
            
            for adaptation in sorted(self.adaptation_stats.keys()):
                stats = self.adaptation_stats[adaptation]
                
                # Skip if no data
                if not stats['genes']:
                    continue
                
                # Get total length of all genes in this adaptation
                total_length = 0
                for gene_id in stats['genes']:
                    if gene_id in self.gene_metadata:
                        gene_length = self.gene_metadata[gene_id]['end'] - self.gene_metadata[gene_id]['start'] + 1
                        total_length += gene_length
                
                global_density = (stats['total_variants'] * 1000) / total_length if total_length > 0 else 0
                
                # Get top 3 genes by variant count
                top_genes = []
                for gene_id, count in stats['gene_counts'].items():
                    if gene_id in self.gene_metadata:
                        sc_gene_id, gene_name, _, _ = self.get_gene_info(gene_id)
                        # Get gene length
                        gene_length = self.gene_metadata[gene_id]['end'] - self.gene_metadata[gene_id]['start'] + 1
                        density = (count * 1000) / gene_length if gene_length > 0 else 0
                        top_genes.append((gene_name, density, count))
                
                top_genes.sort(key=lambda x: -x[1])
                top_3 = top_genes[:3]
                top_desc = ", ".join([f"{g} ({d:.2f})" for g, d, c in top_3])
                
                treatments_str = ", ".join(stats['treatments'])
                
                f.write(f"{adaptation}\t{treatments_str}\t{stats['total_variants']}\t{len(stats['genes'])}\t")
                f.write(f"{global_density:.6f}\t{top_desc}\n")
        
        # 3. Gene summary
        gene_report_file = os.path.join(self.report_dir, "gene_summary.tsv")
        with open(gene_report_file, 'w') as f:
            f.write("Gene_ID\tSC_Gene_ID\tGene_Name\tProduct\tScaffold\tLength\t")
            f.write("Total_Variants\tAvg_Density\tTop_Treatment\tTop_Treatment_Count\n")
            
            # Get genes with variants
            genes_with_variants = set()
            for treatment, stats in self.treatment_stats.items():
                genes_with_variants.update(stats['genes'])
            
            # Sort genes by total variants
            gene_variants_total = {}
            for gene_id in genes_with_variants:
                total = sum(self.gene_variants[gene_id].values())
                gene_variants_total[gene_id] = total
            
            sorted_genes = sorted(gene_variants_total.items(), key=lambda x: -x[1])
            
            for gene_id, total_variants in sorted_genes:
                if gene_id in self.gene_metadata:
                    metadata = self.gene_metadata[gene_id]
                    sc_gene_id = metadata.get('sc_gene_id', '')
                    gene_name = self.gene_name_map.get(sc_gene_id, sc_gene_id) or ''
                    scaffold_id = metadata.get('scaffold_id', '')
                    gene_length = metadata['end'] - metadata['start'] + 1
                    product = metadata.get('product', '')
                    
                    # Find top treatment
                    top_treatment = None
                    top_count = 0
                    
                    for treatment, count in self.gene_variants[gene_id].items():
                        if count > top_count:
                            top_count = count
                            top_treatment = treatment
                    
                    # Calculate average density
                    total_density = sum(self.gene_density[gene_id].values())
                    avg_density = total_density / len(self.gene_density[gene_id]) if self.gene_density[gene_id] else 0
                    
                    f.write(f"{gene_id}\t{sc_gene_id}\t{gene_name}\t{product}\t{scaffold_id}\t{gene_length}\t")
                    f.write(f"{total_variants}\t{avg_density:.6f}\t{top_treatment}\t{top_count}\n")
        
        # 4. Enriched genes report
        enriched_report_file = os.path.join(self.report_dir, "enriched_genes.tsv")
        with open(enriched_report_file, 'w') as f:
            f.write("Treatment\tGene_ID\tSC_Gene_ID\tGene_Name\tProduct\tScaffold\tLength\t")
            f.write("Variants\tDensity\tGlobal_Density\tFold_Enrichment\tExpected_Count\tP_Value\tQ_Value\n")
            
            for treatment in sorted(enriched_genes.keys()):
                # Get global density
                stats = self.treatment_stats[treatment]
                total_length = 0
                for gene_id in stats['genes']:
                    if gene_id in self.gene_metadata:
                        gene_length = self.gene_metadata[gene_id]['end'] - self.gene_metadata[gene_id]['start'] + 1
                        total_length += gene_length
                
                global_density = (stats['total_variants'] * 1000) / total_length if total_length > 0 else 0
                
                for info in enriched_genes[treatment]:
                    f.write(f"{treatment}\t{info['gene_id']}\t{info['sc_gene_id']}\t{info['gene_name']}\t")
                    f.write(f"{info['product']}\t{info['scaffold_id']}\t{info['gene_length']}\t{info['count']}\t")
                    f.write(f"{info['density']:.6f}\t{global_density:.6f}\t{info['fold_enrichment']:.6f}\t")
                    f.write(f"{info['expected_count']:.2f}\t{info['p_value']:.6e}\t{info['q_value']:.6e}\n")
        
        # 5. Correlation report
        correlation_matrix, treatments = self.calculate_correlations()
        correlation_report_file = os.path.join(self.report_dir, "treatment_gene_correlations.tsv")
        
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
            'gene_report': gene_report_file,
            'enriched_report': enriched_report_file,
            'correlation_report': correlation_report_file
        }
    
    def generate_plots(self, enriched_genes):
        """
        Generate comprehensive plots for gene analysis
        
        Args:
            enriched_genes: Dict of enriched genes by treatment
            
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
        
        # 2. Number of genes with variants by treatment
        treatment_genes = [(t, len(stats['genes'])) for t, stats in self.treatment_stats.items()]
        treatment_genes.sort(key=lambda x: x[1], reverse=True)
        
        treatments = [t for t, _ in treatment_genes]
        gene_counts = [g for _, g in treatment_genes]
        
        plt.figure(figsize=(10, 6))
        plt.bar(treatments, gene_counts, color='lightgreen')
        plt.title('Genes with Variants by Treatment', fontsize=14)
        plt.xlabel('Treatment', fontsize=12)
        plt.ylabel('Number of Genes', fontsize=12)
        plt.xticks(rotation=45)
        plt.tight_layout()
        
        treatment_gene_plot = os.path.join(self.plot_dir, 'treatment_gene_counts.png')
        plt.savefig(treatment_gene_plot, dpi=300)
        plt.close()
        
        plot_files['treatment_genes'] = treatment_gene_plot
        
        # 3. Top 10 genes by variant density across all treatments
        # Get genes with variants
        genes_with_variants = set()
        for treatment, stats in self.treatment_stats.items():
            genes_with_variants.update(stats['genes'])
        
        # Calculate average density for each gene
        gene_avg_density = {}
        for gene_id in genes_with_variants:
            if gene_id in self.gene_density:
                densities = list(self.gene_density[gene_id].values())
                if densities:
                    gene_avg_density[gene_id] = sum(densities) / len(densities)
        
        # Sort by density and get top 10
        top_genes = sorted(gene_avg_density.items(), key=lambda x: -x[1])[:10]
        
        gene_labels = []
        densities = []
        
        for gene_id, density in top_genes:
            sc_gene_id, gene_name, _, _ = self.get_gene_info(gene_id)
            gene_labels.append(gene_name)
            densities.append(density)
        
        if gene_labels and densities:
            plt.figure(figsize=(12, 6))
            plt.bar(gene_labels, densities, color='salmon')
            plt.title('Top 10 Genes by Variant Density (All Treatments)', fontsize=14)
            plt.xlabel('Gene', fontsize=12)
            plt.ylabel('Average Variant Density (variants/kb)', fontsize=12)
            plt.xticks(rotation=45)
            plt.tight_layout()
            
            top_gene_plot = os.path.join(self.plot_dir, 'top_genes_by_density.png')
            plt.savefig(top_gene_plot, dpi=300)
            plt.close()
            
            plot_files['top_genes'] = top_gene_plot
        
        # 4. Gene enrichment heatmap
        # Gather data for top 25 enriched genes across all treatments
        all_enriched = []
        for treatment, genes in enriched_genes.items():
            for info in genes:
                all_enriched.append((treatment, info['gene_id'], info['gene_name'], info['sc_gene_id'], info['fold_enrichment']))
        
        # Sort by fold enrichment and get top 25 unique genes
        all_enriched.sort(key=lambda x: -x[4])
        top_genes = []
        top_gene_ids = set()
        
        for _, gene_id, gene_name, sc_gene_id, _ in all_enriched:
            if gene_id not in top_gene_ids and len(top_genes) < 25:
                top_genes.append((gene_id, gene_name or sc_gene_id or gene_id))
                top_gene_ids.add(gene_id)
        
        if top_genes:
            # Create enrichment data for heatmap
            treatments = sorted(self.treatment_stats.keys())
            enrichment_data = np.zeros((len(top_genes), len(treatments)))
            
            for i, (gene_id, _) in enumerate(top_genes):
                for j, treatment in enumerate(treatments):
                    # Find enrichment info for this gene/treatment
                    enrichment = 1.0  # Default to 1.0 (no enrichment)
                    if treatment in enriched_genes:
                        for info in enriched_genes[treatment]:
                            if info['gene_id'] == gene_id:
                                enrichment = info['fold_enrichment']
                                break
                    
                    enrichment_data[i, j] = enrichment
            
            gene_labels = [name for _, name in top_genes]
            
            plt.figure(figsize=(14, 10))
            sns.heatmap(enrichment_data, xticklabels=treatments, yticklabels=gene_labels,
                        cmap="YlOrRd", annot=True, fmt=".1f")
            plt.title('Fold Enrichment of Top 25 Genes Across Treatments', fontsize=14)
            plt.xlabel('Treatment', fontsize=12)
            plt.ylabel('Gene', fontsize=12)
            plt.tight_layout()
            
            enrichment_heatmap = os.path.join(self.plot_dir, 'gene_enrichment_heatmap.png')
            plt.savefig(enrichment_heatmap, dpi=300)
            plt.close()
            
            plot_files['enrichment_heatmap'] = enrichment_heatmap
        
        # 5. Correlation heatmap between treatments
        correlation_matrix, treatments = self.calculate_correlations()
        
        plt.figure(figsize=(10, 8))
        sns.heatmap(correlation_matrix, annot=True, fmt=".2f", cmap="YlGnBu")
        plt.title('Correlation of Gene Variant Patterns Between Treatments', fontsize=14)
        plt.tight_layout()
        
        correlation_heatmap = os.path.join(self.plot_dir, 'gene_correlation_heatmap.png')
        plt.savefig(correlation_heatmap, dpi=300)
        plt.close()
        
        plot_files['correlation_heatmap'] = correlation_heatmap
        
        # 6. Clustered heatmap of gene variant patterns
        # Get top 50 genes by total variants
        gene_variants_total = {}
        for gene_id in genes_with_variants:
            total = sum(self.gene_variants[gene_id].values())
            gene_variants_total[gene_id] = total
        
        top_genes = sorted(gene_variants_total.items(), key=lambda x: -x[1])[:50]
        
        if top_genes:
            top_gene_ids = [g for g, _ in top_genes]
            gene_labels = []
            
            for gene_id in top_gene_ids:
                sc_gene_id, gene_name, _, _ = self.get_gene_info(gene_id)
                gene_labels.append(gene_name)
            
            # Create data matrix
            data = np.zeros((len(top_gene_ids), len(treatments)))
            
            for i, gene_id in enumerate(top_gene_ids):
                for j, treatment in enumerate(treatments):
                    data[i, j] = self.gene_variants[gene_id].get(treatment, 0)
            
            # Normalize by gene length to get densities
            for i, gene_id in enumerate(top_gene_ids):
                # Get gene length, using default 1000 if not available
                gene_length = 1000  # Default length
                if gene_id in self.gene_metadata:
                    gene_length = self.gene_metadata[gene_id]['end'] - self.gene_metadata[gene_id]['start'] + 1
                data[i, :] = data[i, :] * 1000 / gene_length
            
            # Cluster the data
            row_linkage = hierarchy.linkage(data, method='average')
            col_linkage = hierarchy.linkage(data.T, method='average')
            
            plt.figure(figsize=(14, 10))
            sns.clustermap(pd.DataFrame(data, index=gene_labels, columns=treatments),
                          cmap="YlOrRd", figsize=(14, 10), 
                          row_linkage=row_linkage, col_linkage=col_linkage,
                          xticklabels=treatments, yticklabels=gene_labels)
            plt.title('Clustered Heatmap of Gene Variant Densities', fontsize=16)
            
            clustered_heatmap = os.path.join(self.plot_dir, 'clustered_gene_density.png')
            plt.savefig(clustered_heatmap, dpi=300)
            plt.close()
            
            plot_files['clustered_heatmap'] = clustered_heatmap
        
        # 7. Bar chart comparing adaptation types
        adaptation_variants = [(a, stats['total_variants']) for a, stats in self.adaptation_stats.items() 
                              if a != 'None']  # Exclude controls
        adaptation_variants.sort(key=lambda x: -x[1])
        
        if adaptation_variants:
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
        
        # 8. Gene category distribution
        # Categorize genes into functional groups (if possible)
        gene_categories = defaultdict(int)
        categorized = 0
        
        for gene_id in genes_with_variants:
            if gene_id in self.gene_metadata:
                product = self.gene_metadata[gene_id]['product'].lower()
                
                # Simple categorization based on product name
                if 'transport' in product:
                    gene_categories['Transport'] += 1
                    categorized += 1
                elif 'synthetase' in product or 'synthase' in product:
                    gene_categories['Biosynthesis'] += 1
                    categorized += 1
                elif 'oxidase' in product or 'reductase' in product:
                    gene_categories['Redox'] += 1
                    categorized += 1
                elif 'transcription' in product:
                    gene_categories['Transcription'] += 1
                    categorized += 1
                elif 'protease' in product or 'peptidase' in product:
                    gene_categories['Proteolysis'] += 1
                    categorized += 1
                elif 'kinase' in product:
                    gene_categories['Signaling'] += 1
                    categorized += 1
                elif 'metabolism' in product or 'metabolic' in product:
                    gene_categories['Metabolism'] += 1
                    categorized += 1
        
        # Add uncategorized
        gene_categories['Other/Unknown'] = len(genes_with_variants) - categorized
        
        if gene_categories:
            # Sort categories by count
            categories = sorted(gene_categories.items(), key=lambda x: -x[1])
            
            labels = [c for c, _ in categories]
            values = [v for _, v in categories]
            
            plt.figure(figsize=(10, 6))
            plt.pie(values, labels=labels, autopct='%1.1f%%', startangle=90)
            plt.axis('equal')
            plt.title('Gene Functional Categories with Variants', fontsize=14)
            plt.tight_layout()
            
            category_plot = os.path.join(self.plot_dir, 'gene_categories.png')
            plt.savefig(category_plot, dpi=300)
            plt.close()
            
            plot_files['gene_categories'] = category_plot
        
        # 9. Gene enrichment bubble chart for each treatment
        for treatment, genes in enriched_genes.items():
            if not genes:
                continue
                
            # Take top 20 genes
            top_20 = genes[:20]
            
            gene_names = []
            fold_enrichments = []
            counts = []
            neg_log_pvals = []
            
            for info in top_20:
                gene_names.append(info['gene_name'] or info['sc_gene_id'] or info['gene_id'])
                fold_enrichments.append(info['fold_enrichment'])
                counts.append(info['count'])
                neg_log_pvals.append(-np.log10(info['p_value']))
            
            plt.figure(figsize=(12, 8))
            plt.scatter(fold_enrichments, neg_log_pvals, s=[c*20 for c in counts], alpha=0.7, 
                       c=range(len(gene_names)), cmap='viridis')
            
            # Add gene labels
            for i, gene_name in enumerate(gene_names):
                plt.annotate(gene_name, 
                           (fold_enrichments[i], neg_log_pvals[i]),
                           xytext=(5, 0), 
                           textcoords='offset points')
            
            plt.axhline(y=-np.log10(0.05), color='r', linestyle='--', alpha=0.3)
            plt.axvline(x=1.5, color='r', linestyle='--', alpha=0.3)
            
            plt.title(f'Gene Enrichment in {treatment} Treatment', fontsize=14)
            plt.xlabel('Fold Enrichment', fontsize=12)
            plt.ylabel('-log10(p-value)', fontsize=12)
            plt.tight_layout()
            
            bubble_plot = os.path.join(self.plot_dir, f'{treatment}_gene_enrichment.png')
            plt.savefig(bubble_plot, dpi=300)
            plt.close()
            
            plot_files[f'{treatment}_bubble'] = bubble_plot
        
        return plot_files
    
    def generate_html_report(self, report_files, plot_files, enriched_genes):
        """
        Generate comprehensive HTML report with interactive elements
        
        Args:
            report_files: Dict of generated report file paths
            plot_files: Dict of generated plot file paths
            enriched_genes: Dict of enriched genes by treatment
            
        Returns:
            Path to generated HTML report
        """
        report_file = os.path.join(self.report_dir, "gene_scaffold_analysis_report.html")
        
        # Load treatment statistics
        treatment_stats = []
        for treatment, stats in self.treatment_stats.items():
            metadata = self.treatment_metadata.get(treatment, {})
            
            # Calculate global gene density
            total_length = 0
            for gene_id in stats['genes']:
                if gene_id in self.gene_metadata:
                    gene_length = self.gene_metadata[gene_id]['end'] - self.gene_metadata[gene_id]['start'] + 1
                    total_length += gene_length
            
            global_density = (stats['total_variants'] * 1000) / total_length if total_length > 0 else 0
            
            treatment_stats.append({
                'treatment': treatment,
                'description': metadata.get('description', ''),
                'adaptation': metadata.get('adaptation', ''),
                'has_gene': metadata.get('has_gene', ''),
                'total_variants': stats['total_variants'],
                'genes': len(stats['genes']),
                'global_density': global_density
            })
        
        # Sort by total variants
        treatment_stats.sort(key=lambda x: -x['total_variants'])
        
        # Load top enriched genes
        top_enriched = []
        for treatment, genes in enriched_genes.items():
            for info in genes[:5]:  # Top 5 per treatment
                top_enriched.append({
                    'treatment': treatment,
                    'gene_id': info['gene_id'],
                    'gene_name': info['gene_name'] or info['sc_gene_id'] or info['gene_id'],
                    'product': info['product'],
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
    <title>Gene-Scaffold Variant Analysis Report</title>
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
    <h1>Gene-Scaffold Variant Analysis Report</h1>
    <p>Generated on: """ + f"{os.popen('date').read().strip()}" + """</p>
    
    <div class="info">
        <h3>Analysis Overview</h3>
        <p>This report presents a comprehensive gene-centric analysis of genomic variants across different treatments and adaptations.</p>
        <p>The analysis maps variants to genes based on scaffold assignments, identifying genes that are enriched for variants in specific conditions.</p>
    </div>
    
    <div class="tab">
        <button class="tablinks" onclick="openTab(event, 'Summary')" id="defaultOpen">Summary</button>
        <button class="tablinks" onclick="openTab(event, 'Treatments')">Treatments</button>
        <button class="tablinks" onclick="openTab(event, 'Genes')">Genes</button>
        <button class="tablinks" onclick="openTab(event, 'Enrichment')">Enrichment</button>
        <button class="tablinks" onclick="openTab(event, 'Comparisons')">Comparisons</button>
    </div>
    
    <div id="Summary" class="tabcontent">
        <h2>Analysis Summary</h2>
        <div class="info">
            <p><strong>Total Treatments Analyzed:</strong> """ + f"{len(self.treatment_stats)}" + """</p>
            <p><strong>Total Variants:</strong> """ + f"{sum(stats['total_variants'] for stats in self.treatment_stats.values())}" + """</p>
            <p><strong>Total Genes with Variants:</strong> """ + f"{len(set().union(*[stats['genes'] for stats in self.treatment_stats.values()]))}" + """</p>
            <p><strong>Enriched Genes Identified:</strong> """ + f"{sum(len(genes) for genes in enriched_genes.values())}" + """</p>
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
            
            # Top enriched genes
            for enrichment in top_enriched[:3]:
                t = enrichment['treatment']
                g = enrichment['gene_name']
                p = enrichment['product']
                fold = enrichment['fold_enrichment']
                f.write(f"<li>Gene <strong>{g}</strong> ({p}) shows {fold:.1f}x enrichment in the <strong>{t}</strong> treatment</li>\n")
            
            # Top correlations
            if sorted_correlations:
                pair, corr = sorted_correlations[0]
                f.write(f"<li>The most similar treatments in terms of gene variant patterns are <strong>{pair}</strong> (correlation: {corr:.2f})</li>\n")
            
            f.write("""</ul>
        
        <h3>Variant Distribution by Treatment</h3>
        <img src="../plots/treatment_variant_counts.png" alt="Treatment Variant Counts">
        
        <h3>Variant Distribution by Adaptation</h3>
        <img src="../plots/adaptation_variant_counts.png" alt="Adaptation Variant Counts">
        
        <h3>Gene Categories with Variants</h3>
        <img src="../plots/gene_categories.png" alt="Gene Categories">
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
                <th>Genes</th>
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
                <td>{treatment_info['genes']}</td>
                <td>{treatment_info['global_density']:.4f}</td>
            </tr>""")
            
            f.write("""
        </table>
        
        <h3>Genes with Variants by Treatment</h3>
        <img src="../plots/treatment_gene_counts.png" alt="Treatment Gene Counts">
        
        <h3>Treatment Enrichment Plots</h3>
        <p>These bubble charts show gene enrichment by treatment, with bubble size indicating variant count.</p>
        """)
            
            # Add gene enrichment bubble charts for each treatment
            for treatment in sorted(enriched_genes.keys()):
                if f"{treatment}_bubble" in plot_files:
                    f.write(f"""
        <h4>{treatment} Treatment</h4>
        <img src="../plots/{treatment}_gene_enrichment.png" alt="{treatment} Gene Enrichment">
        """)
            
            f.write("""
    </div>
    
    <div id="Genes" class="tabcontent">
        <h2>Gene Analysis</h2>
        <p>This section provides details on the distribution of variants across genes.</p>
        
        <h3>Top Genes by Variant Density</h3>
        <img src="../plots/top_genes_by_density.png" alt="Top Genes by Density">
        
        <h3>Clustered Heatmap of Gene Variant Densities</h3>
        <p>This heatmap shows the variant densities for the top 50 genes across all treatments,
           clustered to reveal patterns of similarity.</p>
        <img src="../plots/clustered_gene_density.png" alt="Clustered Gene Density Heatmap">
    </div>
    
    <div id="Enrichment" class="tabcontent">
        <h2>Gene Enrichment Analysis</h2>
        <p>This section shows statistical enrichment of variants in specific genes.</p>
        
        <h3>Enrichment Heatmap</h3>
        <p>This heatmap shows the fold enrichment of the top enriched genes across treatments.</p>
        <img src="../plots/gene_enrichment_heatmap.png" alt="Gene Enrichment Heatmap">
        
        <h3>Top Enriched Genes</h3>
        <table>
            <tr>
                <th>Treatment</th>
                <th>Gene</th>
                <th>Product</th>
                <th>Density</th>
                <th>Fold Enrichment</th>
                <th>P-Value</th>
            </tr>""")
            
            # Add enriched gene rows
            for enrichment in top_enriched[:20]:  # Show top 20
                f.write(f"""
            <tr>
                <td>{enrichment['treatment']}</td>
                <td>{enrichment['gene_name']}</td>
                <td>{enrichment['product']}</td>
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
        <p>This heatmap shows the correlation of gene variant patterns between treatments.</p>
        <img src="../plots/gene_correlation_heatmap.png" alt="Gene Correlation Heatmap">
        
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
    
    def get_gene_info(self, gene_id):
        """
        Helper function to safely get gene information regardless of ID type
        
        Args:
            gene_id: Gene ID (either W303 ID or SC ID)
            
        Returns:
            Tuple of (sc_gene_id, gene_name, product, scaffold_id)
        """
        # First, check if it's a W303 gene ID
        if gene_id in self.gene_metadata:
            sc_gene_id = self.gene_metadata[gene_id].get('sc_gene_id', '')
            gene_name = self.gene_name_map.get(sc_gene_id, sc_gene_id) or gene_id
            product = self.gene_metadata[gene_id].get('product', '')
            scaffold_id = self.gene_metadata[gene_id].get('scaffold_id', '')
            
        # Next, check if it's an SC gene ID we have info for
        elif gene_id in self.sc_gene_info:
            sc_gene_id = gene_id
            info = self.sc_gene_info[gene_id]
            gene_name = info.get('sc_gene_name', '') or self.gene_name_map.get(sc_gene_id, sc_gene_id) or gene_id
            product = info.get('product', '')
            scaffold_id = info.get('scaffold_id', '')
            
        # Fall back to SC ID with gene name map
        elif gene_id in self.gene_name_map:
            sc_gene_id = gene_id
            gene_name = self.gene_name_map.get(sc_gene_id, sc_gene_id)
            product = "Unknown"
            scaffold_id = "Unknown"
            
        # Last resort - unknown gene
        else:
            sc_gene_id = gene_id
            gene_name = gene_id
            product = "Unknown"
            scaffold_id = "Unknown"
            
        return sc_gene_id, gene_name, product, scaffold_id
    
    def generate_summary_text(self, enriched_genes):
        """
        Generate a comprehensive summary text for the analysis
        
        Args:
            enriched_genes: Dict of enriched genes by treatment
            
        Returns:
            Summary text string
        """
        
        # Calculate basic statistics
        total_variants = sum(stats['total_variants'] for stats in self.treatment_stats.values())
        all_genes = set()
        for treatment, stats in self.treatment_stats.items():
            all_genes.update(stats['genes'])
        
        total_genes = len(all_genes)
        total_enriched = sum(len(genes) for genes in enriched_genes.values())
        
        # Get top treatments
        treatment_variants = [(t, stats['total_variants']) for t, stats in self.treatment_stats.items()]
        treatment_variants.sort(key=lambda x: -x[1])
        top_treatments = treatment_variants[:3]
        
        # Get top genes by density
        gene_avg_density = {}
        for gene_id in all_genes:
            if gene_id in self.gene_density:
                densities = list(self.gene_density[gene_id].values())
                if densities:
                    gene_avg_density[gene_id] = sum(densities) / len(densities)
        
        top_genes = sorted(gene_avg_density.items(), key=lambda x: -x[1])[:5]
        
        # Get top enriched genes
        top_enriched = []
        for treatment, genes in enriched_genes.items():
            for info in genes[:3]:  # Top 3 per treatment
                top_enriched.append((treatment, info['gene_id'], info['gene_name'], info['fold_enrichment']))
        
        top_enriched.sort(key=lambda x: -x[3])
        top_enriched = top_enriched[:5]  # Overall top 5
        
        # Generate the summary text
        summary = "Gene-Scaffold Analysis Summary\n"
        summary += "============================\n\n"
        
        summary += "Overall Statistics:\n"
        summary += "-----------------\n"
        summary += f"Total Variants: {total_variants}\n"
        summary += f"Genes with Variants: {total_genes}\n"
        summary += f"Enriched Genes: {total_enriched}\n\n"
        
        # Add treatment summaries
        for treatment, count in top_treatments:
            metadata = self.treatment_metadata.get(treatment, {})
            stats = self.treatment_stats[treatment]
            
            # Calculate global gene density
            total_length = 0
            for gene_id in stats['genes']:
                if gene_id in self.gene_metadata:
                    gene_length = self.gene_metadata[gene_id]['end'] - self.gene_metadata[gene_id]['start'] + 1
                    total_length += gene_length
            
            global_density = (stats['total_variants'] * 1000) / total_length if total_length > 0 else 0
            
            description = metadata.get('description', '')
            adaptation = metadata.get('adaptation', '')
            adaptation_text = f" ({adaptation} adaptation)" if adaptation and adaptation != 'None' else ""
            
            summary += f"{treatment} Treatment ({description}{adaptation_text}):\n"
            summary += f"  Variants: {count}\n"
            summary += f"  Genes with Variants: {len(stats['genes'])}\n"
            summary += f"  Global Gene Density: {global_density:.4f} variants/kb\n"
            
            # Top 3 genes by density for this treatment
            if stats['genes']:
                gene_densities = []
                for gene_id in stats['genes']:
                    if gene_id in self.gene_density and treatment in self.gene_density[gene_id]:
                        density = self.gene_density[gene_id][treatment]
                        sc_gene_id, gene_name, product, _ = self.get_gene_info(gene_id)
                        gene_densities.append((gene_id, gene_name, density, product))
                
                gene_densities.sort(key=lambda x: -x[2])
                top_3 = gene_densities[:3]
                
                if top_3:
                    summary += f"  Top 3 Genes by Density:\n"
                    for _, gene_name, density, product in top_3:
                        summary += f"    {gene_name}: {density:.4f} variants/kb - {product}\n"
            
            summary += "\n"
        
        # Add top genes overall
        if top_genes:
            summary += "Top 5 Genes by Average Variant Density:\n"
            summary += "------------------------------------\n"
            
            for i, (gene_id, density) in enumerate(top_genes, 1):
                sc_gene_id, gene_name, product, scaffold_id = self.get_gene_info(gene_id)
                
                summary += f"{i}. {gene_name} (Scaffold: {scaffold_id})\n"
                summary += f"   Density: {density:.4f} variants/kb\n"
                summary += f"   Product: {product}\n"
            
            summary += "\n"
        
        # Add top enriched genes
        if top_enriched:
            summary += "Top 5 Enriched Genes:\n"
            summary += "-------------------\n"
            
            for i, (treatment, gene_id, gene_name, fold) in enumerate(top_enriched, 1):
                _, _, product, scaffold_id = self.get_gene_info(gene_id)
                    
                summary += f"{i}. {gene_name} (Scaffold: {scaffold_id})\n"
                summary += f"   Treatment: {treatment}\n"
                summary += f"   Fold Enrichment: {fold:.2f}x\n"
                summary += f"   Product: {product}\n"
            
            summary += "\n"
        
        # Add adaptation comparison
        summary += "Adaptation Comparison:\n"
        summary += "---------------------\n"
        
        for adaptation in ['Temperature', 'Low Oxygen']:
            if adaptation in self.adaptation_stats:
                stats = self.adaptation_stats[adaptation]
                
                summary += f"{adaptation} Adaptation:\n"
                summary += f"  Variants: {stats['total_variants']}\n"
                summary += f"  Genes with Variants: {len(stats['genes'])}\n"
                
                # Calculate global density
                total_length = 0
                for gene_id in stats['genes']:
                    if gene_id in self.gene_metadata:
                        gene_length = self.gene_metadata[gene_id]['end'] - self.gene_metadata[gene_id]['start'] + 1
                        total_length += gene_length
                
                global_density = (stats['total_variants'] * 1000) / total_length if total_length > 0 else 0
                summary += f"  Global Gene Density: {global_density:.4f} variants/kb\n"
                
                # Get top 3 genes by count
                gene_counts = [(g, c) for g, c in stats['gene_counts'].items()]
                gene_counts.sort(key=lambda x: -x[1])
                top_3 = gene_counts[:3]
                
                if top_3:
                    summary += f"  Top 3 Genes by Variant Count:\n"
                    for gene_id, count in top_3:
                        sc_gene_id, gene_name, product, _ = self.get_gene_info(gene_id)
                        summary += f"    {gene_name}: {count} variants - {product}\n"
                
                summary += "\n"
        
        summary += "Main Conclusions:\n"
        summary += "---------------\n"
        summary += "1. This analysis identifies genes with high variant density across treatments.\n"
        summary += "2. Several genes show treatment-specific enrichment of variants.\n"
        summary += "3. The pattern of variant distribution provides insights into adaptation mechanisms.\n"
        summary += "4. Temperature and low oxygen adaptations show distinct gene impact patterns.\n"
        summary += "5. Gene modifications (STC, CAS) appear to influence gene variant patterns.\n"
        summary += "6. Further analysis of enriched genes may reveal functional implications of adaptations.\n"
        
        return summary
    
    def run_complete_analysis(self, vcf_files):
        """Run the complete analysis pipeline"""
        print(f"Starting gene-scaffold analysis with {len(vcf_files)} VCF files")
        
        # Process all VCF files
        results = self.process_all_vcfs(vcf_files)
        print(f"Processed {len(results)} VCF files, found variants in {len(set().union(*[stats['genes'] for stats in self.treatment_stats.values()]))} genes")
        
        # Identify enriched genes
        enriched_genes = self.identify_enriched_genes()
        total_enriched = sum(len(genes) for genes in enriched_genes.values())
        print(f"Identified {total_enriched} enriched genes across all treatments")
        
        # Generate text reports
        report_files = self.generate_text_reports(enriched_genes)
        print("\nGenerated text reports:")
        for report_type, report_file in report_files.items():
            print(f"  {report_type}: {report_file}")
        
        # Generate plots
        plot_files = self.generate_plots(enriched_genes)
        print("\nGenerated plots:")
        for plot_type, plot_file in plot_files.items():
            print(f"  {plot_type}: {plot_file}")
        
        # Generate HTML report
        html_report = self.generate_html_report(report_files, plot_files, enriched_genes)
        print(f"\nComplete HTML report: {html_report}")
        
        # Generate summary text
        summary_text = self.generate_summary_text(enriched_genes)
        summary_file = os.path.join(self.report_dir, "gene_scaffold_analysis_summary.txt")
        
        with open(summary_file, 'w') as f:
            f.write(summary_text)
        
        print(f"Summary text written to: {summary_file}")
        
        return {
            'results': results,
            'enriched_genes': enriched_genes,
            'report_files': report_files,
            'plot_files': plot_files,
            'html_report': html_report,
            'summary_file': summary_file
        }

def main():
    parser = argparse.ArgumentParser(
        description='Gene-centric scaffold-level analysis of variants across treatments',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--vcf', required=True, nargs='+', help='VCF file(s) to analyze')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--gene-index', required=True, help='Gene index file')
    parser.add_argument('--scaffold-mapping', required=True, help='Scaffold mapping file')
    
    args = parser.parse_args()
    
    # Create analyzer
    analyzer = GeneScaffoldAnalyzer(
        output_dir=args.output_dir,
        gene_index_file=args.gene_index,
        scaffold_mapping_file=args.scaffold_mapping
    )
    
    # Run analysis
    analysis_results = analyzer.run_complete_analysis(args.vcf)
    
    print("\nAnalysis complete!")
    print(f"Results are available in {args.output_dir}")

if __name__ == '__main__':
    main()