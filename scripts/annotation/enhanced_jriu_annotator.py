#!/usr/bin/env python3
"""
enhanced_jriu_annotator.py - Comprehensive JRIU-based annotation system for yeast variants

This script implements a JRIU-based approach to variant annotation, addressing
the coordinate system mismatch between VCF files and gene annotations.
It provides comprehensive analysis and visualization of variant patterns
in relation to genes of interest.
"""

import os
import sys
import csv
import argparse
import gzip
import json
import math
from collections import defaultdict, Counter
from scipy import stats as scipy_stats  # Import with a different name
from Bio import SeqIO
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

class EnhancedJRIUAnnotator:
    def __init__(self, genes_of_interest_file, output_dir):
        """
        Initialize the enhanced JRIU-based variant annotator.
        
        Args:
            genes_of_interest_file: Path to TSV file with genes of interest
            output_dir: Directory for output files and reports
        """
        self.genes_of_interest_file = genes_of_interest_file
        self.output_dir = output_dir
        
        # Create output directories
        self.vcf_dir = os.path.join(output_dir, "annotated_vcfs")
        self.report_dir = os.path.join(output_dir, "reports")
        self.plot_dir = os.path.join(output_dir, "plots")
        
        os.makedirs(self.vcf_dir, exist_ok=True)
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
        
        # Reverse lookup: sample to treatment
        self.sample_to_treatment = {}
        for treatment, samples in self.treatment_groups.items():
            for sample in samples:
                self.sample_to_treatment[sample] = treatment
        
        # Data containers
        self.jriu_to_genes = defaultdict(list)   # Maps JRIU to gene info
        self.gene_to_jrius = defaultdict(list)   # Maps gene ID to JRIU IDs
        self.gene_info = {}                      # Gene metadata
        self.sc_gene_map = {}                    # SC gene ID to W303 gene ID mapping
        
        # Statistics containers
        self.sample_stats = {}                   # Statistics by sample
        self.gene_stats = defaultdict(dict)      # Statistics by gene
        self.treatment_stats = defaultdict(dict) # Statistics by treatment
        
        # Load gene information
        self.load_genes_of_interest()
    
    def load_genes_of_interest(self):
        """Load and organize information about genes of interest"""
        print(f"Loading genes of interest from {self.genes_of_interest_file}")
        
        with open(self.genes_of_interest_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gene_id = row.get('w303_gene_id')
                sc_gene_id = row.get('sc_gene_id')
                jriu_id = row.get('jriu_id')
                strand = row.get('strand')
                
                if gene_id and jriu_id:
                    # Store gene to JRIU mapping
                    self.gene_to_jrius[gene_id].append(jriu_id)
                    
                    # Store JRIU to gene mapping
                    self.jriu_to_genes[jriu_id].append({
                        'gene_id': gene_id,
                        'sc_gene_id': sc_gene_id,
                        'strand': strand
                    })
                    
                    # Store SC gene ID mapping
                    if sc_gene_id:
                        self.sc_gene_map[sc_gene_id] = gene_id
                    
                    # Store gene information
                    self.gene_info[gene_id] = {
                        'sc_gene_id': sc_gene_id,
                        'jriu_ids': self.gene_to_jrius[gene_id],
                        'strand': strand
                    }
        
        print(f"Loaded {len(self.gene_info)} genes of interest across {len(self.jriu_to_genes)} JRIUs")
    
    def process_vcf_file(self, vcf_file):
        """
        Process a single VCF file to identify variants on JRIUs with genes of interest
        
        Args:
            vcf_file: Path to VCF file
            
        Returns:
            Dict containing statistics and annotated variants
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
        jriu_variants = 0
        gene_jriu_counts = defaultdict(int)
        
        # Store annotated variants
        annotated_variants = []
        
        # Store header lines
        header_lines = []
        
        # Process VCF file
        with opener(vcf_file, mode) as f:
            # Process header
            for line in f:
                if line.startswith('##'):
                    header_lines.append(line.strip())
                    continue
                elif line.startswith('#CHROM'):
                    header_lines.append(line.strip())
                    break
            
            # Process variants
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                fields = line.split('\t')
                total_variants += 1
                
                # Extract variant info
                chrom = fields[0]  # JRIU ID
                pos = int(fields[1])
                variant_id = fields[2]
                ref = fields[3]
                alt = fields[4]
                qual = fields[5]
                filter_val = fields[6]
                info = fields[7]
                
                # Check if this JRIU has genes of interest
                if chrom in self.jriu_to_genes:
                    jriu_variants += 1
                    genes = self.jriu_to_genes[chrom]
                    
                    # Update counts for each gene
                    for gene in genes:
                        gene_id = gene['gene_id']
                        gene_jriu_counts[gene_id] += 1
                    
                    # Add annotation to variant
                    gene_info = ";".join([f"{g['gene_id']}({g['sc_gene_id']})" for g in genes])
                    
                    # Create annotated variant
                    annotated_variants.append({
                        'chrom': chrom,
                        'pos': pos,
                        'id': variant_id,
                        'ref': ref,
                        'alt': alt,
                        'qual': qual,
                        'filter': filter_val,
                        'info': f"{info};GOI=true;GenesOnJRIU={gene_info}",
                        'format': "GT",
                        'sample': "1",
                        'genes': [g['gene_id'] for g in genes]
                    })
        
        # Write annotated VCF
        output_vcf = os.path.join(self.vcf_dir, f"{sample_name}.annotated.vcf")
        with open(output_vcf, 'w') as f:
            # Write header
            for line in header_lines[:-1]:
                f.write(line + '\n')
            
            # Add custom INFO headers
            f.write('##INFO=<ID=GOI,Number=0,Type=Flag,Description="Variant is on a JRIU containing genes of interest">\n')
            f.write('##INFO=<ID=GenesOnJRIU,Number=.,Type=String,Description="Genes of interest on this JRIU">\n')
            
            # Write column names
            f.write(header_lines[-1] + '\n')
            
            # Write variants
            for variant in annotated_variants:
                f.write(f"{variant['chrom']}\t{variant['pos']}\t{variant['id']}\t"
                       f"{variant['ref']}\t{variant['alt']}\t{variant['qual']}\t"
                       f"{variant['filter']}\t{variant['info']}\t{variant['format']}\t{variant['sample']}\n")
        
        # Store results
        result = {
            'sample': sample_name,
            'treatment': treatment,
            'total_variants': total_variants,
            'jriu_variants': jriu_variants,
            'gene_jriu_counts': dict(gene_jriu_counts),
            'output_vcf': output_vcf
        }
        
        print(f"  Found {jriu_variants}/{total_variants} variants on JRIUs with genes of interest")
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
            
            # Store sample statistics
            self.sample_stats[result['sample']] = {
                'treatment': result['treatment'],
                'total_variants': result['total_variants'],
                'jriu_variants': result['jriu_variants'],
                'gene_counts': result['gene_jriu_counts']
            }
            
            # Update gene statistics
            for gene_id, count in result['gene_jriu_counts'].items():
                if gene_id not in self.gene_stats:
                    self.gene_stats[gene_id] = {
                        'total_variants': 0,
                        'sample_counts': defaultdict(int),
                        'treatment_counts': defaultdict(int)
                    }
                
                self.gene_stats[gene_id]['total_variants'] += count
                self.gene_stats[gene_id]['sample_counts'][result['sample']] += count
                self.gene_stats[gene_id]['treatment_counts'][result['treatment']] += count
            
            # Update treatment statistics
            treatment = result['treatment']
            if 'total_variants' not in self.treatment_stats[treatment]:
                self.treatment_stats[treatment] = {
                    'total_variants': 0,
                    'jriu_variants': 0,
                    'sample_count': 0,
                    'gene_counts': defaultdict(int)
                }
            
            self.treatment_stats[treatment]['total_variants'] += result['total_variants']
            self.treatment_stats[treatment]['jriu_variants'] += result['jriu_variants']
            self.treatment_stats[treatment]['sample_count'] += 1
            
            for gene_id, count in result['gene_jriu_counts'].items():
                self.treatment_stats[treatment]['gene_counts'][gene_id] += count
        
        return all_results
    
    def calculate_enrichment_statistics(self):
        """
        Calculate statistical enrichment of variants on gene JRIUs
        
        Returns:
            Dict with enrichment statistics
        """
        from scipy import stats as scipy_stats  # Import with a different name
        
        enrichment_stats = {}
        
        # Calculate baseline rates across all samples
        total_variants = sum(stats['total_variants'] for stats in self.sample_stats.values())
        total_jriu_variants = sum(stats['jriu_variants'] for stats in self.sample_stats.values())
        baseline_rate = total_jriu_variants / total_variants if total_variants > 0 else 0
        
        print(f"\nCalculating enrichment statistics (baseline rate: {baseline_rate:.4f})")
        
        # Calculate enrichment for each treatment
        for treatment, treatment_stats in self.treatment_stats.items():
            treatment_variant_rate = treatment_stats['jriu_variants'] / treatment_stats['total_variants'] if treatment_stats['total_variants'] > 0 else 0
            enrichment_factor = treatment_variant_rate / baseline_rate if baseline_rate > 0 else 0
            
            # Calculate statistical significance
            if baseline_rate > 0 and baseline_rate < 1:
                # Use binomial test - note we use scipy_stats (not stats)
                try:
                    # Try newer SciPy version function
                    p_value = scipy_stats.binomtest(treatment_stats['jriu_variants'], 
                                                treatment_stats['total_variants'], 
                                                baseline_rate).pvalue
                except AttributeError:
                    # Fall back to older version
                    p_value = scipy_stats.binom_test(treatment_stats['jriu_variants'], 
                                                    treatment_stats['total_variants'], 
                                                    baseline_rate)
            else:
                p_value = 1.0
            
            enrichment_stats[treatment] = {
                'baseline_rate': baseline_rate,
                'treatment_rate': treatment_variant_rate,
                'enrichment_factor': enrichment_factor,
                'p_value': p_value,
                'significant': p_value < 0.05
            }
            
            print(f"  {treatment}: Enrichment = {enrichment_factor:.2f}x (p = {p_value:.4f})")
        
        # Calculate gene-specific enrichment
        gene_enrichment = {}
        for gene_id, gene_stats in self.gene_stats.items():
            gene_enrichment[gene_id] = {}
            
            for treatment, treatment_count in gene_stats['treatment_counts'].items():
                # Calculate expected count based on total variants in treatment
                treatment_total = self.treatment_stats[treatment]['total_variants']
                expected_rate = gene_stats['total_variants'] / total_variants if total_variants > 0 else 0
                expected_count = expected_rate * treatment_total
                
                # Calculate enrichment
                enrichment_factor = treatment_count / expected_count if expected_count > 0 else 0
                
                # Calculate statistical significance with Poisson test
                if expected_count > 0:
                    p_value = scipy_stats.poisson.sf(treatment_count - 1, expected_count)
                else:
                    p_value = 1.0
                
                gene_enrichment[gene_id][treatment] = {
                    'actual_count': treatment_count,
                    'expected_count': expected_count,
                    'enrichment_factor': enrichment_factor,
                    'p_value': p_value,
                    'significant': p_value < 0.05
                }
        
        return {
            'treatment_enrichment': enrichment_stats,
            'gene_enrichment': gene_enrichment
        }
    
    def generate_text_reports(self):
        """Generate comprehensive text-based reports"""
        # Sample summary report
        sample_report_file = os.path.join(self.report_dir, "sample_variant_summary.tsv")
        with open(sample_report_file, 'w') as f:
            f.write("Sample\tTreatment\tTotalVariants\tGeneJRIUVariants\tPercentage\n")
            for sample, stats in sorted(self.sample_stats.items()):
                pct = 100 * stats['jriu_variants'] / stats['total_variants'] if stats['total_variants'] > 0 else 0
                f.write(f"{sample}\t{stats['treatment']}\t{stats['total_variants']}\t{stats['jriu_variants']}\t{pct:.2f}\n")
        
        # Treatment summary report
        treatment_report_file = os.path.join(self.report_dir, "treatment_variant_summary.tsv")
        with open(treatment_report_file, 'w') as f:
            f.write("Treatment\tTotalVariants\tGeneJRIUVariants\tPercentage\tSamples\tAvgVariantsPerSample\n")
            for treatment, stats in sorted(self.treatment_stats.items()):
                pct = 100 * stats['jriu_variants'] / stats['total_variants'] if stats['total_variants'] > 0 else 0
                avg = stats['total_variants'] / stats['sample_count'] if stats['sample_count'] > 0 else 0
                f.write(f"{treatment}\t{stats['total_variants']}\t{stats['jriu_variants']}\t{pct:.2f}\t{stats['sample_count']}\t{avg:.2f}\n")
        
        # Gene summary report
        gene_report_file = os.path.join(self.report_dir, "gene_variant_summary.tsv")
        with open(gene_report_file, 'w') as f:
            f.write("GeneID\tSCGeneID\tGeneName\tTotalVariants\n")
            for gene_id, stats in sorted(self.gene_stats.items(), key=lambda x: -x[1]['total_variants']):
                sc_gene_id = self.gene_info.get(gene_id, {}).get('sc_gene_id', "unknown")
                gene_name = sc_gene_id.split('Y')[0] if sc_gene_id.startswith('Y') else sc_gene_id
                f.write(f"{gene_id}\t{sc_gene_id}\t{gene_name}\t{stats['total_variants']}\n")
        
        # Gene by treatment report
        gene_treatment_file = os.path.join(self.report_dir, "gene_treatment_summary.tsv")
        with open(gene_treatment_file, 'w') as f:
            f.write("GeneID\tSCGeneID\tGeneName\tTreatment\tVariants\n")
            for gene_id, stats in sorted(self.gene_stats.items()):
                sc_gene_id = self.gene_info.get(gene_id, {}).get('sc_gene_id', "unknown")
                gene_name = sc_gene_id.split('Y')[0] if sc_gene_id.startswith('Y') else sc_gene_id
                
                for treatment, count in sorted(stats['treatment_counts'].items()):
                    f.write(f"{gene_id}\t{sc_gene_id}\t{gene_name}\t{treatment}\t{count}\n")
        
        # Enrichment statistics report
        enrichment_stats = self.calculate_enrichment_statistics()
        
        enrichment_file = os.path.join(self.report_dir, "enrichment_statistics.tsv")
        with open(enrichment_file, 'w') as f:
            f.write("Type\tID\tComparison\tActualCount\tExpectedCount\tEnrichmentFactor\tPValue\tSignificant\n")
            
            # Treatment enrichment
            for treatment, stats in enrichment_stats['treatment_enrichment'].items():
                f.write(f"Treatment\t{treatment}\tAll\t")
                f.write(f"{self.treatment_stats[treatment]['jriu_variants']}\t")
                expected = self.treatment_stats[treatment]['total_variants'] * stats['baseline_rate']
                f.write(f"{expected:.2f}\t{stats['enrichment_factor']:.2f}\t{stats['p_value']:.4e}\t{stats['significant']}\n")
            
            # Gene enrichment by treatment
            for gene_id, treatments in enrichment_stats['gene_enrichment'].items():
                sc_gene_id = self.gene_info.get(gene_id, {}).get('sc_gene_id', "unknown")
                gene_label = f"{gene_id}({sc_gene_id})"
                
                for treatment, stats in treatments.items():
                    f.write(f"Gene\t{gene_label}\t{treatment}\t")
                    f.write(f"{stats['actual_count']}\t{stats['expected_count']:.2f}\t")
                    f.write(f"{stats['enrichment_factor']:.2f}\t{stats['p_value']:.4e}\t{stats['significant']}\n")
        
        return {
            'sample_report': sample_report_file,
            'treatment_report': treatment_report_file,
            'gene_report': gene_report_file,
            'gene_treatment_report': gene_treatment_file,
            'enrichment_report': enrichment_file
        }
    
    def generate_plots(self):
        """Generate visualization plots for the analysis"""
        plot_files = {}
        
        # Set visualization style
        plt.style.use('seaborn-v0_8-whitegrid')
        
        # 1. Treatment variant counts
        treatment_counts = [(t, stats['jriu_variants']) for t, stats in self.treatment_stats.items()]
        treatment_counts.sort(key=lambda x: x[1], reverse=True)
        
        treatments = [t for t, _ in treatment_counts]
        counts = [c for _, c in treatment_counts]
        
        plt.figure(figsize=(10, 6))
        ax = plt.bar(treatments, counts, color='skyblue')
        plt.title('Variants on Gene JRIUs by Treatment', fontsize=14)
        plt.xlabel('Treatment', fontsize=12)
        plt.ylabel('Number of Variants', fontsize=12)
        plt.xticks(rotation=45)
        plt.tight_layout()
        
        treatment_plot = os.path.join(self.plot_dir, 'treatment_variant_counts.png')
        plt.savefig(treatment_plot, dpi=300)
        plt.close()
        
        plot_files['treatment_counts'] = treatment_plot
        
        # 2. Gene variant counts
        gene_counts = [(gene_id, stats['total_variants']) for gene_id, stats in self.gene_stats.items()]
        gene_counts.sort(key=lambda x: x[1], reverse=True)

        # Map of SC gene IDs to common names
        GENE_NAME_MAP = {
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
        
        gene_ids = []
        for gene_id, _ in gene_counts:
            sc_gene_id = self.gene_info.get(gene_id, {}).get('sc_gene_id', "unknown")
            # Instead of splitting at 'Y', use the whole SC gene ID as the label
            # This will show YHR190W instead of an empty string
            #gene_name = sc_gene_id if sc_gene_id != "unknown" else gene_id
            # Then use this in your gene label creation
            gene_name = GENE_NAME_MAP.get(sc_gene_id, sc_gene_id) if sc_gene_id != "unknown" else gene_id
            gene_ids.append(gene_name)

        counts = [c for _, c in gene_counts]

        plt.figure(figsize=(12, 6))
        ax = plt.bar(gene_ids, counts, color='lightgreen')
        plt.title('Variants Associated with Genes of Interest', fontsize=14)
        plt.xlabel('Gene', fontsize=12)
        plt.ylabel('Number of Variants', fontsize=12)
        plt.xticks(rotation=45)
        plt.tight_layout()

        gene_plot = os.path.join(self.plot_dir, 'gene_variant_counts.png')
        plt.savefig(gene_plot, dpi=300)
        plt.close()

        plot_files['gene_counts'] = gene_plot
        
        # 3. Gene by treatment heatmap
        # Prepare data for heatmap
        gene_names = []
        for gene_id, _ in gene_counts:
            sc_gene_id = self.gene_info.get(gene_id, {}).get('sc_gene_id', "unknown")
            gene_name = sc_gene_id.split('Y')[0] if sc_gene_id.startswith('Y') else sc_gene_id
            gene_names.append(gene_name)
        
        treatment_names = sorted(self.treatment_stats.keys())
        
        # Create data matrix
        heatmap_data = np.zeros((len(gene_names), len(treatment_names)))
        
        for i, (gene_id, _) in enumerate(gene_counts):
            for j, treatment in enumerate(treatment_names):
                count = self.gene_stats[gene_id]['treatment_counts'].get(treatment, 0)
                heatmap_data[i, j] = count
        
        plt.figure(figsize=(12, 8))
        sns.heatmap(heatmap_data, annot=True, fmt="d", cmap="YlGnBu",
                    xticklabels=treatment_names, yticklabels=gene_names)
        plt.title('Gene Variants by Treatment', fontsize=14)
        plt.xlabel('Treatment', fontsize=12)
        plt.ylabel('Gene', fontsize=12)
        plt.tight_layout()
        
        heatmap_plot = os.path.join(self.plot_dir, 'gene_treatment_heatmap.png')
        plt.savefig(heatmap_plot, dpi=300)
        plt.close()
        
        plot_files['heatmap'] = heatmap_plot
        
        # 4. Enrichment plot
        enrichment_stats = self.calculate_enrichment_statistics()
        
        # Prepare data
        gene_labels = []
        treatment_data = defaultdict(list)
        
        for gene_id, treatments in enrichment_stats['gene_enrichment'].items():
            sc_gene_id = self.gene_info.get(gene_id, {}).get('sc_gene_id', "unknown")
            gene_name = sc_gene_id.split('Y')[0] if sc_gene_id.startswith('Y') else sc_gene_id
            gene_labels.append(gene_name)
            
            for treatment, stats in treatments.items():
                treatment_data[treatment].append(stats['enrichment_factor'])
        
        plt.figure(figsize=(12, 6))
        
        # Set width of bars
        bar_width = 0.15
        index = np.arange(len(gene_labels))
        
        # Plot bars for each treatment
        for i, (treatment, values) in enumerate(sorted(treatment_data.items())):
            # Make sure values match the length of gene_labels
            values_array = np.array(values)
            if len(values_array) != len(gene_labels):
                # Pad or truncate values to match gene_labels length
                if len(values_array) < len(gene_labels):
                    values_array = np.pad(values_array, (0, len(gene_labels) - len(values_array)), 'constant')
                else:
                    values_array = values_array[:len(gene_labels)]
            plt.bar(index + i*bar_width, values_array, bar_width, label=treatment, alpha=0.7)
        
        plt.xlabel('Gene', fontsize=12)
        plt.ylabel('Enrichment Factor', fontsize=12)
        plt.title('Gene Variant Enrichment by Treatment', fontsize=14)
        plt.xticks(index + bar_width * (len(treatment_data) - 1) / 2, gene_labels, rotation=45)
        plt.legend()
        plt.tight_layout()
        
        enrichment_plot = os.path.join(self.plot_dir, 'gene_enrichment.png')
        plt.savefig(enrichment_plot, dpi=300)
        plt.close()
        
        plot_files['enrichment'] = enrichment_plot
        
        # 5. Treatment percentage plot (stacked bar)
        # Calculate percentages
        treatment_percentages = {}
        for treatment, stats in self.treatment_stats.items():
            if stats['total_variants'] > 0:
                treatment_percentages[treatment] = 100 * stats['jriu_variants'] / stats['total_variants']
            else:
                treatment_percentages[treatment] = 0
        
        treatments = sorted(treatment_percentages.keys())
        percentages = [treatment_percentages[t] for t in treatments]
        remaining = [100 - p for p in percentages]
        
        plt.figure(figsize=(10, 6))
        plt.bar(treatments, percentages, label='Gene JRIU Variants', color='skyblue')
        plt.bar(treatments, remaining, bottom=percentages, label='Other Variants', color='lightgray')
        plt.title('Percentage of Variants on Gene JRIUs by Treatment', fontsize=14)
        plt.xlabel('Treatment', fontsize=12)
        plt.ylabel('Percentage', fontsize=12)
        plt.xticks(rotation=45)
        plt.legend()
        plt.tight_layout()
        
        percentage_plot = os.path.join(self.plot_dir, 'treatment_percentages.png')
        plt.savefig(percentage_plot, dpi=300)
        plt.close()
        
        plot_files['percentages'] = percentage_plot
        
        return plot_files
    
    def generate_html_report(self, plot_files):
        """Generate comprehensive HTML report with interactive elements"""
        report_file = os.path.join(self.report_dir, "annotation_report.html")
        
        # Prepare data for charts
        treatment_data = {t: stats['jriu_variants'] for t, stats in self.treatment_stats.items()}
        
        gene_data = {}
        for gene_id, stats in self.gene_stats.items():
            sc_gene_id = self.gene_info.get(gene_id, {}).get('sc_gene_id', "unknown")
            gene_name = sc_gene_id.split('Y')[0] if sc_gene_id.startswith('Y') else sc_gene_id
            gene_data[gene_name] = stats['total_variants']
        
        # Prepare gene treatment data
        gene_treatment_data = defaultdict(dict)
        for gene_id, stats in self.gene_stats.items():
            sc_gene_id = self.gene_info.get(gene_id, {}).get('sc_gene_id', "unknown")
            gene_name = sc_gene_id.split('Y')[0] if sc_gene_id.startswith('Y') else sc_gene_id
            
            for treatment, count in stats['treatment_counts'].items():
                gene_treatment_data[gene_name][treatment] = count
        
        # Generate HTML
        with open(report_file, 'w') as f:
            f.write("""<!DOCTYPE html>
<html>
<head>
    <title>JRIU-Based Variant Annotation Report</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/chartjs-plugin-datalabels@2.0.0"></script>
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
    <h1>JRIU-Based Variant Annotation Report</h1>
    <p>Generated on: """ + f"{os.popen('date').read().strip()}" + """</p>
    
    <div class="warning">
        <h3>Important Note on Analysis Approach</h3>
        <p>This analysis uses a JRIU-based approach rather than position-based matching, due to coordinate system differences between the VCF variants and gene annotations. Variants are associated with genes if they occur on the same JRIU ID (genomic segment).</p>
        <p>While this approach cannot determine the precise impact of variants on gene function, it allows for meaningful comparative analysis between treatment groups and identification of segments with high variant density.</p>
    </div>
    
    <div class="tab">
        <button class="tablinks" onclick="openTab(event, 'Summary')" id="defaultOpen">Summary</button>
        <button class="tablinks" onclick="openTab(event, 'Treatments')">Treatments</button>
        <button class="tablinks" onclick="openTab(event, 'Genes')">Genes</button>
        <button class="tablinks" onclick="openTab(event, 'Enrichment')">Enrichment</button>
        <button class="tablinks" onclick="openTab(event, 'Plots')">Plots</button>
    </div>
    
    <div id="Summary" class="tabcontent">
        <h2>Analysis Summary</h2>
        <div class="info">
            <p><strong>Total Files Analyzed:</strong> """ + f"{len(self.sample_stats)}" + """</p>
            <p><strong>Total Genes of Interest:</strong> """ + f"{len(self.gene_info)}" + """</p>
            <p><strong>Total JRIUs with Genes:</strong> """ + f"{len(self.jriu_to_genes)}" + """</p>
            <p><strong>Total Variants Analyzed:</strong> """ + f"{sum(stats['total_variants'] for stats in self.sample_stats.values())}" + """</p>
            <p><strong>Variants on Gene JRIUs:</strong> """ + f"{sum(stats['jriu_variants'] for stats in self.sample_stats.values())}" + """</p>
        </div>
        
        <h3>Key Findings</h3>
        <ul>""")
            
            # Add key findings
            top_genes = sorted(self.gene_stats.items(), key=lambda x: -x[1]['total_variants'])[:3]
            for gene_id, stats in top_genes:
                sc_gene_id = self.gene_info.get(gene_id, {}).get('sc_gene_id', "unknown")
                f.write(f"<li><strong>{sc_gene_id}</strong> has the highest number of associated variants ({stats['total_variants']})</li>\n")
            
            # Top treatments
            top_treatments = sorted(self.treatment_stats.items(), key=lambda x: -x[1]['jriu_variants'])[:3]
            for treatment, stats in top_treatments:
                f.write(f"<li><strong>{treatment}</strong> treatment has {stats['jriu_variants']} variants on gene JRIUs</li>\n")
            
            # Enrichment stats
            enrichment_stats = self.calculate_enrichment_statistics()
            significant_enrichments = []
            for gene_id, treatments in enrichment_stats['gene_enrichment'].items():
                for treatment, stats in treatments.items():
                    if stats['significant'] and stats['enrichment_factor'] > 1.5:
                        sc_gene_id = self.gene_info.get(gene_id, {}).get('sc_gene_id', "unknown")
                        significant_enrichments.append((gene_id, sc_gene_id, treatment, stats['enrichment_factor'], stats['p_value']))
            
            if significant_enrichments:
                top_enrichments = sorted(significant_enrichments, key=lambda x: -x[3])[:3]
                for _, sc_gene_id, treatment, ef, p in top_enrichments:
                    f.write(f"<li><strong>{sc_gene_id}</strong> shows significant enrichment ({ef:.2f}x, p={p:.1e}) in the <strong>{treatment}</strong> treatment</li>\n")
            
            f.write("""</ul>
        
        <h3>Data Distribution</h3>
        <div class="chart-container">
            <canvas id="summaryChart"></canvas>
        </div>
    </div>
    
    <div id="Treatments" class="tabcontent">
        <h2>Treatment Analysis</h2>
        <p>This section shows the distribution of variants across different treatments.</p>
        
        <h3>Treatment Statistics</h3>
        <table>
            <tr>
                <th>Treatment</th>
                <th>Samples</th>
                <th>Total Variants</th>
                <th>Gene JRIU Variants</th>
                <th>Percentage</th>
                <th>Avg. per Sample</th>
            </tr>""")
            
            # Add treatment rows
            for treatment, stats in sorted(self.treatment_stats.items()):
                pct = 100 * stats['jriu_variants'] / stats['total_variants'] if stats['total_variants'] > 0 else 0
                avg = stats['total_variants'] / stats['sample_count'] if stats['sample_count'] > 0 else 0
                f.write(f"""
            <tr>
                <td>{treatment}</td>
                <td>{stats['sample_count']}</td>
                <td>{stats['total_variants']}</td>
                <td>{stats['jriu_variants']}</td>
                <td>{pct:.2f}%</td>
                <td>{avg:.2f}</td>
            </tr>""")
            
            f.write("""
        </table>
        
        <h3>Treatment Variant Distribution</h3>
        <div class="chart-container">
            <canvas id="treatmentChart"></canvas>
        </div>
        
        <h3>Treatment Comparison</h3>
        <div class="chart-container">
            <canvas id="treatmentComparisonChart"></canvas>
        </div>
    </div>
    
    <div id="Genes" class="tabcontent">
        <h2>Gene Analysis</h2>
        <p>This section shows the distribution of variants associated with genes of interest.</p>
        
        <h3>Gene Statistics</h3>
        <table>
            <tr>
                <th>Gene ID</th>
                <th>SC Gene ID</th>
                <th>Name</th>
                <th>Total Variants</th>
                <th>Most Abundant Treatment</th>
            </tr>""")
            
            # Add gene rows
            for gene_id, stats in sorted(self.gene_stats.items(), key=lambda x: -x[1]['total_variants']):
                sc_gene_id = self.gene_info.get(gene_id, {}).get('sc_gene_id', "unknown")
                gene_name = sc_gene_id.split('Y')[0] if sc_gene_id.startswith('Y') else sc_gene_id
                
                # Find most abundant treatment
                top_treatment = max(stats['treatment_counts'].items(), key=lambda x: x[1], default=("None", 0))
                
                f.write(f"""
            <tr>
                <td>{gene_id}</td>
                <td>{sc_gene_id}</td>
                <td>{gene_name}</td>
                <td>{stats['total_variants']}</td>
                <td>{top_treatment[0]} ({top_treatment[1]} variants)</td>
            </tr>""")
            
            f.write("""
        </table>
        
        <h3>Gene Variant Distribution</h3>
        <div class="chart-container">
            <canvas id="geneChart"></canvas>
        </div>
        
        <h3>Gene by Treatment Heatmap</h3>
        <img src="../plots/gene_treatment_heatmap.png" alt="Gene by Treatment Heatmap">
    </div>
    
    <div id="Enrichment" class="tabcontent">
        <h2>Statistical Enrichment Analysis</h2>
        <p>This section shows statistical enrichment of variants on JRIUs containing genes of interest.</p>
        
        <h3>Treatment Enrichment</h3>
        <table>
            <tr>
                <th>Treatment</th>
                <th>Variants on Gene JRIUs</th>
                <th>Expected Count</th>
                <th>Enrichment Factor</th>
                <th>P-Value</th>
                <th>Significant</th>
            </tr>""")
            
            # Add treatment enrichment rows
            enrichment_stats = self.calculate_enrichment_statistics()
            for treatment, stats in sorted(enrichment_stats['treatment_enrichment'].items()):
                actual = self.treatment_stats[treatment]['jriu_variants']
                expected = self.treatment_stats[treatment]['total_variants'] * stats['baseline_rate']
                
                f.write(f"""
            <tr>
                <td>{treatment}</td>
                <td>{actual}</td>
                <td>{expected:.2f}</td>
                <td>{stats['enrichment_factor']:.2f}x</td>
                <td>{stats['p_value']:.4e}</td>
                <td>{"Yes" if stats['significant'] else "No"}</td>
            </tr>""")
            
            f.write("""
        </table>
        
        <h3>Gene Enrichment by Treatment</h3>
        <p>This table shows which genes are enriched for variants in specific treatments.</p>
        
        <table>
            <tr>
                <th>Gene</th>
                <th>Treatment</th>
                <th>Actual Count</th>
                <th>Expected Count</th>
                <th>Enrichment</th>
                <th>P-Value</th>
                <th>Significant</th>
            </tr>""")
            
            # Add gene enrichment rows
            significant_enrichments = []
            for gene_id, treatments in enrichment_stats['gene_enrichment'].items():
                sc_gene_id = self.gene_info.get(gene_id, {}).get('sc_gene_id', "unknown")
                gene_name = sc_gene_id.split('Y')[0] if sc_gene_id.startswith('Y') else sc_gene_id
                
                for treatment, stats in treatments.items():
                    if stats['enrichment_factor'] > 1.1 or stats['significant']:
                        significant_enrichments.append((gene_id, sc_gene_id, gene_name, treatment, stats))
            
            # Sort by significance and enrichment factor
            significant_enrichments.sort(key=lambda x: (-int(x[4]['significant']), -x[4]['enrichment_factor']))
            
            for gene_id, sc_gene_id, gene_name, treatment, stats in significant_enrichments:
                f.write(f"""
            <tr>
                <td>{gene_name} ({sc_gene_id})</td>
                <td>{treatment}</td>
                <td>{stats['actual_count']}</td>
                <td>{stats['expected_count']:.2f}</td>
                <td>{stats['enrichment_factor']:.2f}x</td>
                <td>{stats['p_value']:.4e}</td>
                <td>{"Yes" if stats['significant'] else "No"}</td>
            </tr>""")
            
            f.write("""
        </table>
        
        <h3>Enrichment Visualization</h3>
        <img src="../plots/gene_enrichment.png" alt="Gene Enrichment by Treatment">
    </div>
    
    <div id="Plots" class="tabcontent">
        <h2>Visualization Plots</h2>
        <p>This section contains all visualization plots generated from the analysis.</p>
        
        <h3>Treatment Variant Counts</h3>
        <img src="../plots/treatment_variant_counts.png" alt="Treatment Variant Counts">
        
        <h3>Gene Variant Counts</h3>
        <img src="../plots/gene_variant_counts.png" alt="Gene Variant Counts">
        
        <h3>Gene Treatment Heatmap</h3>
        <img src="../plots/gene_treatment_heatmap.png" alt="Gene Treatment Heatmap">
        
        <h3>Gene Enrichment</h3>
        <img src="../plots/gene_enrichment.png" alt="Gene Enrichment">
        
        <h3>Treatment Percentages</h3>
        <img src="../plots/treatment_percentages.png" alt="Treatment Percentages">
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
        
        // Summary chart
        const summaryCtx = document.getElementById('summaryChart').getContext('2d');
        const summaryChart = new Chart(summaryCtx, {
            type: 'pie',
            data: {
                labels: ['Variants on Gene JRIUs', 'Other Variants'],
                datasets: [{
                    data: [""" + str(sum(stats['jriu_variants'] for stats in self.sample_stats.values())) + ", " + 
                          str(sum(stats['total_variants'] - stats['jriu_variants'] for stats in self.sample_stats.values())) + """],
                    backgroundColor: ['rgba(54, 162, 235, 0.7)', 'rgba(200, 200, 200, 0.7)'],
                    borderWidth: 1
                }]
            }
        });
        
        // Treatment chart
        const treatmentCtx = document.getElementById('treatmentChart').getContext('2d');
        const treatmentChart = new Chart(treatmentCtx, {
            type: 'bar',
            data: {
                labels: """ + str(list(treatment_data.keys())) + """,
                datasets: [{
                    label: 'Variants on Gene JRIUs',
                    data: """ + str(list(treatment_data.values())) + """,
                    backgroundColor: 'rgba(54, 162, 235, 0.7)',
                    borderColor: 'rgba(54, 162, 235, 1)',
                    borderWidth: 1
                }]
            },
            options: {
                scales: {
                    y: {
                        beginAtZero: true,
                        title: {
                            display: true,
                            text: 'Number of Variants'
                        }
                    },
                    x: {
                        title: {
                            display: true,
                            text: 'Treatment'
                        }
                    }
                }
            }
        });
        
        // Treatment comparison chart
        const treatmentComparisonCtx = document.getElementById('treatmentComparisonChart').getContext('2d');
        const treatmentComparisonChart = new Chart(treatmentComparisonCtx, {
            type: 'bar',
            data: {
                labels: """ + str(list(self.treatment_stats.keys())) + """,
                datasets: [{
                    label: 'Gene JRIU Variants',
                    data: """ + str([stats['jriu_variants'] for _, stats in sorted(self.treatment_stats.items())]) + """,
                    backgroundColor: 'rgba(54, 162, 235, 0.7)',
                    borderColor: 'rgba(54, 162, 235, 1)',
                    borderWidth: 1
                }, {
                    label: 'Other Variants',
                    data: """ + str([stats['total_variants'] - stats['jriu_variants'] for _, stats in sorted(self.treatment_stats.items())]) + """,
                    backgroundColor: 'rgba(200, 200, 200, 0.7)',
                    borderColor: 'rgba(200, 200, 200, 1)',
                    borderWidth: 1
                }]
            },
            options: {
                scales: {
                    x: {
                        stacked: true,
                    },
                    y: {
                        stacked: true,
                        beginAtZero: true
                    }
                }
            }
        });
        
        // Gene chart
        const geneCtx = document.getElementById('geneChart').getContext('2d');
        const geneChart = new Chart(geneCtx, {
            type: 'bar',
            data: {
                labels: """ + str(list(gene_data.keys())) + """,
                datasets: [{
                    label: 'Number of Variants',
                    data: """ + str(list(gene_data.values())) + """,
                    backgroundColor: 'rgba(75, 192, 192, 0.7)',
                    borderColor: 'rgba(75, 192, 192, 1)',
                    borderWidth: 1
                }]
            },
            options: {
                scales: {
                    y: {
                        beginAtZero: true,
                        title: {
                            display: true,
                            text: 'Number of Variants'
                        }
                    },
                    x: {
                        title: {
                            display: true,
                            text: 'Gene'
                        }
                    }
                }
            }
        });
    </script>
</body>
</html>""")
        
        return report_file
    
    def run_complete_analysis(self, vcf_files):
        """Run the complete analysis pipeline"""
        # Process all VCF files
        results = self.process_all_vcfs(vcf_files)
        
        # Generate text reports
        report_files = self.generate_text_reports()
        print("\nGenerated text reports:")
        for report_type, report_file in report_files.items():
            print(f"  {report_type}: {report_file}")
        
        # Generate plots
        plot_files = self.generate_plots()
        print("\nGenerated plots:")
        for plot_type, plot_file in plot_files.items():
            print(f"  {plot_type}: {plot_file}")
        
        # Generate HTML report
        html_report = self.generate_html_report(plot_files)
        print(f"\nComplete HTML report: {html_report}")
        
        return {
            'results': results,
            'report_files': report_files,
            'plot_files': plot_files,
            'html_report': html_report
        }

def main():
    parser = argparse.ArgumentParser(
        description='Enhanced JRIU-based annotation system for yeast variants',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--vcf', required=True, nargs='+', help='VCF file(s) to annotate')
    parser.add_argument('--genes-of-interest', required=True, help='Genes of interest file')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    
    args = parser.parse_args()
    
    # Create annotator
    annotator = EnhancedJRIUAnnotator(
        genes_of_interest_file=args.genes_of_interest,
        output_dir=args.output_dir
    )
    
    # Run analysis
    analysis_results = annotator.run_complete_analysis(args.vcf)
    
    print("\nAnalysis complete!")
    print(f"Results are available in {args.output_dir}")

if __name__ == '__main__':
    main()