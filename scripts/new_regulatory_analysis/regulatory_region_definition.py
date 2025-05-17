#!/usr/bin/env python3
# /Users/zakiralibhai/Documents/GitHub/Yeast_MSA/scripts/new_regulatory_analysis/regulatory_region_definition.py

"""
Enhanced Regulatory Region Definition Framework for the Yeast MSA project.

This script defines regulatory regions with higher precision based on yeast genomic features,
incorporating known regulatory elements from literature and databases, and
creating a comprehensive map of potential regulatory regions around genes.

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

class RegulatoryRegionDefinitionFramework:
    """Framework for defining and analyzing regulatory regions with high precision"""
    
    def __init__(self, regulatory_features_file, output_dir, gene_mapping_file, erg_genes_file, osh_gene_file=None):
        """Initialize the framework with necessary files and configurations"""
        # Setup directories
        self.output_dir = self._ensure_dir(output_dir)
        self.data_dir = self._ensure_dir(os.path.join(output_dir, 'data'))
        self.plot_dir = self._ensure_dir(os.path.join(output_dir, 'plots'))
        
        # Load input files
        print("Loading input files...")
        
        # Load regulatory features
        try:
            with open(regulatory_features_file, 'r') as f:
                self.regulatory_features = json.load(f)
            print(f"Loaded regulatory features with {len(self.regulatory_features)} categories")
        except Exception as e:
            print(f"Error loading regulatory features file: {e}")
            sys.exit(1)
        
        # Load gene mapping and ERG genes
        self.gene_mapping = self._load_file(gene_mapping_file, "gene mapping")
        self.erg_genes = self._load_file(erg_genes_file, "ERG genes")
        
        # Load OSH genes if file provided
        self.osh_genes = None
        if osh_gene_file:
            self.osh_genes = self._load_file(osh_gene_file, "OSH genes")
        else:
            print("OSH gene file not provided, attempting to load from results/osh_analysis/osh_gene_summary.tsv")
            osh_default_path = "/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/osh_analysis/osh_gene_summary.tsv"
            if os.path.exists(osh_default_path):
                self.osh_genes = self._load_file(osh_default_path, "OSH genes")
                print(f"Successfully loaded OSH genes with {len(self.osh_genes)} entries")
            else:
                print("Default OSH gene file not found, proceeding without OSH gene information")
        
        # Extract key region information from loaded features
        self._extract_region_definitions()
        
        # Store analysis results
        self.regions_by_gene = {}
        self.gene_regulatory_map = {}
        self.region_statistics = {}
    
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
            if description in ["gene mapping", "regulatory features"]:
                print("Critical file missing. Exiting.")
                sys.exit(1)
            return pd.DataFrame()
    
    def _extract_region_definitions(self):
        """Extract region definitions from loaded regulatory features"""
        # Initialize region dictionaries
        self.upstream_regions = {}
        self.downstream_regions = {}
        self.special_elements = {}
        self.conservation_zones = {}
        
        # Extract region definitions if available
        if 'upstream_regions' in self.regulatory_features:
            for region, props in self.regulatory_features['upstream_regions'].items():
                if 'range' in props:
                    self.upstream_regions[region] = (props['range'][0], props['range'][1])
        
        if 'downstream_regions' in self.regulatory_features:
            for region, props in self.regulatory_features['downstream_regions'].items():
                if 'range' in props:
                    self.downstream_regions[region] = (props['range'][0], props['range'][1])
        
        if 'special_elements' in self.regulatory_features:
            for element, props in self.regulatory_features['special_elements'].items():
                if 'region' in props and 'consensus' in props:
                    self.special_elements[element] = {
                        'region': (props['region'][0], props['region'][1]),
                        'consensus': props['consensus']
                    }
        
        if 'conservation_zones' in self.regulatory_features:
            for zone, props in self.regulatory_features['conservation_zones'].items():
                if 'range' in props:
                    self.conservation_zones[zone] = (props['range'][0], props['range'][1])
        
        # Log the extracted regions
        print(f"\nExtracted {len(self.upstream_regions)} upstream regions")
        print(f"Extracted {len(self.downstream_regions)} downstream regions")
        print(f"Extracted {len(self.special_elements)} special elements")
        print(f"Extracted {len(self.conservation_zones)} conservation zones")
    
    def create_gene_regulatory_maps(self):
        """Create detailed regulatory maps for all genes"""
        print("\nCreating detailed regulatory maps for genes...")
        
        # Get gene coordinates
        gene_coords = {}
        for _, row in self.gene_mapping.iterrows():
            if 'w303_gene_id' in row and 'start' in row and 'end' in row and 'strand' in row:
                gene_id = row['w303_gene_id']
                scaffold = row.get('w303_scaffold', '')
                start = row['start']
                end = row['end']
                strand = row['strand']
                
                gene_coords[gene_id] = {
                    'scaffold': scaffold,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'length': end - start + 1
                }
        
        print(f"Mapped coordinates for {len(gene_coords)} genes")
        
        # Extract ERG genes if available
        erg_gene_ids = set()
        if 'w303_gene_id' in self.erg_genes.columns:
            erg_gene_ids = set(self.erg_genes['w303_gene_id'])
            print(f"Identified {len(erg_gene_ids)} ERG genes")
        
        # Extract OSH genes if available
        osh_gene_ids = set()
        if self.osh_genes is not None and 'W303_Gene_ID' in self.osh_genes.columns:
            osh_gene_ids = set(self.osh_genes['W303_Gene_ID'])
            print(f"Identified {len(osh_gene_ids)} OSH genes")
        
        # Create detailed regulatory maps
        gene_regulatory_map = {}
        
        for gene_id, coords in gene_coords.items():
            # Basic gene info
            gene_map = {
                'gene_id': gene_id,
                'scaffold': coords['scaffold'],
                'start': coords['start'],
                'end': coords['end'],
                'strand': coords['strand'],
                'length': coords['length'],
                'is_erg': gene_id in erg_gene_ids,
                'is_osh': gene_id in osh_gene_ids,
                'regulatory_regions': {},
                'conservation_zones': {}
            }
            
            # Define reference position based on strand (for regulatory distances)
            if coords['strand'] == '+':
                ref_pos = coords['start']  # Use start (5' end) for + strand
            else:
                ref_pos = coords['end']    # Use end (5' end) for - strand
            
            # Map upstream regions
            for region, (min_dist, max_dist) in self.upstream_regions.items():
                if coords['strand'] == '+':
                    # For + strand, upstream is before the start
                    start_pos = ref_pos - max_dist
                    end_pos = ref_pos - min_dist
                else:
                    # For - strand, upstream is after the end
                    start_pos = ref_pos + min_dist
                    end_pos = ref_pos + max_dist
                
                gene_map['regulatory_regions'][region] = {
                    'start': start_pos,
                    'end': end_pos,
                    'length': end_pos - start_pos + 1,
                    'type': 'upstream'
                }
            
            # Map downstream regions
            for region, (min_dist, max_dist) in self.downstream_regions.items():
                if coords['strand'] == '+':
                    # For + strand, downstream is after the end
                    start_pos = coords['end'] + min_dist
                    end_pos = coords['end'] + max_dist
                else:
                    # For - strand, downstream is before the start
                    start_pos = coords['start'] - max_dist
                    end_pos = coords['start'] - min_dist
                
                gene_map['regulatory_regions'][region] = {
                    'start': start_pos,
                    'end': end_pos,
                    'length': end_pos - start_pos + 1,
                    'type': 'downstream'
                }
            
            # Map conservation zones for ERG and OSH genes
            if gene_id in erg_gene_ids or gene_id in osh_gene_ids:
                for zone, (min_dist, max_dist) in self.conservation_zones.items():
                    if zone == 'core_zone':
                        # Core zone is the gene itself
                        gene_map['conservation_zones'][zone] = {
                            'start': coords['start'],
                            'end': coords['end'],
                            'length': coords['length'],
                            'type': 'conservation'
                        }
                    else:
                        # Other zones are defined by distance from gene boundaries
                        # For buffer, intermediate, and satellite zones
                        # We calculate both upstream and downstream regions
                        
                        # Upstream of gene (5' end)
                        if coords['strand'] == '+':
                            up_start = coords['start'] - max_dist
                            up_end = coords['start'] - min_dist
                        else:
                            up_start = coords['end'] + min_dist
                            up_end = coords['end'] + max_dist
                        
                        # Downstream of gene (3' end)
                        if coords['strand'] == '+':
                            down_start = coords['end'] + min_dist
                            down_end = coords['end'] + max_dist
                        else:
                            down_start = coords['start'] - max_dist
                            down_end = coords['start'] - min_dist
                        
                        gene_map['conservation_zones'][f"{zone}_upstream"] = {
                            'start': up_start,
                            'end': up_end,
                            'length': up_end - up_start + 1,
                            'type': 'conservation'
                        }
                        
                        gene_map['conservation_zones'][f"{zone}_downstream"] = {
                            'start': down_start,
                            'end': down_end,
                            'length': down_end - down_start + 1,
                            'type': 'conservation'
                        }
            
            # Store the gene map
            gene_regulatory_map[gene_id] = gene_map
        
        # Save the gene regulatory map
        self.gene_regulatory_map = gene_regulatory_map
        
        # Save to JSON
        map_file = os.path.join(self.data_dir, 'gene_regulatory_map.json')
        with open(map_file, 'w') as f:
            # Convert map to a more JSON-friendly format
            json_map = {}
            for gene_id, gene_info in gene_regulatory_map.items():
                json_map[gene_id] = {
                    k: (v if not isinstance(v, np.int64) else int(v)) 
                    for k, v in gene_info.items()
                    if k not in ['regulatory_regions', 'conservation_zones']
                }
                
                # Handle nested dictionaries
                json_map[gene_id]['regulatory_regions'] = {
                    region: {
                        k: (v if not isinstance(v, np.int64) else int(v))
                        for k, v in region_info.items()
                    }
                    for region, region_info in gene_info.get('regulatory_regions', {}).items()
                }
                
                json_map[gene_id]['conservation_zones'] = {
                    zone: {
                        k: (v if not isinstance(v, np.int64) else int(v))
                        for k, v in zone_info.items()
                    }
                    for zone, zone_info in gene_info.get('conservation_zones', {}).items()
                }
            
            json.dump(json_map, f, indent=2)
        
        print(f"Created regulatory maps for {len(gene_regulatory_map)} genes")
        print(f"Saved gene regulatory map to {map_file}")
        
        return gene_regulatory_map
    
    def analyze_region_statistics(self):
        """Analyze statistics of various regulatory regions"""
        print("\nAnalyzing statistics of regulatory regions...")
        
        if not self.gene_regulatory_map:
            print("No gene regulatory map available. Run create_gene_regulatory_maps() first.")
            return
        
        # Count genes by type
        gene_types = Counter()
        gene_types['total'] = len(self.gene_regulatory_map)
        gene_types['erg'] = sum(1 for info in self.gene_regulatory_map.values() if info.get('is_erg', False))
        gene_types['osh'] = sum(1 for info in self.gene_regulatory_map.values() if info.get('is_osh', False))
        gene_types['other'] = gene_types['total'] - gene_types['erg'] - gene_types['osh']
        
        # Region counts
        region_counts = defaultdict(int)
        for gene_info in self.gene_regulatory_map.values():
            for region in gene_info.get('regulatory_regions', {}).keys():
                region_counts[region] += 1
        
        # Conservation zone counts
        zone_counts = defaultdict(int)
        for gene_info in self.gene_regulatory_map.values():
            for zone in gene_info.get('conservation_zones', {}).keys():
                zone_counts[zone] += 1
        
        # Region size statistics
        region_sizes = defaultdict(list)
        for gene_info in self.gene_regulatory_map.values():
            for region, region_info in gene_info.get('regulatory_regions', {}).items():
                region_sizes[region].append(region_info.get('length', 0))
        
        # Calculate region size statistics
        region_stats = {}
        for region, sizes in region_sizes.items():
            if sizes:
                region_stats[region] = {
                    'count': len(sizes),
                    'mean_size': np.mean(sizes),
                    'median_size': np.median(sizes),
                    'min_size': np.min(sizes),
                    'max_size': np.max(sizes)
                }
        
        # Store statistics
        statistics = {
            'gene_types': dict(gene_types),
            'region_counts': dict(region_counts),
            'zone_counts': dict(zone_counts),
            'region_stats': region_stats
        }
        
        self.region_statistics = statistics
        
        # Save statistics to JSON
        stats_file = os.path.join(self.data_dir, 'regulatory_region_statistics.json')
        with open(stats_file, 'w') as f:
            # Convert statistics to JSON-serializable format
            json_stats = {}
            for category, cat_stats in statistics.items():
                if isinstance(cat_stats, dict):
                    # Handle nested dictionaries
                    if category == 'region_stats':
                        json_stats[category] = {
                            region: {
                                k: (float(v) if isinstance(v, np.number) else v)
                                for k, v in region_info.items()
                            }
                            for region, region_info in cat_stats.items()
                        }
                    else:
                        json_stats[category] = {
                            k: (int(v) if isinstance(v, np.number) else v)
                            for k, v in cat_stats.items()
                        }
                else:
                    json_stats[category] = cat_stats
            
            json.dump(json_stats, f, indent=2)
        
        print(f"Analyzed {gene_types['total']} genes")
        print(f"Found {gene_types['erg']} ERG genes and {gene_types['osh']} OSH genes")
        print(f"Analyzed {len(region_counts)} regulatory region types")
        print(f"Analyzed {len(zone_counts)} conservation zone types")
        print(f"Saved region statistics to {stats_file}")
        
        # Print summary statistics
        print("\nSummary of regulatory region sizes:")
        for region, stats in sorted(region_stats.items(), key=lambda x: x[1]['count'], reverse=True)[:5]:
            print(f"  {region}: {stats['count']} instances, mean size = {stats['mean_size']:.1f}bp")
        
        return statistics
    
    def visualize_regulatory_regions(self):
        """Generate visualizations of regulatory regions"""
        print("\nGenerating visualizations of regulatory regions...")
        
        if not self.region_statistics:
            print("No region statistics available. Run analyze_region_statistics() first.")
            return
        
        # Set style
        sns.set(style="whitegrid")
        
        # 1. Gene Type Distribution
        gene_types = self.region_statistics['gene_types']
        plt.figure(figsize=(10, 6))
        plt.bar(gene_types.keys(), gene_types.values(), color=['blue', 'green', 'red', 'gray'])
        plt.title('Gene Type Distribution', fontsize=14)
        plt.xlabel('Gene Type', fontsize=12)
        plt.ylabel('Count', fontsize=12)
        plt.tight_layout()
        plt.savefig(os.path.join(self.plot_dir, 'gene_type_distribution.png'), dpi=300)
        plt.close()
        
        # 2. Regulatory Region Count Distribution
        region_counts = self.region_statistics['region_counts']
        plt.figure(figsize=(12, 6))
        plt.bar(region_counts.keys(), region_counts.values(), color='skyblue')
        plt.title('Regulatory Region Distribution', fontsize=14)
        plt.xlabel('Region Type', fontsize=12)
        plt.ylabel('Count', fontsize=12)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(os.path.join(self.plot_dir, 'regulatory_region_counts.png'), dpi=300)
        plt.close()
        
        # 3. Conservation Zone Distribution
        zone_counts = self.region_statistics['zone_counts']
        if zone_counts:
            plt.figure(figsize=(12, 6))
            plt.bar(zone_counts.keys(), zone_counts.values(), color='lightgreen')
            plt.title('Conservation Zone Distribution', fontsize=14)
            plt.xlabel('Zone Type', fontsize=12)
            plt.ylabel('Count', fontsize=12)
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(os.path.join(self.plot_dir, 'conservation_zone_counts.png'), dpi=300)
            plt.close()
        
        # 4. Region Size Distribution
        region_stats = self.region_statistics['region_stats']
        region_types = []
        mean_sizes = []
        min_sizes = []
        max_sizes = []
        
        for region, stats in region_stats.items():
            region_types.append(region)
            mean_sizes.append(stats['mean_size'])
            min_sizes.append(stats['min_size'])
            max_sizes.append(stats['max_size'])
        
        # Convert to DataFrame for plotting
        df = pd.DataFrame({
            'Region': region_types,
            'Mean Size': mean_sizes,
            'Min Size': min_sizes,
            'Max Size': max_sizes
        })
        
        # Sort by mean size
        df = df.sort_values('Mean Size', ascending=False)
        
        plt.figure(figsize=(14, 8))
        sns.barplot(x='Region', y='Mean Size', data=df, color='darkblue')
        plt.title('Mean Size of Regulatory Regions', fontsize=14)
        plt.xlabel('Region Type', fontsize=12)
        plt.ylabel('Mean Size (bp)', fontsize=12)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(os.path.join(self.plot_dir, 'regulatory_region_sizes.png'), dpi=300)
        plt.close()
        
        # 5. ERG vs OSH Gene Region Comparison
        if self.region_statistics['gene_types']['erg'] > 0 and self.region_statistics['gene_types']['osh'] > 0:
            # Extract ERG and OSH genes
            erg_genes = [info for gene_id, info in self.gene_regulatory_map.items() if info.get('is_erg', False)]
            osh_genes = [info for gene_id, info in self.gene_regulatory_map.items() if info.get('is_osh', False)]
            
            # Count regions in ERG and OSH genes
            erg_regions = Counter()
            osh_regions = Counter()
            
            for gene in erg_genes:
                for region in gene.get('regulatory_regions', {}).keys():
                    erg_regions[region] += 1
            
            for gene in osh_genes:
                for region in gene.get('regulatory_regions', {}).keys():
                    osh_regions[region] += 1
            
            # Prepare data for plotting
            all_regions = sorted(set(list(erg_regions.keys()) + list(osh_regions.keys())))
            erg_counts = [erg_regions.get(region, 0) for region in all_regions]
            osh_counts = [osh_regions.get(region, 0) for region in all_regions]
            
            # Plot comparison
            plt.figure(figsize=(14, 8))
            width = 0.35
            x = np.arange(len(all_regions))
            
            plt.bar(x - width/2, erg_counts, width, label='ERG Genes', color='green')
            plt.bar(x + width/2, osh_counts, width, label='OSH Genes', color='purple')
            
            plt.title('Regulatory Region Comparison: ERG vs OSH Genes', fontsize=14)
            plt.xlabel('Region Type', fontsize=12)
            plt.ylabel('Count', fontsize=12)
            plt.xticks(x, all_regions, rotation=45, ha='right')
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(self.plot_dir, 'erg_osh_region_comparison.png'), dpi=300)
            plt.close()
        
        print(f"Generated visualizations in {self.plot_dir}")
    
    def generate_region_definition_report(self):
        """Generate a comprehensive report of the regulatory region definitions"""
        print("\nGenerating comprehensive regulatory region report...")
        
        # Ensure we have necessary data
        if not self.gene_regulatory_map:
            print("No gene regulatory map available. Run create_gene_regulatory_maps() first.")
            return
        
        if not self.region_statistics:
            print("No region statistics available. Run analyze_region_statistics() first.")
            return
        
        # Create report content
        report = ["# Enhanced Regulatory Region Definition Report\n"]
        
        # 1. Overview section
        report.append("## Overview\n")
        report.append(f"Total genes analyzed: {self.region_statistics['gene_types']['total']}")
        report.append(f"ERG pathway genes: {self.region_statistics['gene_types']['erg']}")
        report.append(f"OSH family genes: {self.region_statistics['gene_types']['osh']}")
        report.append(f"Other genes: {self.region_statistics['gene_types']['other']}\n")
        
        # 2. Regulatory Region Definitions
        report.append("## Regulatory Region Definitions\n")
        
        # Upstream regions
        report.append("### Upstream Regulatory Regions\n")
        report.append("| Region | Distance Range (bp) | Description |")
        report.append("|--------|-------------------|-------------|")
        for region, (min_dist, max_dist) in sorted(self.upstream_regions.items()):
            desc = self.regulatory_features.get('upstream_regions', {}).get(region, {}).get('description', '')
            report.append(f"| {region} | {min_dist} - {max_dist} | {desc} |")
        report.append("")
        
        # Downstream regions
        report.append("### Downstream Regulatory Regions\n")
        report.append("| Region | Distance Range (bp) | Description |")
        report.append("|--------|-------------------|-------------|")
        for region, (min_dist, max_dist) in sorted(self.downstream_regions.items()):
            desc = self.regulatory_features.get('downstream_regions', {}).get(region, {}).get('description', '')
            report.append(f"| {region} | {min_dist} - {max_dist} | {desc} |")
        report.append("")
        
        # Special elements
        report.append("### Special Regulatory Elements\n")
        report.append("| Element | Distance Range (bp) | Consensus Sequence | Description |")
        report.append("|---------|-------------------|-------------------|-------------|")
        for element, props in sorted(self.special_elements.items()):
            min_dist, max_dist = props['region']
            consensus = props['consensus']
            desc = self.regulatory_features.get('special_elements', {}).get(element, {}).get('description', '')
            report.append(f"| {element} | {min_dist} - {max_dist} | {consensus} | {desc} |")
        report.append("")
        
        # Conservation zones
        report.append("### Conservation Zones\n")
        report.append("| Zone | Distance Range (bp) | Description |")
        report.append("|------|-------------------|-------------|")
        for zone, (min_dist, max_dist) in sorted(self.conservation_zones.items()):
            desc = self.regulatory_features.get('conservation_zones', {}).get(zone, {}).get('description', '')
            report.append(f"| {zone} | {min_dist} - {max_dist} | {desc} |")
        report.append("")
        
        # 3. Region Statistics
        report.append("## Region Statistics\n")
        
        # Region count statistics
        report.append("### Regulatory Region Counts\n")
        report.append("| Region | Count | Percentage |")
        report.append("|--------|-------|------------|")
        total_regions = sum(self.region_statistics['region_counts'].values())
        for region, count in sorted(self.region_statistics['region_counts'].items(), key=lambda x: x[1], reverse=True):
            pct = (count / total_regions) * 100 if total_regions > 0 else 0
            report.append(f"| {region} | {count} | {pct:.1f}% |")
        report.append("")
        
        # Region size statistics
        report.append("### Regulatory Region Size Statistics\n")
        report.append("| Region | Count | Mean Size (bp) | Median Size (bp) | Min Size (bp) | Max Size (bp) |")
        report.append("|--------|-------|---------------|-----------------|--------------|---------------|")
        for region, stats in sorted(self.region_statistics.get('region_stats', {}).items(), key=lambda x: x[1]['count'], reverse=True):
            report.append(f"| {region} | {stats['count']} | {stats['mean_size']:.1f} | {stats['median_size']:.1f} | {stats['min_size']} | {stats['max_size']} |")
        report.append("")
        
        # 4. ERG and OSH Gene Analysis
        report.append("## ERG and OSH Gene Analysis\n")
        
        # ERG gene regions
        erg_genes = [info for gene_id, info in self.gene_regulatory_map.items() if info.get('is_erg', False)]
        if erg_genes:
            report.append("### ERG Gene Regulatory Profile\n")
            report.append(f"Total ERG genes: {len(erg_genes)}\n")
            
            # Count regions in ERG genes
            erg_regions = Counter()
            for gene in erg_genes:
                for region in gene.get('regulatory_regions', {}).keys():
                    erg_regions[region] += 1
            
            report.append("#### Regulatory Region Distribution in ERG Genes\n")
            report.append("| Region | Count | Percentage |")
            report.append("|--------|-------|------------|")
            total_erg_regions = sum(erg_regions.values())
            for region, count in sorted(erg_regions.items(), key=lambda x: x[1], reverse=True):
                pct = (count / total_erg_regions) * 100 if total_erg_regions > 0 else 0
                report.append(f"| {region} | {count} | {pct:.1f}% |")
            report.append("")
        
        # OSH gene regions
        osh_genes = [info for gene_id, info in self.gene_regulatory_map.items() if info.get('is_osh', False)]
        if osh_genes:
            report.append("### OSH Gene Regulatory Profile\n")
            report.append(f"Total OSH genes: {len(osh_genes)}\n")
            
            # Count regions in OSH genes
            osh_regions = Counter()
            for gene in osh_genes:
                for region in gene.get('regulatory_regions', {}).keys():
                    osh_regions[region] += 1
            
            report.append("#### Regulatory Region Distribution in OSH Genes\n")
            report.append("| Region | Count | Percentage |")
            report.append("|--------|-------|------------|")
            total_osh_regions = sum(osh_regions.values())
            for region, count in sorted(osh_regions.items(), key=lambda x: x[1], reverse=True):
                pct = (count / total_osh_regions) * 100 if total_osh_regions > 0 else 0
                report.append(f"| {region} | {count} | {pct:.1f}% |")
            report.append("")
        
        # 5. Biological Interpretations
        report.append("## Biological Interpretations\n")
        
        # Add interpretations based on the findings
        report.append("### Yeast Promoter Architecture\n")
        report.append("The yeast promoter architecture is characterized by:")
        report.append("- **Core promoter regions** (0-150bp upstream): Essential for basic transcription machinery assembly")
        report.append("- **TATA box regions** (40-120bp upstream): Present in ~20% of yeast genes, associated with regulated genes")
        report.append("- **Upstream Activating Sequences (UAS)**: Enhancer-like elements located 150-1500bp upstream")
        report.append("- **Stress Response Elements (STRE)**: Often found in the UAS regions of stress-responsive genes")
        report.append("- **Short 5' and 3' UTRs**: Typical of yeast genes, with important roles in translation regulation\n")
        
        report.append("### The Four-Zone Conservation Architecture\n")
        report.append("The genomic organization around ERG and OSH genes reveals a hierarchical conservation pattern:")
        report.append("1. **Core Zone**: The genes themselves show absolute conservation, reflecting their essential functions")
        report.append("2. **Buffer Zone** (0-5kb): Minimal variation, primarily regulatory adjustments")
        report.append("3. **Intermediate Zone** (5-50kb): Moderate variation, including pathway modulators")
        report.append("4. **Satellite Zone** (50-100kb): Higher variation, enabling adaptive flexibility while maintaining core functions\n")
        
        report.append("### ERG-OSH Regulatory Relationships\n")
        report.append("The regulatory architecture of ERG (ergosterol biosynthesis) and OSH (oxysterol binding homology) genes suggests:")
        report.append("- Shared regulatory elements controlling both sterol synthesis and transport")
        report.append("- Similar conservation patterns indicating coordinated evolution")
        report.append("- Specialized regulatory elements for oxygen sensing and sterol homeostasis")
        report.append("- Potential co-regulation through chromatin organization\n")
        
        # 6. Implications for Adaptation
        report.append("## Implications for Adaptation\n")
        report.append("The regulatory architecture identified has important implications for adaptation mechanisms:")
        report.append("- **Regulatory flexibility with functional conservation**: Adaptation through changes in gene expression rather than protein structure")
        report.append("- **Zone-specific adaptation**: Different adaptation strategies depending on genomic distance from essential genes")
        report.append("- **Coordinated regulation**: Changes in sterol synthesis and transport systems are likely coordinated")
        report.append("- **Environmental responsiveness**: Specialized regulatory elements for different environmental stresses")
        report.append("- **Transcription factor binding site variation**: Subtle changes in TF binding affinity rather than complete gain/loss of binding sites")
        
        # Write report to file
        report_path = os.path.join(self.output_dir, 'regulatory_region_definition_report.md')
        with open(report_path, 'w') as f:
            f.write('\n'.join(report))
        
        print(f"Generated comprehensive report: {report_path}")
        
        return report
    
    def run_analysis(self):
        """Run the complete regulatory region definition pipeline"""
        print("\n=== Starting Enhanced Regulatory Region Definition Analysis ===\n")
        
        # Step 1: Create gene regulatory maps
        self.create_gene_regulatory_maps()
        
        # Step 2: Analyze region statistics
        self.analyze_region_statistics()
        
        # Step 3: Visualize regulatory regions
        self.visualize_regulatory_regions()
        
        # Step 4: Generate comprehensive report
        self.generate_region_definition_report()
        
        print("\n=== Enhanced Regulatory Region Definition Analysis Complete ===")
        
        # Return the gene regulatory map for potential use by other scripts
        return self.gene_regulatory_map


def main():
    """Main function to parse arguments and run the analysis"""
    parser = argparse.ArgumentParser(description='Enhanced Regulatory Region Definition Framework')
    
    # Required arguments
    parser.add_argument('--output-dir', default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/new_regulatory_analysis',
                        help='Output directory for results')
    parser.add_argument('--regulatory-features', default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/scripts/new_regulatory_analysis/data/regulatory_features.json',
                        help='JSON file containing regulatory feature definitions')
    parser.add_argument('--gene-mapping', default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/reference/gene_mapping_full.tsv',
                        help='Gene mapping file with coordinates')
    parser.add_argument('--erg-genes', default='/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/reference/genes_of_interest_mapping.tsv',
                        help='ERG genes mapping file')
    
    # Optional arguments
    parser.add_argument('--osh-genes', default=None,
                        help='OSH genes mapping file (optional, will try to use default if not provided)')
    
    args = parser.parse_args()
    
    # Initialize the framework
    framework = RegulatoryRegionDefinitionFramework(
        args.regulatory_features,
        args.output_dir,
        args.gene_mapping,
        args.erg_genes,
        args.osh_genes
    )
    
    # Run the analysis
    framework.run_analysis()


if __name__ == "__main__":
    main()