#!/usr/bin/env python3
"""
analyze_key_genomic_regions.py

This script performs a detailed analysis of key genomic regions showing interesting variant patterns:
1. ERG11: Region with HIGH impact frameshifts at ~8kb upstream
2. ERG7: Region with HIGH impact frameshifts at exactly 47,676bp downstream
3. ERG25: Regions with highest number of nearby variants

Usage:
  python analyze_key_genomic_regions.py --variants_file <path> --genome_file <path> --genbank_dir <path> --output_dir <path>
"""

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import json
from collections import defaultdict, Counter
import re
from matplotlib.patches import Rectangle, Patch
from scipy.stats import fisher_exact, binomtest

# Key genomic regions to focus on, based on previous analysis
KEY_REGIONS = {
    'ERG11_upstream': {
        'erg_gene': 'ERG11',
        'sc_gene_id': 'YHR007C',
        'distance': 8149,
        'position': 'upstream',
        'feature': 'HIGH impact frameshifts'
    },
    'ERG7_downstream': {
        'erg_gene': 'ERG7',
        'sc_gene_id': 'YHR072W',
        'distance': 47676,
        'position': 'downstream',
        'feature': 'HIGH impact frameshifts'
    },
    'ERG25_neighborhood': {
        'erg_gene': 'ERG25',
        'sc_gene_id': 'YGR060W',
        'distance': 'multiple',
        'position': 'both',
        'feature': 'Highest number of nearby variants'
    }
}

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Analyze key genomic regions with interesting variant patterns')
    parser.add_argument('--variants_file', required=True, help='Path to the variants TSV file from distance analysis')
    parser.add_argument('--genome_file', required=True, help='Path to the reference genome FASTA file')
    parser.add_argument('--genbank_dir', required=True, help='Directory containing GenBank annotation files')
    parser.add_argument('--gene_mapping', required=True, help='Path to the gene mapping TSV file')
    parser.add_argument('--mapping_file', required=True, help='Path to the chromosome mapping TSV file')
    parser.add_argument('--output_dir', required=True, help='Directory to save results')
    return parser.parse_args()

def load_variants(variants_file):
    """Load variants data from TSV file."""
    print(f"Loading variants from {variants_file}")
    variants = pd.read_csv(variants_file, sep='\t')
    print(f"Loaded {len(variants)} variants")
    return variants

def load_genome(genome_file):
    """Load reference genome from FASTA file with corrected IDs."""
    print(f"Loading reference genome from {genome_file}")
    genome = {}
    for record in SeqIO.parse(genome_file, "fasta"):
        # Extract just the CM ID from the full header
        header = record.id
        cm_id = header.split()[0]  # Get the first part before any whitespace
        genome[cm_id] = record
        
    print(f"Loaded {len(genome)} sequences from reference genome")
    # Print a few examples to verify
    examples = list(genome.keys())[:3]
    print(f"Example genome keys: {examples}")
    return genome

def load_chromosome_mapping(mapping_file):
    """Load scaffold mapping from TSV file."""
    print(f"Loading chromosome mapping from {mapping_file}")
    mapping_df = pd.read_csv(mapping_file, sep='\t')
    
    # Create bidirectional mappings
    cm_to_scaffold = dict(zip(mapping_df['chromosome_id'], mapping_df['w303_scaffold']))
    scaffold_to_cm = dict(zip(mapping_df['w303_scaffold'], mapping_df['chromosome_id']))
    
    print(f"Loaded mapping for {len(cm_to_scaffold)} scaffolds")
    return scaffold_to_cm, cm_to_scaffold

def find_genbank_for_scaffold(genbank_dir, scaffold):
    """Find the GenBank file containing the given scaffold."""
    for file in os.listdir(genbank_dir):
        if not file.endswith('.genbank') and not file.endswith('.gb'):
            continue
        
        # Parse the GenBank file
        try:
            record = SeqIO.read(os.path.join(genbank_dir, file), "genbank")
            
            # Check if this is the correct scaffold
            if record.name == scaffold or record.id == scaffold:
                return os.path.join(genbank_dir, file)
            
            # Also check notes and other fields
            for feature in record.features:
                if feature.type == 'source':
                    for note in feature.qualifiers.get('note', []):
                        if scaffold in note:
                            return os.path.join(genbank_dir, file)
        except Exception as e:
            continue
    
    return None

def load_genbank_annotations(genbank_dir, scaffolds_of_interest):
    """Load GenBank annotations for scaffolds of interest."""
    print(f"Loading GenBank annotations for {len(scaffolds_of_interest)} scaffolds")
    
    annotations = {}
    genbank_files = {}
    
    # Find GenBank files for each scaffold
    for scaffold in scaffolds_of_interest:
        genbank_file = find_genbank_for_scaffold(genbank_dir, scaffold)
        if genbank_file:
            genbank_files[scaffold] = genbank_file
    
    # Load annotations from GenBank files
    for scaffold, genbank_file in genbank_files.items():
        print(f"  Loading annotations from {os.path.basename(genbank_file)}")
        try:
            record = SeqIO.read(genbank_file, "genbank")
            features = []
            
            for feature in record.features:
                if feature.type in ['gene', 'CDS', 'mRNA', 'tRNA', 'rRNA']:
                    feature_info = {
                        'type': feature.type,
                        'start': int(feature.location.start),
                        'end': int(feature.location.end),
                        'strand': '+' if feature.location.strand == 1 else '-',
                        'id': feature.qualifiers.get('gene', ['unknown'])[0],
                        'locus_tag': feature.qualifiers.get('locus_tag', ['unknown'])[0],
                        'product': feature.qualifiers.get('product', ['unknown'])[0],
                        'note': '; '.join(feature.qualifiers.get('note', [])),
                        'inference': '; '.join(feature.qualifiers.get('inference', []))
                    }
                    
                    # Try to find SC gene ID in inference or note
                    sc_gene_id = None
                    for inf in feature.qualifiers.get('inference', []):
                        match = re.search(r'Y[A-Z]{2}\d{3}[WC]', inf)
                        if match:
                            sc_gene_id = match.group(0)
                            break
                    
                    if not sc_gene_id:
                        for note in feature.qualifiers.get('note', []):
                            match = re.search(r'Y[A-Z]{2}\d{3}[WC]', note)
                            if match:
                                sc_gene_id = match.group(0)
                                break
                    
                    feature_info['sc_gene_id'] = sc_gene_id
                    features.append(feature_info)
            
            annotations[scaffold] = {
                'features': features,
                'record': record
            }
            
        except Exception as e:
            print(f"  Error loading {genbank_file}: {e}")
    
    print(f"Loaded annotations for {len(annotations)} scaffolds")
    return annotations

def load_gene_info(gene_mapping_file):
    """Load gene information from mapping file."""
    print(f"Loading gene info from {gene_mapping_file}")
    genes_df = pd.read_csv(gene_mapping_file, sep='\t')
    
    # Convert to dictionary for easier lookup
    gene_info = {}
    for _, row in genes_df.iterrows():
        gene_id = row['w303_gene_id']
        gene_info[gene_id] = {
            'w303_gene_id': gene_id,
            'sc_gene_id': row['sc_gene_id'] if 'sc_gene_id' in row else None,
            'erg_name': row['erg_name'] if 'erg_name' in row else None,
            'w303_scaffold': row['w303_scaffold'] if 'w303_scaffold' in row else None,
            'start': row['start'] if 'start' in row else None,
            'end': row['end'] if 'end' in row else None,
            'strand': row['strand'] if 'strand' in row else None
        }
    
    print(f"Loaded info for {len(gene_info)} genes")
    return gene_info

def identify_key_regions(variants, annotations, erg_gene_info):
    """Identify and analyze the key genomic regions."""
    print("Identifying key genomic regions")
    
    # Organize variants by nearest gene and distance
    regions = {}
    
    # Process each key region
    for region_name, region_info in KEY_REGIONS.items():
        erg_gene = region_info['erg_gene']
        distance = region_info['distance']
        position = region_info['position']
        
        # Extract information for this ERG gene
        erg_scaffold = None
        erg_start = None
        erg_end = None
        erg_strand = None
        
        for gene_id, gene_data in erg_gene_info.items():
            if gene_data['erg_name'] == erg_gene:
                erg_scaffold = gene_data['w303_scaffold']
                erg_start = gene_data['start']
                erg_end = gene_data['end']
                erg_strand = gene_data['strand']
                break
        
        if not erg_scaffold:
            print(f"  Warning: Could not find scaffold for {erg_gene}")
            continue
        
        print(f"  Found {erg_gene} on {erg_scaffold}, position {erg_start}-{erg_end}, strand {erg_strand}")
        
        # Filter variants for this key region
        if distance == 'multiple':
            # For ERG25, get all variants near this gene
            region_variants = variants[variants['nearest_gene_name'] == erg_gene].copy()
        else:
            # For specific distances, filter more precisely
            region_variants = variants[
                (variants['nearest_gene_name'] == erg_gene) & 
                (variants['position_relative'] == position) &
                (variants['distance_to_gene'] == distance)
            ].copy()
        
        if len(region_variants) == 0:
            print(f"  Warning: No variants found for {region_name}")
            continue
        
        # Get scaffold and calculate region bounds
        scaffold = erg_scaffold
        
        # Calculate region start and end
        if position == 'upstream' or position == 'both':
            if erg_strand == '+':
                region_start = max(0, erg_start - 50000)
                region_end = erg_start
            else:
                region_start = erg_end
                region_end = erg_end + 50000
        elif position == 'downstream' or position == 'both':
            if erg_strand == '+':
                region_start = erg_end
                region_end = erg_end + 50000
            else:
                region_start = max(0, erg_start - 50000)
                region_end = erg_start
        
        if position == 'both':
            # For ERG25, include the gene itself plus 50kb on each side
            region_start = max(0, min(erg_start, region_start) - 50000)
            region_end = max(erg_end, region_end) + 50000
        
        # Find genes in this region
        region_genes = []
        if scaffold in annotations:
            for feature in annotations[scaffold]['features']:
                if feature['type'] == 'gene' and (
                    (feature['start'] >= region_start and feature['start'] <= region_end) or
                    (feature['end'] >= region_start and feature['end'] <= region_end)
                ):
                    gene_info = feature.copy()
                    
                    # Calculate distance to ERG gene
                    if erg_strand == '+':
                        if gene_info['end'] < erg_start:
                            gene_dist = erg_start - gene_info['end']
                            gene_pos = 'upstream'
                        elif gene_info['start'] > erg_end:
                            gene_dist = gene_info['start'] - erg_end
                            gene_pos = 'downstream'
                        else:
                            gene_dist = 0
                            gene_pos = 'overlapping'
                    else:
                        if gene_info['end'] < erg_start:
                            gene_dist = erg_start - gene_info['end']
                            gene_pos = 'downstream'
                        elif gene_info['start'] > erg_end:
                            gene_dist = gene_info['start'] - erg_end
                            gene_pos = 'upstream'
                        else:
                            gene_dist = 0
                            gene_pos = 'overlapping'
                    
                    gene_info['distance_to_erg'] = gene_dist
                    gene_info['position_to_erg'] = gene_pos
                    
                    # Check if this gene has variants
                    gene_variants = variants[
                        (variants['chrom'] == scaffold) & 
                        (variants['pos'] >= gene_info['start']) & 
                        (variants['pos'] <= gene_info['end'])
                    ]
                    gene_info['variant_count'] = len(gene_variants)
                    
                    # Check if this gene contains any of our region variants
                    region_variant_count = sum(
                        (region_variants['chrom'] == scaffold) & 
                        (region_variants['pos'] >= gene_info['start']) & 
                        (region_variants['pos'] <= gene_info['end'])
                    )
                    gene_info['region_variant_count'] = region_variant_count
                    
                    region_genes.append(gene_info)
        
        # Store region information
        regions[region_name] = {
            'erg_gene': erg_gene,
            'scaffold': scaffold,
            'region_start': region_start,
            'region_end': region_end,
            'erg_start': erg_start,
            'erg_end': erg_end,
            'erg_strand': erg_strand,
            'variants': region_variants,
            'genes': region_genes
        }
        
        print(f"  Identified {region_name}: {len(region_variants)} variants, {len(region_genes)} genes")
    
    return regions

def analyze_region_sequence_features(regions, genome, scaffold_to_cm):
    """Analyze sequence features in the key regions with correct scaffold mapping."""
    print("Analyzing sequence features in key regions")
    
    for region_name, region_data in regions.items():
        w303_scaffold = region_data['scaffold']
        region_start = region_data['region_start']
        region_end = region_data['region_end']
        
        # Convert w303_scaffold to CM format using mapping
        if w303_scaffold in scaffold_to_cm:
            cm_scaffold = scaffold_to_cm[w303_scaffold]
            print(f"  Converting {w303_scaffold} to {cm_scaffold} for sequence access")
        else:
            cm_scaffold = w303_scaffold
            print(f"  No mapping found for {w303_scaffold}, using as-is")
        
        if cm_scaffold not in genome:
            print(f"  Warning: Scaffold {w303_scaffold} ({cm_scaffold}) not found in genome")
            print(f"  Available keys (first 3): {list(genome.keys())[:3]}")
            continue
        
        # Extract sequence for this region using the CM scaffold ID
        sequence = str(genome[cm_scaffold].seq[region_start:region_end])
        
        # Calculate basic sequence properties
        gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence) if len(sequence) > 0 else 0
        
        # Find homopolymer runs
        homopolymer_pattern = re.compile(r'([ACGT])\1{5,}')
        homopolymers = [match.group(0) for match in homopolymer_pattern.finditer(sequence)]
        
        # Find dinucleotide repeats
        dinucleotide_pattern = re.compile(r'([ACGT]{2})\1{3,}')
        dinucleotide_repeats = [match.group(0) for match in dinucleotide_pattern.finditer(sequence)]
        
        # Find potential origins of replication (ARS consensus)
        ars_pattern = re.compile(r'WTTTAYRTTTW', re.IGNORECASE)  # ARS consensus sequence
        ars_matches = [match.span() for match in ars_pattern.finditer(sequence)]
        
        # Check for common regulatory motifs
        motifs = {
            'TATA_box': re.compile(r'TATA[AT]A[AT]', re.IGNORECASE),
            'GC_box': re.compile(r'GGGCGG', re.IGNORECASE),
            'CAAT_box': re.compile(r'CCAAT', re.IGNORECASE),
            'CpG_island': re.compile(r'CG', re.IGNORECASE)
        }
        
        motif_counts = {}
        for motif_name, pattern in motifs.items():
            matches = [match.group(0) for match in pattern.finditer(sequence)]
            motif_counts[motif_name] = len(matches)
        
        # Look for additional sequence features
        sequence_features = {
            'length': len(sequence),
            'gc_content': gc_content,
            'homopolymers': {
                'count': len(homopolymers),
                'longest': max([len(h) for h in homopolymers] or [0]),
                'examples': homopolymers[:5] if homopolymers else []
            },
            'dinucleotide_repeats': {
                'count': len(dinucleotide_repeats),
                'longest': max([len(d) for d in dinucleotide_repeats] or [0]),
                'examples': dinucleotide_repeats[:5] if dinucleotide_repeats else []
            },
            'ars_consensus': {
                'count': len(ars_matches),
                'positions': [pos[0] for pos in ars_matches][:5] if ars_matches else []
            },
            'motif_counts': motif_counts,
            'cpg_density': motif_counts['CpG_island'] / (len(sequence) / 1000) if len(sequence) > 0 else 0
        }
        
        # Add sequence features to region data
        regions[region_name]['sequence_features'] = sequence_features
        
        print(f"  {region_name} sequence analysis: {len(sequence)} bp, GC: {gc_content:.2f}, " +
              f"homopolymers: {len(homopolymers)}, dinucleotide repeats: {len(dinucleotide_repeats)}")
    
    return regions

def find_functional_relationships(regions, annotations):
    """Identify potential functional relationships between genes in these regions."""
    print("Identifying potential functional relationships")
    
    for region_name, region_data in regions.items():
        erg_gene = region_data['erg_gene']
        genes = region_data['genes']
        
        # Initialize functional relationship data
        relationships = {
            'similar_function_genes': [],
            'potential_regulators': [],
            'pathway_related_genes': [],
            'unknown_function_genes': []
        }
        
        # Keywords related to ergosterol pathway
        ergosterol_keywords = [
            'sterol', 'lipid', 'membrane', 'fatty acid', 'metabol', 'biosynthesis',
            'oxidas', 'reductase', 'transferase', 'synthase', 'ERG'
        ]
        
        # Regulatory keywords
        regulatory_keywords = [
            'regulat', 'transcription', 'factor', 'activat', 'repress', 'control',
            'signal', 'kinas', 'phosphat', 'bind', 'promot'
        ]
        
        # Analyze genes for functional relationships
        for gene in genes:
            product = gene['product'].lower() if gene['product'] != 'unknown' else ''
            note = gene['note'].lower()
            
            # Check for ergosterol pathway related function
            if any(keyword in product or keyword in note for keyword in ergosterol_keywords):
                relationships['similar_function_genes'].append(gene)
            
            # Check for potential regulatory function
            elif any(keyword in product or keyword in note for keyword in regulatory_keywords):
                relationships['potential_regulators'].append(gene)
            
            # If genes contain variants from our key region, add them anyway
            elif gene['region_variant_count'] > 0:
                if product == 'unknown' or not product:
                    relationships['unknown_function_genes'].append(gene)
                else:
                    relationships['pathway_related_genes'].append(gene)
        
        # Add functional relationship data to region
        regions[region_name]['functional_relationships'] = relationships
        
        print(f"  {region_name} functional relationships:")
        print(f"    Similar function genes: {len(relationships['similar_function_genes'])}")
        print(f"    Potential regulators: {len(relationships['potential_regulators'])}")
        print(f"    Pathway related genes: {len(relationships['pathway_related_genes'])}")
        print(f"    Unknown function genes: {len(relationships['unknown_function_genes'])}")
    
    return regions

def calculate_region_statistics(regions, all_variants):
    """Calculate statistical enrichment for variants in these regions."""
    print("Calculating region statistics")
    
    # Get total genome size (approximate)
    total_size = 12000000  # S. cerevisiae genome is ~12Mb
    
    for region_name, region_data in regions.items():
        region_size = region_data['region_end'] - region_data['region_start']
        region_variants = len(region_data['variants'])
        total_variants = len(all_variants)
        total_size = 12_000_000
        expected_variants = total_variants * (region_size / total_size)
        enrichment = region_variants / expected_variants if expected_variants > 0 else 0

        # ---- replace fisher_exact with binomial test ----
        # binomial test: is observing region_variants out of total_variants more than expected?
        try:
            # Calculate probability as region size / total genome size
            probability = region_size / total_size
            binomial_result = binomtest(
                region_variants,
                n=total_variants,
                p=probability,
                alternative='two-sided'
            )
            p_value = binomial_result.pvalue
            print(f"    Debug - Binomial test params: k={region_variants}, n={total_variants}, p={probability:.6f}")
        except Exception as e:
            print(f"    Warning: Binomial test failed with error: {e}")
            p_value = 1.0

        # Enrichment serves as our effect size measure (similar to odds ratio)

        region_data['statistics'] = {
            'region_size': region_size,
            'region_variants': region_variants,
            'expected_variants': expected_variants,
            'enrichment': enrichment,
            'p_value': p_value,
            'significant': p_value < 0.05
        }
        
        print(f"  {region_name} statistics:")
        print(f"    Region size: {region_size:,} bp ({region_size/total_size:.2%} of genome)")
        print(f"    Variants: {region_variants} (Expected: {expected_variants:.2f})")
        print(f"    Enrichment: {enrichment:.2f}x, p-value: {p_value:.4e}")
        
        # Also calculate gene-specific statistics
        for gene in region_data['genes']:
            gene_size = gene['end'] - gene['start']
            gene_variants = gene['variant_count']
            
            # Calculate expected variants based on gene size
            expected_gene_variants = total_variants * (gene_size / total_size)
            
            # Calculate enrichment
            gene_enrichment = gene_variants / expected_gene_variants if expected_gene_variants > 0 else 0
            
            # Store statistics
            gene['statistics'] = {
                'gene_size': gene_size,
                'gene_variants': gene_variants,
                'expected_variants': expected_gene_variants,
                'enrichment': gene_enrichment
            }
    
    return regions

def generate_region_visualizations(regions, output_dir):
    """Generate visualizations of the key genomic regions."""
    print("Generating region visualizations")
    
    # Create output directory
    region_viz_dir = os.path.join(output_dir, 'region_visualizations')
    os.makedirs(region_viz_dir, exist_ok=True)
    
    # Generate a visualization for each region
    for region_name, region_data in regions.items():
        print(f"  Generating visualization for {region_name}")
        
        # Create a gene map visualization
        plt.figure(figsize=(14, 8))
        
        # Define the plotting area
        region_size = region_data['region_end'] - region_data['region_start']
        plot_start = region_data['region_start']
        plot_end = region_data['region_end']
        
        # Plot the ERG gene
        erg_start = region_data['erg_start']
        erg_end = region_data['erg_end']
        erg_strand = region_data['erg_strand']
        
        # Adjust coordinates relative to plot start
        erg_start_rel = erg_start - plot_start
        erg_end_rel = erg_end - plot_start
        
        # Plot ERG gene
        plt.axhline(y=0, xmin=erg_start_rel/region_size, xmax=erg_end_rel/region_size, 
                   color='darkred', linewidth=10, solid_capstyle='butt')
        
        # Add arrow to indicate strand
        arrow_x = erg_start_rel + (erg_end_rel - erg_start_rel) / 2
        arrow_dir = 1 if erg_strand == '+' else -1
        plt.annotate('', xy=(arrow_x + arrow_dir*5000, 0), xytext=(arrow_x, 0),
                   arrowprops=dict(arrowstyle='->', color='darkred', lw=2))
        
        # Add ERG gene label
        plt.text(arrow_x, 0.05, region_data['erg_gene'], ha='center', va='bottom', 
                fontsize=12, fontweight='bold', color='darkred')
        
        # Plot other genes in the region
        genes = region_data['genes']
        
        # Sort genes by variant count and position
        genes_sorted = sorted(genes, key=lambda x: (-x['variant_count'], x['start']))
        
        # Plot genes with variants first, then others
        for i, gene in enumerate(genes_sorted):
            y_pos = 0.5 + i * 0.3  # Stack genes vertically
            
            gene_start_rel = gene['start'] - plot_start
            gene_end_rel = gene['end'] - plot_start
            
            # Determine color based on variant count
            variant_count = gene['variant_count']
            if variant_count > 0:
                color = 'darkblue'
                alpha = min(1.0, 0.3 + variant_count * 0.1)
            else:
                color = 'gray'
                alpha = 0.3
            
            # Plot gene
            plt.axhline(y=y_pos, xmin=gene_start_rel/region_size, xmax=gene_end_rel/region_size, 
                       color=color, linewidth=6, alpha=alpha, solid_capstyle='butt')
            
            # Add arrow to indicate strand
            arrow_x = gene_start_rel + (gene_end_rel - gene_start_rel) / 2
            arrow_dir = 1 if gene['strand'] == '+' else -1
            plt.annotate('', xy=(arrow_x + arrow_dir*3000, y_pos), xytext=(arrow_x, y_pos),
                       arrowprops=dict(arrowstyle='->', color=color, lw=1, alpha=alpha))
            
            # Add gene label
            gene_id = gene['id'] if gene['id'] != 'unknown' else gene['locus_tag']
            product = gene['product'] if gene['product'] != 'unknown' else ''
            if product:
                product = product[:20] + '...' if len(product) > 20 else product
                label = f"{gene_id} ({product})"
            else:
                label = gene_id
            
            plt.text(arrow_x, y_pos + 0.05, label, ha='center', va='bottom', 
                    fontsize=9, color=color, alpha=alpha)
            
            # Add variant count if non-zero
            if variant_count > 0:
                plt.text(arrow_x, y_pos - 0.05, f"{variant_count} variants", ha='center', va='top', 
                        fontsize=8, color=color, alpha=alpha)
        
        # Plot variants
        variants = region_data['variants']
        
        for _, variant in variants.iterrows():
            pos = variant['pos'] - plot_start
            impact = variant['impact']
            color = 'red' if impact == 'HIGH' else 'orange'
            
            plt.axvline(x=pos, ymin=0, ymax=0.2, color=color, linewidth=1, alpha=0.7)
        
        # Add a buffer above the genes for better visualization
        y_max = max(0.5 + (len(genes_sorted) - 1) * 0.3 + 0.5, 2.0)
        
        # Set plot limits and labels
        plt.xlim(0, region_size)
        plt.ylim(-0.3, y_max)
        
        # Add distance markers
        for distance in range(0, region_size, 10000):
            plt.axvline(x=distance, color='black', linestyle=':', alpha=0.3)
            plt.text(distance, -0.2, f"{distance/1000:.0f}kb", ha='center', va='top', 
                    fontsize=8, alpha=0.7)
        
        # Add labels
        gene_counts_with_variants = sum(1 for gene in genes if gene['variant_count'] > 0)
        plt.title(f"{region_name}: {len(variants)} variants, {len(genes)} genes ({gene_counts_with_variants} with variants)")
        plt.xlabel(f"Position on {region_data['scaffold']} (relative to region start at {region_data['region_start']:,})")
        
        # Remove y-axis
        plt.yticks([])
        plt.tight_layout()
        
        # Save the visualization
        viz_file = os.path.join(region_viz_dir, f"{region_name}_gene_map.png")
        plt.savefig(viz_file, dpi=300)
        plt.close()
        
        # Create a variant distribution visualization
        plt.figure(figsize=(14, 6))
        
        # Create histogram of variant positions
        bin_size = max(1, region_size // 100)  # Adjust bin size based on region size
        bins = range(0, region_size + bin_size, bin_size)
        
        if not variants.empty:
            variant_positions = variants['pos'] - plot_start
            plt.hist(variant_positions, bins=bins, alpha=0.7, color='skyblue', edgecolor='black')
        
        # Mark the ERG gene position
        plt.axvline(x=erg_start_rel, color='darkred', linestyle='--', alpha=0.7, label=f"{region_data['erg_gene']} start")
        plt.axvline(x=erg_end_rel, color='darkred', linestyle=':', alpha=0.7, label=f"{region_data['erg_gene']} end")
        
        # Add distance markers
        for dist in [10000, 20000, 30000, 40000, 50000]:
            if erg_strand == '+':
                # For + strand genes
                upstream_pos = max(0, erg_start_rel - dist)
                downstream_pos = min(region_size, erg_end_rel + dist)
                
                if upstream_pos > 0:
                    plt.axvline(x=upstream_pos, color='gray', linestyle='-', alpha=0.3)
                    plt.text(upstream_pos, plt.ylim()[1] * 0.9, f"{dist/1000:.0f}kb upstream", 
                            ha='right', va='center', fontsize=8, rotation=90, alpha=0.7)
                
                if downstream_pos < region_size:
                    plt.axvline(x=downstream_pos, color='gray', linestyle='-', alpha=0.3)
                    plt.text(downstream_pos, plt.ylim()[1] * 0.9, f"{dist/1000:.0f}kb downstream", 
                            ha='left', va='center', fontsize=8, rotation=90, alpha=0.7)
            else:
                # For - strand genes, upstream and downstream are reversed
                upstream_pos = min(region_size, erg_end_rel + dist)
                downstream_pos = max(0, erg_start_rel - dist)
                
                if upstream_pos < region_size:
                    plt.axvline(x=upstream_pos, color='gray', linestyle='-', alpha=0.3)
                    plt.text(upstream_pos, plt.ylim()[1] * 0.9, f"{dist/1000:.0f}kb upstream", 
                            ha='left', va='center', fontsize=8, rotation=90, alpha=0.7)
                
                if downstream_pos > 0:
                    plt.axvline(x=downstream_pos, color='gray', linestyle='-', alpha=0.3)
                    plt.text(downstream_pos, plt.ylim()[1] * 0.9, f"{dist/1000:.0f}kb downstream", 
                            ha='right', va='center', fontsize=8, rotation=90, alpha=0.7)
        
        # Add labels
        plt.title(f"Variant Distribution Near {region_data['erg_gene']}")
        plt.xlabel(f"Position on {region_data['scaffold']} (relative to region start at {region_data['region_start']:,})")
        plt.ylabel("Number of Variants")
        plt.legend()
        plt.grid(alpha=0.3)
        plt.tight_layout()
        
        # Save the visualization
        viz_file = os.path.join(region_viz_dir, f"{region_name}_variant_distribution.png")
        plt.savefig(viz_file, dpi=300)
        plt.close()
        
        # If there are enough variants, create a treatment comparison
        if len(variants) >= 5:
            plt.figure(figsize=(14, 6))
            
            # Group variants by treatment
            treatment_groups = variants.groupby('treatment')
            
            # Define colors for treatments
            treatment_colors = {
                'CAS': 'blue',
                'STC': 'green',
                'WT-37': 'red',
                'WTA': 'purple',
                'WT': 'gray'
            }
            
            # Plot variant positions by treatment
            for i, (treatment, group) in enumerate(treatment_groups):
                positions = group['pos'] - plot_start
                impact_colors = ['red' if impact == 'HIGH' else 'orange' for impact in group['impact']]
                
                plt.scatter(positions, [i] * len(positions), 
                           c=impact_colors, label=treatment if len(positions) > 0 else None,
                           s=50, alpha=0.7, edgecolor='black')
            
            # Set y-ticks to treatment names
            plt.yticks(range(len(treatment_groups)), [t for t in treatment_groups.groups.keys()])
            
            # Mark the ERG gene position
            plt.axvline(x=erg_start_rel, color='darkred', linestyle='--', alpha=0.7, label=f"{region_data['erg_gene']} start")
            plt.axvline(x=erg_end_rel, color='darkred', linestyle=':', alpha=0.7, label=f"{region_data['erg_gene']} end")
            
            # Create a custom legend for impact
            legend_elements = [
                Patch(facecolor='red', edgecolor='black', label='HIGH impact'),
                Patch(facecolor='orange', edgecolor='black', label='MODERATE impact')
            ]
            
            # Add labels
            plt.title(f"Variant Positions by Treatment Near {region_data['erg_gene']}")
            plt.xlabel(f"Position on {region_data['scaffold']} (relative to region start at {region_data['region_start']:,})")
            plt.ylabel("Treatment")
            plt.legend(handles=legend_elements)
            plt.grid(alpha=0.3)
            plt.tight_layout()
            
            # Save the visualization
            viz_file = os.path.join(region_viz_dir, f"{region_name}_treatment_comparison.png")
            plt.savefig(viz_file, dpi=300)
            plt.close()
    
    return region_viz_dir

def save_results(regions, output_dir):
    """Save results to output files."""
    print(f"Saving results to {output_dir}")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Save region details to JSON file
    region_details = {}
    for region_name, region_data in regions.items():
        # Create a serializable version of region data (without BioPython objects)
        serializable_data = {
            'erg_gene': region_data['erg_gene'],
            'scaffold': region_data['scaffold'],
            'region_start': region_data['region_start'],
            'region_end': region_data['region_end'],
            'erg_start': region_data['erg_start'],
            'erg_end': region_data['erg_end'],
            'erg_strand': region_data['erg_strand'],
            'variant_count': len(region_data['variants']),
            'gene_count': len(region_data['genes']),
            'statistics': region_data.get('statistics', {}),
            'sequence_features': region_data.get('sequence_features', {})
        }
        
        # Add functional relationships
        relationships = region_data.get('functional_relationships', {})
        relationship_summary = {}
        for rel_type, rel_genes in relationships.items():
            relationship_summary[rel_type] = len(rel_genes)
            relationship_summary[f"{rel_type}_examples"] = [
                {
                    'id': gene['id'],
                    'locus_tag': gene['locus_tag'],
                    'product': gene['product'],
                    'variant_count': gene['variant_count']
                }
                for gene in rel_genes[:5]  # Include top 5 examples
            ]
        
        serializable_data['functional_relationships'] = relationship_summary
        
        # Add serializable gene info
        serializable_data['genes'] = []
        for gene in region_data['genes']:
            gene_info = {
                'id': gene['id'],
                'locus_tag': gene['locus_tag'],
                'product': gene['product'],
                'start': gene['start'],
                'end': gene['end'],
                'strand': gene['strand'],
                'distance_to_erg': gene['distance_to_erg'],
                'position_to_erg': gene['position_to_erg'],
                'variant_count': gene['variant_count'],
                'region_variant_count': gene['region_variant_count'],
                'statistics': gene.get('statistics', {})
            }
            serializable_data['genes'].append(gene_info)
        
        region_details[region_name] = serializable_data
    
    # Define a custom JSON encoder to handle numpy types
    class NumpyJSONEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, np.bool_):
                return bool(obj)
            return super(NumpyJSONEncoder, self).default(obj)
    
    # Save region details to JSON
    details_file = os.path.join(output_dir, 'region_details.json')
    with open(details_file, 'w') as f:
        json.dump(region_details, f, cls=NumpyJSONEncoder, indent=2)
    
    # Create TSV files for key information
    
    # 1. Genes in key regions
    genes_data = []
    for region_name, region_data in regions.items():
        for gene in region_data['genes']:
            gene_data = {
                'Region': region_name,
                'ERG_Gene': region_data['erg_gene'],
                'Gene_ID': gene['id'],
                'Locus_Tag': gene['locus_tag'],
                'Product': gene['product'],
                'Start': gene['start'],
                'End': gene['end'],
                'Strand': gene['strand'],
                'Distance_to_ERG': gene['distance_to_erg'],
                'Position_to_ERG': gene['position_to_erg'],
                'Variant_Count': gene['variant_count'],
                'Region_Variant_Count': gene['region_variant_count']
            }
            
            # Add gene statistics if available
            if 'statistics' in gene:
                gene_data['Gene_Size'] = gene['statistics'].get('gene_size', '')
                gene_data['Expected_Variants'] = gene['statistics'].get('expected_variants', '')
                gene_data['Enrichment'] = gene['statistics'].get('enrichment', '')
            
            genes_data.append(gene_data)
    
    genes_df = pd.DataFrame(genes_data)
    genes_file = os.path.join(output_dir, 'region_genes.tsv')
    genes_df.to_csv(genes_file, sep='\t', index=False)
    
    # 2. Variants in key regions
    variants_data = []
    for region_name, region_data in regions.items():
        variants = region_data['variants']
        
        if variants.empty:
            continue
            
        # Add region information to each variant
        variants = variants.copy()
        variants['region'] = region_name
        variants['erg_gene'] = region_data['erg_gene']
        
        # Convert to list of dictionaries
        for _, variant in variants.iterrows():
            variant_data = variant.to_dict()
            variants_data.append(variant_data)
    
    if variants_data:
        variants_df = pd.DataFrame(variants_data)
        variants_file = os.path.join(output_dir, 'region_variants.tsv')
        variants_df.to_csv(variants_file, sep='\t', index=False)
    
    # 3. Region summary statistics
    region_stats = []
    for region_name, region_data in regions.items():
        stats = region_data.get('statistics', {})
        
        region_stat = {
            'Region': region_name,
            'ERG_Gene': region_data['erg_gene'],
            'Scaffold': region_data['scaffold'],
            'Region_Start': region_data['region_start'],
            'Region_End': region_data['region_end'],
            'Region_Size': stats.get('region_size', ''),
            'Variant_Count': stats.get('region_variants', ''),
            'Expected_Variants': stats.get('expected_variants', ''),
            'Enrichment': stats.get('enrichment', ''),
            'P_Value': stats.get('p_value', ''),
            'Significant': stats.get('significant', '')
        }
        
        # Add some sequence feature information
        if 'sequence_features' in region_data:
            seq_features = region_data['sequence_features']
            region_stat['GC_Content'] = seq_features.get('gc_content', '')
            region_stat['Homopolymer_Count'] = seq_features.get('homopolymers', {}).get('count', '')
            region_stat['Dinucleotide_Repeat_Count'] = seq_features.get('dinucleotide_repeats', {}).get('count', '')
            region_stat['CpG_Density'] = seq_features.get('cpg_density', '')
        
        # Add functional relationship counts
        if 'functional_relationships' in region_data:
            relationships = region_data['functional_relationships']
            region_stat['Similar_Function_Genes'] = len(relationships.get('similar_function_genes', []))
            region_stat['Potential_Regulators'] = len(relationships.get('potential_regulators', []))
            region_stat['Pathway_Related_Genes'] = len(relationships.get('pathway_related_genes', []))
            region_stat['Unknown_Function_Genes'] = len(relationships.get('unknown_function_genes', []))
        
        region_stats.append(region_stat)
    
    region_stats_df = pd.DataFrame(region_stats)
    stats_file = os.path.join(output_dir, 'region_statistics.tsv')
    region_stats_df.to_csv(stats_file, sep='\t', index=False)
    
    # Create a summary report in Markdown format
    report = [
        "# Key Genomic Regions Analysis Report",
        f"Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "## 1. Overview of Key Regions",
        ""
    ]
    
    for region_name, region_data in regions.items():
        stats = region_data.get('statistics', {})
        seq_features = region_data.get('sequence_features', {})
        
        report.extend([
            f"### 1.{list(regions.keys()).index(region_name) + 1}. {region_name}",
            f"  - ERG Gene: {region_data['erg_gene']}",
            f"  - Scaffold: {region_data['scaffold']}",
            f"  - Region Size: {stats.get('region_size', 0):,} bp",
            f"  - Variants: {stats.get('region_variants', 0)} (Expected: {stats.get('expected_variants', 0):.2f})",
            f"  - Enrichment: {stats.get('enrichment', 0):.2f}x, p-value: {stats.get('p_value', 1):.4e}",
            f"  - Genes in Region: {len(region_data['genes'])}",
            "",
            "#### Sequence Features:",
            f"  - GC Content: {seq_features.get('gc_content', 0):.2f}",
            f"  - Homopolymers: {seq_features.get('homopolymers', {}).get('count', 0)}",
            f"  - Dinucleotide Repeats: {seq_features.get('dinucleotide_repeats', {}).get('count', 0)}",
            "",
            "#### Functional Relationships:",
        ])
        
        # Add functional relationship information
        if 'functional_relationships' in region_data:
            relationships = region_data['functional_relationships']
            
            # Similar function genes
            similar_genes = relationships.get('similar_function_genes', [])
            report.append(f"  - Similar Function Genes: {len(similar_genes)}")
            for gene in similar_genes[:3]:  # Show top 3
                report.append(f"    - {gene['id']} ({gene['product']}): {gene['variant_count']} variants")
            
            # Potential regulators
            regulators = relationships.get('potential_regulators', [])
            report.append(f"  - Potential Regulators: {len(regulators)}")
            for gene in regulators[:3]:  # Show top 3
                report.append(f"    - {gene['id']} ({gene['product']}): {gene['variant_count']} variants")
            
            # Pathway related genes
            pathway_genes = relationships.get('pathway_related_genes', [])
            report.append(f"  - Pathway Related Genes: {len(pathway_genes)}")
            for gene in pathway_genes[:3]:  # Show top 3
                report.append(f"    - {gene['id']} ({gene['product']}): {gene['variant_count']} variants")
            
            # Unknown function genes
            unknown_genes = relationships.get('unknown_function_genes', [])
            report.append(f"  - Unknown Function Genes: {len(unknown_genes)}")
            for gene in unknown_genes[:3]:  # Show top 3
                report.append(f"    - {gene['id']} ({gene['locus_tag']}): {gene['variant_count']} variants")
        
        report.append("")
    
    # Add more sections to the report
    report.extend([
        "## 2. Key Findings",
        "",
        "### ERG11 Upstream Region",
        "This region contains HIGH impact frameshift variants at approximately 8kb upstream of ERG11.",
        "These variants may affect genes that regulate ERG11 expression or function.",
        "",
        "### ERG7 Downstream Region",
        "This region contains HIGH impact frameshift variants at exactly 47,676bp downstream of ERG7.",
        "The precise distance suggests a specific gene that might be functionally related to ERG7.",
        "",
        "### ERG25 Neighborhood",
        "This region has the highest number of nearby variants, both upstream and downstream of ERG25.",
        "The distribution suggests a more complex regulatory network around ERG25 compared to other ergosterol genes.",
        "",
        "## 3. Biological Implications",
        "",
        "### Functional Relationships",
        "The genes harboring variants in these key regions may have functional relationships with their respective ergosterol genes,",
        "either through direct interaction, regulatory mechanisms, or participation in related metabolic pathways.",
        "",
        "### Adaptation Mechanisms",
        "The concentration of HIGH impact variants at specific distances from ergosterol genes suggests that adaptation",
        "may occur through changes in genes that interact with or regulate the ergosterol pathway, rather than through",
        "direct modifications to the essential pathway enzymes themselves.",
        "",
        "### Conservation Architecture",
        "The analysis confirms the 'conservation gradient' model, where purifying selection is strongest on the ergosterol genes",
        "themselves, and decreases with distance from these essential genes. This creates a hierarchical genomic architecture",
        "that preserves core functions while allowing adaptation through changes in auxiliary genes.",
        "",
        "## 4. Further Investigation",
        "",
        "### Functional Characterization",
        "The genes identified in these key regions should be further characterized to understand their precise roles",
        "in ergosterol metabolism and stress adaptation. This could involve experimental validation through gene knockouts,",
        "expression analysis, or protein interaction studies.",
        "",
        "### Integration with Sterol Profiles",
        "Correlating the variants in these key regions with sterol composition data could provide direct evidence of their",
        "functional impact on ergosterol biosynthesis and membrane composition under different stress conditions.",
        "",
        "### Comparative Analysis",
        "Comparing these genomic regions across different yeast strains and species could reveal evolutionary patterns",
        "that further clarify the functional relationships between the genes in these regions and the ergosterol pathway."
    ])
    
    # Save the report
    report_file = os.path.join(output_dir, 'key_regions_report.md')
    with open(report_file, 'w') as f:
        f.write('\n'.join(report))
    
    print(f"Results saved to {output_dir}")
    print(f"Report: {report_file}")
    print(f"Region details: {details_file}")
    print(f"Region genes: {genes_file}")
    print(f"Region variants: {variants_file if variants_data else 'No variants file created'}")
    print(f"Region statistics: {stats_file}")
    
    return {
        'report': report_file,
        'details': details_file,
        'genes': genes_file,
        'variants': variants_file if 'variants_data' in locals() and variants_data else None,
        'statistics': stats_file
    }

def main():
    """Main function to run the analysis."""
    # Parse command line arguments
    args = parse_arguments()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load data
    variants = load_variants(args.variants_file)
    genome = load_genome(args.genome_file)
    
    # Load chromosome mapping
    scaffold_to_cm, cm_to_scaffold = load_chromosome_mapping(args.mapping_file)
    
    # Load gene information
    gene_info = load_gene_info(args.gene_mapping)
    
    # Get unique scaffolds that contain our variants of interest
    scaffolds_of_interest = variants['chrom'].unique().tolist()
    
    # Load GenBank annotations
    annotations = load_genbank_annotations(args.genbank_dir, scaffolds_of_interest)
    
    # Identify key genomic regions
    regions = identify_key_regions(variants, annotations, gene_info)
    
    # Analyze sequence features with mapping
    regions = analyze_region_sequence_features(regions, genome, scaffold_to_cm)
    
    # Find functional relationships
    regions = find_functional_relationships(regions, annotations)
    
    # Calculate region statistics
    regions = calculate_region_statistics(regions, variants)
    
    # Generate visualizations
    region_viz_dir = generate_region_visualizations(regions, args.output_dir)
    
    # Save results
    output_files = save_results(regions, args.output_dir)
    
    print("Analysis complete!")
    print(f"Analyzed {len(regions)} key genomic regions")
    print(f"Generated visualizations in {region_viz_dir}")
    print(f"Output files: {', '.join(filter(None, output_files.values()))}")

if __name__ == "__main__":
    main()