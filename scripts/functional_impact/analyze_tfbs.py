#!/usr/bin/env python3
"""
analyze_tfbs.py - Modified version with scaffold name mapping

This script analyzes how upstream variants might affect transcription factor binding sites 
in the ergosterol pathway genes, with proper scaffold name conversion.
"""

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
import re
import json
from collections import defaultdict, Counter
import subprocess
import tempfile

# Known transcription factors involved in ergosterol pathway regulation
ERGOSTEROL_TFS = {
    'UPC2': 'TCGTATA',     # Sterol regulatory element binding protein
    'ECM22': 'TCGTATA',    # Similar binding site to UPC2
    'HAP1': 'CGGNNNTAA',   # Heme-responsive TF, active under aerobic conditions
    'ROX1': 'YSYATTGTT',   # Repressor under aerobic conditions
    'MOT3': 'AAGGKA',      # Repressor of hypoxic genes
    'TUP1': None,          # General repressor, no specific motif
    'RAP1': 'ACACCCATACATT', # Pleiotropic regulator
    'OAF1': 'CGGN{3}TNAN{9}CCG', # Involved in fatty acid metabolism
    'PDR1': 'TCCGCGGA',    # Pleiotropic drug resistance
    'PDR3': 'TCCGCGGA',    # Pleiotropic drug resistance
    'SUT1': 'GATA',        # Sterol uptake
    'SUT2': 'GATA',        # Sterol uptake
    'TBP': 'TATA[AT]A',    # TATA-binding protein
    'GCN4': 'RTGACTCAY',   # General control of amino acid synthesis
    'INO2': 'CATGTGAAAT',  # Phospholipid synthesis, interacts with ergosterol pathway
    'INO4': 'CATGTGAAAT'   # Phospholipid synthesis, interacts with ergosterol pathway
}

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Analyze transcription factor binding sites affected by variants')
    parser.add_argument('--variant_file', required=True, help='Path to the variants with TSS distance TSV file')
    parser.add_argument('--genome_file', required=True, help='Path to the reference genome FASTA file')
    parser.add_argument('--mapping_file', required=True, help='Path to the chromosome mapping TSV file')
    parser.add_argument('--output_dir', required=True, help='Directory to save results')
    parser.add_argument('--context_size', type=int, default=30, 
                      help='Size of sequence context to extract around each variant (default: 30bp)')
    parser.add_argument('--prediction_method', choices=['jaspar', 'motif_match', 'both'], default='both',
                      help='Method to use for TFBS prediction (default: both)')
    return parser.parse_args()

def load_variants(variant_file):
    """Load variant data from TSV file."""
    print(f"Loading variants from {variant_file}")
    variants = pd.read_csv(variant_file, sep='\t')
    variants.columns = [col.lower() for col in variants.columns]
    
    # Group variants by unique position to avoid redundant analysis
    grouped_variants = variants.groupby(['scaffold', 'position', 'ref', 'alt', 'gene_id', 'erg_name']).first().reset_index()
    
    print(f"Loaded {len(variants)} variants, grouped into {len(grouped_variants)} unique variant positions")
    return grouped_variants

def load_genome(genome_file):
    """Load reference genome from FASTA file."""
    print(f"Loading reference genome from {genome_file}")
    genome = {}
    for record in SeqIO.parse(genome_file, "fasta"):
        # Extract the first part of the ID (before any spaces)
        record_id = record.id.split()[0]
        genome[record_id] = str(record.seq)
    
    print(f"Loaded {len(genome)} scaffolds from reference genome")
    return genome

def load_chromosome_mapping(mapping_file):
    """Load chromosome mapping from TSV file."""
    print(f"Loading chromosome mapping from {mapping_file}")
    mapping = pd.read_csv(mapping_file, sep='\t')
    
    # Create bidirectional mappings
    w303_to_cm = dict(zip(mapping['w303_scaffold'], mapping['chromosome_id']))
    cm_to_w303 = dict(zip(mapping['chromosome_id'], mapping['w303_scaffold']))
    
    print(f"Loaded mapping for {len(w303_to_cm)} scaffolds")
    return w303_to_cm, cm_to_w303

def extract_sequence_context(variant, genome, chromosome_mapping, context_size=30):
    """
    Extract sequence context around a variant, handling different scaffold naming conventions.
    
    Args:
        variant: A pandas Series containing variant information
        genome: Dictionary mapping scaffold IDs to sequences
        chromosome_mapping: Dictionary mapping w303_scaffold_X to CM00XXXX.1
        context_size: Number of base pairs on each side of the variant
    
    Returns:
        tuple: (scaffold, position, reference sequence, variant sequence, variant info)
    """
    scaffold = variant['scaffold']
    position = variant['position']
    ref = variant['ref']
    alt = variant['alt']
    gene_id = variant['gene_id']
    erg_name = variant['erg_name']
    
    # Convert scaffold name if needed
    chromosome_id = chromosome_mapping.get(scaffold)
    if not chromosome_id:
        print(f"Warning: No mapping found for scaffold {scaffold}")
        return None
        
    if chromosome_id not in genome:
        print(f"Warning: Chromosome ID {chromosome_id} not found in reference genome")
        print(f"Available IDs: {list(genome.keys())[:5]}...")
        return None
    
    # Extract sequence context
    start = max(0, position - context_size - 1)  # -1 for 0-based indexing
    end = min(len(genome[chromosome_id]), position + context_size)
    
    seq_before = genome[chromosome_id][start:position-1]
    seq_after = genome[chromosome_id][position-1+len(ref):end]
    
    # Create reference and variant sequences
    ref_seq = seq_before + ref + seq_after
    var_seq = seq_before + alt + seq_after
    
    variant_info = {
        'scaffold': scaffold,
        'chromosome_id': chromosome_id,
        'position': position,
        'ref': ref,
        'alt': alt,
        'gene_id': gene_id,
        'erg_name': erg_name,
        'context_start': start + 1,  # Convert to 1-based
        'context_end': end,
        'tss_distance': variant['tss_distance'] if 'tss_distance' in variant else None
    }
    
    return (scaffold, position, ref_seq, var_seq, variant_info)

def predict_tfbs_motif_match(ref_seq, var_seq, tf_motifs=ERGOSTEROL_TFS):
    """
    Predict TFBS using direct motif matching.
    
    Args:
        ref_seq: Reference sequence
        var_seq: Variant sequence
        tf_motifs: Dictionary mapping TF names to motif patterns
        
    Returns:
        dict: Results for each TF with matches in ref and var sequences
    """
    results = {}
    
    for tf, motif in tf_motifs.items():
        if not motif:
            continue
            
        # Search for motif in reference sequence
        ref_matches = []
        if isinstance(motif, str):
            for match in re.finditer(motif, ref_seq, re.IGNORECASE):
                ref_matches.append((match.start(), match.end(), match.group()))
        
        # Search for motif in variant sequence
        var_matches = []
        if isinstance(motif, str):
            for match in re.finditer(motif, var_seq, re.IGNORECASE):
                var_matches.append((match.start(), match.end(), match.group()))
        
        # Only include TFs with at least one match
        if ref_matches or var_matches:
            results[tf] = {
                'motif': motif,
                'ref_matches': ref_matches,
                'var_matches': var_matches,
                'ref_count': len(ref_matches),
                'var_count': len(var_matches),
                'diff': len(var_matches) - len(ref_matches)
            }
    
    return results

def format_jaspar_output(jaspar_result):
    """Format JASPAR output into a more usable structure."""
    formatted = {}
    
    lines = jaspar_result.strip().split('\n')
    current_tf = None
    
    for line in lines:
        if line.startswith('#'):
            continue
        
        if line.startswith('>'):
            current_tf = line[1:].strip()
            formatted[current_tf] = {'sites': []}
        elif current_tf:
            parts = line.strip().split()
            if len(parts) >= 4:  # Position, strand, score, site
                position = int(parts[0])
                strand = parts[1]
                score = float(parts[2])
                site = parts[3]
                
                formatted[current_tf]['sites'].append({
                    'position': position,
                    'strand': strand,
                    'score': score,
                    'site': site
                })
    
    # Calculate summary stats for each TF
    for tf in formatted:
        if formatted[tf]['sites']:
            scores = [site['score'] for site in formatted[tf]['sites']]
            formatted[tf]['max_score'] = max(scores)
            formatted[tf]['site_count'] = len(scores)
    
    return formatted

def predict_tfbs_jaspar(ref_seq, var_seq, temp_dir):
    """
    Predict TFBS using JASPAR database (if available).
    
    Args:
        ref_seq: Reference sequence
        var_seq: Variant sequence
        temp_dir: Temporary directory for files
        
    Returns:
        tuple: (ref_predictions, var_predictions)
    """
    # Check if JASPAR tools are available
    try:
        # Write sequences to temporary files
        ref_file = os.path.join(temp_dir, "ref_seq.fa")
        var_file = os.path.join(temp_dir, "var_seq.fa")
        
        with open(ref_file, "w") as f:
            f.write(">ref\n" + ref_seq + "\n")
        
        with open(var_file, "w") as f:
            f.write(">var\n" + var_seq + "\n")
        
        # Since JASPAR tool access varies, we'll use a Python-based approach
        # This is a simplified example - in practice, you might use a proper TFBS prediction tool
        
        # For now, return empty predictions as a placeholder
        # In a real implementation, you would call the appropriate JASPAR tool here
        return {}, {}
        
    except Exception as e:
        print(f"Warning: JASPAR prediction failed: {e}")
        return {}, {}

def compare_tfbs_predictions(ref_predictions, var_predictions):
    """
    Compare TFBS predictions between reference and variant sequences.
    
    Args:
        ref_predictions: Predictions for reference sequence
        var_predictions: Predictions for variant sequence
        
    Returns:
        dict: Comparison results
    """
    comparison = {
        'gained_sites': [],
        'lost_sites': [],
        'modified_sites': []
    }
    
    # Identify TFs in both predictions
    all_tfs = set(ref_predictions.keys()) | set(var_predictions.keys())
    
    for tf in all_tfs:
        ref_count = ref_predictions.get(tf, {}).get('ref_count', 0)
        var_count = var_predictions.get(tf, {}).get('var_count', 0)
        
        if ref_count == 0 and var_count > 0:
            comparison['gained_sites'].append({
                'tf': tf,
                'var_count': var_count,
                'var_matches': var_predictions[tf]['var_matches']
            })
        elif ref_count > 0 and var_count == 0:
            comparison['lost_sites'].append({
                'tf': tf,
                'ref_count': ref_count,
                'ref_matches': ref_predictions[tf]['ref_matches']
            })
        elif ref_count != var_count:
            comparison['modified_sites'].append({
                'tf': tf,
                'ref_count': ref_count,
                'var_count': var_count,
                'diff': var_count - ref_count,
                'ref_matches': ref_predictions[tf]['ref_matches'] if tf in ref_predictions else [],
                'var_matches': var_predictions[tf]['var_matches'] if tf in var_predictions else []
            })
    
    comparison['summary'] = {
        'gained_count': len(comparison['gained_sites']),
        'lost_count': len(comparison['lost_sites']),
        'modified_count': len(comparison['modified_sites']),
        'total_changes': len(comparison['gained_sites']) + len(comparison['lost_sites']) + len(comparison['modified_sites'])
    }
    
    return comparison

def analyze_tfbs_changes(variant_contexts, output_dir, prediction_method='both'):
    """
    Analyze TFBS changes for each variant context.
    
    Args:
        variant_contexts: List of (scaffold, position, ref_seq, var_seq, variant_info) tuples
        output_dir: Directory to save results
        prediction_method: Method to use for TFBS prediction
        
    Returns:
        dict: Results for all variants
    """
    print(f"Analyzing TFBS changes for {len(variant_contexts)} variant contexts")
    
    # Create temporary directory for files
    with tempfile.TemporaryDirectory() as temp_dir:
        
        all_results = []
        
        for context in variant_contexts:
            if not context:
                continue
                
            scaffold, position, ref_seq, var_seq, variant_info = context
            
            result = {
                'variant_info': variant_info,
                'ref_seq': ref_seq,
                'var_seq': var_seq,
                'tfbs_changes': {}
            }
            
            # Predict TFBS using direct motif matching
            if prediction_method in ['motif_match', 'both']:
                ref_predictions = predict_tfbs_motif_match(ref_seq, var_seq)
                var_predictions = predict_tfbs_motif_match(var_seq, var_seq)
                motif_comparison = compare_tfbs_predictions(ref_predictions, var_predictions)
                result['tfbs_changes']['motif_match'] = {
                    'ref_predictions': ref_predictions,
                    'var_predictions': var_predictions,
                    'comparison': motif_comparison
                }
            
            # Predict TFBS using JASPAR (if available)
            if prediction_method in ['jaspar', 'both']:
                ref_jaspar, var_jaspar = predict_tfbs_jaspar(ref_seq, var_seq, temp_dir)
                result['tfbs_changes']['jaspar'] = {
                    'ref_predictions': ref_jaspar,
                    'var_predictions': var_jaspar
                }
            
            all_results.append(result)
        
        print(f"Completed TFBS analysis for {len(all_results)} variants")
        return all_results

def summarize_tfbs_results(results):
    """
    Create a summary of TFBS analysis results.
    
    Args:
        results: List of TFBS analysis results
        
    Returns:
        dict: Summary statistics
    """
    summary = {
        'total_variants': len(results),
        'variants_by_gene': defaultdict(int),
        'variants_with_tfbs_changes': 0,
        'variants_by_tss_region': defaultdict(int),
        'tf_changes': defaultdict(int),
        'gained_sites': defaultdict(int),
        'lost_sites': defaultdict(int)
    }
    
    for result in results:
        gene = result['variant_info']['erg_name']
        summary['variants_by_gene'][gene] += 1
        
        # Categorize by TSS region
        tss_dist = result['variant_info'].get('tss_distance')
        if tss_dist is not None:
            if tss_dist < -500:
                region = "Far Upstream (>500bp)"
            elif tss_dist < -250:
                region = "Upstream Activating Sequence (250-500bp)"
            elif tss_dist < -50:
                region = "TATA Box Region (50-250bp)"
            elif tss_dist < 0:
                region = "Core Promoter (0-50bp)"
            else:
                region = "Downstream"
            
            summary['variants_by_tss_region'][region] += 1
        
        # Count variants with TFBS changes
        changes = result['tfbs_changes'].get('motif_match', {}).get('comparison', {})
        if changes.get('summary', {}).get('total_changes', 0) > 0:
            summary['variants_with_tfbs_changes'] += 1
            
            # Count changes by TF
            for gained in changes.get('gained_sites', []):
                summary['gained_sites'][gained['tf']] += 1
                summary['tf_changes'][gained['tf']] += 1
            
            for lost in changes.get('lost_sites', []):
                summary['lost_sites'][lost['tf']] += 1
                summary['tf_changes'][lost['tf']] += 1
            
            for modified in changes.get('modified_sites', []):
                summary['tf_changes'][modified['tf']] += 1
    
    # Convert defaultdicts to regular dicts for JSON serialization
    for key in summary:
        if isinstance(summary[key], defaultdict):
            summary[key] = dict(summary[key])
    
    return summary

def generate_tss_tfbs_plot(results, output_dir):
    """
    Generate a plot showing TFBS changes by TSS distance.
    
    Args:
        results: List of TFBS analysis results
        output_dir: Directory to save plots
        
    Returns:
        str: Path to saved plot
    """
    print(f"Generating TSS-TFBS plot")
    
    # Extract data for plotting
    distances = []
    change_types = []
    genes = []
    
    for result in results:
        tss_dist = result['variant_info'].get('tss_distance')
        gene = result['variant_info']['erg_name']
        
        if tss_dist is not None:
            changes = result['tfbs_changes'].get('motif_match', {}).get('comparison', {})
            change_summary = changes.get('summary', {})
            
            if change_summary.get('total_changes', 0) > 0:
                change_type = 'TFBS Changes'
            else:
                change_type = 'No TFBS Changes'
            
            distances.append(tss_dist)
            change_types.append(change_type)
            genes.append(gene)
    
    # Create DataFrame for plotting
    df = pd.DataFrame({
        'TSS Distance': distances,
        'Change Type': change_types,
        'Gene': genes
    })
    
    # Create plot
    if df.empty:
        print("No data for TSS-TFBS plot")
        return None
        
    plt.figure(figsize=(12, 6))
    
    # Plot scatter points
    sns.scatterplot(data=df, x='TSS Distance', y='Gene', hue='Change Type', style='Change Type', s=100)
    
    # Add vertical lines for typical promoter elements
    plt.axvspan(-150, -50, alpha=0.2, color='lightgreen', label='TATA Box Region')
    plt.axvspan(-50, 0, alpha=0.2, color='lightblue', label='Core Promoter')
    
    plt.title('TFBS Changes by Distance from Transcription Start Site')
    plt.xlabel('Distance from TSS (bp) - Negative values are upstream')
    plt.ylabel('Gene')
    plt.legend(title='TFBS Effect')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    
    # Save plot
    plot_path = os.path.join(output_dir, 'tfbs_tss_distribution.png')
    plt.savefig(plot_path, dpi=300)
    plt.close()
    
    return plot_path

def generate_tf_heatmap(results, output_dir):
    """
    Generate a heatmap of TF binding changes by gene.
    
    Args:
        results: List of TFBS analysis results
        output_dir: Directory to save plots
        
    Returns:
        str: Path to saved plot
    """
    print(f"Generating TF heatmap")
    
    # Collect data for the heatmap
    tf_changes = defaultdict(lambda: defaultdict(int))
    genes = set()
    tfs = set()
    
    for result in results:
        gene = result['variant_info']['erg_name']
        genes.add(gene)
        
        changes = result['tfbs_changes'].get('motif_match', {}).get('comparison', {})
        
        for gained in changes.get('gained_sites', []):
            tf = gained['tf']
            tfs.add(tf)
            tf_changes[gene][tf] += 1
        
        for lost in changes.get('lost_sites', []):
            tf = lost['tf']
            tfs.add(tf)
            tf_changes[gene][tf] -= 1
    
    # Convert to DataFrame
    heatmap_data = []
    for gene in genes:
        for tf in tfs:
            heatmap_data.append({
                'Gene': gene,
                'Transcription Factor': tf,
                'Net Change': tf_changes[gene][tf]
            })
    
    df = pd.DataFrame(heatmap_data)
    if df.empty:
        print("No TF binding changes found for heatmap")
        return None
    
    # Pivot for heatmap
    pivot_df = df.pivot(index='Gene', columns='Transcription Factor', values='Net Change').fillna(0)
    
    # Create heatmap
    plt.figure(figsize=(12, 8))
    
    cmap = sns.diverging_palette(240, 10, as_cmap=True)
    sns.heatmap(pivot_df, cmap=cmap, center=0, annot=True, fmt='.0f')
    
    plt.title('Net Changes in Transcription Factor Binding Sites by Gene')
    plt.tight_layout()
    
    # Save plot
    plot_path = os.path.join(output_dir, 'tf_gene_heatmap.png')
    plt.savefig(plot_path, dpi=300)
    plt.close()
    
    return plot_path

def generate_sequence_visualizations(results, output_dir):
    """
    Generate visualizations of sequence changes and their effects on TFBS.
    
    Args:
        results: List of TFBS analysis results
        output_dir: Directory to save visualizations
        
    Returns:
        list: Paths to saved visualizations
    """
    print(f"Generating sequence visualizations")
    
    # Create output directory for sequence visualizations
    seq_vis_dir = os.path.join(output_dir, 'sequence_visualizations')
    os.makedirs(seq_vis_dir, exist_ok=True)
    
    viz_paths = []
    
    for i, result in enumerate(results):
        variant_info = result['variant_info']
        scaffold = variant_info['scaffold']
        position = variant_info['position']
        gene = variant_info['erg_name']
        
        # Check if there are any TFBS changes
        changes = result['tfbs_changes'].get('motif_match', {}).get('comparison', {})
        if changes.get('summary', {}).get('total_changes', 0) == 0:
            continue
        
        # Create sequence visualization
        ref_seq = result['ref_seq']
        var_seq = result['var_seq']
        
        # Find the central index where the variant occurs
        central_index = len(ref_seq) // 2
        
        # Create figure with two sequences
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 6))
        
        # Draw reference sequence
        ax1.text(0.01, 0.5, "Reference:", fontsize=12, ha='left', va='center')
        for j, base in enumerate(ref_seq):
            color = 'black'
            weight = 'normal'
            bg_color = 'white'
            
            # Highlight the reference allele
            if abs(j - central_index) < len(variant_info['ref']):
                color = 'red'
                weight = 'bold'
                bg_color = 'lightyellow'
            
            # Check if this position is part of a TFBS
            for tf, pred in result['tfbs_changes'].get('motif_match', {}).get('ref_predictions', {}).items():
                for start, end, _ in pred.get('ref_matches', []):
                    if start <= j < end:
                        bg_color = 'lightblue'
            
            ax1.text(0.1 + 0.02 * j, 0.5, base, fontsize=10, ha='center', va='center', 
                   color=color, weight=weight, bbox=dict(facecolor=bg_color, alpha=0.5))
        
        # Draw variant sequence
        ax2.text(0.01, 0.5, "Variant:", fontsize=12, ha='left', va='center')
        for j, base in enumerate(var_seq):
            color = 'black'
            weight = 'normal'
            bg_color = 'white'
            
            # Highlight the variant allele
            if abs(j - central_index) < len(variant_info['alt']):
                color = 'red'
                weight = 'bold'
                bg_color = 'lightyellow'
            
            # Check if this position is part of a TFBS
            for tf, pred in result['tfbs_changes'].get('motif_match', {}).get('var_predictions', {}).items():
                for start, end, _ in pred.get('var_matches', []):
                    if start <= j < end:
                        bg_color = 'lightgreen'
            
            ax2.text(0.1 + 0.02 * j, 0.5, base, fontsize=10, ha='center', va='center', 
                   color=color, weight=weight, bbox=dict(facecolor=bg_color, alpha=0.5))
        
        # Add legend
        ax1.axis('off')
        ax2.axis('off')
        
        # Add title with variant information
        title = f"Variant: {scaffold}:{position} {variant_info['ref']}>{variant_info['alt']} in {gene}"
        if variant_info.get('tss_distance') is not None:
            title += f" (TSS Distance: {variant_info['tss_distance']}bp)"
        plt.suptitle(title, fontsize=14)
        
        # Add TFBS changes information
        changes_text = "TFBS Changes:\n"
        for gained in changes.get('gained_sites', []):
            changes_text += f"Gained: {gained['tf']} (+{gained['var_count']})\n"
        for lost in changes.get('lost_sites', []):
            changes_text += f"Lost: {lost['tf']} (-{lost['ref_count']})\n"
        for modified in changes.get('modified_sites', []):
            sign = "+" if modified['diff'] > 0 else ""
            changes_text += f"Modified: {modified['tf']} ({sign}{modified['diff']})\n"
        
        plt.figtext(0.01, 0.02, changes_text, fontsize=10, ha='left', va='bottom', 
                   bbox=dict(facecolor='lightgray', alpha=0.5))
        
        # Save figure
        viz_path = os.path.join(seq_vis_dir, f"{gene}_{scaffold}_{position}_tfbs_viz.png")
        plt.savefig(viz_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        viz_paths.append(viz_path)
    
    print(f"Generated {len(viz_paths)} sequence visualizations")
    return viz_paths

def save_results(results, summary, output_dir):
    """
    Save analysis results to output directory.
    
    Args:
        results: List of TFBS analysis results
        summary: Summary statistics
        output_dir: Directory to save results
        
    Returns:
        dict: Paths to saved files
    """
    print(f"Saving results to {output_dir}")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Save full results
    results_path = os.path.join(output_dir, 'tfbs_analysis_results.json')
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    # Save summary
    summary_path = os.path.join(output_dir, 'tfbs_analysis_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Create a TSV report of variants with TFBS changes
    report_data = []
    for result in results:
        variant_info = result['variant_info']
        changes = result['tfbs_changes'].get('motif_match', {}).get('comparison', {})
        
        if changes.get('summary', {}).get('total_changes', 0) > 0:
            gained_tfs = [item['tf'] for item in changes.get('gained_sites', [])]
            lost_tfs = [item['tf'] for item in changes.get('lost_sites', [])]
            modified_tfs = [item['tf'] for item in changes.get('modified_sites', [])]
            
            report_data.append({
                'Gene': variant_info['erg_name'],
                'Scaffold': variant_info['scaffold'],
                'Chromosome_ID': variant_info['chromosome_id'],
                'Position': variant_info['position'],
                'Ref': variant_info['ref'],
                'Alt': variant_info['alt'],
                'TSS_Distance': variant_info.get('tss_distance', 'Unknown'),
                'TF_Gained': ', '.join(gained_tfs) or 'None',
                'TF_Lost': ', '.join(lost_tfs) or 'None',
                'TF_Modified': ', '.join(modified_tfs) or 'None',
                'Total_Changes': changes.get('summary', {}).get('total_changes', 0)
            })
    
    if report_data:
        report_df = pd.DataFrame(report_data)
        report_path = os.path.join(output_dir, 'tfbs_changes_report.tsv')
        report_df.to_csv(report_path, sep='\t', index=False)
    else:
        report_path = None
        with open(os.path.join(output_dir, 'tfbs_changes_report.txt'), 'w') as f:
            f.write("No TFBS changes detected in any variants")
    
    # Generate a summary table by gene
    gene_summary_data = []
    for gene, count in summary.get('variants_by_gene', {}).items():
        gene_summary_data.append({
            'Gene': gene,
            'Total_Variants': count,
            'Variants_With_TFBS_Changes': len([r for r in results if r['variant_info']['erg_name'] == gene and 
                                          r['tfbs_changes'].get('motif_match', {}).get('comparison', {}).get('summary', {}).get('total_changes', 0) > 0]),
            'Most_Common_Gained_TF': max([(tf, summary.get('gained_sites', {}).get(tf, 0)) 
                                       for tf in summary.get('gained_sites', {}) 
                                       if any(r['variant_info']['erg_name'] == gene and 
                                             tf in [g['tf'] for g in r['tfbs_changes'].get('motif_match', {}).get('comparison', {}).get('gained_sites', [])])],
                                    key=lambda x: x[1], default=('None', 0))[0],
            'Most_Common_Lost_TF': max([(tf, summary.get('lost_sites', {}).get(tf, 0)) 
                                     for tf in summary.get('lost_sites', {}) 
                                     if any(r['variant_info']['erg_name'] == gene and 
                                           tf in [l['tf'] for l in r['tfbs_changes'].get('motif_match', {}).get('comparison', {}).get('lost_sites', [])])],
                                  key=lambda x: x[1], default=('None', 0))[0]
        })
    
    if gene_summary_data:
        gene_summary_df = pd.DataFrame(gene_summary_data)
        gene_summary_path = os.path.join(output_dir, 'tfbs_gene_summary.tsv')
        gene_summary_df.to_csv(gene_summary_path, sep='\t', index=False)
    else:
        gene_summary_path = None
        with open(os.path.join(output_dir, 'tfbs_gene_summary.txt'), 'w') as f:
            f.write("No gene-specific TFBS changes to summarize")
    
    return {
        'results': results_path,
        'summary': summary_path,
        'report': report_path,
        'gene_summary': gene_summary_path
    }

def main():
    """Main function to run the analysis."""
    # Parse command line arguments
    args = parse_arguments()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load data
    variants = load_variants(args.variant_file)
    genome = load_genome(args.genome_file)
    w303_to_cm, cm_to_w303 = load_chromosome_mapping(args.mapping_file)
    
    # Extract sequence contexts for each variant
    print("Extracting sequence contexts for variants")
    variant_contexts = []
    for _, variant in variants.iterrows():
        context = extract_sequence_context(variant, genome, w303_to_cm, args.context_size)
        variant_contexts.append(context)
    
    # Analyze TFBS changes
    results = analyze_tfbs_changes(variant_contexts, args.output_dir, args.prediction_method)
    
    # Generate summary statistics
    summary = summarize_tfbs_results(results)
    
    # Generate visualizations
    tss_plot = generate_tss_tfbs_plot(results, args.output_dir)
    tf_heatmap = generate_tf_heatmap(results, args.output_dir)
    sequence_viz = generate_sequence_visualizations(results, args.output_dir)
    
    # Save results
    file_paths = save_results(results, summary, args.output_dir)
    
    print("Analysis complete!")
    print(f"Analyzed {len(results)} variant contexts")
    print(f"Found {summary['variants_with_tfbs_changes']} variants with TFBS changes")
    print(f"Results saved to {args.output_dir}")

if __name__ == "__main__":
    main()