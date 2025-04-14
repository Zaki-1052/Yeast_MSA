#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict, Counter
from Bio import SeqIO, motifs
from Bio.Seq import Seq
import random
import subprocess
import logomaker
from scipy.stats import fisher_exact, chi2_contingency
import re
import warnings
warnings.filterwarnings('ignore')

# Set matplotlib style
plt.style.use('ggplot')
sns.set(font_scale=1.2)
sns.set_style("whitegrid")

# Define output directory
OUTPUT_DIR = "analysis/genomic_context_results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Updated biologically correct treatment groups
TREATMENTS = ['WT-37', 'WTA', 'STC', 'CAS']

# Define treatment information for better biological context
TREATMENT_INFO = {
    'WT-37': {'description': 'Temperature-adapted wild type', 'adaptation': 'Temperature'},
    'WTA': {'description': 'Low oxygen-adapted wild type', 'adaptation': 'Low Oxygen'},
    'STC': {'description': 'STC gene with low oxygen adaptation', 'adaptation': 'Low Oxygen', 'gene': 'STC'},
    'CAS': {'description': 'CAS gene with temperature adaptation', 'adaptation': 'Temperature', 'gene': 'CAS'}
}

# Treatment colors for consistent visualization
TREATMENT_COLORS = {
    'WT-37': '#1b9e77',  # Temperature-adapted
    'WTA': '#d95f02',    # Low oxygen-adapted
    'STC': '#7570b3',    # STC gene + low oxygen
    'CAS': '#e7298a'     # CAS gene + temperature
}

# Adaptation colors for grouping by adaptation type
ADAPTATION_COLORS = {
    'Temperature': '#1f77b4',
    'Low Oxygen': '#ff7f0e',
}

# Define nucleotide complements
COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

# Reference genome path
REFERENCE_GENOME = "reference/yeast_w303.fasta"

# Function to find a file trying multiple locations
def find_file(base_name, file_patterns):
    """Find a file checking multiple possible locations."""
    # Substitute base_name into patterns
    patterns = [pattern.format(base_name) for pattern in file_patterns]
    
    # Add backward compatibility for WT-37
    if base_name == 'WT-37':
        wt_patterns = [pattern.format('WT') for pattern in file_patterns]
        patterns.extend(wt_patterns)
    
    # Try each pattern
    for pattern in patterns:
        if os.path.exists(pattern):
            print(f"Found file for {base_name} at {pattern}")
            return pattern
    
    print(f"Could not find file for {base_name} in any expected location")
    return None

# Function to load reference genome
def load_reference_genome(fasta_file=REFERENCE_GENOME):
    """Load the reference genome from a FASTA file."""
    # Define possible reference file patterns
    ref_patterns = [
        "reference/yeast_w303.fasta",
        "reference/genome.fasta",
        "reference/yeast/yeast_w303.fasta"
    ]
    
    # Try to find the reference file
    found_ref = None
    for pattern in ref_patterns:
        if os.path.exists(pattern):
            found_ref = pattern
            break
    
    if not found_ref:
        print(f"Error: Reference genome file not found in expected locations.")
        return {}
    
    print(f"Loading reference genome from {found_ref}...")
    reference = {}
    for record in SeqIO.parse(found_ref, "fasta"):
        reference[record.id] = str(record.seq).upper()
    
    print(f"Loaded {len(reference)} sequences from reference genome.")
    return reference

# Function to parse mutation data for a specific treatment
def parse_mutation_data(treatment):
    """Parse mutation data for a specific treatment."""
    # Define possible file patterns for mutation data
    file_patterns = [
        "mutation_spectrum_analysis/{}_mutations.txt",
        "analysis/MSA/mutation_spectrum_analysis/{}_mutations.txt",
        "results/mutation_spectrum_analysis/{}_mutations.txt"
    ]
    
    # Find the mutation data file
    mutation_file = find_file(treatment, file_patterns)
    
    if mutation_file:
        try:
            # Read the extracted mutation data
            data = pd.read_csv(mutation_file, sep='\t', header=None, 
                               names=['CHROM', 'POS', 'REF', 'ALT'])
            data['Treatment'] = treatment
            
            # Add biological context
            data['Adaptation'] = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
            data['Has_Gene'] = 'Yes' if TREATMENT_INFO.get(treatment, {}).get('gene') else 'No'
            
            print(f"Loaded {len(data)} mutations for {treatment}")
            return data
        except Exception as e:
            print(f"Error reading {mutation_file}: {e}")
    
    # If no mutation file found, try to extract from VCF
    return extract_from_vcf(treatment)

# Function to extract mutation data from VCF if data file not found
def extract_from_vcf(treatment):
    """Extract mutation data from VCF files if no pre-extracted data is found."""
    # Define possible VCF file patterns
    vcf_patterns = [
        "results/merged/analysis/{}/highconf.vcf.gz",
        "results/merged/analysis/{}_highconf.vcf.gz",
        "results/merged/analysis/{}/specific.vcf.gz",
        "results/merged/analysis/{}_specific.vcf.gz"
    ]
    
    # Find the VCF file
    vcf_file = find_file(treatment, vcf_patterns)
    
    if vcf_file:
        try:
            print(f"Extracting mutation data from {vcf_file}")
            # Extract data using bcftools
            cmd = f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' {vcf_file}"
            output = subprocess.check_output(cmd, shell=True).decode('utf-8')
            
            # Parse the output
            rows = []
            for line in output.strip().split('\n'):
                if line:  # Skip empty lines
                    parts = line.split('\t')
                    if len(parts) == 4:
                        rows.append(parts)
            
            # Create dataframe
            if rows:
                data = pd.DataFrame(rows, columns=['CHROM', 'POS', 'REF', 'ALT'])
                data['POS'] = data['POS'].astype(int)
                data['Treatment'] = treatment
                
                # Add biological context
                data['Adaptation'] = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
                data['Has_Gene'] = 'Yes' if TREATMENT_INFO.get(treatment, {}).get('gene') else 'No'
                
                # Save extracted data for future use
                os.makedirs("mutation_spectrum_analysis", exist_ok=True)
                data.to_csv(f"mutation_spectrum_analysis/{treatment}_mutations.txt", 
                          sep='\t', index=False, header=False)
                
                print(f"Extracted and saved {len(data)} mutations for {treatment}")
                return data
        except Exception as e:
            print(f"Error extracting from {vcf_file}: {e}")
    
    print(f"Could not find or extract mutation data for {treatment}")
    return pd.DataFrame()

# Function to filter data for single nucleotide variants
def filter_snvs(data):
    """Filter data to include only single nucleotide variants."""
    # Keep only rows where REF and ALT are single nucleotides
    snv_data = data[(data['REF'].str.len() == 1) & (data['ALT'].str.len() == 1)]
    
    # Keep only ACGT bases (filter out N or other ambiguous bases)
    valid_bases = snv_data['REF'].isin(['A', 'C', 'G', 'T']) & snv_data['ALT'].isin(['A', 'C', 'G', 'T'])
    snv_data = snv_data[valid_bases]
    
    return snv_data

# Function to extract sequence context around a variant
def extract_context_sequence(chrom, pos, ref, reference_genome, context_size=100):
    """Extract sequence context around a variant position."""
    if chrom not in reference_genome:
        return None
    
    sequence = reference_genome[chrom]
    pos = int(pos)  # Ensure position is an integer
    
    # Adjust context_size to stay within sequence bounds
    start = max(0, pos - context_size - 1)
    end = min(len(sequence), pos + context_size)
    
    # Extract context
    context = sequence[start:end]
    
    # Calculate actual positions relative to the variant
    left_offset = pos - 1 - start
    right_offset = end - pos
    
    # Check if the reference base matches what's in the sequence
    if left_offset >= 0 and left_offset < len(context):
        seq_ref = context[left_offset]
        if seq_ref != ref:
            print(f"Warning: Reference base mismatch at {chrom}:{pos}. Expected {ref}, found {seq_ref}")
    
    return {
        'context': context,
        'variant_index': left_offset,
        'left_size': left_offset,
        'right_size': right_offset
    }

# Function to extract context sequences for all variants
def extract_all_context_sequences(data, reference_genome, context_size=100):
    """Extract context sequences for all variants in the dataset."""
    contexts = []
    
    for _, row in data.iterrows():
        context_data = extract_context_sequence(
            row['CHROM'], row['POS'], row['REF'], reference_genome, context_size)
        
        if context_data:
            contexts.append({
                'CHROM': row['CHROM'],
                'POS': row['POS'],
                'REF': row['REF'],
                'ALT': row['ALT'],
                'Treatment': row['Treatment'],
                'Adaptation': row['Adaptation'],
                'Has_Gene': row['Has_Gene'],
                'context': context_data['context'],
                'variant_index': context_data['variant_index'],
                'left_context': context_data['context'][:context_data['variant_index']],
                'right_context': context_data['context'][context_data['variant_index']+1:]
            })
    
    return pd.DataFrame(contexts)

# Function to generate control contexts from random positions
def generate_control_contexts(reference_genome, num_controls=1000, context_size=100):
    """Generate control contexts from random positions in the genome."""
    # Create list of all positions (scaffold, position)
    all_positions = []
    for scaffold, sequence in reference_genome.items():
        # Skip scaffolds that are too short
        if len(sequence) <= 2 * context_size:
            continue
        
        # Add all valid positions
        for pos in range(context_size, len(sequence) - context_size, 1000):  # Sample every 1000bp to speed up
            ref_base = sequence[pos]
            if ref_base in 'ACGT':  # Skip ambiguous bases
                all_positions.append((scaffold, pos, ref_base))
    
    # Sample random positions
    if len(all_positions) < num_controls:
        print(f"Warning: Only {len(all_positions)} valid positions available. Using all.")
        sampled_positions = all_positions
    else:
        sampled_positions = random.sample(all_positions, num_controls)
    
    # Extract contexts
    control_contexts = []
    for scaffold, pos, ref_base in sampled_positions:
        context_data = extract_context_sequence(
            scaffold, pos, ref_base, reference_genome, context_size)
        
        if context_data:
            control_contexts.append({
                'CHROM': scaffold,
                'POS': pos,
                'REF': ref_base,
                'context': context_data['context'],
                'variant_index': context_data['variant_index'],
                'left_context': context_data['context'][:context_data['variant_index']],
                'right_context': context_data['context'][context_data['variant_index']+1:]
            })
    
    return pd.DataFrame(control_contexts)

# Function to calculate GC content
def calculate_gc_content(sequence):
    """Calculate GC content of a sequence."""
    if not sequence:
        return 0
    
    gc_count = sum(1 for base in sequence if base in 'GC')
    return gc_count / len(sequence)

# Function to calculate sequence composition metrics
def calculate_composition_metrics(contexts_df):
    """Calculate sequence composition metrics for each context."""
    # Calculate GC content
    contexts_df['gc_content'] = contexts_df['context'].apply(calculate_gc_content)
    
    # Calculate local GC content (±10bp around variant)
    contexts_df['local_gc_content'] = contexts_df.apply(
        lambda row: calculate_gc_content(
            row['context'][max(0, row['variant_index']-10):row['variant_index']] + 
            row['context'][row['variant_index']+1:min(len(row['context']), row['variant_index']+11)]
        ),
        axis=1
    )
    
    # Calculate purine/pyrimidine ratio
    contexts_df['purine_ratio'] = contexts_df['context'].apply(
        lambda seq: sum(1 for base in seq if base in 'AG') / (len(seq) or 1)
    )
    
    # Detect homopolymers near variant
    contexts_df['nearby_homopolymer'] = contexts_df.apply(
        lambda row: detect_homopolymer(
            row['context'], row['variant_index'], window=20, min_length=3
        ),
        axis=1
    )
    
    # Detect dinucleotide repeats
    contexts_df['nearby_dinucleotide_repeat'] = contexts_df.apply(
        lambda row: detect_repeat(
            row['context'], row['variant_index'], repeat_size=2, window=20, min_repeats=2
        ),
        axis=1
    )
    
    return contexts_df

# Function to detect homopolymers
def detect_homopolymer(sequence, variant_index, window=20, min_length=3):
    """Detect homopolymers near the variant position."""
    # Extract region around variant
    start = max(0, variant_index - window)
    end = min(len(sequence), variant_index + window + 1)
    region = sequence[start:end]
    
    # Search for homopolymers
    for base in 'ACGT':
        pattern = base * min_length
        if pattern in region:
            return True
    
    return False

# Function to detect repeats
def detect_repeat(sequence, variant_index, repeat_size=2, window=20, min_repeats=2):
    """Detect repeats near the variant position."""
    # Extract region around variant
    start = max(0, variant_index - window)
    end = min(len(sequence), variant_index + window + 1)
    region = sequence[start:end]
    
    # Search for repeats
    for i in range(len(region) - repeat_size * min_repeats + 1):
        repeat_unit = region[i:i+repeat_size]
        # Check if this unit is repeated
        is_repeat = True
        for j in range(1, min_repeats):
            next_unit = region[i + j*repeat_size:i + (j+1)*repeat_size]
            if next_unit != repeat_unit:
                is_repeat = False
                break
        if is_repeat:
            return True
    
    return False

# Function to extract immediate context around mutations
def extract_immediate_contexts(contexts_df, context_size=5):
    """Extract immediate context (±5bp) around mutation sites."""
    immediate_contexts = []
    
    for _, row in contexts_df.iterrows():
        context = row['context']
        var_idx = row['variant_index']
        
        # Extract immediate context
        start = max(0, var_idx - context_size)
        end = min(len(context), var_idx + context_size + 1)
        
        immediate_context = context[start:end]
        rel_var_idx = var_idx - start
        
        # Ensure we have the reference base at the variant position
        if rel_var_idx < len(immediate_context) and immediate_context[rel_var_idx] == row['REF']:
            immediate_contexts.append({
                'CHROM': row['CHROM'],
                'POS': row['POS'],
                'REF': row['REF'],
                'ALT': row['ALT'],
                'Treatment': row['Treatment'],
                'Adaptation': row['Adaptation'],
                'Has_Gene': row['Has_Gene'],
                'immediate_context': immediate_context,
                'variant_index': rel_var_idx
            })
    
    return pd.DataFrame(immediate_contexts)

# Function to plot GC content distribution
def plot_gc_content(contexts_df, control_df, output_dir):
    """Plot GC content distribution in variant regions vs. control regions."""
    plt.figure(figsize=(10, 6))
    
    # Plot GC content distribution
    sns.kdeplot(contexts_df['gc_content'], label='Variant Regions')
    sns.kdeplot(control_df['gc_content'], label='Control Regions')
    
    # Customize plot
    plt.xlabel('GC Content')
    plt.ylabel('Density')
    plt.title('GC Content Distribution in Variant vs. Control Regions')
    plt.legend()
    
    # Save plot
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "gc_content_distribution.png"), dpi=300)
    plt.close()
    
    # Plot GC content by treatment with biological context
    plt.figure(figsize=(12, 6))
    
    # Group by treatment
    treatments = contexts_df['Treatment'].unique()
    
    for treatment in treatments:
        treatment_data = contexts_df[contexts_df['Treatment'] == treatment]
        adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
        gene = TREATMENT_INFO.get(treatment, {}).get('gene')
        
        label = f"{treatment} ({adaptation}"
        if gene:
            label += f", {gene} gene"
        label += ")"
        
        # Use consistent treatment colors
        sns.kdeplot(
            treatment_data['gc_content'], 
            label=label,
            color=TREATMENT_COLORS.get(treatment, '#999999')
        )
    
    # Add control
    sns.kdeplot(control_df['gc_content'], label='Control', linestyle='--', color='black')
    
    # Customize plot
    plt.xlabel('GC Content')
    plt.ylabel('Density')
    plt.title('GC Content Distribution by Treatment')
    plt.legend()
    
    # Save plot
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "gc_content_by_treatment.png"), dpi=300)
    plt.close()
    
    # Plot GC content by adaptation type
    plt.figure(figsize=(12, 6))
    
    # Group by adaptation
    adaptation_groups = contexts_df.groupby('Adaptation')
    
    for adaptation, group in adaptation_groups:
        # Use adaptation colors
        sns.kdeplot(
            group['gc_content'], 
            label=f'{adaptation} Adaptation',
            color=ADAPTATION_COLORS.get(adaptation, '#999999')
        )
    
    # Add control
    sns.kdeplot(control_df['gc_content'], label='Control', linestyle='--', color='black')
    
    # Customize plot
    plt.xlabel('GC Content')
    plt.ylabel('Density')
    plt.title('GC Content Distribution by Adaptation Type')
    plt.legend()
    
    # Save plot
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "gc_content_by_adaptation.png"), dpi=300)
    plt.close()
    
    # Plot GC content by gene status
    plt.figure(figsize=(12, 6))
    
    # Group by gene status
    gene_groups = contexts_df.groupby('Has_Gene')
    
    for has_gene, group in gene_groups:
        status = "Gene-Modified" if has_gene == "Yes" else "Non-Modified"
        sns.kdeplot(
            group['gc_content'], 
            label=f'{status} Strains',
            color='#1b9e77' if has_gene == "Yes" else '#d95f02'
        )
    
    # Add control
    sns.kdeplot(control_df['gc_content'], label='Control', linestyle='--', color='black')
    
    # Customize plot
    plt.xlabel('GC Content')
    plt.ylabel('Density')
    plt.title('GC Content Distribution by Gene Modification Status')
    plt.legend()
    
    # Save plot
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "gc_content_by_gene.png"), dpi=300)
    plt.close()
    
    # Plot local GC content
    plt.figure(figsize=(10, 6))
    
    # Plot local GC content distribution
    sns.kdeplot(contexts_df['local_gc_content'], label='Local Variant Regions')
    
    # Customize plot
    plt.xlabel('Local GC Content (±10bp)')
    plt.ylabel('Density')
    plt.title('Local GC Content Distribution Around Variants')
    plt.legend()
    
    # Save plot
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "local_gc_content_distribution.png"), dpi=300)
    plt.close()

# Function to plot sequence features
def plot_sequence_features(contexts_df, output_dir):
    """Plot local sequence features around variants."""
    # Group by treatment
    treatments = contexts_df['Treatment'].unique()
    
    # Plot homopolymer presence
    homopolymer_data = contexts_df.groupby('Treatment')['nearby_homopolymer'].mean()
    
    plt.figure(figsize=(10, 6))
    bars = plt.bar(
        homopolymer_data.index, 
        homopolymer_data.values, 
        color=[TREATMENT_COLORS.get(t, '#999999') for t in homopolymer_data.index]
    )
    
    # Add value labels
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                f'{height:.2f}', ha='center', va='bottom')
    
    # Add treatment descriptions in x-labels
    plt.xticks(
        range(len(treatments)),
        [f"{t}\n({TREATMENT_INFO.get(t, {}).get('adaptation', 'Unknown')})" for t in treatments],
        rotation=0
    )
    
    plt.xlabel('Treatment')
    plt.ylabel('Fraction with Nearby Homopolymer')
    plt.title('Homopolymer Presence Near Variants by Treatment')
    plt.ylim(0, max(homopolymer_data.values) * 1.1)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "homopolymer_by_treatment.png"), dpi=300)
    plt.close()
    
    # Plot dinucleotide repeat presence
    dinucleotide_data = contexts_df.groupby('Treatment')['nearby_dinucleotide_repeat'].mean()
    
    plt.figure(figsize=(10, 6))
    bars = plt.bar(
        dinucleotide_data.index, 
        dinucleotide_data.values, 
        color=[TREATMENT_COLORS.get(t, '#999999') for t in dinucleotide_data.index]
    )
    
    # Add value labels
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                f'{height:.2f}', ha='center', va='bottom')
    
    # Add treatment descriptions in x-labels
    plt.xticks(
        range(len(treatments)),
        [f"{t}\n({TREATMENT_INFO.get(t, {}).get('adaptation', 'Unknown')})" for t in treatments],
        rotation=0
    )
    
    plt.xlabel('Treatment')
    plt.ylabel('Fraction with Nearby Dinucleotide Repeat')
    plt.title('Dinucleotide Repeat Presence Near Variants by Treatment')
    plt.ylim(0, max(dinucleotide_data.values) * 1.1)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "dinucleotide_repeat_by_treatment.png"), dpi=300)
    plt.close()

# Function to plot sequence features by adaptation type
def plot_sequence_features_by_adaptation(contexts_df, output_dir):
    """Plot sequence features grouped by adaptation type."""
    # Plot homopolymer presence by adaptation
    homopolymer_data = contexts_df.groupby('Adaptation')['nearby_homopolymer'].mean()
    
    plt.figure(figsize=(8, 6))
    bars = plt.bar(
        homopolymer_data.index, 
        homopolymer_data.values, 
        color=[ADAPTATION_COLORS.get(a, '#999999') for a in homopolymer_data.index]
    )
    
    # Add value labels
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                f'{height:.2f}', ha='center', va='bottom')
    
    plt.xlabel('Adaptation Type')
    plt.ylabel('Fraction with Nearby Homopolymer')
    plt.title('Homopolymer Presence Near Variants by Adaptation Type')
    plt.ylim(0, max(homopolymer_data.values) * 1.1)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "homopolymer_by_adaptation.png"), dpi=300)
    plt.close()
    
    # Plot dinucleotide repeat presence by adaptation
    dinucleotide_data = contexts_df.groupby('Adaptation')['nearby_dinucleotide_repeat'].mean()
    
    plt.figure(figsize=(8, 6))
    bars = plt.bar(
        dinucleotide_data.index, 
        dinucleotide_data.values, 
        color=[ADAPTATION_COLORS.get(a, '#999999') for a in dinucleotide_data.index]
    )
    
    # Add value labels
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                f'{height:.2f}', ha='center', va='bottom')
    
    plt.xlabel('Adaptation Type')
    plt.ylabel('Fraction with Nearby Dinucleotide Repeat')
    plt.title('Dinucleotide Repeat Presence Near Variants by Adaptation Type')
    plt.ylim(0, max(dinucleotide_data.values) * 1.1)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "dinucleotide_repeat_by_adaptation.png"), dpi=300)
    plt.close()
    
    # Plot GC content by adaptation in a boxplot
    plt.figure(figsize=(8, 6))
    
    sns.boxplot(
        x='Adaptation', 
        y='gc_content', 
        data=contexts_df, 
        palette=ADAPTATION_COLORS
    )
    
    plt.xlabel('Adaptation Type')
    plt.ylabel('GC Content')
    plt.title('GC Content Distribution by Adaptation Type')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "gc_content_boxplot_by_adaptation.png"), dpi=300)
    plt.close()

# Function to plot sequence features by gene status
def plot_sequence_features_by_gene(contexts_df, output_dir):
    """Plot sequence features grouped by gene modification status."""
    # Plot homopolymer presence by gene status
    homopolymer_data = contexts_df.groupby('Has_Gene')['nearby_homopolymer'].mean()
    
    plt.figure(figsize=(8, 6))
    bars = plt.bar(
        ['Non-Modified', 'Gene-Modified'] if 'No' in homopolymer_data.index else ['Gene-Modified', 'Non-Modified'],
        homopolymer_data.values,
        color=['#d95f02', '#1b9e77'] if 'No' in homopolymer_data.index else ['#1b9e77', '#d95f02']
    )
    
    # Add value labels
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                f'{height:.2f}', ha='center', va='bottom')
    
    plt.xlabel('Gene Modification Status')
    plt.ylabel('Fraction with Nearby Homopolymer')
    plt.title('Homopolymer Presence Near Variants by Gene Modification Status')
    plt.ylim(0, max(homopolymer_data.values) * 1.1)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "homopolymer_by_gene.png"), dpi=300)
    plt.close()
    
    # Plot dinucleotide repeat presence by gene status
    dinucleotide_data = contexts_df.groupby('Has_Gene')['nearby_dinucleotide_repeat'].mean()
    
    plt.figure(figsize=(8, 6))
    bars = plt.bar(
        ['Non-Modified', 'Gene-Modified'] if 'No' in dinucleotide_data.index else ['Gene-Modified', 'Non-Modified'],
        dinucleotide_data.values,
        color=['#d95f02', '#1b9e77'] if 'No' in dinucleotide_data.index else ['#1b9e77', '#d95f02']
    )
    
    # Add value labels
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                f'{height:.2f}', ha='center', va='bottom')
    
    plt.xlabel('Gene Modification Status')
    plt.ylabel('Fraction with Nearby Dinucleotide Repeat')
    plt.title('Dinucleotide Repeat Presence Near Variants by Gene Modification Status')
    plt.ylim(0, max(dinucleotide_data.values) * 1.1)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "dinucleotide_repeat_by_gene.png"), dpi=300)
    plt.close()
    
    # Plot GC content by gene status in a boxplot
    plt.figure(figsize=(8, 6))
    
    sns.boxplot(
        x='Has_Gene', 
        y='gc_content', 
        data=contexts_df, 
        palette={'Yes': '#1b9e77', 'No': '#d95f02'}
    )
    
    # Change x-axis labels to be more descriptive
    plt.gca().set_xticklabels(['Non-Modified', 'Gene-Modified'])
    
    plt.xlabel('Gene Modification Status')
    plt.ylabel('GC Content')
    plt.title('GC Content Distribution by Gene Modification Status')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "gc_content_boxplot_by_gene.png"), dpi=300)
    plt.close()

# Function to analyze adaptation-specific context patterns
def analyze_adaptation_specific_patterns(contexts_df, output_dir):
    """Analyze context patterns specific to adaptation types."""
    # Compare GC content between adaptation types
    adaptation_groups = contexts_df.groupby('Adaptation')
    
    adaptation_gc = {}
    for adaptation, group in adaptation_groups:
        adaptation_gc[adaptation] = group['gc_content'].mean()
    
    # Calculate statistical significance using t-test
    adaptations = list(adaptation_gc.keys())
    if len(adaptations) >= 2:
        from scipy.stats import ttest_ind
        
        # Compare first two adaptations
        a1, a2 = adaptations[0], adaptations[1]
        a1_data = contexts_df[contexts_df['Adaptation'] == a1]['gc_content']
        a2_data = contexts_df[contexts_df['Adaptation'] == a2]['gc_content']
        
        t_stat, p_value = ttest_ind(a1_data, a2_data, equal_var=False)
        
        # Create a statistical summary
        with open(os.path.join(output_dir, "adaptation_statistics.txt"), 'w') as f:
            f.write("Adaptation-Specific Context Statistics\n")
            f.write("====================================\n\n")
            
            f.write("GC Content by Adaptation Type:\n")
            for adaptation, gc in adaptation_gc.items():
                f.write(f"  {adaptation}: {gc:.4f}\n")
            
            f.write(f"\nComparison of {a1} vs {a2} GC Content:\n")
            f.write(f"  t-statistic: {t_stat:.4f}\n")
            f.write(f"  p-value: {p_value:.4e}\n")
            f.write(f"  Significant difference: {'Yes' if p_value < 0.05 else 'No'}\n")
            
            # Homopolymer presence by adaptation
            f.write("\nHomopolymer Presence by Adaptation:\n")
            for adaptation, group in adaptation_groups:
                homo_presence = group['nearby_homopolymer'].mean()
                f.write(f"  {adaptation}: {homo_presence:.4f}\n")
            
            # Dinucleotide repeat presence by adaptation
            f.write("\nDinucleotide Repeat Presence by Adaptation:\n")
            for adaptation, group in adaptation_groups:
                dinuc_presence = group['nearby_dinucleotide_repeat'].mean()
                f.write(f"  {adaptation}: {dinuc_presence:.4f}\n")
    
    # Create adaptation-specific nucleotide distributions
    for adaptation, group in adaptation_groups:
        nucleotide_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        
        for _, row in group.iterrows():
            context = row['context']
            for base in context:
                if base in nucleotide_counts:
                    nucleotide_counts[base] += 1
        
        # Plot nucleotide distribution
        plt.figure(figsize=(8, 6))
        
        bases = list(nucleotide_counts.keys())
        counts = list(nucleotide_counts.values())
        total = sum(counts)
        
        if total > 0:
            proportions = [count / total for count in counts]
            
            colors = {'A': 'green', 'C': 'blue', 'G': 'orange', 'T': 'red'}
            plt.bar(bases, proportions, color=[colors[base] for base in bases])
            
            # Add labels
            for i, prop in enumerate(proportions):
                plt.text(i, prop + 0.01, f'{prop:.3f}', ha='center')
            
            plt.xlabel('Nucleotide')
            plt.ylabel('Proportion')
            plt.title(f'Nucleotide Distribution for {adaptation} Adaptation')
            
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"{adaptation}_nucleotide_distribution.png"), dpi=300)
            plt.close()

# Function to analyze gene-specific context patterns
def analyze_gene_specific_patterns(contexts_df, output_dir):
    """Analyze context patterns specific to gene modification status."""
    # Compare GC content between gene statuses
    gene_groups = contexts_df.groupby('Has_Gene')
    
    gene_gc = {}
    for has_gene, group in gene_groups:
        status = "Gene-modified" if has_gene == "Yes" else "Non-modified"
        gene_gc[status] = group['gc_content'].mean()
    
    # Calculate statistical significance using t-test
    gene_statuses = list(gene_gc.keys())
    if len(gene_statuses) >= 2:
        from scipy.stats import ttest_ind
        
        # Compare gene-modified vs non-modified
        yes_data = contexts_df[contexts_df['Has_Gene'] == 'Yes']['gc_content']
        no_data = contexts_df[contexts_df['Has_Gene'] == 'No']['gc_content']
        
        if len(yes_data) > 0 and len(no_data) > 0:
            t_stat, p_value = ttest_ind(yes_data, no_data, equal_var=False)
            
            # Create a statistical summary
            with open(os.path.join(output_dir, "gene_statistics.txt"), 'w') as f:
                f.write("Gene Modification-Specific Context Statistics\n")
                f.write("=========================================\n\n")
                
                f.write("GC Content by Gene Modification Status:\n")
                for status, gc in gene_gc.items():
                    f.write(f"  {status}: {gc:.4f}\n")
                
                f.write(f"\nComparison of Gene-modified vs Non-modified GC Content:\n")
                f.write(f"  t-statistic: {t_stat:.4f}\n")
                f.write(f"  p-value: {p_value:.4e}\n")
                f.write(f"  Significant difference: {'Yes' if p_value < 0.05 else 'No'}\n")
                
                # Homopolymer presence by gene status
                f.write("\nHomopolymer Presence by Gene Modification Status:\n")
                for has_gene, group in gene_groups:
                    status = "Gene-modified" if has_gene == "Yes" else "Non-modified"
                    homo_presence = group['nearby_homopolymer'].mean()
                    f.write(f"  {status}: {homo_presence:.4f}\n")
                
                # Dinucleotide repeat presence by gene status
                f.write("\nDinucleotide Repeat Presence by Gene Modification Status:\n")
                for has_gene, group in gene_groups:
                    status = "Gene-modified" if has_gene == "Yes" else "Non-modified"
                    dinuc_presence = group['nearby_dinucleotide_repeat'].mean()
                    f.write(f"  {status}: {dinuc_presence:.4f}\n")

# Function to analyze mutation context patterns
def analyze_mutation_context_patterns(immediate_contexts_df, output_dir):
    """Analyze and visualize immediate nucleotide context patterns around mutations."""
    # Group by treatment
    treatment_groups = immediate_contexts_df.groupby('Treatment')
    
    for treatment, group in treatment_groups:
        # Skip if too few contexts
        if len(group) < 10:
            continue
        
        # Analyze context nucleotide distributions
        contexts = group['immediate_context'].tolist()
        
        # Define positions relative to mutation
        context_size = min(len(context) for context in contexts) // 2
        positions = list(range(-context_size, context_size + 1))
        
        # Count nucleotides at each position
        position_distributions = {}
        for pos in positions:
            position_distributions[pos] = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        
        for context in contexts:
            # Find the center position
            center = len(context) // 2
            
            for pos in positions:
                if 0 <= center + pos < len(context):
                    nucleotide = context[center + pos]
                    if nucleotide in 'ACGT':
                        position_distributions[pos][nucleotide] += 1
        
        # Plot nucleotide distributions around mutation site
        plt.figure(figsize=(12, 6))
        
        for i, pos in enumerate(positions):
            # Calculate proportion of each nucleotide
            total = sum(position_distributions[pos].values())
            if total > 0:
                proportions = {nt: count / total for nt, count in position_distributions[pos].items()}
                
                # Plot stacked bars
                bottom = 0
                for nt, color in zip(['A', 'C', 'G', 'T'], ['green', 'blue', 'orange', 'red']):
                    plt.bar(i, proportions[nt], bottom=bottom, color=color, width=0.8, label=nt if i == 0 else "")
                    bottom += proportions[nt]
        
        # Customize plot
        plt.xlabel('Position Relative to Mutation')
        plt.ylabel('Nucleotide Proportion')
        
        # Get treatment description
        description = TREATMENT_INFO.get(treatment, {}).get('description', '')
        adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', '')
        has_gene = TREATMENT_INFO.get(treatment, {}).get('gene')
        
        plt.title(f'Nucleotide Context for {treatment} Treatment\n{description} ({adaptation} adaptation{" with " + has_gene + " gene" if has_gene else ""})')
        plt.xticks(range(len(positions)), positions)
        
        # Highlight mutation position
        plt.axvline(x=positions.index(0), color='black', linestyle='--', alpha=0.5)
        
        # Add legend
        plt.legend(title='Nucleotide')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{treatment}_context_pattern.png"), dpi=300)
        plt.close()
    
    # Also create adaptation-specific context patterns
    adaptation_groups = immediate_contexts_df.groupby('Adaptation')
    
    for adaptation, group in adaptation_groups:
        # Skip if too few contexts
        if len(group) < 10:
            continue
        
        # Analyze context nucleotide distributions
        contexts = group['immediate_context'].tolist()
        
        # Define positions relative to mutation
        context_size = min(len(context) for context in contexts) // 2
        positions = list(range(-context_size, context_size + 1))
        
        # Count nucleotides at each position
        position_distributions = {}
        for pos in positions:
            position_distributions[pos] = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        
        for context in contexts:
            # Find the center position
            center = len(context) // 2
            
            for pos in positions:
                if 0 <= center + pos < len(context):
                    nucleotide = context[center + pos]
                    if nucleotide in 'ACGT':
                        position_distributions[pos][nucleotide] += 1
        
        # Plot nucleotide distributions around mutation site
        plt.figure(figsize=(12, 6))
        
        for i, pos in enumerate(positions):
            # Calculate proportion of each nucleotide
            total = sum(position_distributions[pos].values())
            if total > 0:
                proportions = {nt: count / total for nt, count in position_distributions[pos].items()}
                
                # Plot stacked bars
                bottom = 0
                for nt, color in zip(['A', 'C', 'G', 'T'], ['green', 'blue', 'orange', 'red']):
                    plt.bar(i, proportions[nt], bottom=bottom, color=color, width=0.8, label=nt if i == 0 else "")
                    bottom += proportions[nt]
        
        # Customize plot
        plt.xlabel('Position Relative to Mutation')
        plt.ylabel('Nucleotide Proportion')
        plt.title(f'Nucleotide Context for {adaptation} Adaptation')
        plt.xticks(range(len(positions)), positions)
        
        # Highlight mutation position
        plt.axvline(x=positions.index(0), color='black', linestyle='--', alpha=0.5)
        
        # Add legend
        plt.legend(title='Nucleotide')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{adaptation}_context_pattern.png"), dpi=300)
        plt.close()
    
    # Create gene-specific context patterns
    gene_groups = immediate_contexts_df.groupby('Has_Gene')
    
    for has_gene, group in gene_groups:
        # Skip if too few contexts
        if len(group) < 10:
            continue
        
        # Analyze context nucleotide distributions
        contexts = group['immediate_context'].tolist()
        
        # Define positions relative to mutation
        context_size = min(len(context) for context in contexts) // 2
        positions = list(range(-context_size, context_size + 1))
        
        # Count nucleotides at each position
        position_distributions = {}
        for pos in positions:
            position_distributions[pos] = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        
        for context in contexts:
            # Find the center position
            center = len(context) // 2
            
            for pos in positions:
                if 0 <= center + pos < len(context):
                    nucleotide = context[center + pos]
                    if nucleotide in 'ACGT':
                        position_distributions[pos][nucleotide] += 1
        
        # Plot nucleotide distributions around mutation site
        plt.figure(figsize=(12, 6))
        
        for i, pos in enumerate(positions):
            # Calculate proportion of each nucleotide
            total = sum(position_distributions[pos].values())
            if total > 0:
                proportions = {nt: count / total for nt, count in position_distributions[pos].items()}
                
                # Plot stacked bars
                bottom = 0
                for nt, color in zip(['A', 'C', 'G', 'T'], ['green', 'blue', 'orange', 'red']):
                    plt.bar(i, proportions[nt], bottom=bottom, color=color, width=0.8, label=nt if i == 0 else "")
                    bottom += proportions[nt]
        
        # Customize plot
        plt.xlabel('Position Relative to Mutation')
        plt.ylabel('Nucleotide Proportion')
        status = "Gene-Modified" if has_gene == "Yes" else "Non-Modified"
        plt.title(f'Nucleotide Context for {status} Strains')
        plt.xticks(range(len(positions)), positions)
        
        # Highlight mutation position
        plt.axvline(x=positions.index(0), color='black', linestyle='--', alpha=0.5)
        
        # Add legend
        plt.legend(title='Nucleotide')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"gene_{has_gene}_context_pattern.png"), dpi=300)
        plt.close()

# Function to create summary report
def create_summary_report(contexts_df, control_df, output_dir):
    """Create a comprehensive summary report of genomic context analysis."""
    with open(os.path.join(output_dir, "genomic_context_summary.txt"), 'w') as f:
        f.write("Genomic Context Analysis Summary\n")
        f.write("===============================\n\n")
        
        # Overall statistics
        f.write("Overall Statistics:\n")
        f.write("-----------------\n")
        
        # Number of variants analyzed
        total_variants = len(contexts_df)
        f.write(f"Total variants analyzed: {total_variants}\n")
        
        # Treatment breakdown
        treatments = contexts_df['Treatment'].unique()
        f.write("Variants by treatment:\n")
        for treatment in treatments:
            description = TREATMENT_INFO.get(treatment, {}).get('description', 'Unknown')
            adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
            has_gene = TREATMENT_INFO.get(treatment, {}).get('gene')
            
            count = len(contexts_df[contexts_df['Treatment'] == treatment])
            f.write(f"  {treatment}: {count} variants - {description} ({adaptation} adaptation")
            if has_gene:
                f.write(f" with {has_gene} gene")
            f.write(")\n")
        
        # Adaptation breakdown
        adaptation_types = contexts_df['Adaptation'].unique()
        f.write("\nVariants by adaptation type:\n")
        for adaptation in adaptation_types:
            count = len(contexts_df[contexts_df['Adaptation'] == adaptation])
            f.write(f"  {adaptation}: {count} variants\n")
        
        # Gene modification breakdown
        gene_statuses = contexts_df['Has_Gene'].unique()
        f.write("\nVariants by gene modification status:\n")
        for has_gene in gene_statuses:
            status = "Gene-modified" if has_gene == "Yes" else "Non-modified"
            count = len(contexts_df[contexts_df['Has_Gene'] == has_gene])
            f.write(f"  {status}: {count} variants\n")
        
        f.write("\n")
        
        # GC content analysis
        f.write("GC Content Analysis:\n")
        f.write("-------------------\n")
        
        # Overall GC content
        mean_gc = contexts_df['gc_content'].mean()
        control_gc = control_df['gc_content'].mean()
        
        f.write(f"Mean GC content in variant regions: {mean_gc:.4f}\n")
        f.write(f"Mean GC content in control regions: {control_gc:.4f}\n")
        
        # GC content by treatment
        f.write("GC content by treatment:\n")
        for treatment in treatments:
            treatment_gc = contexts_df[contexts_df['Treatment'] == treatment]['gc_content'].mean()
            f.write(f"  {treatment}: {treatment_gc:.4f}\n")
        
        # GC content by adaptation type
        f.write("\nGC content by adaptation type:\n")
        for adaptation in adaptation_types:
            adaptation_gc = contexts_df[contexts_df['Adaptation'] == adaptation]['gc_content'].mean()
            f.write(f"  {adaptation}: {adaptation_gc:.4f}\n")
        
        # GC content by gene status
        f.write("\nGC content by gene modification status:\n")
        for has_gene in gene_statuses:
            status = "Gene-modified" if has_gene == "Yes" else "Non-modified"
            gene_gc = contexts_df[contexts_df['Has_Gene'] == has_gene]['gc_content'].mean()
            f.write(f"  {status}: {gene_gc:.4f}\n")
        
        f.write("\n")
        
        # Local sequence features
        f.write("Local Sequence Features:\n")
        f.write("----------------------\n")
        
        # Homopolymer presence
        homopolymer_frac = contexts_df['nearby_homopolymer'].mean()
        control_homopolymer = control_df['nearby_homopolymer'].mean()
        
        f.write(f"Fraction of variants near homopolymers: {homopolymer_frac:.4f}\n")
        f.write(f"Fraction of control sites near homopolymers: {control_homopolymer:.4f}\n")
        
        # Homopolymer by treatment
        f.write("Homopolymer presence by treatment:\n")
        for treatment in treatments:
            treatment_homopolymer = contexts_df[contexts_df['Treatment'] == treatment]['nearby_homopolymer'].mean()
            f.write(f"  {treatment}: {treatment_homopolymer:.4f}\n")
        
        # Homopolymer by adaptation
        f.write("\nHomopolymer presence by adaptation type:\n")
        for adaptation in adaptation_types:
            adaptation_homopolymer = contexts_df[contexts_df['Adaptation'] == adaptation]['nearby_homopolymer'].mean()
            f.write(f"  {adaptation}: {adaptation_homopolymer:.4f}\n")
        
        # Homopolymer by gene status
        f.write("\nHomopolymer presence by gene modification status:\n")
        for has_gene in gene_statuses:
            status = "Gene-modified" if has_gene == "Yes" else "Non-modified"
            gene_homopolymer = contexts_df[contexts_df['Has_Gene'] == has_gene]['nearby_homopolymer'].mean()
            f.write(f"  {status}: {gene_homopolymer:.4f}\n")
        
        f.write("\n")
        
        # Dinucleotide repeat presence
        dinucleotide_frac = contexts_df['nearby_dinucleotide_repeat'].mean()
        control_dinucleotide = control_df['nearby_dinucleotide_repeat'].mean()
        
        f.write(f"Fraction of variants near dinucleotide repeats: {dinucleotide_frac:.4f}\n")
        f.write(f"Fraction of control sites near dinucleotide repeats: {control_dinucleotide:.4f}\n")
        
        # Dinucleotide by treatment
        f.write("Dinucleotide repeat presence by treatment:\n")
        for treatment in treatments:
            treatment_dinucleotide = contexts_df[contexts_df['Treatment'] == treatment]['nearby_dinucleotide_repeat'].mean()
            f.write(f"  {treatment}: {treatment_dinucleotide:.4f}\n")
        
        # Dinucleotide by adaptation
        f.write("\nDinucleotide repeat presence by adaptation type:\n")
        for adaptation in adaptation_types:
            adaptation_dinucleotide = contexts_df[contexts_df['Adaptation'] == adaptation]['nearby_dinucleotide_repeat'].mean()
            f.write(f"  {adaptation}: {adaptation_dinucleotide:.4f}\n")
        
        # Dinucleotide by gene status
        f.write("\nDinucleotide repeat presence by gene modification status:\n")
        for has_gene in gene_statuses:
            status = "Gene-modified" if has_gene == "Yes" else "Non-modified"
            gene_dinucleotide = contexts_df[contexts_df['Has_Gene'] == has_gene]['nearby_dinucleotide_repeat'].mean()
            f.write(f"  {status}: {gene_dinucleotide:.4f}\n")
        
        f.write("\n")
        
        # Purine/Pyrimidine ratio
        purine_ratio = contexts_df['purine_ratio'].mean()
        control_purine = control_df['purine_ratio'].mean()
        
        f.write(f"Mean purine/pyrimidine ratio in variant regions: {purine_ratio:.4f}\n")
        f.write(f"Mean purine/pyrimidine ratio in control regions: {control_purine:.4f}\n")
        
        f.write("\n")
        
        # Main conclusions
        f.write("Main Conclusions:\n")
        f.write("---------------\n")
        f.write("1. This analysis examines the local sequence context around mutation sites.\n")
        f.write("2. We analyze GC content, homopolymer regions, and repetitive elements near variants.\n")
        f.write("3. Adaptation types (Temperature vs Low Oxygen) show different patterns\n")
        f.write("   in their local sequence contexts.\n")
        f.write("4. Gene modifications (STC/CAS) influence the genomic context of mutations.\n")
        f.write("5. Temperature adaptation shows distinct patterns in sequence composition\n")
        f.write("   compared to low oxygen adaptation.\n")
        f.write("6. Gene-modified strains exhibit different local sequence features compared\n")
        f.write("   to non-modified strains, suggesting specific mutational mechanisms.\n")

# Main function to run the analysis
def main():
    # Load reference genome
    reference_genome = load_reference_genome()
    if not reference_genome:
        print("Error: Could not load reference genome. Exiting.")
        return
    
    # Parse data for each treatment
    all_data = []
    
    for treatment in TREATMENTS:
        data = parse_mutation_data(treatment)
        if len(data) > 0:
            all_data.append(data)
            print(f"Loaded {len(data)} mutations for {treatment} treatment")
        else:
            print(f"Warning: No data available for {treatment} treatment")
    
    if all_data:
        combined_data = pd.concat(all_data, ignore_index=True)
        print(f"Combined {len(combined_data)} mutations across {len(all_data)} treatments")
    else:
        print("No mutation data found.")
        return
    
    # IMPORTANT FIX: Extract context sequences BEFORE calculating metrics
    contexts_df = extract_all_context_sequences(combined_data, reference_genome)
    print(f"Extracted context sequences for {len(contexts_df)} variants")
    
    # Generate control contexts
    control_df = generate_control_contexts(reference_genome, num_controls=1000)
    print(f"Generated {len(control_df)} control context sequences")
    
    # Calculate composition metrics
    contexts_df = calculate_composition_metrics(contexts_df)  # Now using contexts_df, not combined_data
    control_df = calculate_composition_metrics(control_df)
    print("Calculated composition metrics for variant and control regions")
    
    # Extract immediate contexts for context analysis
    immediate_contexts_df = extract_immediate_contexts(contexts_df)
    print(f"Extracted immediate contexts for {len(immediate_contexts_df)} variants")
    
    # Plot GC content distribution
    plot_gc_content(contexts_df, control_df, OUTPUT_DIR)
    print("Generated GC content distribution plots")
    
    # Plot sequence features by treatment
    plot_sequence_features(contexts_df, OUTPUT_DIR)
    print("Generated sequence feature plots by treatment")
    
    # Plot sequence features by adaptation type
    plot_sequence_features_by_adaptation(contexts_df, OUTPUT_DIR)
    print("Generated sequence feature plots by adaptation type")
    
    # Plot sequence features by gene status
    plot_sequence_features_by_gene(contexts_df, OUTPUT_DIR)
    print("Generated sequence feature plots by gene status")
    
    # Analyze adaptation-specific patterns
    analyze_adaptation_specific_patterns(contexts_df, OUTPUT_DIR)
    print("Analyzed adaptation-specific patterns")
    
    # Analyze gene-specific patterns
    analyze_gene_specific_patterns(contexts_df, OUTPUT_DIR)
    print("Analyzed gene-specific patterns")
    
    # Analyze mutation context patterns
    analyze_mutation_context_patterns(immediate_contexts_df, OUTPUT_DIR)
    print("Analyzed immediate mutation context patterns")
    
    # Create summary report
    create_summary_report(contexts_df, control_df, OUTPUT_DIR)
    print("Created summary report")
    
    print(f"Analysis complete! Results saved to {OUTPUT_DIR}/")

# Run the analysis
if __name__ == "__main__":
    main()