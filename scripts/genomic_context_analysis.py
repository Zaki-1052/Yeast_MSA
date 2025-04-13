#!/usr/bin/env python3

import os
import re
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

# Set matplotlib style
plt.style.use('ggplot')
sns.set(font_scale=1.2)
sns.set_style("whitegrid")

# Define output directory
OUTPUT_DIR = "genomic_context_results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Define nucleotide complements
COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

# Reference genome path - adjust as needed
REFERENCE_GENOME = "reference/yeast_w303.fasta"

# Function to load reference genome
def load_reference_genome(fasta_file=REFERENCE_GENOME):
    """Load the reference genome from a FASTA file."""
    if not os.path.exists(fasta_file):
        print(f"Error: Reference genome file {fasta_file} not found.")
        return {}
    
    print(f"Loading reference genome from {fasta_file}...")
    reference = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        reference[record.id] = str(record.seq).upper()
    
    print(f"Loaded {len(reference)} sequences from reference genome.")
    return reference

# Function to parse mutation data for a specific treatment
def parse_mutation_data(treatment):
    """Parse mutation data for a specific treatment."""
    # Check for extracted data from previous analysis
    extracted_file = f"mutation_spectrum_analysis/{treatment}_mutations.txt"
    if os.path.exists(extracted_file):
        # Read the already extracted mutation data
        data = pd.read_csv(extracted_file, sep='\t', header=None, 
                           names=['CHROM', 'POS', 'REF', 'ALT'])
        data['Treatment'] = treatment
        return data
    
    # If extracted file doesn't exist, try to extract from VCF
    vcf_file = f"results/merged/analysis/{treatment}_highconf.vcf.gz"
    if not os.path.exists(vcf_file):
        print(f"Warning: File {vcf_file} not found")
        return pd.DataFrame()
    
    # Extract data from VCF using bcftools
    try:
        output = subprocess.check_output(
            f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' {vcf_file}", 
            shell=True).decode('utf-8')
        
        # Parse the output
        rows = []
        for line in output.strip().split('\n'):
            if line:  # Skip empty lines
                parts = line.split('\t')
                if len(parts) == 4:
                    rows.append(parts)
        
        # Create dataframe
        data = pd.DataFrame(rows, columns=['CHROM', 'POS', 'REF', 'ALT'])
        data['POS'] = data['POS'].astype(int)
        data['Treatment'] = treatment
        
        return data
    except:
        print(f"Error extracting data from {vcf_file}")
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

# Function to extract context sequence around a variant
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
        for pos in range(context_size, len(sequence) - context_size):
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
                'immediate_context': immediate_context,
                'variant_index': rel_var_idx
            })
    
    return pd.DataFrame(immediate_contexts)

# Function to standardize mutation context
def standardize_mutation_context(context, var_idx, ref, alt):
    """Standardize mutation context to use pyrimidine reference."""
    if ref in 'CT':  # Reference is already a pyrimidine
        return context, var_idx, ref, alt
    
    # Reference is a purine, need to complement
    comp_context = ''.join(COMPLEMENT.get(base, 'N') for base in context)
    comp_ref = COMPLEMENT.get(ref, 'N')
    comp_alt = COMPLEMENT.get(alt, 'N')
    
    return comp_context, var_idx, comp_ref, comp_alt

# Function to create position-specific context matrix
def create_context_matrix(contexts_df, window=5):
    """Create a position-specific nucleotide frequency matrix."""
    # Initialize matrix dimensions (positions x nucleotides)
    positions = range(-window, window + 1)
    nucleotides = ['A', 'C', 'G', 'T']
    
    # Create empty matrix
    matrix = {pos: {nt: 0 for nt in nucleotides} for pos in positions}
    
    # Group by mutation type
    grouped = contexts_df.groupby(['REF', 'ALT'])
    
    # Process each mutation type separately
    context_matrices = {}
    
    for (ref, alt), group in grouped:
        # Skip if too few contexts
        if len(group) < 5:
            continue
        
        # Reset matrix counts
        for pos in positions:
            for nt in nucleotides:
                matrix[pos][nt] = 0
        
        # Count nucleotides at each position
        for _, row in group.iterrows():
            ctx = row['immediate_context']
            var_idx = row['variant_index']
            
            # Standardize context if needed
            std_ctx, std_idx, std_ref, std_alt = standardize_mutation_context(
                ctx, var_idx, row['REF'], row['ALT'])
            
            # Count nucleotides at each position relative to variant
            for i, pos in enumerate(range(-window, window + 1)):
                abs_pos = std_idx + pos
                if 0 <= abs_pos < len(std_ctx):
                    nt = std_ctx[abs_pos]
                    if nt in nucleotides:
                        matrix[pos][nt] += 1
        
        # Create a copy of the matrix for this mutation type
        mutation_key = f"{ref}>{alt}"
        context_matrices[mutation_key] = {pos: dict(nts) for pos, nts in matrix.items()}
    
    return context_matrices

# Function to plot sequence logos
def plot_sequence_logos(context_matrices, output_dir):
    """Generate sequence logo plots for mutation contexts."""
    for mutation_type, matrix in context_matrices.items():
        # Convert matrix to pandas DataFrame
        df = pd.DataFrame(matrix).T  # Transpose so positions are rows
        
        # Normalize columns to get probabilities
        for pos in df.index:
            total = df.loc[pos].sum()
            if total > 0:
                df.loc[pos] = df.loc[pos] / total
        
        # Create sequence logo
        logo = logomaker.Logo(df, color_scheme='classic')
        
        # Customize plot
        logo.ax.set_xlabel('Position Relative to Mutation')
        logo.ax.set_ylabel('Probability')
        logo.ax.set_title(f'Sequence Context for {mutation_type} Mutations')
        
        # Mark mutation position
        logo.ax.axvline(x=0, color='red', linestyle='--', alpha=0.5)
        
        # Save plot
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"logo_{mutation_type.replace('>', '_to_')}.png"), dpi=300)
        plt.close()

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
    
    # Plot GC content by treatment
    plt.figure(figsize=(12, 6))
    
    # Group by treatment
    treatments = contexts_df['Treatment'].unique()
    
    for treatment in treatments:
        treatment_data = contexts_df[contexts_df['Treatment'] == treatment]
        sns.kdeplot(treatment_data['gc_content'], label=f'{treatment} Treatment')
    
    # Add control
    sns.kdeplot(control_df['gc_content'], label='Control', linestyle='--')
    
    # Customize plot
    plt.xlabel('GC Content')
    plt.ylabel('Density')
    plt.title('GC Content Distribution by Treatment')
    plt.legend()
    
    # Save plot
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "gc_content_by_treatment.png"), dpi=300)
    plt.close()

# Function to plot local sequence features
def plot_sequence_features(contexts_df, output_dir):
    """Plot local sequence features around variants."""
    # Group by treatment
    treatments = contexts_df['Treatment'].unique()
    
    # Plot homopolymer presence
    homopolymer_data = contexts_df.groupby('Treatment')['nearby_homopolymer'].mean()
    
    plt.figure(figsize=(8, 6))
    bars = plt.bar(homopolymer_data.index, homopolymer_data.values, color='skyblue')
    
    # Add value labels
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                f'{height:.2f}', ha='center', va='bottom')
    
    plt.xlabel('Treatment')
    plt.ylabel('Fraction with Nearby Homopolymer')
    plt.title('Homopolymer Presence Near Variants by Treatment')
    plt.ylim(0, max(homopolymer_data.values) * 1.1)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "homopolymer_by_treatment.png"), dpi=300)
    plt.close()
    
    # Plot dinucleotide repeat presence
    dinucleotide_data = contexts_df.groupby('Treatment')['nearby_dinucleotide_repeat'].mean()
    
    plt.figure(figsize=(8, 6))
    bars = plt.bar(dinucleotide_data.index, dinucleotide_data.values, color='lightgreen')
    
    # Add value labels
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                f'{height:.2f}', ha='center', va='bottom')
    
    plt.xlabel('Treatment')
    plt.ylabel('Fraction with Nearby Dinucleotide Repeat')
    plt.title('Dinucleotide Repeat Presence Near Variants by Treatment')
    plt.ylim(0, max(dinucleotide_data.values) * 1.1)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "dinucleotide_repeat_by_treatment.png"), dpi=300)
    plt.close()

# Function to analyze immediate mutation context patterns
def analyze_mutation_context_patterns(contexts_df, output_dir):
    """Analyze and visualize immediate nucleotide context patterns around mutations."""
    # Group by reference and alternate base
    grouped = contexts_df.groupby(['REF', 'ALT'])
    
    # Process each mutation type
    for (ref, alt), group in grouped:
        # Skip if too few contexts
        if len(group) < 10:
            continue
        
        mutation_type = f"{ref}>{alt}"
        
        # Analyze -1 position
        minus1_counts = Counter()
        
        for _, row in group.iterrows():
            ctx = row['immediate_context']
            var_idx = row['variant_index']
            
            # Get -1 position
            if var_idx > 0:
                minus1_base = ctx[var_idx - 1]
                if minus1_base in 'ACGT':
                    minus1_counts[minus1_base] += 1
        
        # Analyze +1 position
        plus1_counts = Counter()
        
        for _, row in group.iterrows():
            ctx = row['immediate_context']
            var_idx = row['variant_index']
            
            # Get +1 position
            if var_idx < len(ctx) - 1:
                plus1_base = ctx[var_idx + 1]
                if plus1_base in 'ACGT':
                    plus1_counts[plus1_base] += 1
        
        # Plot -1 position distribution
        plt.figure(figsize=(12, 5))
        
        plt.subplot(1, 2, 1)
        bases = ['A', 'C', 'G', 'T']
        values = [minus1_counts.get(base, 0) for base in bases]
        plt.bar(bases, values, color=['green', 'blue', 'orange', 'red'])
        plt.xlabel('Nucleotide at -1 Position')
        plt.ylabel('Count')
        plt.title(f'{mutation_type}: Base Before Mutation')
        
        # Plot +1 position distribution
        plt.subplot(1, 2, 2)
        values = [plus1_counts.get(base, 0) for base in bases]
        plt.bar(bases, values, color=['green', 'blue', 'orange', 'red'])
        plt.xlabel('Nucleotide at +1 Position')
        plt.ylabel('Count')
        plt.title(f'{mutation_type}: Base After Mutation')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"adjacent_bases_{ref}_to_{alt}.png"), dpi=300)
        plt.close()

# Function to analyze mutation context by treatment
def analyze_context_by_treatment(contexts_df, output_dir):
    """Analyze and compare sequence contexts between treatments."""
    # Group by treatment and mutation type
    grouped = contexts_df.groupby(['Treatment', 'REF', 'ALT'])
    
    # Count each mutation type by treatment
    counts = {}
    for (treatment, ref, alt), group in grouped:
        mutation_type = f"{ref}>{alt}"
        if treatment not in counts:
            counts[treatment] = {}
        counts[treatment][mutation_type] = len(group)
    
    # Create dataframe for heatmap
    treatments = sorted(contexts_df['Treatment'].unique())
    mutation_types = sorted(set(f"{ref}>{alt}" for ref, alt in 
                            zip(contexts_df['REF'], contexts_df['ALT'])))
    
    data = []
    for mt in mutation_types:
        row = {'Mutation': mt}
        for treatment in treatments:
            row[treatment] = counts.get(treatment, {}).get(mt, 0)
        data.append(row)
    
    df = pd.DataFrame(data)
    df.set_index('Mutation', inplace=True)
    
    # Create heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(df, annot=True, fmt='d', cmap='YlGnBu')
    plt.title('Mutation Types by Treatment')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "mutation_type_by_treatment_heatmap.png"), dpi=300)
    plt.close()
    
    # Analyze GC content by treatment
    plt.figure(figsize=(8, 6))
    
    # Prepare data for boxplot
    data = []
    for treatment in treatments:
        treatment_data = contexts_df[contexts_df['Treatment'] == treatment]
        data.append(treatment_data['gc_content'].values)
    
    plt.boxplot(data, labels=treatments)
    plt.xlabel('Treatment')
    plt.ylabel('GC Content')
    plt.title('GC Content Distribution by Treatment')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "gc_content_boxplot_by_treatment.png"), dpi=300)
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
            count = len(contexts_df[contexts_df['Treatment'] == treatment])
            f.write(f"  {treatment}: {count} variants\n")
        
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
        
        f.write("\n")
        
        # Main conclusions
        f.write("Main Conclusions:\n")
        f.write("---------------\n")
        f.write("1. This analysis examines the local sequence context around mutation sites.\n")
        f.write("2. We analyze GC content, homopolymer regions, and repetitive elements near variants.\n")
        f.write("3. The sequence logos reveal nucleotide preferences around different mutation types.\n")
        f.write("4. Treatment-specific context patterns may provide insights into mutation mechanisms.\n")

# Main function to run the analysis
def main():
    # Define treatments
    treatments = ['WT', 'STC', 'CAS', 'WTA']
    
    # Load reference genome
    reference_genome = load_reference_genome()
    if not reference_genome:
        print("Error: Could not load reference genome. Exiting.")
        return
    
    # Calculate genome-wide GC content
    genome_gc = 0
    genome_length = 0
    for scaffold, sequence in reference_genome.items():
        genome_gc += sum(1 for base in sequence if base in 'GC')
        genome_length += len(sequence)
    
    genome_gc_content = genome_gc / genome_length if genome_length > 0 else 0
    print(f"Genome-wide GC content: {genome_gc_content:.4f}")
    
    # Parse data for each treatment
    all_data = pd.DataFrame()
    for treatment in treatments:
        data = parse_mutation_data(treatment)
        if len(data) > 0:
            all_data = pd.concat([all_data, data])
            print(f"Loaded {len(data)} variants for {treatment} treatment")
        else:
            print(f"Warning: No data available for {treatment} treatment")
    
    # Filter for SNVs
    snv_data = filter_snvs(all_data)
    print(f"Filtered to {len(snv_data)} single nucleotide variants")
    
    # Extract context sequences
    contexts_df = extract_all_context_sequences(snv_data, reference_genome)
    print(f"Extracted context sequences for {len(contexts_df)} variants")
    
    # Generate control contexts
    control_df = generate_control_contexts(reference_genome, num_controls=1000)
    print(f"Generated {len(control_df)} control context sequences")
    
    # Calculate composition metrics
    contexts_df = calculate_composition_metrics(contexts_df)
    control_df = calculate_composition_metrics(control_df)
    print("Calculated composition metrics for variant and control regions")
    
    # Extract immediate contexts for mutation context analysis
    immediate_contexts_df = extract_immediate_contexts(contexts_df)
    print(f"Extracted immediate contexts for {len(immediate_contexts_df)} variants")
    
    # Create context matrices for sequence logos
    context_matrices = create_context_matrix(immediate_contexts_df)
    print(f"Created context matrices for {len(context_matrices)} mutation types")
    
    # Plot sequence logos
    plot_sequence_logos(context_matrices, OUTPUT_DIR)
    print("Generated sequence logo plots")
    
    # Plot GC content distribution
    plot_gc_content(contexts_df, control_df, OUTPUT_DIR)
    print("Generated GC content distribution plots")
    
    # Plot sequence features
    plot_sequence_features(contexts_df, OUTPUT_DIR)
    print("Generated sequence feature plots")
    
    # Analyze mutation context patterns
    analyze_mutation_context_patterns(immediate_contexts_df, OUTPUT_DIR)
    print("Analyzed immediate mutation context patterns")
    
    # Analyze context by treatment
    analyze_context_by_treatment(contexts_df, OUTPUT_DIR)
    print("Analyzed contexts by treatment")
    
    # Create summary report
    create_summary_report(contexts_df, control_df, OUTPUT_DIR)
    print("Created summary report")
    
    print(f"Analysis complete! Results saved to {OUTPUT_DIR}/")

# Run the analysis
if __name__ == "__main__":
    main()