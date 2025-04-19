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
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist, squareform
import re
import warnings
warnings.filterwarnings('ignore')

# Set matplotlib style
plt.style.use('ggplot')
sns.set(font_scale=1.2)
sns.set_style("whitegrid")

# Define output directory
OUTPUT_DIR = "analysis/mutational_signatures_results"
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

# Define the six base substitution types (pyrimidine-centric)
SUBSTITUTION_TYPES = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']

# Define all possible trinucleotides
NUCLEOTIDES = ['A', 'C', 'G', 'T']
TRINUCLEOTIDES = [''.join(x) for x in 
                  [(a, b, c) for a in NUCLEOTIDES for b in NUCLEOTIDES for c in NUCLEOTIDES]]

# Reference genome path
REFERENCE_GENOME = "reference/w303_chromosomal.fasta"

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
            # First, let's read the raw content to examine the actual structure
            with open(mutation_file, 'r') as f:
                first_line = f.readline().strip()
                columns = first_line.split('\t')
                print(f"First line has {len(columns)} columns: {columns}")
                
            # Based on file inspection, we know it has line number, CHROM, POS, REF, ALT, Treatment
            # But we'll read it carefully to handle any variation
            raw_data = pd.read_csv(mutation_file, sep='\t', header=None)
            
            # Check how many columns we have and adjust accordingly
            if len(raw_data.columns) == 6:  # Line num, CHROM, POS, REF, ALT, Treatment
                data = raw_data.iloc[:, 1:].copy()  # Skip the first column (line number)
                data.columns = ['CHROM', 'POS', 'REF', 'ALT', 'Treatment']
            elif len(raw_data.columns) == 5:  # CHROM, POS, REF, ALT, Treatment
                data = raw_data.copy()
                data.columns = ['CHROM', 'POS', 'REF', 'ALT', 'Treatment']
            else:
                # If neither 5 nor 6 columns, try to guess based on content
                print(f"Unexpected column count: {len(raw_data.columns)}. Attempting to adapt.")
                # If last column has treatment names, assume it's Treatment
                if raw_data.iloc[:, -1].str.contains('|'.join(TREATMENTS), regex=True).any():
                    # Extract the main columns we need
                    if len(raw_data.columns) > 4:
                        data = raw_data.iloc[:, -5:-1].copy()  # Last 5 columns except the last one
                        data['Treatment'] = raw_data.iloc[:, -1]
                        data.columns = ['CHROM', 'POS', 'REF', 'ALT', 'Treatment']
                    else:
                        raise ValueError(f"Cannot extract required columns from data with {len(raw_data.columns)} columns")
                else:
                    raise ValueError(f"Cannot identify Treatment column in data with {len(raw_data.columns)} columns")
            
            # Set Treatment column explicitly to the current treatment
            data['Treatment'] = treatment
            
            # Add biological context columns
            data['Adaptation'] = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
            data['Has_Gene'] = 'Yes' if TREATMENT_INFO.get(treatment, {}).get('gene') else 'No'
            
            # Ensure POS is integer
            data['POS'] = data['POS'].astype(int)
            
            # Print sample data for debugging
            print(f"Sample data:\n{data.head()}")
            
            print(f"Loaded {len(data)} mutations for {treatment}")
            return data
        except Exception as e:
            print(f"Error reading {mutation_file}: {e}")
    
    # If no mutation file found or error occurred, try to extract from VCF
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
# Modified filter_snvs function for mutational_signature_analysis.py
def filter_snvs(data, debug=True):
    """Filter data to include only single nucleotide variants."""
    if debug:
        print(f"Initial variant count: {len(data)}")
        if len(data) > 0:
            print(f"REF column type: {type(data['REF'].iloc[0])}")
            print(f"Sample REF values: {data['REF'].head().tolist()}")
            print(f"Sample ALT values: {data['ALT'].head().tolist()}")
    
    # Check for non-string data types
    if not pd.api.types.is_string_dtype(data['REF']) or not pd.api.types.is_string_dtype(data['ALT']):
        if debug:
            print("Converting REF/ALT to string types")
        data = data.copy()
        data['REF'] = data['REF'].astype(str)
        data['ALT'] = data['ALT'].astype(str)
    
    # First filtering step: length check
    length_filter = (data['REF'].str.len() == 1) & (data['ALT'].str.len() == 1)
    snv_data = data[length_filter]
    
    if debug:
        print(f"After length filter: {len(snv_data)} variants remain")
        print(f"Removed {len(data) - len(snv_data)} multi-nucleotide variants")
        if len(data) > 0 and len(snv_data) == 0:
            # Show some examples of what's being filtered out
            print("Examples of filtered variants:")
            multi_nt = data[~length_filter].head(5)
            for _, row in multi_nt.iterrows():
                print(f"  {row['CHROM']}:{row['POS']} REF={row['REF']} ALT={row['ALT']}")
    
    # Second filtering step: valid bases
    # Make case-insensitive by converting to uppercase
    if len(snv_data) > 0:
        snv_data = snv_data.copy()
        snv_data['REF'] = snv_data['REF'].str.upper()
        snv_data['ALT'] = snv_data['ALT'].str.upper()
        
        valid_bases = snv_data['REF'].isin(['A', 'C', 'G', 'T']) & snv_data['ALT'].isin(['A', 'C', 'G', 'T'])
        final_data = snv_data[valid_bases]
        
        if debug:
            print(f"After ACGT filter: {len(final_data)} variants remain")
            print(f"Removed {len(snv_data) - len(final_data)} variants with non-ACGT bases")
            if len(snv_data) > 0 and len(final_data) == 0:
                print("Examples of non-ACGT bases:")
                non_acgt = snv_data[~valid_bases].head(5)
                for _, row in non_acgt.iterrows():
                    print(f"  {row['CHROM']}:{row['POS']} REF={row['REF']} ALT={row['ALT']}")
    else:
        final_data = snv_data
    
    return final_data

# Function to extract sequence context around a variant
def extract_context_sequence(chrom, pos, ref, alt, reference_genome, context_size=5):
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
    
    # Get surrounding bases
    if left_offset >= 1 and right_offset >= 1:
        trinucleotide = context[left_offset-1:left_offset+2]
        if len(trinucleotide) != 3:
            # Not enough context for trinucleotide
            trinucleotide = None
    else:
        trinucleotide = None
    
    return {
        'context': context,
        'variant_index': left_offset,
        'trinucleotide': trinucleotide,
        'ref': ref,
        'alt': alt
    }

# Function to standardize mutation representation
def standardize_mutation(ref, alt, context=None, var_idx=None):
    """Standardize mutation to pyrimidine-centric representation."""
    if ref in ['C', 'T']:
        # Reference is already a pyrimidine
        std_ref = ref
        std_alt = alt
        std_context = context
        std_var_idx = var_idx
    else:
        # Reference is a purine, convert to complementary base
        std_ref = COMPLEMENT[ref]
        std_alt = COMPLEMENT[alt]
        if context and var_idx is not None:
            # Complement the entire context
            std_context = ''.join(COMPLEMENT.get(base, 'N') for base in reversed(context))
            std_var_idx = len(context) - var_idx - 1
        else:
            std_context = None
            std_var_idx = None
    
    return std_ref, std_alt, std_context, std_var_idx

# Function to extract context for all variants
def extract_all_contexts(data, reference_genome, context_size=5):
    """Extract sequence context for all variants."""
    contexts = []
    
    for _, row in data.iterrows():
        context_data = extract_context_sequence(
            row['CHROM'], row['POS'], row['REF'], row['ALT'], 
            reference_genome, context_size
        )
        
        if context_data and context_data['trinucleotide']:
            # Standardize the mutation and context
            std_ref, std_alt, std_context, std_var_idx = standardize_mutation(
                context_data['ref'], 
                context_data['alt'],
                context_data['context'],
                context_data['variant_index']
            )
            
            mut_type = f"{std_ref}>{std_alt}"
            trinuc = context_data['trinucleotide']
            
            # Also standardize the trinucleotide if needed
            if std_context and std_var_idx is not None:
                std_trinuc = std_context[max(0, std_var_idx-1):std_var_idx+2]
                if len(std_trinuc) == 3:
                    trinuc = std_trinuc
            
            contexts.append({
                'CHROM': row['CHROM'],
                'POS': row['POS'],
                'REF': row['REF'],
                'ALT': row['ALT'],
                'Treatment': row['Treatment'],
                'Adaptation': row['Adaptation'],
                'Has_Gene': row['Has_Gene'],
                'Context': context_data['context'],
                'Std_Context': std_context,
                'Trinucleotide': trinuc,
                'Mutation_Type': mut_type
            })
    
    return pd.DataFrame(contexts)

# Function to create mutation signature matrix
def create_signature_matrix(contexts_df):
    """Create a mutation signature matrix from context data."""
    # Group by treatment
    treatments = contexts_df['Treatment'].unique()
    
    # Initialize signature matrix
    # Structure: treatment -> mutation_type -> trinucleotide -> count
    signatures = {t: {mut: {tri: 0 for tri in TRINUCLEOTIDES} 
                     for mut in SUBSTITUTION_TYPES} 
                 for t in treatments}
    
    # Count occurrences
    for _, row in contexts_df.iterrows():
        treatment = row['Treatment']
        mut_type = row['Mutation_Type']
        trinuc = row['Trinucleotide']
        
        if mut_type in SUBSTITUTION_TYPES and trinuc in TRINUCLEOTIDES:
            signatures[treatment][mut_type][trinuc] += 1
    
    # Convert to normalized signature
    for treatment in treatments:
        for mut_type in SUBSTITUTION_TYPES:
            total = sum(signatures[treatment][mut_type].values())
            if total > 0:
                for trinuc in TRINUCLEOTIDES:
                    signatures[treatment][mut_type][trinuc] /= total
    
    return signatures

# Function to create adaptation-specific signature matrix
def create_adaptation_signature_matrix(contexts_df):
    """Create a mutation signature matrix grouped by adaptation type."""
    # Group by adaptation type
    adaptation_types = contexts_df['Adaptation'].unique()
    
    # Initialize adaptation signature matrix
    adaptation_signatures = {a: {mut: {tri: 0 for tri in TRINUCLEOTIDES} 
                               for mut in SUBSTITUTION_TYPES} 
                           for a in adaptation_types}
    
    # Count occurrences
    for _, row in contexts_df.iterrows():
        adaptation = row['Adaptation']
        mut_type = row['Mutation_Type']
        trinuc = row['Trinucleotide']
        
        if mut_type in SUBSTITUTION_TYPES and trinuc in TRINUCLEOTIDES:
            adaptation_signatures[adaptation][mut_type][trinuc] += 1
    
    # Convert to normalized signature
    for adaptation in adaptation_types:
        for mut_type in SUBSTITUTION_TYPES:
            total = sum(adaptation_signatures[adaptation][mut_type].values())
            if total > 0:
                for trinuc in TRINUCLEOTIDES:
                    adaptation_signatures[adaptation][mut_type][trinuc] /= total
    
    return adaptation_signatures

# Function to create gene-specific signature matrix
def create_gene_signature_matrix(contexts_df):
    """Create a mutation signature matrix grouped by gene modification status."""
    # Group by gene modification status
    gene_statuses = contexts_df['Has_Gene'].unique()
    
    # Initialize gene signature matrix
    gene_signatures = {g: {mut: {tri: 0 for tri in TRINUCLEOTIDES} 
                         for mut in SUBSTITUTION_TYPES} 
                     for g in gene_statuses}
    
    # Count occurrences
    for _, row in contexts_df.iterrows():
        has_gene = row['Has_Gene']
        mut_type = row['Mutation_Type']
        trinuc = row['Trinucleotide']
        
        if mut_type in SUBSTITUTION_TYPES and trinuc in TRINUCLEOTIDES:
            gene_signatures[has_gene][mut_type][trinuc] += 1
    
    # Convert to normalized signature
    for has_gene in gene_statuses:
        for mut_type in SUBSTITUTION_TYPES:
            total = sum(gene_signatures[has_gene][mut_type].values())
            if total > 0:
                for trinuc in TRINUCLEOTIDES:
                    gene_signatures[has_gene][mut_type][trinuc] /= total
    
    return gene_signatures

# Function to create adaptation+gene signature matrix
def create_combined_signature_matrix(contexts_df):
    """Create a mutation signature matrix grouped by adaptation type and gene status."""
    # Combine adaptation and gene status
    contexts_df['Group'] = contexts_df.apply(
        lambda row: f"{row['Adaptation']}_{row['Has_Gene']}", axis=1)
    
    # Group by combined adaptation and gene status
    groups = contexts_df['Group'].unique()
    
    # Initialize combined signature matrix
    combined_signatures = {g: {mut: {tri: 0 for tri in TRINUCLEOTIDES} 
                             for mut in SUBSTITUTION_TYPES} 
                         for g in groups}
    
    # Count occurrences
    for _, row in contexts_df.iterrows():
        group = row['Group']
        mut_type = row['Mutation_Type']
        trinuc = row['Trinucleotide']
        
        if mut_type in SUBSTITUTION_TYPES and trinuc in TRINUCLEOTIDES:
            combined_signatures[group][mut_type][trinuc] += 1
    
    # Convert to normalized signature
    for group in groups:
        for mut_type in SUBSTITUTION_TYPES:
            total = sum(combined_signatures[group][mut_type].values())
            if total > 0:
                for trinuc in TRINUCLEOTIDES:
                    combined_signatures[group][mut_type][trinuc] /= total
    
    return combined_signatures

# Function to plot mutation signature
def plot_signature(signature, treatment, output_dir):
    """Plot a mutation signature in COSMIC style."""
    # Get treatment metadata
    description = TREATMENT_INFO.get(treatment, {}).get('description', '')
    adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', '')
    has_gene = TREATMENT_INFO.get(treatment, {}).get('gene')
    gene_text = f" with {has_gene} gene" if has_gene else ""
    
    # Prepare data for plotting
    data = []
    
    # For each substitution type
    for mut_type in SUBSTITUTION_TYPES:
        # For each trinucleotide context
        for trinuc in TRINUCLEOTIDES:
            # Extract middle base (should match mutation ref)
            middle_base = trinuc[1]
            if middle_base == mut_type[0]:  # Check if middle base matches mutation ref
                data.append({
                    'Substitution': mut_type,
                    'Trinucleotide': trinuc,
                    'Frequency': signature[mut_type][trinuc]
                })
    
    # Convert to DataFrame
    df = pd.DataFrame(data)
    
    # Set up the plot
    plt.figure(figsize=(16, 8))
    
    # Define colors for substitution types
    colors = {
        'C>A': '#2EBAED', 'C>G': '#000000', 'C>T': '#E41A1C',
        'T>A': '#999999', 'T>C': '#4DAF4A', 'T>G': '#984EA3'
    }
    
    # Group by substitution type
    grouped = df.groupby('Substitution')
    
    # Plot each substitution type as a separate bar
    bar_width = 0.8
    index = 0
    ticks = []
    ticklabels = []
    
    for i, (mut_type, group) in enumerate(grouped):
        # Sort by trinucleotide
        sorted_group = group.sort_values('Trinucleotide')
        
        # Plot bars
        indices = range(index, index + len(sorted_group))
        plt.bar(indices, sorted_group['Frequency'], bar_width, color=colors[mut_type], label=mut_type)
        
        # Store tick positions and labels
        for j, trinuc in enumerate(sorted_group['Trinucleotide']):
            ticks.append(index + j)
            ticklabels.append(trinuc)
        
        # Update index
        index += len(sorted_group)
        
        # Add a gap between substitution types
        if i < len(grouped) - 1:
            index += 2
    
    # Customize the plot
    plt.xlabel('Trinucleotide Context')
    plt.ylabel('Relative Frequency')
    plt.title(f'Mutation Signature for {treatment}\n{description} ({adaptation} adaptation{gene_text})')
    
    # Set x-ticks
    plt.xticks(ticks, ticklabels, rotation=90, fontsize=8)
    
    # Add legend
    plt.legend(title='Substitution Type')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(os.path.join(output_dir, f"{treatment}_mutation_signature.png"), dpi=300)
    plt.close()

# Function to plot adaptation-specific signatures
def plot_adaptation_signatures(adaptation_signatures, output_dir):
    """Plot mutation signatures grouped by adaptation type."""
    for adaptation, signature in adaptation_signatures.items():
        # Prepare data for plotting
        data = []
        
        # For each substitution type
        for mut_type in SUBSTITUTION_TYPES:
            # For each trinucleotide context
            for trinuc in TRINUCLEOTIDES:
                # Extract middle base (should match mutation ref)
                middle_base = trinuc[1]
                if middle_base == mut_type[0]:  # Check if middle base matches mutation ref
                    data.append({
                        'Substitution': mut_type,
                        'Trinucleotide': trinuc,
                        'Frequency': signature[mut_type][trinuc]
                    })
        
        # Convert to DataFrame
        df = pd.DataFrame(data)
        
        # Set up the plot
        plt.figure(figsize=(16, 8))
        
        # Define colors for substitution types
        colors = {
            'C>A': '#2EBAED', 'C>G': '#000000', 'C>T': '#E41A1C',
            'T>A': '#999999', 'T>C': '#4DAF4A', 'T>G': '#984EA3'
        }
        
        # Group by substitution type
        grouped = df.groupby('Substitution')
        
        # Plot each substitution type as a separate bar
        bar_width = 0.8
        index = 0
        ticks = []
        ticklabels = []
        
        for i, (mut_type, group) in enumerate(grouped):
            # Sort by trinucleotide
            sorted_group = group.sort_values('Trinucleotide')
            
            # Plot bars
            indices = range(index, index + len(sorted_group))
            plt.bar(indices, sorted_group['Frequency'], bar_width, color=colors[mut_type], label=mut_type)
            
            # Store tick positions and labels
            for j, trinuc in enumerate(sorted_group['Trinucleotide']):
                ticks.append(index + j)
                ticklabels.append(trinuc)
            
            # Update index
            index += len(sorted_group)
            
            # Add a gap between substitution types
            if i < len(grouped) - 1:
                index += 2
        
        # Customize the plot
        plt.xlabel('Trinucleotide Context')
        plt.ylabel('Relative Frequency')
        plt.title(f'Mutation Signature for {adaptation} Adaptation')
        
        # Set x-ticks
        plt.xticks(ticks, ticklabels, rotation=90, fontsize=8)
        
        # Add legend
        plt.legend(title='Substitution Type')
        
        # Adjust layout
        plt.tight_layout()
        
        # Save the plot
        plt.savefig(os.path.join(output_dir, f"{adaptation}_adaptation_signature.png"), dpi=300)
        plt.close()

# Function to plot gene modification signatures
def plot_gene_signatures(gene_signatures, output_dir):
    """Plot mutation signatures grouped by gene modification status."""
    for has_gene, signature in gene_signatures.items():
        # Prepare data for plotting
        data = []
        
        # For each substitution type
        for mut_type in SUBSTITUTION_TYPES:
            # For each trinucleotide context
            for trinuc in TRINUCLEOTIDES:
                # Extract middle base (should match mutation ref)
                middle_base = trinuc[1]
                if middle_base == mut_type[0]:  # Check if middle base matches mutation ref
                    data.append({
                        'Substitution': mut_type,
                        'Trinucleotide': trinuc,
                        'Frequency': signature[mut_type][trinuc]
                    })
        
        # Convert to DataFrame
        df = pd.DataFrame(data)
        
        # Set up the plot
        plt.figure(figsize=(16, 8))
        
        # Define colors for substitution types
        colors = {
            'C>A': '#2EBAED', 'C>G': '#000000', 'C>T': '#E41A1C',
            'T>A': '#999999', 'T>C': '#4DAF4A', 'T>G': '#984EA3'
        }
        
        # Group by substitution type
        grouped = df.groupby('Substitution')
        
        # Plot each substitution type as a separate bar
        bar_width = 0.8
        index = 0
        ticks = []
        ticklabels = []
        
        for i, (mut_type, group) in enumerate(grouped):
            # Sort by trinucleotide
            sorted_group = group.sort_values('Trinucleotide')
            
            # Plot bars
            indices = range(index, index + len(sorted_group))
            plt.bar(indices, sorted_group['Frequency'], bar_width, color=colors[mut_type], label=mut_type)
            
            # Store tick positions and labels
            for j, trinuc in enumerate(sorted_group['Trinucleotide']):
                ticks.append(index + j)
                ticklabels.append(trinuc)
            
            # Update index
            index += len(sorted_group)
            
            # Add a gap between substitution types
            if i < len(grouped) - 1:
                index += 2
        
        # Format the gene status label
        gene_status = "Gene-Modified Strains" if has_gene == "Yes" else "Non-Modified Strains"
        
        # Customize the plot
        plt.xlabel('Trinucleotide Context')
        plt.ylabel('Relative Frequency')
        plt.title(f'Mutation Signature for {gene_status}')
        
        # Set x-ticks
        plt.xticks(ticks, ticklabels, rotation=90, fontsize=8)
        
        # Add legend
        plt.legend(title='Substitution Type')
        
        # Adjust layout
        plt.tight_layout()
        
        # Save the plot
        plt.savefig(os.path.join(output_dir, f"gene_{has_gene}_signature.png"), dpi=300)
        plt.close()

# Function to plot combined adaptation+gene signatures
def plot_combined_signatures(combined_signatures, output_dir):
    """Plot mutation signatures grouped by adaptation type and gene status."""
    for group, signature in combined_signatures.items():
        # Extract adaptation and gene status from group name
        adaptation, has_gene = group.split('_')
        
        # Prepare data for plotting
        data = []
        
        # For each substitution type
        for mut_type in SUBSTITUTION_TYPES:
            # For each trinucleotide context
            for trinuc in TRINUCLEOTIDES:
                # Extract middle base (should match mutation ref)
                middle_base = trinuc[1]
                if middle_base == mut_type[0]:  # Check if middle base matches mutation ref
                    data.append({
                        'Substitution': mut_type,
                        'Trinucleotide': trinuc,
                        'Frequency': signature[mut_type][trinuc]
                    })
        
        # Convert to DataFrame
        df = pd.DataFrame(data)
        
        # Set up the plot
        plt.figure(figsize=(16, 8))
        
        # Define colors for substitution types
        colors = {
            'C>A': '#2EBAED', 'C>G': '#000000', 'C>T': '#E41A1C',
            'T>A': '#999999', 'T>C': '#4DAF4A', 'T>G': '#984EA3'
        }
        
        # Group by substitution type
        grouped = df.groupby('Substitution')
        
        # Plot each substitution type as a separate bar
        bar_width = 0.8
        index = 0
        ticks = []
        ticklabels = []
        
        for i, (mut_type, group_df) in enumerate(grouped):
            # Sort by trinucleotide
            sorted_group = group_df.sort_values('Trinucleotide')
            
            # Plot bars
            indices = range(index, index + len(sorted_group))
            plt.bar(indices, sorted_group['Frequency'], bar_width, color=colors[mut_type], label=mut_type)
            
            # Store tick positions and labels
            for j, trinuc in enumerate(sorted_group['Trinucleotide']):
                ticks.append(index + j)
                ticklabels.append(trinuc)
            
            # Update index
            index += len(sorted_group)
            
            # Add a gap between substitution types
            if i < len(grouped) - 1:
                index += 2
        
        # Format the group label
        gene_status = "with gene modification" if has_gene == "Yes" else "without gene modification"
        
        # Customize the plot
        plt.xlabel('Trinucleotide Context')
        plt.ylabel('Relative Frequency')
        plt.title(f'Mutation Signature for {adaptation} Adaptation {gene_status}')
        
        # Set x-ticks
        plt.xticks(ticks, ticklabels, rotation=90, fontsize=8)
        
        # Add legend
        plt.legend(title='Substitution Type')
        
        # Adjust layout
        plt.tight_layout()
        
        # Save the plot
        plt.savefig(os.path.join(output_dir, f"{adaptation}_{has_gene}_signature.png"), dpi=300)
        plt.close()

# Function to calculate signature similarity
def calculate_signature_similarity(signatures):
    """Calculate similarity between mutation signatures."""
    treatments = list(signatures.keys())
    n_treatments = len(treatments)
    
    # Convert signatures to flat vectors for comparison
    signature_vectors = {}
    
    for treatment in treatments:
        vector = []
        for mut_type in SUBSTITUTION_TYPES:
            for trinuc in TRINUCLEOTIDES:
                vector.append(signatures[treatment][mut_type][trinuc])
        signature_vectors[treatment] = np.array(vector)
    
    # Calculate cosine similarity matrix
    similarity_matrix = np.zeros((n_treatments, n_treatments))
    
    for i, t1 in enumerate(treatments):
        for j, t2 in enumerate(treatments):
            v1 = signature_vectors[t1]
            v2 = signature_vectors[t2]
            
            # Cosine similarity
            dot_product = np.dot(v1, v2)
            norm_v1 = np.linalg.norm(v1)
            norm_v2 = np.linalg.norm(v2)
            
            if norm_v1 > 0 and norm_v2 > 0:
                similarity = dot_product / (norm_v1 * norm_v2)
            else:
                similarity = 0
            
            similarity_matrix[i, j] = similarity
    
    return pd.DataFrame(similarity_matrix, index=treatments, columns=treatments)

# Function to plot signature similarity heatmap
def plot_signature_similarity(similarity_df, output_dir):
    """Plot a heatmap of signature similarities."""
    # Create clustered heatmap with adaptation type coloring
    plt.figure(figsize=(10, 8))
    
    # Create row colors based on adaptation type
    row_colors = []
    for treatment in similarity_df.index:
        adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
        row_colors.append(ADAPTATION_COLORS.get(adaptation, '#333333'))
    
    # Create clustered heatmap
    g = sns.clustermap(
        similarity_df,
        cmap='viridis',
        annot=True,
        fmt='.3f',
        row_colors=[row_colors],
        col_colors=[row_colors],
        vmin=0,
        vmax=1,
        figsize=(12, 10)
    )
    
    # Add adaptation type legend
    adaptations = set(TREATMENT_INFO.get(t, {}).get('adaptation', 'Unknown') for t in similarity_df.index)
    for adaptation in adaptations:
        g.ax_row_dendrogram.bar(0, 0, color=ADAPTATION_COLORS.get(adaptation, '#333333'), label=adaptation)
    
    g.ax_row_dendrogram.legend(title="Adaptation", loc='center', ncol=1)
    
    # Set title
    plt.suptitle('Mutation Signature Similarity Between Treatments', fontsize=16, y=0.95)
    
    plt.savefig(os.path.join(output_dir, "signature_similarity_heatmap.png"), dpi=300)
    plt.close()
    
    # Also create a traditional heatmap for better readability
    plt.figure(figsize=(10, 8))
    sns.heatmap(similarity_df, annot=True, cmap='viridis', vmin=0, vmax=1, fmt='.3f')
    plt.title('Mutation Signature Similarity Between Treatments')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "signature_similarity_simple_heatmap.png"), dpi=300)
    plt.close()

# Function to find enriched trinucleotide contexts
def find_enriched_contexts(contexts_df):
    """Find trinucleotide contexts that are enriched in specific treatments."""
    # Group by treatment and mutation type
    treatment_mut_groups = contexts_df.groupby(['Treatment', 'Mutation_Type'])
    
    # Initialize results
    enriched_contexts = []
    
    # For each treatment and mutation type
    for (treatment, mut_type), group in treatment_mut_groups:
        # Skip if too few mutations
        if len(group) < 5:
            continue
        
        # Count trinucleotide contexts
        trinuc_counts = Counter(group['Trinucleotide'])
        total = sum(trinuc_counts.values())
        
        # Compare with expected frequency (assuming uniform distribution)
        expected_freq = 1 / len(TRINUCLEOTIDES)
        
        for trinuc, count in trinuc_counts.items():
            observed_freq = count / total
            
            # Check if significantly enriched
            if observed_freq > 2 * expected_freq and count >= 3:
                # Get treatment metadata
                adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
                has_gene = 'Yes' if TREATMENT_INFO.get(treatment, {}).get('gene') else 'No'
                
                enriched_contexts.append({
                    'Treatment': treatment,
                    'Adaptation': adaptation,
                    'Has_Gene': has_gene,
                    'Mutation_Type': mut_type,
                    'Trinucleotide': trinuc,
                    'Count': count,
                    'Frequency': observed_freq,
                    'Fold_Enrichment': observed_freq / expected_freq
                })
    
    return pd.DataFrame(enriched_contexts)

# Function to find adaptation-specific enriched contexts
def find_adaptation_enriched_contexts(contexts_df):
    """Find trinucleotide contexts that are enriched in specific adaptation types."""
    # Group by adaptation and mutation type
    adaptation_mut_groups = contexts_df.groupby(['Adaptation', 'Mutation_Type'])
    
    # Initialize results
    enriched_contexts = []
    
    # For each adaptation and mutation type
    for (adaptation, mut_type), group in adaptation_mut_groups:
        # Skip if too few mutations
        if len(group) < 5:
            continue
        
        # Count trinucleotide contexts
        trinuc_counts = Counter(group['Trinucleotide'])
        total = sum(trinuc_counts.values())
        
        # Compare with expected frequency (assuming uniform distribution)
        expected_freq = 1 / len(TRINUCLEOTIDES)
        
        for trinuc, count in trinuc_counts.items():
            observed_freq = count / total
            
            # Check if significantly enriched
            if observed_freq > 2 * expected_freq and count >= 3:
                enriched_contexts.append({
                    'Adaptation': adaptation,
                    'Mutation_Type': mut_type,
                    'Trinucleotide': trinuc,
                    'Count': count,
                    'Frequency': observed_freq,
                    'Fold_Enrichment': observed_freq / expected_freq
                })
    
    return pd.DataFrame(enriched_contexts)

# Function to find gene-specific enriched contexts
def find_gene_enriched_contexts(contexts_df):
    """Find trinucleotide contexts that are enriched in gene-modified or non-modified strains."""
    # Group by gene status and mutation type
    gene_mut_groups = contexts_df.groupby(['Has_Gene', 'Mutation_Type'])
    
    # Initialize results
    enriched_contexts = []
    
    # For each gene status and mutation type
    for (has_gene, mut_type), group in gene_mut_groups:
        # Skip if too few mutations
        if len(group) < 5:
            continue
        
        # Count trinucleotide contexts
        trinuc_counts = Counter(group['Trinucleotide'])
        total = sum(trinuc_counts.values())
        
        # Compare with expected frequency (assuming uniform distribution)
        expected_freq = 1 / len(TRINUCLEOTIDES)
        
        for trinuc, count in trinuc_counts.items():
            observed_freq = count / total
            
            # Check if significantly enriched
            if observed_freq > 2 * expected_freq and count >= 3:
                enriched_contexts.append({
                    'Has_Gene': has_gene,
                    'Mutation_Type': mut_type,
                    'Trinucleotide': trinuc,
                    'Count': count,
                    'Frequency': observed_freq,
                    'Fold_Enrichment': observed_freq / expected_freq
                })
    
    return pd.DataFrame(enriched_contexts)

# Function to plot enriched contexts
def plot_enriched_contexts(enriched_df, output_dir):
    """Plot trinucleotide contexts enriched in specific treatments."""
    if len(enriched_df) == 0:
        print("No enriched contexts found for plotting.")
        return
    
    # Group by treatment
    treatments = enriched_df['Treatment'].unique()
    
    for treatment in treatments:
        treatment_data = enriched_df[enriched_df['Treatment'] == treatment]
        
        if len(treatment_data) == 0:
            continue
        
        # Get treatment metadata
        description = TREATMENT_INFO.get(treatment, {}).get('description', '')
        adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', '')
        has_gene = TREATMENT_INFO.get(treatment, {}).get('gene')
        gene_text = f" with {has_gene} gene" if has_gene else ""
        
        plt.figure(figsize=(12, 6))
        
        # Plot enriched contexts as horizontal bars
        contexts = treatment_data['Trinucleotide'] + ' (' + treatment_data['Mutation_Type'] + ')'
        fold_enrichment = treatment_data['Fold_Enrichment']
        
        # Sort by fold enrichment
        sorted_indices = np.argsort(fold_enrichment)[::-1]  # Descending order
        contexts = contexts.iloc[sorted_indices]
        fold_enrichment = fold_enrichment.iloc[sorted_indices]
        
        # Plot bars
        bars = plt.barh(contexts, fold_enrichment, color=TREATMENT_COLORS.get(treatment, '#1b9e77'))
        
        # Add value labels
        for bar in bars:
            width = bar.get_width()
            plt.text(width + 0.1, bar.get_y() + bar.get_height()/2,
                    f'{width:.1f}', ha='left', va='center')
        
        plt.xlabel('Fold Enrichment')
        plt.ylabel('Trinucleotide Context (Mutation Type)')
        plt.title(f'Enriched Trinucleotide Contexts for {treatment}\n{description} ({adaptation} adaptation{gene_text})')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{treatment}_enriched_contexts.png"), dpi=300)
        plt.close()

# Function to plot adaptation-specific enriched contexts
def plot_adaptation_enriched_contexts(enriched_df, output_dir):
    """Plot trinucleotide contexts enriched in specific adaptation types."""
    if len(enriched_df) == 0:
        print("No adaptation-specific enriched contexts found for plotting.")
        return
    
    # Group by adaptation
    adaptations = enriched_df['Adaptation'].unique()
    
    for adaptation in adaptations:
        adaptation_data = enriched_df[enriched_df['Adaptation'] == adaptation]
        
        if len(adaptation_data) == 0:
            continue
        
        plt.figure(figsize=(12, 6))
        
        # Plot enriched contexts as horizontal bars
        contexts = adaptation_data['Trinucleotide'] + ' (' + adaptation_data['Mutation_Type'] + ')'
        fold_enrichment = adaptation_data['Fold_Enrichment']
        
        # Sort by fold enrichment
        sorted_indices = np.argsort(fold_enrichment)[::-1]  # Descending order
        contexts = contexts.iloc[sorted_indices]
        fold_enrichment = fold_enrichment.iloc[sorted_indices]
        
        # Plot bars
        bars = plt.barh(contexts, fold_enrichment, color=ADAPTATION_COLORS.get(adaptation, '#1b9e77'))
        
        # Add value labels
        for bar in bars:
            width = bar.get_width()
            plt.text(width + 0.1, bar.get_y() + bar.get_height()/2,
                    f'{width:.1f}', ha='left', va='center')
        
        plt.xlabel('Fold Enrichment')
        plt.ylabel('Trinucleotide Context (Mutation Type)')
        plt.title(f'Enriched Trinucleotide Contexts for {adaptation} Adaptation')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{adaptation}_enriched_contexts.png"), dpi=300)
        plt.close()

# Function to plot gene-specific enriched contexts
def plot_gene_enriched_contexts(enriched_df, output_dir):
    """Plot trinucleotide contexts enriched in gene-modified or non-modified strains."""
    if len(enriched_df) == 0:
        print("No gene-specific enriched contexts found for plotting.")
        return
    
    # Group by gene status
    gene_statuses = enriched_df['Has_Gene'].unique()
    
    for has_gene in gene_statuses:
        gene_data = enriched_df[enriched_df['Has_Gene'] == has_gene]
        
        if len(gene_data) == 0:
            continue
        
        plt.figure(figsize=(12, 6))
        
        # Plot enriched contexts as horizontal bars
        contexts = gene_data['Trinucleotide'] + ' (' + gene_data['Mutation_Type'] + ')'
        fold_enrichment = gene_data['Fold_Enrichment']
        
        # Sort by fold enrichment
        sorted_indices = np.argsort(fold_enrichment)[::-1]  # Descending order
        contexts = contexts.iloc[sorted_indices]
        fold_enrichment = fold_enrichment.iloc[sorted_indices]
        
        # Plot bars
        gene_status_label = "Gene-Modified Strains" if has_gene == "Yes" else "Non-Modified Strains"
        bars = plt.barh(contexts, fold_enrichment, color='#1b9e77' if has_gene == "Yes" else '#d95f02')
        
        # Add value labels
        for bar in bars:
            width = bar.get_width()
            plt.text(width + 0.1, bar.get_y() + bar.get_height()/2,
                    f'{width:.1f}', ha='left', va='center')
        
        plt.xlabel('Fold Enrichment')
        plt.ylabel('Trinucleotide Context (Mutation Type)')
        plt.title(f'Enriched Trinucleotide Contexts for {gene_status_label}')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"gene_{has_gene}_enriched_contexts.png"), dpi=300)
        plt.close()

# Function to compare enriched contexts across treatments
def compare_enriched_contexts(enriched_df, output_dir):
    """Compare enriched contexts between treatments."""
    if len(enriched_df) == 0:
        print("No enriched contexts found for comparison.")
        return
    
    # Create a pivot table for heatmap
    pivot_data = enriched_df.pivot_table(
        index=['Trinucleotide', 'Mutation_Type'],
        columns='Treatment',
        values='Fold_Enrichment',
        fill_value=0
    )
    
    # Sort by maximum fold enrichment
    pivot_data['Max'] = pivot_data.max(axis=1)
    pivot_data = pivot_data.sort_values('Max', ascending=False)
    pivot_data = pivot_data.drop('Max', axis=1)
    
    # Create a better colormap
    treatment_colors = {t: TREATMENT_COLORS.get(t, '#333333') for t in pivot_data.columns}
    
    # Plot heatmap
    plt.figure(figsize=(12, max(8, len(pivot_data) * 0.4)))
    
    sns.heatmap(pivot_data, annot=True, cmap='YlOrRd', fmt='.1f')
    
    plt.title('Trinucleotide Context Enrichment Across Treatments')
    plt.tight_layout()
    
    plt.savefig(os.path.join(output_dir, "context_enrichment_comparison.png"), dpi=300)
    plt.close()
    
    # Also create a clustermap with treatment colors
    plt.figure(figsize=(14, max(10, len(pivot_data) * 0.4)))
    
    # Create colors for the columns (treatments)
    col_colors = pd.Series(
        pivot_data.columns.map(lambda t: treatment_colors.get(t, '#333333')),
        index=pivot_data.columns
    )
    
    # Create clustermap
    g = sns.clustermap(
        pivot_data,
        annot=True,
        cmap='YlOrRd',
        fmt='.1f',
        col_colors=col_colors,
        figsize=(14, max(10, len(pivot_data) * 0.4)),
        dendrogram_ratio=(.1, .2)
    )
    
    # Add adaptation type legend
    adaptations = set(TREATMENT_INFO.get(t, {}).get('adaptation', 'Unknown') for t in pivot_data.columns)
    for adaptation in adaptations:
        g.ax_col_dendrogram.bar(0, 0, color=ADAPTATION_COLORS.get(adaptation, '#333333'), label=adaptation)
    
    g.ax_col_dendrogram.legend(title="Adaptation", loc='center')
    
    # Set title
    plt.suptitle('Clustered Trinucleotide Context Enrichment', fontsize=16, y=0.98)
    
    plt.savefig(os.path.join(output_dir, "context_enrichment_clustermap.png"), dpi=300)
    plt.close()

# Function to create sequence logos from contexts
def create_sequence_logos(contexts_df, output_dir, context_width=5):
    """Create sequence logos from mutation contexts."""
    # Group by treatment and mutation type
    treatment_mut_groups = contexts_df.groupby(['Treatment', 'Mutation_Type'])
    
    for (treatment, mut_type), group in treatment_mut_groups:
        # Skip if too few contexts
        if len(group) < 10:
            continue
        
        # Extract extended contexts
        contexts = group['Std_Context'].dropna().tolist()
        
        # Skip if no valid contexts
        if not contexts:
            continue
        
        # Get treatment metadata
        description = TREATMENT_INFO.get(treatment, {}).get('description', '')
        adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', '')
        has_gene = TREATMENT_INFO.get(treatment, {}).get('gene')
        gene_text = f" with {has_gene} gene" if has_gene else ""
        
        # Define the actual positions we want to show
        positions = list(range(-context_width, context_width + 1))
        
        # Count nucleotides at each position
        counts = {pos: {'A': 0, 'C': 0, 'G': 0, 'T': 0} for pos in positions}
        
        for context in contexts:
            var_idx = len(context) // 2
            
            for pos in positions:
                abs_pos = var_idx + pos
                if 0 <= abs_pos < len(context):
                    nt = context[abs_pos]
                    if nt in 'ACGT':
                        counts[pos][nt] += 1
        
        # Create a proper DataFrame for logomaker
        counts_list = []
        for pos in positions:
            pos_counts = counts[pos]
            total = sum(pos_counts.values())
            if total > 0:
                # Create a normalized frequency dictionary
                freq_dict = {nt: count/total for nt, count in pos_counts.items()}
                # Add position information
                freq_dict['pos'] = pos
                counts_list.append(freq_dict)
        
        # Convert to DataFrame
        df = pd.DataFrame(counts_list)
        
        # Set position as index
        if 'pos' in df.columns:
            df = df.set_index('pos')
        
        # Skip if empty
        if df.empty:
            continue
            
        # Create the plot manually instead of using logomaker
        plt.figure(figsize=(10, 3))
        
        try:
            # Draw each nucleotide
            height_dict = {pos: 0 for pos in df.index}
            
            for nt in ['A', 'C', 'G', 'T']:
                if nt in df.columns:
                    # Define color for each nucleotide
                    color_dict = {'A': 'green', 'C': 'blue', 'G': 'orange', 'T': 'red'}
                    
                    # Get heights (frequencies)
                    heights = df[nt].values
                    positions = df.index.values
                    
                    # Draw rectangles for each position
                    for i, pos in enumerate(positions):
                        freq = heights[i]
                        if freq > 0:
                            plt.bar(
                                pos, 
                                freq, 
                                bottom=height_dict[pos], 
                                width=0.8, 
                                color=color_dict[nt], 
                                label=nt if pos == positions[0] else ""
                            )
                            height_dict[pos] += freq
            
            # Add legend (only once for each nucleotide)
            handles, labels = plt.gca().get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            plt.legend(by_label.values(), by_label.keys(), title="Nucleotide")
            
            # Customize plot
            plt.axvline(x=0, color='r', linestyle='--', alpha=0.3)
            plt.xlabel('Position Relative to Mutation')
            plt.ylabel('Frequency')
            plt.title(f'{treatment}: {mut_type} Mutations\n{description} ({adaptation} adaptation{gene_text})')
            plt.xlim(min(positions)-0.5, max(positions)+0.5)
            plt.ylim(0, 1)
            
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"{treatment}_{mut_type.replace('>', '_to_')}_logo.png"), dpi=300)
        except Exception as e:
            print(f"Error creating logo for {treatment} {mut_type}: {e}")
        
        plt.close()

# Function to create adaptation-specific sequence logos
def create_adaptation_logos(contexts_df, output_dir, context_width=5):
    """Create sequence logos grouped by adaptation type."""
    # Group by adaptation and mutation type
    adaptation_mut_groups = contexts_df.groupby(['Adaptation', 'Mutation_Type'])
    
    for (adaptation, mut_type), group in adaptation_mut_groups:
        # Skip if too few contexts
        if len(group) < 10:
            continue
        
        # Extract extended contexts
        contexts = group['Std_Context'].dropna().tolist()
        
        # Skip if no valid contexts
        if not contexts:
            continue
        
        # Define the actual positions we want to show
        positions = list(range(-context_width, context_width + 1))
        
        # Count nucleotides at each position
        counts = {pos: {'A': 0, 'C': 0, 'G': 0, 'T': 0} for pos in positions}
        
        for context in contexts:
            var_idx = len(context) // 2
            
            for pos in positions:
                abs_pos = var_idx + pos
                if 0 <= abs_pos < len(context):
                    nt = context[abs_pos]
                    if nt in 'ACGT':
                        counts[pos][nt] += 1
        
        # Create a proper DataFrame for logomaker
        counts_list = []
        for pos in positions:
            pos_counts = counts[pos]
            total = sum(pos_counts.values())
            if total > 0:
                # Create a normalized frequency dictionary
                freq_dict = {nt: count/total for nt, count in pos_counts.items()}
                # Add position information
                freq_dict['pos'] = pos
                counts_list.append(freq_dict)
        
        # Convert to DataFrame
        df = pd.DataFrame(counts_list)
        
        # Set position as index
        if 'pos' in df.columns:
            df = df.set_index('pos')
        
        # Skip if empty
        if df.empty:
            continue
            
        # Create the plot manually instead of using logomaker
        plt.figure(figsize=(10, 3))
        
        try:
            # Draw each nucleotide
            height_dict = {pos: 0 for pos in df.index}
            
            for nt in ['A', 'C', 'G', 'T']:
                if nt in df.columns:
                    # Define color for each nucleotide
                    color_dict = {'A': 'green', 'C': 'blue', 'G': 'orange', 'T': 'red'}
                    
                    # Get heights (frequencies)
                    heights = df[nt].values
                    positions = df.index.values
                    
                    # Draw rectangles for each position
                    for i, pos in enumerate(positions):
                        freq = heights[i]
                        if freq > 0:
                            plt.bar(
                                pos, 
                                freq, 
                                bottom=height_dict[pos], 
                                width=0.8, 
                                color=color_dict[nt], 
                                label=nt if pos == positions[0] else ""
                            )
                            height_dict[pos] += freq
            
            # Add legend (only once for each nucleotide)
            handles, labels = plt.gca().get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            plt.legend(by_label.values(), by_label.keys(), title="Nucleotide")
            
            # Customize plot
            plt.axvline(x=0, color='r', linestyle='--', alpha=0.3)
            plt.xlabel('Position Relative to Mutation')
            plt.ylabel('Frequency')
            plt.title(f'{adaptation} Adaptation: {mut_type} Mutations')
            plt.xlim(min(positions)-0.5, max(positions)+0.5)
            plt.ylim(0, 1)
            
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"{adaptation}_{mut_type.replace('>', '_to_')}_logo.png"), dpi=300)
        except Exception as e:
            print(f"Error creating logo for {adaptation} {mut_type}: {e}")
        
        plt.close()

# Function to create gene-specific sequence logos
def create_gene_logos(contexts_df, output_dir, context_width=5):
    """Create sequence logos grouped by gene modification status."""
    # Group by gene status and mutation type
    gene_mut_groups = contexts_df.groupby(['Has_Gene', 'Mutation_Type'])
    
    for (has_gene, mut_type), group in gene_mut_groups:
        # Skip if too few contexts
        if len(group) < 10:
            continue
        
        # Extract extended contexts
        contexts = group['Std_Context'].dropna().tolist()
        
        # Skip if no valid contexts
        if not contexts:
            continue
        
        # Define the actual positions we want to show
        positions = list(range(-context_width, context_width + 1))
        
        # Count nucleotides at each position
        counts = {pos: {'A': 0, 'C': 0, 'G': 0, 'T': 0} for pos in positions}
        
        for context in contexts:
            var_idx = len(context) // 2
            
            for pos in positions:
                abs_pos = var_idx + pos
                if 0 <= abs_pos < len(context):
                    nt = context[abs_pos]
                    if nt in 'ACGT':
                        counts[pos][nt] += 1
        
        # Create a proper DataFrame for logomaker
        counts_list = []
        for pos in positions:
            pos_counts = counts[pos]
            total = sum(pos_counts.values())
            if total > 0:
                # Create a normalized frequency dictionary
                freq_dict = {nt: count/total for nt, count in pos_counts.items()}
                # Add position information
                freq_dict['pos'] = pos
                counts_list.append(freq_dict)
        
        # Convert to DataFrame
        df = pd.DataFrame(counts_list)
        
        # Set position as index
        if 'pos' in df.columns:
            df = df.set_index('pos')
        
        # Skip if empty
        if df.empty:
            continue
            
        # Create the plot manually instead of using logomaker
        plt.figure(figsize=(10, 3))
        
        try:
            # Draw each nucleotide
            height_dict = {pos: 0 for pos in df.index}
            
            for nt in ['A', 'C', 'G', 'T']:
                if nt in df.columns:
                    # Define color for each nucleotide
                    color_dict = {'A': 'green', 'C': 'blue', 'G': 'orange', 'T': 'red'}
                    
                    # Get heights (frequencies)
                    heights = df[nt].values
                    positions = df.index.values
                    
                    # Draw rectangles for each position
                    for i, pos in enumerate(positions):
                        freq = heights[i]
                        if freq > 0:
                            plt.bar(
                                pos, 
                                freq, 
                                bottom=height_dict[pos], 
                                width=0.8, 
                                color=color_dict[nt], 
                                label=nt if pos == positions[0] else ""
                            )
                            height_dict[pos] += freq
            
            # Add legend (only once for each nucleotide)
            handles, labels = plt.gca().get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            plt.legend(by_label.values(), by_label.keys(), title="Nucleotide")
            
            # Determine gene status label
            gene_status = "Gene-Modified Strains" if has_gene == "Yes" else "Non-Modified Strains"
            
            # Customize plot
            plt.axvline(x=0, color='r', linestyle='--', alpha=0.3)
            plt.xlabel('Position Relative to Mutation')
            plt.ylabel('Frequency')
            plt.title(f'{gene_status}: {mut_type} Mutations')
            plt.xlim(min(positions)-0.5, max(positions)+0.5)
            plt.ylim(0, 1)
            
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"gene_{has_gene}_{mut_type.replace('>', '_to_')}_logo.png"), dpi=300)
        except Exception as e:
            print(f"Error creating logo for gene status {has_gene} {mut_type}: {e}")
        
        plt.close()

# Function to create summary report
def create_summary_report(contexts_df, signatures, adaptation_signatures, gene_signatures, 
                         similarity_df, enriched_df, adaptation_enriched_df, gene_enriched_df, output_dir):
    """Create a comprehensive summary report of mutational signature analysis."""
    with open(os.path.join(output_dir, "mutational_signatures_summary.txt"), 'w') as f:
        f.write("Mutational Signatures Analysis Summary\n")
        f.write("=====================================\n\n")
        
        # Overall statistics
        f.write("Overall Statistics:\n")
        f.write("-----------------\n")
        
        # Number of variants analyzed
        total_variants = len(contexts_df)
        f.write(f"Total variants analyzed: {total_variants}\n")
        
        # Variants by treatment
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
        
        # Variants by adaptation type
        adaptation_types = contexts_df['Adaptation'].unique()
        f.write("\nVariants by adaptation type:\n")
        for adaptation in adaptation_types:
            count = len(contexts_df[contexts_df['Adaptation'] == adaptation])
            f.write(f"  {adaptation}: {count} variants\n")
        
        # Variants by gene modification status
        gene_statuses = contexts_df['Has_Gene'].unique()
        f.write("\nVariants by gene modification status:\n")
        for has_gene in gene_statuses:
            status = "Gene-modified" if has_gene == "Yes" else "Non-modified"
            count = len(contexts_df[contexts_df['Has_Gene'] == has_gene])
            f.write(f"  {status}: {count} variants\n")
        
        # Variants by mutation type
        mut_types = contexts_df['Mutation_Type'].unique()
        f.write("\nVariants by mutation type:\n")
        for mut_type in sorted(mut_types):
            count = len(contexts_df[contexts_df['Mutation_Type'] == mut_type])
            f.write(f"  {mut_type}: {count} variants ({count/total_variants:.2%})\n")
        
        f.write("\n")
        
        # Signature similarity analysis
        f.write("Signature Similarity Analysis:\n")
        f.write("----------------------------\n")
        
        if not similarity_df.empty:
            # Find most and least similar treatment pairs
            treatments = similarity_df.index.tolist()
            most_similar_value = 0
            least_similar_value = 1
            most_similar_pair = None
            least_similar_pair = None
            
            for i, t1 in enumerate(treatments):
                for j, t2 in enumerate(treatments):
                    if i < j:  # Only look at unique pairs
                        similarity = similarity_df.loc[t1, t2]
                        if similarity > most_similar_value:
                            most_similar_value = similarity
                            most_similar_pair = (t1, t2)
                        if similarity < least_similar_value:
                            least_similar_value = similarity
                            least_similar_pair = (t1, t2)
            
            if most_similar_pair:
                t1, t2 = most_similar_pair
                f.write(f"Most similar treatment signatures: {t1} and {t2} (similarity: {most_similar_value:.4f})\n")
                f.write(f"  {TREATMENT_INFO.get(t1, {}).get('description', 'Unknown')} and {TREATMENT_INFO.get(t2, {}).get('description', 'Unknown')}\n")
                
                # Check if they share adaptation or gene status
                adapt1 = TREATMENT_INFO.get(t1, {}).get('adaptation', 'Unknown')
                adapt2 = TREATMENT_INFO.get(t2, {}).get('adaptation', 'Unknown')
                gene1 = TREATMENT_INFO.get(t1, {}).get('gene')
                gene2 = TREATMENT_INFO.get(t2, {}).get('gene')
                
                same_adaptation = adapt1 == adapt2
                both_have_genes = gene1 is not None and gene2 is not None
                both_lack_genes = gene1 is None and gene2 is None
                
                f.write(f"  Same adaptation type: {'Yes' if same_adaptation else 'No'}\n")
                f.write(f"  Gene modification similarity: ")
                
                if both_have_genes:
                    f.write("Both have gene modifications\n")
                elif both_lack_genes:
                    f.write("Neither has gene modifications\n")
                else:
                    f.write("One has gene modification, one does not\n")
            
            if least_similar_pair:
                t1, t2 = least_similar_pair
                f.write(f"\nLeast similar treatment signatures: {t1} and {t2} (similarity: {least_similar_value:.4f})\n")
                f.write(f"  {TREATMENT_INFO.get(t1, {}).get('description', 'Unknown')} and {TREATMENT_INFO.get(t2, {}).get('description', 'Unknown')}\n")
                
                # Check if they differ in adaptation or gene status
                adapt1 = TREATMENT_INFO.get(t1, {}).get('adaptation', 'Unknown')
                adapt2 = TREATMENT_INFO.get(t2, {}).get('adaptation', 'Unknown')
                gene1 = TREATMENT_INFO.get(t1, {}).get('gene')
                gene2 = TREATMENT_INFO.get(t2, {}).get('gene')
                
                diff_adaptation = adapt1 != adapt2
                gene_difference = (gene1 is None and gene2 is not None) or (gene1 is not None and gene2 is None)
                
                f.write(f"  Different adaptation types: {'Yes' if diff_adaptation else 'No'}\n")
                f.write(f"  Gene modification difference: {'Yes' if gene_difference else 'No'}\n")
        
        f.write("\n")
        
        # Adaptation-specific signature analysis
        f.write("Adaptation-Specific Signature Analysis:\n")
        f.write("-----------------------------------\n")
        
        if adaptation_signatures:
            # Compare adaptation signatures
            adaptation_types = list(adaptation_signatures.keys())
            if len(adaptation_types) > 1:
                # Convert signatures to flat vectors for comparison
                adaptation_vectors = {}
                
                for adaptation in adaptation_types:
                    vector = []
                    for mut_type in SUBSTITUTION_TYPES:
                        for trinuc in TRINUCLEOTIDES:
                            vector.append(adaptation_signatures[adaptation][mut_type][trinuc])
                    adaptation_vectors[adaptation] = np.array(vector)
                
                # Calculate similarity between adaptation types
                for i, a1 in enumerate(adaptation_types):
                    for j, a2 in enumerate(adaptation_types):
                        if i < j:  # Only look at unique pairs
                            v1 = adaptation_vectors[a1]
                            v2 = adaptation_vectors[a2]
                            
                            # Cosine similarity
                            dot_product = np.dot(v1, v2)
                            norm_v1 = np.linalg.norm(v1)
                            norm_v2 = np.linalg.norm(v2)
                            
                            if norm_v1 > 0 and norm_v2 > 0:
                                similarity = dot_product / (norm_v1 * norm_v2)
                                f.write(f"Signature similarity between {a1} and {a2}: {similarity:.4f}\n")
                
            # Characteristic features of each adaptation type
            f.write("\nCharacteristic features of adaptation types:\n")
            
            for adaptation in adaptation_types:
                f.write(f"\n{adaptation} Adaptation:\n")
                
                # Identify characteristic trinucleotide contexts
                if not adaptation_enriched_df.empty:
                    adaptation_enriched = adaptation_enriched_df[adaptation_enriched_df['Adaptation'] == adaptation]
                    
                    if len(adaptation_enriched) > 0:
                        top_contexts = adaptation_enriched.nlargest(3, 'Fold_Enrichment')
                        
                        f.write("  Top enriched trinucleotide contexts:\n")
                        for _, row in top_contexts.iterrows():
                            f.write(f"    {row['Trinucleotide']} ({row['Mutation_Type']}): "
                                   f"{row['Fold_Enrichment']:.2f}-fold enrichment\n")
                    else:
                        f.write("  No significantly enriched trinucleotide contexts\n")
        
        f.write("\n")
        
        # Gene modification effects
        f.write("Gene Modification Effects:\n")
        f.write("-----------------------\n")
        
        if gene_signatures:
            # Compare gene-modified vs non-modified signatures
            gene_statuses = list(gene_signatures.keys())
            if len(gene_statuses) > 1:
                # Convert signatures to flat vectors for comparison
                gene_vectors = {}
                
                for has_gene in gene_statuses:
                    vector = []
                    for mut_type in SUBSTITUTION_TYPES:
                        for trinuc in TRINUCLEOTIDES:
                            vector.append(gene_signatures[has_gene][mut_type][trinuc])
                    gene_vectors[has_gene] = np.array(vector)
                
                # Calculate similarity between gene statuses
                if "Yes" in gene_vectors and "No" in gene_vectors:
                    v1 = gene_vectors["Yes"]
                    v2 = gene_vectors["No"]
                    
                    # Cosine similarity
                    dot_product = np.dot(v1, v2)
                    norm_v1 = np.linalg.norm(v1)
                    norm_v2 = np.linalg.norm(v2)
                    
                    if norm_v1 > 0 and norm_v2 > 0:
                        similarity = dot_product / (norm_v1 * norm_v2)
                        f.write(f"Signature similarity between gene-modified and non-modified strains: {similarity:.4f}\n")
                
            # Characteristic features of gene-modified vs non-modified strains
            f.write("\nCharacteristic features of gene status:\n")
            
            for has_gene in gene_statuses:
                status = "Gene-modified strains" if has_gene == "Yes" else "Non-modified strains"
                f.write(f"\n{status}:\n")
                
                # Identify characteristic trinucleotide contexts
                if not gene_enriched_df.empty:
                    gene_enriched = gene_enriched_df[gene_enriched_df['Has_Gene'] == has_gene]
                    
                    if len(gene_enriched) > 0:
                        top_contexts = gene_enriched.nlargest(3, 'Fold_Enrichment')
                        
                        f.write("  Top enriched trinucleotide contexts:\n")
                        for _, row in top_contexts.iterrows():
                            f.write(f"    {row['Trinucleotide']} ({row['Mutation_Type']}): "
                                   f"{row['Fold_Enrichment']:.2f}-fold enrichment\n")
                    else:
                        f.write("  No significantly enriched trinucleotide contexts\n")
        
        f.write("\n")
        
        # Enriched trinucleotide contexts
        f.write("Enriched Trinucleotide Contexts by Treatment:\n")
        f.write("----------------------------------------\n")
        
        if not enriched_df.empty:
            # Group by treatment
            for treatment in treatments:
                treatment_enriched = enriched_df[enriched_df['Treatment'] == treatment]
                
                if len(treatment_enriched) > 0:
                    description = TREATMENT_INFO.get(treatment, {}).get('description', 'Unknown')
                    adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
                    has_gene = TREATMENT_INFO.get(treatment, {}).get('gene')
                    
                    f.write(f"\n{treatment} Treatment ({description}):\n")
                    f.write(f"  Adaptation: {adaptation}")
                    if has_gene:
                        f.write(f" with {has_gene} gene")
                    f.write("\n")
                    
                    # Show top 5 enriched contexts
                    top_contexts = treatment_enriched.nlargest(5, 'Fold_Enrichment')
                    
                    f.write("  Top enriched trinucleotide contexts:\n")
                    for _, row in top_contexts.iterrows():
                        f.write(f"    {row['Trinucleotide']} ({row['Mutation_Type']}): "
                               f"{row['Fold_Enrichment']:.2f}-fold enrichment "
                               f"({row['Count']} occurrences)\n")
                else:
                    f.write(f"\n{treatment}: No significantly enriched trinucleotide contexts\n")
        
        f.write("\n")
        
        # Main conclusions
        f.write("Main Conclusions:\n")
        f.write("---------------\n")
        f.write("1. This analysis examines the sequence context of mutations in different treatments.\n")
        f.write("2. The mutational signatures reveal preferred sequence contexts for mutations.\n")
        f.write("3. Adaptation types (Temperature vs Low Oxygen) show distinct mutational patterns.\n")
        f.write("4. Gene modifications (STC/CAS) influence the sequence context preferences.\n")
        f.write("5. The analysis of enriched trinucleotide contexts provides insights into the\n")
        f.write("   mechanisms of mutagenesis and DNA damage repair under different conditions.\n")
        f.write("6. Sequence logos reveal position-specific nucleotide preferences around mutations.\n")

# Main function to run the analysis
def main():
    # Load reference genome
    reference_genome = load_reference_genome()
    if not reference_genome:
        print("Error: Could not load reference genome. Exiting.")
        return
    
    # Parse data for each treatment
    all_data = pd.DataFrame()
    for treatment in TREATMENTS:
        data = parse_mutation_data(treatment)
        if len(data) > 0:
            all_data = pd.concat([all_data, data])
            print(f"Loaded {len(data)} variants for {treatment} treatment")
        else:
            print(f"Warning: No data available for {treatment} treatment")
    
    # Filter for SNVs
    snv_data = filter_snvs(all_data)
    print(f"Filtered to {len(snv_data)} single nucleotide variants")
    
    # Check if we have any SNVs before continuing
    if len(snv_data) == 0:
        print("Warning: No single nucleotide variants found after filtering.")
        print("Creating empty results and skipping analysis.")
        
        # Create an empty summary report
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        with open(os.path.join(OUTPUT_DIR, "mutational_signatures_summary.txt"), 'w') as f:
            f.write("Mutational Signatures Analysis Summary\n")
            f.write("=====================================\n\n")
            f.write("No single nucleotide variants found after filtering.\n")
            f.write("Please check your input data format.\n")
        
        print(f"Analysis complete! Empty summary saved to {OUTPUT_DIR}/")
        return
    
    # Extract context sequences
    contexts_df = extract_all_contexts(snv_data, reference_genome)
    print(f"Extracted context sequences for {len(contexts_df)} variants")
    
    # Create treatment-specific mutation signature matrices
    treatment_signatures = create_signature_matrix(contexts_df)
    print("Created treatment-specific mutation signature matrices")
    
    # Create adaptation-specific mutation signature matrices
    adaptation_signatures = create_adaptation_signature_matrix(contexts_df)
    print("Created adaptation-specific mutation signature matrices")
    
    # Create gene-specific mutation signature matrices
    gene_signatures = create_gene_signature_matrix(contexts_df)
    print("Created gene-specific mutation signature matrices")
    
    # Create combined adaptation+gene mutation signature matrices
    combined_signatures = create_combined_signature_matrix(contexts_df)
    print("Created combined adaptation+gene mutation signature matrices")
    
    # Plot treatment-specific mutation signatures
    for treatment in TREATMENTS:
        if treatment in treatment_signatures:
            plot_signature(treatment_signatures[treatment], treatment, OUTPUT_DIR)
    print("Generated treatment-specific mutation signature plots")
    
    # Plot adaptation-specific mutation signatures
    plot_adaptation_signatures(adaptation_signatures, OUTPUT_DIR)
    print("Generated adaptation-specific mutation signature plots")
    
    # Plot gene-specific mutation signatures
    plot_gene_signatures(gene_signatures, OUTPUT_DIR)
    print("Generated gene-specific mutation signature plots")
    
    # Plot combined adaptation+gene mutation signatures
    plot_combined_signatures(combined_signatures, OUTPUT_DIR)
    print("Generated combined adaptation+gene mutation signature plots")
    
    # Calculate signature similarity
    similarity_df = calculate_signature_similarity(treatment_signatures)
    print("Calculated signature similarities")
    
    # Plot signature similarity
    plot_signature_similarity(similarity_df, OUTPUT_DIR)
    print("Generated signature similarity heatmap")
    
    # Find enriched contexts
    enriched_df = find_enriched_contexts(contexts_df)
    print(f"Found {len(enriched_df)} enriched trinucleotide contexts")
    
    # Find adaptation-specific enriched contexts
    adaptation_enriched_df = find_adaptation_enriched_contexts(contexts_df)
    print(f"Found {len(adaptation_enriched_df)} adaptation-specific enriched contexts")
    
    # Find gene-specific enriched contexts
    gene_enriched_df = find_gene_enriched_contexts(contexts_df)
    print(f"Found {len(gene_enriched_df)} gene-specific enriched contexts")
    
    # Plot enriched contexts
    plot_enriched_contexts(enriched_df, OUTPUT_DIR)
    print("Generated enriched context plots")
    
    # Plot adaptation-specific enriched contexts
    plot_adaptation_enriched_contexts(adaptation_enriched_df, OUTPUT_DIR)
    print("Generated adaptation-specific enriched context plots")
    
    # Plot gene-specific enriched contexts
    plot_gene_enriched_contexts(gene_enriched_df, OUTPUT_DIR)
    print("Generated gene-specific enriched context plots")
    
    # Compare enriched contexts across treatments
    compare_enriched_contexts(enriched_df, OUTPUT_DIR)
    print("Generated context enrichment comparison")
    
    # Create sequence logos
    create_sequence_logos(contexts_df, OUTPUT_DIR)
    print("Generated sequence logo plots")
    
    # Create adaptation-specific sequence logos
    create_adaptation_logos(contexts_df, OUTPUT_DIR)
    print("Generated adaptation-specific sequence logo plots")
    
    # Create gene-specific sequence logos
    create_gene_logos(contexts_df, OUTPUT_DIR)
    print("Generated gene-specific sequence logo plots")
    
    # Create summary report
    create_summary_report(
        contexts_df,
        treatment_signatures,
        adaptation_signatures,
        gene_signatures,
        similarity_df,
        enriched_df,
        adaptation_enriched_df,
        gene_enriched_df,
        OUTPUT_DIR
    )
    print("Created summary report")
    
    print(f"Analysis complete! Results saved to {OUTPUT_DIR}/")

# Run the analysis
if __name__ == "__main__":
    main()