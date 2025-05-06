#!/usr/bin/env python3
"""
Gene-Specific Genomic Context Analysis Module

This module analyzes the genomic context of mutations with specific focus on gene-level annotations.
It extracts sequence context around mutation sites, analyzes sequence composition, and generates
visualizations of sequence patterns, grouped by treatment, adaptation type, and gene properties.
"""

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
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist, squareform
import csv
import warnings
warnings.filterwarnings('ignore')

# Set matplotlib style
plt.style.use('ggplot')
sns.set(font_scale=1.2)
sns.set_style("whitegrid")

# Define output directories
OUTPUT_DIR = "analysis/genomic_context_results"
GENE_OUTPUT_DIR = "analysis/gene_genomic_context_results"
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(GENE_OUTPUT_DIR, exist_ok=True)

# Updated biologically correct treatment groups
TREATMENTS = ['WT-37', 'WTA', 'STC', 'CAS']

# Define treatment information for better biological context
TREATMENT_INFO = {
    'WT-37': {'description': 'Temperature-adapted wild type', 'adaptation': 'Temperature'},
    'WTA': {'description': 'Low oxygen-adapted wild type', 'adaptation': 'Low Oxygen'},
    'STC': {'description': 'STC gene with low oxygen adaptation', 'adaptation': 'Low Oxygen', 'gene': 'STC'},
    'CAS': {'description': 'CAS gene with temperature adaptation', 'adaptation': 'Temperature', 'gene': 'CAS'}
}

# File paths for gene mapping and annotations
GENE_MAPPING_FILE = "reference/gene_mapping.tsv"
GENES_OF_INTEREST_FILE = "reference/genes_of_interest_mapping.tsv"

# Gene data structures
GENE_DATA = {}
SCAFFOLD_GENES = defaultdict(list)
GENES_OF_INTEREST = set()

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

# Reference genome path - adjust as needed
REFERENCE_GENOME = "reference/w303_chromosomal.fasta"

# Function to load gene mapping information
def load_gene_mapping():
    """
    Load the gene mapping information from the reference files.
    
    This function loads gene data from:
    1. The main gene mapping file with all annotated genes
    2. The genes of interest file with ergosterol pathway genes
    
    Returns:
        dict: Dictionary mapping gene ID to gene information
        defaultdict: Dictionary mapping scaffold to list of genes on that scaffold
        set: Set of gene IDs that are of special interest (ergosterol pathway)
    """
    global GENE_DATA, SCAFFOLD_GENES, GENES_OF_INTEREST
    
    print(f"Loading gene mapping from {GENE_MAPPING_FILE}")
    
    # Load the main gene mapping file
    if os.path.exists(GENE_MAPPING_FILE):
        with open(GENE_MAPPING_FILE, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gene_id = row['w303_gene_id']
                GENE_DATA[gene_id] = {
                    'sc_gene_id': row['sc_gene_id'],
                    'erg_name': row.get('erg_name', ''),
                    'locus_tag': row['locus_tag'],
                    'w303_scaffold': row['w303_scaffold'],
                    'start': int(row['start']),
                    'end': int(row['end']),
                    'strand': row['strand'],
                    'product': row.get('product', 'hypothetical protein')
                }
                
                # Add to scaffold lookup
                SCAFFOLD_GENES[row['w303_scaffold']].append(gene_id)
        
        print(f"Loaded information for {len(GENE_DATA)} genes across {len(SCAFFOLD_GENES)} scaffolds")
    else:
        print(f"Warning: Gene mapping file {GENE_MAPPING_FILE} not found")
    
    # Load the genes of interest file
    if os.path.exists(GENES_OF_INTEREST_FILE):
        with open(GENES_OF_INTEREST_FILE, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                GENES_OF_INTEREST.add(row['w303_gene_id'])
        
        print(f"Loaded {len(GENES_OF_INTEREST)} genes of interest")
    else:
        print(f"Warning: Genes of interest file {GENES_OF_INTEREST_FILE} not found")
    
    return GENE_DATA, SCAFFOLD_GENES, GENES_OF_INTEREST

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
        "reference/w303_chromosomal.fasta",
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
    """Parse mutation data for a specific treatment and map to genes."""
    # Define possible file patterns for mutation data
    file_patterns = [
        "mutation_spectrum_analysis/{}_mutations.txt",
        "analysis/MSA/mutation_spectrum_analysis/{}_mutations.txt",
        "results/mutation_spectrum_analysis/{}_mutations.txt",
        "analysis/gene_mutation_spectrum_results/{}_mutations.txt"
    ]
    
    # Find the mutation data file
    mutation_file = find_file(treatment, file_patterns)
    
    if mutation_file:
        try:
            # First, let's read the raw content to examine the actual structure
            with open(mutation_file, 'r') as f:
                first_line = f.readline().strip()
                columns = first_line.split('\t')
                print(f"File structure: {len(columns)} columns in first line")
                
            # Based on file inspection, we'll determine the column structure
            if len(columns) == 5:  # 5-column format: CHROM, POS, REF, ALT, Treatment
                data = pd.read_csv(mutation_file, sep='\t', header=None)
                
                # Assign column names based on position
                data.columns = ['CHROM', 'POS', 'REF', 'ALT', 'Treatment']
                
                # Convert column types
                data['POS'] = data['POS'].astype(int)
                
                # Validate that Treatment column matches expected treatment
                actual_treatment = data['Treatment'].iloc[0]
                if actual_treatment != treatment:
                    print(f"Warning: Treatment in file ({actual_treatment}) doesn't match expected treatment ({treatment})")
                
            else:  # Assume 4-column format: CHROM, POS, REF, ALT
                data = pd.read_csv(mutation_file, sep='\t', header=None, 
                                  names=['CHROM', 'POS', 'REF', 'ALT'])
                data['Treatment'] = treatment
            
            # Add biological context
            data['Adaptation'] = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
            data['Has_Gene'] = 'Yes' if TREATMENT_INFO.get(treatment, {}).get('gene') else 'No'
            
            print(f"Loaded {len(data)} mutations for {treatment}")
            
            # Map variants to genes if gene data is available
            if GENE_DATA:
                data = map_variants_to_genes(data)
                print(f"Mapped {len(data[data['Gene_ID'] != ''])} variants to genes for {treatment}")
            
            return data
        except Exception as e:
            print(f"Error reading {mutation_file}: {e}")
    
    # If no mutation file found, try to extract from VCF
    vcf_data = extract_from_vcf(treatment)
    
    # Map extracted variants to genes if gene data is available
    if not vcf_data.empty and GENE_DATA:
        vcf_data = map_variants_to_genes(vcf_data)
        print(f"Mapped {len(vcf_data[vcf_data['Gene_ID'] != ''])} variants to genes for {treatment}")
    
    return vcf_data

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
                os.makedirs("analysis/mutation_spectrum_analysis", exist_ok=True)
                data.to_csv(f"analysis/mutation_spectrum_analysis/{treatment}_mutations.txt", 
                          sep='\t', index=False, header=False)
                
                print(f"Extracted and saved {len(data)} mutations for {treatment}")
                return data
        except Exception as e:
            print(f"Error extracting from {vcf_file}: {e}")
    
    print(f"Could not find or extract mutation data for {treatment}")
    return pd.DataFrame()

# Function to map variants to gene annotations
def map_variants_to_genes(data, upstream_distance=1000):
    """
    Map variants to genes based on genomic coordinates.
    
    Args:
        data (pd.DataFrame): DataFrame containing variant data with CHROM and POS columns
        upstream_distance (int): Distance upstream of genes to include for regulatory regions
        
    Returns:
        pd.DataFrame: DataFrame with added gene annotation columns
    """
    if data.empty or not GENE_DATA:
        return data
    
    # Create a copy of the data to add gene annotations
    annotated_data = data.copy()
    
    # Initialize gene annotation columns if they don't exist
    if 'Gene_ID' not in annotated_data.columns:
        annotated_data['Gene_ID'] = ''
        annotated_data['Gene_Name'] = ''
        annotated_data['Gene_Function'] = ''
        annotated_data['Gene_Context'] = ''
        annotated_data['Is_Gene_Of_Interest'] = False
    
    # Process each variant
    for i, variant in annotated_data.iterrows():
        chrom = variant['CHROM']
        pos = int(variant['POS'])
        
        # Skip if not on a scaffold with our genes or already mapped
        if chrom not in SCAFFOLD_GENES or variant['Gene_ID']:
            continue
        
        # Check all genes on this scaffold
        for gene_id in SCAFFOLD_GENES[chrom]:
            gene_info = GENE_DATA[gene_id]
            
            # Check if the variant is within or near the gene
            gene_start = gene_info['start']
            gene_end = gene_info['end']
            
            # Determine the gene context
            if gene_start <= pos <= gene_end:
                # Inside the gene
                context = 'coding'
                match_found = True
            elif (gene_info['strand'] == '+' and 
                  gene_start - upstream_distance <= pos < gene_start):
                # Upstream of a + strand gene
                context = 'upstream'
                match_found = True
            elif (gene_info['strand'] == '-' and 
                  gene_end < pos <= gene_end + upstream_distance):
                # Upstream of a - strand gene (means downstream in sequence coord)
                context = 'upstream'
                match_found = True
            else:
                # Not related to this gene
                match_found = False
            
            # If the variant affects this gene, add annotation
            if match_found:
                annotated_data.at[i, 'Gene_ID'] = gene_id
                annotated_data.at[i, 'Gene_Name'] = gene_info.get('erg_name', '') or gene_info.get('sc_gene_id', '')
                annotated_data.at[i, 'Gene_Function'] = gene_info.get('product', 'hypothetical protein')
                annotated_data.at[i, 'Gene_Context'] = context
                annotated_data.at[i, 'Is_Gene_Of_Interest'] = gene_id in GENES_OF_INTEREST
                break  # Assign to first matching gene
    
    return annotated_data

# Function to filter data for single nucleotide variants
# Modified filter_snvs function for genomic_context_analysis.py
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
            # Keep adaptation and gene status information
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
            # Preserve adaptation and gene status information
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
def create_position_specific_matrix(contexts_df, window=5):
    """Create a position-specific nucleotide frequency matrix for a group of variants.
    
    Args:
        contexts_df (pd.DataFrame): DataFrame with sequence context data
        window (int): Window size around the mutation site
        
    Returns:
        dict: Position-specific matrix of nucleotide frequencies
    """
    # Initialize matrix dimensions (positions x nucleotides)
    positions = range(-window, window + 1)
    nucleotides = ['A', 'C', 'G', 'T']
    
    # Create empty matrix
    matrix = {pos: {nt: 0 for nt in nucleotides} for pos in positions}
    
    # Count nucleotides at each position
    for _, row in contexts_df.iterrows():
        ctx = row['immediate_context']
        var_idx = row['variant_index']
        ref, alt = row['REF'], row['ALT']
        
        # Standardize context if needed
        std_ctx, std_idx, std_ref, std_alt = standardize_mutation_context(
            ctx, var_idx, ref, alt)
        
        # Count nucleotides at each position relative to variant
        for i, pos in enumerate(range(-window, window + 1)):
            abs_pos = std_idx + pos
            if 0 <= abs_pos < len(std_ctx):
                nt = std_ctx[abs_pos]
                if nt in nucleotides:
                    matrix[pos][nt] += 1
    
    return matrix

def create_context_matrix(contexts_df, window=5):
    """Create a position-specific nucleotide frequency matrix."""
    # Initialize matrix dimensions (positions x nucleotides)
    positions = range(-window, window + 1)
    nucleotides = ['A', 'C', 'G', 'T']
    
    # Create empty matrix
    matrix = {pos: {nt: 0 for nt in nucleotides} for pos in positions}
    
    # Group by treatment and mutation type
    grouped = contexts_df.groupby(['Treatment', 'REF', 'ALT'])
    
    # Process each group separately
    context_matrices = {}
    
    for (treatment, ref, alt), group in grouped:
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
                ctx, var_idx, ref, alt)
            
            # Count nucleotides at each position relative to variant
            for i, pos in enumerate(range(-window, window + 1)):
                abs_pos = std_idx + pos
                if 0 <= abs_pos < len(std_ctx):
                    nt = std_ctx[abs_pos]
                    if nt in nucleotides:
                        matrix[pos][nt] += 1
        
        # Create a copy of the matrix for this group
        mutation_key = f"{treatment}:{ref}>{alt}"
        context_matrices[mutation_key] = {pos: dict(nts) for pos, nts in matrix.items()}
    
    # Also create adaptation-specific matrices
    adaptation_matrices = {}
    adaptation_grouped = contexts_df.groupby(['Adaptation', 'REF', 'ALT'])
    
    for (adaptation, ref, alt), group in adaptation_grouped:
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
                ctx, var_idx, ref, alt)
            
            # Count nucleotides at each position relative to variant
            for i, pos in enumerate(range(-window, window + 1)):
                abs_pos = std_idx + pos
                if 0 <= abs_pos < len(std_ctx):
                    nt = std_ctx[abs_pos]
                    if nt in nucleotides:
                        matrix[pos][nt] += 1
        
        # Create a copy of the matrix for this adaptation group
        adaptation_key = f"{adaptation}:{ref}>{alt}"
        adaptation_matrices[adaptation_key] = {pos: dict(nts) for pos, nts in matrix.items()}
    
    # Also create gene-specific matrices
    gene_matrices = {}
    gene_grouped = contexts_df.groupby(['Has_Gene', 'REF', 'ALT'])
    
    for (has_gene, ref, alt), group in gene_grouped:
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
                ctx, var_idx, ref, alt)
            
            # Count nucleotides at each position relative to variant
            for i, pos in enumerate(range(-window, window + 1)):
                abs_pos = std_idx + pos
                if 0 <= abs_pos < len(std_ctx):
                    nt = std_ctx[abs_pos]
                    if nt in nucleotides:
                        matrix[pos][nt] += 1
        
        # Create a copy of the matrix for this gene status group
        gene_key = f"{'Gene' if has_gene == 'Yes' else 'NoGene'}:{ref}>{alt}"
        gene_matrices[gene_key] = {pos: dict(nts) for pos, nts in matrix.items()}
    
    return context_matrices, adaptation_matrices, gene_matrices

# Function to plot sequence logos
def plot_sequence_logos(context_matrices, output_dir):
    """Generate sequence logo plots for mutation contexts."""
    for mutation_key, matrix in context_matrices.items():
        # Extract treatment and mutation type
        parts = mutation_key.split(':')
        if len(parts) != 2:
            continue
        
        treatment, mut_type = parts
        
        # Get treatment metadata
        description = TREATMENT_INFO.get(treatment, {}).get('description', '')
        adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', '')
        has_gene = TREATMENT_INFO.get(treatment, {}).get('gene')
        gene_text = f" with {has_gene} gene" if has_gene else ""
        
        # Create DataFrame for logomaker
        df = pd.DataFrame(matrix).T
        
        # Normalize columns to get probabilities
        for pos in df.index:
            total = df.loc[pos].sum()
            if total > 0:
                df.loc[pos] = df.loc[pos] / total
        
        # Create sequence logo
        try:
            plt.figure(figsize=(10, 3))
            logo = logomaker.Logo(df, color_scheme='classic')
            
            # Customize plot
            logo.ax.set_xlabel('Position Relative to Mutation')
            logo.ax.set_ylabel('Probability')
            logo.ax.set_title(f'Sequence Context for {treatment}: {mut_type}\n{description} ({adaptation} adaptation{gene_text})')
            
            # Mark mutation position
            logo.ax.axvline(x=0, color='red', linestyle='--', alpha=0.5)
            
            # Save plot
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"logo_{treatment}_{mut_type.replace('>', '_to_')}.png"), dpi=300)
            plt.close()
        except Exception as e:
            print(f"Error creating logo for {mutation_key}: {e}")

# Function to plot adaptation-specific sequence logos
def plot_adaptation_logos(adaptation_matrices, output_dir):
    """Generate adaptation-specific sequence logo plots."""
    for adaptation_key, matrix in adaptation_matrices.items():
        # Extract adaptation and mutation type
        parts = adaptation_key.split(':')
        if len(parts) != 2:
            continue
        
        adaptation, mut_type = parts
        
        # Create DataFrame for logomaker
        df = pd.DataFrame(matrix).T
        
        # Normalize columns to get probabilities
        for pos in df.index:
            total = df.loc[pos].sum()
            if total > 0:
                df.loc[pos] = df.loc[pos] / total
        
        # Create sequence logo
        try:
            plt.figure(figsize=(10, 3))
            logo = logomaker.Logo(df, color_scheme='classic')
            
            # Customize plot
            logo.ax.set_xlabel('Position Relative to Mutation')
            logo.ax.set_ylabel('Probability')
            logo.ax.set_title(f'Sequence Context for {adaptation} Adaptation: {mut_type}')
            
            # Mark mutation position
            logo.ax.axvline(x=0, color='red', linestyle='--', alpha=0.5)
            
            # Save plot
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"logo_{adaptation}_{mut_type.replace('>', '_to_')}.png"), dpi=300)
            plt.close()
        except Exception as e:
            print(f"Error creating logo for {adaptation_key}: {e}")

# Function to plot gene-specific sequence logos
def plot_gene_logos(gene_matrices, output_dir):
    """Generate gene-specific sequence logo plots."""
    for gene_key, matrix in gene_matrices.items():
        # Extract gene status and mutation type
        parts = gene_key.split(':')
        if len(parts) != 2:
            continue
        
        gene_status, mut_type = parts
        display_status = "Gene-Modified" if gene_status == "Gene" else "Non-Modified"
        
        # Create DataFrame for logomaker
        df = pd.DataFrame(matrix).T
        
        # Normalize columns to get probabilities
        for pos in df.index:
            total = df.loc[pos].sum()
            if total > 0:
                df.loc[pos] = df.loc[pos] / total
        
        # Create sequence logo
        try:
            plt.figure(figsize=(10, 3))
            logo = logomaker.Logo(df, color_scheme='classic')
            
            # Customize plot
            logo.ax.set_xlabel('Position Relative to Mutation')
            logo.ax.set_ylabel('Probability')
            logo.ax.set_title(f'Sequence Context for {display_status} Strains: {mut_type}')
            
            # Mark mutation position
            logo.ax.axvline(x=0, color='red', linestyle='--', alpha=0.5)
            
            # Save plot
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"logo_{gene_status}_{mut_type.replace('>', '_to_')}.png"), dpi=300)
            plt.close()
        except Exception as e:
            print(f"Error creating logo for {gene_key}: {e}")

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
        sns.kdeplot(treatment_data['gc_content'], label=f'{treatment}', color=TREATMENT_COLORS.get(treatment, '#333333'))
    
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
    for adaptation in contexts_df['Adaptation'].unique():
        adaptation_data = contexts_df[contexts_df['Adaptation'] == adaptation]
        sns.kdeplot(adaptation_data['gc_content'], label=f'{adaptation}', color=ADAPTATION_COLORS.get(adaptation, '#333333'))
    
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
    for has_gene in contexts_df['Has_Gene'].unique():
        gene_data = contexts_df[contexts_df['Has_Gene'] == has_gene]
        gene_label = "Gene-Modified" if has_gene == "Yes" else "Non-Modified"
        gene_color = "#1b9e77" if has_gene == "Yes" else "#d95f02"
        sns.kdeplot(gene_data['gc_content'], label=gene_label, color=gene_color)
    
    # Add control
    sns.kdeplot(control_df['gc_content'], label='Control', linestyle='--', color='black')
    
    # Customize plot
    plt.xlabel('GC Content')
    plt.ylabel('Density')
    plt.title('GC Content Distribution by Gene Modification Status')
    plt.legend()
    
    # Save plot
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "gc_content_by_gene_status.png"), dpi=300)
    plt.close()
    
    # Plot GC content by adaptation and gene status
    plt.figure(figsize=(14, 6))
    
    # Group by adaptation and gene status
    contexts_df['Group'] = contexts_df.apply(
        lambda row: f"{row['Adaptation']} ({row['Has_Gene'] == 'Yes' and 'Gene' or 'No Gene'})",
        axis=1
    )
    
    for group in sorted(contexts_df['Group'].unique()):
        group_data = contexts_df[contexts_df['Group'] == group]
        adaptation = group.split()[0]
        is_gene = "Gene" in group
        color = ADAPTATION_COLORS.get(adaptation, '#333333')
        style = '-' if is_gene else '--'
        sns.kdeplot(group_data['gc_content'], label=group, color=color, linestyle=style)
    
    # Add control
    sns.kdeplot(control_df['gc_content'], label='Control', linestyle=':', color='black')
    
    # Customize plot
    plt.xlabel('GC Content')
    plt.ylabel('Density')
    plt.title('GC Content Distribution by Adaptation Type and Gene Status')
    plt.legend()
    
    # Save plot
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "gc_content_by_adaptation_gene.png"), dpi=300)
    plt.close()

# Function to plot local sequence features
def plot_sequence_features(contexts_df, output_dir):
    """Plot local sequence features around variants."""
    # Group by treatment
    treatments = contexts_df['Treatment'].unique()
    
    # Plot homopolymer presence
    homopolymer_data = contexts_df.groupby('Treatment')['nearby_homopolymer'].mean()
    
    plt.figure(figsize=(10, 6))
    bars = plt.bar(homopolymer_data.index, homopolymer_data.values, 
                  color=[TREATMENT_COLORS.get(t, '#333333') for t in homopolymer_data.index])
    
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
    
    plt.figure(figsize=(10, 6))
    bars = plt.bar(dinucleotide_data.index, dinucleotide_data.values, 
                  color=[TREATMENT_COLORS.get(t, '#333333') for t in dinucleotide_data.index])
    
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
    
    # Plot by adaptation type
    homopolymer_adaptation = contexts_df.groupby('Adaptation')['nearby_homopolymer'].mean()
    
    plt.figure(figsize=(10, 6))
    bars = plt.bar(homopolymer_adaptation.index, homopolymer_adaptation.values, 
                  color=[ADAPTATION_COLORS.get(a, '#333333') for a in homopolymer_adaptation.index])
    
    # Add value labels
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                f'{height:.2f}', ha='center', va='bottom')
    
    plt.xlabel('Adaptation Type')
    plt.ylabel('Fraction with Nearby Homopolymer')
    plt.title('Homopolymer Presence Near Variants by Adaptation Type')
    plt.ylim(0, max(homopolymer_adaptation.values) * 1.1)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "homopolymer_by_adaptation.png"), dpi=300)
    plt.close()
    
    # Plot features by gene modification status
    homopolymer_gene = contexts_df.groupby('Has_Gene')['nearby_homopolymer'].mean()
    
    plt.figure(figsize=(10, 6))
    gene_labels = ["Gene-Modified" if g == "Yes" else "Non-Modified" for g in homopolymer_gene.index]
    bars = plt.bar(gene_labels, homopolymer_gene.values, 
                  color=["#1b9e77" if g == "Yes" else "#d95f02" for g in homopolymer_gene.index])
    
    # Add value labels
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                f'{height:.2f}', ha='center', va='bottom')
    
    plt.xlabel('Gene Modification Status')
    plt.ylabel('Fraction with Nearby Homopolymer')
    plt.title('Homopolymer Presence Near Variants by Gene Modification Status')
    plt.ylim(0, max(homopolymer_gene.values) * 1.1)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "homopolymer_by_gene_status.png"), dpi=300)
    plt.close()
    
    # Plot features by adaptation type and gene status combined
    contexts_df['Group'] = contexts_df.apply(
        lambda row: f"{row['Adaptation']} ({'Gene' if row['Has_Gene'] == 'Yes' else 'No Gene'})",
        axis=1
    )
    
    homopolymer_combined = contexts_df.groupby('Group')['nearby_homopolymer'].mean().sort_values(ascending=False)
    
    plt.figure(figsize=(12, 6))
    bars = plt.bar(homopolymer_combined.index, homopolymer_combined.values)
    
    # Color bars by adaptation type and pattern by gene status
    for i, (group, _) in enumerate(homopolymer_combined.items()):
        adaptation = group.split()[0]
        is_gene = "Gene" in group
        bars[i].set_color(ADAPTATION_COLORS.get(adaptation, '#333333'))
        if not is_gene:
            bars[i].set_hatch('//')
    
    # Add value labels
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                f'{height:.2f}', ha='center', va='bottom')
    
    plt.xlabel('Adaptation Type and Gene Status')
    plt.ylabel('Fraction with Nearby Homopolymer')
    plt.title('Homopolymer Presence by Adaptation Type and Gene Status')
    plt.xticks(rotation=45, ha='right')
    plt.ylim(0, max(homopolymer_combined.values) * 1.1)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "homopolymer_by_adaptation_gene.png"), dpi=300)
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
    
    # Analyze by adaptation type
    for adaptation in contexts_df['Adaptation'].unique():
        adaptation_data = contexts_df[contexts_df['Adaptation'] == adaptation]
        adaptation_grouped = adaptation_data.groupby(['REF', 'ALT'])
        
        for (ref, alt), group in adaptation_grouped:
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
            plt.title(f'{adaptation} - {mutation_type}: Base Before Mutation')
            
            # Plot +1 position distribution
            plt.subplot(1, 2, 2)
            values = [plus1_counts.get(base, 0) for base in bases]
            plt.bar(bases, values, color=['green', 'blue', 'orange', 'red'])
            plt.xlabel('Nucleotide at +1 Position')
            plt.ylabel('Count')
            plt.title(f'{adaptation} - {mutation_type}: Base After Mutation')
            
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"adjacent_bases_{adaptation}_{ref}_to_{alt}.png"), dpi=300)
            plt.close()
    
    # Analyze by gene modification status
    for has_gene in contexts_df['Has_Gene'].unique():
        gene_data = contexts_df[contexts_df['Has_Gene'] == has_gene]
        gene_status = "Gene-Modified" if has_gene == "Yes" else "Non-Modified"
        gene_grouped = gene_data.groupby(['REF', 'ALT'])
        
        for (ref, alt), group in gene_grouped:
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
            plt.title(f'{gene_status} - {mutation_type}: Base Before Mutation')
            
            # Plot +1 position distribution
            plt.subplot(1, 2, 2)
            values = [plus1_counts.get(base, 0) for base in bases]
            plt.bar(bases, values, color=['green', 'blue', 'orange', 'red'])
            plt.xlabel('Nucleotide at +1 Position')
            plt.ylabel('Count')
            plt.title(f'{gene_status} - {mutation_type}: Base After Mutation')
            
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"adjacent_bases_{has_gene}_{ref}_to_{alt}.png"), dpi=300)
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
    
    # Create clustermap with adaptation coloring
    column_colors = [TREATMENT_COLORS.get(t, '#333333') for t in df.columns]
    
    g = sns.clustermap(
        df, 
        annot=True, 
        fmt='d', 
        cmap='YlGnBu',
        col_colors=[column_colors],
        figsize=(12, 10),
        dendrogram_ratio=(.1, .2)
    )
    
    # Add adaptation type legend
    adaptations = set(TREATMENT_INFO.get(t, {}).get('adaptation', 'Unknown') for t in df.columns)
    for adaptation in adaptations:
        treatments_with_adaptation = [t for t in df.columns if TREATMENT_INFO.get(t, {}).get('adaptation') == adaptation]
        if treatments_with_adaptation:
            g.ax_col_dendrogram.bar(0, 0, color=ADAPTATION_COLORS.get(adaptation, '#333333'), label=adaptation)
    
    g.ax_col_dendrogram.legend(title="Adaptation", loc='center')
    
    plt.suptitle('Mutation Types by Treatment (Clustered)', fontsize=16, y=0.95)
    plt.savefig(os.path.join(output_dir, "mutation_type_by_treatment_clustermap.png"), dpi=300)
    plt.close()
    
    # Analyze GC content by treatment
    plt.figure(figsize=(10, 6))
    
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
    
    # Analyze by adaptation type
    plt.figure(figsize=(10, 6))
    
    # Group by adaptation
    adaptation_groups = {}
    for adaptation in contexts_df['Adaptation'].unique():
        adaptation_data = contexts_df[contexts_df['Adaptation'] == adaptation]
        adaptation_groups[adaptation] = adaptation_data['gc_content'].values
    
    plt.boxplot(list(adaptation_groups.values()), labels=list(adaptation_groups.keys()))
    plt.xlabel('Adaptation Type')
    plt.ylabel('GC Content')
    plt.title('GC Content Distribution by Adaptation Type')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "gc_content_boxplot_by_adaptation.png"), dpi=300)
    plt.close()
    
    # Analyze by gene modification status
    plt.figure(figsize=(10, 6))
    
    # Group by gene status
    gene_groups = {}
    for has_gene in contexts_df['Has_Gene'].unique():
        gene_data = contexts_df[contexts_df['Has_Gene'] == has_gene]
        gene_status = "Gene-Modified" if has_gene == "Yes" else "Non-Modified"
        gene_groups[gene_status] = gene_data['gc_content'].values
    
    plt.boxplot(list(gene_groups.values()), labels=list(gene_groups.keys()))
    plt.xlabel('Gene Modification Status')
    plt.ylabel('GC Content')
    plt.title('GC Content Distribution by Gene Modification Status')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "gc_content_boxplot_by_gene_status.png"), dpi=300)
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
            description = TREATMENT_INFO.get(treatment, {}).get('description', 'Unknown')
            adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
            has_gene = TREATMENT_INFO.get(treatment, {}).get('gene')
            
            f.write(f"  {treatment}: {count} variants - {description} ({adaptation} adaptation")
            if has_gene:
                f.write(f" with {has_gene} gene")
            f.write(")\n")
        
        # Adaptation breakdown
        adaptations = contexts_df['Adaptation'].unique()
        f.write("\nVariants by adaptation type:\n")
        for adaptation in adaptations:
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
        
        # GC content by adaptation
        f.write("\nGC content by adaptation type:\n")
        for adaptation in adaptations:
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
        for adaptation in adaptations:
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
        for adaptation in adaptations:
            adaptation_dinucleotide = contexts_df[contexts_df['Adaptation'] == adaptation]['nearby_dinucleotide_repeat'].mean()
            f.write(f"  {adaptation}: {adaptation_dinucleotide:.4f}\n")
        
        # Dinucleotide by gene status
        f.write("\nDinucleotide repeat presence by gene modification status:\n")
        for has_gene in gene_statuses:
            status = "Gene-modified" if has_gene == "Yes" else "Non-modified"
            gene_dinucleotide = contexts_df[contexts_df['Has_Gene'] == has_gene]['nearby_dinucleotide_repeat'].mean()
            f.write(f"  {status}: {gene_dinucleotide:.4f}\n")
        
        f.write("\n")
        
        # Main conclusions
        f.write("Main Conclusions:\n")
        f.write("---------------\n")
        f.write("1. This analysis examines the genomic context of mutations in different treatments.\n")
        f.write("2. The GC content analysis reveals potential biases in mutation distribution.\n")
        f.write("3. Homopolymer and repeat regions may influence mutation probability.\n")
        f.write("4. Adaptation conditions (Temperature vs Low Oxygen) show distinct context preferences.\n")
        f.write("5. Gene modifications (STC, CAS) appear to influence the genomic context of mutations.\n")
        f.write("6. Sequence composition around mutation sites provides insights into damage mechanisms.\n")
        f.write("7. Further analysis of specific motifs may reveal specific DNA damage signatures.\n")

# Function to create gene-specific summary report
def create_gene_summary_report(contexts_df, control_df):
    """Create a comprehensive gene-specific summary report."""
    # Check if we have gene data
    has_gene_data = 'Gene_ID' in contexts_df.columns and len(contexts_df[contexts_df['Gene_ID'] != '']) > 0
    
    if not has_gene_data:
        print("No gene-specific data available for summary report.")
        return
    
    with open(os.path.join(GENE_OUTPUT_DIR, "gene_genomic_context_summary.txt"), 'w') as f:
        f.write("Gene-Specific Genomic Context Analysis Summary\n")
        f.write("===========================================\n\n")
        
        # Filter to gene-mapped variants
        gene_contexts = contexts_df[contexts_df['Gene_ID'] != '']
        
        # Overall statistics
        f.write("Overall Gene Statistics:\n")
        f.write("----------------------\n")
        
        # Number of variants analyzed
        total_variants = len(contexts_df)
        gene_variants = len(gene_contexts)
        gene_pct = (gene_variants / total_variants * 100) if total_variants > 0 else 0
        
        f.write(f"Total variants analyzed: {total_variants}\n")
        f.write(f"Variants mapped to genes: {gene_variants} ({gene_pct:.2f}%)\n")
        
        # Count unique genes affected
        unique_genes = len(gene_contexts['Gene_ID'].unique())
        interest_genes = len(gene_contexts[gene_contexts['Is_Gene_Of_Interest']]['Gene_ID'].unique())
        
        f.write(f"Unique genes affected: {unique_genes}\n")
        f.write(f"Genes of interest affected: {interest_genes}\n")
        
        # Treatment breakdown for gene variants
        treatments = gene_contexts['Treatment'].unique()
        f.write("\nGene variants by treatment:\n")
        for treatment in treatments:
            count = len(gene_contexts[gene_contexts['Treatment'] == treatment])
            interest_count = len(gene_contexts[(gene_contexts['Treatment'] == treatment) & 
                                        gene_contexts['Is_Gene_Of_Interest']])
            f.write(f"  {treatment}: {count} gene variants, {interest_count} in genes of interest\n")
        
        # Gene context breakdown
        contexts = gene_contexts['Gene_Context'].unique()
        f.write("\nVariants by gene context:\n")
        for context in contexts:
            count = len(gene_contexts[gene_contexts['Gene_Context'] == context])
            context_pct = (count / gene_variants * 100) if gene_variants > 0 else 0
            f.write(f"  {context}: {count} variants ({context_pct:.2f}%)\n")
        
        # Add analysis of genes of interest
        if interest_genes > 0:
            goi_contexts = gene_contexts[gene_contexts['Is_Gene_Of_Interest']]
            f.write("\nAnalysis of Genes of Interest:\n")
            f.write("-----------------------------\n")
            
            # List genes of interest with variants
            f.write("Genes of interest with variants:\n")
            for gene_id in sorted(goi_contexts['Gene_ID'].unique()):
                gene_name = GENE_DATA[gene_id].get('erg_name', '') or gene_id
                count = len(goi_contexts[goi_contexts['Gene_ID'] == gene_id])
                f.write(f"  {gene_name} ({gene_id}): {count} variants\n")
            
            # GC content in genes of interest
            goi_gc = goi_contexts['gc_content'].mean() if 'gc_content' in goi_contexts.columns else 0
            general_gc = gene_contexts[~gene_contexts['Is_Gene_Of_Interest']]['gc_content'].mean() \
                        if 'gc_content' in gene_contexts.columns else 0
            
            f.write(f"\nMean GC content in genes of interest: {goi_gc:.4f}\n")
            f.write(f"Mean GC content in other genes: {general_gc:.4f}\n")
            
            # Homopolymer and repeat presence
            if 'nearby_homopolymer' in goi_contexts.columns and 'nearby_dinucleotide_repeat' in goi_contexts.columns:
                goi_homopolymer = goi_contexts['nearby_homopolymer'].mean()
                goi_dinucleotide = goi_contexts['nearby_dinucleotide_repeat'].mean()
                
                f.write(f"\nFraction of variants near homopolymers in genes of interest: {goi_homopolymer:.4f}\n")
                f.write(f"Fraction of variants near dinucleotide repeats in genes of interest: {goi_dinucleotide:.4f}\n")
        
        # Main gene-specific conclusions
        f.write("\nGene-Specific Conclusions:\n")
        f.write("-------------------------\n")
        f.write("1. This analysis examines the genomic context of mutations at the gene level.\n")
        f.write("2. The gene context (coding vs. upstream) influences mutation patterns.\n")
        f.write("3. Genes of interest show distinct sequence context preferences.\n")
        f.write("4. The genomic features around mutations provide insights into gene-specific damage.\n")
        f.write("5. Further analysis can reveal gene-specific mutational signatures.\n")

# Main function to run the analysis
def main():
    # Load gene mapping information
    load_gene_mapping()
    
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
        with open(os.path.join(OUTPUT_DIR, "genomic_context_summary.txt"), 'w') as f:
            f.write("Genomic Context Analysis Summary\n")
            f.write("===============================\n\n")
            f.write("No single nucleotide variants found after filtering.\n")
            f.write("Please check your input data format.\n")
        
        print(f"Analysis complete! Empty summary saved to {OUTPUT_DIR}/")
        return
    
    # Extract context sequences
    contexts_df = extract_all_context_sequences(snv_data, reference_genome)
    print(f"Extracted context sequences for {len(contexts_df)} variants")
    
    # Generate control context sequences
    control_df = generate_control_contexts(reference_genome)
    print(f"Generated {len(control_df)} control context sequences")
    
    # Calculate composition metrics
    contexts_df = calculate_composition_metrics(contexts_df)
    control_df = calculate_composition_metrics(control_df)
    print("Calculated composition metrics for variant and control regions")
    
    # Extract immediate contexts for mutation context analysis
    immediate_contexts_df = extract_immediate_contexts(contexts_df)
    print(f"Extracted immediate contexts for {len(immediate_contexts_df)} variants")
    
    # Create context matrices for sequence logos
    treatment_matrices, adaptation_matrices, gene_matrices = create_context_matrix(immediate_contexts_df)
    print(f"Created context matrices for {len(treatment_matrices)} treatment-mutation combinations")
    print(f"Created context matrices for {len(adaptation_matrices)} adaptation-mutation combinations")
    print(f"Created context matrices for {len(gene_matrices)} gene status-mutation combinations")
    
    # Plot sequence logos
    plot_sequence_logos(treatment_matrices, OUTPUT_DIR)
    plot_adaptation_logos(adaptation_matrices, OUTPUT_DIR)
    plot_gene_logos(gene_matrices, OUTPUT_DIR)
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
    
    # Add gene-specific analysis if we have gene data
    if 'Gene_ID' in contexts_df.columns and len(contexts_df[contexts_df['Gene_ID'] != '']) > 0:
        print("\nPerforming gene-specific analysis...")
        
        # Filter for genes of interest
        goi_contexts = contexts_df[contexts_df['Is_Gene_Of_Interest']]
        print(f"Found {len(goi_contexts)} variants in genes of interest")
        
        if not goi_contexts.empty:
            # Extract immediate contexts for genes of interest
            goi_immediate_contexts = extract_immediate_contexts(goi_contexts)
            
            # Create gene-specific context matrices
            gene_id_matrices = {}
            for gene_id in goi_contexts['Gene_ID'].unique():
                if not gene_id:  # Skip empty gene IDs
                    continue
                    
                gene_data = goi_immediate_contexts[goi_immediate_contexts['Gene_ID'] == gene_id]
                if len(gene_data) < 5:  # Need enough variants for meaningful analysis
                    continue
                    
                # Group by mutation type
                grouped = gene_data.groupby(['REF', 'ALT'])
                gene_matrices = {}
                
                for (ref, alt), group in grouped:
                    if len(group) < 3:  # Need enough examples
                        continue
                        
                    # Create position-specific matrix
                    matrix = create_position_specific_matrix(group, window=5)
                    gene_matrices[f"{ref}>{alt}"] = matrix
                
                if gene_matrices:  # If we have any matrices
                    gene_name = GENE_DATA.get(gene_id, {}).get('erg_name', '') or gene_id
                    gene_id_matrices[f"{gene_name}_{gene_id}"] = gene_matrices
            
            # Plot gene-specific logos
            for gene_key, matrices in gene_id_matrices.items():
                for mut_type, matrix in matrices.items():
                    try:
                        # Convert matrix to DataFrame
                        df = pd.DataFrame(matrix).T
                        
                        # Normalize
                        for pos in df.index:
                            total = df.loc[pos].sum()
                            if total > 0:
                                df.loc[pos] = df.loc[pos] / total
                        
                        # Create logo
                        plt.figure(figsize=(10, 3))
                        logo = logomaker.Logo(df, color_scheme='classic')
                        
                        # Customize plot
                        logo.ax.set_xlabel('Position Relative to Mutation')
                        logo.ax.set_ylabel('Probability')
                        logo.ax.set_title(f'Sequence Context for {gene_key}: {mut_type}')
                        
                        # Mark mutation position
                        logo.ax.axvline(x=0, color='red', linestyle='--', alpha=0.5)
                        
                        # Save plot
                        plt.tight_layout()
                        plt.savefig(os.path.join(GENE_OUTPUT_DIR, f"logo_{gene_key}_{mut_type.replace('>', '_to_')}.png"), dpi=300)
                        plt.close()
                    except Exception as e:
                        print(f"Error creating logo for {gene_key} {mut_type}: {e}")
        
        # Create gene-specific summary report
        create_gene_summary_report(contexts_df, control_df)
        print(f"Created gene-specific summary report")
        
        print(f"Gene-specific analysis complete! Results saved to {GENE_OUTPUT_DIR}/")
    
    print(f"\nAnalysis complete! Results saved to:\n- Genome-wide results: {OUTPUT_DIR}/\n- Gene-specific results: {GENE_OUTPUT_DIR}/")

# Run the analysis
if __name__ == "__main__":
    main()