#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict, Counter
from Bio import SeqIO
import subprocess
import logomaker
from scipy.stats import chi2_contingency, fisher_exact
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
OUTPUT_DIR = "mutational_signatures_results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Define nucleotide complements
COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

# Define the six base substitution types (pyrimidine-centric)
SUBSTITUTION_TYPES = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']

# Define all possible trinucleotides
NUCLEOTIDES = ['A', 'C', 'G', 'T']
TRINUCLEOTIDES = [''.join(x) for x in 
                  [(a, b, c) for a in NUCLEOTIDES for b in NUCLEOTIDES for c in NUCLEOTIDES]]

# Reference genome path
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

# Function to plot mutation signature
def plot_signature(signature, treatment, output_dir):
    """Plot a mutation signature in COSMIC style."""
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
    plt.title(f'Mutation Signature for {treatment} Treatment')
    
    # Set x-ticks
    plt.xticks(ticks, ticklabels, rotation=90, fontsize=8)
    
    # Add legend
    plt.legend(title='Substitution Type')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(os.path.join(output_dir, f"{treatment}_mutation_signature.png"), dpi=300)
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
    plt.figure(figsize=(10, 8))
    
    # Create heatmap
    sns.heatmap(similarity_df, annot=True, cmap='viridis', vmin=0, vmax=1)
    
    plt.title('Mutation Signature Similarity Between Treatments')
    plt.tight_layout()
    
    plt.savefig(os.path.join(output_dir, "signature_similarity_heatmap.png"), dpi=300)
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
                enriched_contexts.append({
                    'Treatment': treatment,
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
        
        plt.figure(figsize=(12, 6))
        
        # Plot enriched contexts as horizontal bars
        contexts = treatment_data['Trinucleotide'] + ' (' + treatment_data['Mutation_Type'] + ')'
        fold_enrichment = treatment_data['Fold_Enrichment']
        
        # Sort by fold enrichment
        sorted_indices = np.argsort(fold_enrichment)
        contexts = contexts.iloc[sorted_indices]
        fold_enrichment = fold_enrichment.iloc[sorted_indices]
        
        # Plot bars
        bars = plt.barh(contexts, fold_enrichment, color='skyblue')
        
        # Add value labels
        for bar in bars:
            width = bar.get_width()
            plt.text(width + 0.1, bar.get_y() + bar.get_height()/2,
                    f'{width:.1f}', ha='left', va='center')
        
        plt.xlabel('Fold Enrichment')
        plt.ylabel('Trinucleotide Context (Mutation Type)')
        plt.title(f'Enriched Trinucleotide Contexts for {treatment} Treatment')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{treatment}_enriched_contexts.png"), dpi=300)
        plt.close()

# Function to compare enriched contexts across treatments
def compare_enriched_contexts(enriched_df, output_dir):
    """Compare enriched contexts across treatments."""
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
    
    # Plot heatmap
    plt.figure(figsize=(12, max(8, len(pivot_data) * 0.4)))
    
    sns.heatmap(pivot_data, annot=True, cmap='YlOrRd', fmt='.1f')
    
    plt.title('Trinucleotide Context Enrichment Across Treatments')
    plt.tight_layout()
    
    plt.savefig(os.path.join(output_dir, "context_enrichment_comparison.png"), dpi=300)
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
            plt.title(f'{treatment} Treatment: {mut_type} Mutations')
            plt.xlim(min(positions)-0.5, max(positions)+0.5)
            plt.ylim(0, 1)
            
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"{treatment}_{mut_type.replace('>', '_to_')}_logo.png"), dpi=300)
        except Exception as e:
            print(f"Error creating logo for {treatment} {mut_type}: {e}")
        
        plt.close()

# Function to create summary report
def create_summary_report(contexts_df, signatures, similarity_df, enriched_df, output_dir):
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
            count = len(contexts_df[contexts_df['Treatment'] == treatment])
            f.write(f"  {treatment}: {count} variants\n")
        
        # Variants by mutation type
        mut_types = contexts_df['Mutation_Type'].unique()
        f.write("\nVariants by mutation type:\n")
        for mut_type in sorted(mut_types):
            count = len(contexts_df[contexts_df['Mutation_Type'] == mut_type])
            f.write(f"  {mut_type}: {count} variants ({count/total_variants:.2%})\n")
        
        f.write("\nSignature Similarity Analysis:\n")
        f.write("----------------------------\n")
        
        # Find most and least similar treatment pairs
        max_similarity = 0
        min_similarity = 1
        max_pair = None
        min_pair = None
        
        for i, t1 in enumerate(treatments):
            for j, t2 in enumerate(treatments):
                if i < j:  # Only look at unique pairs
                    similarity = similarity_df.loc[t1, t2]
                    
                    if similarity > max_similarity:
                        max_similarity = similarity
                        max_pair = (t1, t2)
                    
                    if similarity < min_similarity:
                        min_similarity = similarity
                        min_pair = (t1, t2)
        
        if max_pair:
            f.write(f"Most similar treatment signatures: {max_pair[0]} and {max_pair[1]} (similarity: {max_similarity:.4f})\n")
        
        if min_pair:
            f.write(f"Least similar treatment signatures: {min_pair[0]} and {min_pair[1]} (similarity: {min_similarity:.4f})\n")
        
        f.write("\n")
        
        # Enriched contexts
        f.write("Enriched Trinucleotide Contexts:\n")
        f.write("-----------------------------\n")
        
        if len(enriched_df) > 0:
            # Group by treatment
            for treatment in treatments:
                treatment_enriched = enriched_df[enriched_df['Treatment'] == treatment]
                
                if len(treatment_enriched) > 0:
                    f.write(f"\n{treatment} Treatment:\n")
                    
                    # Sort by fold enrichment
                    sorted_enriched = treatment_enriched.sort_values('Fold_Enrichment', ascending=False)
                    
                    for i, row in sorted_enriched.head(5).iterrows():
                        f.write(f"  {row['Trinucleotide']} ({row['Mutation_Type']}): "
                                f"{row['Fold_Enrichment']:.2f}-fold enrichment "
                                f"({row['Count']} occurrences)\n")
        else:
            f.write("No significantly enriched trinucleotide contexts found.\n")
        
        f.write("\n")
        
        # Main conclusions
        f.write("Main Conclusions:\n")
        f.write("---------------\n")
        f.write("1. This analysis examines the sequence context of mutations in different treatments.\n")
        f.write("2. Treatment-specific mutational signatures reveal preferred sequence contexts for mutations.\n")
        f.write("3. The signature similarity analysis shows relationships between mutational mechanisms.\n")
        f.write("4. Enriched trinucleotide contexts provide insights into potential damage mechanisms.\n")

# Main function to run the analysis
def main():
    # Define treatments
    treatments = ['WT', 'STC', 'CAS', 'WTA']
    
    # Load reference genome
    reference_genome = load_reference_genome()
    if not reference_genome:
        print("Error: Could not load reference genome. Exiting.")
        return
    
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
    contexts_df = extract_all_contexts(snv_data, reference_genome)
    print(f"Extracted context sequences for {len(contexts_df)} variants")
    
    # Create mutation signature matrix
    signatures = create_signature_matrix(contexts_df)
    print("Created mutation signature matrices")
    
    # Plot mutation signatures
    for treatment in treatments:
        if treatment in signatures:
            plot_signature(signatures[treatment], treatment, OUTPUT_DIR)
    print("Generated mutation signature plots")
    
    # Calculate signature similarity
    similarity_df = calculate_signature_similarity(signatures)
    print("Calculated signature similarities")
    
    # Plot signature similarity
    plot_signature_similarity(similarity_df, OUTPUT_DIR)
    print("Generated signature similarity heatmap")
    
    # Find enriched contexts
    enriched_df = find_enriched_contexts(contexts_df)
    print(f"Found {len(enriched_df)} enriched trinucleotide contexts")
    
    # Plot enriched contexts
    plot_enriched_contexts(enriched_df, OUTPUT_DIR)
    print("Generated enriched context plots")
    
    # Compare enriched contexts across treatments
    compare_enriched_contexts(enriched_df, OUTPUT_DIR)
    print("Generated context enrichment comparison")
    
    # Create sequence logos
    create_sequence_logos(contexts_df, OUTPUT_DIR)
    print("Generated sequence logo plots")
    
    # Create summary report
    create_summary_report(contexts_df, signatures, similarity_df, enriched_df, OUTPUT_DIR)
    print("Created summary report")
    
    print(f"Analysis complete! Results saved to {OUTPUT_DIR}/")

# Run the analysis
if __name__ == "__main__":
    main()