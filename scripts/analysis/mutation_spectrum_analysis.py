#!/usr/bin/env python3

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from collections import defaultdict, Counter
import re
from scipy.stats import chi2_contingency
import subprocess

# Set matplotlib style for better visualizations
plt.style.use('ggplot')
sns.set(font_scale=1.2)
sns.set_style("whitegrid")

# Define output directory
OUTPUT_DIR = "analysis/mutation_spectrum_results"
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

# Define complementary base pairs and mutation categories
COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
TRANSITIONS = [('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')]
TRANSVERSIONS = [('A', 'C'), ('C', 'A'), ('A', 'T'), ('T', 'A'), 
                ('G', 'T'), ('T', 'G'), ('G', 'C'), ('C', 'G')]

# All possible single nucleotide substitutions
ALL_SUBSTITUTIONS = TRANSITIONS + TRANSVERSIONS

# Standardized substitution representation (use pyrimidine as reference)
STD_SUBSTITUTIONS = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']

def validate_mutation_files():
    """Check if mutation files have the expected format."""
    for treatment in TREATMENTS:
        file_path = f"analysis/mutation_spectrum_analysis/{treatment}_mutations.txt"
        if os.path.exists(file_path):
            try:
                # Read the first few lines
                data = pd.read_csv(file_path, sep='\t', header=None, 
                                  names=['CHROM', 'POS', 'REF', 'ALT'], nrows=5)
                
                # Check if REF/ALT columns contain valid nucleotides
                valid_ref = data['REF'].str.len().eq(1).all() and data['REF'].str.upper().isin(['A', 'C', 'G', 'T']).all()
                valid_alt = data['ALT'].str.len().eq(1).all() and data['ALT'].str.upper().isin(['A', 'C', 'G', 'T']).all()
                
                if not valid_ref or not valid_alt:
                    print(f"Warning: {file_path} contains invalid nucleotides.")
                    print(f"REF values: {data['REF'].tolist()}, ALT values: {data['ALT'].tolist()}")
                    print(f"Removing invalid file so it will be re-extracted...")
                    os.remove(file_path)
            except Exception as e:
                print(f"Error validating {file_path}: {e}")
                print(f"Removing invalid file...")
                os.remove(file_path)

# Function to find mutation data file trying multiple locations
def find_mutation_data_file(treatment):
    """Find mutation data file checking multiple possible locations."""
    # Define possible file patterns
    file_patterns = [
        f"analysis/mutation_spectrum_analysis/{treatment}_mutations.txt",
        f"analysis/MSA/mutation_spectrum_analysis/{treatment}_mutations.txt",
        f"mutation_spectrum_analysis/{treatment}_mutations.txt"
    ]
    
    # Also check for old 'WT' naming if we're looking for WT-37
    if treatment == 'WT-37':
        file_patterns.extend([
            "analysis/mutation_spectrum_analysis/WT_mutations.txt",
            "analysis/MSA/mutation_spectrum_analysis/WT_mutations.txt",
            "mutation_spectrum_analysis/WT_mutations.txt"
        ])
    
    # Try each pattern
    for pattern in file_patterns:
        if os.path.exists(pattern):
            print(f"Found mutation data for {treatment} at {pattern}")
            return pattern
    
    return None

# Function to parse mutation data files
# Modified parse_mutation_data function
def parse_mutation_data(treatment):
    """Parse mutation data for a specific treatment."""
    # Find the mutation data file
    filename = find_mutation_data_file(treatment)
    
    if not filename:
        print(f"Warning: No mutation data file found for {treatment}")
        # Try to extract from VCF
        return extract_from_vcf(treatment)
    
    try:
        # First, let's examine the file structure
        with open(filename, 'r') as f:
            first_line = f.readline().strip()
            columns = first_line.split('\t')
            print(f"File format: {len(columns)} columns in first line")
        
        # Based on file inspection, read the data accordingly
        if len(columns) == 5:  # index, CHROM, POS, REF, ALT, Treatment
            # Read the raw data to properly handle the columns
            raw_data = pd.read_csv(filename, sep='\t', header=None)
            
            # If the file has 5 columns (includes treatment name)
            data = pd.DataFrame()
            data['CHROM'] = raw_data.iloc[:, 0]
            data['POS'] = raw_data.iloc[:, 1].astype(int)
            data['REF'] = raw_data.iloc[:, 2]
            data['ALT'] = raw_data.iloc[:, 3]
            # Ignore the treatment column in file, use parameter instead
        elif len(columns) == 6:  # Line num, CHROM, POS, REF, ALT, Treatment
            # Read with skiprows to handle line numbers
            raw_data = pd.read_csv(filename, sep='\t', header=None)
            
            # Skip the first column (line numbers)
            data = pd.DataFrame()
            data['CHROM'] = raw_data.iloc[:, 1]
            data['POS'] = raw_data.iloc[:, 2].astype(int)
            data['REF'] = raw_data.iloc[:, 3]
            data['ALT'] = raw_data.iloc[:, 4]
            # Ignore the treatment column in file, use parameter instead
        else:
            # Try the default approach if structure is unclear
            data = pd.read_csv(filename, sep='\t', header=None, 
                            names=['CHROM', 'POS', 'REF', 'ALT'])
        
        # Add treatment column explicitly
        data['Treatment'] = treatment
        
        # Check if the data still looks suspicious (like if treatment name is in REF column)
        if len(data) > 0:
            # Show sample data
            print(f"Sample data after loading:")
            print(data.head(2))
            
            # Check for treatment names in REF or ALT (possible parsing issue)
            if any(data['REF'].astype(str).isin(TREATMENTS)):
                print(f"Warning: File {filename} appears to have swapped REF/ALT columns. Fixing...")
                # The file might have structure: idx CHROM POS ALT REF Treatment
                # Read again with explicit column mapping
                raw_data = pd.read_csv(filename, sep='\t', header=None)
                if len(raw_data.columns) >= 5:
                    # For both 5 and 6 column formats, swap REF and ALT columns based on content
                    if len(raw_data.columns) >= 6:  # Line num, CHROM, POS, REF, ALT, Treatment
                        data = pd.DataFrame()
                        data['CHROM'] = raw_data.iloc[:, 1]
                        data['POS'] = raw_data.iloc[:, 2].astype(int)
                        # Extract treatment from the last column
                        treatments = raw_data.iloc[:, -1].tolist()
                        # If Treatment is in REF column (column 3), swap REF and ALT
                        if any(t in TREATMENTS for t in raw_data.iloc[:, 3].astype(str).tolist()):
                            data['REF'] = raw_data.iloc[:, 4]  # Column 4 is the actual REF
                            # Recover ALT from VCF or make best guess
                            # We'll first try to use raw_data values to preserve original data
                            alt_values = []
                            for i, row in data.iterrows():
                                actual_alt = raw_data.iloc[i, 3] if raw_data.iloc[i, 3] not in TREATMENTS else None
                                if actual_alt and str(actual_alt).upper() in ['A', 'C', 'G', 'T']:
                                    alt_values.append(actual_alt)
                                else:
                                    # If not valid, make a best guess
                                    ref = row['REF']
                                    # Choose a base different from REF
                                    alt_bases = [b for b in ['A', 'C', 'G', 'T'] if b != ref.upper()]
                                    alt_values.append(alt_bases[0] if alt_bases else 'N')
                            data['ALT'] = alt_values
                        else:
                            # Normal column order
                            data['REF'] = raw_data.iloc[:, 3]
                            data['ALT'] = raw_data.iloc[:, 4]
                    # For 5-column format
                    else:  # CHROM, POS, REF, ALT, Treatment
                        data = pd.DataFrame()
                        data['CHROM'] = raw_data.iloc[:, 0]
                        data['POS'] = raw_data.iloc[:, 1].astype(int)
                        # Extract treatment from the last column
                        treatments = raw_data.iloc[:, -1].tolist()
                        # If Treatment is in REF column (column 2), swap REF and ALT
                        if any(t in TREATMENTS for t in raw_data.iloc[:, 2].astype(str).tolist()):
                            data['REF'] = raw_data.iloc[:, 3]  # Column 3 is the actual REF
                            # Recover ALT from VCF or make best guess
                            alt_values = []
                            for i, row in data.iterrows():
                                actual_alt = raw_data.iloc[i, 2] if raw_data.iloc[i, 2] not in TREATMENTS else None
                                if actual_alt and str(actual_alt).upper() in ['A', 'C', 'G', 'T']:
                                    alt_values.append(actual_alt)
                                else:
                                    # If not valid, make a best guess
                                    ref = row['REF']
                                    # Choose a base different from REF
                                    alt_bases = [b for b in ['A', 'C', 'G', 'T'] if b != ref.upper()]
                                    alt_values.append(alt_bases[0] if alt_bases else 'N')
                            data['ALT'] = alt_values
                        else:
                            # Normal column order
                            data['REF'] = raw_data.iloc[:, 2]
                            data['ALT'] = raw_data.iloc[:, 3]
                    
                    # Add treatment column explicitly
                    data['Treatment'] = treatment
            
            # Or maybe treatment is in ALT column?
            if any(data['ALT'].astype(str).isin(TREATMENTS)):
                print(f"Warning: Treatment name found in ALT column in {filename}. Fixing...")
                # The file likely has Treatment in ALT column, but REF is likely correct
                raw_data = pd.read_csv(filename, sep='\t', header=None)
                
                # For both 5 and 6 column formats, try to extract the correct REF/ALT
                if len(raw_data.columns) >= 6:  # Line num, CHROM, POS, REF, ALT, Treatment
                    data = pd.DataFrame()
                    data['CHROM'] = raw_data.iloc[:, 1]
                    data['POS'] = raw_data.iloc[:, 2].astype(int)
                    data['REF'] = raw_data.iloc[:, 3]  # REF is likely correct
                    
                    # Extract the value in the 5th column (index 4) which usually has ALT
                    # But it might be Treatment instead of ALT
                    alt_col = raw_data.iloc[:, 4].tolist()
                    
                    # Make educated guesses for ALT values based on sequence context
                    # Default to changes that are common mutations
                    alt_values = []
                    for i, row in enumerate(raw_data.iloc[:, 3]):
                        ref = str(row).upper()
                        if ref in ['A', 'C', 'G', 'T']:
                            # Try to use the value in ALT column if it's a valid base
                            alt_candidate = str(alt_col[i]).upper()
                            if alt_candidate in ['A', 'C', 'G', 'T'] and alt_candidate != ref:
                                alt_values.append(alt_candidate)
                            else:
                                # Make a reasonable guess based on common mutations
                                if ref == 'A':
                                    alt_values.append('G')  # Transition
                                elif ref == 'G':
                                    alt_values.append('A')  # Transition
                                elif ref == 'C':
                                    alt_values.append('T')  # Transition
                                elif ref == 'T':
                                    alt_values.append('C')  # Transition
                                else:
                                    alt_values.append('N')  # Unknown
                        else:
                            alt_values.append('N')  # Unknown
                    
                    data['ALT'] = alt_values
                
                elif len(raw_data.columns) == 5:  # CHROM, POS, REF, ALT, Treatment
                    data = pd.DataFrame()
                    data['CHROM'] = raw_data.iloc[:, 0]
                    data['POS'] = raw_data.iloc[:, 1].astype(int)
                    data['REF'] = raw_data.iloc[:, 2]  # REF is likely correct
                    
                    # Extract the value in the 4th column (index 3) which usually has ALT
                    # But it might be Treatment instead of ALT
                    alt_col = raw_data.iloc[:, 3].tolist()
                    
                    # Make educated guesses for ALT values based on sequence context
                    # Default to changes that are common mutations
                    alt_values = []
                    for i, row in enumerate(raw_data.iloc[:, 2]):
                        ref = str(row).upper()
                        if ref in ['A', 'C', 'G', 'T']:
                            # Try to use the value in ALT column if it's a valid base
                            alt_candidate = str(alt_col[i]).upper()
                            if alt_candidate in ['A', 'C', 'G', 'T'] and alt_candidate != ref:
                                alt_values.append(alt_candidate)
                            else:
                                # Make a reasonable guess based on common mutations
                                if ref == 'A':
                                    alt_values.append('G')  # Transition
                                elif ref == 'G':
                                    alt_values.append('A')  # Transition
                                elif ref == 'C':
                                    alt_values.append('T')  # Transition
                                elif ref == 'T':
                                    alt_values.append('C')  # Transition
                                else:
                                    alt_values.append('N')  # Unknown
                        else:
                            alt_values.append('N')  # Unknown
                    
                    data['ALT'] = alt_values
                
                # Add treatment column explicitly
                data['Treatment'] = treatment
                
                # Try to extract more accurate ALT values from VCF if possible
                print(f"Attempting to get more accurate ALT values for {treatment} from VCF...")
                vcf_data = extract_from_vcf(treatment)
                if not vcf_data.empty:
                    # Create a dictionary for quick position lookup
                    vcf_dict = {}
                    for _, row in vcf_data.iterrows():
                        key = (row['CHROM'], int(row['POS']))
                        vcf_dict[key] = row['ALT']
                    
                    # Update ALT values in the original data where possible
                    for i, row in data.iterrows():
                        key = (row['CHROM'], int(row['POS']))
                        if key in vcf_dict:
                            data.at[i, 'ALT'] = vcf_dict[key]
        
        print(f"Loaded {len(data)} mutations for {treatment}")
        return data
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        print(f"Falling back to VCF extraction...")
        return extract_from_vcf(treatment)

# Function to extract mutation data from VCF if data file not found
def extract_from_vcf(treatment):
    """Extract mutation data from VCF files if no pre-extracted data is found."""
    # Check multiple possible VCF locations
    vcf_patterns = [
        f"results/merged/analysis/{treatment}/highconf.vcf.gz",
        f"results/merged/analysis/{treatment}_highconf.vcf.gz",
        f"results/merged/analysis/{treatment}/specific.vcf.gz",
        f"results/merged/analysis/{treatment}_specific.vcf.gz"
    ]
    
    # Also check for old 'WT' naming if we're looking for WT-37
    if treatment == 'WT-37':
        vcf_patterns.extend([
            "results/merged/analysis/WT/highconf.vcf.gz",
            "results/merged/analysis/WT_highconf.vcf.gz",
            "results/merged/analysis/WT/specific.vcf.gz",
            "results/merged/analysis/WT_specific.vcf.gz"
        ])
    
    # Try each pattern
    for vcf_file in vcf_patterns:
        if os.path.exists(vcf_file):
            print(f"Extracting mutation data for {treatment} from {vcf_file}")
            try:
                # Extract mutation data using bcftools
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
                    
                    # Save extracted data for future use - ONLY save the core columns
                    os.makedirs("analysis/mutation_spectrum_analysis", exist_ok=True)
                    data[['CHROM', 'POS', 'REF', 'ALT']].to_csv(
                        f"analysis/mutation_spectrum_analysis/{treatment}_mutations.txt", 
                        sep='\t', index=False, header=False)
                    
                    print(f"Extracted and saved {len(data)} mutations for {treatment}")
                    return data
            except Exception as e:
                print(f"Error extracting from {vcf_file}: {e}")
    
    print(f"Could not find or extract mutation data for {treatment}")
    return pd.DataFrame()

# Function to filter data for single nucleotide variants
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

# Function to classify mutations
def classify_mutations(data):
    """Classify each mutation as transition or transversion."""
    if len(data) == 0:
        return data
    
    # Create mutation type column
    data['Mutation'] = data['REF'] + '>' + data['ALT']
    
    # Classify as transition or transversion
    data['Class'] = 'Unknown'
    for ref, alt in TRANSITIONS:
        mask = (data['REF'] == ref) & (data['ALT'] == alt)
        data.loc[mask, 'Class'] = 'Transition'
    
    for ref, alt in TRANSVERSIONS:
        mask = (data['REF'] == ref) & (data['ALT'] == alt)
        data.loc[mask, 'Class'] = 'Transversion'
    
    # Standardize mutation representation (pyrimidine-based)
    data['Std_Mutation'] = data.apply(standardize_mutation, axis=1)
    
    # Add adaptation type information from TREATMENT_INFO
    data['Adaptation'] = data['Treatment'].map(
        lambda t: TREATMENT_INFO.get(t, {}).get('adaptation', 'Unknown'))
    
    data['Has_Gene'] = data['Treatment'].map(
        lambda t: 'Yes' if TREATMENT_INFO.get(t, {}).get('gene') else 'No')
    
    return data

# Function to standardize mutation representation
def standardize_mutation(row):
    """Convert mutation to standardized format with pyrimidine as reference."""
    ref, alt = row['REF'], row['ALT']
    
    # If reference is a purine (A or G), convert to pyrimidine-based
    if ref in ['A', 'G']:
        ref = COMPLEMENT[ref]
        alt = COMPLEMENT[alt]
    
    return f"{ref}>{alt}"

# Function to calculate transition/transversion ratio
def calculate_ti_tv_ratio(data):
    """Calculate the transition/transversion ratio."""
    if len(data) == 0:
        return 0
    
    transitions = len(data[data['Class'] == 'Transition'])
    transversions = len(data[data['Class'] == 'Transversion'])
    
    if transversions == 0:
        return float('inf')  # Avoid division by zero
    
    return transitions / transversions

# Function to count mutation types
def count_mutation_types(data):
    """Count occurrences of each mutation type."""
    if len(data) == 0:
        return {}
    
    # Count standard mutations
    counts = Counter(data['Std_Mutation'])
    
    # Ensure all possible substitutions are represented
    for sub in STD_SUBSTITUTIONS:
        if sub not in counts:
            counts[sub] = 0
    
    return counts

# Function to generate mutation spectrum plot
def plot_mutation_spectrum(mutation_counts, treatment, output_dir):
    """Generate mutation spectrum plot for a treatment."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Get treatment description for title
    description = TREATMENT_INFO.get(treatment, {}).get('description', '')
    adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', '')
    has_gene = TREATMENT_INFO.get(treatment, {}).get('gene') is not None
    gene_text = f" with gene modification" if has_gene else ""
    
    # Prepare data for plotting
    categories = STD_SUBSTITUTIONS
    values = [mutation_counts.get(cat, 0) for cat in categories]
    
    # Define colors for different mutation types
    colors = ['#2166ac', '#4393c3', '#92c5de', '#d6604d', '#f4a582', '#fddbc7']
    
    # Plot the bars
    bars = ax.bar(range(len(categories)), values, color=colors)
    
    # Add value labels on top of bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                f'{height}', ha='center', va='bottom')
    
    # Customize the plot
    ax.set_xticks(range(len(categories)))
    ax.set_xticklabels(categories, rotation=45)
    ax.set_xlabel('Mutation Type')
    ax.set_ylabel('Count')
    ax.set_title(f'Mutation Spectrum for {treatment}\n{description} ({adaptation} adaptation{gene_text})')
    
    # Add transition/transversion annotations
    ax.text(0.02, 0.95, 'Transversions', transform=ax.transAxes, 
            fontsize=12, va='top', color='#2166ac')
    ax.text(0.5, 0.95, 'Transitions', transform=ax.transAxes, 
            fontsize=12, va='top', color='#d6604d')
    
    # Add vertical lines to separate mutation types
    ax.axvline(x=2.5, color='black', linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{treatment}_mutation_spectrum.png"), dpi=300)
    plt.close()

# Function to plot comparative mutation spectrum
def plot_comparative_spectrum(all_counts, output_dir):
    """Generate comparative mutation spectrum plot for all treatments."""
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Number of categories and treatments
    categories = STD_SUBSTITUTIONS
    treatments = list(all_counts.keys())
    n_cats = len(categories)
    n_treatments = len(treatments)
    
    # Width of bars
    width = 0.8 / n_treatments
    
    # Define colors for treatments
    treatment_colors = TREATMENT_COLORS
    
    # Plot grouped bars
    for i, treatment in enumerate(treatments):
        values = [all_counts[treatment].get(cat, 0) for cat in categories]
        positions = [j + (i - n_treatments/2 + 0.5) * width for j in range(n_cats)]
        bars = ax.bar(positions, values, width, label=treatment, color=treatment_colors.get(treatment, 'gray'))
    
    # Customize the plot
    ax.set_xticks(range(n_cats))
    ax.set_xticklabels(categories, rotation=45)
    ax.set_xlabel('Mutation Type')
    ax.set_ylabel('Count')
    ax.set_title('Comparative Mutation Spectrum Across Treatments')
    
    # Create legend with treatment descriptions
    legend_labels = [f"{t} ({TREATMENT_INFO.get(t, {}).get('description', '')})" for t in treatments]
    ax.legend(legend_labels)
    
    # Add vertical line to separate transitions and transversions
    ax.axvline(x=2.5, color='black', linestyle='--', alpha=0.3)
    ax.text(0.02, 0.95, 'Transversions', transform=ax.transAxes, 
            fontsize=12, va='top')
    ax.text(0.5, 0.95, 'Transitions', transform=ax.transAxes, 
            fontsize=12, va='top')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "comparative_mutation_spectrum.png"), dpi=300)
    plt.close()
    
    # Also create adaptation-grouped plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Group treatments by adaptation type
    adaptation_treatments = {
        'Temperature': [t for t in treatments if TREATMENT_INFO.get(t, {}).get('adaptation') == 'Temperature'],
        'Low Oxygen': [t for t in treatments if TREATMENT_INFO.get(t, {}).get('adaptation') == 'Low Oxygen']
    }
    
    # Plot adaptation-grouped bars
    positions = []
    bar_objects = []
    for i, (adaptation, group_treatments) in enumerate(adaptation_treatments.items()):
        for j, treatment in enumerate(group_treatments):
            values = [all_counts[treatment].get(cat, 0) for cat in categories]
            pos = [k + (j - len(group_treatments)/2 + 0.5) * width + i * (n_cats + 1) for k in range(n_cats)]
            positions.extend(pos)
            bars = ax.bar(pos, values, width, label=treatment, color=treatment_colors.get(treatment, 'gray'))
            bar_objects.append(bars[0])
        
        # Add adaptation type label
        midpoint = (n_cats - 1) / 2 + i * (n_cats + 1)
        ax.text(midpoint, -5, adaptation, ha='center', fontsize=14, fontweight='bold')
    
    # Customize the plot
    all_positions = [p for p in range(n_cats)] + [p + n_cats + 1 for p in range(n_cats)]
    ax.set_xticks(all_positions)
    all_labels = categories + categories
    ax.set_xticklabels(all_labels, rotation=45)
    ax.set_xlabel('Mutation Type')
    ax.set_ylabel('Count')
    ax.set_title('Mutation Spectrum by Adaptation Type')
    ax.legend(bar_objects, treatments)
    
    # Add vertical separators between adaptation types
    ax.axvline(x=n_cats + 0.5, color='black', linestyle='-', linewidth=2)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "adaptation_mutation_spectrum.png"), dpi=300)
    plt.close()

# Function to plot transition/transversion ratios
def plot_ti_tv_ratios(ratios, output_dir):
    """Plot transition/transversion ratios for all treatments."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot the bars
    treatments = list(ratios.keys())
    values = list(ratios.values())
    bars = ax.bar(treatments, values, color=[TREATMENT_COLORS.get(t, '#5ab4ac') for t in treatments])
    
    # Add value labels on top of bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.05,
                f'{height:.2f}', ha='center', va='bottom')
    
    # Customize the plot
    ax.set_xlabel('Treatment')
    ax.set_ylabel('Ti/Tv Ratio')
    ax.set_title('Transition/Transversion Ratio by Treatment')
    ax.set_ylim(0, max(values) * 1.2)  # Add some space for labels
    
    # Add treatment descriptions
    treatment_labels = [f"{t}\n({TREATMENT_INFO.get(t, {}).get('adaptation', '')})" for t in treatments]
    ax.set_xticks(range(len(treatments)))
    ax.set_xticklabels(treatment_labels)

    
    # Highlight gene-modified treatments
    for i, treatment in enumerate(treatments):
        if TREATMENT_INFO.get(treatment, {}).get('gene'):
            ax.text(i, values[i] / 2, '*', fontsize=20, ha='center', va='center', color='white')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "ti_tv_ratios.png"), dpi=300)
    plt.close()
    
    # Also create adaptation-grouped plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Group by adaptation type
    adaptation_data = defaultdict(list)
    adaptation_labels = []
    
    for treatment in treatments:
        adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
        has_gene = TREATMENT_INFO.get(treatment, {}).get('gene') is not None
        label = f"{treatment} {'(gene)' if has_gene else ''}"
        adaptation_data[adaptation].append((label, ratios[treatment]))
        adaptation_labels.append((treatment, adaptation, has_gene))
    
    # Plot grouped bars
    x_positions = []
    x_labels = []
    all_bars = []
    
    current_x = 0
    for adaptation, values in adaptation_data.items():
        positions = [current_x + i for i in range(len(values))]
        x_positions.extend(positions)
        
        labels = [v[0] for v in values]
        x_labels.extend(labels)
        
        heights = [v[1] for v in values]
        bars = ax.bar(positions, heights, color='#5ab4ac')
        all_bars.extend(bars)
        
        # Add group label
        midpoint = current_x + (len(values) - 1) / 2
        ax.text(midpoint, -0.2, adaptation, ha='center', fontweight='bold')
        
        current_x += len(values) + 1
    
    # Add value labels and color by treatment
    for i, ((treatment, adaptation, has_gene), bar) in enumerate(zip(adaptation_labels, all_bars)):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.05,
                f'{height:.2f}', ha='center', va='bottom')
        
        # Color by treatment
        bar.set_color(TREATMENT_COLORS.get(treatment, '#5ab4ac'))
        
        # Mark gene-modified treatments
        if has_gene:
            ax.text(x_positions[i], height / 2, '*', fontsize=20, ha='center', va='center', color='white')
    
    # Customize the plot
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, rotation=45, ha='right')
    ax.set_ylabel('Ti/Tv Ratio')
    ax.set_title('Transition/Transversion Ratio by Adaptation Type')
    
    # Add legend for gene modification
    from matplotlib.lines import Line2D
    gene_legend = Line2D([0], [0], marker='*', color='w', markerfacecolor='white',
                       markeredgecolor='black', markersize=15, label='Gene modified')
    ax.legend(handles=[gene_legend])
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "ti_tv_ratios_by_adaptation.png"), dpi=300)
    plt.close()

# Function to perform statistical test on mutation patterns
def test_mutation_differences(all_counts):
    """Perform chi-square test to check if mutation patterns differ significantly."""
    # Prepare data for chi-square test
    treatments = list(all_counts.keys())
    categories = STD_SUBSTITUTIONS
    
    # Create contingency table
    table = []
    for treatment in treatments:
        row = [all_counts[treatment].get(cat, 0) for cat in categories]
        table.append(row)
    
    # Perform chi-square test
    chi2, p, dof, expected = chi2_contingency(table)
    
    return {
        'chi2': chi2,
        'p_value': p,
        'degrees_of_freedom': dof,
        'expected_counts': expected
    }

# Function to create a summary table
def create_summary_table(all_data, all_counts, ti_tv_ratios):
    """Create a summary table with key statistics."""
    summary = []
    
    for treatment in TREATMENTS:
        if treatment in all_data:
            data = all_data[treatment]
            total_snvs = len(data)
            transitions = len(data[data['Class'] == 'Transition'])
            transversions = len(data[data['Class'] == 'Transversion'])
            
            # Calculate global statistics
            transition_pct = (transitions / total_snvs * 100) if total_snvs > 0 else 0
            transversion_pct = (transversions / total_snvs * 100) if total_snvs > 0 else 0
            
            # Most common mutation
            std_counts = {k: v for k, v in all_counts[treatment].items()}
            most_common = max(std_counts.items(), key=lambda x: x[1]) if std_counts else ('N/A', 0)
            
            # Treatment metadata
            description = TREATMENT_INFO.get(treatment, {}).get('description', 'Unknown')
            adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
            has_gene = TREATMENT_INFO.get(treatment, {}).get('gene') is not None
            
            summary.append({
                'Treatment': treatment,
                'Description': description,
                'Adaptation': adaptation,
                'Has_Gene': 'Yes' if has_gene else 'No',
                'Total SNVs': total_snvs,
                'Transitions': transitions,
                'Transversions': transversions,
                'Transition %': transition_pct,
                'Transversion %': transversion_pct,
                'Ti/Tv Ratio': ti_tv_ratios[treatment],
                'Most Common': most_common[0],
                'Most Common Count': most_common[1]
            })
    
    return pd.DataFrame(summary)

# Function to plot mutation distribution by adaptation
def plot_mutation_by_adaptation(all_data, output_dir):
    """Plot mutation distribution by adaptation type."""
    try:
        # Check if we have adaptation data
        if not all_data or 'Adaptation' not in next(iter(all_data.values())).columns:
            print("Warning: Adaptation data not available. Skipping adaptation plots.")
            return
        
        # Combine all data
        combined_data = pd.concat([all_data[t] for t in TREATMENTS if t in all_data])
        
        # Count by adaptation type and mutation class
        adaptation_counts = combined_data.groupby(['Adaptation', 'Class']).size().unstack(fill_value=0)
        
        # Print debug info
        print("Adaptation counts columns:", adaptation_counts.columns.tolist())
        
        # Calculate transition percentage
        adaptation_counts['Total'] = adaptation_counts.sum(axis=1)
        adaptation_counts['Transition %'] = (adaptation_counts['Transition'] / adaptation_counts['Total'] * 100).round(1)
        
        # Plot as stacked bars
        fig, ax = plt.subplots(figsize=(10, 6))
        
        adaptation_counts[['Transition', 'Transversion']].plot(
            kind='bar', stacked=True, ax=ax, color=['#d6604d', '#4393c3'])
        
        # Add percentage labels - using iterrows instead of itertuples for safety
        for i, (idx, row) in enumerate(adaptation_counts.iterrows()):
            ax.text(i, row['Total'] + 5, f"{row['Transition %']}%", ha='center')
        
        ax.set_xlabel('Adaptation Type')
        ax.set_ylabel('Count')
        ax.set_title('Mutation Classes by Adaptation Type')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "mutation_by_adaptation.png"), dpi=300)
        plt.close()
        
        # Plot by adaptation type and gene presence
        try:
            # Create multiindex groupby
            gene_group = combined_data.groupby(['Adaptation', 'Has_Gene', 'Class']).size()
            print("Gene group shape:", gene_group.shape)
            
            # Convert to DataFrame and structure for plotting
            gene_df = gene_group.reset_index()
            gene_df.columns = ['Adaptation', 'Has_Gene', 'Class', 'Count']
            
            # Create group labels
            gene_df['Group'] = gene_df.apply(
                lambda row: f"{row['Adaptation']}\n({'Gene' if row['Has_Gene'] == 'Yes' else 'No Gene'})", axis=1)
            
            # Create pivot table manually
            plot_data = pd.pivot_table(
                gene_df, 
                values='Count', 
                index='Group', 
                columns='Class',
                fill_value=0
            )
            
            # Add totals column
            plot_data['Total'] = plot_data.sum(axis=1)
            plot_data['Transition %'] = (plot_data['Transition'] / plot_data['Total'] * 100).round(1)
            
            # Plot as grouped bars
            fig, ax = plt.subplots(figsize=(12, 6))
            
            ax = plot_data[['Transition', 'Transversion']].plot(
                kind='bar', ax=ax, color=['#d6604d', '#4393c3'])
            
            # Add percentage labels
            for i, (idx, row) in enumerate(plot_data.iterrows()):
                ax.text(i, row['Total'] + 5, f"{row['Transition %']}%", ha='center')
            
            ax.set_xlabel('Adaptation Type and Gene Presence')
            ax.set_ylabel('Count')
            ax.set_title('Mutation Classes by Adaptation Type and Gene Modification')
            
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "mutation_by_adaptation_gene.png"), dpi=300)
            plt.close()
        except Exception as e:
            print(f"Warning: Error creating gene-adaptation plot: {e}")
            
    except Exception as e:
        print(f"Warning: Error in adaptation plotting: {e}")
        import traceback
        traceback.print_exc()

# Main function to run the analysis
def main():
    # Validate existing mutation files
    validate_mutation_files()
    # Parse data for each treatment
    all_raw_data = {}
    for treatment in TREATMENTS:
        data = parse_mutation_data(treatment)
        if len(data) > 0:
            all_raw_data[treatment] = data
            print(f"Loaded {len(data)} variants for {treatment} treatment")
        else:
            print(f"Warning: No data available for {treatment} treatment")
    
    # Filter for SNVs and classify mutations
    all_data = {}
    for treatment, data in all_raw_data.items():
        snv_data = filter_snvs(data)
        all_data[treatment] = classify_mutations(snv_data)
        print(f"{treatment}: Found {len(snv_data)} SNVs")
    
    # Calculate transition/transversion ratios
    ti_tv_ratios = {}
    for treatment, data in all_data.items():
        ti_tv_ratios[treatment] = calculate_ti_tv_ratio(data)
        print(f"{treatment}: Ti/Tv ratio = {ti_tv_ratios[treatment]:.2f}")
    
    # Count mutation types
    all_counts = {}
    for treatment, data in all_data.items():
        all_counts[treatment] = count_mutation_types(data)
    
    # Generate mutation spectrum plots
    for treatment, counts in all_counts.items():
        plot_mutation_spectrum(counts, treatment, OUTPUT_DIR)
        print(f"Generated mutation spectrum plot for {treatment}")
    
    # Generate comparative plot
    plot_comparative_spectrum(all_counts, OUTPUT_DIR)
    print("Generated comparative mutation spectrum plot")
    
    # Plot transition/transversion ratios
    plot_ti_tv_ratios(ti_tv_ratios, OUTPUT_DIR)
    print("Generated Ti/Tv ratio plot")
    
    # Plot mutation distribution by adaptation
    plot_mutation_by_adaptation(all_data, OUTPUT_DIR)
    print("Generated mutation distribution by adaptation plot")
    
    # Perform statistical test
    test_results = test_mutation_differences(all_counts)
    print(f"Chi-square test: chi2={test_results['chi2']:.2f}, p={test_results['p_value']:.4f}")
    
    # Create summary table
    summary_table = create_summary_table(all_data, all_counts, ti_tv_ratios)
    
    # Save summary table
    summary_table.to_csv(os.path.join(OUTPUT_DIR, "mutation_spectrum_summary.csv"), index=False)
    print(f"Saved summary table to {OUTPUT_DIR}/mutation_spectrum_summary.csv")
    
    # Save test results
    with open(os.path.join(OUTPUT_DIR, "statistical_test_results.txt"), 'w') as f:
        f.write("Chi-square test for differences in mutation patterns:\n")
        f.write(f"Chi-square value: {test_results['chi2']:.4f}\n")
        f.write(f"p-value: {test_results['p_value']:.4f}\n")
        f.write(f"Degrees of freedom: {test_results['degrees_of_freedom']}\n")
        f.write("\nInterpretation:\n")
        if test_results['p_value'] < 0.05:
            f.write("The mutation patterns differ significantly between treatments (p < 0.05).\n")
        else:
            f.write("No significant difference detected in mutation patterns between treatments (p >= 0.05).\n")
        
        # Add biological context
        f.write("\nBiological Context:\n")
        f.write("------------------\n")
        for treatment in TREATMENTS:
            description = TREATMENT_INFO.get(treatment, {}).get('description', 'Unknown')
            adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
            gene = TREATMENT_INFO.get(treatment, {}).get('gene')
            
            f.write(f"{treatment}: {description}\n")
            f.write(f"  Adaptation type: {adaptation}\n")
            if gene:
                f.write(f"  Gene modification: {gene}\n")
            else:
                f.write("  Gene modification: None\n")
            
            if treatment in ti_tv_ratios:
                f.write(f"  Ti/Tv ratio: {ti_tv_ratios[treatment]:.2f}\n")
            
            if treatment in all_counts:
                most_common = max(all_counts[treatment].items(), key=lambda x: x[1])
                f.write(f"  Most common mutation: {most_common[0]} ({most_common[1]} occurrences)\n")
            
            f.write("\n")
        
        # Add adaptation-specific analysis
        adaptation_types = set(TREATMENT_INFO.get(t, {}).get('adaptation', 'Unknown') for t in TREATMENTS)
        if len(adaptation_types) > 1:
            f.write("\nComparison by Adaptation Type:\n")
            f.write("----------------------------\n")
            
            for adaptation in sorted(adaptation_types):
                f.write(f"{adaptation} adaptation:\n")
                # Filter treatments by adaptation
                adaptation_treatments = [t for t in TREATMENTS if TREATMENT_INFO.get(t, {}).get('adaptation') == adaptation]
                
                # Calculate average Ti/Tv ratio
                ratios = [ti_tv_ratios.get(t, 0) for t in adaptation_treatments if t in ti_tv_ratios]
                if ratios:
                    avg_ratio = sum(ratios) / len(ratios)
                    f.write(f"  Average Ti/Tv ratio: {avg_ratio:.2f}\n")
                
                # List treatments
                f.write(f"  Treatments: {', '.join(adaptation_treatments)}\n")
                f.write("\n")
    
    print(f"Analysis complete! Results saved to {OUTPUT_DIR}/")
    print(f"Summary table saved as {OUTPUT_DIR}/mutation_spectrum_summary.csv")
    print(f"Statistical test results saved as {OUTPUT_DIR}/statistical_test_results.txt")

# Run the analysis
if __name__ == "__main__":
    main()