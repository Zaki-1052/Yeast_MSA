#!/usr/bin/env python3
"""
Mutation Spectrum Analysis Module

This module analyzes single nucleotide mutation patterns across different treatment conditions
in yeast adaptation experiments. It calculates mutation spectra, transition/transversion ratios,
and performs statistical comparisons between treatments.

Functions:
- Data loading and validation
- Mutation classification and standardization
- Spectra visualization and comparison
- Statistical analysis of mutation patterns
"""

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

# Treatment groups
TREATMENTS = ['WT-37', 'WTA', 'STC', 'CAS']

# Treatment information for biological context
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

# Nucleotide definitions
COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
TRANSITIONS = [('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')]
TRANSVERSIONS = [('A', 'C'), ('C', 'A'), ('A', 'T'), ('T', 'A'), 
                ('G', 'T'), ('T', 'G'), ('G', 'C'), ('C', 'G')]

# All possible single nucleotide substitutions
ALL_SUBSTITUTIONS = TRANSITIONS + TRANSVERSIONS

# Standardized substitution representation (pyrimidine-based)
STD_SUBSTITUTIONS = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']

def validate_mutation_files():
    """
    Check if mutation files have the expected format with valid nucleotides.
    
    Validates each treatment's mutation file to ensure it contains properly formatted
    single nucleotides for REF and ALT columns. Invalid files are removed to be
    regenerated from VCF data.
    """
    for treatment in TREATMENTS:
        file_path = f"analysis/mutation_spectrum_analysis/{treatment}_mutations.txt"
        if not os.path.exists(file_path):
            continue
            
        try:
            # Read the first few lines
            data = pd.read_csv(
                file_path, 
                sep='\t', 
                header=None, 
                names=['CHROM', 'POS', 'REF', 'ALT'], 
                nrows=5
            )
            
            # Check if REF/ALT columns contain valid nucleotides
            valid_ref = (data['REF'].str.len().eq(1).all() and 
                        data['REF'].str.upper().isin(['A', 'C', 'G', 'T']).all())
            valid_alt = (data['ALT'].str.len().eq(1).all() and 
                        data['ALT'].str.upper().isin(['A', 'C', 'G', 'T']).all())
            
            if not valid_ref or not valid_alt:
                print(f"Warning: {file_path} contains invalid nucleotides.")
                print(f"REF values: {data['REF'].tolist()}, ALT values: {data['ALT'].tolist()}")
                print(f"Removing invalid file so it will be re-extracted...")
                os.remove(file_path)
        except Exception as e:
            print(f"Error validating {file_path}: {e}")
            print(f"Removing invalid file...")
            os.remove(file_path)

def find_mutation_data_file(treatment):
    """
    Find mutation data file by checking multiple possible locations.
    
    Args:
        treatment (str): Treatment identifier (e.g., 'WT-37', 'WTA')
        
    Returns:
        str or None: Path to the mutation data file if found, None otherwise
    """
    # Define possible file patterns
    file_patterns = [
        f"analysis/mutation_spectrum_analysis/{treatment}_mutations.txt",
        f"analysis/MSA/mutation_spectrum_analysis/{treatment}_mutations.txt",
        f"mutation_spectrum_analysis/{treatment}_mutations.txt"
    ]
    
    # Handle backward compatibility with old 'WT' naming for WT-37
    if treatment == 'WT-37':
        file_patterns.extend([
            "analysis/mutation_spectrum_analysis/WT_mutations.txt",
            "analysis/MSA/mutation_spectrum_analysis/WT_mutations.txt",
            "mutation_spectrum_analysis/WT_mutations.txt"
        ])
    
    # Check each location
    for pattern in file_patterns:
        if os.path.exists(pattern):
            print(f"Found mutation data for {treatment} at {pattern}")
            return pattern
    
    return None

def parse_mutation_data(treatment):
    """
    Parse mutation data for a specific treatment, handling various file formats.
    
    This function attempts to load mutation data from text files, handles different
    file formats, and fixes common formatting issues. If no valid file is found
    or if parsing fails, it extracts data directly from VCF files.
    
    Args:
        treatment (str): Treatment identifier (e.g., 'WT-37', 'WTA')
        
    Returns:
        pd.DataFrame: DataFrame containing mutation data with columns:
            CHROM, POS, REF, ALT, Treatment
    """
    # Find the mutation data file
    filename = find_mutation_data_file(treatment)
    
    if not filename:
        print(f"Warning: No mutation data file found for {treatment}")
        # Fall back to VCF extraction
        return extract_from_vcf(treatment)
    
    try:
        # Examine the file structure
        with open(filename, 'r') as f:
            first_line = f.readline().strip()
            columns = first_line.split('\t')
            print(f"File format: {len(columns)} columns in first line")
        
        # Read data according to the detected format
        raw_data = pd.read_csv(filename, sep='\t', header=None)
        data = pd.DataFrame()
        
        # Process based on column count
        if len(columns) == 5:  # CHROM, POS, REF, ALT, Treatment
            data['CHROM'] = raw_data.iloc[:, 0]
            data['POS'] = raw_data.iloc[:, 1].astype(int)
            data['REF'] = raw_data.iloc[:, 2]
            data['ALT'] = raw_data.iloc[:, 3]
        elif len(columns) == 6:  # Line num, CHROM, POS, REF, ALT, Treatment
            data['CHROM'] = raw_data.iloc[:, 1]
            data['POS'] = raw_data.iloc[:, 2].astype(int)
            data['REF'] = raw_data.iloc[:, 3]
            data['ALT'] = raw_data.iloc[:, 4]
        else:
            # Default approach for unexpected formats
            data = pd.read_csv(
                filename, 
                sep='\t', 
                header=None, 
                names=['CHROM', 'POS', 'REF', 'ALT']
            )
        
        # Always use the provided treatment parameter
        data['Treatment'] = treatment
        
        # Validation and correction logic
        if len(data) > 0:
            # Debug output
            print(f"Sample data after loading:")
            print(data.head(2))
            
            # Fix issues with treatment names in REF column
            if any(data['REF'].astype(str).isin(TREATMENTS)):
                print(f"Warning: File {filename} appears to have swapped REF/ALT columns. Fixing...")
                _fix_ref_column_issues(data, raw_data, treatment)
            
            # Fix issues with treatment names in ALT column
            if any(data['ALT'].astype(str).isin(TREATMENTS)):
                print(f"Warning: Treatment name found in ALT column in {filename}. Fixing...")
                data = _fix_alt_column_issues(data, raw_data, treatment)
                
                # Try to get more accurate ALT values from VCF
                vcf_data = extract_from_vcf(treatment)
                if not vcf_data.empty:
                    data = _update_alt_from_vcf(data, vcf_data)
        
        print(f"Loaded {len(data)} mutations for {treatment}")
        return data
        
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        print(f"Falling back to VCF extraction...")
        return extract_from_vcf(treatment)

def _fix_ref_column_issues(data, raw_data, treatment):
    """
    Fix issues where REF column contains treatment names.
    
    Args:
        data (pd.DataFrame): DataFrame to fix
        raw_data (pd.DataFrame): Raw data from file
        treatment (str): Current treatment
        
    Returns:
        None (modifies data in-place)
    """
    if len(raw_data.columns) >= 6:  # Line num, CHROM, POS, REF, ALT, Treatment
        data['CHROM'] = raw_data.iloc[:, 1]
        data['POS'] = raw_data.iloc[:, 2].astype(int)
        
        # Check if treatment appears in what should be the REF column
        if any(t in TREATMENTS for t in raw_data.iloc[:, 3].astype(str).tolist()):
            data['REF'] = raw_data.iloc[:, 4]  # Use column 4 as REF
            
            # Create best-guess ALT values
            alt_values = []
            for i, row in data.iterrows():
                actual_alt = raw_data.iloc[i, 3]
                if actual_alt not in TREATMENTS and str(actual_alt).upper() in ['A', 'C', 'G', 'T']:
                    alt_values.append(actual_alt)
                else:
                    # Choose a nucleotide different from REF
                    ref = str(row['REF']).upper()
                    alt_bases = [b for b in ['A', 'C', 'G', 'T'] if b != ref]
                    alt_values.append(alt_bases[0] if alt_bases else 'N')
            data['ALT'] = alt_values
        else:
            # Normal column order
            data['REF'] = raw_data.iloc[:, 3]
            data['ALT'] = raw_data.iloc[:, 4]
    
    elif len(raw_data.columns) == 5:  # CHROM, POS, REF, ALT, Treatment
        data['CHROM'] = raw_data.iloc[:, 0]
        data['POS'] = raw_data.iloc[:, 1].astype(int)
        
        # Check if treatment appears in what should be the REF column
        if any(t in TREATMENTS for t in raw_data.iloc[:, 2].astype(str).tolist()):
            data['REF'] = raw_data.iloc[:, 3]  # Use column 3 as REF
            
            # Create best-guess ALT values
            alt_values = []
            for i, row in data.iterrows():
                actual_alt = raw_data.iloc[i, 2]
                if actual_alt not in TREATMENTS and str(actual_alt).upper() in ['A', 'C', 'G', 'T']:
                    alt_values.append(actual_alt)
                else:
                    # Choose a nucleotide different from REF
                    ref = str(row['REF']).upper()
                    alt_bases = [b for b in ['A', 'C', 'G', 'T'] if b != ref]
                    alt_values.append(alt_bases[0] if alt_bases else 'N')
            data['ALT'] = alt_values
        else:
            # Normal column order
            data['REF'] = raw_data.iloc[:, 2]
            data['ALT'] = raw_data.iloc[:, 3]
    
    # Reset treatment
    data['Treatment'] = treatment

def _fix_alt_column_issues(data, raw_data, treatment):
    """
    Fix issues where ALT column contains treatment names.
    
    Args:
        data (pd.DataFrame): DataFrame to fix
        raw_data (pd.DataFrame): Raw data from file
        treatment (str): Current treatment
        
    Returns:
        pd.DataFrame: Fixed DataFrame
    """
    new_data = pd.DataFrame()
    
    if len(raw_data.columns) >= 6:  # Line num, CHROM, POS, REF, ALT, Treatment
        new_data['CHROM'] = raw_data.iloc[:, 1]
        new_data['POS'] = raw_data.iloc[:, 2].astype(int)
        new_data['REF'] = raw_data.iloc[:, 3]  # REF is usually correct
        
        # Create best-guess ALT values based on common mutations
        alt_values = []
        for i, ref in enumerate(raw_data.iloc[:, 3]):
            ref = str(ref).upper()
            if ref in ['A', 'C', 'G', 'T']:
                # Try using the value in ALT column if it's valid
                alt_candidate = str(raw_data.iloc[i, 4]).upper()
                if alt_candidate in ['A', 'C', 'G', 'T'] and alt_candidate != ref:
                    alt_values.append(alt_candidate)
                else:
                    # Make best guess based on common transitions
                    if ref == 'A': alt_values.append('G')
                    elif ref == 'G': alt_values.append('A')
                    elif ref == 'C': alt_values.append('T')
                    elif ref == 'T': alt_values.append('C')
                    else: alt_values.append('N')
            else:
                alt_values.append('N')
        
        new_data['ALT'] = alt_values
        
    elif len(raw_data.columns) == 5:  # CHROM, POS, REF, ALT, Treatment
        new_data['CHROM'] = raw_data.iloc[:, 0]
        new_data['POS'] = raw_data.iloc[:, 1].astype(int)
        new_data['REF'] = raw_data.iloc[:, 2]  # REF is usually correct
        
        # Create best-guess ALT values based on common mutations
        alt_values = []
        for i, ref in enumerate(raw_data.iloc[:, 2]):
            ref = str(ref).upper()
            if ref in ['A', 'C', 'G', 'T']:
                # Try using the value in ALT column if it's valid
                alt_candidate = str(raw_data.iloc[i, 3]).upper()
                if alt_candidate in ['A', 'C', 'G', 'T'] and alt_candidate != ref:
                    alt_values.append(alt_candidate)
                else:
                    # Make best guess based on common transitions
                    if ref == 'A': alt_values.append('G')
                    elif ref == 'G': alt_values.append('A')
                    elif ref == 'C': alt_values.append('T')
                    elif ref == 'T': alt_values.append('C')
                    else: alt_values.append('N')
            else:
                alt_values.append('N')
        
        new_data['ALT'] = alt_values
    
    # Always set treatment
    new_data['Treatment'] = treatment
    return new_data

def _update_alt_from_vcf(data, vcf_data):
    """
    Update ALT values from VCF data where possible.
    
    Args:
        data (pd.DataFrame): DataFrame with potentially incorrect ALT values
        vcf_data (pd.DataFrame): Data extracted from VCF
        
    Returns:
        pd.DataFrame: Updated DataFrame
    """
    # Create a dictionary for quick position lookup
    vcf_dict = {}
    for _, row in vcf_data.iterrows():
        key = (row['CHROM'], int(row['POS']))
        vcf_dict[key] = row['ALT']
    
    # Update ALT values from VCF where available
    for i, row in data.iterrows():
        key = (row['CHROM'], int(row['POS']))
        if key in vcf_dict:
            data.at[i, 'ALT'] = vcf_dict[key]
    
    return data

def extract_from_vcf(treatment):
    """
    Extract mutation data directly from VCF files.
    
    This function is used when no pre-extracted mutation data file is found,
    or if the existing file has invalid format.
    
    Args:
        treatment (str): Treatment identifier (e.g., 'WT-37', 'WTA')
        
    Returns:
        pd.DataFrame: DataFrame containing mutation data with columns:
            CHROM, POS, REF, ALT, Treatment
    """
    # Define possible VCF locations to check
    vcf_patterns = [
        f"results/merged/analysis/{treatment}/highconf.vcf.gz",
        f"results/merged/analysis/{treatment}_highconf.vcf.gz",
        f"results/merged/analysis/{treatment}/specific.vcf.gz",
        f"results/merged/analysis/{treatment}_specific.vcf.gz"
    ]
    
    # Handle backward compatibility with old 'WT' naming for WT-37
    if treatment == 'WT-37':
        vcf_patterns.extend([
            "results/merged/analysis/WT/highconf.vcf.gz",
            "results/merged/analysis/WT_highconf.vcf.gz",
            "results/merged/analysis/WT/specific.vcf.gz",
            "results/merged/analysis/WT_specific.vcf.gz"
        ])
    
    # Try each VCF location
    for vcf_file in vcf_patterns:
        if not os.path.exists(vcf_file):
            continue
            
        print(f"Extracting mutation data for {treatment} from {vcf_file}")
        try:
            # Extract mutation data using bcftools
            cmd = f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' {vcf_file}"
            output = subprocess.check_output(cmd, shell=True).decode('utf-8')
            
            # Parse the output
            rows = []
            for line in output.strip().split('\n'):
                if not line:  # Skip empty lines
                    continue
                    
                parts = line.split('\t')
                if len(parts) == 4:
                    rows.append(parts)
            
            # Create dataframe if we have rows
            if not rows:
                continue
                
            data = pd.DataFrame(rows, columns=['CHROM', 'POS', 'REF', 'ALT'])
            data['POS'] = data['POS'].astype(int)
            data['Treatment'] = treatment
            
            # Save extracted data for future use (only core columns)
            os.makedirs("analysis/mutation_spectrum_analysis", exist_ok=True)
            data[['CHROM', 'POS', 'REF', 'ALT']].to_csv(
                f"analysis/mutation_spectrum_analysis/{treatment}_mutations.txt", 
                sep='\t', index=False, header=False
            )
            
            print(f"Extracted and saved {len(data)} mutations for {treatment}")
            return data
            
        except Exception as e:
            print(f"Error extracting from {vcf_file}: {e}")
    
    print(f"Could not find or extract mutation data for {treatment}")
    return pd.DataFrame()

def filter_snvs(data, debug=True, include_mnvs=False):
    """
    Filter data to include only single nucleotide variants.
    
    This function filters the input data to include only single nucleotide variants
    (SNVs) with valid nucleotide bases. It can optionally decompose multi-nucleotide
    variants (MNVs) into individual SNVs.
    
    Args:
        data (pd.DataFrame): DataFrame containing variant data
        debug (bool): Whether to print debug information
        include_mnvs (bool): Whether to include decomposed multi-nucleotide variants
        
    Returns:
        pd.DataFrame: Filtered DataFrame containing only SNVs
    """
    if data.empty:
        return data
        
    if debug:
        print(f"Initial variant count: {len(data)}")
        print(f"REF column type: {type(data['REF'].iloc[0])}")
        print(f"Sample REF values: {data['REF'].head().tolist()}")
        print(f"Sample ALT values: {data['ALT'].head().tolist()}")
    
    # Ensure REF/ALT are string type
    processed_data = data.copy()
    if not pd.api.types.is_string_dtype(processed_data['REF']) or not pd.api.types.is_string_dtype(processed_data['ALT']):
        if debug:
            print("Converting REF/ALT to string types")
        processed_data['REF'] = processed_data['REF'].astype(str)
        processed_data['ALT'] = processed_data['ALT'].astype(str)
    
    # First separate SNVs and MNVs
    length_filter = (processed_data['REF'].str.len() == 1) & (processed_data['ALT'].str.len() == 1)
    snv_data = processed_data[length_filter].copy()
    mnv_data = processed_data[~length_filter].copy()
    
    if debug:
        print(f"Found {len(snv_data)} single-nucleotide variants")
        print(f"Found {len(mnv_data)} multi-nucleotide variants")
    
    # Handle MNVs decomposition if requested
    decomposed_mnvs = pd.DataFrame()
    if include_mnvs and not mnv_data.empty:
        decomposed_mnvs = _decompose_mnvs(mnv_data, debug)
    
    # Combine SNVs and decomposed MNVs if needed
    if include_mnvs and not decomposed_mnvs.empty:
        combined_data = pd.concat([snv_data, decomposed_mnvs])
    else:
        combined_data = snv_data
    
    if debug:
        print(f"After length filter: {len(combined_data)} variants remain")
        print(f"Skipped {len(data) - len(combined_data)} variants that couldn't be processed as SNVs")
        
        if not mnv_data.empty:
            print("Examples of multi-nucleotide variants:")
            for _, row in mnv_data.head(5).iterrows():
                print(f"  {row['CHROM']}:{row['POS']} REF={row['REF']} ALT={row['ALT']}")
    
    # Filter for valid ACGT bases
    if not combined_data.empty:
        final_data = _filter_valid_bases(combined_data, debug)
    else:
        final_data = combined_data
    
    print(f"Filtered to {len(final_data)} single nucleotide variants")
    return final_data

def _decompose_mnvs(mnv_data, debug=False):
    """
    Decompose multi-nucleotide variants (MNVs) into single nucleotide variants.
    
    Args:
        mnv_data (pd.DataFrame): DataFrame containing MNV data
        debug (bool): Whether to print debug information
        
    Returns:
        pd.DataFrame: DataFrame with decomposed MNVs
    """
    decomposed_rows = []
    
    for _, row in mnv_data.iterrows():
        ref = row['REF']
        alt = row['ALT']
        
        # Skip indels (large length differences)
        if abs(len(ref) - len(alt)) > 2:
            continue
        
        # Process substitution-type MNVs (e.g., AT>AG)
        if len(ref) == len(alt) and len(ref) > 1:
            # Find positions where bases differ
            diff_positions = [i for i in range(len(ref)) 
                             if i < len(alt) and ref[i] != alt[i]]
            
            # Create a separate SNV for each position that differs
            for pos in diff_positions:
                new_row = row.copy()
                new_row['REF'] = ref[pos]
                new_row['ALT'] = alt[pos]
                new_row['POS'] = int(row['POS']) + pos  # Adjust position
                new_row['MNV_Source'] = f"{ref}>{alt}"  # Track the source MNV
                decomposed_rows.append(new_row)
        
        # Skip deletions and insertions
        # (These could be handled here if needed in the future)
    
    if decomposed_rows:
        result = pd.DataFrame(decomposed_rows)
        if debug:
            print(f"Decomposed {len(result)} SNVs from multi-nucleotide variants")
        return result
    
    return pd.DataFrame()

def _filter_valid_bases(data, debug=False):
    """
    Filter variants to include only those with valid ACGT nucleotide bases.
    
    Args:
        data (pd.DataFrame): DataFrame to filter
        debug (bool): Whether to print debug information
        
    Returns:
        pd.DataFrame: Filtered DataFrame
    """
    # Standardize to uppercase
    data['REF'] = data['REF'].str.upper()
    data['ALT'] = data['ALT'].str.upper()
    
    # Keep only variants with valid nucleotide bases
    valid_bases = (data['REF'].isin(['A', 'C', 'G', 'T']) & 
                  data['ALT'].isin(['A', 'C', 'G', 'T']))
    filtered_data = data[valid_bases].copy()
    
    if debug:
        print(f"After ACGT filter: {len(filtered_data)} variants remain")
        print(f"Removed {len(data) - len(filtered_data)} variants with non-ACGT bases")
        
        if len(data) > 0 and len(filtered_data) == 0:
            print("Examples of non-ACGT bases:")
            non_acgt = data[~valid_bases].head(5)
            for _, row in non_acgt.iterrows():
                print(f"  {row['CHROM']}:{row['POS']} REF={row['REF']} ALT={row['ALT']}")
    
    return filtered_data

def classify_mutations(data):
    """
    Classify each mutation as transition or transversion and add annotation.
    
    This function classifies each mutation as either a transition or transversion,
    standardizes mutation representation to be pyrimidine-based, and adds metadata
    about adaptation type and gene modification status.
    
    Args:
        data (pd.DataFrame): DataFrame containing mutation data
        
    Returns:
        pd.DataFrame: DataFrame with added classification and annotation columns
    """
    if data.empty:
        return data
    
    # Create a working copy of the data
    result = data.copy()
    
    # Create mutation type column (REF>ALT format)
    result['Mutation'] = result['REF'] + '>' + result['ALT']
    
    # Classify mutations as transitions or transversions
    result['Class'] = 'Unknown'
    
    # Mark transitions (purine↔purine or pyrimidine↔pyrimidine)
    for ref, alt in TRANSITIONS:
        mask = (result['REF'] == ref) & (result['ALT'] == alt)
        result.loc[mask, 'Class'] = 'Transition'
    
    # Mark transversions (purine↔pyrimidine)
    for ref, alt in TRANSVERSIONS:
        mask = (result['REF'] == ref) & (result['ALT'] == alt)
        result.loc[mask, 'Class'] = 'Transversion'
    
    # Standardize mutation representation (pyrimidine-based)
    result['Std_Mutation'] = result.apply(standardize_mutation, axis=1)
    
    # Add metadata based on treatment
    result['Adaptation'] = result['Treatment'].map(
        lambda t: TREATMENT_INFO.get(t, {}).get('adaptation', 'Unknown'))
    
    result['Has_Gene'] = result['Treatment'].map(
        lambda t: 'Yes' if TREATMENT_INFO.get(t, {}).get('gene') else 'No')
    
    return result

def standardize_mutation(row):
    """
    Convert mutation to standardized format with pyrimidine as reference.
    
    This function standardizes mutation representation by ensuring the reference
    base is a pyrimidine (C or T). If the reference is a purine (A or G), both
    the reference and alternate bases are converted to their complements.
    
    Args:
        row (pd.Series): DataFrame row containing 'REF' and 'ALT' columns
        
    Returns:
        str: Standardized mutation in the format "REF>ALT"
    """
    ref, alt = row['REF'], row['ALT']
    
    # If reference is a purine (A or G), convert to pyrimidine-based
    if ref in ['A', 'G']:
        ref = COMPLEMENT[ref]  # Convert A→T, G→C
        alt = COMPLEMENT[alt]  # Convert the alternate base too
    
    return f"{ref}>{alt}"

def calculate_ti_tv_ratio(data):
    """
    Calculate the transition/transversion ratio.
    
    Computes the ratio of transition mutations to transversion mutations
    in the dataset. Returns 0 for empty datasets and infinity if there
    are no transversions.
    
    Args:
        data (pd.DataFrame): DataFrame with mutations classified as 
                            'Transition' or 'Transversion'
                            
    Returns:
        float: Ratio of transitions to transversions
    """
    if data.empty:
        return 0
    
    transitions = len(data[data['Class'] == 'Transition'])
    transversions = len(data[data['Class'] == 'Transversion'])
    
    if transversions == 0:
        return float('inf')  # Avoid division by zero
    
    return transitions / transversions

def count_mutation_types(data):
    """
    Count occurrences of each standardized mutation type.
    
    Counts the frequency of each standardized mutation type (C>A, C>G, etc.)
    and ensures all six possible types are represented in the results.
    
    Args:
        data (pd.DataFrame): DataFrame containing mutations with 'Std_Mutation' column
        
    Returns:
        dict: Dictionary with mutation types as keys and counts as values
    """
    if data.empty:
        return {sub: 0 for sub in STD_SUBSTITUTIONS}
    
    # Count standardized mutations
    counts = Counter(data['Std_Mutation'])
    
    # Ensure all possible substitutions are represented
    for sub in STD_SUBSTITUTIONS:
        if sub not in counts:
            counts[sub] = 0
    
    return counts

def plot_mutation_spectrum(mutation_counts, treatment, output_dir):
    """
    Generate mutation spectrum plot for a specific treatment.
    
    Creates a bar chart showing the distribution of different mutation types
    (C>A, C>G, etc.) for the given treatment, with annotations for 
    transitions/transversions and treatment metadata.
    
    Args:
        mutation_counts (dict): Dictionary with mutation types and their counts
        treatment (str): Treatment identifier (e.g., 'WT-37', 'WTA')
        output_dir (str): Directory to save the output plot
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Get treatment metadata for title
    metadata = TREATMENT_INFO.get(treatment, {})
    description = metadata.get('description', '')
    adaptation = metadata.get('adaptation', '')
    has_gene = metadata.get('gene') is not None
    gene_text = " with gene modification" if has_gene else ""
    
    # Prepare plot data
    categories = STD_SUBSTITUTIONS
    values = [mutation_counts.get(cat, 0) for cat in categories]
    
    # Color scheme: blues for transversions, reds for transitions
    colors = ['#2166ac', '#4393c3', '#92c5de', '#d6604d', '#f4a582', '#fddbc7']
    
    # Create the bar plot
    bars = ax.bar(range(len(categories)), values, color=colors)
    
    # Add value labels on top of each bar
    for bar in bars:
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width()/2., 
            height + 0.1,
            f'{height}', 
            ha='center', 
            va='bottom'
        )
    
    # Customize plot appearance
    ax.set_xticks(range(len(categories)))
    ax.set_xticklabels(categories, rotation=45)
    ax.set_xlabel('Mutation Type')
    ax.set_ylabel('Count')
    ax.set_title(
        f'Mutation Spectrum for {treatment}\n'
        f'{description} ({adaptation} adaptation{gene_text})'
    )
    
    # Add annotations for transition/transversion categories
    ax.text(0.02, 0.95, 'Transversions', transform=ax.transAxes, 
            fontsize=12, va='top', color='#2166ac')
    ax.text(0.5, 0.95, 'Transitions', transform=ax.transAxes, 
            fontsize=12, va='top', color='#d6604d')
    
    # Add divider between transition and transversion sections
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

def main():
    """
    Main function to run the mutation spectrum analysis pipeline.
    
    This function orchestrates the entire analysis workflow:
    1. Load and validate mutation data for each treatment
    2. Filter for SNVs and classify mutations
    3. Calculate transition/transversion ratios
    4. Generate visualizations
    5. Perform statistical tests
    6. Create and save summary reports
    """
    print("Starting mutation spectrum analysis...")
    
    # Step 1: Validate existing mutation files and load data
    validate_mutation_files()
    all_raw_data = _load_treatment_data()
    
    # Step 2: Filter for SNVs and classify mutations
    all_data = _process_mutation_data(all_raw_data)
    
    # Step 3: Calculate statistics
    ti_tv_ratios = _calculate_ti_tv_ratios(all_data)
    all_counts = _count_mutation_types(all_data)
    
    # Step 4: Generate visualizations
    _generate_plots(all_counts, all_data, ti_tv_ratios)
    
    # Step 5: Perform statistical tests
    test_results = test_mutation_differences(all_counts)
    print(f"Chi-square test: chi2={test_results['chi2']:.2f}, p={test_results['p_value']:.4f}")
    
    # Step 6: Create and save summary reports
    _save_summary_reports(all_raw_data, all_data, all_counts, ti_tv_ratios, test_results)
    
    print(f"Analysis complete! Results saved to {OUTPUT_DIR}/")
    print(f"Summary table saved as {OUTPUT_DIR}/mutation_spectrum_summary.csv")
    print(f"Statistical test results saved as {OUTPUT_DIR}/statistical_test_results.txt")

def _load_treatment_data():
    """Load mutation data for all treatments."""
    all_raw_data = {}
    for treatment in TREATMENTS:
        data = parse_mutation_data(treatment)
        if not data.empty:
            all_raw_data[treatment] = data
            print(f"Loaded {len(data)} variants for {treatment} treatment")
        else:
            print(f"Warning: No data available for {treatment} treatment")
    return all_raw_data

def _process_mutation_data(all_raw_data):
    """Filter for SNVs and classify mutations."""
    all_data = {}
    for treatment, data in all_raw_data.items():
        # Include decomposed multi-nucleotide variants
        snv_data = filter_snvs(data, include_mnvs=True)
        all_data[treatment] = classify_mutations(snv_data)
        print(f"{treatment}: Found {len(snv_data)} SNVs")
    return all_data

def _calculate_ti_tv_ratios(all_data):
    """Calculate transition/transversion ratios for all treatments."""
    ti_tv_ratios = {}
    for treatment, data in all_data.items():
        ti_tv_ratios[treatment] = calculate_ti_tv_ratio(data)
        print(f"{treatment}: Ti/Tv ratio = {ti_tv_ratios[treatment]:.2f}")
    return ti_tv_ratios

def _count_mutation_types(all_data):
    """Count occurrences of each mutation type for all treatments."""
    all_counts = {}
    for treatment, data in all_data.items():
        all_counts[treatment] = count_mutation_types(data)
    return all_counts

def _generate_plots(all_counts, all_data, ti_tv_ratios):
    """Generate all visualization plots."""
    # Individual mutation spectrum plots
    for treatment, counts in all_counts.items():
        plot_mutation_spectrum(counts, treatment, OUTPUT_DIR)
        print(f"Generated mutation spectrum plot for {treatment}")
    
    # Comparative plots
    plot_comparative_spectrum(all_counts, OUTPUT_DIR)
    print("Generated comparative mutation spectrum plot")
    
    # Ti/Tv ratio plots
    plot_ti_tv_ratios(ti_tv_ratios, OUTPUT_DIR)
    print("Generated Ti/Tv ratio plot")
    
    # Adaptation-based plots
    plot_mutation_by_adaptation(all_data, OUTPUT_DIR)
    print("Generated mutation distribution by adaptation plot")

def _save_summary_reports(all_raw_data, all_data, all_counts, ti_tv_ratios, test_results):
    """Create and save summary reports and statistics."""
    # Create main summary table
    summary_table = create_summary_table(all_data, all_counts, ti_tv_ratios)
    
    # Calculate filtering statistics
    filtering_stats = _calculate_filtering_stats(all_raw_data, all_data)
    
    # Save tables to CSV
    summary_table.to_csv(os.path.join(OUTPUT_DIR, "mutation_spectrum_summary.csv"), index=False)
    filtering_stats.to_csv(os.path.join(OUTPUT_DIR, "variant_filtering_stats.csv"), index=False)
    
    # Generate and save statistical test report
    _save_statistical_report(test_results, ti_tv_ratios, all_counts)

def _calculate_filtering_stats(all_raw_data, all_data):
    """Calculate statistics on variant filtering."""
    original_counts = {t: len(all_raw_data.get(t, pd.DataFrame())) 
                      for t in TREATMENTS if t in all_raw_data}
    final_counts = {t: len(all_data.get(t, pd.DataFrame())) 
                   for t in TREATMENTS if t in all_data}
    
    filtered_info = []
    for treatment in TREATMENTS:
        if treatment in original_counts and treatment in final_counts:
            filtered_info.append({
                'Treatment': treatment,
                'Total_Variants': original_counts[treatment],
                'SNVs_After_Filtering': final_counts[treatment],
                'MNVs_and_Other_Filtered': original_counts[treatment] - final_counts[treatment],
                'Percent_Kept': round(final_counts[treatment] / original_counts[treatment] * 100, 1)
            })
    
    return pd.DataFrame(filtered_info)

def _save_statistical_report(test_results, ti_tv_ratios, all_counts):
    """Generate and save statistical test report with biological context."""
    report_path = os.path.join(OUTPUT_DIR, "statistical_test_results.txt")
    
    with open(report_path, 'w') as f:
        # Chi-square test results
        f.write("Chi-square test for differences in mutation patterns:\n")
        f.write(f"Chi-square value: {test_results['chi2']:.4f}\n")
        f.write(f"p-value: {test_results['p_value']:.4f}\n")
        f.write(f"Degrees of freedom: {test_results['degrees_of_freedom']}\n")
        
        # Interpretation
        f.write("\nInterpretation:\n")
        if test_results['p_value'] < 0.05:
            f.write("The mutation patterns differ significantly between treatments (p < 0.05).\n")
        else:
            f.write("No significant difference detected in mutation patterns between treatments (p >= 0.05).\n")
        
        # Biological context section
        f.write("\nBiological Context:\n")
        f.write("------------------\n")
        
        # Treatment details
        for treatment in TREATMENTS:
            metadata = TREATMENT_INFO.get(treatment, {})
            description = metadata.get('description', 'Unknown')
            adaptation = metadata.get('adaptation', 'Unknown')
            gene = metadata.get('gene')
            
            f.write(f"{treatment}: {description}\n")
            f.write(f"  Adaptation type: {adaptation}\n")
            f.write(f"  Gene modification: {gene if gene else 'None'}\n")
            
            if treatment in ti_tv_ratios:
                f.write(f"  Ti/Tv ratio: {ti_tv_ratios[treatment]:.2f}\n")
            
            if treatment in all_counts:
                most_common = max(all_counts[treatment].items(), key=lambda x: x[1])
                f.write(f"  Most common mutation: {most_common[0]} ({most_common[1]} occurrences)\n")
            
            f.write("\n")
        
        # Adaptation-specific analysis
        adaptation_types = set(TREATMENT_INFO.get(t, {}).get('adaptation', 'Unknown') 
                              for t in TREATMENTS)
        
        if len(adaptation_types) > 1:
            f.write("\nComparison by Adaptation Type:\n")
            f.write("----------------------------\n")
            
            for adaptation in sorted(adaptation_types):
                f.write(f"{adaptation} adaptation:\n")
                
                # Filter treatments by adaptation
                adaptation_treatments = [t for t in TREATMENTS 
                                        if TREATMENT_INFO.get(t, {}).get('adaptation') == adaptation]
                
                # Calculate average Ti/Tv ratio
                ratios = [ti_tv_ratios.get(t, 0) for t in adaptation_treatments 
                         if t in ti_tv_ratios]
                
                if ratios:
                    avg_ratio = sum(ratios) / len(ratios)
                    f.write(f"  Average Ti/Tv ratio: {avg_ratio:.2f}\n")
                
                # List treatments
                f.write(f"  Treatments: {', '.join(adaptation_treatments)}\n")
                f.write("\n")

# Run the analysis
if __name__ == "__main__":
    main()