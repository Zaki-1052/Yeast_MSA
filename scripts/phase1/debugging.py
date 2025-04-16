#!/usr/bin/env python3

import os
import pandas as pd
import subprocess

def inspect_mutation_file(filename):
    """Examine the contents of a mutation file to identify issues."""
    print(f"Inspecting file: {filename}")
    
    if not os.path.exists(filename):
        print(f"  File does not exist")
        return None
    
    try:
        # Try to read the file with tab separator
        data = pd.read_csv(filename, sep='\t', header=None, 
                          names=['CHROM', 'POS', 'REF', 'ALT'])
        
        print(f"  File shape: {data.shape}")
        print(f"  First few rows:")
        print(data.head())
        
        # Check REF/ALT columns for valid nucleotides
        ref_valid = data['REF'].str.len().eq(1).all() and data['REF'].str.upper().isin(['A', 'C', 'G', 'T']).all()
        alt_valid = data['ALT'].str.len().eq(1).all() and data['ALT'].str.upper().isin(['A', 'C', 'G', 'T']).all()
        
        print(f"  REF column contains valid nucleotides: {ref_valid}")
        print(f"  ALT column contains valid nucleotides: {alt_valid}")
        
        if not ref_valid or not alt_valid:
            print(f"  Sample REF values: {data['REF'].head().tolist()}")
            print(f"  Sample ALT values: {data['ALT'].head().tolist()}")
        
        return data
    except Exception as e:
        print(f"  Error reading file: {e}")
        return None

def find_vcf_file(treatment):
    """Find a valid VCF file for a treatment."""
    # Define possible VCF file patterns
    vcf_patterns = [
        f"results/merged/analysis/{treatment}/highconf.vcf.gz",
        f"results/merged/analysis/{treatment}_highconf.vcf.gz",
        f"results/merged/analysis/{treatment}/specific.vcf.gz",
        f"results/merged/analysis/{treatment}_specific.vcf.gz"
    ]
    
    # Handle WT-37 backward compatibility
    if treatment == 'WT-37':
        vcf_patterns.extend([
            "results/merged/analysis/WT/highconf.vcf.gz",
            "results/merged/analysis/WT_highconf.vcf.gz",
            "results/merged/analysis/WT/specific.vcf.gz",
            "results/merged/analysis/WT_specific.vcf.gz"
        ])
    
    # Try each pattern
    for pattern in vcf_patterns:
        if os.path.exists(pattern):
            print(f"Found VCF file for {treatment}: {pattern}")
            return pattern
    
    print(f"No VCF file found for {treatment}")
    return None

def extract_from_vcf(treatment, vcf_file, output_path):
    """Extract mutation data from a VCF file and save to output file."""
    if not vcf_file or not os.path.exists(vcf_file):
        print(f"VCF file not found: {vcf_file}")
        return False
    
    try:
        print(f"Extracting mutations from {vcf_file}")
        # Extract data using bcftools
        cmd = f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' {vcf_file}"
        output = subprocess.check_output(cmd, shell=True).decode('utf-8')
        
        # Parse the output
        rows = []
        for line in output.strip().split('\n'):
            if line:  # Skip empty lines
                parts = line.split('\t')
                # Ensure we have all required parts
                if len(parts) == 4:
                    # Validate that REF and ALT are valid nucleotides
                    ref, alt = parts[2], parts[3]
                    if (len(ref) == 1 and ref.upper() in ['A', 'C', 'G', 'T'] and
                        len(alt) == 1 and alt.upper() in ['A', 'C', 'G', 'T']):
                        rows.append(parts)
        
        if not rows:
            print(f"No valid mutations found in {vcf_file}")
            return False
        
        # Write to output file
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'w') as f:
            for row in rows:
                f.write('\t'.join(row) + '\n')
        
        print(f"Successfully extracted {len(rows)} mutations for {treatment}")
        return True
    except Exception as e:
        print(f"Error extracting from {vcf_file}: {e}")
        return False

def regenerate_all_mutation_files():
    """Regenerate all mutation files from VCF sources."""
    treatments = ['WT-37', 'WTA', 'STC', 'CAS']
    
    # First, check current files
    print("Checking current mutation files...")
    for treatment in treatments:
        file_path = f"mutation_spectrum_analysis/{treatment}_mutations.txt"
        inspect_mutation_file(file_path)
    
    # Now regenerate each file
    print("\nRegenerating mutation files...")
    for treatment in treatments:
        output_path = f"mutation_spectrum_analysis/{treatment}_mutations.txt"
        
        # Remove corrupted file if it exists
        if os.path.exists(output_path):
            print(f"Removing existing file: {output_path}")
            os.remove(output_path)
        
        # Find a valid VCF file
        vcf_file = find_vcf_file(treatment)
        
        if vcf_file:
            # Extract mutations
            success = extract_from_vcf(treatment, vcf_file, output_path)
            
            if success:
                # Verify the regenerated file
                print(f"Verifying regenerated file for {treatment}:")
                inspect_mutation_file(output_path)
            else:
                print(f"Failed to regenerate mutation file for {treatment}")
        else:
            print(f"Could not find a VCF file for {treatment}")
    
    print("\nMutation file regeneration complete!")

# Run the regeneration
if __name__ == "__main__":
    regenerate_all_mutation_files()