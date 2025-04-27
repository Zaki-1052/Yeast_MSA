#!/usr/bin/env python3
"""
validate_chromosome_mapping.py

This script validates that our chromosome mapping covers all chromosomes in the VCF files.

Usage:
    python validate_chromosome_mapping.py --vcf_dir <vcf_directory> --mapping_file <mapping_file>
"""

import os
import argparse
import pandas as pd
import gzip
import re

def extract_chromosomes_from_vcf(vcf_dir):
    """
    Extract all chromosome identifiers from VCF files.
    
    Args:
        vcf_dir (str): Directory containing VCF files
        
    Returns:
        set: Set of unique chromosome identifiers
    """
    chromosomes = set()
    
    # List VCF files
    vcf_files = [f for f in os.listdir(vcf_dir) if f.endswith('.vcf.gz') or f.endswith('.vcf')]
    
    if not vcf_files:
        print(f"No VCF files found in {vcf_dir}")
        return chromosomes
    
    print(f"Found {len(vcf_files)} VCF files to process")
    
    # Check first 10 lines of each file for chromosomes
    for vcf_file in vcf_files:
        file_path = os.path.join(vcf_dir, vcf_file)
        print(f"Examining {vcf_file}...")
        
        try:
            # Open gzipped or regular VCF file
            opener = gzip.open if vcf_file.endswith('.gz') else open
            mode = 'rt' if vcf_file.endswith('.gz') else 'r'
            
            with opener(file_path, mode) as vcf:
                # Check header lines for contig information
                for line in vcf:
                    if line.startswith('#'):
                        # Look for contig header lines
                        if line.startswith('##contig='):
                            match = re.search(r'ID=([^,]+)', line)
                            if match:
                                chrom = match.group(1)
                                chromosomes.add(chrom)
                    else:
                        # Process first few data lines to find chromosomes
                        chrom = line.split('\t')[0]
                        chromosomes.add(chrom)
                        
                        # Limit to checking a few data lines
                        if len(chromosomes) >= 30:
                            break
        
        except Exception as e:
            print(f"  ERROR examining {vcf_file}: {str(e)}")
    
    print(f"Found {len(chromosomes)} unique chromosomes in VCF files")
    return chromosomes

def validate_mapping(vcf_chromosomes, mapping_file):
    """
    Validate that all VCF chromosomes are in the mapping.
    
    Args:
        vcf_chromosomes (set): Set of chromosome identifiers from VCF files
        mapping_file (str): Path to chromosome mapping file
        
    Returns:
        tuple: (missing chromosomes, mapping coverage percentage)
    """
    # Load mapping file
    mapping_df = pd.read_csv(mapping_file, sep='\t')
    
    # Get list of chromosomes in mapping
    mapped_chromosomes = set(mapping_df['cm_identifier'])
    
    # Find missing chromosomes
    missing = vcf_chromosomes - mapped_chromosomes
    
    # Calculate coverage
    coverage_pct = 100 * (len(vcf_chromosomes) - len(missing)) / len(vcf_chromosomes) if vcf_chromosomes else 0
    
    return missing, coverage_pct

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Validate chromosome mapping coverage")
    parser.add_argument("--vcf_dir", required=True, help="Directory containing VCF files")
    parser.add_argument("--mapping_file", required=True, help="Chromosome mapping file (tsv)")
    args = parser.parse_args()
    
    # Extract chromosomes from VCF files
    vcf_chromosomes = extract_chromosomes_from_vcf(args.vcf_dir)
    
    # Validate mapping
    missing, coverage_pct = validate_mapping(vcf_chromosomes, args.mapping_file)
    
    # Print results
    print("\nValidation Results:")
    print(f"  Total chromosomes in VCF files: {len(vcf_chromosomes)}")
    print(f"  Mapping coverage: {coverage_pct:.2f}%")
    
    if missing:
        print("\nWARNING: The following chromosomes are in VCF files but missing from mapping:")
        for chrom in sorted(missing):
            print(f"  - {chrom}")
    else:
        print("\nSuccess: All chromosomes in VCF files are covered by the mapping")

if __name__ == "__main__":
    main()