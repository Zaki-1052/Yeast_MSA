#!/usr/bin/env python3
"""
vcf_preparation.py - Prepare VCF files for SnpEff annotation
"""

import os
import gzip
import sys
import subprocess
from collections import defaultdict

# Configuration
VCF_DIR = "vcf/merged/filtered"
OUTPUT_DIR = "annotation/vcf_prepared"
REPORT_FILE = "annotation/vcf_preparation_report.txt"

def main():
    """Main function to prepare VCF files for annotation"""
    
    # Create output directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Get list of VCF files
    vcf_files = [f for f in os.listdir(VCF_DIR) 
                if f.endswith('.vcf.gz') or f.endswith('.vcf')]
    
    if not vcf_files:
        print(f"No VCF files found in {VCF_DIR}")
        sys.exit(1)

    print(f"Found {len(vcf_files)} VCF files for processing")
    
    # Initialize report
    with open(REPORT_FILE, 'w') as report:
        report.write("VCF Preparation Report\n")
        report.write("=====================\n\n")
        report.write(f"Total VCF files found: {len(vcf_files)}\n\n")
    
    # Process each file
    chrom_counts = defaultdict(int)
    sample_names = set()
    
    for vcf_file in vcf_files:
        vcf_path = os.path.join(VCF_DIR, vcf_file)
        print(f"Analyzing {vcf_file}...")
        
        # Check if file is gzipped
        is_gzipped = vcf_file.endswith('.gz')
        opener = gzip.open if is_gzipped else open
        
        # Extract information from VCF
        with opener(vcf_path, 'rt') as f:
            # Extract sample info and chromosome names
            sample_name = None
            chromosomes = set()
            
            for line in f:
                # Skip empty lines
                if not line.strip():
                    continue
                
                # Parse header lines
                if line.startswith('#'):
                    if line.startswith('#CHROM'):
                        # Get sample name from header
                        parts = line.strip().split('\t')
                        if len(parts) > 9:  # Has at least one sample
                            sample_name = parts[9]
                            sample_names.add(sample_name)
                    continue
                
                # Parse variant lines
                parts = line.strip().split('\t')
                if len(parts) > 0:
                    chrom = parts[0]
                    chromosomes.add(chrom)
                    chrom_counts[chrom] += 1
                
                # Only read a few lines to get chromosome info
                if len(chromosomes) >= 5:
                    break
            
            # Write file info to report
            with open(REPORT_FILE, 'a') as report:
                report.write(f"File: {vcf_file}\n")
                report.write(f"  Sample name: {sample_name}\n")
                report.write(f"  Chromosome naming examples: {', '.join(list(chromosomes)[:5])}\n")
                report.write("\n")
    
    # Summarize chromosome naming
    with open(REPORT_FILE, 'a') as report:
        report.write("\nChromosome Naming Summary\n")
        report.write("========================\n")
        for chrom, count in sorted(chrom_counts.items()):
            report.write(f"{chrom}: {count} occurrences\n")
        
        report.write("\nSample Names Summary\n")
        report.write("===================\n")
        for sample in sorted(sample_names):
            report.write(f"{sample}\n")
    
    print(f"Analysis complete. See {REPORT_FILE} for details.")

if __name__ == "__main__":
    main()