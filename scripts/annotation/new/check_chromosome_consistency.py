#!/usr/bin/env python3
"""
check_chromosome_consistency.py - Verify chromosome names match between VCF and SnpEff
"""

import os
import gzip
import sys
import subprocess

# Configuration
VCF_DIR = "vcf/merged/filtered"
SNPEFF_DIR = "~/snpEff"

def main():
    # Get chromosome names from VCF files
    vcf_chroms = set()
    vcf_files = [f for f in os.listdir(VCF_DIR) 
                if f.endswith('.vcf.gz') or f.endswith('.vcf')]
    
    for vcf_file in vcf_files[:1]:  # Just check first file
        vcf_path = os.path.join(VCF_DIR, vcf_file)
        is_gzipped = vcf_file.endswith('.gz')
        opener = gzip.open if is_gzipped else open
        
        with opener(vcf_path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if parts:
                    vcf_chroms.add(parts[0])
                if len(vcf_chroms) >= 20:  # Get reasonable sample
                    break
    
    print(f"Found {len(vcf_chroms)} chromosome names in VCF files:")
    for chrom in sorted(vcf_chroms):
        print(f"  {chrom}")
    
    # Check SnpEff chromosomes
    print("\nChecking SnpEff database chromosomes...")
    result = subprocess.run(
        f"java -jar {SNPEFF_DIR}/snpEff.jar dump -v w303 genome | grep 'chromosome'", 
        shell=True, 
        capture_output=True, 
        text=True
    )
    
    snpeff_output = result.stdout
    print("\nSnpEff chromosome information:")
    print(snpeff_output)
    
    # Compare and report
    if result.returncode != 0:
        print("WARNING: Could not retrieve SnpEff chromosome information")
    else:
        print("\nManual verification needed:")
        print("Compare the VCF chromosome names listed above with the SnpEff chromosome information")
        print("They should match exactly for annotation to work correctly")

if __name__ == "__main__":
    main()