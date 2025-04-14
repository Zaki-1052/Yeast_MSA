#!/usr/bin/env python3

# File: scripts/annotation/05_fix_chromosome_mapping.py
# Purpose: Fix chromosome naming mismatches between VCF files and SnpEff database

import os
import gzip
import re
import pandas as pd
from datetime import datetime
import subprocess
import shutil

def main():
    print("=== Fixing Chromosome Mapping for SnpEff Annotation ===")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("")
    
    # Define directories
    vcf_source = "annotation/vcf_ready"
    snpeff_dir = "/Users/zakiralibhai/snpEff"
    data_dir = os.path.join(snpeff_dir, "data", "w303")
    fixed_dir = "annotation/vcf_fixed_chrom"
    debug_dir = "annotation/debug"
    
    # Create output directories
    os.makedirs(fixed_dir, exist_ok=True)
    os.makedirs(debug_dir, exist_ok=True)
    
    # Step 1: Extract chromosome names from VCF files
    print("Extracting chromosome names from VCF files...")
    vcf_chroms = set()
    vcf_files = []
    
    for root, dirs, files in os.walk(vcf_source):
        for file in files:
            if file.endswith(".sorted.vcf.gz"):
                vcf_files.append(os.path.join(root, file))
                
                # Extract chromosome names
                sample = file.replace(".sorted.vcf.gz", "")
                print(f"- Examining {sample}...")
                
                # Use bcftools to list chromosomes
                try:
                    output = subprocess.check_output(
                        f"bcftools view -h {os.path.join(root, file)} | grep '##contig='",
                        shell=True, text=True
                    )
                    
                    # Extract contig IDs
                    for line in output.splitlines():
                        match = re.search(r'ID=([^,]+)', line)
                        if match:
                            vcf_chroms.add(match.group(1))
                except subprocess.CalledProcessError:
                    print(f"  Error extracting contigs from {file}")
    
    vcf_chrom_list = sorted(list(vcf_chroms))
    print(f"Found {len(vcf_chrom_list)} unique chromosomes in VCF files")
    
    # Save to file
    with open(f"{debug_dir}/vcf_chromosomes.txt", 'w') as f:
        for chrom in vcf_chrom_list:
            f.write(f"{chrom}\n")
    
    print(f"First 5 chromosome names: {', '.join(vcf_chrom_list[:5])}")
    print("")
    
    # Step 2: Extract chromosome names from SnpEff database
    print("Extracting chromosome names from SnpEff database...")
    snpeff_chroms = set()
    
    # Method 1: Check genome file
    genome_file = os.path.join(data_dir, "sequences.fa")
    if os.path.exists(genome_file):
        print(f"- Examining {genome_file}...")
        try:
            output = subprocess.check_output(
                f"grep '>' {genome_file} | head -10",
                shell=True, text=True
            )
            
            for line in output.splitlines():
                if line.startswith('>'):
                    chrom = line[1:].split()[0]  # Remove '>' and take first word
                    snpeff_chroms.add(chrom)
                    
            print(f"  Found chromosome names in genome file")
        except subprocess.CalledProcessError:
            print(f"  Error examining genome file")
    
    # Method 2: Check genes file
    genes_file = os.path.join(data_dir, "genes.gbk")
    if os.path.exists(genes_file):
        print(f"- Examining {genes_file}...")
        try:
            # Use grep to extract LOCUS lines which contain sequence IDs
            output = subprocess.check_output(
                f"grep 'LOCUS' {genes_file} | head -10",
                shell=True, text=True
            )
            
            for line in output.splitlines():
                if line.startswith('LOCUS'):
                    parts = line.split()
                    if len(parts) > 1:
                        chrom = parts[1]
                        snpeff_chroms.add(chrom)
                        
            print(f"  Found chromosome names in genes file")
        except subprocess.CalledProcessError:
            print(f"  Error examining genes file")
    
    snpeff_chrom_list = sorted(list(snpeff_chroms))
    print(f"Found {len(snpeff_chrom_list)} unique chromosomes in SnpEff database")
    
    # Save to file
    with open(f"{debug_dir}/snpeff_chromosomes.txt", 'w') as f:
        for chrom in snpeff_chrom_list:
            f.write(f"{chrom}\n")
    
    if snpeff_chrom_list:
        print(f"First 5 chromosome names: {', '.join(snpeff_chrom_list[:5])}")
    else:
        print("No chromosome names found in SnpEff database files")
    print("")
    
    # Step 3: Create a mapping between VCF and SnpEff chromosomes
    print("Creating chromosome mapping...")
    
    # First attempt: Check if there's a predictable pattern
    if vcf_chrom_list and snpeff_chrom_list:
        print("Analyzing chromosome naming patterns...")
        
        # Count chromosomes with JRIU prefix in VCF
        jriu_count = sum(1 for c in vcf_chrom_list if c.startswith("JRIU"))
        if jriu_count > 0:
            print(f"Found {jriu_count} chromosomes with JRIU prefix in VCF files")
            
        # Check for chromosome naming patterns in SnpEff
        roman_count = sum(1 for c in snpeff_chrom_list if c in 
                         ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"])
        
        if roman_count > 0:
            print(f"Found {roman_count} chromosomes with Roman numerals in SnpEff database")
            
        # Check if chromosomes in SnpEff are like 'chrI', 'chrII', etc.
        chr_count = sum(1 for c in snpeff_chrom_list if c.startswith("chr"))
        if chr_count > 0:
            print(f"Found {chr_count} chromosomes with 'chr' prefix in SnpEff database")
    
    # For now, let's assume the VCF scaffold IDs need to be mapped to standard yeast chromosomes
    # We'll create a chromosome mapping file as a placeholder for now
    print("Creating a placeholder chromosome mapping file...")
    placeholder_mapping = []
    
    # Add a few example mappings (this will be incomplete and needs to be manually reviewed)
    if vcf_chrom_list and len(vcf_chrom_list) > 0:
        for i, chrom in enumerate(vcf_chrom_list[:min(16, len(vcf_chrom_list))]):
            if i < len(snpeff_chrom_list):
                placeholder_mapping.append({
                    "VCF_Chromosome": chrom,
                    "SnpEff_Chromosome": snpeff_chrom_list[i],
                    "Notes": "Auto-generated mapping (NEEDS VERIFICATION)"
                })
    
    mapping_df = pd.DataFrame(placeholder_mapping)
    mapping_file = f"{debug_dir}/chromosome_mapping.tsv"
    mapping_df.to_csv(mapping_file, sep='\t', index=False)
    
    print(f"Created placeholder mapping file: {mapping_file}")
    print("⚠️  IMPORTANT: This mapping file is just a starting point!")
    print("⚠️  You MUST manually verify and correct the mappings!")
    print("")
    
    # Step 4: Generate report on chromosome issues
    print("Generating summary report...")
    
    with open(f"{debug_dir}/chromosome_report.txt", 'w') as report:
        report.write("Chromosome Mapping Analysis Report\n")
        report.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        report.write("==============================================\n\n")
        
        report.write("VCF Chromosomes:\n")
        for chrom in vcf_chrom_list[:20]:
            report.write(f"- {chrom}\n")
        if len(vcf_chrom_list) > 20:
            report.write(f"... and {len(vcf_chrom_list) - 20} more\n")
        report.write("\n")
        
        report.write("SnpEff Database Chromosomes:\n")
        for chrom in snpeff_chrom_list[:20]:
            report.write(f"- {chrom}\n")
        if len(snpeff_chrom_list) > 20:
            report.write(f"... and {len(snpeff_chrom_list) - 20} more\n")
        report.write("\n")
        
        report.write("Observations:\n")
        report.write("1. The chromosome names in VCF files don't match those in the SnpEff database\n")
        report.write("2. This mismatch prevents SnpEff from properly annotating variants\n")
        report.write("3. To fix this, we need to either:\n")
        report.write("   a. Create a chromosome mapping file for SnpEff\n")
        report.write("   b. Rename the chromosomes in the VCF files\n")
        report.write("\n")
        
        report.write("Recommended Actions:\n")
        report.write("1. Examine the genes.gbk file to understand the exact chromosome naming convention\n")
        report.write("2. Examine the GenBank files from Israel's annotation pipeline\n")
        report.write("3. Create a complete and accurate chromosome mapping file\n")
        report.write("4. Either configure SnpEff to use this mapping or modify the VCF files\n")
    
    print(f"Summary report saved to: {debug_dir}/chromosome_report.txt")
    print("")
    
    print("=== Chromosome Mapping Analysis complete ===")

if __name__ == "__main__":
    main()
