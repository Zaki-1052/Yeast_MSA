#!/usr/bin/env python3

# File: scripts/annotation/06_create_chromosome_mapping.py
# Purpose: Create mapping between NCBI accessions and W303 scaffold names

import os
import re
import pandas as pd
from datetime import datetime
import subprocess
import glob

def main():
    print("=== Creating Chromosome Mapping ===")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("")
    
    # Define directories
    gbk_dir = "annotation/reference/w303_scaffolds"  # Directory with GenBank files from Israel
    snpeff_dir = "/Users/zakiralibhai/snpEff/data/w303"
    mapping_dir = "annotation/chromosome_mapping"
    
    # Create output directory
    os.makedirs(mapping_dir, exist_ok=True)
    
    # Check if we have the required GenBank files
    print("Checking for GenBank files...")
    if not os.path.exists(gbk_dir):
        print(f"ERROR: GenBank directory not found: {gbk_dir}")
        print("Please create this directory and add the GenBank files from Israel's annotation pipeline")
        return
    
    # Look for files with .genbank extension
    gbk_files = glob.glob(f"{gbk_dir}/*.genbank")
    if not gbk_files:
        print(f"ERROR: No GenBank files (.genbank) found in {gbk_dir}")
        print("Please add the GenBank files from Israel's annotation pipeline")
        return
    
    print(f"Found {len(gbk_files)} GenBank files")
    print("")
    
    # Extract mapping information from GenBank files
    print("Analyzing GenBank files to extract mapping information...")
    mapping_data = []
    
    for gbk_file in gbk_files[:50]:  # Process a subset for initial analysis
        gbk_name = os.path.basename(gbk_file)
        print(f"- Analyzing {gbk_name}...")
        
        # Extract the w303 scaffold name from the filename
        w303_scaffold = gbk_name.replace('.genbank', '')
        
        # Parse the GenBank file to find the original NCBI accession
        ncbi_accession = None
        version = None
        definition = None
        
        with open(gbk_file, 'r', errors='ignore') as f:
            for line in f:
                # Look for ACCESSION line
                if line.startswith('ACCESSION'):
                    parts = line.strip().split()
                    if len(parts) > 1:
                        ncbi_accession = parts[1]
                
                # Look for VERSION line (might include accession.version)
                elif line.startswith('VERSION'):
                    version = line.strip()
                    # Extract accession.version format
                    match = re.search(r'([A-Z0-9]+\.[0-9]+)', version)
                    if match:
                        ncbi_accession = match.group(1)
                
                # Look for DEFINITION line as fallback
                elif line.startswith('DEFINITION'):
                    definition = line.strip()
                    
                # We found what we need, break the loop
                if ncbi_accession:
                    break
        
        # Add to mapping data
        mapping_data.append({
            'W303_Scaffold': w303_scaffold,
            'NCBI_Accession': ncbi_accession if ncbi_accession else "UNKNOWN",
            'Version_Info': version,
            'Definition': definition
        })
    
    # Create mapping DataFrame
    print("Creating mapping DataFrame...")
    mapping_df = pd.DataFrame(mapping_data)
    
    # Save the raw mapping
    raw_mapping_file = f"{mapping_dir}/raw_scaffold_mapping.tsv"
    mapping_df.to_csv(raw_mapping_file, sep='\t', index=False)
    print(f"Raw mapping saved to: {raw_mapping_file}")
    
    # Create the final mapping file
    final_mapping = []
    for _, row in mapping_df.iterrows():
        w303_scaffold = row['W303_Scaffold']
        ncbi_accession = row['NCBI_Accession']
        
        if ncbi_accession != "UNKNOWN":
            final_mapping.append({
                'VCF_Chromosome': ncbi_accession,
                'SnpEff_Chromosome': w303_scaffold,
                'Notes': 'Extracted from GenBank'
            })
    
    final_df = pd.DataFrame(final_mapping)
    final_mapping_file = f"{mapping_dir}/chromosome_mapping.tsv"
    final_df.to_csv(final_mapping_file, sep='\t', index=False)
    print(f"Final mapping saved to: {final_mapping_file}")
    print("")
    
    # Generate summary statistics
    print("Generating summary statistics...")
    total_mappings = len(final_mapping)
    unique_vcf_chroms = len(final_df['VCF_Chromosome'].unique())
    unique_snpeff_chroms = len(final_df['SnpEff_Chromosome'].unique())
    
    print(f"Total mappings: {total_mappings}")
    print(f"Unique VCF chromosomes: {unique_vcf_chroms}")
    print(f"Unique SnpEff chromosomes: {unique_snpeff_chroms}")
    
    # Generate report
    report_file = f"{mapping_dir}/mapping_report.txt"
    with open(report_file, 'w') as f:
        f.write("Chromosome Mapping Report\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("==============================================\n\n")
        
        f.write(f"Total GenBank files analyzed: {len(gbk_files)}\n")
        f.write(f"Total mappings established: {total_mappings}\n")
        f.write(f"Unique VCF chromosomes: {unique_vcf_chroms}\n")
        f.write(f"Unique SnpEff chromosomes: {unique_snpeff_chroms}\n\n")
        
        f.write("First 10 mappings:\n")
        for i, row in final_df.head(10).iterrows():
            f.write(f"  {row['VCF_Chromosome']} → {row['SnpEff_Chromosome']}\n")
        
        f.write("\nRecommended next steps:\n")
        f.write("1. Review the mapping file to ensure accuracy\n")
        f.write("2. Use the mapping to modify VCF files or configure SnpEff\n")
        f.write("3. Re-run the annotation with the fixed chromosome names\n")
    
    print(f"Report saved to: {report_file}")
    print("")
    
    print("=== Chromosome Mapping Creation complete ===")

if __name__ == "__main__":
    main()
