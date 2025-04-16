#!/usr/bin/env python3

# File: scripts/annotation/26_extract_target_genes_simple.py
# Purpose: Simple extraction of target gene references from GenBank file

import os
import subprocess
import re
from datetime import datetime

def main():
    print("=== Simple Target Gene Extraction ===")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("")
    
    # Define directories
    snpeff_dir = "/Users/zakiralibhai/snpEff"
    gene_dir = "annotation/genes_of_interest"
    output_dir = "annotation/gene_results_simple"
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Read target genes
    target_genes_file = f"{gene_dir}/target_genes.txt"
    with open(target_genes_file, 'r') as f:
        target_genes = [line.strip() for line in f if line.strip()]
    
    print(f"Looking for {len(target_genes)} target genes")
    print("")
    
    # Step 1: Extract references to target genes from the GenBank file
    print("Step 1: Extracting references to target genes...")
    
    # Save all output to a single file for analysis
    output_file = f"{output_dir}/gene_references.txt"
    with open(output_file, 'w') as output:
        output.write(f"Target Gene Extraction Results\n")
        output.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        output.write("==============================================\n\n")
        
        # Loop through each target gene
        for gene in target_genes:
            print(f"Searching for {gene}...")
            output.write(f"Gene: {gene}\n")
            output.write("-----------------------\n")
            
            # Try multiple search patterns
            found = False
            
            # Pattern 1: Exact match
            try:
                cmd = f"cd {snpeff_dir} && grep -n -A 5 -B 5 '\\b{gene}\\b' data/w303/genes.gbk | head -30"
                result = subprocess.check_output(cmd, shell=True, text=True)
                if result:
                    print(f"  ✓ Found exact match for {gene}")
                    output.write("Exact match:\n")
                    output.write(result)
                    output.write("\n")
                    found = True
            except subprocess.CalledProcessError:
                pass
            
            # Pattern 2: Without Y prefix
            if gene.startswith('Y'):
                gene_no_y = gene[1:]
                try:
                    cmd = f"cd {snpeff_dir} && grep -n -A 5 -B 5 '\\b{gene_no_y}\\b' data/w303/genes.gbk | head -30"
                    result = subprocess.check_output(cmd, shell=True, text=True)
                    if result:
                        print(f"  ✓ Found match without Y prefix for {gene}")
                        output.write("Match without Y prefix:\n")
                        output.write(result)
                        output.write("\n")
                        found = True
                except subprocess.CalledProcessError:
                    pass
            
            # Pattern 3: Case insensitive
            try:
                cmd = f"cd {snpeff_dir} && grep -n -i -A 5 -B 5 '{gene}' data/w303/genes.gbk | head -30"
                result = subprocess.check_output(cmd, shell=True, text=True)
                if result:
                    print(f"  ✓ Found case insensitive match for {gene}")
                    output.write("Case insensitive match:\n")
                    output.write(result)
                    output.write("\n")
                    found = True
            except subprocess.CalledProcessError:
                pass
            
            if not found:
                print(f"  ✗ No matches found for {gene}")
                output.write("No matches found\n\n")
            
            output.write("\n-----------------------------------------\n\n")
    
    print(f"\nAll target gene references saved to: {output_file}")
    
    # Step 2: Extract complete gene entries to understand format
    print("\nStep 2: Extracting complete gene entries for reference...")
    
    gene_entries_file = f"{output_dir}/sample_gene_entries.txt"
    with open(gene_entries_file, 'w') as f:
        f.write("Sample Gene Entries from SnpEff Database\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("==============================================\n\n")
        
        # Extract a complete gene entry
        try:
            cmd = f"cd {snpeff_dir} && grep -A 50 'gene ' data/w303/genes.gbk | head -100"
            result = subprocess.check_output(cmd, shell=True, text=True)
            f.write("Sample gene entry:\n")
            f.write(result)
            f.write("\n")
        except subprocess.CalledProcessError:
            f.write("Error extracting sample gene entry\n")
    
    print(f"Sample gene entries saved to: {gene_entries_file}")
    
    # Step 3: Get statistics on the genes.gbk file
    print("\nStep 3: Analyzing genes.gbk file...")
    
    stats_file = f"{output_dir}/gbk_stats.txt"
    with open(stats_file, 'w') as f:
        f.write("GenBank File Statistics\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("==============================================\n\n")
        
        # Get file size
        try:
            cmd = f"ls -lh {snpeff_dir}/data/w303/genes.gbk"
            result = subprocess.check_output(cmd, shell=True, text=True)
            f.write("File size information:\n")
            f.write(result)
            f.write("\n")
        except subprocess.CalledProcessError:
            f.write("Error getting file size\n")
        
        # Count genes
        try:
            cmd = f"cd {snpeff_dir} && grep -c 'gene' data/w303/genes.gbk"
            result = subprocess.check_output(cmd, shell=True, text=True)
            f.write(f"Number of 'gene' occurrences: {result}")
            f.write("\n")
        except subprocess.CalledProcessError:
            f.write("Error counting genes\n")
        
        # Count SGD gene references
        try:
            cmd = f"cd {snpeff_dir} && grep -c 'Y[A-Z][A-Z][0-9][0-9][0-9][WC]' data/w303/genes.gbk"
            result = subprocess.check_output(cmd, shell=True, text=True)
            f.write(f"Number of SGD gene references: {result}")
            f.write("\n")
        except subprocess.CalledProcessError:
            f.write("Error counting SGD gene references\n")
    
    print(f"GenBank file statistics saved to: {stats_file}")
    
    # Step 4: Extract SnpEff database information
    print("\nStep 4: Extracting SnpEff database information...")
    
    db_info_file = f"{output_dir}/snpeff_db_info.txt"
    with open(db_info_file, 'w') as f:
        f.write("SnpEff Database Information\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("==============================================\n\n")
        
        # Run SnpEff dump
        try:
            cmd = f"cd {snpeff_dir} && java -jar snpEff.jar dump w303 | head -100"
            result = subprocess.check_output(cmd, shell=True, text=True)
            f.write("SnpEff database dump:\n")
            f.write(result)
            f.write("\n")
        except subprocess.CalledProcessError:
            f.write("Error running SnpEff dump\n")
    
    print(f"SnpEff database information saved to: {db_info_file}")
    
    print("\n=== Simple Target Gene Extraction Complete ===")
    print("\nAll results have been saved to the following files:")
    print(f"- Target gene references: {output_file}")
    print(f"- Sample gene entries: {gene_entries_file}")
    print(f"- GenBank file statistics: {stats_file}")
    print(f"- SnpEff database information: {db_info_file}")
    print("\nPlease examine these files to understand the gene format and find locations")

if __name__ == "__main__":
    main()
