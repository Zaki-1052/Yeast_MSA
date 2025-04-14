#!/usr/bin/env python3

# File: scripts/annotation/29_locate_and_analyze_vcfs.py
# Purpose: Find all VCF files and search for target genes

import os
import subprocess
import pandas as pd
from datetime import datetime
import re

def main():
    print("=== Locating and Analyzing VCF Files ===")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("")
    
    # Define target genes
    gene_mapping = [
        {"SGD_Gene": "YHR190W", "W303_ID": "W303_0BY00190", "Gene": "W3030BY00190", 
         "Description": "similar to ERG9"},
        {"SGD_Gene": "YGR175C", "W303_ID": "W303_0AJ00440", "Gene": "W3030AJ00440", 
         "Description": "similar to ERG1"},
        {"SGD_Gene": "YHR072W", "W303_ID": "W303_0CB00150", "Gene": "W3030CB00150", 
         "Description": "similar to ERG7"},
        {"SGD_Gene": "YHR007C", "W303_ID": "W303_0EI00110", "Gene": "W3030EI00110", 
         "Description": "similar to ERG11"},
        {"SGD_Gene": "YNL280C", "W303_ID": "W303_0O00140", "Gene": "W3030O00140", 
         "Description": "similar to ERG24"},
        {"SGD_Gene": "YGR060W", "W303_ID": "W303_0W00270", "Gene": "W3030W00270", 
         "Description": "similar to ERG25"},
        {"SGD_Gene": "YML008C", "W303_ID": "W303_0S00270", "Gene": "W3030S00270", 
         "Description": "similar to ERG6"},
        {"SGD_Gene": "YMR202W", "W303_ID": "W303_0AD00260", "Gene": "W3030AD00260", 
         "Description": "similar to ERG2"},
        {"SGD_Gene": "YLR056W", "W303_ID": "W303_0E01010", "Gene": "W3030E01010", 
         "Description": "similar to ERG3"},
        {"SGD_Gene": "YMR015C", "W303_ID": "W303_0S00490", "Gene": "W3030S00490", 
         "Description": "similar to ERG5"},
        {"SGD_Gene": "YGL012W", "W303_ID": "W303_0Y00390", "Gene": "W3030Y00390", 
         "Description": "similar to ERG4"}
    ]
    
    # Define output directory
    output_dir = "annotation/gene_results_final"
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(f"{output_dir}/by_gene", exist_ok=True)
    
    # Step 1: Find all VCF files in the project
    print("Step 1: Locating all VCF files in the project...")
    cmd = "find . -name '*.vcf.gz'"
    try:
        all_vcfs = subprocess.check_output(cmd, shell=True, text=True).splitlines()
        print(f"Found {len(all_vcfs)} VCF files")
        print("First 5 files:")
        for i, vcf in enumerate(all_vcfs[:5]):
            print(f"  {i+1}. {vcf}")
        
        # Save to file
        with open(f"{output_dir}/all_vcf_files.txt", 'w') as f:
            for vcf in all_vcfs:
                f.write(f"{vcf}\n")
    except subprocess.CalledProcessError as e:
        print(f"Error finding VCF files: {e}")
        all_vcfs = []
    
    if not all_vcfs:
        print("No VCF files found. Exiting.")
        return
    
    # Step 2: Test access to VCF files
    print("\nStep 2: Testing access to VCF files...")
    
    # Test with different commands
    test_vcf = all_vcfs[0]
    print(f"Testing access to: {test_vcf}")
    
    methods = [
        {"name": "zcat", "cmd": f"zcat {test_vcf} | head -5"},
        {"name": "gunzip -c", "cmd": f"gunzip -c {test_vcf} | head -5"},
        {"name": "bcftools", "cmd": f"bcftools view -h {test_vcf} | head -5"},
        {"name": "file", "cmd": f"file {test_vcf}"}
    ]
    
    working_method = None
    
    for method in methods:
        try:
            print(f"  Testing {method['name']}...")
            result = subprocess.check_output(method['cmd'], shell=True, text=True)
            print(f"  ✓ {method['name']} works!")
            print(f"  Sample output: {result.splitlines()[0]}")
            working_method = method
            break
        except subprocess.CalledProcessError as e:
            print(f"  ✗ {method['name']} failed: {e}")
    
    if not working_method:
        print("Could not find a working method to access VCF files. Exiting.")
        return
    
    print(f"Will use {working_method['name']} for accessing VCF files")
    
    # Step 3: Analyze VCF files for target genes
    # Step 3: Create BED file for target genes (Placeholder - Requires correct scaffold mapping)
    # Skipping accurate BED creation for now as the W303_ID to scaffold mapping needs verification.
    # Instead, we will proceed directly to searching the VCFs.
    print("\nStep 3: Skipping BED file creation (requires scaffold mapping verification).")


    # Step 4: Search VCFs directly for target gene identifiers
    print("\nStep 4: Searching VCFs for target gene identifiers...")

    all_found_variants = []
    access_cmd = ""
    if working_method['name'] == 'zcat':
        access_cmd = "zcat"
    elif working_method['name'] == 'gunzip -c':
        access_cmd = "gunzip -c"
    elif working_method['name'] == 'bcftools':
        # bcftools view requires specific syntax for streaming, zcat/gunzip easier with grep
        print("  Using zcat as fallback for bcftools with grep.")
        access_cmd = "zcat"
    else:
        access_cmd = "zcat" # Default fallback

    print(f"Using command '{access_cmd}' to access VCF content.")

    for vcf_file in all_vcfs:
        # Extract a clean sample name
        sample_basename = os.path.basename(vcf_file)
        sample = sample_basename
        common_suffixes = [
            '.sorted.vcf.gz', '.fixed.vcf.gz', '.renamed.vcf.gz',
            '.snpeff.vcf.gz', '.vcf.gz'
        ]
        for suffix in common_suffixes:
            if sample.endswith(suffix):
                sample = sample[:-len(suffix)]
                break

        print(f"  Processing {sample} ({vcf_file})...")

        found_in_sample_count = 0
        gene_matches_in_sample = {}

        for gene_info in gene_mapping:
            # Create a comprehensive search pattern including all known identifiers
            search_terms = list(filter(None, [
                gene_info['SGD_Gene'],          # e.g., YHR190W
                gene_info['W303_ID'],           # e.g., W303_0BY00190
                gene_info['Gene'],              # e.g., W3030BY00190
                gene_info['Description'].split()[-1] # e.g., ERG9
            ]))
            # Escape special characters for grep -E
            escaped_terms = [re.escape(term) for term in search_terms]
            search_pattern = '|'.join(escaped_terms)

            if not search_pattern:
                print(f"    Skipping search for {gene_info['SGD_Gene']} - no valid identifiers.")
                continue

            # Use grep -E to search within the decompressed VCF content
            # Limit matches per file to avoid excessive output for common terms if any
            cmd = f"{access_cmd} {vcf_file} | grep -v '^#' | grep -E '{search_pattern}' | head -50"
            try:
                result = subprocess.check_output(cmd, shell=True, text=True, stderr=subprocess.PIPE)
                if result:
                    lines = result.strip().split('\n')
                    gene_match_count = 0
                    for line in lines:
                        if not line: continue
                        fields = line.split('\t')
                        if len(fields) >= 8:
                            chrom = fields[0]
                            pos = fields[1]
                            ref = fields[3]
                            alt = fields[4]
                            info = fields[7]
                            all_found_variants.append({
                                'Sample': sample,
                                'SGD_Gene': gene_info['SGD_Gene'],
                                'VCF_File': vcf_file,
                                'Chromosome': chrom,
                                'Position': pos,
                                'Ref': ref,
                                'Alt': alt,
                                'Info': info
                            })
                            found_in_sample_count += 1
                            gene_match_count +=1
                    if gene_match_count > 0:
                        gene_matches_in_sample[gene_info['SGD_Gene']] = gene_match_count
                        print(f"    ✓ Found {gene_match_count} potential variant lines for {gene_info['SGD_Gene']}")

            except subprocess.CalledProcessError:
                # No matches found for this gene in this file
                pass
            except Exception as e:
                print(f"      Error processing {gene_info['SGD_Gene']} in {sample}: {e}")

        if found_in_sample_count > 0:
            print(f"    Found {found_in_sample_count} total potential variant lines in {sample}")
        else:
                print(f"    No potential variants found for target genes in {sample}")


    # Save all found variants
    if all_found_variants:
        variants_df = pd.DataFrame(all_found_variants)
        variants_file = f"{output_dir}/potential_target_variants.tsv"
        variants_df.to_csv(variants_file, sep='\t', index=False)
        print(f"\nSaved {len(all_found_variants)} potential variant lines to {variants_file}")
        print("Please review this file to confirm if these variants are within the target genes.")
    else:
        print("\nNo potential variants found containing target gene identifiers across all files.")


    print("\n=== VCF Analysis Complete ===")

if __name__ == "__main__":
    main()