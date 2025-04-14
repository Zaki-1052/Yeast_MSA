#!/usr/bin/env python3

# File: scripts/annotation/20_debug_vcf_annotations.py
# Purpose: Debug VCF file access and gene annotations

import os
import subprocess
import glob
import pandas as pd
from datetime import datetime

def main():
    print("=== Debugging VCF File Access and Annotations ===")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("")
    
    # Define directories
    results_dir = "annotation/results"
    debug_dir = "annotation/debug_annotations"
    mapped_dir = "annotation/mapped_genes"
    
    # Ensure directories exist
    os.makedirs(debug_dir, exist_ok=True)
    os.makedirs(mapped_dir, exist_ok=True)
    
    # Step 1: Check VCF file accessibility
    print("Step 1: Checking VCF file accessibility...")
    
    # List files in results directory
    if os.path.exists(results_dir):
        files = os.listdir(results_dir)
        vcf_files = [f for f in files if f.endswith('.vcf.gz')]
        print(f"Found {len(vcf_files)} .vcf.gz files in {results_dir}")
        
        if vcf_files:
            sample_file = os.path.join(results_dir, vcf_files[0])
            print(f"Sample file: {sample_file}")
            
            # Check if file is accessible with zcat
            try:
                print("Testing zcat access...")
                zcat_output = subprocess.check_output(f"zcat {sample_file} | head -5", shell=True, text=True)
                print("zcat successful, first 5 lines:")
                print(zcat_output)
            except subprocess.CalledProcessError as e:
                print(f"Error with zcat: {e}")
                print("Trying gunzip -c instead...")
                try:
                    gunzip_output = subprocess.check_output(f"gunzip -c {sample_file} | head -5", shell=True, text=True)
                    print("gunzip successful, first 5 lines:")
                    print(gunzip_output)
                except subprocess.CalledProcessError as e:
                    print(f"Error with gunzip: {e}")
                    print("Trying other options...")
                    try:
                        # List file details
                        file_info = subprocess.check_output(f"ls -la {sample_file}", shell=True, text=True)
                        print(f"File info: {file_info}")
                        
                        # Try file command
                        file_type = subprocess.check_output(f"file {sample_file}", shell=True, text=True)
                        print(f"File type: {file_type}")
                        
                        # Try bcftools
                        bcftools_output = subprocess.check_output(f"bcftools view -h {sample_file} | head -5", shell=True, text=True)
                        print("bcftools successful, first 5 header lines:")
                        print(bcftools_output)
                    except subprocess.CalledProcessError as e:
                        print(f"Other file access methods failed: {e}")
    else:
        print(f"Directory not found: {results_dir}")
    
    # Step 2: Check for our mapped genes in annotations
    print("\nStep 2: Checking for mapped genes in annotations...")
    
    # Load existing mappings
    mapping_file = f"{mapped_dir}/simple_gene_mapping.tsv"
    if os.path.exists(mapping_file):
        try:
            mapping_df = pd.read_csv(mapping_file, sep='\t')
            print(f"Loaded mapping file with {len(mapping_df)} entries")
            
            # Extract valid mappings
            valid_mappings = mapping_df[mapping_df['W303_ID'] != 'NOT_FOUND']
            print(f"Found {len(valid_mappings)} valid mappings:")
            print(valid_mappings)
            
            # Check for these genes in VCF annotations
            if len(valid_mappings) > 0 and len(vcf_files) > 0:
                print("\nSearching for mapped genes in annotations...")
                sample_file = os.path.join(results_dir, vcf_files[0])
                
                for _, row in valid_mappings.iterrows():
                    sgd_gene = row['SGD_Gene']
                    w303_id = row['W303_ID']
                    
                    print(f"Looking for {sgd_gene} ({w303_id}) in annotations...")
                    
                    # Try using bcftools to extract annotations
                    try:
                        # Use grep to find lines with the gene
                        grep_cmd = f"bcftools view {sample_file} | grep -w '{w303_id}' | head -1"
                        grep_output = subprocess.check_output(grep_cmd, shell=True, text=True, stderr=subprocess.PIPE)
                        
                        if grep_output.strip():
                            print(f"  ✓ Found {w303_id} in annotations")
                            print(f"    Sample line: {grep_output.strip()}")
                        else:
                            print(f"  ✗ Could not find {w303_id} in annotations")
                            
                            # Try looser matching
                            print("    Trying broader pattern matching...")
                            # For w303_scaffold_15, try both w303_scaffold_15 and scaffold_15
                            if "scaffold" in w303_id:
                                # Extract just the number
                                scaffold_num = w303_id.split("_")[-1]
                                loose_pattern = f"scaffold.{scaffold_num}"
                                grep_cmd = f"bcftools view {sample_file} | grep '{loose_pattern}' | head -1"
                                try:
                                    loose_output = subprocess.check_output(grep_cmd, shell=True, text=True, stderr=subprocess.PIPE)
                                    if loose_output.strip():
                                        print(f"    ✓ Found with looser pattern '{loose_pattern}'")
                                        print(f"      Sample line: {loose_output.strip()}")
                                except subprocess.CalledProcessError:
                                    print(f"    ✗ Not found with looser pattern either")
                    except subprocess.CalledProcessError:
                        print(f"  ✗ Error searching for {w303_id}")
        except Exception as e:
            print(f"Error loading mapping file: {e}")
    else:
        print(f"Mapping file not found: {mapping_file}")
    
    # Step 3: Look more systematically at annotation formats
    print("\nStep 3: Examining annotation formats...")
    
    if len(vcf_files) > 0:
        sample_file = os.path.join(results_dir, vcf_files[0])
        
        # Extract some annotations to look for patterns
        try:
            # Use bcftools to extract all ANN fields
            ann_cmd = f"bcftools view {sample_file} | grep 'ANN=' | head -10 > {debug_dir}/ann_examples.txt"
            subprocess.run(ann_cmd, shell=True)
            print(f"Saved annotation examples to {debug_dir}/ann_examples.txt")
            
            # Count total variants and annotations
            try:
                total_variants = int(subprocess.check_output(
                    f"bcftools view {sample_file} | grep -v '^#' | wc -l",
                    shell=True, text=True
                ).strip())
                
                total_anns = int(subprocess.check_output(
                    f"bcftools view {sample_file} | grep 'ANN=' | wc -l",
                    shell=True, text=True
                ).strip())
                
                print(f"File contains {total_variants} variants, {total_anns} with annotations")
                
                # Extract representative gene IDs
                genes_cmd = f"bcftools view {sample_file} | grep 'ANN=' | grep -o 'ANN=[^;]*' | cut -d'|' -f4 | sort | uniq -c | sort -nr | head -20 > {debug_dir}/common_genes.txt"
                subprocess.run(genes_cmd, shell=True)
                
                print(f"Saved most common gene IDs to {debug_dir}/common_genes.txt")
                
                # Look at the common genes
                common_genes = subprocess.check_output(
                    f"cat {debug_dir}/common_genes.txt",
                    shell=True, text=True
                ).strip()
                
                print("\nMost common gene IDs in annotations:")
                print(common_genes)
                
                # Extract naming patterns from these genes
                gene_patterns = set()
                for line in common_genes.splitlines():
                    parts = line.strip().split()
                    if len(parts) > 1:
                        gene_id = parts[-1]
                        if gene_id:
                            # Extract pattern (prefix)
                            prefix = ''.join([c for c in gene_id if not c.isdigit()])
                            gene_patterns.add(prefix)
                
                print("\nGene ID patterns found:")
                for pattern in gene_patterns:
                    print(f"- {pattern}")
            except subprocess.CalledProcessError as e:
                print(f"Error counting variants: {e}")
        except subprocess.CalledProcessError as e:
            print(f"Error extracting annotations: {e}")
    
    # Step 4: Generate mapping recommendations
    print("\nStep 4: Generating mapping recommendations...")
    
    # Based on observed patterns, suggest mappings for remaining genes
    with open(f"{debug_dir}/mapping_recommendations.txt", 'w') as f:
        f.write("Gene Mapping Recommendations\n")
        f.write("==============================================\n\n")
        
        f.write("Based on the debug analysis, here are the key findings:\n\n")
        
        f.write("1. VCF File Access:\n")
        if 'zcat_output' in locals():
            f.write("   - VCF files are accessible with zcat\n")
        elif 'gunzip_output' in locals():
            f.write("   - VCF files are accessible with gunzip -c\n")
        elif 'bcftools_output' in locals():
            f.write("   - VCF files are accessible with bcftools\n")
        else:
            f.write("   - VCF files may have access issues\n")
        
        f.write("\n2. Current Gene Mappings:\n")
        if 'valid_mappings' in locals() and len(valid_mappings) > 0:
            for _, row in valid_mappings.iterrows():
                f.write(f"   - {row['SGD_Gene']} → {row['W303_ID']} (Confidence: {row['Confidence']})\n")
        else:
            f.write("   - No valid mappings found\n")
        
        f.write("\n3. Observed Gene ID Patterns:\n")
        if 'gene_patterns' in locals() and len(gene_patterns) > 0:
            for pattern in gene_patterns:
                f.write(f"   - {pattern}\n")
        else:
            f.write("   - No clear gene ID patterns observed\n")
        
        f.write("\n4. Recommendations for next steps:\n")
        f.write("   a. Manually check the GenBank files with grep for each target gene\n")
        f.write("      Example: grep -A 10 -B 10 'YHR190W' /Users/zakiralibhai/snpEff/data/w303/genes.gbk\n")
        f.write("   b. Try extracting the LOCUS line before each gene reference\n")
        f.write("   c. Look for gene name patterns in the SnpEff database dump\n")
        f.write("   d. Consider a more aggressive grep for partial matches\n")
        f.write("   e. For unmapped genes, try using other genes on the same chromosome as guides\n")
    
    print(f"Saved mapping recommendations to {debug_dir}/mapping_recommendations.txt")
    
    print("\n=== Debugging Complete ===")

if __name__ == "__main__":
    main()
