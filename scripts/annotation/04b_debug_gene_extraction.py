#!/usr/bin/env python3

# File: scripts/annotation/04b_debug_gene_extraction.py
# Purpose: Debug why we're not finding variants in target genes

import os
import gzip
import re
from datetime import datetime

def main():
    print("=== Debugging Gene Extraction ===")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("")
    
    # Define directories
    anno_source = "annotation/results"
    gene_dest = "annotation/genes_of_interest"
    target_genes_file = f"{gene_dest}/target_genes.txt"
    debug_dir = "annotation/debug"
    
    # Create debug directory
    os.makedirs(debug_dir, exist_ok=True)
    
    # Read target genes
    print("Loading target genes...")
    with open(target_genes_file, 'r') as f:
        target_genes = [line.strip() for line in f if line.strip()]
    
    print(f"Found {len(target_genes)} target genes")
    print("")
    
    # Find all annotated VCF files
    print("Scanning for annotated VCF files...")
    vcf_files = []
    for root, dirs, files in os.walk(anno_source):
        for file in files:
            if file.endswith(".snpeff.vcf.gz"):
                vcf_files.append(os.path.join(root, file))
    
    if not vcf_files:
        print(f"ERROR: No annotated VCF files found in {anno_source}")
        return
    
    print(f"Found {len(vcf_files)} annotated VCF files")
    print("")
    
    # Select first file for detailed examination
    debug_file = vcf_files[0]
    debug_sample = os.path.basename(debug_file).replace(".snpeff.vcf.gz", "")
    print(f"Using {debug_sample} for detailed examination")
    
    # Examine annotation format
    print("Examining annotation format...")
    annotation_examples = []
    line_count = 0
    
    with gzip.open(debug_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            line_count += 1
            if line_count > 1000:  # Limit to first 1000 non-header lines
                break
                
            if "ANN=" in line:
                fields = line.strip().split('\t')
                info = fields[7]
                ann_match = re.search(r'ANN=([^;]+)', info)
                if ann_match and len(annotation_examples) < 5:
                    annotation_examples.append(ann_match.group(1))
    
    if annotation_examples:
        print("Found annotation examples:")
        for i, example in enumerate(annotation_examples):
            print(f"Example {i+1}:")
            print(example[:200] + "..." if len(example) > 200 else example)
            print("")
    else:
        print("No annotations found in sample file!")
    
    # Save complete examples to debug file
    with open(f"{debug_dir}/annotation_examples.txt", 'w') as f:
        f.write(f"Annotation Examples from {debug_sample}\n")
        f.write("==============================================\n\n")
        for i, example in enumerate(annotation_examples):
            f.write(f"Example {i+1}:\n")
            f.write(example + "\n\n")
    
    # Check for any gene mentions
    print("Checking for any gene-like mentions...")
    gene_pattern = r'\|[A-Z0-9]+\|'
    gene_examples = set()
    
    with gzip.open(debug_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            if "ANN=" in line:
                matches = re.findall(gene_pattern, line)
                for match in matches:
                    gene = match.strip('|')
                    if gene and len(gene_examples) < 20:
                        gene_examples.add(gene)
    
    if gene_examples:
        print(f"Found {len(gene_examples)} potential gene identifiers:")
        for gene in sorted(list(gene_examples))[:10]:
            print(f"- {gene}")
        if len(gene_examples) > 10:
            print(f"... and {len(gene_examples) - 10} more")
    else:
        print("No gene-like identifiers found!")
    
    # Save gene examples to debug file
    with open(f"{debug_dir}/gene_examples.txt", 'w') as f:
        f.write(f"Gene Examples from {debug_sample}\n")
        f.write("==============================================\n\n")
        for gene in sorted(list(gene_examples)):
            f.write(f"{gene}\n")
    
    # Search for target genes with different patterns
    print("\nSearching for target genes with more flexible patterns...")
    
    for gene in target_genes:
        print(f"Searching for {gene}...")
        
        # Try different patterns
        patterns = [
            f"|{gene}|",        # Exact match with delimiters
            gene,               # Exact match anywhere
            gene.replace("Y", ""), # Without Y prefix
            gene[:-1],          # Without last character
            gene.lower(),       # Lowercase
            gene[1:]            # Without first character
        ]
        
        found = False
        for pattern in patterns:
            with gzip.open(debug_file, 'rt') as f:
                for line in f:
                    if pattern in line:
                        print(f"  ✓ Found with pattern '{pattern}'")
                        found = True
                        break
            if found:
                break
        
        if not found:
            print(f"  ✗ Not found with any pattern")
    
    # Examine SnpEff configuration and database
    print("\nExamining directory structure to check SnpEff database setup...")
    
    snpeff_dir = "/Users/zakiralibhai/snpEff"
    snpeff_config = os.path.join(snpeff_dir, "snpEff.config")
    snpeff_data = os.path.join(snpeff_dir, "data", "w303")
    
    if os.path.exists(snpeff_config):
        print(f"SnpEff config exists: {snpeff_config}")
    else:
        print(f"SnpEff config not found: {snpeff_config}")
    
    if os.path.exists(snpeff_data):
        print(f"w303 database exists: {snpeff_data}")
        
        # Check for genes.gbk file
        genes_gbk = os.path.join(snpeff_data, "genes.gbk")
        if os.path.exists(genes_gbk):
            print(f"  genes.gbk exists: {genes_gbk}")
            
            # Check file size
            size = os.path.getsize(genes_gbk) / (1024 * 1024)  # Size in MB
            print(f"  genes.gbk size: {size:.2f} MB")
            
            # Check for target genes in genes.gbk
            print("  Checking for target genes in genes.gbk...")
            found_in_gbk = 0
            
            with open(genes_gbk, 'rb') as f:
                content = f.read()
                for gene in target_genes:
                    if gene.encode() in content:
                        print(f"    ✓ Found {gene} in genes.gbk")
                        found_in_gbk += 1
            
            if found_in_gbk == 0:
                print("    ✗ No target genes found in genes.gbk")
        else:
            print(f"  genes.gbk not found: {genes_gbk}")
    else:
        print(f"w303 database not found: {snpeff_data}")
    
    print("\nDebugging information saved to: {debug_dir}/")
    print("\nRecommended next steps:")
    print("1. Examine the annotation examples to understand the format")
    print("2. Check if genes.gbk contains the target genes in the expected format")
    print("3. Modify the extraction script to match the gene naming convention used by SnpEff")
    print("")
    
    print("=== Debugging complete ===")

if __name__ == "__main__":
    main()
