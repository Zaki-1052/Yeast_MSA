#!/usr/bin/env python3

# File: scripts/annotation/15_find_gene_mappings.py
# Purpose: Find mappings between SGD gene names and W303 identifiers

import os
import re
import subprocess
import pandas as pd
from datetime import datetime
import gzip

def main():
    print("=== Finding Gene Name Mappings ===")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("")
    
    # Define directories
    snpeff_dir = "/Users/zakiralibhai/snpEff"
    gene_dir = "annotation/genes_of_interest"
    mapped_dir = "annotation/mapped_genes"
    anno_source = "annotation/results_renamed"
    stats_dir = "annotation/stats_renamed"
    diag_dir = "annotation/diagnosis"
    
    # Create output directories
    os.makedirs(mapped_dir, exist_ok=True)
    os.makedirs(diag_dir, exist_ok=True)
    
    # Read target genes
    target_genes_file = f"{gene_dir}/target_genes.txt"
    with open(target_genes_file, 'r') as f:
        target_genes = [line.strip() for line in f if line.strip()]
    
    print(f"Looking for mappings for {len(target_genes)} target genes")
    print("")
    
    # Step 1: Examine GenBank file for gene mappings
    print("Step 1: Examining genes.gbk for SGD gene mappings...")
    genes_gbk = f"{snpeff_dir}/data/w303/genes.gbk"
    
    mappings = []
    sgd_patterns = []
    
    # Create search patterns for each target gene
    for gene in target_genes:
        # Regular patterns
        sgd_patterns.append(gene)
        # Without Y prefix
        if gene.startswith('Y'):
            sgd_patterns.append(gene[1:])
    
    # Compile unique patterns
    unique_patterns = list(set(sgd_patterns))
    pattern_str = '|'.join(unique_patterns)
    print(f"Searching for patterns: {pattern_str}")
    
    # Use grep to quickly find matches in the GenBank file
    try:
        grep_cmd = f"grep -A 15 -E '({pattern_str})' {genes_gbk}"
        output = subprocess.check_output(grep_cmd, shell=True, text=True)
        
        # Save output for analysis
        with open(f"{diag_dir}/gbk_matches.txt", 'w') as f:
            f.write(output)
        
        print(f"Found {len(output.splitlines())} lines with potential matches")
        print("Sample of matches:")
        print('\n'.join(output.splitlines()[:20]))
    except subprocess.CalledProcessError:
        print("No exact matches found with grep")
    
    # Step 2: Search for aliases or alternative formats
    print("\nStep 2: Searching for gene aliases or alternative naming formats...")
    
    # Create a more flexible search for SGD gene names in product or note fields
    sgd_mappings = []
    
    try:
        # Extract product and note fields containing target genes
        for gene in target_genes:
            base_gene = gene[1:] if gene.startswith('Y') else gene
            
            # Try multiple patterns
            patterns = [
                f"/product=\".*{gene}.*\"",
                f"/product=\".*{base_gene}.*\"",
                f"/note=\".*{gene}.*\"",
                f"/note=\".*{base_gene}.*\"",
                f"/gene=\"{gene}\"",
                f"/gene=\"{base_gene}\"",
                f"\\b{gene}\\b",
                f"\\b{base_gene}\\b"
            ]
            
            for pattern in patterns:
                try:
                    grep_cmd = f"grep -A 15 -B 5 '{pattern}' {genes_gbk}"
                    output = subprocess.check_output(grep_cmd, shell=True, text=True)
                    
                    # Found match for this gene
                    print(f"Found potential matches for {gene} using pattern '{pattern}'")
                    
                    # Parse out LOCUS and gene information
                    matches = output.split("--")
                    for match in matches:
                        w303_id = None
                        locus_tag = None
                        product = None
                        sgd_id = None
                        
                        # Extract LOCUS name
                        locus_match = re.search(r'LOCUS\s+(\S+)', match)
                        if locus_match:
                            w303_id = locus_match.group(1)
                        
                        # Extract locus_tag
                        tag_match = re.search(r'/locus_tag="([^"]+)"', match)
                        if tag_match:
                            locus_tag = tag_match.group(1)
                        
                        # Extract product
                        product_match = re.search(r'/product="([^"]+)"', match)
                        if product_match:
                            product = product_match.group(1)
                        
                        # Extract SGD ID if present
                        sgd_match = re.search(r'[YA-Z]{2,3}[0-9]{3}[WC](?:-[A-Z])?', match)
                        if sgd_match:
                            sgd_id = sgd_match.group(0)
                        
                        if w303_id and (sgd_id or (product and gene in product)):
                            sgd_mappings.append({
                                'SGD_Gene': gene,
                                'W303_ID': w303_id,
                                'Locus_Tag': locus_tag,
                                'Product': product,
                                'SGD_ID_Found': sgd_id,
                                'Pattern': pattern
                            })
                except subprocess.CalledProcessError:
                    continue
    
    # Remove duplicates based on W303_ID
    except Exception as e:
        print(f"Error processing SGD mappings: {e}")
    if sgd_mappings:
        unique_mappings = []
        seen_w303 = set()
        for mapping in sgd_mappings:
            if mapping['W303_ID'] not in seen_w303:
                unique_mappings.append(mapping)
                seen_w303.add(mapping['W303_ID'])
        
        sgd_mappings = unique_mappings
        
        # Save mappings to file
        mapping_df = pd.DataFrame(sgd_mappings)
        mapping_df.to_csv(f"{mapped_dir}/gene_mappings.tsv", sep='\t', index=False)
        
        print(f"\nFound {len(sgd_mappings)} potential gene mappings")
        print("Sample of mappings:")
        for i, mapping in enumerate(sgd_mappings[:min(5, len(sgd_mappings))]):
            print(f"{i+1}. {mapping['SGD_Gene']} → {mapping['W303_ID']} ({mapping['Product']})")
    else:
        print("No SGD mappings found")
    
    # Step 3: Check stats files for gene information
    print("\nStep 3: Checking stats files for gene names...")
    
    stats_files = [f for f in os.listdir(stats_dir) if f.endswith('.stats.genes.txt')]
    if stats_files:
        # Use first stats file
        stats_file = os.path.join(stats_dir, stats_files[0])
        
        # Extract top genes
        try:
            top_genes = subprocess.check_output(
                f"head -20 {stats_file}",
                shell=True, text=True
            ).strip()
            
            print("Top genes from stats file:")
            print(top_genes)
            
            # Save full gene stats
            with open(f"{diag_dir}/all_genes_stats.txt", 'w') as f:
                f.write(subprocess.check_output(
                    f"cat {stats_file}",
                    shell=True, text=True
                ))
            
            print(f"Full gene stats saved to {diag_dir}/all_genes_stats.txt")
        except subprocess.CalledProcessError:
            print("Error reading stats file")
    else:
        print("No stats files found")
    
    # Step 4: Extract gene IDs from annotated VCF files
    print("\nStep 4: Extracting gene IDs from annotated VCF files...")
    
    vcf_files = [f for f in os.listdir(anno_source) if f.endswith('.snpeff.vcf.gz')]
    if vcf_files:
        # Use first VCF file
        vcf_file = os.path.join(anno_source, vcf_files[0])
        
        # Extract genes from annotations
        genes_in_vcf = set()
        try:
            with gzip.open(vcf_file, 'rt') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    # Look for ANN field
                    if 'ANN=' in line:
                        ann_match = re.search(r'ANN=([^;]+)', line)
                        if ann_match:
                            ann = ann_match.group(1)
                            
                            # Extract gene IDs (4th field in pipe-delimited format)
                            parts = ann.split(',')
                            for part in parts:
                                fields = part.split('|')
                                if len(fields) >= 4 and fields[3]:
                                    genes_in_vcf.add(fields[3])
            
            print(f"Found {len(genes_in_vcf)} unique gene IDs in VCF")
            print("Sample of gene IDs:")
            sample_genes = list(genes_in_vcf)[:min(10, len(genes_in_vcf))]
            for i, gene in enumerate(sample_genes):
                print(f"{i+1}. {gene}")
            
            # Save all gene IDs to file
            with open(f"{diag_dir}/vcf_gene_ids.txt", 'w') as f:
                for gene in sorted(genes_in_vcf):
                    f.write(f"{gene}\n")
            
            print(f"All gene IDs saved to {diag_dir}/vcf_gene_ids.txt")
        except Exception as e:
            print(f"Error extracting genes from VCF: {e}")
    else:
        print("No annotated VCF files found")
    
    # Step 5: Create a manual mapping file template
    print("\nStep 5: Creating manual mapping file template...")
    
    # Create a template with target genes
    manual_mapping = []
    for gene in target_genes:
        manual_mapping.append({
            'SGD_Gene': gene,
            'W303_ID': '',
            'Notes': 'Fill in corresponding W303 ID'
        })
    
    manual_df = pd.DataFrame(manual_mapping)
    manual_file = f"{mapped_dir}/manual_mapping_template.tsv"
    manual_df.to_csv(manual_file, sep='\t', index=False)
    
    print(f"Manual mapping template saved to {manual_file}")
    print("You'll need to fill in the W303_ID column with the appropriate identifiers")
    
    # Step 6: Generate recommendations
    print("\nStep 6: Generating recommendations...")
    
    with open(f"{mapped_dir}/mapping_recommendations.txt", 'w') as f:
        f.write("Gene Mapping Recommendations\n")
        f.write("==============================================\n\n")
        
        f.write("The SnpEff database uses different gene identifiers than the SGD-style YHR190W format.\n")
        f.write("To find variants in your genes of interest, you'll need to:\n\n")
        
        f.write("1. Identify the corresponding W303 identifiers for your target genes\n")
        f.write("   - Look in the genes.gbk file for product or note fields containing your gene names\n")
        f.write("   - Search for gene descriptions or aliases that match your genes\n")
        f.write("   - Fill in the manual_mapping_template.tsv file with the correct mappings\n\n")
        
        f.write("2. Update your gene extraction script to:\n")
        f.write("   - Use the W303 gene identifiers instead of SGD names\n")
        f.write("   - Extract the correct field from the ANN string (appears to be field 4)\n")
        f.write("   - Map back to the SGD names in your final output\n\n")
        
        f.write("3. Alternative approaches:\n")
        f.write("   - Run SnpEff on the VCF files with the -canon option to get canonical gene names\n")
        f.write("   - Use an external database to map between identifier systems\n")
        f.write("   - Create a custom annotation using BED files with your genes of interest\n")
    
    print(f"Recommendations saved to {mapped_dir}/mapping_recommendations.txt")
    print("")
    print("=== Gene Mapping Analysis Complete ===")

if __name__ == "__main__":
    main()
