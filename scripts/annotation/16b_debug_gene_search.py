#!/usr/bin/env python3

# File: scripts/annotation/16b_debug_gene_search.py
# Purpose: Simple grep-based search for target genes in GenBank file

import os
import subprocess
from datetime import datetime

def main():
    print("=== Basic Gene Search Debug ===")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("")
    
    # Define directories
    snpeff_dir = "/Users/zakiralibhai/snpEff"
    gene_dir = "annotation/genes_of_interest"
    debug_dir = "annotation/debug_mappings"
    
    # Create output directory
    os.makedirs(debug_dir, exist_ok=True)
    
    # Read target genes
    target_genes_file = f"{gene_dir}/target_genes.txt"
    with open(target_genes_file, 'r') as f:
        target_genes = [line.strip() for line in f if line.strip()]
    
    print(f"Looking for {len(target_genes)} target genes")
    print("")
    
    # Path to GenBank file
    genes_gbk = f"{snpeff_dir}/data/w303/genes.gbk"
    
    # Check file existence and size
    if os.path.exists(genes_gbk):
        file_size = os.path.getsize(genes_gbk) / (1024 * 1024)  # Size in MB
        print(f"GenBank file exists: {genes_gbk}")
        print(f"File size: {file_size:.2f} MB")
    else:
        print(f"ERROR: GenBank file not found: {genes_gbk}")
        return
    
    # Search for each target gene using basic grep
    all_results = {}
    for gene in target_genes:
        print(f"\nSearching for {gene}...")
        
        # Try different patterns
        patterns = [
            gene,                     # Exact match
            gene.replace('Y', ''),    # Without Y prefix
            gene.lower(),             # Lowercase
            f'"{gene}"',              # Quoted
            f'"{gene.replace("Y", "")}"'  # Quoted without Y
        ]
        
        found = False
        for pattern in patterns:
            cmd = f"grep -n '{pattern}' {genes_gbk}"
            try:
                result = subprocess.check_output(cmd, shell=True, text=True)
                if result:
                    print(f"  ✓ Found using pattern '{pattern}'")
                    first_match = result.split('\n')[0]
                    print(f"    First match (line number): {first_match[:100]}...")
                    
                    # Get more context for the first match
                    line_num = first_match.split(':')[0]
                    context_cmd = f"head -n {int(line_num) + 10} {genes_gbk} | tail -n 20"
                    context = subprocess.check_output(context_cmd, shell=True, text=True)
                    
                    all_results[gene] = {
                        'pattern': pattern,
                        'matches': result,
                        'context': context
                    }
                    
                    found = True
                    break
            except subprocess.CalledProcessError:
                pass
        
        if not found:
            print(f"  ✗ Not found with any pattern")
            all_results[gene] = {
                'pattern': None,
                'matches': None,
                'context': None
            }
    
    # Save detailed results to file
    with open(f"{debug_dir}/gene_search_results.txt", 'w') as f:
        f.write("Gene Search Results\n")
        f.write("==============================================\n\n")
        
        for gene, result in all_results.items():
            f.write(f"Gene: {gene}\n")
            if result['pattern']:
                f.write(f"  Found using pattern: '{result['pattern']}'\n")
                f.write("  Matches:\n")
                for line in result['matches'].split('\n')[:10]:  # First 10 matches
                    if line:
                        f.write(f"    {line}\n")
                if len(result['matches'].split('\n')) > 10:
                    f.write(f"    ... and {len(result['matches'].split('\n')) - 10} more\n")
                
                f.write("\n  Context of first match:\n")
                for line in result['context'].split('\n'):
                    f.write(f"    {line}\n")
            else:
                f.write("  Not found with any pattern\n")
            
            f.write("\n")
    
    print(f"\nDetailed results saved to: {debug_dir}/gene_search_results.txt")
    
    # Try alternative data source - check genes.txt file
    print("\nChecking alternative gene information sources...")
    
    genes_txt = f"{snpeff_dir}/data/w303/genes.txt"
    if os.path.exists(genes_txt):
        print(f"Found genes.txt file: {genes_txt}")
        
        # Look for our target genes
        for gene in target_genes:
            cmd = f"grep '{gene}' {genes_txt}"
            try:
                result = subprocess.check_output(cmd, shell=True, text=True)
                if result:
                    print(f"  ✓ Found {gene} in genes.txt")
                    print(f"    {result[:100]}...")
            except subprocess.CalledProcessError:
                print(f"  ✗ {gene} not found in genes.txt")
    else:
        print(f"genes.txt not found: {genes_txt}")
    
    # Look for gene information in snpEffectPredictor.bin by dumping info
    print("\nTrying to extract gene information from SnpEff database...")
    
    try:
        os.chdir(snpeff_dir)
        dump_cmd = "java -jar snpEff.jar dump w303 | grep -A 5 'Genes '"
        dump_result = subprocess.check_output(dump_cmd, shell=True, text=True)
        print("Database dump result:")
        print(dump_result)
        
        # Return to original directory
        os.chdir(os.path.dirname(os.path.abspath(__file__)) + "/../..")
    except subprocess.CalledProcessError as e:
        print(f"Error dumping gene information: {e}")
    
    # Check if we can use a different approach - search VCF annotations directly
    print("\nExamining VCF annotations directly to find gene identifiers...")
    
    # Find genes present in results
    try:
        cmd = "zcat annotation/results_renamed/*.snpeff.vcf.gz | grep 'ANN=' | head -5"
        ann_sample = subprocess.check_output(cmd, shell=True, text=True)
        print("Sample ANN fields:")
        print(ann_sample)
        
        # Extract gene IDs from ANN field
        print("\nExtracted Gene IDs from sample ANN fields:")
        for line in ann_sample.splitlines():
            if 'ANN=' in line:
                ann_parts = line.split('ANN=')[1].split(';')[0]
                for part in ann_parts.split(','):
                    fields = part.split('|')
                    if len(fields) >= 4:
                        print(f"  Gene ID: {fields[3]}")
        
        # Find common genes in annotations
        cmd = "zcat annotation/results_renamed/*.snpeff.vcf.gz | grep 'ANN=' | cut -d'|' -f4 | sort | uniq -c | sort -nr | head -20"
        common_genes = subprocess.check_output(cmd, shell=True, text=True)
        print("\nMost common genes in VCF annotations:")
        print(common_genes)
        
        # Save to file
        with open(f"{debug_dir}/common_genes_in_vcf.txt", 'w') as f:
            f.write("Most Common Genes in VCF Annotations\n")
            f.write("==============================================\n\n")
            f.write(common_genes)
        
        print(f"Common genes saved to: {debug_dir}/common_genes_in_vcf.txt")
    except subprocess.CalledProcessError as e:
        print(f"Error examining VCF annotations: {e}")
    
    print("\n=== Gene Search Debug Complete ===")

if __name__ == "__main__":
    main()
