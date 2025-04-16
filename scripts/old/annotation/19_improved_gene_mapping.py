#!/usr/bin/env python3

# File: scripts/annotation/19_improved_gene_mapping.py
# Purpose: Improved approach to extract W303 identifiers for target genes

import os
import re
from datetime import datetime
import pandas as pd
import subprocess

def main():
    print("=== Improved Gene Mapping Extraction ===")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("")
    
    # Define directories
    snpeff_dir = "/Users/zakiralibhai/snpEff"
    gene_dir = "annotation/genes_of_interest"
    mapped_dir = "annotation/mapped_genes"
    debug_dir = "annotation/debug_mappings"
    
    # Create output directories
    os.makedirs(mapped_dir, exist_ok=True)
    os.makedirs(debug_dir, exist_ok=True)
    
    # Read target genes
    target_genes_file = f"{gene_dir}/target_genes.txt"
    with open(target_genes_file, 'r') as f:
        target_genes = [line.strip() for line in f if line.strip()]
    
    # Load existing mappings if available
    existing_mappings = {}
    mapping_file = f"{mapped_dir}/simple_gene_mapping.tsv"
    if os.path.exists(mapping_file):
        try:
            mapping_df = pd.read_csv(mapping_file, sep='\t')
            for _, row in mapping_df.iterrows():
                if row['W303_ID'] != 'NOT_FOUND':
                    existing_mappings[row['SGD_Gene']] = row['W303_ID']
            
            print(f"Loaded {len(existing_mappings)} existing mappings")
        except Exception as e:
            print(f"Error loading existing mappings: {e}")
    
    # Identify genes we still need to find
    remaining_genes = [g for g in target_genes if g not in existing_mappings]
    print(f"Looking for mappings for {len(remaining_genes)} remaining genes")
    print("")
    
    # Try approach 1: Extract from test VCF annotations
    print("Approach 1: Creating a test VCF with target genes as comments...")
    
    # Create a test VCF file that includes target genes in INFO field
    test_vcf = f"{debug_dir}/test_genes.vcf"
    with open(test_vcf, 'w') as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##reference=w303\n")
        
        # Add each target gene as a comment
        for gene in remaining_genes:
            f.write(f"##INFO=<ID=TEST_{gene},Number=0,Type=Flag,Description=\"Test flag for {gene}\">\n")
        
        # Add chromosome info for all w303 scaffolds
        for i in range(1, 416):  # Up to 415 scaffolds
            f.write(f"##contig=<ID=w303_scaffold_{i}>\n")
        
        # Header line
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        
        # Add one variant on each chromosome with target gene in info
        for i, gene in enumerate(remaining_genes, 1):
            # Use a different scaffold for each gene
            f.write(f"w303_scaffold_{i}\t100\t.\tA\tT\t100\tPASS\tTEST_{gene}\n")
    
    print(f"Created test VCF: {test_vcf}")
    
    # Run SnpEff on the test VCF
    print("Running SnpEff on test VCF...")
    test_out = f"{debug_dir}/test_genes_annotated.vcf"
    
    cmd = f"java -jar {snpeff_dir}/snpEff.jar -v w303 {test_vcf} > {test_out}"
    try:
        subprocess.run(cmd, shell=True, check=True)
        print("SnpEff annotation completed")
        
        # Look for annotations that might reveal the gene mappings
        print("Checking annotations for gene information...")
        with open(test_out, 'r') as f:
            content = f.read()
            
            # Save for debugging
            with open(f"{debug_dir}/test_genes_annotations.txt", 'w') as debug_f:
                debug_f.write(content)
            
            # Look for gene annotations (ANN= fields)
            ann_matches = re.findall(r'ANN=([^;]+)', content)
            if ann_matches:
                print(f"Found {len(ann_matches)} annotation fields")
                for ann in ann_matches[:5]:  # Show first 5
                    print(f"  {ann[:100]}...")
                
                # Try to extract gene IDs
                for ann in ann_matches:
                    parts = ann.split('|')
                    if len(parts) >= 4 and parts[3]:
                        print(f"Found gene ID: {parts[3]}")
            else:
                print("No annotation fields found")
    except subprocess.CalledProcessError as e:
        print(f"Error running SnpEff: {e}")
    
    # Try approach 2: Extract from stats files
    print("\nApproach 2: Examining stats files for gene information...")
    
    stats_dir = "annotation/stats_renamed"
    if os.path.exists(stats_dir):
        stats_files = [f for f in os.listdir(stats_dir) if f.endswith('.stats.genes.txt')]
        if stats_files:
            # Use first stats file
            stats_file = os.path.join(stats_dir, stats_files[0])
            print(f"Using stats file: {stats_file}")
            
            # Look for gene information that might contain our target genes
            for gene in remaining_genes:
                print(f"Searching for {gene} in stats file...")
                
                # Try different patterns
                patterns = [
                    gene,                     # Exact match
                    gene.replace('Y', ''),    # Without Y prefix
                    gene.lower(),             # Lowercase
                ]
                
                found = False
                for pattern in patterns:
                    try:
                        cmd = f"grep '{pattern}' {stats_file}"
                        result = subprocess.check_output(cmd, shell=True, text=True)
                        if result:
                            print(f"  ✓ Found using pattern '{pattern}'")
                            print(f"    {result.strip()}")
                            found = True
                            break
                    except subprocess.CalledProcessError:
                        pass
                
                if not found:
                    print(f"  ✗ Not found in stats file")
        else:
            print("No stats files found")
    else:
        print(f"Stats directory not found: {stats_dir}")
    
    # Try approach 3: Look at sample VCF annotations
    print("\nApproach 3: Examining sample VCF annotations...")
    
    # Find an annotated VCF file
    results_dir = "annotation/results_renamed"
    vcf_files = [f for f in os.listdir(results_dir) if f.endswith('.snpeff.vcf.gz')]
    
    if vcf_files:
        sample_vcf = os.path.join(results_dir, vcf_files[0])
        print(f"Using sample VCF: {sample_vcf}")
        
        # Extract some annotations to see format
        try:
            cmd = f"zcat {sample_vcf} | grep 'ANN=' | head -5"
            annotations = subprocess.check_output(cmd, shell=True, text=True)
            
            print("Sample annotations:")
            for line in annotations.splitlines():
                if 'ANN=' in line:
                    ann_parts = line.split('ANN=')[1].split(';')[0]
                    for part in ann_parts.split(','):
                        fields = part.split('|')
                        if len(fields) >= 4:
                            print(f"  Gene ID: {fields[3]}")
                            break
                    break
        except subprocess.CalledProcessError as e:
            print(f"Error extracting annotations: {e}")
    else:
        print("No VCF files found")
    
    # Try approach 4: Use direct W303 to SGD conversion based on naming patterns
    print("\nApproach 4: Trying pattern-based conversion...")
    
    # Look at all the identifiers in the database
    w303_identifiers = []
    try:
        os.chdir(snpeff_dir)
        cmd = "java -jar snpEff.jar dump w303 | grep -A 50000 'gene:' | grep 'id:'"
        idents = subprocess.check_output(cmd, shell=True, text=True)
        os.chdir(os.path.dirname(os.path.abspath(__file__)) + "/../..")
        
        # Extract gene IDs
        for line in idents.splitlines():
            if "id:" in line:
                id_match = re.search(r'id: "([^"]+)"', line)
                if id_match:
                    w303_identifiers.append(id_match.group(1))
        
        print(f"Found {len(w303_identifiers)} gene identifiers in database")
        print("Sample gene IDs:")
        for id in w303_identifiers[:5]:
            print(f"  {id}")
    except subprocess.CalledProcessError as e:
        print(f"Error extracting gene identifiers: {e}")
    
    # Now try to find any that might match our target genes patterns
    pattern_mappings = {}
    
    # Try to find patterns like "W3030X00YYY" where X might be a chromosome letter and YYY might relate to position
    for gene in remaining_genes:
        # Extract chromosome and position from SGD name
        # Format is typically YXR012W where:
        # - X is chromosome (A-P)
        # - R is right/left arm (R/L)
        # - 012 is position
        # - W/C is Watson/Crick strand
        match = re.search(r'Y([A-P])([RL])(\d+)([WC])', gene)
        if match:
            chr_letter = match.group(1)
            arm = match.group(2)
            position = match.group(3)
            strand = match.group(4)
            
            print(f"{gene} -> Chr: {chr_letter}, Arm: {arm}, Pos: {position}, Strand: {strand}")
            
            # Look for W303 identifiers that might match this pattern
            possible_matches = []
            for w303_id in w303_identifiers:
                # Try some pattern matching heuristics
                if any(x in w303_id.upper() for x in [chr_letter, position]):
                    possible_matches.append(w303_id)
            
            if possible_matches:
                print(f"  Possible W303 matches for {gene}:")
                for i, match in enumerate(possible_matches[:5]):
                    print(f"    {i+1}. {match}")
                if len(possible_matches) > 5:
                    print(f"    ... and {len(possible_matches) - 5} more")
                    
                # Add the first match as a candidate
                if len(possible_matches) > 0:
                    pattern_mappings[gene] = possible_matches[0]
            else:
                print(f"  No matches found for {gene}")
    
    print(f"\nFound {len(pattern_mappings)} potential mappings using pattern analysis")
    
    # Combine results from all approaches
    new_mappings = {}
    
    # First, carry over existing mappings
    for gene, w303_id in existing_mappings.items():
        new_mappings[gene] = {'W303_ID': w303_id, 'Confidence': 'High'}
    
    # Add any new mappings from pattern analysis
    for gene, w303_id in pattern_mappings.items():
        if gene not in new_mappings:
            new_mappings[gene] = {'W303_ID': w303_id, 'Confidence': 'Low'}
    
    # For any remaining genes, mark as NOT_FOUND
    for gene in target_genes:
        if gene not in new_mappings:
            new_mappings[gene] = {'W303_ID': 'NOT_FOUND', 'Confidence': 'None'}
    
    # Save updated mappings
    updated_mappings = []
    for gene in target_genes:
        if gene in new_mappings:
            updated_mappings.append({
                'SGD_Gene': gene,
                'W303_ID': new_mappings[gene]['W303_ID'],
                'Confidence': new_mappings[gene]['Confidence']
            })
    
    if updated_mappings:
        updated_df = pd.DataFrame(updated_mappings)
        updated_file = f"{mapped_dir}/updated_gene_mapping.tsv"
        updated_df.to_csv(updated_file, sep='\t', index=False)
        
        print(f"\nUpdated mappings saved to: {updated_file}")
    
    print("\n=== Improved Gene Mapping Complete ===")

if __name__ == "__main__":
    main()
