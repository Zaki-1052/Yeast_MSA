#!/usr/bin/env python3

# File: scripts/annotation/25_extract_target_gene_regions.py
# Purpose: Extract target gene locations and find variants directly from VCF files

import os
import subprocess
import re
import pandas as pd
from datetime import datetime

def main():
    print("=== Extracting Target Gene Regions ===")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("")
    
    # Define directories
    snpeff_dir = "/Users/zakiralibhai/snpEff"
    gene_dir = "annotation/genes_of_interest"
    vcf_source = "annotation/vcf_ready"
    output_dir = "annotation/gene_variants_direct"
    debug_dir = "annotation/debug_gene_locations"
    
    # Create output directories
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(f"{output_dir}/by_gene", exist_ok=True)
    os.makedirs(debug_dir, exist_ok=True)
    
    # Read target genes
    target_genes_file = f"{gene_dir}/target_genes.txt"
    with open(target_genes_file, 'r') as f:
        target_genes = [line.strip() for line in f if line.strip()]
    
    print(f"Looking for {len(target_genes)} target genes")
    print("")
    
    # Step 1: Find gene locations in the SnpEff database
    print("Step 1: Finding gene locations in the SnpEff database...")
    
    gene_locations = {}
    
    try:
        # Move to the SnpEff directory
        os.chdir(snpeff_dir)
        
        # Search for our target genes in the genes.gbk file
        gbk_file = "data/w303/genes.gbk"
        for gene in target_genes:
            print(f"Searching for {gene}...")
            
            # Use grep to find gene references
            cmd = f"grep -n -A 10 -B 10 '{gene}' {gbk_file} | head -50"
            try:
                output = subprocess.check_output(cmd, shell=True, text=True)
                if output:
                    print(f"  ✓ Found references to {gene}")
                    
                    # Save to file for inspection
                    with open(f"{debug_dir}/{gene}_references.txt", 'w') as f:
                        f.write(output)
                    
                    # Try to extract the LOCUS and gene boundaries
                    locus = None
                    start = None
                    end = None
                    strand = None
                    
                    # Look for LOCUS line before the gene reference
                    locus_match = re.search(r'LOCUS\s+(\S+)', output)
                    if locus_match:
                        locus = locus_match.group(1)
                    
                    # Look for gene coordinates
                    loc_match = re.search(r'(\d+)\.\.(\d+)', output)
                    if loc_match:
                        start = int(loc_match.group(1))
                        end = int(loc_match.group(2))
                    
                    # Determine strand (complement or not)
                    if "complement(" in output:
                        strand = "-"
                    else:
                        strand = "+"
                    
                    if locus and start and end:
                        gene_locations[gene] = {
                            'scaffold': locus,
                            'start': start,
                            'end': end,
                            'strand': strand
                        }
                        
                        print(f"    Location: {locus}:{start}-{end} ({strand})")
                    else:
                        print(f"    ✗ Could not extract location information")
                else:
                    print(f"  ✗ No references found for {gene}")
            except subprocess.CalledProcessError:
                print(f"  ✗ Error searching for {gene}")
        
        # Return to original directory
        os.chdir(os.path.dirname(os.path.abspath(__file__)) + "/../..")
    except Exception as e:
        print(f"Error finding gene locations: {e}")
        # In case of error, try to return to original directory
        try:
            os.chdir(os.path.dirname(os.path.abspath(__file__)) + "/../..")
        except:
            pass
    
    # Save gene locations to file
    if gene_locations:
        locations_df = pd.DataFrame([
            {
                'Gene': gene,
                'Scaffold': info['scaffold'],
                'Start': info['start'],
                'End': info['end'],
                'Strand': info['strand']
            }
            for gene, info in gene_locations.items()
        ])
        
        locations_df.to_csv(f"{output_dir}/gene_locations.tsv", sep='\t', index=False)
        print(f"\nSaved {len(gene_locations)} gene locations to {output_dir}/gene_locations.tsv")
    else:
        print("\nNo gene locations found")
    
    # Step 2: Extract all possible variants from target regions
    print("\nStep 2: Extracting variants from target gene regions...")
    
    # Find VCF files
    vcf_files = []
    for root, dirs, files in os.walk(vcf_source):
        for file in files:
            if file.endswith('.sorted.vcf.gz'):
                vcf_files.append(os.path.join(root, file))
    
    if not vcf_files:
        print(f"ERROR: No VCF files found in {vcf_source}")
        return
    
    print(f"Found {len(vcf_files)} VCF files")
    
    # Get the JRIU to W303 scaffold mapping
    mapping_file = "annotation/chromosome_mapping/jriu_to_w303_direct_mapping.tsv"
    jriu_to_w303 = {}
    w303_to_jriu = {}
    
    if os.path.exists(mapping_file):
        mapping_df = pd.read_csv(mapping_file, sep='\t')
        jriu_to_w303 = dict(zip(mapping_df['JRIU_Chromosome'], mapping_df['W303_Scaffold']))
        w303_to_jriu = dict(zip(mapping_df['W303_Scaffold'], mapping_df['JRIU_Chromosome']))
        print(f"Loaded {len(jriu_to_w303)} chromosome mappings")
    else:
        print(f"WARNING: Chromosome mapping file not found: {mapping_file}")
    
    # Extract variants for each gene
    all_variants = []
    gene_variants = {gene: [] for gene in target_genes}
    
    for gene, location in gene_locations.items():
        scaffold = location['scaffold']
        start = location['start']
        end = location['end']
        
        print(f"\nExtracting variants for {gene} ({scaffold}:{start}-{end})...")
        
        # Find corresponding JRIU chromosome(s)
        jriu_chroms = []
        if scaffold in w303_to_jriu:
            jriu_chroms.append(w303_to_jriu[scaffold])
        else:
            # Try to match by scaffold number
            match = re.search(r'scaffold_(\d+)', scaffold)
            if match:
                scaffold_num = int(match.group(1))
                
                # Look for JRIU chromosomes with similar numbers
                for jriu_chrom, w303_scaffold in jriu_to_w303.items():
                    if f"scaffold_{scaffold_num}" in w303_scaffold:
                        jriu_chroms.append(jriu_chrom)
        
        if not jriu_chroms:
            print(f"  ✗ Could not find corresponding JRIU chromosome for {scaffold}")
            continue
        
        print(f"  Corresponding JRIU chromosome(s): {', '.join(jriu_chroms)}")
        
        # Process each VCF file
        for vcf_file in vcf_files:
            sample = os.path.basename(vcf_file).replace('.sorted.vcf.gz', '')
            
            variant_count = 0
            for jriu_chrom in jriu_chroms:
                try:
                    # Use bcftools to extract variants in this region
                    cmd = f"bcftools view {vcf_file} {jriu_chrom} | grep -v '^#'"
                    output = subprocess.check_output(cmd, shell=True, text=True, stderr=subprocess.DEVNULL)
                    
                    for line in output.splitlines():
                        fields = line.strip().split('\t')
                        if len(fields) < 8:
                            continue
                        
                        chrom = fields[0]
                        pos = int(fields[1])
                        ref = fields[3]
                        alt = fields[4]
                        info = fields[7]
                        
                        # Create variant record
                        record = {
                            'Sample': sample,
                            'Gene': gene,
                            'Chromosome': chrom,
                            'Position': pos,
                            'W303_Scaffold': scaffold,
                            'W303_Start': start,
                            'W303_End': end,
                            'Ref': ref,
                            'Alt': alt,
                            'Info': info
                        }
                        
                        all_variants.append(record)
                        gene_variants[gene].append(record)
                        variant_count += 1
                except subprocess.CalledProcessError:
                    # No variants for this chromosome
                    pass
            
            if variant_count > 0:
                print(f"  Found {variant_count} variants in {sample}")
    
    # Save results
    all_df = pd.DataFrame(all_variants)
    if not all_df.empty:
        all_df.to_csv(f"{output_dir}/all_variants.tsv", sep='\t', index=False)
        print(f"\nSaved {len(all_df)} total variants to {output_dir}/all_variants.tsv")
        
        # Save variants by gene
        for gene, variants in gene_variants.items():
            if variants:
                gene_df = pd.DataFrame(variants)
                gene_df.to_csv(f"{output_dir}/by_gene/{gene}_variants.tsv", sep='\t', index=False)
                print(f"Saved {len(gene_df)} variants for {gene}")
    else:
        print("\nNo variants found in target gene regions")
    
    # Step 3: Analyze the GenBank file to find diagnostic information
    print("\nStep 3: Extracting diagnostic information from GenBank file...")
    
    # Try a more aggressive search for our target genes
    try:
        os.chdir(snpeff_dir)
        
        # First, try a word boundary search
        for gene in target_genes:
            if gene not in gene_locations:
                print(f"Trying alternative search methods for {gene}...")
                
                # Try word boundary search
                cmd = f"grep -n -B 3 -A 3 '\\b{gene}\\b' data/w303/genes.gbk | head -20"
                try:
                    output = subprocess.check_output(cmd, shell=True, text=True)
                    if output:
                        print(f"  ✓ Found with word boundary search: {gene}")
                        with open(f"{debug_dir}/{gene}_word_boundary.txt", 'w') as f:
                            f.write(output)
                except subprocess.CalledProcessError:
                    pass
                
                # Try case-insensitive search
                cmd = f"grep -n -i -B 3 -A 3 '{gene}' data/w303/genes.gbk | head -20"
                try:
                    output = subprocess.check_output(cmd, shell=True, text=True)
                    if output:
                        print(f"  ✓ Found with case-insensitive search: {gene}")
                        with open(f"{debug_dir}/{gene}_case_insensitive.txt", 'w') as f:
                            f.write(output)
                except subprocess.CalledProcessError:
                    pass
                
                # Try searching without the Y prefix
                if gene.startswith('Y'):
                    gene_no_y = gene[1:]
                    cmd = f"grep -n -B 3 -A 3 '\\b{gene_no_y}\\b' data/w303/genes.gbk | head -20"
                    try:
                        output = subprocess.check_output(cmd, shell=True, text=True)
                        if output:
                            print(f"  ✓ Found without Y prefix: {gene} as {gene_no_y}")
                            with open(f"{debug_dir}/{gene}_no_y_prefix.txt", 'w') as f:
                                f.write(output)
                    except subprocess.CalledProcessError:
                        pass
        
        # Try to extract one complete gene entry
        cmd = "grep -A 50 'gene ' data/w303/genes.gbk | head -50"
        try:
            output = subprocess.check_output(cmd, shell=True, text=True)
            with open(f"{debug_dir}/sample_gene_entry.txt", 'w') as f:
                f.write(output)
            print(f"Saved sample gene entry to {debug_dir}/sample_gene_entry.txt")
        except subprocess.CalledProcessError:
            pass
        
        # Return to original directory
        os.chdir(os.path.dirname(os.path.abspath(__file__)) + "/../..")
    except Exception as e:
        print(f"Error in diagnostic extraction: {e}")
        # Try to return to original directory
        try:
            os.chdir(os.path.dirname(os.path.abspath(__file__)) + "/../..")
        except:
            pass
    
    print("\n=== Target Gene Region Extraction Complete ===")

if __name__ == "__main__":
    main()
