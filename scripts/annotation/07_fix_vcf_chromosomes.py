#!/usr/bin/env python3

# File: scripts/annotation/07_fix_vcf_chromosomes.py
# Purpose: Fix chromosome names in VCF files and re-run SnpEff annotation

import os
import gzip
import re
import subprocess
import glob
from datetime import datetime
import pandas as pd

def extract_scaffold_number(filename):
    """Extract scaffold number from filename like w303_13Apr2025_scaffold_123.genbank"""
    match = re.search(r'scaffold_(\d+)\.genbank$', filename)
    if match:
        return int(match.group(1))
    return None

def extract_jriu_number(chrom):
    """Extract number from JRIU identifier like JRIU01000123.1"""
    match = re.search(r'JRIU\d+0*(\d+)\.1$', chrom)
    if match:
        return int(match.group(1))
    return None

def main():
    print("=== Fixing VCF Chromosome Names and Re-running Annotation ===")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("")
    
    # Define directories
    vcf_source = "annotation/vcf_ready"
    gbk_dir = "annotation/reference/w303_scaffolds"
    fixed_vcf_dir = "annotation/vcf_fixed_chrom"
    annotation_dir = "annotation/results_fixed"
    stats_dir = "annotation/stats_fixed"
    mapping_dir = "annotation/chromosome_mapping"
    
    # Create output directories
    os.makedirs(fixed_vcf_dir, exist_ok=True)
    os.makedirs(annotation_dir, exist_ok=True)
    os.makedirs(stats_dir, exist_ok=True)
    os.makedirs(mapping_dir, exist_ok=True)
    
    # Step 1: Get a list of all GenBank scaffolds
    print("Getting list of all GenBank scaffolds...")
    gbk_files = glob.glob(f"{gbk_dir}/*.genbank")
    
    if not gbk_files:
        print(f"ERROR: No GenBank files found in {gbk_dir}")
        return
    
    # Create a mapping from scaffold number to full scaffold name
    scaffold_dict = {}
    for gbk_file in gbk_files:
        scaffold_name = os.path.basename(gbk_file).replace('.genbank', '')
        scaffold_num = extract_scaffold_number(gbk_file)
        if scaffold_num is not None:
            scaffold_dict[scaffold_num] = scaffold_name
    
    print(f"Found {len(scaffold_dict)} scaffold mappings")
    
    # Step 2: Get a list of all JRIU chromosomes from VCF files
    print("Getting list of chromosomes from VCF files...")
    vcf_files = glob.glob(f"{vcf_source}/*.sorted.vcf.gz")
    
    if not vcf_files:
        print(f"ERROR: No VCF files found in {vcf_source}")
        return
    
    # Extract chromosomes from first VCF file
    vcf_chroms = []
    try:
        output = subprocess.check_output(
            f"bcftools view -h {vcf_files[0]} | grep '##contig=' | head -10",
            shell=True, text=True
        )
        print("Sample of contig lines from VCF:")
        print(output)
        
        full_output = subprocess.check_output(
            f"bcftools view -h {vcf_files[0]} | grep '##contig='",
            shell=True, text=True
        )
        
        # Extract all contig IDs
        for line in full_output.splitlines():
            match = re.search(r'ID=([^,]+)', line)
            if match:
                chrom = match.group(1)
                vcf_chroms.append(chrom)
    except subprocess.CalledProcessError:
        print(f"Error extracting contigs from VCF")
        return
    
    print(f"Found {len(vcf_chroms)} chromosomes in VCF files")
    
    # Step 3: Create a mapping between JRIU and scaffold numbers
    print("Creating chromosome mapping...")
    mapping = []
    
    for chrom in vcf_chroms:
        jriu_num = extract_jriu_number(chrom)
        if jriu_num is not None and jriu_num in scaffold_dict:
            mapping.append({
                'VCF_Chromosome': chrom,
                'SnpEff_Chromosome': scaffold_dict[jriu_num],
                'JRIU_Number': jriu_num,
                'Notes': 'Pattern-based mapping'
            })
    
    # Add any missing mappings in sequential order
    jriu_nums = sorted([extract_jriu_number(c) for c in vcf_chroms if extract_jriu_number(c) is not None])
    scaffold_nums = sorted(scaffold_dict.keys())
    
    # Ensure we have enough scaffold numbers for all JRIU numbers
    if len(jriu_nums) > len(scaffold_nums):
        print(f"WARNING: More JRIU numbers ({len(jriu_nums)}) than scaffold numbers ({len(scaffold_nums)})")
        missing_count = len(jriu_nums) - len(scaffold_nums)
        print(f"Some mappings may be missing ({missing_count})")
    
    # Create mapping for any JRIU numbers that didn't get mapped yet
    mapped_jriu_nums = [m['JRIU_Number'] for m in mapping]
    for jriu_num in jriu_nums:
        if jriu_num not in mapped_jriu_nums and jriu_num <= len(scaffold_nums):
            # For unmapped JRIU, use the scaffold with the same position in the sorted list
            scaffold_name = scaffold_dict[scaffold_nums[jriu_num - 1]]
            mapping.append({
                'VCF_Chromosome': f"JRIU01000{jriu_num:03d}.1",
                'SnpEff_Chromosome': scaffold_name,
                'JRIU_Number': jriu_num,
                'Notes': 'Position-based mapping'
            })
    
    # Create mapping DataFrame
    mapping_df = pd.DataFrame(mapping)
    mapping_file = f"{mapping_dir}/chromosome_mapping_complete.tsv"
    mapping_df.to_csv(mapping_file, sep='\t', index=False)
    
    print(f"Saved {len(mapping)} mappings to: {mapping_file}")
    print("")
    
    # Create lookup dictionary for faster access
    chrom_lookup = {row['VCF_Chromosome']: row['SnpEff_Chromosome'] for _, row in mapping_df.iterrows()}
    
    # Step 4: Modify the VCF files to use the SnpEff-compatible chromosome names
    print("Modifying VCF files...")
    
    for vcf_file in vcf_files:
        sample = os.path.basename(vcf_file).replace('.sorted.vcf.gz', '')
        output_file = f"{fixed_vcf_dir}/{sample}.fixed.vcf"
        
        print(f"- Processing {sample}...")
        
        # Create a new header with modified contig lines
        header_file = f"{fixed_vcf_dir}/{sample}.header.vcf"
        
        with open(header_file, 'w') as out:
            # Get header from original file
            header = subprocess.check_output(
                f"bcftools view -h {vcf_file}",
                shell=True, text=True
            )
            
            # Replace contig IDs in header
            header_lines = header.splitlines()
            modified_header_lines = []
            
            for line in header_lines:
                if line.startswith('##contig='):
                    match = re.search(r'ID=([^,]+)', line)
                    if match:
                        old_id = match.group(1)
                        if old_id in chrom_lookup:
                            new_id = chrom_lookup[old_id]
                            modified_line = line.replace(f'ID={old_id}', f'ID={new_id}')
                            modified_header_lines.append(modified_line)
                else:
                    modified_header_lines.append(line)
            
            # Write modified header
            out.write('\n'.join(modified_header_lines) + '\n')
        
        # Create a modified VCF with the new header and chromosome names
        with gzip.open(vcf_file, 'rt') as vcf_in, open(output_file, 'w') as vcf_out:
            # Skip header in input file (we'll use our modified header)
            for line in vcf_in:
                if line.startswith('#'):
                    continue
                
                # Replace chromosome name in variant lines
                fields = line.strip().split('\t')
                if len(fields) >= 8:
                    old_chrom = fields[0]
                    if old_chrom in chrom_lookup:
                        fields[0] = chrom_lookup[old_chrom]
                        vcf_out.write('\t'.join(fields) + '\n')
        
        # Combine header and variant records
        os.system(f"cat {header_file} >> {output_file}")
        
        # Compress and index
        print(f"  Compressing and indexing {sample}...")
        os.system(f"bgzip -f {output_file}")
        os.system(f"tabix -p vcf {output_file}.gz")
        
        # Clean up temporary header file
        os.remove(header_file)
    
    print("VCF modification complete")
    print("")
    
    # Step 5: Run SnpEff annotation on the fixed VCF files
    print("Running SnpEff annotation on fixed VCF files...")
    
    # Path to SnpEff
    snpeff_path = "/Users/zakiralibhai/snpEff"
    snpeff_jar = f"{snpeff_path}/snpEff.jar"
    
    # Ensure SnpEff jar exists
    if not os.path.exists(snpeff_jar):
        print(f"ERROR: SnpEff jar not found at {snpeff_jar}")
        return
    
    fixed_vcf_files = glob.glob(f"{fixed_vcf_dir}/*.fixed.vcf.gz")
    
    for vcf_file in fixed_vcf_files:
        sample = os.path.basename(vcf_file).replace('.fixed.vcf.gz', '')
        output_file = f"{annotation_dir}/{sample}.snpeff.vcf"
        stats_file = f"{stats_dir}/{sample}.snpeff.stats.html"
        
        print(f"- Annotating {sample}...")
        
        # Run SnpEff
        cmd = f"java -Xmx4g -jar {snpeff_jar} -v -stats {stats_file} w303 {vcf_file} > {output_file}"
        os.system(cmd)
        
        # Compress and index
        os.system(f"bgzip -f {output_file}")
        os.system(f"tabix -p vcf {output_file}.gz")
    
    print("SnpEff annotation complete")
    print("")
    
    print("=== VCF Chromosome Fixing and Annotation complete ===")

if __name__ == "__main__":
    main()
