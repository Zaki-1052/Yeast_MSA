#!/usr/bin/env python3
"""
rename_vcf_chromosomes.py

This script renames chromosomes in VCF files to match the scaffold naming convention 
used in the SnpEff database (w303_scaffold_X).

Usage:
    python rename_vcf_chromosomes.py --vcf_dir <vcf_directory> --mapping_file <mapping_file> --output_dir <output_directory>
"""

import os
import argparse
import pandas as pd
import subprocess
import gzip
import shutil
from concurrent.futures import ProcessPoolExecutor

def read_chromosome_mapping(mapping_file):
    """
    Read chromosome mapping file.
    
    Args:
        mapping_file (str): Path to chromosome mapping file
        
    Returns:
        dict: Dictionary mapping CM007XXX.1 to w303_scaffold_X
    """
    try:
        mapping_df = pd.read_csv(mapping_file, sep='\t')
        # Create mapping dictionary
        mapping_dict = dict(zip(mapping_df['chromosome_id'], mapping_df['w303_scaffold']))
        print(f"Read mapping for {len(mapping_dict)} chromosomes")
        return mapping_dict
    except Exception as e:
        print(f"ERROR reading mapping file: {str(e)}")
        return {}

def rename_chromosomes_in_vcf(vcf_file, output_file, mapping_dict):
    """
    Rename chromosomes in a VCF file.
    
    Args:
        vcf_file (str): Input VCF file
        output_file (str): Output VCF file
        mapping_dict (dict): Chromosome mapping dictionary
        
    Returns:
        bool: Success status
    """
    try:
        # Determine if input is gzipped
        is_gzipped = vcf_file.endswith('.gz')
        opener = gzip.open if is_gzipped else open
        mode = 'rt' if is_gzipped else 'r'
        
        # Determine if output should be gzipped
        out_is_gzipped = output_file.endswith('.gz')
        
        # Create temporary output file (uncompressed)
        temp_output = output_file[:-3] if out_is_gzipped else output_file
        
        header_lines = []
        with opener(vcf_file, mode) as vcf, open(temp_output, 'w') as out:
            # Process header lines
            for line in vcf:
                if line.startswith('#'):
                    # Update contig lines
                    if line.startswith('##contig='):
                        for old_id, new_id in mapping_dict.items():
                            if f'ID={old_id}' in line:
                                line = line.replace(f'ID={old_id}', f'ID={new_id}')
                                break
                    header_lines.append(line)
                    out.write(line)
                else:
                    # First non-header line, start replacing chromosomes
                    break
            
            # Process variant lines
            line_count = 0
            renamed_count = 0
            
            # Process the first non-header line
            if line:
                fields = line.strip().split('\t')
                old_chrom = fields[0]
                
                if old_chrom in mapping_dict:
                    fields[0] = mapping_dict[old_chrom]
                    renamed_count += 1
                
                out.write('\t'.join(fields) + '\n')
                line_count += 1
            
            # Process remaining lines
            for line in vcf:
                fields = line.strip().split('\t')
                old_chrom = fields[0]
                
                if old_chrom in mapping_dict:
                    fields[0] = mapping_dict[old_chrom]
                    renamed_count += 1
                
                out.write('\t'.join(fields) + '\n')
                line_count += 1
        
        # Compress if needed
        if out_is_gzipped:
            with open(temp_output, 'rb') as f_in, gzip.open(output_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(temp_output)
        
        # Create tabix index for bgzipped VCF
        if out_is_gzipped:
            subprocess.run(['tabix', '-p', 'vcf', output_file], check=True)
        
        print(f"{os.path.basename(vcf_file)}: Renamed {renamed_count}/{line_count} variants")
        return True
    
    except Exception as e:
        print(f"ERROR processing {vcf_file}: {str(e)}")
        return False

def process_vcf_files(vcf_dir, output_dir, mapping_dict, max_workers=4):
    """
    Process all VCF files in directory.
    
    Args:
        vcf_dir (str): Directory containing VCF files
        output_dir (str): Output directory
        mapping_dict (dict): Chromosome mapping dictionary
        max_workers (int): Maximum number of parallel processes
        
    Returns:
        int: Number of successfully processed files
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all VCF files
    vcf_files = []
    for file in os.listdir(vcf_dir):
        if file.endswith('.vcf.gz') or file.endswith('.vcf'):
            vcf_files.append(os.path.join(vcf_dir, file))
    
    if not vcf_files:
        print(f"No VCF files found in {vcf_dir}")
        return 0
    
    print(f"Found {len(vcf_files)} VCF files to process")
    
    # Process files in parallel
    success_count = 0
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        
        for vcf_file in vcf_files:
            basename = os.path.basename(vcf_file)
            output_file = os.path.join(output_dir, basename)
            
            # Submit task
            future = executor.submit(rename_chromosomes_in_vcf, vcf_file, output_file, mapping_dict)
            futures.append((future, basename))
        
        # Process results
        for future, basename in futures:
            try:
                success = future.result()
                if success:
                    success_count += 1
            except Exception as e:
                print(f"ERROR processing {basename}: {str(e)}")
    
    return success_count

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Rename chromosomes in VCF files")
    parser.add_argument("--vcf_dir", required=True, help="Directory containing VCF files")
    parser.add_argument("--mapping_file", required=True, help="Chromosome mapping file (tsv)")
    parser.add_argument("--output_dir", required=True, help="Output directory for renamed VCF files")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel processes")
    args = parser.parse_args()
    
    # Read chromosome mapping
    mapping_dict = read_chromosome_mapping(args.mapping_file)
    if not mapping_dict:
        print("ERROR: Failed to read chromosome mapping")
        return
    
    # Process VCF files
    success_count = process_vcf_files(args.vcf_dir, args.output_dir, mapping_dict, args.threads)
    
    print(f"\nProcessing complete: Successfully renamed chromosomes in {success_count} VCF files")
    print(f"Renamed VCF files are in: {args.output_dir}")

if __name__ == "__main__":
    main()