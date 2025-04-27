#!/usr/bin/env python3
"""
rename_all_vcfs.py

This script renames chromosomes in all VCF files using bcftools.
"""

import os
import argparse
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

def rename_vcf(input_vcf, mapping_file, output_dir):
    """Rename chromosomes in a VCF file using bcftools."""
    # Get base filename
    basename = os.path.basename(input_vcf)
    output_vcf = os.path.join(output_dir, basename.replace('.vcf.gz', '.renamed.vcf.gz'))
    
    # Run bcftools rename
    cmd = f"bcftools annotate --rename-chrs {mapping_file} {input_vcf} -Oz -o {output_vcf}"
    
    try:
        subprocess.run(cmd, shell=True, check=True)
        
        # Index the output VCF
        subprocess.run(f"bcftools index {output_vcf}", shell=True, check=True)
        
        return (basename, True)
    
    except subprocess.CalledProcessError as e:
        return (basename, False)

def process_vcf_dir(vcf_dir, mapping_file, output_dir, max_workers=4):
    """Process all VCF files in a directory."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all VCF files
    vcf_files = []
    for file in os.listdir(vcf_dir):
        if file.endswith('.vcf.gz'):
            vcf_files.append(os.path.join(vcf_dir, file))
    
    if not vcf_files:
        print(f"No VCF files found in {vcf_dir}")
        return []
    
    print(f"Found {len(vcf_files)} VCF files to process")
    
    # Process files in parallel
    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(rename_vcf, vcf_file, mapping_file, output_dir): vcf_file
            for vcf_file in vcf_files
        }
        
        for future in as_completed(futures):
            basename, success = future.result()
            results.append((basename, success))
            status = "SUCCESS" if success else "FAILED"
            print(f"{basename}: {status}")
    
    return results

def main():
    parser = argparse.ArgumentParser(description="Rename chromosomes in all VCF files")
    parser.add_argument("--vcf_dir", required=True, help="Directory containing VCF files")
    parser.add_argument("--mapping", required=True, help="bcftools mapping file")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel processes")
    args = parser.parse_args()
    
    # Process all VCF files
    results = process_vcf_dir(args.vcf_dir, args.mapping, args.output_dir, args.threads)
    
    # Summarize results
    successful = sum(1 for _, success in results if success)
    total = len(results)
    
    print(f"\nProcessing complete: {successful}/{total} VCF files successfully renamed")
    print(f"Renamed VCF files are in: {args.output_dir}")

if __name__ == "__main__":
    main()