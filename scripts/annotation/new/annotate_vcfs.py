#!/usr/bin/env python3
"""
annotate_vcfs.py

This script performs batch annotation of VCF files using SnpEff.

Usage:
    python annotate_vcfs.py --vcf_dir <vcf_directory> --output_dir <output_directory> \
                            --snpeff_dir <snpeff_directory> --genome <genome_name> \
                            --threads <number_of_threads>
"""

import os
import argparse
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import time

def annotate_vcf(vcf_file, output_dir, snpeff_dir, genome, stats_dir):
    """
    Annotate a VCF file using SnpEff.
    
    Args:
        vcf_file (str): Path to VCF file
        output_dir (str): Output directory for annotated VCF
        snpeff_dir (str): SnpEff installation directory
        genome (str): Genome name in SnpEff
        stats_dir (str): Directory for annotation statistics
        
    Returns:
        tuple: (basename, success, message)
    """
    basename = os.path.basename(vcf_file)
    sample_name = basename.split('.')[0]
    output_vcf = os.path.join(output_dir, f"{sample_name}.annotated.vcf")
    stats_file = os.path.join(stats_dir, f"{sample_name}.stats.html")
    genes_file = os.path.join(stats_dir, f"{sample_name}.genes.txt")
    
    start_time = time.time()
    
    cmd = [
        "java", "-Xmx4g", "-jar", os.path.join(snpeff_dir, "snpEff.jar"),
        "-v",
        "-stats", stats_file, 
        "-csvStats", os.path.join(stats_dir, f"{sample_name}.stats.csv"),
        genome,
        vcf_file
    ]
    
    cmd_str = " ".join(cmd) + f" > {output_vcf}"
    
    try:
        # Run SnpEff
        process = subprocess.run(cmd_str, shell=True, capture_output=True, text=True)
        
        # Check for errors
        if process.returncode != 0:
            return (basename, False, f"SnpEff returned code {process.returncode}: {process.stderr}")
        
        # Check if output file exists and has content
        if not os.path.exists(output_vcf):
            return (basename, False, "Output file not created")
        
        if os.path.getsize(output_vcf) == 0:
            return (basename, False, "Output file is empty")
        
        # Check for ERROR_CHROMOSOME_NOT_FOUND in output
        with open(output_vcf, 'r') as f:
            if "ERROR_CHROMOSOME_NOT_FOUND" in f.read():
                return (basename, False, "Chromosome not found error in output")
        
        # Success
        duration = time.time() - start_time
        return (basename, True, f"Completed in {duration:.1f} seconds")
    
    except Exception as e:
        return (basename, False, f"Exception: {str(e)}")

def batch_annotate(vcf_dir, output_dir, snpeff_dir, genome, threads=4):
    """
    Annotate multiple VCF files using SnpEff.
    
    Args:
        vcf_dir (str): Directory containing VCF files
        output_dir (str): Output directory for annotated VCFs
        snpeff_dir (str): SnpEff installation directory
        genome (str): Genome name in SnpEff
        threads (int): Number of parallel processes
        
    Returns:
        list: List of (basename, success, message) tuples
    """
    # Create output directories
    os.makedirs(output_dir, exist_ok=True)
    stats_dir = os.path.join(output_dir, "stats")
    os.makedirs(stats_dir, exist_ok=True)
    
    # Find VCF files
    vcf_files = []
    for file in os.listdir(vcf_dir):
        if file.endswith('.renamed.vcf.gz'):
            vcf_files.append(os.path.join(vcf_dir, file))
    
    if not vcf_files:
        print(f"No renamed VCF files found in {vcf_dir}")
        return []
    
    print(f"Found {len(vcf_files)} VCF files to annotate")
    
    # Process files in parallel
    results = []
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(
                annotate_vcf, vcf_file, output_dir, snpeff_dir, genome, stats_dir
            ): vcf_file for vcf_file in vcf_files
        }
        
        for future in as_completed(futures):
            basename, success, message = future.result()
            results.append((basename, success, message))
            status = "SUCCESS" if success else "FAILED"
            print(f"{basename}: {status} - {message}")
    
    return results

def create_annotation_summary(results, output_dir):
    """
    Create a summary report of annotation results.
    
    Args:
        results (list): List of (basename, success, message) tuples
        output_dir (str): Output directory
    """
    summary_file = os.path.join(output_dir, "annotation_summary.txt")
    
    with open(summary_file, 'w') as f:
        f.write("# SnpEff Annotation Summary\n")
        f.write(f"# Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # Overall statistics
        total = len(results)
        successful = sum(1 for _, success, _ in results if success)
        f.write(f"Total VCF files processed: {total}\n")
        f.write(f"Successfully annotated: {successful} ({successful/total*100:.1f}%)\n")
        f.write(f"Failed annotation: {total-successful} ({(total-successful)/total*100:.1f}%)\n\n")
        
        # Detailed results
        f.write("## Detailed Results\n\n")
        f.write("File\tStatus\tMessage\n")
        for basename, success, message in sorted(results):
            status = "SUCCESS" if success else "FAILED"
            f.write(f"{basename}\t{status}\t{message}\n")
    
    print(f"Created annotation summary: {summary_file}")

def main():
    parser = argparse.ArgumentParser(description="Batch annotate VCF files using SnpEff")
    parser.add_argument("--vcf_dir", required=True, help="Directory containing renamed VCF files")
    parser.add_argument("--output_dir", required=True, help="Output directory for annotated VCFs")
    parser.add_argument("--snpeff_dir", required=True, help="SnpEff installation directory")
    parser.add_argument("--genome", default="standard_w303", help="Genome name in SnpEff")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel processes")
    args = parser.parse_args()
    
    print(f"=== Starting Batch Annotation of VCF Files ===")
    print(f"Genome: {args.genome}")
    print(f"Using {args.threads} threads")
    
    start_time = time.time()
    
    # Run batch annotation
    results = batch_annotate(
        args.vcf_dir, args.output_dir, args.snpeff_dir, args.genome, args.threads
    )
    
    # Create summary report
    create_annotation_summary(results, args.output_dir)
    
    duration = time.time() - start_time
    successful = sum(1 for _, success, _ in results if success)
    total = len(results)
    
    print(f"\n=== Annotation Complete ===")
    print(f"Total time: {duration:.1f} seconds")
    print(f"Successfully annotated: {successful}/{total} VCF files")
    print(f"Annotated VCF files are in: {args.output_dir}")
    print(f"Annotation statistics are in: {os.path.join(args.output_dir, 'stats')}")

if __name__ == "__main__":
    main()