#!/usr/bin/env python3
"""
test_bcftools_rename.py

This script tests chromosome renaming using bcftools on a sample VCF file.
"""

import os
import argparse
import subprocess
import shutil

def test_bcftools_rename(input_vcf, mapping_file, output_dir):
    """Test chromosome renaming with bcftools."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Get base filename
    basename = os.path.basename(input_vcf)
    output_vcf = os.path.join(output_dir, basename.replace('.vcf.gz', '.renamed.vcf.gz'))
    
    print(f"Renaming chromosomes in {basename}...")
    
    # Run bcftools rename
    cmd = f"bcftools annotate --rename-chrs {mapping_file} {input_vcf} -Oz -o {output_vcf}"
    print(f"Command: {cmd}")
    
    try:
        subprocess.run(cmd, shell=True, check=True)
        
        # Index the output VCF
        subprocess.run(f"bcftools index {output_vcf}", shell=True, check=True)
        
        print(f"Successfully renamed and indexed: {output_vcf}")
        
        # Return path to renamed VCF
        return output_vcf
    
    except subprocess.CalledProcessError as e:
        print(f"ERROR: bcftools command failed: {e}")
        return None

def test_snpeff_annotation(snpeff_dir, input_vcf, genome_name="standard_w303", output_dir=None):
    """Test SnpEff annotation on the renamed VCF."""
    if output_dir is None:
        output_dir = os.path.dirname(input_vcf)
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Get base filename
    basename = os.path.basename(input_vcf)
    output_vcf = os.path.join(output_dir, basename.replace('.renamed.vcf.gz', '.annotated.vcf'))
    
    print(f"\nTesting SnpEff annotation on renamed VCF...")
    
    # Run SnpEff
    cmd = f"java -Xmx4g -jar {os.path.join(snpeff_dir, 'snpEff.jar')} -v {genome_name} {input_vcf} > {output_vcf}"
    print(f"Command: {cmd}")
    
    try:
        subprocess.run(cmd, shell=True, check=True)
        
        print(f"Successfully created {output_vcf}")
        
        # Check for ERROR_CHROMOSOME_NOT_FOUND
        with open(output_vcf, 'r') as f:
            content = f.read()
            if "ERROR_CHROMOSOME_NOT_FOUND" in content:
                print("ERROR: Still seeing ERROR_CHROMOSOME_NOT_FOUND in output")
                return False
            else:
                print("SUCCESS: No ERROR_CHROMOSOME_NOT_FOUND found in output")
                return True
    
    except subprocess.CalledProcessError as e:
        print(f"ERROR: SnpEff command failed: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Test chromosome renaming using bcftools")
    parser.add_argument("--input", required=True, help="Input VCF file")
    parser.add_argument("--mapping", required=True, help="bcftools mapping file")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--snpeff_dir", required=True, help="SnpEff installation directory")
    parser.add_argument("--genome", default="standard_w303", help="Genome name in SnpEff")
    args = parser.parse_args()
    
    # Test bcftools rename
    renamed_vcf = test_bcftools_rename(args.input, args.mapping, args.output_dir)
    
    if renamed_vcf:
        # Test SnpEff annotation
        test_snpeff_annotation(args.snpeff_dir, renamed_vcf, args.genome, args.output_dir)

if __name__ == "__main__":
    main()