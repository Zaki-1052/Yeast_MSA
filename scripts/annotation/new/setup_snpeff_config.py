#!/usr/bin/env python3
"""
setup_snpeff_config.py

This script configures SnpEff for proper chromosome name handling and tests annotation
on a sample VCF file to verify it works correctly.

Usage:
    python setup_snpeff_config.py --snpeff_dir <snpeff_directory> --synonym_file <synonym_file> --test_vcf <test_vcf_file>
"""

import os
import argparse
import subprocess
import shutil
import pandas as pd

def setup_chromosome_synonyms(snpeff_dir, synonym_file, genome_name="standard_w303"):
    """
    Set up chromosome synonyms file for SnpEff.
    
    Args:
        snpeff_dir (str): SnpEff installation directory
        synonym_file (str): Path to chromosome synonyms file
        genome_name (str): Name of the genome in SnpEff
    
    Returns:
        bool: Success status
    """
    # Check if files exist
    if not os.path.exists(synonym_file):
        print(f"ERROR: Synonym file {synonym_file} not found")
        return False
    
    # Read synonym file
    try:
        synonyms = pd.read_csv(synonym_file, sep='\t', header=None)
        print(f"Read {len(synonyms)} chromosome synonyms")
    except Exception as e:
        print(f"ERROR reading synonym file: {str(e)}")
        return False
    
    # Create target directory if it doesn't exist
    genome_dir = os.path.join(snpeff_dir, "data", genome_name)
    os.makedirs(genome_dir, exist_ok=True)
    
    # Copy synonym file to SnpEff directory
    target_file = os.path.join(genome_dir, "chromosome.synonym.txt")
    try:
        shutil.copy(synonym_file, target_file)
        print(f"Copied synonyms to {target_file}")
    except Exception as e:
        print(f"ERROR copying synonym file: {str(e)}")
        return False
    
    return True

def test_snpeff_annotation(snpeff_dir, test_vcf, genome_name="standard_w303"):
    """
    Test SnpEff annotation on a sample VCF file.
    
    Args:
        snpeff_dir (str): SnpEff installation directory
        test_vcf (str): Path to a sample VCF file
        genome_name (str): Name of the genome in SnpEff
    
    Returns:
        bool: Success status
    """
    if not os.path.exists(test_vcf):
        print(f"ERROR: Test VCF file {test_vcf} not found")
        return False
    
    # Create test output directory
    test_output_dir = os.path.join(os.path.dirname(test_vcf), "snpeff_test")
    os.makedirs(test_output_dir, exist_ok=True)
    
    # Output file path
    output_vcf = os.path.join(test_output_dir, "test_annotated.vcf")
    
    # Run SnpEff command
    cmd = [
        "java", "-Xmx4g", "-jar", os.path.join(snpeff_dir, "snpEff.jar"),
        "-v",
        "-stats", os.path.join(test_output_dir, "snpEff_summary.html"),
        genome_name,
        test_vcf,
        ">", output_vcf
    ]
    
    print("Running SnpEff test annotation...")
    print(f"Command: {' '.join(cmd)}")
    
    try:
        # Execute command
        process = subprocess.run(" ".join(cmd), shell=True, capture_output=True, text=True)
        
        # Check for errors
        if process.returncode != 0:
            print(f"ERROR: SnpEff command failed with return code {process.returncode}")
            print(f"STDERR: {process.stderr}")
            return False
        
        # Print output
        print(f"SnpEff output:\n{process.stdout}")
        
        # Check if output file exists and has content
        if not os.path.exists(output_vcf):
            print(f"ERROR: Output file {output_vcf} not created")
            return False
        
        file_size = os.path.getsize(output_vcf)
        if file_size == 0:
            print(f"ERROR: Output file {output_vcf} is empty")
            return False
        
        print(f"Test annotation successful. Output file: {output_vcf} (size: {file_size} bytes)")
        return True
        
    except Exception as e:
        print(f"ERROR executing SnpEff command: {str(e)}")
        return False

def check_genome_database(snpeff_dir, genome_name="standard_w303"):
    """
    Check if the genome database exists and is properly configured.
    
    Args:
        snpeff_dir (str): SnpEff installation directory
        genome_name (str): Name of the genome in SnpEff
    
    Returns:
        bool: Whether the database exists
    """
    # Check if database directory exists
    database_dir = os.path.join(snpeff_dir, "data", genome_name)
    if not os.path.exists(database_dir):
        print(f"WARNING: Database directory {database_dir} does not exist")
        return False
    
    # Check for essential files
    essential_files = ["genes.gbk", "snpEffectPredictor.bin"]
    missing_files = []
    
    for file in essential_files:
        file_path = os.path.join(database_dir, file)
        if not os.path.exists(file_path):
            missing_files.append(file)
    
    if missing_files:
        print(f"WARNING: Missing essential files: {', '.join(missing_files)}")
        return False
    
    # Verify that the database is in snpEff.config
    config_file = os.path.join(snpeff_dir, "snpEff.config")
    if not os.path.exists(config_file):
        print(f"WARNING: snpEff.config not found at {config_file}")
        return False
    
    # Check if genome is in config
    found_in_config = False
    with open(config_file, 'r') as f:
        for line in f:
            if line.startswith(genome_name + ".genome"):
                found_in_config = True
                break
    
    if not found_in_config:
        print(f"WARNING: Genome {genome_name} not found in snpEff.config")
        print("You may need to add the following line to snpEff.config:")
        print(f"{genome_name}.genome : W303")
        return False
    
    print(f"Genome database {genome_name} exists and appears to be properly configured")
    return True

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Configure SnpEff for proper chromosome name handling")
    parser.add_argument("--snpeff_dir", required=True, help="SnpEff installation directory")
    parser.add_argument("--synonym_file", required=True, help="Path to chromosome synonyms file")
    parser.add_argument("--test_vcf", required=True, help="Path to a sample VCF file for testing")
    parser.add_argument("--genome", default="standard_w303", help="Name of the genome in SnpEff")
    args = parser.parse_args()
    
    # Check if SnpEff directory exists
    if not os.path.exists(args.snpeff_dir):
        print(f"ERROR: SnpEff directory {args.snpeff_dir} not found")
        return
    
    # Check if the genome database exists
    db_exists = check_genome_database(args.snpeff_dir, args.genome)
    if not db_exists:
        print("WARNING: SnpEff database not properly configured")
        print("Please make sure you have built the SnpEff database for this genome")
        print("Example command:")
        print(f"  java -Xmx4g -jar {os.path.join(args.snpeff_dir, 'snpEff.jar')} build -genbank -v {args.genome}")
        return
    
    # Set up chromosome synonyms
    if not setup_chromosome_synonyms(args.snpeff_dir, args.synonym_file, args.genome):
        print("ERROR: Failed to set up chromosome synonyms")
        return
    
    # Test annotation
    test_success = test_snpeff_annotation(args.snpeff_dir, args.test_vcf, args.genome)
    
    if test_success:
        print("\nSuccess! SnpEff is properly configured for annotation.")
        print("\nYou can now use the following command to annotate your VCF files:")
        print(f"  java -Xmx4g -jar {os.path.join(args.snpeff_dir, 'snpEff.jar')} {args.genome} [input.vcf] > [output.vcf]")
    else:
        print("\nWARNING: SnpEff test annotation failed. Please review the errors above.")
        print("You may need to check:")
        print("  1. SnpEff database is built correctly")
        print("  2. Chromosome synonyms are set up properly")
        print("  3. Input VCF file format is valid")

if __name__ == "__main__":
    main()