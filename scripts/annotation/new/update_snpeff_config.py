#!/usr/bin/env python3
"""
update_snpeff_config.py

This script adds explicit chromosome and synonym entries to the SnpEff configuration.
"""

import os
import argparse
import pandas as pd
import shutil

def read_chromosome_mapping(mapping_file):
    """Read chromosome mapping file."""
    df = pd.read_csv(mapping_file, sep='\t')
    return df

def update_snpeff_config(snpeff_dir, mapping_file, genome_name="standard_w303"):
    """Update SnpEff config with explicit chromosome and synonym entries."""
    print(f"Reading mapping from {mapping_file}")
    mapping = read_chromosome_mapping(mapping_file)
    
    # Backup original config
    config_path = os.path.join(snpeff_dir, "snpEff.config")
    backup_path = os.path.join(snpeff_dir, "snpEff.config.backup")
    
    if not os.path.exists(backup_path):
        print(f"Creating backup of {config_path}")
        shutil.copy(config_path, backup_path)
    else:
        print(f"Backup already exists at {backup_path}")
    
    # Read original config
    with open(config_path, 'r') as f:
        config_content = f.read()
    
    # Check if we've already added chromosomes to avoid duplicates
    if f"{genome_name}.chromosomes" in config_content:
        print(f"WARNING: {genome_name}.chromosomes already exists in config file")
        print("Would you like to proceed anyway? (y/n)")
        response = input().strip().lower()
        if response != 'y':
            print("Operation cancelled")
            return False
    
    # Create chromosome list and synonym entries
    chromosome_list = ", ".join(mapping['w303_scaffold'])
    additions = f"\n# Explicit chromosome list for {genome_name}\n{genome_name}.chromosomes = {chromosome_list}\n\n"
    additions += f"# Explicit chromosome synonyms for {genome_name}\n"
    
    for _, row in mapping.iterrows():
        additions += f"{genome_name}.{row['w303_scaffold']}.synonyms = {row['chromosome_id']}\n"
    
    # Append to config file
    print(f"Appending {len(additions.splitlines())} lines to SnpEff config")
    with open(config_path, 'a') as f:
        f.write(additions)
    
    print(f"Successfully updated {config_path}")
    return True

def create_minimal_vcf(input_vcf, output_vcf):
    """Create a minimal test VCF using gunzip for Mac compatibility."""
    print(f"Creating minimal test VCF from {input_vcf}")
    
    if not os.path.exists(input_vcf):
        print(f"ERROR: Input VCF {input_vcf} not found")
        return False
    
    # Determine correct command for macOS
    if input_vcf.endswith('.gz'):
        # Use gunzip -c instead of zcat for better compatibility
        header_cmd = f"gunzip -c {input_vcf} | grep '^#' > {output_vcf}"
        variant_cmd = f"gunzip -c {input_vcf} | grep -v '^#' | head -1 >> {output_vcf}"
    else:
        header_cmd = f"grep '^#' {input_vcf} > {output_vcf}"
        variant_cmd = f"grep -v '^#' {input_vcf} | head -1 >> {output_vcf}"
    
    try:
        print(f"Extracting headers: {header_cmd}")
        os.system(header_cmd)
        
        print(f"Extracting first variant: {variant_cmd}")
        os.system(variant_cmd)
        
        if os.path.exists(output_vcf):
            print(f"Successfully created {output_vcf}")
            
            # Verify content
            with open(output_vcf, 'r') as f:
                lines = f.readlines()
                header_count = sum(1 for line in lines if line.startswith('#'))
                variant_count = len(lines) - header_count
            
            print(f"  {header_count} header lines")
            print(f"  {variant_count} variant lines")
            
            return True
        else:
            print(f"ERROR: Failed to create {output_vcf}")
            return False
    
    except Exception as e:
        print(f"ERROR: {str(e)}")
        return False

def test_updated_config(snpeff_dir, test_vcf, genome_name="standard_w303"):
    """Test the updated configuration."""
    print(f"\nTesting updated configuration with {test_vcf}")
    
    output_vcf = os.path.join(os.path.dirname(test_vcf), "test_annotated.vcf")
    
    cmd = f"java -Xmx4g -jar {os.path.join(snpeff_dir, 'snpEff.jar')} -v {genome_name} {test_vcf} > {output_vcf}"
    print(f"Command: {cmd}")
    
    try:
        os.system(cmd)
        
        if os.path.exists(output_vcf):
            print(f"Successfully created {output_vcf}")
            
            # Check for ERROR_CHROMOSOME_NOT_FOUND
            with open(output_vcf, 'r') as f:
                content = f.read()
                if "ERROR_CHROMOSOME_NOT_FOUND" in content:
                    print("ERROR: Still seeing ERROR_CHROMOSOME_NOT_FOUND in output")
                    print("Configuration update did not resolve the issue")
                    return False
                else:
                    print("SUCCESS: No ERROR_CHROMOSOME_NOT_FOUND found in output")
                    print("Configuration update appears to have resolved the issue")
                    return True
        else:
            print(f"ERROR: Failed to create {output_vcf}")
            return False
    
    except Exception as e:
        print(f"ERROR: {str(e)}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Update SnpEff config with explicit chromosome and synonym entries")
    parser.add_argument("--snpeff_dir", required=True, help="SnpEff installation directory")
    parser.add_argument("--mapping_file", required=True, help="Chromosome mapping file (tsv)")
    parser.add_argument("--test_vcf", required=True, help="VCF file to test")
    parser.add_argument("--genome", default="standard_w303", help="Genome name in SnpEff")
    args = parser.parse_args()
    
    # Update SnpEff config
    success = update_snpeff_config(args.snpeff_dir, args.mapping_file, args.genome)
    
    if not success:
        print("Failed to update SnpEff config, exiting")
        return
    
    # Create minimal test VCF
    test_dir = os.path.dirname(args.test_vcf)
    minimal_vcf = os.path.join(test_dir, "minimal_test.vcf")
    create_minimal_vcf(args.test_vcf, minimal_vcf)
    
    # Test updated configuration
    test_updated_config(args.snpeff_dir, minimal_vcf, args.genome)
    
    print("\nIf testing was successful, you can proceed with annotation")
    print("If not, you can restore the original config:")
    print(f"  cp {os.path.join(args.snpeff_dir, 'snpEff.config.backup')} {os.path.join(args.snpeff_dir, 'snpEff.config')}")

if __name__ == "__main__":
    main()