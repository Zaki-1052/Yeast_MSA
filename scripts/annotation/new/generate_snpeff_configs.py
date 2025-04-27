#!/usr/bin/env python3
"""
generate_snpeff_configs.py

This script generates alternative SnpEff configurations to test chromosome synonym handling.
"""

import os
import argparse
import pandas as pd
import shutil

def read_chromosome_mapping(mapping_file):
    """Read chromosome mapping file."""
    df = pd.read_csv(mapping_file, sep='\t')
    return df

def create_config_options(snpeff_dir, mapping_file, genome_name="standard_w303"):
    """Create different config options to test."""
    print(f"Reading mapping from {mapping_file}")
    mapping = read_chromosome_mapping(mapping_file)
    
    # Backup original config
    config_path = os.path.join(snpeff_dir, "snpEff.config")
    backup_path = os.path.join(snpeff_dir, "snpEff.config.backup")
    
    if not os.path.exists(backup_path):
        print(f"Creating backup of {config_path}")
        shutil.copy(config_path, backup_path)
    
    # Read original config
    with open(config_path, 'r') as f:
        config_content = f.read()
    
    # Create directory for alternative configs
    configs_dir = os.path.join(os.path.dirname(mapping_file), "configs")
    os.makedirs(configs_dir, exist_ok=True)
    
    # Create options
    print("Generating configuration options...")
    
    # Option 1: Add explicit chromosome list
    chromosome_list = ", ".join(mapping['w303_scaffold'])
    option1 = f"\n# Explicit chromosome list for {genome_name}\n{genome_name}.chromosomes = {chromosome_list}\n"
    
    # Option 2: Add explicit synonyms for each chromosome
    option2 = "\n# Explicit chromosome synonyms\n"
    for _, row in mapping.iterrows():
        option2 += f"{genome_name}.{row['w303_scaffold']}.synonyms = {row['chromosome_id']}\n"
    
    # Save configurations
    options = {
        "chromosome_list": option1,
        "explicit_synonyms": option2,
        "combined": option1 + option2
    }
    
    for name, content in options.items():
        option_path = os.path.join(configs_dir, f"snpEff_option_{name}.txt")
        with open(option_path, 'w') as f:
            f.write(content)
        print(f"Created option file: {option_path}")
    
    # Save test configurations that could be appended to snpEff.config
    return options

def create_alternative_synonym_files(mapping_file, output_dir):
    """Create alternative formats for synonym files to test."""
    print(f"\nCreating alternative synonym file formats")
    
    # Read mapping
    mapping = read_chromosome_mapping(mapping_file)
    
    # Create alternatives directory
    alt_dir = os.path.join(output_dir, "alternative_synonyms")
    os.makedirs(alt_dir, exist_ok=True)
    
    # Alternative 1: With header
    alt1_path = os.path.join(alt_dir, "synonyms_with_header.txt")
    with open(alt1_path, 'w') as f:
        f.write("#alias\tname\n")
        for _, row in mapping.iterrows():
            f.write(f"{row['chromosome_id']}\t{row['w303_scaffold']}\n")
    print(f"Created alternative 1 (with header): {alt1_path}")
    
    # Alternative 2: Reversed order
    alt2_path = os.path.join(alt_dir, "synonyms_reversed.txt")
    with open(alt2_path, 'w') as f:
        for _, row in mapping.iterrows():
            f.write(f"{row['w303_scaffold']}\t{row['chromosome_id']}\n")
    print(f"Created alternative 2 (reversed order): {alt2_path}")
    
    # Alternative 3: Multiple synonyms per line
    alt3_path = os.path.join(alt_dir, "synonyms_multiple.txt")
    chrom_dict = {}
    for _, row in mapping.iterrows():
        key = row['w303_scaffold']
        value = row['chromosome_id']
        chrom_dict.setdefault(key, []).append(value)
    
    with open(alt3_path, 'w') as f:
        for key, values in chrom_dict.items():
            f.write(f"{key}\t{', '.join(values)}\n")
    print(f"Created alternative 3 (multiple synonyms): {alt3_path}")
    
    return {
        "with_header": alt1_path,
        "reversed": alt2_path,
        "multiple": alt3_path
    }

def main():
    parser = argparse.ArgumentParser(description="Generate alternative SnpEff configurations")
    parser.add_argument("--snpeff_dir", required=True, help="SnpEff installation directory")
    parser.add_argument("--mapping_file", required=True, help="Chromosome mapping file (tsv)")
    parser.add_argument("--output_dir", required=True, help="Output directory for generated files")
    parser.add_argument("--genome", default="standard_w303", help="Genome name in SnpEff")
    args = parser.parse_args()
    
    # Create config options
    config_options = create_config_options(args.snpeff_dir, args.mapping_file, args.genome)
    
    # Create alternative synonym files
    synonym_options = create_alternative_synonym_files(args.mapping_file, args.output_dir)
    
    print("\n=== Test Implementation Instructions ===")
    print("1. First, test with debug script:")
    print(f"   python debug_snpeff_synonyms.py --snpeff_dir {args.snpeff_dir} --test_vcf <your_vcf> --genome {args.genome}")
    
    print("\n2. Try appending config options to snpEff.config:")
    print(f"   cat {os.path.join(args.output_dir, 'configs/snpEff_option_combined.txt')} >> {os.path.join(args.snpeff_dir, 'snpEff.config')}")
    
    print("\n3. Try alternative synonym files:")
    print(f"   cp {synonym_options['reversed']} {os.path.join(args.snpeff_dir, 'data', args.genome, 'chromosome.synonym.txt')}")
    
    print("\n4. After testing, restore original config:")
    print(f"   cp {os.path.join(args.snpeff_dir, 'snpEff.config.backup')} {os.path.join(args.snpeff_dir, 'snpEff.config')}")

if __name__ == "__main__":
    main()