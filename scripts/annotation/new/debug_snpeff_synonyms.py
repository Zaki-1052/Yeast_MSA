#!/usr/bin/env python3
"""
debug_snpeff_synonyms.py

This script diagnoses issues with SnpEff chromosome synonyms configuration.
"""

import os
import subprocess
import re

def check_synonym_file(snpeff_dir, genome_name="standard_w303"):
    """Check if the synonym file exists and has the correct format."""
    synonym_path = os.path.join(snpeff_dir, "data", genome_name, "chromosome.synonym.txt")
    
    print(f"Checking synonym file: {synonym_path}")
    
    if not os.path.exists(synonym_path):
        print("ERROR: Synonym file not found!")
        return False
    
    # Check file permissions
    if not os.access(synonym_path, os.R_OK):
        print("ERROR: Synonym file exists but is not readable!")
        return False
    
    # Check file format
    with open(synonym_path, 'r') as f:
        lines = f.readlines()
        
    if not lines:
        print("ERROR: Synonym file is empty!")
        return False
    
    print(f"Synonym file contains {len(lines)} lines")
    
    # Check first few lines
    print("First few lines:")
    for i, line in enumerate(lines[:5]):
        print(f"  {i+1}: {line.strip()}")
    
    # Check if format is as expected (two columns, tab-separated)
    for i, line in enumerate(lines):
        parts = line.strip().split('\t')
        if len(parts) != 2:
            print(f"ERROR: Line {i+1} does not have exactly 2 tab-separated columns: {line.strip()}")
            return False
    
    print("Synonym file format appears correct (tab-separated, 2 columns)")
    return True

def check_snpeff_config(snpeff_dir, genome_name="standard_w303"):
    """Check if the SnpEff config file references our genome correctly."""
    config_path = os.path.join(snpeff_dir, "snpEff.config")
    
    print(f"\nChecking SnpEff config: {config_path}")
    
    if not os.path.exists(config_path):
        print("ERROR: SnpEff config file not found!")
        return False
    
    with open(config_path, 'r') as f:
        config = f.read()
    
    # Check for genome entry
    genome_pattern = re.compile(f"{re.escape(genome_name)}\\.genome\\s*:\\s*(\\S+)")
    match = genome_pattern.search(config)
    
    if not match:
        print(f"ERROR: No entry for {genome_name} found in config file!")
        return False
    
    print(f"Found genome entry: {genome_name}.genome : {match.group(1)}")
    
    # Check for chromosome entries
    chrom_pattern = re.compile(f"{re.escape(genome_name)}\\.chromosomes\\s*=")
    if chrom_pattern.search(config):
        print(f"Found chromosomes entry for {genome_name}")
    else:
        print(f"NOTE: No explicit chromosomes list found for {genome_name}")
        print("  This is not necessarily an error, but adding one might help")
    
    # Check for explicit synonyms
    synonym_pattern = re.compile(f"{re.escape(genome_name)}\\.\\S+\\.synonyms\\s*=")
    synonyms_found = [m.group(0) for m in synonym_pattern.finditer(config)]
    
    if synonyms_found:
        print(f"Found {len(synonyms_found)} explicit synonym entries:")
        for s in synonyms_found[:3]:
            print(f"  {s}")
        if len(synonyms_found) > 3:
            print(f"  ... and {len(synonyms_found) - 3} more")
    else:
        print("NOTE: No explicit synonym entries found in config")
        print("  This might be fine if the chromosome.synonym.txt file is properly read")
    
    return True

def run_snpeff_test(snpeff_dir, test_vcf, genome_name="standard_w303"):
    """Run SnpEff with debug flags to see detailed output."""
    print(f"\nRunning SnpEff with debug flags on {test_vcf}")
    
    output_dir = os.path.dirname(test_vcf)
    output_vcf = os.path.join(output_dir, "debug_test.vcf")
    
    # Run SnpEff with debug flags
    cmd = [
        "java", "-Xmx4g", "-jar", os.path.join(snpeff_dir, "snpEff.jar"),
        "-v",  # Verbose
        "-debug",  # Debug output
        "-csvStats", os.path.join(output_dir, "debug_stats.csv"),
        genome_name,
        test_vcf
    ]
    
    print(f"Command: {' '.join(cmd)}")
    print("Running (this may take a moment)...")
    
    try:
        output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True)
        
        # Save the full output to a file
        debug_log = os.path.join(output_dir, "snpeff_debug.log")
        with open(debug_log, 'w') as f:
            f.write(output)
        
        print(f"Full debug output saved to: {debug_log}")
        
        # Look for specific synonym-related log entries
        synonym_related = []
        for line in output.splitlines():
            if "synonym" in line.lower() or "chromosome" in line.lower():
                synonym_related.append(line)
        
        if synonym_related:
            print("\nChromosome/Synonym related log entries:")
            for line in synonym_related[:20]:  # Limit to first 20 entries
                print(f"  {line}")
            if len(synonym_related) > 20:
                print(f"  ... and {len(synonym_related) - 20} more entries")
        else:
            print("\nNo specific chromosome/synonym entries found in the log")
        
        # Check if the error persists
        with open(debug_log, 'r') as f:
            if "ERROR_CHROMOSOME_NOT_FOUND" in f.read():
                print("\nERROR: Chromosome not found issue still persists!")
                print("Let's analyze why...")
            else:
                print("\nSuccess! No 'ERROR_CHROMOSOME_NOT_FOUND' detected")
        
        return True
    
    except subprocess.CalledProcessError as e:
        print(f"ERROR: SnpEff command failed with exit code {e.returncode}")
        print(f"Output: {e.output}")
        return False

def create_minimal_test_vcf(input_vcf, output_vcf):
    """Create a minimal test VCF with just one variant."""
    print(f"\nCreating minimal test VCF from {input_vcf}")
    
    # Determine if input is gzipped
    is_gzipped = input_vcf.endswith('.gz')
    
    # Command to extract header and first variant
    if is_gzipped:
        header_cmd = f"zcat {input_vcf} | grep '^#' > {output_vcf}"
        first_variant_cmd = f"zcat {input_vcf} | grep -v '^#' | head -1 >> {output_vcf}"
    else:
        header_cmd = f"grep '^#' {input_vcf} > {output_vcf}"
        first_variant_cmd = f"grep -v '^#' {input_vcf} | head -1 >> {output_vcf}"
    
    try:
        subprocess.run(header_cmd, shell=True, check=True)
        subprocess.run(first_variant_cmd, shell=True, check=True)
        
        # Verify the file was created
        if os.path.exists(output_vcf):
            print(f"Created minimal test VCF: {output_vcf}")
            
            # Print the content
            with open(output_vcf, 'r') as f:
                first_variant = next((line for line in f if not line.startswith('#')), None)
                
            if first_variant:
                print(f"First variant: {first_variant.strip()}")
                chrom = first_variant.split('\t')[0]
                print(f"Chromosome: {chrom}")
            else:
                print("WARNING: No variant found in the test VCF")
                
            return True
        else:
            print(f"ERROR: Failed to create {output_vcf}")
            return False
    
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Command failed: {e}")
        return False

def suggest_next_steps():
    """Suggest possible next steps based on our findings."""
    print("\n=== Possible Next Steps ===")
    print("1. Add explicit chromosome synonyms to snpEff.config:")
    print("   standard_w303.chromosomes = w303_scaffold_1, w303_scaffold_2, ...")
    print("   standard_w303.w303_scaffold_1.synonyms = CM007964.1")
    print("   standard_w303.w303_scaffold_2.synonyms = CM007965.1")
    print("   ...")
    
    print("\n2. Try a different synonym file format:")
    print("   # Original")
    print("   CM007964.1\tw303_scaffold_1")
    print("   # Alternative 1 (with header)")
    print("   #alias\tname")
    print("   CM007964.1\tw303_scaffold_1")
    print("   # Alternative 2 (reversed order)")
    print("   w303_scaffold_1\tCM007964.1")
    
    print("\n3. If the above fails, consider using bcftools for chromosome renaming:")
    print("   bcftools annotate --rename-chrs reference/chromosome_map.txt input.vcf > output.vcf")

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Debug SnpEff chromosome synonyms configuration")
    parser.add_argument("--snpeff_dir", required=True, help="SnpEff installation directory")
    parser.add_argument("--test_vcf", required=True, help="VCF file to test")
    parser.add_argument("--genome", default="standard_w303", help="Genome name in SnpEff")
    args = parser.parse_args()
    
    print("=== SnpEff Chromosome Synonyms Debugging ===")
    
    # Check synonym file
    check_synonym_file(args.snpeff_dir, args.genome)
    
    # Check SnpEff config
    check_snpeff_config(args.snpeff_dir, args.genome)
    
    # Create minimal test VCF
    test_dir = os.path.dirname(args.test_vcf)
    minimal_vcf = os.path.join(test_dir, "minimal_test.vcf")
    create_minimal_test_vcf(args.test_vcf, minimal_vcf)
    
    # Run SnpEff test
    run_snpeff_test(args.snpeff_dir, minimal_vcf, args.genome)
    
    # Suggest next steps
    suggest_next_steps()

if __name__ == "__main__":
    main()