#!/usr/bin/env python3

# File: scripts/annotation/07b_debug_snpeff_chromosomes.py
# Purpose: Debug SnpEff chromosome naming issues in depth

import os
import subprocess
import re
import glob
from datetime import datetime

def main():
    print("=== Debugging SnpEff Chromosome Naming ===")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("")
    
    # Define directories
    snpeff_dir = "/Users/zakiralibhai/snpEff"
    data_dir = os.path.join(snpeff_dir, "data", "w303")
    debug_dir = "annotation/debug_chromosomes"
    
    # Create debug directory
    os.makedirs(debug_dir, exist_ok=True)
    
    # Step 1: Check if SnpEff database is built correctly
    print("Checking SnpEff database configuration...")
    
    # Check if genome exists in config
    config_file = os.path.join(snpeff_dir, "snpEff.config")
    
    if os.path.exists(config_file):
        print(f"- SnpEff config exists: {config_file}")
        
        # Extract w303 configuration
        try:
            with open(config_file, 'r') as f:
                config_text = f.read()
                
            w303_config = re.search(r'# w303.*?(?=\n#|\Z)', config_text, re.DOTALL)
            if w303_config:
                print("Found w303 configuration:")
                print(w303_config.group(0))
            else:
                print("w303 configuration not found in snpEff.config")
        except Exception as e:
            print(f"Error reading config file: {e}")
    else:
        print(f"SnpEff config not found: {config_file}")
    
    # Extract database information using SnpEff
    print("\nQuerying SnpEff database information...")
    try:
        output = subprocess.check_output(
            f"java -jar {os.path.join(snpeff_dir, 'snpEff.jar')} databases | grep -i w303",
            shell=True, text=True
        )
        print("SnpEff database output:")
        print(output)
    except subprocess.CalledProcessError:
        print("Error querying SnpEff databases")
    
    # Step 2: Check what chromosome names SnpEff is expecting
    print("\nExamining SnpEff chromosome names...")
    
    # Method 1: Look for chr names in sequence file
    sequence_file = os.path.join(data_dir, "sequences.fa")
    if os.path.exists(sequence_file):
        print(f"- Examining {sequence_file}...")
        try:
            output = subprocess.check_output(
                f"grep '>' {sequence_file} | head -10",
                shell=True, text=True
            )
            print("First 10 chromosome names in sequences.fa:")
            print(output)
            
            # Count total chromosomes
            count_output = subprocess.check_output(
                f"grep -c '>' {sequence_file}",
                shell=True, text=True
            )
            print(f"Total chromosomes in sequences.fa: {count_output.strip()}")
            
            # Save all chromosome names to file
            all_chrom_output = subprocess.check_output(
                f"grep '>' {sequence_file}",
                shell=True, text=True
            )
            all_chrom_file = os.path.join(debug_dir, "snpeff_chromosomes.txt")
            with open(all_chrom_file, 'w') as f:
                for line in all_chrom_output.splitlines():
                    chrom = line[1:].split()[0]  # Remove '>' and take first word
                    f.write(f"{chrom}\n")
            print(f"All chromosome names saved to: {all_chrom_file}")
        except subprocess.CalledProcessError:
            print("Error examining sequence file")
    else:
        print(f"Sequence file not found: {sequence_file}")
    
    # Step 3: Create a test VCF with the correct chromosome names
    print("\nCreating a test VCF with the correct chromosome names...")
    test_vcf_file = os.path.join(debug_dir, "test_variants.vcf")
    
    try:
        # Get all chromosome names
        with open(all_chrom_file, 'r') as f:
            chroms = [line.strip() for line in f if line.strip()]
        
        # Create a simple VCF with one variant per chromosome
        with open(test_vcf_file, 'w') as f:
            # Write VCF header
            f.write("##fileformat=VCFv4.2\n")
            f.write("##source=testVCF\n")
            for chrom in chroms:
                f.write(f"##contig=<ID={chrom}>\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            
            # Write one variant per chromosome
            for chrom in chroms:
                f.write(f"{chrom}\t100\t.\tA\tT\t100\tPASS\t.\n")
        
        print(f"Test VCF created: {test_vcf_file}")
        
        # Compress and index
        subprocess.run(f"bgzip -f {test_vcf_file}", shell=True)
        subprocess.run(f"tabix -p vcf {test_vcf_file}.gz", shell=True)
        
        # Run SnpEff on test VCF
        print("\nRunning SnpEff on test VCF...")
        test_output = os.path.join(debug_dir, "test_annotated.vcf")
        test_stats = os.path.join(debug_dir, "test_stats.html")
        
        cmd = f"java -jar {os.path.join(snpeff_dir, 'snpEff.jar')} -v -stats {test_stats} w303 {test_vcf_file}.gz > {test_output}"
        subprocess.run(cmd, shell=True)
        
        print(f"Test annotation complete: {test_output}")
        
        # Check for chromosome errors
        error_count = 0
        with open(test_output, 'r') as f:
            for line in f:
                if "ERROR_CHROMOSOME_NOT_FOUND" in line:
                    error_count += 1
        
        print(f"Found {error_count} chromosome errors in test annotation")
        
        if error_count == 0:
            print("✅ Test annotation successful! The chromosome names in the test VCF match what SnpEff expects.")
        else:
            print("❌ Test annotation failed with chromosome errors.")
    except Exception as e:
        print(f"Error creating or running test VCF: {e}")
    
    # Step 4: Generate SnpEff database rebuild recommendation
    print("\nGenerating recommendation for fixing the chromosome naming issue...")
    
    # Check which GenBank accession format is in the database
    genes_file = os.path.join(data_dir, "genes.gbk")
    jriu_in_genes = False
    
    if os.path.exists(genes_file):
        try:
            # Check first few lines for JRIU pattern
            with open(genes_file, 'r') as f:
                first_lines = ''.join([f.readline() for _ in range(100)])
                jriu_in_genes = bool(re.search(r'JRIU\d+', first_lines))
        except Exception:
            pass
    
    recommendation_file = os.path.join(debug_dir, "chromosome_fix_recommendation.txt")
    with open(recommendation_file, 'w') as f:
        f.write("SnpEff Chromosome Naming Issue Recommendation\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("==============================================\n\n")
        
        f.write("ISSUE DIAGNOSIS:\n")
        f.write("The chromosome names in your VCF files don't match the chromosome names in the SnpEff database.\n\n")
        
        f.write("RECOMMENDED SOLUTION OPTIONS:\n\n")
        
        f.write("OPTION 1: Modify the VCF files to use the chromosome names from the SnpEff database\n")
        f.write("- Extract a complete list of chromosome names from the SnpEff database\n")
        f.write("- Create a comprehensive mapping between your VCF chromosomes and SnpEff chromosomes\n")
        f.write("- Use bcftools to reheader the VCF files with the correct chromosome names\n\n")
        
        f.write("OPTION 2: Rebuild the SnpEff database using the GenBank files with JRIU identifiers\n")
        if jriu_in_genes:
            f.write("- The genes.gbk file already contains JRIU identifiers, but the database might not be properly built\n")
            f.write("- Delete the existing w303 database and rebuild it using the correct GenBank files\n")
        else:
            f.write("- The genes.gbk file doesn't contain JRIU identifiers\n")
            f.write("- Create new GenBank files using the JRIU identifiers from your VCF files\n")
            f.write("- Rebuild the SnpEff database using these modified GenBank files\n")
        f.write("\n")
        
        f.write("OPTION 3: Use a different annotation tool that can handle the chromosome name mismatch\n")
        f.write("- Tools like VEP or ANNOVAR might be more flexible with chromosome name mappings\n")
        f.write("- You could also use a tool like vcfanno with a custom annotation file\n\n")
        
        f.write("NEXT STEPS:\n")
        f.write("1. Run a detailed comparison between the exact chromosome names in your VCF files and the SnpEff database\n")
        f.write("2. If there's a consistent pattern, create a comprehensive mapping file\n")
        f.write("3. Use bcftools reheader to rename chromosomes in your VCF files\n")
        f.write("4. Try a test annotation on a single modified VCF to verify the fix\n")
    
    print(f"Recommendation saved to: {recommendation_file}")
    print("\nFinished debugging SnpEff chromosome naming issues.")
    print("")
    
    print("=== Debugging complete ===")

if __name__ == "__main__":
    main()
