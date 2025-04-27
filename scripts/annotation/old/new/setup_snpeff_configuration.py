# File: scripts/annotation/new/setup_snpeff_configuration.py

import os
import shutil

def setup_snpeff_configuration():
    """Set up a custom SnpEff configuration that properly handles chromosome structure."""
    
    # Directory setup
    snpeff_dir = os.path.expanduser("~/snpEff")
    data_dir = os.path.join(snpeff_dir, "data", "w303_fixed")
    os.makedirs(data_dir, exist_ok=True)
    
    # Read chromosome summary
    chromosome_info = {}
    fragments = []
    circular_chromosomes = []
    
    with open("reference/chromosome_summary.tsv", 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if i == 0:  # Skip header
                continue
            
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                chrom_id, roman_numeral, scaffold, is_fragment, is_circular = parts
                
                chromosome_info[chrom_id] = {
                    'roman_numeral': roman_numeral,
                    'scaffold': scaffold,
                    'is_fragment': is_fragment,
                    'is_circular': is_circular
                }
                
                if is_fragment == "Yes":
                    fragments.append(chrom_id)
                
                if is_circular == "Yes":
                    circular_chromosomes.append(chrom_id)
    
    # Create genome.yaml file
    genome_yaml = os.path.join(data_dir, "genome.yaml")
    with open(genome_yaml, 'w') as f:
        f.write("---\n")
        f.write("name: w303\n")
        f.write("id: w303_fixed\n")
        f.write("circular_chromosomes:\n")
        for chrom in circular_chromosomes:
            f.write(f"  - {chrom}\n")
        f.write("\n")
        f.write("chromosome_alias:\n")
        # Map VCF chromosome names to SnpEff chromosome names
        for chrom_id, info in chromosome_info.items():
            if chrom_id not in fragments:  # Skip fragments
                f.write(f"  {chrom_id}: {chrom_id}\n")
                # Also add alias without .1 suffix if needed
                if '.' in chrom_id:
                    base_id = chrom_id.split('.')[0]
                    f.write(f"  {base_id}: {chrom_id}\n")
    
    # Copy GenBank files (excluding fragments)
    source_dir = os.path.join(snpeff_dir, "data", "w303_cm", "genes")
    genes_dir = os.path.join(data_dir, "genes")
    os.makedirs(genes_dir, exist_ok=True)
    
    for filename in os.listdir(source_dir):
        if filename.endswith('.gbk'):
            # Skip fragment chromosomes
            if any(fragment in filename for fragment in fragments):
                print(f"Skipping fragment file: {filename}")
                continue
            
            source_file = os.path.join(source_dir, filename)
            dest_file = os.path.join(genes_dir, filename)
            shutil.copy2(source_file, dest_file)
            print(f"Copied {filename} to {genes_dir}")
    
    # Create a combined genbank file
    combined_file = os.path.join(data_dir, "genes.gbk")
    with open(combined_file, 'w') as outfile:
        for filename in sorted(os.listdir(genes_dir)):
            if filename.endswith('.gbk'):
                with open(os.path.join(genes_dir, filename), 'r') as infile:
                    outfile.write(infile.read())
                    outfile.write("\n\n")
    print(f"Created combined GenBank file: {combined_file}")
    
    # Update SnpEff config file
    config_file = os.path.join(snpeff_dir, "snpEff.config")
    config_updated = False
    
    # Check if w303_fixed already exists in config
    with open(config_file, 'r') as f:
        config_content = f.read()
        if "w303_fixed.genome" in config_content:
            config_updated = True
    
    if not config_updated:
        with open(config_file, 'a') as f:
            f.write("\n# W303 genome with fixed configuration\n")
            f.write("w303_fixed.genome : Saccharomyces cerevisiae W303 (fixed)\n")
    
    # Create a build script
    build_script = "build_snpeff_w303_fixed.sh"
    with open(build_script, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("cd ~/snpEff\n")
        f.write("java -jar snpEff.jar build -genbank -v w303_fixed\n")
    
    # Make the script executable
    os.chmod(build_script, 0o755)
    
    # Create a test annotation script
    test_script = "test_snpeff_annotation.sh"
    with open(test_script, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("cd ~/snpEff\n")
        f.write("# Test annotation on a sample VCF file\n")
        f.write("java -jar snpEff.jar ann w303_fixed -v ~/Documents/GitHub/Yeast_MSA/vcf/merged/filtered/WT-CTRL.filtered.vcf.gz > ~/Documents/GitHub/Yeast_MSA/vcf/merged/filtered/WT-CTRL.filtered.ann.vcf\n")
    
    # Make the test script executable
    os.chmod(test_script, 0o755)
    
    print("\nSetup complete!")
    print(f"SnpEff configuration created in: {data_dir}")
    print(f"Build script created: {build_script}")
    print(f"Test script created: {test_script}")
    print("\nNext steps:")
    print(f"1. Run the build script: ./{build_script}")
    print(f"2. Test annotation with: ./{test_script}")

if __name__ == "__main__":
    setup_snpeff_configuration()