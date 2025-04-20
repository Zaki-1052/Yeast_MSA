# File: scripts/annotation/new/setup_snpeff_cm_genome.py

import os
import shutil
import subprocess

def setup_snpeff_cm_genome():
    """Set up SnpEff for the w303_cm genome with CM chromosome names."""
    snpeff_dir = os.path.expanduser("~/snpEff")
    
    # 1. Create directories
    data_dir = os.path.join(snpeff_dir, "data", "w303_cm")
    genes_dir = os.path.join(data_dir, "genes")
    
    os.makedirs(genes_dir, exist_ok=True)
    print(f"Created directory structure: {genes_dir}")
    
    # 2. Copy modified GenBank files
    source_dir = "reference/w303_annotations/genbank_cm_names"
    gbk_files = [f for f in os.listdir(source_dir) if f.endswith('.gbk')]
    
    for gbk_file in gbk_files:
        source_path = os.path.join(source_dir, gbk_file)
        dest_path = os.path.join(genes_dir, gbk_file)
        shutil.copy2(source_path, dest_path)
    
    print(f"Copied {len(gbk_files)} GenBank files to {genes_dir}")
    
    # 3. Update SnpEff config
    config_path = os.path.join(snpeff_dir, "snpEff.config")
    
    # Check if w303_cm entry already exists
    config_updated = False
    if os.path.exists(config_path):
        with open(config_path, 'r') as f:
            config_content = f.read()
            if "w303_cm" in config_content:
                config_updated = True
                print("w303_cm genome already exists in SnpEff config")
    
    # Add w303_cm entry if needed
    if not config_updated:
        with open(config_path, 'a') as f:
            f.write("\n# W303 genome with CM chromosome names\n")
            f.write("w303_cm.genome : Saccharomyces cerevisiae W303\n")
        print("Added w303_cm genome to SnpEff config")
    
    # 4. Create a build script
    build_script = "build_snpeff_w303cm.sh"
    with open(build_script, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("cd ~/snpEff\n")
        f.write("java -jar snpEff.jar build -genbank -v w303_cm\n")
    
    # Make the script executable
    os.chmod(build_script, 0o755)
    print(f"Created build script: {build_script}")
    
    print("\nNext steps:")
    print(f"1. Run the build script: ./{build_script}")
    print("2. Once built, test SnpEff annotation with a sample VCF file")

if __name__ == "__main__":
    setup_snpeff_cm_genome()