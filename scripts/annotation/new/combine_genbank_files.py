# File: scripts/annotation/new/combine_genbank_files.py

import os
import shutil

def combine_genbank_files():
    """Combine all GenBank files for the w303_cm genome."""
    snpeff_dir = os.path.expanduser("~/snpEff")
    genes_dir = os.path.join(snpeff_dir, "data", "w303_cm", "genes")
    combined_path = os.path.join(snpeff_dir, "data", "w303_cm", "genes.gbk")
    
    # Check if genes directory exists and contains files
    if not os.path.exists(genes_dir):
        print(f"Error: Genes directory does not exist: {genes_dir}")
        return False
    
    gbk_files = [f for f in os.listdir(genes_dir) if f.endswith('.gbk')]
    if not gbk_files:
        print(f"Error: No GenBank files found in {genes_dir}")
        return False
    
    print(f"Combining {len(gbk_files)} GenBank files into {combined_path}")
    
    # Combine the GenBank files
    with open(combined_path, 'w') as outfile:
        for i, gbk_file in enumerate(sorted(gbk_files)):
            file_path = os.path.join(genes_dir, gbk_file)
            with open(file_path, 'r') as infile:
                content = infile.read()
                outfile.write(content)
                # Only add two newlines between files (not after the last file)
                if i < len(gbk_files) - 1:
                    outfile.write("\n\n")
    
    print(f"Successfully created combined file: {combined_path}")
    
    # Create a build script
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
    return True

if __name__ == "__main__":
    combine_genbank_files()