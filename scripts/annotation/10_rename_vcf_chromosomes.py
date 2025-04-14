#!/usr/bin/env python3

# File: scripts/annotation/10_rename_vcf_chromosomes.py
# Purpose: Rename chromosomes in VCF files to match SnpEff database

import os
import re
import subprocess
import glob
from datetime import datetime
import pandas as pd

def main():
    print("=== Renaming Chromosomes in VCF Files ===")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("")
    
    # Define directories
    vcf_source = "annotation/vcf_ready"
    gbk_dir = "annotation/reference/w303_scaffolds"
    renamed_vcf_dir = "annotation/vcf_renamed"
    mapping_dir = "annotation/chromosome_mapping"
    snpeff_dir = "/Users/zakiralibhai/snpEff"
    
    # Create output directories
    os.makedirs(renamed_vcf_dir, exist_ok=True)
    os.makedirs(mapping_dir, exist_ok=True)
    
    # Step 1: Get SnpEff chromosome names
    print("Extracting SnpEff chromosome names...")
    
    # Use SnpEff dump command to get chromosome info
    try:
        os.chdir(snpeff_dir)
        output = subprocess.check_output(
            "java -jar snpEff.jar dump w303 | grep -A 1000 '# Chromosomes' | grep -v '#'",
            shell=True, text=True
        )
        os.chdir(os.path.dirname(os.path.abspath(__file__)) + "/../..")  # Return to original directory
        
        # Parse chromosome names
        snpeff_chroms = []
        for line in output.splitlines():
            line = line.strip()
            if line:
                # Extract chromosome name (format: 'w303_scaffold_1' 336272 Standard)
                match = re.search(r"'([^']+)'", line)
                if match:
                    snpeff_chroms.append(match.group(1))
        
        print(f"Found {len(snpeff_chroms)} chromosomes in SnpEff database")
        print(f"First 5: {', '.join(snpeff_chroms[:5])}")
    except subprocess.CalledProcessError:
        print("Error getting SnpEff chromosome names")
        return
    
    # Step 2: Get GenBank scaffold numbers
    print("Extracting GenBank scaffold numbers...")
    gbk_files = glob.glob(f"{gbk_dir}/*.genbank")
    
    if not gbk_files:
        print(f"ERROR: No GenBank files found in {gbk_dir}")
        return
    
    # Extract scaffold numbers from filenames
    scaffold_info = []
    for gbk_file in gbk_files:
        filename = os.path.basename(gbk_file)
        
        # Extract scaffold number
        scaffold_match = re.search(r'scaffold_(\d+)\.genbank$', filename)
        if scaffold_match:
            scaffold_num = int(scaffold_match.group(1))
            
            # Read the first few lines to extract sequence info
            sequence_id = None
            sequence_length = None
            
            with open(gbk_file, 'r', errors='ignore') as f:
                for _ in range(50):  # Read first 50 lines
                    line = f.readline()
                    if not line:
                        break
                    
                    # Look for LOCUS line which has sequence ID and length
                    if line.startswith('LOCUS'):
                        parts = line.strip().split()
                        if len(parts) > 1:
                            sequence_id = parts[1]
                        if len(parts) > 2:
                            for part in parts:
                                if part.isdigit():
                                    sequence_length = int(part)
                                    break
                    
                    # Check for VERSION line with accession
                    if line.startswith('VERSION'):
                        version_match = re.search(r'([A-Z0-9]+\.\d+)', line)
                        if version_match:
                            sequence_id = version_match.group(1)
            
            scaffold_info.append({
                'GenBank_File': filename,
                'Scaffold_Number': scaffold_num,
                'Scaffold_Name': f"w303_scaffold_{scaffold_num}",
                'Sequence_ID': sequence_id,
                'Sequence_Length': sequence_length
            })
    
    # Create scaffold DataFrame
    scaffold_df = pd.DataFrame(scaffold_info)
    scaffold_df = scaffold_df.sort_values('Scaffold_Number')
    
    print(f"Found {len(scaffold_df)} scaffolds in GenBank files")
    
    # Step 3: Get VCF chromosome names
    print("Extracting VCF chromosome names...")
    vcf_files = glob.glob(f"{vcf_source}/*.sorted.vcf.gz")
    
    if not vcf_files:
        print(f"ERROR: No VCF files found in {vcf_source}")
        return
    
    # Get chromosomes from first VCF
    vcf_chroms = []
    jriu_nums = []
    
    try:
        output = subprocess.check_output(
            f"bcftools view -h {vcf_files[0]} | grep '##contig='",
            shell=True, text=True
        )
        
        # Extract contig names and parse JRIU numbers
        for line in output.splitlines():
            match = re.search(r'ID=([^,]+)', line)
            if match:
                chrom = match.group(1)
                vcf_chroms.append(chrom)
                
                # Extract JRIU number
                jriu_match = re.search(r'JRIU\d+0*(\d+)\.1$', chrom)
                if jriu_match:
                    jriu_num = int(jriu_match.group(1))
                    jriu_nums.append((chrom, jriu_num))
        
        print(f"Found {len(vcf_chroms)} chromosomes in VCF files")
        print(f"Found {len(jriu_nums)} JRIU numbers")
    except subprocess.CalledProcessError:
        print("Error extracting VCF chromosome names")
        return
    
    # Step 4: Create a mapping between JRIU and w303_scaffold names
    print("Creating chromosome mapping...")
    
    # Sort JRIU identifiers by number
    jriu_nums.sort(key=lambda x: x[1])
    
    # Create mapping by aligning sorted JRIU numbers with sorted scaffold numbers
    # This assumes the numeric ordering corresponds between the two naming systems
    mapping = []
    
    # Track stats
    mapped_count = 0
    unmapped_count = 0
    
    # First try to map by matching JRIU numbers to scaffold numbers
    mapped_jrius = set()
    for jriu_chrom, jriu_num in jriu_nums:
        # Find scaffold with matching number (if exists)
        scaffold_match = scaffold_df[scaffold_df['Scaffold_Number'] == jriu_num]
        
        if not scaffold_match.empty:
            snpeff_chrom = scaffold_match.iloc[0]['Scaffold_Name']
            mapping.append({
                'VCF_Chromosome': jriu_chrom,
                'SnpEff_Chromosome': snpeff_chrom,
                'JRIU_Number': jriu_num,
                'Scaffold_Number': jriu_num,
                'Mapping_Type': 'Number match',
                'GenBank_File': scaffold_match.iloc[0]['GenBank_File']
            })
            mapped_jrius.add(jriu_chrom)
            mapped_count += 1
    
    # For remaining unmapped JRIU chroms, map in numeric order to remaining scaffolds
    remaining_scaffolds = scaffold_df.sort_values('Scaffold_Number')
    scaffold_idx = 0
    
    for jriu_chrom, jriu_num in jriu_nums:
        if jriu_chrom not in mapped_jrius:
            # Find next available scaffold
            while scaffold_idx < len(remaining_scaffolds):
                scaffold_row = remaining_scaffolds.iloc[scaffold_idx]
                scaffold_idx += 1
                
                # Check if this scaffold is already mapped
                if not any(m['SnpEff_Chromosome'] == scaffold_row['Scaffold_Name'] for m in mapping):
                    mapping.append({
                        'VCF_Chromosome': jriu_chrom,
                        'SnpEff_Chromosome': scaffold_row['Scaffold_Name'],
                        'JRIU_Number': jriu_num,
                        'Scaffold_Number': scaffold_row['Scaffold_Number'],
                        'Mapping_Type': 'Sequential',
                        'GenBank_File': scaffold_row['GenBank_File']
                    })
                    mapped_count += 1
                    break
            else:
                # No more scaffolds available
                unmapped_count += 1
    
    # Create mapping DataFrame
    mapping_df = pd.DataFrame(mapping)
    mapping_file = f"{mapping_dir}/chromosome_mapping_final.tsv"
    mapping_df.to_csv(mapping_file, sep='\t', index=False)
    
    print(f"Created mapping with {mapped_count} chromosomes")
    print(f"Unable to map {unmapped_count} chromosomes")
    print(f"Mapping saved to: {mapping_file}")
    
    # Create lookup dictionary for faster access
    chrom_lookup = {row['VCF_Chromosome']: row['SnpEff_Chromosome'] for _, row in mapping_df.iterrows()}
    
    # Step 5: Modify VCF files with new chromosome names
    print("\nModifying VCF files with new chromosome names...")
    
    for vcf_file in vcf_files:
        sample = os.path.basename(vcf_file).replace('.sorted.vcf.gz', '')
        output_file = f"{renamed_vcf_dir}/{sample}.renamed.vcf"
        
        print(f"- Processing {sample}...")
        
        # Create a new header with modified contig lines
        header_file = f"{renamed_vcf_dir}/{sample}.header.vcf"
        
        with open(header_file, 'w') as out:
            # Get header from original file
            header = subprocess.check_output(
                f"bcftools view -h {vcf_file}",
                shell=True, text=True
            )
            
            # Replace contig IDs in header
            header_lines = header.splitlines()
            modified_header_lines = []
            
            for line in header_lines:
                if line.startswith('##contig='):
                    match = re.search(r'ID=([^,]+)', line)
                    if match:
                        old_id = match.group(1)
                        if old_id in chrom_lookup:
                            new_id = chrom_lookup[old_id]
                            modified_line = line.replace(f'ID={old_id}', f'ID={new_id}')
                            modified_header_lines.append(modified_line)
                else:
                    modified_header_lines.append(line)
            
            # Write modified header
            out.write('\n'.join(modified_header_lines) + '\n')
        
        # Extract data lines with modified chromosome names
        data_file = f"{renamed_vcf_dir}/{sample}.data.vcf"
        
        with open(data_file, 'w') as out:
            # Use bcftools to extract data lines
            subprocess.run(
                f"bcftools view -H {vcf_file} > {data_file}",
                shell=True
            )
        
        # Replace chromosome names in data file
        with open(data_file, 'r') as f_in, open(output_file, 'w') as f_out:
            for line in f_in:
                fields = line.strip().split('\t')
                if len(fields) >= 8:
                    old_chrom = fields[0]
                    if old_chrom in chrom_lookup:
                        fields[0] = chrom_lookup[old_chrom]
                    f_out.write('\t'.join(fields) + '\n')
        
        # Combine header and data
        os.system(f"cat {header_file} >> {output_file}.tmp")
        os.system(f"cat {output_file} >> {output_file}.tmp")
        os.system(f"mv {output_file}.tmp {output_file}")
        
        # Compress and index
        print(f"  Compressing and indexing {sample}...")
        os.system(f"bgzip -f {output_file}")
        os.system(f"tabix -p vcf {output_file}.gz")
        
        # Clean up temporary files
        os.remove(header_file)
        os.remove(data_file)
    
    print("All VCF files have been renamed")
    print("")
    
    # Step 6: Test an annotation with a renamed VCF
    print("Testing annotation with renamed VCF...")
    
    test_vcf = f"{renamed_vcf_dir}/{os.path.basename(vcf_files[0]).replace('.sorted.vcf.gz', '.renamed.vcf.gz')}"
    test_output = f"{renamed_vcf_dir}/test_annotation.vcf"
    
    cmd = f"java -jar {snpeff_dir}/snpEff.jar -v w303 {test_vcf} > {test_output}"
    os.system(cmd)
    
    # Check for errors
    error_count = 0
    with open(test_output, 'r') as f:
        for line in f:
            if "ERROR_CHROMOSOME_NOT_FOUND" in line:
                error_count += 1
    
    if error_count > 0:
        print(f"❌ Test failed: Found {error_count} chromosome errors")
        print("First few errors:")
        os.system(f"grep 'ERROR_CHROMOSOME_NOT_FOUND' {test_output} | head -5")
    else:
        print("✅ Test successful! No chromosome errors found")
        print("Proceeding to annotate all VCF files")
        
        # Create scripts for full annotation
        annotation_script = "scripts/annotation/11_annotate_all_renamed_vcfs.sh"
        with open(annotation_script, 'w') as f:
            f.write("#!/bin/bash\n\n")
            f.write("# File: scripts/annotation/11_annotate_all_renamed_vcfs.sh\n")
            f.write("# Purpose: Annotate all renamed VCF files\n\n")
            f.write("echo \"=== Annotating All Renamed VCF Files ===\"\n")
            f.write("echo \"Date: $(date)\"\n")
            f.write("echo \"\"\n\n")
            f.write("# Define directories\n")
            f.write("RENAMED_VCF_DIR=\"annotation/vcf_renamed\"\n")
            f.write("RESULTS_DIR=\"annotation/results_renamed\"\n")
            f.write("STATS_DIR=\"annotation/stats_renamed\"\n")
            f.write("SNPEFF_JAR=\"/Users/zakiralibhai/snpEff/snpEff.jar\"\n\n")
            f.write("# Create output directories\n")
            f.write("mkdir -p \"$RESULTS_DIR\"\n")
            f.write("mkdir -p \"$STATS_DIR\"\n\n")
            f.write("# Find all renamed VCF files\n")
            f.write("VCF_FILES=$(find \"$RENAMED_VCF_DIR\" -name \"*.renamed.vcf.gz\")\n\n")
            f.write("# Process each file\n")
            f.write("for VCF_FILE in $VCF_FILES; do\n")
            f.write("    SAMPLE=$(basename \"$VCF_FILE\" .renamed.vcf.gz)\n")
            f.write("    echo \"Processing $SAMPLE...\"\n")
            f.write("    \n")
            f.write("    # Run SnpEff\n")
            f.write("    java -Xmx4g -jar \"$SNPEFF_JAR\" \\\n")
            f.write("        -v \\\n")
            f.write("        -stats \"$STATS_DIR/${SAMPLE}.stats.html\" \\\n")
            f.write("        w303 \\\n")
            f.write("        \"$VCF_FILE\" \\\n")
            f.write("        > \"$RESULTS_DIR/${SAMPLE}.snpeff.vcf\"\n")
            f.write("    \n")
            f.write("    # Compress and index\n")
            f.write("    bgzip -f \"$RESULTS_DIR/${SAMPLE}.snpeff.vcf\"\n")
            f.write("    tabix -p vcf \"$RESULTS_DIR/${SAMPLE}.snpeff.vcf.gz\"\n")
            f.write("    \n")
            f.write("    echo \"Done processing $SAMPLE\"\n")
            f.write("    echo \"\"\n")
            f.write("done\n\n")
            f.write("echo \"All VCF files have been annotated\"\n")
            f.write("echo \"\"\n")
            f.write("echo \"=== Annotation Complete ===\"\n")
        
        os.system(f"chmod +x {annotation_script}")
        print(f"Created annotation script: {annotation_script}")
        print("Run it to annotate all VCF files")
    
    print("")
    print("=== Chromosome Renaming Complete ===")

if __name__ == "__main__":
    main()
