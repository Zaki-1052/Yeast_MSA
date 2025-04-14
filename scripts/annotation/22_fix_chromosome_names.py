#!/usr/bin/env python3

# File: scripts/annotation/22_fix_chromosome_names.py
# Purpose: Create direct mapping between JRIU chromosomes and w303 scaffolds

import os
import subprocess
import pandas as pd
import re
from datetime import datetime

def main():
    print("=== Creating Direct JRIU to W303 Chromosome Mapping ===")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("")
    
    # Define directories
    vcf_source = "annotation/vcf_ready"
    snpeff_dir = "/Users/zakiralibhai/snpEff"
    mapping_dir = "annotation/chromosome_mapping"
    fixed_vcf_dir = "annotation/vcf_fixed_direct"
    debug_dir = "annotation/debug_chroms"
    
    # Create output directories
    os.makedirs(mapping_dir, exist_ok=True)
    os.makedirs(fixed_vcf_dir, exist_ok=True)
    os.makedirs(debug_dir, exist_ok=True)
    
    # Step 1: Extract JRIU chromosome IDs from VCF
    print("Step 1: Extracting JRIU chromosome IDs from VCF...")
    
    vcf_files = os.listdir(vcf_source)
    vcf_files = [f for f in vcf_files if f.endswith('.sorted.vcf.gz')]
    
    if not vcf_files:
        print("ERROR: No sorted VCF files found")
        return
    
    sample_vcf = os.path.join(vcf_source, vcf_files[0])
    jriu_chroms = []
    
    try:
        cmd = f"bcftools view -h {sample_vcf} | grep '##contig=' | cut -d'=' -f3 | cut -d',' -f1"
        output = subprocess.check_output(cmd, shell=True, text=True)
        jriu_chroms = [line.strip() for line in output.splitlines() if line.strip()]
        
        # Extract numbers from JRIU IDs
        jriu_nums = []
        for chrom in jriu_chroms:
            # Remove the ID= prefix
            chrom = chrom.replace('ID=', '')
            # Match JRIU\d+(\d+)\.1 pattern
            match = re.search(r'JRIU\d+0*(\d+)\.1', chrom)
            if match:
                jriu_nums.append(int(match.group(1)))
            else:
                jriu_nums.append(None)
        
        print(f"Found {len(jriu_chroms)} JRIU chromosomes")
        print(f"First few: {', '.join(jriu_chroms[:5])}")
        
        # Save to file
        with open(f"{debug_dir}/jriu_chroms.txt", 'w') as f:
            for i, chrom in enumerate(jriu_chroms):
                num = jriu_nums[i] if jriu_nums[i] is not None else "NA"
                f.write(f"{chrom}\t{num}\n")
        
        print(f"Saved JRIU chromosomes to {debug_dir}/jriu_chroms.txt")
    except subprocess.CalledProcessError as e:
        print(f"Error extracting JRIU chromosomes: {e}")
        return
    
    # Step 2: Extract w303 scaffold IDs from SnpEff database
    print("\nStep 2: Extracting w303 scaffold IDs from SnpEff database...")
    
    w303_scaffolds = []
    
    try:
        os.chdir(snpeff_dir)
        cmd = "java -jar snpEff.jar dump w303 | grep -A500 'Chromosomes' | grep 'w303_scaffold_' | cut -d\"'\" -f2"
        output = subprocess.check_output(cmd, shell=True, text=True)
        os.chdir(os.path.dirname(os.path.abspath(__file__)) + "/../..")  # Return to original directory
        
        w303_scaffolds = [line.strip() for line in output.splitlines() if line.strip()]
        
        # Extract numbers from w303 scaffolds
        w303_nums = []
        for scaffold in w303_scaffolds:
            match = re.search(r'w303_scaffold_(\d+)', scaffold)
            if match:
                w303_nums.append(int(match.group(1)))
            else:
                w303_nums.append(None)
        
        print(f"Found {len(w303_scaffolds)} w303 scaffolds")
        print(f"First few: {', '.join(w303_scaffolds[:5])}")
        
        # Save to file
        with open(f"{debug_dir}/w303_scaffolds.txt", 'w') as f:
            for i, scaffold in enumerate(w303_scaffolds):
                num = w303_nums[i] if w303_nums[i] is not None else "NA"
                f.write(f"{scaffold}\t{num}\n")
        
        print(f"Saved w303 scaffolds to {debug_dir}/w303_scaffolds.txt")
    except subprocess.CalledProcessError as e:
        print(f"Error extracting w303 scaffolds: {e}")
        return
    
    # Step 3: Create direct mapping based on numeric ordering
    print("\nStep 3: Creating direct mapping between JRIU and w303 scaffolds...")
    
    # Sort JRIU chromosomes by number
    jriu_sorted = sorted(zip(jriu_chroms, jriu_nums), key=lambda x: x[1] if x[1] is not None else float('inf'))
    
    # Sort w303 scaffolds by number
    w303_sorted = sorted(zip(w303_scaffolds, w303_nums), key=lambda x: x[1] if x[1] is not None else float('inf'))
    
    # Create direct mapping
    direct_mapping = []
    
    # Only map up to the smaller of the two lists
    map_count = min(len(jriu_sorted), len(w303_sorted))
    
    for i in range(map_count):
        jriu_chrom, jriu_num = jriu_sorted[i]
        w303_scaffold, w303_num = w303_sorted[i]
        
        direct_mapping.append({
            'JRIU_Chromosome': jriu_chrom,
            'JRIU_Number': jriu_num,
            'W303_Scaffold': w303_scaffold,
            'W303_Number': w303_num,
            'Mapping_Type': 'Numeric_Order'
        })
    
    # Create mapping DataFrame
    mapping_df = pd.DataFrame(direct_mapping)
    mapping_file = f"{mapping_dir}/jriu_to_w303_direct_mapping.tsv"
    mapping_df.to_csv(mapping_file, sep='\t', index=False)
    
    print(f"Created mapping between {map_count} chromosomes")
    print(f"Saved mapping to: {mapping_file}")
    print(f"First few mappings:")
    for i, row in mapping_df.head(5).iterrows():
        print(f"  {row['JRIU_Chromosome']} → {row['W303_Scaffold']}")
    
    # Step 4: Modify VCF files using the direct mapping
    print("\nStep 4: Modifying VCF files with direct chromosome mapping...")
    
    # Create a lookup dictionary for faster access
    chrom_lookup = dict(zip(mapping_df['JRIU_Chromosome'], mapping_df['W303_Scaffold']))
    
    for vcf_file in vcf_files:
        sample = vcf_file.replace('.sorted.vcf.gz', '')
        input_path = os.path.join(vcf_source, vcf_file)
        output_path = os.path.join(fixed_vcf_dir, f"{sample}.fixed.vcf")
        
        print(f"- Processing {sample}...")
        
        # Extract and modify header
        header_file = f"{fixed_vcf_dir}/{sample}.header.txt"
        try:
            cmd = f"bcftools view -h {input_path} > {header_file}"
            subprocess.run(cmd, shell=True, check=True)
            
            # Modify the header
            with open(header_file, 'r') as f:
                header_lines = f.readlines()
            
            with open(f"{fixed_vcf_dir}/{sample}.new_header.txt", 'w') as f:
                for line in header_lines:
                    if line.startswith('##contig='):
                        # Extract the JRIU ID
                        match = re.search(r'ID=([^,]+)', line)
                        if match:
                            jriu_id = match.group(1)
                            if jriu_id in chrom_lookup:
                                # Replace with w303 scaffold
                                new_line = line.replace(f'ID={jriu_id}', f'ID={chrom_lookup[jriu_id]}')
                                f.write(new_line)
                            else:
                                # Keep original if not mapped
                                f.write(line)
                    else:
                        f.write(line)
            
            # Extract and modify data
            data_file = f"{fixed_vcf_dir}/{sample}.data.txt"
            cmd = f"bcftools view -H {input_path} > {data_file}"
            subprocess.run(cmd, shell=True, check=True)
            
            # Process the data file
            with open(data_file, 'r') as f:
                data_lines = f.readlines()
            
            with open(f"{fixed_vcf_dir}/{sample}.new_data.txt", 'w') as f:
                for line in data_lines:
                    fields = line.strip().split('\t')
                    if len(fields) > 0:
                        jriu_id = fields[0]
                        if jriu_id in chrom_lookup:
                            # Replace chromosome name
                            fields[0] = chrom_lookup[jriu_id]
                            f.write('\t'.join(fields) + '\n')
                        else:
                            # Keep original if not mapped
                            f.write(line)
                    else:
                        f.write(line)
            
            # Combine header and data
            cmd = f"cat {fixed_vcf_dir}/{sample}.new_header.txt {fixed_vcf_dir}/{sample}.new_data.txt > {output_path}"
            subprocess.run(cmd, shell=True, check=True)
            
            # Compress and index
            cmd = f"bgzip -f {output_path}"
            subprocess.run(cmd, shell=True, check=True)
            
            cmd = f"tabix -p vcf {output_path}.gz"
            subprocess.run(cmd, shell=True, check=True)
            
            # Clean up temporary files
            os.remove(header_file)
            os.remove(f"{fixed_vcf_dir}/{sample}.new_header.txt")
            os.remove(data_file)
            os.remove(f"{fixed_vcf_dir}/{sample}.new_data.txt")
            
            print(f"  ✓ Successfully processed {sample}")
        except subprocess.CalledProcessError as e:
            print(f"  ✗ Error processing {sample}: {e}")
    
    # Step 5: Run SnpEff on a sample fixed VCF to verify
    print("\nStep 5: Testing SnpEff annotation on fixed VCF...")
    
    fixed_vcfs = os.listdir(fixed_vcf_dir)
    fixed_vcfs = [f for f in fixed_vcfs if f.endswith('.fixed.vcf.gz')]
    
    if fixed_vcfs:
        sample_vcf = os.path.join(fixed_vcf_dir, fixed_vcfs[0])
        sample = fixed_vcfs[0].replace('.fixed.vcf.gz', '')
        test_output = os.path.join(debug_dir, f"{sample}.test_annotated.vcf")
        
        print(f"Testing annotation on {sample}...")
        
        cmd = f"java -jar {snpeff_dir}/snpEff.jar -v w303 {sample_vcf} > {test_output}"
        try:
            subprocess.run(cmd, shell=True, check=True)
            
            # Check if annotation worked
            cmd = f"grep -c 'ERROR_CHROMOSOME_NOT_FOUND' {test_output}"
            try:
                error_count = int(subprocess.check_output(cmd, shell=True, text=True).strip())
                if error_count > 0:
                    print(f"  ✗ Still found {error_count} chromosome errors")
                    print("  Annotation still not working correctly")
                else:
                    print("  ✓ No chromosome errors found!")
                    print("  Annotation appears to be working correctly")
                    
                    # Check for gene names
                    cmd = f"grep -c 'ANN=' {test_output}"
                    ann_count = int(subprocess.check_output(cmd, shell=True, text=True).strip())
                    
                    # Extract a sample annotation
                    cmd = f"grep 'ANN=' {test_output} | head -1"
                    sample_ann = subprocess.check_output(cmd, shell=True, text=True).strip()
                    
                    print(f"  Found {ann_count} annotations")
                    print(f"  Sample annotation: {sample_ann[:200]}...")
                    
                    # Create a script to run SnpEff on all fixed VCF files
                    run_script = os.path.join("scripts/annotation", "23_run_snpeff_fixed.sh")
                    with open(run_script, 'w') as f:
                        f.write("""#!/bin/bash

# File: scripts/annotation/23_run_snpeff_fixed.sh
# Purpose: Run SnpEff on all fixed VCF files

echo "=== Running SnpEff on Fixed VCF Files ==="
echo "Date: $(date)"
echo ""

# Define directories
SNPEFF_DIR="/Users/zakiralibhai/snpEff"
VCF_DIR="annotation/vcf_fixed_direct"
RESULTS_DIR="annotation/results_fixed"
STATS_DIR="annotation/stats_fixed"

# Create output directories
mkdir -p "$RESULTS_DIR"
mkdir -p "$STATS_DIR"

# Find all fixed VCF files
VCF_FILES=$(find "$VCF_DIR" -name "*.fixed.vcf.gz")
VCF_COUNT=$(echo "$VCF_FILES" | wc -l)

echo "Found $VCF_COUNT fixed VCF files to annotate"
echo ""

# Process each file
for VCF_FILE in $VCF_FILES; do
    SAMPLE=$(basename "$VCF_FILE" .fixed.vcf.gz)
    
    echo "Processing $SAMPLE..."
    
    # Run SnpEff
    java -Xmx4g -jar "$SNPEFF_DIR/snpEff.jar" \\
        -v \\
        -stats "$STATS_DIR/${SAMPLE}.stats.html" \\
        w303 \\
        "$VCF_FILE" \\
        > "$RESULTS_DIR/${SAMPLE}.snpeff.vcf"
    
    # Check for errors
    ERROR_COUNT=$(grep -c "ERROR_CHROMOSOME_NOT_FOUND" "$RESULTS_DIR/${SAMPLE}.snpeff.vcf")
    if [ "$ERROR_COUNT" -gt 0 ]; then
        echo "  ✗ Found $ERROR_COUNT chromosome errors"
    else
        echo "  ✓ No chromosome errors found"
    fi
    
    # Compress and index
    bgzip -f "$RESULTS_DIR/${SAMPLE}.snpeff.vcf"
    tabix -p vcf "$RESULTS_DIR/${SAMPLE}.snpeff.vcf.gz"
    
    echo "Done processing $SAMPLE"
    echo ""
done

echo "All VCF files have been annotated"
echo ""
echo "=== SnpEff Annotation Complete ==="
""")
                    
                    os.chmod(run_script, 0o755)
                    print(f"\nCreated script to run SnpEff on all fixed VCF files: {run_script}")
                    print("You can run this script to annotate all fixed VCF files")
                    
                    # Also create a script to extract target genes
                    extract_script = os.path.join("scripts/annotation", "24_extract_target_genes_fixed.py")
                    with open(extract_script, 'w') as f:
                        f.write("""#!/usr/bin/env python3

# File: scripts/annotation/24_extract_target_genes_fixed.py
# Purpose: Extract variants in target genes from fixed annotated VCF files

import os
import gzip
import re
import pandas as pd
from datetime import datetime

def main():
    print("=== Extracting Target Genes from Fixed Annotated VCFs ===")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("")
    
    # Define directories
    fixed_dir = "annotation/results_fixed"
    gene_dir = "annotation/genes_of_interest"
    output_dir = "annotation/gene_variants_fixed"
    
    # Create output directories
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(f"{output_dir}/by_gene", exist_ok=True)
    
    # Read target genes
    target_genes_file = f"{gene_dir}/target_genes.txt"
    with open(target_genes_file, 'r') as f:
        target_genes = [line.strip() for line in f if line.strip()]
    
    print(f"Looking for {len(target_genes)} target genes")
    
    # Find VCF files
    vcf_files = []
    for root, dirs, files in os.walk(fixed_dir):
        for file in files:
            if file.endswith('.snpeff.vcf.gz'):
                vcf_files.append(os.path.join(root, file))
    
    if not vcf_files:
        print(f"ERROR: No annotated VCF files found in {fixed_dir}")
        return
    
    print(f"Found {len(vcf_files)} annotated VCF files")
    print("")
    
    # Create a dictionary to store variants by gene
    gene_variants = {gene: [] for gene in target_genes}
    all_variants = []
    
    # Process each VCF file
    for vcf_file in vcf_files:
        sample = os.path.basename(vcf_file).replace('.snpeff.vcf.gz', '')
        print(f"Processing {sample}...")
        
        # Count variants found per gene
        gene_counts = {gene: 0 for gene in target_genes}
        
        # Parse the VCF file
        try:
            with gzip.open(vcf_file, 'rt') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    # Check if this variant affects any of our target genes
                    for gene in target_genes:
                        if gene in line:
                            # Parse the line
                            fields = line.strip().split('\\t')
                            if len(fields) < 8:
                                continue
                            
                            chrom = fields[0]
                            pos = fields[1]
                            ref = fields[3]
                            alt = fields[4]
                            info = fields[7]
                            
                            # Extract the annotation
                            ann_match = re.search(r'ANN=([^;]+)', info)
                            if not ann_match:
                                continue
                            
                            ann_field = ann_match.group(1)
                            
                            # Find annotations for this gene
                            gene_anns = [a for a in ann_field.split(',') if gene in a]
                            
                            for ann in gene_anns:
                                parts = ann.split('|')
                                if len(parts) < 10:
                                    continue
                                
                                # Get annotation details
                                effect = parts[1] if len(parts) > 1 else ""
                                impact = parts[2] if len(parts) > 2 else ""
                                gene_name = parts[3] if len(parts) > 3 else ""
                                protein_change = parts[10] if len(parts) > 10 else ""
                                
                                # Create a record
                                record = {
                                    "Sample": sample,
                                    "Gene": gene,
                                    "Chromosome": chrom,
                                    "Position": pos,
                                    "Ref": ref,
                                    "Alt": alt,
                                    "Effect": effect,
                                    "Impact": impact,
                                    "Gene_Name": gene_name,
                                    "Protein_Change": protein_change
                                }
                                
                                # Add to results
                                gene_variants[gene].append(record)
                                all_variants.append(record)
                                
                                # Increment count
                                gene_counts[gene] += 1
                                
                                # Only process the first matching annotation for this gene in this variant
                                break
        except Exception as e:
            print(f"  Error processing {vcf_file}: {e}")
            continue
        
        # Report counts
        for gene, count in gene_counts.items():
            if count > 0:
                print(f"  Found {count} variants affecting {gene}")
    
    # Save results
    print("\nSaving results...")
    
    # Save all variants
    all_df = pd.DataFrame(all_variants)
    if not all_df.empty:
        all_df.to_csv(f"{output_dir}/all_target_variants.tsv", sep='\\t', index=False)
        print(f"Saved {len(all_df)} variants to {output_dir}/all_target_variants.tsv")
    else:
        print("No variants found for target genes")
    
    # Save variants by gene
    for gene, variants in gene_variants.items():
        gene_df = pd.DataFrame(variants)
        if not gene_df.empty:
            gene_df.to_csv(f"{output_dir}/by_gene/{gene}_variants.tsv", sep='\\t', index=False)
            print(f"Saved {len(gene_df)} variants for {gene}")
    
    # Create a summary by gene
    gene_summary = []
    for gene in target_genes:
        df = pd.DataFrame(gene_variants[gene])
        if df.empty:
            gene_summary.append({
                "Gene": gene,
                "Total_Variants": 0,
                "High_Impact": 0,
                "Moderate_Impact": 0,
                "Low_Impact": 0,
                "Modifier_Impact": 0,
                "Samples_With_Variants": 0
            })
        else:
            unique_samples = len(df['Sample'].unique())
            gene_summary.append({
                "Gene": gene,
                "Total_Variants": len(df),
                "High_Impact": len(df[df["Impact"] == "HIGH"]),
                "Moderate_Impact": len(df[df["Impact"] == "MODERATE"]),
                "Low_Impact": len(df[df["Impact"] == "LOW"]),
                "Modifier_Impact": len(df[df["Impact"] == "MODIFIER"]),
                "Samples_With_Variants": unique_samples
            })
    
    # Save gene summary
    summary_df = pd.DataFrame(gene_summary)
    summary_df.to_csv(f"{output_dir}/gene_summary.tsv", sep='\\t', index=False)
    print(f"Saved gene summary to {output_dir}/gene_summary.tsv")
    
    # Create a summary by sample
    sample_summary = []
    for sample in set([v["Sample"] for v in all_variants]):
        sample_variants = [v for v in all_variants if v["Sample"] == sample]
        sample_df = pd.DataFrame(sample_variants)
        
        unique_genes = len(sample_df['Gene'].unique())
        sample_summary.append({
            "Sample": sample,
            "Total_Variants": len(sample_df),
            "High_Impact": len(sample_df[sample_df["Impact"] == "HIGH"]),
            "Moderate_Impact": len(sample_df[sample_df["Impact"] == "MODERATE"]),
            "Low_Impact": len(sample_df[sample_df["Impact"] == "LOW"]),
            "Modifier_Impact": len(sample_df[sample_df["Impact"] == "MODIFIER"]),
            "Genes_With_Variants": unique_genes
        })
    
    if sample_summary:
        # Save sample summary
        sample_df = pd.DataFrame(sample_summary)
        sample_df.to_csv(f"{output_dir}/sample_summary.tsv", sep='\\t', index=False)
        print(f"Saved sample summary to {output_dir}/sample_summary.tsv")
    
    print("\n=== Target Gene Extraction Complete ===")

if __name__ == "__main__":
    main()
""")
                    
                    os.chmod(extract_script, 0o755)
                    print(f"Created script to extract target genes: {extract_script}")
                    print("Run the annotation script first, then the extraction script")
            except subprocess.CalledProcessError:
                print("  No errors found, annotation may be successful")
        except subprocess.CalledProcessError as e:
            print(f"  ✗ Error running SnpEff: {e}")
    else:
        print("No fixed VCF files to test")
    
    print("\n=== Direct Chromosome Mapping Complete ===")

if __name__ == "__main__":
    main()
