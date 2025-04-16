#!/usr/bin/env python3

# File: scripts/annotation/21_examine_annotations.py
# Purpose: Examine the exact format of annotation fields

import os
import subprocess
import re
from datetime import datetime

def main():
    print("=== Examining Annotation Formats in Detail ===")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("")
    
    # Define directories
    results_dir = "annotation/results"
    debug_dir = "annotation/debug_annotations"
    os.makedirs(debug_dir, exist_ok=True)
    
    # Step 1: Examine the VCF header format for ANN field
    print("Step 1: Examining VCF header for annotation format...")
    
    vcf_files = [f for f in os.listdir(results_dir) if f.endswith('.vcf.gz')]
    if vcf_files:
        sample_file = os.path.join(results_dir, vcf_files[0])
        print(f"Using sample file: {sample_file}")
        
        # Try alternative methods to extract header
        try:
            # Try bcftools
            header_cmd = f"bcftools view -h {sample_file} > {debug_dir}/header.txt"
            subprocess.run(header_cmd, shell=True, check=True)
            print(f"Saved header to {debug_dir}/header.txt")
            
            # Look for ANN format in header
            ann_format_cmd = f"grep -A 5 'ID=ANN' {debug_dir}/header.txt > {debug_dir}/ann_format.txt"
            subprocess.run(ann_format_cmd, shell=True)
            
            # Display ANN format
            if os.path.exists(f"{debug_dir}/ann_format.txt") and os.path.getsize(f"{debug_dir}/ann_format.txt") > 0:
                with open(f"{debug_dir}/ann_format.txt", 'r') as f:
                    ann_format = f.read()
                print("\nAnnotation format from header:")
                print(ann_format)
            else:
                print("\nCould not find ANN format in header")
        except subprocess.CalledProcessError:
            print("Error extracting header with bcftools")
    else:
        print("No VCF files found")
    
    # Step 2: Extract sample annotations
    print("\nStep 2: Extracting sample annotations...")
    
    if vcf_files:
        sample_file = os.path.join(results_dir, vcf_files[0])
        
        try:
            # Extract first few variants with annotations (not just header)
            variants_cmd = f"bcftools view {sample_file} | grep -v '^#' | head -5 > {debug_dir}/variants.txt"
            subprocess.run(variants_cmd, shell=True, check=True)
            
            # Display variants
            with open(f"{debug_dir}/variants.txt", 'r') as f:
                variants = f.read()
            print("\nSample variants:")
            print(variants)
            
            # Try to extract ANN fields from these variants
            print("\nParsing ANN fields...")
            with open(f"{debug_dir}/variants.txt", 'r') as f:
                for i, line in enumerate(f):
                    print(f"\nVariant {i+1}:")
                    fields = line.strip().split('\t')
                    if len(fields) >= 8:
                        info = fields[7]
                        print(f"INFO field: {info}")
                        
                        # Look for ANN field
                        ann_match = re.search(r'ANN=([^;]+)', info)
                        if ann_match:
                            ann = ann_match.group(1)
                            print(f"ANN field: {ann}")
                            
                            # Split by pipe to see format
                            parts = ann.split('|')
                            print(f"Number of fields: {len(parts)}")
                            print("Fields:")
                            for j, part in enumerate(parts[:min(len(parts), 15)]):
                                print(f"  {j}: '{part}'")
                        else:
                            print("No ANN field found")
                    else:
                        print("Invalid variant line format")
        except subprocess.CalledProcessError:
            print("Error extracting variants")
    
    # Step 3: Try direct grep search in VCF for target genes
    print("\nStep 3: Searching directly for target genes in VCF...")
    
    target_genes = ["YHR072W", "YHR007C", "YNL280C", "YHR190W", "YGR175C", "YML008C"]
    mapped_genes = ["w303_scaffold_80", "w303_scaffold_139", "w303_scaffold_15"]
    
    if vcf_files:
        sample_file = os.path.join(results_dir, vcf_files[0])
        
        # Check for target genes by SGD name
        print("\nSearching for SGD gene names:")
        for gene in target_genes:
            try:
                cmd = f"bcftools view {sample_file} | grep -v '^#' | grep -w '{gene}' | wc -l"
                count = subprocess.check_output(cmd, shell=True, text=True).strip()
                if int(count) > 0:
                    print(f"  ✓ Found {count} occurrences of {gene}")
                else:
                    print(f"  ✗ No occurrences of {gene}")
            except subprocess.CalledProcessError:
                print(f"  Error searching for {gene}")
        
        # Check for mapped genes by w303 scaffold ID
        print("\nSearching for W303 scaffold IDs:")
        for gene in mapped_genes:
            try:
                cmd = f"bcftools view {sample_file} | grep -v '^#' | grep -w '{gene}' | wc -l"
                count = subprocess.check_output(cmd, shell=True, text=True).strip()
                if int(count) > 0:
                    print(f"  ✓ Found {count} occurrences of {gene}")
                else:
                    print(f"  ✗ No occurrences of {gene}")
            except subprocess.CalledProcessError:
                print(f"  Error searching for {gene}")
    
    # Step 4: Check gene IDs in a different format
    print("\nStep 4: Checking for alternative gene ID formats...")
    
    if vcf_files:
        sample_file = os.path.join(results_dir, vcf_files[0])
        
        # Extract all unique identifiers from field 4 of ANN (presumed gene ID)
        try:
            # This command extracts the pipe-delimited field that should contain gene IDs
            cmd = f"bcftools view {sample_file} | grep -v '^#' | grep 'ANN=' | head -100 | sed 's/.*ANN=//' | sed 's/;.*//' | tr ',' '\\n' | cut -d'|' -f4 | sort | uniq > {debug_dir}/gene_ids.txt"
            subprocess.run(cmd, shell=True, check=True)
            
            # Display the gene IDs
            with open(f"{debug_dir}/gene_ids.txt", 'r') as f:
                gene_ids = f.read()
            
            print("\nUnique values in field 4 of ANN (potential gene IDs):")
            print(gene_ids)
            
            # Count the number of unique gene IDs
            count = len([line for line in gene_ids.splitlines() if line.strip()])
            print(f"Found {count} unique potential gene identifiers")
            
            # Check for common prefixes
            prefixes = {}
            for line in gene_ids.splitlines():
                if line.strip():
                    prefix = ''.join([c for c in line if not c.isdigit()]).strip()
                    prefixes[prefix] = prefixes.get(prefix, 0) + 1
            
            print("\nCommon gene ID prefixes:")
            for prefix, count in sorted(prefixes.items(), key=lambda x: x[1], reverse=True):
                if prefix:  # Skip empty prefix
                    print(f"  {prefix}: {count} occurrences")
        except subprocess.CalledProcessError:
            print("Error extracting gene IDs")
    
    # Step 5: Check one random full annotation
    print("\nStep 5: Examining one complete annotation entry...")
    
    if vcf_files:
        sample_file = os.path.join(results_dir, vcf_files[0])
        
        try:
            # Extract one random annotation with context
            cmd = f"bcftools view {sample_file} | grep -v '^#' | grep 'ANN=' | head -1 > {debug_dir}/full_annotation.txt"
            subprocess.run(cmd, shell=True, check=True)
            
            with open(f"{debug_dir}/full_annotation.txt", 'r') as f:
                full_anno = f.read()
            
            print("\nComplete annotation entry:")
            print(full_anno)
            
            # Analyze the annotation format
            if 'ANN=' in full_anno:
                # Extract and parse the ANN field
                ann_match = re.search(r'ANN=([^;]+)', full_anno)
                if ann_match:
                    ann_field = ann_match.group(1)
                    annotations = ann_field.split(',')
                    
                    print(f"\nNumber of comma-separated annotations: {len(annotations)}")
                    print("First annotation parsed:")
                    
                    parts = annotations[0].split('|')
                    field_names = ["Allele", "Annotation", "Impact", "Gene_Name", "Gene_ID", 
                                  "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", 
                                  "HGVS.c", "HGVS.p", "cDNA.pos", "CDS.pos", "AA.pos", 
                                  "Distance", "Errors"]
                    
                    for i, part in enumerate(parts):
                        field_name = field_names[i] if i < len(field_names) else f"Field_{i}"
                        print(f"  {field_name}: '{part}'")
        except subprocess.CalledProcessError:
            print("Error extracting full annotation")
    
    print("\n=== Annotation Examination Complete ===")

if __name__ == "__main__":
    main()
