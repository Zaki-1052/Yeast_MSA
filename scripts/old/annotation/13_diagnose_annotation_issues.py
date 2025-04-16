#!/usr/bin/env python3

# File: scripts/annotation/13_diagnose_annotation_issues.py
# Purpose: Diagnose why we're not finding variants in target genes

import os
import gzip
import re
import subprocess
import glob
from datetime import datetime

def main():
    print("=== Diagnosing Annotation Issues ===")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("")
    
    # Define directories
    anno_source = "annotation/results_renamed"
    diag_dir = "annotation/diagnosis"
    snpeff_dir = "/Users/zakiralibhai/snpEff"
    gene_dir = "annotation/genes_of_interest"
    
    # Create output directory
    os.makedirs(diag_dir, exist_ok=True)
    
    # Make sure target genes file exists
    target_genes_file = f"{gene_dir}/target_genes.txt"
    if not os.path.exists(target_genes_file):
        with open(target_genes_file, 'w') as f:
            f.write("YHR190W\n")
            f.write("YGR175C\n")
            f.write("YHR072W\n")
            f.write("YHR007C\n")
            f.write("YNL280C\n")
            f.write("YGR060W\n")
            f.write("YML008C\n")
            f.write("YMR202W\n")
            f.write("YLR056W\n")
            f.write("YMR015C\n")
            f.write("YGL012W\n")
    
    # Step 1: Check if the annotated VCF files have ANY annotations
    print("Step 1: Checking for ANY annotations in VCF files...")
    
    vcf_files = []
    for root, dirs, files in os.walk(anno_source):
        for file in files:
            if file.endswith(".snpeff.vcf.gz"):
                vcf_files.append(os.path.join(root, file))
    
    if not vcf_files:
        print(f"ERROR: No annotated VCF files found in {anno_source}")
        return
    
    print(f"Found {len(vcf_files)} annotated VCF files")
    
    # Check first file for annotations
    test_file = vcf_files[0]
    sample = os.path.basename(test_file).replace(".snpeff.vcf.gz", "")
    print(f"Examining {sample} for annotations...")
    
    # Count ANN entries
    try:
        ann_count = int(subprocess.check_output(
            f"zcat {test_file} | grep -c 'ANN='",
            shell=True, text=True
        ).strip())
        
        print(f"Found {ann_count} entries with ANN= field")
        
        if ann_count > 0:
            # Extract a few examples
            print("Sample annotations:")
            sample_anns = subprocess.check_output(
                f"zcat {test_file} | grep 'ANN=' | head -5",
                shell=True, text=True
            ).strip()
            
            for i, line in enumerate(sample_anns.splitlines()):
                parts = line.split('\t')
                if len(parts) >= 8:
                    info = parts[7]
                    ann_match = re.search(r'ANN=([^;]+)', info)
                    if ann_match:
                        print(f"Example {i+1}: {ann_match.group(1)[:200]}...")
            
            # Save all annotations to file for examination
            with open(f"{diag_dir}/sample_annotations.txt", 'w') as f:
                f.write(f"Annotations from {sample}\n")
                f.write("==============================================\n\n")
                f.write(sample_anns)
        else:
            print("⚠️ No annotations found in the VCF file!")
    except subprocess.CalledProcessError:
        print("Error extracting annotations")
    
    # Step 2: Check if our target genes exist in the SnpEff database
    print("\nStep 2: Checking if target genes exist in SnpEff database...")
    
    # Read target genes
    with open(target_genes_file, 'r') as f:
        target_genes = [line.strip() for line in f if line.strip()]
    
    print(f"Looking for {len(target_genes)} target genes in database...")
    
    # Use snpEff dump to check for genes
    try:
        # Change to SnpEff directory
        os.chdir(snpeff_dir)
        
        # Get all genes in the database
        gene_output = subprocess.check_output(
            "java -jar snpEff.jar dump w303 | grep -A 50000 'Genes:' | grep -v '#'",
            shell=True, text=True
        )
        
        # Change back to original directory
        os.chdir(os.path.dirname(os.path.abspath(__file__)) + "/../..")
        
        # Save genes to file
        with open(f"{diag_dir}/all_genes.txt", 'w') as f:
            f.write(gene_output)
        
        # Count genes
        gene_count = len(gene_output.splitlines())
        print(f"Found {gene_count} genes in the database")
        
        # Check for target genes
        found_genes = []
        potential_matches = []
        
        for gene in target_genes:
            # Exact match
            if re.search(rf'\b{gene}\b', gene_output):
                found_genes.append(gene)
            else:
                # Look for similar names (case insensitive, with/without Y prefix)
                if re.search(rf'\b{gene[1:]}\b', gene_output, re.IGNORECASE):
                    potential_matches.append((gene, f"{gene[1:]} (without Y prefix)"))
                elif re.search(rf'\b{gene}\b', gene_output, re.IGNORECASE):
                    potential_matches.append((gene, f"{gene} (different case)"))
        
        print(f"Found {len(found_genes)} exact matches for target genes")
        if found_genes:
            print("Matched genes:")
            for gene in found_genes:
                print(f"- {gene}")
        
        if potential_matches:
            print(f"\nFound {len(potential_matches)} potential matches:")
            for gene, note in potential_matches:
                print(f"- {gene}: {note}")
        
        if not found_genes and not potential_matches:
            print("⚠️ None of the target genes were found in the database!")
        
        # Save a detailed comparison to file
        with open(f"{diag_dir}/gene_comparison.txt", 'w') as f:
            f.write("Target Genes vs Database Genes\n")
            f.write("==============================================\n\n")
            
            for gene in target_genes:
                f.write(f"Gene: {gene}\n")
                
                # Collect lines with this gene or potential matches
                matches = []
                for line in gene_output.splitlines():
                    if gene in line:
                        matches.append(line)
                    elif gene[1:] in line:  # Without Y prefix
                        matches.append(f"{line} (potential match without Y prefix)")
                
                if matches:
                    f.write("  Matches in database:\n")
                    for match in matches[:10]:  # Limit to 10 matches
                        f.write(f"  - {match}\n")
                    if len(matches) > 10:
                        f.write(f"  ... and {len(matches) - 10} more\n")
                else:
                    f.write("  No matches found in database\n")
                
                f.write("\n")
    except subprocess.CalledProcessError as e:
        print(f"Error checking genes in database: {e}")
    
    # Step 3: Check annotation statistics
    print("\nStep 3: Examining annotation statistics...")
    
    stats_dir = "annotation/stats_renamed"
    stats_files = glob.glob(f"{stats_dir}/*.stats.genes.txt")
    
    if stats_files:
        print(f"Found {len(stats_files)} stats files")
        
        # Check first stats file
        stats_file = stats_files[0]
        sample = os.path.basename(stats_file).replace(".stats.genes.txt", "")
        print(f"Examining stats for {sample}...")
        
        try:
            # Check top genes with variants
            top_genes = subprocess.check_output(
                f"cat {stats_file} | head -20",
                shell=True, text=True
            ).strip()
            
            print("Top genes with variants:")
            print(top_genes)
            
            # Check if any target genes are in the stats
            target_in_stats = []
            for gene in target_genes:
                try:
                    matches = subprocess.check_output(
                        f"grep -i {gene} {stats_file}",
                        shell=True, text=True
                    ).strip()
                    
                    if matches:
                        target_in_stats.append((gene, matches))
                except subprocess.CalledProcessError:
                    pass
            
            if target_in_stats:
                print("\nTarget genes found in stats:")
                for gene, matches in target_in_stats:
                    print(f"- {gene}: {matches}")
            else:
                print("\n⚠️ No target genes found in annotation statistics!")
        except subprocess.CalledProcessError:
            print("Error examining stats file")
    else:
        print("No stats files found")
    
    # Step 4: Check raw VCF for annotations
    print("\nStep 4: Examining annotation format in raw VCF...")
    
    try:
        # Get a few annotated variants
        variants = subprocess.check_output(
            f"zcat {test_file} | grep -v '^#' | head -10",
            shell=True, text=True
        ).strip()
        
        # Save to file
        with open(f"{diag_dir}/sample_variants.txt", 'w') as f:
            f.write(f"Sample variants from {sample}\n")
            f.write("==============================================\n\n")
            f.write(variants)
        
        # Extract and analyze ANN fields
        ann_fields = []
        for line in variants.splitlines():
            if "ANN=" in line:
                parts = line.split('\t')
                if len(parts) >= 8:
                    info = parts[7]
                    ann_match = re.search(r'ANN=([^;]+)', info)
                    if ann_match:
                        ann_fields.append(ann_match.group(1))
        
        if ann_fields:
            print(f"Found {len(ann_fields)} annotations in sample variants")
            print("Analyzing annotation format...")
            
            # Analyze first annotation
            ann = ann_fields[0]
            parts = ann.split('|')
            
            print(f"Annotation has {len(parts)} fields separated by '|'")
            
            # Check which field would have gene name
            if len(parts) >= 4:
                print(f"Potential gene name field (4th): '{parts[3]}'")
                
                # Is it empty?
                if not parts[3]:
                    print("⚠️ Gene name field is empty!")
            else:
                print("⚠️ Annotation format doesn't have enough fields for gene name!")
            
            # Save detailed analysis
            with open(f"{diag_dir}/annotation_analysis.txt", 'w') as f:
                f.write("Annotation Format Analysis\n")
                f.write("==============================================\n\n")
                
                for i, ann in enumerate(ann_fields):
                    f.write(f"Annotation {i+1}:\n")
                    parts = ann.split('|')
                    
                    for j, part in enumerate(parts):
                        f.write(f"  Field {j+1}: '{part}'\n")
                    
                    f.write("\n")
        else:
            print("No annotations found in sample variants")
    except subprocess.CalledProcessError:
        print("Error examining raw VCF")
    
    # Step 5: Generate recommendations
    print("\nStep 5: Generating recommendations...")
    
    with open(f"{diag_dir}/recommendations.txt", 'w') as f:
        f.write("Annotation Issue Recommendations\n")
        f.write("==============================================\n\n")
        
        f.write("Based on the diagnostic analysis, consider the following:\n\n")
        
        f.write("1. Gene Name Format:\n")
        f.write("   - Check if the gene names in the database match what we're searching for\n")
        f.write("   - Try searching with different formats (without 'Y' prefix, different case, etc.)\n\n")
        
        f.write("2. Database Issues:\n")
        f.write("   - Verify the w303 SnpEff database has proper gene annotations\n")
        f.write("   - Consider rebuilding the database with explicit gene annotations\n\n")
        
        f.write("3. Mapping Issues:\n")
        f.write("   - Ensure the chromosome mapping is correct for genes of interest\n")
        f.write("   - Verify which scaffolds contain our target genes\n\n")
        
        f.write("4. Annotation Format:\n")
        f.write("   - Check if the annotation format includes gene names as expected\n")
        f.write("   - Modify extraction script to handle the actual annotation format\n\n")
        
        f.write("5. Alternative Approaches:\n")
        f.write("   - Consider using a different annotation tool (e.g., VEP, ANNOVAR)\n")
        f.write("   - Create a custom annotation based on gene coordinates\n")
    
    print(f"Diagnostic results saved to: {diag_dir}/")
    print("")
    print("=== Diagnosis Complete ===")

if __name__ == "__main__":
    main()
