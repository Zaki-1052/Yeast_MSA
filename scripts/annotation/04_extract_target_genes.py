#!/usr/bin/env python3

# File: scripts/annotation/04_extract_target_genes.py
# Purpose: Extract variants affecting target genes from annotated VCFs

import os
import gzip
import re
import pandas as pd
from datetime import datetime

def parse_ann_field(ann_field, target_gene):
    """Parse the ANN field and extract information for the target gene."""
    # Remove ANN= prefix
    ann = ann_field[4:] if ann_field.startswith("ANN=") else ann_field
    
    # Split into individual annotations
    annotations = ann.split(",")
    
    # Find annotations that match our target gene
    gene_anns = [a for a in annotations if f"|{target_gene}|" in a]
    
    if not gene_anns:
        return None
    
    # Take the first annotation that matches our gene
    fields = gene_anns[0].split("|")
    
    # Check if we have enough fields
    if len(fields) < 11:
        return {
            "effect": fields[1] if len(fields) > 1 else "",
            "impact": fields[2] if len(fields) > 2 else "",
            "gene_name": fields[3] if len(fields) > 3 else "",
            "protein_change": ""
        }
    
    return {
        "effect": fields[1],
        "impact": fields[2],
        "gene_name": fields[3],
        "protein_change": fields[10]
    }

def main():
    print("=== Extracting Variants in Target Genes ===")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("")
    
    # Define directories
    anno_source = "annotation/results"
    gene_dest = "annotation/genes_of_interest"
    target_genes_file = f"{gene_dest}/target_genes.txt"
    
    # Create output directory
    os.makedirs(f"{gene_dest}/variants", exist_ok=True)
    
    # Read target genes
    print("Loading target genes...")
    with open(target_genes_file, 'r') as f:
        target_genes = [line.strip() for line in f if line.strip()]
    
    print(f"Found {len(target_genes)} target genes")
    print("")
    
    # Create a combined results dataframe
    combined_results = []
    
    # Create gene-specific dataframes
    gene_variants = {gene: [] for gene in target_genes}
    
    # Find all annotated VCF files
    print("Scanning for annotated VCF files...")
    vcf_files = []
    for root, dirs, files in os.walk(anno_source):
        for file in files:
            if file.endswith(".snpeff.vcf.gz"):
                vcf_files.append(os.path.join(root, file))
    
    if not vcf_files:
        print(f"ERROR: No annotated VCF files found in {anno_source}")
        return
    
    print(f"Found {len(vcf_files)} annotated VCF files")
    print("")
    
    # Create a log file
    log_file = f"{gene_dest}/gene_extraction.log"
    with open(log_file, 'w') as log:
        log.write("Gene Extraction Log\n")
        log.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        log.write("==============================================\n\n")
        
        # Process each annotated VCF file
        for vcf_file in vcf_files:
            filename = os.path.basename(vcf_file)
            sample = filename.replace(".snpeff.vcf.gz", "")
            
            print(f"Processing {sample}...")
            log.write(f"[{sample}] Processing started at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            
            # Process each target gene
            for gene in target_genes:
                print(f"- Extracting variants for gene {gene}...")
                log.write(f"[{sample}] Extracting variants for gene {gene}\n")
                
                variant_count = 0
                
                # Parse the VCF file
                with gzip.open(vcf_file, 'rt') as f:
                    for line in f:
                        if line.startswith('#'):
                            continue
                        
                        # Check if this variant affects our gene
                        if gene not in line:
                            continue
                        
                        # Parse the line
                        fields = line.strip().split('\t')
                        if len(fields) < 8:
                            continue
                        
                        chrom = fields[0]
                        pos = fields[1]
                        ref = fields[3]
                        alt = fields[4]
                        info = fields[7]
                        
                        # Find the ANN field
                        ann_match = re.search(r'ANN=([^;]+)', info)
                        if not ann_match:
                            continue
                        
                        ann_field = ann_match.group(0)
                        
                        # Process annotation to get effect, impact, etc.
                        ann_data = parse_ann_field(ann_field, gene)
                        if not ann_data:
                            continue
                        
                        variant_count += 1
                        
                        # Create a record
                        record = {
                            "Sample": sample,
                            "Gene": gene,
                            "Chromosome": chrom,
                            "Position": pos,
                            "Ref": ref,
                            "Alt": alt,
                            "Effect": ann_data["effect"],
                            "Impact": ann_data["impact"],
                            "Gene_Name": ann_data["gene_name"],
                            "Protein_Change": ann_data["protein_change"]
                        }
                        
                        # Add to combined results
                        combined_results.append(record)
                        
                        # Add to gene-specific results (without the Gene column)
                        gene_record = record.copy()
                        del gene_record["Gene"]
                        gene_variants[gene].append(gene_record)
                
                print(f"  Found {variant_count} variants affecting {gene}")
                log.write(f"[{sample}] Found {variant_count} variants affecting {gene}\n")
            
            print(f"Done processing {sample}")
            log.write(f"[{sample}] Processing completed at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            print("")
    
    # Convert to dataframes and save results
    print("Saving results...")
    
    # Save combined results
    combined_df = pd.DataFrame(combined_results)
    if not combined_df.empty:
        combined_df.to_csv(f"{gene_dest}/all_target_variants.tsv", sep='\t', index=False)
    
    # Save gene-specific results
    for gene, variants in gene_variants.items():
        gene_df = pd.DataFrame(variants)
        if not gene_df.empty:
            gene_df.to_csv(f"{gene_dest}/variants/{gene}_variants.tsv", sep='\t', index=False)
    
    # Generate summary statistics
    print("Generating summary statistics...")
    
    # Summary by gene
    gene_summary = []
    for gene in target_genes:
        gene_df = pd.DataFrame(gene_variants[gene])
        if gene_df.empty:
            gene_summary.append({
                "Gene": gene,
                "Total_Variants": 0,
                "High_Impact": 0,
                "Moderate_Impact": 0,
                "Low_Impact": 0,
                "Modifier_Impact": 0
            })
        else:
            gene_summary.append({
                "Gene": gene,
                "Total_Variants": len(gene_df),
                "High_Impact": len(gene_df[gene_df["Impact"] == "HIGH"]),
                "Moderate_Impact": len(gene_df[gene_df["Impact"] == "MODERATE"]),
                "Low_Impact": len(gene_df[gene_df["Impact"] == "LOW"]),
                "Modifier_Impact": len(gene_df[gene_df["Impact"] == "MODIFIER"])
            })
    
    gene_summary_df = pd.DataFrame(gene_summary)
    gene_summary_df.to_csv(f"{gene_dest}/gene_summary.tsv", sep='\t', index=False)
    
    # Summary by sample
    if not combined_df.empty:
        sample_summary = []
        for sample in combined_df["Sample"].unique():
            sample_df = combined_df[combined_df["Sample"] == sample]
            sample_summary.append({
                "Sample": sample,
                "Total_Variants": len(sample_df),
                "High_Impact": len(sample_df[sample_df["Impact"] == "HIGH"]),
                "Moderate_Impact": len(sample_df[sample_df["Impact"] == "MODERATE"]),
                "Low_Impact": len(sample_df[sample_df["Impact"] == "LOW"]),
                "Modifier_Impact": len(sample_df[sample_df["Impact"] == "MODIFIER"])
            })
        
        sample_summary_df = pd.DataFrame(sample_summary)
        sample_summary_df.to_csv(f"{gene_dest}/sample_summary.tsv", sep='\t', index=False)
    
    print("Summary statistics saved to:")
    print(f"- {gene_dest}/gene_summary.tsv")
    print(f"- {gene_dest}/sample_summary.tsv")
    print("")
    
    print("=== Target Gene Extraction complete ===")

if __name__ == "__main__":
    main()
