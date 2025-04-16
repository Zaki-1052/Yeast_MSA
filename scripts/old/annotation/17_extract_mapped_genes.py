#!/usr/bin/env python3

# File: scripts/annotation/17_extract_mapped_genes.py
# Purpose: Extract variants affecting target genes using W303 identifiers

import os
import gzip
import re
import pandas as pd
from datetime import datetime

def parse_ann_field(ann_field, target_w303_ids):
    """Parse the ANN field and extract information for the target W303 gene IDs."""
    # Remove ANN= prefix
    ann = ann_field[4:] if ann_field.startswith("ANN=") else ann_field
    
    # Split into individual annotations
    annotations = ann.split(",")
    
    # Find annotations that match our target gene
    matched_anns = []
    for a in annotations:
        fields = a.split("|")
        if len(fields) < 4:
            continue
        
        gene_id = fields[3]
        if gene_id in target_w303_ids:
            matched_anns.append({
                "effect": fields[1] if len(fields) > 1 else "",
                "impact": fields[2] if len(fields) > 2 else "",
                "gene_id": gene_id,
                "protein_change": fields[10] if len(fields) > 10 else "",
                "w303_id": gene_id
            })
    
    return matched_anns

def main():
    print("=== Extracting Variants in Target Genes (Using W303 Mappings) ===")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("")
    
    # Define directories
    anno_source = "annotation/results_renamed"
    gene_dest = "annotation/genes_of_interest"
    mapped_dir = "annotation/mapped_genes"
    
    # Create output directory
    os.makedirs(f"{gene_dest}/variants", exist_ok=True)
    
    # Load gene mappings
    mapping_file = f"{mapped_dir}/simple_gene_mapping.tsv"
    if not os.path.exists(mapping_file):
        print(f"ERROR: Mapping file not found: {mapping_file}")
        return
    
    mappings_df = pd.read_csv(mapping_file, sep='\t')
    
    # Create a lookup from W303 to SGD
    w303_to_sgd = {}
    target_w303_ids = []
    
    for _, row in mappings_df.iterrows():
        if row['W303_ID'] != 'NOT_FOUND':
            w303_to_sgd[row['W303_ID']] = row['SGD_Gene']
            target_w303_ids.append(row['W303_ID'])
    
    print(f"Loaded {len(target_w303_ids)} mapped W303 identifiers")
    print("")
    
    # Create a combined results dataframe
    combined_results = []
    
    # Create gene-specific dataframes
    gene_variants = {sgd_gene: [] for sgd_gene in mappings_df['SGD_Gene']}
    
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
    log_file = f"{gene_dest}/gene_extraction_mapped.log"
    with open(log_file, 'w') as log:
        log.write("Gene Extraction Log (Using W303 Mappings)\n")
        log.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        log.write("==============================================\n\n")
        
        # Process each annotated VCF file
        for vcf_file in vcf_files:
            filename = os.path.basename(vcf_file)
            sample = filename.replace(".snpeff.vcf.gz", "")
            
            print(f"Processing {sample}...")
            log.write(f"[{sample}] Processing started at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            
            # Track variants per gene for this sample
            gene_variant_counts = {w303_id: 0 for w303_id in target_w303_ids}
            
            # Parse the VCF file
            with gzip.open(vcf_file, 'rt') as f:
                for line in f:
                    if line.startswith('#'):
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
                    ann_data_list = parse_ann_field(ann_field, target_w303_ids)
                    
                    for ann_data in ann_data_list:
                        w303_id = ann_data["w303_id"]
                        sgd_gene = w303_to_sgd.get(w303_id, "Unknown")
                        
                        gene_variant_counts[w303_id] += 1
                        
                        # Create a record
                        record = {
                            "Sample": sample,
                            "SGD_Gene": sgd_gene,
                            "W303_ID": w303_id,
                            "Chromosome": chrom,
                            "Position": pos,
                            "Ref": ref,
                            "Alt": alt,
                            "Effect": ann_data["effect"],
                            "Impact": ann_data["impact"],
                            "Protein_Change": ann_data["protein_change"]
                        }
                        
                        # Add to combined results
                        combined_results.append(record)
                        
                        # Add to gene-specific results
                        gene_variants[sgd_gene].append(record)
            
            # Report variants found per gene
            for w303_id, count in gene_variant_counts.items():
                if count > 0:
                    sgd_gene = w303_to_sgd.get(w303_id, "Unknown")
                    print(f"  Found {count} variants affecting {sgd_gene} (W303: {w303_id})")
                    log.write(f"[{sample}] Found {count} variants affecting {sgd_gene} (W303: {w303_id})\n")
            
            print(f"Done processing {sample}")
            log.write(f"[{sample}] Processing completed at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            print("")
    
    # Convert to dataframes and save results
    print("Saving results...")
    
    # Save combined results
    combined_df = pd.DataFrame(combined_results)
    if not combined_df.empty:
        combined_df.to_csv(f"{gene_dest}/all_target_variants_mapped.tsv", sep='\t', index=False)
        print(f"Saved {len(combined_df)} variants to all_target_variants_mapped.tsv")
    else:
        print("No variants found for target genes")
    
    # Save gene-specific results
    for sgd_gene, variants in gene_variants.items():
        gene_df = pd.DataFrame(variants)
        if not gene_df.empty:
            gene_df.to_csv(f"{gene_dest}/variants/{sgd_gene}_variants_mapped.tsv", sep='\t', index=False)
            print(f"Saved {len(gene_df)} variants for {sgd_gene}")
    
    # Generate summary statistics
    print("Generating summary statistics...")
    
    # Summary by gene
    gene_summary = []
    for sgd_gene in mappings_df['SGD_Gene']:
        gene_df = pd.DataFrame(gene_variants[sgd_gene])
        if gene_df.empty:
            gene_summary.append({
                "SGD_Gene": sgd_gene,
                "W303_ID": next((m['W303_ID'] for _, m in mappings_df.iterrows() if m['SGD_Gene'] == sgd_gene), "NOT_FOUND"),
                "Total_Variants": 0,
                "High_Impact": 0,
                "Moderate_Impact": 0,
                "Low_Impact": 0,
                "Modifier_Impact": 0,
                "Samples_With_Variants": 0
            })
        else:
            w303_id = gene_df['W303_ID'].iloc[0]
            unique_samples = len(gene_df['Sample'].unique())
            gene_summary.append({
                "SGD_Gene": sgd_gene,
                "W303_ID": w303_id,
                "Total_Variants": len(gene_df),
                "High_Impact": len(gene_df[gene_df["Impact"] == "HIGH"]),
                "Moderate_Impact": len(gene_df[gene_df["Impact"] == "MODERATE"]),
                "Low_Impact": len(gene_df[gene_df["Impact"] == "LOW"]),
                "Modifier_Impact": len(gene_df[gene_df["Impact"] == "MODIFIER"]),
                "Samples_With_Variants": unique_samples
            })
    
    gene_summary_df = pd.DataFrame(gene_summary)
    gene_summary_df.to_csv(f"{gene_dest}/gene_summary_mapped.tsv", sep='\t', index=False)
    print(f"Generated gene summary with {len(gene_summary_df)} genes")
    
    # Summary by sample
    if not combined_df.empty:
        sample_summary = []
        for sample in combined_df["Sample"].unique():
            sample_df = combined_df[combined_df["Sample"] == sample]
            unique_genes = len(sample_df['SGD_Gene'].unique())
            sample_summary.append({
                "Sample": sample,
                "Total_Variants": len(sample_df),
                "High_Impact": len(sample_df[sample_df["Impact"] == "HIGH"]),
                "Moderate_Impact": len(sample_df[sample_df["Impact"] == "MODERATE"]),
                "Low_Impact": len(sample_df[sample_df["Impact"] == "LOW"]),
                "Modifier_Impact": len(sample_df[sample_df["Impact"] == "MODIFIER"]),
                "Genes_With_Variants": unique_genes
            })
        
        sample_summary_df = pd.DataFrame(sample_summary)
        sample_summary_df.to_csv(f"{gene_dest}/sample_summary_mapped.tsv", sep='\t', index=False)
        print(f"Generated sample summary with {len(sample_summary_df)} samples")
    
    # Generate report of common variants
    if not combined_df.empty:
        # Group variants by position and gene
        grouped = combined_df.groupby(['SGD_Gene', 'Chromosome', 'Position', 'Ref', 'Alt'])
        common_variants = []
        
        for (gene, chrom, pos, ref, alt), group in grouped:
            if len(group) > 1:  # Found in multiple samples
                samples = ', '.join(group['Sample'].unique())
                effect = group.iloc[0]['Effect']
                impact = group.iloc[0]['Impact']
                protein = group.iloc[0]['Protein_Change']
                
                common_variants.append({
                    'SGD_Gene': gene,
                    'W303_ID': group.iloc[0]['W303_ID'],
                    'Chromosome': chrom,
                    'Position': pos,
                    'Ref': ref,
                    'Alt': alt,
                    'Effect': effect,
                    'Impact': impact,
                    'Protein_Change': protein,
                    'Sample_Count': len(group['Sample'].unique()),
                    'Samples': samples
                })
        
        if common_variants:
            common_df = pd.DataFrame(common_variants)
            common_df = common_df.sort_values(['Sample_Count', 'SGD_Gene'], ascending=[False, True])
            common_df.to_csv(f"{gene_dest}/common_variants_mapped.tsv", sep='\t', index=False)
            print(f"Generated report of {len(common_df)} common variants")
    
    print("")
    print("Summary statistics saved to:")
    print(f"- {gene_dest}/gene_summary_mapped.tsv")
    print(f"- {gene_dest}/sample_summary_mapped.tsv")
    print(f"- {gene_dest}/common_variants_mapped.tsv")
    print("")
    
    print("=== Target Gene Extraction (Using W303 Mappings) complete ===")

if __name__ == "__main__":
    main()
