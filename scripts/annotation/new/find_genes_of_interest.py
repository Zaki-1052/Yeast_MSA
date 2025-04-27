#!/usr/bin/env python3
"""
find_genes_of_interest.py

This script identifies the 11 target genes in the ergosterol biosynthesis pathway
within the GenBank annotation files and extracts their information.

Usage:
    python find_genes_of_interest.py --genbank_dir <genbank_directory> --output_dir <output_directory>
"""

import os
import argparse
from Bio import SeqIO
import pandas as pd
import re

# Define the genes of interest
GENES_OF_INTEREST = [
    "YHR190W",  # ERG9 - Squalene synthase
    "YGR175C",  # ERG1 - Squalene epoxidase
    "YHR072W",  # ERG7 - Lanosterol synthase
    "YHR007C",  # ERG11 - Lanosterol 14-alpha-demethylase
    "YNL280C",  # ERG24 - C-14 sterol reductase
    "YGR060W",  # ERG25 - C-4 methyl sterol oxidase
    "YML008C",  # ERG6 - Delta(24)-sterol C-methyltransferase
    "YMR202W",  # ERG2 - C-8 sterol isomerase
    "YLR056W",  # ERG3 - C-5 sterol desaturase
    "YMR015C",  # ERG5 - C-22 sterol desaturase
    "YGL012W",  # ERG4 - C-24(28) sterol reductase
]

# Common ergosterol gene names mapping
ERG_NAMES = {
    "YHR190W": "ERG9",
    "YGR175C": "ERG1",
    "YHR072W": "ERG7",
    "YHR007C": "ERG11",
    "YNL280C": "ERG24",
    "YGR060W": "ERG25",
    "YML008C": "ERG6",
    "YMR202W": "ERG2",
    "YLR056W": "ERG3",
    "YMR015C": "ERG5",
    "YGL012W": "ERG4",
}

def find_genes_in_genbank(genbank_dir):
    """
    Parse GenBank files to find genes of interest.
    
    Args:
        genbank_dir (str): Directory containing GenBank files
        
    Returns:
        list: List of dictionaries with gene information
    """
    gene_data = []
    
    # List all GenBank files
    genbank_files = [f for f in os.listdir(genbank_dir) if f.endswith('.genbank')]
    
    if not genbank_files:
        print(f"No GenBank files found in {genbank_dir}")
        return gene_data
    
    print(f"Searching for genes of interest in {len(genbank_files)} GenBank files")
    print(f"Target genes: {', '.join(GENES_OF_INTEREST)}")
    
    # Process each GenBank file
    for genbank_file in sorted(genbank_files):
        file_path = os.path.join(genbank_dir, genbank_file)
        print(f"Processing {genbank_file}...")
        
        try:
            # Parse GenBank file
            record = SeqIO.read(file_path, "genbank")
            scaffold_name = record.name  # w303_scaffold_X
            
            # Extract chromosome identifier
            chromosome_id = None
            for feature in record.features:
                if feature.type == "source":
                    notes = feature.qualifiers.get("note", [])
                    for note in notes:
                        # Look for CM007XXX.1 or LYZE01XXXXXX.1 pattern
                        if (note.startswith("CM007") and "." in note) or (note.startswith("LYZE") and "." in note):
                            chromosome_id = note.strip()
                            break
                    if chromosome_id:
                        break
            
            # Look for genes of interest
            for feature in record.features:
                if feature.type == "gene":
                    # Extract information from gene feature
                    gene_id = feature.qualifiers.get("gene", [""])[0]
                    
                    # Skip if this doesn't look like a W303 gene ID
                    if not gene_id.startswith("W303"):
                        continue
                    
                    # Check CDS feature for more information
                    for next_feature in record.features:
                        if next_feature.type == "CDS" and next_feature.qualifiers.get("gene", [""])[0] == gene_id:
                            
                            # Check inference and note fields for genes of interest
                            inferences = next_feature.qualifiers.get("inference", [])
                            notes = next_feature.qualifiers.get("note", [])
                            sc_gene_match = None
                            erg_name = None
                            
                            # Look for matches in inference field
                            for inference in inferences:
                                for target_gene in GENES_OF_INTEREST:
                                    if target_gene in inference:
                                        sc_gene_match = target_gene
                                        erg_name = ERG_NAMES.get(target_gene)
                                        break
                                if sc_gene_match:
                                    break
                            
                            # If not found in inference, try note field
                            if not sc_gene_match:
                                for note in notes:
                                    for target_gene in GENES_OF_INTEREST:
                                        if target_gene in note:
                                            sc_gene_match = target_gene
                                            erg_name = ERG_NAMES.get(target_gene)
                                            break
                                    if sc_gene_match:
                                        break
                            
                            # If we found a match, extract full gene information
                            if sc_gene_match:
                                # Get gene location and orientation
                                location = feature.location
                                start = int(location.start)
                                end = int(location.end)
                                strand = "+" if location.strand == 1 else "-"
                                
                                # Get gene product if available
                                product = next_feature.qualifiers.get("product", ["unknown"])[0]
                                
                                # Combine annotations from note field
                                annotation = "; ".join(notes)
                                
                                # Create gene entry
                                entry = {
                                    "w303_gene_id": gene_id,
                                    "locus_tag": next_feature.qualifiers.get("locus_tag", [""])[0],
                                    "sc_gene_id": sc_gene_match,
                                    "erg_name": erg_name,
                                    "w303_scaffold": scaffold_name,
                                    "chromosome_id": chromosome_id,
                                    "start": start,
                                    "end": end,
                                    "length": end - start,
                                    "strand": strand,
                                    "product": product,
                                    "annotation": annotation
                                }
                                
                                gene_data.append(entry)
                                print(f"  Found {sc_gene_match} ({erg_name}) on {scaffold_name} ({chromosome_id})")
                                break
        
        except Exception as e:
            print(f"  ERROR processing {genbank_file}: {str(e)}")
    
    return gene_data

def create_gene_mapping_files(gene_data, output_dir):
    """
    Create gene mapping files from the extracted data.
    
    Args:
        gene_data (list): List of dictionaries with gene information
        output_dir (str): Directory to save output files
    """
    if not gene_data:
        print("No genes of interest found")
        return
    
    # Create DataFrame
    df = pd.DataFrame(gene_data)
    
    # Sort by SC gene ID
    df = df.sort_values('sc_gene_id')
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Create gene mapping file
    gene_mapping_file = os.path.join(output_dir, "gene_mapping.tsv")
    df.to_csv(gene_mapping_file, sep='\t', index=False)
    print(f"Created gene mapping file: {gene_mapping_file}")
    
    # Create gene of interest mapping
    goi_mapping_file = os.path.join(output_dir, "genes_of_interest_mapping.tsv")
    goi_columns = [
        'sc_gene_id', 'erg_name', 'w303_gene_id', 'locus_tag',
        'w303_scaffold', 'chromosome_id', 'start', 'end', 'strand', 'product'
    ]
    
    # Ensure we have all columns
    for col in goi_columns:
        if col not in df.columns:
            df[col] = ""
    
    df[goi_columns].to_csv(goi_mapping_file, sep='\t', index=False)
    print(f"Created genes of interest mapping file: {goi_mapping_file}")
    
    # Print summary
    found_genes = set(df['sc_gene_id'].tolist())
    missing_genes = set(GENES_OF_INTEREST) - found_genes
    
    print(f"\nSummary: Found {len(found_genes)} out of {len(GENES_OF_INTEREST)} target genes")
    
    if missing_genes:
        print(f"\nMissing genes: {', '.join(sorted(missing_genes))}")
    
    print("\nFound genes:")
    for gene_id in sorted(found_genes):
        gene_rows = df[df['sc_gene_id'] == gene_id]
        for _, row in gene_rows.iterrows():
            print(f"  {gene_id} ({row['erg_name']}): {row['w303_gene_id']} on {row['w303_scaffold']} ({row['chromosome_id']})")

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Find genes of interest in GenBank annotations")
    parser.add_argument("--genbank_dir", required=True, help="Directory containing GenBank files")
    parser.add_argument("--output_dir", required=True, help="Directory to save output files")
    args = parser.parse_args()
    
    # Find genes of interest
    gene_data = find_genes_in_genbank(args.genbank_dir)
    
    # Create gene mapping files
    create_gene_mapping_files(gene_data, args.output_dir)

if __name__ == "__main__":
    main()