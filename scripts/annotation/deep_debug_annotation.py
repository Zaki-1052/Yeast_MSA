#!/usr/bin/env python3
"""
deep_debug_annotation.py - Deep debugging of variant-gene annotation issues
"""

import os
import sys
import csv
import argparse
import gzip
from Bio import SeqIO
from collections import defaultdict

def deep_debug(vcf_file, genbank_dir, genes_of_interest_file, output_file):
    """
    Perform deep debugging of variant and gene annotation matching
    """
    # Step 1: Load genes of interest with detailed logging
    print("\n=== Step 1: Loading Genes of Interest ===")
    genes_of_interest = []
    
    with open(genes_of_interest_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            genes_of_interest.append(row)
            print(f"Gene: {row.get('w303_gene_id')}, SC Gene: {row.get('sc_gene_id')}")
            print(f"  Location: {row.get('scaffold_id')} at {row.get('start')}-{row.get('end')} ({row.get('strand')})")
            print(f"  JRIU ID: {row.get('jriu_id')}")
    
    print(f"Loaded {len(genes_of_interest)} genes of interest")
    
    # Step 2: Find a specific GenBank file for a gene of interest
    print("\n=== Step 2: Finding GenBank File for a Specific Gene ===")
    sample_gene = genes_of_interest[0]  # Use the first gene as an example
    sample_gene_id = sample_gene.get('w303_gene_id')
    sample_jriu = sample_gene.get('jriu_id')
    
    print(f"Looking for GenBank file containing JRIU ID: {sample_jriu}")
    
    genbank_file_for_sample = None
    for filename in os.listdir(genbank_dir):
        if not filename.endswith(('.genbank', '.gb', '.gbk')):
            continue
        
        filepath = os.path.join(genbank_dir, filename)
        try:
            record = SeqIO.read(filepath, "genbank")
            
            for feature in record.features:
                if feature.type == "source":
                    for note in feature.qualifiers.get('note', []):
                        if note == sample_jriu:
                            genbank_file_for_sample = filepath
                            print(f"Found matching GenBank file: {filename}")
                            break
            
            if genbank_file_for_sample:
                break
        except Exception as e:
            print(f"Error processing {filename}: {str(e)}")
    
    if not genbank_file_for_sample:
        print(f"ERROR: Could not find GenBank file for JRIU {sample_jriu}")
        return
    
    # Step 3: Parse the GenBank file to understand structure and gene locations
    print("\n=== Step 3: Analyzing GenBank File Structure ===")
    genbank_record = None
    try:
        genbank_record = SeqIO.read(genbank_file_for_sample, "genbank")
        
        print(f"GenBank Record ID: {genbank_record.id}")
        print(f"Length: {len(genbank_record.seq)} bp")
        print(f"Total features: {len(genbank_record.features)}")
        
        # Count feature types
        feature_types = defaultdict(int)
        for feature in genbank_record.features:
            feature_types[feature.type] += 1
        
        print("Feature types:")
        for ftype, count in feature_types.items():
            print(f"  {ftype}: {count}")
        
        # Find gene features
        gene_features = [f for f in genbank_record.features if f.type == "gene"]
        print(f"\nFound {len(gene_features)} gene features")
        
        # Print details of a few gene features
        print("Sample gene features:")
        for i, feature in enumerate(gene_features[:3]):
            gene_id = feature.qualifiers.get('gene', ['unknown'])[0]
            locus_tag = feature.qualifiers.get('locus_tag', ['unknown'])[0]
            location = feature.location
            print(f"  Gene {i+1}: {gene_id} / {locus_tag} at {location}")
            
            # Check notes for SC gene IDs
            if 'note' in feature.qualifiers:
                for note in feature.qualifiers['note']:
                    if 'similar to Saccharomyces cerevisiae' in note:
                        print(f"    SC note: {note}")
        
        # Look for our specific gene of interest
        print(f"\nLooking for gene {sample_gene_id} in GenBank file:")
        found_specific_gene = False
        
        for feature in gene_features:
            gene_id = feature.qualifiers.get('gene', ['unknown'])[0]
            locus_tag = feature.qualifiers.get('locus_tag', ['unknown'])[0]
            
            # Check direct matches
            if gene_id == sample_gene_id or locus_tag == sample_gene_id:
                found_specific_gene = True
                print(f"  FOUND by direct match: {locus_tag} at {feature.location}")
            
            # Check notes for SC gene ID
            elif 'note' in feature.qualifiers:
                sc_gene_id = sample_gene.get('sc_gene_id')
                for note in feature.qualifiers['note']:
                    if sc_gene_id in note:
                        found_specific_gene = True
                        print(f"  FOUND by SC gene ID in note: {locus_tag} at {feature.location}")
                        print(f"  Note: {note}")
        
        if not found_specific_gene:
            print(f"  WARNING: Could not find gene {sample_gene_id} in the GenBank file")
            print("  Showing all gene locus tags in this file:")
            for feature in gene_features[:10]:  # Show first 10
                locus_tag = feature.qualifiers.get('locus_tag', ['unknown'])[0]
                print(f"    {locus_tag}")
    
    except Exception as e:
        print(f"Error analyzing GenBank file: {str(e)}")
        return
    
    # Step 4: Check VCF to understand structure and coordinates
    print("\n=== Step 4: Analyzing VCF Structure and Coordinates ===")
    
    # Determine if the file is gzipped
    is_gzipped = vcf_file.endswith('.gz')
    opener = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'
    
    vcf_header = []
    vcf_variant_count = 0
    vcf_jriu_counts = defaultdict(int)
    
    print(f"Analyzing VCF file: {vcf_file}")
    
    with opener(vcf_file, mode) as f:
        # Process header lines
        for line in f:
            line = line.strip()
            if line.startswith('##'):
                vcf_header.append(line)
                # Check for reference information
                if 'reference=' in line:
                    print(f"Reference info: {line}")
            elif line.startswith('#CHROM'):
                vcf_header.append(line)
                break
        
        # Process variant lines
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            vcf_variant_count += 1
            fields = line.split('\t')
            
            # Extract JRIU ID
            chrom = fields[0]
            pos = int(fields[1])
            
            vcf_jriu_counts[chrom] += 1
            
            # Check if this variant is on our sample JRIU
            if chrom == sample_jriu:
                # Print details of a few variants on our sample JRIU
                if vcf_jriu_counts[chrom] <= 5:
                    print(f"Variant on {chrom} at position {pos}: {fields[3]} → {fields[4]}")
                    
                    # Check if this variant overlaps with our sample gene
                    gene_start = int(sample_gene.get('start', 0))
                    gene_end = int(sample_gene.get('end', 0))
                    
                    if gene_start <= pos <= gene_end:
                        print(f"  MATCH! Position {pos} is within gene {sample_gene_id} ({gene_start}-{gene_end})")
                    else:
                        print(f"  No overlap with gene {sample_gene_id} ({gene_start}-{gene_end})")
    
    print(f"\nProcessed {vcf_variant_count} variants on {len(vcf_jriu_counts)} unique JRIU IDs")
    print(f"Variants on sample JRIU {sample_jriu}: {vcf_jriu_counts.get(sample_jriu, 0)}")
    
    # Step 5: Cross-reference coordinates to check for potential issues
    print("\n=== Step 5: Cross-Referencing Coordinates ===")
    
    # Check if there are variants on the JRIU of our sample gene
    if sample_jriu in vcf_jriu_counts:
        print(f"Found {vcf_jriu_counts[sample_jriu]} variants on JRIU {sample_jriu}")
        
        # Parse the GenBank file again to check coordinate ranges
        if genbank_record:
            genbank_length = len(genbank_record.seq)
            print(f"GenBank record length for {sample_jriu}: {genbank_length} bp")
            
            # Check if any variants fall outside GenBank sequence length
            if vcf_jriu_counts[sample_jriu] > 0:
                outside_count = 0
                with opener(vcf_file, mode) as f:
                    for line in f:
                        if line.startswith('#'):
                            continue
                        
                        fields = line.strip().split('\t')
                        chrom = fields[0]
                        pos = int(fields[1])
                        
                        if chrom == sample_jriu and pos > genbank_length:
                            outside_count += 1
                
                if outside_count > 0:
                    print(f"WARNING: {outside_count} variants on {sample_jriu} have positions beyond the GenBank sequence length")
                    print("This suggests a reference version mismatch or coordinate system incompatibility")
                else:
                    print(f"All variants on {sample_jriu} are within the GenBank sequence length")
                    
                    # At this point, check the exact position ranges of variants vs. genes
                    variants_on_sample_jriu = []
                    with opener(vcf_file, mode) as f:
                        for line in f:
                            if line.startswith('#'):
                                continue
                            
                            fields = line.strip().split('\t')
                            chrom = fields[0]
                            pos = int(fields[1])
                            
                            if chrom == sample_jriu:
                                variants_on_sample_jriu.append(pos)
                    
                    # Sort positions for easier comparison
                    variants_on_sample_jriu.sort()
                    
                    # Get gene positions from GenBank
                    gene_positions = []
                    for feature in genbank_record.features:
                        if feature.type == "gene":
                            start = int(feature.location.start)
                            end = int(feature.location.end)
                            locus_tag = feature.qualifiers.get('locus_tag', ['unknown'])[0]
                            gene_positions.append((start, end, locus_tag))
                    
                    # Sort gene positions
                    gene_positions.sort()
                    
                    # Compare position ranges
                    print("\nPosition range comparison:")
                    print(f"Variant position range: {min(variants_on_sample_jriu)} - {max(variants_on_sample_jriu)}")
                    if gene_positions:
                        print(f"Gene position range: {gene_positions[0][0]} - {gene_positions[-1][1]}")
                        
                        # Check for overlap between variant and gene ranges
                        variant_min = min(variants_on_sample_jriu)
                        variant_max = max(variants_on_sample_jriu)
                        gene_min = min(gene_positions, key=lambda x: x[0])[0]
                        gene_max = max(gene_positions, key=lambda x: x[1])[1]
                        
                        if variant_max < gene_min or variant_min > gene_max:
                            print("CRITICAL ISSUE: No overlap between variant position range and gene position range")
                            print("This suggests a fundamental coordinate system incompatibility")
                        else:
                            print("Variant and gene ranges overlap, but specific positions might not match")
                            print("Try detailed position comparison between variants and genes")
    else:
        print(f"No variants found on JRIU {sample_jriu}")
        print("Check if the VCF file contains variants on JRIU IDs of other genes of interest")
        
        # Check other genes of interest
        found_any = False
        for gene in genes_of_interest:
            jriu_id = gene.get('jriu_id')
            if jriu_id in vcf_jriu_counts:
                found_any = True
                print(f"Found {vcf_jriu_counts[jriu_id]} variants on JRIU {jriu_id} for gene {gene.get('w303_gene_id')}")
        
        if not found_any:
            print("CRITICAL ISSUE: No variants found on any JRIU IDs of genes of interest")
            print("This suggests a complete mismatch between VCF and gene annotation references")
    
    # Step 6: Diagnose the issue and suggest solutions
    print("\n=== Step 6: Issue Diagnosis ===")
    
    # Check for JRIU and scaffold notations in VCF
    vcf_format_issue = True
    sample_jrius = list(vcf_jriu_counts.keys())[:5]
    for jriu in sample_jrius:
        if jriu.startswith('JRIU'):
            vcf_format_issue = False
            break
    
    if vcf_format_issue:
        print("ISSUE: VCF files might be using scaffold IDs instead of JRIU IDs")
        print(f"Sample chromosome IDs in VCF: {sample_jrius}")
        print("SOLUTION: Update VCF files to use JRIU IDs or modify annotation code to map scaffold IDs")
    
    # Alternatively, check for reference version mismatch
    if genbank_record and sample_jriu in vcf_jriu_counts:
        if any(gene_positions) and any(variants_on_sample_jriu):
            if max(variants_on_sample_jriu) > genbank_length:
                print("ISSUE: Position coordinates in VCF exceed GenBank sequence length")
                print("This indicates a reference version mismatch")
                print("SOLUTION: Ensure VCF and GenBank files are based on the same reference version")
            
            # Check for potential 0-based vs 1-based coordinate issue
            elif min(variants_on_sample_jriu) > 0 and max(variants_on_sample_jriu) <= genbank_length:
                if min(gene_positions, key=lambda x: x[0])[0] > max(variants_on_sample_jriu) or \
                   max(gene_positions, key=lambda x: x[1])[1] < min(variants_on_sample_jriu):
                    print("ISSUE: No overlap between variant and gene coordinate ranges despite both being within GenBank length")
                    print("This suggests a possible assembly version difference or coordinate system mismatch")
                    print("SOLUTION: Consider using a liftover tool to convert coordinates between versions")
    
    # Output summary to file
    with open(output_file, 'w') as f:
        f.write("=== Variant-Gene Annotation Debug Summary ===\n\n")
        f.write(f"VCF file: {vcf_file}\n")
        f.write(f"GenBank directory: {genbank_dir}\n")
        f.write(f"Genes of interest file: {genes_of_interest_file}\n\n")
        
        f.write(f"Total genes of interest: {len(genes_of_interest)}\n")
        f.write(f"Total variants in VCF: {vcf_variant_count}\n")
        f.write(f"Unique JRIU IDs in VCF: {len(vcf_jriu_counts)}\n\n")
        
        f.write("Key findings:\n")
        if vcf_format_issue:
            f.write("- VCF files may be using scaffold IDs instead of JRIU IDs\n")
        
        if genbank_record and sample_jriu in vcf_jriu_counts:
            if max(variants_on_sample_jriu) > genbank_length:
                f.write("- Position coordinates in VCF exceed GenBank sequence length, indicating reference version mismatch\n")
            
            if min(gene_positions, key=lambda x: x[0])[0] > max(variants_on_sample_jriu) or \
               max(gene_positions, key=lambda x: x[1])[1] < min(variants_on_sample_jriu):
                f.write("- No overlap between variant and gene coordinate ranges despite both being within GenBank length\n")
        
        if not any(jriu in vcf_jriu_counts for jriu in [gene.get('jriu_id') for gene in genes_of_interest]):
            f.write("- No variants found on any JRIU IDs of genes of interest, suggesting complete mismatch between references\n")
    
    print(f"\nDebug summary written to {output_file}")
    return

def main():
    parser = argparse.ArgumentParser(description='Deep debug variant gene annotation')
    parser.add_argument('--vcf', required=True, help='VCF file to analyze')
    parser.add_argument('--genbank-dir', required=True, help='Directory containing GenBank files')
    parser.add_argument('--genes-of-interest', required=True, help='Genes of interest file')
    parser.add_argument('--output', required=True, help='Output file for debug results')
    
    args = parser.parse_args()
    deep_debug(args.vcf, args.genbank_dir, args.genes_of_interest, args.output)

if __name__ == '__main__':
    main()