#!/usr/bin/env python3
# /Users/zakiralibhai/Documents/GitHub/Yeast_MSA/scripts/regulatory_analysis/debug_variants.py

"""
Debug script to examine specific variant characteristics in detail
"""

import pandas as pd
import sys
import os

def main():
    # Load the variant data
    variants_file = '/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/gene_variants_expanded/all_gene_variants.tsv'
    
    try:
        variants = pd.read_csv(variants_file, sep='\t')
        print(f"Loaded {len(variants)} variants from {variants_file}")
    except Exception as e:
        print(f"Error loading variants: {e}")
        sys.exit(1)
    
    # Get unique variants only (removing duplicates across samples)
    cols_for_uniqueness = ['Scaffold', 'Position', 'Ref', 'Alt', 'Gene_ID', 'Effect', 'Distance']
    unique_variants = variants.drop_duplicates(subset=cols_for_uniqueness)
    
    print(f"\nFound {len(unique_variants)} unique variants after removing duplicates")
    
    # Show the actual unique variants
    print("\nUnique variants:")
    for _, row in unique_variants.iterrows():
        print(f"Gene: {row['Gene_ID']}, Distance: {row['Distance']}bp, Effect: {row['Effect']}")
    
    # Analyze distance distribution
    print("\nDistance analysis:")
    distances = unique_variants['Distance'].dropna()
    print(f"Range: {distances.min()} to {distances.max()} bp")
    print(f"Unique distance values: {sorted(distances.unique())}")
    
    # Manual classification
    print("\nManual regulatory region classification:")
    
    # Define yeast regulatory regions
    regulatory_regions = {
        'core_promoter': (0, 150),        # Core promoter
        'proximal_UAS': (150, 500),       # Proximal UAS
        'distal_UAS': (500, 1500),        # Distal UAS
        'far_upstream': (1500, float('inf'))  # Far upstream
    }
    
    # Classify each unique variant
    for _, row in unique_variants.iterrows():
        distance = row['Distance']
        gene_id = row['Gene_ID']
        
        region = "unknown"
        for reg_region, (min_dist, max_dist) in regulatory_regions.items():
            if min_dist <= distance < max_dist:
                region = reg_region
                break
        
        print(f"Gene: {gene_id}, Distance: {distance}bp -> {region}")
    
    # Count by region
    region_counts = {}
    
    for _, row in unique_variants.iterrows():
        distance = row['Distance']
        
        region = "unknown"
        for reg_region, (min_dist, max_dist) in regulatory_regions.items():
            if min_dist <= distance < max_dist:
                region = reg_region
                break
        
        region_counts[region] = region_counts.get(region, 0) + 1
    
    print("\nRegion counts for unique variants:")
    for region, count in region_counts.items():
        percent = (count / len(unique_variants) * 100)
        print(f"{region}: {count} variants ({percent:.1f}%)")
    
    # Now let's check how this translates to the full dataset with duplicates
    all_regions = []
    for _, row in variants.iterrows():
        distance = row['Distance']
        
        region = "unknown"
        for reg_region, (min_dist, max_dist) in regulatory_regions.items():
            if min_dist <= distance < max_dist:
                region = reg_region
                break
        
        all_regions.append(region)
    
    variants['Region'] = all_regions
    
    # Calculate full dataset counts
    full_region_counts = variants['Region'].value_counts()
    
    print("\nRegion counts for full dataset (including duplicates):")
    for region, count in full_region_counts.items():
        percent = (count / len(variants) * 100)
        print(f"{region}: {count} variants ({percent:.1f}%)")
    
    # Check distribution across treatments
    print("\nRegion distribution by treatment:")
    treatment_region = pd.crosstab(
        variants['Treatment'],
        variants['Region'],
        normalize='index'
    ) * 100  # Convert to percentages
    
    print(treatment_region.round(1))
    
    # Output the specific variant details (for verification)
    output_dir = '/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/regulatory_analysis'
    output_file = os.path.join(output_dir, 'unique_variants_debug.tsv')
    
    unique_variants.to_csv(output_file, sep='\t', index=False)
    print(f"\nSaved unique variants to {output_file}")

if __name__ == "__main__":
    main()