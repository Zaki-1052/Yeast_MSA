#!/usr/bin/env python3
# /Users/zakiralibhai/Documents/GitHub/Yeast_MSA/scripts/specific_regulatory_analysis/fix_variant_annotations.py

"""
Fix variant regulatory annotations by merging data from the original treatment-specific variants file
with regulatory annotations, ensuring we analyze the complete set of treatment-specific variants.
"""

import os
import sys
import json
import pandas as pd
import numpy as np
from collections import defaultdict

def main():
    """Main function to fix variant annotations"""
    print("Starting to fix variant regulatory annotations...")
    
    # Define file paths
    scaffold_variants_file = '/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/gene_variants/all_gene_variants.tsv'
    gene_regulatory_map_file = '/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/new_regulatory_analysis/data/gene_regulatory_map.json'
    old_annotations_file = '/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/specific_regulatory_analysis/data/variant_regulatory_annotations.tsv'
    output_file = '/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/specific_regulatory_analysis/data/fixed_variant_regulatory_annotations.tsv'
    output_dir = '/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/specific_regulatory_analysis/data'
    
    # Make sure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Load original scaffolds variants file
    print("Loading scaffold variants file...")
    try:
        scaffold_variants = pd.read_csv(scaffold_variants_file, sep='\t')
        print(f"Loaded {len(scaffold_variants)} scaffold variants")
    except Exception as e:
        print(f"Error loading scaffold variants file: {e}")
        sys.exit(1)
    
    # Load gene regulatory map if available
    print("Loading gene regulatory map...")
    gene_map = None
    try:
        with open(gene_regulatory_map_file, 'r') as f:
            gene_map = json.load(f)
        print(f"Loaded gene regulatory map with {len(gene_map)} genes")
    except Exception as e:
        print(f"Warning: Could not load gene regulatory map: {e}")
        print("Will continue without gene regulatory information")
    
    # Load existing regulatory annotations (if available)
    print("Loading existing regulatory annotations...")
    existing_annotations = None
    try:
        existing_annotations = pd.read_csv(old_annotations_file, sep='\t')
        print(f"Loaded {len(existing_annotations)} existing regulatory annotations")
    except Exception as e:
        print(f"Warning: Could not load existing annotations: {e}")
        print("Will create new annotations from scratch")
    
    # Create a variant key for matching
    scaffold_variants['variant_key'] = scaffold_variants.apply(
        lambda row: f"{row['Scaffold']}:{row['Position']}:{row['Ref']}:{row['Alt']}", axis=1
    )
    
    # If we have existing annotations, create a lookup dictionary
    annotation_lookup = {}
    if existing_annotations is not None:
        existing_annotations['variant_key'] = existing_annotations.apply(
            lambda row: f"{row['Scaffold']}:{row['Position']}:{row['Ref']}:{row['Alt']}", axis=1
        )
        
        print("Creating annotation lookup from existing data...")
        for _, row in existing_annotations.iterrows():
            key = row['variant_key']
            annotation_lookup[key] = {
                'Regulatory_Region': row.get('Regulatory_Region', 'unknown'),
                'Regulatory_Type': row.get('Regulatory_Type', 'unknown'),
                'Conservation_Zone': row.get('Conservation_Zone', 'unknown'),
                'TFBS': row.get('TFBS', ''),
                'TFBS_Impact': row.get('TFBS_Impact', ''),
                'Motif_Disruption_Score': row.get('Motif_Disruption_Score', '')
            }
    
    # Function to determine regulatory region based on distance and effect
    def classify_regulatory_region(distance, effect=None):
        if effect and 'upstream_gene_variant' in effect:
            if distance <= 250:
                return 'core_promoter'
            elif distance <= 500:
                return 'UAS_proximal'
            else:
                return 'distal_regulatory'
        elif effect and 'downstream_gene_variant' in effect:
            if distance <= 250:
                return 'terminator'
            elif distance <= 500:
                return 'downstream_regulatory'
            else:
                return 'distal_regulatory'
        elif effect and 'missense_variant' in effect:
            return 'coding_missense'
        elif effect and 'synonymous_variant' in effect:
            return 'coding_synonymous'
        elif effect and 'intron_variant' in effect:
            return 'intronic'
        else:
            # Default based on distance only
            if distance <= 250:
                return 'core_promoter'
            elif distance <= 500:
                return 'UAS_proximal'
            else:
                return 'distal_regulatory'
    
    # Function to classify conservation zone
    def classify_conservation_zone(distance, position_relative):
        if position_relative == 'upstream':
            prefix = 'buffer_zone_'
        else:
            prefix = 'buffer_zone_'
            
        if 0 < distance <= 5000:
            return f"{prefix}upstream"
        elif 5000 < distance <= 50000:
            return 'intermediate_zone'
        elif 50000 < distance <= 100000:
            return 'satellite_zone'
        else:
            return 'distant_zone'
    
    # Add regulatory information to the scaffold variants
    print("Adding regulatory information to scaffold variants...")
    
    # Initialize new columns
    scaffold_variants['Regulatory_Region'] = 'unknown'
    scaffold_variants['Regulatory_Type'] = 'unknown'
    scaffold_variants['Conservation_Zone'] = 'unknown'
    scaffold_variants['TFBS'] = ''
    scaffold_variants['TFBS_Impact'] = ''
    scaffold_variants['Motif_Disruption_Score'] = ''
    
    # Update with existing annotations or infer new ones
    for i, row in scaffold_variants.iterrows():
        key = row['variant_key']
        
        # Check if we have existing annotation
        if key in annotation_lookup:
            annotation = annotation_lookup[key]
            scaffold_variants.at[i, 'Regulatory_Region'] = annotation['Regulatory_Region']
            scaffold_variants.at[i, 'Regulatory_Type'] = annotation['Regulatory_Type']
            scaffold_variants.at[i, 'Conservation_Zone'] = annotation['Conservation_Zone']
            scaffold_variants.at[i, 'TFBS'] = annotation['TFBS']
            scaffold_variants.at[i, 'TFBS_Impact'] = annotation['TFBS_Impact']
            scaffold_variants.at[i, 'Motif_Disruption_Score'] = annotation['Motif_Disruption_Score']
        else:
            # Infer regulatory information from distance and position
            if 'Distance' in row:
                distance = row['Distance']
                
                # For gene variants data, we assume 'upstream' since these are upstream variants
                # In gene variants, the Effect column tells us if it's upstream
                position_relative = 'upstream' if 'upstream' in row.get('Effect', '') else 'downstream'
                
                # Classify based on distance and effect
                if pd.notna(distance):
                    effect = row.get('Effect', '')
                    scaffold_variants.at[i, 'Regulatory_Region'] = classify_regulatory_region(distance, effect)
                    scaffold_variants.at[i, 'Conservation_Zone'] = classify_conservation_zone(distance, position_relative)
                    scaffold_variants.at[i, 'Regulatory_Type'] = position_relative
    
    # Remove the temporary variant_key column
    scaffold_variants = scaffold_variants.drop('variant_key', axis=1)
    
    # Save the result
    print(f"Saving {len(scaffold_variants)} fixed variant annotations to {output_file}...")
    scaffold_variants.to_csv(output_file, sep='\t', index=False)
    
    # Summary statistics
    control_count = len(scaffold_variants[scaffold_variants['Treatment'].str.endswith('-CTRL')])
    treatment_count = len(scaffold_variants) - control_count
    
    regulatory_counts = scaffold_variants['Regulatory_Region'].value_counts()
    
    print("\nSummary of fixed variant annotations:")
    print(f"Total variants: {len(scaffold_variants)}")
    print(f"Control variants: {control_count}")
    print(f"Treatment variants: {treatment_count}")
    print("\nRegulatory region distribution:")
    for region, count in regulatory_counts.items():
        print(f"  {region}: {count} ({count/len(scaffold_variants)*100:.1f}%)")
    
    print("\nFix completed successfully!")
    
if __name__ == "__main__":
    main()