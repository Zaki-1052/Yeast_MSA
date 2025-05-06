#!/usr/bin/env python3

"""
Filter variants to only include those on chromosomes with gene mapping data.

This script extracts variants from the VCF files and keeps only those variants on 
chromosomes that have gene mapping data available (CM007970.1, CM007971.1, CM007975.1, 
CM007976.1, CM007977.1).
"""

import os
import subprocess
import pandas as pd
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("filter_analysis.log", mode='w'),
        logging.StreamHandler()
    ]
)

# Define the chromosomes that have gene mapping data
GENE_MAPPED_CHROMOSOMES = {
    'CM007970.1', 'CM007971.1', 'CM007975.1', 'CM007976.1', 'CM007977.1'
}

def find_vcf_file(base_name, possible_locations):
    """Find a VCF file by trying multiple possible locations and formats."""
    for location in possible_locations:
        if os.path.exists(location.format(base_name)):
            return location.format(base_name)
    return None

def extract_variants_from_vcf(vcf_file):
    """Extract variants from a VCF file."""
    if not vcf_file or not os.path.exists(vcf_file):
        logging.error(f"Error: VCF file not found: {vcf_file}")
        return None
    
    try:
        # Extract variant information from VCF
        cmd = f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' {vcf_file}"
        output = subprocess.check_output(cmd, shell=True, universal_newlines=True)
        
        # Parse the output into a DataFrame
        lines = output.strip().split('\n')
        if not lines or not lines[0]:
            return pd.DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT'])
            
        data = [line.split('\t') for line in lines if line]
        variant_df = pd.DataFrame(data, columns=['CHROM', 'POS', 'REF', 'ALT'])
        
        # Convert position to numeric
        variant_df['POS'] = pd.to_numeric(variant_df['POS'])
        
        # Add variant ID column
        variant_df['Variant_ID'] = variant_df.apply(
            lambda row: f"{row['CHROM']}_{row['POS']}_{row['REF']}_{row['ALT']}", axis=1
        )
        
        # Count variants by chromosome
        chrom_counts = variant_df['CHROM'].value_counts().to_dict()
        logging.info(f"Extracted {len(variant_df)} variants from {vcf_file}")
        
        # Filter to only include chromosomes with gene mapping data
        filtered_df = variant_df[variant_df['CHROM'].isin(GENE_MAPPED_CHROMOSOMES)]
        
        # Log the results
        filtered_counts = filtered_df['CHROM'].value_counts().to_dict()
        logging.info(f"Filtered to {len(filtered_df)} variants on gene-mapped chromosomes")
        
        # Print chromosome statistics
        logging.info(f"Chromosome statistics:")
        for chrom in GENE_MAPPED_CHROMOSOMES:
            count = filtered_counts.get(chrom, 0)
            total = chrom_counts.get(chrom, 0)
            logging.info(f"  {chrom}: {count}/{total} variants (retained/total)")
        
        return filtered_df
    
    except Exception as e:
        logging.error(f"Error extracting variants from {vcf_file}: {e}")
        return None

def main():
    """Main function to filter variants."""
    logging.info("Starting variant filtering for gene-mapped chromosomes")
    
    # Define treatments and their controls
    treatments = {
        'WT-37-55-1': {'treatment': 'WT-37-55-1', 'control': 'WT-CTRL'},
        'WT-37-55-2': {'treatment': 'WT-37-55-2', 'control': 'WT-CTRL'},
        'WT-37-55-3': {'treatment': 'WT-37-55-3', 'control': 'WT-CTRL'},
        'WTA-55-1': {'treatment': 'WTA-55-1', 'control': 'WT-CTRL'},
        'WTA-55-2': {'treatment': 'WTA-55-2', 'control': 'WT-CTRL'},
        'WTA-55-3': {'treatment': 'WTA-55-3', 'control': 'WT-CTRL'},
        'STC-55-1': {'treatment': 'STC-55-1', 'control': 'WT-CTRL'},
        'STC-55-2': {'treatment': 'STC-55-2', 'control': 'WT-CTRL'},
        'STC-55-3': {'treatment': 'STC-55-3', 'control': 'WT-CTRL'},
        'CAS-55-1': {'treatment': 'CAS-55-1', 'control': 'WT-CTRL'},
        'CAS-55-2': {'treatment': 'CAS-55-2', 'control': 'WT-CTRL'},
        'CAS-55-3': {'treatment': 'CAS-55-3', 'control': 'WT-CTRL'},
    }
    
    # Define possible VCF locations
    vcf_locations = [
        "vcf/filtered/{}.filtered.vcf.gz",
        "vcf/annotated/{}.annotated.vcf",
        "vcf/individual/{}.vcf.gz",
        "vcf/individual/{}.vcf",
        "vcf/individual/{}.norm.vcf",
        "vcf/renamed/{}.filtered.renamed.vcf.gz",
    ]
    
    # Create output directory
    os.makedirs('filtered_variants', exist_ok=True)
    
    # Process each treatment
    for name, info in treatments.items():
        logging.info(f"Processing treatment: {name}")
        
        # Find treatment VCF
        treatment_file = find_vcf_file(info['treatment'], vcf_locations)
        if not treatment_file:
            treatment_file = find_vcf_file(name, vcf_locations)
        
        # Find control VCF
        control_file = find_vcf_file(info['control'], vcf_locations)
        
        logging.info(f"Treatment VCF: {treatment_file}")
        logging.info(f"Control VCF: {control_file}")
        
        # Extract and filter variants
        if treatment_file:
            treatment_variants = extract_variants_from_vcf(treatment_file)
            if treatment_variants is not None and len(treatment_variants) > 0:
                output_file = f"filtered_variants/{name}_treatment_filtered.csv"
                treatment_variants.to_csv(output_file, index=False)
                logging.info(f"Saved {len(treatment_variants)} filtered treatment variants to {output_file}")
        
        if control_file:
            control_variants = extract_variants_from_vcf(control_file)
            if control_variants is not None and len(control_variants) > 0:
                output_file = f"filtered_variants/{name}_control_filtered.csv"
                control_variants.to_csv(output_file, index=False)
                logging.info(f"Saved {len(control_variants)} filtered control variants to {output_file}")
    
    logging.info("Variant filtering complete")

if __name__ == "__main__":
    main()