#!/usr/bin/env python3
"""
analyze_genes_of_interest.py

This script extracts and analyzes variants in the 11 target genes from annotated VCF files.

Usage:
    python analyze_genes_of_interest.py --annotated_dir <annotated_vcf_directory> \
                                       --gene_mapping <gene_mapping_file> \
                                       --output_dir <output_directory>
"""

import os
import argparse
import pandas as pd
import re
import subprocess
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns

# Define sample to treatment mapping
SAMPLE_TREATMENTS = {
    "CAS-55-1": "CAS",
    "CAS-55-2": "CAS",
    "CAS-55-3": "CAS",
    "CAS-CTRL": "CAS-CTRL",
    "STC-55-1": "STC",
    "STC-55-2": "STC",
    "STC-55-3": "STC",
    "STC-CTRL": "STC-CTRL",
    "WT-37-55-1": "WT-37",
    "WT-37-55-2": "WT-37",
    "WT-37-55-3": "WT-37",
    "WT-CTRL": "WT-CTRL",
    "WTA-55-1": "WTA",
    "WTA-55-2": "WTA",
    "WTA-55-3": "WTA"
}

# Define treatment groups for comparison
TREATMENT_GROUPS = {
    "Temperature Adaptation": ["WT-37", "WT-CTRL"],
    "Low Oxygen Adaptation": ["WTA", "WT-CTRL"],
    "STC Effect": ["STC", "WTA"],
    "CAS Effect": ["CAS", "WT-37"]
}

# Impact severity levels
IMPACT_LEVELS = {
    "HIGH": 4,
    "MODERATE": 3,
    "LOW": 2,
    "MODIFIER": 1
}

def read_gene_mapping(mapping_file):
    """Read gene mapping file to get genes of interest."""
    df = pd.read_csv(mapping_file, sep='\t')
    # Create a mapping of W303 gene ID to standard gene ID and name
    gene_mapping = {}
    for _, row in df.iterrows():
        w303_id = row['w303_gene_id']
        sc_id = row['sc_gene_id']
        erg_name = row['erg_name'] if 'erg_name' in row else sc_id
        gene_mapping[w303_id] = {
            'sc_id': sc_id,
            'erg_name': erg_name,
            'scaffold': row['w303_scaffold'],
            'start': row['start'],
            'end': row['end'],
            'strand': row['strand']
        }
    
    print(f"Read mapping for {len(gene_mapping)} genes of interest")
    return gene_mapping

def extract_variants_in_genes(vcf_file, gene_mapping):
    """Extract variants in genes of interest from an annotated VCF file."""
    sample_name = os.path.basename(vcf_file).split('.')[0]
    print(f"Processing {sample_name}...")
    
    # Collect variants by gene
    variants_by_gene = defaultdict(list)
    gene_pattern = re.compile(r'gene=([^|;]+)')
    impact_pattern = re.compile(r'IMPACT=([^|;]+)')
    effect_pattern = re.compile(r'ANN=.*?\|.*?\|([^|]+)\|')
    
    # Open VCF file
    try:
        # Check if VCF is gzipped
        is_gzipped = vcf_file.endswith('.gz')
        command = f"{'gunzip -c' if is_gzipped else 'cat'} {vcf_file} | grep -v '^#'"
        
        # Run command to get non-header lines
        process = subprocess.run(command, shell=True, capture_output=True, text=True)
        
        if process.returncode != 0:
            print(f"  Error reading VCF file: {process.stderr}")
            return None
        
        # Process each variant line
        for line in process.stdout.splitlines():
            if not line.strip():
                continue
                
            # Extract basic variant information
            fields = line.split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            info = fields[7]
            
            # Extract gene information
            gene_match = gene_pattern.search(info)
            if not gene_match:
                continue
                
            gene_id = gene_match.group(1)
            
            # Check if this is one of our genes of interest
            if gene_id not in gene_mapping:
                continue
            
            # Extract impact
            impact_match = impact_pattern.search(info)
            impact = impact_match.group(1) if impact_match else "UNKNOWN"
            
            # Extract effect
            effect_match = effect_pattern.search(info)
            effect = effect_match.group(1) if effect_match else "UNKNOWN"
            
            # Add variant to list
            variant = {
                'sample': sample_name,
                'treatment': SAMPLE_TREATMENTS.get(sample_name, "UNKNOWN"),
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'gene_id': gene_id,
                'sc_id': gene_mapping[gene_id]['sc_id'],
                'erg_name': gene_mapping[gene_id]['erg_name'],
                'impact': impact,
                'effect': effect,
                'impact_value': IMPACT_LEVELS.get(impact, 0)
            }
            
            variants_by_gene[gene_id].append(variant)
        
        # Count variants by gene
        gene_counts = {gene: len(variants) for gene, variants in variants_by_gene.items()}
        total_variants = sum(gene_counts.values())
        
        print(f"  Found {total_variants} variants in genes of interest")
        for gene, count in gene_counts.items():
            sc_id = gene_mapping[gene]['sc_id']
            erg_name = gene_mapping[gene]['erg_name']
            print(f"    {sc_id} ({erg_name}): {count} variants")
        
        return variants_by_gene
        
    except Exception as e:
        print(f"  Error processing {vcf_file}: {str(e)}")
        return None

def combine_all_variants(vcf_dir, gene_mapping):
    """Process all VCF files and combine variants from genes of interest."""
    # Find annotated VCF files
    vcf_files = []
    for file in os.listdir(vcf_dir):
        if file.endswith('.annotated.vcf') or file.endswith('.annotated.vcf.gz'):
            vcf_files.append(os.path.join(vcf_dir, file))
    
    if not vcf_files:
        print(f"No annotated VCF files found in {vcf_dir}")
        return None
    
    print(f"Found {len(vcf_files)} annotated VCF files")
    
    # Process each VCF file
    all_variants = []
    for vcf_file in sorted(vcf_files):
        variants_by_gene = extract_variants_in_genes(vcf_file, gene_mapping)
        if variants_by_gene:
            for gene, variants in variants_by_gene.items():
                all_variants.extend(variants)
    
    # Convert to DataFrame
    if all_variants:
        df = pd.DataFrame(all_variants)
        print(f"Combined dataset contains {len(df)} variants across {df['gene_id'].nunique()} genes")
        return df
    else:
        print("No variants found in genes of interest")
        return None

def create_summary_tables(variants_df, output_dir):
    """Create summary tables for variants in genes of interest."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Table 1: Variants by gene and sample
    gene_sample_counts = variants_df.groupby(['sc_id', 'erg_name', 'sample']).size().reset_index(name='variant_count')
    gene_sample_counts['treatment'] = gene_sample_counts['sample'].map(SAMPLE_TREATMENTS)
    gene_sample_table = os.path.join(output_dir, "gene_sample_variants.tsv")
    gene_sample_counts.to_csv(gene_sample_table, sep='\t', index=False)
    print(f"Created gene-sample variant count table: {gene_sample_table}")
    
    # Table 2: Variants by gene and treatment
    gene_treatment_counts = variants_df.groupby(['sc_id', 'erg_name', 'treatment']).size().reset_index(name='variant_count')
    gene_treatment_table = os.path.join(output_dir, "gene_treatment_variants.tsv")
    gene_treatment_counts.to_csv(gene_treatment_table, sep='\t', index=False)
    print(f"Created gene-treatment variant count table: {gene_treatment_table}")
    
    # Table 3: Impact distribution by gene
    impact_counts = variants_df.groupby(['sc_id', 'erg_name', 'impact']).size().reset_index(name='count')
    impact_table = os.path.join(output_dir, "gene_impact_distribution.tsv")
    impact_counts.to_csv(impact_table, sep='\t', index=False)
    print(f"Created impact distribution table: {impact_table}")
    
    # Table 4: Effect types by gene
    effect_counts = variants_df.groupby(['sc_id', 'erg_name', 'effect']).size().reset_index(name='count')
    effect_table = os.path.join(output_dir, "gene_effect_distribution.tsv")
    effect_counts.to_csv(effect_table, sep='\t', index=False)
    print(f"Created effect distribution table: {effect_table}")
    
    # Table 5: Treatment comparison
    treatment_comparison = {}
    for comparison_name, groups in TREATMENT_GROUPS.items():
        if len(groups) >= 2:
            group1, group2 = groups[:2]
            group1_variants = variants_df[variants_df['treatment'] == group1]
            group2_variants = variants_df[variants_df['treatment'] == group2]
            
            # Count unique variants by gene in each group
            g1_counts = group1_variants.groupby('sc_id').size()
            g2_counts = group2_variants.groupby('sc_id').size()
            
            # Compare
            for gene in set(g1_counts.index) | set(g2_counts.index):
                count1 = g1_counts.get(gene, 0)
                count2 = g2_counts.get(gene, 0)
                if gene not in treatment_comparison:
                    treatment_comparison[gene] = {}
                treatment_comparison[gene][f"{group1}_count"] = count1
                treatment_comparison[gene][f"{group2}_count"] = count2
                treatment_comparison[gene][f"{comparison_name}_ratio"] = count1 / count2 if count2 > 0 else float('inf')
    
    # Convert to DataFrame
    if treatment_comparison:
        treatment_df = pd.DataFrame.from_dict(treatment_comparison, orient='index')
        treatment_df.index.name = 'sc_id'
        treatment_df = treatment_df.reset_index()
        
        # Add gene names
        gene_names = variants_df[['sc_id', 'erg_name']].drop_duplicates().set_index('sc_id')['erg_name']
        treatment_df['erg_name'] = treatment_df['sc_id'].map(gene_names)
        
        # Save to file
        treatment_table = os.path.join(output_dir, "treatment_comparison.tsv")
        treatment_df.to_csv(treatment_table, sep='\t', index=False)
        print(f"Created treatment comparison table: {treatment_table}")
    
    # Create comprehensive gene summary
    gene_summary = variants_df.groupby(['sc_id', 'erg_name']).agg(
        variant_count=('pos', 'size'),
        samples=('sample', 'nunique'),
        treatments=('treatment', 'nunique'),
        high_impact=('impact', lambda x: sum(x == 'HIGH')),
        moderate_impact=('impact', lambda x: sum(x == 'MODERATE')),
        low_impact=('impact', lambda x: sum(x == 'LOW')),
        modifier_impact=('impact', lambda x: sum(x == 'MODIFIER'))
    ).reset_index()
    
    gene_summary_table = os.path.join(output_dir, "gene_summary.tsv")
    gene_summary.to_csv(gene_summary_table, sep='\t', index=False)
    print(f"Created gene summary table: {gene_summary_table}")
    
    return {
        'gene_sample': gene_sample_counts,
        'gene_treatment': gene_treatment_counts,
        'impact': impact_counts,
        'effect': effect_counts,
        'treatment_comparison': treatment_df if treatment_comparison else None,
        'gene_summary': gene_summary
    }

def create_visualizations(summary_tables, output_dir):
    """Create visualizations for variants in genes of interest."""
    viz_dir = os.path.join(output_dir, "visualizations")
    os.makedirs(viz_dir, exist_ok=True)
    
    # Set plot style
    plt.style.use('seaborn-v0_8-whitegrid')
    
    # 1. Heatmap of variant counts by gene and treatment
    if 'gene_treatment' in summary_tables:
        gene_treatment = summary_tables['gene_treatment']
        pivot_data = gene_treatment.pivot_table(
            index='erg_name', 
            columns='treatment', 
            values='variant_count',
            fill_value=0
        )
        
        plt.figure(figsize=(10, 8))
        sns.heatmap(pivot_data, annot=True, cmap='YlGnBu', fmt='d')
        plt.title('Variant Counts by Gene and Treatment')
        plt.xlabel('Treatment')
        plt.ylabel('Gene')
        plt.tight_layout()
        heatmap_file = os.path.join(viz_dir, "gene_treatment_heatmap.png")
        plt.savefig(heatmap_file, dpi=300)
        plt.close()
        print(f"Created gene-treatment heatmap: {heatmap_file}")
    
    # 2. Bar chart of variant counts by gene
    if 'gene_summary' in summary_tables:
        gene_summary = summary_tables['gene_summary']
        plt.figure(figsize=(12, 6))
        bars = plt.bar(gene_summary['erg_name'], gene_summary['variant_count'])
        
        # Add count labels on top of bars
        for bar in bars:
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                    f'{int(height)}', ha='center', va='bottom')
        
        plt.title('Total Variants by Gene')
        plt.xlabel('Gene')
        plt.ylabel('Variant Count')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        bar_file = os.path.join(viz_dir, "gene_variant_counts.png")
        plt.savefig(bar_file, dpi=300)
        plt.close()
        print(f"Created gene variant count bar chart: {bar_file}")
    
    # 3. Impact distribution by gene
    if 'impact' in summary_tables:
        impact = summary_tables['impact']
        
        # Order impacts by severity
        impact_order = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']
        impact['impact'] = pd.Categorical(impact['impact'], categories=impact_order, ordered=True)
        
        plt.figure(figsize=(12, 8))
        pivot_data = impact.pivot_table(
            index='erg_name', 
            columns='impact', 
            values='count',
            fill_value=0
        )
        
        # Plot a stacked bar chart
        pivot_data.plot(kind='bar', stacked=True, colormap='viridis', figsize=(12, 8))
        plt.title('Impact Distribution by Gene')
        plt.xlabel('Gene')
        plt.ylabel('Variant Count')
        plt.xticks(rotation=45, ha='right')
        plt.legend(title='Impact')
        plt.tight_layout()
        impact_file = os.path.join(viz_dir, "gene_impact_distribution.png")
        plt.savefig(impact_file, dpi=300)
        plt.close()
        print(f"Created impact distribution chart: {impact_file}")
    
    # 4. Treatment comparison
    if 'treatment_comparison' in summary_tables and summary_tables['treatment_comparison'] is not None:
        treatment_comp = summary_tables['treatment_comparison']
        
        # Plot ratio comparisons
        ratio_cols = [col for col in treatment_comp.columns if col.endswith('_ratio')]
        
        if ratio_cols:
            plt.figure(figsize=(12, 8))
            for col in ratio_cols:
                comparison_name = col.replace('_ratio', '')
                plt.figure(figsize=(10, 6))
                
                # Get data and order by ratio
                plot_data = treatment_comp[['erg_name', col]].sort_values(col)
                
                # Plot horizontal bars
                bars = plt.barh(plot_data['erg_name'], plot_data[col])
                
                # Add labels
                for i, bar in enumerate(bars):
                    width = bar.get_width()
                    if width < 10:  # Only label reasonably sized bars
                        plt.text(width + 0.1, bar.get_y() + bar.get_height()/2,
                                f'{width:.2f}', ha='left', va='center')
                
                plt.title(f'Variant Ratio: {comparison_name}')
                plt.xlabel('Ratio')
                plt.ylabel('Gene')
                plt.axvline(x=1, color='red', linestyle='--', alpha=0.5)  # Mark ratio=1 line
                plt.tight_layout()
                
                ratio_file = os.path.join(viz_dir, f"{comparison_name}_ratio.png")
                plt.savefig(ratio_file, dpi=300)
                plt.close()
                print(f"Created {comparison_name} ratio chart: {ratio_file}")

def main():
    parser = argparse.ArgumentParser(description="Analyze variants in genes of interest")
    parser.add_argument("--annotated_dir", required=True, help="Directory containing annotated VCF files")
    parser.add_argument("--gene_mapping", required=True, help="File containing gene mapping information")
    parser.add_argument("--output_dir", required=True, help="Output directory for results")
    args = parser.parse_args()
    
    # Read gene mapping
    gene_mapping = read_gene_mapping(args.gene_mapping)
    
    # Extract and combine variants from all samples
    variants_df = combine_all_variants(args.annotated_dir, gene_mapping)
    
    if variants_df is not None:
        # Save the combined dataset
        os.makedirs(args.output_dir, exist_ok=True)
        variants_file = os.path.join(args.output_dir, "all_variants_in_genes.tsv")
        variants_df.to_csv(variants_file, sep='\t', index=False)
        print(f"Saved all variants to: {variants_file}")
        
        # Create summary tables
        summary_tables = create_summary_tables(variants_df, args.output_dir)
        
        # Create visualizations
        create_visualizations(summary_tables, args.output_dir)
        
        print("\nAnalysis complete!")
        print(f"Results available in: {args.output_dir}")
    else:
        print("Analysis failed: No variants found in genes of interest")

if __name__ == "__main__":
    main()